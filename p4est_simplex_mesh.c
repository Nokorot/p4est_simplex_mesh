#include "p4est_base.h"
#include "sc.h"
#include <sc_containers.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "sc_mpi.h"
#include "utils.c"

#ifdef VIM_LS
// #include <p4est_to_p8est.h>
#endif

#ifndef P4_TO_P8
#include "p4est_ghost.h"
#include <p4est_iterate.h>
#include "p4est_simplex_mesh.h"
#else
#include <p8est_iterate.h>
#include "p8est_simplex_mesh.h"
#include "p8est_connectivity.h"
#include "p8est_simplex_mesh.h"
#endif


// Local declarations


typedef struct {
  p4est_quadrant_t element; // We assume p.which_tree is set to treeid
  p4est_locidx_t quadid; 

  
  p4est_locidx_t volume_node;

  p4est_locidx_t face_nodes[ P4EST_FACES ];
#ifdef P4_TO_P8
  p4est_locidx_t edge_nodes[ P8EST_EDGES ];
#endif
  p4est_locidx_t corner_nodes[ P4EST_CHILDREN ];
}
element_nodes_t;

typedef struct
{
  sc_MPI_Comm mpicomm;
  int mpirank;

  p4est_t *p4est;
  p4est_ghost_t *ghost;
  p4est_geometry_t *geom;

  sc_array_t *tree_offsets;

  sc_array_t *owner;
  sc_hash_array_t *element_hash;
  sc_array_t *vertices;
}
p4est_simplex_nodes_data_t;


element_nodes_t *
p4est_find_or_insert_new_element_node(
    sc_hash_array_t *element_hash,
    p4est_quadrant_t quad, 
    p4est_locidx_t quadid,
    p4est_topidx_t treeid)
{
  size_t position;
  element_nodes_t *elem_node;
  quad.p.which_tree = treeid;

  elem_node = sc_hash_array_insert_unique(element_hash, &quad, &position);
  if (elem_node) {
    memset(elem_node, 0, sizeof(element_nodes_t));
    elem_node->element = quad;
    elem_node->quadid = quadid;
  } 
  else
    elem_node = sc_array_index(&element_hash->a, position);

  return elem_node;
}

// ITERATORS
void
iter_volume(p4est_iter_volume_info_t *info, void *user_data) {
  p4est_simplex_nodes_data_t *d = user_data;

  sc_hash_array_t *element_hash = d->element_hash;

  size_t position;
  element_nodes_t *elem_node;
  elem_node = p4est_find_or_insert_new_element_node(
      element_hash, *info->quad, info->quadid, info->treeid);

  // Add the quad centre node
  p4est_locidx_t cc = d->vertices->elem_count;
  elem_node->volume_node = cc;

  double abc[3];
  double *vx = sc_array_push(d->vertices);
  p4est_quad_center_coords(info->quad, abc);
  d->geom->X(d->geom, info->treeid, abc, vx);
}

void
iter_face(p4est_iter_face_info_t *info, void *user_data) {
  p4est_simplex_nodes_data_t *d = user_data;

  p4est_iter_face_side_t *side, *side_full, *side_hanging;
  p4est_quadrant_t *elem;
  element_nodes_t *elem_node;

  side_hanging = NULL;
  for (size_t zz=0; zz < info->sides.elem_count; ++zz) {
    side = sc_array_index(&info->sides, zz);

    if (side->is_hanging)
      side_hanging = side;
    else
      side_full = side;
  }

  if (!side_hanging)
    return;

  p4est_locidx_t cf;

  if (!side_full->is.full.is_ghost) {
    elem_node = p4est_find_or_insert_new_element_node(d->element_hash, 
        *side_full->is.full.quad, side_full->is.full.quadid, side_full->treeid);
  
    cf = elem_node->face_nodes[side_full->face];
  } else {
    cf = 0;
  }

  // cf != 0; Might happen if the node was added by a hanging edge
  if (!cf) { 
    // Add the face centre node
    cf = d->vertices->elem_count;

    double abc[3];
    double *vx = sc_array_push(d->vertices);

    p4est_quad_face_center_coords(side_full->is.full.quad, side_full->face, abc);
    d->geom->X(d->geom, side_full->treeid, abc, vx);
  }
 
  if (!side_full->is.full.is_ghost) {
    elem_node->face_nodes[side_full->face] = cf;
  }

  p4est_quadrant_t *quad;
  int face = side_hanging->face;

  for (size_t sz = 0; sz < P4EST_HALF; ++sz) {
    if (side_hanging->is.hanging.is_ghost[sz])
      continue;

    quad = side_hanging->is.hanging.quad[sz];
    p4est_locidx_t quadid = side_hanging->is.hanging.quadid[sz];

    elem_node = p4est_find_or_insert_new_element_node(d->element_hash, 
        *quad, quadid, side_hanging->treeid);

    int cid = p4est_quadrant_child_id(quad);
    int cfc = p4est_corner_face_corners[cid][face];
#ifndef P4_TO_P8
    int c = p4est_face_corners[face][cfc ^ 1];
#else
    int c = p4est_face_corners[face][cfc ^ 3];
#endif

    elem_node->corner_nodes[c] = cf;
  }
}

#ifdef P4_TO_P8
void
iter_edge(p8est_iter_edge_info_t *info, void *user_data) {
  p4est_simplex_nodes_data_t *d;

  p8est_iter_edge_side_t *side, *side_hanging, *side_full;
  element_nodes_t *elem_node;

  p4est_topidx_t treeid;
  p4est_locidx_t quadid;
  p4est_quadrant_t *quad;
  
  double abc[3], *vx;

  p4est_locidx_t ce;
  size_t zz, position;


  d = user_data;

  // Find a full and a hanging face
  side_full = NULL;
  side_hanging = NULL;
  for (zz = 0; zz < info->sides.elem_count; ++zz) {
    side = sc_array_index(&info->sides, zz);
    if (side->is_hanging) {
      side_hanging = side;

      if (side_full)
        break;

    } else {
      side_full = side;

      if (side_hanging)
        break;
    }
  }
  
  // If no side is hanging we will not add an edge node
  if (!side_hanging)
    return;


  // Add the edge node
  ce = d->vertices->elem_count;
  quad = side_full->is.full.quad;

  vx = sc_array_push(d->vertices);

  p8est_quad_edge_center_coords(quad, side_full->edge, abc);
  d->geom->X(d->geom, side_full->treeid, abc, vx);

  // We record the index in the element_nodes for each side
  for (size_t zz=0; zz < info->sides.elem_count; ++zz) {
    side = sc_array_index(&info->sides, zz);

    int edge = side->edge;

    if (!side->is_hanging) {
      if (side->is.full.is_ghost)
        continue;

      quad = side->is.full.quad;
      quadid = side->is.full.quadid;
      treeid = side->treeid;

      elem_node = p4est_find_or_insert_new_element_node(
          d->element_hash, *quad, quadid, side->treeid);
      
      elem_node->edge_nodes[side->edge] = ce;
  
      for (size_t efz = 0; efz < 2; ++efz) {
        int face = p8est_edge_faces[edge][efz];
      
        int cf;
        if (!elem_node->face_nodes[face]) {
          // Add face node
          cf = d->vertices->elem_count;

          vx = sc_array_push(d->vertices);

          p4est_quad_face_center_coords(quad, face, abc);
          d->geom->X(d->geom, treeid, abc, vx);

          elem_node->face_nodes[face] = cf;
        }

      }
    }
    else {

      for(size_t sz = 0; sz < 2; ++sz) {
        if (side->is.hanging.is_ghost[sz])
          continue;
    
        quad = side->is.hanging.quad[sz];
        quadid = side->is.hanging.quadid[sz];
  
        elem_node = p4est_find_or_insert_new_element_node(
            d->element_hash, *quad, quadid, side->treeid);

        int cid = p4est_quadrant_child_id(quad);
        int cfc = p8est_corner_edge_corners[cid][edge];
        int c = p8est_edge_corners[edge][cfc ^ 1];
  
        elem_node->corner_nodes[c] = ce;
      }
    }
  }
}
#endif

void
iter_corner(p4est_iter_corner_info_t *info, void *user_data) {
  p4est_simplex_nodes_data_t *d;

  p4est_iter_corner_side_t *side;

  sc_hash_array_t *element_hash;
  p4est_quadrant_t *elem;
  element_nodes_t *elem_node;
  size_t zz, position;

  d = user_data;

  zz = 0;
  side = sc_array_index(&info->sides, zz);

  // Add the quad corner node
  p4est_locidx_t cc = d->vertices->elem_count;

  double abc[3];
  double *vx = sc_array_push(d->vertices);

  p4est_quad_corner_coords(side->quad, side->corner, abc);
  d->geom->X(d->geom, side->treeid, abc, vx);

  // Record this corner_node for each of the adjacent elements
  element_hash = d->element_hash;
  for (zz = 0; zz < info->sides.elem_count; ++zz) {
    side = sc_array_index(&info->sides, zz);
    
    p4est_locidx_t quadid = side->quadid;

    elem_node = p4est_find_or_insert_new_element_node(
        d->element_hash, *side->quad, quadid, side->treeid);

    elem_node->corner_nodes[side->corner] = cc;
  }
}

void
p4est_new_simplex_mesh_nodes(
    p4est_t *p4est, 
    p4est_ghost_t *ghost,
    p4est_geometry_t *geom,

    p4est_simplex_nodes_data_t *d
    )
{
  // p4est_simplex_nodes_data_t d;
  
  size_t tz, t_off;
  p4est_tree_t *tree;
  
  d->p4est = p4est;
  d->ghost = ghost;
  d->geom = geom;
  d->mpicomm = sc_MPI_COMM_WORLD;
  d->mpirank = p4est->mpirank;

  t_off = 0;
  d->tree_offsets = sc_array_new_count(sizeof(size_t), p4est->trees->elem_count + 1);
  for (tz = 0; tz < p4est->trees->elem_count; ++tz) {
    *((size_t *) sc_array_index(d->tree_offsets, tz)) = t_off;
    tree = sc_array_index(p4est->trees, tz);
    t_off += tree->quadrants.elem_count;
  }

  d->vertices = sc_array_new(3*sizeof(double));


  // int vnodes = d->vnodes = P4EST_DIM_POW(3);
  p4est_locidx_t num_quads = p4est->local_num_quadrants;

  // This does not need to be a hash_array,
  //  it may simply be an array indexed by  element_id = tree_offsets[treeid] + quadid
  // d->element_hash = sc_array_new_count(sizeof(p4est_locidx_t), num_quads * vnodes);
  d->element_hash = sc_hash_array_new(
      sizeof(element_nodes_t),
      p4est_node_hash_piggy_fn,
      p4est_node_equal_piggy_fn,
      NULL);

  p4est_iterate(p4est, ghost, d,
      iter_volume, iter_face,
#ifdef P4_TO_P8
      iter_edge,
#endif
      iter_corner);
}

// Implementations

p4est_simplex_mesh_t *
p4est_new_simplex_mesh(
    p4est_t *p4est,
    p4est_geometry_t *geom,
    p4est_ghost_t *ghost)
{
  p4est_simplex_nodes_data_t d;

  // sc_hash_array_t *element_hash;
  sc_array_t *simplicies;
  
  size_t zz, fi, ei;
  p4est_locidx_t elemid;
  p4est_topidx_t t;
  p4est_tree_t *tree;
  sc_array_t *quadrants;
  p4est_quadrant_t qid, *quad;

  int face;
  int fc0, fc1; // fc2, fc3

  p4est_locidx_t cf, cc, c0, c1;

#ifdef P4_TO_P8
  int edge;
  // int face_has_hanging_edge;

  p4est_locidx_t ce, c2, c3;

  // p4est_locidx_t edge_noes[P8EST_EDGES];
#endif

  element_nodes_t *elem_node;
  p4est_locidx_t *simplex;
  double abc[3], *vx;
  size_t position;
  int found;

  
  p4est_new_simplex_mesh_nodes(p4est, ghost, geom, &d);
  // d.element_hash;

  // mesh = P4EST_ALLOC(p4est_simplex_mesh_t, 1);

  // mesh->simplicies = sc_array_new((P4EST_DIM + 1) *sizeof(p4est_locidx_t));
  simplicies = sc_array_new((P4EST_DIM + 1)*sizeof(p4est_locidx_t));

  // We iterate through the local quads to make the simplicies
  elemid = 0;
  for (t=p4est->first_local_tree;
      t<=p4est->last_local_tree; ++t) {

    tree = p4est_tree_array_index(p4est->trees, t);

    quadrants = &tree->quadrants;
    for (zz = 0; zz < quadrants->elem_count; ++zz, ++elemid) {
      quad = p4est_quadrant_array_index(quadrants, zz);

      qid = *quad;
      qid.p.which_tree = t;
    
      found = sc_hash_array_lookup(d.element_hash, &qid, &position);

      SC_CHECK_ABORT(found, "Element was not found in the element_hash");

      elem_node = sc_array_index(&d.element_hash->a, position);

      cc = elem_node->volume_node;

      for (fi = 0; fi < P4EST_FACES; ++fi) {
        cf = elem_node->face_nodes[fi];

        if (cf) // TODO:  cf = 0, should refer to the first node. Should do cf < 0.
                //        Would be fixed if we use a bitmap
        {
#ifndef P4_TO_P8
          c0 = elem_node->corner_nodes[ p4est_face_corners[fi][0] ];
          c1 = elem_node->corner_nodes[ p4est_face_corners[fi][1] ];

          simplex = sc_array_push(simplicies);
          simplex[0] = cc;
          simplex[1] = c0;
          simplex[2] = cf;

          simplex = sc_array_push(simplicies);
          simplex[0] = cc;
          simplex[1] = cf;
          simplex[2] = c1;
#else
          for (ei=0; ei<4; ++ei) {
            edge = p8est_face_edges[fi][ei];

            c0 = elem_node->corner_nodes[ p8est_edge_corners[edge][0] ];
            c1 = elem_node->corner_nodes[ p8est_edge_corners[edge][1] ];

            ce = elem_node->edge_nodes[edge];

            if (!ce)
            {
              simplex = sc_array_push(simplicies);
              simplex[0] = cc;
              simplex[1] = cf;
              simplex[2] = c0;
              simplex[3] = c1;
            }
            else
            {
              // edge center
              simplex = sc_array_push(simplicies);
              simplex[0] = cc;
              simplex[1] = cf;
              simplex[2] = c0;
              simplex[3] = ce;

              simplex = sc_array_push(simplicies);
              simplex[0] = cc;
              simplex[1] = cf;
              simplex[2] = ce;
              simplex[3] = c1;
            }
          }
#endif

        }
        else
        {  // Split as usual
          c0 = elem_node->corner_nodes[ p4est_face_corners[fi][0] ];
          c1 = elem_node->corner_nodes[ p4est_face_corners[fi][1] ];

#ifndef P4_TO_P8

          simplex = sc_array_push(simplicies);
          simplex[0] = cc;
          simplex[1] = c0;
          simplex[2] = c1;
#else
          c2 = elem_node->corner_nodes[ p4est_face_corners[fi][2] ];
          c3 = elem_node->corner_nodes[ p4est_face_corners[fi][3] ];

          int rot_corners = 0;
          if (quad->level == 0) {
            // The face split most be consistent on neighbouring faces.
            // This may be achieved by determining the diagonal, based on the global index
            // of the vertices, since neighbouring faces has the same vertices.
            // We do this by always splitting the diagonal with the vertex with the smallest
            // global index.
            //
            // NOTE: We observe that for level = 0, all the vertices are independent in p4est_nodes_t
            // so we may assume they have a global index

            
            // TODO: Determine this by global node nummer



#if 0
            p4est_locidx_t first = c0;

            int rank = first < nodes->offset_owned_indeps
                          ? nodes->nonlocal_ranks[first]
                          : p4est->mpirank;


            // THIS IS SO HACKY!
            int rank1;
#define cmp_node(NODE) \
            rank1 = NODE < nodes->offset_owned_indeps  \
                          ? nodes->nonlocal_ranks[NODE] \
                          : p4est->mpirank; \
            if (rank1 < rank || (rank1 == rank && NODE < first)) { \
              first = NODE; \
              rank = rank1; \
            }

            cmp_node(c1);
            cmp_node(c2);
            cmp_node(c3);
#undef cmp_node

            rot_corners = (first == c1 || first == c2);
#endif

            rot_corners = 0;

          } else {
            // We split the face, so that the parent face is split by an X:
            //    +---+---+
            //    | \ | / |
            //    +---+---+
            //    | / | \ |
            //    +---+---+
            // This way, the split is the same for all face-orientations (rotation, mirroring)
            // and thus the split will always be the same on any neighbouring face.
            //
            int dim = (fi >> 1);
            int c = p4est_quadrant_child_id(quad);
            int fchild =
                (dim == 0 ? c >> 1               : 0)
              | (dim == 1 ? c & 1 | (c >> 1) & 2 : 0)
              | (dim == 2 ? c & 3                : 0);

            rot_corners = (fchild == 1 || fchild == 2);
          }

          // The face is (by default) split along the c0-c3 diagonal:
          //   c2--c3
          //   |  / |
          //   | /  |
          //   c0--c1
          //
          // In order to split along the opposite diagonal, we rotate the
          // variables once anti-clockwise
          //
          //   c3--c1
          //   |  / |
          //   | /  |
          //   c2--c0
          //
          if (rot_corners) {
            p4est_locidx_t tmp;
            tmp = c1;
            c1 = c0;
            c0 = c2;
            c2 = c3;
            c3 = tmp;
          }

          simplex = sc_array_push(simplicies);
          simplex[0] = cc;
          simplex[1] = c0;
          simplex[2] = c1;
          simplex[3] = c3;

          simplex = sc_array_push(simplicies);
          simplex[0] = cc;
          simplex[1] = c3;
          simplex[2] = c2;
          simplex[3] = c0;

#endif
        } // if fc == 0
          
      }  // for (fi = 0..P4EST_FACES)
    } // for (zz = 0..|quadtrees|)
  } // for (t = first-..last-localtrees)

  sc_hash_array_destroy(d.element_hash);
  sc_array_destroy(d.tree_offsets);

  p4est_simplex_mesh_t *mesh;

  // TODO: Finish mesh, by makeing lnodes structure

  mesh = P4EST_ALLOC(p4est_simplex_mesh_t, 1);
  mesh->vertices = d.vertices;
  mesh->simplicies = simplicies;

  return mesh;
}

void
p4est_simplex_mesh_destroy(
    p4est_simplex_mesh_t *smesh)
{
  sc_array_destroy(smesh->vertices);
  sc_array_destroy(smesh->simplicies);
  P4EST_FREE(smesh);
}


void
p4est_quad_center_coords(
    p4est_quadrant_t *quad,
    double abc[3])
{
  double frac = 1.0 / (double) P4EST_ROOT_LEN;
  double half_len = P4EST_QUADRANT_LEN(quad->level + 1);

  abc[0] = (double) ( quad->x + half_len ) * frac;
  abc[1] = (double) ( quad->y + half_len ) * frac;
#ifdef P4_TO_P8
  abc[2] = (double) ( quad->z + half_len ) * frac;
#else
  abc[2] = 0;
#endif
}

void
p4est_quad_corner_coords(
    p4est_quadrant_t *quad,
    int c,
    double abc[3])
{
  double frac = 1.0 / (double) P4EST_ROOT_LEN;
  double quad_len = P4EST_QUADRANT_LEN(quad->level);

  abc[0] = frac * ( (double) quad->x
                    + ( (c & 1) ? quad_len : 0 ) );

  abc[1] = frac * ( (double) quad->y
                    + ( (c & 2) ? quad_len : 0 ) );

#ifdef P4_TO_P8
  abc[2] = frac * ( (double) quad->z
                    + ( (c & 4) ? quad_len : 0 ) );
#else
  abc[2] = 0;
#endif
}

void
p4est_quad_face_center_coords(
    p4est_quadrant_t *quad,
    int face,
    double abc[3])
{
  double frac = 1.0 / (double) P4EST_ROOT_LEN;
  double quad_len = P4EST_QUADRANT_LEN(quad->level);
  double half_len = P4EST_QUADRANT_LEN(quad->level + 1);

  abc[0] = frac * ( quad->x
           + ( !(face ^ 1) ? quad_len : 0)  // For face 1
           + (  (face >> 1) ? half_len : 0) // For all not x-orth faces
           );                             // (ie. fi != 0, 1)

  abc[1] = frac * ( quad->y
           + ( !(face ^ 3) ? quad_len : 0) // For face 3
           + ( !(face & 2) ? half_len : 0) // For all not y-orth faces
           );                            // (ie. fi != 2, 3)

#ifdef P4_TO_P8
  abc[2] = frac * ( quad->z
           + ( !(face ^ 5) ? quad_len : 0)  // For face 5
           + ( !(face >> 2) ? half_len : 0) // For all not z-orth faces
           );                             // (ie. fi != 4, 5)
#else
  abc[2] = 0;
#endif
}

#ifdef P4_TO_P8
void
p8est_quad_edge_center_coords(
    p4est_quadrant_t *quad,
    int edge,
    double abc[3])
{
  double frac = 1.0 / (double) P4EST_ROOT_LEN;
  double len = P4EST_QUADRANT_LEN(quad->level);
  double half_len = 0.5 * len;

  int axis = (edge >> 2) + (edge >> 4);
  int num = (edge >> 2*axis) & 3;

  abc[0] = frac * ( quad->x
      + ( axis == 0 ? half_len : 0 )  // x-parallel
      + ( (axis != 0 && edge & 1) ? len : 0 )
      );
  abc[1] = frac * ( quad->y
      + ( axis == 1 ? half_len : 0 )  // y-parallel
      + ( (axis == 0 && edge & 1) || (axis == 2 && edge & 2) ? len : 0 )
      );
  abc[2] = frac * ( quad->z
      + ( axis == 2 ? half_len : 0 )  // z-parallel
      + ( (axis != 2 && edge & 2) ? len : 0 )
      );
}
#endif


/*** Export as gmsh/vtk file **/

void p4est_simplex_mesh_write_gmsh_file(const char *filename, p4est_simplex_mesh_t *mesh) {
    FILE *file = fopen(filename, "w");
    if (!file) {
        perror("Error opening file");
        return;
    }

    // Write the header for Gmsh 2.2 format
    fprintf(file, "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n");

    // Write the nodes section
    size_t num_vertices = mesh->vertices->elem_count;
    fprintf(file, "$Nodes\n%zu\n", num_vertices);
    double *vertices = (double *) mesh->vertices->array;
    for (size_t i = 0; i < num_vertices; i++) {
        fprintf(file, "%zu %.16g %.16g %.16g\n", i + 1, vertices[3 * i], vertices[3 * i + 1], vertices[3 * i + 2]);
    }
    fprintf(file, "$EndNodes\n");

    // Write the elements section
    size_t num_simplicies = mesh->simplicies->elem_count;
    fprintf(file, "$Elements\n%zu\n", num_simplicies);
    p4est_locidx_t *simplicies = (p4est_locidx_t *) mesh->simplicies->array;
    for (size_t i = 0; i < num_simplicies; i++) {
#ifdef P4_TO_P8
        fprintf(file, "%zu 4 0 %u %u %u %u\n", i + 1, simplicies[4 * i] + 1, simplicies[4 * i + 1] + 1, simplicies[4 * i + 2] + 1, simplicies[4 * i + 3] + 1);
#else
        fprintf(file, "%zu 2 0 %u %u %u\n", i + 1, simplicies[3 * i] + 1, simplicies[3 * i + 1] + 1, simplicies[3 * i + 2] + 1);
#endif
    }
    fprintf(file, "$EndElements\n");

    fclose(file);
}

void p4est_simplex_mesh_write_vtk_file(const char *filename, p4est_simplex_mesh_t *mesh) {
    FILE *file = fopen(filename, "w");
    if (!file) {
        perror("Error opening file");
        return;
    }

    // Write the header for VTK file
    fprintf(file, "# vtk DataFile Version 3.0\n");
    fprintf(file, "Simplex Mesh\n");
    fprintf(file, "ASCII\n");
    fprintf(file, "DATASET UNSTRUCTURED_GRID\n");

    // Write the points section
    size_t num_vertices = mesh->vertices->elem_count;
    fprintf(file, "POINTS %zu double\n", num_vertices);
    double *vertices = (double *) mesh->vertices->array;
    for (size_t i = 0; i < num_vertices; i++) {
        fprintf(file, "%.16g %.16g %.16g\n", vertices[3 * i], vertices[3 * i + 1], vertices[3 * i + 2]);
    }

    // Write the cells section
    size_t num_simplicies = mesh->simplicies->elem_count;
#ifdef P4_TO_P8
    fprintf(file, "CELLS %zu %zu\n", num_simplicies, num_simplicies * 5);
#else
    fprintf(file, "CELLS %zu %zu\n", num_simplicies, num_simplicies * 4);
#endif
    p4est_locidx_t *simplicies = (p4est_locidx_t *) mesh->simplicies->array;
    for (size_t i = 0; i < num_simplicies; i++) {
#ifdef P4_TO_P8
        fprintf(file, "4 %u %u %u %u\n", simplicies[4 * i], simplicies[4 * i + 1], simplicies[4 * i + 2], simplicies[4 * i + 3]);
#else
        fprintf(file, "3 %u %u %u\n", simplicies[3 * i], simplicies[3 * i + 1], simplicies[3 * i + 2]);
#endif
    }

    // Write the cell types section
    fprintf(file, "CELL_TYPES %zu\n", num_simplicies);
    for (size_t i = 0; i < num_simplicies; i++) {
#ifdef P4_TO_P8
        fprintf(file, "10\n");  // VTK_TETRA type
#else
        fprintf(file, "5\n");   // VTK_TRIANGLE type
#endif
    }

    fclose(file);
}

