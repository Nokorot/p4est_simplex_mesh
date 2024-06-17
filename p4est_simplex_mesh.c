
#include <sc.h>
#include <sc_containers.h>
#include <sc_mpi.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "utils.c"

#ifdef VIM_LS
// #include <p4est_to_p8est.h>
#endif

#include <p4est_base.h>

#ifndef P4_TO_P8
#include "p4est_lnodes.h"
#include "p4est_ghost.h"
#include <p4est_iterate.h>
#include "p4est_simplex_mesh.h"
#include "p4est_communication.h"
#include "p4est_geometry.h"
#else
#include "p8est_lnodes.h"
#include <p8est_iterate.h>
#include "p8est_simplex_mesh.h"
#include "p8est_simplex_mesh.h"
#include "p8est_communication.h"
#include "p8est_geometry.h"
#endif


// Local declarations


typedef struct {
  p4est_locidx_t volume_node;

  p4est_locidx_t face_nodes[ P4EST_FACES ];
#ifdef P4_TO_P8
  p4est_locidx_t edge_nodes[ P8EST_EDGES ];
#endif
  p4est_locidx_t corner_nodes[ P4EST_CHILDREN ];
}
element_node_t;

typedef struct
{
  sc_MPI_Comm mpicomm;
  int mpirank;

  p4est_t *p4est;
  p4est_ghost_t *ghost;
  p4est_geometry_t *geom;

  sc_array_t *tree_offsets;

  sc_array_t *owned_count;

  sc_array_t *owners;
  sc_array_t *element_nodes;

  sc_array_t *vertices;
}
p4est_simplex_nodes_data_t;


element_node_t *
p4est_element_node_get(
    p4est_simplex_nodes_data_t *d,
    p4est_topidx_t treeid,
    p4est_locidx_t quadid
    )
{
  p4est_locidx_t *toff, elemid;

  toff = sc_array_index(d->tree_offsets, treeid);
  elemid = *toff + quadid;

  return sc_array_index(d->element_nodes, elemid);
}

p4est_locidx_t
p4est_simplex_push_node(
    p4est_simplex_nodes_data_t *d,
    p4est_topidx_t treeid,
    double abc[3],
    int owner
    )
{
  p4est_locidx_t node_id;
  double *vx;
  int *o;
  size_t *ocount;

  node_id = d->vertices->elem_count;

  vx = sc_array_push(d->vertices);
  d->geom->X(d->geom, treeid, abc, vx);

  o = sc_array_push(d->owners);
  *o = owner;

  ocount = sc_array_index(d->owned_count, owner);
  (*ocount)++;

  return node_id;
}

int
p4est_find_ghost_owner(
    p4est_t *p4est,
    p4est_ghost_t *ghost,
    p4est_locidx_t ghostid)
{
  p4est_quadrant_t *ghost_quad;
  p4est_topidx_t treeid;


  ghost_quad = sc_array_index(&ghost->ghosts, ghostid);
  treeid = ghost_quad->p.piggy3.which_tree;

  return p4est_comm_find_owner(
            p4est,
            treeid,
            ghost_quad,
            p4est->mpirank);
}

// ITERATORS
void
iter_volume(p4est_iter_volume_info_t *info, void *user_data)
{
  p4est_simplex_nodes_data_t *d;
  element_node_t *elem_node;

  double abc[3];
  int owner;
  p4est_locidx_t cc;

  d = user_data;
  elem_node = p4est_element_node_get(d, info->treeid, info->quadid);

  // Add the quad centre node
  p4est_quad_center_coords(info->quad, abc);

  owner = d->mpirank;

  cc = p4est_simplex_push_node(d, info->treeid, abc, owner);
  elem_node->volume_node = cc;
}

void
iter_face(p4est_iter_face_info_t *info, void *user_data)
{
  p4est_simplex_nodes_data_t *d;

  p4est_iter_face_side_t *side, *side_full, *side_hanging;
  p4est_quadrant_t *quad;
  element_node_t *elem_node;

  double abc[3];
  int owner;

  p4est_locidx_t cf;
  size_t sz, zz;
  int face;

  d = user_data;

  side_hanging = NULL;
  for (zz = 0; zz < info->sides.elem_count; ++zz) {
    side = sc_array_index(&info->sides, zz);

    if (side->is_hanging)
      side_hanging = side;
    else
      side_full = side;
  }

  if (!side_hanging)
    return;

  if (!side_full->is.full.is_ghost) {
    elem_node = p4est_element_node_get(d,
        side_full->treeid, side_full->is.full.quadid);

    cf = elem_node->face_nodes[side_full->face];
  } else {
    cf = -1;
  }

  // cf != 0; Might happen if the node was added by a hanging edge
  if (cf < 0) {
    // Add the face centre node
    p4est_quad_face_center_coords(side_full->is.full.quad,
        side_full->face, abc);

    // The full_side is always the owner
    if (side_full->is.full.is_ghost) {
      int ghostid = side_full->is.full.quadid;
      owner = p4est_find_ghost_owner(d->p4est, d->ghost, ghostid);
    } else {
      owner = d->mpirank;
    }

    cf = p4est_simplex_push_node(d, side_full->treeid, abc, owner);
  }

  if (!side_full->is.full.is_ghost) {
    elem_node->face_nodes[side_full->face] = cf;
  }

  face = side_hanging->face;

  for (sz = 0; sz < P4EST_HALF; ++sz) {
    if (side_hanging->is.hanging.is_ghost[sz])
      continue;

    quad = side_hanging->is.hanging.quad[sz];
    p4est_locidx_t quadid = side_hanging->is.hanging.quadid[sz];

    elem_node = p4est_element_node_get(d, side_hanging->treeid, quadid);

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
iter_edge(p8est_iter_edge_info_t *info, void *user_data)
{
  p4est_simplex_nodes_data_t *d;

  p8est_iter_edge_side_t *side, *side_hanging, *side_full;
  element_node_t *elem_node;

  p4est_topidx_t treeid;
  p4est_locidx_t quadid;
  p4est_quadrant_t *quad;

  double abc[3];
  int owner;

  p4est_locidx_t ce, cf;
  size_t sz, zz, efz;
  int face, edge;

  d = user_data;

  // Find a full and a hanging face
  side_full = NULL;
  side_hanging = NULL;
  for (zz = 0; zz < info->sides.elem_count; ++zz) {
    side = sc_array_index(&info->sides, zz);
    if (side->is_hanging)
      side_hanging = side;
    else {
      side_full = side;

      // Only a full side may be the owner
      owner = !(side->is.full.is_ghost)
                  ? d->mpirank
                  : p4est_find_ghost_owner(
                      d->p4est, d->ghost, side->is.full.quadid);
    }

    // If only need to find one full and one hanging side
    // if (side_hanging && side_full)
    //   break;
  }

  // If no side is hanging we will not add an edge node
  if (!side_hanging) {
    return;
  }

  SC_CHECK_ABORT(side_full, "Got edge with only hanging sides");

  // Add the edge node
  quad = side_full->is.full.quad;
  p8est_quad_edge_center_coords(quad, side_full->edge, abc);
  ce = p4est_simplex_push_node(d, side_full->treeid, abc, owner);

  // We record the index in the element_nodes for each side
  for (zz = 0; zz < info->sides.elem_count; ++zz) {
    side = sc_array_index(&info->sides, zz);
    edge = side->edge;

    if (!side->is_hanging) {
      if (side->is.full.is_ghost)
        continue;

      quad = side->is.full.quad;
      quadid = side->is.full.quadid;
      treeid = side->treeid;

      elem_node = p4est_element_node_get(d, side->treeid, quadid);
      elem_node->edge_nodes[side->edge] = ce;

      for (efz = 0; efz < 2; ++efz) {
        face = p8est_edge_faces[edge][efz];

        if (elem_node->face_nodes[face] < 0) {
          // Add face node
          p4est_quad_face_center_coords(quad, face, abc);

          cf = p4est_simplex_push_node(d, treeid, abc, owner);
          elem_node->face_nodes[face] = cf;
        }
      }
    }
    else {
      for(sz = 0; sz < 2; ++sz) {
        if (side->is.hanging.is_ghost[sz])
          continue;

        quad = side->is.hanging.quad[sz];
        quadid = side->is.hanging.quadid[sz];

        int cid = p4est_quadrant_child_id(quad);
        int cfc = p8est_corner_edge_corners[cid][edge];
        int c = p8est_edge_corners[edge][cfc ^ 1];

        elem_node = p4est_element_node_get(d, side->treeid, quadid);
        elem_node->corner_nodes[c] = ce;
      }
    }
  }
}
#endif

void
iter_corner(p4est_iter_corner_info_t *info, void *user_data)
{
  p4est_simplex_nodes_data_t *d;

  p4est_iter_corner_side_t *side;

  element_node_t *elem_node;
  p4est_locidx_t cc;
  size_t zz;
  int owner, ghost_owner;

  d = user_data;

  owner = d->mpirank;
  for (zz = 0; zz < info->sides.elem_count; ++zz) {
    side = sc_array_index(&info->sides, zz);

    if (!side->is_ghost)
      continue;

    ghost_owner = p4est_find_ghost_owner(d->p4est, d->ghost, side->quadid);

    if (ghost_owner < owner)
      owner = ghost_owner;
  }

  cc = -1;
  // Record this corner_node for each of the adjacent elements
  for (zz = 0; zz < info->sides.elem_count; ++zz) {
    side = sc_array_index(&info->sides, zz);

    if (side->is_ghost)
      continue;

    if (cc < 0) {
      // Add the quad corner node
      double abc[3];
      p4est_quad_corner_coords(side->quad, side->corner, abc);
      cc = p4est_simplex_push_node(d, side->treeid, abc, owner);
    }

    p4est_locidx_t quadid = side->quadid;

    elem_node = p4est_element_node_get(d, side->treeid, quadid);
    elem_node->corner_nodes[side->corner] = cc;
  }
}

void
p4est_new_simplex_mesh_nodes(
    p4est_simplex_nodes_data_t *d
    )
{
  p4est_t *p4est;
  p4est_ghost_t *ghost;
  p4est_tree_t *tree;

  p4est_locidx_t num_local_elements;
  size_t tz, t_off;

  p4est = d->p4est;
  ghost = d->ghost;

  d->mpicomm = sc_MPI_COMM_WORLD;
  d->mpirank = p4est->mpirank;

  t_off = 0;
  d->tree_offsets = sc_array_new_count(sizeof(size_t), p4est->trees->elem_count + 1);
  for (tz = 0; tz < p4est->trees->elem_count; ++tz) {
    *((size_t *) sc_array_index(d->tree_offsets, tz)) = t_off;
    tree = sc_array_index(p4est->trees, tz);
    t_off += tree->quadrants.elem_count;
  }

  d->vertices     = sc_array_new(3*sizeof(double));
  d->owners       = sc_array_new(sizeof(int));
  d->owned_count  = sc_array_new_count(sizeof(size_t), p4est->mpisize);

  for (size_t zz=0; zz < p4est->mpisize; ++zz) {
    size_t *owned = sc_array_index(d->owned_count, zz);
    *owned = 0;
  }

  num_local_elements = p4est->local_num_quadrants;
  d->element_nodes = sc_array_new_count(sizeof(element_node_t),
                                        num_local_elements);

  // Initialize all nodes with -1
  memset(d->element_nodes->array, 0xff, num_local_elements*sizeof(element_node_t));
  // memset(d->element_nodes->array, 0, num_quads*sizeof(element_nodes_t));

  p4est_iterate(p4est, ghost, d,
      iter_volume, iter_face,
#ifdef P4_TO_P8
      iter_edge,
#endif
     iter_corner);

}

void
p4est_simplex_mesh_reoder_nodes(
    p4est_simplex_nodes_data_t *d)
{
  p4est_t *p4est;

  int mpisize, mpirank;

  p4est_locidx_t num_local_elements,  num_local_nodes;

  p4est_locidx_t *proc_offsets, offset;
  size_t *owned;
  int *owner;
  size_t zz, elem, cc_off;
  int vnodes;

  p4est_locidx_t *new_local_nodes, *proc_node_count, new_cc;
  p4est_locidx_t *elem_nodes, *new_elem_nodes, cc;
  sc_array_t *new_vertices;
  sc_array_t *new_element_nodes;

  p4est = d->p4est;

  mpisize = p4est->mpisize;
  mpirank = p4est->mpirank;

  num_local_elements = p4est->local_num_quadrants;

  // We begin by computing where to put the nodes owned
  // by the respective processes

  proc_offsets = P4EST_ALLOC(p4est_locidx_t, mpisize + 1);

  owned = sc_array_index(d->owned_count, mpirank);
  offset = *owned;

  for (int zz=0; zz < mpisize; ++zz) {
    if (zz == mpirank) { // We put all locally owned nodes first
      proc_offsets[zz] = 0;
      continue;
    }

    proc_offsets[zz] = offset;

    owned = sc_array_index(d->owned_count, zz);
    offset += *owned;
  }

  num_local_nodes = proc_offsets[mpisize] = offset;

  // Now we move the vertices into the new order and compute
  // the new location for each local node

  proc_node_count = P4EST_ALLOC_ZERO(p4est_locidx_t, mpisize);
  new_local_nodes = P4EST_ALLOC(p4est_locidx_t, num_local_nodes);
  new_vertices = sc_array_new_count(3*sizeof(double), num_local_nodes);

  for (zz = 0; zz < num_local_nodes; ++zz) {
    owner = sc_array_index(d->owners, zz);

    new_cc = proc_offsets[*owner] + (proc_node_count[*owner]++);
    new_local_nodes[zz] = new_cc;

    // printf("[%d] A %ld -> %d\n", mpirank, zz, new_cc);
    memcpy(sc_array_index(new_vertices, new_cc),
           sc_array_index( d->vertices, zz),
           3*sizeof(double));
  }

  P4EST_FREE(proc_offsets);
  P4EST_FREE(proc_node_count);

  sc_array_destroy(d->vertices);
  sc_array_destroy(d->owners);
  d->vertices = new_vertices;

  // Finally we remap the element_nodes to refer to the new node indices

  vnodes = P4EST_DIM_POW(3); // = sizeof(element_node_t) / sizeof(p4est_locidx_t);

  new_element_nodes = sc_array_new_count(sizeof(element_node_t),
                                         num_local_elements);

  for (elem = 0; elem < num_local_elements; ++elem) {
    elem_nodes = sc_array_index(d->element_nodes, elem);
    new_elem_nodes = sc_array_index(new_element_nodes, elem);

    for (cc_off = 0; cc_off < vnodes; ++cc_off) {
      cc = elem_nodes[cc_off];
      new_elem_nodes[cc_off] = (cc < 0) ? -1 : new_local_nodes[cc];
    }
  }

  P4EST_FREE(new_local_nodes);

  sc_array_destroy(d->element_nodes);
  d->element_nodes = new_element_nodes;

  // p4est_lnodes_t lnode;
  // lnode.owned_count = local_count;
  // lnode.num_local_nodes = num_local_nodes;

  // lnode.element_nodes = new_elem_nodes->array;

  // lnode.global_offset: Needs communication
  // lnode.global_owned_count: Needs communication

  // lnode.sharers : Might need to be constructed in the iterator

  sc_array_destroy(d->owned_count);
}

// Implementations

p4est_simplex_mesh_t *
p4est_new_simplex_mesh(
    p4est_t *p4est,
    p4est_geometry_t *geom,
    p4est_ghost_t *ghost)
{
  p4est_simplex_nodes_data_t d;
  sc_array_t *simplicies;

  size_t zz, fi;
  p4est_locidx_t elemid;
  p4est_topidx_t t;
  p4est_tree_t *tree;
  sc_array_t *quadrants;
  p4est_quadrant_t qid, *quad;

  p4est_locidx_t cf, cc, c0, c1;

#ifdef P4_TO_P8
  size_t ei;
  int edge;
  // int face_has_hanging_edge;

  p4est_locidx_t ce, c2, c3;
#endif

  element_node_t *elem_node;
  p4est_locidx_t *simplex;

  d.p4est = p4est;
  d.ghost = ghost;
  d.geom = geom;

  p4est_new_simplex_mesh_nodes(&d);
  p4est_simplex_mesh_reoder_nodes(&d);

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

      elem_node = sc_array_index(d.element_nodes, elemid);
      cc = elem_node->volume_node;

      for (fi = 0; fi < P4EST_FACES; ++fi) {
        cf = elem_node->face_nodes[fi];

        if (0 <= cf)
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
          for (ei = 0; ei < 4; ++ei) {
            edge = p8est_face_edges[fi][ei];

            c0 = elem_node->corner_nodes[ p8est_edge_corners[edge][0] ];
            c1 = elem_node->corner_nodes[ p8est_edge_corners[edge][1] ];

            ce = elem_node->edge_nodes[edge];

            if (ce < 0)
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

            // TODO: Determine this by global node number
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

  sc_array_destroy(d.element_nodes);
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

