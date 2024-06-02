// #include <p4est_to_p8est.h>

#include "p4est_base.h"
#include "sc.h"
#include "sc_containers.h"
#include <sc_options.h>

#ifndef P4_TO_P8
#include "p4est_search.h"
#include "p4est_nodes.h"
#include <p4est_extended.h>
#include "p4est_vtk.h"
#else
#include "p8est_search.h"
#include "p8est_nodes.h"
#include <p8est_extended.h>
#include <p8est_vtk.h>

#include "simplex3.h" // This is only for my vim-language-server
#endif

#include "simplex2.h"

#include "p4est_simplex_to_gmsh.h"

p4est_vtk_context_t *
p4est_vtk_write_header_points(
    p4est_vtk_context_t *cont,
    sc_array_t *vertices);

#ifdef P4_TO_P8
#define QUAD_COORDS(quad) \
            (quad).x / (double) P4EST_ROOT_LEN, \
            (quad).y / (double) P4EST_ROOT_LEN, \
            (quad).z / (double) P4EST_ROOT_LEN

#else
#define QUAD_COORDS(quad) \
            (quad).x / (double) P4EST_ROOT_LEN, \
            (quad).y / (double) P4EST_ROOT_LEN
#endif

void 
p4est_face_center_coords(
    p4est_quadrant_t *quad,
    int face,
    double abc[3]
    )
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

unsigned int
p4est_quad_face_equal_fn(
    const p4est_quad_face_node_t *v1,
    const p4est_quad_face_node_t *v2,
    const void *used_data)
{
  return p4est_quadrant_equal_fn(&v1->quad, &v2->quad, NULL) \
    && v1->which_tree == v2->which_tree
    && v1->face == v2->face;
}

unsigned int
p4est_quad_face_hash_fn(const p4est_quad_face_node_t *v, const void *user_data)
{
  unsigned int qhash = p4est_quadrant_hash_fn(&v->quad, user_data);
  int b, c;
  b = v->which_tree;
  c = v->face;

  sc_hash_mix(qhash, b, c);
  sc_hash_final(qhash, b, c);

  return (unsigned) c;
}


/** Refines the p4est, producing hanging-nodes
 */
static int refine_fn (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t *quadrant)
{
  context_t *g = (context_t *) p4est->user_pointer;
  options_t *opts = g->opts;

  int refine_level = g->opts->maxlevel;
  if ((int) quadrant->level >= (refine_level - (int) (which_tree % 3))) {
    return 1;
  }
  if (quadrant->level == 1 && p4est_quadrant_child_id (quadrant) == 3) {
    // return 0;
    return 1;
  }
  if (quadrant->x == P4EST_LAST_OFFSET (2) &&
      quadrant->y == P4EST_LAST_OFFSET (2)) {
    return 1;
  }
#ifndef P4_TO_P8
  if (quadrant->x >= P4EST_QUADRANT_LEN (2)) {
    return 0;
  }
#else
  if (quadrant->z >= P8EST_QUADRANT_LEN (2)) {
    return 0;
  }
#endif

  return 1;
}

int
create_context_from_opts(context_t *g, options_t *opts)
{
  p4est_t *p4est;
  p4est_connectivity_t *conn;
  p4est_ghost_t *ghost_layer;
  p4est_geometry_t *geom = NULL;

  /* Create p4est */
#ifndef P4_TO_P8
  if (strcmp(opts->conn, "sphere2d") == 0) {
    conn = p4est_connectivity_new_cubed();
  }
  else
#endif
  {
    conn = p4est_connectivity_new_byname (opts->conn);
  }

  if (!conn) {
      printf("ERROR: HEY\n");
  }

  g->conn = conn;
  p4est = p4est_new_ext (g->mpicomm, conn, 0, opts->minlevel, 1, 0, NULL, NULL);
  g->p4est = p4est;
  p4est->user_pointer = (void *) g;
  p4est_refine_ext (p4est, 1, opts->maxlevel, refine_fn, NULL, NULL);
  p4est_partition (p4est, 0, NULL);
  p4est_balance (p4est, P4EST_CONNECT_FULL, NULL);

  /* Create ghost layer */
  ghost_layer = p4est_ghost_new(p4est, P4EST_CONNECT_FULL);
  g->ghost_layer = ghost_layer;

  /* Create geometry */
#ifndef P4_TO_P8
  if (opts->conn && !strcmp(opts->conn, "disk2d")) {
    geom = p4est_geometry_new_disk2d (conn, 1.0, 2.0);
  } else if (opts->conn && !strcmp(opts->conn, "icosahedron")) {
    geom = p4est_geometry_new_icosahedron (conn, 1.0);
  } else if (opts->conn && !strcmp(opts->conn, "shell2d")) {
    geom = p4est_geometry_new_shell2d (conn, 2.0, 1.0);
  } else if (opts->conn && !strcmp(opts->conn, "sphere2d")) {
    geom = p4est_geometry_new_sphere2d(conn, 1.0);
  }
#else
  if (opts->conn && !strcmp(opts->conn, "shell")) {
    geom = p8est_geometry_new_shell (conn, 2.0, 1.0);
  } else if (opts->conn && !strcmp(opts->conn, "sphere")) {
    geom = p8est_geometry_new_sphere (conn, 3.0, 2.0, 1.0);
  } else if (opts->conn && !strcmp(opts->conn, "torus")) {
    geom = p8est_geometry_new_torus (conn, 1.0, 2.0, 4.0);
  }
#endif
  else {
    geom = p4est_geometry_new_connectivity(conn);
  }
  g->geom = geom;

  return 0;
}

void
destroy_context(context_t *g)
{
  if (g->p4est)
    p4est_destroy(g->p4est);
  if (g->conn)
    p4est_connectivity_destroy(g->conn);
  if (g->ghost_layer)
    p4est_ghost_destroy(g->ghost_layer);
  if (g->geom)
    p4est_geometry_destroy(g->geom);
}

void
p4est_simplex_mesh_destroy(
    p4est_simplex_mesh_t *smesh)
{
  sc_array_destroy(smesh->vertices);
  sc_array_destroy(smesh->simplicies);
  P4EST_FREE(smesh);
}

typedef struct {
  context_t *g;
  p4est_nodes_t *nodes;

  sc_hash_array_t *hash_record;
//   sc_hash_array_t *hanging_face_nodes;
// #ifdef P4_TO_P8
//   sc_hash_array_t *hanging_edge_nodes;
// #endif

  sc_array_t *vertices;

  size_t num_indep;
  size_t num_fnodes;

} p4est_record_hanging_data_t;


void p4est_record_per_quad(
    p4est_record_hanging_data_t *d,
    p4est_topidx_t t,
    p4est_quadrant_t *quad,
    size_t elem,
    int is_ghost)
{
  p4est_connectivity_t *conn;
  p4est_ghost_t *ghost;

  sc_hash_array_t *hash_record;

  size_t num_indep, num_fnodes;

  int fi, face, nface;
  int c, fc, opp_corner0;

#ifdef P4_TO_P8
  int edge;
  int ec, opp_ec;

  p4est_locidx_t ce;
#endif

  double *vx, abc[3];

  p4est_topidx_t rt;
  p4est_quadrant_t parent, n;
  p4est_locidx_t c0;

#ifndef P4_TO_P8
  p4est_quad_face_node_t *added, qid;
  hash_record = d->hash_record;
#else
  p8est_quad_edge_node_t *added, qid;

  sc_array_t *quads, *treeids, *nedges;
#endif

  conn = d->g->conn;
  ghost = d->g->ghost_layer;

  hash_record = d->hash_record;

  num_indep = d->num_indep;
  num_fnodes = d->num_fnodes;

  // if (quad->level == 0) { // This case need to be considered differently
  // TODO: Need to make sure faces are split the same from both sides
  //  Only relevant in 3D
  // }

  c = p4est_quadrant_child_id(quad);

  // Suffices to do this for child == 0, 7
  // if (c != 0 && c != P4EST_CHILDREN - 1)
  //   return;


  // Loop through the faces that are on the boundary of the parent,
  // and thus may be hanging
  for (fi = 0; fi < P4EST_DIM; ++fi) {
    face = p4est_corner_faces[c][fi];
    p4est_quadrant_parent(quad, &parent);

    if (is_ghost) {
        rt = p4est_quadrant_face_neighbor_extra(
                  &parent, t, face, &n, &nface, conn);
        if (rt < 0) // neighbour is outside conn
          continue;

        n.p.piggy3.which_tree = rt;

        ssize_t lnid; // p4est_quadrant_compare_piggy
        lnid = sc_array_bsearch (&(ghost->mirrors),
                                 &n,
                                 p4est_quadrant_compare_piggy);
        if (lnid < 0)  // neighbour \r n is not a local quad.
          continue;
    }
    else {
      fc = p4est_corner_face_corners[c][face];
      opp_corner0 = p4est_face_corners[face][(P4EST_HALF - 1) ^ fc];

      c0 = d->nodes->local_nodes[P4EST_CHILDREN * elem + opp_corner0];

      if (c0 < num_indep || c0 >= num_indep + num_fnodes)
        continue;

      rt = p4est_quadrant_face_neighbor_extra(
                &parent, t, face, &n, &nface, conn);

      if (rt < 0)
        continue;
    }

    // Hanging face
#ifndef P4_TO_P8
    qid.which_tree = rt;
    qid.quad = n;
    qid.face = nface % P4EST_FACES;

    added = sc_hash_array_insert_unique(hash_record, &qid, NULL);
    if (added) {

      *added = qid;
      if (is_ghost) {
        vx = sc_array_push(d->vertices);

        p4est_face_center_coords(&qid.quad, qid.face, abc);
        d->g->geom->X(d->g->geom, qid.which_tree, abc, vx);

        added->node_index = -d->vertices->elem_count;
      }
      else {
        added->node_index = c0;
      }
    }
#else
    for (size_t ei = 0; ei < 3; ++ei) {
      int edge = p8est_corner_edges[c][ei];
      int efe = p8est_edge_face_edges[edge][face];

      if (efe < 0)
        continue;

      int nedge;

      if (t == rt) {
        nedge = p8est_face_edges[nface][efe];
      } else {
        int nf = nface % P4EST_FACES, o = nface / P4EST_FACES;

        int nfedge =
            p8est_connectivity_face_neighbor_face_edge(efe, face, nf, o);

        nedge = p8est_face_edges[nf][nfedge];
      }

      qid.quad = n;
      qid.which_tree = rt;
      qid.edge = nedge;

      added = sc_hash_array_insert_unique(hash_record, &qid, NULL);
      if (added) {
        *added = qid;
  
        if (is_ghost) {
          vx = sc_array_push(d->vertices);
          ce = - d->vertices->elem_count;

          p8est_edge_center_coords(&qid.quad, qid.edge, abc);
          d->g->geom->X(d->g->geom, qid.which_tree, abc, vx);
        } else {
          ec = p8est_corner_edge_corners[c][edge];
          opp_ec = p8est_edge_corners[edge][1 ^ ec];

          ce = d->nodes->local_nodes[P4EST_CHILDREN * elem + opp_ec];
        }

        added->node_index = ce;
      }
    }
#endif
  }

#ifdef P4_TO_P8
  // Record hanging edges

  quads = sc_array_new(sizeof(p4est_quadrant_t));
  treeids = sc_array_new(sizeof(p4est_topidx_t));
  nedges = sc_array_new(sizeof(int));

  for (fi = 0; fi < P4EST_DIM; ++fi) {
    edge = p8est_corner_edges[c][fi];
    ec = p8est_edge_corners[edge][0];

    p4est_locidx_t c0 = d->nodes->local_nodes[P4EST_CHILDREN * elem + ec];

    if (c0 >= d->num_indep + d->num_fnodes) // Hanging edge
    {
      sc_array_reset(quads);
      sc_array_reset(treeids);
      sc_array_reset(nedges);

      p4est_quadrant_parent(quad, &parent);
      p8est_quadrant_edge_neighbor_extra(&parent, t, edge, quads, treeids,
                                         nedges, d->g->conn);

      for (int jj = 0; jj < quads->elem_count; ++jj) {
        qid.quad = *((p8est_quadrant_t *)sc_array_index(quads, jj));
        qid.which_tree = *((p4est_locidx_t *)sc_array_index(treeids, jj));
        qid.edge = *((int *)sc_array_index(nedges, jj));

        added = sc_hash_array_insert_unique(hash_record, &qid, NULL);
        if (added) {
          *added = qid;
          added->node_index = c0;
        }
      }
    }
  }

  sc_array_destroy(quads);
  sc_array_destroy(treeids);
  sc_array_destroy(nedges);
#endif
}

void
p4est_record_hanging_faces(
    context_t *g,
    sc_array_t *vertices,
    p4est_nodes_t *nodes,
    sc_hash_array_t **hash_record
    ) {
  p4est_record_hanging_data_t d;

  p4est_topidx_t t, rt;
  p4est_tree_t *tree;
  p4est_quadrant_t *quad;

  sc_array_t *quadrants;
  size_t zz;

  d.g = g;
  d.vertices = vertices;
  d.nodes = nodes;
  d.num_indep = nodes->indep_nodes.elem_count;
  d.num_fnodes = nodes->face_hangings.elem_count;

#ifndef P4_TO_P8
  d.hash_record = sc_hash_array_new(sizeof(p4est_quad_face_node_t),
                                        (sc_hash_function_t) p4est_quad_face_hash_fn,
                                        (sc_equal_function_t) p4est_quad_face_equal_fn,
                                        NULL);
#else
  d.hash_record = sc_hash_array_new(sizeof(p8est_quad_edge_node_t),
                                        (sc_hash_function_t)  p4est_quad_edge_hash_fn,
                                        (sc_equal_function_t) p4est_quad_edge_equal_fn,
                                        NULL);
#endif
  *hash_record = d.hash_record;

  p4est_locidx_t elem = 0;
  for (t=g->p4est->first_local_tree;
      t<=g->p4est->last_local_tree; ++t) {

    tree = p4est_tree_array_index(g->p4est->trees, t);
    quadrants = &tree->quadrants;

    for (zz = 0; zz < quadrants->elem_count; ++zz, ++elem) {
      quad = p4est_quadrant_array_index(quadrants, zz);
      p4est_record_per_quad(&d, t, quad, elem, 0);
    }
  }

  for (zz=0; zz < g->ghost_layer->ghosts.elem_count; ++zz) {
    quad = sc_array_index(&(g->ghost_layer->ghosts), zz);

    p4est_record_per_quad(&d, quad->p.which_tree, quad, 0, 1);
  }
}

p4est_simplex_mesh_t *
create_simplex_mesh(context_t *g)
{
  p4est_simplex_mesh_t *mesh;
  p4est_nodes_t *nodes;
  sc_hash_array_t *hash_record;

  p4est_indep_t *cnode;

  size_t num_verticies;
  size_t num_elems, num_cnodes, num_fnodes;
  size_t off_elem_nodes, off_cnodes, off_fnodes;

#ifdef P4_TO_P8
  p8est_hang4_t *fnode;
  p8est_hang2_t *enode;

  size_t num_enodes;
  size_t off_enodes;
#else
  p4est_hang2_t *fnode;
#endif

#ifndef P4_TO_P8
  p4est_quad_face_node_t *found, qid;
#else
  p8est_quad_edge_node_t *found, qid;
#endif

  size_t zz, vi, fi, ei;
  p4est_locidx_t elem;
  p4est_topidx_t t, rt;
  p4est_tree_t *tree;
  sc_array_t *quadrants;
  p4est_quadrant_t *quad;
  
  int face;

  p4est_locidx_t cf, cc, c0, c1;

#ifdef P4_TO_P8
  int edge;
  int face_has_hanging_edge;

  p4est_locidx_t ce, c2, c3;

  p4est_locidx_t edge_noes[P8EST_EDGES];
#endif

  p4est_locidx_t *simplex;
  double abc[3], *vx;
  size_t position;

  // Generate p4est_nodes structure
  nodes = p4est_nodes_new(g->p4est, g->ghost_layer);

  mesh = P4EST_ALLOC(p4est_simplex_mesh_t, 1);


  mesh->vertices = sc_array_new(3*sizeof(double));

  p4est_record_hanging_faces(
      g,
      mesh->vertices,
      nodes,
      &hash_record
      );

    /// The vertex order is
    ///    ghost_nodes
    ///    cnodes
    ///    fnodes
    ///    enodes
    ///    elem_nodes
    ///    additional_face_nodes

  num_elems = nodes->num_local_quadrants;
  num_cnodes = nodes->indep_nodes.elem_count;
  num_fnodes = nodes->face_hangings.elem_count;

  off_cnodes = mesh->vertices->elem_count;
  off_fnodes = off_cnodes + num_cnodes;

#ifdef P4_TO_P8
  num_enodes = nodes->edge_hangings.elem_count;
  off_enodes = off_fnodes + num_fnodes;

  off_elem_nodes = off_enodes + num_enodes;

  num_verticies = num_elems + num_cnodes + num_fnodes + num_enodes;
#else
  off_elem_nodes = off_fnodes + num_fnodes;
  num_verticies = num_elems + num_cnodes + num_fnodes;
#endif

  sc_array_push_count(mesh->vertices, num_verticies);
  mesh->simplicies = sc_array_new((P4EST_DIM + 1) *sizeof(p4est_locidx_t));

  // Add all the vertices coming from independent nodes (indep_nodes)
  for (vi=0; vi < num_cnodes; ++vi) {
    cnode = sc_array_index(&nodes->indep_nodes, vi);
    vx = sc_array_index(mesh->vertices, off_cnodes + vi);

    abc[0] = (double) cnode->x / P4EST_ROOT_LEN;
    abc[1] = (double) cnode->y / P4EST_ROOT_LEN;
#ifdef P4_TO_P8
    abc[2] = (double) cnode->z / P4EST_ROOT_LEN;
#else
    abc[2] = 0;
#endif

    g->geom->X(g->geom, cnode->p.which_tree, abc, vx);
  }

  // Add all the vertices coming from face hangings nods (face_hangings)
  for (vi=0; vi < num_fnodes; ++vi) {
    fnode = sc_array_index(&nodes->face_hangings, vi);
    vx = sc_array_index(mesh->vertices, off_fnodes + vi);

    abc[0] = (double) fnode->x / P4EST_ROOT_LEN;
    abc[1] = (double) fnode->y / P4EST_ROOT_LEN;
#ifdef P4_TO_P8
    abc[2] = (double) fnode->z / P4EST_ROOT_LEN;
#else
    abc[2] = 0;
#endif

    g->geom->X(g->geom, fnode->p.which_tree, abc, vx);
  }

#ifdef P4_TO_P8
  // Add all the vertices coming from edge hangings nods (edge_hangings)
  for (vi=0; vi < num_enodes; ++vi) {
    enode = sc_array_index(&nodes->edge_hangings, vi);
    vx = sc_array_index(mesh->vertices, off_enodes + vi);

    abc[0] = (double) enode->x / P4EST_ROOT_LEN;
    abc[1] = (double) enode->y / P4EST_ROOT_LEN;
    abc[2] = (double) enode->z / P4EST_ROOT_LEN;

    g->geom->X(g->geom, enode->p.which_tree, abc, vx);
  }
#endif

  // We iterate through the local quads to make the simplices
  elem = 0;
  for (t=g->p4est->first_local_tree;
      t<=g->p4est->last_local_tree; ++t) {

    tree = p4est_tree_array_index(g->p4est->trees, t);
    qid.which_tree = t;

    quadrants = &tree->quadrants;
    for (zz = 0; zz < quadrants->elem_count; ++zz, ++elem) {
      quad = p4est_quadrant_array_index(quadrants, zz);
      qid.quad = *quad;

      // Add the quad centre node
      vx = sc_array_index(mesh->vertices, off_elem_nodes + elem);

      abc[0] = (double) ( quad->x + P4EST_QUADRANT_LEN(quad->level + 1) ) / P4EST_ROOT_LEN;
      abc[1] = (double) ( quad->y + P4EST_QUADRANT_LEN(quad->level + 1) ) / P4EST_ROOT_LEN;
#ifdef P4_TO_P8
      abc[2] = (double) ( quad->z + P4EST_QUADRANT_LEN(quad->level + 1) ) / P4EST_ROOT_LEN;
#else
      abc[2] = 0;
#endif
      g->geom->X(g->geom, t, abc, vx);
      cc = off_elem_nodes + elem;

#ifdef P4_TO_P8
      //  This should be in the record
      //  Could also be per face? This would avoid the first loop below
      //   However this way, there is only one lookup per quad in stead of per face

      for (ei=0; ei<P8EST_EDGES; ++ei) {
        qid.edge = ei;

        if (sc_hash_array_lookup(hash_record, &qid, &position)) {
          found = sc_array_index(&hash_record->a, position);

          if (found->node_index < 0)
            edge_noes[ei] = ~found->node_index;
          else 
            edge_noes[ei] = off_cnodes + found->node_index;
        }
        else
        {
          edge_noes[ei] = -1;
        }
      }
#endif

      for (fi = 0; fi < P4EST_FACES; ++fi) {
#ifndef P4_TO_P8
        qid.face = fi;

        c0 = off_cnodes +
          nodes->local_nodes[P4EST_CHILDREN*elem + p4est_face_corners[fi][0]];
        c1 = off_cnodes +
          nodes->local_nodes[P4EST_CHILDREN*elem + p4est_face_corners[fi][1]];

        // Check if the face is hanging
        if (sc_hash_array_lookup(hash_record, &qid, &position))
        {
          found = sc_array_index(&hash_record->a, position);
          if (found->node_index < 0)
            cf = ~found->node_index;
          else 
            cf = off_cnodes + found->node_index;

          simplex = sc_array_push(mesh->simplicies);
          simplex[0] = cc;
          simplex[1] = c0;
          simplex[2] = cf;

          simplex = sc_array_push(mesh->simplicies);
          simplex[0] = cc;
          simplex[1] = cf;
          simplex[2] = c1;
        }
        else
        {  // Split as usual
          simplex = sc_array_push(mesh->simplicies);
          simplex[0] = cc;
          simplex[1] = c0;
          simplex[2] = c1;
        }
#else
        // Could skip this first loop if we know no edges are hanging
        face_has_hanging_edge = 0;
        for (ei = 0; ei < P4EST_HALF; ++ei) {
          edge = p8est_face_edges[fi][ei];
          if (0 <= edge_noes[edge])
          {
            face_has_hanging_edge = 1;
            break;
          }
        }

        if (face_has_hanging_edge)
        {
          cf = mesh->vertices->elem_count;
          vx = sc_array_push(mesh->vertices);

          // Compute face center vertex
          P4EST_ASSERT(0 <= fi && fi < 6);

          p4est_face_center_coords(quad, fi, abc);

          g->geom->X(g->geom, t, abc, vx);

          for (ei=0; ei<4; ++ei) {
            edge = p8est_face_edges[fi][ei];

            // Check these indices
            c0 = off_cnodes +
              nodes->local_nodes[P4EST_CHILDREN*elem + p8est_edge_corners[edge][0]];
            c1 = off_cnodes +
              nodes->local_nodes[P4EST_CHILDREN*elem + p8est_edge_corners[edge][1]];

            ce = edge_noes[edge];

            if (ce < 0)
            {
              simplex = sc_array_push(mesh->simplicies);
              simplex[0] = cc;
              simplex[1] = cf;
              simplex[2] = c0;
              simplex[3] = c1;
            }
            else
            {
              // edge center
              simplex = sc_array_push(mesh->simplicies);
              simplex[0] = cc;
              simplex[1] = cf;
              simplex[2] = c0;
              simplex[3] = ce;

              simplex = sc_array_push(mesh->simplicies);
              simplex[0] = cc;
              simplex[1] = cf;
              simplex[2] = ce;
              simplex[3] = c1;
            }
          }
        }
        else
        {
          c0 = off_cnodes +
            nodes->local_nodes[P4EST_CHILDREN*elem + p4est_face_corners[fi][0]];
          c1 = off_cnodes +
            nodes->local_nodes[P4EST_CHILDREN*elem + p4est_face_corners[fi][1]];
          c2 = off_cnodes +
            nodes->local_nodes[P4EST_CHILDREN*elem + p4est_face_corners[fi][2]];
          c3 = off_cnodes +
            nodes->local_nodes[P4EST_CHILDREN*elem + p4est_face_corners[fi][3]];

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
    
            p4est_locidx_t first = c0;

            int rank = first < nodes->offset_owned_indeps 
                          ? nodes->nonlocal_ranks[first] 
                          : g->mpirank;
  

            // THIS IS SO HACKY!
            int rank1;
#define cmp_node(NODE) \
            rank1 = NODE < nodes->offset_owned_indeps  \
                          ? nodes->nonlocal_ranks[NODE] \
                          : g->mpirank; \
            if (rank1 < rank || (rank1 == rank && NODE < first)) { \
              first = NODE; \
              rank = rank1; \
            }

            cmp_node(c1);
            cmp_node(c2);
            cmp_node(c3);
#undef cmp_node

            rot_corners = (first == c1 || first == c2);
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

          simplex = sc_array_push(mesh->simplicies);
          simplex[0] = cc;
          simplex[1] = c0;
          simplex[2] = c1;
          simplex[3] = c3;

          simplex = sc_array_push(mesh->simplicies);
          simplex[0] = cc;
          simplex[1] = c3;
          simplex[2] = c2;
          simplex[3] = c0;
        }
#endif
      }  // for (fi = 0..P4EST_FACES)
    } // for (zz = 0..|quadtrees|)
  } // for (t = first-..last-localtrees)

  sc_hash_array_destroy(hash_record);
  p4est_nodes_destroy(nodes);

  return mesh;
}

int
main(int argc, char **argv) {
  int mpiret;
  int first_argc;
  int exitcode;
  sc_options_t *sc_opts;
  options_t options, *opts = &options;
  context_t context = {0}, *g = &context;
  p4est_simplex_mesh_t *smesh;
  char filepath[1024];


  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  g->mpicomm = sc_MPI_COMM_WORLD;
  mpiret = sc_MPI_Comm_size (g->mpicomm, &g->mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (g->mpicomm, &g->mpirank);
  SC_CHECK_MPI (mpiret);

  /* Does this really need to be before the parsing? */
  sc_init (g->mpicomm, 1, 1, NULL, SC_LP_DEFAULT);

  // User info flags
  int help = 0;
  int list_conn = 0;

  smesh = NULL;

  /* Initialize default option values */
  opts->minlevel = 2;
  opts->maxlevel = 5;
  opts->conn = "unit";
  opts->mshpath = NULL;  // Defaults to "{opts->conn}.msh"

  sc_opts = sc_options_new (argv[0]);
  sc_options_add_int (sc_opts, 'l', "minlevel", &opts->minlevel, opts->minlevel, "Lowest level");
  sc_options_add_int (sc_opts, 'L', "maxlevel", &opts->maxlevel, opts->maxlevel, "Highest level");
  // sc_options_add_string (sc_opts, 'v', "vtk", &opts->vtk, opts->vtk, "VTK basename");
  sc_options_add_string (sc_opts, 'm', "msh", &opts->mshpath, opts->mshpath, "MSH base filename");
  sc_options_add_string (sc_opts, 'c', "conn", &opts->conn, opts->conn, "Name of the connectivity");


  sc_options_add_bool (sc_opts, 'h', "help", &help, help, "Print help message");
  sc_options_add_bool (sc_opts, 0, "list-conn", &list_conn, list_conn, "Lists the available connectivity 'conn' options");

  exitcode = 0;
  do {
    first_argc = sc_options_parse(p4est_package_id, SC_LP_DEFAULT, sc_opts, argc, argv);
    if (first_argc < 0) {
      g->errmsg = "Failed to parse option format";
      sc_options_print_usage (p4est_package_id, SC_LP_ERROR, sc_opts, NULL);
      exitcode = 1;
      break;
    }

    if (first_argc != argc) {
      g->errmsg = "Too many arguments";
      sc_options_print_usage (p4est_package_id, SC_LP_ERROR, sc_opts, NULL);
      exitcode = 1;
      break;
    }

    if (help) {
      sc_options_print_usage (p4est_package_id, SC_LP_INFO, sc_opts, NULL);
      exitcode = 0;
      break;
    }

    if (list_conn) {
      char **conn_name = p4est_connectivity_names;
      for (;*conn_name; ++conn_name)
        printf("  %s\n", *conn_name);
      exitcode = 0;
      break;
    }

    g->opts = opts;
    // sc_init (g->mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
    p4est_init (NULL, SC_LP_DEFAULT);

    if (create_context_from_opts(g, opts)) {
      exitcode = 1;
      break;
    }

    smesh = create_simplex_mesh(g);

#ifndef P4_TO_P8
    sprintf(filepath, "out2d/box_%s_%d-%d", opts->conn, opts->minlevel, opts->maxlevel);
#else
    sprintf(filepath, "out3d/box_%s_%d-%d", opts->conn, opts->minlevel, opts->maxlevel);
#endif

    p4est_vtk_write_file(g->p4est, g->geom, filepath);

    // sprintf(filepath, "out/sx_%s.msh", g->opts->conn);
    // write_gmsh_file(filepath, smesh);

#ifndef P4_TO_P8
    sprintf(filepath, "out2d/sx_%s_%d-%d_%d.vtk", opts->conn, opts->minlevel, opts->maxlevel, g->mpirank);
#else
    sprintf(filepath, "out3d/sx_%s_%d-%d_%d.vtk", opts->conn, opts->minlevel, opts->maxlevel, g->mpirank);
#endif
    write_vtk_file(filepath, smesh);

  } while (0);
  if (g->errmsg) {
    P4EST_GLOBAL_LERROR (g->errmsg);
  }

  if (smesh) {
    p4est_simplex_mesh_destroy(smesh);
  }

  /*** clean up and exit ***/
  destroy_context(g);

  sc_options_destroy (sc_opts);
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return exitcode ? EXIT_FAILURE : EXIT_SUCCESS;
}

