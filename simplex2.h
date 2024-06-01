

#ifndef P4_TO_P8
#include <p4est_bits.h>
#include <p4est_ghost.h>
#include <p4est_geometry.h>
#else
#include <p8est_bits.h>
#include <p8est_ghost.h>
#include <p8est_geometry.h>
#endif

#include "p4est_simplex_to_gmsh.h"


static char *p4est_connectivity_names[] = {
#ifndef P4_TO_P8
  "brick23","corner","cubed","disk","icosahedron","moebius",
  "periodic","pillow","rotwrap","star","shell2d","disk2d",
  "bowtie","unit","sphere2d",
#else
  "brick235","periodic","rotcubes","rotwrap","shell","sphere",
  "twocubes","twowrap","unit",
#endif
  NULL
};


typedef struct
{
  int minlevel, maxlevel;
  const char *conn;
  const char *mshpath;
}
options_t;

typedef struct
{
  sc_MPI_Comm         mpicomm;
  int                 mpisize, mpirank;

  options_t *opts;

  p4est_t *p4est;
  p4est_connectivity_t *conn;
  p4est_ghost_t *ghost_layer;
  p4est_geometry_t *geom;

  char *errmsg;
}
context_t;

typedef struct {
    p4est_quadrant_t quad;
    p4est_locidx_t which_tree;
    int8_t face;

    p4est_locidx_t node_index;
}
p4est_quad_face_node_t;

unsigned int
p4est_quad_face_hash_fn(const p4est_quad_face_node_t *v, const void *user_data);

unsigned int
p4est_quad_face_equal_fn(
    const p4est_quad_face_node_t *v1,
    const p4est_quad_face_node_t *v2,
    const void *used_data);


#ifdef P4_TO_P8
typedef struct {
    p8est_quadrant_t quad;
    p4est_locidx_t which_tree;
    int8_t edge;

    p4est_locidx_t node_index;
}
p4est_quad_edge_node_t;

unsigned int
p4est_quad_edge_hash_fn(const p4est_quad_edge_node_t *v, const void *user_data);
#endif



p4est_simplex_mesh_t *
create_simplex_mesh(context_t *g);

