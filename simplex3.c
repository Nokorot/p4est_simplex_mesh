#include "p4est_to_p8est.h"

#include "p8est_bits.h"

#include "simplex3.h"

#include "simplex2.c"

unsigned int
p4est_quad_edge_hash_fn(const p8est_quad_edge_node_t *v, const void *user_data)
{
  unsigned int qhash = p4est_quadrant_hash_fn(&v->quad, user_data);
  int b, c;
  b = v->which_tree;
  c = v->edge;

  sc_hash_mix(qhash, b, c);
  sc_hash_final(qhash, b, c);

  return (unsigned) c;
}

unsigned int
p4est_quad_edge_equal_fn(
    const p8est_quad_edge_node_t *v1,
    const p8est_quad_edge_node_t *v2,
    const void *used_data)
{
  return p4est_quadrant_equal_fn(&v1->quad, &v2->quad, NULL) \
    && v1->which_tree == v2->which_tree
    && v1->edge == v2->edge;
}

void
p8est_edge_center_coords(
    p4est_quadrant_t *quad,
    int edge,
    double abc[3]
    )
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


