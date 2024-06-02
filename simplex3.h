#ifndef SIMPLEX3_H_
#define SIMPLEX3_H_

#include <p8est_bits.h>

typedef struct {
    p8est_quadrant_t quad;
    p4est_locidx_t which_tree;
    int8_t edge;

    p4est_locidx_t node_index;
}
p8est_quad_edge_node_t;


unsigned int
p4est_quad_edge_hash_fn(const p8est_quad_edge_node_t *v, const void *user_data);

unsigned int
p4est_quad_edge_equal_fn(
    const p8est_quad_edge_node_t *v1,
    const p8est_quad_edge_node_t *v2,
    const void *used_data);

void
p8est_edge_coords(
    p8est_quadrant_t *quad,
    int edge,
    double abc[3]
    );

#endif
