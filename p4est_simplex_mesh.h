#ifndef P4EST_SIMPLEX_MESH_H_
#define P4EST_SIMPLEX_MESH_H_

#ifndef P4_TO_P8
#include <p4est_tnodes.h>
#include <p4est_bits.h>
#include <p4est_lnodes.h>
#include <p4est_ghost.h>
#include <p4est_geometry.h>
#else
#include <p8est_tnodes.h>
#include <p8est_lnodes.h>
#include <p8est_ghost.h>
#include <p8est_geometry.h>
#endif

typedef struct p4est_snodes {
  union {
    struct {
      // lnode
      sc_MPI_Comm         mpicomm;
      p4est_locidx_t      num_local_nodes;
      p4est_locidx_t      owned_count;
      p4est_gloidx_t      global_offset;
      p4est_gloidx_t     *nonlocal_nodes;
      sc_array_t         *sharers;
      p4est_locidx_t     *global_owned_count;

      int                 degree, vnodes; // Not defined
      p4est_locidx_t      num_local_elements;
      p4est_lnodes_code_t *face_code; //  Not defined
      p4est_locidx_t      *element_nodes;
      // end lnode
    };
    p4est_lnodes_t lnodes;
  };
}
p4est_snodes_t;

p4est_tnodes_t *
p4est_new_simplex_mesh(
    p4est_t *p4est,
    p4est_geometry_t *geometry,
    p4est_ghost_t *ghost);

#endif
