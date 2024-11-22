#ifndef P4EST_SIMPLEX_MESH_H_
#define P4EST_SIMPLEX_MESH_H_

#ifdef VIM_LS
// #include <p4est_to_p8est.h>
#define P4EST_SIMPLEX_DEBUG
// #define TIMINGS
#endif


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

#ifdef TIMINGS
#include "statistics.h"
#endif

p4est_tnodes_t *
p4est_new_BECK_mesh(
    p4est_t *p4est,
    p4est_geometry_t *geometry,
    p4est_ghost_t *ghost
#ifdef TIMINGS
    , statistics_t *ss
#endif
    );

#endif
