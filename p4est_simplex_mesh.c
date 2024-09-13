
#include <bits/types/locale_t.h>
#include <sc.h>
#include <sc_containers.h>
#include <sc_mpi.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "utils.c"

#ifdef VIM_LS
#include <p4est_to_p8est.h>
#define P4EST_SIMPLEX_DEBUG
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



// #ifndef VIM_LS
/** Shortcut to access full face information */
// typedef struct p4est_iter_face_side_full
// {
//   int8_t              is_ghost; /**< boolean: local (0) or ghost */
//   p4est_quadrant_t   *quad;     /**< the actual quadrant */
//   p4est_locidx_t      quadid;   /**< index in tree or ghost array */
//
// }
// p4est_iter_face_side_full_t;
//
// /** Shortcut to access hanging face information */
// typedef struct p4est_iter_face_side_hanging
// {
//   int8_t              is_ghost[2];      /**< boolean: local (0) or ghost */
//   p4est_quadrant_t   *quad[2];          /**< the actual quadrant */
//   p4est_locidx_t      quadid[2];        /**< index in tree or ghost array */
// }
// p4est_iter_face_side_hanging_t;


#ifdef P4_TO_P8
/** Shortcut to access full face information */
// typedef struct p8est_iter_face_side_full
// {
//   int8_t              is_ghost; /**< boolean: local (0) or ghost */
//   p8est_quadrant_t   *quad;     /**< the actual quadrant */
//   p4est_locidx_t      quadid;   /**< index in tree or ghost array */
//
// }
// p8est_iter_face_side_full_t;
//
// /** Shortcut to access hanging face information */
// typedef struct p8est_iter_face_side_hanging
// {
//   int8_t              is_ghost[4];      /**< boolean: local (0) or ghost */
//   p8est_quadrant_t   *quad[4];          /**< the actual quadrant */
//   p4est_locidx_t      quadid[4];        /**< index in tree or ghost array */
// }
// p8est_iter_face_side_hanging_t;

/** Shortcut to access full edge information */
typedef struct p8est_iter_edge_side_full
{
  int8_t              is_ghost; /**< boolean: local (0) or ghost */
  p8est_quadrant_t   *quad;     /**< the actual quadrant */
  p4est_locidx_t      quadid;   /**< index in tree or ghost array */

}
p8est_iter_edge_side_full_t;

/** Shortcut to access hanging edge information */
typedef struct p8est_iter_edge_side_hanging
{
  int8_t              is_ghost[2];      /**< boolean: local (0) or ghost */
  p8est_quadrant_t   *quad[2];          /**< the actual quadrant */
  p4est_locidx_t      quadid[2];        /**< index in tree or ghost array */
}
p8est_iter_edge_side_hanging_t;

#define p4est_iter_face_side_full_t     p8est_iter_face_side_full_t
#define p4est_iter_face_side_hanging_t  p8est_iter_face_side_hanging_t

#define p4est_iter_face_side_full_t     p8est_iter_face_side_full_t
#define p4est_iter_face_side_hanging_t  p8est_iter_face_side_hanging_t
#endif

// #endif





// Local declarations


static int P4EST_COMM_SNODES_QUERY = 50;
static int P4EST_COMM_SNODES_REPLY = 51;


typedef struct {
  p4est_locidx_t volume_node;

  p4est_locidx_t face_nodes[ P4EST_FACES ];
#ifdef P4_TO_P8
  p4est_locidx_t edge_nodes[ P8EST_EDGES ];
#endif
  p4est_locidx_t corner_nodes[ P4EST_CHILDREN ];
}
element_node_t;

#define SNODES_VNODES P4EST_DIM_POW(3)

/***
 *
 *    5 --- 4 --- 6
 *    |           |
 *    1     0     2
 *    |           |
 *    7 --- 3 --- 8
 */


/* cube center */
static const int n_center = 0;

#ifdef P4_TO_P8
/* face midpoints  */
static const int n_face[6] = { 1, 2,
                               3, 4,
                               5, 6 };

/* edge midpoints  */
static const int n_edge[12] = {  7,  8,  9, 10,
                                11, 12, 13, 14,
                                15, 16, 17, 18 };

/* cube corners */
static const int n_corner[8] = { 19, 20, 21, 22,
                                 23, 24, 25, 26 };

#ifdef P4EST_ENABLE_DEBUG
static const int    alwaysowned[27] =
  { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
#endif

#else
/* face midpoints  */
static const int n_face[4] = { 1, 2, 3, 4 };

/* cube corners */
static const int n_corner[4] = { 5, 6, 7, 8 };

#ifdef P4EST_ENABLE_DEBUG
static const int    alwaysowned[9] =
  { 1, 0, 0, 0, 0, 0, 0, 0, 0 };
#endif
#endif


#define XYZ(i, j, k)  { i , j, k }
// #define XYZ(i, j, k)  (k << 4 | j << 2 | i)

static const int8_t elem_coords_n[ SNODES_VNODES ][3] = {
  // Center
  XYZ(1,1,1),
#ifdef P4_TO_P8
  // Faces
  XYZ(0,1,1), XYZ(2,1,1), // x-orth
  XYZ(1,0,1), XYZ(1,2,1), // y-orth
  XYZ(1,1,0), XYZ(1,1,2), // z-orth
  // Edges
  XYZ(1,0,0), XYZ(1,2,0), XYZ(1,0,2), XYZ(1,2,2), // x-parl
  XYZ(0,1,0), XYZ(2,1,0), XYZ(0,1,2), XYZ(2,1,2), // y-parl
  XYZ(0,0,1), XYZ(2,0,1), XYZ(0,2,1), XYZ(2,2,1), // z-parl
  // Corners
  XYZ(0,0,0), XYZ(2,0,0), XYZ(0,2,0), XYZ(2,2,0),
  XYZ(0,0,2), XYZ(2,0,2), XYZ(0,2,2), XYZ(2,2,2)
#else
  // Faces
  XYZ(0,1,0), XYZ(2,1,0), // x-orth
  XYZ(1,0,0), XYZ(1,2,0), // y-orth
  // Corners
  XYZ(0,0,0), XYZ(2,0,0), XYZ(0,2,0), XYZ(2,2,0)
#endif
};
#undef XYZ


typedef enum {
  NODE_TYPE_VOLUME = 0,
  NODE_TYPE_FACE   = 1,
#ifdef P4_TO_P8
  NODE_TYPE_EDGE   = 1 + P4EST_FACES,
  NODE_TYPE_CORNER = 1 + P4EST_FACES + P8EST_EDGES,
#else
  NODE_TYPE_CORNER = 1 + P4EST_FACES,
#endif
}
snode_connect_type_t;


typedef struct
{
  sc_MPI_Comm mpicomm;
  int mpirank, mpisize;
  int *ghost_rank;
  int emptypeers;
  int locsharer;         /**< Index of local sharer in sharers */
  uint8_t *chilev;

  p4est_locidx_t *element_nodes;
  sc_array_t construct;  /**< Collect nodes during traversal */
  sc_array_t ownsort;    /**< Sorted owned nodes of local process */

  p4est_locidx_t lenum;              // Element number
  p4est_locidx_t num_owned;
  p4est_locidx_t num_owned_shared;           /**< Nodes we both own and share */
  p4est_locidx_t num_notowned_shared;        /**< Nodes we share, owned or not */
  p4est_locidx_t num_shared;                 /**< Nodes we share, owned or not */

  p4est_gloidx_t *goffset;          /**< Global offsets for owned nodes */


  // element configuration might be useful

  // sc_array_t *tree_offsets;
  // sc_array_t *global_owned_count;
  // sc_array_t *owners;

#ifdef P4EST_ENABLE_MPI
  int *proc_peer;
  sc_array_t          sortp;        /**< Sorted array pointing to peers */
  sc_array_t          peers;        /**< Unsorted peer storage */
  sc_array_t          peerreq;        /**< Requests for storage */
#endif

  p4est_snodes_t *snodes;
  // p4est_locidx_t num_local_elements;

  p4est_t *p4est;
  p4est_ghost_t *ghost;
}
snodes_meta_t;


typedef struct snodes_contr {
  int nodene;           /**< Relative to element. */
  int rank;             /**< The referring process. */
  p4est_locidx_t le;    /**< Element/ghost number. */
} snodes_contr_t;

typedef struct snodes_cnode {
  p4est_locidx_t runid;         /**< Running count of node. */
  snodes_contr_t *owner;         /**< Owning contributor. */
  sc_array_t contr;             /**< Contributing processes. */
} snodes_cnode_t;

/** Record one communication partner and/or node sharer. */
typedef struct snodes_peer
{
  int                 rank;         /**< Rank of the peer process */
  int                 done;
  int                 sharind;      /**< Index of corresponding sharer */
  int                 passive;      /**< Number of passively shared nodes */
  p4est_locidx_t      lastadd;      /**< Most recently added node number */
  p4est_locidx_t      bufcount;     /**< Number items in message buffer */
  p4est_locidx_t      shacumul;     /**< Number owned nodes before peer */
  sc_array_t          sharedno;     /**< Remember local node with query */
  sc_array_t          querypos;     /**< Send/receive buffer for messages */
  sc_array_t          remosort;     /**< Pointer array to sort peer nodes */
}
snodes_peer_t;

typedef struct {
  union {
    p4est_locidx_t node_index;
    p4est_locidx_t node_count;
  };

#ifdef P4EST_SIMPLEX_DEBUG
  double x, y, z;
#endif
}
shared_node_comm_data_t;



/**** static function declarations ******/

static void
iter_volume(p4est_iter_volume_info_t *info, void *user_data);

static void
iter_face(p4est_iter_face_info_t *info, void *user_data);

#ifdef P4_TO_P8
static void
iter_edge(p8est_iter_edge_info_t *info, void *user_data);
#endif

static void
iter_corner(p4est_iter_corner_info_t *info, void *user_data);

// void
// p4est_snodes_new(snodes_meta_t *d);
p4est_snodes_t *
p4est_snodes_new (p4est_t *p4est, p4est_ghost_t *ghost);



// Implementations
static void
check_node (snodes_meta_t * me, p4est_locidx_t lni)
{
#ifdef P4EST_ENABLE_DEBUG
  snodes_cnode_t     *cnode;
  snodes_contr_t     *contr;
  int                 owner_rank;
  size_t              zz, siz;

  cnode = (snodes_cnode_t *) sc_array_index (&me->construct, (size_t) lni);
  P4EST_ASSERT (cnode->runid == lni);
  owner_rank = cnode->owner->rank;
  siz = cnode->contr.elem_count;
  for (zz = 0; zz < siz; ++zz) {
    contr = (snodes_contr_t *) sc_array_index (&cnode->contr, zz);
    P4EST_ASSERT (owner_rank <= contr->rank);
    if (owner_rank == contr->rank) {
      P4EST_ASSERT (contr == cnode->owner);
    }
  }
#endif /* P4EST_ENABLE_DEBUG */
}


static p4est_locidx_t
tree_quad_to_le (p4est_t * p4est, p4est_topidx_t treeid,
                 p4est_locidx_t quadid)
{
  p4est_tree_t       *tree;

  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (p4est->trees != NULL);
  tree = p4est_tree_array_index (p4est->trees, treeid);

  P4EST_ASSERT (0 <= quadid &&
                quadid < (p4est_locidx_t) tree->quadrants.elem_count);
  return tree->quadrants_offset + quadid;
}


/**
 *
 * me : meta_data,
 * lni : index pointer,
 * rank : rank of owner,
 * le :     local element index,
 * nodene : element node index,
 */
void
node_register(
    snodes_meta_t *me, p4est_locidx_t *lni,
    int rank, p4est_locidx_t le, int nodene
    )
{
  p4est_locidx_t lnis;
  snodes_cnode_t *cnode;
  snodes_contr_t *contr;
  int owner_rank;
  size_t size, zz;

  /* basic checks */
  P4EST_ASSERT (me != NULL);
  P4EST_ASSERT(0 <= le && le < me->snodes->num_local_elements);
  P4EST_ASSERT(0 <= nodene && nodene < SNODES_VNODES);


  /* a new node is to be created or an existing one is passed in */
  if (lni == NULL) {
    lnis = -1;
    lni = &lnis;
  }

  P4EST_ASSERT (-1 == *lni ||
                *lni < (p4est_locidx_t) me->construct.elem_count);

  /* abbreviate local rank */
  P4EST_ASSERT (rank >= -1 && rank != me->mpirank);
  if (rank == -1) {
    rank = me->mpirank;
  }

  P4EST_ASSERT (0 <= rank && rank < me->mpisize);

  /* check remaining arguments */
  P4EST_ASSERT (0 <= le && le < (p4est_locidx_t)
                (me->p4est->global_first_quadrant[rank + 1] -
                 me->p4est->global_first_quadrant[rank]));

  /* if a node is already registered we don't need to do anything */
  if (rank == me->mpirank && me->element_nodes[le * SNODES_VNODES + nodene] != -1) {
    P4EST_ASSERT (*lni == -1 ||
                  *lni == me->element_nodes[le * SNODES_VNODES + nodene]);

    *lni = me->element_nodes[le * SNODES_VNODES + nodene];
    return;
  }


  if (*lni == -1) {
    /* create a new node with one instance */
    *lni = (p4est_locidx_t) me->construct.elem_count;
    cnode = (snodes_cnode_t *) sc_array_push (&me->construct);
    cnode->runid = *lni;
    cnode->owner = NULL;
    sc_array_init (&cnode->contr, sizeof (snodes_contr_t));
  }
  else {
    /* create a new instance of an existing node */
    cnode = (snodes_cnode_t *) sc_array_index (&me->construct, (size_t) *lni);

    P4EST_ASSERT (cnode->runid == *lni);
    P4EST_ASSERT (cnode->contr.elem_size == sizeof (snodes_contr_t));
    P4EST_ASSERT (cnode->contr.elem_count > 0);
    P4EST_ASSERT (cnode->owner != NULL);
    check_node (me, *lni);
  }

  /* assign node to the local element position */
  if (rank == me->mpirank) {
    P4EST_ASSERT (me->element_nodes[le * SNODES_VNODES + nodene] == -1);
    me->element_nodes[le * SNODES_VNODES + nodene] = *lni;
  }

  /* iterate through instances to find matching process */
  size = cnode->contr.elem_count;
  P4EST_ASSERT (size == 0 || cnode->owner != NULL);
  for (zz = 0; zz < size; ++zz) {
    contr = (snodes_contr_t *) sc_array_index (&cnode->contr, zz);
    P4EST_ASSERT (cnode->owner->rank <= contr->rank);
    if (contr->rank == rank) {
      /* rank is found and we remember the smallest node position */
      if (le < contr->le || (le == contr->le && nodene < contr->nodene)) {
        contr->nodene = nodene;
        contr->le = le;
      }
      check_node (me, *lni);
      return;
    }
  }

  /* add new node process to the list for this node */
  owner_rank = cnode->owner == NULL ? -1 : cnode->owner->rank;
  contr = (snodes_contr_t *) sc_array_push (&cnode->contr);
  contr->nodene = nodene;
  contr->rank = rank;
  contr->le = le;
  if (cnode->owner == NULL || rank <= owner_rank) {
    /* the array was empty before or we know to set a new owner */
    P4EST_ASSERT (rank != owner_rank);
    cnode->owner = contr;
  }
  else {
    // TODO: Consider using an offset instead
    /* pushing to the array has invalidated the previous owner pointer */
    cnode->owner = (snodes_contr_t *) sc_array_index (&cnode->contr, 0);
    size = cnode->contr.elem_count;
    for (zz = 1; zz < size - 1; ++zz) {
      contr = (snodes_contr_t *) sc_array_index (&cnode->contr, zz);
      if (contr->rank < cnode->owner->rank) {
        cnode->owner = contr;
      }
    }
  }

  check_node(me, *lni);
}

void
node_lregister(
    snodes_meta_t *me, p4est_locidx_t *lni,
    p4est_locidx_t le, int nodene
    )
{
  node_register(me, lni, -1, le, nodene);
}

void
node_gregister(
    snodes_meta_t *me, p4est_locidx_t *lni,
    p4est_locidx_t ghostid, int nodene
    )
{
  p4est_quadrant_t *gquad;

  P4EST_ASSERT (me != NULL);
  // P4EST_ASSERT (me->sm != NULL);
  P4EST_ASSERT (0 <= nodene && nodene < SNODES_VNODES);

  if (me->ghost != NULL) {
    P4EST_ASSERT (me->ghost_rank != NULL);
    P4EST_ASSERT (0 <= ghostid &&
                  ghostid < (p4est_locidx_t) me->ghost->ghosts.elem_count);

    /* extract remote element number from ghost quadrant */
    gquad = (p4est_quadrant_t *) sc_array_index (&me->ghost->ghosts, ghostid);
    node_register (me, lni, me->ghost_rank[ghostid],
                   gquad->p.piggy3.local_num, nodene);
  }
}


// ITERATORS
void
iter_volume(p4est_iter_volume_info_t *vi, void *user_data)
{
  snodes_meta_t *me = (snodes_meta_t *) user_data;

  p4est_locidx_t le;
  p4est_locidx_t childid;
  int8_t level;

  le = me->lenum++;
  level = vi->quad->level;
  childid = p4est_quadrant_child_id (vi->quad);
  me->chilev[le] = (((uint8_t) level) << 3) | ((uint8_t) childid);

  node_lregister(me, NULL, le, n_center);
}

static void
iter_face(p4est_iter_face_info_t *fi, void *user_data)
{
  snodes_meta_t *me = (snodes_meta_t *) user_data;

  int nodene;
  p4est_locidx_t le;              /**< local element number */
  p4est_locidx_t lni;             /**< local node number */

  p4est_iter_face_side_t *side, *sf, *sh;
  p4est_iter_face_side_full_t *sfd;
  p4est_iter_face_side_hanging_t *shd;

  size_t sz, zz;
  int face;

  /* initial checks  */
  P4EST_ASSERT (fi->p4est == me->p4est);

  /* We check if there is a hanging side */
  sh = NULL;
  for (zz = 0; zz < fi->sides.elem_count; ++zz) {
    side = sc_array_index(&fi->sides, zz);

    if (side->is_hanging) {
      sh = side;
      shd = (p4est_iter_face_side_hanging_t *) &sh->is.hanging;
    }
    else {
      sf = side;
      sfd = (p4est_iter_face_side_full_t *) &sf->is.full;
    }
  }

  /* we only add a node if one of the sides are hanging */
  if (!sh)
    return;

  P4EST_ASSERT(sh && shd && sf && shd);


  lni = -1;
  nodene = n_face[sf->face];

  if (!sfd->is_ghost) {
    le = tree_quad_to_le(me->p4est, sf->treeid, sfd->quadid);
    lni = me->element_nodes[le*SNODES_VNODES + nodene];
    // if (lni >= 0)
    //   return;

    node_lregister(me, &lni, le, n_face[sf->face]);
  } else {
    node_gregister(me, &lni, sfd->quadid, n_face[sf->face]);
  }


  face = sh->face;

  for (sz = 0; sz < P4EST_HALF; ++sz) {
    nodene = n_corner[p4est_face_corners[face][sz ^ (P4EST_HALF - 1)]];
    if (!shd->is_ghost[sz]) {
      le = tree_quad_to_le(me->p4est, sh->treeid, shd->quadid[sz]);
      node_lregister(me, &lni, le, nodene);
    }
    else {
      node_gregister(me, &lni, shd->quadid[sz], nodene);
    }

    // // Compute the corner that is on the centre of the parent face
    // int cid = p4est_quadrant_child_id(quad);
    // int cfc = p4est_corner_face_corners[cid][face];
    // int c = p4est_face_corners[face][cfc ^ (P4EST_HALF - 1)];

    // p4est_simplex_record_node(me,
    //     sh->treeid,
    //     shd->quadid[sz],
    //     shd->is_ghost[sz],
    //     NODE_TYPE_CORNER + c,
    //     cf);
  }
}

#ifdef P4_TO_P8
void
iter_edge(p8est_iter_edge_info_t *info, void *user_data)
{
  snodes_meta_t *me = (snodes_meta_t *) user_data;

  p4est_locidx_t le;
  p4est_locidx_t lni;

  int nodene;

  p8est_iter_edge_side_t *side, *sh, *sf;
  p8est_iter_edge_side_full_t *sfd;
  p8est_iter_edge_side_hanging_t *shd;

  // p4est_topidx_t treeid;
  // p4est_locidx_t quadid;
  p4est_quadrant_t *quad;

  size_t sz, zz, efz;
  int face, edge;
  int faces_count;

  // We determine the number of adjacent faces
  faces_count = 0;
  for (zz = 0; zz < info->sides.elem_count; ++zz) {
    side = sc_array_index(&info->sides, zz);
    faces_count = SC_MAX(faces_count,
                    1 + SC_MAX(side->faces[0], side->faces[1]));
  }

  p4est_locidx_t *face_nodes;
  face_nodes = P4EST_ALLOC_ZERO(p4est_locidx_t, faces_count);
  // memset(face_nodes, 0xff, sizeof(p4est_locidx_t) * faces_count);

  // Find a full and a hanging face
  sf = NULL;
  sh = NULL;
  for (zz = 0; zz < info->sides.elem_count; ++zz) {
    side = sc_array_index(&info->sides, zz);
    if (side->is_hanging) {
      sh = side;
    }
    else {
      sf = side;

      if (!side->is.full.is_ghost) {
        face_nodes[side->faces[0]] = 1;
        face_nodes[side->faces[1]] = 1;
      }
    }
  }

  // If no side is hanging we will not add an edge node
  if (!sh) {
    P4EST_FREE(face_nodes);
    return;
  }

  SC_CHECK_ABORT(sf, "Got edge with only hanging sides");

  // Add the edge node
  // quad = side_full->is.full.quad;
  // p8est_quad_edge_center_coords(quad, side_full->edge, abc);
  // ce = p4est_simplex_push_node(me, side_full->treeid, abc, owner);

  lni = -1;

  // We record the index in the element_nodes for each side
  for (zz = 0; zz < info->sides.elem_count; ++zz) {
    side = sc_array_index(&info->sides, zz);
    edge = side->edge;

    if (!side->is_hanging) {
      sfd = (p8est_iter_edge_side_full_t *) &side->is.full;

      nodene = n_edge[side->edge];
      if (!sfd->is_ghost) {
        le = tree_quad_to_le(me->p4est, side->treeid, sfd->quadid);
        node_lregister(me, &lni, le, nodene);
      } else {
        node_gregister(me, &lni, sfd->quadid, nodene);
      }

      // For each face adjacent to a hanging edge, we also add a face node
      for (efz = 0; efz < 2; ++efz) {
        if (!face_nodes[ side->faces[efz] ])
          continue;

        nodene = n_face[ p8est_edge_faces[edge][efz] ];

        if (!sfd->is_ghost) {
          le = tree_quad_to_le(me->p4est, side->treeid, sfd->quadid);
          node_lregister(me, NULL, le, nodene);
        } else {
          node_gregister(me, NULL, sfd->quadid, nodene);
        }

      }
    } else {
      shd = (p8est_iter_edge_side_hanging_t *) &side->is.hanging;

      for (sz = 0; sz < 2; ++sz) {
        quad = shd->quad[sz];

        int cid = p4est_quadrant_child_id(quad);
        int cfc = p8est_corner_edge_corners[cid][edge];
        // int c = p8est_edge_corners[edge][cfc ^ 1];

        nodene = n_corner[p8est_edge_corners[edge][cfc ^ 1]];
        if (!shd->is_ghost[sz]) {
          le = tree_quad_to_le(me->p4est, side->treeid, shd->quadid[sz]);
          node_lregister(me, &lni, le, nodene);
        }
        else {
          node_gregister(me, &lni, shd->quadid[sz], nodene);
        }
      }
    }
  }

  P4EST_FREE(face_nodes);
}
#endif

void
iter_corner(p4est_iter_corner_info_t *info, void *user_data)
{
  snodes_meta_t *me = (snodes_meta_t *) user_data;

  p4est_iter_corner_side_t *side;

  p4est_locidx_t le;
  p4est_locidx_t lni;
  int nodene;

  size_t zz;

  lni = -1;
  // Record this corner_node for each of the adjacent elements
  for (zz = 0; zz < info->sides.elem_count; ++zz) {
    side = sc_array_index(&info->sides, zz);


    nodene = n_corner[side->corner];

    if (!side->is_ghost) {
      le = tree_quad_to_le(me->p4est, side->treeid, side->quadid);
      node_lregister(me, &lni, le, nodene);
    }
    else {
      node_gregister(me, &lni, side->quadid, nodene);
    }
  }
}

static int
cnode_compare (const void *v1, const void *v2)
{
  const snodes_cnode_t **cc1 = (const snodes_cnode_t **) v1;
  const snodes_cnode_t **cc2 = (const snodes_cnode_t **) v2;
  const snodes_contr_t *o1, *o2;
  p4est_locidx_t      ldiff;

  P4EST_ASSERT (cc1 != NULL && *cc1 != NULL);
  P4EST_ASSERT (cc2 != NULL && *cc2 != NULL);

  P4EST_ASSERT ((*cc1)->contr.elem_size == sizeof (snodes_contr_t));
  P4EST_ASSERT ((*cc2)->contr.elem_size == sizeof (snodes_contr_t));

  /* we sort within the same owner process */
  o1 = (*cc1)->owner;
  o2 = (*cc2)->owner;
  P4EST_ASSERT (o1->rank == o2->rank);

  /* if the elements are the same, compare node position */
  if (o1->le == o2->le) {
    return o1->nodene - o2->nodene;
  }

  /* nodes are sorted according to their element number */
  ldiff = o1->le - o2->le;
  return ldiff == 0 ? 0 : ldiff < 0 ? -1 : 1;
}

#ifdef P4EST_ENABLE_MPI

static int
peer_compare (const void *v1, const void *v2)
{
  const snodes_peer_t **p1 = (const snodes_peer_t **) v1;
  const snodes_peer_t **p2 = (const snodes_peer_t **) v2;
  return (*p1)->rank - (*p2)->rank;
}

static int
rnode_compare (const void *v1, const void *v2)
{
  const snodes_cnode_t **cc1 = (const snodes_cnode_t **) v1;
  const snodes_cnode_t **cc2 = (const snodes_cnode_t **) v2;
#ifdef P4EST_ENABLE_DEBUG
  const snodes_contr_t *o1, *o2;
#endif
  p4est_locidx_t      ldiff;

  P4EST_ASSERT (cc1 != NULL && *cc1 != NULL);
  P4EST_ASSERT (cc2 != NULL && *cc2 != NULL);

  P4EST_ASSERT ((*cc1)->contr.elem_size == sizeof (snodes_contr_t));
  P4EST_ASSERT ((*cc2)->contr.elem_size == sizeof (snodes_contr_t));

#ifdef P4EST_ENABLE_DEBUG
  /* we sort within the same owner process */
  o1 = (*cc1)->owner;
  o2 = (*cc2)->owner;
  P4EST_ASSERT (o1->rank == o2->rank);
#endif

  /* nodes are sorted according to their runid member */
  P4EST_ASSERT ((*cc1)->runid >= 0);
  P4EST_ASSERT ((*cc2)->runid >= 0);
  ldiff = (*cc1)->runid - (*cc2)->runid;
  return ldiff == 0 ? 0 : ldiff < 0 ? -1 : 1;
}


static void
push_sharer (snodes_meta_t * me, int *sindex, int rank)
{
  p4est_lnodes_rank_t *sharer;
  // p4est_lnodes_t     *ln = me->tm->lnodes;
  p4est_snodes_t *snodes = me->snodes;

  P4EST_ASSERT (me != NULL);
  P4EST_ASSERT (sindex != NULL);
  P4EST_ASSERT (0 <= rank && rank < me->mpisize);

  /* push empty sharer structure */
  *sindex = (int) snodes->sharers->elem_count;
  sharer = (p4est_lnodes_rank_t *) sc_array_push (snodes->sharers);
  memset (sharer, -1, sizeof (p4est_lnodes_rank_t));
  sharer->rank = rank;
  sc_array_init (&sharer->shared_nodes, sizeof (p4est_locidx_t));
}
#endif

static void
clean_construct (snodes_meta_t * me)
{
  snodes_cnode_t     *cnode;
  size_t              zz;

  for (zz = 0; zz < me->construct.elem_count; ++zz) {
    cnode = (snodes_cnode_t *) sc_array_index (&me->construct, zz);
    sc_array_reset (&cnode->contr);
  }
}

#ifdef  P4EST_ENABLE_MPI
static p4est_lnodes_rank_t *
peer_sharer (snodes_meta_t * me, int q)
{
  int                 pi;
  snodes_peer_t      *peer;

  P4EST_ASSERT (me != NULL);
  P4EST_ASSERT (me->ghost != NULL);
  P4EST_ASSERT (me->ghost_rank != NULL);
  P4EST_ASSERT (me->proc_peer != NULL);
  P4EST_ASSERT (0 <= q && q < me->mpisize);

  /* currently we do not store a peer for the local process */
  P4EST_ASSERT (q != me->mpirank);

  pi = me->proc_peer[q];
  P4EST_ASSERT (0 < pi && pi <= me->mpisize);
  peer = (snodes_peer_t *) sc_array_index_int (&me->peers, pi - 1);
  P4EST_ASSERT (peer->rank == q);
  return (p4est_lnodes_rank_t *) sc_array_index_int (me->snodes->sharers,
                                                     peer->sharind);
}

static snodes_peer_t *
peer_access (snodes_meta_t * me, int q)
{
  int                 pi;
  snodes_peer_t      *peer;

  P4EST_ASSERT (me != NULL);
  P4EST_ASSERT (me->ghost != NULL);
  P4EST_ASSERT (me->ghost_rank != NULL);
  P4EST_ASSERT (me->proc_peer != NULL);
  P4EST_ASSERT (0 <= q && q < me->mpisize);

  /* currently we do not store a peer for the local process */
  P4EST_ASSERT (q != me->mpirank);

  if ((pi = me->proc_peer[q]) == 0) {
    peer = (snodes_peer_t *) sc_array_push (&me->peers);
    me->proc_peer[q] = (int) me->peers.elem_count;
    peer->rank = q;
    peer->done = 0;
    peer->sharind = -1;
    peer->passive = 0;
    peer->lastadd = -1;
    peer->bufcount = 0;
    sc_array_init (&peer->sharedno, sizeof (p4est_locidx_t));
    sc_array_init (&peer->querypos, sizeof (p4est_locidx_t));
    sc_array_init (&peer->remosort, sizeof (snodes_cnode_t *));
  }
  else {
    P4EST_ASSERT (0 < pi && pi <= me->mpisize);
    peer = (snodes_peer_t *) sc_array_index_int (&me->peers, pi - 1);
    P4EST_ASSERT (peer->rank == q);
  }
  return peer;
}

#endif /* P4EST_ENABLE_MPI */

/** The local owner process will receive a query for a node number. */
static void
peer_add_reply (snodes_meta_t * me, snodes_peer_t * peer, p4est_locidx_t lni)
{
  P4EST_ASSERT (me != NULL);
  P4EST_ASSERT (peer != NULL);
  P4EST_ASSERT (peer->rank > me->mpirank);
  P4EST_ASSERT (peer->lastadd < lni);
  P4EST_ASSERT (0 <= lni && lni < (p4est_locidx_t) me->construct.elem_count);

  ++peer->bufcount;
  peer->lastadd = lni;
}

/** The local process queries a remote owner for its node number. */
static void
peer_add_query (snodes_meta_t * me, snodes_peer_t * peer,
                p4est_locidx_t lni, p4est_locidx_t epos)
{
#ifdef P4EST_ENABLE_DEBUG
  p4est_gloidx_t      gdiff;

  P4EST_ASSERT (me != NULL);
  P4EST_ASSERT (peer != NULL);
  P4EST_ASSERT (peer->rank < me->mpirank);
  P4EST_ASSERT (peer->lastadd < lni);
  P4EST_ASSERT (0 <= lni && lni < (p4est_locidx_t) me->construct.elem_count);
  gdiff = (me->p4est->global_first_quadrant[peer->rank + 1] -
           me->p4est->global_first_quadrant[peer->rank]);
  P4EST_ASSERT (0 <= epos && epos < SNODES_VNODES * (p4est_locidx_t) gdiff);
#endif

  ++peer->bufcount;
  *(p4est_locidx_t *) sc_array_push (&peer->querypos) = epos;
  *(p4est_locidx_t *) sc_array_push (&peer->sharedno) = peer->lastadd = lni;
}

static void
owned_query_reply (snodes_meta_t *me)
{
  snodes_cnode_t     *cnode, **ccn;
  snodes_contr_t     *owner;
#ifdef P4EST_ENABLE_MPI
  snodes_contr_t     *contr;
  snodes_peer_t      *peer;
  int                 withloc;
  size_t              zc, sic;
#endif
  size_t              zz, size;

  P4EST_ASSERT (me->num_owned == 0);
  P4EST_ASSERT (me->num_owned_shared == 0);
  P4EST_ASSERT (me->num_notowned_shared == 0);
  P4EST_ASSERT (me->num_shared == 0);

  size = me->construct.elem_count;
  for (zz = 0; zz < size; ++zz) {
    cnode = (snodes_cnode_t *) sc_array_index (&me->construct, zz);
    P4EST_ASSERT (cnode->runid == (p4est_locidx_t) zz);
    check_node (me, cnode->runid);

    owner = cnode->owner;
    if (owner->rank == me->mpirank) {
      ccn = (snodes_cnode_t **) sc_array_push(&me->ownsort);
      *ccn = cnode;
      ++me->num_owned;

#ifdef P4EST_ENABLE_MPI
      /* post replies for all queries to self */
      sic = cnode->contr.elem_count;
      for (zc = 0; zc < sic; ++zc) {
        contr = (snodes_contr_t *) sc_array_index (&cnode->contr, zc);
        if (contr->rank != me->mpirank) {
          P4EST_ASSERT (contr->rank > me->mpirank);
          peer = peer_access (me, contr->rank);
          peer_add_reply (me, peer, cnode->runid);
        }
        else {
          P4EST_ASSERT (owner == contr);
        }
      }
      if (sic > 1) {
        ++me->num_owned_shared;
        // ++me->num_shared;
      }
#endif
    }
    else {
#ifdef P4EST_ENABLE_MPI
      /* weed out remote-only nodes */
      withloc = 0;
      sic = cnode->contr.elem_count;
      for (zc = 0; zc < sic; ++zc) {
        contr = (snodes_contr_t *) sc_array_index (&cnode->contr, zc);
        if (contr->rank == me->mpirank) {
          withloc = 1;
          break;
        }
      }
      if (!withloc) {
        P4EST_ASSERT( 0 ); // Don't expect this to exists
        cnode->runid = -1;
        continue;
      }
      P4EST_ASSERT (owner->rank < me->mpirank);

      /* check for passively shared nodes */
      for (zc = 0; zc < sic; ++zc) {
        contr = (snodes_contr_t *) sc_array_index (&cnode->contr, zc);
        if (contr->rank != me->mpirank && contr->rank != owner->rank) {
          /* passively share a remotely owned node */
          P4EST_ASSERT (contr->rank > owner->rank);
          peer = peer_access (me, contr->rank);
          ++peer->passive;
        }
      }

      /* post query to remote owner */
      peer = peer_access (me, owner->rank);
      peer_add_query (me, peer, cnode->runid,
                      owner->le * SNODES_VNODES + owner->nodene);
      ccn = (snodes_cnode_t **) sc_array_push (&peer->remosort);
      *ccn = cnode;
      ++me->num_notowned_shared;
#else
      SC_ABORT_NOT_REACHED ();
#endif
    }


    /* the running id will be replaced by the owner's node number */
    cnode->runid = -1;
  }

  // me->num_notowned_shared = me->num_shared - me->num_owned_shared;
  me->num_shared = me->num_owned_shared + me->num_notowned_shared;
}

static void
post_query_reply (snodes_meta_t * me)
{
#ifndef P4EST_ENABLE_MPI
  P4EST_ASSERT (me != NULL);
  P4EST_ASSERT (me->num_all_shared == 0);
#else
  int                 mpiret;
  size_t              peer_count, iz;
  sc_MPI_Request     *preq;
  snodes_peer_t      *peer;

  /* explicitly do nothing without a ghost layer */
  peer_count = me->peers.elem_count;
  if (me->ghost == NULL || peer_count == 0) {
    P4EST_ASSERT (me->num_shared == 0);
    return;
  }
  P4EST_ASSERT (me->num_shared >= 0);

  /* go through peers (unsorted) and post messages */
  P4EST_ASSERT (me->emptypeers == 0);
  sc_array_resize (&me->peerreq, peer_count);
  for (iz = 0; iz < peer_count; ++iz) {
    peer = (snodes_peer_t *) sc_array_index (&me->peers, iz);
    preq = (sc_MPI_Request *) sc_array_index (&me->peerreq, iz);
    if (peer->bufcount == 0) {
      /* purely passive peers do not send messages */
      P4EST_ASSERT (peer->passive > 0);
      *preq = sc_MPI_REQUEST_NULL;
      ++me->emptypeers;
      continue;
    }
    if (peer->rank > me->mpirank) {
      /* expect query from higher rank */
      P4EST_ASSERT (peer->querypos.elem_count == 0);
      sc_array_resize (&peer->querypos, peer->bufcount);
      mpiret = sc_MPI_Irecv (sc_array_index (&peer->querypos, 0),
                             peer->bufcount, P4EST_MPI_LOCIDX, peer->rank,
                             P4EST_COMM_SNODES_QUERY, me->mpicomm, preq);
      SC_CHECK_MPI (mpiret);
      peer->done = 1;
    }
    else {
      /* address query to lower rank */
      P4EST_ASSERT (peer->rank < me->mpirank);
      P4EST_ASSERT (peer->bufcount ==
                    (p4est_locidx_t) peer->querypos.elem_count);
      mpiret =
        sc_MPI_Isend (sc_array_index (&peer->querypos, 0), peer->bufcount,
                      P4EST_MPI_LOCIDX, peer->rank, P4EST_COMM_SNODES_QUERY,
                      me->mpicomm, preq);
      SC_CHECK_MPI (mpiret);
      peer->done = 3;
    }
  }
#endif /* P4EST_ENABLE_MPI */
}

static void
wait_query_reply (snodes_meta_t * me)
{
#ifndef P4EST_ENABLE_MPI
  P4EST_ASSERT (me != NULL);
  P4EST_ASSERT (me->num_shared == 0);
#else
  int                 i, j;
  int                 mpiret;
  int                 nwalloc;
  int                 nwtotal;
  int                 nwaited;
  int                *waitind;
  sc_MPI_Request     *preq;
  p4est_locidx_t      lbc, lcl, lni, nonloc;
  p4est_locidx_t      epos, oind;
  p4est_gloidx_t      gof, gni;
  p4est_snodes_t     *snodes = me->snodes;
  snodes_cnode_t     *cnode;
  snodes_peer_t      *peer;

  /* explicitly do nothing without a ghost layer */
  nwalloc = (int) me->peers.elem_count;
  if (me->ghost == NULL || nwalloc == 0) {
    P4EST_ASSERT (me->num_shared == 0);
    return;
  }
  P4EST_ASSERT (me->num_shared >= 0);

  /* currently the local process does not count as peer */
  nwtotal = nwalloc - me->emptypeers;
  P4EST_ASSERT (nwtotal > 0);
  waitind = P4EST_ALLOC (int, nwalloc);
  while (nwtotal > 0) {
    mpiret = sc_MPI_Waitsome
      (nwalloc, (sc_MPI_Request *) sc_array_index (&me->peerreq, 0),
       &nwaited, waitind, sc_MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);
    SC_CHECK_ABORT (nwaited > 0, "Invalid count after MPI_Waitsome");
    for (i = 0; i < nwaited; ++i) {
      j = waitind[i];
      peer = (snodes_peer_t *) sc_array_index (&me->peers, j);
      P4EST_ASSERT (peer->rank != me->mpirank);
      preq = (sc_MPI_Request *) sc_array_index (&me->peerreq, j);
      P4EST_ASSERT (*preq == sc_MPI_REQUEST_NULL);
      if (peer->rank > me->mpirank) {
        P4EST_ASSERT (peer->shacumul == (me->num_shared - me->num_owned_shared));
        P4EST_ASSERT (peer->sharedno.elem_count == 0);
        if (peer->done == 1) {
          /* we have received a request and shall send a reply */
          lbc = peer->bufcount;
          for (lcl = 0; lcl < lbc; ++lcl) {
            epos = *(p4est_locidx_t *) sc_array_index (&peer->querypos, lcl);
#if 0
            P4EST_LDEBUGF ("Got %d gquad %d pos %d\n from %d\n", lcl,
                           epos / ln->vnodes, epos % ln->vnodes, peer->rank);
#endif
            P4EST_ASSERT (0 <= epos && epos < SNODES_VNODES * snodes->owned_count);
// #ifndef P4_TO_P8
//             P4EST_ASSERT (!alwaysowned[epos % ln->vnodes]);
// #else
//             /* extend this as further progress is made */
//             P4EST_ASSERT (epos % SNODES_VNODES < 8);
// #endif
            lni = snodes->element_nodes[epos];
            P4EST_ASSERT (0 <= lni && lni <
                          (p4est_locidx_t) me->construct.elem_count);
            cnode = (snodes_cnode_t *) sc_array_index (&me->construct, lni);
            oind = cnode->runid;
            P4EST_ASSERT (0 <= oind && oind < snodes->owned_count);

            /* send back the number of node owned by the local process */
            *(p4est_locidx_t *) sc_array_index (&peer->querypos, lcl) = oind;
          }
          mpiret = sc_MPI_Isend (sc_array_index (&peer->querypos, 0),
                                 peer->bufcount, P4EST_MPI_LOCIDX, peer->rank,
                                 P4EST_COMM_SNODES_REPLY, me->mpicomm, preq);
          SC_CHECK_MPI (mpiret);
          peer->done = 2;
        }
        else {
          /* our reply has been received */
          P4EST_ASSERT (peer->done == 2);
          peer->done = 0;
          --nwtotal;
        }
      }
      else {
        P4EST_ASSERT (peer->rank < me->mpirank);
        if (peer->done == 3) {
          /* our request has been sent and we await the reply */
          mpiret = sc_MPI_Irecv (sc_array_index (&peer->querypos, 0),
                                 peer->bufcount, P4EST_MPI_LOCIDX, peer->rank,
                                 P4EST_COMM_SNODES_REPLY, me->mpicomm, preq);
          SC_CHECK_MPI (mpiret);
          peer->done = 4;
        }
        else {
          /* process owner's node numbers in reply received */
          P4EST_ASSERT (peer->done == 4);
          lbc = peer->bufcount;
          for (lcl = 0; lcl < lbc; ++lcl) {
            oind = *(p4est_locidx_t *) sc_array_index (&peer->querypos, lcl);
            lni = *(p4est_locidx_t *) sc_array_index (&peer->sharedno, lcl);
            cnode = (snodes_cnode_t *) sc_array_index (&me->construct, lni);
            P4EST_ASSERT (cnode->owner->rank == peer->rank);
            P4EST_ASSERT (0 <= oind &&
                          oind < snodes->global_owned_count[peer->rank]);
            cnode->runid = oind;
          }
          sc_array_sort (&peer->remosort, rnode_compare);

          /* store shared node's global index */
          gof = me->goffset[peer->rank];
          for (lcl = 0; lcl < lbc; ++lcl) {
            cnode =
              *(snodes_cnode_t **) sc_array_index (&peer->remosort, lcl);
            P4EST_ASSERT (cnode->owner->rank == peer->rank);
            P4EST_ASSERT (lcl <= cnode->runid);
            nonloc = peer->shacumul + lcl;
            P4EST_ASSERT (nonloc < (me->num_shared - me->num_owned_shared));
            gni = gof + cnode->runid;
            P4EST_ASSERT (me->goffset[peer->rank] <= gni &&
                          gni < me->goffset[peer->rank + 1]);
            snodes->nonlocal_nodes[nonloc] = gni;

            /* now the runid of each node is the local node number */
            cnode->runid = me->num_owned + nonloc;
          }
          peer->done = 0;
          --nwtotal;
        }
      }
    }
  }
  P4EST_FREE (waitind);
#ifdef P4EST_ENABLE_DEBUG
  gof = -1;
  for (lcl = 0; lcl < me->num_notowned_shared; ++lcl) {
    gni = snodes->nonlocal_nodes[lcl];
    P4EST_ASSERT (0 <= gni && gni < me->goffset[me->mpisize]);
    P4EST_ASSERT (gni < me->goffset[me->mpirank] ||
                  gni >= me->goffset[me->mpirank + 1]);
    P4EST_ASSERT (gni > gof);
    gof = gni;
  }
#endif /* P4EST_ENABLE_DEBUG */
#endif /* P4EST_ENABLE_MPI */
}

static void
populate_sharers (snodes_meta_t * me)
{
#ifdef P4EST_ENABLE_MPI
  int                 i;
  int                 num_peers;
  size_t              zz;
  p4est_locidx_t      lni;
  p4est_locidx_t      lbc, lcl;
  p4est_snodes_t     *snodes = me->snodes;
  p4est_lnodes_rank_t *sharer, *locshare;
  snodes_peer_t      *tp;
  snodes_cnode_t     *cnode;
  snodes_contr_t     *contr;

  /* populate sharers array */
  num_peers = (int) me->peers.elem_count;
  if (me->ghost == NULL || num_peers == 0) {
    P4EST_ASSERT (me->num_shared == 0);
    return;
  }
  P4EST_ASSERT (me->num_shared >= 0);
  P4EST_ASSERT (num_peers + 1 == (p4est_locidx_t) snodes->sharers->elem_count);
  locshare =
    (p4est_lnodes_rank_t *) sc_array_index_int (snodes->sharers, me->locsharer);

  /* first iterate through owned nodes in order */
  lbc = (p4est_locidx_t) me->ownsort.elem_count;
  for (lcl = 0; lcl < lbc; ++lcl) {
    cnode = *(snodes_cnode_t **) sc_array_index (&me->ownsort, lcl);
    P4EST_ASSERT (cnode->owner->rank == me->mpirank);
    P4EST_ASSERT (lcl == cnode->runid);
    if (cnode->contr.elem_count == 1) {
      /* this node is purely local */
      continue;
    }

    /* this node has sharers: iterate through all of them */
    for (zz = 0; zz < cnode->contr.elem_count; ++zz) {
      contr = (snodes_contr_t *) sc_array_index (&cnode->contr, zz);
      if (contr->rank == me->mpirank) {
        /* local process is owner */
        P4EST_ASSERT (cnode->owner == contr);
        sharer = locshare;
      }
      else {
        /* remote process is sharer */
        sharer = peer_sharer (me, contr->rank);
      }
      P4EST_ASSERT (sharer->rank == contr->rank);
      *(p4est_locidx_t *) sc_array_push (&sharer->shared_nodes) = lcl;
    }
  }
  P4EST_ASSERT (me->num_owned == (p4est_locidx_t) me->ownsort.elem_count);
  P4EST_ASSERT (me->num_owned_shared ==
                (p4est_locidx_t) locshare->shared_nodes.elem_count);

  /* determine the sharer offset and count variables */
  locshare->shared_mine_offset = locshare->owned_offset = 0;
  locshare->shared_mine_count = me->num_owned_shared;
  locshare->owned_count = me->num_owned;
  for (i = 0; i < num_peers; ++i) {
    tp = *(snodes_peer_t **) sc_array_index_int (&me->sortp, i);
    sharer =
      (p4est_lnodes_rank_t *) sc_array_index_int (snodes->sharers, tp->sharind);
    P4EST_ASSERT (tp->rank == sharer->rank);
    sharer->shared_mine_offset = 0;
    sharer->shared_mine_count =
      (p4est_locidx_t) sharer->shared_nodes.elem_count;
    sharer->owned_offset = me->num_owned + tp->shacumul;
    if (tp->rank < me->mpirank) {
      P4EST_ASSERT (tp->bufcount > 0 || tp->passive > 0);
      sharer->owned_count = tp->bufcount;
    }
    else {
      P4EST_ASSERT (tp->rank > me->mpirank);
      sharer->owned_count = 0;
    }
  }

  /* iterate through the remote local nodes in order */
  lni = me->num_owned;
  for (i = 0; i < num_peers; ++i) {
    tp = *(snodes_peer_t **) sc_array_index_int (&me->sortp, i);
    if (tp->rank < me->mpirank) {
      lbc = tp->bufcount;
      P4EST_ASSERT (lbc == (p4est_locidx_t) tp->remosort.elem_count);
      for (lcl = 0; lcl < lbc; ++lcl, ++lni) {
        cnode = *(snodes_cnode_t **) sc_array_index (&tp->remosort, lcl);
        P4EST_ASSERT (cnode->owner->rank == tp->rank);
        P4EST_ASSERT (cnode->runid == lni);

        /* this node has sharers: iterate through all of them */
        P4EST_ASSERT (cnode->contr.elem_count > 1);
        for (zz = 0; zz < cnode->contr.elem_count; ++zz) {
          contr = (snodes_contr_t *) sc_array_index (&cnode->contr, zz);
          if (contr->rank == me->mpirank) {
            /* local process is sharer */
            P4EST_ASSERT (cnode->owner != contr);
            sharer = locshare;
          }
          else {
            /* remote process is owner or not */
            sharer = peer_sharer (me, contr->rank);
          }
        }
        P4EST_ASSERT (sharer->rank == contr->rank);
        *(p4est_locidx_t *) sc_array_push (&sharer->shared_nodes) = lni;
      }
    }
  }
  P4EST_ASSERT (lni == snodes->num_local_nodes);
#endif /* P4EST_ENABLE_MPI */
}

static void
sort_allgather (snodes_meta_t * me)
{
  snodes_cnode_t    **ccn;
  p4est_snodes_t    *snodes = me->snodes;
  // p4est_locidx_t      lel, le, lc;
  p4est_locidx_t     *localboth, lb[2];
  p4est_gloidx_t      gc;
  const int           s = me->mpisize;
  int                 mpiret;
#ifndef P4_TO_P8
  int                 cind, lookup;
#endif
  int                 q;
  size_t              zz;

  /* sort local node list */
  sc_array_sort (&me->ownsort, cnode_compare);
  for (zz = 0; zz < me->ownsort.elem_count; ++zz) {
    ccn = (snodes_cnode_t **) sc_array_index (&me->ownsort, zz);
    (*ccn)->runid = zz;
  }

  /* share owned count */
  snodes->owned_count = me->num_owned;
  snodes->num_local_nodes = me->num_owned + me->num_notowned_shared;
  snodes->nonlocal_nodes = P4EST_ALLOC (p4est_gloidx_t, me->num_notowned_shared);;
  snodes->global_owned_count = P4EST_ALLOC (p4est_locidx_t, s);
  lb[0] = snodes->owned_count;

  // TODO:
  lb[1] = 0; // Number of owned local triangles

  //
//   /* establish local triangle count */
//   lel = snodes->num_local_elements;
//   lc = me->local_element_offset[0] = 0;
//   for (le = 0; le < lel; ++le) {
// #ifndef P4_TO_P8
//     cind = config_cind (me->configuration[le]);
//     lookup = p4est_tnodes_config_lookup[cind];
//     P4EST_ASSERT (0 <= lookup && lookup < 6);
//     lc = tm->local_element_offset[le + 1] =
//       lc + p4est_tnodes_lookup_counts[lookup][2];
// #else
//     /* extend this as further progress is made */
//     lc = tm->local_element_offset[le + 1] = lc + 0;
// #endif
//   }
//   lb[1] = me->num_triangles = lc;

  /* parallel sharing of owned node and element counts */
  localboth = P4EST_ALLOC (p4est_locidx_t, 2 * me->mpisize);
  mpiret = sc_MPI_Allgather (lb, 2, P4EST_MPI_LOCIDX,
                             localboth, 2, P4EST_MPI_LOCIDX,
                             me->p4est->mpicomm);
  SC_CHECK_MPI (mpiret);
  me->goffset = P4EST_ALLOC (p4est_gloidx_t, s + 1);
  gc = me->goffset[0] = 0;
  // P4EST_ASSERT (me->num_global_triangles == 0);
  for (q = 0; q < s; ++q) {
    gc = me->goffset[q + 1] =
      gc + (snodes->global_owned_count[q] = localboth[2 * q + 0]);
    // if (q == me->mpirank) {
    //   tm->global_toffset = me->num_global_triangles;
    // }
    // me->num_global_triangles += (tm->local_tcount[q] = localboth[2 * q + 1]);
  }
  snodes->global_offset = me->goffset[me->mpirank];
  P4EST_FREE (localboth);
}

static void
sort_peers (snodes_meta_t * me)
{
#ifndef P4EST_ENABLE_MPI
  P4EST_ASSERT (me != NULL);
  P4EST_ASSERT (me->num_all_shared == 0);
#else
  int                 i;
  int                 num_peers;
  p4est_locidx_t      nonlofs;
  snodes_peer_t      *tp;

  /* explicitly do nothing without a ghost layer */
  num_peers = (int) me->peers.elem_count;
  if (me->ghost == NULL || num_peers == 0) {
    P4EST_ASSERT (me->num_shared == 0);
    return;
  }
  P4EST_ASSERT (me->num_shared > 0);

  /* make it possible to iterate through peers in rank order */
  sc_array_resize (&me->sortp, num_peers);
  for (i = 0; i < num_peers; ++i) {
    *(snodes_peer_t **) sc_array_index_int (&me->sortp, i) =
      (snodes_peer_t *) sc_array_index_int (&me->peers, i);
  }
  sc_array_sort (&me->sortp, peer_compare);
  nonlofs = 0;
  for (i = 0; i < num_peers; ++i) {
    tp = *(snodes_peer_t **) sc_array_index_int (&me->sortp, i);
    tp->shacumul = nonlofs;
    if (tp->rank < me->mpirank) {
      nonlofs += tp->bufcount;
    }
  }

  P4EST_ASSERT (nonlofs == me->num_notowned_shared);

  /* initialize sharers array */
  for (i = 0; i < num_peers; ++i) {
    tp = *(snodes_peer_t **) sc_array_index_int (&me->sortp, i);
    P4EST_ASSERT (tp->rank != me->mpirank);
    if (tp->rank > me->mpirank) {
      break;
    }
    push_sharer (me, &tp->sharind, tp->rank);
  }
  push_sharer (me, &me->locsharer, me->mpirank);
  for (; i < num_peers; ++i) {
    tp = *(snodes_peer_t **) sc_array_index_int (&me->sortp, i);
    P4EST_ASSERT (tp->rank > me->mpirank);
    push_sharer (me, &tp->sharind, tp->rank);
  }
  P4EST_ASSERT (num_peers + 1 == (int) me->snodes->sharers->elem_count);
  P4EST_ASSERT (me->locsharer >= 0);
#endif /* P4EST_ENABLE_MPI */
}

p4est_snodes_t *
p4est_snodes_new (
    p4est_t *p4est, p4est_ghost_t *ghost)
{
  p4est_snodes_t *snodes;
  snodes_meta_t snmeta, *me = &snmeta;

  int mpisize;
  p4est_locidx_t lel;
#ifdef P4EST_ENABLE_MPI
  size_t              nz, zi;
  snodes_peer_t      *peer;
#endif

  memset (me, 0, sizeof (*me));
  //
  snodes = me->snodes = P4EST_ALLOC_ZERO(p4est_snodes_t, 1);
  lel = p4est->local_num_quadrants;
  snodes->num_local_elements = lel;
  // me->snodes->num_local_elements = lel;
  me->chilev = P4EST_ALLOC_ZERO(uint8_t, lel);
  snodes->element_nodes = me->element_nodes =
                P4EST_ALLOC(p4est_locidx_t, lel * SNODES_VNODES);
  memset (me->element_nodes, -1, lel * SNODES_VNODES * sizeof(p4est_locidx_t));

  me->mpicomm = p4est->mpicomm;
  me->mpirank = p4est->mpirank;
  me->mpisize = mpisize = p4est->mpisize;
  me->p4est = p4est;
  me->ghost = ghost;

  p4est_locidx_t lg, ng;
  int q;
  if (ghost != NULL) {
    P4EST_ASSERT (ghost->proc_offsets[0] == 0);
    P4EST_ASSERT (ghost->proc_offsets[me->mpisize] ==
                  (p4est_locidx_t) ghost->ghosts.elem_count);
    me->ghost_rank = P4EST_ALLOC (int, ghost->ghosts.elem_count);
    lg = 0;
    for (q = 0; q < me->mpisize; ++q) {
      ng = ghost->proc_offsets[q + 1];
      for (; lg < ng; ++lg) {
        me->ghost_rank[lg] = q;
      }
    }
    P4EST_ASSERT (lg == (p4est_locidx_t) ghost->ghosts.elem_count);
#ifdef P4EST_ENABLE_MPI
    me->proc_peer = P4EST_ALLOC_ZERO (int, me->mpisize);
    sc_array_init (&me->sortp, sizeof (snodes_peer_t *));
    sc_array_init (&me->peers, sizeof (snodes_peer_t));
    sc_array_init (&me->peerreq, sizeof (sc_MPI_Request));
#endif
  }

  sc_array_init (&me->construct, sizeof (snodes_peer_t));
  sc_array_init (&me->ownsort, sizeof (sc_MPI_Request));

  snodes->mpicomm = p4est->mpicomm;
  snodes->sharers = sc_array_new(sizeof(p4est_lnodes_rank_t));
  snodes->degree = 0;
  snodes->vnodes = SNODES_VNODES;

  // Find all local nodes
  p4est_iterate (p4est, ghost, me,
      iter_volume, iter_face,
#ifdef P4_TO_P8
      iter_edge,
#endif
     iter_corner);

  P4EST_ASSERT(me->lenum == lel);

  owned_query_reply (me);
  P4EST_INFOF ("p4est_snodes_new: nodes owned %ld shared %ld\n",
               (long) me->num_owned, (long) me->num_notowned_shared);

  /* post first round of messages */
  post_query_reply (me);

  /* sort local nodes and allgather owned counts */
  sort_allgather (me);
  // P4EST_INFOF ("p4est_snodes_new: triangles owned %ld\n",
  //              (long) me->num_triangles);
  // P4EST_GLOBAL_PRODUCTIONF
  //   ("p4est_snodes_new: global triangles %lld nodes %lld\n",
  //    (long long) me->num_global_triangles, (long long) me->goffset[mpisize]);
  P4EST_GLOBAL_PRODUCTIONF
    ("p4est_snodes_new: nodes %lld\n", (long long) me->goffset[mpisize]);

  /* sort peers by process */
  sort_peers (me);

  /* receive query messages and send replies */
  wait_query_reply (me);

  // TODO: Check that the indices in sharers agree with those in element_nodes

  /* finalize element node assignment */
// #ifndef P4_TO_P8
//   assign_element_nodes (me);
// #endif /* P4_TO_P8 */

  /* finalize sharer information */
  populate_sharers (me);

  /* cleanup */
  P4EST_FREE (me->goffset);
  if (me->ghost) {
#ifdef P4EST_ENABLE_MPI
    nz = me->peers.elem_count;
    for (zi = 0; zi < nz; ++zi) {
      peer = (snodes_peer_t *) sc_array_index (&me->peers, zi);
      P4EST_ASSERT (!peer->done);
      sc_array_reset (&peer->sharedno);
      sc_array_reset (&peer->querypos);
      sc_array_reset (&peer->remosort);
    }

    P4EST_FREE(me->proc_peer);
    sc_array_reset(&me->sortp);
    sc_array_reset(&me->peers);
    sc_array_reset(&me->peerreq);
#endif
    P4EST_FREE(me->ghost_rank);
  }

  clean_construct(me);
  sc_array_reset(&me->construct);
  sc_array_reset(&me->ownsort);
  P4EST_FREE(me->chilev);

  return snodes;
}

void
p4est_snodes_destory(p4est_snodes_t *snodes)
{
  size_t count, zz;
  p4est_lnodes_rank_t *lrank;

  P4EST_ASSERT(snodes);
  P4EST_ASSERT(snodes->element_nodes);

  P4EST_FREE (snodes->element_nodes);

  if (snodes->nonlocal_nodes)
    P4EST_FREE (snodes->nonlocal_nodes);

  if (snodes->global_owned_count)
    P4EST_FREE (snodes->global_owned_count);


  count = snodes->sharers->elem_count;
  for (zz = 0; zz < count; zz++) {
    lrank = p4est_lnodes_rank_array_index (snodes->sharers, zz);
    sc_array_reset (&(lrank->shared_nodes));
  }
  sc_array_destroy (snodes->sharers);

  P4EST_FREE (snodes);
}

void
p4est_geometry_coordinates_snodes
  (p4est_t *p4est, p4est_snodes_t *snodes,
   p4est_geometry_t *geom,
   sc_array_t *coordinates)
{
  p4est_topidx_t tt;
  size_t qz, quad_count;
  size_t k;

  int quad_len;
  static const double irlen = 1. / P4EST_ROOT_LEN;

  p4est_locidx_t lni;
  p4est_locidx_t *enodes;

  double abc[3], *xyz, *ret;

  p4est_tree_t *tree;
  p4est_quadrant_t *quad;
  char *mask;

  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (p4est->connectivity != NULL);
  P4EST_ASSERT (geom == NULL || geom->X != NULL);
  P4EST_ASSERT (snodes != NULL);
  P4EST_ASSERT (coordinates != NULL &&
                coordinates->elem_size == 3 * sizeof (double));

  sc_array_resize(coordinates, snodes->num_local_nodes);

  mask = P4EST_ALLOC_ZERO(char, snodes->num_local_nodes);

  enodes = snodes->element_nodes;

  for (tt = p4est->first_local_tree; tt <= p4est->last_local_tree; ++tt) {

    tree = p4est_tree_array_index(p4est->trees, tt);
    quad_count = tree->quadrants.elem_count;
    for (qz = 0; qz < quad_count; ++qz, enodes += SNODES_VNODES) {
      quad = p4est_quadrant_array_index(&tree->quadrants, qz);

      for (k = 0; k < SNODES_VNODES; ++k) {
        lni = enodes[k];
        P4EST_ASSERT(lni < snodes->num_local_nodes);

        if (0 <= lni && !mask[lni]) {
          mask[lni] = 1;
          xyz = sc_array_index(coordinates, lni);

          ret = (geom == NULL) ? xyz : abc;
          quad_len = P4EST_QUADRANT_LEN(quad->level) / 2.0;
          ret[0] = (quad->x + elem_coords_n[k][0] * quad_len) * irlen;
          ret[1] = (quad->y + elem_coords_n[k][1] * quad_len) * irlen;
#ifndef P4_TO_P8
          ret[2] = 0.;
#else
          ret[2] = (quad->z + elem_coords_n[k][2] * quad_len) * irlen;
#endif
          if (geom != NULL) {
            /* transform tree reference into geometry image */
            geom->X (geom, tt, abc, xyz);
          }

        } // if mask
      } // for |vnodes|
    } // for |quads|
  } //  for |local_trees|

  P4EST_FREE(mask);

}


p4est_simplex_mesh_t *
p4est_new_simplex_mesh(
    p4est_t *p4est,
    p4est_geometry_t *geom,
    p4est_ghost_t *ghost)
{
  size_t zz, fi;
  p4est_topidx_t t;
  p4est_tree_t *tree;
  p4est_locidx_t le;
  p4est_locidx_t *le_nodes;
  sc_array_t *quadrants;
  p4est_quadrant_t qid, *quad;

  p4est_locidx_t cf, cc, c0, c1;

#ifdef P4_TO_P8
  size_t ei;
  int edge;
  // int face_has_hanging_edge;

  p4est_locidx_t ce, c2, c3;
#endif

  // element_node_t *elem_node;
  p4est_locidx_t *simplex;


  p4est_snodes_t *snodes;
  p4est_simplex_mesh_t *smesh;

  smesh = P4EST_ALLOC_ZERO(p4est_simplex_mesh_t, 1);

  snodes = smesh->snodes = p4est_snodes_new(p4est, ghost);

  // Compute coordinates
  smesh->coordinates = sc_array_new(3 * sizeof(double));
  p4est_geometry_coordinates_snodes
    (p4est, snodes, geom, smesh->coordinates);

  // Simplices
  smesh->simplicies = sc_array_new((P4EST_DIM + 1)*sizeof(p4est_locidx_t));

  // return smesh;

  // We iterate through the local quads to make the simplices
  le = 0;
  for (t=p4est->first_local_tree;
          t<=p4est->last_local_tree; ++t) {

    tree = p4est_tree_array_index(p4est->trees, t);

    quadrants = &tree->quadrants;
    for (zz = 0; zz < quadrants->elem_count; ++zz, ++le) {
      quad = p4est_quadrant_array_index(quadrants, zz);

      qid = *quad;
      qid.p.which_tree = t;

      le_nodes = &snodes->element_nodes[le*SNODES_VNODES];

      // elem_node = sc_array_index(me->sa_element_nodes, elemid);
      // cc = elem_node->volume_node;
      cc = le_nodes[n_center];

      for (fi = 0; fi < P4EST_FACES; ++fi) {
        // cf = elem_node->face_nodes[fi];
        cf = le_nodes[n_face[fi]];

        if (0 <= cf)
        {
#ifndef P4_TO_P8
          // c0 = elem_node->corner_nodes[ p4est_face_corners[fi][0] ];
          // c1 = elem_node->corner_nodes[ p4est_face_corners[fi][1] ];

          c0 = le_nodes[ n_corner[p4est_face_corners[fi][0]] ];
          c1 = le_nodes[ n_corner[p4est_face_corners[fi][1]] ];

          simplex = sc_array_push(smesh->simplicies);
          simplex[0] = cc;
          simplex[1] = c0;
          simplex[2] = cf;

          simplex = sc_array_push(smesh->simplicies);
          simplex[0] = cc;
          simplex[1] = cf;
          simplex[2] = c1;
#else
          for (ei = 0; ei < 4; ++ei) {
            edge = p8est_face_edges[fi][ei];

            // c0 = elem_node->corner_nodes[ p8est_edge_corners[edge][0] ];
            // c1 = elem_node->corner_nodes[ p8est_edge_corners[edge][1] ];
            //
            // ce = elem_node->edge_nodes[edge];

            c0 = le_nodes[ n_corner[p8est_edge_corners[edge][0]] ];
            c1 = le_nodes[ n_corner[p8est_edge_corners[edge][1]] ];

            ce = le_nodes[ n_edge[edge] ];

            if (ce < 0)
            {
              simplex = sc_array_push(smesh->simplicies);
              simplex[0] = cc;
              simplex[1] = cf;
              simplex[2] = c0;
              simplex[3] = c1;
            }
            else
            {
              // edge center
              simplex = sc_array_push(smesh->simplicies);
              simplex[0] = cc;
              simplex[1] = cf;
              simplex[2] = c0;
              simplex[3] = ce;

              simplex = sc_array_push(smesh->simplicies);
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
          // c0 = elem_node->corner_nodes[ p4est_face_corners[fi][0] ];
          // c1 = elem_node->corner_nodes[ p4est_face_corners[fi][1] ];

          c0 = le_nodes[ n_corner[p4est_face_corners[fi][0]] ];
          c1 = le_nodes[ n_corner[p4est_face_corners[fi][1]] ];

#ifndef P4_TO_P8

          simplex = sc_array_push(smesh->simplicies);
          simplex[0] = cc;
          simplex[1] = c0;
          simplex[2] = c1;
#else
          // c2 = elem_node->corner_nodes[ p4est_face_corners[fi][2] ];
          // c3 = elem_node->corner_nodes[ p4est_face_corners[fi][3] ];

          c2 = le_nodes[ n_corner[p4est_face_corners[fi][2]] ];
          c3 = le_nodes[ n_corner[p4est_face_corners[fi][3]] ];

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

          simplex = sc_array_push(smesh->simplicies);
          simplex[0] = cc;
          simplex[1] = c0;
          simplex[2] = c1;
          simplex[3] = c3;

          simplex = sc_array_push(smesh->simplicies);
          simplex[0] = cc;
          simplex[1] = c3;
          simplex[2] = c2;
          simplex[3] = c0;

#endif
        } // if fc == 0

      }  // for (fi = 0..P4EST_FACES)
    } // for (zz = 0..|quadrants|)
  } // for (t = first-..last-localtrees)


  // sc_array_destroy(me.sa_element_nodes);
  // sc_array_destroy(me.tree_offsets);

  return smesh;
}


void
p4est_simplex_mesh_destroy(
    p4est_simplex_mesh_t *smesh)
{
  P4EST_ASSERT(smesh);

  if (smesh->snodes)
    p4est_snodes_destory(smesh->snodes);

  if (smesh->coordinates)
    sc_array_destroy(smesh->coordinates);

  if (smesh->simplicies)
    sc_array_destroy(smesh->simplicies);

  P4EST_FREE(smesh);
}


// Some utility functions for computing the coordinates of the nodes

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
    // size_t num_vertices = mesh->vertices->elem_count;
    size_t num_vertices = mesh->coordinates->elem_count;
    fprintf(file, "$Nodes\n%zu\n", num_vertices);
    double *vertices = (double *) mesh->coordinates->array;
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
    size_t num_vertices = mesh->coordinates->elem_count;
    fprintf(file, "POINTS %zu double\n", num_vertices);
    double *vertices = (double *) mesh->coordinates->array;
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
