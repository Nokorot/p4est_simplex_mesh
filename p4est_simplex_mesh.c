
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
typedef struct p4est_iter_face_side_full
{
  int8_t              is_ghost; /**< boolean: local (0) or ghost */
  p4est_quadrant_t   *quad;     /**< the actual quadrant */
  p4est_locidx_t      quadid;   /**< index in tree or ghost array */

}
p4est_iter_face_side_full_t;

/** Shortcut to access hanging face information */
typedef struct p4est_iter_face_side_hanging
{
  int8_t              is_ghost[2];      /**< boolean: local (0) or ghost */
  p4est_quadrant_t   *quad[2];          /**< the actual quadrant */
  p4est_locidx_t      quadid[2];        /**< index in tree or ghost array */
}
p4est_iter_face_side_hanging_t;


#ifdef P4_TO_P8
/** Shortcut to access full face information */
typedef struct p8est_iter_face_side_full
{
  int8_t              is_ghost; /**< boolean: local (0) or ghost */
  p8est_quadrant_t   *quad;     /**< the actual quadrant */
  p4est_locidx_t      quadid;   /**< index in tree or ghost array */

}
p8est_iter_face_side_full_t;

/** Shortcut to access hanging face information */
typedef struct p8est_iter_face_side_hanging
{
  int8_t              is_ghost[4];      /**< boolean: local (0) or ghost */
  p8est_quadrant_t   *quad[4];          /**< the actual quadrant */
  p4est_locidx_t      quadid[4];        /**< index in tree or ghost array */
}
p8est_iter_face_side_hanging_t;

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


static int P4EST_COMM_SIMPLEX_PASS = 50;


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


#else
/* face midpoints  */
static const int n_face[4] = { 1, 2, 3, 4 };

/* cube corners */
static const int n_corner[4] = { 5, 6, 7, 8 };
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
  uint8_t *chilev;


  p4est_locidx_t lenum;              // Element number
  p4est_locidx_t num_owned;
  p4est_locidx_t num_owned_shared;  /**< Nodes we both own and share */
  p4est_locidx_t num_shared;        /**< Nodes we share, owned or not */


  sc_array_t construct;
  p4est_locidx_t *element_nodes;

  // element configuration might be useful


  // sc_array_t *tree_offsets;
  // sc_array_t *global_owned_count;
  // sc_array_t *owners;

  // sc_array_t **shared_element_nodes;

  // sc_array_t *sa_element_nodes;
  // sc_array_t *vertices;

  p4est_locidx_t num_local_elements;

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
  snodes_contr_t *owner;         /**< Codimension of node. */
  // snode_connect_type_t bcon;    /**< Owning contributor. */
  sc_array_t contr;             /**< Contributing processes. */
} snodes_cnode_t;

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
  size_t siz, zz;

  /* basic checks */
  P4EST_ASSERT (me != NULL);
  P4EST_ASSERT(0 <= le && le < me->num_local_elements);
  P4EST_ASSERT(0 <= nodene && nodene < SNODES_VNODES);


  /* a new node is to be created or an existing one is passed in */
  if (lni == NULL) {
    lnis = -1;
    lni = &lnis;
  }

  /* if a node is already registered we don't need to do anything */
  if (me->element_nodes[le * SNODES_VNODES + nodene] != -1) {
    P4EST_ASSERT (*lni == -1 ||
                  *lni == me->element_nodes[le * SNODES_VNODES + nodene]);

    *lni = me->element_nodes[le * SNODES_VNODES + nodene];

    return;
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

  if (*lni == -1) {
    /* create a new node with one instance */
    *lni = (p4est_locidx_t) me->construct.elem_count;
    cnode = (snodes_cnode_t *) sc_array_push (&me->construct);
    cnode->runid = *lni;

    // cnode->bcon = bcon;
    cnode->owner = NULL;
    sc_array_init (&cnode->contr, sizeof (snodes_contr_t));
  }
  else {
    /* create a new instance of an existing node */
    cnode = (snodes_cnode_t *) sc_array_index (&me->construct, (size_t) *lni);

    P4EST_ASSERT (cnode->runid == *lni);
    // P4EST_ASSERT (cnode->bcon == bcon);
    P4EST_ASSERT (cnode->contr.elem_size == sizeof (snodes_contr_t));
    P4EST_ASSERT (cnode->contr.elem_count > 0);
    P4EST_ASSERT (cnode->owner != NULL);
    // TODO: check_node (me, *lni);
  }

  /* assign node to the local element position */
  if (rank == me->mpirank) {
    P4EST_ASSERT (me->element_nodes[le * SNODES_VNODES + nodene] == -1);
    me->element_nodes[le * SNODES_VNODES + nodene] = *lni;
  }

  /* iterate through instances to find matching process */
  siz = cnode->contr.elem_count;
  P4EST_ASSERT (siz == 0 || cnode->owner != NULL);
  for (zz = 0; zz < siz; ++zz) {
    contr = (snodes_contr_t *) sc_array_index (&cnode->contr, zz);
    P4EST_ASSERT (cnode->owner->rank <= contr->rank);
    if (contr->rank == rank) {
      /* rank is found and we remember the smallest node position */
      if (le < contr->le || (le == contr->le && nodene < contr->nodene)) {
        contr->nodene = nodene;
        contr->le = le;
      }
      // TODO: check_node (me, *lni);
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
    siz = cnode->contr.elem_count;
    for (zz = 1; zz < siz - 1; ++zz) {
      contr = (snodes_contr_t *) sc_array_index (&cnode->contr, zz);
      if (contr->rank < cnode->owner->rank) {
        cnode->owner = contr;
      }
    }
  }
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


p4est_snodes_t *
p4est_snodes_new (
    p4est_t *p4est, p4est_ghost_t *ghost)
{
  p4est_snodes_t *snodes;
  snodes_meta_t snmeta, *me = &snmeta;

  p4est_locidx_t lel;

  memset (me, 0, sizeof (*me));
  sc_array_init(&me->construct, sizeof(snodes_cnode_t));
  //
  lel = p4est->local_num_quadrants;
  me->num_local_elements = lel;
  me->chilev = P4EST_ALLOC_ZERO(uint8_t, lel);
  me->element_nodes = P4EST_ALLOC(p4est_locidx_t, lel * SNODES_VNODES);
  memset (me->element_nodes, -1, lel * SNODES_VNODES * sizeof(p4est_locidx_t));

  me->mpicomm = p4est->mpicomm;
  me->mpirank = p4est->mpirank;
  me->mpisize = p4est->mpisize;
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
// #ifdef P4EST_ENABLE_MPI
//     me->proc_peer = P4EST_ALLOC_ZERO (int, s);
//     sc_array_init (&me->sortp, sizeof (tnodes_peer_t *));
//     sc_array_init (&me->peers, sizeof (tnodes_peer_t));
//     sc_array_init (&me->pereq, sizeof (sc_MPI_Request));
// #endif
  }


  // Find all local nodes
  p4est_iterate(p4est, ghost, me,
      iter_volume, iter_face,
#ifdef P4_TO_P8
      iter_edge,
#endif
     iter_corner);


  snodes = P4EST_ALLOC_ZERO(p4est_snodes_t, 1);
  snodes->mpicomm = me->mpicomm;

  snodes->num_local_nodes = me->construct.elem_count;
  // snodes->owned_count;
  // snodes->global_offset;
  // snodes->nonlocal_nodes;
  // snodes->sharers;
  // snodes->global_owned_count;
  snodes->num_local_elements = p4est->local_num_quadrants;
  snodes->element_nodes = me->element_nodes;

  if (me->ghost) {
    P4EST_FREE(me->ghost_rank);
  }

  clean_construct(me);
  sc_array_reset(&me->construct);
  P4EST_FREE(me->chilev);

  return snodes;

#if 0
  // *** We sort the local nodes according to process rank, with the
  // current process nodes first.
  //
  // We thus reorder the vertices array and remap element_nodes

  p4est_locidx_t num_local_nodes, num_owned_nodes;

  p4est_locidx_t *proc_offsets, offset;
  size_t *owned;
  int *owner;
  size_t zz, elem, cc_off;
  int vnodes;

  p4est_locidx_t *new_local_nodes, *proc_node_count, new_cc;
  p4est_locidx_t *elem_nodes, cc;
  sc_array_t *new_vertices;

  lel = p4est->local_num_quadrants;

  // We begin by computing where to put the nodes owned
  // by the respective processes

  proc_offsets = P4EST_ALLOC(p4est_locidx_t, me->mpisize + 1);

  // owned = sc_array_index(me->global_owned_count, me->mpirank);
  // offset = *owned;

  for (int zz=0; zz < me->mpisize; ++zz) {
    if (zz == me->mpirank) { // We put all locally owned nodes first
      proc_offsets[zz] = 0;
      continue;
    }

    proc_offsets[zz] = offset;

    owned = sc_array_index(me->global_owned_count, zz);
    offset += *owned;
  }

  owned = sc_array_index(me->global_owned_count, me->mpirank);
  num_owned_nodes = *owned;

  num_local_nodes = proc_offsets[me->mpisize] = offset;

  // Now we move the vertices into the new order and compute
  // the new location for each local node

  proc_node_count = P4EST_ALLOC_ZERO(p4est_locidx_t, me->mpisize);
  new_local_nodes = P4EST_ALLOC(p4est_locidx_t, num_local_nodes);
  new_vertices = sc_array_new_count(3*sizeof(double), num_local_nodes);

  for (zz = 0; zz < num_local_nodes; ++zz) {
    owner = sc_array_index(me->owners, zz);

    new_cc = proc_offsets[*owner] + (proc_node_count[*owner]++);

    new_local_nodes[zz] = new_cc;
  }

  P4EST_FREE(proc_node_count);

  // Finally we remap the element_nodes to refer to the new node indices
  vnodes = SNODES_VNODES; // = sizeof(element_node_t) / sizeof(p4est_locidx_t);

  for (elem = 0; elem < lel; ++elem) {
    elem_nodes = sc_array_index(me->sa_element_nodes, elem);

    for (cc_off = 0; cc_off < vnodes; ++cc_off) {
      cc = elem_nodes[cc_off];
      elem_nodes[cc_off] = (cc < 0) ? -1 : new_local_nodes[cc];
    }
  }


  int owned_count = *( (int *) sc_array_index(me->global_owned_count, me->mpirank) );

  sc_array_destroy(me->owners);
  sc_array_destroy(me->global_owned_count);


// Send node indices

  int mpiret;
  shared_node_comm_data_t *lp;
  sc_array_t *send, *send_buf, *send_requests;

  int num_send_procs;
  int send_count, total_sent;
  int i, proc;
  double *vx;

  send_buf = P4EST_ALLOC (sc_array_t, me->mpisize);

  printf("[%d] Num local nodes %ld\n", p4est->mpirank, me->vertices->elem_count);
  printf("[%d] Owned count %d\n", p4est->mpirank, owned_count);

  for (i = 0; i < me->mpisize; i++) {
    sc_array_init (&(send_buf[i]), sizeof (shared_node_comm_data_t));
  }

  send_requests = sc_array_new(sizeof(sc_MPI_Request));

  num_send_procs = 0;
  total_sent = 0;

  shared_node_t *snode;
  for (proc=0; proc < me->mpisize; ++proc) {
    if (proc == me->mpirank)
      continue;

    send = &(send_buf[proc]);

    lp = sc_array_push(send);
    lp->node_count = num_owned_nodes;

    sc_array_t *shared_nodes = me->shared_element_nodes[proc];
    if (shared_nodes) {
      printf("[%d] rr %d\n", me->mpirank, proc);
      for (size_t iz=0; iz<shared_nodes->elem_count; ++iz) {
        snode = sc_array_index(shared_nodes, iz);


        // Remap shared_nodes, might want to do in separate loop
        cc = snode->node_index;
        snode->node_index = new_local_nodes[cc];

        vx = sc_array_index(me->vertices, snode->node_index);

        if (snode->node_index < num_owned_nodes) {
          lp = sc_array_push(send);
          lp->node_index = snode->node_index;

#ifdef P4EST_SIMPLEX_DEBUG
          lp->x = vx[0];
          lp->y = vx[1];
          lp->z = vx[2];
#endif

           printf(" * ");
        } else {
           printf("   ");
        }


#ifdef P4EST_SIMPLEX_DEBUG
        printf("[%d] SA %2d -> %2d -> %2d,    (%.4f, %.4f, %.4f)\n", me->mpirank,
            cc,
            snode->node_index,
            snode->vnode_offset,
            vx[0], vx[1], vx[2]
            );
#endif

      }
    }

    num_send_procs++;
    send_count = send->elem_count;
    sc_MPI_Request *send_request = sc_array_push(send_requests);
    mpiret = sc_MPI_Isend(
                send->array,
                (int) (send_count * sizeof(shared_node_comm_data_t)),
                sc_MPI_BYTE, proc, P4EST_COMM_SIMPLEX_PASS,
                p4est->mpicomm, send_request);
    SC_CHECK_MPI(mpiret);

    total_sent += (send_count * sizeof (shared_node_comm_data_t));

    printf("[%d] Sending %d indices to %d\n", p4est->mpirank, send_count, proc);
  }


  P4EST_INFOF ("Total of %llu bytes sent to %d processes\n",
                  (unsigned long long) total_sent, num_send_procs);


// Receive node indices

  p4est_locidx_t *global_owned_count; // , *shared_nodes;
  global_owned_count = P4EST_ALLOC(p4est_locidx_t, me->mpisize);

  global_owned_count[me->mpirank] = owned_count;

  p4est_gloidx_t *nonlocal_nodes = P4EST_ALLOC(p4est_gloidx_t, num_local_nodes - owned_count);


  sc_array_t *recv, *recv_buf;
  sc_MPI_Status probe_status, recv_status;

  int byte_count, elem_count;
  int node_count;

  recv_buf = P4EST_ALLOC (sc_array_t, me->mpisize);
  for (i = 0; i < me->mpisize; i++) {
    sc_array_init (&(recv_buf[i]), sizeof (shared_node_comm_data_t));
  }

  for (i = 0; i < me->mpisize - 1; i++) {
    printf("[%d] Probing\n", p4est->mpirank);

    mpiret = sc_MPI_Probe (sc_MPI_ANY_SOURCE, P4EST_COMM_SIMPLEX_PASS,
                           p4est->mpicomm, &probe_status);

    SC_CHECK_MPI (mpiret);
    proc = probe_status.MPI_SOURCE;
    P4EST_ASSERT (proc != p4est->mpirank);
    recv = &(recv_buf[proc]);
    mpiret = sc_MPI_Get_count (&probe_status, sc_MPI_BYTE, &byte_count);
    SC_CHECK_MPI (mpiret);
    printf("[%d] FOUND: %d bytes\n", p4est->mpirank, byte_count);
    P4EST_ASSERT (byte_count % ((int) sizeof (shared_node_comm_data_t)) == 0);
    elem_count = ((size_t) byte_count) / sizeof (shared_node_comm_data_t);
    sc_array_resize (recv, elem_count);
    mpiret = sc_MPI_Recv (recv->array, byte_count, sc_MPI_BYTE, proc,
                          P4EST_COMM_SIMPLEX_PASS, p4est->mpicomm,
                          &recv_status);
    SC_CHECK_MPI (mpiret);


    printf("[%d] Received %d indices from %d\n", p4est->mpirank, elem_count, proc);

    lp = sc_array_index(recv, 0);
    global_owned_count[proc] = lp->node_count;

    printf("[%d] R process %d has %d nodes\n", me->mpirank, proc, lp->node_count);

    node_count = 0;

    sc_array_t *shared_nodes = me->shared_element_nodes[proc];
    if (shared_nodes) {
      p4est_locidx_t next_proc_off = proc_offsets[proc+1]
                                ? proc_offsets[proc+1]
                                : proc_offsets[proc+2];

      printf("[%d] Rr[%d, %d] %d\n", me->mpirank,
              proc_offsets[proc],
              next_proc_off,
              proc);

      for (size_t iz=0; iz<shared_nodes->elem_count; ++iz) {
        snode = sc_array_index(shared_nodes, iz);

        cc = snode->node_index;

        // Check that the node is owned by process \ref proc
        if (proc_offsets[proc] <= cc && cc < next_proc_off) {
          lp = sc_array_index(recv, 1 + node_count);

          nonlocal_nodes[ snode->node_index - owned_count ] = lp->node_index;

          vx = sc_array_index(me->vertices, snode->node_index);

#ifdef P4EST_SIMPLEX_DEBUG
          printf("[%d] RA  %2d, %2d -> %2d,    (%.4f, %.4f, %.4f) -> (%.4f, %.4f, %.4f)\n", p4est->mpirank,
              proc_offsets[proc],
              snode->node_index,
              lp->node_index,
              vx[0], vx[1], vx[2],
              lp->x, lp->y, lp->z
              );

          if (vx[0] != lp->x && vx[1] != lp->y && vx[2] != lp->z) {
            SC_ABORT1("[%d] Send cords are incorrect\n", p4est->mpirank);
          }
#endif

          node_count++;
        }
      }

    }

    P4EST_ASSERT (1 + node_count == elem_count);
  }


  p4est_gloidx_t global_node_count, global_offset;
  p4est_gloidx_t proc_global_offset;
  p4est_locidx_t next_proc_off;

  proc_global_offset = 0;
  proc = -1;
  next_proc_off = -1;

  for (int zz = 0; zz < num_local_nodes - owned_count; zz++) {

    while (zz > next_proc_off) {
      proc_global_offset += global_owned_count[ proc ];
      if (++proc == me->mpirank) {
        global_offset = proc_global_offset;
        continue;
      }

      next_proc_off =  ( proc_offsets[proc+1]
                            ? proc_offsets[proc+1]
                            : proc_offsets[proc+2] ) - owned_count;

      printf("[%d] B Off(%d) %ld,  %d\n", p4est->mpirank, proc, proc_global_offset, next_proc_off);
    }

    printf("[%d]  zz %d:  %d -> %ld -> %ld\n", p4est->mpirank, zz,
        owned_count + zz,
        nonlocal_nodes[ zz ], nonlocal_nodes[ zz ] + proc_global_offset);
    nonlocal_nodes[ zz ] += proc_global_offset;
  }

  for (; proc < me->mpisize; proc++) {
    if (proc == me->mpirank)
      global_offset = proc_global_offset;
    proc_global_offset += global_owned_count[ proc ];
  }

  global_node_count = proc_global_offset;


  printf("[%d] There are a total of %ld nodes, and my offset is %ld\n", p4est->mpirank,
      global_node_count, global_offset);


  for (int zz = 0; zz < owned_count; zz++) {
    vx = sc_array_index(me->vertices, zz);

    printf("[%d] LocNode  %2d, -> %2ld,   (%.4f, %.4f, %.4f)\n", p4est->mpirank,
        zz,
        global_offset + zz,
        vx[0], vx[1], vx[2]);
  }


  // if (send_requests->elem_count > 0) {
  //   mpiret = sc_MPI_Waitall ((int) send_requests->elem_count,
  //                            (sc_MPI_Request *) send_requests->array,
  //                            sc_MPI_STATUSES_IGNORE);
  //   SC_CHECK_MPI (mpiret);
  // }

  sc_array_destroy (send_requests);

  for (i = 0; i < me->mpisize; i++) {
    sc_array_reset (&(send_buf[i]));
    sc_array_reset (&(recv_buf[i]));
  }

  P4EST_FREE (send_buf);
  P4EST_FREE (recv_buf);

  P4EST_FREE (proc_offsets);
  P4EST_FREE (new_local_nodes);

  for (int proc=0; proc < me->mpisize; ++proc) {
    sc_array_t *shared_nodes = me->shared_element_nodes[proc];
    if (!shared_nodes)
      continue;

    sc_array_destroy(shared_nodes);
  }


  P4EST_FREE (global_owned_count);
  P4EST_FREE (nonlocal_nodes);
  P4EST_FREE (me->shared_element_nodes);



  return;


  p4est_lnodes_t lnode;
  // lnode.owned_count = local_count;
  // lnode.num_local_nodes = num_local_nodes;

  // lnode.element_nodes = new_elem_nodes->array;

  // lnode.global_offset: Needs communication
  // lnode.global_owned_count: Needs communication

  // lnode.sharers : Might need to be constructed in the iterator

#endif
}

void
p4est_snodes_destory(p4est_snodes_t *snodes)
{
  P4EST_ASSERT(snodes);
  P4EST_ASSERT(snodes->element_nodes);

  P4EST_FREE (snodes->element_nodes);

  if (snodes->nonlocal_nodes)
    P4EST_FREE (snodes->nonlocal_nodes);

  if (snodes->global_owned_count)
    P4EST_FREE (snodes->global_owned_count);


  // count = lnodes->sharers->elem_count;
  // for (zz = 0; zz < count; zz++) {
  //   lrank = p4est_lnodes_rank_array_index (lnodes->sharers, zz);
  //   sc_array_reset (&(lrank->shared_nodes));
  // }
  // sc_array_destroy (snodes->sharers);

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
