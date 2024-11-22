
#include <bits/types/locale_t.h>
#include <limits.h>
#include <linux/limits.h>
#include <math.h>
#include <sc.h>
#include <sc_containers.h>
#include <sc_mpi.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "utils.h"

#ifdef VIM_LS
// #include <p4est_to_p8est.h>
#define P4EST_SIMPLEX_DEBUG
// #define TIMINGS
#endif

#include "statistics.h"

#ifndef P4_TO_P8
#include <p4est_bits.h>
#include "p4est_tnodes.h"

#include "p4est_lnodes.h"
#include "p4est_ghost.h"
#include <p4est_iterate.h>
#include "p4est_simplex_mesh.h"
#include "p4est_communication.h"
#include "p4est_geometry.h"
#else
#include <p8est_bits.h>
#include "p8est_tnodes.h"

#include "p8est_lnodes.h"
#include <p8est_iterate.h>
#include "p8est_simplex_mesh.h"
#include "p8est_communication.h"
#include "p8est_geometry.h"
#endif




#ifdef P4_TO_P8

// This code implementations the
// BECK: Big Element Centre Knut triangulation of a quadforest




/************************************************************
 * The following code belongs in other p4est header files,
 * but to avoid modifying the p4est code I have added them here.
 ************************************************************/

// The following structures are missing in p8est_iterate.h

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

// These should really be in the enum p4est_comm_tag in p4est_base.h :
static int P4EST_COMM_SNODES_QUERY = 50;
static int P4EST_COMM_SNODES_REPLY = 51;


/************************************************************
 * End of code belongs in other p4est header files.
 ************************************************************/

/************************************************************
 * The following code declares constants pertaining to the BECK implementation.
 ************************************************************/

#define SNODES_SOURCES (P4EST_DIM + 1)
#define SNODES_VNODES P4EST_DIM_POW(3)

typedef uint32_t p4est_BECK_element_config_t;
typedef uint32_t p4est_BECK_face_config_t;


/***
 *
 *    5 --- 4 --- 6
 *    |           |
 *    1     0     2
 *    |           |
 *    7 --- 3 --- 8
 */

// The n_* constants define the indexing of the nodes relative to an element.
// f.exp.  The centre node has index 0.

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



static int8_t n_face_edges_initialized = 0;
static int n_face_edges[6][4];
inline
static void
_n_face_edges_initialize()
{
  if (n_face_edges_initialized)
    return;
  for (int face=0; face < 6; ++face)
    for (int edge=0; edge < 4; ++edge) {
      n_face_edges[face][edge] = n_edge[ p8est_face_edges[face][edge] ];
    }
}

//fconfig is 5 bits.
// The most significant bit is for the face node, and the
// remaining 4 bits indicate whether or not an each edge has an edge node
//
// This defines the number of tets touching a face depending on the fconfig
static const int fconfig_tets[32] = {
   2, -1, -1, -1,     -1, -1, -1, -1,  // Without a face node we
  -1, -1, -1, -1,     -1, -1, -1, -1,  // do not allow any edge nodes

  // 000, 001, 010, 011, 100, 101, 110, 111,
     4,   5,   5,   6,   5,   6,   6,   7,
     5,   6,   6,   7,   6,   7,   7,   8
};

// Compute the fconfig for the face of an element
//
static int
config_face_config(p4est_BECK_element_config_t config, int8_t face)
{
  _n_face_edges_initialize();

  p4est_BECK_face_config_t fconf;
  fconf  = ((config >> n_face[face]) & 1) << 4;

  fconf |= ((config >> n_face_edges[face][3]) & 1) << 3;
  fconf |= ((config >> n_face_edges[face][2]) & 1) << 2;
  fconf |= ((config >> n_face_edges[face][1]) & 1) << 1;
  fconf |= ((config >> n_face_edges[face][0]) & 1) << 0;
  return fconf;
}


#else
/* face midpoints  */
static const int n_face[4] = { 1, 2, 3, 4 };

/* cube corners */
static const int n_corner[4] = { 5, 6, 7, 8 };

// fconfig is 1 bits, that is set if there is a face node
//
// This defines the number of tets touching a face depending on the fconfig
static const int fconfig_tets[2] = { 1, 2 };

// Compute the fconfig for the face of an element
static int
config_face_config(p4est_BECK_element_config_t config, int8_t face)
{
  return (config >> n_face[face]) & 1;
}


#endif

// *elem_coords_n* define the element local coordinates of the nodes
// in the range [0,2]^dim

static const int8_t elem_coords_n[ SNODES_VNODES ][3] = {
  // Center
  {1,1,1},
#ifdef P4_TO_P8
  // Faces
  {0,1,1}, {2,1,1}, // x-orth
  {1,0,1}, {1,2,1}, // y-orth
  {1,1,0}, {1,1,2}, // z-orth
  // Edges
  {1,0,0}, {1,2,0}, {1,0,2}, {1,2,2}, // x-parl
  {0,1,0}, {2,1,0}, {0,1,2}, {2,1,2}, // y-parl
  {0,0,1}, {2,0,1}, {0,2,1}, {2,2,1}, // z-parl
  // Corners
  {0,0,0}, {2,0,0}, {0,2,0}, {2,2,0},
  {0,0,2}, {2,0,2}, {0,2,2}, {2,2,2}
#else
  // Faces
  {0,1,0}, {2,1,0}, // x-orth
  {1,0,0}, {1,2,0}, // y-orth
  // Corners
  {0,0,0}, {2,0,0}, {0,2,0}, {2,2,0}
#endif
};

/************************************************************
 * End of constants for BECK
 ************************************************************/


/************************************************************
 * Static structures for BECK
 ************************************************************/

// The data required for the construction of the nodes. A pointer to an instance
// of this structure is passed as user_data in iteration step.

typedef struct
{
  sc_MPI_Comm mpicomm;
  int mpisize, mpirank;

  int *ghost_rank;                     /**< Mpi rank of the owners of the ghost elements */
  int emptypeers;                      /**< Number of empty peers,  */
  int locsharer;                       /**< Index of local sharer in sharers */
  uint8_t *chilev;                     /**< Level of the elements */

  p4est_BECK_element_config_t *configuration; /**< The element node configurations */
  /** Offsets into local triangles per element and one beyond. */
  // p4est_locidx_t     *local_element_offset;


  p4est_locidx_t *element_nodes;       /**< The node index relative to an element */
  p4est_locidx_t *ghost_element_nodes; /**< Same as element_nodes, for the ghost elements,
                                            necessary to check whether some nodes has been added before */
  sc_array_t construct;                /**< Collect nodes during traversal */
  sc_array_t ownsort;                  /**< Sorted owned nodes of local process */

  p4est_locidx_t lenum;                /**< Local element counter, incremented at volume iteration */
  p4est_locidx_t num_owned;            /**< Nodes we own */
  p4est_locidx_t num_owned_shared;     /**< Nodes we both own and share */
  p4est_locidx_t num_notowned_shared;  /**< Nodes we share, but not own */
  p4est_locidx_t num_shared;           /**< Nodes we share, owned or not */

  p4est_gloidx_t *goffset;             /**< Global offsets for owned nodes */

  // TODO: Array of element configurations using a bitmap.

#ifdef P4EST_ENABLE_MPI

  // The peers, refer to the processes that the current process share nodes with
  int *proc_peer;                      /**< Maps mpi rank to (index+1) of the
                                            corresponding peer in in peers
                                            array (0 if not in peers) */
  sc_array_t          sortp;           /**< Sorted array pointing to peers */
  sc_array_t          peers;           /**< Unsorted peer storage */
  sc_array_t          peerreq;         /**< Requests for storage */
#endif

  // Input structures
  p4est_t *p4est;
  p4est_ghost_t *ghost;
  p4est_geometry_t *geom;

  // Output structures
  p4est_lnodes_t *ln;
  p4est_tnodes_t *tn;

#ifdef TIMINGS
  statistics_t *ss;
#endif

}
BECK_nodes_meta_t;


// A 'contributor' to an element.
// Per node, there is one contributor per process touching the node.
typedef struct snodes_contr {
  int nodene;           /**< The node index relative to an element touching the node. */
  int rank;             /**< The referring process. */
  p4est_locidx_t le;    /**< Element/ghost number. */
} snodes_contr_t;

// The data we record per 'node' during the iteration step
typedef struct snodes_cnode {
  p4est_locidx_t runid;         /**< Running count of node. */
  snodes_contr_t *owner;        /**< Owning contributor, pointing into contr array. */
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

/************************************************************
 * End of Static structures for BECK
 ************************************************************/


/************************************************************
 * Static function declarations for BECK
 ************************************************************/

/* The following functions are used in the iteration step */

// Check that cnode at index lni in me->construct is valid.
// Purely a debug functionality
static void
check_node (BECK_nodes_meta_t * me, p4est_locidx_t lni);

// Compute the local element index of element quad in quadforest.
static p4est_locidx_t
tree_quad_to_le (p4est_t * p4est, p4est_topidx_t treeid,
                 p4est_locidx_t quadid);

/* Register a node contribution.
 *  (This function is an attoptation of the function in p4est(feature-triangle)
 *  by Carsten Burstedde)
 *
 * *lni should be the index of a cnode in the me->construct, and is used to
 * identify the various node contributions pertaining to a specific node.
 * For if *lni = -1, there will be constructed a new cnode and *lni will be set to
 * to the corresponding index.
 * It is important to call this function for all contributions to a specific node and
 * *lni most be have the same value as it was set to by the fist call.
 *
 * For the case, when there is only one contributor (namely the centre volume node)
 * lni may be set to NULL.
 *
 * rank is the mpi rank of the contributing element.
 * If the element is local, then rank should be -1 and le should be the local
 * element index. Otherwise, le should be the index of the element in the ghost
 * element array.
 *
 * nodene should be the element local, node index of the node.
 */
static void
node_register(
    BECK_nodes_meta_t *me, p4est_locidx_t *lni,
    int rank, p4est_locidx_t le, int nodene);

// Register a local contribution.
// le should be the index of a local element
void
node_lregister(
    BECK_nodes_meta_t *me, p4est_locidx_t *lni,
    p4est_locidx_t le, int nodene);


// Register a global contribution.
// ghostid should be the index of a ghost element
void
node_gregister(
    BECK_nodes_meta_t *me, p4est_locidx_t *lni,
    p4est_locidx_t ghostid, int nodene);


// Volume iteration
static void
iter_volume(p4est_iter_volume_info_t *info, void *user_data);

// Face iteration
static void
iter_face(p4est_iter_face_info_t *info, void *user_data);

#ifdef P4_TO_P8
// Edge iteration
static void
iter_edge(p8est_iter_edge_info_t *info, void *user_data);
#endif

// Corner iteration
static void
iter_corner(p4est_iter_corner_info_t *info, void *user_data);


/* The following functions are used to establish the sharers and
 * prolate the lnodes structure based on the cnodes constructed in
 * the iteration step
 * (The functions are adapted from p4est(feature-triangle)
 *  by Carsten Burstedde)
 */

// Sort the owned nodes.
// Order is based on the element id of the contributing
// element with the lowest index.
// If this is the same element for both nodes, then the order is that of
// the element local index.
static int
cnode_compare (const void *v1, const void *v2);

// Sort the peers according to mpi rank.
static int
peer_compare (const void *v1, const void *v2);

/** Sort remote nodes pertaining to a specific peer
 * The order is determined by the runid, which at the time of sorting is
 * peer process local index of the node received through communication.
 */
static int
rnode_compare (const void *v1, const void *v2);

/** Push sharers into lnodes sharers. */
static void
push_sharer (BECK_nodes_meta_t * me, int *sindex, int rank);




#ifdef  P4EST_ENABLE_MPI
/** Retrieve sharer data corresponding to process q, using me->proc_peer */
static p4est_lnodes_rank_t *
peer_sharer (BECK_nodes_meta_t * me, int q);

/** Retrieve peer data corresponding to process q, using me->proc_peer */
static snodes_peer_t *
peer_access (BECK_nodes_meta_t * me, int q);

#endif // P4EST_ENABLE_MPI

/** The local owner process will receive a query for a node number. */
static void
peer_add_reply (BECK_nodes_meta_t * me, snodes_peer_t * peer, p4est_locidx_t lni);

/** The local process queries a remote owner for its node number. */
static void
peer_add_query (BECK_nodes_meta_t * me, snodes_peer_t * peer,
                p4est_locidx_t lni, p4est_locidx_t epos);


/** Process the me->construct data accumulated during the iteration step.
 * The owner of a node is determined the contributing process with the lowest
 * mpirank. This is determined when we register the nodes.
 * Here we collect, count and sort the nodes owned by the peers and the current
 * process. We also construct the messages we will send to the peers requesting
 * the node indexes of the nodes owned by the peer.
 */
static void
owned_query_reply (BECK_nodes_meta_t *me);

/** Handle the request messages that we constructed in owned_query_reply.
 * Since peers with a lower mpi rank will own the shared nodes nodes,
 * we need to send our request to these peers.
 * For peers with a higher mpi rank, the current process is the owner and we
 * need to receive the requests.
 */
static void
post_query_reply (BECK_nodes_meta_t * me);

/** Sort the owned nodes and allGather owned node count, to determine the
 * global node offsets and global node count.
 *
 * TODO: allGather the owned tets counts in parallel.
 */
static void
allgather_counts (BECK_nodes_meta_t * me);

/** Sort the peers and initialize the lnodes shares array */
static void
sort_peers (BECK_nodes_meta_t * me);

/** We now wait for the send and receive request in post_query_reply
 * to finalize.
 * For the messages:
 *  - we have sent we, start waiting for a reply.
 *  - we have received, we send a reply
 * When we have sent a reply we are done, when we revive a reply we need to
 * handed the reply and we are done.
 */
static void
wait_query_reply (BECK_nodes_meta_t * me);

/** Finalize the element_nodes.
 * At the moment the values in element_nodes are indexing the
 * me->construct array. Instead we want it to be the local node
 * index. Following the wait_query_reply call, the cnodes->runid is set
 * to the correct node index. Therefore, what we need to do is, update
 * the element_nodes array as follows. and update the value at index i, to be
 * element_nodes[i] := me->construct[ element_nodes[i] ]->runid.
 */
static void
assign_element_nodes (BECK_nodes_meta_t * me);

/** Finalize the lnodes->sharers */
static void
populate_sharers (BECK_nodes_meta_t * me);

/** Clean up the cnodes in me->construct, constructed at the iteration step */
static void
clean_construct (BECK_nodes_meta_t * me);


/** Simplex construction functions */

/** Append simplex */
static void
snodes_push_simplex(
    p4est_tnodes_t *tnodes,
    const p4est_locidx_t *enodes,
    const int eindex[SNODES_SOURCES]);

static void
p4est_geometry_coordinates_BECK
  (p4est_t *p4est, p4est_lnodes_t *ln,
   p4est_geometry_t *geom,
   sc_array_t *coordinates);

static void
collect_local_tets(
    p4est_t *p4est,
    p4est_tnodes_t *tnodes,
    p4est_lnodes_t *ln);

/************************************************************
 * End of static function declarations for BECK
 ************************************************************/

/************************************************************
 * Static function implementation for BECK
 ************************************************************/

static void
check_node (BECK_nodes_meta_t * me, p4est_locidx_t lni)
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

static void
node_register(
    BECK_nodes_meta_t *me, p4est_locidx_t *lni,
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
    P4EST_ASSERT(0 <= le && le < me->ln->num_local_elements);

    P4EST_ASSERT (me->element_nodes[le * SNODES_VNODES + nodene] == -1);
    me->element_nodes[le * SNODES_VNODES + nodene] = *lni;
    me->configuration[le] |= (1 << nodene);
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
    BECK_nodes_meta_t *me, p4est_locidx_t *lni,
    p4est_locidx_t le, int nodene
    )
{
  P4EST_ASSERT (0 <= nodene && nodene < SNODES_VNODES);
  P4EST_ASSERT (0 <= le && le < me->p4est->local_num_quadrants);

  P4EST_ASSERT (lni == NULL
                || *lni == -1
                || me->element_nodes[le * SNODES_VNODES + nodene] == -1
                || *lni == me->element_nodes[le * SNODES_VNODES + nodene]);

  /* Check if this node is already registered we don't need to do anything */
  if ( me->element_nodes[le * SNODES_VNODES + nodene] != -1 ) {
    P4EST_ASSERT (lni != NULL );

    *lni = me->element_nodes[le * SNODES_VNODES + nodene];
    return;
  }

  node_register(me, lni, -1, le, nodene);
}

void
node_gregister(
    BECK_nodes_meta_t *me, p4est_locidx_t *lni,
    p4est_locidx_t ghostid, int nodene
    )
{
  p4est_locidx_t lnis;
  p4est_quadrant_t *gquad;

  P4EST_ASSERT (me != NULL);
  P4EST_ASSERT (0 <= nodene && nodene < SNODES_VNODES);

  if (me->ghost != NULL) {
    P4EST_ASSERT (me->ghost_rank != NULL);
    P4EST_ASSERT (0 <= ghostid &&
                  ghostid < (p4est_locidx_t) me->ghost->ghosts.elem_count);

    /* Check if this node is already registered we don't need to do anything */
    if ( me->ghost_element_nodes[ghostid * SNODES_VNODES + nodene] != -1 ) {
      P4EST_ASSERT (lni != NULL );

      *lni = me->ghost_element_nodes[ghostid * SNODES_VNODES + nodene];
      return;
    }

    if (lni == NULL) {
      lnis = -1;
      lni = &lnis;
    }

    /* extract remote element number from ghost quadrant */
    gquad = (p4est_quadrant_t *) sc_array_index (&me->ghost->ghosts, ghostid);
    node_register (me, lni, me->ghost_rank[ghostid],
                   gquad->p.piggy3.local_num, nodene);

    me->ghost_element_nodes[ghostid * SNODES_VNODES + nodene] = *lni;
  }
}

// ITERATORS
static void
iter_volume(p4est_iter_volume_info_t *vi, void *user_data)
{
  BECK_nodes_meta_t *me = (BECK_nodes_meta_t *) user_data;
  p4est_locidx_t le;
  p4est_locidx_t childid;
  int8_t level;

  TSTAT_BEGIN(me->ss, TIMINGS_ITER_VOLUMES);

  le = me->lenum++;
  level = vi->quad->level;
  childid = p4est_quadrant_child_id (vi->quad);
  me->chilev[le] = (((uint8_t) level) << 3) | ((uint8_t) childid);

  node_lregister(me, NULL, le, n_center);

  TSTAT_END(me->ss, TIMINGS_ITER_VOLUMES);
}

static void
iter_face(p4est_iter_face_info_t *fi, void *user_data)
{
  BECK_nodes_meta_t *me = (BECK_nodes_meta_t *) user_data;

  int nodene;
  p4est_locidx_t le;              /**< local element number */
  p4est_locidx_t lni;             /**< local node number */

  TSTAT_BEGIN(me->ss, TIMINGS_ITER_FACES);

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

#ifdef TIMINGS
  if (!sfd->is_ghost)
    counter_inc(me->ss, 1, COUNTER_HANGING_FACE);
#endif



  P4EST_ASSERT(sh && shd && sf && shd);

  lni = -1;
  nodene = n_face[sf->face];

  if (!sfd->is_ghost) {
    le = tree_quad_to_le(me->p4est, sf->treeid, sfd->quadid);
    lni = me->element_nodes[le*SNODES_VNODES + nodene];
    if (lni >= 0)
      return;

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
  }

  TSTAT_END(me->ss, TIMINGS_ITER_FACES);
}

#ifdef P4_TO_P8
static void
iter_edge(p8est_iter_edge_info_t *info, void *user_data)
{
  BECK_nodes_meta_t *me = (BECK_nodes_meta_t *) user_data;

  size_t sz, zz, efz;
  int face, edge, nodene;
  int faces_count;
  int hanging;

  p4est_locidx_t le;
  p4est_locidx_t lni, *face_nodes;

  p8est_iter_edge_side_t *side;
  p8est_iter_edge_side_full_t *sfd;
  p8est_iter_edge_side_hanging_t *shd;

  p4est_quadrant_t *quad;

  TSTAT_BEGIN(me->ss, TIMINGS_ITER_EDGES);
#ifdef TIMINGS
  int orank = INT_MAX, srank;
#endif

  hanging = 0;

  // Find a full and a hanging face
  for (zz = 0; zz < info->sides.elem_count; ++zz) {
    side = sc_array_index(&info->sides, zz);
    hanging |= side->is_hanging;

    faces_count = SC_MAX(faces_count,
                    1 + SC_MAX(side->faces[0], side->faces[1]));
  }

  // If no side is hanging we will not add an edge node
  if (!hanging) {
    return;
  }

  face_nodes = MY__ALLOC_N(p4est_locidx_t, faces_count, "face_nodes");
  memset(face_nodes, -1, sizeof(p4est_locidx_t) * faces_count);

  lni = -1;

  // We record the index in the element_nodes for each side
  for (zz = 0; zz < info->sides.elem_count; ++zz) {
    side = sc_array_index(&info->sides, zz);
    edge = side->edge;

    if (!side->is_hanging) {
      sfd = (p8est_iter_edge_side_full_t *) &side->is.full;

 #ifdef TIMINGS
      srank = side->is.full.is_ghost ? me->ghost_rank[side->is.full.quadid] : me->mpirank;
      if (srank < orank)
        orank = srank;
 #endif

      nodene = n_edge[side->edge];
      if (!sfd->is_ghost) {
        le = tree_quad_to_le(me->p4est, side->treeid, sfd->quadid);
        node_lregister(me, &lni, le, nodene);
      } else {
        node_gregister(me, &lni, sfd->quadid, nodene);
      }

      // For each face adjacent to a hanging edge, we also add a face node
      for (efz = 0; efz < 2; ++efz) {
        face = p8est_edge_faces[edge][efz];
        nodene = n_face[ face ];

        if (!sfd->is_ghost) {
          le = tree_quad_to_le(me->p4est, side->treeid, sfd->quadid);

          node_lregister(me, &face_nodes[ side->faces[efz] ], le, nodene);
        } else {
          node_gregister(me, &face_nodes[ side->faces[efz] ], sfd->quadid, nodene);
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

  MY__FREE(face_nodes);

#ifdef TIMINGS
  if (orank == me->mpirank)
    counter_inc(me->ss, 1, COUNTER_HANGING_EDGE);
#endif

  TSTAT_END(me->ss, TIMINGS_ITER_EDGES);
}
#endif

static void
iter_corner(p4est_iter_corner_info_t *info, void *user_data)
{
  BECK_nodes_meta_t *me = (BECK_nodes_meta_t *) user_data;

  p4est_iter_corner_side_t *side;

  p4est_locidx_t le;
  p4est_locidx_t lni;
  int nodene;

  size_t zz;

  TSTAT_BEGIN(me->ss, TIMINGS_ITER_CORNERS);

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

  TSTAT_END(me->ss, TIMINGS_ITER_CORNERS);
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
push_sharer (BECK_nodes_meta_t * me, int *sindex, int rank)
{
  p4est_lnodes_rank_t *sharer;
  p4est_lnodes_t *ln = me->ln;

  P4EST_ASSERT (me != NULL);
  P4EST_ASSERT (sindex != NULL);
  P4EST_ASSERT (0 <= rank && rank < me->mpisize);

  /* push empty sharer structure */
  *sindex = (int) ln->sharers->elem_count;
  sharer = (p4est_lnodes_rank_t *) sc_array_push (ln->sharers);
  memset (sharer, -1, sizeof (p4est_lnodes_rank_t));
  sharer->rank = rank;
  sc_array_init (&sharer->shared_nodes, sizeof (p4est_locidx_t));
}
#endif

static void
clean_construct (BECK_nodes_meta_t * me)
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
peer_sharer (BECK_nodes_meta_t * me, int q)
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
  return (p4est_lnodes_rank_t *) sc_array_index_int (me->ln->sharers,
                                                     peer->sharind);
}

static snodes_peer_t *
peer_access (BECK_nodes_meta_t * me, int q)
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
peer_add_reply (BECK_nodes_meta_t * me, snodes_peer_t * peer, p4est_locidx_t lni)
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
peer_add_query (BECK_nodes_meta_t * me, snodes_peer_t * peer,
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
populate_ghost_rank(BECK_nodes_meta_t *me)
{
  p4est_locidx_t lg, ng;
  int q;

  p4est_ghost_t *ghost = me->ghost;

  P4EST_ASSERT(ghost != NULL);
  P4EST_ASSERT (ghost->proc_offsets[0] == 0);
  P4EST_ASSERT (ghost->proc_offsets[me->mpisize] ==
                (p4est_locidx_t) ghost->ghosts.elem_count);

  me->ghost_rank = MY__ALLOC (int, ghost->ghosts.elem_count);
  lg = 0;
  for (q = 0; q < me->mpisize; ++q) {
    ng = ghost->proc_offsets[q + 1];
    for (; lg < ng; ++lg) {
      me->ghost_rank[lg] = q;
    }
  }
  P4EST_ASSERT (lg == (p4est_locidx_t) ghost->ghosts.elem_count);
  me->ghost_element_nodes = MY__ALLOC(p4est_locidx_t, lg * SNODES_VNODES);
  memset (me->ghost_element_nodes, -1, lg * SNODES_VNODES * sizeof(p4est_locidx_t));
#ifdef P4EST_ENABLE_MPI
  me->proc_peer = MY__ALLOC_ZERO (int, me->mpisize);
  sc_array_init (&me->sortp, sizeof (snodes_peer_t *));
  sc_array_init (&me->peers, sizeof (snodes_peer_t));
  sc_array_init (&me->peerreq, sizeof (sc_MPI_Request));
#endif
}

static void
owned_query_reply (BECK_nodes_meta_t *me)
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
post_query_reply (BECK_nodes_meta_t * me)
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
wait_query_reply (BECK_nodes_meta_t * me)
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
  p4est_lnodes_t     *ln = me->ln;
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
  waitind = MY__ALLOC (int, nwalloc);
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
            P4EST_ASSERT (0 <= epos && epos < SNODES_VNODES * ln->owned_count);

            lni = ln->element_nodes[epos];
            P4EST_ASSERT (0 <= lni && lni <
                          (p4est_locidx_t) me->construct.elem_count);
            cnode = (snodes_cnode_t *) sc_array_index (&me->construct, lni);
            oind = cnode->runid;

            P4EST_ASSERT (0 <= oind && oind < ln->owned_count);

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
            //  This check is not possible because
            //  global_owned_count is not populated yet
            // P4EST_ASSERT (0 <= oind &&
            //               oind < ln->global_owned_count[peer->rank]);
            cnode->runid = oind;
          }
          sc_array_sort (&peer->remosort, rnode_compare);

          /* store shared node's global index */
          for (lcl = 0; lcl < lbc; ++lcl) {
            cnode =
              *(snodes_cnode_t **) sc_array_index (&peer->remosort, lcl);
            P4EST_ASSERT (cnode->owner->rank == peer->rank);
            P4EST_ASSERT (lcl <= cnode->runid);
            nonloc = peer->shacumul + lcl;
            P4EST_ASSERT (nonloc < me->num_shared);

            /* we temperately set nonlocal_nodes to be the peers local index,
             * since we still do not have access to the global offsets */
            ln->nonlocal_nodes[nonloc] = cnode->runid;

            /* now the runid of each node is the local node number */
            cnode->runid = me->num_owned + nonloc;
          }
          peer->done = 0;
          --nwtotal;
        }
      }
    }
  }
  MY__FREE (waitind);
// #ifdef P4EST_ENABLE_DEBUG
//   gof = -1;
//   for (lcl = 0; lcl < me->num_notowned_shared; ++lcl) {
//     gni = ln->nonlocal_nodes[lcl];
//     P4EST_ASSERT (0 <= gni && gni < me->goffset[me->mpisize]);
//     P4EST_ASSERT (gni < me->goffset[me->mpirank] ||
//                   gni >= me->goffset[me->mpirank + 1]);
//     P4EST_ASSERT (gni > gof);
//     gof = gni;
//   }
// #endif /* P4EST_ENABLE_DEBUG */
#endif /* P4EST_ENABLE_MPI */
}

static void
assign_element_nodes (BECK_nodes_meta_t * me)
{
  int                 nodene;
  p4est_locidx_t      le, lel;
  p4est_locidx_t      lni, runid;
  snodes_cnode_t     *cnode;
  p4est_lnodes_t     *ln = me->ln;

  lel = me->ln->num_local_elements;
  for (le = 0; le < lel; ++le) {
    for (nodene = 0; nodene < SNODES_VNODES; ++nodene) {
      if ((lni = me->element_nodes[le * SNODES_VNODES + nodene]) < 0)
        continue;
      P4EST_ASSERT (lni < (p4est_locidx_t) me->construct.elem_count);

      cnode = (snodes_cnode_t *) sc_array_index (&me->construct, lni);
      runid = cnode->runid;
// #ifdef P4EST_ENABLE_DEBUG
//       if (runid >= me->num_owned) {
//         lni = (p4est_locidx_t) (ln->nonlocal_nodes[runid - me->num_owned] -
//                                 me->goffset[cnode->owner->rank]);
//         P4EST_ASSERT (0 <= lni &&
//                       lni < ln->global_owned_count[cnode->owner->rank]);
//       }
// #endif
      ln->element_nodes[le * ln->vnodes + nodene] = runid;
    }
  }
}

static void
populate_sharers (BECK_nodes_meta_t * me)
{
#ifdef P4EST_ENABLE_MPI
  int                 i;
  p4est_gloidx_t      gof, gni;
  int                 num_peers;
  size_t              zz;
  p4est_locidx_t      nonloc, lni;
  p4est_locidx_t      lbc, lcl;
  p4est_lnodes_t     *ln = me->ln;
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
  P4EST_ASSERT (num_peers + 1 == (p4est_locidx_t) ln->sharers->elem_count);
  locshare =
    (p4est_lnodes_rank_t *) sc_array_index_int (ln->sharers, me->locsharer);

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
      (p4est_lnodes_rank_t *) sc_array_index_int (ln->sharers, tp->sharind);
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

      gof = me->goffset[tp->rank];
      for (lcl = 0; lcl < lbc; ++lcl, ++lni) {
        cnode = *(snodes_cnode_t **) sc_array_index (&tp->remosort, lcl);
        P4EST_ASSERT (cnode->owner->rank == tp->rank);
        P4EST_ASSERT (cnode->runid == lni);

        /* populate nonlocal_nodes */
        P4EST_ASSERT (me->num_notowned_shared
                      == (me->num_shared - me->num_owned_shared));
        P4EST_ASSERT (lni < ln->num_local_nodes);

        nonloc = lni - me->num_owned;
        gni = gof + ln->nonlocal_nodes[nonloc]; // cnode->runid;
                                             //
        if (! (me->goffset[tp->rank] <= gni &&
                      gni < me->goffset[tp->rank + 1]))
        {
            printf("H\n");
        }
                                             //
        P4EST_ASSERT (me->goffset[tp->rank] <= gni &&
                      gni < me->goffset[tp->rank + 1]);
        ln->nonlocal_nodes[nonloc] = gni;

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
  P4EST_ASSERT (lni == ln->num_local_nodes);
#endif /* P4EST_ENABLE_MPI */
}

static void
allgather_counts (BECK_nodes_meta_t * me)
{
  int                mpiret;
  int                q, s = me->mpisize;
  p4est_locidx_t    *localboth, lb[2];
  p4est_gloidx_t     gc, gt;

  p4est_lnodes_t    *ln = me->ln;
  p4est_tnodes_t    *tn = me->tn;

  /* share owned count */
  lb[0] = ln->owned_count;
  lb[1] = tn->simplices->elem_count;

  // /* parallel sharing of owned node and element counts */
  localboth = MY__ALLOC (p4est_locidx_t, 2 * me->mpisize);
  mpiret = sc_MPI_Allgather (lb, 2, P4EST_MPI_LOCIDX,
                             localboth, 2, P4EST_MPI_LOCIDX,
                             me->mpicomm);
  SC_CHECK_MPI (mpiret);

  /* compute offsets */
  me->goffset = MY__ALLOC (p4est_gloidx_t, s + 1);
  ln->global_owned_count = MY__ALLOC (p4est_locidx_t, s);
  tn->local_tcount = MY__ALLOC (p4est_locidx_t, s);

  gc = me->goffset[0] = 0;
  gt = 0;

  for (q = 0; q < s; ++q) {
    if (q == me->mpirank) {
      ln->global_offset = gc;
      tn->global_toffset = gt;
    }

    gc = me->goffset[q + 1] =
      gc + (ln->global_owned_count[q] = localboth[2 * q + 0]);
    gt += (tn->local_tcount[q] = localboth[2 * q + 1]);
  }
  ln->global_offset = me->goffset[me->mpirank];
  MY__FREE (localboth);
}

static void
sort_peers (BECK_nodes_meta_t * me)
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
  P4EST_ASSERT (num_peers + 1 == (int) me->ln->sharers->elem_count);
  P4EST_ASSERT (me->locsharer >= 0);
#endif /* P4EST_ENABLE_MPI */
}

/************************************************************
 * End of static function implementation for BECK
 ************************************************************/



/** This is the function that generates the
 *
 *  Big Element Centre Knut triangulation
 *
 **/

static void
p4est_geometry_coordinates_BECK
  (p4est_t *p4est, p4est_lnodes_t *ln,
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
  P4EST_ASSERT (ln != NULL);
  P4EST_ASSERT (coordinates != NULL &&
                coordinates->elem_size == 3 * sizeof (double));

  sc_array_resize(coordinates, ln->num_local_nodes);
  mask = MY__ALLOC_ZERO(char, ln->num_local_nodes);
  enodes = ln->element_nodes;

  for (tt = p4est->first_local_tree; tt <= p4est->last_local_tree; ++tt) {

    tree = p4est_tree_array_index(p4est->trees, tt);
    quad_count = tree->quadrants.elem_count;
    for (qz = 0; qz < quad_count; ++qz, enodes += SNODES_VNODES) {
      quad = p4est_quadrant_array_index(&tree->quadrants, qz);

      for (k = 0; k < SNODES_VNODES; ++k) {
        lni = enodes[k];
        P4EST_ASSERT(lni < ln->num_local_nodes);

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

  MY__FREE(mask);

}

static void
p4est_tnodes_icoord_arrow (const int a[P4EST_DIM], const int b[P4EST_DIM],
                           int r[P4EST_DIM])
{
  int                 j;

  /* compute b - a */
  for (j = 0; j < P4EST_DIM; ++j) {
    r[j] = b[j] - a[j];
  }
}

static void
p4est_tnodes_icoord_cross (const int a[P4EST_DIM], const int b[P4EST_DIM],
                           int *r)
{
  /* compute cross product */
#ifndef P4_TO_P8
  *r = a[0] * b[1] - a[1] * b[0];
#else
  r[0] = a[1] * b[2] - a[2] * b[1];
  r[1] = a[2] * b[0] - a[0] * b[2];
  r[2] = a[0] * b[1] - a[1] * b[0];
#endif
}

#ifdef P4_TO_P8

static int
p4est_tnodes_icoord_inner (const int a[P4EST_DIM], const int b[P4EST_DIM])
{
  int                 j;
  int                 r;

  /* compute inner product */
  r = 0;
  for (j = 0; j < P4EST_DIM; ++j) {
    r += a[j] * b[j];
  }
  return r;
}

#endif

static void
snodes_push_simplex(
    p4est_tnodes_t *tnodes,
    const p4est_locidx_t *enodes,
    const int eindex[SNODES_SOURCES]
    )
{
  int i, j;
  int icoord[SNODES_SOURCES][P4EST_DIM];
  int taxes[P4EST_DIM][P4EST_DIM];
  int product;

#ifndef P4_TO_P8
  int                *cross = &product;
#else
  int                 cross[P4EST_DIM];
#endif

  p4est_locidx_t *simplex, windex;

#ifdef P4EST_ENABLE_DEBUG
  /* verify range of node indices */
  for (i = 0; i < SNODES_SOURCES; ++i) {
    P4EST_ASSERT (0 <= eindex[i] && eindex[i] < P4EST_INSUL);
    for (j = 0; j < i; ++j) {
      P4EST_ASSERT (eindex[j] != eindex[i]);
    }
  }
#endif

  simplex = (p4est_locidx_t *) sc_array_push(tnodes->simplices);

  /* local nodes are identified with coordinates */
  for (i = 0; i < SNODES_SOURCES; ++i) {
    P4EST_ASSERT(enodes[eindex[i]] >= 0);

    simplex[i] = enodes[eindex[i]];
  }

  return;

  // TODO /* ensure right-handed orientation of simplex */
  for (i = 0; i < SNODES_SOURCES; ++i) {
    icoord[i][0] = elem_coords_n[ eindex[i] ][0];
    icoord[i][1] = elem_coords_n[ eindex[i] ][1];
#ifdef P4_TO_P8
    icoord[i][2] = elem_coords_n[ eindex[i] ][2];
#endif
  }
  for (j = 0; j < P4EST_DIM; ++j) {
    p4est_tnodes_icoord_arrow (icoord[0], icoord[j + 1], taxes[j]);
  }
  p4est_tnodes_icoord_cross (taxes[0], taxes[1], cross);
#ifdef P4_TO_P8
  product = p4est_tnodes_icoord_inner (cross, taxes[2]);
#endif
  P4EST_ASSERT (product != 0);

  if (product < 0) {
    P4EST_LDEBUG("Flipped simplex");

    windex = simplex[1];
    simplex[1] = simplex[2];
    simplex[2] = windex;
  }
}

static void
collect_local_tets(
    p4est_t *p4est,
    p4est_tnodes_t *tnodes,
    p4est_lnodes_t *ln
    )
{
  size_t zz, fi;
  int eindex[SNODES_VNODES];
  int cf, c0, c1;
#ifdef P4_TO_P8
  size_t ei;
  int edge;
  int ce, c2, c3;
#endif

  p4est_topidx_t t, nt, lft, llt;
  p4est_locidx_t le, lel, nte;
  p4est_locidx_t *le_nodes;
  p4est_tree_t *tree;
  sc_array_t *quadrants;
  p4est_quadrant_t qid, *quad;

  tnodes->simplices = sc_array_new((P4EST_DIM + 1)*sizeof(p4est_locidx_t));

  // The last index is always the center
  eindex[P4EST_DIM] = n_center;

  lel = ln->num_local_elements;

  tnodes->local_element_offset = MY__ALLOC (p4est_locidx_t, lel + 1);
  tnodes->local_element_offset[0] = 0;

  tnodes->local_first_tree = lft = p4est->first_local_tree;
  tnodes->local_last_tree = llt = p4est->last_local_tree;

  nt = llt - lft;
  tnodes->local_tree_offset = MY__ALLOC (p4est_locidx_t, nt + 1);

  // We iterate through the local quads to make the simplices
  le = 0;
  for (t = lft; t <= llt; ++t) {

    tree = p4est_tree_array_index(p4est->trees, t);
    quadrants = &tree->quadrants;
    nte = quadrants->elem_count;

    for (zz = 0; zz < nte; ++zz, ++le) {
      quad = p4est_quadrant_array_index(quadrants, zz);

      qid = *quad;
      qid.p.which_tree = t;

      le_nodes = &ln->element_nodes[le*SNODES_VNODES];

      for (fi = 0; fi < P4EST_FACES; ++fi) {
        // Check of there is a face node
        cf = n_face[fi];
        if (0 <= le_nodes[ cf ])
        {
          eindex[0] = n_face[fi];

#ifndef P4_TO_P8
          eindex[1] = n_corner[p4est_face_corners[fi][0]];
          snodes_push_simplex(tnodes, le_nodes, eindex);

          eindex[1] = n_corner[p4est_face_corners[fi][1]];
          snodes_push_simplex(tnodes, le_nodes, eindex);
#else
          for (ei = 0; ei < 4; ++ei) {
            edge = p8est_face_edges[fi][ei];
            ce = n_edge[edge];

            if (le_nodes[ ce ] < 0)
            {
              eindex[1] = n_corner[p8est_edge_corners[edge][0]];
              eindex[2] = n_corner[p8est_edge_corners[edge][1]];
              snodes_push_simplex(tnodes, le_nodes, eindex);
            }
            else
            {
              eindex[1] = ce;

              eindex[2] = n_corner[p8est_edge_corners[edge][0]];
              snodes_push_simplex(tnodes, le_nodes, eindex);

              eindex[2] = n_corner[p8est_edge_corners[edge][1]];
              snodes_push_simplex(tnodes, le_nodes, eindex);
            }
          }
#endif
        }
        else
        {  // Split as usual
          c0 = n_corner[p4est_face_corners[fi][0]];
          c1 = n_corner[p4est_face_corners[fi][1]];
#ifndef P4_TO_P8
          eindex[0] = c0;
          eindex[1] = c1;
          snodes_push_simplex(tnodes, le_nodes, eindex);
#else
          c2 = n_corner[p4est_face_corners[fi][2]];
          c3 = n_corner[p4est_face_corners[fi][3]];

          // When splitting the face we most ensure that the split
          // is consistent on both sides of the faces.
          int offd_split = 0;
          if (quad->level == 0) {
            // At level zero we achieve this by splitting along the diagonal touching the vertex
            // with the lowest global index, which is independent of face-orientation

            p4est_gloidx_t g0 = p4est_lnodes_global_index(tnodes->lnodes, le_nodes[c0]);
            p4est_gloidx_t g1 = p4est_lnodes_global_index(tnodes->lnodes, le_nodes[c1]);
            p4est_gloidx_t g2 = p4est_lnodes_global_index(tnodes->lnodes, le_nodes[c2]);
            p4est_gloidx_t g3 = p4est_lnodes_global_index(tnodes->lnodes, le_nodes[c3]);

            offd_split = SC_MIN(g1, g2) < SC_MIN(g0, g3);
          } else {
            // At higher levels we split the face, so that the parent face is split by an X:
            //    +---+---+
            //    | \ | / |
            //    +---+---+
            //    | / | \ |
            //    +---+---+
            // This way, the split is the same for all face-orientations (rotation, mirroring)
            // and thus the split will always be the same on any neighbouring face.
            int dim = (fi >> 1);
            int c = p4est_quadrant_child_id(quad);
            int fchild =
                (dim == 0 ? c >> 1               : 0)
              | (dim == 1 ? c & 1 | (c >> 1) & 2 : 0)
              | (dim == 2 ? c & 3                : 0);

            offd_split = (fchild == 1 || fchild == 2);
          }

          // Normally (not \ref offd_split) the face is split along the c0-c3 diagonal:
          //   c2--c3
          //   |  / |
          //   | /  |
          //   c0--c1
          //
          // Otherwise the face is split along the c1-c2 diagonal:
          //   c2--c3
          //   | \  |
          //   |  \ |
          //   c0--c1

          if (!offd_split) {
            eindex[0] = c0;
            eindex[2] = c3;

            eindex[1] = c1;
            snodes_push_simplex(tnodes, le_nodes, eindex);

            eindex[1] = c2;
            snodes_push_simplex(tnodes, le_nodes, eindex);
          } else
          {
            eindex[0] = c1;
            eindex[2] = c2;

            eindex[1] = c0;
            snodes_push_simplex(tnodes, le_nodes, eindex);

            eindex[1] = c3;
            snodes_push_simplex(tnodes, le_nodes, eindex);
          }
#endif
        } // if fc == 0

      }  // for (fi = 0..P4EST_FACES)

      /* update element simplex offset list */
      tnodes->local_element_offset[le + 1] =
        (p4est_locidx_t) tnodes->simplices->elem_count;


    } // for (zz = 0..|quadrants|)

    tnodes->local_tree_offset[t - lft + 1] =
      (p4est_locidx_t) tnodes->simplices->elem_count;

  } // for (t = first-..last-localtrees)

  P4EST_ASSERT(le == lel);
  P4EST_ASSERT(nt < 0 || tnodes->local_tree_offset[llt - lft + 1]
                              == tnodes->simplices->elem_count);

  P4EST_ASSERT(tnodes->local_element_offset[lel]
                              == tnodes->simplices->elem_count);

  P4EST_INFOF ("Created %ld local simplices\n",
               tnodes->simplices->elem_count );
}

p4est_tnodes_t *
p4est_new_BECK_mesh(
    p4est_t *p4est,
    p4est_geometry_t *geom,
    p4est_ghost_t *ghost
#ifdef TIMINGS
    , statistics_t *ss
#endif
    )
{
  int mpisize, mpirank;
  p4est_locidx_t lel;
  size_t zz;
  snodes_cnode_t **ccn;
#ifdef P4EST_ENABLE_MPI
  size_t              nz, zi;
  snodes_peer_t      *peer;
#endif

  BECK_nodes_meta_t snmeta, *me = &snmeta;
  p4est_lnodes_t *ln;
  p4est_tnodes_t *tnodes;

  // Basic checks
  P4EST_ASSERT(p4est != NULL);
  P4EST_ASSERT(geom != NULL);
  P4EST_ASSERT(ghost != NULL);  /// Might remove

  // shorthand
  mpirank = p4est->mpirank;
  mpisize = p4est->mpisize;
  lel = p4est->local_num_quadrants;

  // Initialize meta structure
  memset (me, 0, sizeof (*me));

  me->mpicomm = p4est->mpicomm;
  me->mpisize = p4est->mpisize;
  me->mpirank = p4est->mpirank;
  me->p4est = p4est;
  me->geom = geom;
  me->ghost = ghost;

#ifdef TIMINGS
  P4EST_ASSERT( ss != NULL);
  me->ss = ss;
#endif

  me->configuration = MY__ALLOC_ZERO_N(p4est_BECK_face_config_t, lel, "me->configuration");
  // me->local_element_offset = MY__ALLOC (p4est_locidx_t, lel + 1);
  me->chilev = MY__ALLOC_ZERO_N(uint8_t, lel, "me->chilev");
  me->element_nodes = MY__ALLOC_N(p4est_locidx_t, lel * SNODES_VNODES, "me->element_nodes");
  memset (me->element_nodes, -1, lel * SNODES_VNODES * sizeof(p4est_locidx_t));

  sc_array_init (&me->construct, sizeof (snodes_peer_t));
  sc_array_init (&me->ownsort, sizeof (sc_MPI_Request));

  if (ghost != NULL)
    populate_ghost_rank(me);


  // Initialize lnodes
  ln = me->ln = MY__ALLOC_ZERO(p4est_lnodes_t, 1);
  ln->num_local_elements = lel;
  ln->element_nodes = me->element_nodes;
  ln->mpicomm = me->mpicomm;
  ln->sharers = sc_array_new(sizeof(p4est_lnodes_rank_t));
  ln->degree = 0;
  ln->vnodes = SNODES_VNODES;

  // Initialize tnodes structure
  tnodes = me->tn = MY__ALLOC_ZERO(p4est_tnodes_t, 1);
  tnodes->lnodes = ln;
  tnodes->lnodes_owned = 0;
  tnodes->local_first_child = -1;

  tnodes->coord_to_lnode = NULL;
  tnodes->local_element_level = MY__ALLOC_N(int8_t, 1, "local_elem");
  tnodes->coordinates = sc_array_new(3 * sizeof(double));

  /**** Iteration Step ******/

  // Find all local nodes

  TSTAT_BEGIN(me->ss, TIMINGS_ITET_STEP);

  p4est_iterate (p4est, me->ghost, me,
      iter_volume, iter_face,
#ifdef P4_TO_P8
      iter_edge,
#endif
     iter_corner);

  TSTAT_END(me->ss, TIMINGS_ITET_STEP);
  P4EST_ASSERT(me->lenum == lel);

  /* perp peer communication */

  owned_query_reply (me);
  P4EST_INFOF ("p4est_snodes_new: nodes owned %ld shared %ld\n",
               (long) me->num_owned, (long) me->num_notowned_shared);

  ln->owned_count = me->num_owned;
  ln->num_local_nodes = me->num_owned + me->num_notowned_shared;

#ifdef TIMINGS
  // Count the number of nodes that are shared between processes
  counter_inc(me->ss, me->num_owned_shared, COUNTER_SHARED_NODES);
#endif


  /* post first round of messages */
  post_query_reply (me);

  /* sort local node list */
  sc_array_sort (&me->ownsort, cnode_compare);
  for (zz = 0; zz < me->ownsort.elem_count; ++zz) {
    ccn = (snodes_cnode_t **) sc_array_index (&me->ownsort, zz);
    (*ccn)->runid = zz;
  }

  /* sort peers by process */
  sort_peers(me);

  ln->nonlocal_nodes = MY__ALLOC (p4est_gloidx_t, me->num_notowned_shared);;

  /* receive query messages and send replies */
  wait_query_reply (me);

  /* finalize element node assignment */
  assign_element_nodes (me);


  /* Now that element_nodes is finalized we may construct the tets */
  collect_local_tets(p4est, tnodes, ln);

  // if (flag & COORDS)
  p4est_geometry_coordinates_BECK
    (p4est, ln, geom, tnodes->coordinates);

  /* allgather owned node and tet counts */
  allgather_counts (me);
  P4EST_INFOF ("p4est_snodes_new: triangles owned %ld\n",
               (long) tnodes->local_tcount[mpisize]);
  P4EST_GLOBAL_PRODUCTIONF
    ("p4est_snodes_new: global triangles %lld nodes %lld\n",
     (long long) tnodes->global_tcount, (long long) me->goffset[mpisize]);

  /* finalize sharer information */
  populate_sharers (me);

  /* cleanup */
  MY__FREE (me->goffset);
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

    MY__FREE(me->proc_peer);
    sc_array_reset(&me->sortp);
    sc_array_reset(&me->peers);
    sc_array_reset(&me->peerreq);
#endif
    MY__FREE(me->ghost_rank);
    MY__FREE(me->ghost_element_nodes);
  }

  clean_construct(me);
  sc_array_reset(&me->construct);
  sc_array_reset(&me->ownsort);
  MY__FREE(me->chilev);
  MY__FREE(me->configuration);

  return tnodes;
}
