#include <string.h>
#ifdef VIM_LS
// #include <p4est_to_p8est.h>
#endif

#include <sc_options.h>

#ifndef P4_TO_P8
#include <p4est_bits.h>
#include <p4est_extended.h>
#include <p4est_vtk.h>

#include "p4est_simplex_mesh.h"
#else
#include <p8est_bits.h>
#include <p8est_extended.h>
#include <p8est_vtk.h>

#include "p8est_simplex_mesh.h"
#endif

#include "utils.c"

typedef struct
{
  int minlevel, maxlevel;
  const char *conn;
  const char *mshpath;

  int hangine_node_level;
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

/** Refines the p4est, producing hanging-nodes
 */
static int refine_fn (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t *quadrant)
{
  context_t *g = (context_t *) p4est->user_pointer;

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

  /* Initialize default option values */
  opts->minlevel = 2;
  opts->maxlevel = 5;
  opts->conn = "unit";
  opts->mshpath = NULL;  // Defaults to "{opts->conn}.msh"

  sc_opts = sc_options_new (argv[0]);
  sc_options_add_int (sc_opts, 'l', "minlevel", &opts->minlevel, opts->minlevel, "Lowest level");
  sc_options_add_int (sc_opts, 'L', "maxlevel", &opts->maxlevel, opts->maxlevel, "Highest level");

  sc_options_add_int (sc_opts, 'h', "hanging-node-level", &opts->hangine_node_level, opts->hangine_node_level,
      "This determines what hanging nodes to include. h=0,1,2\n"
      // " h=-1: Do not add element centre nodes either,\n"
      " h=0: Include no hanging nodes,\n"
      " h=1: Include hanging face nodes, but not edge nodes,\n"
      " h=2: Include all hanging nodes. (default)");
  // sc_options_add_string (sc_opts, 'v', "vtk", &opts->vtk, opts->vtk, "VTK basename");
  sc_options_add_string (sc_opts, 'm', "msh", &opts->mshpath, opts->mshpath, "MSH base filename");
  sc_options_add_string (sc_opts, 'c', "conn", &opts->conn, opts->conn, "Name of the connectivity");


  smesh = NULL;

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

    g->opts = opts;
    // sc_init (g->mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
    p4est_init (NULL, SC_LP_DEFAULT);

    if (create_context_from_opts(g, opts)) {
      exitcode = 1;
      break;
    }

    smesh = p4est_new_simplex_mesh(
        g->p4est,
        g->geom,
        g->ghost_layer);

    if (!smesh)
      break;

#ifndef P4_TO_P8
    sprintf(filepath, "out/box_2d_%s_%d-%d", opts->conn, opts->minlevel, opts->maxlevel);
#else
    sprintf(filepath, "out/box_3d_%s_%d-%d", opts->conn, opts->minlevel, opts->maxlevel);
#endif

    p4est_vtk_write_file(g->p4est, g->geom, filepath);

    if (opts->mshpath == NULL) {
#ifndef P4_TO_P8
      sprintf(filepath, "out/sx_2d_%s_%d-%d_%d.vtk", opts->conn, opts->minlevel, opts->maxlevel, g->mpirank);
#else
      sprintf(filepath, "out/sx_3d_%s_%d-%d_%d.vtk", opts->conn, opts->minlevel, opts->maxlevel, g->mpirank);
#endif
    } else {
      strcpy(filepath, opts->mshpath);
    }

    p4est_simplex_mesh_write_vtk_file(filepath, smesh);

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
