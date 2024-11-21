#include "sc.h"
#include "sc_statistics.h"
#include <stdint.h>
#ifdef VIM_LS
// #include <p4est_to_p8est.h>
#define TIMINGS
#endif

#include <string.h>
#include <sc_options.h>

#ifndef P4_TO_P8
#include "p4est_lnodes.h"
#include "p4est_tnodes.h"
#include <p4est_bits.h>
#include <p4est_extended.h>
#include <p4est_vtk.h>

#include "p4est_simplex_mesh.h"
#else
#include "p8est_lnodes.h"
#include "p8est_tnodes.h"
#include <p8est_bits.h>
#include <p8est_extended.h>
#include <p8est_vtk.h>

#include "p8est_simplex_mesh.h"
#endif

#ifdef TIMINGS
#include "statistics.h"
#endif

#include "utils.c"

typedef struct
{
  int minlevel, maxlevel;
  const char *conn;
  const char *mshpath;
  const char *outdir;

  int output_vtk;

  double vtk_scale;
}
options_t;

typedef struct
{
  sc_MPI_Comm         mpicomm;
  int                 mpisize, mpirank;

  options_t *opts;

  p4est_t *p4est;
  p4est_connectivity_t *conn;
  p4est_ghost_t *ghost;
  p4est_geometry_t *geom;

  char *errmsg;

#ifdef TIMINGS
  statistics_t *ss;
#endif


}
context_t;

/** Refines the p4est, producing hanging-nodes
 */
static int
refine_fn (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t *quadrant)
{
  context_t *g = (context_t *) p4est->user_pointer;


  if (which_tree == 1)
    return 1;
  return 0;


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

static int
create_context_from_opts(context_t *g, options_t *opts)
{
  p4est_t *p4est;
  p4est_connectivity_t *conn;
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
  g->ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FULL);

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

static void
destroy_context(context_t *g)
{
  if (g->p4est)
    p4est_destroy(g->p4est);
  if (g->conn)
    p4est_connectivity_destroy(g->conn);
  if (g->ghost)
    p4est_ghost_destroy(g->ghost);
  if (g->geom)
    p4est_geometry_destroy(g->geom);
}

static void
p4est_write_vtk(context_t *g)
{
  char filepath[1024];

  sprintf(filepath, "out/p4est_box_d%d_%s_%d-%d",
      P4EST_DIM, g->opts->conn, g->opts->minlevel, g->opts->maxlevel);

  p4est_vtk_write_file(g->p4est, g->geom, filepath);
}

static void
tnodes_run(context_t *g)
{
  p4est_tnodes_t *tnodes;
  char filepath[1024];

  int retval;
  p4est_vtk_context_t *cont;

  P4EST_ASSERT(g->p4est);
  P4EST_ASSERT(g->ghost);

#ifdef TIMINGS
  g->ss = SC_ALLOC(statistics_t, 1);
  timing_init(g->ss);

  timing_interval_begin(g->ss, TIMINGS_ALL);
#endif

  tnodes = p4est_new_BECK_mesh(
        g->p4est,
        g->geom,
        g->ghost
#ifdef TIMINGS
        , g->ss
#endif
        );

#ifdef TIMINGS
  timing_interval_end(g->ss, TIMINGS_ALL);
  statistics_finelize(g->ss);

  sc_stats_compute(g->mpicomm, NUM_COUNTERS, g->ss->counters);
  sc_stats_compute(g->mpicomm, TIMINGS_NUM_STATS, g->ss->timings);

  if (g->mpirank == 0) {
    size_t stat;
    for (stat = 0; stat < TIMINGS_NUM_STATS; ++stat) {
      P4EST_INFOF("TimeTital[%s]: %f\n", TIMEING_STAT_NAMES[stat],
          g->ss->timings[stat].sum_values);
    }

    for (stat = 0; stat < NUM_COUNTERS; ++stat) {
      P4EST_INFOF("Counter[%s]: %lu\n", COUNTER_NAMES[stat],
          (uint64_t) (g->ss->counters[stat].sum_values));
    }

  }

  timing_reset(g->ss);
  SC_FREE(g->ss);
#endif


  if (g->opts->output_vtk) {
    sprintf(filepath, "%s/" P4EST_STRING "_snodes_simplices_%s_%d-%d",
        g->opts->outdir, g->opts->conn,
        g->opts->minlevel, g->opts->maxlevel);

    /* write VTK output */
    /* the geometry was passed to the tnodes already, don't use it here */
    cont = p4est_vtk_context_new (g->p4est, filepath);
    SC_CHECK_ABORT (cont != NULL, "Open VTK context");
    // p4est_vtk_context_set_geom (cont, g->geom);
    p4est_vtk_context_set_continuous (cont, g->opts->vtk_scale >= 1.0);

    /* beware: values < 1. cause a lot more mesh nodes */
    p4est_vtk_context_set_scale (cont, g->opts->vtk_scale);

    cont = p4est_vtk_write_header_tnodes (cont, tnodes);
    SC_CHECK_ABORT (cont != NULL, "Write tnodes VTK header");
    // TODO: Add levels to snodes
    cont = p4est_vtk_write_cell_dataf (cont, 1, 0, 1, 0,
                                        0, 0, cont);

    SC_CHECK_ABORT (cont != NULL, "Write tnodes VTK cells");
    retval = p4est_vtk_write_footer (cont);
    SC_CHECK_ABORT (!retval, "Close VTK context");
  }

  /* free triangle mesh */
  p4est_lnodes_destroy (tnodes->lnodes);
  p4est_tnodes_destroy (tnodes);
}

int
main(int argc, char **argv) {
  int mpiret;
  int first_argc;
  int exitcode;
  sc_options_t *sc_opts;
  options_t options, *opts = &options;
  context_t context = {0}, *g = &context;
  // p4est_simplex_mesh_t *smesh;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  g->mpicomm = sc_MPI_COMM_WORLD;
  mpiret = sc_MPI_Comm_size (g->mpicomm, &g->mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (g->mpicomm, &g->mpirank);
  SC_CHECK_MPI (mpiret);

  /* Does this really need to be before the parsing? */
  sc_init (g->mpicomm, 1, 1, NULL, SC_LP_DEFAULT);

  sc_opts = sc_options_new (argv[0]);
  sc_options_add_int (sc_opts, 'l', "minlevel",
              &opts->minlevel, 2, "Lowest level");
  sc_options_add_int (sc_opts, 'L', "maxlevel",
              &opts->maxlevel, 5, "Highest level");
  sc_options_add_double (sc_opts, 's', "scale",
              &opts->vtk_scale, .95, "Parameter to shrink the simplices");

  sc_options_add_string (sc_opts, 'c', "conn",
            &opts->conn, "unit", "Name of the connectivity");
  sc_options_add_string (sc_opts, 'o', "outdir",
            &opts->outdir, "out",
            "Directory in which to output output the vtk files");


  sc_options_add_bool (sc_opts, 'v', "",
            &opts->output_vtk, 0,
            "If true, the resulting mesh is written to a vtk file");

  exitcode = 0;
  do {
    first_argc = sc_options_parse(p4est_package_id,
                    SC_LP_DEFAULT, sc_opts, argc, argv);
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
    p4est_init (NULL, SC_LP_DEFAULT);

    if (create_context_from_opts(g, opts)) {
      exitcode = 1;
      break;
    }

    tnodes_run(g);

    if (opts->output_vtk)
      p4est_write_vtk(g);


  } while (0);
  if (g->errmsg) {
    P4EST_GLOBAL_LERROR (g->errmsg);
  }

  /*** clean up and exit ***/
  destroy_context(g);

  sc_options_destroy (sc_opts);
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return exitcode ? EXIT_FAILURE : EXIT_SUCCESS;
}
