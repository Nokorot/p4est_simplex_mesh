#include "p4est_base.h"
#include "sc.h"
#include "sc_statistics.h"
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
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

#include "utils.h"

static const char *log_dirpath = "local_logs";


typedef struct
{
  int minlevel, maxlevel;
  const char *conn;
  const char *mshpath;
  const char *outdir;

  int output_vtk;
  int log_to_file;

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
MY_p4est_lnodes_destroy (p4est_lnodes_t * lnodes)
{
  size_t              zz, count;
  p4est_lnodes_rank_t *lrank;

  MY__FREE (lnodes->element_nodes);
  MY__FREE (lnodes->nonlocal_nodes);
  MY__FREE (lnodes->global_owned_count);
  MY__FREE (lnodes->face_code);

  count = lnodes->sharers->elem_count;
  for (zz = 0; zz < count; zz++) {
    lrank = p4est_lnodes_rank_array_index (lnodes->sharers, zz);
    sc_array_reset (&(lrank->shared_nodes));
  }
  sc_array_destroy (lnodes->sharers);

  MY__FREE (lnodes);
}

static void
MY_p4est_tnodes_destroy (p4est_tnodes_t * tm)
{
  P4EST_ASSERT (tm != NULL);
  P4EST_ASSERT (tm->lnodes != NULL);

  if (tm->lnodes_owned) {
    MY_p4est_lnodes_destroy (tm->lnodes);
  }
  if (tm->coordinates != NULL) {
    sc_array_destroy (tm->coordinates);
  }
  if (tm->coord_to_lnode != NULL) {
    sc_array_destroy (tm->coord_to_lnode);
  }
  if (tm->simplices != NULL) {
    sc_array_destroy (tm->simplices);
  }
  if (tm->simplex_level != NULL) {
    sc_array_destroy (tm->simplex_level);
  }
  MY__FREE (tm->local_element_offset);
  MY__FREE (tm->local_element_level);
  MY__FREE (tm->local_tree_offset);
  MY__FREE (tm->local_tcount);
  MY__FREE (tm);
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
#endif


  TSTAT_BEGIN(g->ss, TIMINGS_ALL);
  tnodes = p4est_new_BECK_mesh(
        g->p4est,
        g->geom,
        g->ghost
#ifdef TIMINGS
        , g->ss
#endif
        );
  TSTAT_END(g->ss, TIMINGS_ALL);


#ifdef TIMINGS
  // statistics_print_tots(g->ss);

  statistics_finelize(g->ss);

  sc_stats_compute(g->mpicomm, NUM_COUNTERS, g->ss->counters);
  sc_stats_compute(g->mpicomm, TIMINGS_NUM_STATS, g->ss->timings);

  // sc_stats_print(p4est_package_id, SC_LP_INFO, TIMINGS_NUM_STATS,
  //                 g->ss->timings, 0, 0);


  sc_statistics_t *stats = g->ss->stats;
  // sc_statistics_print(stats, p4est_package_id, SC_LP_STATISTICS, 1, 0);

  size_t zz;
  sc_statinfo_t *stat;
  for (zz = 0; zz < stats->sarray->elem_count; ++zz) {
    stat = sc_array_index( stats->sarray, zz);
    P4EST_GLOBAL_STATISTICSF(
        "Total time for [%s]: %lf\n", stat->variable, stat->sum_values);
  }

  for (zz = 0; zz < TIMINGS_NUM_STATS; ++zz) {
    P4EST_GLOBAL_STATISTICSF("Total time for [%s]: %f\n", TIMEING_STAT_NAMES[zz],
        g->ss->timings[zz].sum_values);
  }

  for (zz = 0; zz < NUM_COUNTERS; ++zz) {
    P4EST_GLOBAL_STATISTICSF("Total count for [%s]: %lu\n", COUNTER_NAMES[zz],
        (uint64_t) (g->ss->counters[zz].sum_values));
  }


  if (g->mpirank == 0) {
    statistics_print_tots(g->ss);
  }



  sc_statistics_destroy(g->ss->stats);

  // timing_reset(g->ss);
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
  MY_p4est_tnodes_destroy (tnodes);
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

#ifdef MEM_COUNT
  mem_counter_init(&simplex_mem_counter);
#endif

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

  sc_options_add_bool (sc_opts, 'l', "",
            &opts->log_to_file, 0,
            "Log to file");

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
    // p4est_init (NULL, SC_LP_DEFAULT);

    char log_filepath[1024];
    sprintf(log_filepath, "%s/simplex_d%d_n%d_%s_l%d_L%d.log",
        log_dirpath, P4EST_DIM, g->mpisize, opts->conn, opts->minlevel, opts->maxlevel);

    FILE *log_file = fopen(log_filepath, "w");
    if (log_file == NULL) {
        fprintf(stderr, "Error opening log file\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    sc_set_log_defaults(log_file, NULL, SC_LP_DEFAULT);

    if (create_context_from_opts(g, opts)) {
      exitcode = 1;
      break;
    }

    tnodes_run(g);

    if (opts->output_vtk)
      p4est_write_vtk(g);

    fclose(log_file);


  } while (0);
  if (g->errmsg) {
    P4EST_GLOBAL_LERROR (g->errmsg);
  }

  /*** clean up and exit ***/
  destroy_context(g);

#ifdef MEM_COUNT
  mem_counter_print_unalloced(&simplex_mem_counter);
  mem_counter_deinit(&simplex_mem_counter);
#endif

  sc_options_destroy (sc_opts);
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return exitcode ? EXIT_FAILURE : EXIT_SUCCESS;
}
