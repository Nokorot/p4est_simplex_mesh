

#ifndef TIMINGS
#define TIMINGS
#include "p4est_base.h"
#endif

#include "sc.h"
#include "sc_flops.h"
#include "sc_statistics.h"

#include "p4est.h"
#include "statistics.h"
#include <alloca.h>
#include <string.h>

const char * TIMEING_STAT_NAMES[TIMINGS_NUM_STATS];
const char * COUNTER_NAMES[NUM_COUNTERS];

void
timing_init(statistics_t *ss)
{
  memset(ss->counters, 0, NUM_COUNTERS*sizeof(sc_statinfo_t));
  memset(ss->timings, 0, TIMINGS_NUM_STATS*sizeof(sc_statinfo_t));

  TIMEING_STAT_NAMES[TIMINGS_ALL]           = "Total Time";
  TIMEING_STAT_NAMES[TIMINGS_ITET_STEP]     = "Iteration Step";
  TIMEING_STAT_NAMES[TIMINGS_ITER_VOLUMES]  = "Iter Volume";
  TIMEING_STAT_NAMES[TIMINGS_ITER_FACES]    = "Iter Face";
  TIMEING_STAT_NAMES[TIMINGS_ITER_EDGES]    = "Iter Edge";
  TIMEING_STAT_NAMES[TIMINGS_ITER_CORNERS]  = "Iter Corner";

  COUNTER_NAMES[COUNTER_HANGING_FACE] = "Hanging Faces";
  COUNTER_NAMES[COUNTER_HANGING_EDGE] = "Hanging Edges";
  COUNTER_NAMES[COUNTER_SHARED_NODES] = "Shared Nodes";

  sc_flops_start(&ss->fi);
}

void
counter_inc(statistics_t *ss, int count, int counter)
{
  ss->counters[counter].sum_values += count;
}

void
statistics_finelize(statistics_t *ss)
{
  size_t stat;
  for (stat = 0; stat < TIMINGS_NUM_STATS; ++stat) {
    sc_stats_set1(&ss->timings[stat], ss->timings[stat].sum_values,
        TIMEING_STAT_NAMES[stat]);
  }

  for (stat = 0; stat < NUM_COUNTERS; ++stat) {
    sc_stats_set1(&ss->counters[stat], ss->counters[stat].sum_values,
        COUNTER_NAMES[stat]);
  }
}


void
statistics_print_tots(statistics_t *ss)
{
  size_t stat;
  for (stat = 0; stat < TIMINGS_NUM_STATS; ++stat) {
    P4EST_GLOBAL_STATISTICSF("Tital time for  %s : %f\n", TIMEING_STAT_NAMES[stat],
        ss->timings[stat].sum_values);
  }

  for (stat = 0; stat < NUM_COUNTERS; ++stat) {
    P4EST_INFOF("Counter[%s]: %lu\n", COUNTER_NAMES[stat],
        (uint64_t) (ss->counters[stat].sum_values));
  }
}

void
timing_reset(statistics_t *ss)
{
  SC_NOOP();
}

inline void
timing_interval_begin(statistics_t *ss, int stat)
{
  sc_flops_snap(&ss->fi, &ss->snapshots[stat]);
}

inline void
timing_interval_end(statistics_t *ss, int stat)
{
  sc_flopinfo_t *snap = &ss->snapshots[stat];
  sc_flops_shot(&ss->fi, snap);
  ss->timings[stat].sum_values += snap->iwtime;

  // sc_stats_accumulate(&ss->stats[stat], snap->iwtime);
  //
}
