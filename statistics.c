
#include "sc.h"
#include "sc_flops.h"
#include "sc_statistics.h"

#include "statistics.h"
#include <alloca.h>
#include <string.h>

const char * TIMEING_STAT_NAMES[TIMINGS_NUM_STATS] = {
  "Total Time",
  "iter Volume",
  "iter Face",
  "iter Edge",
  "iter Corner",
};

const char * COUNTER_NAMES[NUM_COUNTERS] = {
  "Hanging Faces",
  "Hanging Edges",
  "Shared Nodes",
};


void
timing_init(statistics_t *ss)
{
  memset(ss->counters, 0, NUM_COUNTERS*sizeof(sc_statinfo_t));
  memset(ss->timings, 0, TIMINGS_NUM_STATS*sizeof(sc_statinfo_t));

  // sc_stats_init(&ss->stats[TIMINGS_ALL], "ALL");

  // sc_stats_init(&ss->stats[TIMINGS_ITER_CORNERS], "Itter Corners");
  // sc_stats_init(&ss->stats[TIMINGS_ITER_VOLUMES], "Itter Volume");
  // sc_stats_init(&ss->stats[TIMINGS_ITER_FACES], "Itter Faces");
  // sc_stats_init(&ss->stats[TIMINGS_ITER_EDGES], "Itter Edges");
  ss->fi = SC_ALLOC(sc_flopinfo_t, 1);
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

  for (stat = 0; stat < NUM_COUNTERS; ++stat) {
    sc_stats_set1(&ss->counters[stat],
                  ss->counters[stat].sum_values,
                  COUNTER_NAMES[stat]);
  }

  for (stat = 0; stat < TIMINGS_NUM_STATS; ++stat) {
    sc_stats_set1(&ss->timings[stat],
                  ss->timings[stat].sum_values,
                  TIMEING_STAT_NAMES[stat]);
  }
}

void
timing_reset(statistics_t *ss)
{
  SC_FREE(ss->fi);
}


inline void
timing_interval_begin(statistics_t *ss, int stat)
{
  sc_flops_snap(ss->fi, &ss->snapshots[stat]);
}

inline void
timing_interval_end(statistics_t *ss, int stat)
{
  sc_flopinfo_t *snap = &ss->snapshots[stat];
  sc_flops_shot(ss->fi, snap);
  ss->timings[stat].sum_values += snap->iwtime;

  // sc_stats_accumulate(&ss->stats[stat], snap->iwtime);
}
