#ifndef TIMINGS_H_
#define TIMINGS_H_

#ifndef TIMINGS

#define TSTAT_BEGIN(ss, stat)
#define TSTAT_END(ss, stat)

// #define TIME_CALL(ss, name, cmd) cmd

#else

#define TSTAT_BEGIN(ss, stat) timing_interval_begin(ss, stat)
#define TSTAT_END(ss, stat) timing_interval_end(ss, stat)

// #define TIME_CALL(ss, name, cmd) cmd

#include <sc_statistics.h>
#include <sc_flops.h>
#include <time.h>

enum {
  TIMINGS_ALL,
  TIMINGS_ITET_STEP,
  TIMINGS_ITER_VOLUMES,
  TIMINGS_ITER_FACES,
  TIMINGS_ITER_EDGES,
  TIMINGS_ITER_CORNERS,

  TIMINGS_NUM_STATS,
};

enum {
  COUNTER_HANGING_FACE,
  COUNTER_HANGING_EDGE,
  COUNTER_SHARED_NODES,
  NUM_COUNTERS,
};

extern const char * TIMEING_STAT_NAMES[TIMINGS_NUM_STATS];
extern const char * COUNTER_NAMES[NUM_COUNTERS];

typedef struct timing_state
{
  int mpirank;

  sc_statinfo_t counters[NUM_COUNTERS];
  sc_statinfo_t timings[TIMINGS_NUM_STATS];
  sc_flopinfo_t fi, snap, snapshots[TIMINGS_NUM_STATS];
}
statistics_t;


void counter_inc(statistics_t *ss, int count, int couner);


void timing_init(statistics_t *ss);
void statistics_finelize(statistics_t *ss);

void statistics_print_tots(statistics_t *ss);
void timing_reset(statistics_t *ss);

void timing_interval_begin(statistics_t *ss, int stat);
void timing_interval_end(statistics_t *ss, int stat);

#endif // TIMINGS
#endif // TIMINGS_H_
