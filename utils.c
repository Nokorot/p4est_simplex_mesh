
#ifndef MEM_COUNT
#define MEM_COUNT
#endif

#include <sc_containers.h>
#include "utils.h"

mem_counter_t simplex_mem_counter;


static unsigned int
ptr_pair_hash(const void *a, const void *u)
{
  long key = (long) ((ptr_pair_t *) a)->ptr;

  key = (~key) + (key << 18); // key = (key << 18) - key - 1;
  key = key ^ (key >> 31);
  key = key * 21; // key = (key + (key << 2)) + (key << 4);
  key = key ^ (key >> 11);
  key = key + (key << 6);
  key = key ^ (key >> 22);

  return key;
}

static int
ptr_pair_equal(const void *a, const void *b, const void *u)
{
  return ((ptr_pair_t *) a)->ptr == ((ptr_pair_t *) b)->ptr;
}

void
mem_counter_init(mem_counter_t *mc)
{
  mc->mempool = sc_mempool_new(sizeof(ptr_pair_t));
  mc->mem_hash = sc_hash_new(ptr_pair_hash, ptr_pair_equal, NULL, NULL);
}

void
mem_counter_deinit(mem_counter_t *mc)
{
  sc_hash_destroy(mc->mem_hash);
  sc_mempool_destroy(mc->mempool);
}

void *
mem_counter_add(mem_counter_t *mc, void *ptr, const char *name)
{
  int added;
  ptr_pair_t **new, ins = {ptr, name};

  added = sc_hash_insert_unique(mc->mem_hash, &ins, (void ***) &new);
  SC_ASSERT(added);

  // SC_INFOF("[%s] (%p) was loaded\n", name, ptr);

  *new = sc_mempool_alloc(mc->mempool);

  (*new)->ptr = ptr;
  (*new)->name = name;

  return ptr;
}


void *
mem_counter_remove(mem_counter_t *mc, void *ptr)
{
  ptr_pair_t *pair, ins = {ptr, NULL};
  sc_hash_remove(mc->mem_hash, &ins, (void **) &pair);
  // SC_INFOF("[%s] (%p) was unloaded\n", pair->name, pair->ptr, pair);

  return ptr;
}


static int
ptr_pair_foreach(void **v, const void *u)
{
  ptr_pair_t *pair = *v;

  SC_LERRORF("[%s] (%p) is loaded\n", pair->name, pair->ptr);
  return 1;
}

void
mem_counter_print_unalloced(mem_counter_t *mc)
{
  sc_hash_foreach(mc->mem_hash, ptr_pair_foreach);
}

