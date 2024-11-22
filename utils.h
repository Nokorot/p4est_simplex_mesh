
#ifdef MEM_COUNT

#include "sc_containers.h"

typedef struct mem_counter {
  sc_mempool_t *mempool;
  sc_hash_t *mem_hash;
}
mem_counter_t;

extern mem_counter_t simplex_mem_counter;


typedef struct ptr_pair {
    const void *ptr;
    const char *name;
} ptr_pair_t;

void mem_counter_init(mem_counter_t *mc);
void mem_counter_deinit(mem_counter_t *mc);

void * mem_counter_add(mem_counter_t *mc, void *ptr, const char *name);

void * mem_counter_remove(mem_counter_t *mc, void *ptr);

void mem_counter_print_unalloced(mem_counter_t *mc);

// #ifdef P4EST_ENABLE_DEBUG


#define MY__ALLOC_N(t,n,name) (t *) mem_counter_add(&simplex_mem_counter, P4EST_ALLOC(t,n), (name) )
#define MY__ALLOC_ZERO_N(t,n,name) (t *) mem_counter_add(&simplex_mem_counter, P4EST_ALLOC_ZERO(t,n), name )


// #else
// #endif
//
#else

#define MY__ALLOC_N(t,n,name)       P4EST_ALLOC(t, n)
#define MY__ALLOC_ZERO_N(t,n,name)  P4EST_ALLOC_ZERO(t, n)

#define MY__FREE(ptr)               P4EST_FREE(ptr)

#endif

#define MY__ALLOC(t,n) MY__ALLOC_N(t, n, #t )
#define MY__ALLOC_ZERO(t,n) MY__ALLOC_ZERO_N(t,n, #t )
