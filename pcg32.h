#ifndef pcg32_H
#define pcg32_H

#define _POSIX_C_SOURCE 200809L
#include <stdint.h>

// State structure
typedef struct {
  uint64_t state;
  uint64_t inc;
} pcg32_random_t;

// rng seed
void pcg32Seed(pcg32_random_t *rng, uint64_t seed);
// generate random number in interval <0, bound) (i.e., 0-inclusive)
uint32_t pcg32Rand0Int(pcg32_random_t *rng, uint32_t bound);
// Uniform random integer in [a, b] (i.e., a- and b-inclusive)
int32_t pcg32RandIntInt(pcg32_random_t *rng, int32_t a, int32_t b);

#endif
