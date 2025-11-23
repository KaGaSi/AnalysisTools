#include "pcg32.h"

// Syrex... TODO: explain

// Step function
uint32_t pcg32_random_r(pcg32_random_t *rng) {
  uint64_t oldstate = rng->state;
  rng->state = oldstate * 6364136223846793005ULL + (rng->inc | 1);
  uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
  uint32_t rot = oldstate >> 59u;
  return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}
// Seed
void pcg32Seed(pcg32_random_t *rng, uint64_t seed) {
  rng->state = 0U;
  // init sequence
  rng->inc = (seed << 1u) | 1u;
  pcg32_random_r(rng);
  // init state - flip bits via XOR operator
  rng->state += seed ^ 0xda3e39cb94b95bdbULL;
  pcg32_random_r(rng);
}

// Return uniform random int in [0, bound)
uint32_t pcg32Rand0Int(pcg32_random_t *rng, uint32_t bound) {
  uint32_t threshold = -bound % bound; // smallest number >= 2^32 % bound
  while (1) {
    uint32_t r = pcg32_random_r(rng);
    if (r >= threshold)
      return r % bound;
  }
}

// Uniform random integer in [a, b] inclusive
int32_t pcg32RandIntInt(pcg32_random_t *rng, int32_t a, int32_t b) {
  if (b < a) { // swap if user gives reversed bounds
    int32_t tmp = a;
    a = b;
    b = tmp;
  }
  uint32_t span = (uint32_t)(b - a) + 1u;
  uint32_t r = pcg32Rand0Int(rng, span);
  return a + (int32_t)r;
}
