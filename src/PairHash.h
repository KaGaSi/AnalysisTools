#ifndef PAIRHASH_H
#define PAIRHASH_H

#define _POSIX_C_SOURCE 200809L

#include <stdint.h>
#include "khash.h"

/*
 * Sparse hash map for counting intermolecular contacts (instead of the naive
 * approach of 2D array); stores only the non-zero contact counts in a hash map.
 *
 * Key encoding:
 * A pair (i, j) with i > j is packed into a single uint64_t:
 *   key = ((uint64_t) i << 32) | (uint64_t) j
 *
 * The convention i > j matches the triangular layout used in the rest of
 * Aggregates. Storing the two 32-bit indices side-by-side in one 64-bit word is
 * cheap to encode and decode and produces well-distributed hash keys.
 *
 * khash (klib, MIT licence): a header-only open-addressing hash table that is
 * widely used in bioinformatics / HPC code.  It is vendored as src/khash.h.
 * KHASH_MAP_INIT_INT64 registers a specialisation for uint64_t keys; the
 * generated functions are all static inline so including this header in
 * multiple translation units is safe.
 */

// Register the khash specialisation for uint64_t -> uint8_t (-c opt cap at 255)
KHASH_MAP_INIT_INT64(contacts, uint8_t)

/* Public type alias — callers only need to see "PairHash *". */
typedef khash_t(contacts) PairHash;

/*
 * PairHashKey() encodes a molecule pair as a hash key.
 *
 * Always call with i > j (larger index first) so that each unordered pair
 * has exactly one canonical key. The cast via unsigned int prevents sign
 * extension if mol indices are ever negative (shouldn't happen, but safe).
 */
static inline uint64_t PairHashKey(int i, int j) {
  return ((uint64_t)(unsigned int)i << 32) | (uint64_t)(unsigned int)j;
}

/*
 * PairHash_mol_i() / PairHash_mol_j()
 *
 * Decode the larger (i) and smaller (j) molecule indices back from a key.
 * Used in EvaluateContacts when iterating over all stored contacts.
 */
static inline int PairHash_mol_i(uint64_t key) {
  return (int)(uint32_t)(key >> 32);
}
static inline int PairHash_mol_j(uint64_t key) {
  return (int)(uint32_t)(key & 0xFFFFFFFFu);
}

// Functions implemented in PairHash.c
// Allocate and initialise an empty hash table
PairHash *PairHashAlloc(void);
/*
 * Increment the contact count for the pair (i, j), i > j. If the pair is not
 * yet in the table it is inserted with count 1. The count is capped at 255 (-c
 * option cap) to prevent overflow
 */
void PairHashIncrement(PairHash *h, int i, int j);
void PairHashFree(PairHash *h);

#endif /* PAIRHASH_H */
