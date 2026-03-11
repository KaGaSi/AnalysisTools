#include "PairHash.h"
#include "Errors.h"   /* err_msg / ErrorAlloc conventions used elsewhere */

#include <stdlib.h>

/*
 * PairHashAlloc()
 *
 * kh_init() allocates the khash bookkeeping structure (bucket arrays, etc.)
 * and returns a pointer.  The table starts completely empty; khash will
 * resize it automatically as entries are added.
 */
PairHash *PairHashAlloc(void) {
  PairHash *h = kh_init(contacts);
  if (!h) {
    /*
     * kh_init() returns NULL only on a failed malloc.  Use the same
     * error-reporting style as the rest of the codebase.
     */
    err_msg("PairHashAlloc: allocation failed");
    exit(1);
  }
  return h;
}

/*
 * PairHashIncrement()
 *
 * Look up the pair (i, j) in the table:
 *   - if it is absent, insert it with count = 1;
 *   - if it is present, increment its count by 1 (capped at UINT8_MAX=255).
 *
 * khash API recap:
 *   kh_get(name, h, key) ... returns an iterator (== kh_end if absent)
 *   kh_put(name, h, key, &ret) ... insert key, return iterator;
 *                                  *ret: 1 = inserted, 0 = already present,
 *                                       -1 = error
 *   kh_val(h, iter) ... lvalue for the value at iterator position
 */
void PairHashIncrement(PairHash *h, int i, int j) {
  uint64_t key = PairHashKey(i, j);

  /* Look for the pair in the table. */
  khiter_t it = kh_get(contacts, h, key);

  if (it == kh_end(h)) {
    // Pair not yet seen -> insert with count 1
    int absent; // kh_put sets this to 1 (new), 0 (existed), or -1 (error)
    it = kh_put(contacts, h, key, &absent);
    if (absent < 0) {
      err_msg("PairHashIncrement: kh_put failed (out of memory?)");
      exit(1);
    }
    kh_val(h, it) = 1;
  } else {
    // Pair already recorded -> increment, but don't overflow uint8_t.
    if (kh_val(h, it) < 255) { // cap at 255 as max value (-c option)
      kh_val(h, it)++;
    }
  }
}

// kh_destroy() frees the bucket arrays and the khash_t struct itself
void PairHashFree(PairHash *h) {
  kh_destroy(contacts, h);
}
