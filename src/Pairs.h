#ifndef PAIRS_H
#define PAIRS_H

#define _POSIX_C_SOURCE 200809L

#include "Structs.h"

// callback function to call for each pair of particles
// i & j ... indices in System.BeadCoor array
typedef void (*pair_cb_t)(int i, int j, const SYSTEM System, void *userdata);
// callback function for checking if bead/mol should be used
// i ... index in System.BeadCoor array
typedef bool (*check_cb_t)(int i, const SYSTEM System, void *userdata);
// use linked list to traverse bead pairs
void TraverseLinkedListPairs(const SYSTEM System, const double cell_size,
                             pair_cb_t pair_callback, void *pair_ud,
                             check_cb_t check_callback, void *check_ud);
// use brute way to traverse bead pairs
void TraverseBrutePairs(const SYSTEM System,
                        pair_cb_t pair_callback, void *pair_ud,
                        check_cb_t check_callback, void *check_ud);
// use linked list traversal if cell_size is not too large, brute otherwise
void TraversePairs(const SYSTEM System, const double cell_size,
                   pair_cb_t pair_callback, void *pair_ud,
                   check_cb_t check_callback, void *check_ud);
#endif
