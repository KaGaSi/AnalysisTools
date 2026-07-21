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
/*
 * 2D variant: cells are built only in the two axes perpendicular to
 * norm_axis (0/1/2), so beads any distance apart along norm_axis still meet
 * if their in-plane separation is within cell_size. Use for slit-like (2D)
 * systems where distances are measured in-plane only.
 */
void TraverseLinkedListPairs2D(const SYSTEM System, const double cell_size,
                               const int norm_axis,
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
// ditto, but with 2D cells (see TraverseLinkedListPairs2D)
void TraversePairs2D(const SYSTEM System, const double cell_size,
                     const int norm_axis,
                     pair_cb_t pair_callback, void *pair_ud,
                     check_cb_t check_callback, void *check_ud);
#endif
