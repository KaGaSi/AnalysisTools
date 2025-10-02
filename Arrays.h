#ifndef ARRAYS_H
#define ARRAYS_H

#define _POSIX_C_SOURCE 200809L

#include <stdlib.h>
#include <stdio.h>

// TODO: comment & explain

// generic N-D array
typedef struct {
  size_t ndim;     // number of dimensions
  size_t *shape;   // array of dimension sizes
  size_t *stride;  // array of strides
  double *data;    // flat storage
} ArrND;
// allocate new N-D array
ArrND NewArrND(size_t ndim, const size_t *shape);
// free N-D array
void FreeArrND(ArrND arr);
// 2D helper
static inline size_t idx2d(const ArrND arr, size_t i, size_t j) {
  return i * arr.stride[0] + j * arr.stride[1];
}
// 3D helper
static inline size_t idx3d(const ArrND arr, size_t i, size_t j, size_t k) {
  return i * arr.stride[0] + j * arr.stride[1] + k * arr.stride[2];
}
// 4D helper
static inline size_t idx4d(const ArrND arr, size_t i, size_t j,
                           size_t k, size_t l) {
  return i * arr.stride[0] + j * arr.stride[1] +
         k * arr.stride[2] + l * arr.stride[3];
}

#endif
