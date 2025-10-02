#include "Arrays.h"
// TODO: comment & explain

// general ND indexer
static inline size_t idND(const ArrND arr, const size_t *indices);

// // allocate new N-D array //{{{
// ArrND *new_ArrND(size_t ndim, const size_t *shape) {
//   ArrND *arr = malloc(sizeof(ArrND));
//   if (!arr) {
//     return NULL;
//   }
//
//   arr->ndim = ndim;
//   arr->shape = calloc(ndim, sizeof(size_t));
//   arr->stride = calloc(ndim, sizeof(size_t));
//   if (!arr->shape || !arr->stride) {
//     free(arr->shape);
//     free(arr->stride);
//     free(arr);
//     return NULL;
//   }
//
//   size_t total = 1;
//   for (size_t i = 0; i < ndim; i++) {
//     arr->shape[i] = shape[i];
//     total *= shape[i];
//   }
//
//   // compute strides (row-major)
//   arr->stride[ndim - 1] = 1;
//   for (size_t i = ndim - 1; i-- > 0; )
//     arr->stride[i] = arr->stride[i + 1] * arr->shape[i + 1];
//
//   arr->data = calloc(total, sizeof(double));
//   if (!arr->data) {
//     free(arr->shape);
//     free(arr->stride);
//     free(arr);
//     return NULL;
//   }
//
//   return arr;
// } //}}}
// allocate new N-D array //{{{
ArrND NewArrND(size_t ndim, const size_t *shape) {
  ArrND arr;

  arr.ndim = ndim;
  arr.shape = calloc(ndim, sizeof(size_t));
  arr.stride = calloc(ndim, sizeof(size_t));
  if (!arr.shape || !arr.stride) {
    free(arr.shape);
    free(arr.stride);
    fprintf(stderr, "ERROR CREATING N-D ARRAY (.shape/.stride bit)!");
    exit(1);
  }

  size_t total = 1;
  for (size_t i = 0; i < ndim; i++) {
    arr.shape[i] = shape[i];
    total *= shape[i];
  }

  // compute strides (row-major)
  arr.stride[ndim - 1] = 1;
  for (size_t i = ndim - 1; i-- > 0; )
    arr.stride[i] = arr.stride[i + 1] * arr.shape[i + 1];

  arr.data = calloc(total, sizeof(double));
  if (!arr.data) {
    free(arr.shape);
    free(arr.stride);
    fprintf(stderr, "ERROR CREATING N-D ARRAY (.data bit)!");
    exit(1);
  }

  return arr;
} //}}}
// free N-D array //{{{
void FreeArrND(ArrND arr) {
  free(arr.data);
  free(arr.shape);
  free(arr.stride);
} //}}}
// general ND indexer //{{{
static inline size_t idND(const ArrND arr, const size_t *indices) {
  size_t offset = 0;
  for (size_t d = 0; d < arr.ndim; d++)
    offset += indices[d] * arr.stride[d];
  return offset;
} //}}}
