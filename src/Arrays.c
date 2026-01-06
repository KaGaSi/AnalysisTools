#include "Arrays.h"
#include <string.h>

// TODO: EXPLAIN!!!

// Helpers //{{{
static size_t *CalcArrNDStride(size_t ndim, const size_t *shape) {
  size_t *stride = malloc(ndim * sizeof(size_t));
  if (!stride) {
    return NULL;
  }
  stride[ndim - 1] = 1;
  for (ssize_t d = ndim - 2; d >= 0; d--) {
    stride[d] = stride[d+1] * shape[d+1];
  }
  return stride;
}
static size_t CalcArrNDTotalSize(size_t ndim, const size_t *shape) {
  size_t total = 1;
  for (size_t d = 0; d < ndim; d++) {
    total *= shape[d];
  }
  return total;
}
static size_t CalcArrNDOffset(size_t ndim, const size_t *idx,
                              const size_t *stride) {
  size_t off = 0;
  for (size_t d = 0; d < ndim; d++) {
    off += idx[d] * stride[d];
  }
  return off;
}
static int InitNDBase(size_t ndim, const size_t *shape,
                      size_t **shape_out, size_t **stride_out) {
  *shape_out = malloc(ndim * sizeof(size_t));
  if (!*shape_out) {
    return -1;
  }
  for (size_t d = 0; d < ndim; d++) {
    (*shape_out)[d] = shape[d];
  }
  *stride_out = CalcArrNDStride(ndim, shape);
  if (!*stride_out) {
    free(*shape_out);
    return -1;
  }
  return 0;
} //}}}
// Constructors //{{{
ArrNDd *CreateArrNDd(size_t ndim, const size_t *shape) {
  ArrNDd *a = malloc(sizeof *a);
  if (!a) {
    return NULL;
  }
  a->ndim = ndim;
  if (InitNDBase(ndim, shape, &a->shape, &a->stride) != 0) {
    free(a);
    return NULL;
  }
  size_t total = CalcArrNDTotalSize(ndim, shape);
  a->d = calloc(total, sizeof(double));
  if (!a->d) {
    free(a->shape);
    free(a->stride);
    free(a);
    return NULL;
  }
  return a;
}
ArrNDld *CreateArrNDld(size_t ndim, const size_t *shape) {
  ArrNDld *a = malloc(sizeof *a);
  if (!a) {
    return NULL;
  }
  a->ndim = ndim;
  if (InitNDBase(ndim, shape, &a->shape, &a->stride) != 0) {
    free(a);
    return NULL;
  }
  size_t total = CalcArrNDTotalSize(ndim, shape);
  a->d = calloc(total, sizeof(long double));
  if (!a->d) {
    free(a->shape);
    free(a->stride);
    free(a);
    return NULL;
  }
  return a;
}
ArrNDi *CreateArrNDi(size_t ndim, const size_t *shape) {
  ArrNDi *a = malloc(sizeof *a);
  if (!a) {
    return NULL;
  }
  a->ndim = ndim;
  if (InitNDBase(ndim, shape, &a->shape, &a->stride) != 0) {
    free(a);
    return NULL;
  }
  size_t total = CalcArrNDTotalSize(ndim, shape);
  a->d = calloc(total, sizeof(int));
  if (!a->d) {
    free(a->shape);
    free(a->stride);
    free(a);
    return NULL;
  }
  return a;
}
ArrNDli *CreateArrNDli(size_t ndim, const size_t *shape) {
  ArrNDli *a = malloc(sizeof *a);
  if (!a) {
    return NULL;
  }
  a->ndim = ndim;
  if (InitNDBase(ndim, shape, &a->shape, &a->stride) != 0) {
    free(a);
    return NULL;
  }
  size_t total = CalcArrNDTotalSize(ndim, shape);
  a->d = calloc(total, sizeof(long double));
  if (!a->d) {
    free(a->shape);
    free(a->stride);
    free(a);
    return NULL;
  }
  return a;
}
ArrNDb *CreateArrNDb(size_t ndim, const size_t *shape) {
  ArrNDb *a = malloc(sizeof *a);
  if (!a) {
    return NULL;
  }
  a->ndim = ndim;
  if (InitNDBase(ndim, shape, &a->shape, &a->stride) != 0) {
    free(a);
    return NULL;
  }
  size_t total = CalcArrNDTotalSize(ndim, shape);
  a->d = calloc(total, sizeof(bool));
  if (!a->d) {
    free(a->shape);
    free(a->stride);
    free(a);
    return NULL;
  }
  return a;
}
//}}}
// Destructors //{{{
static void FreeBaseND(void *d, size_t *shape, size_t *stride) {
  free(shape);
  free(stride);
  free(d);
}
void FreeArrNDd(ArrNDd *a) {
  if (!a) {
    return;
  }
  FreeBaseND(a->d, a->shape, a->stride);
  free(a);
}
void FreeArrNDld(ArrNDld *a) {
  if (!a) {
    return;
  }
  FreeBaseND(a->d, a->shape, a->stride);
  free(a);
}
void FreeArrNDi(ArrNDi *a) {
  if (!a) {
    return;
  }
  FreeBaseND(a->d, a->shape, a->stride);
  free(a);
}
void FreeArrNDli(ArrNDli *a) {
  if (!a) {
    return;
  }
  FreeBaseND(a->d, a->shape, a->stride);
  free(a);
}
void FreeArrNDb(ArrNDb *a) {
  if (!a) {
    return;
  }
  FreeBaseND(a->d, a->shape, a->stride);
  free(a);
} //}}}
// Fillers - fill all array elements by given value //{{{
void ArrND_double_fill(ArrNDd *a, double v) {
  size_t total = CalcArrNDTotalSize(a->ndim, a->shape);
  for (size_t i = 0; i < total; i++) {
    a->d[i] = v;
  }
}
void ArrND_longdouble_fill(ArrNDld *a, long double v) {
  size_t total = CalcArrNDTotalSize(a->ndim, a->shape);
  for (size_t i = 0; i < total; i++) {
    a->d[i] = v;
  }
}
void ArrND_int_fill(ArrNDi *a, int v) {
  size_t total = CalcArrNDTotalSize(a->ndim, a->shape);
  for (size_t i = 0; i < total; i++) {
    a->d[i] = v;
  }
}
void ArrND_longint_fill(ArrNDli *a, long int v) {
  size_t total = CalcArrNDTotalSize(a->ndim, a->shape);
  for (size_t i = 0; i < total; i++) {
    a->d[i] = v;
  }
}
void ArrND_bool_fill(ArrNDb *a, bool v) {
  size_t total = CalcArrNDTotalSize(a->ndim, a->shape);
  for (size_t i = 0; i < total; i++) {
    a->d[i] = v;
  }
} //}}}
// Setters - set one element to given value //{{{
void ArrND_double_set(ArrNDd *a, const size_t *idx, double v) {
  size_t off = CalcArrNDOffset(a->ndim, idx, a->stride);
  a->d[off] = v;
}
void ArrND_longdouble_set(ArrNDld *a, const size_t *idx, long double v) {
  size_t off = CalcArrNDOffset(a->ndim, idx, a->stride);
  a->d[off] = v;
}
void ArrND_int_set(ArrNDi *a, const size_t *idx, int v) {
  size_t off = CalcArrNDOffset(a->ndim, idx, a->stride);
  a->d[off] = v;
}
void ArrND_longint_set(ArrNDli *a, const size_t *idx, long int v) {
  size_t off = CalcArrNDOffset(a->ndim, idx, a->stride);
  a->d[off] = v;
}
void ArrND_bool_set(ArrNDb *a, const size_t *idx, bool v) {
  size_t off = CalcArrNDOffset(a->ndim, idx, a->stride);
  a->d[off] = v;
} //}}}
// Adders - increment array element by given value //{{{
void ArrND_double_add(ArrNDd *a, const size_t *idx, double v) {
  size_t off = CalcArrNDOffset(a->ndim, idx, a->stride);
  a->d[off] += v;
}
void ArrND_longdouble_add(ArrNDld *a, const size_t *idx, long double v) {
  size_t off = CalcArrNDOffset(a->ndim, idx, a->stride);
  a->d[off] += v;
}
void ArrND_int_add(ArrNDi *a, const size_t *idx, int v) {
  size_t off = CalcArrNDOffset(a->ndim, idx, a->stride);
  a->d[off] += v;
}
void ArrND_longint_add(ArrNDli *a, const size_t *idx, long int v) {
  size_t off = CalcArrNDOffset(a->ndim, idx, a->stride);
  a->d[off] += v;
} //}}}
// Getters - get one array element //{{{
double ArrND_double_get(const ArrNDd *a, const size_t *idx) {
  size_t off = CalcArrNDOffset(a->ndim, idx, a->stride);
  return a->d[off];
}
long double ArrND_longdouble_get(const ArrNDld *a, const size_t *idx) {
  size_t off = CalcArrNDOffset(a->ndim, idx, a->stride);
  return a->d[off];
}
int ArrND_int_get(const ArrNDi *a, const size_t *idx) {
  size_t off = CalcArrNDOffset(a->ndim, idx, a->stride);
  return a->d[off];
}
long int ArrND_longint_get(const ArrNDli *a, const size_t *idx) {
  size_t off = CalcArrNDOffset(a->ndim, idx, a->stride);
  return a->d[off];
}
bool ArrND_bool_get(const ArrNDb *a, const size_t *idx) {
  size_t off = CalcArrNDOffset(a->ndim, idx, a->stride);
  return a->d[off];
}
 //}}}
