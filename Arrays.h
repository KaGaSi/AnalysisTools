#ifndef ARRAYS_H
#define ARRAYS_H

#define _POSIX_C_SOURCE 200809L

#include <stdlib.h>
#include <stddef.h>
#include <unistd.h>

// TODO: EXPLAIN!!!

// ArrND structs //{{{
typedef struct {
  size_t ndim;
  size_t *shape;
  size_t *stride;
  double *d;
} ArrNDd;
typedef struct {
  size_t ndim;
  size_t *shape;
  size_t *stride;
  long double *d;
} ArrNDld;
typedef struct {
  size_t ndim;
  size_t *shape;
  size_t *stride;
  int *d;
} ArrNDi;
typedef struct {
  size_t ndim;
  size_t *shape;
  size_t *stride;
  long int *d;
} ArrNDli; //}}}
// constructors //{{{
// general
ArrNDd  *CreateArrNDd(size_t ndim, const size_t *shape);
ArrNDld *CreateArrNDld(size_t ndim, const size_t *shape);
ArrNDi  *CreateArrNDi(size_t ndim, const size_t *shape);
ArrNDli *CreateArrNDli(size_t ndim, const size_t *shape);
// 2-dimensional constructors
static inline ArrNDd *CreateArr2Dd(const size_t x, const size_t y) {
  size_t shape2D[] = {x, y};
  return CreateArrNDd(2, shape2D);
}
static inline const ArrNDld *CreateArr2Dld(const size_t x, const size_t y) {
  size_t shape2D[] = {x, y};
  return CreateArrNDld(2, shape2D);
}
static inline ArrNDi *CreateArr2Di(const size_t x, const size_t y) {
  size_t shape2D[] = {x, y};
  return CreateArrNDi(2, shape2D);
}
static inline ArrNDli *CreateArr2Dli(const size_t x, const size_t y) {
  size_t shape2D[] = {x, y};
  return CreateArrNDli(2, shape2D);
}
// 3-dimensional constructors
static inline ArrNDd *CreateArr3Dd(const size_t x, const size_t y, const size_t z) {
  size_t shape3D[] = {x, y, z};
  return CreateArrNDd(3, shape3D);
}
static inline const ArrNDld *CreateArr3Dld(const size_t x, const size_t y,
                                           const size_t z) {
  size_t shape3D[] = {x, y, z};
  return CreateArrNDld(3, shape3D);
}
static inline ArrNDi *CreateArr3Di(const size_t x, const size_t y,
                                   const size_t z) {
  size_t shape3D[] = {x, y, z};
  return CreateArrNDi(3, shape3D);
}
static inline ArrNDli *CreateArr3Dli(const size_t x, const size_t y,
                                     const size_t z) {
  size_t shape3D[] = {x, y, z};
  return CreateArrNDli(3, shape3D);
} //}}}
// destructors //{{{
void FreeArrNDd(ArrNDd *a);
void FreeArrNDld(ArrNDld *a);
void FreeArrNDi(ArrNDi *a);
void FreeArrNDli(ArrNDli *a);
#define FreeArrND(a) \
  _Generic((a), \
           ArrNDd*: FreeArrNDd, \
           ArrNDld*: FreeArrNDld, \
           ArrNDi*: FreeArrNDi, \
           ArrNDli*: FreeArrNDli)(a)
 //}}}
// fillers //{{{
void ArrND_double_fill(ArrNDd *a, double val);
void ArrND_longdouble_fill(ArrNDld *a, long double val);
void ArrND_int_fill(ArrNDi *a, int val);
void ArrND_longint_fill(ArrNDli *a, long int val);
#define FillArrND(a, val) \
  _Generic((a), \
           ArrNDd*: ArrND_double_fill, \
           ArrNDld*: ArrND_longdouble_fill, \
           ArrNDi*: ArrND_int_fill, \
           ArrNDli*: ArrND_longint_fill)(a, val)
 //}}}
// general setters //{{{
void ArrND_double_set(ArrNDd *a, const size_t *idx, double val);
void ArrND_longdouble_set(ArrNDld *a, const size_t *idx, long double val);
void ArrND_int_set(ArrNDi *a, const size_t *idx, int val);
void ArrND_longint_set(ArrNDli *a, const size_t *idx, long int val);
#define SetArrND(a, idx, val) \
  _Generic((a), \
           ArrNDd*: ArrND_double_set, \
           ArrNDld*: ArrND_longdouble_set, \
           ArrNDi*: ArrND_int_set, \
           ArrNDli*: ArrND_longint_set)(a, idx, val)
 //}}}
// 2-dimensional setters //{{{
static inline void Arr2D_double_set(ArrNDd *a, const size_t i,
                                    const size_t j, double val) {
  size_t index[] = {i, j};
  ArrND_double_set(a, index, val);
}
static inline void Arr2D_longdouble_set(ArrNDld *a, const size_t i,
                                        const size_t j, long double val) {
  size_t index[] = {i, j};
  ArrND_longdouble_set(a, index, val);
}
static inline void Arr2D_int_set(ArrNDi *a, const size_t i,
                                 const size_t j, int val) {
  size_t index[] = {i, j};
  ArrND_int_set(a, index, val);
}
static inline void Arr2D_longint_set(ArrNDli *a, const size_t i,
                                     const size_t j, long int val) {
  size_t index[] = {i, j};
  ArrND_longint_set(a, index, val);
}
#define SetArr2D(a, i, j, val) \
  _Generic((a), \
           ArrNDd*: Arr2D_double_set, \
           ArrNDld*: Arr2D_longdouble_set, \
           ArrNDi*: Arr2D_int_set, \
           ArrNDli*: Arr2D_longint_set)(a, i, j, val)
 //}}}
// 3-dimensional setters //{{{
static inline void Arr3D_double_set(ArrNDd *a, const size_t i, const size_t j,
                                    const size_t k, double val) {
  size_t index[] = {i, j, k};
  ArrND_double_set(a, index, val);
}
static inline void Arr3D_longdouble_set(ArrNDld *a, const size_t i, const size_t j,
                                        const size_t k, long double val) {
  size_t index[] = {i, j, k};
  ArrND_longdouble_set(a, index, val);
}
static inline void Arr3D_int_set(ArrNDi *a, const size_t i, const size_t j,
                                 const size_t k, int val) {
  size_t index[] = {i, j, k};
  ArrND_int_set(a, index, val);
}
static inline void Arr3D_longint_set(ArrNDli *a, const size_t i, const size_t j,
                                     const size_t k, long int val) {
  size_t index[] = {i, j, k};
  ArrND_longint_set(a, index, val);
}
#define SetArr3D(a, i, j, k, val) \
  _Generic((a), \
           ArrNDd*: Arr3D_double_set, \
           ArrNDld*: Arr3D_longdouble_set, \
           ArrNDi*: Arr3D_int_set, \
           ArrNDli*: Arr3D_longint_set)(a, i, j, k, val)
 //}}}
// general adders //{{{
void ArrND_double_add(ArrNDd *a, const size_t *idx, double val);
void ArrND_longdouble_add(ArrNDld *a, const size_t *idx, long double val);
void ArrND_int_add(ArrNDi *a, const size_t *idx, int val);
void ArrND_longint_add(ArrNDli *a, const size_t *idx, long int val);
#define AddArrND(a, idx, val) \
  _Generic((a), \
           ArrNDd*: ArrND_double_add, \
           ArrNDld*: ArrND_longdouble_add, \
           ArrNDi*: ArrND_int_add, \
           ArrNDli*: ArrND_longint_add)(a, idx, val)
 //}}}
// 2-dimensional adders //{{{
static inline void Arr2D_double_add(ArrNDd *a, const size_t i,
                                    const size_t j, double val) {
  size_t index[] = {i, j};
  ArrND_double_add(a, index, val);
}
static inline void Arr2D_longdouble_add(ArrNDld *a, const size_t i,
                                        const size_t j, long double val) {
  size_t index[] = {i, j};
  ArrND_longdouble_add(a, index, val);
}
static inline void Arr2D_int_add(ArrNDi *a, const size_t i,
                                 const size_t j, int val) {
  size_t index[] = {i, j};
  ArrND_int_add(a, index, val);
}
static inline void Arr2D_longint_add(ArrNDli *a, const size_t i,
                                     const size_t j, long int val) {
  size_t index[] = {i, j};
  ArrND_longint_add(a, index, val);
}
#define AddArr2D(a, i, j, val) \
  _Generic((a), \
           ArrNDd*: Arr2D_double_add, \
           ArrNDld*: Arr2D_longdouble_add, \
           ArrNDi*: Arr2D_int_add, \
           ArrNDli*: Arr2D_longint_add)(a, i, j, val)
 //}}}
// 3-dimensional adders //{{{
static inline void Arr3D_double_add(ArrNDd *a, const size_t i, const size_t j,
                                    const size_t k, double val) {
  size_t index[] = {i, j, k};
  ArrND_double_add(a, index, val);
}
static inline void Arr3D_longdouble_add(ArrNDld *a, const size_t i,
                                        const size_t j, const size_t k,
                                        long double val) {
  size_t index[] = {i, j, k};
  ArrND_longdouble_add(a, index, val);
}
static inline void Arr3D_int_add(ArrNDi *a, const size_t i, const size_t j,
                                 const size_t k, int val) {
  size_t index[] = {i, j, k};
  ArrND_int_add(a, index, val);
}
static inline void Arr3D_longint_add(ArrNDli *a, const size_t i, const size_t j,
                                     const size_t k, long int val) {
  size_t index[] = {i, j, k};
  ArrND_longint_add(a, index, val);
}
#define AddArr3D(a, i, j, k, val) \
  _Generic((a), \
           ArrNDd*: Arr3D_double_add, \
           ArrNDld*: Arr3D_longdouble_add, \
           ArrNDi*: Arr3D_int_add, \
           ArrNDli*: Arr3D_longint_add)(a, i, j, k, val)
 //}}}
// general getters //{{{
double ArrND_double_get(const ArrNDd *a, const size_t *idx);
long double ArrND_longdouble_get(const ArrNDld *a, const size_t *idx);
int ArrND_int_get(const ArrNDi *a, const size_t *idx);
long int ArrND_longint_get(const ArrNDli *a, const size_t *idx);
#define GetArrND(a, idx) \
  _Generic((a), \
           ArrNDd*: ArrND_double_get, \
           const ArrNDd*: ArrND_double_get, \
           ArrNDld*: ArrND_longdouble_get, \
           const ArrNDld*: ArrND_longdouble_get, \
           ArrNDi*: ArrND_int_get, \
           const ArrNDi*: ArrND_int_get, \
           ArrNDli*: ArrND_longint_get, \
           const ArrNDli*: ArrND_longint_get)(a, idx)
 //}}}
// 2-dimensional getters //{{{
static inline double Arr2D_double_get(const ArrNDd *a,
                                      const size_t i, const size_t j) {
  size_t index[] = {i, j};
  return ArrND_double_get(a, index);
}
static inline long double Arr2D_longdouble_get(const ArrNDld *a,
                                               const size_t i, const size_t j) {
  size_t index[] = {i, j};
  return ArrND_longdouble_get(a, index);
}
static inline int Arr2D_int_get(const ArrNDi *a,
                                const size_t i, const size_t j) {
  size_t index[] = {i, j};
  return ArrND_int_get(a, index);
}
static inline long int Arr2D_longint_get(const ArrNDli *a,
                                         const size_t i, const size_t j) {
  size_t index[] = {i, j};
  return ArrND_longint_get(a, index);
}
#define GetArr2D(a, i, j) \
  _Generic((a), \
           ArrNDd*: Arr2D_double_get, \
           const ArrNDd*: Arr2D_double_get, \
           ArrNDld*: Arr2D_longdouble_get, \
           const ArrNDld*: Arr2D_longdouble_get, \
           ArrNDi*: Arr2D_int_get, \
           const ArrNDi*: Arr2D_int_get, \
           ArrNDli*: Arr2D_longint_get, \
           const ArrNDli*: Arr2D_longint_get)(a, i, j)
 //}}}
// 3-dimensional getters //{{{
static inline double Arr3D_double_get(const ArrNDd *a, const size_t i,
                                      const size_t j, const size_t k) {
  size_t index[] = {i, j, k};
  return ArrND_double_get(a, index);
}
static inline long double Arr3D_longdouble_get(const ArrNDld *a, const size_t i,
                                               const size_t j, const size_t k) {
  size_t index[] = {i, j, k};
  return ArrND_longdouble_get(a, index);
}
static inline int Arr3D_int_get(const ArrNDi *a, const size_t i,
                                const size_t j, const size_t k) {
  size_t index[] = {i, j, k};
  return ArrND_int_get(a, index);
}
static inline long int Arr3D_longint_get(const ArrNDli *a, const size_t i,
                                         const size_t j, const size_t k) {
  size_t index[] = {i, j, k};
  return ArrND_longint_get(a, index);
}
#define GetArr3D(a, i, j, k) \
  _Generic((a), \
           ArrNDd*: Arr3D_double_get, \
           const ArrNDd*: Arr3D_double_get, \
           ArrNDld*: Arr3D_longdouble_get, \
           const ArrNDld*: Arr3D_longdouble_get, \
           ArrNDi*: Arr3D_int_get, \
           const ArrNDi*: Arr3D_int_get, \
           ArrNDli*: Arr3D_longint_get, \
           const ArrNDli*: Arr3D_longint_get)(a, i, j, k)
 //}}}

#endif
