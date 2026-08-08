/*
 * ArrND container tests.
 *
 * Every utility that accumulates data (density profiles, pair correlations,
 * per-type tables) stores it in an ArrND, so a stride bug here corrupts
 * results everywhere without crashing. The shapes used are deliberately
 * non-cubic and non-square: a transposed or mis-computed stride still passes
 * on a cubic shape, because the strides happen to be equal.
 *
 * The layout is asserted twice over - once through the public accessors, and
 * once against the flat backing buffer - which pins row-major order as part
 * of the contract rather than an accident of the implementation.
 */
#include "../src/AnalysisTools.h"
#include "test_util.h"
#include <stdlib.h>

// ---- layout ---------------------------------------------------------------

// 2D: every element addressable, and the flat layout is row-major //{{{
static void test_2d_layout(void) {
  const size_t nx = 3, ny = 7; // deliberately different extents
  ArrNDd *a = CreateArr2Dd(nx, ny);
  CHECK(a != nullptr);
  CHECK(a->ndim == 2);
  CHECK(a->shape[0] == nx && a->shape[1] == ny);
  // row-major: the last axis is contiguous
  CHECK(a->stride[1] == 1);
  CHECK(a->stride[0] == ny);
  // write a value that encodes its own index
  for (size_t i = 0; i < nx; i++) {
    for (size_t j = 0; j < ny; j++) {
      SetArr2D(a, i, j, (double)(i * 100 + j));
    }
  }
  // read every one back through the accessor
  for (size_t i = 0; i < nx; i++) {
    for (size_t j = 0; j < ny; j++) {
      CHECK_CLOSE(GetArr2D(a, i, j), (double)(i * 100 + j), 0);
    }
  }
  // and directly against the flat buffer
  for (size_t i = 0; i < nx; i++) {
    for (size_t j = 0; j < ny; j++) {
      CHECK_CLOSE(a->d[i * ny + j], (double)(i * 100 + j), 0);
    }
  }
  FreeArrND(a);
} //}}}
// 3D: same, with three distinct extents //{{{
static void test_3d_layout(void) {
  const size_t nx = 2, ny = 3, nz = 5;
  ArrNDi *a = CreateArr3Di(nx, ny, nz);
  CHECK(a != nullptr);
  CHECK(a->ndim == 3);
  CHECK(a->stride[2] == 1);
  CHECK(a->stride[1] == nz);
  CHECK(a->stride[0] == ny * nz);
  for (size_t i = 0; i < nx; i++) {
    for (size_t j = 0; j < ny; j++) {
      for (size_t k = 0; k < nz; k++) {
        SetArr3D(a, i, j, k, (int)(i * 100 + j * 10 + k));
      }
    }
  }
  for (size_t i = 0; i < nx; i++) {
    for (size_t j = 0; j < ny; j++) {
      for (size_t k = 0; k < nz; k++) {
        CHECK(GetArr3D(a, i, j, k) == (int)(i * 100 + j * 10 + k));
        CHECK(a->d[i * ny * nz + j * nz + k] == (int)(i * 100 + j * 10 + k));
      }
    }
  }
  FreeArrND(a);
} //}}}
// writing one element must not disturb its neighbours //{{{
/*
 * A stride that is too large or too small usually still round-trips the
 * element just written; what it breaks is everything next to it.
 */
static void test_no_neighbour_corruption(void) {
  const size_t nx = 4, ny = 6, nz = 3;
  ArrNDd *a = CreateArr3Dd(nx, ny, nz);
  FillArrND(a, -1.0);
  // set exactly one element
  SetArr3D(a, 2, 4, 1, 42.0);
  size_t total = nx * ny * nz;
  size_t hits = 0;
  for (size_t f = 0; f < total; f++) {
    if (a->d[f] != -1.0) {
      hits++;
      CHECK_CLOSE(a->d[f], 42.0, 0);
      // it must sit exactly where row-major says it should
      CHECK(f == 2 * ny * nz + 4 * nz + 1);
    }
  }
  CHECK(hits == 1);
  FreeArrND(a);
} //}}}
// arbitrary rank through the general accessors //{{{
static void test_nd_layout(void) {
  size_t shape[4] = {2, 3, 2, 4};
  ArrNDli *a = CreateArrNDli(4, shape);
  CHECK(a != nullptr);
  CHECK(a->stride[3] == 1);
  CHECK(a->stride[2] == 4);
  CHECK(a->stride[1] == 2 * 4);
  CHECK(a->stride[0] == 3 * 2 * 4);
  for (size_t i = 0; i < shape[0]; i++) {
    for (size_t j = 0; j < shape[1]; j++) {
      for (size_t k = 0; k < shape[2]; k++) {
        for (size_t l = 0; l < shape[3]; l++) {
          size_t idx[4] = {i, j, k, l};
          SetArrND(a, idx, (long)(((i * 3 + j) * 2 + k) * 4 + l));
        }
      }
    }
  }
  // the encoded values are exactly 0..total-1 in flat order
  size_t total = shape[0] * shape[1] * shape[2] * shape[3];
  for (size_t f = 0; f < total; f++) {
    CHECK(a->d[f] == (long)f);
  }
  FreeArrND(a);
} //}}}
// degenerate shapes: rank 1, and an axis of extent 1 //{{{
static void test_edge_shapes(void) {
  // one-dimensional
  size_t shape1[1] = {5};
  ArrNDd *a = CreateArrNDd(1, shape1);
  CHECK(a != nullptr);
  CHECK(a->ndim == 1);
  CHECK(a->stride[0] == 1);
  for (size_t i = 0; i < 5; i++) {
    size_t idx[1] = {i};
    SetArrND(a, idx, (double)i * 2);
  }
  for (size_t i = 0; i < 5; i++) {
    size_t idx[1] = {i};
    CHECK_CLOSE(GetArrND(a, idx), (double)i * 2, 0);
  }
  FreeArrND(a);

  // an axis of extent 1 collapses a stride; the others must still be right
  ArrNDd *b = CreateArr3Dd(1, 4, 3);
  CHECK(b->stride[0] == 4 * 3);
  CHECK(b->stride[1] == 3);
  CHECK(b->stride[2] == 1);
  SetArr3D(b, 0, 3, 2, 9.0);
  CHECK_CLOSE(GetArr3D(b, 0, 3, 2), 9.0, 0);
  CHECK_CLOSE(b->d[3 * 3 + 2], 9.0, 0);
  FreeArrND(b);

  ArrNDd *c = CreateArr3Dd(4, 1, 3);
  CHECK(c->stride[0] == 1 * 3);
  CHECK(c->stride[1] == 3);
  CHECK(c->stride[2] == 1);
  SetArr3D(c, 3, 0, 2, 8.0);
  CHECK_CLOSE(c->d[3 * 3 + 2], 8.0, 0);
  FreeArrND(c);
} //}}}

// ---- element operations ---------------------------------------------------

// a fresh array is zeroed (the accumulating utilities rely on it) //{{{
static void test_zero_initialised(void) {
  ArrNDd *d = CreateArr2Dd(3, 4);
  ArrNDi *i = CreateArr2Di(3, 4);
  ArrNDli *li = CreateArr2Dli(3, 4);
  ArrNDb *b = CreateArr2Db(3, 4);
  ArrNDld *ld = CreateArr2Dld(3, 4);
  for (size_t x = 0; x < 3; x++) {
    for (size_t y = 0; y < 4; y++) {
      CHECK_CLOSE(GetArr2D(d, x, y), 0, 0);
      CHECK(Arr2D_int_get(i, x, y) == 0);
      CHECK(Arr2D_longint_get(li, x, y) == 0);
      CHECK(Arr2D_bool_get(b, x, y) == false);
      CHECK((double)Arr2D_longdouble_get(ld, x, y) == 0);
    }
  }
  FreeArrND(d);
  FreeArrND(i);
  FreeArrND(li);
  FreeArrND(b);
  FreeArrND(ld);
} //}}}
// fill touches every element and nothing beyond //{{{
static void test_fill(void) {
  const size_t nx = 5, ny = 2, nz = 3;
  ArrNDd *a = CreateArr3Dd(nx, ny, nz);
  FillArrND(a, 7.5);
  size_t total = nx * ny * nz;
  for (size_t f = 0; f < total; f++) {
    CHECK_CLOSE(a->d[f], 7.5, 0);
  }
  // refilling overwrites rather than accumulates
  FillArrND(a, -2.0);
  for (size_t f = 0; f < total; f++) {
    CHECK_CLOSE(a->d[f], -2.0, 0);
  }
  FreeArrND(a);

  // bool arrays have their own filler
  ArrNDb *b = CreateArr2Db(3, 3);
  FillArrND(b, true);
  for (size_t x = 0; x < 3; x++) {
    for (size_t y = 0; y < 3; y++) {
      CHECK(Arr2D_bool_get(b, x, y) == true);
    }
  }
  FreeArrND(b);
} //}}}
// add accumulates in place, which is how the histograms are built //{{{
static void test_add_accumulates(void) {
  ArrNDd *a = CreateArr2Dd(3, 4);
  for (int rep = 0; rep < 5; rep++) {
    AddArr2D(a, 1, 2, 1.5);
  }
  CHECK_CLOSE(GetArr2D(a, 1, 2), 7.5, 1e-12);
  // neighbours untouched
  CHECK_CLOSE(GetArr2D(a, 1, 1), 0, 0);
  CHECK_CLOSE(GetArr2D(a, 1, 3), 0, 0);
  CHECK_CLOSE(GetArr2D(a, 0, 2), 0, 0);
  CHECK_CLOSE(GetArr2D(a, 2, 2), 0, 0);
  // adding a negative value subtracts
  AddArr2D(a, 1, 2, -7.5);
  CHECK_CLOSE(GetArr2D(a, 1, 2), 0, 1e-12);
  FreeArrND(a);

  ArrNDli *li = CreateArr3Dli(2, 2, 2);
  for (int rep = 0; rep < 3; rep++) {
    AddArr3D(li, 1, 0, 1, 1000000L);
  }
  CHECK(Arr3D_longint_get(li, 1, 0, 1) == 3000000L);
  FreeArrND(li);
} //}}}
// every element type stores its own range without truncation //{{{
/*
 * The _Generic dispatch picks the accessor from the pointer type; a wrong
 * entry in the macro would silently narrow the value.
 */
static void test_type_ranges(void) {
  ArrNDd *d = CreateArr2Dd(2, 2);
  SetArr2D(d, 0, 0, 0.1);
  CHECK_CLOSE(GetArr2D(d, 0, 0), 0.1, 1e-15); // not truncated to an int
  FreeArrND(d);

  ArrNDli *li = CreateArr2Dli(2, 2);
  long big = 4000000000L; // beyond a 32-bit int
  Arr2D_longint_set(li, 1, 1, big);
  CHECK(Arr2D_longint_get(li, 1, 1) == big);
  FreeArrND(li);

  ArrNDi *i = CreateArr2Di(2, 2);
  Arr2D_int_set(i, 0, 1, -12345);
  CHECK(Arr2D_int_get(i, 0, 1) == -12345);
  FreeArrND(i);

  ArrNDb *b = CreateArr2Db(2, 2);
  Arr2D_bool_set(b, 0, 1, true);
  CHECK(Arr2D_bool_get(b, 0, 1) == true);
  CHECK(Arr2D_bool_get(b, 0, 0) == false);
  Arr2D_bool_set(b, 0, 1, false);
  CHECK(Arr2D_bool_get(b, 0, 1) == false);
  FreeArrND(b);
} //}}}
// the getters accept a const pointer (the _Generic table lists both) //{{{
static void test_const_getters(void) {
  ArrNDd *a = CreateArr2Dd(2, 3);
  SetArr2D(a, 1, 2, 3.25);
  const ArrNDd *ca = a;
  CHECK_CLOSE(GetArr2D(ca, 1, 2), 3.25, 0);
  size_t idx[2] = {1, 2};
  CHECK_CLOSE(GetArrND(ca, idx), 3.25, 0);
  FreeArrND(a);
} //}}}
// freeing a null pointer is a no-op rather than a crash //{{{
static void test_free_null(void) {
  ArrNDd *d = nullptr;
  ArrNDi *i = nullptr;
  ArrNDli *li = nullptr;
  ArrNDb *b = nullptr;
  ArrNDld *ld = nullptr;
  FreeArrND(d);
  FreeArrND(i);
  FreeArrND(li);
  FreeArrND(b);
  FreeArrND(ld);
  CHECK(true); // reaching here without a segfault is the assertion
} //}}}
// a large allocation is fully addressable (ASan catches an undersized calloc) //{{{
static void test_large_allocation(void) {
  const size_t nx = 40, ny = 30, nz = 20;
  ArrNDd *a = CreateArr3Dd(nx, ny, nz);
  CHECK(a != nullptr);
  // touch the corners, including the very last element
  SetArr3D(a, 0, 0, 0, 1.0);
  SetArr3D(a, nx - 1, ny - 1, nz - 1, 2.0);
  CHECK_CLOSE(GetArr3D(a, 0, 0, 0), 1.0, 0);
  CHECK_CLOSE(GetArr3D(a, nx - 1, ny - 1, nz - 1), 2.0, 0);
  CHECK_CLOSE(a->d[nx * ny * nz - 1], 2.0, 0);
  FreeArrND(a);
} //}}}

int main(void) {
  RUN(test_2d_layout);
  RUN(test_3d_layout);
  RUN(test_no_neighbour_corruption);
  RUN(test_nd_layout);
  RUN(test_edge_shapes);
  RUN(test_zero_initialised);
  RUN(test_fill);
  RUN(test_add_accumulates);
  RUN(test_type_ranges);
  RUN(test_const_getters);
  RUN(test_free_null);
  RUN(test_large_allocation);
  return test_main_end();
}
