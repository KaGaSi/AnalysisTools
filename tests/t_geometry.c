/*
 * Geometry / box unit tests.
 *
 * Pins the pure numerical routines that the pair traversal and every PBC
 * calculation depend on: CalculateBoxData (both modes), the transform/inverse
 * matrices, RestorePBC, DistancePBC, and the fractional round-trip.
 */
#include "../src/AnalysisTools.h"
#include "test_util.h"
#include <stdlib.h>

static pcg32_random_t rng;
// uniform double in [0,1)
static double rand01(void) {
  return pcg32Rand0Int(&rng, 1u << 24) / (double)(1u << 24);
}
// uniform double in [lo,hi)
static double randr(double lo, double hi) {
  return lo + (hi - lo) * rand01();
}

// build a box from Length + angles (mode 0)
static BOX box_mode0(vec3d L, double al, double be, double ga) {
  BOX b = InitBox;
  b.Length = L;
  b.alpha = al;
  b.beta = be;
  b.gamma = ga;
  CHECK(CalculateBoxData(&b, 0));
  return b;
}

// mode 0 -> mode 1 must recover the original cell //{{{
static void test_box_roundtrip(void) {
  const struct { vec3d L; double al, be, ga; } cases[] = {
    { { .v = {10, 10, 10} }, 90, 90, 90 },   // orthogonal
    { { .v = {10, 12,  8} }, 90, 90, 60 },   // monoclinic (gamma tilt)
    { { .v = {13, 11,  9} }, 80, 95, 105 },  // fully triclinic
  };
  for (size_t t = 0; t < sizeof cases / sizeof *cases; t++) {
    BOX b = box_mode0(cases[t].L, cases[t].al, cases[t].be, cases[t].ga);
    // reconstruct from OrthoLength + the three tilt factors (mode 1)
    BOX b2 = InitBox;
    b2.OrthoLength = b.OrthoLength;
    b2.transform[0][1] = b.transform[0][1]; // xy
    b2.transform[0][2] = b.transform[0][2]; // xz
    b2.transform[1][2] = b.transform[1][2]; // yz
    CHECK(CalculateBoxData(&b2, 1));
    for (int dd = 0; dd < 3; dd++) {
      CHECK_CLOSE(b2.Length.v[dd], b.Length.v[dd], 1e-9);
    }
    CHECK_CLOSE(b2.alpha, b.alpha, 1e-7);
    CHECK_CLOSE(b2.beta,  b.beta,  1e-7);
    CHECK_CLOSE(b2.gamma, b.gamma, 1e-7);
    CHECK_CLOSE(b2.Volume, b.Volume, 1e-7 * b.Volume);
  }
} //}}}

// transform * inverse == I, and Volume == det(transform) //{{{
static void test_transform_inverse(void) {
  BOX b = box_mode0((vec3d){ .v = {13, 11, 9} }, 80, 95, 105);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      double s = 0;
      for (int k = 0; k < 3; k++) {
        s += b.transform[i][k] * b.inverse[k][j];
      }
      double expected = 0.0;
      if (i == j) {
        expected = 1.0;
      }
      CHECK_CLOSE(s, expected, 1e-9);
    }
  }
  // transform is upper-triangular, so det == product of the diagonal
  double det = b.transform[0][0] * b.transform[1][1] * b.transform[2][2];
  CHECK_CLOSE(det, b.Volume, 1e-7 * b.Volume);
} //}}}

// RestorePBC: in-range, idempotent, invariant under integer box shifts //{{{
static void test_restore_pbc(void) {
  vec3d L = { .v = {10, 12, 8} };
  for (int n = 0; n < 2000; n++) {
    vec3d c = { .v = { randr(-30, 30), randr(-30, 30), randr(-30, 30) } };
    vec3d w = RestorePBC(c, L);
    for (int dd = 0; dd < 3; dd++) {
      CHECK(w.v[dd] >= 0.0 && w.v[dd] < L.v[dd]);
    }
    // idempotent
    vec3d w2 = RestorePBC(w, L);
    for (int dd = 0; dd < 3; dd++) {
      CHECK_CLOSE(w2.v[dd], w.v[dd], 1e-12);
    }
    // shifting the input by an integer number of box lengths is invariant
    int k[3] = { pcg32RandIntInt(&rng, -3, 3), pcg32RandIntInt(&rng, -3, 3),
                 pcg32RandIntInt(&rng, -3, 3) };
    vec3d shifted = { .v = { c.v[0] + k[0] * L.v[0],
                             c.v[1] + k[1] * L.v[1],
                             c.v[2] + k[2] * L.v[2] } };
    vec3d ws = RestorePBC(shifted, L);
    for (int dd = 0; dd < 3; dd++) {
      CHECK_CLOSE(ws.v[dd], w.v[dd], 1e-8);
    }
  }
} //}}}

// DistancePBC on an orthogonal box == plain Distance; components <= L/2 //{{{
static void test_distance_orthogonal(void) {
  BOX b = box_mode0((vec3d){ .v = {10, 12, 8} }, 90, 90, 90);
  for (int n = 0; n < 2000; n++) {
    vec3d r1 = { .v = { randr(0, 10), randr(0, 12), randr(0, 8) } };
    vec3d r2 = { .v = { randr(0, 10), randr(0, 12), randr(0, 8) } };
    vec3d d  = DistancePBC(r1, r2, &b);
    vec3d dr = Distance(r1, r2, b.OrthoLength);
    for (int dd = 0; dd < 3; dd++) {
      CHECK_CLOSE(d.v[dd], dr.v[dd], 1e-12);
      CHECK(fabs(d.v[dd]) <= b.OrthoLength.v[dd] / 2 + 1e-9);
    }
    // symmetric in magnitude
    CHECK_CLOSE(VectLength(d), VectLength(DistancePBC(r2, r1, &b)), 1e-9);
  }
} //}}}

// DistancePBC minimum-image: invariant under lattice translations //{{{
static void test_distance_triclinic_minimum_image(void) {
  BOX b = box_mode0((vec3d){ .v = {13, 11, 9} }, 80, 95, 105);
  // lattice edge vectors are the columns of the transform matrix
  vec3d edge[3];
  for (int j = 0; j < 3; j++) {
    for (int i = 0; i < 3; i++) {
      edge[j].v[i] = b.transform[i][j];
    }
  }
  for (int n = 0; n < 2000; n++) {
    vec3d r1 = { .v = { randr(0, 13), randr(0, 11), randr(0, 9) } };
    vec3d r2 = { .v = { randr(0, 13), randr(0, 11), randr(0, 9) } };
    double ref = VectLength(DistancePBC(r1, r2, &b));
    // translate r2 by a random integer combination of lattice vectors
    int m[3] = { pcg32RandIntInt(&rng, -2, 2), pcg32RandIntInt(&rng, -2, 2),
                 pcg32RandIntInt(&rng, -2, 2) };
    vec3d r2t = r2;
    for (int dd = 0; dd < 3; dd++) {
      for (int j = 0; j < 3; j++) {
        r2t.v[dd] += m[j] * edge[j].v[dd];
      }
    }
    CHECK_CLOSE(VectLength(DistancePBC(r1, r2t, &b)), ref, 1e-7);
    // symmetric
    CHECK_CLOSE(VectLength(DistancePBC(r2, r1, &b)), ref, 1e-9);
  }
} //}}}

// CoorToFractional followed by CoorFromFractional is the identity //{{{
static void test_fractional_roundtrip(void) {
  BOX b = box_mode0((vec3d){ .v = {13, 11, 9} }, 80, 95, 105);
  int N = 300;
  SYSTEM s = (SYSTEM){0};
  s.Box = b;
  s.Count = InitCount;
  s.Count.Bead = N;
  s.Count.BeadCoor = N;
  s.Bead = malloc(N * sizeof *s.Bead);
  s.BeadCoor = malloc(N * sizeof *s.BeadCoor);
  CHECK(s.Bead && s.BeadCoor);
  vec3d *saved = malloc(N * sizeof *saved);
  for (int i = 0; i < N; i++) {
    InitBead(&s.Bead[i]);
    // place inside the parallelepiped: r = transform * s_frac
    double sf[3] = { rand01(), rand01(), rand01() };
    for (int dd = 0; dd < 3; dd++) {
      s.Bead[i].Position.v[dd] = b.transform[dd][0] * sf[0] +
                                 b.transform[dd][1] * sf[1] +
                                 b.transform[dd][2] * sf[2];
    }
    s.Bead[i].InTimestep = true;
    s.BeadCoor[i] = i;
    saved[i] = s.Bead[i].Position;
  }
  CoorToFractional(&s);
  CoorFromFractional(&s);
  for (int i = 0; i < N; i++) {
    for (int dd = 0; dd < 3; dd++) {
      CHECK_CLOSE(s.Bead[i].Position.v[dd], saved[i].v[dd], 1e-9);
    }
  }
  free(saved);
  free(s.Bead);
  free(s.BeadCoor);
} //}}}

// CosAngle / AngleDegrees: known angles, clamp, and degenerate zero vector //{{{
static void test_cos_angle(void) {
  vec3d x = { .v = {1, 0, 0} };
  vec3d y = { .v = {0, 1, 0} };
  vec3d z = { .v = {0, 0, 0} };
  // known angles
  CHECK_CLOSE(CosAngle(x, x),  1.0, 1e-12); // parallel
  CHECK_CLOSE(CosAngle(x, y),  0.0, 1e-12); // perpendicular
  CHECK_CLOSE(CosAngle(x, (vec3d){ .v = {-1, 0, 0} }), -1.0, 1e-12); // opposite
  CHECK_CLOSE(AngleDegrees(x, y), 90.0, 1e-9);
  CHECK_CLOSE(AngleDegrees(x, x),  0.0, 1e-9);
  // a zero-length vector has no defined angle -> NaN, not a div-by-zero inf
  CHECK(isnan(CosAngle(x, z)));
  CHECK(isnan(CosAngle(z, z)));
  CHECK(isnan(AngleDegrees(x, z)));
  // The cosine must stay within acos()'s domain even when fp round-off pushes
  // a (near-)parallel self-dot just past 1 (sqrt(s)*sqrt(s) < s happens for a
  // large fraction of vectors); without the clamp AngleDegrees(a, a) is NaN.
  for (int n = 0; n < 5000; n++) {
    vec3d a = { .v = { randr(-5, 5), randr(-5, 5), randr(-5, 5) } };
    // parallel / antiparallel: exactly the two edges of the domain
    CHECK(CosAngle(a, a) <= 1.0);
    CHECK(!isnan(AngleDegrees(a, a)));
    vec3d na = { .v = { -a.x, -a.y, -a.z } };
    CHECK(CosAngle(a, na) >= -1.0);
    CHECK(!isnan(AngleDegrees(a, na)));
    // an arbitrary second vector must land in range too
    vec3d b = { .v = { randr(-5, 5), randr(-5, 5), randr(-5, 5) } };
    double c = CosAngle(a, b);
    CHECK(c >= -1.0 && c <= 1.0);
    CHECK(!isnan(AngleDegrees(a, b)));
  }
} //}}}

int main(void) {
  pcg32Seed(&rng, 20260721u);
  RUN(test_cos_angle);
  RUN(test_box_roundtrip);
  RUN(test_transform_inverse);
  RUN(test_restore_pbc);
  RUN(test_distance_orthogonal);
  RUN(test_distance_triclinic_minimum_image);
  RUN(test_fractional_roundtrip);
  return test_main_end();
}
