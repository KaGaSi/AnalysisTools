/*
 * Differential test for the pair traversal.
 *
 * The cell-linked-list traversal is only an acceleration structure: for any
 * system it must enumerate exactly the same set of within-cutoff pairs as the
 * O(N^2) brute force, with no misses and no double counting. We check that
 * equivalence over random orthogonal AND triclinic boxes (3D), plus the 2D
 * slit path. The triclinic case is the regression guard for the fractional
 * cell-binning fix in Pairs.c.
 */
#include "../src/AnalysisTools.h"
#include "test_util.h"
#include <stdlib.h>
#include <string.h>

static pcg32_random_t rng;
static double rand01(void) {
  return pcg32Rand0Int(&rng, 1u << 24) / (double)(1u << 24);
}

// shared state for the recording callback //{{{
struct rec {
  const BOX *box;
  double cutoff;
  int N;         // number of beads (== Count.BeadCoor)
  int norm_axis; // -1 = full 3D, else the non-binned (slit-normal) axis
  long *count;   // N*N matrix of hit counts, indexed [lo*N + hi], lo < hi
};

// in-plane minimum-image distance for the 2D (slit) path (orthogonal box)
static double inplane_dist(vec3d a, vec3d b, const BOX *box, int norm) {
  double d2 = 0;
  for (int dd = 0; dd < 3; dd++) {
    if (dd == norm) {
      continue;
    }
    double delta = a.v[dd] - b.v[dd];
    double L = box->Length.v[dd];
    delta -= L * round(delta / L);
    d2 += delta * delta;
  }
  return sqrt(d2);
}

static void record(int i, int j, const SYSTEM System, void *ud) {
  struct rec *r = (struct rec *)ud;
  int a = System.BeadCoor[i], b = System.BeadCoor[j];
  double dist;
  if (r->norm_axis < 0) {
    dist = VectLength(DistancePBC(System.Bead[a].Position,
                                  System.Bead[b].Position, r->box));
  } else {
    dist = inplane_dist(System.Bead[a].Position, System.Bead[b].Position,
                        r->box, r->norm_axis);
  }
  if (dist < r->cutoff) {
    int lo = i < j ? i : j, hi = i < j ? j : i;
    r->count[lo * r->N + hi]++;
  }
}

static bool always(int i, const SYSTEM System, void *ud) {
  (void)i; (void)System; (void)ud;
  return true;
} //}}}

// build an in-memory system with N beads placed inside the cell //{{{
static SYSTEM make_system(int N, BOX box) {
  SYSTEM s = (SYSTEM){0};
  s.Box = box;
  s.Count = InitCount;
  s.Count.Bead = N;
  s.Count.BeadCoor = N;
  s.Bead = malloc(N * sizeof *s.Bead);
  s.BeadCoor = malloc(N * sizeof *s.BeadCoor);
  CHECK(s.Bead && s.BeadCoor);
  for (int i = 0; i < N; i++) {
    InitBead(&s.Bead[i]);
    // r = transform * s_frac places the bead in the primary cell
    double sf[3] = { rand01(), rand01(), rand01() };
    for (int dd = 0; dd < 3; dd++) {
      s.Bead[i].Position.v[dd] = box.transform[dd][0] * sf[0] +
                                 box.transform[dd][1] * sf[1] +
                                 box.transform[dd][2] * sf[2];
    }
    s.Bead[i].InTimestep = true;
    s.BeadCoor[i] = i;
  }
  return s;
}
static void free_system(SYSTEM *s) {
  free(s->Bead);
  free(s->BeadCoor);
}
static BOX box_mode0(vec3d L, double al, double be, double ga) {
  BOX b = InitBox;
  b.Length = L;
  b.alpha = al;
  b.beta = be;
  b.gamma = ga;
  CHECK(CalculateBoxData(&b, 0));
  return b;
} //}}}

// core differential check: brute vs linked list produce identical pair sets //{{{
static void diff_check(BOX box, int N, double cutoff, int norm_axis) {
  SYSTEM s = make_system(N, box);
  long *brute = calloc((size_t)N * N, sizeof *brute);
  long *cells = calloc((size_t)N * N, sizeof *cells);
  CHECK(brute && cells);

  struct rec rb = { &box, cutoff, N, norm_axis, brute };
  struct rec rc = { &box, cutoff, N, norm_axis, cells };

  // ground truth: brute force enumerates every pair, callback applies cutoff
  rb.count = brute;
  TraverseBrutePairs(s, record, &rb, always, nullptr);
  // acceleration structure under test
  if (norm_axis < 0) {
    TraverseLinkedListPairs(s, cutoff, record, &rc, always, nullptr);
  } else {
    TraverseLinkedListPairs2D(s, cutoff, norm_axis, record, &rc, always, nullptr);
  }

  long n_pairs = 0, n_dup = 0, n_mismatch = 0;
  for (long k = 0; k < (long)N * N; k++) {
    if (brute[k] > 1 || cells[k] > 1) {
      n_dup++; // a pair was reported more than once by some traversal
    }
    if (brute[k] != cells[k]) {
      n_mismatch++;
    }
    n_pairs += brute[k];
  }
  CHECK(n_dup == 0);
  CHECK(n_mismatch == 0);
  CHECK(n_pairs > 0); // sanity: the config actually has contacts to compare

  free(brute);
  free(cells);
  free_system(&s);
} //}}}

static void test_diff_orthogonal_3d(void) {
  diff_check(box_mode0((vec3d){ .v = {24, 26, 22} }, 90, 90, 90), 400, 3.0, -1);
}
static void test_diff_triclinic_3d(void) {
  // mild and strong tilts; both must match brute exactly
  diff_check(box_mode0((vec3d){ .v = {24, 26, 22} }, 85, 92, 97), 400, 3.0, -1);
  diff_check(box_mode0((vec3d){ .v = {24, 26, 22} }, 70, 80, 115), 400, 2.5, -1);
}
static void test_diff_2d_slit(void) {
  // orthogonal slit: bin only in the plane perpendicular to each axis
  for (int norm = 0; norm < 3; norm++) {
    diff_check(box_mode0((vec3d){ .v = {24, 26, 22} }, 90, 90, 90), 400, 3.0, norm);
  }
}

int main(void) {
  pcg32Seed(&rng, 424242u);
  RUN(test_diff_orthogonal_3d);
  RUN(test_diff_triclinic_3d);
  RUN(test_diff_2d_slit);
  return test_main_end();
}
