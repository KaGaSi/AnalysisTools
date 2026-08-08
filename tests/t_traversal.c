/*
 * Differential test for the pair traversal.
 *
 * The cell-linked-list traversal is only an acceleration structure: for any
 * system it must enumerate exactly the same set of within-cutoff pairs as the
 * O(N^2) brute force, with no misses and no double counting. We check that
 * equivalence over random orthogonal AND triclinic boxes (3D), plus the 2D
 * slit path. The triclinic case is the regression guard for the fractional
 * cell-binning fix in Pairs.c.
 *
 * The same equivalence is checked through three further paths: with a check
 * callback that filters beads out (the early-outs in both traversals), and
 * through TraversePairs/TraversePairs2D, the dispatchers the utilities
 * actually call, on both sides of the cell-size threshold at which they fall
 * back to brute force.
 */
#include "../src/AnalysisTools.h"
#include "test_util.h"
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/wait.h>

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
  long *count;   // N*N matrix of hit counts, indexed [lo*N+hi], lo < hi
  const bool *use; // if set, a bead the filter rejected must never show up
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
  if (r->use) {
    // a rejected bead must not reach the pair callback at all
    CHECK(r->use[a]);
    CHECK(r->use[b]);
  }
  double dist;
  if (r->norm_axis < 0) {
    dist = VectLength(DistancePBC(System.Bead[a].Position,
                                  System.Bead[b].Position, r->box));
  } else {
    dist = inplane_dist(System.Bead[a].Position, System.Bead[b].Position,
                        r->box, r->norm_axis);
  }
  if (dist < r->cutoff) {
    int lo = j,
        hi = i;
    if (i < j) {
      lo = i;
      hi = j;
    }
    r->count[lo*r->N+hi]++;
  }
}

static bool always(int i, const SYSTEM System, void *ud) {
  (void)i;
  (void)System;
  (void)ud;
  return true;
}

// filter used to exercise the check-callback early-outs //{{{
/*
 * The callback is handed an index into BeadCoor, not a bead id, which is the
 * detail worth pinning: both traversals must agree on that convention, or the
 * two would filter different beads and the pair sets would diverge.
 */
static bool use_flagged(int i, const SYSTEM System, void *ud) {
  const bool *use = (const bool *)ud;
  return use[System.BeadCoor[i]];
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

// which entry point enumerates the pairs //{{{
enum {
  TRAV_DIRECT,  // TraverseLinkedListPairs{,2D} - the concrete implementation
  TRAV_DISPATCH // TraversePairs{,2D} - what the utilities actually call
}; //}}}
// core differential check: brute vs linked list produce identical pair sets //{{{
/*
 * Both traversals get the same check callback, so filtering is part of the
 * equivalence rather than something tested separately.
 */
static void diff_check_full(BOX box, int N, double cutoff, int norm_axis,
                            int mode, check_cb_t check, void *check_ud,
                            const bool *use) {
  SYSTEM s = make_system(N, box);
  long *brute = calloc((size_t)N * N, sizeof *brute);
  long *cells = calloc((size_t)N * N, sizeof *cells);
  CHECK(brute && cells);

  struct rec rb = { &box, cutoff, N, norm_axis, brute, use };
  struct rec rc = { &box, cutoff, N, norm_axis, cells, use };

  // ground truth: brute force enumerates every pair, callback applies cutoff
  TraverseBrutePairs(s, record, &rb, check, check_ud);
  // acceleration structure under test
  if (mode == TRAV_DIRECT) {
    if (norm_axis < 0) {
      TraverseLinkedListPairs(s, cutoff, record, &rc, check, check_ud);
    } else {
      TraverseLinkedListPairs2D(s, cutoff, norm_axis, record, &rc,
                                check, check_ud);
    }
  } else {
    if (norm_axis < 0) {
      TraversePairs(s, cutoff, record, &rc, check, check_ud);
    } else {
      TraversePairs2D(s, cutoff, norm_axis, record, &rc, check, check_ud);
    }
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
}
static void diff_check(BOX box, int N, double cutoff, int norm_axis) {
  diff_check_full(box, N, cutoff, norm_axis, TRAV_DIRECT, always, nullptr,
                  nullptr);
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

// the check callback must filter identically in both traversals //{{{
/*
 * Nothing else exercises the early-outs in TraverseLL and TraverseBrutePairs:
 * every other test passes a callback that accepts everything, leaving those
 * branches dead. A filter that rejects a random subset must leave the two
 * traversals still enumerating exactly the same pairs, and must keep the
 * rejected beads out of the pair callback entirely (asserted in record()).
 */
static void test_check_callback_filter(void) {
  const int N = 400;
  bool *use = malloc(N * sizeof *use);
  CHECK(use != nullptr);
  int kept = 0;
  for (int i = 0; i < N; i++) {
    use[i] = (pcg32Rand0Int(&rng, 100) < 60); // keep roughly 60%
    if (use[i]) {
      kept++;
    }
  }
  CHECK(kept > 0 && kept < N); // the filter must actually filter

  BOX ortho = box_mode0((vec3d){ .v = {24, 26, 22} }, 90, 90, 90);
  BOX tric = box_mode0((vec3d){ .v = {24, 26, 22} }, 85, 92, 97);
  diff_check_full(ortho, N, 3.0, -1, TRAV_DIRECT, use_flagged, use, use);
  diff_check_full(tric, N, 3.0, -1, TRAV_DIRECT, use_flagged, use, use);
  // and through the 2D path, for every slit normal
  for (int norm = 0; norm < 3; norm++) {
    diff_check_full(ortho, N, 3.0, norm, TRAV_DIRECT, use_flagged, use, use);
  }
  free(use);
}
// a filter that rejects everything must produce no pairs at all //{{{
static bool never(int i, const SYSTEM System, void *ud) {
  (void)i;
  (void)System;
  (void)ud;
  return false;
}
static void count_calls(int i, int j, const SYSTEM System, void *ud) {
  (void)i;
  (void)j;
  (void)System;
  (*(long *)ud)++;
}
static void test_check_callback_rejects_all(void) {
  BOX box = box_mode0((vec3d){ .v = {24, 26, 22} }, 90, 90, 90);
  SYSTEM s = make_system(200, box);
  long calls = 0;
  TraverseBrutePairs(s, count_calls, &calls, never, nullptr);
  CHECK(calls == 0);
  TraverseLinkedListPairs(s, 3.0, count_calls, &calls, never, nullptr);
  CHECK(calls == 0);
  TraverseLinkedListPairs2D(s, 3.0, 2, count_calls, &calls, never, nullptr);
  CHECK(calls == 0);
  TraversePairs(s, 3.0, count_calls, &calls, never, nullptr);
  CHECK(calls == 0);
  TraversePairs2D(s, 3.0, 2, count_calls, &calls, never, nullptr);
  CHECK(calls == 0);
  free_system(&s);
} //}}}
 //}}}
// the dispatchers must agree with brute force on both sides of the cut //{{{
/*
 * TraversePairs falls back to brute force once cell_size exceeds a third of
 * the smallest perpendicular cell width; TraversePairs2D applies the same
 * rule to the two in-plane axes only. Both branches have to enumerate the
 * same pairs, so the cut-offs below straddle the threshold deliberately: for
 * a 24x26x22 box the 3D threshold is 22/3 = 7.33, and the 2D one is at worst
 * 22/3 for a normal along x or y and 24/3 = 8 for one along z.
 */
static void test_dispatch_3d(void) {
  BOX ortho = box_mode0((vec3d){ .v = {24, 26, 22} }, 90, 90, 90);
  /*
   * For an orthogonal box the perpendicular widths are just the box lengths,
   * so the switch-over point can be computed here and the two cut-offs below
   * asserted to sit on opposite sides of it. Without this the test could
   * silently drift onto one branch after a change to the heuristic.
   */
  double threshold = Min3(ortho.Length.x, ortho.Length.y, ortho.Length.z) / 3;
  CHECK(3.0 < threshold);
  CHECK(9.0 > threshold);
  // below the threshold: the dispatcher takes the linked-list branch
  diff_check_full(ortho, 400, 3.0, -1, TRAV_DISPATCH, always, nullptr, nullptr);
  // above it: the dispatcher falls back to brute force
  diff_check_full(ortho, 400, 9.0, -1, TRAV_DISPATCH, always, nullptr, nullptr);
  // triclinic, where the widths come from the lattice rather than Length
  BOX tric = box_mode0((vec3d){ .v = {24, 26, 22} }, 85, 92, 97);
  diff_check_full(tric, 400, 3.0, -1, TRAV_DISPATCH, always, nullptr, nullptr);
  diff_check_full(tric, 400, 12.0, -1, TRAV_DISPATCH, always, nullptr, nullptr);
}
static void test_dispatch_2d(void) {
  BOX ortho = box_mode0((vec3d){ .v = {24, 26, 22} }, 90, 90, 90);
  for (int norm = 0; norm < 3; norm++) {
    // the 2D rule looks only at the two in-plane axes
    double in_plane = -1;
    for (int dd = 0; dd < 3; dd++) {
      if (dd != norm && (in_plane < 0 || ortho.Length.v[dd] < in_plane)) {
        in_plane = ortho.Length.v[dd];
      }
    }
    CHECK(3.0 < in_plane / 3);
    CHECK(10.0 > in_plane / 3);
    diff_check_full(ortho, 400, 3.0, norm, TRAV_DISPATCH, always, nullptr,
                    nullptr);
    diff_check_full(ortho, 400, 10.0, norm, TRAV_DISPATCH, always, nullptr,
                    nullptr);
  }
} //}}}
// TraversePairs2D rejects an out-of-range normal axis instead of indexing //{{{
static void child_bad_axis(int axis) {
  TestChildNoLeakReports(); // the validation exits; see test_util.h
  int devnull = open("/dev/null", O_WRONLY);
  if (devnull >= 0) {
    dup2(devnull, STDERR_FILENO);
    dup2(devnull, STDOUT_FILENO);
  }
  BOX box = box_mode0((vec3d){ .v = {24, 26, 22} }, 90, 90, 90);
  SYSTEM s = make_system(50, box);
  long calls = 0;
  TraversePairs2D(s, 3.0, axis, count_calls, &calls, always, nullptr);
  _exit(0); // reached only if the bad axis was accepted
}
static void expect_axis_rejected(int axis) {
  fflush(nullptr);
  pid_t pid = fork();
  if (pid == 0) {
    child_bad_axis(axis);
  }
  int status = 0;
  waitpid(pid, &status, 0);
  CHECK(WIFEXITED(status));            // a clean error, not a crash
  if (WIFEXITED(status)) {
    CHECK(WEXITSTATUS(status) != 0);   // and it must report failure
  }
}
static void test_dispatch_2d_bad_axis(void) {
  expect_axis_rejected(-1);
  expect_axis_rejected(3);
  expect_axis_rejected(100);
} //}}}

int main(void) {
  pcg32Seed(&rng, 424242u);
  RUN(test_diff_orthogonal_3d);
  RUN(test_diff_triclinic_3d);
  RUN(test_diff_2d_slit);
  RUN(test_check_callback_filter);
  RUN(test_check_callback_rejects_all);
  RUN(test_dispatch_3d);
  RUN(test_dispatch_2d);
  RUN(test_dispatch_2d_bad_axis);
  return test_main_end();
}
