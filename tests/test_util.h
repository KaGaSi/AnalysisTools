#ifndef TEST_UTIL_H
#define TEST_UTIL_H

/*
 * Minimal assert-based test harness (no external framework).
 *
 * Each test executable defines test functions and runs them from main() via
 * RUN(). Checks accumulate into globals; TEST_MAIN_END() returns non-zero if
 * anything failed, so CTest reports the executable as failed.
 */

#include <stdio.h>
#include <math.h>

static int g_failures = 0;
static int g_checks = 0;

// boolean check
#define CHECK(cond)                                                           \
  do {                                                                        \
    g_checks++;                                                               \
    if (!(cond)) {                                                            \
      g_failures++;                                                           \
      fprintf(stderr, "  FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond);       \
    }                                                                         \
  } while (0)

// floating-point closeness check
#define CHECK_CLOSE(a, b, tol)                                                \
  do {                                                                        \
    g_checks++;                                                               \
    double _a = (a), _b = (b), _t = (tol);                                    \
    if (!(fabs(_a - _b) <= _t)) {                                             \
      g_failures++;                                                           \
      fprintf(stderr, "  FAIL %s:%d: |%.12g - %.12g| = %.3g > %.3g  (%s ~ %s)\n", \
              __FILE__, __LINE__, _a, _b, fabs(_a - _b), _t, #a, #b);         \
    }                                                                         \
  } while (0)

// run a single test function, reporting per-test pass/fail
#define RUN(test)                                                             \
  do {                                                                        \
    int _before = g_failures;                                                 \
    test();                                                                   \
    if (g_failures == _before) {                                              \
      fprintf(stderr, "[  OK  ] %s\n", #test);                                \
    } else {                                                                  \
      fprintf(stderr, "[ FAIL ] %s\n", #test);                               \
    }                                                                         \
  } while (0)

// print summary and yield process exit status
static inline int test_main_end(void) {
  fprintf(stderr, "\n%d checks, %d failure(s)\n", g_checks, g_failures);
  if (g_failures) {
    return 1;
  } else {
    return 0;
  }
}

#endif
