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

// detect an AddressSanitizer build (GCC and clang spell it differently) //{{{
#if defined(__SANITIZE_ADDRESS__)
#define TEST_HAVE_ASAN 1
#elif defined(__has_feature)
#if __has_feature(address_sanitizer)
#define TEST_HAVE_ASAN 1
#endif
#endif
#ifdef TEST_HAVE_ASAN
#include <sanitizer/lsan_interface.h>
#endif //}}}

// silence leak reports in a forked child that is expected to exit() //{{{
/*
 * The library exits from deep inside on bad input without unwinding, so a
 * child that exercises an error path leaves allocations behind. That is
 * harmless in a process about to die, but LeakSanitizer reports it at exit
 * and the ASAN_OPTIONS=abort_on_error=1 that CTest sets turns the report into
 * a SIGABRT - which a test classifying "died by signal" as a crash would count
 * as a failure. Call this first thing in such a child: it suppresses only the
 * leak report, leaving ASan's memory-error checking (the thing the malformed
 * inputs are actually fuzzing for) fully active, and leaves the parent and
 * every other test executable with leak detection intact.
 */
static inline void TestChildNoLeakReports(void) {
#ifdef TEST_HAVE_ASAN
  __lsan_disable();
#endif
} //}}}

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
