/*
 * General utility and option-parsing tests.
 *
 * These are the small functions every utility funnels its command line and
 * every reader funnels its input through, so a change in their behaviour is
 * felt everywhere at once and usually as wrong results rather than a crash.
 *
 * Several of the assertions below pin behaviour that is surprising rather
 * than obviously correct - strtol's base-0 parsing making "010" octal, a
 * partially numeric string being accepted, the two option parsers disagreeing
 * about negative numbers. Those are marked where they appear. The point is to
 * make a future change to any of them visible, not to bless them.
 *
 * Option parsers exit() on malformed input, so those cases run in a forked
 * child and are classified by how the child terminated.
 */
#include "../src/AnalysisTools.h"
#include "test_util.h"
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/wait.h>
#include <limits.h>

// ---- number parsers -------------------------------------------------------

// IsRealNumber / IsPosRealNumber //{{{
static void test_real_numbers(void) {
  double v;
  // plain values
  CHECK(IsRealNumber("2.5", &v) && v == 2.5);
  CHECK(IsRealNumber("-3.5", &v) && v == -3.5);
  CHECK(IsRealNumber("0", &v) && v == 0);
  CHECK(IsRealNumber("1e3", &v) && v == 1000);
  CHECK(IsRealNumber(".5", &v) && v == 0.5);
  // leading whitespace is skipped by strtod
  CHECK(IsRealNumber("  5", &v) && v == 5);
  /*
   * Documented in General.c: conversion stops at the first illegal character,
   * so a string that merely starts with a number is accepted. This is what
   * makes "5x" a valid option argument.
   */
  CHECK(IsRealNumber("02.2x", &v) && v == 2.2);
  CHECK(IsRealNumber("5 6", &v) && v == 5);
  // nothing numeric at the front is a failure
  CHECK(!IsRealNumber("x2.2", &v));
  CHECK(!IsRealNumber("", &v));
  CHECK(!IsRealNumber("abc", &v));
  // strtod also accepts these spellings, so they reach the utilities
  CHECK(IsRealNumber("inf", &v) && isinf(v));
  CHECK(IsRealNumber("nan", &v) && isnan(v));

  // positive-only variant
  CHECK(IsPosRealNumber("0.001", &v) && v == 0.001);
  CHECK(!IsPosRealNumber("0", &v));    // zero is not positive
  CHECK(!IsPosRealNumber("-1", &v));
  CHECK(!IsPosRealNumber("abc", &v));
} //}}}
// IsIntegerNumber / IsNaturalNumber / IsWholeNumber //{{{
static void test_integer_numbers(void) {
  long v;
  CHECK(IsIntegerNumber("42", &v) && v == 42);
  CHECK(IsIntegerNumber("-7", &v) && v == -7);
  CHECK(IsIntegerNumber("0", &v) && v == 0);
  CHECK(IsIntegerNumber("12abc", &v) && v == 12); // stops at 'a'
  CHECK(!IsIntegerNumber("abc", &v));
  CHECK(!IsIntegerNumber("", &v));
  /*
   * strtol is called with base 0, so C integer-literal prefixes apply: a
   * leading 0 means octal and 0x means hex. "-sk 010" is therefore 8, not 10.
   */
  CHECK(IsIntegerNumber("0x10", &v) && v == 16);
  CHECK(IsIntegerNumber("010", &v) && v == 8);
  CHECK(IsIntegerNumber("08", &v) && v == 0); // invalid octal digit: stops at 8

  // natural: strictly positive
  CHECK(IsNaturalNumber("1", &v) && v == 1);
  CHECK(!IsNaturalNumber("0", &v));
  CHECK(!IsNaturalNumber("-1", &v));
  // whole: non-negative
  CHECK(IsWholeNumber("0", &v) && v == 0);
  CHECK(IsWholeNumber("3", &v) && v == 3);
  CHECK(!IsWholeNumber("-1", &v));
} //}}}

// ---- line splitting -------------------------------------------------------

// SplitLine: counts, delimiters, and the max_str cap //{{{
static void test_split_line(void) {
  char *out[8];
  char buf[LINE];

  // ordinary case
  s_strcpy(buf, "one two three", LINE);
  CHECK(SplitLine(8, out, buf, " ") == 3);
  CHECK(strcmp(out[0], "one") == 0);
  CHECK(strcmp(out[1], "two") == 0);
  CHECK(strcmp(out[2], "three") == 0);

  // runs of delimiters, and mixed delimiter characters, collapse
  s_strcpy(buf, "  a\t\tb \t c  ", LINE);
  CHECK(SplitLine(8, out, buf, " \t") == 3);
  CHECK(strcmp(out[0], "a") == 0);
  CHECK(strcmp(out[2], "c") == 0);

  // trailing newline is a delimiter in the readers' usual call
  s_strcpy(buf, "a b\n", LINE);
  CHECK(SplitLine(8, out, buf, " \t\n") == 2);
  CHECK(strcmp(out[1], "b") == 0);

  // degenerate inputs
  s_strcpy(buf, "", LINE);
  CHECK(SplitLine(8, out, buf, " ") == 0);
  s_strcpy(buf, "     ", LINE);
  CHECK(SplitLine(8, out, buf, " ") == 0);
  s_strcpy(buf, "single", LINE);
  CHECK(SplitLine(8, out, buf, " ") == 1);
} //}}}
// SplitLine must never write past out[max_str-1] //{{{
/*
 * The cap is what keeps a long input line from running off the end of the
 * caller's array, so check it with sentinels rather than trusting the count.
 */
static void test_split_line_cap(void) {
  const int max_str = 4;
  char *out[6];
  char *canary = (char *)0x1;
  for (int i = 0; i < 6; i++) {
    out[i] = canary;
  }
  char buf[LINE];
  s_strcpy(buf, "a b c d e f g", LINE); // 7 words into a 4-slot array
  int n = SplitLine(max_str, out, buf, " ");
  CHECK(n == max_str - 1); // one slot is reserved for the terminating lookahead
  CHECK(strcmp(out[0], "a") == 0);
  CHECK(strcmp(out[max_str - 1], "d") == 0);
  // everything past the array the caller declared must be untouched
  CHECK(out[max_str] == canary);
  CHECK(out[max_str + 1] == canary);

  // exactly max_str-1 words fits without truncation
  for (int i = 0; i < 6; i++) {
    out[i] = canary;
  }
  s_strcpy(buf, "a b c", LINE);
  CHECK(SplitLine(max_str, out, buf, " ") == 3);
  CHECK(out[3] == nullptr); // the lookahead that ended the loop
  CHECK(out[4] == canary);
} //}}}

// ---- line reading ---------------------------------------------------------

// write content to a temp file and return the path (static buffer) //{{{
static const char *temp_with(const char *content) {
  static char path[64];
  strcpy(path, "/tmp/general_test_XXXXXX");
  int fd = mkstemp(path);
  CHECK(fd >= 0);
  size_t len = strlen(content);
  CHECK(write(fd, content, len) == (ssize_t)len);
  close(fd);
  return path;
} //}}}
// ReadLine: normal lines, EOF, and a missing final newline //{{{
static void test_read_line(void) {
  char buf[LINE];
  const char *path = temp_with("first\nsecond\nthird");
  FILE *fr = fopen(path, "r");
  CHECK(fr != nullptr);
  CHECK(ReadLine(fr, buf));
  CHECK(strcmp(buf, "first\n") == 0);
  CHECK(ReadLine(fr, buf));
  CHECK(strcmp(buf, "second\n") == 0);
  // a final line without a newline is still returned
  CHECK(ReadLine(fr, buf));
  CHECK(strcmp(buf, "third") == 0);
  CHECK(!ReadLine(fr, buf)); // EOF
  fclose(fr);
  unlink(path);

  // an empty file yields nothing at all
  const char *empty = temp_with("");
  FILE *fe = fopen(empty, "r");
  CHECK(!ReadLine(fe, buf));
  fclose(fe);
  unlink(empty);
} //}}}
// a line longer than LINE is truncated, and its tail is discarded //{{{
/*
 * Without the drain loop the remainder of an over-long line would come back
 * as the next line, silently shifting everything a reader sees after it.
 */
static void test_read_line_overlong(void) {
  char big[LINE * 2];
  for (size_t i = 0; i < sizeof big - 1; i++) {
    big[i] = 'x';
  }
  big[sizeof big - 1] = '\0';
  char content[LINE * 3];
  snprintf(content, sizeof content, "%s\nnext line\n", big);

  const char *path = temp_with(content);
  FILE *fr = fopen(path, "r");
  CHECK(fr != nullptr);
  char buf[LINE];
  CHECK(ReadLine(fr, buf));
  CHECK(strlen(buf) == LINE - 1); // filled to capacity
  // the next read must be the following line, not the tail of the long one
  CHECK(ReadLine(fr, buf));
  CHECK(strcmp(buf, "next line\n") == 0);
  CHECK(!ReadLine(fr, buf));
  fclose(fr);
  unlink(path);
} //}}}
// ReadAndSplitLine drives the globals the readers use //{{{
static void test_read_and_split_line(void) {
  const char *path = temp_with("alpha beta\n\ngamma\n");
  FILE *fr = fopen(path, "r");
  CHECK(fr != nullptr);
  CHECK(ReadAndSplitLine(fr, SPL_STR, " \t\n"));
  CHECK(words == 2);
  CHECK(strcmp(split[0], "alpha") == 0);
  CHECK(strcmp(split[1], "beta") == 0);
  // a blank line reads successfully but splits into nothing
  CHECK(ReadAndSplitLine(fr, SPL_STR, " \t\n"));
  CHECK(words == 0);
  CHECK(ReadAndSplitLine(fr, SPL_STR, " \t\n"));
  CHECK(words == 1);
  CHECK(strcmp(split[0], "gamma") == 0);
  CHECK(!ReadAndSplitLine(fr, SPL_STR, " \t\n")); // EOF
  fclose(fr);
  unlink(path);
} //}}}
// SkipLine consumes exactly one line //{{{
static void test_skip_line(void) {
  const char *path = temp_with("one\ntwo\nthree\n");
  FILE *fr = fopen(path, "r");
  CHECK(fr != nullptr);
  SkipLine(fr);
  char buf[LINE];
  CHECK(ReadLine(fr, buf));
  CHECK(strcmp(buf, "two\n") == 0);
  SkipLine(fr);
  CHECK(!ReadLine(fr, buf)); // both remaining lines consumed
  fclose(fr);
  unlink(path);
} //}}}

// ---- string and array helpers ---------------------------------------------

// s_strcpy always terminates and never overruns //{{{
static void test_s_strcpy(void) {
  char dst[8];
  // exact fit
  memset(dst, 'Z', sizeof dst);
  s_strcpy(dst, "abcdefg", 8); // 7 chars + NUL
  CHECK(strcmp(dst, "abcdefg") == 0);
  // truncation keeps dest_size-1 characters and a NUL
  memset(dst, 'Z', sizeof dst);
  s_strcpy(dst, "abcdefghij", 8);
  CHECK(strlen(dst) == 7);
  CHECK(strcmp(dst, "abcdefg") == 0);
  // a one-byte destination becomes the empty string
  char one[1] = {'Z'};
  s_strcpy(one, "abc", 1);
  CHECK(one[0] == '\0');
  // empty source
  s_strcpy(dst, "", 8);
  CHECK(dst[0] == '\0');
} //}}}
// StripPath returns the basename //{{{
static void test_strip_path(void) {
  CHECK(strcmp(StripPath("/usr/bin/Info"), "Info") == 0);
  CHECK(strcmp(StripPath("Info"), "Info") == 0);
  CHECK(strcmp(StripPath("./Info"), "Info") == 0);
  CHECK(strcmp(StripPath("a/b/c/d"), "d") == 0);
  CHECK(strcmp(StripPath(""), "") == 0);
  // a trailing slash leaves nothing after it
  CHECK(strcmp(StripPath("dir/"), "") == 0);
} //}}}
// the array initialisers, including the row-major 2D ones //{{{
static void test_init_arrays(void) {
  int ai[5];
  InitIntArray(ai, 5, -3);
  for (int i = 0; i < 5; i++) {
    CHECK(ai[i] == -3);
  }
  double ad[4];
  InitDoubleArray(ad, 4, 1.5);
  for (int i = 0; i < 4; i++) {
    CHECK(ad[i] == 1.5);
  }
  bool ab[3];
  InitBoolArray(ab, 3, true);
  for (int i = 0; i < 3; i++) {
    CHECK(ab[i] == true);
  }
  // n == 0 must touch nothing
  int guard[2] = {7, 7};
  InitIntArray(guard, 0, 99);
  CHECK(guard[0] == 7 && guard[1] == 7);

  // 2D initialisers index row-major into a flat buffer
  int flat[3 * 4];
  for (int i = 0; i < 12; i++) {
    flat[i] = 0;
  }
  InitInt2DArray(flat, 3, 4, 5);
  for (int i = 0; i < 12; i++) {
    CHECK(flat[i] == 5);
  }
  double dflat[2 * 3];
  InitDouble2DArray(dflat, 2, 3, 0.25);
  for (int i = 0; i < 6; i++) {
    CHECK(dflat[i] == 0.25);
  }
} //}}}
// SameArrayInt //{{{
static void test_same_array_int(void) {
  int a[4] = {1, 2, 3, 4};
  int b[4] = {1, 2, 3, 4};
  int c[4] = {1, 2, 9, 4};
  CHECK(SameArrayInt(a, b, 4));
  CHECK(!SameArrayInt(a, c, 4));
  // a shorter comparison stops before the difference
  CHECK(SameArrayInt(a, c, 2));
  // an empty comparison is vacuously true
  CHECK(SameArrayInt(a, c, 0));
} //}}}

// ---- timestep selection ---------------------------------------------------

// UseStep: the -st / -e / -sk boundary matrix //{{{
/*
 * These three options decide which frames a utility actually processes, so an
 * off-by-one changes every reported average without any visible error.
 */
static void test_use_step(void) {
  COMMON_OPT o = {0};
  // defaults after CommonOptions with no options given: start 1, end -1, skip 1
  o.start = 1;
  o.end = -1;
  o.skip = 1;
  CHECK(!UseStep(o, 0));  // before the start
  CHECK(UseStep(o, 1));
  CHECK(UseStep(o, 2));
  CHECK(UseStep(o, 1000)); // end == -1 means unbounded

  // -st 3
  o.start = 3;
  CHECK(!UseStep(o, 2));
  CHECK(UseStep(o, 3)); // the start step itself is used
  CHECK(UseStep(o, 4));

  // -e 5 with -st 3: the end step is inclusive
  o.end = 5;
  CHECK(UseStep(o, 5));
  CHECK(!UseStep(o, 6));

  // -sk 1 becomes skip == 2: every other step from the start
  o.start = 1;
  o.end = -1;
  o.skip = 2;
  CHECK(UseStep(o, 1));
  CHECK(!UseStep(o, 2));
  CHECK(UseStep(o, 3));
  CHECK(!UseStep(o, 4));
  CHECK(UseStep(o, 5));

  // the skip pattern is anchored on start, not on step 1
  o.start = 4;
  o.skip = 3;
  CHECK(UseStep(o, 4));
  CHECK(!UseStep(o, 5));
  CHECK(!UseStep(o, 6));
  CHECK(UseStep(o, 7));
} //}}}

// ---- option parsing -------------------------------------------------------

// BoolOption //{{{
static void test_bool_option(void) {
  char *argv[] = {(char *)"prog", (char *)"file", (char *)"--verbose",
                  (char *)"-x", nullptr};
  int argc = 4;
  CHECK(BoolOption(argc, argv, "--verbose"));
  CHECK(BoolOption(argc, argv, "-x"));
  CHECK(!BoolOption(argc, argv, "--silent"));
  // argv[0] is never considered an option
  CHECK(!BoolOption(argc, argv, "prog"));
} //}}}
// NumbersOption and its one/two/three wrappers //{{{
static void test_numbers_option(void) {
  {
    char *argv[] = {(char *)"prog", (char *)"-n", (char *)"5", nullptr};
    int val = 0;
    CHECK(OneNumberOption(3, argv, "-n", &val, 'i'));
    CHECK(val == 5);
    // an absent option leaves the value alone and returns false
    int untouched = 99;
    CHECK(!OneNumberOption(3, argv, "-q", &untouched, 'i'));
    CHECK(untouched == 99);
  }
  {
    char *argv[] = {(char *)"prog", (char *)"-n", (char *)"2", (char *)"8",
                    nullptr};
    int val[2] = {0, 0};
    CHECK(TwoNumbersOption(4, argv, "-n", val, 'i'));
    CHECK(val[0] == 2 && val[1] == 8);
  }
  {
    char *argv[] = {(char *)"prog", (char *)"-off", (char *)"1.5",
                    (char *)"-2.5", (char *)"0.25", nullptr};
    double val[3] = {0, 0, 0};
    CHECK(ThreeNumbersOption(5, argv, "-off", val, 'd'));
    CHECK(val[0] == 1.5);
    // NumbersOption parses numbers with IsIntegerNumber/IsRealNumber, so a
    // negative value is accepted even though it starts with '-'
    CHECK(val[1] == -2.5);
    CHECK(val[2] == 0.25);
  }
  {
    // parsing stops at the first non-numeric argument
    char *argv[] = {(char *)"prog", (char *)"-n", (char *)"3",
                    (char *)"--verbose", nullptr};
    int count = 0;
    int val[4] = {0};
    CHECK(NumbersOption(4, argv, 4, "-n", &count, val, 'i'));
    CHECK(count == 1);
    CHECK(val[0] == 3);
  }
  {
    // more arguments than max: the count is clamped and the first max kept
    char *argv[] = {(char *)"prog", (char *)"-n", (char *)"1", (char *)"2",
                    (char *)"3", nullptr};
    int count = 0;
    int val[2] = {0, 0};
    CHECK(NumbersOption(5, argv, 2, "-n", &count, val, 'i'));
    CHECK(count == 2);
    CHECK(val[0] == 1 && val[1] == 2);
  }
} //}}}
// FileOption / FileNumbersOption //{{{
static void test_file_options(void) {
  {
    char *argv[] = {(char *)"prog", (char *)"-o", (char *)"out.txt", nullptr};
    char file[LINE] = "";
    CHECK(FileOption(3, argv, "-o", file));
    CHECK(strcmp(file, "out.txt") == 0);
    char none[LINE] = "";
    CHECK(!FileOption(3, argv, "-q", none));
    CHECK(none[0] == '\0'); // cleared even when the option is absent
  }
  {
    char *argv[] = {(char *)"prog", (char *)"-d", (char *)"data.txt",
                    (char *)"3", (char *)"7", nullptr};
    char file[LINE] = "";
    int val[2] = {0, 0};
    int count = 0;
    CHECK(FileNumbersOption(5, argv, 0, 2, "-d", val, &count, file, 'i'));
    CHECK(strcmp(file, "data.txt") == 0);
    CHECK(count == 2);
    CHECK(val[0] == 3 && val[1] == 7);
  }
  {
    /*
     * Unlike NumbersOption, this one stops at any argument beginning with '-',
     * so a negative number cannot be passed here. Pinned because the two
     * parsers disagreeing is a trap for anyone adding an option.
     */
    char *argv[] = {(char *)"prog", (char *)"-d", (char *)"data.txt",
                    (char *)"-5", nullptr};
    char file[LINE] = "";
    int val[2] = {0, 0};
    int count = 0;
    CHECK(FileNumbersOption(4, argv, 0, 2, "-d", val, &count, file, 'i'));
    CHECK(strcmp(file, "data.txt") == 0);
    CHECK(count == 0); // "-5" was treated as the next option, not a number
  }
} //}}}

// run a callback in a child process and return how it terminated //{{{
struct child_exit { bool exited; int code; bool signaled; };
static struct child_exit run_in_child(void (*fn)(void)) {
  fflush(nullptr);
  pid_t pid = fork();
  if (pid == 0) {
    TestChildNoLeakReports(); // the option parsers exit(); see test_util.h
    int devnull = open("/dev/null", O_WRONLY);
    if (devnull >= 0) {
      dup2(devnull, STDERR_FILENO);
      dup2(devnull, STDOUT_FILENO);
    }
    fn();
    _exit(0); // reached only if the parser accepted the input
  }
  int status = 0;
  waitpid(pid, &status, 0);
  struct child_exit c = {0};
  if (WIFEXITED(status)) {
    c.exited = true;
    c.code = WEXITSTATUS(status);
  }
  if (WIFSIGNALED(status)) {
    c.signaled = true;
  }
  return c;
} //}}}
// malformed option arguments must be rejected, not silently accepted //{{{
static void child_missing_number(void) {
  char *argv[] = {(char *)"prog", (char *)"-n", nullptr};
  int count = 0;
  int val[2] = {0, 0};
  NumbersOption(2, argv, 2, "-n", &count, val, 'i');
}
static void child_wrong_count(void) {
  // -n given one number where two are required
  char *argv[] = {(char *)"prog", (char *)"-n", (char *)"5", nullptr};
  int val[2] = {0, 0};
  TwoNumbersOption(3, argv, "-n", val, 'i');
}
static void child_non_numeric_file_arg(void) {
  char *argv[] = {(char *)"prog", (char *)"-d", (char *)"f.txt",
                  (char *)"abc", nullptr};
  char file[LINE] = "";
  int val[2] = {0, 0};
  int count = 0;
  FileNumbersOption(4, argv, 1, 2, "-d", val, &count, file, 'i');
}
static void child_too_few_numbers(void) {
  // min == 2 but only one number supplied
  char *argv[] = {(char *)"prog", (char *)"-d", (char *)"f.txt",
                  (char *)"1", nullptr};
  char file[LINE] = "";
  int val[2] = {0, 0};
  int count = 0;
  FileNumbersOption(4, argv, 2, 2, "-d", val, &count, file, 'i');
}
static void expect_option_rejected(const char *label, void (*fn)(void)) {
  struct child_exit c = run_in_child(fn);
  CHECK(!c.signaled);            // a clean error, not a crash
  CHECK(c.exited && c.code != 0);
  if (c.signaled || (c.exited && c.code == 0)) {
    fprintf(stderr, "  (note: '%s' was not rejected)\n", label);
  }
}
static void test_option_errors(void) {
  expect_option_rejected("-n with no argument", child_missing_number);
  expect_option_rejected("-n with one of two numbers", child_wrong_count);
  expect_option_rejected("file option, non-numeric", child_non_numeric_file_arg);
  expect_option_rejected("file option, too few numbers", child_too_few_numbers);
} //}}}

int main(void) {
  RUN(test_real_numbers);
  RUN(test_integer_numbers);
  RUN(test_split_line);
  RUN(test_split_line_cap);
  RUN(test_read_line);
  RUN(test_read_line_overlong);
  RUN(test_read_and_split_line);
  RUN(test_skip_line);
  RUN(test_s_strcpy);
  RUN(test_strip_path);
  RUN(test_init_arrays);
  RUN(test_same_array_int);
  RUN(test_use_step);
  RUN(test_bool_option);
  RUN(test_numbers_option);
  RUN(test_file_options);
  RUN(test_option_errors);
  return test_main_end();
}
