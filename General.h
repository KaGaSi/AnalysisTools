#ifndef GENERAL_H
#define GENERAL_H

#define _POSIX_C_SOURCE 200809L

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <unistd.h>
#include <stdarg.h>
#include <stdint.h>
#include "Globals.h"
#include "MathUtils.h"

// swap values //{{{
static inline void SwapInt(int *a, int *b) {
  int swap = *a;
  *a = *b;
  *b = swap;
}
static inline void SwapDouble(double *a, double *b) {
  double swap = *a;
  *a = *b;
  *b = swap;
}
static inline void SwapBool(bool *a, bool *b) {
  bool swap = *a;
  *a = *b;
  *b = swap;
} //}}}
// minimum/maximum from three numbers //{{{
static inline double Min3(double x, double y, double z) {
  double min = y;
  if (x < y) {
    min = x;
  }
  if (min < z) {
    return min;
  } else {
    return z;
  }
}
static inline double Max3(double x, double y, double z) {
  double max = y;
  if (x > y) {
    max = x;
  }
  if (max > z) {
    return max;
  } else {
    return z;
  }
} //}}}
// changing the text colour (and making it bold) for cli output //{{{
static inline const char *Colour(FILE *f, const char colour[]) {
  if (isatty(fileno(f))) {
    return colour;
  } else {
    return "";
  }
}
// colours for stderr
static inline const char *ErrRed() {
  return Colour(stderr, RED);
}
static inline const char *ErrCyan() {
  return Colour(stderr, CYAN);
}
static inline const char *ErrYellow() {
  return Colour(stderr, YELLOW);
}
static inline const char *ErrColourReset() {
  return Colour(stderr, C_RESET);
}
// colours for stdout
static inline const char *Red() {
  return Colour(stdout, RED);
}
static inline const char *Cyan() {
  return Colour(stdout, CYAN);
}
static inline const char *Yellow() {
  return Colour(stdout, YELLOW);
}
static inline const char *Magenta() {
  return Colour(stdout, MAGENTA);
}
static inline const char *Green() {
  return Colour(stdout, GREEN);
}
static inline const char *White() {
  return Colour(stdout, WHITE);
}
static inline const char *ColourReset() {
  return Colour(stdout, C_RESET);
}
void ColourChange(const int a, const char *colour);
 //}}}
// convert string into number if possible
bool IsRealNumber(const char *str, double *val);
bool IsPosRealNumber(const char *str, double *val);
bool IsIntegerNumber(const char *str, long *val);
bool IsNaturalNumber(const char *str, long *val);
bool IsWholeNumber(const char *str, long *val);
// bubble sort int/double array ascendingly/descendingly
void SortArray(void *array, const int length, const int mode, const char type);
// line reading and splitting
bool ReadLine(FILE *fr, char *line);
int SplitLine(const int max_str, char **out, char *line, const char *delim);
bool ReadAndSplitLine(FILE *fr, const int max_strings, const char *delim);
void WriteSplitLine(FILE *f);
// strip path from command
const char * StripPath(const char cmd[]);
void PrintCommand(FILE *ptr, const int argc, char **argv);
// wrapper around fopen
FILE * OpenFile(const char *file, char *mode);
// initialize arrays to specified value
void InitDoubleArray (double *array, const int n, const double val);
void InitIntArray (int *array, const int n, const int val);
void InitBoolArray (bool *array, const int n, const bool val);
void InitLong2DArray (long **arr, const int m, const int n, const long val);
void InitDouble2DArray (double **arr, const int m, const int n,
                        const double val);
void InitInt2DArray (int *arr, const int m, const int n, const int val);
void InitBool2DArray (bool **arr, const int m, const int n, const bool val);
// test whether two arrays are the same
bool SameArrayInt(const int arr_1[], const int arr_2[], const int n);
// safe function alternatives
void s_strcpy(char *dest, const char *src, size_t dest_size);
void * s_realloc(void *ptr, size_t new_size);

#endif
