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

#define PI 3.14159265358979323846264338327950
#define LINE 1024 // maximum length of an array for strings
#define SPL_STR 32 // maximum number of split strings
#define SPL_LEN 64 // maximum length of split strings
#define SQR(x) ((x)*(x)) // macro for algebraic square
#define CUBE(x) ((x)*(x)*(x)) // macro for algebraic cube
#define VECTORLENGTH(x) (sqrt(SQR(x[0])+SQR(x[1])+SQR(x[2])))
#define ID2D(x,y,size) ((x)*(size)[0]+(y))

#define BLACK   "\033[1;30m"
#define RED     "\033[1;31m"
#define GREEN   "\033[1;32m"
#define YELLOW  "\033[1;33m"
#define BLUE    "\033[1;34m"
#define MAGENTA "\033[1;35m"
#define CYAN    "\033[1;36m"
#define WHITE   "\033[1;37m"
#define C_RESET "\033[0m"

extern char line[LINE], *split[SPL_STR];
extern int words;

// deffine vector structures //{{{
typedef struct Vector {
  double x, y, z;
} VECTOR;
typedef struct LongVector {
  long double x, y, z;
} LONGVECTOR;
typedef struct IntVector {
  int x, y, z;
} INTVECTOR;
typedef struct LongIntVector {
  long int x, y, z;
} LONGINTVECTOR;
//}}}
// convert string into number if possible
bool IsRealNumber(char *str, double *val);
bool IsPosRealNumber(char *str, double *val);
bool IsIntegerNumber(char *str, long *val);
bool IsNaturalNumber(char *str, long *val);
bool IsWholeNumber(char *str, long *val);
// minimum/maximum from three numbers
double Min3(double x, double y, double z);
double Max3(double x, double y, double z);
// swap values
void SwapInt(int *a, int *b);
void SwapDouble(double *a, double *b);
void SwapBool(bool *a, bool *b);

void SortArrayInt(int *array, int length, int mode);
void SortArrayDouble(double *array, int length, int mode);

int SplitLine(int max_str, char *out[], char *line, const char *delim);
bool ReadAndSplitLine(FILE *fr, int max_strings, const char *delim);

char * BareCommand(char cmd[]);
void PrintCommand(FILE *ptr, int argc, char *argv[]);

// changing the text colour (and making it bold) for cli output //{{{
char *Colour(FILE *f, char *colour);
// colours for stderr
char *ErrRed();
char *ErrCyan();
char *ErrYellow();
char *ErrColourReset();
// colours for stdout
char *Red();
char *Cyan();
char *Yellow();
char *Magenta();
char *Green();
char *White();
char *ColourReset();
void ColourChange(int a, char *colour);
 //}}}

FILE *OpenFile(char *file, char *mode);

// initialize arrays to specified value
void InitDoubleArray (double array[], int n, double val);
void InitIntArray (int array[], int n, int val);
void InitBoolArray (bool array[], int n, bool val);
void InitVecArray (VECTOR array[], int n, bool val);
void InitLong2DArray (long *array[], int m, int n, long val);
void InitDouble2DArray (double *array[], int m, int n, double val);

// test whether two arrays are the same
bool SameArrayInt(const int arr_1[], const int arr_2[], int n);

// safe function alternatives
void s_strcpy(char *dest, const char *src, size_t dest_size);
void s_snprintf(char *buffer, size_t size, const char *format, ...);
void* s_realloc(void *ptr, size_t new_size);

#endif
