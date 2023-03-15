#ifndef _GENERAL_H_
#define _GENERAL_H_

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <unistd.h>

#define PI 3.141593 // value of pi
#define LINE 1024 // maximum length of an array for strings
#define SPL_STR 32 // maximum number of split strings
#define SPL_LEN 64 // maximum length of split strings
#define SQR(x) ((x)*(x)) // macro for algebraic square
#define CUBE(x) ((x)*(x)*(x)) // macro for algebraic cube

#define BLACK   "\033[1;30m"
#define RED     "\033[1;31m"
#define GREEN   "\033[1;32m"
#define YELLOW  "\033[1;33m"
#define BLUE    "\033[1;34m"
#define MAGENTA "\033[1;35m"
#define CYAN    "\033[1;36m"
#define WHITE   "\033[1;37m"
#define C_RESET "\033[0m"

// needs to be here so I can use fileno() from a standard library
int fileno(const FILE *stream);

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
//}}}

double VectorLength(VECTOR a);

// convert string into number if possible //{{{
bool IsRealNumber(char *str, double *val);
bool IsPosRealNumber(char *str, double *val);
bool IsIntegerNumber(char *str, long *val);
bool IsNaturalNumber(char *str, long *val);
bool IsWholeNumber(char *str, long *val);
 //}}}

double Min3(double x, double y, double z);
double Max3(double x, double y, double z);

// swapping functions //{{{
void SwapInt(int *a, int *b);
void SwapDouble(double *a, double *b);
void SwapBool(bool *a, bool *b);
// }}}

void SortArray(int *array, int length, int mode);
VECTOR SortVector(VECTOR in);

bool ReadLine(FILE *fr, int max_char, char *line);
int SplitLine(int max_str, char *out[], char *line, const char *delim);
bool ReadAndSplitLine(FILE *fr, int max_char, char *line, int *words,
                      char *out[], int max_strings, const char *delim);

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
void InitIntArray (int array[], int n, int val);
void InitBoolArray (int array[], int n, bool val);
void InitLong2DArray (long *array[], int m, int n, long val);
void InitDouble2DArray (double *array[], int m, int n, double val);

// test whether two arrays are the same
bool SameArray(int arr_1[], int arr_2[], int n);

#endif
