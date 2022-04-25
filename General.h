/**
 * \file
 * \brief Functions independent of the analysis utilities.
 */

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

/*
// tell gcc to ignore certain warnings
// helper macro
#define DO_PRAGMA(x) _Pragma(#x)

// macro to ignore warning; x='-W...'
#define P_IGNORE(x) \
  _Pragma("GCC diagnostic push") \
  DO_PRAGMA(GCC diagnostic ignored #x) \

// macro to pop back to the 'pushed' diagnostic state
#define P_POP _Pragma("GCC diagnostic pop")
*/

int fileno(const FILE *stream); // needs to be here for some reason

// struct Vector //{{{
/**
 * \brief 3D vector of floats.
 */
typedef struct Vector {
  double x, y, z;
} VECTOR; //}}}

// struct LongVector //{{{
/**
 * \brief 3D vector of floats.
 */
typedef struct LongVector {
  long double x, y, z;
} LONGVECTOR; //}}}

// struct IntVector //{{{
/**
 * \brief 3D vector of integers.
 */
typedef struct IntVector {
  int x, y, z;
} INTVECTOR; //}}}

// Length() //{{{
/*
 * \brief Function to calculate vector length.
 *
 * \param [in] a    vector
 * \return a's length
 */
double Length(VECTOR a); //}}}

// IsReal() //{{{
/*
 * \brief Function to test if a string is a real number.
 *
 * \param [in] a   string to test
 * \return 'true' if a is double, 'false' otherwise
 */
bool IsReal(char *a); //}}}

// IsReal2() //{{{
/*
 * \brief Function to test if a string is a real number.
 *
 * \param [in] a   string to test
 * \return 'true' if a is double, 'false' otherwise
 */
bool IsReal2(char *str, double *val); //}}}

// IsPosReal() //{{{
/*
 * \brief Function to test if a string is a non-negative real number.
 *
 * \param [in] a   string to test
 * \return 'true' if a is non-negative double, 'false' otherwise
 */
bool IsPosReal(char *a); //}}}

// IsInteger() //{{{
/*
 * \brief Function to test if a string is a non-negative whole number.
 *
 * \param [in] a   string to test
 * \return 'true' if a is integer, 'false' otherwise
 */
bool IsInteger(char *a); //}}}

bool IsInteger2(char *str, long *val);

// IsNatural()  //{{{
bool IsNatural(char *a); //}}}

// Min3() //{{{
/**
 * \brief Function returning the lowest number from three floats.
 *
 * \param [in] x   first double precision number
 * \param [in] y   second double precision number
 * \param [in] z   third double precision number
 * \return lowest of the supplied numbers
 */
double Min3(double x, double y, double z); //}}}

// Max3() //{{{
/**
 * \brief Function returning the highest number from three floats.
 *
 * \param [in] x   first double precision number
 * \param [in] y   second double precision number
 * \param [in] z   third double precision number
 * \return highest of the supplied numbers
 */
double Max3(double x, double y, double z); //}}}

// Sort3() //{{{
/**
 * \brief Function returning sorted numbers x < y < z.
 *
 * \param [in] in   first double precision number
 * \return sorted vector
 */
VECTOR Sort3(VECTOR in); //}}}

// SwapInt() //{{{
/**
 * \brief Function to swap two integers.
 *
 * \param [in] a   first integer to swap
 * \param [in] b   second integer to swap
 */
void SwapInt(int *a, int *b);
// }}}

// SwapDouble() //{{{
/**
 * \brief Function to swap two doubles.
 *
 * \param [in] a   first integer to swap
 * \param [in] b   second integer to swap
 */
void SwapDouble(double *a, double *b);
// }}}

// SwapBool() //{{{
/**
 * \brief Function to swap two booleans.
 *
 * \param [in] a   first integer to swap
 * \param [in] b   second integer to swap
 */
void SwapBool(bool *a, bool *b);
// }}}

// SortArray() //{{{
/**
 * \brief Function to sort an integer array.
 *
 * \param [out] array   integer array to sort
 * \param [in]  length  array length
 * \param [in]  mode    0 for ascending order, 1 for descending order
 */
void SortArray(int *array, int length, int mode); //}}}

// ReadAndSplitLine  //{{{
bool ReadAndSplitLine(FILE *fr, int *words, char split[SPL_STR][SPL_LEN]); //}}}

// ReadLine() //{{{
bool ReadLine(FILE *fr, char line[LINE]); //}}}

bool ReadAndSplitLine2(FILE *fr, int *words, char *split[SPL_STR]);

// SplitLine() //{{{
/*
 * \brief Function to split provided line.
 *
 * \param [out] out    array of strings
 * \param [in]  line   string to split
 * \return number of strings in the line
 */
int SplitLine(char out[SPL_STR][SPL_LEN], char *line, const char *delim); //}}}

// SplitLine2() //{{{
/*
 * \brief Function to split provided line.
 *
 * \param [out] out    array of strings
 * \param [in]  line   string to split
 * \return number of strings in the line
 */
int SplitLine2(char *out[SPL_STR], int strings, char *line, const char *delim);
 //}}}

// TrimLine() //{{{
/**
 * \brief Function to trim whitespace from
 * the beginning and end of a string.
 *
 * \param line [in]   string to trim
 *
 * \return trimmed string
 */
char* TrimLine(char *line); //}}}

void PrintCommand(FILE *ptr, int argc, char *argv[]);

// changing colour of the text (only for cli output)
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
char *ColourReset();

void SafeStrcat(char **out, char *in, int initial_size);

FILE *OpenFile(char *file, char *mode);

// TODO remove //{{{
void ColourChange(int a, char *colour);
void ColourReset_old(int a);
 //}}}
#endif
