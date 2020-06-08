/**
 * \file
 * \brief Functions independent of the analysis utilities.
 */

#ifndef _GENERAL_H_
#define _GENERAL_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>

#define PI 3.141593 ///< value of pi
#define LINE 1024 ///< maximum length of an array (for strings)
#define SQR(x) ((x)*(x)) ///< macro for algebraic square
#define CUBE(x) ((x)*(x)*(x)) ///< macro for algebraic cube

// struct Vector //{{{
/**
 * \brief 3D vector of floats.
 */
typedef struct Vector {
  double x, y, z;
} Vector; //}}}

// struct LongVector //{{{
/**
 * \brief 3D vector of floats.
 */
typedef struct LongVector {
  long double x, y, z;
} LongVector; //}}}

// struct IntVector //{{{
/**
 * \brief 3D vector of integers.
 */
typedef struct IntVector {
  int x, y, z;
} IntVector; //}}}

// Length() //{{{
/*
 * \brief Function to calculate vector length.
 *
 * \param [in] a    vector
 * \return a's length
 */
double Length(Vector a); //}}}

// IsDouble() //{{{
/*
 * \brief Function to test if a string is a real number.
 *
 * \param [in] a   string to test
 * \return 'true' if a is double, 'false' otherwise
 */
bool IsDouble(char *a); //}}}

// IsPosDouble() //{{{
/*
 * \brief Function to test if a string is a non-negative real number.
 *
 * \param [in] a   string to test
 * \return 'true' if a is non-negative double, 'false' otherwise
 */
bool IsPosDouble(char *a); //}}}

// IsInteger() //{{{
/*
 * \brief Function to test if a string is a non-negative whole number.
 *
 * \param [in] a   string to test
 * \return 'true' if a is integer, 'false' otherwise
 */
bool IsInteger(char *a); //}}}

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
Vector Sort3(Vector in); //}}}

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
void SortArray(int **array, int length, int mode); //}}}

// SplitLine() //{{{
/*
 * \brief Function to split provided line.
 *
 * \param [out] out    array of strings
 * \param [in]  line   string to split
 * \return number of strings in the line
 */
int SplitLine(char out[30][100], char *line, char delim[8]); //}}}

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
#endif
