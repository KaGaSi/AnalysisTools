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

// Length() //{{{
/*
 * Function to calculate vector's Euclidian length.
 *
 * [in] a .. vector
 * return: length of the vector
 */
double Length(VECTOR a); //}}}

// convert string into number if possible //{{{
/*
 * Functions to test and convert strings into numbers
 *
 * Parameters for all of them:
 * [in]  str .. string to test
 * [out] val .. output number of the proper type
 * return: 'true' if str start with the proper number, 'false' otherwise
 */
bool IsReal(char *str, double *val);
bool IsPosReal(char *str, double *val);
bool IsInteger(char *str, long *val);
bool IsPosInteger(char *str, long *val);
bool IsNatural(char *str, long *val);
 //}}}

// Min3() //{{{
/**
 * Function returning the lowest number from three floats.
 *
 * Parameters:
 * [in] x .. first double precision number
 * [in] y .. second double precision number
 * [in] z .. third double precision number
 * return: lowest of the supplied numbers
 */
double Min3(double x, double y, double z); //}}}

// Max3() //{{{
/**
 * Function returning the highest number from three floats.
 *
 * Parameters:
 * [in] x .. first double precision number
 * [in] y .. second double precision number
 * [in] z .. third double precision number
 * return: highest of the supplied numbers
 */
double Max3(double x, double y, double z); //}}}

// Sort3() //{{{
/**
 * Function returning sorted numbers x < y < z.
 *
 * Parameters:
 * [in] in .. vector to sort
 * return: sorted vector
 */
VECTOR SortVector(VECTOR in); //}}}

// swapping functions //{{{
/**
 * Functions to swap two numbers of given type.
 *
 * Parameters for all of them:
 * [in & out] a .. first number to swap
 * [in & out] b .. second number to swap
 */
void SwapInt(int *a, int *b);
void SwapDouble(double *a, double *b);
void SwapBool(bool *a, bool *b);
// }}}

// SortArray() //{{{
/**
 * Function to sort an integer array.
 *
 * Parameters:
 * [in & out] array .. integer array to sort
 *       [in] length .. array length
 *       [in] mode .. 0 for ascending order, 1 for descending order
 */
void SortArray(int *array, int length, int mode); //}}}

// ReadLine() //{{{
bool ReadLine(FILE *fr, int max_char, char *line); //}}}

// SplitLine() //{{{
/*
 * Function to split a string based on provided delimiters.
 *
 * Parameters:
 * [out] out .. array of strings
 *  [in] line .. string to split
 * return: number of created strings
 */
int SplitLine(int max_str, char *out[], char *line, const char *delim);
 //}}}

bool ReadAndSplitLine(FILE *fr, int max_char, char *line, int *words,
                      char *out[], int max_strings, const char *delim);

void PrintCommand(FILE *ptr, int argc, char *argv[]);

// changing colour the text for cli output //{{{
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
 //}}}

FILE *OpenFile(char *file, char *mode);

// TODO remove //{{{
void ColourChange(int a, char *colour);
void ColourReset_old(int a);
// IsReal_old() //{{{
/*
 * \brief Function to test if a string is a real number.
 *
 * \param [in] a   string to test
 * \return 'true' if a is double, 'false' otherwise
 */
bool IsReal_old(char *a); //}}}
// IsInteger_old() //{{{
/*
 * \brief Function to test if a string is a non-negative whole number.
 *
 * \param [in] a   string to test
 * \return 'true' if a is integer, 'false' otherwise
 */
bool IsInteger_old(char *a); //}}}
// IsPosReal_old() //{{{
/*
 * \brief Function to test if a string is a non-negative real number.
 *
 * \param [in] a   string to test
 * \return 'true' if a is non-negative double, 'false' otherwise
 */
bool IsPosReal_old(char *a); //}}}
// IsNatural_old()  //{{{
bool IsNatural_old(char *a); //}}}
// SplitLine_old() //{{{
/*
 * \brief Function to split provided line.
 *
 * \param [out] out    array of strings
 * \param [in]  line   string to split
 * \return number of strings in the line
 */
int SplitLine_old(char out[SPL_STR][SPL_LEN], char *line, const char *delim); //}}}
// TrimLine() //{{{
/**
 * Function to trim white space from the beginning and end of a string.
 *
 * Parameters:
 * [in] line .. string to trim
 * return: string without preceding or trailing white space
 */
char* TrimLine(char *line); //}}}
 //}}}

#endif
