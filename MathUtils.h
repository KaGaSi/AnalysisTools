#ifndef MATHUTILS_H
#define MATHUTILS_H

#include <math.h>

static inline double Square(const double x) {
  return x * x;
}
static inline double Cube(const double x) {
  return x * x * x;
}
static inline double SqVectLength(const double x[3]) {
  return Square(x[0]) + Square(x[1]) + Square(x[2]);
}
static inline double VectLength(const double x[3]) {
  return sqrt(SqVectLength(x));
}

#endif
