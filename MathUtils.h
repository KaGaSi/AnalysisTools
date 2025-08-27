#ifndef MATHUTILS_H
#define MATHUTILS_H

#include <math.h>
#include "Globals.h"

// x^2
static inline double Square(const double x) {
  return x * x;
}
// x^3
static inline double Cube(const double x) {
  return x * x * x;
}
// squared Euclidean vector of 3D vector: x^2+y^2+z^2
static inline double SqVectLength(const vec3 a) {
  return Square(a.v[0]) + Square(a.v[1]) + Square(a.v[2]);
}
// Euclidean length of 3D vector: (x^2+y^2+z^2)^0.5
static inline double VectLength(const vec3 a) {
  return sqrt(SqVectLength(a));
}
// dot product of 3D vectors
static inline double Dot(const vec3 a, const vec3 b) {
  return a.v[0] * b.v[0] + a.v[1] * b.v[1] + a.v[2] * b.v[2];
}
// cosine of angle between two 3D vectors
static inline double CosAngle(const vec3 a, const vec3 b) {
  return Dot(a, b) / (VectLength(a) * VectLength(b));
}
static inline vec3 Vector(const vec3 a, const vec3 b) {
  vec3 c;
  for (int dd = 0; dd < 3; dd++) {
    c.v[dd] = a.v[dd] - b.v[dd];
  }
  return c;
}
static inline double AngleDegrees(const vec3 a, const vec3 b) {
  return acos(CosAngle(a, b)) * 180 / PI;
}

#endif
