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
static inline double SqVectLength(const vec3d a) {
  return Square(a.x) + Square(a.y) + Square(a.z);
}
// Euclidean length of 3D vector: (x^2+y^2+z^2)^0.5
static inline double VectLength(const vec3d a) {
  return sqrt(SqVectLength(a));
}
// dot product of 3D vectors
static inline double Dot(const vec3d a, const vec3d b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}
// cosine of angle between two 3D vectors
static inline double CosAngle(const vec3d a, const vec3d b) {
  return Dot(a, b) / (VectLength(a) * VectLength(b));
}
static inline vec3d Vector(const vec3d a, const vec3d b) {
  vec3d c;
  for (int dd = 0; dd < 3; dd++) {
    c.v[dd] = a.v[dd] - b.v[dd];
  }
  return c;
}
static inline double AngleDegrees(const vec3d a, const vec3d b) {
  return acos(CosAngle(a, b)) * 180 / PI;
}

#endif
