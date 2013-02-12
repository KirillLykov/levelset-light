//  (C) Copyright Kirill Lykov 2012.
//
// Distributed under the FreeBSD Software License (See accompanying file license.txt)

#ifndef RAW_MATH_VECTOR_H_
#define RAW_MATH_VECTOR_H_

#include <math.h>
#include <memory.h>

namespace ls
{
namespace geometry_utils
{
  /**
   * Operations on vectors applied to block of 3 doubles.
   * To be used when usage of MathVector is not efficient.
   * Code of operation is self-explanatory. In order to avoid name conflicts,
   * it is highly recommended to use these functions with the namespace name.
   * Always assumed that double* is a pointer to the head of 3 element array
   */
  namespace raw_math_vector
  {
    static const double linearTolerance = 1e-6; //TODO use tolerance from tolerance class

    inline bool isEqualLinear(double left, double right)
    {
      return fabs(left - right) < linearTolerance;
    }

    inline void copy(double* targetVector, const double* sourceVector)
    {
      memcpy(targetVector, sourceVector, 3 * sizeof(targetVector[0]));
    }

    inline void zero(double* targetVector)
    {
      memset(targetVector, 0, 3 * sizeof(targetVector[0]));
    }

    inline void add(double* targetVector, const double* sourceVector)
    {
      for (size_t i = 0; i < 3; ++i)
        targetVector[i] += sourceVector[i];
    }

    inline void add(double* resVector, const double* firstVector, const double* secondVector)
    {
      copy(resVector, firstVector);
      add(resVector, secondVector);
    }

    inline void substract(double* targetVector, const double* sourceVector)
    {
      for (size_t i = 0; i < 3; ++i)
        targetVector[i] -= sourceVector[i];
    }

    inline void substract(double* resVector, const double* firstVector, const double* secondVector)
    {
      copy(resVector, firstVector);
      substract(resVector, secondVector);
    }

    inline double length(const double* vector)
    {
      return sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);
    }

    inline double squaredLength(const double* vector)
    {
      return vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2];
    }

    inline double dotProduct(const double* firstVector, const double* secondVector)
    {
      double res = 0;
      for (size_t i = 0; i < 3; ++i)
        res += firstVector[i] * secondVector[i];
      return res;
    }

    inline void crossProduct(double* resVector, const double* v1, const double* v2)
    {
      resVector[0] = v1[1] * v2[2] - v1[2] * v2[1];
      resVector[1] = v1[2] * v2[0] - v1[0] * v2[2];
      resVector[2] = v1[0] * v2[1] - v1[1] * v2[0];
    }

    inline void scale(double s, double* vector)
    {
      vector[0] *= s;
      vector[1] *= s;
      vector[2] *= s;
    }

    inline void normalize(double* vector)
    {
      double scaleCoeff = 1.0 / length(vector);
      scale(scaleCoeff, vector);
    }
  }
}
}

#endif /* RAW_MATH_VECTOR_H_ */
