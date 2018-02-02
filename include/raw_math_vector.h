//  (C) Copyright Kirill Lykov 2012.
//
// Distributed under the FreeBSD Software License (See accompanying file license.txt)

#ifndef RAW_MATH_VECTOR_H_
#define RAW_MATH_VECTOR_H_

#include <math.h>
#include <memory.h>
#include <cstring>
#include "tolerance.h"

namespace ls
{
namespace geometry_utils
{
  /**
   * Operations on vectors applied to block of 3 Ts/floats.
   * To be used when usage of MathVector is not efficient.
   * Code of operation is self-explanatory. In order to avoid name conflicts,
   * it is highly recommended to use these functions with the namespace name.
   * Always assumed that T* is a pointer to the head of 3 element array
   */
  namespace raw_math_vector
  {
  	template<typename T>
    inline bool isEqualLinear(T left, T right)
    {
      return Tolerance::close(left, right);
    }

    template<typename T>
	inline void copy(T* targetVector, const T* sourceVector)
    {
      memcpy(targetVector, sourceVector, 3 * sizeof(targetVector[0]));
    }

	template<typename T>
    inline void zero(T* targetVector)
    {
      memset(targetVector, 0, 3 * sizeof(targetVector[0]));
    }

	template<typename T>
    inline void add(T* targetVector, const T* sourceVector)
    {
      for (size_t i = 0; i < 3; ++i)
        targetVector[i] += sourceVector[i];
    }

	template<typename T>
    inline void add(T* resVector, const T* firstVector, const T* secondVector)
    {
      copy(resVector, firstVector);
      add(resVector, secondVector);
    }

	template<typename T>
    inline void substract(T* targetVector, const T* sourceVector)
    {
      for (size_t i = 0; i < 3; ++i)
        targetVector[i] -= sourceVector[i];
    }

	template<typename T>
    inline void substract(T* resVector, const T* firstVector, const T* secondVector)
    {
      copy(resVector, firstVector);
      substract(resVector, secondVector);
    }

	template<typename T>
    inline T length(const T* vector)
    {
      return sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);
    }

	template<typename T>
    inline T squaredLength(const T* vector)
    {
      return vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2];
    }

	template<typename T>
    inline T dotProduct(const T* firstVector, const T* secondVector)
    {
      T res = 0;
      for (size_t i = 0; i < 3; ++i)
        res += firstVector[i] * secondVector[i];
      return res;
    }

	template<typename T>
    inline void crossProduct(T* resVector, const T* v1, const T* v2)
    {
      resVector[0] = v1[1] * v2[2] - v1[2] * v2[1];
      resVector[1] = v1[2] * v2[0] - v1[0] * v2[2];
      resVector[2] = v1[0] * v2[1] - v1[1] * v2[0];
    }

	template<typename T>
    inline void scale(T s, T* vector)
    {
      vector[0] *= s;
      vector[1] *= s;
      vector[2] *= s;
    }

	template<typename T>
    inline void normalize(T* vector)
    {
      T scaleCoeff = 1.0 / length(vector);
      scale(scaleCoeff, vector);
    }
  }
}
}

#endif /* RAW_MATH_VECTOR_H_ */
