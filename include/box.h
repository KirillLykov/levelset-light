//  (C) Copyright Kirill Lykov 2012.
//
// Distributed under the FreeBSD Software License (See accompanying file license.txt)

#ifndef BOX_H_
#define BOX_H_

#include <cstddef>
#include <stdexcept>
#include <math.h>
#include <assert.h>
#include "raw_math_vector.h"
#include "tolerance.h"

namespace ls
{
namespace geometry_utils
{
  class Box
  {
    double low[3];
    double top[3];
  public:

    Box()
    {
      raw_math_vector::zero(low);
      raw_math_vector::zero(top);
    }

    /**
     * Cube with center at origin
     */
    Box(double domainSz)
    {
      double half = domainSz / 2.0;
      for (size_t i = 0; i < 3; ++i) {
        low[i] = -half;
        top[i] = half;
      }
    }

    /**
     * Cuboid with center at origin
     */
    Box(const double* domainSz)
    {
      for (size_t i = 0; i < 3; ++i) {
        double half = domainSz[i] / 2.0;
        low[i] = -half;
        top[i] = half;
      }
    }

    Box(double domainSzX, double domainSzY, double domainSzZ)
    {
      double domainSz[] = {domainSzX, domainSzY, domainSzZ};
      for (size_t i = 0; i < 3; ++i) {
        double half = domainSz[i] / 2.0;
        low[i] = -half;
        top[i] = half;
      }
    }

    Box(const double* low, const double* top)
    {
      if (low[0] > top[0] || low[1] >= top[1] || low[2] >= top[2])
        throw std::logic_error("low must be bottomLeft point, while top - upper right");
      raw_math_vector::copy(this->low, low);
      raw_math_vector::copy(this->top, top);
    }

    Box(const Box& anotherBox)
    {
      raw_math_vector::copy(this->low, anotherBox.low);
      raw_math_vector::copy(this->top, anotherBox.top);
    }

    Box& operator= (const Box& anotherBox)
    {
      if (this == &anotherBox)
        return *this;

      raw_math_vector::copy(this->low, anotherBox.low);
      raw_math_vector::copy(this->top, anotherBox.top);

      return *this;
    }

    double getIthSize(size_t i) const
    {
      assert(i < 3);
      return fabs(top[i] - low[i]);
    }

    double getSizeX() const
    {
      return getIthSize(0);
    }

    double getSizeY() const
    {
      return getIthSize(1);
    }

    double getSizeZ() const
    {
      return getIthSize(2);
    }

    bool inside(const double* point) const
    {
      if (point[0] >= low[0] && point[0] <= top[0]
       && point[1] >= low[1] && point[1] <= top[1]
       && point[2] >= low[2] && point[2] <= top[2])
        return true;
      return false;
    }

    void getCenter(double* center) const
    {
      for (size_t i =0 ; i < 3; ++i)
        center[i] = 0.5 * (top[i] + low[i]);
    }

    double getVolume() const
    {
      return getSizeX() * getSizeY() * getSizeZ();
    }

    void shift(double* newOrigin)
    {
      raw_math_vector::add(low, newOrigin);
      raw_math_vector::add(top, newOrigin);
    }

    double getMaxSize() const
    {
      return std::max(getSizeX(), std::max(getSizeY(), getSizeZ()));
    }

    /**
     * return true if box volume is 0
     */
    bool isTrivial() const
    {
      //TODO use tolerance
      return getVolume() == 0.0;
    }

    //unsafe functions because a user can try to modify internal data
    const double* getLow() const
    {
      return low;
    }

    const double* getTop() const
    {
      return top;
    }

    bool operator==(const Box& anotherBox) const
    {
      double diff[3];
      raw_math_vector::substract(diff, this->low, anotherBox.low);
      double difflow = raw_math_vector::length(diff);

      raw_math_vector::substract(diff, this->top, anotherBox.top);
      double difftop = raw_math_vector::length(diff);

      if (difflow < Tolerance::globalTolerance && difftop < Tolerance::globalTolerance)
        return true;
      return false;
    }

    bool operator!=(const Box& anotherBox) const
    {
      return !(*this == anotherBox);
    }
  };
}
}

#endif /* BOX_H_ */
