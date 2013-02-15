//  (C) Copyright Kirill Lykov 2012.
//
// Distributed under the FreeBSD Software License (See accompanying file license.txt)

#ifndef BOX_H_
#define BOX_H_

#include <cstddef>
#include <stdexcept>
#include <math.h>
#include <assert.h>
#include "math_vector.h"
#include "tolerance.h"

namespace ls
{
namespace geometry_utils
{
  class Box
  {
    MathVector3D low;
    MathVector3D top;
  public:

    Box()
    {
    }

    /**
     * Cube with center at origin
     */
    Box(double domainSz)
    {
      double half = domainSz / 2.0;
      for (size_t i = 0; i < 3; ++i) {
        low.setCoord(i, -half);
        top.setCoord(i, half);
      }
    }

    /**
     * Cuboid with center at origin
     */
    Box(const double* domainSz)
    {
      for (size_t i = 0; i < 3; ++i) {
        double half = domainSz[i] / 2.0;
        low.setCoord(i, -half);
        top.setCoord(i, half);
      }
    }

    Box(double domainSzX, double domainSzY, double domainSzZ)
    {
      double domainSz[] = {domainSzX, domainSzY, domainSzZ};
      for (size_t i = 0; i < 3; ++i) {
        double half = domainSz[i] / 2.0;
        low.setCoord(i, -half);
        top.setCoord(i, half);
      }
    }

    Box(const MathVector3D& low, const MathVector3D& top)
    {
      if (low.getX() > top.getX() || low.getY() >= top.getY() || low.getZ() >= top.getZ())
        throw std::logic_error("low must be bottomLeft point, while top - upper right");
      this->low = low;
      this->top = top;
    }

    Box(const Box& anotherBox)
    {
      this->low = anotherBox.low;
      this->top = anotherBox.top;
    }

    Box& operator= (const Box& anotherBox)
    {
      if (this == &anotherBox)
        return *this;

      this->low = anotherBox.low;
      this->top = anotherBox.top;

      return *this;
    }

    double getIthSize(size_t i) const
    {
      assert(i < 3);
      return fabs(top.getCoord(i) - low.getCoord(i));
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
      if (point[0] >= low.getX() && point[0] <= top.getX()
       && point[1] >= low.getY() && point[1] <= top.getY()
       && point[2] >= low.getZ() && point[2] <= top.getZ())
        return true;
      return false;
    }


    bool inside(const MathVector3D& point) const
    {
      if (point.getX() >= low.getX() && point.getX() <= top.getX()
       && point.getY() >= low.getY() && point.getY() <= top.getY()
       && point.getZ() >= low.getZ() && point.getZ() <= top.getZ())
        return true;
      return false;
    }

    void getCenter(MathVector3D& center) const
    {
      for (size_t i =0 ; i < 3; ++i)
        center.setCoord(i, 0.5 * (top.getCoord(i) + low.getCoord(i)));
    }

    double getVolume() const
    {
      return getSizeX() * getSizeY() * getSizeZ();
    }

    void shift(const MathVector3D& newOrigin)
    {
      low += newOrigin;
      top += newOrigin;
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
      return Tolerance::close(getVolume(), 0.0);
    }

    //unsafe functions because a user can try to modify internal data
    MathVector3D getLow() const
    {
      return low;
    }

    MathVector3D getTop() const
    {
      return top;
    }

    bool operator==(const Box& anotherBox) const
    {
      MathVector3D diff = this->low - anotherBox.low;
      double difflow = diff.getLength();

      diff = this->top - anotherBox.top;
      double difftop = diff.getLength();

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
