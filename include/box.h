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
  template<class T, size_t Dim>
  class Box
  {
    MathVector<T, Dim> low;
    MathVector<T, Dim> top;
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
      for (size_t i = 0; i < Dim; ++i) {
        low.setCoord(i, -half);
        top.setCoord(i, half);
      }
    }

    /**
     * Cuboid with center at origin
     */
    Box(const double* domainSz)
    {
      for (size_t i = 0; i < Dim; ++i) {
        double half = domainSz[i] / 2.0;
        low.setCoord(i, -half);
        top.setCoord(i, half);
      }
    }

    Box(double domainSzX, double domainSzY, double domainSzZ = 0.0)
    {
      double domainSz[] = {domainSzX, domainSzY, domainSzZ};
      for (size_t i = 0; i < Dim; ++i) {
        double half = domainSz[i] / 2.0;
        low.setCoord(i, -half);
        top.setCoord(i, half);
      }
    }

    Box(const MathVector<T, Dim>& low, const MathVector<T, Dim>& top)
    {
      for (size_t i = 0; i < Dim; ++i) {
        if (low.getCoord(i) >= top.getCoord(i))
          throw std::logic_error("low must be bottomLeft point, while top - upper right");
      }
      this->low = low;
      this->top = top;
    }

    Box(const Box<T, Dim>& anotherBox)
    {
      this->low = anotherBox.low;
      this->top = anotherBox.top;
    }

    Box<T, Dim>& operator= (const Box<T, Dim>& anotherBox)
    {
      if (this == &anotherBox)
        return *this;

      this->low = anotherBox.low;
      this->top = anotherBox.top;

      return *this;
    }

    double getIthSize(size_t i) const
    {
      assert(i < Dim);
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

    //TODO instead of if use templates like in MathVector
    bool inside(const MathVector<T, Dim>& point) const
    {
      if (Dim == 3) {
        if (point.getX() >= low.getX() - Tolerance::globalTolerance
         && point.getX() <= top.getX() + Tolerance::globalTolerance
         && point.getY() >= low.getY() - Tolerance::globalTolerance
         && point.getY() <= top.getY() + Tolerance::globalTolerance
         && point.getZ() >= low.getZ() - Tolerance::globalTolerance
         && point.getZ() <= top.getZ() + Tolerance::globalTolerance)
          return true;
      } else {
        if (point.getX() >= low.getX() - Tolerance::globalTolerance
         && point.getX() <= top.getX() + Tolerance::globalTolerance
         && point.getY() >= low.getY() - Tolerance::globalTolerance
         && point.getY() <= top.getY() + Tolerance::globalTolerance)
          return true;
      }
      return false;
    }

    void getCenter(MathVector<T, Dim>& center) const
    {
      for (size_t i =0 ; i < Dim; ++i)
        center.setCoord(i, 0.5 * (top.getCoord(i) + low.getCoord(i)));
    }

    double getVolume() const
    {
      return getSizeX() * getSizeY() * getSizeZ();
    }

    void shift(const MathVector<T, Dim>& newOrigin)
    {
      low += newOrigin;
      top += newOrigin;
    }

    double getMaxSize() const
    {
      return Dim == 3 ? std::max(getSizeX(), std::max(getSizeY(), getSizeZ())) :
          std::max(getSizeX(), getSizeY());
    }

    /**
     * return true if box volume is 0
     */
    bool isTrivial() const
    {
      return Tolerance::close(getVolume(), 0.0);
    }

    //unsafe functions because a user can try to modify internal data
    MathVector<T, Dim> getLow() const
    {
      return low;
    }

    MathVector<T, Dim> getTop() const
    {
      return top;
    }

    bool operator==(const Box<T, Dim>& anotherBox) const
    {
      MathVector<T, Dim> diff = this->low - anotherBox.low;
      double difflow = diff.getLength();

      diff = this->top - anotherBox.top;
      double difftop = diff.getLength();

      if (difflow < Tolerance::globalTolerance && difftop < Tolerance::globalTolerance)
        return true;
      return false;
    }

    bool operator!=(const Box<T, Dim>& anotherBox) const
    {
      return !(*this == anotherBox);
    }
  };

  typedef Box<double, 2> Box2D;
  typedef Box<double, 3> Box3D;
}
}

#endif /* BOX_H_ */
