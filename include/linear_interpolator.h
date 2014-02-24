//  (C) Copyright Kirill Lykov 2012.
//
// Distributed under the FreeBSD Software License (See accompanying file license.txt)

#ifndef LINEAR_INTERPOLATOR_H_
#define LINEAR_INTERPOLATOR_H_

#include <cstddef>
#include <algorithm>
#include <grid.h>
#include "box.h"
#include "raw_math_vector.h"

namespace ls
{
  template < typename T, template<typename X>  class AccessStrategy >
  class LinearInterpolator : public AccessStrategy<T>
  {
    typedef AccessStrategy<T> _AS;
    typedef typename _AS::_Grid _Grid;

    geometry_utils::Box3D m_bbox; // bounding box for grid
    double h[3];

  public:
    LinearInterpolator(const _Grid& grid)
    : _AS(grid), m_bbox(grid.getBoundingBox())
    {
      for (size_t i = 0; i < 3; ++i)
        h[i] = m_bbox.getIthSize(i) / (_AS::m_grid.size(i) - 1.0);
    }

    double compute(double x, double y, double z) const
    {
      geometry_utils::MathVector3D point(x, y, z);
      return compute(point);
    }

    void computeIndex(const geometry_utils::MathVector3D& relativePosition, size_t* index) const
    {
      assert(relativePosition.getX() >= 0.0 && relativePosition.getY() >= 0.0 && relativePosition.getZ() >= 0.0);
      // index_x = floor( p_x / h_x )
      for (size_t i = 0; i < 3; ++i) {
        // coef := 1 / h
        double coef = (_AS::m_grid.size(i) - 1.0) / static_cast<double>(m_bbox.getIthSize(i));
        index[i] = _AS::mapIndex( static_cast<int>(relativePosition.getCoord(i) * coef), i );
      }
    }

    double compute(const geometry_utils::MathVector3D& point) const
    {
      assert(m_bbox.inside(point));

      // work with Cartesian with origin in left bottom point of the domain
      // thus shift the input point. Then fin index of the cell where the point is.
      // After that interpolate function value for the point.
      geometry_utils::MathVector3D relativePosition = point - m_bbox.getLow();

      size_t index[3];
      computeIndex(relativePosition, index);

      return trilin_interp(relativePosition, index);
    }

  private:

     /**
      * Interpolate using trilinear interpolation method
      * names of variables are from http://en.wikipedia.org/wiki/Trilinear_interpolation
      */
     double trilin_interp(const geometry_utils::MathVector3D& inputPoint, const size_t* index) const
     {
       geometry_utils::MathVector3D x0(index[0] * h[0], index[1] * h[1], index[2] * h[2]);
       geometry_utils::MathVector3D xd = inputPoint - x0;
       double c[2][2];
       for (size_t i = 0; i < 2; ++i) {
         for (size_t j = 0; j < 2; ++j) {
           c[i][j] = _AS::getValue(index[0], index[1] + i, index[2] + j) * (h[0] - xd.getX()) +
               _AS::getValue(index[0] + 1, index[1] + i, index[2] + j) * xd.getX();
         }
       }

       double c0 = c[0][0] * (h[1] - xd.getY()) + c[1][0] * xd.getY();
       double c1 = c[0][1] * (h[1] - xd.getY()) + c[1][1] * xd.getY();

       double res = 1.0 / h[0] / h[1] / h[2] * (c0 * (h[2] - xd.getZ()) + c1 * xd.getZ());
       return res;
     }
  };
}

#endif /* LINEAR_INTERPOLATOR_H_ */
