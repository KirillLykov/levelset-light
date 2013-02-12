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

    geometry_utils::Box m_bbox; // bounding box for grid
    double h[3];

  public:
    LinearInterpolator(const geometry_utils::Box& box, const _Grid& grid)
    : _AS(grid), m_bbox(box)
    {
      for (size_t i = 0; i < 3; ++i)
        h[i] = m_bbox.getIthSize(i) / (_AS::m_grid.size(i) - 1.0);
    }

    double run(double x, double y, double z) const
    {
      double point[3] = {x, y, z};
      return run(point);
    }

    void computeIndex(const double* relativePosition, size_t* index) const
    {
      assert(relativePosition[0] >= 0.0 && relativePosition[1] >= 0.0 && relativePosition[2] >= 0.0);
      // index_x = floor( p_x / h_x )
      for (size_t i = 0; i < 3; ++i) {
        // coef := 1 / h
        double coef = (_AS::m_grid.size(i) - 1.0) / static_cast<double>(m_bbox.getIthSize(i));
        index[i] = _AS::mapIndex( static_cast<int>(relativePosition[i] * coef), i );
      }
    }

    double run(const double* point) const
    {
      // work with Cartesian with origin in left bottom point of the domain
      // thus shift the input point. Then fin index of the cell where the point is.
      // After that interpolate function value for the point.
      double relativePosition[3];
      geometry_utils::raw_math_vector::copy(relativePosition, point);
      geometry_utils::raw_math_vector::substract(relativePosition, m_bbox.getLow());

      size_t index[3];
      computeIndex(relativePosition, index);

      return trilin_interp(relativePosition, index);
    }

  private:

     /**
      * Interpolate using trilinear interpolation method
      * names of variables are from http://en.wikipedia.org/wiki/Trilinear_interpolation
      */
     double trilin_interp(const double* inputPoint, const size_t* index) const
     {
       double x0[] = {index[0] * h[0], index[1] * h[1], index[2] * h[2]};
       double xd[3];
       geometry_utils::raw_math_vector::substract(xd, inputPoint, x0);
       double c[2][2];
       for (size_t i = 0; i < 2; ++i) {
         for (size_t j = 0; j < 2; ++j) {
           c[i][j] = _AS::getValue(index[0], index[1] + i, index[2] + j) * (h[0] - xd[0]) +
               _AS::getValue(index[0] + 1, index[1] + i, index[2] + j) * xd[0];
         }
       }

       double c0 = c[0][0] * (h[1] - xd[1]) + c[1][0] * xd[1];
       double c1 = c[0][1] * (h[1] - xd[1]) + c[1][1] * xd[1];

       double res = 1.0 / h[0] / h[1] / h[2] * (c0 * (h[2] - xd[2]) + c1 * xd[2]);
       return res;
     }
  };
}

#endif /* LINEAR_INTERPOLATOR_H_ */
