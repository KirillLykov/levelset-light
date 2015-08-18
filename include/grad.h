//  (C) Copyright Kirill Lykov 2015.
//
// Distributed under the FreeBSD Software License (See accompanying file license.txt)

#ifndef GRAD_H_
#define GRAD_H_

#include <cstddef>
#include <algorithm>
#include <grid.h>
#include "box.h"
#include "raw_math_vector.h"

namespace ls
{
  template < typename T, template<typename X>  class AccessStrategy >
  class Grad : public AccessStrategy<T>
  {
    typedef AccessStrategy<T> _AS;
    typedef typename _AS::_Grid _Grid;

    geometry_utils::Box3D m_bbox; // bounding box for grid
    double h[3];

  public:
    Grad(const _Grid& grid)
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

    geometry_utils::MathVector3D compute_forward(const geometry_utils::MathVector3D& point) const
    {
      assert(m_bbox.inside(point));

      // work with Cartesian with origin in left bottom point of the domain
      // thus shift the input point. Then fin index of the cell where the point is.
      geometry_utils::MathVector3D relativePosition = point - m_bbox.getLow();

      size_t index[3];
      computeIndex(relativePosition, index);

      double val = _AS::getValue(index[0], index[1], index[2]);
      geometry_utils::MathVector3D grad;
      grad.setCoord(0, (_AS::getValue(index[0] + 1, index[1],     index[2])     - val)/h[0]);
      grad.setCoord(1, (_AS::getValue(index[0],     index[1] + 1, index[2])     - val)/h[1]);
      grad.setCoord(2, (_AS::getValue(index[0],     index[1],     index[2] + 1) - val)/h[2]);
      return grad;
    }

    geometry_utils::MathVector3D compute_central(const geometry_utils::MathVector3D& point) const
    {
      assert(m_bbox.inside(point));

      // work with Cartesian with origin in left bottom point of the domain
      // thus shift the input point. Then fin index of the cell where the point is.
      geometry_utils::MathVector3D relativePosition = point - m_bbox.getLow();

      size_t index[3];
      computeIndex(relativePosition, index);

      geometry_utils::MathVector3D grad;
      grad.setCoord(0, _AS::getValue(index[0] + 1, index[1],     index[2])     - _AS::getValue(index[0] - 1, index[1], index[2]));
      grad.setCoord(1, _AS::getValue(index[0],     index[1] + 1, index[2])     - _AS::getValue(index[0], index[1] - 1, index[2]));
      grad.setCoord(2, _AS::getValue(index[0],     index[1],     index[2] + 1) - _AS::getValue(index[0], index[1], index[2] - 1));

      grad.normalize();
      return grad;
    }
  };
}

#endif /* GRID_H_ */
