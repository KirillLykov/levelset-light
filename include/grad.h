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

  //template < typename T, template<typename X>  class AccessStrategy >
  //class Grad : public AccessStrategy<T>
  template <typename T, template<typename X> class AccessStrategy>
  class Grad : public ls::LinearInterpolator<T,AccessStrategy>
  {
    typedef AccessStrategy<T> _AS;
    typedef ls::LinearInterpolator<T, AccessStrategy> _LI;
    typedef typename _LI::_Grid _Grid;

    geometry_utils::Box3D m_bbox; // bounding box for grid
    T h[3];

  public:
    Grad(const _Grid& grid)
    : _LI(grid), m_bbox(grid.getBoundingBox())
    {
      for (size_t i = 0; i < 3; ++i)
        h[i] = m_bbox.getIthSize(i) / (_LI::m_grid.size(i) - 1.0);
    }

    //T compute(T x, T y, T z) const
    //{
    //  geometry_utils::MathVector3D point(x, y, z);
    //  return compute(point);
    //}

    void computeIndex(const geometry_utils::MathVector3D& relativePosition, size_t* index) const
    {
      assert(relativePosition.getX() >= 0.0 && relativePosition.getY() >= 0.0 && relativePosition.getZ() >= 0.0);
      // index_x = floor( p_x / h_x )
      for (size_t i = 0; i < 3; ++i) {
        // coef := 1 / h
        T coef = (_LI::m_grid.size(i) - 1.0) / static_cast<T>(m_bbox.getIthSize(i));
        index[i] = _LI::mapIndex( static_cast<int>(relativePosition.getCoord(i) * coef), i );
      }
    }

    template<size_t i, size_t j, size_t k>
    T df(size_t index[3], T val) const
    {
      return _LI::getValue(index[0] + i, index[1] + j, index[2] + k) - val;
    }

    template<size_t i, size_t j, size_t k>
    T db(size_t index[3], T val) const
    {
      return val - _LI::getValue(index[0] - i, index[1] - j, index[2] - k);
    }

    geometry_utils::MathVector3D compute_biased(const geometry_utils::MathVector3D& point,
                                                const geometry_utils::MathVector3D& vel) const
    {
      assert(m_bbox.inside(point));

      // work with Cartesian with origin in left bottom point of the domain
      // thus shift the input point. Then fin index of the cell where the point is.
      geometry_utils::MathVector3D relativePosition = point - m_bbox.getLow();

      size_t index[3];
      computeIndex(relativePosition, index);

      T val = _LI::getValue(index[0], index[1], index[2]);

      return geometry_utils::MathVector3D(
                      (vel.getX() > 0.0 ? df<1,0,0>(index, val) : db<1,0,0>(index, val))/h[0] ,
                      (vel.getY() > 0.0 ? df<0,1,0>(index, val) : db<0,1,0>(index, val))/h[1] ,
                      (vel.getZ() > 0.0 ? df<0,0,1>(index, val) : db<0,0,1>(index, val))/h[2] );
    }

    geometry_utils::MathVector3D compute_backward(const geometry_utils::MathVector3D& point) const
    {
      assert(m_bbox.inside(point));

      // work with Cartesian with origin in left bottom point of the domain
      // thus shift the input point. Then fin index of the cell where the point is.
      geometry_utils::MathVector3D relativePosition = point - m_bbox.getLow();

      size_t index[3];
      computeIndex(relativePosition, index);

      T val = _LI::getValue(index[0], index[1], index[2]);
      return geometry_utils::MathVector3D(db<1,0,0>(index, val)/h[0],
                                          db<0,1,0>(index, val)/h[1],
                                          db<0,0,1>(index, val)/h[2]);
    }

    geometry_utils::MathVector3D compute_forward(const geometry_utils::MathVector3D& point) const
    {
      assert(m_bbox.inside(point));

      // work with Cartesian with origin in left bottom point of the domain
      // thus shift the input point. Then fin index of the cell where the point is.
      geometry_utils::MathVector3D relativePosition = point - m_bbox.getLow();

      size_t index[3];
      computeIndex(relativePosition, index);

      T val = _LI::getValue(index[0], index[1], index[2]);
      return geometry_utils::MathVector3D(df<1,0,0>(index, val)/h[0],
                                   df<0,1,0>(index, val)/h[1],
                                   df<0,0,1>(index, val)/h[2]);
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
      grad.setCoord(0, _LI::getValue(index[0] + 1, index[1],     index[2])     - _LI::getValue(index[0] - 1, index[1], index[2]));
      grad.setCoord(1, _LI::getValue(index[0],     index[1] + 1, index[2])     - _LI::getValue(index[0], index[1] - 1, index[2]));
      grad.setCoord(2, _LI::getValue(index[0],     index[1],     index[2] + 1) - _LI::getValue(index[0], index[1], index[2] - 1));

      grad.normalize();
      return grad;
    }

    // expensive way to compute gradient, when |computeFGrad()| > 1
    //from the paper "Adaptively Sampled Distance Fields : A General Representation of Shape for Computer Graphics"
    geometry_utils::MathVector3D compute_precise(const geometry_utils::MathVector3D& point) const
    {
      Box3D bb = _AS::m_grid.getBoundingBox();
      assert(bb.inside(point));

      // work with Cartesian with origin in left bottom point of the domain
      // thus shift the input point. Then fin index of the cell where the point is.
      geometry_utils::MathVector3D relativePosition = point - bb.getLow();
      geometry_utils::MathVector3D grad;
      for (size_t i = 0; i < 3; ++i) {
        T h =  static_cast<T>(bb.getIthSize(i)) / (_LI::m_grid.size(i) - T(1.0));
        T index = _LI::mapIndex( static_cast<int>(relativePosition.getCoord(i) / h ), i );

        geometry_utils::MathVector3D p1 = point;
        p1.setCoord(i, index*h +  bb.getLow().getCoord(i));
        //bb.applyPBC(p1);

        geometry_utils::MathVector3D p2 = point;
        p2.setCoord(i, (index+1)*h  +  bb.getLow().getCoord(i));
        //bb.applyPBC(p2);

        grad.setCoord(i, (_LI::compute(p2) - _LI::compute(p1)) / h);
      }

      return grad;
    }
  };
}

#endif /* GRID_H_ */
