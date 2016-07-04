//  (C) Copyright Kirill Lykov 2015.
//
// Distributed under the FreeBSD Software License (See accompanying file license.txt)

#ifndef GRAD_H_
#define GRAD_H_

#include <cstddef>
#include <algorithm>
#include <grid.h>
#include <box.h>
#include <math_vector.h>
#include <linear_interpolator.h>

namespace ls
{
  template <typename T, template<typename X> class AccessStrategy>
  class Grad : public ls::LinearInterpolator<T,AccessStrategy>
  {
    typedef AccessStrategy<T> _AS;
    typedef ls::LinearInterpolator<T, AccessStrategy> _LI;
    typedef typename _LI::_Grid _Grid;
  public:
    typedef geometry_utils::MathVector<T,3> MV;

    Grad(const _Grid& grid)
    : _LI(grid)
    {
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

    MV compute_biased(const MV& point, const MV& vel) const
    {
      //assert(m_bbox.inside(point));

      // work with Cartesian with origin in left bottom point of the domain
      // thus shift the input point. Then fin index of the cell where the point is.
      MV relativePosition = _AS::getRelativePosition(point);

      size_t index[3];
      _LI::computeIndex(relativePosition, index);

      T val = _LI::getValue(index[0], index[1], index[2]);

      return MV(
              (vel.getX() > 0.0 ? df<1,0,0>(index, val) : db<1,0,0>(index, val))/_LI::h[0] ,
              (vel.getY() > 0.0 ? df<0,1,0>(index, val) : db<0,1,0>(index, val))/_LI::h[1] ,
              (vel.getZ() > 0.0 ? df<0,0,1>(index, val) : db<0,0,1>(index, val))/_LI::h[2] );
    }

    MV compute_backward(const MV& point) const
    {
      //assert(m_bbox.inside(point));

      // work with Cartesian with origin in left bottom point of the domain
      // thus shift the input point. Then fin index of the cell where the point is.
      MV relativePosition = _AS::getRelativePosition(point);

      size_t index[3];
      _LI::computeIndex(relativePosition, index);

      T val = _LI::getValue(index[0], index[1], index[2]);
      return MV(db<1,0,0>(index, val)/_LI::h[0],
                                          db<0,1,0>(index, val)/_LI::h[1],
                                          db<0,0,1>(index, val)/_LI::h[2]);
    }

    MV compute_forward(const MV& point) const
    {
      //assert(m_bbox.inside(point));

      // work with Cartesian with origin in left bottom point of the domain
      // thus shift the input point. Then fin index of the cell where the point is.
      MV relativePosition = _AS::getRelativePosition(point);

      size_t index[3];
      _LI::computeIndex(relativePosition, index);

      T val = _LI::getValue(index[0], index[1], index[2]);
      return MV(df<1,0,0>(index, val)/_LI::h[0],
                                   df<0,1,0>(index, val)/_LI::h[1],
                                   df<0,0,1>(index, val)/_LI::h[2]);
    }

    MV compute_central(const MV& point) const
    {
      //assert(m_bbox.inside(point));

      // work with Cartesian with origin in left bottom point of the domain
      // thus shift the input point. Then fin index of the cell where the point is.
      MV relativePosition = _AS::getRelativePosition(point);

      size_t index[3];
      _LI::computeIndex(relativePosition, index);

      MV grad;
      grad.setCoord(0, 0.5*(_LI::getValue(index[0] + 1, index[1],     index[2])     - _LI::getValue(index[0] - 1, index[1], index[2]))/_LI::h[0]);
      grad.setCoord(1, 0.5*(_LI::getValue(index[0],     index[1] + 1, index[2])     - _LI::getValue(index[0], index[1] - 1, index[2]))/_LI::h[0]);
      grad.setCoord(2, 0.5*(_LI::getValue(index[0],     index[1],     index[2] + 1) - _LI::getValue(index[0], index[1], index[2] - 1))/_LI::h[0]);

      //grad.normalize();
      return grad;
    }

    // expensive way to compute gradient, when |computeFGrad()| > 1
    //from the paper "Adaptively Sampled Distance Fields : A General Representation of Shape for Computer Graphics"
    MV compute_precise(const MV& point) const
    {
      geometry_utils::Box3D bb = _AS::m_grid.getBoundingBox();
      //assert(bb.inside(point));

      // work with Cartesian with origin in left bottom point of the domain
      // thus shift the input point. Then fin index of the cell where the point is.
      MV relativePosition = _AS::getRelativePosition(point);
      MV grad;
      for (size_t i = 0; i < 3; ++i) {
        T h =  static_cast<T>(bb.getIthSize(i)) / (_LI::m_grid.size(i) - T(1.0));
        T index = _LI::mapIndex( static_cast<int>(relativePosition.getCoord(i) / h ), i );

        MV p1 = point;
        p1.setCoord(i, index*h +  bb.getLow().getCoord(i));
        //bb.applyPBC(p1);

        MV p2 = point;
        p2.setCoord(i, (index+1)*h  +  bb.getLow().getCoord(i));
        //bb.applyPBC(p2);

        grad.setCoord(i, (_LI::compute(p2) - _LI::compute(p1)) / h);
      }

      return grad;
    }
  };
}

#endif /* GRID_H_ */
