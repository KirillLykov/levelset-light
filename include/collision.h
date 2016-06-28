//  (C) Copyright Kirill Lykov 2016.
//
// Distributed under the FreeBSD Software License (See accompanying file license.txt)

#ifndef INCLUDE_COLLISION_H_
#define INCLUDE_COLLISION_H_

#include <grid.h>
#include <grad.h>
#include <linear_interpolator.h>
#include <basic_access_strategy.h>
#include <math_vector.h>
#include <typeinfo>
#include <limits>
#include <cassert>

namespace ls
{
using namespace ls::geometry_utils;
class Collision
{
  typedef ls::Grid3D<double> _LSGrid;
  typedef ls::Grad<double, ls::BasicReadAccessStrategy > _BasicGrad;
  typedef ls::LinearInterpolator<double, ls::BasicReadAccessStrategy > _Interpolator;

  _LSGrid* m_grid; // grid is stored somewhere else, don't clean it
  _Interpolator m_interpolator;
  _BasicGrad m_grad;

  MathVector3D computeFGrad(const MathVector3D& point) const
  {
    MathVector3D grad;
    if (!m_grid->getBoundingBox().inside(point))
      grad = MathVector3D(-std::numeric_limits<double>::infinity());
    else {
      grad = m_grad.compute_forward(point);
    }
    return grad;
  }

public:

  Collision(_LSGrid* grid)
  : m_grid(grid), m_interpolator(*m_grid), m_grad(*m_grid)
  {
  }

  // Just value in the left bottom corner of the cell where the point is
  double computeCheapSDF(const MathVector3D& point) const
  {
    double dist;
    if (!m_grid->getBoundingBox().inside(point))
     dist = -std::numeric_limits<double>::infinity();
    else {
     ls::geometry_utils::MathVector3D relativePosition = point - m_grid->getBoundingBox().getLow();

     size_t index[3];
     m_interpolator.computeIndex(relativePosition, index);
     assert(index[0] < m_grid->size(0) && index[1] < m_grid->size(1) && index[2] < m_grid->size(2));
     dist = (*m_grid)(index[0], index[1], index[2]);
    }
    return dist;
  }

  double computeSDF(const MathVector3D& point) const
  {
    double dist;
    if (!m_grid->getBoundingBox().inside(point))
      dist = -std::numeric_limits<double>::infinity();
    else
      dist = m_interpolator.compute(point);
    return dist;
  }

  void bounceBack(double currsdf, double dt, MathVector3D& pos, MathVector3D& vel)
  {
    using namespace ls::geometry_utils::raw_math_vector;
    const double tolerance = 1e-4;
    double subdt = dt;
    MathVector3D posOld;
    size_t nmultipleReflections = 0; // to tackle multiple reflections
    do
    {
      assert(currsdf >= 0.0);
      posOld = pos - dt * vel;

      assert(computeSDF(posOld) < 0.0);

      MathVector3D xstar = pos;
      double xstarSdf = currsdf;
      subdt = dt;
      // iterations of newton method t^(n+1)=t^n - phi(t^n)/phi'(t^n)
      for (size_t i = 0; i < 5; ++i)
      {
        MathVector3D grad = computeFGrad(xstar);

        assert(grad.getLength() < 1.1);
        const double DphiDt = std::max(tolerance, grad * vel);

        assert(DphiDt > 0);

        subdt = std::min(dt, std::max(0.0, subdt - xstarSdf / DphiDt * (1.0 + tolerance)));

        MathVector3D xstarNew = posOld + subdt * vel;
        double diff = (xstar - xstarNew).getLength();
        if (diff < tolerance)
          break;

        xstar = xstarNew;
        xstarSdf = computeSDF(xstar);
      }

      const double lambda = 2.0 * subdt - dt;

      pos = posOld + lambda * vel;

      vel *= -1.0;
      dt -= subdt;
      currsdf = computeSDF(pos);
      ++nmultipleReflections;
    } while (currsdf >= 0 && nmultipleReflections < 5);

    // to many reflections, could not resolve
    if (currsdf >= 0)
    {
      pos = posOld;
      assert(computeSDF(pos) < 0);
    }

    return;
  }

};


}

#endif /* INCLUDE_COLLISION_H_ */
