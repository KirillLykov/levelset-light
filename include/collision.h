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
#include <iomanip>
#include <cassert>
#include <iostream>

namespace ls
{
using namespace ls::geometry_utils;

template <typename T, template<typename X> class AccessStrategy>
class Collision : private ls::LinearInterpolator<T,AccessStrategy>
{
  typedef ls::LinearInterpolator<T, AccessStrategy> _LI;
  typedef ls::Grid3D<T> _LSGrid;
  typedef ls::Grad<T, AccessStrategy> _BasicGrad;

  _LSGrid* m_grid; // grid is stored somewhere else, don't clean it
  _BasicGrad m_grad;
  const T m_tolerance;
  T m_diagLenght;

public:
  MathVector<T, 3> computeGrad(const MathVector<T, 3>& point) const
  {
    MathVector<T, 3> grad = m_grad.compute_forward(point);
    if ( fabs(grad*grad - T(1.0)) > T(0.2)) {
      grad = m_grad.compute_precise(point);
      grad.normalize();
    }
    return grad;
  }

  Collision(_LSGrid* grid)
  : _LI(*grid), m_grid(grid), m_grad(*m_grid), m_tolerance(1e-4)
  {
    T res = 0.0;
    for (size_t i = 0; i < 3; ++i)
      res += pow(_LI::h[i], 2);
    m_diagLenght = sqrt(res);
  }

  // Just value in the left bottom corner of the cell where the point is
  T computeCheapSDF(const MathVector<T, 3>& point) const
  {
    ls::geometry_utils::MathVector<T, 3> relativePosition = _LI::getRelativePosition(point);;

    size_t index[3];
    _LI::computeIndex(relativePosition, index);
    assert(index[0] < m_grid->size(0) && index[1] < m_grid->size(1) && index[2] < m_grid->size(2));
    return  (*m_grid)(index[0], index[1], index[2]);
  }

  // if it is true than it might be that the point require bounce back
  bool cheapOutsideCheck(const MathVector<T, 3>& point) const
  {
    return computeCheapSDF(point) >= -m_diagLenght;
  }

  T computeSDF(const MathVector<T, 3>& point) const
  {
    return _LI::compute(point);
  }

  // it should be called rarely, happens due to numerical reasons
  void rescueParticle(T currsdf, MathVector<T, 3>& pos) const
  {
    MathVector<T, 3> grad = computeGrad(pos);
    for (size_t i = 0; i < 5; ++i) {
      T stepsize = std::max(m_tolerance, fabs(currsdf));
      pos -= grad*stepsize;
      currsdf = computeSDF(pos);
      if (currsdf < T(-2.0)*m_tolerance)
        break;
    }
  }

  // sometimes particle is on the border, we want to shift it inside to avoid numerical problems
  // it is close to the interface so use gradient computed earlier for the intersection point xstar
  void shiftInside(T currsdf, const MathVector<T, 3>& grad, MathVector<T, 3>& pos) const
  {
    T stepsize = T(8.0)*std::max(m_tolerance, fabs(currsdf));
    pos -= grad*stepsize;
  }

  void bounceBack(T currsdf, T dt, MathVector<T, 3>& pos, MathVector<T, 3>& vel)
  {
    //MathVector<T, 3> origpos = pos; MathVector<T, 3> origvel = vel; // to debug
    using namespace ls::geometry_utils::raw_math_vector;
    T subdt = dt;
    MathVector<T, 3> posOld, grad;
    size_t nmultipleReflections = 0; // to tackle multiple reflections
    do
    {
      assert(currsdf >= 0.0);

      posOld = pos - dt * vel;

      if (computeSDF(posOld) > 0.0) {
        rescueParticle(currsdf, posOld);
        pos = posOld;
        vel *= -1.0;
        return;
      }

      assert(computeSDF(posOld) <= 0.0);

      MathVector<T, 3> xstar = pos;
      T xstarSdf = currsdf;
      subdt = dt;
      // iterations of newton method t^(n+1)=t^n - phi(t^n)/phi'(t^n)
      for (size_t i = 0; i < 5; ++i)
      {
        grad = computeGrad(xstar);

        const T DphiDt = std::max(m_tolerance, grad * vel);

        assert(DphiDt > 0);

        subdt = std::min(dt, std::max(T(0.0), subdt - xstarSdf / DphiDt * (1.0 + m_tolerance)));

        MathVector<T, 3> xstarNew = posOld + subdt * vel;
        MathVector<T, 3> diffXstar = xstar - xstarNew;
        T diff2 = diffXstar * diffXstar;
        if (diff2 < m_tolerance*m_tolerance)
          break;

        xstar = xstarNew;
        xstarSdf = computeSDF(xstar);
      }

      const T lambda = 2.0 * subdt - dt;

      pos = posOld + lambda * vel;

      vel *= -1.0;
      dt -= subdt;
      currsdf = computeSDF(pos);
      ++nmultipleReflections;
    } while (currsdf >= 0.0 && nmultipleReflections < 5);

    if (currsdf > -m_tolerance && currsdf < 0.0)
    {
      shiftInside(currsdf, grad, pos);
      assert(computeSDF(pos) < 0);
    }

    // could not resolve
    if (currsdf > 0)
    {
      std::cout << "Could not bb, nm" << nmultipleReflections << "\n";
      pos = posOld;
      assert(computeSDF(pos) < 0);
    }

    return;
  }

  // assumed to be called only for those edges which are in proximity to the interface
  // idea is from Fuhrmann et al. "Distance Fields for Rapid Collision Detection in Physically Based Modeling"
  // left and right vtx forming an edge
  void bounceBackEdge(T dt, MathVector<T, 3>& posleft, MathVector<T, 3>& velleft, MathVector<T, 3>& posright, MathVector<T, 3>& velright)
  {
    assert(computeSDF(posleft) < 0.0 && computeSDF(posright) < 0.0);
    const MathVector<T, 3> posmiddle = 0.5*(posleft + posright);
    MathVector<T, 3> velmiddle = 0.5*(velleft + velright);

    if (cheapOutsideCheck(posmiddle)) {
      T dist = computeSDF(posmiddle);
      if (dist >= 0.0) {
        MathVector<T, 3> bbpos = posmiddle;
        bounceBack(dist, dt, bbpos, velmiddle);
        MathVector<T, 3> shift = bbpos - posmiddle;
        posleft += shift;
        velleft *= -1.0;
        posright += shift;
        velright *= -1.0;
      }
    }
  }

};
typedef Collision<double, ls::BasicReadAccessStrategy>  CollisionBasic;
typedef Collision<double, ls::PeriodicReadAS>  CollisionPeriodic;

}

#endif /* INCLUDE_COLLISION_H_ */
