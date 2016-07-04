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
  MathVector3D computeGrad(const MathVector3D& point) const
  {
    MathVector3D grad = m_grad.compute_forward(point);
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

  /*
  void applyPBC(double coord[3]) const
  {
    m_grid->getBoundingBox().applyPBC(coord);
  }

  void applyPBC(MathVector3D& coord) const
  {
    m_grid->getBoundingBox().applyPBC(coord);
  }
  */

  // Just value in the left bottom corner of the cell where the point is
  double computeCheapSDF(const MathVector3D& point) const
  {
    T dist;
    if (!m_grid->getBoundingBox().inside(point))
     dist = -std::numeric_limits<double>::infinity();
    else {
     ls::geometry_utils::MathVector3D relativePosition = point - m_grid->getBoundingBox().getLow();

     size_t index[3];
     _LI::computeIndex(relativePosition, index);
     assert(index[0] < m_grid->size(0) && index[1] < m_grid->size(1) && index[2] < m_grid->size(2));
     dist = (*m_grid)(index[0], index[1], index[2]);
    }
    return dist;
  }

  // if it is true than it might be that the point require bounce back
  bool cheapOutsideCheck(const MathVector3D& point) const
  {
    return computeCheapSDF(point) >= -m_diagLenght;
  }

  T computeSDF(const MathVector3D& point) const
  {
    return _LI::compute(point);
    //return fabs(dist) < m_tolerance ? 0.0 : dist;
  }

  // it should be called rarely, happens due to numerical reasons
  void rescueParticle(T currsdf, MathVector3D& pos) const
  {
    MathVector3D grad = computeGrad(pos);
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
  void shiftInside(T currsdf, const MathVector3D& grad, MathVector3D& pos) const
  {
    T stepsize = T(8.0)*std::max(m_tolerance, fabs(currsdf));
    pos -= grad*stepsize;
  }

  void bounceBack(T currsdf, double dt, MathVector3D& pos, MathVector3D& vel)
  {
    MathVector3D origpos = pos; MathVector3D origvel = vel; // to debug
    using namespace ls::geometry_utils::raw_math_vector;
    T subdt = dt;
    MathVector3D posOld, grad;
    size_t nmultipleReflections = 0; // to tackle multiple reflections
    do
    {
      if (currsdf < 0) {
        std::cout << "currsdf < 0| " << currsdf << " " << dt << " " << nmultipleReflections << "\n";
      }

      assert(currsdf >= 0.0);

      posOld = pos - dt * vel;

      if (computeSDF(posOld) > 0.0) {
        std::cout << "computeSDF(posOld) >= 0| " << computeSDF(posOld) << " " << dt << " " << nmultipleReflections << "\n";

        rescueParticle(currsdf, posOld);
        pos = posOld;
        vel *= -1.0;
        return;
      }

      assert(computeSDF(posOld) <= 0.0);

      MathVector3D xstar = pos;
      double xstarSdf = currsdf;
      subdt = dt;
      // iterations of newton method t^(n+1)=t^n - phi(t^n)/phi'(t^n)
      for (size_t i = 0; i < 5; ++i)
      {
        grad = computeGrad(xstar);

        const T DphiDt = std::max(m_tolerance, grad * vel);

        assert(DphiDt > 0);

        subdt = std::min(dt, std::max(T(0.0), subdt - xstarSdf / DphiDt * (1.0 + m_tolerance)));

        MathVector3D xstarNew = posOld + subdt * vel;
        MathVector3D diffXstar = xstar - xstarNew;
        //applyPBC(diffXstar);
        T diff2 = diffXstar * diffXstar;
        if (diff2 < m_tolerance*m_tolerance)
          break;

        xstar = xstarNew;
        //applyPBC(xstar);
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
      //std::cout << "fabs(currsdf) <= m_tolerance " << currsdf << " " << nmultipleReflections << "\n";
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
  void bounceBackEdge(double dt, MathVector3D& posleft, MathVector3D& velleft, MathVector3D& posright, MathVector3D& velright)
  {
    assert(computeSDF(posleft) < 0.0 && computeSDF(posright) < 0.0);
    const MathVector3D posmiddle = 0.5*(posleft + posright);
    MathVector3D velmiddle = 0.5*(velleft + velright);

    if (cheapOutsideCheck(posmiddle)) {
      MathVector3D bbpos = posmiddle;
      bounceBack(computeSDF(bbpos), dt, bbpos, velmiddle);
      MathVector3D shift = bbpos - posmiddle;
      posleft += shift;
      velleft *= -1.0;
      posright += shift;
      velright *= -1.0;
    }
  }

};
typedef Collision<double, ls::BasicReadAccessStrategy>  CollisionBasic;
typedef Collision<double, ls::PeriodicReadAS>  CollisionPeriodic;

}

#endif /* INCLUDE_COLLISION_H_ */
