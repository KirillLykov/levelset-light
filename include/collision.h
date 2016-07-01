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

namespace ls
{
using namespace ls::geometry_utils;
typedef ls::LinearInterpolator<double, ls::BasicReadAccessStrategy > _Interpolator;
class Collision : private _Interpolator
{
  typedef ls::Grid3D<double> _LSGrid;
  typedef ls::Grad<double, ls::BasicReadAccessStrategy > _BasicGrad;

  _LSGrid* m_grid; // grid is stored somewhere else, don't clean it
  //_Interpolator m_interpolator;
  _BasicGrad m_grad;
  const double m_tolerance;

public:
  MathVector3D computeCGrad(const MathVector3D& point) const
  {
    MathVector3D grad;
    if (!m_grid->getBoundingBox().inside(point))
      grad = MathVector3D(-std::numeric_limits<double>::infinity());
    else {
      grad = m_grad.compute_central(point);
    }
    return grad;
  }

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

  // expensive way to compute gradient, when |computeFGrad()| > 1
  //from the paper "Adaptively Sampled Distance Fields : A General Representation of Shape for Computer Graphics"
  geometry_utils::MathVector3D computePreciseGrad(const geometry_utils::MathVector3D& point) const
  {
    Box3D bb = m_grid->getBoundingBox();
    assert(bb.inside(point));

    // work with Cartesian with origin in left bottom point of the domain
    // thus shift the input point. Then fin index of the cell where the point is.
    geometry_utils::MathVector3D relativePosition = point - bb.getLow();
    geometry_utils::MathVector3D grad;
    for (size_t i = 0; i < 3; ++i) {
      double h =  static_cast<double>(bb.getIthSize(i)) / (m_grid->size(i) - 1.0);
      double index = mapIndex( static_cast<int>(relativePosition.getCoord(i) / h ), i );

      geometry_utils::MathVector3D p1 = point;
      p1.setCoord(i, index*h +  bb.getLow().getCoord(i));
      bb.applyPBC(p1);

      geometry_utils::MathVector3D p2 = point;
      p2.setCoord(i, (index+1)*h  +  bb.getLow().getCoord(i));
      bb.applyPBC(p2);

      grad.setCoord(i, (compute(p2) - compute(p1)) / h);
    }

    return grad;
  }

  MathVector3D computeGrad(const MathVector3D& point) const
  {
    MathVector3D grad = computeFGrad(point);
    if (grad*grad > 1.1*1.1) { // gradCheap.getLength() > 1.1
      grad = computePreciseGrad(point);
    }
    grad.normalize();
    assert( fabs(grad.getLength() - 1.0) < 1e-6);
    return grad;
  }

  Collision(_LSGrid* grid)
  : _Interpolator(*grid), m_grid(grid), m_grad(*m_grid), m_tolerance(1e-4)
  {
  }

  void applyPBC(double coord[3]) const
  {
    m_grid->getBoundingBox().applyPBC(coord);
  }

  void applyPBC(MathVector3D& coord) const
  {
    m_grid->getBoundingBox().applyPBC(coord);
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
     computeIndex(relativePosition, index);
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
      dist = compute(point);

    return dist;
    //return fabs(dist) < m_tolerance ? 0.0 : dist;
  }

  // it should be called rarely, happens due to numerical reasons
  void rescueParticle(double currsdf, MathVector3D& pos) const
  {
    MathVector3D grad = computeGrad(pos);
    for (size_t i = 0; i < 5; ++i) {
      double stepsize = std::max(m_tolerance, fabs(currsdf));
      pos -= grad*stepsize;
      currsdf = computeSDF(pos);
      if (currsdf < -2.0*m_tolerance)
        break;
    }
  }

  // sometimes particle is on the border, we want to shift it inside to avoid numerical problems
  // it is close to the interface so use gradient computed earlier for the intersection point xstar
  void shiftInside(double currsdf, const MathVector3D& grad, MathVector3D& pos) const
  {
    double stepsize = 8.0*std::max(m_tolerance, fabs(currsdf));
    pos -= grad*stepsize;
  }

  void bounceBack(double currsdf, double dt, MathVector3D& pos, MathVector3D& vel)
  {
    MathVector3D origpos = pos; MathVector3D origvel = vel;
    using namespace ls::geometry_utils::raw_math_vector;
    double subdt = dt;
    MathVector3D posOld, grad;
    size_t nmultipleReflections = 0; // to tackle multiple reflections
    do
    {
      if (currsdf < 0) {
        std::cout << "currsdf < 0| " << currsdf << " " << dt << " " << nmultipleReflections << "\n";
        std::cout << origpos.getX() << " " << origpos.getY() << " " << origpos.getZ() <<"\n";
        std::cout << origvel.getX() << " " << origvel.getY() << " " << origvel.getZ() <<"\n";
      }

      assert(currsdf >= 0.0);

      posOld = pos - dt * vel;

      if (computeSDF(posOld) > 0) {
        std::cout << "computeSDF(posOld) >= 0| " << computeSDF(posOld) << " " << dt << " " << nmultipleReflections << "\n";
        std::cout << origpos.getX() << " " << origpos.getY() << " " << origpos.getZ() <<"\n";
        std::cout << origvel.getX() << " " << origvel.getY() << " " << origvel.getZ() <<"\n";

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

        const double DphiDt = std::max(m_tolerance, grad * vel);

        assert(DphiDt > 0);

        subdt = std::min(dt, std::max(0.0, subdt - xstarSdf / DphiDt * (1.0 + m_tolerance)));

        MathVector3D xstarNew = posOld + subdt * vel;
        MathVector3D diffXstar = xstar - xstarNew;
        applyPBC(diffXstar);
        double diff = diffXstar.getLength();
        if (diff < m_tolerance)
          break;

        xstar = xstarNew;
        applyPBC(xstar);
        xstarSdf = computeSDF(xstar);
      }

      const double lambda = 2.0 * subdt - dt;

      pos = posOld + lambda * vel;

      vel *= -1.0;
      dt -= subdt;
      currsdf = computeSDF(pos);
      ++nmultipleReflections;
    } while (currsdf >= 0 && nmultipleReflections < 5);

    if (currsdf > -m_tolerance && currsdf < 0.0)
    {
      //std::cout << "fabs(currsdf) <= m_tolerance " << currsdf << " " << nmultipleReflections << "\n";
      shiftInside(currsdf, grad, pos);
      //std::cout << std::setprecision(16) << origpos.getX() << " " << origpos.getY() << " " << origpos.getZ() <<"\n";
      //std::cout << std::setprecision(16) << origvel.getX() << " " << origvel.getY() << " " << origvel.getZ() <<"\n";
      assert(computeSDF(pos) < 0);
    }

    // could not resolve
    if (currsdf > 0)
    {
      std::cout << "currsdf >= 0" << nmultipleReflections << "\n";
      pos = posOld;
      assert(computeSDF(pos) < 0);
    }

    return;
  }

};


}

#endif /* INCLUDE_COLLISION_H_ */
