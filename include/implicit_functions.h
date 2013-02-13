//  (C) Copyright Kirill Lykov 2012.
//
// Distributed under the FreeBSD Software License (See accompanying file license.txt)

#ifndef IMPLICITFUNCTIONS_H_
#define IMPLICITFUNCTIONS_H_

#include <stddef.h> //for size_t
#include <vector>
#include <boost/smart_ptr.hpp>
#include "box.h"
#include "math_vector.h"

/**
 * various implicit functions and operations on them
 */
namespace ls
{
using geometry_utils::MathVector;

enum Dimension
{
  xDim,
  yDim,
  zDim
};

template<typename T>
class IImplicitFunction
{
public:
  IImplicitFunction() {}
  virtual double compute(const MathVector<T, 3>& p) const = 0;
  virtual ~IImplicitFunction() {}

  typedef boost::shared_ptr< const IImplicitFunction<T> > Ptr;
private:
  IImplicitFunction(const IImplicitFunction&);
  IImplicitFunction& operator= (const IImplicitFunction&);
};

template<typename T>
class AxialCylinder : public IImplicitFunction<T>
{
  size_t m_i, m_j;
  T m_radius;
public:
  AxialCylinder(Dimension d, T radius);
  T compute(const MathVector<T, 3>& p) const;
};

template<typename T>
class NonAxialCylinder : public IImplicitFunction<T>
{
  MathVector<T, 3> m_pointOnLine;
  MathVector<T, 3> m_direction;
  T m_r, m_length;
public:
  NonAxialCylinder(const MathVector<T, 3>& pointOnLine, const MathVector<T, 3>& direction, T r, T length);
  T compute(const MathVector<T, 3>& point) const;
};

template<typename T>
class Sphere : public IImplicitFunction<T>
{
  MathVector<T, 3> m_center;
  T m_r;
public:
  Sphere(const MathVector<T, 3>& center, T r);
  T compute(const MathVector<T, 3>& point) const;
};

template<typename T>
class AxialTorus : public IImplicitFunction<T>
{
  size_t m_i, m_j, m_k; //indexes used for having one code for 3 axial-central toruses
  ls::geometry_utils::MathVector3D m_center; //center of the torus
  T m_torusRadius; // distance from the center of the tube to the center of the torus
  T m_tubeRadius; // radius of the tube
  T m_stretchCoef; // stretching coefficient for the torus in the 2nd dimension
public:
  AxialTorus(Dimension dim, const MathVector<T, 3>& center, T torusRadius, T tubeRadius, T stretchCoef = 1.0);
  T compute(const MathVector<T, 3>& point) const;
};

/**
 * cut the figure by plane. if a point is above the plane, the value will not be computed
 */
template<typename T>
class CutByPlane : public IImplicitFunction<T>
{
  typedef typename IImplicitFunction<T>::Ptr IImplicitFunctionPtr;
  const IImplicitFunctionPtr m_func;
  MathVector<T, 3> m_normal, m_pointOnPlane;
public:
  CutByPlane(const IImplicitFunctionPtr func, const MathVector<T, 3>& normal, const MathVector<T, 3>& point);
  T compute(const MathVector<T, 3>& point) const;
};

template<typename T>
class Union : public IImplicitFunction<T>
{
  typedef typename IImplicitFunction<T>::Ptr IImplicitFunctionPtr;
  typedef std::vector<IImplicitFunctionPtr> Functions;
  typedef typename std::vector<IImplicitFunctionPtr>::const_iterator FIterator;
  Functions m_functions;
public:
  Union(const Functions& functions);
  T compute(const MathVector<T, 3>& point) const;
};

template<typename T>
class Intersection : public IImplicitFunction<T>
{
  typedef typename IImplicitFunction<T>::Ptr IImplicitFunctionPtr;
  typedef std::vector<IImplicitFunctionPtr> Functions;
  typedef typename std::vector<IImplicitFunctionPtr>::const_iterator FIterator;
  Functions m_functions;
public:
  Intersection(const Functions& functions);
  T compute(const MathVector<T, 3>& point) const;
};

template<typename T>
class Difference : public IImplicitFunction<T>
{
  typedef typename IImplicitFunction<T>::Ptr IImplicitFunctionPtr;
  typedef std::vector<IImplicitFunctionPtr> Functions;
  typedef typename std::vector<IImplicitFunctionPtr>::const_iterator FIterator;
  Functions m_functions;
public:
  Difference(const Functions& functions);
  T compute(const MathVector<T, 3>& point) const;
};

// Methods bodies
template<typename T>
AxialCylinder<T>::AxialCylinder(Dimension d, T radius)
: m_i(0), m_j(0), m_radius(radius)
{
  if (d == xDim) {
    m_i = 1;
    m_j = 2;
  }
  else if (d == yDim) {
    m_i = 0;
    m_j = 2;
  }
  else {
    m_i = 0;
    m_j = 1;
  }
}

template<typename T>
T AxialCylinder<T>::compute(const MathVector<T, 3>& point) const
{
  T res = sqrt(pow(point.getCoord(m_i), 2) + pow(point.getCoord(m_j), 2) ) - m_radius;
  return res;
}

template<typename T>
NonAxialCylinder<T>::NonAxialCylinder(const MathVector<T, 3>& pointOnLine, const MathVector<T, 3>& direction, T r, T length)
: m_pointOnLine(pointOnLine), m_direction(direction), m_r(r), m_length(length)
{
  m_direction.normalize();
}

template<typename T>
T NonAxialCylinder<T>::compute(const MathVector<T, 3>& point) const
{
  // find projection of x on the line: xp = x0 + <x - x0, u>u
  // it is a length of the segment [x0, xint]
  T coeff = (point - m_pointOnLine) * m_direction;
  MathVector<T, 3> xp = m_pointOnLine + coeff * m_direction;

  // vector from point to projection on line
  MathVector<T, 3> diff = point - xp;

  T distToAxe = diff.getLength();

  T res = distToAxe - m_r;

  // take into account that cylinder is not infinite
  // l/2 because pointOnLine is in the center of the cylinder
  if (fabs(coeff) > m_length/2.) {

    T lextra = fabs(coeff) - m_length/2.;

    if (res < 0.0)
      res = lextra;
    else
      res = sqrt( res * res + lextra * lextra );
  }

  return res;
}

template<typename T>
Sphere<T>::Sphere(const MathVector<T, 3>& center, T r)
  : m_center(center), m_r(r)
{
}

template<typename T>
T Sphere<T>::compute(const MathVector<T, 3>& point) const
{
  T scaledNorm2 = 0.0;
  for (size_t i = 0; i < 3; ++i) {
    scaledNorm2 +=  pow(point.getCoord(i) - m_center.getCoord(i), 2);
  }

  return sqrt(scaledNorm2) -  m_r;
}

template<typename T>
AxialTorus<T>::AxialTorus(Dimension dim, const MathVector<T, 3>& center, T torusRadius, T tubeRadius, T stretchCoef)
: m_center(center), m_torusRadius(torusRadius), m_tubeRadius(tubeRadius), m_stretchCoef(stretchCoef)
{
  if (dim == xDim) {
    // plane yz
    m_i = 2;
    m_j = 1;
    m_k = 0;
  }
  else if (dim == yDim) {
    // plane xz
    m_i = 2;
    m_j = 0;
    m_k = 1;
  }
  else {
    // plane xy
    m_i = 1;
    m_j = 0;
    m_k = 2;
  }
}

template<typename T>
T AxialTorus<T>::compute(const MathVector<T, 3>& point) const
{
  // get distance from tube-wall
  T scoeff = (point.getCoord(m_i) - m_center.getCoord(m_i)) / m_stretchCoef;
  T dist_plane = fabs( sqrt( pow(point.getCoord(m_j) - m_center.getCoord(m_j), 2) + scoeff * scoeff) - m_torusRadius);
  T dist_other = point.getCoord(m_k) - m_center.getCoord(m_k); //here +- coef will be slope
  T res = sqrt(dist_other * dist_other + dist_plane * dist_plane) - m_tubeRadius;
  return res;
}

template<typename T>
CutByPlane<T>::CutByPlane(const IImplicitFunctionPtr func,
    const MathVector<T, 3>& normal,
    const MathVector<T, 3>& point)
: m_func(func), m_normal(normal), m_pointOnPlane(point)
{
}

template<typename T>
T CutByPlane<T>::compute(const MathVector<T, 3>& point) const
{
  T distToPlane = (point - m_pointOnPlane) * m_normal;
  if (distToPlane <= 0.0) {
    return m_func->compute(point);
  }

  // Projection( point ) = point - dist * normal
  MathVector<T, 3> projection = point - distToPlane * m_normal;
  // func(Projection( point ))^2 + d^2
  return sqrt(pow(m_func->compute(projection), 2) + distToPlane * distToPlane);
}

template<typename T>
Union<T>::Union(const Functions& functions)
: m_functions(functions)
{
}

template<typename T>
T Union<T>::compute(const MathVector<T, 3>& point) const
{
  T res = std::numeric_limits<T>::max();
  for (FIterator it = m_functions.begin(); it != m_functions.end(); ++it) {
    res = std::min(res, (*it)->compute(point));
  }
  return res;
}

template<typename T>
Intersection<T>::Intersection(const Functions& functions)
: m_functions(functions)
{
}

template<typename T>
T Intersection<T>::compute(const MathVector<T, 3>& point) const
{
  T res = std::numeric_limits<T>::min();
  for (FIterator it = m_functions.begin(); it != m_functions.end(); ++it) {
    res = std::max(res, (*it)->compute(point));
  }
  return res;
}

template<typename T>
Difference<T>::Difference(const Functions& functions)
: m_functions(functions)
{
}

template<typename T>
T Difference<T>::compute(const MathVector<T, 3>& point) const
{
  T res = std::numeric_limits<T>::min();
  for (FIterator it = m_functions.begin(); it != m_functions.end(); ++it) {
    res = std::max(res, -(*it)->compute(point));
  }
  return res;
}


#define DEF_IMPLICIT_FUNCTIONS(type,postfix) \
typedef IImplicitFunction<type> IImplicitFunction##postfix; \
typedef AxialCylinder<type> AxialCylinder##postfix; \
typedef NonAxialCylinder<type> NonAxialCylinder##postfix; \
typedef Sphere<type> Sphere##postfix; \
typedef AxialTorus<type> AxialTorus##postfix; \
typedef CutByPlane<type> CutByPlane##postfix; \
typedef Union<type> Union##postfix; \
typedef Intersection<type> Intersection##postfix; \
typedef Difference<type> Difference##postfix;

DEF_IMPLICIT_FUNCTIONS(double,D);

typedef boost::shared_ptr<const IImplicitFunctionD> IImplicitFunctionDPtr;
typedef std::vector<IImplicitFunctionDPtr> FunctionsD;
typedef typename FunctionsD::const_iterator FIteratorD;

DEF_IMPLICIT_FUNCTIONS(float,F);

typedef boost::shared_ptr<const IImplicitFunctionF> IImplicitFunctionFPtr;
typedef std::vector<IImplicitFunctionFPtr> FunctionsF;
typedef typename FunctionsF::const_iterator FIteratorF;

}

#endif /* IMPLICITFUNCTIONS_H_ */

