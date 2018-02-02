//  (C) Copyright Kirill Lykov 2012.
//
// Distributed under the FreeBSD Software License (See accompanying file license.txt)

#ifndef TOLERANCE_H_
#define TOLERANCE_H_

#include <math.h>

namespace ls
{
/**
 * class to check if two floating point numbers are the same.
 * For instance, if tolerance is 10e-11, then 1.2e-11 and 1.4e-11 cannot be distinguished
 */
template<typename T>
class Close_absolut
{
  const T m_tolerance;
public:
  explicit Close_absolut(T tolerance)
  : m_tolerance(tolerance)
  {
  }

  bool operator() (T left, T right) const
  {
    T diff = fabs(left  - right);
    return diff < m_tolerance;
  }
};

/**
 * class to check if two floating point numbers are the same.
 * For instance, if tolerance is 10e-11, then 1.2e-11 and 1.4e-11 will be distingueshed
 */
template<typename T>
class Close_relative
{
  const T m_tolerance;
public:
  explicit Close_relative(T tolerance)
  : m_tolerance(tolerance)
  {
  }

  bool operator() (T left, T right) const
  {
    T diff = fabs(left  - right);
    return (diff < m_tolerance * fabs(left)) && (diff < m_tolerance * fabs(right));
  }
};

//to avoid creating cpp file
namespace {
  /**
   * wrapper around Close_at_tolerance used for sugar
   * TODO think about concurrency issues
   */
  template<typename T>
  struct Tolerance;
  
  template<>
#ifdef SINGLE_PRECISION
  struct Tolerance<float>
  {
  	typedef float T;
#else
  struct Tolerance<double>
  {
	typedef double T;
#endif
    static const T globalTolerance;
    static const Close_absolut<T> close_at_tol;

    static bool close(T left, T right)
    {
      return close_at_tol(left, right);
    }
  };

#ifdef SINGLE_PRECISION
  const float Tolerance<float>::globalTolerance = 1e-6f;
  const Close_absolut<float> Tolerance<float>::close_at_tol(Tolerance::globalTolerance);
#else
  const double Tolerance<double>::globalTolerance = 1e-9;
  const Close_absolut<double> Tolerance<double>::close_at_tol(Tolerance::globalTolerance);
#endif

}

}

#endif /* TOLERANCE_H_ */
