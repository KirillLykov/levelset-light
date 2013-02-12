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
class Close_absolut
{
  const double m_tolerance;
public:
  explicit Close_absolut(double tolerance)
  : m_tolerance(tolerance)
  {
  }

  bool operator() (double left, double right) const
  {
    double diff = fabs(left  - right);
    return diff < m_tolerance;
  }
};

/**
 * class to check if two floating point numbers are the same.
 * For instance, if tolerance is 10e-11, then 1.2e-11 and 1.4e-11 will be distingueshed
 */
class Close_relative
{
  const double m_tolerance;
public:
  explicit Close_relative(double tolerance)
  : m_tolerance(tolerance)
  {
  }

  bool operator() (double left, double right) const
  {
    double diff = fabs(left  - right);
    return (diff < m_tolerance * fabs(left)) && (diff < m_tolerance * fabs(right));
  }
};

//to avoid creating cpp file
namespace {
  /**
   * wrapper around Close_at_tolerance used for sugar
   * TODO think about concurrency issues
   */
  struct Tolerance
  {
    static const double globalTolerance;
    static const Close_absolut close_at_tol;

    static bool close(double left, double right)
    {
      return close_at_tol(left, right);
    }
  };

  const double Tolerance::globalTolerance = 10e-10;
  const Close_absolut Tolerance::close_at_tol(Tolerance::globalTolerance);
}

}

#endif /* TOLERANCE_H_ */
