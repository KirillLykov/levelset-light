//  (C) Copyright Kirill Lykov 2012.
//
// Distributed under the FreeBSD Software License (See accompanying file license.txt)

#ifndef BASIC_ACCESS_STRATEGY_H_
#define BASIC_ACCESS_STRATEGY_H_

#include "grid.h"

namespace ls
{
  /**
   * Access strategy is used in algorithms working with Grid
   * such as interpolation and deserializing.
   *
   * @param inputIndex
   *  may be positive or negative and may be greater that grid size
   * @param dimInd
   *  Index of dim, dim <- {0,2}
   * @return valid index in range [0, grid_ith_size)
   * size_t mapIndex(int inputIndex, size_t dimInd) const;
   *
   * @param i
   * @param j
   * @param k
   *  i,j,k - indexes in grid space, supposed to be less than grid size
   * T getValue(size_t i, size_t j, size_t k) const
   *
   * BasicReadAccessStrategy is the most basic case, when grid consists of objects
   * of primitive types.
   */

  template<class T>
  class BasicReadAccessStrategy
  {
  protected:

    typedef Grid3D<T> _Grid;
    const _Grid& m_grid;

    explicit BasicReadAccessStrategy(const _Grid& grid)
    : m_grid(grid)
    {
    }

    geometry_utils::MathVector<T, 3> getRelativePosition(const geometry_utils::MathVector<T, 3>& point) const
    {
      return point - m_grid.getBoundingBox().getLow();
    }

    size_t mapIndex(int inputIndex, size_t dimInd) const
    {
      assert(inputIndex >= 0);
      return static_cast<size_t>(inputIndex);
    }

    T getValue(size_t i, size_t j, size_t k) const
    {
      // handle this condition to correctly compute value on the border
      if (i >= m_grid.size(0) || j >= m_grid.size(1) || k >= m_grid.size(2)) {
        return T(0.0);
      }
      return m_grid(i, j, k);
    }

	bool inside(const geometry_utils::MathVector<T, 3>& point) const 
	{
		return m_grid.getBoundingBox().inside(point);
	}
  };

  template<class T>
  class BasicWriteAccessStrategy
  {
  protected:
    typedef Grid3D<T> _Grid;
    _Grid& m_grid;

    explicit BasicWriteAccessStrategy(_Grid& grid)
    : m_grid(grid)
    {
    }

    void setValue(size_t i, size_t j, size_t k, const T* newValues)
    {
      assert(i < m_grid.size(0) && j < m_grid.size(1) && k < m_grid.size(2));
      m_grid(i, j, k) = *newValues;
    }
  };


  /**
   * AccessStrategy for periodic domain. Points may be out of the domain cube
   * so index at getValue is also out of range
   */
  template<class T>
  class PeriodicReadAS
  {
  protected:
    typedef ls::Grid3D<T> _Grid;
    const _Grid& m_grid;
    geometry_utils::Box3 m_bbox;

    explicit PeriodicReadAS(const _Grid& grid)
    : m_grid(grid)
    {
      m_bbox = m_grid.getBoundingBox();
    }

    geometry_utils::MathVector<T, 3> getRelativePosition(const geometry_utils::MathVector<T, 3>& point) const
    {
      geometry_utils::MathVector<T, 3> res = point;
      m_bbox.applyPBC(res);
      res -= m_bbox.getLow();
      res.max(T(0.0));
      return res;
    }

    size_t mapIndex(int inputIndex, size_t dimInd) const
    {
      // it is needed for the interpolation
      int gsize = static_cast<int>(m_grid.size(dimInd));
      inputIndex %= gsize;
      assert(inputIndex >= 0);
      //while (inputIndex < 0)
      //  inputIndex += gsize; //gsize - 1;
      return static_cast<size_t>(inputIndex);
    }

    T getValue(size_t i, size_t j, size_t k) const
    {
      // it can be greater, just remap assert (i < m_grid.size(0) && j < m_grid.size(1) && k < m_grid.size(2));
      return m_grid( mapIndex(i, 0), mapIndex(j, 1), mapIndex(k, 2) );
    }

	bool inside(const geometry_utils::MathVector<T, 3>& point) const 
	{
	  return true;			    
	}
  };

}

#endif /* BASIC_ACCESS_STRATEGY_H_ */
