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
}

#endif /* BASIC_ACCESS_STRATEGY_H_ */
