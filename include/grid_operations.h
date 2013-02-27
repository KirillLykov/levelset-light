//  (C) Copyright Kirill Lykov 2012.
//
// Distributed under the FreeBSD Software License (See accompanying file license.txt)

#ifndef GRID_OPERATIONS_H_
#define GRID_OPERATIONS_H_


#include "grid.h"
#include "box.h"
#include "implicit_functions.h"
#include "linear_interpolator.h"
#include "basic_access_strategy.h"

namespace ls
{
  /**
   * Operations on grids, grid and domain are passed into the run method while algorithm-specific
   * parameters shall be passed to the constructor
   */
  class IGridOperation
  {
  public:
    virtual void run(Grid3D<double>& grid, geometry_utils::Box& domain) const = 0;
    virtual ~IGridOperation() {}
  };

  /**
   * reflects existing grid regarding the plane z=0
   */
  class ReflectGrid : public IGridOperation
  {
  public:
    void run(Grid3D<double>& grid, geometry_utils::Box& domain) const;
  };

  /**
   * Fills in a grid by values generated by implicit function
   */
  class FillInGrid : public IGridOperation
  {
    IImplicitFunctionDPtr m_lsfunc;
  public:

    FillInGrid(const IImplicitFunctionDPtr lsfunc)
    : m_lsfunc(lsfunc)
    {}

    void run(Grid3D<double>& grid, geometry_utils::Box& domain) const;
  };

  /**
   * Reduces the grid resolution
   */
  class CoarsenGrid : public IGridOperation
  {
    size_t m_dim[3];
  public:
    CoarsenGrid(size_t n, size_t m, size_t w)
    {
      m_dim[0] = n;
      m_dim[1] = m;
      m_dim[2] = w;
    }

    void run(Grid3D<double>& grid, geometry_utils::Box& domain) const;
  };

  // methods bodies
  void ReflectGrid::run(Grid3D<double>& grid, geometry_utils::Box& domain) const
  {
    Grid3D<double> twiceBiggerGrid(grid.size(0), grid.size(1), 2 * grid.size(2));
    geometry_utils::Box newDomain(domain.getSizeX(), domain.getSizeY(), 2.0 * domain.getSizeZ());

    size_t center = grid.size(2);
    for (size_t iz = 0; iz < twiceBiggerGrid.size(2); ++iz)
      for (size_t iy = 0; iy < twiceBiggerGrid.size(1); ++iy)
        for (size_t ix = 0; ix < twiceBiggerGrid.size(0); ++ix)
        {
          if (iz >= center)
            twiceBiggerGrid(ix, iy, iz) = grid(ix, iy, twiceBiggerGrid.size(2) - iz - 1);
          else
            twiceBiggerGrid(ix, iy, iz) = grid(ix, iy, iz);
        }

    std::swap(twiceBiggerGrid, grid);
    std::swap(domain, newDomain);
  }

  void FillInGrid::run(Grid3D<double>& grid, geometry_utils::Box& domain) const
  {
    size_t n = grid.size(0), m = grid.size(1), w = grid.size(2);

    double h[] = {domain.getSizeX() / (n - 1.0), domain.getSizeY() / (m - 1.0), domain.getSizeZ() / (w - 1.0)};

    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < m; ++j) {
        for (size_t k = 0; k < w; ++k) {
          geometry_utils::MathVector3D p(i * h[0], j * h[1], k * h[2]);
          p += domain.getLow();
          assert(domain.inside(p));
          grid(i, j, k) = m_lsfunc->compute(p);
        }
      }
    }
  }

  void CoarsenGrid::run(Grid3D<double>& grid, geometry_utils::Box& domain) const
  {
    if (grid.size(0) < m_dim[0] || grid.size(1) < m_dim[1] || grid.size(2) < m_dim[2]) {
      throw std::logic_error("at least one of the input grid's dimensions are smaller than output grid dimensions");
    }

    Grid3D<double> outGrid(m_dim[0], m_dim[1], m_dim[2]);
    ls::LinearInterpolator<double, ls::BasicReadAccessStrategy > li(domain, grid);

    double h[] = {domain.getSizeX() / (m_dim[0] - 1.0), domain.getSizeY() / (m_dim[1] - 1.0), domain.getSizeZ() / (m_dim[2] - 1.0)};

    for (size_t i = 0; i < m_dim[0]; ++i) {
      for (size_t j = 0; j < m_dim[1]; ++j) {
        for (size_t k = 0; k < m_dim[2]; ++k) {
          geometry_utils::MathVector3D p(i * h[0], j * h[1], k * h[2]);
          p += domain.getLow();
          assert(domain.inside(p));
          outGrid(i, j, k) = li.compute(p);
        }
      }
    }

    std::swap(outGrid, grid);
  }

}

#endif /* GRID_OPERATIONS_H_ */