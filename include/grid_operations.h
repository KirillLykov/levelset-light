//  (C) Copyright Kirill Lykov 2012.
//
// Distributed under the FreeBSD Software License (See accompanying file license.txt)

#ifndef GRID_OPERATIONS_H_
#define GRID_OPERATIONS_H_


#include "grid.h"
#include "box.h"
#include "implicit_functions.h"

namespace ls
{
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

  class FillInGrid : public IGridOperation
  {
    IImplicitFunctionDPtr m_lsfunc;
  public:

    FillInGrid(const IImplicitFunctionDPtr lsfunc)
    : m_lsfunc(lsfunc)
    {}

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

}

#endif /* GRID_OPERATIONS_H_ */
