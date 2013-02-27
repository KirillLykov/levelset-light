//  (C) Copyright Kirill Lykov 2012.
//
// Distributed under the FreeBSD Software License (See accompanying file license.txt)

#include <gtest/gtest.h>
#include <random> //c++11
#include "tolerance.h"
#include "grid_operations.h"
#include "implicit_functions.h"
#include "linear_interpolator.h"
#include "basic_access_strategy.h"

using namespace ls;
using namespace geometry_utils;

TEST(GridOperations, ReflectGrid)
{
  size_t n = 8, m = 8, w = 16;
  Grid3D<double> originalGrid(n, m, w);
  Box originalDomain(10.0, 20.0, 30.0);
  double h = 1.0 / (n - 1);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      for (size_t k = 0; k < w; ++k) {
        MathVector3D point( i * h - 0.5, j * h - 0.5, k * h - 0.5);
        originalGrid(i, j, k) = point.getX() + point.getY() * point.getZ() + 6.0 * point.getZ();
      }
    }
  }

  Grid3D<double> reflectedGrid = originalGrid;
  Box reflecteDomain = originalDomain;

  ReflectGrid rg;
  rg.run(reflectedGrid, reflecteDomain);

  EXPECT_EQ(reflectedGrid.size(0), originalGrid.size(0));
  EXPECT_EQ(reflectedGrid.size(1), originalGrid.size(1));
  EXPECT_EQ(reflectedGrid.size(2), 2 * originalGrid.size(2));
  EXPECT_EQ(reflecteDomain.getIthSize(2), 2.0 * originalDomain.getIthSize(2));

  size_t center = originalGrid.size(2);
  for (size_t iz = 0; iz < reflectedGrid.size(2); ++iz)
    for (size_t iy = 0; iy < reflectedGrid.size(1); ++iy)
      for (size_t ix = 0; ix < reflectedGrid.size(0); ++ix)
      {
        if (iz >= center)
          EXPECT_EQ(reflectedGrid(ix, iy, iz), originalGrid(ix, iy, reflectedGrid.size(2) - iz - 1));
        else
          EXPECT_EQ(reflectedGrid(ix, iy, iz), originalGrid(ix, iy, iz));
      }
}

TEST(GridOperations, FillInGrid)
{
  size_t n = 8, m = 8, w = 16;
  IImplicitFunctionDPtr func( new SphereD(MathVector3D(0.0, 0.0, 0.0), 1.0) );
  FillInGrid fill(func);

  Grid3D<double> grid(n, m, w);
  Box domain(10.0, 20.0, 30.0);
  fill.run(grid, domain);

  double h[] = {domain.getSizeX() / (n - 1.0), domain.getSizeY() / (m - 1.0), domain.getSizeZ() / (w - 1.0)};

  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      for (size_t k = 0; k < w; ++k) {
        MathVector3D point(i * h[0], j * h[1], k * h[2]);
        point += domain.getLow();
        EXPECT_EQ(grid(i, j, k), func->compute(point));
      }
    }
  }
}

TEST(GridOperations, CoarsenGrid)
{
  Close_absolut close_at_absolut(10e-12);
  size_t rndPointsCount = 10;

  size_t n = 20, m = 25, w = 30;
  IImplicitFunctionDPtr func( new SphereD(MathVector3D(0.0, 0.0, 0.0), 1.0) );
  FillInGrid fill(func);

  Grid3D<double> grid(n, m, w);
  Box domain(10.0, 20.0, 30.0);
  fill.run(grid, domain);

  CoarsenGrid cg(20, 20, 20);

  Grid3D<double> outGrid = grid;
  cg.run(outGrid, domain);

  ls::LinearInterpolator<double, ls::BasicReadAccessStrategy > liOriginal(domain, grid);
  ls::LinearInterpolator<double, ls::BasicReadAccessStrategy > liOut(domain, outGrid);

  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution1(domain.getLow().getX(), domain.getTop().getX());
  std::uniform_real_distribution<double> distribution2(-domain.getLow().getY(), domain.getTop().getY());
  std::uniform_real_distribution<double> distribution3(-domain.getLow().getZ(), domain.getTop().getZ());

  for (size_t i = 0; i < rndPointsCount; ++i) {
    MathVector3D point(distribution1(generator), distribution2(generator), distribution3(generator));
    EXPECT_TRUE( close_at_absolut(liOriginal.compute(point), liOut.compute(point)) );
  }
}
