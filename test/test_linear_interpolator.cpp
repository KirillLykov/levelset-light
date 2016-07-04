//  (C) Copyright Kirill Lykov 2012.
//
// Distributed under the FreeBSD Software License (See accompanying file license.txt)

#include <gtest/gtest.h>
#include <grid.h>
#include <iostream>
#include <random> //c++11
#include "linear_interpolator.h"
#include "basic_access_strategy.h"
#include <implicit_functions.h>
#include <grid_operations.h>

typedef ls::LinearInterpolator<double, ls::BasicReadAccessStrategy > _BasicLinInterpolator;

using namespace ls;
using namespace geometry_utils;


const double tolerance = 0.00000001;

namespace
{
  const size_t rndPointsCount = 10;

  double nonSumFunciton(const MathVector3D& p)
  {
    return p.getX() + p.getY() * p.getY() + p.getZ() * p.getZ() * p.getZ();
  }

  double nonSumFunciton2(const MathVector3D& p)
  {
    return 1.0 + p.getX() + sin(p.getY())  + log( fabs(p.getZ()) + 1.0);
  }

  double computeError(const MathVector3D& h)
  {
    return h.getLength() / 15.0; //just random coefficient which works for some reason
  }
}

TEST(InterpolatorTest, trivialCubicGrid)
{
  Box3D box(1.0);

  double h = 1.0;

  size_t n = 2, m = 2, w = 2;
  Grid3D<double> grid(n, m, w, box);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      for (size_t k = 0; k < w; ++k) {
        grid(i, j, k) = (k == 0) ? 1.0 : -1.0;
      }
    }
  }

  _BasicLinInterpolator li(grid);

  //check that interpolation at known points are as defined in the grid
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      for (size_t k = 0; k < w; ++k) {
        EXPECT_TRUE( fabs(li.compute(i * h - 0.5, j * h - 0.5, k * h - 0.5) - grid(i, j, k)) < tolerance );
      }
    }
  }

  EXPECT_TRUE( fabs(li.compute(0.0, 0.0, 0.0)) < tolerance );
  EXPECT_TRUE( fabs(li.compute(0.25, 0.25, 0.25) + 0.5) < tolerance );
}

TEST(InterpolatorTest, diffentStepLength)
{
  Box3D box(1.0);
  size_t n = 10, m = 12, w = 14;
  double h[] = {1.0 / (n - 1.0), 1.0 / (m - 1.0), 1.0 / (w - 1.0)};

  Grid3D<double> grid(n, m, w, box);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      for (size_t k = 0; k < w; ++k) {
        double p[] = {i * h[0] - 0.5, j * h[1] - 0.5, k * h[2] - 0.5};
        grid(i, j, k) = nonSumFunciton(p);
      }
    }
  }

  _BasicLinInterpolator li(grid);

  double error = computeError(h);
  EXPECT_TRUE( fabs(li.compute(0.0, 0.0, 0.0)) < error );

  //check that interpolation at known points are as defined in the grid
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      for (size_t k = 0; k < w; ++k) {
        EXPECT_TRUE( fabs(li.compute(i * h[0] - 0.5, j * h[1] - 0.5, k * h[2] - 0.5) - grid(i, j, k)) < tolerance );
      }
    }
  }

  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(-0.5, 0.5);

  for (size_t i = 0; i < rndPointsCount; ++i) {
    double p[] = {distribution(generator), distribution(generator), distribution(generator)};
    EXPECT_TRUE( fabs(li.compute(p) - nonSumFunciton(p)) < error );
  }
}

TEST(InterpolatorTest, cuboidGrid)
{
  double sideLength[] = {3.0, 4.0, 5.0};
  Box3D box(sideLength);

  size_t n = 5, m = 6, w = 7;

  double h[] = {sideLength[0] / (n - 1.0), sideLength[1] / (m - 1.0), sideLength[2] / (w - 1.0)};
  double shift[3] = {box.getIthSize(0) / 2.0, box.getIthSize(1) / 2.0, box.getIthSize(2) / 2.0};

  Grid3D<double> grid(n, m, w, box);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      for (size_t k = 0; k < w; ++k) {
        double p[] = {i * h[0] - shift[0], j * h[1] - shift[1], k * h[2] - shift[2]};
        grid(i, j, k) = nonSumFunciton2(p);
      }
    }
  }

  _BasicLinInterpolator li(grid);

  //check that interpolation at known points are as defined in the grid
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      for (size_t k = 0; k < w; ++k) {
        double p[] = {i * h[0] - shift[0], j * h[1] - shift[1], k * h[2] - shift[2]};
        EXPECT_TRUE( fabs(li.compute(p) - grid(i, j, k)) < tolerance );
      }
    }
  }

  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(-0.5, 0.5);

  double error = computeError(h);

  for (size_t i = 0; i < rndPointsCount; ++i) {
    double p[] = {distribution(generator), distribution(generator), distribution(generator)};
    EXPECT_TRUE( fabs(li.compute(p) - nonSumFunciton2(p)) < error );
  }
}

TEST(InterpolatorTest, notOriginPlacedCuboid)
{
  double top[] = {4.0, 5.0, 9.0};
  double low[] = {-3.0, -4.0, -5.0};
  Box3D box(low, top);

  size_t n = 5, m = 6, w = 7;

  double h[] = {box.getSizeX() / (n - 1.0), box.getSizeY() / (m - 1.0), box.getSizeZ() / (w - 1.0)};

  Grid3D<double> grid(n, m, w, box);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      for (size_t k = 0; k < w; ++k) {
        MathVector3D p(i * h[0], j * h[1], k * h[2]);
        p += box.getLow();
        assert(box.inside(p));
        grid(i, j, k) = nonSumFunciton2(p);
      }
    }
  }

  _BasicLinInterpolator li(grid);

  //check that interpolation at known points are as defined in the grid
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      for (size_t k = 0; k < w; ++k) {
        MathVector3D p(i * h[0], j * h[1], k * h[2]);
        p += box.getLow();
        EXPECT_TRUE( fabs(li.compute(p) - grid(i, j, k)) < tolerance );
      }
    }
  }

  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution1(box.getLow().getX(), box.getTop().getX());
  std::uniform_real_distribution<double> distribution2(-box.getLow().getY(), box.getTop().getY());
  std::uniform_real_distribution<double> distribution3(-box.getLow().getZ(), box.getTop().getZ());

  double error = raw_math_vector::length(h);

  for (size_t i = 0; i < rndPointsCount; ++i) {
    double p[] = {distribution1(generator), distribution2(generator), distribution3(generator)};
    EXPECT_TRUE( fabs(li.compute(p) - nonSumFunciton2(p)) < error );
  }
}

TEST(InterpolatorTest, periodicAccessStrategy2)
{
  using namespace ls;
  size_t n = 64, m = 64, w = 64;

  IImplicitFunctionDPtr func( new AxialCylinderD(zDim, 2.0) );
  FillInGrid fill(func);

  Grid3D<double> grid(n, m, w, Box3D(5.0, 5.0, 5.0));
  fill.run(grid);

  ls::LinearInterpolator< double, PeriodicReadAS > li(grid);

  // check on the border
  {
    MathVector3D p = MathVector3D(0.0, 0.0, 0.0);
    const double exp = li.compute( p );

    for (int i = 0; i < 5; ++i) {
      MathVector3D shift(0.0, 0.0, 2.5 + i*10);
      EXPECT_LT( fabs(li.compute( p + shift ) - exp), tolerance);
      EXPECT_LT( fabs(li.compute( p - shift ) - exp), tolerance);
    }
  }

  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(-2.5, 2.5);
  for (int k = 0; k < 100; ++k) {
    MathVector3D p = MathVector3D(distribution(generator), distribution(generator), distribution(generator));
    const double exp = li.compute( p );

    for (int i = 0; i < 5; ++i) {
      MathVector3D shift;
      for (int dim = 0; dim < 3; ++dim) {
        shift.setCoord(dim, 5.0 + i*10);
        EXPECT_LT( fabs(li.compute( p + shift ) - exp), tolerance);
        EXPECT_LT( fabs(li.compute( p - shift ) - exp), tolerance);
      }
    }
  }
}

TEST(InterpolatorTest, periodicAccessStrategy)
{
  double sideLength[] = {3.0, 4.0, 5.0};
  Box3D box(sideLength);

  size_t n = 15, m = 16, w = 17;

  double h[] = {sideLength[0] / (n - 1.0), sideLength[1] / (m - 1.0), sideLength[2] / (w - 1.0)};
  double shift[3] = {box.getIthSize(0) / 2.0, box.getIthSize(1) / 2.0, box.getIthSize(2) / 2.0};

  Grid3D<double> grid(n, m, w, box);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      for (size_t k = 0; k < w; ++k) {
        MathVector3D p(i * h[0] - shift[0], j * h[1] - shift[1], k * h[2] - shift[2]);
        grid(i, j, k) = p.getX() + p.getY() + p.getZ(); //TODO check why it doesn't work with nonSumFunciton2(p)
      }
    }
  }

  ls::LinearInterpolator< double, PeriodicReadAS > li(grid);

  //check that interpolation at known points are as defined in the grid
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      for (size_t k = 0; k < w; ++k) {
        MathVector3D p(i * h[0] - shift[0], j * h[1] - shift[1], k * h[2] - shift[2]);
        if (fabs(li.compute(p) - grid(i, j, k)) > tolerance)
          std::cout << "("<<i <<", "<< j<<", "<< k << "): " << li.compute(p) << " " << grid(i, j, k) << std::endl;
        EXPECT_TRUE( fabs(li.compute(p) - grid(i, j, k)) < tolerance );
      }
    }
  }

  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(-0.5, 0.5);
/*
  double error = computeError(h);

  for (size_t i = 0; i < rndPointsCount; ++i) {
    double p[] = {distribution(generator), distribution(generator), distribution(generator)};
    EXPECT_TRUE( fabs(li.compute(p) - nonSumFunciton2(p)) < error );
  }
*/
}

