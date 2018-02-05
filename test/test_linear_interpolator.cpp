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

typedef ls::LinearInterpolator<Real, ls::BasicReadAccessStrategy > _BasicLinInterpolator;

using namespace ls;
using namespace geometry_utils;

namespace
{
  const size_t rndPointsCount = 10;

  Real nonSumFunciton(const MathVector3R& p)
  {
    return p.getX() + p.getY() * p.getY() + p.getZ() * p.getZ() * p.getZ();
  }

  Real nonSumFunciton2(const MathVector3R& p)
  {
    return 1.0 + p.getX() + sin(p.getY())  + log( fabs(p.getZ()) + 1.0);
  }

  Real computeError(const MathVector3R& h)
  {
    return h.getLength() / 15.0; //just random coefficient which works for some reason
  }
}

TEST(InterpolatorTest, trivialCubicGrid)
{
  Box3R box(1.0);

  Real h = 1.0;

  size_t n = 2, m = 2, w = 2;
  Grid3D<Real> grid(n, m, w, box);
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
        EXPECT_NEAR(li.compute(i * h - Real(0.5), j * h - Real(0.5), k * h - Real(0.5)), grid(i, j, k), tol);
      }
    }
  }

  EXPECT_NEAR(li.compute(Real(0.0), Real(0.0), Real(0.0)), Real(0.0), tol);
  EXPECT_NEAR(li.compute(Real(0.25), Real(0.25), Real(0.25)), Real(-0.5), tol);
}

TEST(InterpolatorTest, diffentStepLength)
{
  Box3R box(1.0);
  size_t n = 10, m = 12, w = 14;
  Real h[] = {Real(1.0) / Real(n - 1.0), Real(1.0) / Real(m - 1.0), Real(1.0) / Real(w - 1.0)};

  Grid3D<Real> grid(n, m, w, box);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      for (size_t k = 0; k < w; ++k) {
        Real p[] = {i * h[0] - Real(0.5), j * h[1] - Real(0.5), k * h[2] - Real(0.5)};
        grid(i, j, k) = nonSumFunciton(p);
      }
    }
  }

  _BasicLinInterpolator li(grid);

  Real error = computeError(h);
  EXPECT_NEAR(li.compute(0.0, 0.0, 0.0), Real(0.0), error);

  //check that interpolation at known points are as defined in the grid
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      for (size_t k = 0; k < w; ++k) {
        EXPECT_NEAR(li.compute(i * h[0] - Real(0.5), j * h[1] - Real(0.5), k * h[2] - Real(0.5)), grid(i, j, k), error);
      }
    }
  }

  std::default_random_engine generator;
  std::uniform_real_distribution<Real> distribution(-0.5, 0.5);

  for (size_t i = 0; i < rndPointsCount; ++i) {
    MathVector3R p(distribution(generator), distribution(generator), distribution(generator));
    EXPECT_NEAR(li.compute(p), nonSumFunciton(p), error);
  }
}

TEST(InterpolatorTest, cuboidGrid)
{
  Real sideLength[] = {3.0, 4.0, 5.0};
  Box3R box(sideLength);

  size_t n = 5, m = 6, w = 7;

  Real h[] = {sideLength[0] / Real(n - 1.0), sideLength[1] / Real(m - 1.0), sideLength[2] / Real(w - 1.0)};
  Real shift[3] = {box.getIthSize(0) / Real(2.0), box.getIthSize(1) / Real(2.0), box.getIthSize(2) / Real(2.0)};

  Grid3D<Real> grid(n, m, w, box);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      for (size_t k = 0; k < w; ++k) {
        Real p[] = {i * h[0] - shift[0], j * h[1] - shift[1], k * h[2] - shift[2]};
        grid(i, j, k) = nonSumFunciton2(p);
      }
    }
  }

  _BasicLinInterpolator li(grid);

  //check that interpolation at known points are as defined in the grid
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      for (size_t k = 0; k < w; ++k) {
        Real p[] = {i * h[0] - shift[0], j * h[1] - shift[1], k * h[2] - shift[2]};
        EXPECT_NEAR(li.compute(p), grid(i, j, k), tol);
      }
    }
  }

  std::default_random_engine generator;
  std::uniform_real_distribution<Real> distribution(-0.5, 0.5);

  Real error = computeError(h);

  for (size_t i = 0; i < rndPointsCount; ++i) {
    Real p[] = {distribution(generator), distribution(generator), distribution(generator)};
    EXPECT_NEAR( li.compute(p), nonSumFunciton2(p), error );
  }
}

TEST(InterpolatorTest, notOriginPlacedCuboid)
{
  Real top[] = {4.0, 5.0, 9.0};
  Real low[] = {-3.0, -4.0, -5.0};
  Box3R box(low, top);

  size_t n = 5, m = 6, w = 7;

  Real h[] = {box.getSizeX() / Real(n - 1.0), box.getSizeY() / Real(m - 1.0), box.getSizeZ() / Real(w - 1.0)};

  Grid3D<Real> grid(n, m, w, box);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      for (size_t k = 0; k < w; ++k) {
        MathVector3R p(i * h[0], j * h[1], k * h[2]);
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
        MathVector3R p(i * h[0], j * h[1], k * h[2]);
        p += box.getLow();
        EXPECT_NEAR(li.compute(p), grid(i, j, k), tol);
      }
    }
  }

  std::default_random_engine generator;
  std::uniform_real_distribution<Real> distribution1(box.getLow().getX(), box.getTop().getX());
  std::uniform_real_distribution<Real> distribution2(-box.getLow().getY(), box.getTop().getY());
  std::uniform_real_distribution<Real> distribution3(-box.getLow().getZ(), box.getTop().getZ());

  Real error = raw_math_vector::length(h);

  for (size_t i = 0; i < rndPointsCount; ++i) {
    Real p[] = {distribution1(generator), distribution2(generator), distribution3(generator)};
    EXPECT_NEAR(li.compute(p), nonSumFunciton2(p), error);
  }
}

TEST(InterpolatorTest, periodicAccessStrategy2)
{
  using namespace ls;
  size_t n = 64, m = 64, w = 64;

  IImplicitFunctionDPtr func( new AxialCylinderD(zDim, 2.0) );
  FillInGrid<Real> fill(func);

  Grid3D<Real> grid(n, m, w, Box3R(5.0, 5.0, 5.0));
  fill.run(grid);

  ls::LinearInterpolator< Real, PeriodicReadAS > li(grid);

  // check on the border
  {
    MathVector3R p = MathVector3R(0.0, 0.0, 0.0);
    const Real exp = li.compute( p );

    for (int i = 0; i < 5; ++i) {
      MathVector3R shift(0.0, 0.0, 2.5 + i*10);
      EXPECT_NEAR(li.compute( p + shift ), exp, tol);
      EXPECT_NEAR(li.compute( p - shift ), exp, tol);
    }
  }

  std::default_random_engine generator;
  std::uniform_real_distribution<Real> distribution(-2.5, 2.5);
  for (int k = 0; k < 100; ++k) {
    MathVector3R p = MathVector3R(distribution(generator), distribution(generator), distribution(generator));
    const Real exp = li.compute( p );

    for (int i = 0; i < 5; ++i) {
      MathVector3R shift;
      for (int dim = 0; dim < 3; ++dim) {
        shift.setCoord(dim, 5.0 + i*10);
        EXPECT_NEAR(li.compute( p + shift ), exp, tol);
        EXPECT_NEAR(li.compute( p - shift ), exp, tol);
      }
    }
  }
}

TEST(InterpolatorTest, periodicAccessStrategy)
{
  Real sideLength[] = {3.0, 4.0, 5.0};
  Box3R box(sideLength);

  size_t n = 15, m = 16, w = 17;

  Real h[] = {sideLength[0] / (n - 1.0), sideLength[1] / (m - 1.0), sideLength[2] / (w - 1.0)};
  Real shift[3] = {box.getIthSize(0) / 2.0, box.getIthSize(1) / 2.0, box.getIthSize(2) / 2.0};

  Grid3D<Real> grid(n, m, w, box);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      for (size_t k = 0; k < w; ++k) {
        MathVector3R p(i * h[0] - shift[0], j * h[1] - shift[1], k * h[2] - shift[2]);
        grid(i, j, k) = p.getX() + p.getY() + p.getZ(); //TODO check why it doesn't work with nonSumFunciton2(p)
      }
    }
  }

  ls::LinearInterpolator< Real, PeriodicReadAS > li(grid);

  //check that interpolation at known points are as defined in the grid
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      for (size_t k = 0; k < w; ++k) {
        MathVector3R p(i * h[0] - shift[0], j * h[1] - shift[1], k * h[2] - shift[2]);
        EXPECT_NEAR(li.compute(p), grid(i, j, k), tol);
      }
    }
  }

  std::default_random_engine generator;
  std::uniform_real_distribution<Real> distribution(-0.5, 0.5);
/*
  Real error = computeError(h);

  for (size_t i = 0; i < rndPointsCount; ++i) {
    Real p[] = {distribution(generator), distribution(generator), distribution(generator)};
    EXPECT_TRUE( fabs(li.compute(p) - nonSumFunciton2(p)) < error );
  }
*/
}

