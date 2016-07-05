//  (C) Copyright Kirill Lykov 2012.
//
// Distributed under the FreeBSD Software License (See accompanying file license.txt)

#include <gtest/gtest.h>
#include <grid.h>
#include <iostream>
#include <random> //c++11
#include "Grad.h"
#include "basic_access_strategy.h"
typedef ls::Grad<double, ls::BasicReadAccessStrategy > _BasicGrad;
typedef ls::Grad<double, ls::PeriodicReadAS > _PeriodicGrad;

using namespace ls;
using namespace geometry_utils;

TEST(GradTest, sphere)
{
  Box3D box(16.0);
  MathVector3D low = box.getLow();

  size_t n = 32, m = 32, w = 32;
  double h[] = {box.getSizeX()/n, box.getSizeY()/m, box.getSizeZ()/w};
  Grid3D<double> grid(n, m, w, box);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      for (size_t k = 0; k < w; ++k) {
        MathVector3D p = MathVector3D(i * h[0], j * h[1] , k * h[2]) + low;
        grid(i, j, k) = p.getLength() - 6.0; // sphere levelset
      }
    }
  }

  _BasicGrad grad(grid);

  double ferror = -1e6;
  double cerror = -1e6;
  //check that interpolation at known points are as defined in the grid
  for (size_t i = 2; i < n-2; ++i) {
    for (size_t j = 2; j < m-2; ++j) {
      for (size_t k = 2; k < w-2; ++k) {
        //skip discontinous points
        if ((i >= n/2 - 3 && i <= n/2 + 3) &&
            (j >= m/2 - 3 && j <= m/2 + 3) &&
            (k >= w/2 - 3 && k <= w/2 + 3))
              continue;

        MathVector3D p = MathVector3D(i * h[0], j * h[1] , k * h[2]) + low;
        MathVector3D fgrad = grad.compute_forward(p);
        MathVector3D cgrad = grad.compute_central(p);
        MathVector3D groundTruth = p;
        groundTruth.normalize();

        ferror = std::max(ferror, (fgrad - groundTruth).getLength());
        cerror = std::max(cerror, (cgrad - groundTruth).getLength());
      }
    }
  }
  EXPECT_TRUE( ferror < 0.5 );
  EXPECT_TRUE( cerror < 0.5 );
}

TEST(GradTest, bback_grad)
{
  using namespace ls;
  size_t n = 32, m = 32, w = 32;

  IImplicitFunctionDPtr func( new SphereD(MathVector3D(0.0, 0.0, 0.0), 2.0) );
  FillInGrid fill(func);

  Grid3D<double> grid(n, m, w, Box3D(5.0, 5.0, 5.0));
  fill.run(grid);

  const double tolerance = 5.0/n * sqrt(3);

  _BasicGrad collision(grid);
  //
  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(-2.49, 2.49);

  int i = 0, j = 0;

  while (i < 1000) {
    MathVector3D pos(distribution(generator), distribution(generator), distribution(generator));
    const MathVector3D grad_exact = pos / pos.getLength();
    MathVector3D grad1 = collision.compute_forward(pos);
    EXPECT_LT( fabs(grad1.getLength() - 1.0), 0.2);
    double diff1 = (grad1 - grad_exact).getLength();
    EXPECT_TRUE(diff1 < 2*tolerance);
    MathVector3D grad2 = collision.compute_precise(pos);
    EXPECT_LT( fabs(grad2.getLength() - 1.0), 0.2);
    double diff2 = (grad2 - grad_exact).getLength();
    EXPECT_TRUE(diff2 < tolerance);

    if (diff2 > diff1)
      ++j;
    ++i;
  }

  EXPECT_LT( j/1000.0, 0.5 ); // usually precise scheme is better
}

TEST(GradTest, grad_discont)
{
  using namespace ls;
  const size_t n = 128, m = 128, w = 128;

  double r = 5.0;
  IImplicitFunctionDPtr func( new RecPipe<double>(r) );
  FillInGrid fill(func);

  Grid3D<double> grid(n, m, w, Box3D(15.0, 15.0, 15.0));
  fill.run(grid);

  const double tolerance = 15.0/n * sqrt(3);

  _PeriodicGrad collision(grid);

  const MathVector3D top(r, 0.0, 0.0);
  double dist = collision.compute(top);
  EXPECT_LT(dist, tolerance);

  // now add pertrubations
  {
    MathVector3D grad_exact(1.0, 0.0, 0.0); // also strictly speaking it does not exist
    MathVector3D p = top + MathVector3D(1.0/n, 1.0/n, 1.0/n);
    MathVector3D grad1 = collision.compute_forward(p);
    double diff1 = (grad1 - grad_exact).getLength();
    EXPECT_LT(diff1, 1e-12);
    MathVector3D grad2 = collision.compute_precise(p);
    double diff2 = (grad2 - grad_exact).getLength();
    EXPECT_LT(diff2, 1e-12);
  }

  {
    MathVector3D grad_exact(1.0, 1.0, 0.0);
    grad_exact.normalize();
    MathVector3D p = top + MathVector3D(-1e-1, 1e-1, 0.0);
    MathVector3D grad1 = collision.compute_precise(p);
    grad1.normalize();
    EXPECT_LT((grad1 - grad_exact).getLength(), 1e-12);

    MathVector3D grad2 = collision.compute_forward(p);
    EXPECT_LT((grad2 - MathVector3D(1.0, 1.0, 0.0)).getLength(), 1e-12);
    MathVector3D grad3 = collision.compute_backward(p);
    EXPECT_LT((grad3 - MathVector3D(1.0, 0.0, 0.0)).getLength(), 1e-12);

    MathVector3D vel(0.0, 1.0, 0.0);
    MathVector3D grad4 = collision.compute_biased(p, vel);
    EXPECT_LT((grad4 - MathVector3D(1.0, 1.0, 0.0)).getLength(), 1e-12);

    vel = MathVector3D(0.0, -1.0, 0.0);
    MathVector3D grad5 = collision.compute_biased(p, vel);
    EXPECT_LT((grad5 - MathVector3D(1.0, 0.0, 0.0)).getLength(), 1e-12);
  }
}

