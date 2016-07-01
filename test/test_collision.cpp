//  (C) Copyright Kirill Lykov 2016.
//
// Distributed under the FreeBSD Software License (See accompanying file license.txt)

#include <gtest/gtest.h>
#include <math_vector.h>
#include <implicit_functions.h>
#include <grid.h>
#include <grid_operations.h>
#include <collision.h>
#include <box.h>

#include <random>
#include <iostream>

namespace {
  template<typename T>
  class Square : public IImplicitFunction<T>
  {
    T m_d;
  public:
    Square(T d) : m_d(d) {}
    T compute(const MathVector<T, 3>& point) const
    {
      return fabs(point.getX()) + fabs(point.getY()) + fabs(point.getZ()) - m_d;
    }
  };

  template<typename T>
  class TwoPlanceInY : public IImplicitFunction<T>
  {
    T m_d;
  public:
    TwoPlanceInY(T d) : m_d(d) {}
    T compute(const MathVector<T, 3>& point) const
    {
      return fabs(point.getY()) - m_d;
    }
  };
}

TEST(CollisionTest, bback_grad)
{
  using namespace ls;
  size_t n = 32, m = 32, w = 32;

  IImplicitFunctionDPtr func( new SphereD(MathVector3D(0.0, 0.0, 0.0), 2.0) );
  FillInGrid fill(func);

  Grid3D<double> grid(n, m, w, Box3D(5.0, 5.0, 5.0));
  fill.run(grid);

  const double tolerance = 5.0/n * sqrt(3);

  Collision collision(&grid);
  //
  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(-2.49, 2.49);

  int i = 0, j = 0;

  while (i < 1000) {
    MathVector3D pos(distribution(generator), distribution(generator), distribution(generator));
    const MathVector3D grad_exact = pos / pos.getLength();
    MathVector3D grad1 = collision.computeFGrad(pos);
    EXPECT_LT( fabs(grad1.getLength() - 1.0), 0.2);
    double diff1 = (grad1 - grad_exact).getLength();
    EXPECT_TRUE(diff1 < 2*tolerance);
    MathVector3D grad2 = collision.computePreciseGrad(pos);
    EXPECT_LT( fabs(grad2.getLength() - 1.0), 0.2);
    double diff2 = (grad2 - grad_exact).getLength();
    EXPECT_TRUE(diff2 < tolerance);

    if (diff2 > diff1)
      ++j;
    ++i;
  }

  EXPECT_LT( j/1000.0, 0.5 ); // usually precise scheme is better
}

TEST(CollisionTest, bback_grad_discont)
{
  using namespace ls;
  const size_t n = 128, m = 128, w = 128;

  double r = 5.0;
  IImplicitFunctionDPtr func( new Square<double>(r) );
  FillInGrid fill(func);

  Grid3D<double> grid(n, m, w, Box3D(15.0, 15.0, 15.0));
  fill.run(grid);

  const double tolerance = 15.0/n * sqrt(3);

  Collision collision(&grid);

  const MathVector3D top(r, 0.0, 0.0);
  double dist = collision.computeSDF(top);
  EXPECT_LT(dist, tolerance);

  // now add pertrubations
  {
    MathVector3D grad_exact(1.0, 0.0, 0.0); // also strictly speaking it does not exist
    MathVector3D p = top + MathVector3D(1.0/n, 1.0/n, 1.0/n);
    MathVector3D grad1 = collision.computeFGrad(p);
    double diff1 = (grad1 - grad_exact).getLength();
    if (diff1 > 1e-12)
      std::cout << "AAA";
    EXPECT_LT(diff1, 1e-12);
    MathVector3D grad2 = collision.computePreciseGrad(p);
    double diff2 = (grad2 - grad_exact).getLength();
    EXPECT_LT(diff2, 1e-12);
  }

  {
    MathVector3D grad_exact(1.0, 1.0, 1.0);
    grad_exact.normalize();
    MathVector3D p = top + MathVector3D(0, 10.0/n, 10.0/n);
    MathVector3D grad1 = collision.computeGrad(p); // since grad will be spoiled
    double diff1 = (grad1 - grad_exact).getLength();
    if (diff1 > 1e-12)
      std::cout << "AAA";
    EXPECT_LT(diff1, 1e-12);
  }
}

TEST(CollisionTest, bback_cylinder)
{
  using namespace ls;
  size_t n = 64, m = 64, w = 64;
  const double tolerance = 1.0/n;
  IImplicitFunctionDPtr func( new AxialCylinderD(zDim, 1.0) );
  FillInGrid fill(func);

  Grid3D<double> grid(n, m, w, Box3D(10.0, 10.0, 10.0));
  fill.run(grid);

  Collision collision(&grid);

  MathVector3D pos(0.0, 1.25, 0.0);
  MathVector3D vel(0.0, 1.0, 0.0);
  double dt = 1.0;
  double currsdf = collision.computeSDF(pos);
  collision.bounceBack(currsdf, dt, pos, vel);

  EXPECT_TRUE( (pos - MathVector3D(0.0, 0.75, 0.0)).getLength() < tolerance );
}

class BBCylinderTest : public testing::Test {
 public:
  virtual void SetUp()
  {
    using namespace ls;
    IImplicitFunctionDPtr func( new AxialCylinderD(zDim, 6) );
    FillInGrid fill(func);

    grid.reset( new Grid3D<double>(n, m, w, Box3D(16.0, 16.0, 16.0)) );
    fill.run(*grid);

    collision.reset( new Collision(grid.get()) );
  }
  virtual void TearDown() {
    grid = nullptr;
    collision = nullptr;
  }

  std::shared_ptr< Grid3D<double> > grid;
  std::shared_ptr< Collision > collision;
  const size_t n = 32, m = 32, w = 32;
  const double tolerance = 1.0/n;
};

TEST_F(BBCylinderTest, bback_cylinder_pbc1)
{
  MathVector3D pos(-4.05265, -4.42603, 8.00028);
  MathVector3D vel(-0.759647, -0.528299, 0.615282);
  double dt = 1.0;

  collision->applyPBC(pos);
  double currsdf = collision->computeSDF(pos);
  collision->bounceBack(currsdf, dt, pos, vel);

  MathVector3D diff = pos - MathVector3D(-4.0492098, -4.4236375, 16-8.0025063);
  collision->applyPBC(diff);
  EXPECT_TRUE( diff.getLength() < tolerance );
}

TEST_F(BBCylinderTest, bback_cylinder_pbc2)
{
  double dt = 1.0;
  MathVector3D pos(0.0, 6.25, 8.1);
  MathVector3D vel(0.0, 1.0, 0.0);
  collision->applyPBC(pos);
  double currsdf = collision->computeSDF(pos);
  collision->bounceBack(currsdf, dt, pos, vel);
  MathVector3D diff = pos - MathVector3D(0.0, 5.7473117, 8.1);
  collision->applyPBC(diff);
  EXPECT_TRUE( diff.getLength() < tolerance );
}

TEST_F(BBCylinderTest, bback_cylinder_pbc3)
{
  // too close to the border so it is a bit inside
  double dt = 1e-3;
  MathVector3D pos(3.85635, 4.59372, 4.82354);
  MathVector3D posOrig = pos;
  MathVector3D vel(-2.52048, 2.13566, -1.08796);
  double  currsdf = collision->computeSDF(pos);
  collision->bounceBack(currsdf, dt, pos, vel);
  double length = (pos - posOrig).getLength();
  double newcurrsdf = collision->computeSDF(pos);
  EXPECT_TRUE( length < tolerance && newcurrsdf < -2.0*1e-4);
}

TEST_F(BBCylinderTest, bback_cylinder_pbc4)
{
  double dt = 1e-3;
  // the point is inside but too close to the interface
  MathVector3D pos(1.192889499240467, 5.876849556783832, 6.565664278186738);
  MathVector3D posOrig = pos;
  MathVector3D vel(-0.381401772476884, 0.1058015442625392, -0.2673344284414041);
  double currsdf = collision->computeSDF(pos);
  collision->bounceBack(currsdf, dt, pos, vel);
  double length = (pos - posOrig).getLength();
  double newcurrsdf = collision->computeSDF(pos);
  EXPECT_TRUE( length < tolerance && newcurrsdf < -2.0*1e-4);
}

TEST_F(BBCylinderTest, bback_cylinder_pbc5)
{
  double dt = 1e-3;
  MathVector3D pos(3.5603906528217, 4.823928346832147, 4.021767564215939);
  MathVector3D posOrig = pos;
  MathVector3D vel(1.968563824141501, -0.002774591556346406, -1.484421351051873);
  double currsdf = collision->computeSDF(pos);
  collision->bounceBack(currsdf, dt, pos, vel);
  double length = (pos - posOrig).getLength();
  double newcurrsdf = collision->computeSDF(pos);
  EXPECT_TRUE( length < 1e-3 && newcurrsdf < 0 && fabs(newcurrsdf) < 1e-3);
}

TEST_F(BBCylinderTest, bback_cylinder_pbc6)
{
  double dt = 1e-3;
  MathVector3D pos(-4.5745, 3.87941, 0.801322);
  MathVector3D vel(1.39311, 1.68965, -0.0923976);
  double currsdf = collision->computeSDF(pos);
  collision->bounceBack(currsdf, dt, pos, vel);
  MathVector3D diff = pos - MathVector3D(-4.5754562441755819, 3.8773540850173442, 0.8014143976);
  collision->applyPBC(diff);
  EXPECT_TRUE( diff.getLength() < tolerance );
}

TEST(CollisionTest, bback_two_planes)
{
  using namespace ls;
  const double tolerance = 1e-4;
  size_t n = 32, m = 32, w = 32;
  IImplicitFunctionDPtr func( new TwoPlanceInY<double>(1.0) );
  FillInGrid fill(func);

  Grid3D<double> grid(n, m, w, Box3D(10.0, 10.0, 10.0));
  fill.run(grid);

  Collision collision(&grid);

  // perpendicular to plane
  double dt = 1.0;
  MathVector3D pos(0.0, 1.25, 0.0);
  double currsdf = collision.computeSDF(pos);
  MathVector3D vel(0.0, 1.0, 0.0);
  collision.bounceBack(currsdf, dt, pos, vel);

  EXPECT_TRUE( (pos - MathVector3D(0.0, 0.75, 0.0)).getLength() < tolerance );

  // angle to plane
  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(0.0,1.0);

  int i = 0;
  while (i < 100) {
    pos = MathVector3D(distribution(generator) + 0.5, distribution(generator) + 1.0, 0.0); //;(1.0, 1.25, 0.0);
    currsdf = collision.computeSDF(pos);

    double angle = distribution(generator) * M_PI/2.0;
    vel = MathVector3D(cos(angle), sin(angle), 0.0);

    const MathVector3D posOld = pos - dt * vel;
    MathVector3D exp = pos - posOld;
    if (collision.computeSDF(posOld) >= 0.0)
      continue;

    exp *= (2.0*(1.0-posOld.getY())/exp.getY() - 1.0);
    exp += posOld;

    collision.bounceBack(currsdf, dt, pos, vel);
    EXPECT_TRUE( (pos - exp).getLength() < tolerance );
    ++i;
  }
}

TEST(CollisionTest, bback_multiplereflections)
{
  using namespace ls;
  const double tolerance = 1e-4;
  size_t n = 64, m = 64, w = 64;
  IImplicitFunctionDPtr func( new TwoPlanceInY<double>(1.0) );
  FillInGrid fill(func);

  Grid3D<double> grid(n, m, w, Box3D(20.0, 20.0, 20.0));
  fill.run(grid);

  Collision collision(&grid);

  // perpendicular to plane, two reflections
  double dt = 1.0;
  MathVector3D pos(0.0, 3.5, 0.0);
  double currsdf = collision.computeSDF(pos);
  MathVector3D vel(0.0, 4.0, 0.0);
  collision.bounceBack(currsdf, dt, pos, vel);

  EXPECT_TRUE( (pos - MathVector3D(0.0, -0.5, 0.0)).getLength() < tolerance );

  // many reflections
  pos = MathVector3D(0.0, 7.998, 0.0);
  currsdf = collision.computeSDF(pos);
  vel = MathVector3D(0.0, 8.0, 0.0);
  collision.bounceBack(currsdf, dt, pos, vel);

  EXPECT_TRUE( (pos - MathVector3D(0.0, -0.002, 0.0)).getLength() < tolerance );
}

TEST(CollisionTest, bback_square)
{
  using namespace ls;
  const double tolerance = 1e-4;
  size_t n = 64, m = 64, w = 64;
  IImplicitFunctionDPtr func( new Square<double>(5.0) );
  FillInGrid fill(func);

  Grid3D<double> grid(n, m, w, Box3D(15.0, 15.0, 15.0));
  fill.run(grid);

  BasicSerializerHDF5 writer(grid, "test-aux/square", "data");
  writer.run();

  Collision collision(&grid);
}

/*
template<typename T>
class ThinLayer : public IImplicitFunction<T>
{
  T m_d;
public:
  ThinLayer(T d) : m_d(d) {}
  T compute(const MathVector<T, 3>& point) const
  {
    return -(fabs(point.getY()) - m_d);
  }
};

TEST(CollisionTest, bback_gothrough)
{
  using namespace ls;
  const double tolerance = 1e-4;
  size_t n = 32, m = 32, w = 32;
  IImplicitFunctionDPtr func( new ThinLayer<double>(0.25) );
  FillInGrid fill(func);

  Grid3D<double> grid(n, m, w, Box3D(10.0, 10.0, 10.0));
  fill.run(grid);

  Collision collision(&grid);

  // perpendicular to plane
  double dt = 1.0;
  MathVector3D pos(0.0, 0.5, 0.0);
  double currsdf = collision.computeSDF(pos);
  MathVector3D vel(0.0, 1.0, 0.0);
  bounceBack(collision, currsdf, dt, pos, vel);

  // It will always go throw very thin SDFs, don't know how to detect that
  //EXPECT_TRUE( (pos - MathVector3D(0.0, 0.75, 0.0)).getY() < tolerance );
}*/
