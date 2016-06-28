//  (C) Copyright Kirill Lykov 2016.
//
// Distributed under the FreeBSD Software License (See accompanying file license.txt)

#include <gtest/gtest.h>
#include <math_vector.h>
#include <collision.h>
#include <box.h>

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
