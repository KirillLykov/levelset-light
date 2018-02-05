//  (C) Copyright Kirill Lykov 2012.
//
// Distributed under the FreeBSD Software License (See accompanying file license.txt)

#include <gtest/gtest.h>
#include <math_vector.h>
#include <box.h>

TEST(VectorTest, vector3D)
{
  using namespace geometry_utils;

  MathVector3R v1(0.0, 0.0, 0.0);
  MathVector3R v2(1.0, 2.0, 3.0);
  MathVector3R v3 = v1 + v2;
  EXPECT_EQ(v2, v3);

  v3 *= Real(0.0);
  EXPECT_EQ(v1, v3);

  v3 += MathVector3R(1.0, 0.0, 1.0);
  EXPECT_EQ(v3 * v2, 4.0);
}

TEST(VectorTest, vector2D)
{
  using namespace geometry_utils;

  MathVector2R v1(0.0, 0.0);
  MathVector2R v2(1.0, 2.0);
  MathVector2R v3 = v1 + v2;
  EXPECT_EQ(v2, v3);

  v3 *= Real(0.0);
  EXPECT_EQ(v1, v3);

  v3 += MathVector2R(1.0, 0.0);
  EXPECT_EQ(v3 * v2, 1.0);
}

TEST(VectorTest, vectorSwap)
{
  using namespace geometry_utils;

  MathVector3R v1(-1.0, -2.0, -3.0);
  MathVector3R v2(1.0, 2.0, 3.0);
  std::swap(v1, v2);
  EXPECT_EQ(1.0, v1.getX());
  EXPECT_EQ(2.0, v1.getY());
  EXPECT_EQ(3.0, v1.getZ());
  EXPECT_EQ(-1.0, v2.getX());
  EXPECT_EQ(-2.0, v2.getY());
  EXPECT_EQ(-3.0, v2.getZ());
}

TEST(BoxTest, box3D)
{
  using namespace geometry_utils;

  Box3R box(5.6, 10.0, 21.0);
  EXPECT_NEAR(5.6, box.getSizeX(), tol);
  EXPECT_NEAR(10.0, box.getSizeY(), tol);
  EXPECT_EQ(21.0, box.getSizeZ());
  MathVector3R center;
  box.getCenter(center);
  EXPECT_EQ(0.0, center.getX());
  EXPECT_EQ(0.0, center.getY());
  EXPECT_EQ(0.0, center.getZ());
}

TEST(BoxTest, box3dPointOnBorder)
{
  using namespace geometry_utils;

  MathVector3R low(-42.5, -24.80000000000000, -14);
  MathVector3R top(42.5, 24.80000000000000, 14);
  Box3R box(low, top);
  MathVector3R p(-13.886428108937016, -24.800000000026049, 10.568595969105855);
  EXPECT_TRUE(box.inside(p));
}

TEST(BoxTest, box2D)
{
  using namespace geometry_utils;

  Box2R box(5.6, 10.0);
  EXPECT_NEAR(5.6, box.getSizeX(), tol);
  EXPECT_NEAR(10.0, box.getSizeY(), tol);
  MathVector2R center;
  box.getCenter(center);
  EXPECT_EQ(0.0, center.getX());
  EXPECT_EQ(0.0, center.getY());
}

TEST(BoxTest, box3Dswap)
{
  using namespace geometry_utils;

  Box3R box1(5.6, 10.0, 21.0);
  Box3R box2;
  std::swap(box1, box2);
  EXPECT_NEAR(5.6, box2.getSizeX(), tol);
  EXPECT_NEAR(10.0, box2.getSizeY(), tol);
  EXPECT_NEAR(21.0, box2.getSizeZ(), tol);
  EXPECT_NEAR(0.0, box1.getSizeX(), tol);
  EXPECT_NEAR(0.0, box1.getSizeY(), tol);
  EXPECT_NEAR(0.0, box1.getSizeZ(), tol);
}
