//  (C) Copyright Kirill Lykov 2012.
//
// Distributed under the FreeBSD Software License (See accompanying file license.txt)

#include <gtest/gtest.h>
#include <math_vector.h>
#include <box.h>

//TODO write more tests

TEST(VectorTest, vector3D)
{
  using namespace geometry_utils;

  MathVector3D v1(0.0, 0.0, 0.0);
  MathVector3D v2(1.0, 2.0, 3.0);
  MathVector3D v3 = v1 + v2;
  EXPECT_EQ(v2, v3);

  v3 *= 0.0;
  EXPECT_EQ(v1, v3);

  v3 += MathVector3D(1.0, 0.0, 1.0);
  EXPECT_EQ(v3 * v2, 4.0);
}

TEST(VectorTest, vector2D)
{
  using namespace geometry_utils;

  MathVector2D v1(0.0, 0.0);
  MathVector2D v2(1.0, 2.0);
  MathVector2D v3 = v1 + v2;
  EXPECT_EQ(v2, v3);

  v3 *= 0.0;
  EXPECT_EQ(v1, v3);

  v3 += MathVector2D(1.0, 0.0);
  EXPECT_EQ(v3 * v2, 1.0);
}

TEST(BoxTest, box)
{
  using namespace geometry_utils;

  Box box(5.6, 10.0, 21.0);
  EXPECT_EQ(5.6, box.getSizeX());
  EXPECT_EQ(10.0, box.getSizeY());
  EXPECT_EQ(21.0, box.getSizeZ());
  MathVector3D center;
  box.getCenter(center);
  EXPECT_EQ(0.0, center.getX());
  EXPECT_EQ(0.0, center.getY());
  EXPECT_EQ(0.0, center.getZ());
}
