//  (C) Copyright Kirill Lykov 2012.
//
// Distributed under the FreeBSD Software License (See accompanying file license.txt)

#include <gtest/gtest.h>
#include <implicit_functions.h>
#include <tolerance.h>

using namespace ls;
using namespace geometry_utils;

TEST(ImpilictFunctionsTest, AxialCylinder)
{
  ls::AxialCylinderD axialCylinder(ls::zDim, 1.0);
  EXPECT_TRUE( Tolerance<Real>::close_at_tol(axialCylinder.compute(MathVector3D(0.0, 0.0, 0.0)), -1.0) );
  EXPECT_TRUE( Tolerance<Real>::close_at_tol(axialCylinder.compute(MathVector3D(1.0, 0.0, 0.0)), 0.0) );
  EXPECT_TRUE( Tolerance<Real>::close_at_tol(axialCylinder.compute(MathVector3D(2.0, 0.0, 0.0)), 1.0) );
}

TEST(ImpilictFunctionsTest, NonAxialCylinder)
{
  ls::NonAxialCylinderD nonAxialCylinder(MathVector3D(0.0, 0.0, 0.0), MathVector3D(0.0, 0.0, 1.0), 1.0, 1.0);
  EXPECT_TRUE( Tolerance<Real>::close_at_tol(nonAxialCylinder.compute(MathVector3D(0.0, 0.0, 0.0)), -1.0) );
  EXPECT_TRUE( Tolerance<Real>::close_at_tol(nonAxialCylinder.compute(MathVector3D(1.0, 0.0, 0.0)), 0.0) );
  EXPECT_TRUE( Tolerance<Real>::close_at_tol(nonAxialCylinder.compute(MathVector3D(2.0, 0.0, 0.0)), 1.0) );
  EXPECT_TRUE( Tolerance<Real>::close_at_tol(nonAxialCylinder.compute(MathVector3D(0.0, 0.0, 2.0)), 1.5) );
  EXPECT_TRUE( Tolerance<Real>::close_at_tol(nonAxialCylinder.compute(MathVector3D(0.0, 0.0, -2.0)), 1.5) );
}

TEST(ImpilictFunctionsTest, Sphere)
{
  ls::SphereD sphere(MathVector3D(0.0, 0.0, 0.0), 1.0);
  EXPECT_TRUE( Tolerance<Real>::close_at_tol(sphere.compute(MathVector3D(0.0, 0.0, 0.0)), -1.0) );
  EXPECT_TRUE( Tolerance<Real>::close_at_tol(sphere.compute(MathVector3D(1.0, 0.0, 0.0)), 0.0) );
  EXPECT_TRUE( Tolerance<Real>::close_at_tol(sphere.compute(MathVector3D(2.0, 0.0, 0.0)), 1.0) );
}

TEST(ImpilictFunctionsTest, AxialTorus)
{
  ls::AxialTorusD axialTorus(ls::zDim, MathVector3D(0.0, 0.0, 0.0), 1.5, 0.5);
  EXPECT_TRUE( Tolerance<Real>::close_at_tol(axialTorus.compute(MathVector3D(0.0, 0.0, 0.0)), 1.0) );
  EXPECT_TRUE( Tolerance<Real>::close_at_tol(axialTorus.compute(MathVector3D(1.5, 0.0, 0.0)), -0.5) );
  EXPECT_TRUE( Tolerance<Real>::close_at_tol(axialTorus.compute(MathVector3D(2.0, 0.0, 0.0)), 0.0) );
  EXPECT_TRUE( Tolerance<Real>::close_at_tol(axialTorus.compute(MathVector3D(3.0, 0.0, 0.0)), 1.0) );
}

TEST(ImpilictFunctionsTest, CutByPlane)
{
  std::shared_ptr<ls::SphereD> pSphere( new ls::SphereD(MathVector3D(0.0, 0.0, 0.0), 1.0) );
  ls::CutByPlaneD cut(pSphere, MathVector3D(0.0, 0.0, 1.0), MathVector3D(0.0, 0.0, 0.0));
  EXPECT_TRUE( Tolerance<Real>::close_at_tol(cut.compute(MathVector3D(0.0, 0.0, 0.0)), -1.0) );
  EXPECT_TRUE( Tolerance<Real>::close_at_tol(cut.compute(MathVector3D(1.0, 0.0, 0.0)), 0.0) );
  EXPECT_TRUE( Tolerance<Real>::close_at_tol(cut.compute(MathVector3D(2.0, 0.0, 0.0)), 1.0) );
  EXPECT_TRUE( Tolerance<Real>::close_at_tol(cut.compute(MathVector3D(0.0, 0.0, 1.0)), 1.4142135623730951) );
}

TEST(ImpilictFunctionsTest, Operations)
{
  IImplicitFunctionDPtr sphere( new ls::SphereD(MathVector3D(0.0, 0.0, 0.0), 1.0) );
  IImplicitFunctionDPtr cut( new ls::CutByPlaneD(sphere, MathVector3D(0.0, 0.0, 1.0), MathVector3D(0.0, 0.0, 0.0)) );
  IImplicitFunctionDPtr axialTorus( new ls::AxialTorusD(ls::zDim, MathVector3D(0.0, 0.0, 0.0), 1.5, 0.5) );
  IImplicitFunctionDPtr nonAxialCylinder( new ls::NonAxialCylinderD(MathVector3D(0.0, 0.0, 0.0), MathVector3D(0.0, 0.0, 1.0), 1.0, 1.0) );
  //due to simplicity of this operation, I don't do any math check
  FunctionsD functions;
  functions.push_back(cut);
  functions.push_back(axialTorus);
  functions.push_back(nonAxialCylinder);

  ls::UnionD un(functions);
  ls::UnionD intersect(functions);
  ls::UnionD diff(functions);
}

