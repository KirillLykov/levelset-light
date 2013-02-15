//  (C) Copyright Kirill Lykov 2012.
//
// Distributed under the FreeBSD Software License (See accompanying file license.txt)

#include <vector>
#include <algorithm>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <math.h>

#include <box.h>
#include <tolerance.h>

#include <deserializerVTK.h>
#include <serializerVTK.h>

using namespace ls;
using namespace io;
using namespace geometry_utils;

#include <gtest/gtest.h>
/*
 * declared in test_HDF5
namespace
{
  double testFunction1(const double* point) {
    double r = 0.5;
    double res = sqrt(pow(point[0], 2) + pow(point[1], 2) ) - r;
    return res;
  }

  double testFunction2(const double* point) {
    double res = sin(point[0]) + cos(point[1]) - log(fabs(point[3])) / (point[3] + 1.0);
    return res;
  }
}
*/

TEST(ToleranceTest, tol)
{
  double x[] = {1.2e-11, 1.1e-11};
  EXPECT_TRUE(ls::Tolerance::close(x[0], x[1]));
  EXPECT_FALSE(ls::Tolerance::close(0.00001, 0.01));
  EXPECT_FALSE(ls::Tolerance::close(10e10, 1.0));
}

TEST(VTKTest, writeAndRead1)
{
  Close_relative close_at_tol(10e-8); // for whatever reason vtk writes with error

  size_t n = 8, m = 8, w = 8;
  Box domainWrite(1.0, 1.0, 1.0);
  Grid3D<double> grid(n, m, w);
  double h = 1.0 / (n - 1);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      for (size_t k = 0; k < w; ++k) {
        double point[] = { i * h - 0.5, j * h - 0.5, k * h - 0.5};
        grid(i, j, k) = testFunction1(point);
      }
    }
  }

  BasicSerializerVTK writer(grid, "test-aux/writeAndRead1", domainWrite);
  EXPECT_TRUE(writer.run());

  Grid3D<double> readGrid; // at this point I don't know the size
  Box domainRead;
  BasicDeserializerVTK reader(readGrid, "test-aux/writeAndRead1", domainRead);
  EXPECT_TRUE(reader.run());

  EXPECT_EQ(domainWrite, domainRead);

  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      for (size_t k = 0; k < w; ++k) {
        if (grid(i, j, k) != readGrid(i, j, k)) {
          EXPECT_TRUE(close_at_tol(grid(i, j, k), readGrid(i, j, k)));
        }
      }
    }
  }
}
/*
typedef GridSerializerVTK< ls::Grid3D<double>, StructuredPointsTraits > BasicSerializerVTK_SPoint;
// The list of types we want to test.
typedef testing::Types<BasicSerializerVTK, BasicSerializerVTK_SPoint> Implementations;

TYPED_TEST_CASE(VTKTest, BasicSerializerVTK);
*/
TEST(VTKTest, writeAndRead3)
{
  Close_relative close_at_tol(10e-8); // for whatever reason vtk writes with error

  size_t n = 12, m = 11, w = 14;
  Grid3D<double> grid(n, m, w);

  double top[] = {4.0, 5.0, 9.0};
  double low[] = {-3.0, -4.0, -5.0};
  Box box(low, top);
  double h[] = {box.getSizeX() / (n - 1.0), box.getSizeY() / (m - 1.0), box.getSizeZ() / (w - 1.0)};

  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      for (size_t k = 0; k < w; ++k) {
        MathVector3D p(i * h[0], j * h[1], k * h[2]);
        p += box.getLow();
        assert(box.inside(p));
        grid(i, j, k) = testFunction2(p);
      }
    }
  }

  BasicSerializerVTK writer(grid, "test-aux/writeAndRead3", box);
  writer.run();

  Grid3D<double> readGrid; // at this point I don't know the size
  Box domainRead;
  BasicDeserializerVTK reader(readGrid, "test-aux/writeAndRead3", domainRead);
  reader.run();

  EXPECT_EQ(box, domainRead);

  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      for (size_t k = 0; k < w; ++k) {
        if (grid(i, j, k) != readGrid(i, j, k)) {
          EXPECT_TRUE(close_at_tol(grid(i, j, k), readGrid(i, j, k)));
        }
      }
    }
  }
}

typedef GridSerializerVTK< ls::Grid3D<double>, ls::BasicReadAccessStrategy<double>, StructuredPointsTraits > BasicSerializerVTK_SPoint;
typedef GridDeserializerVTK< ls::Grid3D<double>, ls::BasicWriteAccessStrategy<double>, StructuredPointsTraits > BasicDeserializerVTK_SPoint;

TEST(VTKTest, writeAndRead_SPoints)
{
  Close_relative close_at_tol(10e-6); // for whatever reason structured point have this tolerance

  size_t n = 8, m = 8, w = 8;
  Box domainWrite(1.0, 1.0, 1.0);
  Grid3D<double> grid(n, m, w);
  double h = 1.0 / (n - 1);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      for (size_t k = 0; k < w; ++k) {
        double point[] = { i * h - 0.5, j * h - 0.5, k * h - 0.5};
        grid(i, j, k) = testFunction1(point);
      }
    }
  }

  BasicSerializerVTK_SPoint writer(grid, "test-aux/writeAndRead3", domainWrite);
  EXPECT_TRUE(writer.run());

  Grid3D<double> readGrid; // at this point I don't know the size
  Box domainRead;
  BasicDeserializerVTK_SPoint reader(readGrid, "test-aux/writeAndRead3", domainRead);
  EXPECT_TRUE(reader.run());

  //EXPECT_EQ(domainWrite, domainRead); check with tolerance

  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      for (size_t k = 0; k < w; ++k) {
        if (grid(i, j, k) != readGrid(i, j, k)) {
          EXPECT_TRUE(close_at_tol(grid(i, j, k), readGrid(i, j, k)));
        }
      }
    }
  }
}

