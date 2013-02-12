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

#include <deserializerHDF5.h>
#include <serializerHDF5.h>

using namespace ls;
using namespace io;
using namespace geometry_utils;


#include <gtest/gtest.h>

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

TEST(HDF5Test, writeAndRead1)
{
  size_t n = 8, m = 8, w = 16;
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

  BasicSerializerHDF5 writer(grid, "test-aux/writeAndRead1", "data");
  writer.run();

  Grid3D<double> readGrid; // at this point I don't know the size
  BasicDeserializerHDF5 reader(readGrid, "test-aux/writeAndRead1", "data");
  reader.run();

  EXPECT_EQ(grid, readGrid);
}

TEST(HDF5Test, writeAndRead2)
{
  size_t n = 10, m = 12, w = 14;
  Grid3D<double> grid(n, m, w);
  double h = 10.0 / (n - 1);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      for (size_t k = 0; k < w; ++k) {
        double point[] = { i * h - 5.0, j * h - 5.0, k * h - 5.0};
        grid(i, j, k) = testFunction2(point);
      }
    }
  }

  BasicSerializerHDF5 writer(grid, "test-aux/writeAndRead2", "data");
  writer.run();

  Grid3D<double> readGrid; // at this point I don't know the size
  BasicDeserializerHDF5 reader(readGrid, "test-aux/writeAndRead2", "data");
  reader.run();

  EXPECT_EQ(grid, readGrid);
}

//TODO improve test - check xdmf metadata
TEST(HDF5Test, writeAndRead3)
{
  size_t n = 12, m = 11, w = 14;
  Grid3D<double> grid(n, m, w);

  double top[] = {4.0, 5.0, 9.0};
  double low[] = {-3.0, -4.0, -5.0};
  Box box(low, top);
  double h[] = {box.getSizeX() / (n - 1.0), box.getSizeY() / (m - 1.0), box.getSizeZ() / (w - 1.0)};

  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      for (size_t k = 0; k < w; ++k) {
        double p[] = {i * h[0], j * h[1], k * h[2]};
        raw_math_vector::add(p, box.getLow());
        assert(box.inside(p));
        grid(i, j, k) = testFunction2(p);
      }
    }
  }

  BasicSerializerHDF5 writer(grid, "test-aux/writeAndRead3", "data", box);
  writer.run();

  Grid3D<double> readGrid; // at this point I don't know the size
  BasicDeserializerHDF5 reader(readGrid, "test-aux/writeAndRead3", "data");
  reader.run();

  EXPECT_EQ(grid, readGrid);
}
