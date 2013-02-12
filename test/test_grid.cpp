//  (C) Copyright Kirill Lykov 2012.
//
// Distributed under the FreeBSD Software License (See accompanying file license.txt)

#include <gtest/gtest.h>

#define THROW_ON_INDEX_ERRORS
#include <exceptions.h>
#include <grid.h>

using namespace ls;

TEST(TwoDGridTest, construction)
{
  try {
    Grid2D<double> grid3(0, 0);
  }
  catch (const array_size_error& er)
  {
    return;
  }
  EXPECT_TRUE(false);
}

TEST(TwoDGridTest, size)
{
  size_t n = 3, m = 4;
  Grid2D<double> grid(n, m);
  EXPECT_EQ(grid.size(0), n);
  EXPECT_EQ(grid.size(1), m);
  try {
    grid.size(99);
  }
  catch (const std::out_of_range& er)
  {
    return;
  }
  EXPECT_TRUE(false);
}

TEST(TwoDGridTest, copyConstruction)
{
  Grid2D<double> grid1(1, 2);
  Grid2D<double> grid2(grid1);
  EXPECT_EQ(grid1, grid2);
}

TEST(TwoDGridTest, assignment)
{
  Grid2D<double> grid1(1, 2);
  Grid2D<double> grid2(1, 2);
  grid2 = grid1;
  EXPECT_EQ(grid1, grid2);
}

TEST(TwoDGridTest, assignmentBrake)
{
  Grid2D<double> grid1(1, 2);
  Grid2D<double> grid2(2, 3);
  try {
    grid2 = grid1;
  }
  catch (const array_size_error& er)
  {
    return;
  }
  EXPECT_TRUE(false);
}

TEST(TwoDGridTest, access)
{
  size_t n = 2, m = 4;
  Grid2D<double> grid(n, m);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      grid(i, j) = i * j;
    }
  }

  const Grid2D<double> grid2 = grid;
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      EXPECT_EQ(i * j, grid2(i, j));
    }
  }
}

TEST(TwoDGridTest, swap)
{
  size_t n1 = 3, m1 = 2;
  Grid2D<double> grid1(n1, m1);
  for (size_t i = 0; i < n1; ++i) {
    for (size_t j = 0; j < m1; ++j) {
      grid1(i, j) = i * j;
    }
  }

  size_t n2 = 5, m2 = 3;
  Grid2D<double> grid2(n2, m2);
  for (size_t i = 0; i < n2; ++i) {
    for (size_t j = 0; j < m2; ++j) {
      grid2(i, j) = i + j;
    }
  }

  std::swap(grid1, grid2);

  for (size_t i = 0; i < n1; ++i) {
    for (size_t j = 0; j < m1; ++j) {
      EXPECT_EQ(i * j, grid2(i, j));
    }
  }

  for (size_t i = 0; i < n2; ++i) {
    for (size_t j = 0; j < m2; ++j) {
      EXPECT_EQ(i + j, grid1(i, j));
    }
  }
}

TEST(TwoDGridTest, defaultConstructorAndResize)
{
  Grid2D<double> grid;
  size_t n = 2, m = 4;
  grid.resize(n, m);

  EXPECT_EQ(n, grid.size(0));
  EXPECT_EQ(m, grid.size(1));
}

TEST(TwoDGridTest, resize)
{
  size_t n = 2, m = 4;
  Grid2D<double> grid(n, m);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      grid(i, j) = i * j;
    }
  }

  n = 4;
  m = 6;
  grid.resize(n, m);

  EXPECT_EQ(n, grid.size(0));
  EXPECT_EQ(m, grid.size(1));

  // to check that memory is allocated and there is no access violation
  grid(3, 5) = 15.0;

  // it is not guaranteed up to this moment
  //for (size_t i = 0; i < n; ++i) {
  //  for (size_t j = 0; j < m; ++j) {
  //    EXPECT_EQ(copyOfGrid(i, j), i * j);
  //  }
  //}
}

TEST(ThreeDGridTest, construction)
{
  Grid3D<double> grid(1, 2, 3);
}

TEST(ThreeDGridTest, size)
{
  size_t n = 3, m = 4, w = 7;
  Grid3D<double> grid(n, m, w);
  EXPECT_EQ(grid.size(0), n);
  EXPECT_EQ(grid.size(1), m);
  EXPECT_EQ(grid.size(2), w);
  try {
    grid.size(99);
  }
  catch (const std::out_of_range& er)
  {
    return;
  }
  EXPECT_TRUE(false);
}

TEST(ThreeDGridTest, copyConstruction)
{
  Grid3D<double> grid1(1, 2, 3);
  Grid3D<double> grid2(grid1);
  EXPECT_EQ(grid1, grid2);
}

TEST(ThreeDGridTest, assignment)
{
  Grid3D<double> grid1(2, 2, 3);
  Grid3D<double> grid2(2, 2, 3);
  grid2 = grid1;
  EXPECT_EQ(grid1, grid2);
}

TEST(ThreeDGridTest, assignmentBrake)
{
  Grid3D<double> grid1(1, 2, 5);
  Grid3D<double> grid2(2, 3, 4);
  try {
    grid2 = grid1;
  }
  catch (const array_size_error& er)
  {
    return;
  }
  EXPECT_TRUE(false);
}

TEST(ThreeDGridTest, access)
{
  size_t n = 3, m = 4, w = 7;
  Grid3D<double> grid(n, m, w);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      for (size_t k = 0; k < w; ++k) {
        grid(i, j, k) = i * j * k;
      }
    }
  }

  const Grid3D<double> grid2 = grid;
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      for (size_t k = 0; k < w; ++k) {
        EXPECT_EQ(grid2(i, j, k), i * j * k);
      }
    }
  }
}


TEST(ThreeDGridTest, swap)
{
  size_t n1 = 3, m1 = 2, w1 = 4;
  Grid3D<double> grid1(n1, m1, w1);
  for (size_t i = 0; i < n1; ++i) {
    for (size_t j = 0; j < m1; ++j) {
      for (size_t k = 0; k < w1; ++k) {
        grid1(i, j, k) = i * j * k;
      }
    }
  }

  size_t n2 = 5, m2 = 3, w2 = 5;
  Grid3D<double> grid2(n2, m2, w2);
  for (size_t i = 0; i < n2; ++i) {
    for (size_t j = 0; j < m2; ++j) {
      for (size_t k = 0; k < w2; ++k) {
        grid2(i, j, k) = i + j + k;
      }
    }
  }

  std::swap(grid1, grid2);

  for (size_t i = 0; i < n1; ++i) {
    for (size_t j = 0; j < m1; ++j) {
      for (size_t k = 0; k < w1; ++k) {
        EXPECT_EQ(grid2(i, j, k), i * j * k);
      }
    }
  }

  for (size_t i = 0; i < n2; ++i) {
    for (size_t j = 0; j < m2; ++j) {
      for (size_t k = 0; k < w2; ++k) {
        EXPECT_EQ(grid1(i, j, k), i + j + k);
      }
    }
  }
}

