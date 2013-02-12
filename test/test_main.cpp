//  (C) Copyright Kirill Lykov 2012.
//
// Distributed under the FreeBSD Software License (See accompanying file license.txt)

#include <vector>
#include <algorithm>
#include <assert.h>
#include <iostream>

#include <gtest/gtest.h>

// to run all tests from one point
// corresponding cpp files must be excluded from build
// otherwise "duplicate symbol" linker error happens

#include "test_grid.cpp"

#include "test_hdf5.cpp"

#include "test_vtk.cpp"

#include "test_geometry_utils.cpp"

#include "test_linear_interpolator.cpp"

#include "test_implicit_functions.cpp"

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
