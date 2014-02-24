//  (C) Copyright Kirill Lykov 2013.
//
// Distributed under the FreeBSD Software License (See accompanying file license.txt)

#include <box.h>
#include <tolerance.h>

#include <deserializerVDB.h>
//#include <serializerVDB.h>

using namespace ls;
using namespace io;
using namespace geometry_utils;

#include <gtest/gtest.h>

//TODO develop tests as soon as serializer will be written
TEST(VDBTest, mockup)
{
  ls::Grid3D<float> outputGrid;
  ls::geometry_utils::Box3D lsOutBox;
  ls::io::GridDeserializerVDB< ls::Grid3D<float>, ls::BasicWriteAccessStrategy<float> > vdbReader(outputGrid, "venusstatue", "ls_venus_statue", lsOutBox, 0.5f);
  //EXPECT_TRUE(vdbReader.run());
}

