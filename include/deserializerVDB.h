//  (C) Copyright Kirill Lykov 2013.
//
// Distributed under the FreeBSD Software License (See accompanying file license.txt)

#ifndef DESERIALIZERVDB_H_
#define DESERIALIZERVDB_H_

#include <vector>
#include <string>
#include <stdexcept>
#include <memory.h>

#include <openvdb/openvdb.h>
#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/GridSampling.h>

#include <grid.h>
#include <basic_access_strategy.h>

namespace ls
{
namespace io
{

#ifndef USE_VDB
template<class Grid, class AccessStrategy = BasicWriteAccessStrategy<double> >
class GridDeserializerVDB
{
public:
  GridDeserializerVDB(Grid& grid, const std::string& fileName)
  {}

  GridDeserializerVDB(Grid& grid, const std::string& fileName, const std::string& datasetName, geometry_utils::Box& domain, double coarsenCoeff = 1.0)
  {}
  bool run() throw()
  {
    return false;
  }
};
#else

  template<class Grid, class AccessStrategy = BasicWriteAccessStrategy<double> >
  class GridDeserializerVDB : public AccessStrategy
  {
    typedef AccessStrategy _AS;
    typedef typename _AS::_Grid::_DataType _DataType;

    const std::string m_fileName;
    const std::string m_datasetName;

    geometry_utils::Box& m_domain;
    _DataType m_coarsenCoeff; // used because usually vdb meshes are too big when transformed to ls

  public:

    GridDeserializerVDB(Grid& grid, const std::string& fileName)
    : _AS(grid), m_fileName(fileName + std::string(".vdb")), m_datasetName("data"), m_domain(1.0) ,m_coarsenCoeff(1.0)
    {
    }

    GridDeserializerVDB(Grid& grid, const std::string& fileName, const std::string& datasetName, geometry_utils::Box& domain, _DataType coarsenCoeff = 1.0)
    : _AS(grid), m_fileName(fileName + std::string(".vdb")), m_datasetName(datasetName), m_domain(domain), m_coarsenCoeff(coarsenCoeff)
    {
    }

    bool run() throw()
    {
      assert(m_coarsenCoeff != 0.0);
      try
      {
        // load vdb grid
        openvdb::initialize();

        openvdb::io::File inputFile(m_fileName);
        inputFile.open();
        openvdb::GridBase::Ptr baseGrid = inputFile.readGrid(m_datasetName);
        inputFile.close();

        openvdb::FloatGrid::Ptr inputGrid = openvdb::gridPtrCast<openvdb::FloatGrid>(baseGrid);
        openvdb::CoordBBox inputIndexBb = inputGrid->evalActiveVoxelBoundingBox();
        openvdb::Coord dims = inputGrid->evalActiveVoxelDim();
        openvdb::Coord dimsOut(static_cast<openvdb::Int32>(dims.x() * m_coarsenCoeff),
                      static_cast<openvdb::Int32>(dims.y() * m_coarsenCoeff),
                      static_cast<openvdb::Int32>(dims.z() * m_coarsenCoeff));

        openvdb::Vec3d h = inputGrid->voxelSize();

        openvdb::tools::GridSampler<openvdb::FloatTree, openvdb::tools::BoxSampler> interpolator(inputGrid->constTree(), inputGrid->transform());

        openvdb::Vec3d top(inputIndexBb.max().x() * h.x(),
                  inputIndexBb.max().y() * h.y(),
                  inputIndexBb.max().z() * h.z());
        openvdb::Vec3d low(inputIndexBb.min().x() * h.x(),
                  inputIndexBb.min().y() * h.y(),
                  inputIndexBb.min().z() * h.z());
        openvdb::BBoxd box(low, top);

        // outputGrid voxel sizes
        openvdb::Vec3d hhout = 1.0f/m_coarsenCoeff * h;

        // set in output grid size and bounding box
        _AS::m_grid.resize(dimsOut.x(), dimsOut.y(), dimsOut.z());
        m_domain = vdbBox2lsBox(box);

        for (size_t i = 0; i < dimsOut.x(); ++i) {
          for (size_t j = 0; j < dimsOut.y(); ++j) {
            for (size_t k = 0; k < dimsOut.z(); ++k) {
              openvdb::Vec3d p(i * hhout[0], j * hhout[1], k * hhout[2]);
              p += low;
              _DataType interpolatedValue = interpolator.wsSample(p);

              _AS::setValue(i, j, k, &interpolatedValue);
            }
          }
        }
      }
      catch(const std::exception& error)
      {
        return false;
      }
      catch(...)
      {
        return false;
      }

      return true;
    }
  private:
    ls::geometry_utils::MathVector3D vdbVector2lsVector(const openvdb::Vec3d& inputVector)
    {
      return ls::geometry_utils::MathVector3D(inputVector.asPointer());
    }

    ls::geometry_utils::Box vdbBox2lsBox(const openvdb::BBoxd& inputBox)
    {
      return ls::geometry_utils::Box(vdbVector2lsVector(inputBox.min()), vdbVector2lsVector(inputBox.max()));
    }
  };
#endif
  typedef GridDeserializerVDB<ls::Grid3D<double> > BasicDeserializerVDB;
}
}

#endif /* HDF5DATASET_H_ */
