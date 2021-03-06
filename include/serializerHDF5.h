//  (C) Copyright Kirill Lykov 2012.
//
// Distributed under the FreeBSD Software License (See accompanying file license.txt)

#ifndef HDF5SERIALIZER_H_
#define HDF5SERIALIZER_H_

#include <vector>
#include <string>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <memory.h>

#ifdef USE_HDF5
#include <H5Cpp.h>
#endif

#include <grid.h>
#include <basic_access_strategy.h>
#include <box.h>

namespace ls
{
namespace io
{

#ifndef USE_HDF5
template<typename T, typename Grid = ls::Grid3D<T>, class AccessStrategy = ls::BasicReadAccessStrategy<T> >
class GridSerializerHDF5
{
public:
  GridSerializerHDF5(Grid& grid, const std::string& fileName, const std::string& datasetName)
  {}
  bool run() const throw()
  {
    return false;
  }
};
#else

#ifndef H5_NO_NAMESPACE
  using namespace H5;
#endif

  // TODO up to the moment work only with Grid3D<double>, extend it
  template<typename T, typename Grid = ls::Grid3D<T>, class AccessStrategy = ls::BasicReadAccessStrategy<T> >
  class GridSerializerHDF5 : public AccessStrategy
  {
    typedef AccessStrategy _AS;

    const std::string m_fullFileName;
    const std::string m_metadataFileName; //used by paraview
    const std::string m_datasetName;
    const geometry_utils::Box3 m_bbox; // bounding box for the grid, by default it is origin centerd 1x1x1

    static const hsize_t m_rank = 3;

  public:

    GridSerializerHDF5(const Grid& grid, const std::string& fullFileName,
        const std::string& datasetName)
    : _AS(grid), m_fullFileName(fullFileName + std::string(".h5")),
      m_metadataFileName(fullFileName + ".xmf"), m_datasetName(datasetName),
      m_bbox(grid.getBoundingBox())
    {
    }

    bool run() const throw()
    {
      try
      {
        H5::Exception::dontPrint();
        H5File file(m_fullFileName, H5F_ACC_TRUNC);

        DSetCreatPropList plist;
        hsize_t dims[] = {_AS::m_grid.size(0), _AS::m_grid.size(1), _AS::m_grid.size(2)};
        DataSpace fspace(m_rank, dims);
#ifdef SINGLE_PRECISION
		auto PT = PredType::NATIVE_FLOAT;
#else
		auto PT = PredType::NATIVE_DOUBLE;
#endif
		DataSet dataset(file.createDataSet(m_datasetName, PT, fspace, plist));
        
        // Box's attributes
        hsize_t six = 6;
        ArrayType dataType(PT, 1, &six);
        H5::DataSpace dataSpace = H5::DataSpace();
        Attribute attrBoxSize(dataset.createAttribute("boxSize", dataType, dataSpace)); 
        T data[] = {m_bbox.getLow().getX(), m_bbox.getLow().getY(), m_bbox.getLow().getZ(),
            m_bbox.getTop().getX(), m_bbox.getTop().getY(), m_bbox.getTop().getZ()};
        attrBoxSize.write(dataType, data); 

        hsize_t offset[m_rank]; // hyperslab offset in the file
        memset(offset, 0, m_rank * sizeof(hsize_t));
        fspace.selectHyperslab(H5S_SELECT_SET, dims, offset);

        std::vector<T> buffer( computeBufferSize(dims) );
        fillBufferFromGrid(dims, buffer);

        dataset.write(&buffer[0], PT);

        writeMetadata(dims);
      }
      // catch failure caused by the H5File operations
      catch (const FileIException& error) {
        return false;
      }
      // catch failure caused by the DataSet operations
      catch (const DataSetIException& error) {
        return false;
      }
      // catch failure caused by the DataSpace operations
      catch (const DataSpaceIException& error) {
        return false;
      }
      // catch failure caused by the DataSpace operations
      catch (const DataTypeIException& error ) {
        return false;
      }
      catch (...) {
        return false;
      }

      return true;
    }


  private:

    hsize_t computeBufferSize(const hsize_t* dims) const
    {
      hsize_t res = 1;
      for (size_t i= 0; i < m_rank; ++i)
        res *= dims[i];
      return res;
    }

    /**
     * writing 3D data only
     */
    void fillBufferFromGrid(const hsize_t* dims, std::vector<T>& buffer) const
    {
      for (size_t iz = 0; iz < _AS::m_grid.size(2); ++iz)
        for (size_t iy = 0; iy < _AS::m_grid.size(1); ++iy)
          for (size_t ix = 0; ix < _AS::m_grid.size(0); ++ix)
          {
            T data = _AS::getValue(ix, iy, iz);
            buffer[iz + _AS::m_grid.size(2) * (iy + _AS::m_grid.size(1) * ix)] = data;
          }
    }

    void writeMetadata(const hsize_t* dims) const
    {
      std::stringstream ss;
      geometry_utils::MathVector3D bboxCenter;
      m_bbox.getCenter(bboxCenter);
#ifdef SINGLE_PRECISION
	  int precision = 4;
#else
	  int precision = 8;
#endif

      //Note, that z and x are swapped for xdmf
      ss << "<?xml version=\"1.0\" ?>\n" <<
          "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n" <<
          "<Xdmf Version=\"2.0\">\n" <<
          " <Domain>\n" <<
          "   <Grid GridType=\"Uniform\">\n" <<
          "     <Time Value=\"0\"/>\n" <<
          "     <Topology TopologyType=\"3DCORECTMesh\" Dimensions=\"" <<
                  dims[0] << " "<< dims[1] << " " << dims[2] << "\"/>\n" <<
          "     <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n"  <<
          "       <DataItem Name=\"Origin\" Dimensions=\"3\" NumberType=\"Float\" Precision=\"" << precision << "\" Format=\"XML\">\n"  <<
          "        " << bboxCenter.getX() << " " << bboxCenter.getY() << " " << bboxCenter.getZ() << "\n" <<
          "       </DataItem>\n"  <<
          "       <DataItem Name=\"Spacing\" Dimensions=\"3\" NumberType=\"Float\" Precision=\"" << precision << "\" Format=\"XML\">\n"  <<
          "        " << m_bbox.getIthSize(0) / dims[0] << " " << m_bbox.getIthSize(1) / dims[1] << " " << m_bbox.getIthSize(2) / dims[2]  << "\n" <<
          "       </DataItem>\n"  <<
          "     </Geometry>\n"  <<
          "     <Attribute Name=\"data\" AttributeType=\"\" Center=\"Node\">\n"
          "       <DataItem Dimensions=\"" << dims[0] << " " << dims[1] << " " << dims[2] << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n" <<
          "        " << extractFilename(m_fullFileName) << ":/data\n" <<
          "       </DataItem>\n"  <<
          "     </Attribute>\n"  <<
          "   </Grid>\n"  <<
          " </Domain>\n"  <<
          "</Xdmf>\n";
      std::ofstream metadataFile(m_metadataFileName.c_str());
      metadataFile << ss.str();
      metadataFile.close();
    }

    static std::string extractFilename(const std::string& fullFileName)
    {
      size_t index = fullFileName.find_last_of( '\\' );
      if (index == std::string::npos)
        index = fullFileName.find_last_of( '/' );
      return fullFileName.substr(index + 1);
    }
  };
#endif

#ifdef SINGLE_PRECISION
  typedef GridSerializerHDF5<float> BasicSerializerHDF5;
#else
  typedef GridSerializerHDF5<double> BasicSerializerHDF5;
#endif
}
}

#endif /* HDF5SERIALIZER_H_ */
