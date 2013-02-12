//  (C) Copyright Kirill Lykov 2012.
//
// Distributed under the FreeBSD Software License (See accompanying file license.txt)

#ifndef HDF5DESERIALIZER_H_
#define HDF5DESERIALIZER_H_

#include <vector>
#include <string>
#include <stdexcept>
#include <memory.h>

#include <H5Cpp.h>

#include <grid.h>
#include <basic_access_strategy.h>

namespace ls
{
namespace io
{

#ifndef USE_HDF5
template<class Grid, class AccessStrategy = BasicWriteAccessStrategy<double> >
class GridDeserializerHDF5
{
public:
  GridDeserializerHDF5(Grid& grid, const std::string& fileName, const std::string& datasetName)
  {}
  bool run() throw()
  {
    return false;
  }
};
#else

#ifndef H5_NO_NAMESPACE
  using namespace H5;
#endif

  //TODO up to this moment it works only with Grid3D<double>, think how to extend it
  //TODO add reading xml metadata using XDMF library
  template<class Grid, class AccessStrategy = BasicWriteAccessStrategy<double> >
  class GridDeserializerHDF5 : public AccessStrategy
  {
    typedef AccessStrategy _AS;
    typedef std::vector<hsize_t> _szArray;

    const std::string m_fileName;
    const std::string m_datasetName;

  public:

    GridDeserializerHDF5(Grid& grid, const std::string& fileName)
    : _AS(grid), m_fileName(fileName + std::string(".h5")), m_datasetName("data")
    {
    }

    GridDeserializerHDF5(Grid& grid, const std::string& fileName, const std::string& datasetName)
    : _AS(grid), m_fileName(fileName + std::string(".h5")), m_datasetName(datasetName)
    {
    }

    bool run() throw()
    {
      try
      {
        H5File file(m_fileName, H5F_ACC_RDONLY);
        DataSet dataset = file.openDataSet(m_datasetName);
        H5T_class_t typeClass = dataset.getTypeClass();

        if (typeClass != H5T_FLOAT) //TODO change exception type
          throw std::logic_error("only double is supported by HDF5Deserializer");

        //validateRank();
        DataSpace dataspace = dataset.getSpace();
        int rank = dataspace.getSimpleExtentNdims();
        if (rank < 3) //TODO change exception type
          throw std::logic_error("rank of the data in the hdf5 file is not 3");

        _szArray dims_out(rank);
        dataspace.getSimpleExtentDims(&dims_out[0]);

        _szArray offset(rank); // hyperslab offset in the file
        memset(&offset[0], 0, rank * sizeof(hsize_t));

        dataspace.selectHyperslab(H5S_SELECT_SET, &dims_out[0], &offset[0]);

        DataSpace memspace(rank, &dims_out[0]);
        memspace.selectHyperslab(H5S_SELECT_SET, &dims_out[0], &offset[0]);

        std::vector<double> buffer( computeBufferSize(dims_out) );

        dataset.read(&buffer[0], PredType::NATIVE_DOUBLE, memspace, dataspace);

        fillGridFromBuffer(dims_out, buffer);
      }
      // catch failure caused by the H5File operations
      catch(const FileIException& error) {
        error.printError();
        return false;
      }
      // catch failure caused by the DataSet operations
      catch(const DataSetIException& error) {
        error.printError();
        return false;
      }
      // catch failure caused by the DataSpace operations
      catch(const DataSpaceIException& error) {
        error.printError();
        return false;
      }
      // catch failure caused by the DataSpace operations
      catch(const DataTypeIException& error ) {
        error.printError();
        return false;
      }
      catch(...)
      {
        return false;
      }

      return true;
    }


  private:

    hsize_t computeBufferSize(const _szArray& dims) const
    {
      hsize_t res = 1;
      for (size_t i= 0; i < dims.size(); ++i)
        res *= dims[i];
      return res;
    }

    /**
     * Although dims may be any size, it is assumed that it is 3 or 4
     *  but in this case assumed that dim[3] == 1.
     */
    void fillGridFromBuffer(const _szArray& dims, const std::vector<double>& buffer)
    {
      if (dims.size() == 4 && dims[3] != 1)
        throw std::logic_error("datasets with 4th dimension not equal to 1 are not supported");

      _AS::m_grid.resize(dims[0], dims[1], dims[2]);
      for (size_t iz = 0; iz < _AS::m_grid.size(2); ++iz)
        for (size_t iy = 0; iy < _AS::m_grid.size(1); ++iy)
          for (size_t ix = 0; ix < _AS::m_grid.size(0); ++ix)
          {
            //TODO should be const T*
            const double* data = &buffer[iz + _AS::m_grid.size(2) * (iy + _AS::m_grid.size(1) * ix)];
            _AS::setValue(ix, iy, iz, data);
          }
    }
  };
#endif
  typedef GridDeserializerHDF5<ls::Grid3D<double> > BasicDeserializerHDF5;
}
}

#endif /* HDF5DATASET_H_ */
