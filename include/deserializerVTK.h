//  (C) Copyright Kirill Lykov 2012.
//
// Distributed under the FreeBSD Software License (See accompanying file license.txt)


#ifndef DESERIALIZERVTK_H_
#define DESERIALIZERVTK_H_

#include <vector>
#include <string>
#include <stdexcept>
#include <memory.h>

#include <grid.h>
#include <basic_access_strategy.h>

#ifdef USE_VTK
#include <vtkSmartPointer.h>
#endif

#include "traits_VTK.h"

namespace ls
{
namespace io
{
#ifndef USE_VTK
template<class Grid, class AccessStrategy = BasicWriteAccessStrategy<double> >
class GridDeserializerVTK
{
public:
  GridDeserializerVTK(Grid& grid, const std::string& fileName)
  {}
  bool run() throw()
  {
    return false;
  }
};
#else

  template<class Grid, class AccessStrategy = BasicWriteAccessStrategy<double>, typename DeserializerTraits = ImageTraits >
  class GridDeserializerVTK : public AccessStrategy
  {
    typedef typename DeserializerTraits::Reader Reader;
    typedef typename DeserializerTraits::Representation Representation;

    typedef AccessStrategy _AS;
    typedef typename _AS::_Grid::_DataType _DataType;

    const std::string m_fileName;

  public:

    GridDeserializerVTK(Grid& grid, const std::string& fileName)
    : _AS(grid), m_fileName(fileName + std::string(".vti"))
    {
    }

    bool run() throw()
    {
      try
      {
        vtkSmartPointer<Reader> reader = vtkSmartPointer<Reader>::New();
        reader->SetFileName(m_fileName.c_str());
        reader->Update();

        vtkSmartPointer<Representation> strPoints = reader->GetOutput();

        int extent[6];
        strPoints->GetExtent(extent);
        _AS::m_grid.resize(extent[1] - extent[0] + 1, extent[3] - extent[2] + 1, extent[5] - extent[4] + 1);

        double spacing[3];
        strPoints->GetSpacing(spacing);
        double center[3];
        strPoints->GetOrigin(center);

        double bottom[3], top[3];
        for (size_t i = 0; i < 3; ++i) {
          bottom[i] = center[i] - spacing[i] * (_AS::m_grid.size(i) - 1) / 2.0;
          top[i] = center[i] + spacing[i] * (_AS::m_grid.size(i) - 1) / 2.0;
        }

        geometry_utils::Box3D newDomain(bottom, top);
        _AS::m_grid.setBoundingBox(newDomain);

        for (size_t iz = 0; iz < _AS::m_grid.size(2); ++iz)
          for (size_t iy = 0; iy < _AS::m_grid.size(1); ++iy)
            for (size_t ix = 0; ix < _AS::m_grid.size(0); ++ix)
            {
              double data = strPoints->GetScalarComponentAsDouble(ix, iy, iz, 0);
              _AS::setValue(ix, iy, iz, &data);
            }
      }
      catch(...)
      {
        return false;
      }

      return true;
    }
  };
#endif
  typedef GridDeserializerVTK<ls::Grid3D<double> > BasicDeserializerVTK;
}
}


#endif
