//  (C) Copyright Kirill Lykov 2012.
//
// Distributed under the FreeBSD Software License (See accompanying file license.txt)

#ifndef SERIALIZERVTK_H_
#define SERIALIZERVTK_H_

#include <vector>
#include <string>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <memory.h>

#include <grid.h>
#include <basic_access_strategy.h>
#include <box.h>
#include <math_vector.h>

#ifdef USE_VTK
#include <vtkSmartPointer.h>
#endif

#include "traits_VTK.h"

namespace ls
{
namespace io
{
#ifndef USE_VTK
template<class Grid, class AccessStrategy = ls::BasicReadAccessStrategy<double> >
class GridSerializerVTK
{
public:
  GridSerializerVTK(Grid& grid, const std::string& fileName, const std::string& datasetName)
  {}
  bool run() const throw()
  {
    return false;
  }
};
#else

  // TODO up to the moment work only with Grid3D<double>, extend it
  template<class Grid, class AccessStrategy = ls::BasicReadAccessStrategy<double>, typename SerializerTraits = ImageTraits >
  class GridSerializerVTK : public AccessStrategy
  {
    typedef typename SerializerTraits::Writer Writer;
    typedef typename SerializerTraits::Representation Representation;

    typedef AccessStrategy _AS;

    const std::string m_fullFileName;
    const geometry_utils::Box m_bbox; // bounding box for the grid

  public:

    GridSerializerVTK(const Grid& grid, const std::string& fullFileName,
        const geometry_utils::Box& box = geometry_utils::Box(1.0))
    : _AS(grid), m_fullFileName(fullFileName + std::string(".vti")), m_bbox(box)
    {
    }

    bool run() const throw()
    {
      try
      {
        vtkSmartPointer<Representation> img = Representation::New();

        img->SetExtent(0, _AS::m_grid.size(0) - 1, 0, _AS::m_grid.size(1) - 1, 0, _AS::m_grid.size(2) - 1);
        img->SetScalarTypeToDouble();
        img->SetNumberOfScalarComponents(1);

        img->AllocateScalars();
        img->SetSpacing(m_bbox.getSizeX()/(_AS::m_grid.size(0) - 1), m_bbox.getSizeY()/(_AS::m_grid.size(1) - 1),
            m_bbox.getSizeZ()/(_AS::m_grid.size(2) - 1));
        geometry_utils::MathVector3D origin;
        m_bbox.getCenter(origin);
        double rawOrigin[] = {origin.getX(), origin.getY(), origin.getZ()};
        img->SetOrigin(rawOrigin);

        for (size_t iz = 0; iz < _AS::m_grid.size(2); ++iz)
          for (size_t iy = 0; iy < _AS::m_grid.size(1); ++iy)
            for (size_t ix = 0; ix < _AS::m_grid.size(0); ++ix)
            {
              double data = _AS::getValue(ix, iy, iz);
              img->SetScalarComponentFromFloat(ix, iy, iz, 0, data);
            }

        vtkSmartPointer<Writer> writer = Writer::New();
        writer->SetFileName(m_fullFileName.c_str());
        writer->SetInput(img);
        writer->Write();
      }
      catch(...) {
        return false;
      }
      return true;
    }
  };
#endif

  typedef GridSerializerVTK<ls::Grid3D<double> > BasicSerializerVTK;
}
}

#endif /* SERIALIZERVTK_H_ */
