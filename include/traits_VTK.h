//  (C) Copyright Kirill Lykov 2012.
//
// Distributed under the FreeBSD Software License (See accompanying file license.txt)

#ifndef TRAITS_VTK_H_
#define TRAITS_VTK_H_

#ifdef USE_VTK
#include <vtkStructuredPointsReader.h>
#include <vtkStructuredPointsWriter.h>
#include <vtkStructuredPoints.h>

#include <vtkXMLImageDataReader.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkImageData.h>
#endif

namespace ls
{
namespace io
{
#ifdef USE_VTK
  struct StructuredPointsTraits
  {
    typedef vtkStructuredPointsReader Reader;
    typedef vtkStructuredPointsWriter Writer;
    typedef vtkStructuredPoints Representation;
  };

  struct ImageTraits
  {
    typedef vtkXMLImageDataReader Reader;
    typedef vtkXMLImageDataWriter Writer;
    typedef vtkImageData Representation;
  };
#endif
}
}

#endif /* TRAITS_VTK_H_ */
