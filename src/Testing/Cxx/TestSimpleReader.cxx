#include "vtkAutoInit.h" 
VTK_MODULE_INIT(vtkRenderingOpenGL2); // VTK was built with vtkRenderingOpenGL2
VTK_MODULE_INIT(vtkInteractionStyle);

#include "vtkActor.h"
#include "vtkCompositeDataDisplayAttributes.h"
#include "vtkCompositeDataGeometryFilter.h"
#include "vtkCompositeDataPipeline.h"
#include "vtkCompositeDataSet.h"
#include "vtkCompositePolyDataMapper2.h"
#include "vtkDataSet.h"
#include "vtkDataSetMapper.h"
#include "vtkDataSetSurfaceFilter.h"
#include "vtkDataSetReader.h"
#include "vtkGeometryFilter.h"
#include "vtkInformation.h"
#include "vtkLookupTable.h"
#include "vtkPointData.h"
#include "vtkPolyDataMapper.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkSmartPointer.h"
#include "vtkUnstructuredGrid.h"

#include <vtksys/SystemTools.hxx>
#include <vtksys/CommandLineArguments.hxx>

#define VTK_CREATE(type, var) \
  vtkSmartPointer<type> var = vtkSmartPointer<type>::New();



#define WITH_GRAPHICS 1
int
main(int argc, char **argv)
{
  vtkDataSetReader *reader = vtkDataSetReader::New();
  reader->DebugOff();
  reader->SetFileName("/tmp/foo.vtk");
  reader->Update();

  double range[2];
  reader->GetOutput()->GetPointData()->SetActiveScalars("Density");
  reader->GetOutput()->GetPointData()->GetArray(0)->GetRange(range);
  cerr << ": scalar range = [" << range[0] << ", " << range[1] << "]\n";

  //cout << *reader;


#ifdef WITH_GRAPHICS

  VTK_CREATE(vtkLookupTable, lut);
  lut->SetHueRange(0.66,0.0);
  lut->SetNumberOfTableValues(256);
  lut->SetScaleToLog10();
  lut->Build();

  lut->SetTableRange(range[0], range[1]);
  lut->Build();

  VTK_CREATE(vtkGeometryFilter, geom);
  geom->SetInputConnection(reader->GetOutputPort(0));

  VTK_CREATE(vtkPolyDataMapper, mapper1);
  mapper1->SetInputConnection(geom->GetOutputPort(0));
  mapper1->ScalarVisibilityOn();
  mapper1->SetScalarModeToUsePointFieldData();

  mapper1->SelectColorArray("Density");
  mapper1->SetLookupTable(lut);
  mapper1->UseLookupTableScalarRangeOn();

  VTK_CREATE(vtkActor, actor1);
  actor1->SetMapper(mapper1);

  VTK_CREATE(vtkRenderer, ren);
  VTK_CREATE(vtkRenderWindow, renWin);
  VTK_CREATE(vtkRenderWindowInteractor, iren);

  iren->SetRenderWindow(renWin);
  renWin->AddRenderer(ren);
  ren->AddActor(actor1);

  renWin->SetSize(512, 512);
  renWin->Render();
  ren->ResetCamera();

  renWin->Render();

  iren->Start();
#endif
  reader->Delete();
  }
