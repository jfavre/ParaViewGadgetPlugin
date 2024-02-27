#include "vtkAutoInit.h" 
VTK_MODULE_INIT(vtkRenderingOpenGL2); // VTK was built with vtkRenderingOpenGL2
VTK_MODULE_INIT(vtkInteractionStyle);

#include "vtkActor.h"
#include "vtkDataSet.h"
#include "vtkDataSetWriter.h"
#include "vtkGadgetReader.h"
#include "vtkGeometryFilter.h"
#include "vtkInformation.h"
#include "vtkLookupTable.h"
#include "vtkPartitionedDataSetCollection.h"
#include "vtkPartitionedDataSet.h"
#include "vtkNew.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
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

int
main(int argc, char **argv)
{
  std::string filein;
  std::string varname;
  bool vis;

  double TimeStep = 0.0;
  int k, BlockIndex = 0;

  vtksys::CommandLineArguments args;
  args.Initialize(argc, argv);
  args.AddArgument(
    "-f", vtksys::CommandLineArguments::SPACE_ARGUMENT, &filein, "(the names of the Gadget (HDF5) files to read)");
  args.AddArgument(
    "-var", vtksys::CommandLineArguments::SPACE_ARGUMENT, &varname, "(the name of the SCALAR variable to display)");
  args.AddArgument(
    "-vis", vtksys::CommandLineArguments::NO_ARGUMENT, &vis, "(optional vis(ualization) with display)");

  if ( !args.Parse() || argc == 1 || filein.empty())
    {
    cerr << "\nTestSimpleGadgetReader: Written by Jean M. Favre\n"
         << "options are:\n";
    cerr << args.GetHelp() << "\n";
    exit(1);
    }

  if(!vtksys::SystemTools::FileExists(filein.c_str()))
    {
    cerr << "\nFile " << filein.c_str() << " does not exist\n\n";
    exit(1);
    }

  vtkGadgetReader *reader = vtkGadgetReader::New();
  reader->DebugOff();
  reader->SetDirectoryName(filein.c_str());
  reader->UpdateInformation();
  reader->SetCellType(2);
  for(auto i=0; i < reader->GetNumberOfPointArrays(); i++)
    cout << "found array (" << i << ") = " << reader->GetPointArrayName(i) << endl;
  reader->DisableAllPointArrays();
  reader->SetPointArrayStatus(varname.c_str(), 1);

  reader->UpdateTimeStep(TimeStep); // time value
  reader->Update();

  double range[2];
  vtkDataSet *FirstBlock = static_cast<vtkDataSet *>(reader->GetOutput()->GetPartitionedDataSet(0)->GetPartition(0));
  if(varname.size())
    {
    FirstBlock->GetPointData()->GetArray(0)->GetRange(range);
    cerr << varname.c_str() << ": scalar range = [" << range[0] << ", " << range[1] << "]\n";
    }
  //cout << *reader;
  /*
  VTK_CREATE(vtkDataSetWriter, writer);
  writer->SetInputData(FirstBlock);
  writer->SetFileTypeToBinary();
  writer->SetFileName("/tmp/foo.vtk");
  writer->Write();
  */
  
if(vis)
  {
  VTK_CREATE(vtkLookupTable, lut);
  lut->SetHueRange(0.66,0.0);
  lut->SetNumberOfTableValues(256);
  lut->SetScaleToLog10();
  lut->Build();

  if(varname.size())
    {
    lut->SetTableRange(range[0], range[1]);
    lut->Build();
    }

  VTK_CREATE(vtkGeometryFilter, geom1);
  geom1->SetInputData(FirstBlock);
  geom1->Update();
  
  VTK_CREATE(vtkPolyDataMapper, mapper1);
  mapper1->SetInputConnection(geom1->GetOutputPort(0));
  mapper1->ScalarVisibilityOn();
  mapper1->SetScalarModeToUsePointFieldData();

  if(varname.size())
    {
    mapper1->SelectColorArray(&varname[varname.find("/") + 1]);
    mapper1->SetLookupTable(lut);
    mapper1->UseLookupTableScalarRangeOn();
    }
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
  }
  reader->Delete();
}
