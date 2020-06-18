/*=========================================================================

  Program:    ParaView Plugin for version 5.8
  Module:     vtkGadgetReader.h
  Written by: Jean M. Favre
              Senior Visualization Software Engineer
              Swiss National Supercomputing Center
              CH-6900 Lugano
  Tested:     Thu 18 Jun 2020 08:58:34 AM CEST
=========================================================================*/

#ifndef vtkGadgetReader_h
#define vtkGadgetReader_h

#include <string>
#include <vector>

#define ALL_TYPES 1
#define OUTPUT_UG 1

#ifdef ALL_TYPES
#include "vtkMultiBlockDataSetAlgorithm.h"
static std::vector<std::string> ParticleTypes = {"PartType0", "PartType1", "PartType2", "PartType3", "PartType4", "PartType5"};
#else
#ifdef OUTPUT_UG
#include "vtkUnstructuredGridAlgorithm.h"
#else
#include "vtkPolyDataAlgorithm.h"
#endif
static std::vector<std::string> ParticleTypes = {"PartType0"};
#endif

#include <map>
#include <sstream>

class vtkDataArraySelection;
class vtkStdString;
class vtkMultiProcessController;

enum  ParticleType : int {Gas=0, Halo=1, Disk=2, Bulge=3, Stars=4, Bndry=5};
enum  CellTypes : int {None=0, Vertex=1, PolyVertex=2};

#ifdef ALL_TYPES
class vtkGadgetReader : public vtkMultiBlockDataSetAlgorithm
#else
#ifdef OUTPUT_UG
class vtkGadgetReader : public vtkUnstructuredGridAlgorithm
#else
class vtkGadgetReader : public vtkPolyDataAlgorithm
#endif
#endif
{
public:
  static vtkGadgetReader *New();
#ifdef ALL_TYPES
  vtkTypeMacro(vtkGadgetReader,vtkMultiBlockDataSetAlgorithm);
#else
#ifdef OUTPUT_UG
  vtkTypeMacro(vtkGadgetReader,vtkUnstructuredGridAlgorithm);
#else
  vtkTypeMacro(vtkGadgetReader,vtkPolyDataAlgorithm);
#endif
#endif
  void PrintSelf(ostream& os, vtkIndent indent);

  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);
  
  void SetDirectoryName(const char* dn);
  vtkGetStringMacro(DirectoryName);

  // Description:
  // When using ParaView, cell generation is recommended, without them
  // many filter operations are unavailable
  // Note that Point Gaussian Rendering does not require cells. 
  // When CellType == Vertex, the reader will generate one vertex cell
  // for each point/particle read.
  // When CellType == PolyVertex, the reader will generate a single cell
  // for each ParticleType
  vtkSetMacro(CellType, int);
  vtkGetMacro(CellType, int);

  int         GetNumberOfPointArrays();
  const char* GetPointArrayName(int index);
  int         GetPointArrayStatus(const char* name);
  void        SetPointArrayStatus(const char* name, int status);
  void        DisableAllPointArrays();
  void        EnableAllPointArrays();
  //
  int         GetNumberOfPointArrayStatusArrays() { return GetNumberOfPointArrays(); }
  const char* GetPointArrayStatusArrayName(int index) { return GetPointArrayName(index); }
  int         GetPointArrayStatusArrayStatus(const char* name) { return GetPointArrayStatus(name); }
  void        SetPointArrayStatusArrayStatus(const char* name, int status) { SetPointArrayStatus(name, status); }

#ifdef PARAVIEW_USE_MPI

    // Description:
    // Set/Get the controller use in compositing (set to
    // the global controller by default)
    // If not using the default, this must be called before any
    // other methods.
    virtual void SetController(vtkMultiProcessController* controller);
    vtkGetObjectMacro(Controller, vtkMultiProcessController);

#endif

protected:
   vtkGadgetReader();
  ~vtkGadgetReader();
  //
  int   RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  int   RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  int   OpenFile();
  void  CloseFile();

  //
  // Internal Variables
  //
  char          *FileName;
  char          *DirectoryName;
  int           CellType;
  vtkTimeStamp  FileModifiedTime;
  vtkTimeStamp  FileOpenedTime;
  int           UpdatePiece;
  int           UpdateNumPieces;
  std::map<std::string, bool> PartTypes; // to decide if to load that particular dataset
  int           NumPart_Total[6];
  int           offsets[7] = {0,0,0,0,0,0, -1};
  int           NumFilesPerSnapshot;
  double        Time;
  typedef std::vector<std::string>  stringlist;
  struct field {
    std::string name;
    int type; // H5T_FLOAT or H5T_INTEGER
    int size; // 1 for scalar, 3 for vector
    };
  typedef std::vector<field> hdf5varlist;
  std::map<std::string, hdf5varlist> FieldArrays;
  stringlist                        GadgetFileNames;
  // To allow paraview gui to enable/disable scalar reading
  vtkDataArraySelection* PointDataArraySelection;

  #ifdef PARAVIEW_USE_MPI
  vtkMultiProcessController* Controller;
  #endif

private:
  vtkGadgetReader(const vtkGadgetReader&)  = delete;
  void operator=(const vtkGadgetReader&)  = delete;
};

#endif
