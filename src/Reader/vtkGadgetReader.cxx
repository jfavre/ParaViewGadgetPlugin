/*=========================================================================

  Program:    ParaView Plugin for version 5.8
  Module:     vtkGadgetReader.cxx
  Written by: Jean M. Favre
              Senior Visualization Software Engineer
              Swiss National Supercomputing Center
              CH-6900 Lugano
  Tested:     Thu 18 Jun 2020 08:58:34 AM CEST
  
  A Data snapshot distributed over multiple files can be read in parallel.
  In this case, use a number of pvservers that divide evently the # of files
  such that we spread the load evenly. For example, 16 files and 4 pvservers.
  
  A Time-series of many single-file snapshots can also be read. Parallel reading
  is not supported in this case.
  
  Compile Instructions:
  1) Compile you own ParaView from source
  2) Assume the source code has been placed in a place of your choice,
     where you should have:
     the asscociated CMakeLists.txt file, and
     the src directory
  3) mkdir build; cd build
  4) cmake .. && make
=========================================================================*/

#include "vtkGadgetReader.h"

#include "vtkCellType.h"
#include "vtkCellArray.h"
#include "vtkDataArray.h"
#include "vtkDataArraySelection.h"
#include "vtkDirectory.h"
#include "vtkFloatArray.h"
#include "vtkIdTypeArray.h"
#include "vtkIdList.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#ifdef ALL_TYPES
#include "vtkMultiBlockDataSet.h"
#endif
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkUnstructuredGrid.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include <vtksys/RegularExpression.hxx>
#include <vtksys/SystemTools.hxx>

#ifdef PARAVIEW_USE_MPI
#include "vtkMultiProcessController.h"
vtkCxxSetObjectMacro(vtkGadgetReader, Controller, vtkMultiProcessController);
#endif

#include <algorithm>
#include <functional>
#include <map>
#include <vtk_hdf5.h>


static int ReadHDF5INT64Dataset(const char *name, hid_t mesh_id, long long *data)
{
  hid_t coords_id, filespace, attr1;
  herr_t  status;
  hsize_t dimsf[2]={0,0};

  coords_id = H5Dopen(mesh_id, name, H5P_DEFAULT);

  if (H5Dread(coords_id, H5T_NATIVE_LONG, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT, data) < 0)
    return -1;

  H5Dclose(coords_id);
  return dimsf[0];
}

static int ReadHDF5Dataset(const char *name, hid_t mesh_id, float *data)
{
  hid_t coords_id, filespace, attr1;
  herr_t  status;
  hsize_t dimsf[2]={0,0};
  //double CGSConversionFactor;
  //float aexpScaleExponent, hScaleExponent;

  coords_id = H5Dopen(mesh_id, name, H5P_DEFAULT);
  
  if(coords_id >= 0)
    {
    filespace = H5Dget_space(coords_id);
    H5Sget_simple_extent_dims(filespace, dimsf, NULL);
    H5Sclose(filespace);
  
    if (H5Dread(coords_id, H5T_NATIVE_FLOAT, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT, data) < 0)
      {
      std::cerr << "Error reading dataset " << name << std::endl;
      return -1;
      }
    }
  else
    std::cerr << "Error opening dataset " << name << std::endl;
/*
// read attributes of given dataset
        attr1 = H5Aopen_name(coords_id, "CGSConversionFactor");
        if (H5Aread(attr1, H5T_NATIVE_DOUBLE, &CGSConversionFactor) < 0)
    return -1;
        status = H5Aclose(attr1);

        attr1 = H5Aopen_name(coords_id, "aexp-scale-exponent");
        if (H5Aread(attr1, H5T_NATIVE_FLOAT, &aexpScaleExponent) < 0)
    return -1;
        status = H5Aclose(attr1);

        attr1 = H5Aopen_name(coords_id, "h-scale-exponent");
        if (H5Aread(attr1, H5T_NATIVE_FLOAT, &hScaleExponent) < 0)
    return -1;
        status = H5Aclose(attr1);
*/
// end of attributes read
// TODO Do coordinates have to be scaled with given attributes. Don't know yet.
    H5Dclose(coords_id);
    return dimsf[0];
}

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkGadgetReader);
//----------------------------------------------------------------------------
vtkGadgetReader::vtkGadgetReader()
{
  this->SetNumberOfInputPorts(0);
  this->Time = 0.0;
  this->DirectoryName            = nullptr;
  this->FileName                 = nullptr;
  this->NumFilesPerSnapshot      = 0;
  this->CellType                 = 0;
  this->UpdatePiece              = 0;
  this->UpdateNumPieces          = 0;
  this->PointDataArraySelection  = vtkDataArraySelection::New();
  for (int i=0; i< 6; i++)
    this->NumPart_Total[i] = 0;
  for(auto const &it: ParticleTypes)
    PartTypes[it] = false; // do not construct the dataset to add to the multi-block container

#ifdef PARAVIEW_USE_MPI
  this->Controller = nullptr;
  this->SetController(vtkMultiProcessController::GetGlobalController());
#endif
}
//----------------------------------------------------------------------------
vtkGadgetReader::~vtkGadgetReader()
{
  this->CloseFile();

  if(this->DirectoryName)
    free(this->DirectoryName);
  if(this->FileName)
    free(this->FileName);
  this->PointDataArraySelection->Delete();
  this->PointDataArraySelection = nullptr;

#ifdef PARAVIEW_USE_MPI
  this->SetController(NULL);
#endif
}
//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
void vtkGadgetReader::CloseFile()
{

}
//----------------------------------------------------------------------------
int vtkGadgetReader::OpenFile()
{
  if (!this->DirectoryName)
    {
    vtkErrorMacro(<<"DirectoryName must be specified.");
    return 0;
    }

  if (FileModifiedTime>FileOpenedTime)
    {
    this->CloseFile();
    }

  return 1;
}

void vtkGadgetReader::SetDirectoryName(const char* dn)
{
  this->DirectoryName = strdup(vtksys::SystemTools::GetParentDirectory(dn).c_str());
  this->FileName = strdup(dn);
}

//----------------------------------------------------------------------------
int vtkGadgetReader::RequestInformation(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  outInfo->Set(CAN_HANDLE_PIECE_REQUEST(), 1);

#ifdef PARAVIEW_USE_MPI
  if (this->Controller)
    {
    this->UpdatePiece = this->Controller->GetLocalProcessId();
    this->UpdateNumPieces = this->Controller->GetNumberOfProcesses();
    }
#else
  this->UpdatePiece = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
#endif

    vtkDirectory* dir = vtkDirectory::New();
    int opened = dir->Open(this->DirectoryName);
    if (!opened)
    {
      vtkErrorMacro("Couldn't open " << this->DirectoryName);
      dir->Delete();
      return 0;
    }
    vtkIdType numFiles = dir->GetNumberOfFiles();
    
    this->GadgetFileNames.clear();
    // if opening a directory where a snapshot has been split into multiple files,
    // then the GadgetFileNames should be filled up with all sub-files
    // otherwise we are in the presence of a "normal directory" where
    // multiple timesteps have been stored. In that case, we should ignore
    // the other files found in the directory
 
// open only file0 and look what's inside
  hid_t    file_id = H5Fopen(this->FileName, H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t    root_id, mesh_id, attr1, d_id, filespace;
  herr_t  status;
  root_id = H5Gopen(file_id, "/", H5P_DEFAULT);
  int j=0;

  H5E_auto_t func;
  void *client_data;
  H5Eget_auto2(H5E_DEFAULT, &func, &client_data);
  
  int J; // how many field variable are we finding for the given PartType
  for (auto const i: ParticleTypes) // {"PartType0", "PartType1", "PartType2", "PartType3", "PartType4", "PartType5"}
    {
    J=0;
    if(H5Lexists(root_id, i.c_str(), H5P_DEFAULT))
      {
      mesh_id = H5Gopen(root_id, i.c_str(), H5P_DEFAULT);
      H5G_info_t group_info;
      herr_t err = H5Gget_info(mesh_id, &group_info);
      //std::cerr << __LINE__ << "opening group = " << i.c_str() << " has " << group_info.nlinks << " datasets" << std::endl;
      hdf5varlist group_of_datasets;
      for(hsize_t idx=0; idx < group_info.nlinks; idx++)
        {
        char name[64];
        H5Gget_objname_by_idx(mesh_id, idx, name, 64 );
        if (strcmp(name,"Coordinates"))
          {
          d_id = H5Dopen(mesh_id, name, H5P_DEFAULT);
          if(d_id >= 0)
            {
            group_of_datasets.push_back(field()); // default constructor
        
            hsize_t dimsf[2]={0,0};
            hid_t mytype = H5Dget_type(d_id);
            filespace = H5Dget_space(d_id);
            H5Sget_simple_extent_dims(filespace, dimsf, NULL);
            H5Sclose(filespace);
        
            group_of_datasets[J].type = H5Tget_class(mytype);
            group_of_datasets[J].size = dimsf[1] ? dimsf[1] : 1;
            group_of_datasets[J].name = name;
            H5Tclose(mytype);
            H5Dclose(d_id);
       
            this->PointDataArraySelection->AddArray((i + "/" + name).c_str());
            //std::cerr << __LINE__ << ": " << J << ": adding " << i + "/" + name << endl;
            J++;
            }
          }
        }
      this->FieldArrays[i] = group_of_datasets;

      PartTypes[i] = true;
      H5Gclose(mesh_id);
      }
    j++;
    this->offsets[j] = this->offsets[j-1] + J; // we will deduct that offset later when searching for the allocated arrays
    }
  H5Gclose(root_id);
  
  root_id = H5Gopen(file_id, "/Header", H5P_DEFAULT);

  attr1 = H5Aopen_name(root_id, "NumPart_Total");
  if (H5Aread(attr1, H5T_NATIVE_INT, &this->NumPart_Total) < 0)
    {
    vtkErrorMacro( << "cannot find the NumPart_Total");
    }
  else
     status = H5Aclose(attr1);
  
  attr1 = H5Aopen_name(root_id, "NumFilesPerSnapshot");
  if (H5Aread(attr1, H5T_NATIVE_INT, &this->NumFilesPerSnapshot) < 0)
    {
    vtkErrorMacro( << "cannot find the NumFilesPerSnapshot");
    }
  else
    status = H5Aclose(attr1);
// turn off errors for the moment since we don't have data with Time
  H5Eset_auto2(H5E_DEFAULT, NULL, NULL);
  
  attr1 = H5Aopen_name(root_id, "Time");
  if (H5Aread(attr1, H5T_NATIVE_DOUBLE, &this->Time) < 0)
    {
    vtkErrorMacro( << "cannot read attribute Time");
    }
  else
    status = H5Aclose(attr1);
  H5Eset_auto2(H5E_DEFAULT, func, client_data);
    
  H5Gclose(root_id);
  H5Fclose(file_id);

  double timeRange[2] = {this->Time, this->Time};

  outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), &this->Time, 1);
  outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(), timeRange, 2);
  
  if(this->NumFilesPerSnapshot > 1)
    for (vtkIdType i = 0; i < numFiles; i++)
    {
      if (strcmp(dir->GetFile(i), ".") == 0 ||
          strcmp(dir->GetFile(i), "..") == 0)
      {
        continue;
      }

      std::string fileString = this->DirectoryName;
      fileString += "/";
      fileString += dir->GetFile(i);

      if(fileString.find("snap") != std::string::npos)
        {
        this->GadgetFileNames.push_back(fileString);
        }
    }
  else
    this->GadgetFileNames.push_back(this->FileName);
  
  dir->Delete();
  return 1;
}

int vtkGadgetReader::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkDataObject* doOutput = outInfo->Get(vtkDataObject::DATA_OBJECT());
  
  int length = outInfo->Length( vtkStreamingDemandDrivenPipeline::TIME_STEPS() );
  double* steps = outInfo->Get( vtkStreamingDemandDrivenPipeline::TIME_STEPS() );
  std::cerr << "Got Time = " << steps[0] << std::endl;
  doOutput->GetInformation()->Set( vtkDataObject::DATA_TIME_STEP(), steps[0] );
  
#ifdef ALL_TYPES
  vtkMultiBlockDataSet* mb = vtkMultiBlockDataSet::SafeDownCast(doOutput);
  if (!mb)
    {
    return 0;
    }
#else
#ifdef OUTPUT_UG
  vtkUnstructuredGrid* output = vtkUnstructuredGrid::SafeDownCast(doOutput);
#else
  vtkPolyData* output = vtkPolyData::SafeDownCast(doOutput);
#endif
#endif

#ifdef PARAVIEW_USE_MPI
  if (this->Controller &&
      (this->UpdatePiece != this->Controller->GetLocalProcessId() ||
       this->UpdateNumPieces != this->Controller->GetNumberOfProcesses()))
  {
    vtkDebugMacro(<< "Parallel failure, Id's not right (ignore)");
    this->UpdatePiece = this->Controller->GetLocalProcessId();
    this->UpdateNumPieces = this->Controller->GetNumberOfProcesses();
  }
#else
  this->UpdatePiece = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  this->UpdateNumPieces = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());
#endif

#define PARALLEL_DEBUG 1
#ifdef PARALLEL_DEBUG
  std::ostringstream fname;
  fname << "/tmp/out." << this->UpdatePiece << ".txt" << ends;
  ofstream errs;
  errs.open(fname.str().c_str(),ios::trunc);
  errs << "piece " << this->UpdatePiece << " out of " << this->UpdateNumPieces << endl;

#endif
    
  int load, MyNumber_of_Files;
  int nb_of_Files = this->GadgetFileNames.size();
  if(this->UpdateNumPieces == 1)
    {
    load = MyNumber_of_Files = nb_of_Files;
    }
  else
    {
    load = nb_of_Files / this->UpdateNumPieces;
    if (this->UpdatePiece < (this->UpdateNumPieces-1))
      {
      MyNumber_of_Files = load;
      }
    else
      {
      MyNumber_of_Files = nb_of_Files - (this->UpdateNumPieces-1) * load;
      }
    }

  hid_t   file_id, root_id, mesh_id, coords_id, dataset_id, filespace, attr1;
  herr_t  status;
  hsize_t dimsf[2];
  int LoadPart_Total[6] ={0,0,0,0,0,0};
  int lpT[6];
  for(int myFile = load*this->UpdatePiece; myFile< load*this->UpdatePiece + MyNumber_of_Files; myFile++)
    {
    errs << __LINE__ << ": opening " << this->GadgetFileNames[myFile] << endl;
    file_id = H5Fopen(this->GadgetFileNames[myFile].c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    root_id = H5Gopen(file_id, "/Header", H5P_DEFAULT);

    attr1 = H5Aopen_name(root_id, "NumPart_ThisFile");
    if (H5Aread(attr1, H5T_NATIVE_INT, &lpT) < 0)
      vtkErrorMacro( << "cannot find the NumPart_ThisFile");

    for (int i=0; i< 6; i++)
      LoadPart_Total[i] += lpT[i];

    status = H5Aclose(attr1);

    H5Gclose(root_id);
    H5Fclose(file_id);
    }

#ifdef ALL_TYPES
#ifdef OUTPUT_UG
  vtkUnstructuredGrid *output;
#else
  vtkPolyData *output;
#endif
#endif
// for each particle type, we now double-check if any variable has been selected for reading
// if none were selected, we completely skip constructing the dataset
// this could be made optional with a GUI checkbox.
  for (auto const &it: this->FieldArrays)
    {
      bool ReallyLoad = false;
      for(auto it2 = it.second.begin(); it2 != it.second.end(); ++it2)
        if(this->GetPointArrayStatus((it.first + std::string("/") + (*it2).name).c_str() )) // if variable name enabled to be read
          {
          ReallyLoad = true;
          }
      PartTypes[it.first] = ReallyLoad;
    }
  
  vtkFloatArray  *data;
  vtkIdTypeArray *uidata;
  int validPart=0; // index into the MultiBlock container
  int myType=Gas;
  for(auto const &it: PartTypes)
    {
    if(it.second)
      {
#ifdef PARALLEL_DEBUG
      errs << __LINE__ << ": creating PolyData for PartType " << myType << " with " << LoadPart_Total[myType] << " points\n";
#endif

#ifdef ALL_TYPES
#ifdef OUTPUT_UG
      output = vtkUnstructuredGrid::New();
#else
      output = vtkPolyData::New();
#endif

      mb->SetBlock(validPart, output);
      mb->GetMetaData(validPart)->Set(vtkCompositeDataSet::NAME(), it.first);
      output->Delete();
#endif
      vtkFloatArray *coords = vtkFloatArray::New();
      coords->SetNumberOfComponents(3);
      coords->SetNumberOfTuples(LoadPart_Total[myType]);
      coords->SetName("coords");

      vtkPoints *points = vtkPoints::New();
      points->SetData(coords);
      output->SetPoints(points);
      coords->Delete();
      points->Delete();
// full varnames have been added in alphabetical order, first all related to PartType0, then PartType1, etc...
      for(int i = 0 ; i < this->GetNumberOfPointArrays(); i++)
        {
        if(this->GetPointArrayStatus(this->GetPointArrayName(i))) // if variable name enabled to be read
          {
          if(!strncmp(this->GetPointArrayName(i), it.first.c_str(), 9)) // first 9 characters are equal
          {
          const char *name = &this->GetPointArrayName(i)[10];
#ifdef PARALLEL_DEBUG
      errs << __LINE__ << ": " << this->GetPointArrayName(i) << ": creating data array \"" << name << "\" for PartType " << myType << " with " << LoadPart_Total[myType] << " points\n";
#endif
          int parti = i - offsets[myType];
 std::cerr << " creating data array \"" << this->FieldArrays[it.first][parti].name << "\"(" << myType << "),"
            << this->FieldArrays[it.first][parti].type << ","
            << this->FieldArrays[it.first][parti].size << ", "
            << "\n";

          if(this->FieldArrays[it.first][parti].type == H5T_FLOAT)
            {
            data = vtkFloatArray::New();
            data->SetNumberOfComponents(this->FieldArrays[it.first][parti].size);
            data->SetNumberOfTuples(LoadPart_Total[myType]);
            data->SetName(name); // use the mapped names; HDF5 reading will need the original
            output->GetPointData()->AddArray(data);
            data->Delete();
            }
          else if(this->FieldArrays[it.first][parti].type == H5T_INTEGER)
            {
            uidata = vtkIdTypeArray::New();
            uidata->SetNumberOfComponents(1);
            uidata->SetNumberOfTuples(LoadPart_Total[myType]);
            uidata->SetName(name);
            if(!strcmp(name, "ParticleIDs"))
              {
              output->GetPointData()->SetActiveGlobalIds(name);
              std::cout << " Setting Active GlobalIds = " << name << "\n";
              }
            output->GetPointData()->AddArray(uidata);
            uidata->Delete();
            }
          }
          }
        }

      if (this->CellType == CellTypes::Vertex)
        {
        vtkCellArray *vertices =  vtkCellArray::New();
#ifdef VTK_CELL_ARRAY_V2
        vtkIdTypeArray* ca = vtkIdTypeArray::New();
        ca->SetNumberOfValues(2*LoadPart_Total[myType]);
        vertices->AllocateExact(LoadPart_Total[myType], LoadPart_Total[myType]);
        vtkIdType* cells = ca->GetPointer(0);
        for (auto id = 0; id < LoadPart_Total[myType]; id++)
          {
          *cells++ = 1;  // number of vertices in cell
          *cells++ = id; // internal ID of vertex
          }
        vertices->ImportLegacyFormat(ca);
#else
        vtkIdType* cells = vertices->WritePointer(LoadPart_Total[myType], 2 * LoadPart_Total[myType]);

        for (vtkIdType i = 0; i < LoadPart_Total[myType]; ++i)
          {
          cells[2 * i] = 1;
          cells[2 * i + 1] = i;
          }
#endif
       output->SetCells(VTK_VERTEX, vertices);
       vertices->Delete();
       }
     else if (this->CellType == CellTypes::PolyVertex)
        {
        vtkIdList *list = vtkIdList::New();
        list->SetNumberOfIds(LoadPart_Total[myType]);
        for(unsigned int i=0; i < LoadPart_Total[myType]; i++)
          list->SetId(i, i);
        output->Allocate(1);
        output->InsertNextCell(VTK_POLY_VERTEX, list);
        list->Delete();
        }
      validPart++;
      }
       std::cerr <<"================\n";
      myType++;
    }

  int offset, fileOffsetNodes[6] = {0,0,0,0,0,0};

  for(int myFile = load*this->UpdatePiece; myFile< load*this->UpdatePiece + MyNumber_of_Files; myFile++)
    {
    file_id = H5Fopen(this->GadgetFileNames[myFile].c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    root_id = H5Gopen(file_id, "/", H5P_DEFAULT);
    validPart=0;
    myType=Gas;
    for(auto const &it: PartTypes)
      {
      if(it.second)
        {
        char name[64];
        sprintf(name,"PartType%1d", myType);
#ifdef ALL_TYPES
#ifdef OUTPUT_UG
        output = static_cast<vtkUnstructuredGrid*>(mb->GetBlock(validPart));
#else
        output = static_cast<vtkPolyData*>(mb->GetBlock(validPart));
#endif
#endif
        mesh_id = H5Gopen(root_id, name, H5P_DEFAULT);

// insert coordinates read
        float *dptr = static_cast<vtkFloatArray *>(output->GetPoints()->GetData())->GetPointer(fileOffsetNodes[myType]*3);
        offset = ReadHDF5Dataset("Coordinates", mesh_id, dptr);
// end of coordinates read

// insert PointData here
      for(int i = 0 ; i < this->GetNumberOfPointArrays(); i++)
        {
        if(this->GetPointArrayStatus(this->GetPointArrayName(i)))
          {
          if(!strncmp(this->GetPointArrayName(i), it.first.c_str(), 9)) // first 9 characters are equal
            {
            const char *name = &this->GetPointArrayName(i)[10];
#ifdef PARALLEL_DEBUG
      errs << this->GadgetFileNames[myFile] << ": reading data array " << name << " for PartType " << myType << " with " << this->NumPart_Total[myType] << " points\n";
#endif
          if(!strcmp(name, "ParticleIDs") || !strcmp(name, "ParticleChildIDsNumber") || !strcmp(name, "ParticleIDGenerationNumber"))
            {
            uidata = static_cast<vtkIdTypeArray *>(output->GetPointData()->GetArray(name));
            vtkTypeInt64 *lptr = uidata->GetPointer( uidata->GetNumberOfComponents() * fileOffsetNodes[myType] );
            ReadHDF5INT64Dataset(name, mesh_id, lptr);
            }
          else
            {
            data = static_cast<vtkFloatArray *>(output->GetPointData()->GetArray(name));
            dptr = data->GetPointer( data->GetNumberOfComponents() * fileOffsetNodes[myType] );
            ReadHDF5Dataset(name, mesh_id, dptr);
            // split by component
            if(!strcmp(name, "velocity"))
              {
              double tuple[3];
              int NbTuples = data->GetNumberOfTuples();
              vtkFloatArray *vx = vtkFloatArray::New();
              vx->SetName("vx");
              vx->SetNumberOfTuples(NbTuples);
              output->GetPointData()->AddArray(vx);
              vx->Delete();

              vtkFloatArray *vy = vtkFloatArray::New();
              vy->SetName("vy");
              vy->SetNumberOfTuples(NbTuples);
              output->GetPointData()->AddArray(vy);
              vy->Delete();

              vtkFloatArray *vz = vtkFloatArray::New();
              vz->SetName("vz");
              vz->SetNumberOfTuples(NbTuples);
              output->GetPointData()->AddArray(vz);
              vz->Delete();
              for(vtkIdType i=0; i < NbTuples; i++)
                {
                data->GetTuple(i, tuple);
                vx->SetTuple1(i, tuple[0]);
                vy->SetTuple1(i, tuple[1]);
                vz->SetTuple1(i, tuple[2]);
                }
              }
            }
            }
          }
        }
// end of PointData read

        H5Gclose(mesh_id);
        validPart++;
        }
      fileOffsetNodes[myType] += offset;
      myType++;
    } // for all part types
    H5Gclose(root_id);
    H5Fclose(file_id);
    } // for all files in my load
  this->CloseFile();

#ifdef PARALLEL_DEBUG
  errs.close();
#endif
  return 1;
}


//----------------------------------------------------------------------------
const char* vtkGadgetReader::GetPointArrayName(int index)
{
  return this->PointDataArraySelection->GetArrayName(index);
}
//----------------------------------------------------------------------------
int vtkGadgetReader::GetPointArrayStatus(const char* name)
{
  return this->PointDataArraySelection->ArrayIsEnabled(name);
}
//----------------------------------------------------------------------------
void vtkGadgetReader::SetPointArrayStatus(const char* name, int status)
{
  if (status!=this->GetPointArrayStatus(name))
    {
    if (status)
      {
      this->PointDataArraySelection->EnableArray(name);
      }
    else
      {
      this->PointDataArraySelection->DisableArray(name);
      }
    this->Modified();
    }
}

void vtkGadgetReader::EnableAllPointArrays()
{
    this->PointDataArraySelection->EnableAllArrays();
}

//----------------------------------------------------------------------------
void vtkGadgetReader::DisableAllPointArrays()
{
    this->PointDataArraySelection->DisableAllArrays();
}

//----------------------------------------------------------------------------
int vtkGadgetReader::GetNumberOfPointArrays()
{
  return this->PointDataArraySelection->GetNumberOfArrays();
}
//----------------------------------------------------------------------------
void vtkGadgetReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "Directory Name: " <<
    (this->DirectoryName ? this->DirectoryName : "(none)") << "\n";
}
