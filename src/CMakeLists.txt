
paraview_add_plugin(pvGadgetReader
  REQUIRED_ON_SERVER
  VERSION "1.0"
  MODULES GadgetReader
  MODULE_FILES      "${CMAKE_CURRENT_SOURCE_DIR}/Reader/vtk.module")


