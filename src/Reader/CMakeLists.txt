set(classes
  vtkGadgetReader)

set(private_headers
  vtkGadgetReader.h)

vtk_module_add_module(GadgetReader
  CLASSES ${classes}
  PRIVATE_HEADERS ${private_headers})

paraview_add_server_manager_xmls(
  XMLS  Gadget_Server.xml)

ADD_EXECUTABLE(TestSimpleGadgetReader TestSimpleGadgetReader.cxx)
TARGET_LINK_LIBRARIES(TestSimpleGadgetReader GadgetReader)

ADD_EXECUTABLE(TestSimpleReader TestSimpleReader.cxx)
TARGET_LINK_LIBRARIES(TestSimpleReader GadgetReader)
