set(classes
  vtkGadgetReader)

set(private_headers
  vtkGadgetReader.h)

vtk_module_add_module(GadgetReader
  CLASSES ${classes}
  PRIVATE_HEADERS ${private_headers})

paraview_add_server_manager_xmls(
  XMLS  readers_iogadget.xml)
