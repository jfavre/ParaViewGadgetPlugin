ADD_EXECUTABLE(TestSimpleGadgetReader TestSimpleGadgetReader.cxx)
TARGET_LINK_LIBRARIES(TestSimpleGadgetReader
	PUBLIC GadgetReader
	PRIVATE
	  VTK::FiltersGeometry
          VTK::InteractionStyle
          VTK::RenderingCore
          VTK::RenderingOpenGL2
	)

ADD_EXECUTABLE(TestSimpleReader TestSimpleReader.cxx)
TARGET_LINK_LIBRARIES(TestSimpleReader
	PRIVATE
	  VTK::IOLegacy
	  VTK::FiltersGeometry
          VTK::InteractionStyle
          VTK::RenderingCore
          VTK::RenderingOpenGL2
	)

