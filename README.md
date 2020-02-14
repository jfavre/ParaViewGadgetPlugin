# ParaViewGadgetPlugin
a ParaView plugin for GADGET (HDF5) files

# Compilation instructions:
 
cd ParaViewGadgetPlugin
mkdir build; cd build

assuming your new ParaView install is in your PATH

ccmake ..
make

After successful compilation, the directory "./lib/paraview-5.8/plugins/pvGadgetReader" contains the shared libs

Thus,

export PV_PLUGIN_PATH=<your-path-to>/ParaViewGadgetPlugin/build/lib/paraview-5.8/plugins/pvGadgetReader

