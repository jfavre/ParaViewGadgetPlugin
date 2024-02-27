# ParaViewGadgetPlugin
a ParaView plugin for GADGET (HDF5) files for ParaView 5.12

# Compilation instructions:
 
cd ParaViewGadgetPlugin
mkdir build; cd build

assuming your new ParaView install is in your PATH

ccmake ..; make

After successful compilation, the directory "./lib64/paraview-5.12/plugins/pvGadgetReader" contains the shared libs

Thus,

export PV_PLUGIN_PATH=<your-path-to>/ParaViewGadgetPlugin/build/lib64/paraview-5.12/plugins/pvGadgetReader

