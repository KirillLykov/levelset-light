levelset-light
==============

Levelset-light is an open source C++ library for storing and manipulations on voxel data. The developement of the library on the early stage,
so there it misses many features - some of them a listed in the TODO.txt.
Levelset-light doesn't have it's own data format but supports existing data formats such as vtk, hdf5, vdb, field3D (in plans).
If you don't want to use all of these formats but only some of them, define only those you are interested in by -DUSE_<libname>,
where <libname>::=VTK|HDF5|OPENVDB. VTK and hdf5 files might be visualized by Paraview, while openvdb file might be opened by Houdini with 
openvdb_houdini plugin.