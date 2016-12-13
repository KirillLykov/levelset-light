levelset-light
==============

Levelset-light is an open source C++ library for storing and manipulations on voxel data.

It is used for implementation of the interactions between fluid and deformable bodies from one side and solid walls from other side.
Thus solvers for signed-distance specific equations are not here (they are implemented in openvdb, if needed).
Levelset-light doesn't have it's own data format but supports existing data formats such as vtk, hdf5, vdb.
If you don't want to use all of these formats but only some of them, define only those you are interested in by -DUSE_\<libname\>,
where \<libname\> ::= VTK | HDF5 | OPENVDB. 
VTK and hdf5 files might be visualized by Paraview, while openvdb file might be opened by Houdini with openvdb_houdini plugin.

My typical usage of this code is the following:
1) generate geometry in hdf5 format
2) run particle-based simulation (using different solvers like LAMMPS) in the geomtry described by this Signed distance function: particles which have position with negative sdf are free, others - frozen
This approach was used, in particular, to this simulaiton:

[![](http://lammps.sandia.gov/images/blood_small.png)](http://lammps.sandia.gov/movies/blood.mp4 "RBC simulation")



