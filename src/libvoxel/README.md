##  libvoxel

This library serves as a header-only 3D image manipulation and I/O for my other 
software which work with X-ray computer tomography data.

In addition, voxelImageProcess and voxelToFoam(Par) applications may be included here for convinience; the libvoxel codes (all starting with voxelImage) are independent of these two packages.


The library can read raw data in ascii (.dat) or binary (.raw) formats, in Avizo (.am) formats (only uncompressed and ByteRLE encoded data are supported).  It can also read raw.gz and .tif image formats provided that the  [libz] and [libtiff] libraries are available.

### Usage

This library is used to read 3D image files from other codes, however the standalone app `voxelImageProcess`, which solely acts as an interface to libvoxel, can be used to run basic image processing tasks as well as to print help messages about the keywords supported by libvoxel:

   `voxelImageProcess -h`


### Compilation instructions

See the [../script/README.md](../script/README.md) build instructions.

The code as a library is not currently compiled into binary format. It meant to be header only, but to reduce compile time the library parts are seperated in the voxelImage.cpp file which is compiled and linked statistically to other executables.

The voxelImageProcess utility is pre-compiled into a Win64 executable in pnextract and pnflow codes, in [../../bin.7z](../../bin.7z).

### Download, build and install

###  Licence

The code and executables are provided as is, without any kind of warranty;
use at your own risk.

For further information contact me by email:   a.q.raeini@imperial.ac.uk


### References
See the [Publications on our website], also [Images on our website].

[Publications on our website]: https://www.imperial.ac.uk/earth-science/research/research-groups/pore-scale-modelling/publications/
[Images on our website]: https://www.imperial.ac.uk/earth-science/research/research-groups/pore-scale-modelling/micro-ct-images-and-networks/
[Imperial College - pore-scale consortium]: https://www.imperial.ac.uk/earth-science/research/research-groups/pore-scale-modelling
[libtiff]: https://gitlab.com/libtiff/libtiff
[porefoam]: https://github.com/aliraeini/porefoam
[pnextract]: https://github.com/aliraeini/pnextract
[pnflow]: https://github.com/aliraeini/pnflow
[libz]: https://github.com/madler/zlib
