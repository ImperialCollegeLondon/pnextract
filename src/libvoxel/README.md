##  libvoxel

This library serves as a header-only 3D image manipulation and I/O for my other 
software which work with X-ray computer tomography data.

In addition, voxelImageConvert and voxelToFoam(Par) applications are included here (when needed) for convinience; the libvoxel codes (all starting with voxelImage) are independent of these apps.

voxelImageConvert is also used to generate synthetic images for testing libvoxel and other packages.

The library can read and write 3D image data in ascii (.dat) or binary (.raw) formats, in Avizo (.am) formats (only uncompressed and ByteRLE encoded data are supported).  It can also read raw.gz and .tif image formats provided that the  [libz] and [libtiff] libraries are available.

### Instructions

#### pre-requisites:

To install necassary libraries in Ubuntu (18.08 etc.), run:    
	`sudo apt install libjpeg-dev  #used by libtiff`    
	`sudo apt install liblzma-dev  #used by libtiff`    

### Download, build and install
This library does not to be compiled (its a header-only template library), and a component of my other `apps`: such as [pnextract]/[pnflow] and [porefoam].  

Nevertheless, the individual .cpp files will becompiled into executables with the same base name when running make from the upper directory, along with other other apps.

This library is not developed here, it is however kept synchronized with a personal development branch.  I will be happy to merge and keep record of any pull-requests you send us though.

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
