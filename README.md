##  pnextract -- (Conventional) pore-network extraction

This code extracts a conventional pore network from a micro-CT image.
It is a re-write of the maximal-ball network extraction algorithm by [Dong and Blunt, 2009]. 
It is used as a base for the generalized network extraction code 
[Raeini, Bijeljic and Blunt 2017].  It does not contain the modules
for extracting the generalized network (corner elements) though.   
The network parameters has been calibrated to reproduce single-phase flow 
properties, and water-wet relative permeability curves of a set of sandstone rocks.


### Instructions
A sample input file, Image.mhd, is provided in the doc folder, in ascii 
(text) format. Please use this file together with a 8-bit micro-CT 
image, similar to the [Images on our website].

To extract a pore network, decompress the bin/pnextract.exe.7z and run:   
 pnextract.exe  Image.mhd

###  BUild instructions:
Already compiled to bin/pnextract.exe, a Win64 executable, using mingw compilers.

The compilation can be done in Linux by running './AllMake' bash script.

The './AllMakeMinGW' bash script compiles the code for Windows machines.
Run './AllClean' beforhand, to avoid mixing the intermidiate Linux and Windows object files.

###  Dependencies
The included voxelImage library  is the main prerequisite. 
voxelImage itself can optionally be linked to [libz] and [libtiff] to support
reading .raw.gz and 3D .tif files. 
Stripped down versions of both of these libraries are provided in the 
thirdparty directory for compatibility and ease of compilation.

###  Licence

The code and executables are provided as is, without any kind of warranty;
use at your own risk.

For contact and further information see [Imperial College - pore-scale consortium] website,
or send me an email:   a.qaseminejad-raeini09@imperial.ac.uk


### References
See the [Publications on our website].

[Publications on our website]: http://www.imperial.ac.uk/earth-science/research/research-groups/perm/research/pore-scale-modelling/publications/
[Images on our website]: http://www.imperial.ac.uk/earth-science/research/research-groups/perm/research/pore-scale-modelling/micro-ct-images-and-networks/
[Imperial College - pore-scale consortium]: http://www.imperial.ac.uk/earth-science/research/research-groups/perm/research/pore-scale-modelling
[Raeini, Bijeljic and Blunt 2017]: https://doi.org/10.1103/PhysRevE.97.023308
[Dong and Blunt, 2009]: https://doi.org/10.1103/PhysRevE.80.036307
[libtiff]: https://gitlab.com/libtiff/libtiff
[libz]: https://github.com/madler/zlib
