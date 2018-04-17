##  pnextract -- Conventional pore-network extraction code

This code extracts a conventional pore network from a micro-CT image.
The source code has been developed by me, Ali Qaseminejad Raeini, at
[Imperial College - pore-scale consortium].   
It is a re-write of the maximal-ball network extraction algorithm by Dong and Blunt (2009). 
It is used as a base for the generalized network extraction code 
(Raeini, Bijeljic and Blunt 2017).  It does not contain the modules
for extracting the generalized network (corner elements) though.   
The network parameters has been calibrated to reproduce single-phase flow 
properties, and water-wet relative permeability curves of a set of sandstone rocks.


###  Compiling
Already compiled to pnextract.exe, a Win64 executable, using mingw compilers.


### Instructions
A sample input file, Image.mhd, is provided in the doc folder, in ascii 
(text) format. Please use this file together with a 8-bit micro-CT 
image, similar to the [Images on our website].

To extract a pore network, decompress the pnextract.exe.7z and run:   
 pnextract.exe  Image.mhd

###  Licence
This executable should be used for research purposes only. 

The source code will be placed on public domain
(I am just being extra catuous for now).

The code and executables are provided as is, without any kind of warranty;
use at your own risk.

For contact and further information see [Imperial College - pore-scale consortium] website,
or send me an email:   a.qaseminejad-raeini09@imperial.ac.uk


### References
See the [Publications on our website].

[Publications on our website]: http://www.imperial.ac.uk/earth-science/research/research-groups/perm/research/pore-scale-modelling/publications/
[Images on our website]: http://www.imperial.ac.uk/earth-science/research/research-groups/perm/research/pore-scale-modelling/micro-ct-images-and-networks/
[Imperial College - pore-scale consortium]: http://www.imperial.ac.uk/earth-science/research/research-groups/perm/research/pore-scale-modelling


