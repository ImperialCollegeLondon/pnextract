##  pnextract -- (classical) pore-network extraction

This code extracts a conventional pore network from a micro-CT image.
It is a re-write of the maximal-ball network extraction algorithm by [Dong and Blunt, 2009]. 
It is used as a base for the generalized network extraction code 
[Raeini, Bijeljic and Blunt 2017], sponsored by [TOTAL]. However, here the modules
for extracting the generalized network elements (corners) are not included.   
The network parameters has been calibrated to reproduce single-phase flow 
properties, and water-wet relative permeability for a set of sandstone rocks.
The code uses a new scale-independent shape factor definition to characterize
pores and throats, published in [Bultreys et al. 2018].

----------------------------------------

Note: this is same as [pnflow repository](https://github.com/aliraeini/pnflow) but without the pnflow code.

----------------------------------------

### Instructions

A sample input file, [Image.mhd](../../doc/Image.mhd), is provided
in the [doc](../../doc/) folder, in ascii (text) format. Please use this 
together with a micro-CT image, similar to the [images on our website].

To extract a pore network, decompress the bin.7z file
 and run, in a Windows Command Prompt / PowerShell:     

    PATH\TO\pnextract.exe  Image.mhd

* Replace `PATH\TO\bin\`` with the path to the extracted pnextract.exe.

* To open a command-prompt in Windows, hold the `Shift` key and `right-click`
  into the folder where the input_pnflow.dat is copied and click the `Open Command Window Here` menu.


The instructions are similar in Linux. Additionally, you can source the src/script/bashrc file to set the paths to the compiled binaries:     

     source PATH/TO/src/script/bashrc

See the [pnextract wiki](https://github.com/aliraeini/pnextract/wiki/pnextract-FAQ) for more details.

###  Build instructions

Download and extract [bin.7z](../../bin.7z) for pnextract.exe, 
a Win64 executable compiled using MinGW compilers.

See [../script/README.md](../script/README.md) for build instructions. 
The included [libvoxel](../libvoxel)  is the main dependancy, which itself is (optionally) linked to the included 
[zlib](../../thirdparty/zlib) and [libtiff](../../thirdparty/libtiff) libraries.


###  Licence

The code is release as a free, using a zlib-style licence (TODO: add license).     
The code and executables are provided as is, without any kind of warranty;
use at your own risk.

For contact and further information see [Imperial College - Pore-scale Consortium] website,
raise an `issue` on github, or send me an email:   a.q.raeini@imperial.ac.uk


### References

See the [Publications on our website], in particular:

H. Dong and M. J. Blunt, "Pore-network extraction from micro-computerized-tomography images",  Phys. Rev. E 80, 036307 (2009) 
https://doi.org/10.1103/PhysRevE.80.036307

A Q Raeini, B Bijeljic, and M J Blunt, "Generalized network modeling: Network extraction as a coarse-scale discretization of the void space of porous media", Phys. Rev. E 96, 013312  (2017)
https://doi.org/10.1103/PhysRevE.96.013312

T Bultreys, Q Lin, Y Gao, A Q Raeini, A AlRatrout, B Bijeljic, and M J Blunt . "Validation of model predictions of pore-scale fluid distributions during two-phase flow", Phys. Rev. E 97, 053104 (2018) 
https://link.aps.org/doi/10.1103/PhysRevE.97.053104

[Publications on our website]: https://www.imperial.ac.uk/earth-science/research/research-groups/pore-scale-modelling/publications/
[Imperial College - Pore-scale Consortium]: https://www.imperial.ac.uk/earth-science/research/research-groups/pore-scale-modelling
[images on our website]: https://www.imperial.ac.uk/earth-science/research/research-groups/pore-scale-modelling/micro-ct-images-and-networks/
[Bultreys et al. 2018]: https://link.aps.org/doi/10.1103/PhysRevE.97.053104
[Raeini, Bijeljic and Blunt 2017]: https://doi.org/10.1103/PhysRevE.96.013312
[Dong and Blunt, 2009]: https://doi.org/10.1103/PhysRevE.80.036307
[libtiff]: https://gitlab.com/libtiff/libtiff
[libz]: https://github.com/madler/zlib
[TOTAL]: https://www.total.com
