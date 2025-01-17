



![make and test](https://github.com/aliraeini/pnextract/workflows/make%20and%20test/badge.svg)

##  pnextract - pore-network extraction
**pnextract** is an open-source software written in C++ that extracts pore networks from images or statistical reconstructions of porous materials. These pore networks preserve the topology of the pore space and also calculate other parameters, such as the pore size distribution, interfacial area, and volume needed for analysis and simulations. Network extraction is normally the first step for further analysis and modelling, to predict flow and transport processes in porous materials, and to analyse three-dimensional images. This code provides essential network files for pore network flow model ([**pnflow**](https://github.com/aliraeini/pnflow)). 

**Note: this repository is same as [pnflow repository](https://github.com/aliraeini/pnflow) but without the pnflow code.**


 ----------------------------------------------------------------

## See [src/pnm](src/pnm) and [doc](doc) for details on pnextract code.

## See [src/script/README.md](src/script/README.md) for compile/build instructions.

See also README files of other modules which are located in their own directories:    
[src/libvoxel](src/libvoxel), [src/script](src/script) and in [thirdparty](thirdparty).


 ----------------------------------------------------------------

Download the [bin.7z](bin.7z) for pre-compiled Windows executables. 

### Contact and References ###

For contacts and references please see: 
https://www.imperial.ac.uk/earth-science/research/research-groups/pore-scale-modelling


Alternatively, contact Sajjad Foroughi:
- Email: s.foroughi@imperial.ac.uk
- Additional Email: foroughi.sajad@gmail.com

--

For more in-depth understanding of the pnextract code, users can refer to following papers related to this topic, which may provide additional insights:
- [Pore-network extraction from micro-computerized-tomography images](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.80.036307)
- [Generalized network modeling: Network extraction as a coarse-scale discretization of the void space of porous media](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.96.013312) 
