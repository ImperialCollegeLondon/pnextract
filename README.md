![make and test](https://github.com/aliraeini/pnflow/workflows/make%20and%20test/badge.svg)

##  pnflow - classical pore-network (extraction and) flow simulation

This repository hosts the classical network flow simulation code called pnflow. 
[pnextract network extraction code](https://github.com/aliraeini/pnextract) is also included [here](src/pnm/pnextract) for convenience.


 ----------------------------------------------------------------
 
 **Please see the [src/pnm](src/pnm), [src/pnm/pnextract](src/pnm/pnextract), and [doc](doc) folders for specific details.**
 
 ----------------------------------------------------------------


## General instructions

### Compiling

To compile, open a terminal in the upper most directory and run:

 `make -j`

once everything compiled successfully, to clean the temporary files, type:

 `make clean`

The above command can be run inside most of the subfolders, wherever a 
makefile or Makefile is present.  The libraries, those with a `makefile`,
should be compiled before the apps that contain `Makefile`s.

Compilation requires gnu (Linux) make, cmake, a c++ compiler with -std=c++11
support and an MPI. The compilation is tested using g++ (version 5+) (default)
and using intel-2018 compilers.


### Tests and demos
To test the codes type:

 `make test`

This should copy a series of input files/scripts in a `test` folder and 
run a series of relatively quick test cases (see README.md files in 
subdirectories).  

### Coding conventions

GNU Makefile scripts are used primarily for code compilation and 
running quick tests by code developers.

Automatic tests are written using input files for C++ codes, new C++ 
executables testing internals of the codes and Python scripts. All 
these can be run using `make test` command, which uses 
script/testApp bash script.

All scripts, either for testing or production, which need mathematical 
calculations or plotting and are not performance critical are developed 
using Python.   

All computations which typically take more than few minutes are 
developed using C++.  When interaction with other languages 
(R/Shell/Python) are needed the C++ code shall make use of SiR class to 
achieve this, where possible, using the `_sh_` command to launch other 
scripts/executables.


Bash/Shell scripts are used to launch run Makefile and Python scripts, 
either for testing or in production to simplify the run of openfoam 
solvers.  



### Contact and References ###

For contacts and references please see: 
https://www.imperial.ac.uk/earth-science/research/research-groups/pore-scale-modelling
or contact Ali Q. Raeini, email: a.q.raeini@imperial.ac.uk

More details are given in the apps/doc directory.

