# Compiling codes

To compile, open a terminal in the upper most directory, where `src` and `thirdparty`
 folders are located, and run:    

 `make -j`

To test the compilation, run:    

 `make test`

Once everything compiled successfully, to clean the temporary files, type:

 `make clean`

The above command can be run inside most of the subfolders, wherever a 
makefile or Makefile is present.  The libraries, those with a `makefile`,
will be compiled before the apps that contain `Makefile`s.

Compilation requires gnu (Linux) make, cmake, a c++ compiler with -std=c++11
support and MPI. The compilation is tested using g++ (version 5+) (default)
and using intel-2018 compilers.

For the modules which have `thirdparty/foamx3m` as dependancy, you if you have
any other OpenFOAM you have on your machine, need to temporarily ***deactivate 
your OpenFOAM installation when compiling*** this code to avoid conflict between
the foam3m provided here and your openfoam instalation.
Additionally, you need to install foam3m dependancies, this can be done in Ubuntu
by running:       
`sudo apt install mpi-default-dev flex libscotch-dev`

# Tests and demos
To test the codes type:

 `make test`

This should copy a series of input files/scripts in a `test` folder and 
run a series of relatively quick test cases (see README.md files in 
subdirectories).  


# Technical notes for code developers

GNU Makefile scripts are used primarily for code compilation and 
running quick tests by code developers.

Automatic tests are written using input files for C++ codes, new C++ 
executables testing internals of the codes and Python scripts. All 
these can be run using `make test` command, which uses 
script/testApp bash script.

All scripts, either for testing or production, which need mathematical 
calculations or plotting and are not performance critical are developed 
using Python. We use Python 3, shich should be available as python3 command.
In Ubuntu (18.04-20.04) this can be installed by typing in a terminal:    
 `sudo apt install python3`


Bash/Shell scripts are used to run Makefile and Python scripts, 
either for testing or in production to simplify the run of openfoam 
solvers.  


The initbash is an independent bash script containing utility macros, 
which together with bashrc (that contain installation variabls), is 
required for running of the compiled applications.

---------

## Aims 

These scripts has been released as a separate module, to help 
the code developers with script re-use,and simplification of 
workflow for code compilation, testing, deployment and release. 


---------

##  Caution

The script here use recursive make by running the AllMake and AllClean 
scripts. The script change, add and delete files on your 
computer: the root directory where files are changed, generated or 
deleted is called msRoot which, by default, points to two directories 
upper to the location of these script themselves.  So, if you ever consider
using these scripts for building your applications, make sure these are 
wrapped inside two (sub-)sub-folders, dedicated to source codes.  Here is 
what the directory structure should looks like:


- `apps/ -------------------- msRoot directory`

    - `src/ ----------------- -- source codes`
        * `script/ ---------- -- ** this module, build & install`
        * `include/ --------- -- ** C++ utility codes`
        * `bench/ ----------- -- ** test/example data`
        * `...`
        * `...`
        
    - `thirdparty/ ---------- -- others' source codes`
        * `foamx3m ---------- -- ** a minified openfoam `
        * `svplot ---------- -- ** a modified former svg_plot`
        * `zlib`
        * `libtiff`
        * `...`

    - `bench/ --------------- -- large input files, too large for src/`
    - `test/ ---------------- -- test folder (auto copied from bench/)`
    - `bin/ ----------------- -- executables folder (auto crea/deleted)`
    - `lib/ ----------------- -- library files (auto crea/deleted)`
    - `share/ --------------- -- configuration files (auto crea/deleted)`
    - `build/ --------------- -- temporary files (auto crea/deleted)`
    - `...`


As shown above, this `script` folder is typically located in `src` 
(which holds regularly changed source codes as subdirectories). The 
less frequently changed source code are placed a directory called 
`thirdparty`.  A bench folder is also sometimes included holding 
temporary/client data, as shown above.  A `.git` directory is kept 
outside the msRoot directory and is used to test and release varius 
modules of the code, which it ends up in the msRoot directory in the 
modules published using git.  The `build`, `bin`,`lib`, `test` and 
`share` folders are auto-generated and should not be used to store 
any user data as they are removed by `make distclean` command.



## Makefile.in usage

 In `Makefile`s, define variables `tsts` together with `USE_msTEST=1` 
 and `srcs` with `USE_msMAKE=1`.  If not defined, `srcs` will be set 
 to all `.cpp` files.  Example Makefiles can be found in libvoxel 
 (single source.cpp apps) or pnextract (multiple source.cpp files).



