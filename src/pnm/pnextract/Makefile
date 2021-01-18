
all: msApp 

msSrc ?= $(abspath ../..)


msCFLAGS +=  -I.
msLFLAGS += -static-libgcc -static-libstdc++ 

# keep grouped for auto package
srcs = \
	blockNet.cpp   nextract.cpp medialSurf.cpp \
	writers_vtk.cpp  writers_vxl.cpp 

tsts = test.py





USE_ZLIB=1
USE_TIFF=1
USE_VXIMG=1
USE_OMP=1
USE_msMAKE=1
USE_msTEST=1
include  ${msSrc}/script/Makefile.in



