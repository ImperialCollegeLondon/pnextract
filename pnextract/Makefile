


all: directories Linux
#all: Win
#all: winVS


ifeq "$(OPT)" "dbg"
 COPTS     = -g -fno-inline  -D_debugCompile_  $(psCFLAGSExtraE)
 BuildDir = ../build/pnextractdbg
else
 COPTS     = -O3  $(psCFLAGSExtraE)
endif


myBinDir ?= ../bin
BuildDir ?= ../build/pnextract

directories:
	@mkdir -p $(myBinDir)
	@mkdir -p $(BuildDir)




USE_ZLIB=1
USE_TIFF=1
USE_VXLIMG=1
psDir ?= $(abspath ..)
include  $(psDir)/Makefile.in

psCXXWin 	?= x86_64-w64-mingw32-g++
psCXX 		?= g++


psCFLAGS := $(psCFLAGS) $(COPTS)   -Wall -pedantic -MMD -MP  -std=c++0x  -I../voxelImage  -I.
psLFLAGS := $(psLFLAGS) -static-libgcc -static-libstdc++  


srcs = \
	blockNet.cpp medialSurf.cpp  nextract.cpp \
	writers_vtk.cpp  writers_vxl.cpp


objs = $(addprefix $(BuildDir)/,$(subst /,--,$(srcs:.cpp=.o)))
deps = $(addprefix $(BuildDir)/,$(subst /,--,$(srcs:.cpp=.d)))

# delete the following if nmake (from visual studio) complains, uncomment for linux compillation
Linux: $(objs)
	$(psCXX)  $(objs) -o $(myBinDir)/pnextract$(OPT) $(psLFLAGS)

-include $(deps)


$(BuildDir)/%.o: $(subst $(BuildDir)/,,$(subst --,/,$(subst .o,.cpp,$@)))
	$(psCXX)  $(psCFLAGS)  $(EXTRACFLAGS) -c $(subst $(BuildDir)/,,$(subst --,/,$(subst .o,.cpp,$@))) -o $@



.PHONY: clean $(deps)

clean:
	$(RM) -f $(objs) $(deps) nelmTo*Par.o



Win:
	$(psCXXWin)  $(psCFLAGS)  $(srcs) -o $(myBinDir)/pnextract$(OPT).exe   $(psLFLAGS)


winVS:
	@cl /nologo /MT /EHsc /Ox  $(srcs)  /link -LARGEADDRESSAWARE /out:..\bin\pnextract$(OPT).exe
	@del *.obj





