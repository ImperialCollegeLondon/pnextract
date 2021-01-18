#TODO: implement OPTions: .exe dbg  .so, empty option means Linux executable

srcs := voxelImageConvert.cpp 
tsts := test.py
apps:=$(srcs:.cpp=${OPT})
all: $(apps)

$(warning Apps: ${apps})
$(warning tsts: ${tsts})


msSrc ?= $(abspath ..)
USE_ZLIB=1
USE_TIFF=1
USE_OMP=1
USE_SINGLECPP=1
USE_msTEST=1
include  ${msSrc}/script/Makefile.in

