
msSrc ?= .
all:
	$(warning "running ${msSrc}/script/AllMake on $(shell dirname $(firstword $(MAKEFILE_LIST)))")
	${msSrc}/script/AllMake


clean: 
	$(warning "running ${msSrc}/script/AllClean on $(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))")
	${msSrc}/script/AllClean
