
all:  tiff  zlib 

psDir ?= $(abspath ..)




zlib:
	ls zlib-1.2.11/build/lib/libzxx.a || (cd zlib-1.2.11 && mkdir -p build)
	ls zlib-1.2.11/build/lib/libzxx.a || (cd zlib-1.2.11/build && cmake -DCMAKE_TOOLCHAIN_FILE=$(psDir)/thirdparty/toolchain$(OPT).txt -DCMAKE_INSTALL_PREFIX=. -DBUILD_SHARED_LIBS=False ..)
	ls zlib-1.2.11/build/lib/libzxx.a || (cd zlib-1.2.11/build && make install)



tiff:
	ls libtiff/build/lib/libtiffxx.a || (cd libtiff && mkdir -p build)
	ls libtiff/build/lib/libtiffxx.a || (cd libtiff/build && cmake -DCMAKE_TOOLCHAIN_FILE=$(psDir)/thirdparty/toolchain$(OPT).txt -DCMAKE_INSTALL_PREFIX=. -DBUILD_SHARED_LIBS=False ..)
	ls libtiff/build/lib/libtiffxx.a || (cd libtiff/build && make install)
#	cp libtiff/build/libtiff/libtiff*.so* $(psDir)/bin/


clean: 
	(rm -rf libtiff/build)
	(rm -rf zlib-1.2.11/build)

