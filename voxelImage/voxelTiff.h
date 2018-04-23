
#ifndef VOXELTIFF_H
#define VOXELTIFF_H


#include <stdio.h>
#include <stdlib.h>
#include <typeinfo>
#include <typeindex>
#include <unordered_map>
#include <string.h>
#include <iostream>

#ifdef HAVE_UNISTD_H
# include <unistd.h>
#endif



#include "tiffio.h"
#include "voxelImage.h"


//#define	streq(a,b)	(strcmp(a,b) == 0)
//#define	strneq(a,b,n)	(strncmp(a,b,n) == 0)





//https://stackoverflow.com/questions/24059421/adding-custom-tags-to-a-tiff-file, didn't work properly, also a too dirty solution

template<typename T>
int tifDataType(T)
{
   const std::unordered_map<std::type_index, int> tifTypes
	{
    {std::type_index(typeid(unsigned char)), SAMPLEFORMAT_UINT},
    {std::type_index(typeid(char)),          SAMPLEFORMAT_INT},
    {std::type_index(typeid(int)),           SAMPLEFORMAT_INT},
    {std::type_index(typeid(unsigned int)),  SAMPLEFORMAT_UINT},  
    {std::type_index(typeid(short)),         SAMPLEFORMAT_INT},
    {std::type_index(typeid(unsigned short)), SAMPLEFORMAT_UINT},
    {std::type_index(typeid(float)),         SAMPLEFORMAT_IEEEFP},
    {std::type_index(typeid(double)),        SAMPLEFORMAT_IEEEFP},
	};
	return tifTypes.at(std::type_index(typeid(T)));
}



inline void getTifTags(vec3& X0_, vec3& dx_, TIFF *tif)
{
	float dx,dy,dz;
	TIFFGetField(tif, TIFFTAG_XRESOLUTION, &dx);
	TIFFGetField(tif, TIFFTAG_YRESOLUTION, &dy);
	dz=dy;
	dx_={dx,dy,dz};

	float X0=0.0, Y0=0.0, Z0=-2.1e30; 
	TIFFGetField(tif, TIFFTAG_XPOSITION, &X0);
	TIFFGetField(tif, TIFFTAG_YPOSITION, &Y0);

	uint32 kFrst=0;
	TIFFGetField(tif, TIFFTAG_IMAGEDEPTH, &kFrst);
	Z0=double(kFrst)*dx;

	X0_[0]=X0;	X0_[1]=Y0;	X0_[2]=Z0;


	char * info;
	if(TIFFGetField(tif, TIFFTAG_IMAGEDESCRIPTION, &info))
	{
		std::istringstream instrim;
		//std::cout<<"\" "<<info<<" \""<<std::endl;
		instrim.str(info);
		while(instrim.good())
		{	std::string str;
			instrim>>str;
			if(str=="dx")	instrim >> dx_;
			else if(str=="X0")	instrim >> X0_;
		}
	}
	if(dx_[0]<1.0e-16) {dx_={1.0e-6, 1.0e-6, 1.0e-6}; std::cout<<"\n\n\tError: dx read from tif seems invalid, assuming dx=1.0e-6\n"<<std::endl;}

}

inline void setTifTags(const vec3& X0_, const vec3& dx_, TIFF *tif)
{ // negative numbers not supported in TIFFTAG_XPOSITION :-( for X0, using TIFFTAG_IMAGEDESCRIPTION instead
	float dx=dx_[0], dy=dx_[1];
	TIFFSetField(tif, TIFFTAG_XRESOLUTION, dx);
	TIFFSetField(tif, TIFFTAG_YRESOLUTION, dy);
	//TIFFSetField(tif, TIFFTAG_ZRESOLUTION, dz);, dz=dx_[2]

	
	std::ostringstream ostrim;
	ostrim<<"\t dx "<<dx_;
	if(mag(X0_)>1.0e-16) ostrim<<"\t X0 "<<X0_;

	const char* info = (ostrim.str()+" \0").c_str();//+" \0"
	std::cout<<"\" "<<info<<" \""<<std::endl;
	TIFFSetField(tif, TIFFTAG_IMAGEDESCRIPTION, info);

	TIFFSetField(tif, TIFFTAG_SOFTWARE, "voxelImage+libtiff");
}


template<typename Type>   int readTif(voxelField<Type>&  aa, std::string innam )
{
	uint32 nx, ny;
	//uint16 samplesperpixel;
	//uint16 bitspersample=1;
	//uint16 config;
	//uint16 photometric;
	//uint16* red, green, blue;
	//tsize_t rowsize;
	//register uint32 row;
	//register tsample_t s;
        

	TIFF *tif = (TIFF *) NULL;
	tif = TIFFOpen(innam.c_str(), "r");	 if (tif == NULL)	return (-1);

	voxelImageT<Type>* vxls = dynamic_cast<voxelImageT<Type>*>(&aa);
	if(vxls) 		getTifTags(vxls->X0Ch(),vxls->dxCh(),tif);

	TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &nx);
	TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &ny);
	//TIFFGetField(tif, TIFFTAG_PLANARCONFIG, &config);
	int npages=TIFFNumberOfDirectories(tif);

	aa.reset(nx,ny,npages,0);

	std::cout<<"size:"<<aa.size3()<<" * "<<sizeof(Type)
	<<"  X0:"<<vxls->X0()<<"  dx:"<<vxls->dx()<<std::endl;;


	for(int pn=0;pn<npages;++pn) {

		TIFFReadEncodedStrip( tif, static_cast<tstrip_t>(0), static_cast<void *>(&aa(0,0,pn)), static_cast<tsize_t>(nx * ny) * sizeof(Type) );
		//for (row = 0; row < ny; row++) {		if (TIFFReadScanline(tif, &aa(0,row,pn), row, 0) < 0)		break;	}
		TIFFReadDirectory(tif);
	}
  	TIFFClose(tif);


	return (0);
}

template<typename Type>   int writeTif(const voxelField<Type>&  aa, std::string outnam )
{

	//uint32 rowsperstrip = uint32(-1);
	uint32 nx=aa.size3()[0], ny=aa.size3()[1];
	int pn=0, npages=aa.size3()[2];
	int smplfrmt	 = tifDataType(Type());
	//uint16 samplesperpixel;
	//uint16 bitspersample;
	//uint16 config;
	//uint16 photometric;
	//uint16* red, green, blue;
	//tsize_t rowsize;
	//register uint32 row;
	//register tsample_t s;

	TIFF * tif = (TIFF *) NULL;
	tif = TIFFOpen(outnam.c_str(), "w");		if (tif == NULL)	return (-2);

	const voxelImageT<Type>* vxls = dynamic_cast<const voxelImageT<Type>*>(&aa);
	if(vxls) setTifTags(vxls->X0(),vxls->dx(),tif);
	else std::cout<<"dxXo not set"<<std::endl;


	//std::cout<<aa.size3()<<std::endl;
   for(pn=0;pn<npages;++pn) {
		//TIFFSetField(tif, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);


		TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, nx);
		TIFFSetField(tif, TIFFTAG_IMAGELENGTH, ny);

		TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 8*sizeof(Type));
		TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);
		//TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
		//cpTags(tif, tif);
		TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, smplfrmt);
		TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK); ///. this is needed by stupid paraview!


		if(sizeof(Type)==1) TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_LZW); /// For best compression, choose COMPRESSION_PACKBITS and manually compress as 7z. Other opts:COMPRESSION_NONE  COMPRESSION_DEFLATE COMPRESSION_PACKBITS
		else TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_LZW);

		//TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, TIFFDefaultStripSize(tif, ny));
		TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, ny);


		TIFFSetField(tif, TIFFTAG_PAGENUMBER, pn, npages);

		Type aOrig=aa(nx/2,ny/2,pn);  //! Warn if image modified in libtiff
		TIFFWriteEncodedStrip( tif, 0, const_cast<Type *>(&aa(0,0,pn)), nx * ny * sizeof(Type) ); 
		if(aOrig!=aa(nx/2,ny/2,pn)) std::cout<<"Warning image modified in libtiff"<<std::endl;


		//for (int row = 0; row<ny; row++) {	if (TIFFWriteScanline(tif, &aa(0,row,pn), row, 0) < 0)	break;  }
		TIFFWriteDirectory(tif);
  }


	TIFFClose(tif);
	return (0);
}


template<typename Type>   int writeTif(const voxelField<Type>&  aa, std::string outnam, int iStart,int iEnd , int jStart,int jEnd , int kStart,int kEnd )
{
	uint32 rowsperstrip = uint32(-1);
	TIFF *tif;
	uint32 nx=iEnd-iStart, ny=jEnd-jStart;
	int pn=0, npages=aa.size3()[2];
	int smplfrmt	 = tifDataType(Type());
	//uint16 samplesperpixel;
	//uint16 bitspersample;
	//uint16 config;
	//uint16 photometric;
	//uint16* red, green, blue;
	//tsize_t rowsize;
	register uint32 row;
	//register tsample_t s;
        
  tif = (TIFF *) NULL;

	
	tif = TIFFOpen(outnam.c_str(), "w");		if (tif == NULL)	return (-2);
	
	const voxelImageT<Type>* vxls = dynamic_cast<const voxelImageT<Type>*>(&aa);
	if(vxls) setTifTags(vxls->X0(),vxls->dx(),tif);


	//std::cout<<aa.size3()<<std::endl;
   for(pn=0;pn<npages;++pn) {
		TIFFSetField(tif, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);

		TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, nx);
		TIFFSetField(tif, TIFFTAG_IMAGELENGTH, ny);

		TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 8*sizeof(Type));
		TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);
		//TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);

		TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK); ///. this is needed by stupid paraview!
		TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, smplfrmt);

		if(sizeof(Type)==1) TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_LZW); /// For best compression, choose COMPRESSION_PACKBITS and manually compress as 7z. Other opts:COMPRESSION_NONE  COMPRESSION_DEFLATE COMPRESSION_PACKBITS
		else TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_LZW);


		TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, TIFFDefaultStripSize(tif, rowsperstrip));
		//TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, ny);


		TIFFSetField(tif, TIFFTAG_PAGENUMBER, pn, npages);

		Type aOrig=aa(iStart+nx/2,jStart+ny/2,pn);  //! Warn if image modified in libtiff
		for (row = jStart; row < (uint32)jEnd; row++) {
			if (TIFFWriteScanline(tif, const_cast<Type *>(&aa(iStart,row,kStart+pn)), row, 0) < 0)
				break;
		}
		if(aOrig!=aa(iStart+nx/2,jStart+ny/2,pn)) std::cout<<"Warning image modified in libtiff"<<std::endl;

		TIFFWriteDirectory(tif);
  }

	TIFFClose(tif);
	return (0);
}




inline std::unique_ptr<voxelImageTBase> readTif
(
	std::string innam
)
{
	TIFF *tif = (TIFF *) NULL;
	tif = TIFFOpen(innam.c_str(), "r");	 if (tif == NULL)	return std::unique_ptr<voxelImageTBase>();
 //TIFFTAG_SAMPLEFORMAT, smplfrmt);TIFFTAG_BITSPERSAMPLE
	uint16 frmt=SAMPLEFORMAT_UINT;
	uint16 nbits=8;
	TIFFGetField(tif, TIFFTAG_SAMPLEFORMAT, &frmt);
	TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &nbits);


	std::cout<<"frmt: "<<frmt<<" : "<<nbits<<std::endl;
	switch (frmt) {

		case SAMPLEFORMAT_UINT:
		default:
		  if   (nbits==32)
		  	{ TIFFClose(tif);	return std::unique_ptr<voxelImageTBase>(new voxelImageT<unsigned int>(innam)); }
		  else if(nbits==16)
		  	{ TIFFClose(tif);	return std::unique_ptr<voxelImageTBase>(new voxelImageT<unsigned short>(innam)); }
		  else
			{ TIFFClose(tif);	return std::unique_ptr<voxelImageTBase>(new voxelImageT<unsigned char>(innam)); }
		  break;

		case SAMPLEFORMAT_INT:
		  if   (nbits==32)
			{ TIFFClose(tif);	return std::unique_ptr<voxelImageTBase>(new voxelImageT<int>(innam)); }
		  else if(nbits==16)
			{ TIFFClose(tif);	return std::unique_ptr<voxelImageTBase>(new voxelImageT<short>(innam)); }
		  else
			{ TIFFClose(tif);	return std::unique_ptr<voxelImageTBase>(new voxelImageT<char>(innam)); }
		  break;
		case SAMPLEFORMAT_IEEEFP:
		  if   (nbits==32)
			{ TIFFClose(tif);	return std::unique_ptr<voxelImageTBase>(new voxelImageT<float>(innam)); }
		  else
			{ TIFFClose(tif);	return std::unique_ptr<voxelImageTBase>(new voxelImageT<double>(innam)); }
			break;
		}

	//TIFFDataType type=TIFF_BYTE;
	//TIFFGetField(tif, TIFFTAG_DATATYPE, &type);
	//std::cout<<type<<std::endl;
	//switch (type) {
	//switch (type) {
		//case TIFF_BYTE:
		//case TIFF_ASCII:
		//case TIFF_UNDEFINED:
		//default:
			//{ TIFFClose(tif);	return std::unique_ptr<voxelImageTBase>(new voxelImageT<unsigned char>(innam)); }
		//case TIFF_SBYTE:
			//{ TIFFClose(tif);	return std::unique_ptr<voxelImageTBase>(new voxelImageT<char>(innam)); }
			//break;
		//case TIFF_SHORT:
			//{ TIFFClose(tif);	return std::unique_ptr<voxelImageTBase>(new voxelImageT<unsigned short>(innam)); }
		//case TIFF_SSHORT:
			//{ TIFFClose(tif);	return std::unique_ptr<voxelImageTBase>(new voxelImageT<short>(innam)); }
			//break;
		//case TIFF_LONG:
			//{ TIFFClose(tif);	return std::unique_ptr<voxelImageTBase>(new voxelImageT<unsigned int>(innam)); }
		//case TIFF_SLONG:
			//{ TIFFClose(tif);	return std::unique_ptr<voxelImageTBase>(new voxelImageT<int>(innam)); }
		//case TIFF_IFD:
		//case TIFF_RATIONAL:
		//case TIFF_SRATIONAL:
		//case TIFF_FLOAT:
			//{ TIFFClose(tif);	return std::unique_ptr<voxelImageTBase>(new voxelImageT<float>(innam)); }
			//break;
		//case TIFF_DOUBLE:
			//{ TIFFClose(tif);	return std::unique_ptr<voxelImageTBase>(new voxelImageT<double>(innam)); }
			//break;
		//}
}



#endif
