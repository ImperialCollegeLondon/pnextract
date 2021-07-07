/*-------------------------------------------------------------------------*\
You can redistribute this code and/or modify this code under the
terms of the GNU General Public License (GPL) as published by the 
Free Software Foundation, either version 3 of the License, or (at
your option) any later version. see <http://www.gnu.org/licenses/>.

The code is part of voxLib, developed by Ali Qaseminejad Raeini

For further information please contact me by email:
Ali Q Raeini: a.q.raeini@imperial.ac.uk

\*-------------------------------------------------------------------------*/




#include "voxelImage.h"

using namespace std;


#include "voxelImage.cpp"

int usage()  {
	cout<<"\nvoxelImageProcess:\n utility to read a 3D image\n optionally run some image processing commands on it (cropD, threshold, ...),"
	<<"\n and write it back to disk,\n in the same or in a different format\n .dat suffix is used for ascii files and .raw implies binary format. \n Other image formats supported are .tif and .am  and compressed .raw.gz"
		"\nusages: \n"
		"    voxelImageProcess inputMetaImage.mhd outPutFile.tif \n"
		"    voxelImageProcess inputMetaImage.tif outPutFile.mhd UChar\n"
		"    voxelImageProcess inputMetaImage.tif outPutFile.dat \n"
		"    voxelImageProcess inputMetaImage.tif outPutFile.raw \n"
		"    voxelImageProcess inputMetaImage.mhd outPutFile.tif \n"
		"    voxelImageProcess inputMetaImage.mhd NO_WRITE \n"
		"Header file contents should be like: \n"
		" ObjectType  = Image\n"
		" NDims       = 3\n"
		" ElementType = MET_SHORT\n"
		"\n"
		" DimSize     = 650  650  650\n"
		" ElementSize = 5.   5.   5.\n"
		" Offset      = 0    0    0\n"
		"\n"
		" ElementDataFile = Berea.tif\n"
		"\n"
		"Optional keywords:\n"
		<< VxlKeysHelp("","")<<endl;

		cout<<"\nFor argumrntd of individual keywords run:\n"<<" vxlImageProcess ? <keyword name>\n"<<endl;
	return -1;
}

int main(int argc, char *argv[])  {

	if(argc<2 || argc>4) return usage();

	string header(argv[1]);
	string outputName(argc>2 ? argv[2] : "");
	string extra(argc>3 ?  argv[3] : "" );

	if(header.size()<4) return usage();

	if(header[0]=='?') {
		cout<<outputName<<" "<<VxlKeysHelp(outputName,extra)<<endl;
		return 0; }

	cout<<"//-*-C-*-\\ voxelImageProcess, in: "<<header<<",  out:"<<outputName<<endl;
	if(!outputName.size()) cerr<<"\n\nWarning: no output (2nd argument) \n\n"<<endl;


	unique_ptr<voxelImageTBase> vxlImage = readImage(header);

	if(extra=="UChar")  {
		voxelImage vImgUChar(vxlImage->size3(), vxlImage->dx(), vxlImage->X0(), 0);
		//forAlliii_(vImgUChar)			vImgUChar(iii) = vxlImage->getInt(iii);
		resetFromImageT<unsigned char,short,unsigned short,int,unsigned int,float,double, char>(vImgUChar,vxlImage.get());
		vImgUChar.write(outputName);
	}
	else if(outputName.size()) vxlImage->write(outputName);


	cout<< "end" << endl;

	vxlImage->printInfo();


	return 0;
}


