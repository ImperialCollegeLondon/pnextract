#ifndef CONFIG_H
#define CONFIG_H


#include "inputFileGNE.h"
#include "voxelImage.h"
#include "ElementGNE.h"



#define MICROCTDATA 1
#define BINARYDATA 2
#define ASCIIDATA 3
#define MHDDATA 4


class threshold :
public std::pair<unsigned char,unsigned char>
{
public:
	threshold(unsigned char lower, unsigned char upper) : std::pair<unsigned char,unsigned char>(lower,upper){};
	threshold(std::string nam, unsigned char lower,unsigned char upper)
	  : std::pair<unsigned char,unsigned char>(lower,upper),name(nam){};
	threshold(const threshold& copyt): std::pair<unsigned char,unsigned char>(copyt),name(copyt.name) {};
	bool outside(unsigned char value){return this->first > value  ||  value > this->second; };
	std::string name;
};


class inputDataNE : public InputFileNE
{
public:
 inputDataNE(std::string file, int minvoid, int maxvoid, int outsiderange, bool readCfg) :
	InputFileNE(file,readCfg),  precision(1), X0(0.0,0.0,0.0),  datatype(MHDDATA), invalidSeg{-10000, 255}
	{
		if(!readCfg)
		{
			setKeyword("ElementDataFile",file);
			setKeyword("void_range",toStr(minvoid)+" "+toStr(maxvoid));
			if(outsiderange<256) setKeyword("outside_range", toStr(outsiderange));
		}
		init(true);
		readImage();
		createSegments();
	}
 inputDataNE(std::string file) :
	InputFileNE(file), precision(1), X0(0.0,0.0,0.0),  datatype(BINARYDATA), invalidSeg{-10000, 255}
	{
		init();
		readImage();
		createSegments();
	}
 void init(bool ignor=false)
 {

	cout<< "Reading inputDataNE data:"<<endl;

	std::istringstream inputKeyData;


	Assert(getVar(imgfileName,"ElementDataFile") || getVar(imgfileName,"imageFile") || ignor, "ElementDataFile or imageFile", "keyword not found");
	if(!ignor)
		cout<<" image file: "<<imgfileName<<endl;



	Assert(getData(inputKeyData, "DimSize") ||  getData(inputKeyData, "imageSize") || ignor, "DimSize or imageSize", "keyword not found");
		inputKeyData >> nx >> ny >> nz; ///. not fixed yet
	if(!ignor)
		cout<< " image size: " << nx << " " << ny << " " << nz <<endl;


	Assert(getVar(precision, "ElementSpacing") || getVar(precision, "voxelSize") || ignor, "ElementSpacing or voxelSize", "keyword not found");
	if(!ignor)
		cout<< " voxel size: " <<precision <<endl;


	std::string dataType("binary char");
	if (getVar(dataType,"ElementType") || getVar(dataType,"fileType"))
	{
		if (dataType == "binary")	datatype = BINARYDATA;
		if (dataType == "MET_UCHAR")	datatype = MHDDATA;
		else if	(dataType == "microct")	                        datatype = MICROCTDATA;
		else if	(dataType == "ascii")	                        datatype = ASCIIDATA;
		else cout<<"wrong data type, going with binary"<<endl;
		cout<<"  file type: "<<dataType <<endl;
	}

	if (!getVar(imgfrmt,"DefaultImageFormat")) imgfrmt=".tif";
	if(imgfrmt[0]!='.') imgfrmt="."+imgfrmt;
	suffix(imgfrmt);
	cout<<"DefaultImageFormat: "<<imgfrmt<<endl;
	

	cout<<" voxel indices:"<<endl;
	if (getData(inputKeyData,"voidSpace"))
	{
		std::string rockTypeName;
		inputKeyData>>rockTypeName;
		cout<<"  "<<0<<": rockTypeName = "<<rockTypeName<<endl;
		threshold ithRockType(rockTypeName,0,0);

		_rockTypes.push_back(ithRockType);
	}
	else
	{
		threshold ithRockType("void",0,0);
		_rockTypes.push_back(ithRockType);
		cout<<"  "<<0<<": void voxels "<<endl;
	}

	int nRTypes(0);
	if(getData(inputKeyData,"porousRocks"))
	{
		inputKeyData >> nRTypes;
		for(int i = 1; i <=  nRTypes; ++i)
		{
				std::string rockTypeName;
				inputKeyData>>rockTypeName;
				cout<<"  "<<i<<": porousRocks = "<<rockTypeName<<endl;
				threshold ithRockType(rockTypeName,i,i);
				_rockTypes.push_back(ithRockType);
		}
		cout<< "  number of rock types: "<<_rockTypes.size()<<endl;
	}

	segValues.resize(256, _rockTypes.size());

	for(int i = 0; i <=  nRTypes; ++i)
	{
		if(getData(inputKeyData, _rockTypes[i].name+"_range") || getData(inputKeyData, _rockTypes[i].name+"_thresholds"))
		{
			int lower,upper;
			inputKeyData >> lower >> upper;
			 _rockTypes[i].first = lower;
			 _rockTypes[i].second = upper;
			if (_rockTypes[i].first > _rockTypes[i].second)
				cout<<"  Wrong entries for keyword \""<<_rockTypes[i].name+"_range"<<"\":\n"<<keywordData(_rockTypes[i].name+"_thresholds")<<"\n lower value is higher than upper value"<<endl;
		}

		for(size_t j = _rockTypes[i].first; j <=  _rockTypes[i].second; ++j)
			segValues[j] = i;

		cout<<"  "<< _rockTypes[i].name<<" voxel values: ["<<int(_rockTypes[i].first)<< " "<<int(_rockTypes[i].second)<<"]"<<endl;
	}

	cout<<"  Voxel value indices:";
	for(size_t i = 0; i <=  12; ++i)		cout<<" "<<segValues[i];
	cout<<" ... "<<endl;

 }



 void readImage()
 {
	cout<< "\nLoading voxel data:"<<endl;

	VImage.reset(nx,ny,nz,0);
	if (datatype == MICROCTDATA)	VImage.readMicroCT(imgfileName);
	else if (datatype == BINARYDATA)	{if (!VImage.readBin(imgfileName))   {cout<<"\nError: didn't read binary image!\n"<<endl; exit(-1);}}
	else if (datatype == ASCIIDATA)	{if (!VImage.readAscii(imgfileName)) {cout<<"\nError: didn't read ascii image!\n"<<endl; exit(-1);}}
	else
	{
		VImage.reset(0,0,0,255);
		VImage.readFromHeader(fileName());
		int3 siz=VImage.size3();
		nx= siz[0];  ny= siz[1];  nz= siz[2];
		precision = VImage.dx()[0];
		X0=VImage.X0();
	}
	if(!VImage.size3()[0]) {  cout<<"\nError: no image read!\n"<<endl; exit(-1);}
	VImage.printInfo();

	nInside= (long long)(nx)*ny*nz;
	std::istringstream inputKeyData;
	if(getData(inputKeyData, "outside_range") || getData(inputKeyData,"outside_thresholds"))
	{
		int lower,upper;
		inputKeyData >> lower >> upper;
		forAllcp(VImage) if(lower<=(*cp) && (*cp)<=upper) --nInside;

	}
 }
 void createSegments()
 {
	nVxlVs.resize(_rockTypes.size()+1,0);

	std::vector<segment> segTmp(nx+1);
 	segs_.resize(nz,std::vector<segments>(ny));

 	for (int iz = 0; iz<nz; iz++)
 	 for (int iy = 0; iy<ny; iy++)
 		if(segs_[iz][iy].s != NULL) cout<<"ERROR"<<endl;


	for (int iz = 0; iz<nz; ++iz)
	{
 		for (int iy = 0; iy<ny; ++iy)
 		{
			int cnt = 0;
			int  currentSegValue = 257;
			for (int ix = 0; ix<nx; ++ix)
			{
				unsigned char vV = VImage(ix,iy,iz);
				if (segValues[vV] != currentSegValue)
				{
					currentSegValue = segValues[vV];
					segTmp[cnt].start = ix;
					segTmp[cnt].value = currentSegValue;
					++cnt;
				}
				++(nVxlVs[segValues[vV]]);
			}

			segments & ss = segs_[iz][iy];
			ss.cnt = cnt;
			ss.reSize(cnt+1);
			if (segs_[iz][iy].s == NULL)		{cout<<"\n   iz "<<iz<<" iy "<<iy<<" cnt"<<cnt<<" ERROR XX XX"<<endl;	}

			for (int i = 0; i<cnt; ++i)
			{
				ss.s[i].start = segTmp[i].start;
				ss.s[i].value = segTmp[i].value;
				if (i>0 && ss.s[i].value == ss.s[i-1].value) 	cout<<"\n   ERROR XXX"<<i<<" "<<int(ss.s[i-1].value)<<" "<<int(ss.s[i].value)<<"    "<<int(segValues[3])<<"    "<<int(VImage(ss.s[i-1].start,iy,iz))<<" "<<int(VImage(ss.s[i].start,iy,iz))<<endl;
			}
			ss.s[cnt].start = nx;
			ss.s[cnt].value = 254; ///. Warning: 254 is reserved
			if (cnt>0 && ss.s[cnt].value == 0) 			{cout<<"\n   ERROR XX "<<endl; 	}
 		}
 	}


	cout<< endl;
	for (int i = 0;i<int(_rockTypes.size());i++)
		cout<<" "<<i<< ". " << _rockTypes[i].name<< ": " << nVxlVs[i] << " voxels, " <<  (nVxlVs[i]*(100.0/nx/ny/nz))<< "%"<<endl;
	cout<< endl;

 }


 const segment* segptr(int i, int j, int k) const
 {
 	if (i<0 || j<0 || k<0 || i>= nx || j>= ny || k>= nz)  return &invalidSeg;

 	const segments& s = segs_[k][j];
 	for (int p = 0; p<s.cnt; ++p)
 		if (i >= s.s[p].start && i < s.s[p+1].start)	  return s.s+p;

	cout<<"Error can not find segment at "<<i<<" "<<j<<" "<<k<<" nSegs: "<<s.cnt<<endl;
 	return (s.s+s.cnt);
 }


public:

	int nx, ny, nz;
	double precision;
	vec3 X0;
	int datatype;
	std::string imgfrmt;
	long long nInside;

	std::string imgfileName;

	std::vector< std::vector<segments> > segs_;
	segment invalidSeg;

	std::vector<int> nVxlVs;
	std::vector<int> segValues;
	std::vector<threshold> _rockTypes; ///. vValue-index pairs segValues
	voxelImage VImage;

};

#endif
