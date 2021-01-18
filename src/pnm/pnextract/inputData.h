#ifndef CONFIG_H
#define CONFIG_H


#include "InputFile.h"
#include "voxelImage.h"
#include "ElementGNE.h"


using namespace std;

#define MICROCTDATA 1
#define BINARYDATA 2
#define ASCIIDATA 3
#define MHDDATA 4

#ifndef PORORANGE_H
#define PORORANGE_H
class poroRange :
public std::pair<unsigned char,unsigned char>
{
public:
	poroRange(unsigned char lower, unsigned char upper) : std::pair<unsigned char,unsigned char>(lower,upper){};
	poroRange(std::string nam, unsigned char lower,unsigned char upper)
	  : std::pair<unsigned char,unsigned char>(lower,upper),name(nam){};
	poroRange(const poroRange& copyt): std::pair<unsigned char,unsigned char>(copyt),name(copyt.name) {};
	bool outside(unsigned char value){return this->first > value  ||  value > this->second; };
	std::string name;
};
#endif //PORORANGE_H


inline int createSample_input_nextract(std::string fnam, std::string opts)
{
	if (fnam.empty()) fnam = "vxlImage.mhd";


	ensure(!std::ifstream(fnam),"\n\nFile "+fnam+" exists,\n"
	         +     "  to run simulation: rerun with "+fnam+" as the only argument \n"
	         +     "  to regenerate: delete it and try again,\n or provide a different file name:\n"
	         +     "   pnextract -g input_pnextract.mhd\n",  -1);

	ofstream of(fnam);
	if(opts=="-g")
	{
	 of	<<"ObjectType =  Image\n"
		<<"NDims =       3\n"
		<<"ElementType = MET_UCHAR\n"
		<<"ElementByteOrderMSB = False\n"
		<<"ElementNumberOfChannels = 1\n"
		<<"CompressedData = True\n\n"
		<<"HeaderSize = 0\n"
		<<"DimSize =    	1000	1000	1000\n"
		<<"ElementSize = 	1.6 	1.6 	1.6\n"
		<<"Offset =      	0   	0   	0\n"
		<<"\n"
		<<"ElementDataFile = input_image.raw.gz\n"
		<<"\n\n"
		<<"//! The above keywords are compatible with mhd format and \n"
		<<"//! can be used to open the Image in Fiji/ImageJ or Paraview\n"
		<<"//! The following commands are optional, remove the \"//\" to activate them\n"
		<<"\n"
		<<"//DefaultImageFormat = tif\n"
		<<"\n"
		<<"//!______________  image processing  commands _________________\n"
		<<"\n"
		<<"//! crop image to  [ Nxyz_begin  Nxyz_end )"
		<<"//cropD                0 0 0    300 300 300 \n"
		<<"\n"
		<<"//! flip x direction with y or z \n"
		<<"//direction z\n"
		<<"\n"
		<<"//! manipulate voxel values:\n"
		<<"//!    range   [start...end] -> value \n"
		<<"//replaceRange   0     127       0 \n"
		<<"//replaceRange   128   255       1 \n"
		<<"\n"
		<<"//! threshold image: range -> 0 (void-space), rest->1 (solid) \n"
		<<"//threshold   0  128 \n"
		<<"\n\n\n";
	 of	<<"//!_______________  network extraction keywords __________________\n"
		<<"//!______(should be after image processing  commands above) ______\n"
		<<"\n"
		<<"//title:   output_network"
		<<"\n"
		<<"//write_all:	true; // use `write_all` is a memorable alternative to all visualization keywords "
		<<"// write_radius:	true\n"
		<<"// write_statistics:	true\n"
		<<"// write_elements:	true\n"
		<<"// write_poreMaxBalls:	true\n"
		<<"// write_throatMaxBalls:	true\n"
		<<"// write_throats:	true\n"
		//<<"// write_poroats:	true\n" // leads to seg fault
		<<"// write_hierarchy:	true\n"
		<<"// write_medialSurface:	true\n"
		<<"// write_throatHierarchy:	true\n"
		<<"// write_vtkNetwork:	true\n"
		<<endl;
	}
	else
		cout<<"Error unknown option (first argument)"<<endl;

	of.close();

	cout <<" file "<<fnam<<" generated, edit: set Image size, name etc, and rerun\n";
	return 0;

}




class inputDataNE : public InputFile
{
public:
 inputDataNE(const std::string& fnam)
	: InputFile(fnam,false), nx(0), ny(0), nz(0), vxlSize(1), X0(0.0,0.0,0.0),  datatype(BINARYDATA), invalidSeg{-10000, 255,0}
	{}

 inputDataNE(const std::string& fnam, int minvoid, int maxvoid, int outsiderange, bool readCfg)
	: InputFile(false),  nx(0), ny(0), nz(0), vxlSize(1), X0(0.0,0.0,0.0),  datatype(MHDDATA), invalidSeg{-10000, 255,0}
	{
		if(readCfg) 		read(fnam);
		else
		{
			setKeyword("ElementDataFile",fnam);
			fileName_ = fnam;
		}
		if(maxvoid)          setKeyword("void_range",_s(minvoid)+" "+_s(maxvoid));
		if(outsiderange<256) setKeyword("outside_range", _s(outsiderange));
		setTitle(fnam);
		//echoKeywords(cout);
	}

 void init(bool verbos=true)
 {

	cout<< "Reading inputDataNE data:"<<endl;
	echoKeywords(std::cout);

	std::istringstream iss;


	Assert(getVar(imgfileName,"ElementDataFile") || getVar(imgfileName,"shapeToVoxel") || !verbos, "ElementDataFile or shapeToVoxel", "keyword not found", true);
	if(verbos)
		cout<<" image file: "<<imgfileName<<endl;


	Assert(getData(iss, "DimSize") ||  getData(iss, "imageSize") || !verbos, "DimSize or imageSize", "keyword not found", false);
		iss >> nx >> ny >> nz;
	if(verbos)
		cout<< " image size: " << nx << " " << ny << " " << nz <<endl;


	Assert(getVar(vxlSize, "ElementSize") || getVar(vxlSize,  "ElementSpacing") || getVar(vxlSize, "voxelSize") || !verbos, "ElementSpacing or voxelSize", "keyword not found", false);
	if(verbos)
		cout<< " voxel size: " <<vxlSize <<endl;


	std::string dataType("binary char");
	if (getVar(dataType,"ElementType") || getVar(dataType,"fileType"))
	{
		if (dataType == "binary")        	datatype = BINARYDATA;
		else if (dataType == "MET_UCHAR")	datatype = MHDDATA;
		else if	(dataType == "microct")  	datatype = MICROCTDATA;
		else if	(dataType == "ascii")    	datatype = ASCIIDATA;
		else cout<<"wrong data type, going with binary"<<endl;
		cout<<"  file type: "<<dataType<<", "<<datatype <<endl;
	}
	else if	(fileName_.size()>4 && fileName_.compare(fileName_.size()-4,4,".mhd")==0)	datatype = MHDDATA;

	if (!getVar(imgfrmt,"DefaultImageFormat")) imgfrmt=".raw.gz";
	if(imgfrmt[0]!='.') imgfrmt="."+imgfrmt;
	imgExt(imgfrmt);
	cout<<"DefaultImageFormat: "<<imgfrmt<<endl;


	cout<<" voxel indices:"<<endl;
	if (getData(iss,"voidSpace"))
	{
		std::string rockTypeName;
		iss>>rockTypeName;
		cout<<"  "<<0<<": rockTypeName = "<<rockTypeName<<endl;
		poroRange ithRockType(rockTypeName,0,0);

		_rockTypes.push_back(ithRockType);
	}
	else
	{
		poroRange ithRockType("void",0,0);
		_rockTypes.push_back(ithRockType);
		cout<<"  "<<0<<": void voxels "<<endl;
	}


	segValues.resize(256, _rockTypes.size());

	for_i(_rockTypes)
	{	auto& rt=_rockTypes[i];
		if(getData(iss, rt.name+"_range") || getData(iss, rt.name+"_thresholds"))
		{
			int lower,upper;
			iss >> lower >> upper;
			rt.first = lower;
			rt.second = upper;
			if (rt.first > rt.second)
				cout<<"  Wrong entries for keyword \""<<rt.name+"_range"<<"\":\n"<<keyvals(rt.name+"_thresholds")<<"\n lower value is higher than upper value"<<endl;
		}

		for(size_t j = rt.first; j <=  rt.second; ++j)
			segValues[j] = i;

		cout<<"  "<< rt.name<<" voxel values: ["<<int(rt.first)<< " "<<int(rt.second)<<"]"<<endl;
	}

	cout<<"  Voxel value indices:";
	for(size_t i = 0; i <=  12; ++i)		cout<<" "<<segValues[i];
	cout<<" ... "<<endl;




 }


 void readImage()
 {
	cout<< "\nLoading voxel data, format:"<<datatype <<" fileName:"<<fileName()<<endl;

	VImage.reset(nx,ny,nz,255);
	if (datatype == MICROCTDATA)	VImage.readMicroCT(imgfileName);
	else if (datatype == BINARYDATA)  {if (!VImage.readBin(imgfileName))   {cout<<"\nError: didn't read binary image!\n"<<endl; exit(-1);}}
	else if (datatype == ASCIIDATA)   {if (!VImage.readAscii(imgfileName)) {cout<<"\nError: didn't read ascii image!\n"<<endl; exit(-1);}}
	else
	{
		VImage.reset(0,0,0,255);
		readConvertFromHeader(VImage, fileName());
		int3 siz=VImage.size3();
		nx= siz[0];  ny= siz[1];  nz= siz[2];
		vxlSize = VImage.dx()[0];
		X0=VImage.X0();
		cout<<" siz:"<<siz<<" vxlSize:"<<vxlSize<<" X0:"<<X0<<endl;
		ensure(siz[2],"no image read",2);
	}
	if(!VImage.nx()) {  cout<<"\nError: no image read!\n"<<endl; exit(-1);}
	VImage.printInfo();

	nInside= (long long)(nx)*ny*nz;
	std::istringstream iss;
	if(getData(iss, "outside_range") || getData(iss,"outside_thresholds"))
	{
		int lower,upper;
		iss >> lower >> upper;
		forAllcp(VImage) if(lower<=(*cp) && (*cp)<=upper) --nInside;
		cout<<" outside_range: "<<lower <<" "<<upper<<"; inside_fraction: "<<nInside/(double(nx)*ny*nz)<<endl;

	}
 }

 void createSegments()
 {




	std::vector<size_t> nVxlVs(_rockTypes.size()+1,0);
	segs_.resize(nz,std::vector<segments>(ny));



	#ifdef OpenMP
		#pragma   omp declare reduction(vec_sizet_plus : std::vector<size_t> : \
				  std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<size_t>())) \
				  initializer(omp_priv = omp_orig)
	#endif

	OMPragma("omp parallel for  reduction(vec_sizet_plus : nVxlVs)")
	for (int iz = 0; iz<nz; ++iz)
	{
		std::vector<segment> segTmp(nx+1);
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
			ss.s[cnt].value = 254;
		}
	}


	cout<< endl;
	size_t nVInsids=0.0;
	for_i(_rockTypes)
	{
		nVInsids+=nVxlVs[i];
		cout<<" "<<i<<". "<< _rockTypes[i].name<<": "<<nVxlVs[i]<<" voxels, "<< nVxlVs[i]/(double(0.01*nx)*ny*nz) << "%"<<endl;
	}
	ensure(nVInsids>double(0.01*nx)*ny*nz, "too low porosity, set 'void_range' maybe?", -1);
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


	std::string netsufix() const { return (flowBaseDir.empty() ? "DS0" : (flowBaseDir.back()=='/' ? "DS1": "DS4")); }
public:

	int nx, ny, nz;
	double vxlSize;
	dbl3 X0;
	int datatype;
	std::string imgfrmt;
	long long nInside;

	std::string imgfileName;
	std::string flowBaseDir;

	std::vector< std::vector<segments> > segs_;
	segment invalidSeg;

	std::vector<int> segValues;
	std::vector<poroRange> _rockTypes;
	voxelImage VImage;

};

#endif
