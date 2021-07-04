#ifndef CONFIG_H
#define CONFIG_H


#include "InputFile.h"
#include "voxelImage.h"
#include "ElementGNE.h"


using namespace std;


#ifndef PORORANGE_H
#define PORORANGE_H
class poroRange : public std::pair<unsigned char,unsigned char> {

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
	if(opts=="-g")  {
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
	: InputFile(fnam,false), nx(0), ny(0), nz(0), vxlSize(1), X0(0.,0.,0.),  invalidSeg{-10000, 255,0}
	{}



 void init(bool verbos=true)  {

	if (!giv("DefaultImageFormat", imgfrmt)) imgfrmt=".raw.gz";
	if(imgfrmt[0]!='.') imgfrmt="."+imgfrmt;
	imgExt(imgfrmt);
	cout<<"DefaultImageFormat: "<<imgfrmt<<endl;

	nBP6 = getOr("multiDir",false) ? 6 : 2;

	cout<<" voxel indices:"<<endl;
	_rockTypes.push_back(poroRange("void",0,0));
	cout<<"  "<<0<<": void voxels "<<endl;

	std::istringstream iss;


	segValues.resize(256, _rockTypes.size());

	for_i(_rockTypes)  {
		auto& rt=_rockTypes[i];
		if(giv(rt.name+"_range", iss))  {
			int lower,upper;
			iss >> lower >> upper;
			rt.first = lower;
			rt.second = upper;
			if (rt.first > rt.second)
				cout<<"  Wrong entries for keyword \""<<rt.name+"_range"<<"\":  lower value is higher than upper value"<<endl;
		}

		for(size_t j=rt.first; j<=rt.second; ++j)  segValues[j] = i;

		cout<<"  "<< rt.name<<" voxel values: ["<<int(rt.first)<< " "<<int(rt.second)<<"]"<<endl;
	}

	cout<<"  Voxel value indices:";
	for(size_t i=0; i<=12; ++i)		cout<<" "<<segValues[i];
	cout<<" ... "<<endl;

 }


 void readImage() {


	std::string fnam(fileName());
	if(fnam.empty()) { giv("ElementDataFile",fnam) || giv("read",fnam); }
	//Assert(giv("ElementDataFile") || kwrd("read").size()



	//std::string dataType("mhd");  giv("ElementType", dataType) || giv("fileType", dataType);
	cout<<" Image file: "<<fnam<<endl;

	//VImage.reset(nx,ny,nz,255);
	//if (dataType == "microct")	VImage.readMicroCT(fnam);
	//else if (dataType == "binary")  {if (!VImage.readBin(fnam))   alert("could not read binary image!",-1); }
	//else if (dataType == "ascii")   {if (!VImage.readAscii(fnam)) alert("could not read ascii image!",-1); }
	//else  { VImage.reset(0,0,0,255);  
	readConvertFromHeader(VImage,fnam);  //dataType="mhd"; }

	const std::string& vxlkys=kwrd("VxlPro");
	if(vxlkys.size())  vxlProcess(vxlkys,VImage,"GNE:VxlPro");

	//if(dataType=="mhd") 
	//{
	vxlSize = VImage.dx().x;
	X0=VImage.X0();
	//}
	nx=VImage.nx();  ny=VImage.ny();  nz=VImage.nz();
	ensure(nz,"no image read",2);
	VImage.printInfo();

	nInside= (long long)(nx)*ny*nz;
	cout<<" siz:"<<VImage.size3()<<" vxlSize:"<<vxlSize<<" X0:"<<X0<<endl;
	
	int2 outrange;
	if(giv("outside_range", outrange))  {
		forAllcp(VImage) if(outrange.a<=(*cp) && (*cp)<=outrange.b) --nInside;
		cout<<" outside_range: "<<outrange<<"; inside_fraction: "<<nInside/(double(nx)*ny*nz)<<endl;
	}
 }

 void createSegments() {





	stvec<size_t> nVxlVs(_rockTypes.size()+1,0);
	segs_.resize(nz,stvec<segments>(ny));



	#ifdef OpenMP
		#pragma   omp declare reduction(vec_sizet_plus : std::vector<size_t> : \
				  std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<size_t>())) \
				  initializer(omp_priv = omp_orig)
	#endif

	OMPragma("omp parallel for  reduction(vec_sizet_plus : nVxlVs)")
	for (int iz = 0; iz<nz; ++iz)  {
		stvec<segment> segTmp(nx+1);
		for (int iy = 0; iy<ny; ++iy)  {
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
			if (segs_[iz][iy].s == nullptr)		{cout<<"\n   iz "<<iz<<" iy "<<iy<<" cnt"<<cnt<<" ERROR XX XX"<<endl;	}

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
	size_t nVInsids=0;
	for_i(_rockTypes)  {
		nVInsids+=nVxlVs[i];
		cout<<" "<<i<<". "<< _rockTypes[i].name<<": "<<nVxlVs[i]<<" voxels, "<< nVxlVs[i]/(double(0.01*nx)*ny*nz) << "%"<<endl;
	}
	ensure(nVInsids>double(0.01*nx)*ny*nz, "too low porosity, set 'void_range' maybe?", -1);
	cout<< endl;

 }


 const segment* segptr(int i, int j, int k) const {

	if (i<0 || j<0 || k<0 || i>= nx || j>= ny || k>= nz)  return &invalidSeg;

	const segments& s = segs_[k][j];
	for (int p = 0; p<s.cnt; ++p)
		if (i >= s.s[p].start && i < s.s[p+1].start)	  return s.s+p;

	cout<<"Error can not find segment at "<<i<<" "<<j<<" "<<k<<" nSegs: "<<s.cnt<<endl;
	return (s.s+s.cnt);
 }


	std::string netName() const { return name()+(flowBaseDir.empty() ? "DS0" : (flowBaseDir.back()=='/' ? "DS1": "DS4")); }
public:

	int                    nx, ny, nz;
	int                    nBP6;
	double                 vxlSize;
	dbl3                   X0;
	std::string            imgfrmt;
	long long              nInside;

	std::string            flowBaseDir;

	stvec<stvec<segments>> segs_;
	segment                invalidSeg;

	stvec<int>             segValues;
	stvec<poroRange>      _rockTypes;
	voxelImage             VImage;

};

#endif
