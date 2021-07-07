/*-------------------------------------------------------------------------*\
You can redistribute this code and/or modify this code under the
terms of the GNU General Public License (GPL) as published by the
Free Software Foundation, either version 3 of the License, or (at
your option) any later version. see <http://www.gnu.org/licenses/>.

This file is part of voxelImage library, a C++ template library  
developed by Ali Qaseminejad Raeini for handelling 3D raw images.


Please see our website for relavant literature making use of this code:
https://www.imperial.ac.uk/earth-science/research/research-groups/pore-scale-modelling/

For further information please contact us by email:
Ali Q Raeini: a.q.raeini@imperial.ac.uk

\*-------------------------------------------------------------------------*/
#include <memory>
#include <sstream>
#include "shapeToVoxel.h"
#include "voxelEndian.h"
#include "globals.h"  // ensure...
#include "voxelRegions.h"
#include "InputFile.h"

using namespace std; //cin cout endl string stringstream  istream istringstream regex*

#ifdef LPNG  // PRI_0:
#include "voxelPng.h"
#endif       // PRI_0;


int maxNz = 50000;


								namespace MCTProcessing _begins_


template<typename T> bool ignore( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("{ section to ignore }");
	return 0;
}

template<typename T> bool fillHoles( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("maxHoleSize // fill small isolated features");
	unsigned int maxHoleSize;
	ins>>maxHoleSize;

	cout<<"  fillHoles: eliminating isolated rocks/pores; maxHoleSize:" <<maxHoleSize<<" (default is 2) "<<endl;
	vImg.fillHoles(maxHoleSize);

	vImg.FaceMedian06(1,5);
	//vImg.FaceMedian07(2,5);
	//vImg.FaceMedian07(2,5);
	return 0;
}

template<typename T> bool info( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("// print porosity...");
	vImg.printInfo();
	return 0;
}

template<typename T> bool selectPore( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("thresholdBegin ThresholdEnd;// segment the image");
	cout<<"  Converting to binary (0 and 1):";
	unsigned int  thresholdMin=0,thresholdMax=0;
	ins>>thresholdMin;
	ins>>thresholdMax;

	(cout<<"  pore (=0) <- ["<<int(thresholdMin)<<" "<<int(thresholdMax)<<"] ").flush();
	vImg.threshold101(thresholdMin,std::min(unsigned(maxT(T)),thresholdMax));
	return 0;
}

template<typename T,  enable_if_t<std::is_arithmetic<T>::value, int> = 0> 
bool rescale( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("range_min range_max;// rescale values to be within range");
	(cout<<"  rescaling voxel values to [ ").flush();
	unsigned int  thresholdMin=0,thresholdMax=0;
	ins>>thresholdMin;
	ins>>thresholdMax;

	(cout<<thresholdMin<<", "<<thresholdMax<<" ]    ").flush();
	rescale(vImg,T(thresholdMin),T(thresholdMax));
	(cout<<".").flush();
	return 0;
}

template<typename T> bool growPore( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("valToGrow(0/1) Algorithm(f) ...;// grow labels to adjacent voxels");
	cout<<"  growing voxels: ";
	int vvalTogrow(0); ins>>vvalTogrow;
	char alg('f'); ins>>alg;

	while (ins.good())		  // loop while extraction from file is possible
	{
		if (alg=='f')  {
		cout<<"   "<<vvalTogrow<<" "<<alg<<" ";

			if(vvalTogrow==0)
				vImg.growPore();
			else if (vvalTogrow==1)
				vImg.shrinkPore();
			else
			{
				std::cerr<<"growing is only implemented for binary images: "<<
				"selected voxel value to grow is "<<vvalTogrow << ", which is not acceptable"<<endl;
				return false;
			}
		} else {
			std::cerr<<"selected growing algorithm: "<<alg<<
			" the only implemented algorithm is f which stands for faceGrowing"<<endl;
			return false;
		}

		ins>>vvalTogrow;
		ins>>alg;
	}
	cout<<"."<<endl;
	return 0;
}


template<typename T> bool resampleMean( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("by_resampleFactor // assigning mean val");
	double nResample=1;
	(ins>>nResample, cout<<__FUNCTION__<<" factor: "<<nResample<<" ").flush();
	vImg = resampleMean(vImg,nResample);
	return 0;
}


template<typename T> bool resampleMax( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("by_resampleFactor // assigning max val");
	double nResample=1;
	(ins>>nResample, cout<<__FUNCTION__<<" factor: "<<nResample<<" ").flush();
	vImg = resampleMax(vImg,nResample);
	return 0;
}

template<typename T> bool resliceZ( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("nResample(-1)// reslize z to make dz same as dx and dy");
	double nResample=-1;
	(ins>>nResample, cout<<__FUNCTION__<<" factor: "<<nResample<<" ").flush();
	vImg = resliceZ(vImg,nResample);
	return 0;
}

template<typename T> bool resampleMode( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("nResample //  assigning mode val");
	double nResample=1;
	(ins>>nResample, cout<<__FUNCTION__<<" factor: "<<nResample<<" ").flush();
	vImg = resampleMode(vImg,nResample);
	return 0;
}

template<typename T> bool redirect( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("directiom(y/z) // flip x with y or z axes");
	char axs;
	ins>>axs;
	(cout<<axs<<", swapping x and "<<axs<<" axes").flush();

	vImg.rotate(axs);
	cout<<endl;
	return 0;
}

template<typename T> bool replaceRange( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("FromVal ToValue NewValue");
	int  thresholdMin(0),thresholdMax(0); //. Warning don't use T, uchar wont work
	ins >> thresholdMin >> thresholdMax;

	int  value=(thresholdMin+thresholdMax)/2; //. Warning don't use T, uchar wont work
	ins >> value;

	cout<<"  Replacing range  ["<<thresholdMin<<"  "<<thresholdMax<<"] with "<<value<<", ";
	replaceRange(vImg,T(thresholdMin),T(thresholdMax),T(value));
	(cout<<".").flush();
	return 0;
}


template<typename T> bool cropD( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("cropBegin(0 0 0) cropEnd(nx ny nz)");
	int3 cropBegin(0,0,0), cropEnd=vImg.size3();
	int nLayers(0); int value(1);
	ins>>cropBegin>>cropEnd >> nLayers >> value;  
	(cout<<" "<<cropBegin<<" -- "<<cropEnd<<" ").flush();
	if (nLayers) { cout<<"  + "<<nLayers<<" layers of "<<value<<" "<<endl; } 
	vImg.cropD(cropBegin,cropEnd,nLayers,value,true);
	return 0;
}

template<typename T> bool cropf( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("beginFraction(0 0 0) endFraction(0.5 0.5 0.5)");
	dbl3 bgn(0,0,0), end(1,1,1);
	int nLayers(0); int value(1);
	ins>>bgn>>end >> nLayers >> value;  
	(cout<<" "<<bgn<<" -- "<<end<<" ").flush();
	if (nLayers) { cout<<"  + "<<nLayers<<" layers of "<<value<<" "<<endl; } 
	vImg.cropD(bgn*dbl3(vImg.size3())+0.5,end*dbl3(vImg.size3())+0.5,nLayers,value,true);
	return 0;
}

template<typename T> bool write( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("outputImageName.raw/mhd/tif/am/.raw.gz");
	string outName("dump.tif");    ins >> outName;
	vImg.write(outName);
	(cout<<".").flush();
	return 0;
}

template<typename T> bool write8bit( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("outputImageName_8bit.raw/mhd/tif/am/.raw.gz");
	string outName("dump.tif");     ins >> outName;
	double minv=-0.5, maxv=255.;   ins>>minv>>maxv;
	double delv=255.499999999/(maxv-minv);
	(cout<<minv<<" "<<maxv).flush();
	voxelImageT<unsigned char> voxels(vImg.size3(),vImg.dx(),vImg.X0(),255);
	forAlliii_(voxels) voxels(iii)=std::max(0,std::min(255,int(delv*(vImg(iii)-minv))));
	voxels.write(outName);
	(cout<<".").flush();
	return 0;
}

template<typename T> bool read( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("ImageToRead.mhd/.am/.tif");
	int3 nnn = vImg.size3();
	int processHdr=1;  string fnam;   ins>>fnam>>processHdr;
	cout<<"  reading from  image "<<fnam<<endl;
	if(fnam.size()>4)  {
		if ((nnn[2] && (hasExt(fnam,7,".raw.gz") || hasExt(fnam,4,".raw"))) || hasExt(fnam,4,".tif") )  {
			vImg.reset(nnn,T(0));
			vImg.readBin(fnam);
		}
		else vImg.readFromHeader(fnam,processHdr);
	}
	return 0;
}

template<typename T> bool readAtZ( stringstream& ins, voxelImageT<T>& vImg)  {// used to stitch images
	KeyHint("ImageToReadAndReplacePreviousFromZ.mhd/.am iSlic");
	int3 nnn = vImg.size3();
	size_t iSlic=0;	string fnam;	ins>>fnam>>iSlic;
	cout<<"  reading from  image "<<fnam<<", assigning to slices after "<<iSlic<<endl;
	voxelImageT<T> img(fnam);
	ensure(img.nx()==nnn.x);	ensure(img.ny()==nnn.y);
	std::copy(img.begin(),img.end(),vImg.begin()+iSlic*nnn[0]*nnn[1]);
	return 0;
}



template<typename T> bool medianFilter( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("nIterations");
	int nIterations(1);  ins >> nIterations;
	(cout<<"  median Filter, nIterations: "<<nIterations).flush();
	vImg.growBox(2);
	for (int i=0; i<nIterations; ++i)  vImg=median(vImg);
	vImg.shrinkBox(2);
	(cout<<".").flush();
	return 0;
}

template<typename T> bool modeFilter( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("nIterations(1) nMinNeis(2)");
	int nIterations(1), nMinNeis(2);   ins >> nIterations >> nMinNeis;
	(cout<<"  mode Filter, nIterations: "<<nIterations<<"  nMinNeis"<<nMinNeis).flush();
	vImg.growBox(2);
	for (int i=0; i<nIterations; ++i)  modeNSames(vImg,nMinNeis,true);
	vImg.shrinkBox(2);
	(cout<<".").flush();
	return 0;
}
template<typename T> bool medianX( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("nIterations// only applied in x direction to reduce compressed file size");
	int nIterations(1);   ins >> nIterations;
	(cout<<"  median Filter, nIterations: "<<nIterations).flush();
	for (int i=0; i<nIterations; ++i)
		vImg=medianx(vImg);
	(cout<<".").flush();
	return 0;
}


template<typename T> bool FaceMedian06( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("nAdj0(2), nAdj1(4),  nIterations(1)");
	int nAdj0(2), nAdj1(4),  nIterations(1);     ins >>nAdj0 >>nAdj1 >>nIterations;
	(cout<<"  FaceMedian06: "<<nAdj0<<" "<<nAdj1<<" "<<nIterations<<"     ").flush();
	vImg.growBox(2);
	for (int i=0; i<nIterations; ++i) vImg.FaceMedian06(nAdj0,nAdj1);
	vImg.shrinkBox(2);
	(cout<<".").flush();
	return 0;
}



template<typename T> bool PointMedian032( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("nItrs(1),  nAdj(11), lbl0(0), lbl1(1)");
	int nItrs(1),  nAdj(11), lbl0(0), lbl1(1);
	ins >> nItrs >> nAdj >> lbl0 >> lbl1;
	(cout<<"  PointMedian032, "<<" nItrs:"<<nItrs<< "; nAdjThreshold "<<nAdj<<"  lbl0:"<<lbl0<<"  lbl1;"<<lbl1<<"s \n  PointMedian032 is only applied to the labels  lbl0 and  lbl1").flush();
	//vImg.growBox(2);

	for (int i=0; i<nItrs; ++i)  vImg.PointMedian032(nAdj,nAdj,lbl0,lbl1);

	//vImg.shrinkBox(2);
	(cout<<".").flush();
	return 0;
}


template<typename T> bool faceMedNgrowToFrom( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("nItrs(2),  lblTo(0), lblFrm(1), ndif(-3)");
	int nItrs(2),  ndif(-3); 
	Tint lblTo(0), lblFrm(1);
	ins  >> nItrs >> lblTo >> lblFrm >> ndif;
	(cout<<"{ "<<" nItrs:"<<nItrs<<"; "<<lblFrm<<" --> "<<lblTo<< "; ndif: "<<ndif<<";  ").flush();

	vImg.growBox(2); cout<<endl;
	for (int i=0; i<nItrs; ++i) FaceMedGrowToFrom(vImg,T(lblTo),T(lblFrm),ndif);
	vImg.shrinkBox(2);

	(cout<<"};\n").flush();
	return 0;
}

template<typename T> bool delense032( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("nItrs(2) lbl0(0) lbl1(1) nAdj0(10) nAdj1(6)");
	int nItrs(2),  nAdj0(10),  nAdj1(6);    Tint lbl0(0), lbl1(1);
	ins >> nItrs >> lbl0 >> lbl1 >> nAdj0 >> nAdj1;
	(cout<<"{ "<<" nItrs:"<<nItrs<<"; lbls: "<<lbl0<<" "<<lbl1<< "; nAdjThresholds: "<<nAdj0<<" "<<nAdj1<<";  ").flush();

	vImg.growBox(2); cout<<endl;
	voxelImageT<T> vimgo=vImg;
	for (int i=0; i<nItrs; ++i)   vImg.PointMedian032(25,nAdj1,lbl0,lbl1);
	FaceMedGrowToFrom(vImg,T(lbl1),T(lbl0),1);
	FaceMedGrowToFrom(vImg,T(lbl0),T(lbl1),-1);
	for (int i=0; i<2*nItrs; ++i) { vImg.PointMedian032(nAdj0,25,lbl0,lbl1);	FaceMedGrowToFrom(vImg,T(lbl0),T(lbl1),-1); }
	FaceMedGrowToFrom(vImg,T(lbl0),T(lbl1),-3);
	FaceMedGrowToFrom(vImg,T(lbl0),T(lbl1),-1);
	FaceMedGrowToFrom(vImg,T(lbl0),T(lbl1),-1);
	
	FaceMedGrowToFrom(vimgo,T(lbl1),T(lbl0),2);//41 51 -> lbl1 
	FaceMedGrowToFrom(vimgo,T(lbl1),T(lbl0),2);//41 51 -> lbl1 
	forAlliii_(vimgo) if(vimgo(iii)==lbl1) vImg(iii)=lbl1;
	vImg.shrinkBox(2);

	(cout<<"};\n").flush();
	return 0;
}


template<typename T> bool circleOut( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("Axis(x/y/z) X0(=N/2) Y0(=N/2) R(=N/2) outVal(Max)");

	char d='z';  ins >> d;
	int i = std::max<int>(d-'x',0);
	int X0(vImg.size3()[(i+1)%3]/2), Y0(vImg.size3()[(i+2)%3]/2);
	int R((X0+Y0)/2);
	Tint outVal=maxT(T);

	ins >> X0 >> Y0 >> R >> outVal;
	(cout<<"  circleOut: dir="<<d<<",  X0="<<X0 <<"  Y0="<<Y0  <<"  R="<<R<<"  out="<<outVal ).flush();

	circleOut(vImg,X0,Y0,R,d,T(outVal));

	(cout<<".").flush();
	return 0;
}


template<typename T> bool maskWriteFraction( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("not implemented");
	int maskvv(2);
	Tint minIelm(1), maxIelm=std::numeric_limits<T>::max();
	string maskname, outName("maskWriteFraction.txt");
	ins >> maskname >> outName >> maskvv >> minIelm >> maxIelm;
	(cout<<"  maskWriteFraction:  mask:"<<maskname <<"  outName:"<<outName<<"  maskvv:"<<maskvv  <<"  minIelm:"<<minIelm<<"  maxIelm:"<<maxIelm ).flush();

	//maskWriteFraction(vImg,maskname,outName,maskvv,minIelm,maxIelm);

	(cout<<"Error: not implemented.").flush();
	return 0;
}


template<typename T> bool Offset( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("offset(0.,0.,0.)");
	dbl3 offset(0.);  ins >> offset;
	(cout<<"  Offset:"<<offset<<" " ).flush();
	vImg.X0Ch()=offset;
	(cout<<".").flush();
	return 0;
}


template<typename T>  bool keepLargest0( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint(" // sets smaller isolated regions (of value 0) to 254, computationally expensive");	
	keepLargest0(vImg); //! CtrlF:isolated=254
	(cout<<".").flush();
	return 0;
}

template<typename T>  bool growLabel( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("vvalue(255)  nIters(0) ");
	int  vv(255), nIters(0);  ins >> vv >> nIters;
	(cout<<"  growLabel: "<<vv<<" x"<<nIters ).flush();

	for (int i=0; i<=nIters; ++i)  vImg.growLabel(vv);

	(cout<<".").flush();
	return 0;
}


template<typename T>  bool reset( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("param Value // param can be N X0 dX V VN or NdXX0");
	string param;
	ins >>param;
	while(param.size())  {
		if(param=="maxNz")  { 	ins>>maxNz;	 cout<<"maxNz:"<<maxNz<<" "; } // stop reading after maxNz layers
		else if(param=="N")  { 	int3 N{0,0,0}; 	ins>>N; 	vImg.reset(N,T{0}); cout<<"N:"<<N<<" "; 	}
		else if(param[0]=='X')  { 	dbl3 X0; 	ins>>X0; 	vImg.X0Ch()=X0; cout<<"X0:"<<X0<<" "; 	}
		else if(param[0]=='d')  { 	dbl3 dx(1.,-2e9,1.); 	ins>>dx; 
			if(dx[1]<-1e9) { dx[1]=dx[0]; dx[2]=dx[1]; }
			vImg.dxCh()=dx; cout<<"dX:"<<dx<<" "; 	}
		else if(param[0]=='V') // VN
		{ 	int3 N=vImg.size3();
			Tint vv(0.); 	ins>>vv>>N;
			vImg.reset(N,vv);
		}
		else if(param[0]=='N' && param[1]=='d') // NdX or Nd0
		{ 	dbl3 dx(1.,-2e9,1.), X0(0.,0.,0.); int3 N=vImg.size3(); 	ins>>N>>dx>>X0;
			if(dx[1]<-1e9) { dx[1]=dx[0]; dx[2]=dx[0]; }
			vImg.reset(N); cout<<"N:"<<N<<" ";  vImg.dxCh()=dx; cout<<"dX:"<<dx<<" "; 	vImg.X0Ch()=X0; cout<<"X0:"<<X0<<" "; 	}
		else cout<<"reset does not support "<<param<<endl;
		param="";
		ins >>param;
	}
	(cout<<".").flush();
	return 0;
}


template<typename T>  bool operat( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("operation(+-^...) [img2Nam/number]  [shift]");
	string opr=" ", img2Nam;    ins>>opr>>img2Nam;
	operat(vImg,opr[0],img2Nam,ins);
	return 0;
}

template<typename T> 
bool mapFrom( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("image2name  minv maxv");
	Tint minv(0), maxv(maxT(T)); double scale=0, shift(0.5-T(0.5));
	string image2name;
	ins>>image2name>>minv>>maxv;
	ensure(maxv>=minv);
	cout<<"\n{  mapping from image "<<image2name<<", assigning to values originally in range: ["<<minv<<" "<<maxv<<"], by "; if(scale>1e-16) { cout<<shift<<"+"<<scale<<"*"; } cout<<image2name<<endl;
	voxelImageT<T> image2(image2name);

	mapToFrom(vImg,image2,T(minv),T(maxv), scale, shift);

	cout<<" } //mapFrom "<<endl;
	return 0;
}



template<typename T,  enable_if_t<std::is_arithmetic<T>::value, int> = 0> 
std::unordered_map<string,bool(*)( stringstream&, voxelImageT<T>&)>
 namedProcesses()  {
	typedef bool(*ProcessP)( stringstream&  ins, voxelImageT<T>& vImg);
	return std::unordered_map<string,ProcessP>{
		{  ""             , ProcessP(& ignore)},// ProcessP can be removed if using g++
		{  ";"	          , ProcessP(& ignore )},// TODO delete
		{  "skip"	      , ProcessP(& ignore )},
		{  "fillHoles"    , ProcessP(& fillHoles )},
		{  "reset"        , ProcessP(& reset )},
		{  "info"         , ProcessP(& info )},
		{  "rescale"	  , ProcessP(& rescale )},
		{  "pore" 	      , ProcessP(& selectPore )},
		{  "threshold"    , ProcessP(& selectPore )},
		{  "threshold101" , ProcessP(& selectPore )},
		{  "Offset"       , ProcessP(& Offset )},
		{  "redirect"     , ProcessP(& redirect )},
		{  "direction"    , ProcessP(& redirect )},
		{  "crop"         , ProcessP(& cropD )},
		{  "cropD"        , ProcessP(& cropD )},
		{  "cropf"        , ProcessP(& cropf )},
		{  "resample"     , ProcessP(& resampleMean )},
		{  "resampleMean" , ProcessP(& resampleMean )},
		{  "resampleMax"  , ProcessP(& resampleMax )},
		{  "resampleMode" , ProcessP(& resampleMode )},
		{  "resliceZ"     , ProcessP(& resliceZ )},
		{  "rangeTo"      , ProcessP(& replaceRange )},
		{  "replaceRange" , ProcessP(& replaceRange )},
		{  "write"        , ProcessP(& write )},
		{  "write8bit"    , ProcessP(& write8bit )},
		{  "read"         , ProcessP(& read )},
		{  "readAtZ"      , ProcessP(& readAtZ )},
		{  "modeFilter"   , ProcessP(& modeFilter )},
		{  "medianFilter" , ProcessP(& medianFilter )},
		{  "medianX"      , ProcessP(& medianX )},
		{  "FaceMedian06" , ProcessP(& FaceMedian06 )},
		{  "PointMedian032" , ProcessP(& PointMedian032 )},
		{  "faceMedNgrowToFrom"   , ProcessP(& faceMedNgrowToFrom )},
		{  "delense032"   , ProcessP(& delense032 )},
		{  "circleOut"    , ProcessP(& circleOut )},
		{  "growLabel"    , ProcessP(& growLabel )},
		{  "keepLargest0"    , ProcessP(& keepLargest0 )},
		{  "maskWriteFraction",ProcessP(& maskWriteFraction )},
		{  "mapFrom"      , ProcessP(& mapFrom )},
		{  "Paint"        , ProcessP(& Paint )},
		{  "PaintAdd"     , ProcessP(& PaintAdd )},
		{  "PaintBefore"    ,ProcessP(& PaintBefore )},
		{  "PaintAfter"     ,ProcessP(& PaintAfter )},
		{  "PaintAddBefore" ,ProcessP(& PaintAddBefore )},
		{  "PaintAddAfter"  ,ProcessP(& PaintAddAfter )},
		#ifdef LPNG
		{  "sliceToPng"    , ProcessP(& sliceToPng )},
		{  "sliceToPngBW"  , ProcessP(& sliceToPngBW )},
		#endif
		{  "operation"    , ProcessP(& operat )},
		{  "operat"       , ProcessP(& operat )}
	};
}


template<typename T,  enable_if_t<std::is_class<T>::value, int> = 0> 
std::unordered_map<string,bool(*)( stringstream&, voxelImageT<T>&)> namedProcesses()  {
	typedef bool(*ProcessP)( stringstream&  ins, voxelImageT<T>& vImg);
	return std::unordered_map<string,ProcessP>{
		{  ""      , ProcessP(& ignore)},// ProcessP can be removed if using g++
		{  ";"		   , ProcessP(& ignore )}, // TODO delete
		{  "skip"	   , ProcessP(& ignore )},
		{  "Offset"    , ProcessP(& Offset )},
		{  "redirect"  , ProcessP(& redirect )},
		{  "direction" , ProcessP(& redirect )},
		{  "write"     , ProcessP(& write )},
		{  "read"      , ProcessP(& read )},
		{  "circleOut" , ProcessP(& circleOut )}
		};
}



								_end_of_(namespace MCTProcessing)






template<typename T>
class  voxelplugins
{
	public:
	typedef bool(*ProcessP)( stringstream&  inputs, voxelImageT<T>& vImg);

	std::unordered_map<string,ProcessP> key_funs;
	const std::unordered_map<string,ProcessP>& operator()() const { return key_funs; }
	voxelplugins() {
		using namespace MCTProcessing;
		key_funs = MCTProcessing::namedProcesses<T>();
	};




	int process(const InputFile& inp, voxelImageT<T>& img, string nam="") const  { // nam is ignored here
		if(inp.data().size()>2) std::cout<<std::endl;
		for(const auto& ky:inp.data())  {
			auto paer = key_funs.find(ky.first);
			if (paer!=key_funs.end())  {
				(cout<<" "<<ky.first<<": ").flush();
				stringstream ss(ky.second);
				(*(paer->second))(ss, img);
				if(inp.data().size()>2) cout<<endl;
			}
			else  {
				if(ky.first!="end") { cout<<"  stopped executing "+inp.fileName()+" before \""+ky.first+"\"  :/ "; 
												return -1; }
				break;
			}
		}
		return 0;
	};
	int process(const string&  keystr, voxelImageT<T>& img, string nam)	{  return  process(InputFile(keystr,nam,false), img);   }
	/*int process(      istream& keyins, voxelImageT<T>& img, string nam)	{  return  process(InputFile(keyins,nam,false), img);
		//while (true)  {
			//std::streampos begLine = keyins.tellg();
			//string ky;  keyins>>ky;
			//if (keyins.fail()) 	{cout<<"  @"<<keyins.tellg()<<"  "<<nam<<" done."<<endl;  break; }
			//else if (ky[0]=='{' || ky[0]=='}') { keyins.seekg(int(keyins.tellg())-ky.size()-1); continue; }
			//else if (ky[0]=='#' || ky[0]=='\'' || ky[0]=='/' || ky[0]=='%')  keyins.ignore(10000,'\n');  
			//else  {
				//auto paer = key_funs.find(ky);
				//if (paer!=key_funs.end())  {
					//(cout<<" "<<ky<<": ").flush();
					//stringstream ss;  if(keyins.peek()!='\n') keyins.get (*(ss.rdbuf()));
					//(*(paer->second))(ss,img);
					//cout<<endl;
				//}
				//else  {	
					//cout<<"  stopped processing "<<nam<<" before \""<<ky<<"\" :/ "<<endl; 
					//keyins.clear(); keyins.seekg(begLine);
					//return -1;
				//}
		//}	}
		//return 0;
	};*/

};

template<class InpT, typename T>  //! run voxel plugins
 int vxlProcess(const InpT& ins, voxelImageT<T>& img, string nam) {  return voxelplugins<T>().process(ins,img,nam);  }

template<class InpT, typename First=uint8_t, typename... Rest>
 int vxlProcess(const InpT& ins, voxelImageTBase* imgPtr, string nam)  { //! detect type and run voxel plugins
	if(auto img = dynamic_cast<voxelImageT<First>*>(imgPtr))  
		return vxlProcess<InpT,First>(ins,*img,nam);
	else if(sizeof...(Rest))
		return vxlProcess<InpT,Rest...>(ins, imgPtr, nam);
	cout<<"Unknown image type."<<endl;
	return -1;
}
template int vxlProcess<InputFile, unsigned char,unsigned short,short,int,float>(const InputFile& inp, voxelImageTBase* imgPtr, string nam);
template int vxlProcess<string, unsigned char,unsigned short,short,int,float>(string const& ins, voxelImageTBase* imgPtr, string nam);


string VxlKeysHelp(string keyname, string subkey)  {
	//! Query and print MCTProcessing keyword usage messages
	typedef bool(*ProcessP)( stringstream&  ins, voxelImageT<unsigned char>& vImg);

	std::unordered_map<string,ProcessP> key_funs = voxelplugins<unsigned char>()();

	voxelImage vImg;
	stringstream keys;
	if(keyname.size())  {
		auto paer = key_funs.find(keyname);
		if (paer!=key_funs.end())  {
			stringstream ss(subkey.empty()?  "?" : "? "+subkey);
			try                         {  (*(paer->second))(ss, vImg); }
			catch (std::exception &exc) {  std::cerr <<keyname<<" KeyHelp not implemented:" << exc.what() << endl; }
			catch (...)                 {  std::cerr <<keyname<<" KeyHelp not implemented:" << endl; }
			return ss.str();
		}
		else
			cout<<" Error: no such keyword "<<keyname<<endl; 
		keys<<"//!-*- C -*- keywords:\n";
		for(const auto& proc:key_funs) 	keys<<proc.first<<"\n";
		keys<<" Error: no such keyword "<<keyname<<"\n\n";
	}
	else  {
		std::vector<std::pair<string,ProcessP>> keyfuns(key_funs.begin(), key_funs.end());
		std::sort(keyfuns.begin(), keyfuns.end());
		for(const auto& proc:keyfuns)  if(proc.first.size()>1)  {
			stringstream ss("?");
			try                         {  (*(proc.second))(ss, vImg); }
			catch (std::exception &exc) {  std::cerr <<proc.first<<" KeyHelp not implemented:" << exc.what() << endl; }
			catch (...)                 {  std::cerr <<proc.first<<" KeyHelp not implemented:" << endl; }
			keys<<"\t"<<proc.first<<": \t"<<ss.str()<<"\n\n";
		}
	}
	return keys.str();
}


template<typename T>
void voxelImageT<T>::readFromHeader(istream& hdrFile, const string& hdrNam, int procesKeys)  {
	//! read image from file header, format detected based on image extension
	auto& vImg=*this; string inputName;

	int3 nnn(0,0,0);
	string BinaryData="XXX", flipSigByt="False";
	bool X0read=false, dxread=false, autoUnit=true; //auto unit only applies to .mhd format
	double unit_=1.;
	int nSkipBytes(0);
	if (hasExt(hdrNam,4,".mhd") || hasExt(hdrNam,3,".py"))	{
		cout<<" mhd:"<<hdrNam<<": "<<endl;
		while (true)  {
			std::streampos begLine = hdrFile.tellg();
			string ky, tmp;   hdrFile>>ky>>tmp;
			stringstream ss;  if(hdrFile.peek()!='\n') hdrFile.get (*(ss.rdbuf()));
			if (hdrFile.fail()) break;
			if (ky=="ObjectType")  {  ss>> tmp;  if (tmp != "Image") cout<<"  Warning: ObjectType != Image :="<<tmp<<endl;	}
			else if (ky=="NDims")  {  ss>> tmp;  if (tmp != "3"    ) cout<<"  Warning: NDims != 3 :="<<tmp<<endl;	}
			else if (ky=="ElementType")  { ss>> tmp;  if ((tmp != "MET_UCHAR") && (sizeof(T)==1)) cout<<"  Warning: ElementType != MET_UCHAR :="<<tmp<<endl; 	}
			else if (ky=="Offset")       { ss>> vImg.X0_;   cout<<"  X0: "<<  vImg.X0_<<",  ";  X0read=true; }
			else if (ky=="ElementSize" || ky=="ElementSpacing")  {  ss>> vImg.dx_;  cout<<"  dX: "<<vImg.dx_<<",  ";  dxread=true;	}
			else if (ky=="DimSize")                              {  ss>> nnn;  nnn.z=std::min(nnn.z,maxNz);  cout<<"  Nxyz: "<<nnn<<",  ";	}
			else if (ky=="ElementDataFile")  {  if (inputName.empty()) ss>> inputName;
				size_t is=hdrNam.find_last_of("\\/");
				if (is<hdrNam.size() && inputName[0]!='/' &&  inputName[1]!=':') inputName=hdrNam.substr(0,is+1)+inputName;
				cout<<"  Img: "<<inputName<<",	";
			}
			else if (ky=="BinaryData")  {  ss>> BinaryData;     cout<<"  BinaryData: "<<BinaryData<<"	"<<endl; }
			else if (ky=="Unit")        {  ss>> unit_;  autoUnit=false;   cout<<"  Unit, OneMeter: "<<unit_<<endl; 	}
			else if (ky=="HeaderSize")  {  ss>> nSkipBytes;         cout<<"  Ski pHeaderSize: "<<nSkipBytes<<endl;	}
			else if (ky=="OutputFormat" || ky=="DefaultImageFormat" )  {  if(tmp=="=") ss>> tmp;  cout<<"  OutputFormat: "<<tmp<<", suffix:"<<imgExt(tmp)<<"	"<<endl; }///. sets suffix+format
			else if (ky=="BinaryDataByteOrderMSB" || ky=="ElementByteOrderMSB")  {  ss>> flipSigByt; }
			else if (ky!="CompressedData" &&  ky!="CompressedDataSize" &&  ky!="TransformMatrix" &&
					 ky!="ElementNumberOfChannels" && ky!="CenterOfRotation" && ky!="AnatomicalOrientation" && ky!="AnatomicalOrientation")  {
				hdrFile.clear();  hdrFile.seekg(begLine);
				(cout<<"; ").flush();
				break;
			}
		}
		cout<<endl;

	}
	#ifdef TIFLIB
	else if (hasExt(hdrNam,".tif"))  {  readTif(vImg, hdrNam);  return;  }
	#endif
	else if (hasExt(hdrNam,".am"))	{
		inputName=hdrNam;
		procesKeys=0;
	}
	else if (hasExt(hdrNam,7,".raw.gz") || hasExt(hdrNam,4,".raw") || hasExt(hdrNam,4,".dat"))  { // detect size and voxel size from image name.
		string 
		data=replaceFromTo(replaceFromTo(replaceFromTo(replaceFromTo(replaceFromTo(
									hdrNam,".gz$",""), ".raw$",""), ".dat$",""),"__","\n"),"_"," ");
		data=replaceFromTo(replaceFromTo(replaceFromTo(data,"voxel",""),"size",""),"um","\n");
		data=regex_replace(data,regex("( [0-9][0-9]*)c"), " $1 $1 $1 ", regex_constants::format_first_only);
		data=regex_replace(data,regex("( [0-9][0-9]*)[ x]*([0-9][0-9]*)[ x]*([0-9][0-9]* )"), 
		                                        "\n   reset Nd0 $1 $2 $3 ", regex_constants::format_first_only);
		data=regex_replace(data,regex("^[^\n]*\n"), "", regex_constants::format_first_only);
		data=regex_replace(data,regex("\n|($)"),"\n   read "+hdrNam+"\n", regex_constants::format_first_only);
		for(auto&da:data)  { if(da=='p') da='.'; else if(da=='\n') break; }
		cout<<"  Keywords: {\n"<<data<<"  }"<<endl;
		vxlProcess(data,vImg,hdrNam);
		procesKeys=0;
	}
	else if (hasExt(hdrNam,7,"_header"))  {
		cout<<" (depricated) _header:"<<hdrNam<<","<<endl;

		char tmpc;
		for (int i=0; i<8; ++i)   hdrFile>>tmpc, cout<<tmpc;  //ignore the first 8 characters (ascii 3uc)

		if (hasExt(hdrNam,7,"_header"))  inputName=hdrNam.substr(0,hdrNam.size()-7);
		hdrFile>>nnn >> vImg.dx_ >>	vImg.X0_ ;
		cout<<"\n Nxyz: "<<nnn<<"    dX: "<< vImg.dx_<<"   X0: "<< vImg.X0_ <<" um"<< endl;
		if (!hdrFile)	 { cout<<"   Incomplete/bad header name. Aborting"<<endl; exit(-1); }
	}
	else  alert("Unknown (header) file type: "+hdrNam,-1); // exit

	if(nnn.z) vImg.reset(nnn);
	if( !inputName.empty() && inputName!="NO_READ" && procesKeys!=2)  {
	  if (hasExt(inputName,4,".tif"))  {
			dbl3 dx=vImg.dx_, X0=vImg.X0_;
			bool readingImage = vImg.readBin(inputName);
			assert(readingImage);
			if(X0read) vImg.X0_=X0;
			if(dxread) vImg.dx_=dx;
	  }
	  else if ((hasExt(inputName,4,".raw") && BinaryData!="False") || BinaryData=="True")   {
			bool readingImage = vImg.readBin(inputName, nSkipBytes);
			assert(readingImage);
	  }
	  else if (hasExt(inputName,3,".am"))    {
			int RLECompressed;
			dbl3 dx=vImg.dx_, X0=vImg.X0_;
			getAmiraHeaderSize(inputName, nnn,vImg.dx_,vImg.X0_,nSkipBytes,RLECompressed);
			bool readingImage = vImg.readBin(inputName, nSkipBytes);
			assert(readingImage);
			if(X0read) vImg.X0_=X0;
			if(dxread) vImg.dx_=dx;
	  }
	  else if (hasExt(inputName,7,".raw.gz"))   {
			bool readingImage = vImg.readBin(inputName);
			assert(readingImage);
	  }
	  else   {
		std::ifstream in(inputName);  assert(in);
		if(nSkipBytes) in.ignore(nSkipBytes);
		vImg.voxelField<T>::readAscii(in);
	  }
	}

	if(flipSigByt=="True") {
		cout<<"  flipEndian "<<endl;
		flipEndian(vImg);	}

	if(autoUnit  && vImg.dx_[0]>0.02)	{ //&& dxread
		cout<<"   dx="<<vImg.dx_[0]<<"(>0.02 -> assuming unit is um), ";
		unit_ = 1e-6;
	}
	vImg.dx_*=unit_;
	vImg.X0_*=unit_;
	if(abs(unit_-1.)>epsT(float)) cout<<"  unit= "<<unit_<<" => dx= "<<vImg.dx_<<", X0= "<<vImg.X0_<<endl;


	if (procesKeys) voxelplugins<T>().process(InputFile(hdrFile,hdrNam),vImg);

}




template void voxelImageT<unsigned char>::readFromHeader(istream&,	const string&, int );
template void voxelImageT<unsigned short>::readFromHeader(istream&,	const string&, int );
template void voxelImageT<int>::readFromHeader(istream&,	const string&, int );
template void voxelImageT<float>::readFromHeader(istream&,	const string&, int );
template void voxelImageT<double>::readFromHeader(istream&,	const string&, int );
template void voxelImageT<float3>::readFromHeader(istream&,	const string&, int );






std::unique_ptr<voxelImageTBase> readImage(string hdrNam,	int procesKeys)  {
	//! read or create image
	using namespace std;
	(cout<<"voxelImage \""<<hdrNam<<"\": ").flush();
	if (hasExt(hdrNam,".am"))  {
		string vtype = getAmiraDataType(hdrNam);
		cout<<"reading '"<<vtype<<"'s from .am file"<<endl;

		if (vtype=="int")       return make_unique<voxelImageT<int>>(hdrNam,0);
		if (vtype=="short")     return make_unique<voxelImageT<short>>(hdrNam,0);
		if (vtype=="ushort")    return make_unique<voxelImageT<unsigned short>>(hdrNam,0);
		if (vtype=="byte")      return make_unique<voxelImageT<unsigned char>>(hdrNam,0);

		cout<<"  Error: data type "<<vtype<<" not supported, when reading "<<hdrNam<<endl;
		exit(-1);
	}

	#ifdef TIFLIB
	if (hasExt(hdrNam,".tif"))  return readTif(hdrNam);
	#endif

	string typ;
	std::ifstream hdrFile(hdrNam); // header file
	if(!hdrFile)  
	{
		ensure(hdrNam.size()<4 || hdrNam[hdrNam.size()-4]!='.', "can not open header file '"+hdrNam+"', pwd: "+getpwd(), -1);
		typ = hdrNam; hdrNam="NO_READ";
	}
	else if (hasExt(hdrNam,4,".mhd"))  {
		while (true)  {
			string ky;  hdrFile>>ky;
			stringstream ss;
			if(hdrFile.peek()!='\n') hdrFile.get (*(ss.rdbuf()));
			if (hdrFile.fail()) {  cout<<"\n\n\nWarning: readImage, 'ElementType =' not set in "<<hdrNam<<endl; break; }
			if (ky == "ElementType")  {  ss >> typ >> typ;  break; }
		}
	}
	hdrFile.close();

	if (typ=="MET_UCHAR")        return make_unique<voxelImageT<unsigned char>>(hdrNam, procesKeys);
	if (typ=="MET_CHAR")         return make_unique<voxelImageT<char>>          (hdrNam, procesKeys);
	if (typ=="MET_USHORT")       return make_unique<voxelImageT<unsigned short>>(hdrNam, procesKeys);
	if (typ=="MET_SHORT")        return make_unique<voxelImageT<short>>         (hdrNam, procesKeys);
	if (typ=="MET_UINT")         return make_unique<voxelImageT<unsigned int>>  (hdrNam, procesKeys);
	if (typ=="MET_INT")          return make_unique<voxelImageT<int>>           (hdrNam, procesKeys);
	if (typ=="MET_FLOAT")        return make_unique<voxelImageT<float>>         (hdrNam, procesKeys);
	if (typ=="MET_DOUBLE")       return make_unique<voxelImageT<double>>        (hdrNam, procesKeys);
	//if (typ=="MET_FLOAT_ARRAY")  return make_unique<voxelImageT<float3>>        (hdrNam, procesKeys);
	//if (typ=="MET_DOUBLE_ARRAY") return make_unique<voxelImageT<dbl3>>          (hdrNam, procesKeys);

	return                              make_unique<voxelImage>(hdrNam, procesKeys);

}

