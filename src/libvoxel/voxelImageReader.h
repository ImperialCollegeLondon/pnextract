/*-------------------------------------------------------------------------*\
You can redistribute this code and/or modify this code under the
terms of the GNU General Public License (GPL) as published by the
Free Software Foundation, either version 3 of the License, or (at
your option) any later version. see <http://www.gnu.org/licenses/>.

This file is part of voxelImage library, a C++ template library  
developed by Ali Qaseminejad Raeini for handelling 3D raw images.


Please see our website for relavant literature making use of this code:
http://www3.imperial.ac.uk/earthscienceandengineering/research/perm/porescalemodelling

For further information please contact us by email:
Ali Q Raeini: a.q.raeini@imperial.ac.uk

\*-------------------------------------------------------------------------*/

#include <sstream>
#include "shapeToVoxel.h"

#ifndef GLOBALS_SkipH
#define Dbg_SkipH
#include "globals.h"
#undef  Dbg_SkipH
#endif

namespace MCTProcessing
{

template<typename T> bool ignore( std::stringstream & ins, voxelImageT<T>& vImg)  {
	if (ins.good() && vImg.nx()==-1) std::cout<<" ";
	
	return true;
}

template<typename T> bool fillHoles( std::stringstream & ins, voxelImageT<T>& vImg)  {
	unsigned int maxHoleSize;
	ins>>maxHoleSize;

	std::cout<<"fillHoles: eliminating isolated rocks/pores; maxHoleSize:" <<maxHoleSize<<" (default is 2) "<<std::endl;
	vImg.fillHoles(maxHoleSize);

	vImg.FaceMedian06(1,5);
	//vImg.FaceMedian07(2,5);
	//vImg.FaceMedian07(2,5);
	return true;
}

template<typename T> bool info( std::stringstream & ins, voxelImageT<T>& vImg)  {
	vImg.printInfo();
	return true;
}

template<typename T> bool selectPore( std::stringstream & ins, voxelImageT<T>& vImg)  {
	std::cout<<"  converting to binary (0 and 1):"<<std::endl
		 <<"  selecting pore (->0) with values between:";
	unsigned int  thresholdMin=0,thresholdMax=0;
	ins>>thresholdMin;
	ins>>thresholdMax;

	std::cout<<" "<<int(thresholdMin)<<"  and "<<int(thresholdMax)<<"  inclusive."<<std::endl;
	vImg.threshold101(thresholdMin,thresholdMax);
	return true;
}

template<typename T,  enable_if_t<std::is_arithmetic<T>::value, int> = 0> 
bool rescale( std::stringstream & ins, voxelImageT<T>& vImg)  {
	(std::cout<<"  rescaling voxel values to [ ").flush();
	unsigned int  thresholdMin=0,thresholdMax=0;
	ins>>thresholdMin;
	ins>>thresholdMax;

	(std::cout<<thresholdMin<<", "<<thresholdMax<<" ]    ").flush();
	rescale(vImg,T(thresholdMin),T(thresholdMax));
	(std::cout<<".").flush();
	return true;
}

template<typename T> bool growPore( std::stringstream & ins, voxelImageT<T>& vImg)  {
		std::cout<<"  growing voxels:"<<std::endl;
		int voxelValueTogrow; ins>>voxelValueTogrow;
 		char growingAlgorithm; ins>>growingAlgorithm;

		while (ins.good())		  // loop while extraction from file is possible
		{
			if (growingAlgorithm!='f')
			{
				if(voxelValueTogrow==0)
					vImg.growPore();
				else if (voxelValueTogrow==0)
					vImg.shrinkPore();
				else
				{
					std::cerr<<"growing is only implemented for binary images: "<<
					"selected voxel value to grow is "<<voxelValueTogrow << ", which is not acceptable"<<std::endl;
					return false;
				}
			} else {
				std::cerr<<"selected growing algorithm: "<<growingAlgorithm<<
				" the only implemented algorithm is f which stands for faceGrowing"<<std::endl;
				return false;
			}

			ins>>voxelValueTogrow;
			ins>>growingAlgorithm;
		}
		std::cout<<" done"<<std::endl;
		return true;
}


template<typename T> bool resampleMean( std::stringstream & ins, voxelImageT<T>& vImg)  {
	double nResample=1;
		ins>>nResample, std::cout<<__FUNCTION__<<" factor: "<<nResample<<std::endl;
		vImg = resampleMean(vImg,nResample);
		return true;
}


template<typename T> bool resampleMax( std::stringstream & ins, voxelImageT<T>& vImg)  {
	double nResample=1;
		ins>>nResample, std::cout<<__FUNCTION__<<" factor: "<<nResample<<std::endl;
		vImg = resampleMax(vImg,nResample);
		return true;
}

template<typename T> bool resliceZ( std::stringstream & ins, voxelImageT<T>& vImg)  {
	double nResample=1;
		ins>>nResample, std::cout<<__FUNCTION__<<" factor: "<<nResample<<std::endl;
		vImg = resliceZ(vImg,nResample);
		return true;
}

template<typename T> bool resampleMode( std::stringstream & ins, voxelImageT<T>& vImg)  {
	double nResample=1;
		ins>>nResample, std::cout<<__FUNCTION__<<" factor: "<<nResample<<std::endl;
		vImg = resampleMode(vImg,nResample);
		return true;
}

template<typename T> bool redirect( std::stringstream & ins, voxelImageT<T>& vImg)  {
		char direction;
		ins>>direction;
		(std::cout<<direction<<", swapping x and "<<direction<<" directions").flush();

		vImg.rotate(direction);
		std::cout<<std::endl;
		return true;
}

template<typename T> bool replaceRange( std::stringstream & ins, voxelImageT<T>& vImg)  {
	int  thresholdMin(0),thresholdMax(0); ///. Warning don't use T, uchar wont work
	ins >> thresholdMin >> thresholdMax;

	int  value=(thresholdMin+thresholdMax)/2; ///. Warning don't use T, uchar wont work
	ins >> value;

	std::cout<<" Replacing range  ["<<thresholdMin<<"  "<<thresholdMax<<"] with "<<value<<";   ";
	replaceRange(vImg,T(thresholdMin),T(thresholdMax),T(value));
	(std::cout<<".").flush();
	return true;
}

//template<typename T> bool crop( std::stringstream & ins, voxelImageT<T>& vImg)
//{
	//int cropBegin[3], cropEnd[3];

	//std::cout<<"Crop:   ";
	//for (int i=0; i<3;++i)   ins>>cropBegin[i] >>cropEnd[i],  std::cout<<cropBegin[i]<<' '<<cropEnd[i]<<"	";
	//std::cout<<' '<<std::endl;

	//vImg.cropD(cropBegin,cropEnd);
	//return true;
//}


template<typename T> bool cropD( std::stringstream & ins, voxelImageT<T>& vImg)  {
	int3 cropBegin(0,0,0), cropEnd=vImg.size3();
	int nLayers(0); int value(1);
	std::cout<<"cropD:   ";
	ins>>cropBegin[0] >>cropBegin[1] >>cropBegin[2];  std::cout<<" "<<cropBegin[0] <<" "<<cropBegin[1] <<" "<<cropBegin[2]<<" -- ";
	ins>>cropEnd[0] >>cropEnd[1] >>cropEnd[2];		std::cout<<cropEnd[0] <<" "<<cropEnd[1] <<" "<<cropEnd[2]<<"  +  ";;
	ins >> nLayers >> value;
	std::cout<<nLayers<<" layers of "<<value<<std::endl;
	vImg.cropD(cropBegin,cropEnd,nLayers,value);
	return true;
}

template<typename T> bool write( std::stringstream & ins, voxelImageT<T>& vImg)  {
	std::string outName("dump.tif");
	ins >> outName;
	vImg.write(outName);
	(std::cout<<".").flush();
	return true;
}

template<typename T> bool writeUchar( std::stringstream & ins, voxelImageT<T>& vImg)  {
	std::string outName("dump.tif");
	ins >> outName;
	double minv=-0.5, maxv=255.0;
	ins>>minv>>maxv;
	double delv=255.499999999/(maxv-minv);
	(std::cout<<minv<<" "<<maxv).flush();
	voxelImageT<unsigned char> voxels(vImg.size3(),vImg.dx(),vImg.X0(),255);
	forAlliii_(voxels) voxels(iii)=std::max(0,std::min(255,int(delv*(vImg(iii)-minv))));
	voxels.write(outName);
	(std::cout<<".").flush();
	return true;
}

template<typename T> bool read( std::stringstream & ins, voxelImageT<T>& vImg)  {
	int3 nnn = vImg.size3();
	std::string fnam;
	ins>>fnam;
	std::cout<<"  reading from  image "<<fnam<<std::endl;
	if(fnam.size()>4)
	{
		if ( hasExt(fnam,4,".tif") || hasExt(fnam,7,".raw.gz") || hasExt(fnam,4,".raw") )
		{
			  vImg.reset(nnn,0);
			  vImg.readBin(fnam);
		}
		else vImg.readFromHeader(fnam,0);
	}
	return true;
}
template<typename T> bool readAtZ( std::stringstream & ins, voxelImageT<T>& vImg)  {// used to stitch images
	int3 nnn = vImg.size3();
	size_t iSlic=0;
	std::string fnam;
	ins>>fnam>>iSlic;
	std::cout<<"  reading from  image "<<fnam<<", assigning to slices after "<<iSlic<<std::endl;
	voxelImageT<T> img(fnam);
	ensure(img.nx()==nnn.x);	ensure(img.ny()==nnn.y);
	std::copy(img.begin(),img.end(),vImg.begin()+iSlic*nnn[0]*nnn[1]);
	return true;
}



template<typename T> bool medianFilter( std::stringstream & ins, voxelImageT<T> & vImg)  {
	int nIterations(1);
	ins >> nIterations;
	(std::cout<<"  median Filter, nIterations: "<<nIterations).flush();
	vImg.growBox(2);
	for (int i=0; i<nIterations; ++i)
	{
		vImg=median(vImg);
	}
	vImg.shrinkBox(2);
	(std::cout<<".").flush();
	return true;
}

template<typename T> bool modeFilter( std::stringstream & ins, voxelImageT<T> & vImg)  {
	int nIterations(1), nMinNeis(2);
	ins >> nIterations >> nMinNeis;
	(std::cout<<"  mode Filter, nIterations: "<<nIterations<<"  nMinNeis"<<nMinNeis).flush();
	vImg.growBox(2);
	for (int i=0; i<nIterations; ++i)
	{
		modeNSames(vImg,nMinNeis,true);
	}
	vImg.shrinkBox(2);
	(std::cout<<".").flush();
	return true;
}
template<typename T> bool medianX( std::stringstream & ins, voxelImageT<T> & vImg)  {
	int nIterations(1);
	ins >> nIterations;
	(std::cout<<"  median Filter, nIterations: "<<nIterations).flush();
	for (int i=0; i<nIterations; ++i)
	{
		vImg=medianx(vImg);
	}
	(std::cout<<".").flush();
	return true;
}


template<typename T> bool FaceMedian06( std::stringstream & ins, voxelImageT<T> & vImg)  {
	if(ins.peek()=='?') { ins.str("nAdj0(2), nAdj1(4),  nIterations(1)"); return true; }
	int nAdj0(2), nAdj1(4),  nIterations(1);
	ins >> nAdj0>> nAdj1>> nIterations;
	(std::cout<<"  FaceMedian06: "<<nAdj0<<" "<<nAdj1<<" "<<nIterations<<"     ").flush();
	vImg.growBox(2);
	for (int i=0; i<nIterations; ++i)
	{
		vImg.FaceMedian06(nAdj0,nAdj1);
	}
	vImg.shrinkBox(2);
	(std::cout<<".").flush();
	return true;
}



template<typename T> bool PointMedian032( std::stringstream & ins, voxelImageT<T> & vImg)  {
	int nItrs(1),  nAdj(11), lbl0(0), lbl1(1);
	ins >> nItrs>> nAdj>> lbl0>> lbl1;
	(std::cout<<"  PointMedian032, "<<" nItrs:"<<nItrs<< "; nAdjThreshold "<<nAdj<<"  lbl0:"<<lbl0<<"  lbl1;"<<lbl1<<"s    ").flush();
	//vImg.growBox(2);

	for (int i=0; i<nItrs; ++i)  vImg.PointMedian032(nAdj,nAdj,lbl0,lbl1);

	//vImg.shrinkBox(2);
	(std::cout<<".").flush();
	return true;
}


template<typename T> bool delense032( std::stringstream & ins, voxelImageT<T> & vImg)  {
	int nItrs(2),  nAdj0(10),  nAdj1(6); 
	Tint lbl0(0), lbl1(1);
	ins >> nItrs>> lbl0>> lbl1>> nAdj0>>nAdj1;
	(std::cout<<"{ "<<" nItrs:"<<nItrs<<"; lbls: "<<lbl0<<" "<<lbl1<< "; nAdjThresholds: "<<nAdj0<<" "<<nAdj1<<";  ").flush();

	vImg.growBox(2); std::cout<<endl;
	voxelImageT<T> vimgo=vImg;
	for (int i=0; i<nItrs; ++i)   vImg.PointMedian032(25,nAdj1,lbl0,lbl1);
	FaceMedGrowTo(vImg,T(lbl1),T(lbl0),1);
	FaceMedGrowTo(vImg,T(lbl0),T(lbl1),-1);
	for (int i=0; i<2*nItrs; ++i) { vImg.PointMedian032(nAdj0,25,lbl0,lbl1);	FaceMedGrowTo(vImg,T(lbl0),T(lbl1),-1); }
	FaceMedGrowTo(vImg,T(lbl0),T(lbl1),-3);
	FaceMedGrowTo(vImg,T(lbl0),T(lbl1),-1);
	FaceMedGrowTo(vImg,T(lbl0),T(lbl1),-1);
	
	FaceMedGrowTo(vimgo,T(lbl1),T(lbl0),2);//41 51 -> lbl1 
	FaceMedGrowTo(vimgo,T(lbl1),T(lbl0),2);//41 51 -> lbl1 
	forAlliii_(vimgo) if(vimgo(iii)==lbl1) vImg(iii)=lbl1;
	vImg.shrinkBox(2);

	(std::cout<<"};\n").flush();
	return true;
}


template<typename T> bool circleOut( std::stringstream & ins, voxelImageT<T> & vImg)  {

	char d='z';
	ins >> d;
	int i = std::max<int>(d-'x',0);
	int X0(vImg.size3()[(i+1)%3]/2), Y0(vImg.size3()[(i+2)%3]/2);
	int R((X0+Y0)/2);

	ins >> X0 >> Y0 >> R;
	(std::cout<<"  circleOut: dir="<<d<<",  X0="<<X0 <<"  Y0="<<Y0  <<"  R="<<R ).flush();

	circleOut(vImg,X0,Y0,R,d);

	(std::cout<<".").flush();
	return true;
}


template<typename T> bool maskWriteFraction( std::stringstream & ins, voxelImageT<T> & vImg)  {
	int maskvv(2);
	T minIelm(1), maxIelm=std::numeric_limits<T>::max();
	std::string maskname, outName("maskWriteFraction.txt");
	ins >> maskname >> outName >> maskvv >> minIelm >> maxIelm;
	(std::cout<<"  maskWriteFraction:  mask:"<<maskname <<"  outName:"<<outName<<"  maskvv:"<<maskvv  <<"  minIelm:"<<minIelm<<"  maxIelm:"<<maxIelm ).flush();

	maskWriteFraction(vImg,maskname,outName,maskvv,minIelm,maxIelm);

	(std::cout<<".").flush();
	return true;
}


template<typename T> bool Offset( std::stringstream & ins, voxelImageT<T> & vImg)  {
	dbl3 offset;
	ins >> offset;
	(std::cout<<"  Offset:"<<offset<<" " ).flush();
	vImg.X0Ch()=offset;
	(std::cout<<".").flush();
	return true;
}


template<typename T>  bool growLabel( std::stringstream & ins, voxelImageT<T>& vImg)  {
	int  vv(255),nIters(0);
	ins >> vv>>nIters;
	(std::cout<<"  growLabel: "<<vv<<" x"<<nIters ).flush();


	for (int i=0; i<=nIters;++i)
	 vImg.growLabel(vv);

	(std::cout<<".").flush();
	return true;
}


template<typename T>  bool reset( std::stringstream & ins, voxelImageT<T>& vImg)  {
	std::string param;
	ins >>param;
	while(param.size())
	{
		if(param=="N")
		{ 	int3 N{0,0,0}; 	ins>>N; 	vImg.reset(N); cout<<"N:"<<N<<" "; 	}
		else if(param=="X0")
		{ 	dbl3 X0; 	ins>>X0; 	vImg.X0Ch()=X0; cout<<"X0:"<<X0<<" "; 	}
		else if(param=="dx")
		{ 	dbl3 dx(1.,-1e64,1.); 	ins>>dx; 
			if(dx[1]<0) {  dx[1]=dx[0]; dx[2]=dx[1]; }
			vImg.dxCh()=dx; cout<<"dX:"<<dx<<" "; 	}
		else if(param[0]=='V')
		{ 	int3 N=vImg.size3();
			Tint vv(0.0); 	ins>>vv>>N;
			vImg.reset(N,vv);
		}
		else if(param=="Nd0")
		{ 	dbl3 dx(1.,1.,1.), X0(0.,0.,0.); int3 N=vImg.size3(); 	ins>>N>>dx>>X0;
			vImg.reset(N); cout<<"N:"<<N<<" ";  vImg.dxCh()=dx; cout<<"dX:"<<dx<<" "; 	vImg.X0Ch()=X0; cout<<"X0:"<<X0<<" "; 	}
		else cout<<"reset does not support "<<param<<endl;
		param="";
		ins >>param;
	}
	(std::cout<<".").flush();
	return true;
}

template<typename T > // enable_if_t wont work with the old C++11
  bool operation( stringstream & ins, voxelImageT<T> & vxlImage)
{
	char operate=' ';
	std::string image2name,catchImg;
	ins>>operate>>image2name;
	ins>>catchImg;
	if(image2name.empty())
	{
	 cout<<"\n  image ="<<operate<<"image ";
	 switch (operate) {
		case '!':
			forAllvp_(vxlImage) *vp=!(*vp);	break;
		//case '~':
			//forAllvp_(vxlImage) *vp=~(*vp);	break;
		case '-':
			forAllvp_(vxlImage) *vp=-(*vp);	break;
		default:   cout<<"\n  Image  "<<operate<<"= !!!not supported!!! ";
	 }
	}else{
	 if (isdigit(image2name[0])&&isdigit(image2name.back()))
	 {
		T   ii; ii=strTo<Tint>(image2name);
		double dd; dd=strTo<double>(image2name);
		cout<<"  image "<<operate<<"= "<< ii<<" ";
		switch (operate) {
			case '=':
				forAllvp_(vxlImage) { (*vp)=ii; }	break;
			case '+':
				forAllvp_(vxlImage) { (*vp)+=ii; }	break;
			case '-':
				forAllvp_(vxlImage) { (*vp)-=ii; }	break;
			//case '&':
				//forAllvp_(vxlImage) (*vp)=(*vp)&ii;	break;
			//case '|':
				//forAllvp_(vxlImage) (*vp)=(*vp)|ii;	break;
			//case '%':
				//forAllvp_(vxlImage) (*vp)=(*vp)%ii;	break;
			case '*':
				forAllvp_(vxlImage) { (*vp)*=dd; }	break;
			case '/':
				forAllvp_(vxlImage) { (*vp)/=dd; }	break;
			case 'b':
				forAllvp_(vxlImage) { (*vp)=max(ii,(*vp)); }	break;
			case 'e':
				forAllvp_(vxlImage) { (*vp)=min(ii,(*vp)); }	break;
			default:   cout<<"\n  Image  "<<operate<<"= !!!not supported!!! ";
		}
	 }else{
		cout<<"\n  image  "<<operate<<"= "<<image2name<<" ";
		voxelImageT<T> image2;
		if (hasExt(image2name,4,".mhd"))
		{
			  image2.reset(vxlImage.size3(),0);
			  image2.readBin(image2name);
		}
		else image2.readFromHeader(image2name);

		switch (operate) {
			case '+':
				forAlliii_(vxlImage) { vxlImage(iii)=min(Tint(vxlImage(iii))+image2(iii),Tint(maxT(T))); }	break;
			case '-':
				//forAlliii_(vxlImage) vxlImage(iii)=min(max(int(vxlImage(iii))+128-image2(iii),0),maxT(T));	break;
				forAlliii_(vxlImage) { vxlImage(iii)=min(max(Tint(vxlImage(iii))-image2(iii),Tint(0)),Tint(maxT(T))); }	break; //! Analyse needs this
			//case '&':
				//forAlliii_(vxlImage) vxlImage(iii)=vxlImage(iii)&image2(iii);	break;
			//case '|':
				//forAlliii_(vxlImage) vxlImage(iii)=vxlImage(iii)|image2(iii);	break;
			//case '%':
				//forAlliii_(vxlImage) vxlImage(iii)=vxlImage(iii)%image2(iii);	break;
			case '*':
				forAlliii_(vxlImage) { vxlImage(iii)=min(Tint(vxlImage(iii))*image2(iii),Tint(maxT(T))); }	break;
			//case '/':
				//forAlliii_(vxlImage) vxlImage(iii)/=image2(iii);	break;
			default:   cout<<"\n  Image  "<<operate<<"= !!!not supported!!! ";
		}
	 }
	}


	return true;
}





template<typename T > 
  bool mapFrom( std::stringstream & ins, voxelImageT<T> & vImg)  {
	Tint minv(0),maxv(255);
	std::string image2name,catchImg;
	ins>>image2name>>minv>>maxv;

	std::cout<<"\n{  mapping from image "<<image2name<<", assigning to values originally in range: ["<<minv<<" "<<maxv<<"]"<<std::endl;
	voxelImageT<T> image2(image2name);

	mapToFrom(vImg,image2,T(minv),T(maxv));

	std::cout<<" } //mapFrom "<<std::endl;
	return true;
}

template<typename T> std::unordered_map<std::string,bool(*)( std::stringstream&, voxelImageT<T>&)>
 namedProcesses()
{
	typedef bool(*ProcessP)( std::stringstream&  ins, voxelImageT<T>& vImg);
	return std::unordered_map<std::string,ProcessP>{
		{  ""             , ProcessP(& ignore)},// ProcessP can be removed if using g++
		{  ";"		      , ProcessP(& ignore )},
		{  "fillHoles"    , ProcessP(& fillHoles )},
		{  "reset"	      , ProcessP(& reset )},
		{  "info"	      , ProcessP(& info )},
		{  "rescale"	  , ProcessP(& rescale )},
		{  "pore"		  , ProcessP(& selectPore )},
		{  "threshold"    , ProcessP(& selectPore )},
		{  "threshold101" , ProcessP(& selectPore )},
		{  "Offset"       , ProcessP(& Offset )},
		{  "direction"    , ProcessP(& redirect )},
		{  "crop"	      , ProcessP(& cropD )},
		{  "cropD"	      , ProcessP(& cropD )},
		{  "resampleMean" , ProcessP(& resampleMean )},
		{  "resampleMax"  , ProcessP(& resampleMax )},
		{  "resampleMode" , ProcessP(& resampleMode )},
		{  "resliceZ"     , ProcessP(& resliceZ )},
		{  "replaceRange" , ProcessP(& replaceRange )},
		{  "write"        , ProcessP(& write )},
		{  "writeUchar"   , ProcessP(& writeUchar )},
		{  "read"         , ProcessP(& read )},
		{  "readAtZ"      , ProcessP(& readAtZ )},
		{  "modeFilter" , ProcessP(& modeFilter )},
		{  "medianFilter" , ProcessP(& medianFilter )},
		{  "medianX"      , ProcessP(& medianX )},
		{  "FaceMedian06" , ProcessP(& FaceMedian06 )},
		{  "PointMedian032" , ProcessP(& PointMedian032 )},
		{  "delense032" , ProcessP(& delense032 )},
		{  "circleOut"  , ProcessP(& circleOut )},
		{  "growLabel"  , ProcessP(& growLabel )},
		{  "maskWriteFraction"  ,ProcessP(& maskWriteFraction )},
		{  "mapFrom"  ,ProcessP(& mapFrom )},
		{  "shapeToVoxel"  ,ProcessP(& shapeToVoxel )},
		{  "shapeToVoxelAdd"  ,ProcessP(& shapeToVoxelAdd )},
		{  "operation"  ,ProcessP(& operation )},
	};
}



}


template<typename T>
void VxlKeysProcess
(
	std::istream& keyins,
	voxelImageT<T>& vImg,
	const std::string& hdrNam
)
{
	typedef bool(*ProcessP)( std::stringstream&, voxelImageT<T>&);
	std::unordered_map<std::string,ProcessP> key_funs = MCTProcessing::namedProcesses<T>();


	while (true)
	{
		std::streampos begLine = keyins.tellg();
		std::string tmpStr;
		keyins>>tmpStr;
		//bool validKey=false;
		//cout<<tmpStr<<endl;///. keep me
		if (keyins.fail())
		{std::cout<<" Read "<<hdrNam<<":/  "<<keyins.tellg()<<std::endl;  break; }
		else if (tmpStr[0]=='#' || tmpStr[0]=='\'' || tmpStr[0]=='/' || tmpStr[0]=='%')
		{
			keyins.ignore(10000,'\n');
			//validKey=true;
		}
		else
		{
			auto paer = key_funs.find(tmpStr);
			if (paer!=key_funs.end())	
			{
				(std::cout<<" "<<tmpStr<<": ").flush();
				std::stringstream ss;
				if(keyins.peek()!='\n') keyins.get (*(ss.rdbuf()));
				(*(paer->second))(ss,vImg);
				std::cout<<std::endl;
				//validKey=true;
			}
			else
			{	std::cout<<"  read "<<hdrNam<<" util entry \""<<tmpStr<<"\":/ \n"<<std::endl;
				keyins.clear();
				keyins.seekg(begLine);
				break;
			}
		}
	}
}

inline std::string VxlKeysHelp(std::string keyname="", std::string subkey="")
{
	typedef bool(*ProcessP)( std::stringstream&  ins, voxelImageT<unsigned char>& vImg);
	std::unordered_map<std::string,ProcessP> key_funs = MCTProcessing::namedProcesses<unsigned char>();

	std::stringstream keys;
	if(keyname.size())
	{
		auto paer = key_funs.find(keyname);
		if (paer!=key_funs.end())
		{
			voxelImageT<unsigned char> vImg;
			std::stringstream ss;
			ss.str(subkey.empty()?  "?" : "? "+subkey);
			(*(paer->second))(ss, vImg);
			return ss.str();
		}
		else
			std::cout<<" Error: no such keyword "<<keyname<<std::endl; 
		keys<<"//!-*- C -*- keywords:\n";
		for(const auto& proc:key_funs) 	keys<<proc.first<<"\n";
		keys<<" Error: no such keyword "<<keyname<<"\n\n";
	}
	else
		for(const auto& proc:key_funs) 	keys<<proc.first<<"\n";
	//std::cout<<keys.str();
	return keys.str();
}


template<typename T>
void voxelImageT<T>::readFromHeader
(
	std::istream& headerFile,
	const std::string& hdrNam,
	int procesKeys,
	std::string inputName
)
{
	int3 n(0,0,0);
	std::string BinaryData="XXX";
	bool X0read=false, dxread=false, autoUnit=true; //auto unit only applies to .mhd format
	double unit_=1.0;
	int nSkipBytes(0);
	#ifdef TIFLIB
	if (hasExt(hdrNam,4,".tif"))
	{
		(std::cout<<  " reading tif file "<<hdrNam<<" ").flush();
		readTif(*this, hdrNam);
		std::cout<<  " ."<<std::endl;
		return;
	}
	else
	#endif
	if (hasExt(hdrNam,4,".mhd") || hasExt(hdrNam,3,".py"))
	{
		std::cout<<" mhd:"<<hdrNam<<": "<<std::endl;
		while (true)
		{
			std::string tmpStr;
			std::streampos begLine = headerFile.tellg();
			headerFile>>tmpStr;


			//ObjectType = Image
			//NDims = 3
			//Offset = 0 0 0
			//ElementSize = 8 8 8
			//DimSize = 200 225 153
			//ElementType = MET_UCHAR
			//ElementDataFile = Ketton100.raw
			std::stringstream ss;
			if(headerFile.peek()!='\n') headerFile.get (*(ss.rdbuf()));
			if (headerFile.fail()) break;
			std::string tmp;
			if (tmpStr == "ObjectType")  {
				ss >> tmp;  ss >> tmp;  if (tmp != "Image") std::cout<<" Warning: ObjectType != Image :="<<tmp<<std::endl;
			}
			else if (tmpStr == "NDims")  {
				ss >> tmp;  ss >> tmp;  if (tmp != "3") std::cout<<" Warning: NDims != 3 :="<<tmp<<std::endl;
			}
			else if (tmpStr == "ElementType")  {
				ss >> tmp;  ss >> tmp;  if ((tmp != "MET_UCHAR") && (sizeof(T)==1)) std::cout<<" Warning: ElementType != MET_UCHAR :="<<tmp<<std::endl;
			}
			else if (tmpStr == "Offset")  {
				ss >> tmp;  ss>> X0_;   std::cout<<" X0: "<<  X0_<<",  ";  X0read=true;
			}
			else if (tmpStr == "ElementSize" || tmpStr == "ElementSpacing")  {
				ss >> tmp;  ss>> dx_ ;  std::cout<<" dX: "<< dx_<<",  ";  dxread=true;
			}
			else if (tmpStr == "DimSize")  {
				ss >> tmp;  ss>> n;     std::cout<<" Nxyz: "<<n<<",  ";
			}
			else if (tmpStr == "ElementDataFile")  {
				ss >> tmp; if (inputName.empty()) ss >> inputName;
				size_t is=hdrNam.find_last_of("\\/");
				if (is<hdrNam.size() && inputName[0]!='/' &&  inputName[1]!=':') inputName=hdrNam.substr(0,is+1)+inputName;
				std::cout<<"\n ElementDataFile = "<<inputName<<"	";
			}
			else if (tmpStr == "BinaryData")  {
				ss >> tmp; ss >> BinaryData;  std::cout<<" BinaryData = "<<BinaryData<<"	"<<std::endl;
			}
			else if (tmpStr == "OutputFormat" || tmpStr == "DefaultImageFormat" )  {
				std::string ext;  ss >> ext;  if(ext=="=") ss >> ext;  std::cout<<" OutputFormat = "<<ext<<", suffix:"<<imgExt(ext)<<"	"<<std::endl; ///. sets suffix+format
			}
			else if (tmpStr == "Unit")  {
				ss >> tmp; ss >> unit_;      std::cout<<" Unit, OneMeter = "<<unit_<<std::endl;  autoUnit=false;
			}
			else if (tmpStr == "HeaderSize")  {
				ss >> tmp; ss >> nSkipBytes;  std::cout<<"HeaderSize, nSkipBytes = "<<nSkipBytes<<std::endl;
			}
			else if (tmpStr!="BinaryDataByteOrderMSB" && tmpStr!="ElementByteOrderMSB" && tmpStr!="CompressedData" &&  tmpStr!="CompressedDataSize" &&  tmpStr!="TransformMatrix" &&
					 tmpStr!="ElementNumberOfChannels" && tmpStr!="CenterOfRotation" && tmpStr!="AnatomicalOrientation" && tmpStr!="AnatomicalOrientation")  {
				headerFile.clear();
				headerFile.seekg(begLine);
				(std::cout<<" ; ").flush();
				break;
			}

		}
		std::cout<<std::endl;

	}
	else if (hasExt(hdrNam,3,".am"))
	{
		inputName=hdrNam;
		procesKeys=0;
	}
	else
	{
		std::cout<<" (depricated) _header:"<<hdrNam<<","<<std::endl;

		char tmpc;
		for (int i=0; i<8;++i)   headerFile>>tmpc, std::cout<<tmpc;  //ignore the first 8 characters (ascii 3uc)

		if (hasExt(hdrNam,7,"_header"))  inputName=hdrNam.substr(0,hdrNam.size()-7);
		headerFile>>n[0]>>n[1]>>n[2];						// number of variables (dimension of
		std::cout<<"\n Nxyz: "<<n[0]<<" "<<n[1]<<" "<<n[2]<<"   "; std::cout.flush();
		headerFile>>	dx_[0]>>dx_[1]>>dx_[2] ;
		std::cout<<" dX: "<< dx_[0]<<"  "<<dx_[1]<<"  "<<dx_[2]<<"   "; std::cout.flush();
		headerFile>>	X0_[0]>>X0_[1]>>X0_[2] ;
		std::cout<<" X0: "<<  X0_[0]<<"  "<<X0_[1]<<"   "<<X0_[2] <<" um"<< std::endl;
		if (!headerFile)	 { std::cout<<"  Incomplete/bad header name, aborting"<<std::endl; exit(-1);}
		//if (!headerFile)	 { std::cout<<"  Incomplete/bad hdrNam, continuing anyway"<<std::endl; }

	}


	this->reset(n);
	if( !inputName.empty() && inputName!="NO_READ" && procesKeys!=2 )
	{
	  if (hasExt(inputName,4,".tif"))
	  {
			dbl3 dx=dx_, X0=X0_;
			bool readingImage = this->readBin(inputName);
			assert(readingImage);
			if(X0read) X0_=X0;
			if(dxread) dx_=dx;
	  }
	  else if ((hasExt(inputName,4,".raw") && BinaryData!="False") || BinaryData=="True")
	  {
			bool readingImage = this->readBin(inputName, nSkipBytes);
			assert(readingImage);
	  }
	  else if (hasExt(inputName,3,".am"))
	  {
			int RLECompressed;
			dbl3 dx=dx_, X0=X0_;
			getAmiraHeaderSize(inputName, n,dx_,X0_,nSkipBytes,RLECompressed);
			bool readingImage = this->readBin(inputName, nSkipBytes);
			assert(readingImage);
			if(X0read) X0_=X0;
			if(dxread) dx_=dx;
	  }
	  else if (hasExt(inputName,7,".raw.gz"))
	  {
			bool readingImage = this->readBin(inputName);
			assert(readingImage);
	  }
	  else
	  {
		std::ifstream in(inputName.c_str());
		assert(in);
		if(nSkipBytes) in.ignore(nSkipBytes);
		voxelField<T>::readAscii(in);
	  }
	}

	if(autoUnit  && dx_[0]>0.01) //&& dxread
	{
		std::cout<<" Warning: too large dx (="<<dx_[0]<<"), assuming unit is um, ";
		unit_ = 1.0e-6;
	}
	dx_*=unit_;
	X0_*=unit_;
	if(std::abs(unit_-1.0)>epsT(float)) std::cout<<"  unit= "<<unit_<<" => dx= "<<dx_<<", X0= "<<X0_<<std::endl;


	if (procesKeys) VxlKeysProcess(headerFile,*this,hdrNam);

}









// read or create image
inline std::unique_ptr<voxelImageTBase> readImage
(
	std::string hdrNam, //!< headername or image type
	int procesKeys = 1
)
{

	(std::cout<<"voxelImage \""<<hdrNam<<"\": ").flush();
	if (hasExt(hdrNam,3,".am"))
	{
		std::string vtype = getAmiraDataType(hdrNam);
		if (vtype=="int," || vtype=="int")
		{
			std::cout<<"reading int .am file: "<<hdrNam<<std::endl;
			return std::unique_ptr<voxelImageTBase>(new voxelImageT<int>(hdrNam,0));
		}
		if (vtype=="short," || vtype=="short")
		{
			std::cout<<"reading short .am file: "<<hdrNam<<std::endl;
			return std::unique_ptr<voxelImageTBase>(new voxelImageT<short>(hdrNam,0));
		}
		else if (vtype=="ushort," || vtype=="ushort")
		{
			std::cout<<"reading ushort .am file: "<<hdrNam<<std::endl;
			return std::unique_ptr<voxelImageTBase>(new voxelImageT<unsigned short>(hdrNam,0));
		}
		else if (vtype=="byte," || vtype=="byte")
		{
			std::cout<<"reading unsigned byte .am file: "<<hdrNam<<std::endl;
			return std::unique_ptr<voxelImageTBase>(new voxelImageT<unsigned char>(hdrNam,0));
		}
		else 
		{
			std::cout<<"data type "<<vtype<<" not supported, when reading "<<hdrNam<<std::endl;
			exit(-1);
		}
	}

	#ifdef TIFLIB
	if (hasExt(hdrNam,4,".tif"))  return readTif(hdrNam);
	#endif

	std::string typ;
	std::ifstream headerFile(hdrNam.c_str());
	if(!headerFile)  
	{
		if (hdrNam.size()>4 && hdrNam[hdrNam.size()-4]=='.') std::cout<<"\n\n\nError: can not open hdrNam file, "<<hdrNam<<std::endl<<std::endl;
		typ = hdrNam; hdrNam="NO_READ";
	}
	else if (hasExt(hdrNam,4,".mhd"))
	{
		while (true)
		{
			std::string tmpStr;  headerFile>>tmpStr;
			std::stringstream ss;
			if(headerFile.peek()!='\n') headerFile.get (*(ss.rdbuf()));
			if (headerFile.fail()) {  std::cout<<"\n\n\nWarning: readImage, 'ElementType =' not set in "<<hdrNam<<std::endl; break; }
			if (tmpStr == "ElementType")  {  ss >> typ >> typ;  break; }
		}
	}
	headerFile.close();

	if (typ=="MET_UCHAR")
	 { headerFile.close();	return std::unique_ptr<voxelImageTBase>(new voxelImageT<unsigned char>(hdrNam, procesKeys)); }
	if (typ=="MET_CHAR")
	 { headerFile.close();	return std::unique_ptr<voxelImageTBase>(new voxelImageT<char>(hdrNam, procesKeys)); }
	if (typ=="MET_USHORT")
	 { headerFile.close();	return std::unique_ptr<voxelImageTBase>(new voxelImageT<unsigned short>(hdrNam, procesKeys)); }
	if (typ=="MET_SHORT")
	 { headerFile.close();	return std::unique_ptr<voxelImageTBase>(new voxelImageT<short>(hdrNam, procesKeys)); }
	if (typ=="MET_UINT")
	 { headerFile.close();	return std::unique_ptr<voxelImageTBase>(new voxelImageT<unsigned int>(hdrNam, procesKeys)); }
	if (typ=="MET_INT")
	 { headerFile.close();	return std::unique_ptr<voxelImageTBase>(new voxelImageT<int>(hdrNam, procesKeys)); }
	if (typ=="MET_FLOAT")
	 { headerFile.close();	return std::unique_ptr<voxelImageTBase>(new voxelImageT<float>(hdrNam, procesKeys)); }
	if (typ=="MET_DOUBLE")
	 { headerFile.close();	return std::unique_ptr<voxelImageTBase>(new voxelImageT<double>(hdrNam, procesKeys)); }


	return std::unique_ptr<voxelImageTBase>(new voxelImage(hdrNam, procesKeys));

}




template<typename T>
void readConvertFromHeader
(	voxelImageT<T>& vImg,
	std::string hdrNam,
	int procesKeys = 1
)
{
	std::unique_ptr<voxelImageTBase> vImgUptr = readImage(hdrNam,procesKeys);
	voxelImageTBase* imgPtr = vImgUptr.get();
	bool red = false; //read
	if(!red) {auto img = dynamic_cast<voxelImageT<T>*>            (imgPtr);  if(img) { vImg = *img; red=true; } }
	if(!red) {auto img = dynamic_cast<voxelImageT<char>*>         (imgPtr);  if(img) { vImg.resetFrom(*img); red=true; std::cout<<"read into "<<vImg.size3()<<" chars "; } }
	if(!red) {auto img = dynamic_cast<voxelImageT<unsigned char>*>(imgPtr);  if(img) { vImg.resetFrom(*img); red=true; std::cout<<"read into "<<vImg.size3()<<" uchrs "; } }
	if(!red) {auto img = dynamic_cast<voxelImageT<short>*>        (imgPtr);  if(img) { vImg.resetFrom(*img); red=true; std::cout<<"read into "<<vImg.size3()<<" shrts "; } }
	if(!red) {auto img = dynamic_cast<voxelImageT<unsigned short>*>(imgPtr); if(img) { vImg.resetFrom(*img); red=true; std::cout<<"read into "<<vImg.size3()<<" usrts "; } }
	if(!red) {auto img = dynamic_cast<voxelImageT<int>*>         (imgPtr);   if(img) { vImg.resetFrom(*img); red=true; std::cout<<"read into "<<vImg.size3()<<" intgs "; } }
	if(!red) {auto img = dynamic_cast<voxelImageT<unsigned int>*>(imgPtr);   if(img) { vImg.resetFrom(*img); red=true; std::cout<<"read into "<<vImg.size3()<<" uints "; } }
	if(!red) {auto img = dynamic_cast<voxelImageT<float>*>       (imgPtr);   if(img) { vImg.resetFrom(*img); red=true; std::cout<<"read into "<<vImg.size3()<<" flots "; } }
	if(!red) {auto img = dynamic_cast<voxelImageT<double>*>      (imgPtr);   if(img) { vImg.resetFrom(*img); red=true; std::cout<<"read into "<<vImg.size3()<<" dobls "; } }
	if(!red) std::cout<<"\n\ncan not convert image\n\n"<<std::endl;
}

template<class T, typename First=uint8_t, typename... Rest>
int resetFromImageT(voxelImageT<T>& vImg, voxelImageTBase* vxlImgPtr)  {
	auto img = dynamic_cast<voxelImageT<First>*>(vxlImgPtr); 
	if(img) { vImg.resetFrom(*img); return 0; }
	else if(sizeof...(Rest)) return resetFromImageT<T,Rest...>(vImg, vxlImgPtr);
	std::cout<<"Error in resetFromImageT: Unknown image type."<<std::endl;
	return -1;
}

