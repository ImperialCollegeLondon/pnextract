/*-------------------------------------------------------------------------*\
You can redistribute this code and/or modify this code under the
terms of the GNU General Public License (GPL) as published by the
Free Software Foundation, either version 3 of the License, or (at
your option) any later version. see <http://www.gnu.org/licenses/>.


This file is part of voxelImage library, a small c++ template library  
developed by Ali Qaseminejad Raeini for handelling 3D raw images.


Please see our website for relavant literature making use of this code:
http://www3.imperial.ac.uk/earthscienceandengineering/research/perm/porescalemodelling

For further information please contact us by email:
Ali Q Raeini:	a.qaseminejad-raeini09@imperial.ac.uk

\*-------------------------------------------------------------------------*/

#include <sstream>
//~ #include <streambuf>


namespace MCTProcessing
{
template<typename T> bool ignore( std::stringstream & inputs, voxelImageT<T>& vxlImage)
{
	if (inputs.good() && vxlImage.size3()[0]==-1) std::cout<<" ";
	
	return true;
}

template<typename T> bool fillHoles( std::stringstream & inputs, voxelImageT<T>& vxlImage)
{
	unsigned int maxHoleSize;
	inputs>>maxHoleSize;

		std::cout<<"fillHoles: eliminating isolated rocks/pores; maxHoleSize:" <<maxHoleSize<<" (default is 2) "<<std::endl;
		vxlImage.fillHoles(maxHoleSize);

		vxlImage.FaceMedian06(1,5);
		//~ vxlImage.FaceMedian07(2,5);
		//~ vxlImage.FaceMedian07(2,5);
		return true;
}

template<typename T> bool selectPore( std::stringstream & inputs, voxelImageT<T>& vxlImage)
{
		std::cout<<"  converting to binary (0 and 1):"<<std::endl
			 <<"  selecting pore (->0) with values between:";
		unsigned int  thresholdMin=0,thresholdMax=0;
		inputs>>thresholdMin;
		inputs>>thresholdMax;

		std::cout<<" "<<int(thresholdMin)<<"  and "<<int(thresholdMax)<<"  inclusive."<<std::endl;
		vxlImage.threshold101(thresholdMin,thresholdMax);
		return true;
}

template<typename T> bool growPore( std::stringstream & inputs, voxelImageT<T>& vxlImage)
{
		std::cout<<"  growing voxels:"<<std::endl;
		int voxelValueTogrow; inputs>>voxelValueTogrow;
 		char growingAlgorithm; inputs>>growingAlgorithm;

		while (inputs.good())		  // loop while extraction from file is possible
		{
			if (growingAlgorithm!='f')
			{
				if(voxelValueTogrow==0)
					vxlImage.growPore();
				else if (voxelValueTogrow==0)
					vxlImage.shrinkPore();
				else
				{
					std::cerr<<"growing is only implemented for binary images: "<<
					"selected voxel value to grow is "<<voxelValueTogrow << ", which is not acceptable"<<std::endl;
					return false;//error occurred
				}
			}
			else
			{
				std::cerr<<"selected growing algorithm: "<<growingAlgorithm<<
				" the only implemented algorithm is f which stands for faceGrowing"<<std::endl;
				return false;//error occurred
			}

			inputs>>voxelValueTogrow;
			inputs>>growingAlgorithm;
		}
		std::cout<<" done"<<std::endl;
		return true;
}


template<typename T> bool resample( std::stringstream & inputs, voxelImageT<T>& vxlImage)
{
	double nResample=1;
		inputs>>nResample, std::cout<<"  resampling factor: "<<nResample<<std::endl;
		vxlImage.resample(nResample);
		return true;
}


template<typename T> bool resampleMax( std::stringstream & inputs, voxelImageT<T>& vxlImage)
{
	double nResample=1;
		inputs>>nResample, std::cout<<"  resampling factor: "<<nResample<<std::endl;
		vxlImage.resampleMax(nResample);
		return true;
}

template<typename T> bool redirect( std::stringstream & inputs, voxelImageT<T>& vxlImage)
{
		char direction;
		inputs>>direction, std::cout<<direction<<", swapping x and "<<direction<<" directions"<<std::endl;

		vxlImage.rotate(direction);
		return true;
}

template<typename T> bool replaceRange( std::stringstream & inputs, voxelImageT<T>& vxlImage)
{
	int  thresholdMin(0),thresholdMax(0); ///. Warning don't use T, uchar wont work
	inputs >> thresholdMin >> thresholdMax;

	int  value=(thresholdMin+thresholdMax)/2; ///. Warning don't use T, uchar wont work
	inputs >> value;

	std::cout<<" Replacing range  ["<<thresholdMin<<"  "<<thresholdMax<<"] with "<<value<<";   ";
	replaceRange(vxlImage,T(thresholdMin),T(thresholdMax),T(value));
	(std::cout<<".").flush();
	return true;
}

template<typename T> bool crop( std::stringstream & inputs, voxelImageT<T>& vxlImage)
{
	int cropBegin[3], cropEnd[3];

	std::cout<<"Crop:   ";
	for (int i=0; i<3;++i)   inputs>>cropBegin[i] >>cropEnd[i],  std::cout<<cropBegin[i]<<' '<<cropEnd[i]<<"	";
	std::cout<<' '<<std::endl;

	//cropEnd[0]+=1; cropEnd[1]+=1; cropEnd[2]+=1;
	vxlImage.crop(cropBegin,cropEnd);
	return true;
}


template<typename T> bool cropD( std::stringstream & inputs, voxelImageT<T>& vxlImage)
{
	int3 cropBegin{{0,0,0}}, cropEnd=vxlImage.sizeu3();
	int nLayers(0); int value(1);
	std::cout<<"cropD:   ";
	inputs>>cropBegin[0] >>cropBegin[1] >>cropBegin[2];  std::cout<<" "<<cropBegin[0] <<" "<<cropBegin[1] <<" "<<cropBegin[2]<<" --  ";  
	inputs>>cropEnd[0] >>cropEnd[1] >>cropEnd[2];		std::cout<<" "<<cropEnd[0] <<" "<<cropEnd[1] <<" "<<cropEnd[2]<<"  +  ";;  
	inputs >> nLayers >> value;
	std::cout<<nLayers<<" layers with "<<value<<std::endl;
	vxlImage.cropD(cropBegin,cropEnd,nLayers,value);
	return true;
}

template<typename T> bool write( std::stringstream & inputs, voxelImageT<T>& voximage)
{
	std::string outName("dump.tif");
	inputs >> outName;
	voximage.write(outName);
	(std::cout<<".").flush();
	return true;
}


template<typename T> bool medianFilter( std::stringstream & inputs, voxelImageT<T> & voximage)
{
	int nIterations(1); 
	inputs >> nIterations;
	(std::cout<<"  median Filter, nIterations: "<<nIterations).flush();
	voximage.growBox(2);
	for (int i=0; i<nIterations; ++i)
	{
		voximage=median(voximage);
	}
	voximage.shrinkBox(2);
	(std::cout<<".").flush();
	return true;
}


template<typename T> bool FaceMedian06( std::stringstream & inputs, voxelImageT<T> & voximage)
{
	int nIterations(1),  thereshold0(2), thereshold1(4); 
	inputs >> thereshold0>> thereshold1>> nIterations;
	(std::cout<<"  FaceMedian06: "<<thereshold0<<" "<<thereshold1<<" "<<nIterations<<"     ").flush();
	voximage.growBox(2);
	for (int i=0; i<nIterations; ++i)
	{
		voximage.FaceMedian06(thereshold0,thereshold1);
	}
	voximage.shrinkBox(2);
	(std::cout<<".").flush();
	return true;
}



template<typename T> bool PointMedian026( std::stringstream & inputs, voxelImageT<T> & voximage)
{
	int nIterations(1),  thereshold0(11), thereshold1(15); 
	inputs >> thereshold0>> thereshold1>> nIterations;
	(std::cout<<"  PointMedian026: "<<thereshold0<<" "<<thereshold1<<" "<<nIterations<<"     ").flush();
	voximage.growBox(2);
	for (int i=0; i<nIterations; ++i)
	{
		voximage.PointMedian026(thereshold0,thereshold1);
	}
	voximage.shrinkBox(2);
	(std::cout<<".").flush();
	return true;
}


template<typename T> bool circleOut( std::stringstream & inputs, voxelImageT<T> & voximage)
{

	char d='z';
	inputs >> d;
	int i = std::max<int>(d-'x',0);
	int X0(voximage.size3()[(i+1)%3]/2), Y0(voximage.size3()[(i+2)%3]/2);
	int R((X0+Y0)/2);

	inputs >> X0 >> Y0 >> R;
	(std::cout<<"  circleOut: dir="<<d<<",  X0="<<X0 <<"  Y0="<<Y0  <<"  R="<<R ).flush();

	circleOut(voximage,X0,Y0,R,d);

	(std::cout<<".").flush();
	return true;
}


template<typename T> bool maskWriteFraction( std::stringstream & inputs, voxelImageT<T> & voximage)
{
	int maskvv(2); 
	T minIelm(1), maxIelm=std::numeric_limits<T>::max();
	std::string maskname, outName("maskWriteFraction.txt");
	inputs >> maskname >> outName >> maskvv >> minIelm >> maxIelm;
	(std::cout<<"  maskWriteFraction:  mask:"<<maskname <<"  outName:"<<outName<<"  maskvv:"<<maskvv  <<"  minIelm:"<<minIelm<<"  maxIelm:"<<maxIelm ).flush();

	maskWriteFraction(voximage,maskname,outName,maskvv,minIelm,maxIelm);

	(std::cout<<".").flush();
	return true;
}


template<typename T> bool Offset( std::stringstream & inputs, voxelImageT<T> & voximage)
{
	vec3 offset; 
	inputs >> offset;
	(std::cout<<"  Offset:"<<offset<<" " ).flush();
	voximage.X0Ch()=offset;
	(std::cout<<".").flush();
	return true;
}


 



template<typename T> std::unordered_map<std::string,bool(*)( std::stringstream&, voxelImageT<T>&)> namedProcesses()
{

	typedef bool(*ProcessP)( std::stringstream&  inputs, voxelImageT<T>& vxlImage);
	return std::unordered_map<std::string,ProcessP>{
		{  "",& ignore },
		{  ";"		   , & ignore },
		{  "fillHoles"   , & fillHoles },
		{  "pore"		, & selectPore },
		{  "threshold"   , & selectPore },
		{  "threshold101"   , & selectPore },
		{  "resample"	, & resample },
		{  "Offset"   , & Offset },
		{  "direction"   , & redirect },
		{  "crop"		, & crop },
		{  "cropD"	   , & cropD },
		{  "resampleMax" , & resampleMax },
		{  "replaceRange", & replaceRange },
		{  "write"  ,& write },
		{  "medianFilter"  ,& medianFilter },
		{  "FaceMedian06"  ,& FaceMedian06 },
		{  "PointMedian026"  ,& PointMedian026 },
		{  "circleOut"  ,& circleOut },
		{  "maskWriteFraction"  ,& maskWriteFraction }	};
}


}




template<typename T>
void voxelImageT<T>::readFromHeader
(
	std::ifstream& headerFile,
	std::string header,
	int processKeys,
	std::string inputName
)
{
	int3 n{{0,0,0}};
	std::string BinaryData="XXX";
	bool X0read=false, dxread=false;
	double unit_=1.0;
	#ifdef TIFLIB
	if ((header.size()>4 && header.compare(header.size()-4,4,".tif") == 0))
	{
		(std::cout<<  " reading tif file "<<header<<" ").flush();
		readTif(*this, header);
		std::cout<<  "."<<std::endl;
		return;
	}
	else
	#endif
	if ((header.size()>4 && header.compare(header.size()-4,4,".mhd") == 0) || (header.size()>3 && header.compare(header.size()-3,3,".py") == 0))
	{
		std::cout<<" mhd:"<<header<<": "<<std::endl;
		while (true)
		{
			std::string tmpStr;
			std::streampos begLine = headerFile.tellg();
			headerFile>>tmpStr;

			

			if (headerFile.fail()) break;
			//~ ObjectType = Image
			//~ NDims = 3
			//~ Offset = 0 0 0
			//~ ElementSpacing = 8 8 8
			//~ DimSize = 200 225 153
			//~ ElementType = MET_UCHAR
			//~ ElementDataFile = Ketton100.raw
			std::stringstream keywordData;
			headerFile.get (*(keywordData.rdbuf()));
			std::string tmp;
			if (tmpStr == "ObjectType")
			{
				keywordData >> tmp; keywordData >> tmp;
				if (tmp != "Image") std::cout<<" Warning: ObjectType != Image :="<<tmp<<std::endl;
			}
			else if (tmpStr == "NDims")
			{
				keywordData >> tmp; keywordData >> tmp;
				if (tmp != "3") std::cout<<" Warning: NDims != 3 :="<<tmp<<std::endl;
			}
			else if (tmpStr == "ElementType")
			{
				keywordData >> tmp; keywordData >> tmp;
				if (tmp != "MET_UCHAR") std::cout<<" Warning: ElementType != MET_UCHAR :="<<tmp<<std::endl;
			}
			else if (tmpStr == "Offset")
			{
				keywordData >> tmp; keywordData>>	X0_[0]>>X0_[1]>>X0_[2] ;
				std::cout<<" X0: "<<  X0_[0]<<"  "<<X0_[1]<<"   "<<X0_[2]<<std::endl ;
			}
			else if (tmpStr == "ElementSpacing")
			{
				keywordData >> tmp; keywordData>>	dx_[0]>>dx_[1]>>dx_[2] ;
				std::cout<<" dX: "<< dx_[0]<<"  "<<dx_[1]<<"  "<<dx_[2]<<"   "<<std::endl; 
				if(dx_[0]>0.01)
				{
					std::cout<<"	 Warning: too large dx (="<<dx_[0]<<"), assuming unit is um. "<<std::endl;
					unit_ = 1.0e-6;
				}
			}
			else if (tmpStr == "DimSize")
			{
				keywordData >> tmp; keywordData>>	n[0]>>n[1]>>n[2];
				std::cout<<" Nxyz: "<<n[0]<<" "<<n[1]<<" "<<n[2]<<"   "<<std::endl; 
			}
			else if (tmpStr == "ElementDataFile")
			{
				keywordData >> tmp; if (inputName.empty()) keywordData >> inputName;

				size_t islash=header.find_last_of("\\/");
				if (islash<header.size() && inputName[0]!='/' &&  inputName[1]!=':') inputName=header.substr(0,islash+1)+inputName;
				std::cout<<" ElementDataFile = "<<inputName<<"	"<<std::endl;
			}
			else if (tmpStr == "BinaryData")
			{
				keywordData >> tmp; keywordData >> BinaryData;
				std::cout<<" BinaryData = "<<BinaryData<<"	"<<std::endl;
			}
			else if (tmpStr == "DefaultImageFormat")
			{
				std::string defSuffix;
				keywordData >> tmp; keywordData >> defSuffix;
				std::cout<<" OutputFormat = "<<defSuffix<<", suffix:"<<suffix(defSuffix)<<"	"<<std::endl; ///. sets suffix+format
			}
			else if (tmpStr == "Unit")
			{
				keywordData >> tmp; keywordData >> unit_;
				std::cout<<" Unit, OneMeter = "<<unit_<<std::endl;
			}
			else if (tmpStr!="BinaryDataByteOrderMSB" && tmpStr!="ElementByteOrderMSB" && tmpStr!="CompressedData" &&  tmpStr!="CompressedDataSize" &&  tmpStr!="TransformMatrix" &&
					 tmpStr!="CenterOfRotation" && tmpStr!="AnatomicalOrientation" && tmpStr!="AnatomicalOrientation")
			{
				headerFile.clear();
				headerFile.seekg(begLine);
				std::cout<<std::endl;
				break;
			}

		}


	}
	else
	{
		std::cout<<" (depricated) _header:"<<header<<","<<std::endl;

		char tmpc;
		for (int i=0; i<8;++i)   headerFile>>tmpc, std::cout<<tmpc;  //ignore the first 8 characters (ascii 3uc)

		if (header.size()>7 && header.compare(header.size()-7,7,"_header") == 0)  inputName=header.substr(0,header.size()-7);
		headerFile>>n[0]>>n[1]>>n[2];						// number of variables (dimension of
		std::cout<<"\n Nxyz: "<<n[0]<<" "<<n[1]<<" "<<n[2]<<"   "; std::cout.flush();
		headerFile>>	dx_[0]>>dx_[1]>>dx_[2] ;
		std::cout<<" dX: "<< dx_[0]<<"  "<<dx_[1]<<"  "<<dx_[2]<<"   "; std::cout.flush();
		headerFile>>	X0_[0]>>X0_[1]>>X0_[2] ;
		std::cout<<" X0: "<<  X0_[0]<<"  "<<X0_[1]<<"   "<<X0_[2] <<" um"<< std::endl;
		if (!headerFile)	 { std::cout<<"  Incomplete/bad header, aborting"<<std::endl; exit(-1);}
		//if (!headerFile)	 { std::cout<<"  Incomplete/bad header, continuing anyway"<<std::endl; }
		if(dx_[0]>0.01)
		{
			std::cout<<"Warning: too large dx (="<<dx_[0]<<"), assuming unit is um"<<std::endl;
			unit_ = 1.0e-6;
		}
	}



	dx_*=unit_;
	X0_*=unit_;
	this->reset(n,1);
	if( !inputName.empty() && inputName!="NO_READ" && processKeys!=2 )
	{
	  if (inputName.compare(inputName.size()-4,4,".tif") == 0)
	  {
			vec3 dx=dx_, X0=X0_;
			bool readingImage = this->readBin(inputName);
			assert(readingImage);
			if(X0read) X0_=X0;
			if(dxread) dx_=dx;
	  }
	  else if ((inputName.compare(inputName.size()-4,4,".raw") == 0 && BinaryData!="False") || BinaryData=="True")
	  {
			bool readingImage = this->readBin(inputName);
			assert(readingImage);
	  }
	  else if (inputName.size()>7 && inputName.compare(inputName.size()-7,7,".raw.gz") == 0)
	  {
			bool readingImage = this->readBin(inputName);
			assert(readingImage);
	  }
	  else
	  {
		std::ifstream infile(inputName.c_str());
		assert(infile);
		readAscii(infile);
	  }
	}

	typedef bool(*ProcessP)( std::stringstream&  inputs, voxelImageT<T>& vxlImage);


	std::unordered_map<std::string,ProcessP> name_Processes = MCTProcessing::namedProcesses<T>();

	if (processKeys)
	{


	while (true)
	{
		std::streampos begLine = headerFile.tellg();
		std::string tmpStr;
		headerFile>>tmpStr;
		//bool validKey=false;
		//cout<<tmpStr<<endl;///. keep me
		if (headerFile.fail())
		{std::cout<<" Finished reading "<<header<<":/  "<<headerFile.tellg()<<std::endl;  break; }
		else if (tmpStr[0]=='#' || tmpStr[0]=='\'' || tmpStr[0]=='%')
		{
			headerFile.ignore(10000,'\n');
			//validKey=true;
		}
		else
		{
			auto paer = name_Processes.find(tmpStr);
			if (paer!=name_Processes.end())	
			{
				(std::cout<<" "<<tmpStr<<": ").flush();
				std::stringstream keywordData;
				headerFile.get (*(keywordData.rdbuf()));
				(*(paer->second))(keywordData,*this);
				std::cout<<std::endl;
				//validKey=true;
			}
			else
			{	std::cout<<"  finished reading "<<header<<" before entry \""<<tmpStr<<"\":/ "<<std::endl;
				headerFile.clear();
				headerFile.seekg(begLine);
				break;
			}
		}
	}


	}

}










inline std::unique_ptr<voxelImageTBase> readImage
(
	std::string headerName,
	int processKeys = 1
)
{

	std::cout<<"Openning header file: "<<headerName<<std::endl;
	std::ifstream headerFile(headerName.c_str());
	if(!headerFile)  {std::cout<<"\n\n\nError: can not open header file, "<<headerName<<std::endl<<std::endl; }
	else
	{
	#ifdef TIFLIB
	if (headerName.size()>4 && headerName.compare(headerName.size()-4,4,".tif") == 0)
		{ headerFile.close(); return readTif(headerName); }
	#endif
	if (headerName.size()>4 && headerName.compare(headerName.size()-4,4,".mhd") == 0)
	{
		while (true)
		{
			std::string tmpStr;
			headerFile>>tmpStr;


			if (headerFile.fail()) break;

			std::stringstream keywordData;
			headerFile.get (*(keywordData.rdbuf()));
			std::string tmp;
			if (tmpStr == "ElementType")
			{
				keywordData >> tmp; keywordData >> tmp;
				headerFile.close();
				if (tmp=="MET_UCHAR")
				 { headerFile.close();	return std::unique_ptr<voxelImageTBase>(new voxelImageT<unsigned char>(headerName, processKeys)); }
				if (tmp=="MET_CHAR")
				 { headerFile.close();	return std::unique_ptr<voxelImageTBase>(new voxelImageT<char>(headerName, processKeys)); }
				if (tmp=="MET_USHORT")
				 { headerFile.close();	return std::unique_ptr<voxelImageTBase>(new voxelImageT<unsigned short>(headerName, processKeys)); }
				if (tmp=="MET_SHORT")
				 { headerFile.close();	return std::unique_ptr<voxelImageTBase>(new voxelImageT<short>(headerName, processKeys)); }
				if (tmp=="MET_UINT")
				 { headerFile.close();	return std::unique_ptr<voxelImageTBase>(new voxelImageT<unsigned int>(headerName, processKeys)); }
				if (tmp=="MET_INT")
				 { headerFile.close();	return std::unique_ptr<voxelImageTBase>(new voxelImageT<int>(headerName, processKeys)); }
				if (tmp=="MET_FLOAT")
				 { headerFile.close();	return std::unique_ptr<voxelImageTBase>(new voxelImageT<float>(headerName, processKeys)); }
				if (tmp=="MET_DOUBLE")
				 { headerFile.close();	return std::unique_ptr<voxelImageTBase>(new voxelImageT<double>(headerName, processKeys)); }
				  
			}

		}
	 }
	}
	
	headerFile.close();
	return std::unique_ptr<voxelImageTBase>(new voxelImageT<unsigned char>(headerName, processKeys));

}




template<typename T>
void readConvertFromHeader
(	voxelImageT<T>& vxlImg,
	std::string headerName,
	int processKeys = 1
)
{
	std::unique_ptr<voxelImageTBase> vxlImgTup = readImage(headerName,processKeys);
	voxelImageTBase* vxlImgT = vxlImgTup.get();
	
	bool red = false;
	{auto vxlImage = dynamic_cast<voxelImageT<char>* >(vxlImgT);			if(vxlImage) { vxlImg.resetFrom(*vxlImage); red=true; std::cout<<"read into "<<vxlImg.size3()<<" chars "; } }
	{auto vxlImage = dynamic_cast<voxelImageT<unsigned char>* >(vxlImgT);   if(vxlImage) { vxlImg.resetFrom(*vxlImage); red=true; std::cout<<"read into "<<vxlImg.size3()<<" ucars "; } }
	{auto vxlImage = dynamic_cast<voxelImageT<short>* >(vxlImgT);		   if(vxlImage) { vxlImg.resetFrom(*vxlImage); red=true; std::cout<<"read into "<<vxlImg.size3()<<" shrts "; } }
	{auto vxlImage = dynamic_cast<voxelImageT<unsigned short>* >(vxlImgT);  if(vxlImage) { vxlImg.resetFrom(*vxlImage); red=true; std::cout<<"read into "<<vxlImg.size3()<<" usrts "; } }
	{auto vxlImage = dynamic_cast<voxelImageT<int>* >(vxlImgT);			 if(vxlImage) { vxlImg.resetFrom(*vxlImage); red=true; std::cout<<"read into "<<vxlImg.size3()<<" intgs "; } }
	{auto vxlImage = dynamic_cast<voxelImageT<unsigned int>* >(vxlImgT); 	if(vxlImage) { vxlImg.resetFrom(*vxlImage); red=true; std::cout<<"read into "<<vxlImg.size3()<<" uints "; } }
	{auto vxlImage = dynamic_cast<voxelImageT<float>* >(vxlImgT);		   if(vxlImage) { vxlImg.resetFrom(*vxlImage); red=true; std::cout<<"read into "<<vxlImg.size3()<<" flots "; } }
	{auto vxlImage = dynamic_cast<voxelImageT<double>* >(vxlImgT);		  if(vxlImage) { vxlImg.resetFrom(*vxlImage); red=true; std::cout<<"read into "<<vxlImg.size3()<<" dobls "; } }

	if(!red) std::cout<<"\n\ncan not convert image\n\n"<<std::endl;
	if(!red) std::cerr<<"\n\ncan not convert image\n\n"<<std::endl;
}
