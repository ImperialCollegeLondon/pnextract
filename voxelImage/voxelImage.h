/*-------------------------------------------------------------------------*\


This file is part of voxelImage library, a small c++ template library  
developed by Ali Qaseminejad Raeini for handelling 3D raw images.


Please see our website for relavant literature making use of this code:
http://www3.imperial.ac.uk/earthscienceandengineering/research/perm/porescalemodelling

For further information please contact us by email:
Ali Q Raeini:    a.qaseminejad-raeini09@imperial.ac.uk

\*-------------------------------------------------------------------------*/


#ifndef voxelImageT_H
#define voxelImageT_H

#include <fstream>
#include <iostream>
#include <string.h>
#include <vector>
#include <valarray>
#include <cassert>
#include <sstream>
#include <memory>
#include <algorithm>
#include <numeric>
#include <functional>


#include "vec3.h"

/// class Config { 
inline const std::string& suffix(const std::string& defSuffix="")
{
	#ifdef TIFLIB
	 static std::string defSuffix_=".tif";
	#else
	  #ifdef ZLIB
  	   static std::string defSuffix_=".raw.gz";
	  #else
  	   static std::string defSuffix_=".raw";
	  #endif //ZLIB
	#endif //TIFLIB

	///. set via OutputFormat keyword
	if (defSuffix.size()) {
		if(defSuffix[0]!='.') defSuffix_="."+defSuffix;
		else                  defSuffix_=defSuffix;
		if( defSuffix_!=".tif" && defSuffix_!=".raw.gz" &&
		    defSuffix_!=".raw" && defSuffix_!=".dat" && defSuffix_!=".txt" )
			std::cout<<"\nError: wrong default image format: "<<defSuffix_<<"\n"<<std::endl;
	}

	return defSuffix_;
} /// };// Config:  set/get default suffix, uses static storage.






template <typename T> class voxelField
{
 protected:
 public:
	long long nij_;
	int3  nnn_;
	std::vector<T> data_;
	voxelField(): nij_(0), nnn_{{0,0,0}} {};
	voxelField(int3 n) {  reset(n);  };
	voxelField(int3 n, T value) {  reset(n, value);  };
	voxelField(int n1, int n2, int n3, T value) {  reset(int3{{n1,n2,n3}}, value);  };
	voxelField(const voxelField<T>& vm ): nij_(vm.nij_), nnn_(vm.nnn_), data_(vm.data_) {} ;
	virtual ~voxelField() {};

	void reset(int3 n);
	void reset(int3 n, T value);
	const T& operator()(const int i, const int j, const size_t k) const { return data_[k*nij_+j*nnn_[0]+i]; }
	      T& operator()(const int i, const int j, const size_t k)       { return data_[k*nij_+j*nnn_[0]+i]; }
	const T& operator()(const size_t iii) const { return data_[iii]; }
	      T& operator()(const size_t iii)       { return data_[iii]; }
	      T& operator()(long long ii, const size_t k)       { return data_[k*nij_+ii]; }
	const T& v_i(const int i, const T* vr) const { return *(vr+i); }
	const T& v_j(const int j, const T* vr) const { return *(vr+j*int(nnn_[0])); }
	const T& v_k(const int k, const T* vr) const { return *(vr+k* static_cast<long long>(nij_)); }
	const T* p_i(const int i, const T* vr) const { return (vr+i); }
	const T* p_j(const int j, const T* vr) const { return (vr+j*int(nnn_[0])); }
	const T* p_k(const int k, const T* vr) const { return (vr+k* static_cast<long long>(nij_)); }
	void reset(int n1, int n2, int n3, T value) {  reset(int3{{n1,n2,n3}}, value);  };
	void readMicroCT(std::string);
	void readMicroCTHeader(std::ifstream);
	bool readAscii(std::string);
	void readAscii(std::ifstream& in);
	bool readBin(std::string fileName);
	bool readBin(std::string fileName,int iStart,int iEndp1 , int jStart,int jEndp1 , int kStart,int kEndp1 );
	void writeNoHdr(std::string fileName) const;
	void writeBin(std::string fileName) const;
	void writeBin(std::string fileName,int iStart,int iEndp1 , int jStart,int jEndp1 , int kStart,int kEndp1 ) const;
	void writeAscii(std::string fileName) const;
	void writeAscii(std::string fileName,int iStart,int iEndp1 , int jStart,int jEndp1 , int kStart,int kEndp1) const;
	void writeRotatedXZ(std::ofstream& of) const;
	int3 size3() const;
	int3 sizeu3() const;
	void getSize(int& n1, int& n2, int& n3) const;


	void setLayer(int k, const T* Values);
	void setSlice(char dir, int ijk, T vv);
	void replaceyLayer(int j, int fromj);
	void replacexLayer(int i, int fromi);
	void setBlock(int n1, int n2, int n3, const voxelField<T>&Values);
};




class voxelImageTBase
{
public:
	virtual ~voxelImageTBase() {};
	virtual void write(std::string fileName) const = 0;
	virtual void printInfo() const = 0;
	virtual int3 sizeu3() const = 0;
	virtual int getInt(int i, int j, int k) const = 0;
	virtual double getDbl(int i, int j, int k) const = 0;
	virtual double vv_mp5(double i, double j, double k) const = 0;
	virtual const vec3& dx() const = 0;
	virtual const vec3& X0() const = 0;

};


template <typename T>
class voxelImageT: public voxelField<T>, public voxelImageTBase
{

	vec3	X0_, dx_;

 public:

	voxelImageT():X0_(0.0,0.0,0.0),dx_(1,1,1) {};


	voxelImageT(int n1, int n2, int n3, T value)
	: voxelField<T>( n1,  n2,  n3,  value),  X0_(0.0,0.0,0.0), dx_(1,1,1) {}


	voxelImageT(int3 n, vec3 dx, vec3 xmin, T value)
	: voxelField<T>( n[0],  n[1],  n[2],  value), X0_(xmin),dx_(dx) {}

	voxelImageT(const voxelImageT & vm)
	:  voxelField<T>(vm), X0_(vm.X0_), dx_(vm.dx_) {}


	voxelImageT(std::string headerName, int processKeys=1, std::string fileName="")
	: X0_(0.0,0.0,0.0),dx_(1,1,1)  {readFromHeader(headerName, processKeys,fileName);}
	void readFromHeader(std::string headerName, int processKeys=1, std::string fileName="")
	{	if (!headerName.empty())
		{	std::cout<<"Openning header file: "<<headerName<<std::endl;
			std::ifstream headerFile(headerName.c_str());
			if(!headerFile)  {std::cout<<"\n\n\nError: can not open header file, "<<headerName<<std::endl<<std::endl; }
			else
				readFromHeader(headerFile,headerName,processKeys,fileName);
			headerFile.close();
		}
	}

	void readFromHeader( std::ifstream& headerFile,	std::string headerName, int processKeys=1, std::string fileName="");



	bool readAscii(std::string fileName)
	{	///  overwrite as the parent -voxelField<T>- interprets
		///  numerical values as characters not integers

		std::cout<<  " reading "<<fileName<<std::endl;

		//if ( (fileName.compare(fileName.size()-4,4,".dat")==0) || (fileName.compare(fileName.size()-4,4,".txt") == 0) )
		//{
			std::ifstream in(fileName.c_str());
			assert(in);

			char tmpc[8];
			for ( int i=0; i<8;i++)   in>>tmpc[i];
			if (std::string(tmpc).compare(0,4,"ascii") == 0) //ignore first lines
			{
				int n[3];
				in>>n[2]>>n[0]>>n[1];//ignore first lines
				double  xmin[3],xmax[3];
				in>> xmin[0]>>xmax[0]>>xmin[1]>>xmax[1]>>xmin[2]>>xmax[2] ;
				std::cout<<"Warning: ignoring the header of file "<<fileName<<std::endl;
			}
			else
				in.seekg(0, in.beg);
			readAscii(in);
			in.close();
			return !in.fail();
		//}
		//else
			//this->readBin(fileName);

	}


	void  readAscii(std::ifstream& in)
	{	///  overwrite as the parent -voxelField<T>- interprets
		///  numerical values as characters not integers
		int tmp=0;
		for (auto& vv : voxelField<T>::data_)
		{
			in>>tmp;
			vv=tmp;
		}
	}



	void cropD( int3 cropBegin,  int3 cropEnd,int emptylayers=0, T emptylayersValue=1) ;
	void crop( int cropBegin[3],  int cropEnd[3],int emptylayers=0, T emptylayersValue=1) ;
	void crop(int iStart, int iEnd ,
				 int jStart, int jEnd ,
				 int kStart, int kEnd ,
				 int emptylayers=0,T emptylayersValue=1);

	void writeHeader(std::string fileName) const;
	void writeHeader(std::string fileName, int3 iStart, int3 iEnd) const;


	void erodeLayer(int i);
	void resample(double i);
	void resampleMax(double i);
	void rotate(char direction);
	void PointMedian026(int thereshold0,int thereshold1);
	void FaceMedian06(int thereshold0,int thereshold1);
	void mode(short nNeist);

	void AND(const voxelImageT& data2);
	void NOT(const voxelImageT& data2);
	void OR(const voxelImageT& data2);
	void XOR(const voxelImageT& data2);
	void maxEq(const voxelImageT& data2);
	void minEq(const voxelImageT& data2);


	void fillHoles(int maxHoleRadius);

	void shrinkPore();
	void growPore();


   void threshold101(T theresholdMin,T theresholdMax);


	void writeAConnectedPoreVoxel(std::string fileName) const;
	template<typename T2> void resetFrom(const voxelImageT<T2>&Values);
	void setFrom(const voxelImageT<T>&Values, int n1, int n2, int n3);
	void growBox(int nlyr);
	void shrinkBox(int nlyr)
		{int3 beg={{nlyr,nlyr,nlyr}},
		end={{((*this).size3()[0]-nlyr),((*this).size3()[1]-nlyr),((*this).size3()[2]-nlyr)}};
		cropD(beg,end);};

	void write(std::string fileName) const;
	double volFraction(T vv1,T vv2) const;
	void printInfo() const;

	const vec3& X0() const {return X0_;};
	vec3& X0Ch()    {return X0_;};
	const vec3& dx() const {return dx_;};
	vec3& dxCh()    {return dx_;};
	int3 sizeu3() const {return voxelField<T>::sizeu3();};
	int getInt(int i, int j, int k) const { return (*this)(i,j,k); };
	double getDbl(int i, int j, int k) const { return (*this)(i,j,k); };
	double vv_mp5(double i, double j, double k) const /// set i,j,k -=0.5 before passing them here
	{
		const int i0=i-0.0000001, j0=j-0.0000001, k0=k-0.0000001; ///. note neg fracs round to zero
		//if (i0<0 || i0+1>=voxelField<T>::nnn_[0] || j0<0 || j0+1>=voxelField<T>::nnn_[1] || k0<0 || k0+1>=voxelField<T>::nnn_[2])
			//cout<<"Error out of range: "<<"ijk0: "<<i0<<" "<<j0<<" "<<k0<<endl;
		double dd=i-i0;
		double Imd=1.0-dd;
		const T* vp=&(*this)(i0,j0,k0);
		const double v00=*vp*(Imd) + this->v_i(1,vp)*(dd);
		vp=this->p_j(1,vp);
		const double v10=*vp*(Imd) + this->v_i(1,vp)*(dd);
		vp=this->p_k(1,vp);
		const double v11=*vp*(Imd) + this->v_i(1,vp)*(dd);
		vp=&(*this)(i0,j0,k0+1);
		const double v01=*vp*(Imd) + this->v_i(1,vp)*(dd);
		dd=j-j0;
		Imd=1.0-dd;
		return ( v00*Imd+v10*dd )*(1.0-(k-k0))  + ( v01*Imd+v11*dd )*(k-k0);
	};

 };





typedef voxelImageT<unsigned char> voxelImage;

#include "voxelImageI.h"
#include "voxelImageReader.h"


#endif
