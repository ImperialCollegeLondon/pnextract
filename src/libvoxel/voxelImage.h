/*-------------------------------------------------------------------------*\

This file is part of voxelImage library, a C++ template library  
developed by Ali Qaseminejad Raeini for handelling 3D raw images.


Please see our website for relavant literature making use of this code:
https://www.imperial.ac.uk/earth-science/research/research-groups/pore-scale-modelling/

For further information please contact us by email:
Ali Q Raeini: a.q.raeini@imperial.ac.uk

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
#include <unordered_map>
#include <functional>


#include "typses.h"
#ifdef _STOR_PUB
#include "SiR.h" //_STOR 
#endif //_STOR_PUB 

//! suffix,  set/get default suffix, uses static storage.
inline const std::string& imgExt(const std::string& defSuffix="")
{
	#ifdef ZLIB
  	   static std::string defSuffix_=".raw.gz";
	#else
	  #ifdef TIFLIB
	   static std::string defSuffix_=".tif";
	  #else
  	   static std::string defSuffix_=".raw";
	  #endif //TIFLIB
	#endif //ZLIB

	///. set via OutputFormat keyword
	if (defSuffix.size()) {
		if(defSuffix[0]!='.') defSuffix_="."+defSuffix;
		else                  defSuffix_=defSuffix;
		if( defSuffix_!=".tif" && defSuffix_!=".raw.gz" && defSuffix_!=".am" &&
		    defSuffix_!=".raw" && defSuffix_!=".dat" && defSuffix_!=".txt")
			std::cout<<"\nError: wrong default image format: "<<defSuffix_<<"\n"<<std::endl;
	}

	return defSuffix_;
}


template<typename T>
using intOr = typename std::conditional<(sizeof(T) > sizeof(short)),T,int>::type;
#define Tint  intOr<T>

extern int maxNz;

template <typename T> class voxelField  //! 3D voxel data, on Cartesian uniform grids, with efficient access functions, and I/O into  .tif, .raw.gz .am, .raw, .dat file formats.
{
 protected:
 public:

	long long nij_;
	int3  nnn_;
	std::vector<T> data_;
	

	voxelField(): nij_(0), nnn_(0,0,0) {};
	voxelField(int3 n) {  reset(n);  };
	voxelField(int3 n, T value) {  reset(n, value);  };
	voxelField(int n1, int n2, int n3, T value) {  reset(int3(n1,n2,n3), value);  };
	voxelField(const voxelField<T>& vm ): nij_(vm.nij_), nnn_(vm.nnn_), data_(vm.data_) {} ;
	//voxelField(voxelField<T>&& vm ): nij_(vm.nij_), nnn_(vm.nnn_), data_(std::move(vm.data_))  {} ;

	void reset(int3 n);
	void reset(int3 n, T value);
	void reset(int n1, int n2, int n3, T value) {  reset(int3(n1,n2,n3), value);  };
	void readMicroCT(std::string);
	bool readAscii(std::string);
	void readAscii(std::ifstream& in);
	bool readBin(std::string fileName, int nSkipBytes=0);
	bool readBin(std::string fileName,int iBgn,int iEndp1 , int jBgn,int jEndp1 , int kBgn,int kEndp1, int nSkipBytes=0);
	void writeNoHdr(std::string fileName) const;
	void writeBin(std::string fileName) const;
	void writeBin(std::string fileName,int iBgn,int iEndp1 , int jBgn,int jEndp1 , int kBgn,int kEndp1 ) const;
	void writeAscii(std::string fileName) const;
	void writeAscii(std::string fileName,int iBgn,int iEndp1 , int jBgn,int jEndp1 , int kBgn,int kEndp1) const;
	void writeRotatedXZ(std::ofstream& of) const;

	void writeHeader(std::string fileName, int3 iBgn, int3 iEnd, dbl3 dx, dbl3 X0) const;

	void setLayer(int k, const T* Values);
	void setSlice(char dir, int ijk, T vv);
	void replacezLayer(int j, int fromj);
	void replaceyLayer(int j, int fromj);
	void replacexLayer(int i, int fromi);
	void setBlock(int n1, int n2, int n3, const voxelField<T>&Values);
	void setFrom(const voxelField<T>&Values, int n1, int n2, int n3);



	size_t index(const int i, const int j, const size_t k) const { return k*nij_+j*nnn_.x+i; }
	const T& operator()(const int i, const int j, const size_t k) const { return data_[k*nij_+j*nnn_.x+i]; }
	      T& operator()(const int i, const int j, const size_t k)       { return data_[k*nij_+j*nnn_.x+i]; }
	const T& operator()(const size_t iii) const     { return data_[iii]; }
	      T& operator()(const size_t iii)           { return data_[iii]; }
	      T& operator()(long long ii, const size_t k){return data_[k*nij_+ii]; }
	const T& v_i(const int i, const T* vr) const    { return *(vr+i); }
	const T& v_j(const int j, const T* vr) const    { return *(vr+j*nnn_.x); }
	const T& v_k(const int k, const T* vr) const    { return *(vr+k* nij_); }
	      T& v_i(const int i, T* vr)                { return *(vr+i); }
	      T& v_j(const int j, T* vr)                { return *(vr+j*nnn_.x); }
	      T& v_k(const int k, T* vr)                { return *(vr+k*nij_); }
	const T* p_i(const int i, const T* vr)    const { return (vr+i); }
	const T* p_j(const int j, const T* vr)    const { return (vr+j*nnn_.x); }
	const T* p_k(const int k, const T* vr)    const { return (vr+k*nij_); }
	size_t I_i(const int i, const size_t iii) const { return (iii+i); }
	size_t I_j(const int j, const size_t iii) const { return (iii+j*nnn_.x); }
	size_t I_k(const int k, const size_t iii) const { return (iii+k*nij_); }

	T* begin()              { return &*data_.begin(); };
	T* end()                { return &*data_.end();  };
	const T& back()   const { return data_.back();  };
	const T* cbegin() const { return &*data_.cbegin();};
	const T* cend()   const { return &*data_.cend();};
	//T* operator()()   const { return begin(); };

	const int3& size3() const   {  return nnn_;  };
	int nx() const   {  return nnn_.x;  };
	int ny() const   {  return nnn_.y;  };
	int nz() const   {  return nnn_.z;  };
	long long  nxy() const   {  return nij_;  };
	void getSize(int& n1, int& n2, int& n3) const;
	virtual ~voxelField() {};

};




class voxelImageTBase //! Base class handling different image files with different data types (float, char, int...)
{
 public:
	virtual ~voxelImageTBase() {};
	virtual void write(std::string fileName) const = 0;
	virtual void printInfo() const {};
	virtual std::unique_ptr<voxelImageTBase> copy() const = 0;
	virtual const int3& size3() const = 0;
	virtual const dbl3& dx() const = 0;
	virtual const dbl3& X0() const = 0;
	virtual void readFromHeader(std::istream& headerFile, const std::string& headerName, int processKeys) = 0;
};


template <typename T>
class voxelImageT: public voxelImageTBase, public voxelField<T>   //!  3D image data of different types (T = float, char, int...)
{
	dbl3	X0_;   //!< origin
	dbl3	dx_;   //!< voxel size

 public:
	//string bname_; //!< base name, optional

	voxelImageT():X0_(0.,0.,0.),dx_(1,1,1) {};

	voxelImageT(int n1, int n2, int n3, T value) //do not remove, the following constructor will be misused in old codes!
	: voxelField<T>( n1,  n2,  n3,  value),  X0_(0.,0.,0.), dx_(1.,1.,1.) {}

	voxelImageT(int3 n, dbl3 dx=dbl3(1.,1.,1.), dbl3 xmin=dbl3(0.,0.,0.), T value=0)
	: voxelField<T>( n.x,  n.y,  n.z,  value), X0_(xmin),dx_(dx) {}

	voxelImageT(const voxelImageT & vm)
	:  voxelField<T>(vm), X0_(vm.X0_), dx_(vm.dx_) {}



	voxelImageT(const std::string& headerName, int processKeys=1)
	:	X0_(0.,0.,0.),dx_(1,1,1)  { readFromHeader(headerName, processKeys); }

	void readFromHeader(const std::string& headerName, int processKeys=1)
	{	if (!headerName.empty() && headerName!="NO_READ")
		{	std::cout<<"  Openning header: "<<headerName<<std::endl;
			std::ifstream hdr(headerName);
			if(hdr) this->readFromHeader(hdr,headerName,processKeys);
			else  alert("cannot open header file, "+headerName); 
		}
	}

	void readFromHeader(std::istream& headerFile, const std::string& headerName, int processKeys=1);



	bool readAscii(std::string fileName);
	void readRLE(std::string fileName);

	std::unique_ptr<voxelImageTBase> copy() const { return std::make_unique<voxelImageT<T>>(*this); };


	void cropD( int3 frm,  int3 to,int emptylyrs=0, T eLyrsValue=1, bool verbose=false) ;
	void cropOld(int cropBgn[3],  int cropEnd[3],int emptylyrs=0, T eLyrsValue=1) ;
	void cropOld(int iBgn, int iEnd, int jBgn, int jEnd, int kBgn, int kEnd, int emptylyrs=0,T eLyrsValue=1);

	void writeHeader(std::string fileName) const;
	void writeHeader(std::string fileName, int3 iBgn, int3 iEnd) const { voxelField<T>::writeHeader(fileName,iBgn,iEnd,dx_,X0_); };
	void writeNoHdr(std::string fileName) const;
	void write(std::string fileName) const;

	void erodeLayer(int i);
	void rotate(char direction);
	void PointMedian032(int nAdj0, int nAdj1, T lbl0, T lbl1);
	size_t FaceMedian06(int nAdj0,int nAdj1); // obsolete
	void mode(short nNeist, bool verbose=false);
	void zeroGrad(int nlyr);

	void AND(const voxelImageT& data2);
	void NOT(const voxelImageT& data2);
	void OR(const voxelImageT& data2);
	void XOR(const voxelImageT& data2);
	void maxEq(const voxelImageT& data2);
	void minEq(const voxelImageT& data2);


	void fillHoles(int maxHoleRadius);

	void shrinkPore();
	void growPore();
	void growLabel(T vl);


   void threshold101(T theresholdMin,T theresholdMax);


	void writeAConnectedPoreVoxel(std::string fileName) const;
	template<typename T2> void resetFrom(const voxelImageT<T2>&Values);
	void setFrom(const voxelImageT<T>&Values, int n1, int n2, int n3);
	void growBox(int nlyr);
	void shrinkBox(int nlyr) {  int3 bgn(nlyr,nlyr,nlyr);  cropD(bgn,this->size3()-bgn);  };

	double volFraction(T vv1,T vv2) const;
	void printInfo() const;

	const dbl3& X0  ()  const { return X0_; };
	dbl3&       X0Ch()        { return X0_; };
	const dbl3& dx  ()  const { return dx_; };
	dbl3&       dxCh()        { return dx_; };
	const int3& size3() const { return voxelField<T>::size3(); };

	double vv_mp5(double i, double j, double k) const { /// linear interpolation
		///set i,j,k -=0.5 (=_mp5) before passing them here, assuming vxl centres are at +0.5
		const int i0=std::min(int(i),this->nx()-2), j0=std::min(int(j),this->ny()-2), k0=std::min(int(k),this->nz()-2); ///. note neg fracs round to zero
		double dd=i-i0,    ld=1.-dd;
		const T* 
		vp=&(*this)(i0,j0,k0);  const double v00= *vp*ld + dd*this->v_i(1,vp);
		vp=this->p_j( 1,vp);    const double v10= *vp*ld + dd*this->v_i(1,vp);
		vp=this->p_k( 1,vp);    const double v11= *vp*ld + dd*this->v_i(1,vp);
		vp=this->p_j(-1,vp);    const double v01= *vp*ld + dd*this->v_i(1,vp);
		dd=j-j0;    ld=1.-dd;       k-=k0;
		return (v00*ld + dd*v10) *(1.-k)+k* (v01*ld + dd*v11);
	};
};


std::string VxlKeysHelp(std::string keyname="", std::string subkey="");


template<class InpT, typename T> int vxlProcess(const InpT& inks, voxelImageT<T>& img, std::string nam="");// InpT= string or InputFile

template<class InpT, typename First=uint8_t, typename... Rest> int vxlProcess(const InpT& inks, voxelImageTBase* ptr, std::string nam="");


std::unique_ptr<voxelImageTBase> readImage(std::string hdrNam /*headername or image type*/, int procesKeys = 1 );

template<class T, typename First=uint8_t, typename... Rest>
int resetFromImageT(voxelImageT<T>& vImg, voxelImageTBase* imgPtr) { //! cast to specified vImg type 
	if(auto img = dynamic_cast<voxelImageT<First>*>(imgPtr)) { vImg.resetFrom(*img); return 0; }
	else if(sizeof...(Rest)) return resetFromImageT<T,Rest...>(vImg, imgPtr);
	alert("Unknown image type in resetFromImageT");
	return -1;
}

template<typename T>
void readConvertFromHeader( voxelImageT<T>& vImg, std::string hdrNam, int procesKeys=1)  { //! read image and if needed convert its type
	std::unique_ptr<voxelImageTBase> vImgUptr = readImage(hdrNam,procesKeys);
	voxelImageTBase* imgPtr = vImgUptr.get();
	if  (auto img = dynamic_cast<voxelImageT<T>*>(imgPtr)) { vImg = std::move(*img); } //TODO ensure use swap
	else  ensure(  resetFromImageT(vImg, imgPtr)==0, "can not convert image", -1);
}

template<class T> voxelImageT<T>* vxlCast(voxelImageTBase* imgPtr) { return dynamic_cast<voxelImageT<T>*>(imgPtr); }

template<typename T>
voxelImageT<T> copyOrReadImgT(std::string hdrNam) {
	#ifdef _STOR_PUB 
	if (void* ptr = dbget(_STOR,hdrNam,0)) {
		auto imgPtr = vxlCast<T>(static_cast<voxelImageTBase*>(ptr)); ensure(imgPtr, "wrong image type", -1);
		return *imgPtr;
	}
	else
	#endif
		return  voxelImageT<T>(hdrNam);


//template<typename ImgT>
//std::unique_ptr<ImgT, std::function<void(ImgT*)> getOrReadImgT(std::string hdrNam) {
	//#ifdef _STOR_PUB 
	//if (void* ptr = dbget(_STOR,hdrNam,0)) {
		//auto imgPtr = vxlCast<T>(static_cast<voxelImageTBase*>(ptr)); ensure(imgPtr, "wrong image type", -1);
		//return std::unique_ptr<vxlT, std::function<void(vxlT*)>(imgPtr,[](vxlT* ptr){});
	//}
	//else
	//#endif
		//return  std::unique_ptr<vxlT, std::function<void(vxlT*)>(new ImgT(hdrNam),[](vxlT* ptr){ delete ptr; });
//}
}

#ifdef _STOR_PUB 
	#define _dbgetOrReadImgT(_img,_hdrNam)  \
			voxelImageT<T>* _imgPtr;   std::unique_ptr<voxelImageTBase> _imgRead;  \
			if (void* ptr = dbget(_STOR,_hdrNam,0)) { _imgPtr = vxlCast<T>(static_cast<voxelImageTBase*>(ptr)); }  \
			else       { _imgRead = readImage(_hdrNam,1); _imgPtr=vxlCast<T>(_imgRead.get()); }  \
			ensure(_imgPtr, "wrong image type", -1);  \
			const voxelImageT<T>& _img = *_imgPtr;
#else
	#define _dbgetOrReadImgT(_img,_hdrNam) voxelImageT<T> _img(_hdrNam)
#endif

typedef voxelImageT<unsigned char> voxelImage;   //! default image format

#include "voxelImageI.h"

#endif
