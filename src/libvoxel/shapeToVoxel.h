
/*-------------------------------------------------------------------------*\

This file is part of voxelImage library, a C++ template library  
developed by Ali Qaseminejad Raeini for handelling 3D raw images.


Please see our website for relavant literature making use of this code:
https://www.imperial.ac.uk/earth-science/research/research-groups/pore-scale-modelling/

For further information please contact us by email:
Ali Q Raeini: a.q.raeini@imperial.ac.uk

\*-------------------------------------------------------------------------*/

#include <math.h>
#include <fstream>
#include <assert.h>
#include <algorithm>
#include <vector>
#include <iostream>
#include "voxelImage.h"
	
#include <functional>   // std::mem_fn

namespace MCTProcessing
{
 using namespace std;
 class shape
 {
  public:
	const static int invalidv = 0x80000000; // 2147483648 largest negative
	int insidev, outsidev;
	shape() : insidev(0), outsidev(invalidv) {};
	virtual ~shape() {}; 
	virtual int value(dbl3)=0;
	virtual int isBefore(dbl3) { return false; };
	virtual int isAfter(dbl3) { return false; };
	template<typename T> void setIn(voxelImageT<T> & vImg)
	{	(cout <<__FUNCTION__<<": "<<insidev<<",  ").flush();
		dbl3 xmin = vImg.X0();
		dbl3 dx = vImg.dx();
		forAllkji_(vImg)
		{  int vv = this->value(dbl3(i+0.5,j+0.5,k+0.5)*dx+xmin);  if(vv!=shape::invalidv) vImg(i,j,k)  = vv; }
	}
	template<typename T> void addTo(voxelImageT<T> & vImg)
	{	(cout <<__FUNCTION__<<": "<<insidev<<",  ").flush();
		dbl3 xmin = vImg.X0();
		dbl3 dx = vImg.dx();
		forAllkji_(vImg)
		{  int vv = this->value(dbl3(i+0.5,j+0.5,k+0.5)*dx+xmin);  if(vv!=shape::invalidv) vImg(i,j,k) += vv; }
	}
	template<typename T> void setBefor(voxelImageT<T> & vImg)
	{	(cout <<__FUNCTION__<<": "<<insidev<<",  ").flush();
		dbl3 xmin = vImg.X0();
		dbl3 dx = vImg.dx();
		forAllkji_(vImg)
		{  if(this->isBefore(dbl3(i+0.5,j+0.5,k+0.5)*dx+xmin))  vImg(i,j,k)  = insidev; }
	}
	template<typename T> void setAfter(voxelImageT<T> & vImg)
	{	(cout <<__FUNCTION__<<": "<<insidev<<",  ").flush();
		dbl3 xmin = vImg.X0();
		dbl3 dx = vImg.dx();
		forAllkji_(vImg)
		{  if(this->isAfter (dbl3(i+0.5,j+0.5,k+0.5)*dx+xmin))  vImg(i,j,k)  = insidev; }
	}
	template<typename T> void addBefor(voxelImageT<T> & vImg)
	{	(cout <<__FUNCTION__<<": "<<insidev<<",  ").flush();
		dbl3 xmin = vImg.X0();
		dbl3 dx = vImg.dx();
		forAllkji_(vImg)
		{  if(this->isBefore(dbl3(i+0.5,j+0.5,k+0.5)*dx+xmin))  vImg(i,j,k)  += insidev; }
	}
	template<typename T> void addAfter(voxelImageT<T> & vImg)
	{	(cout <<__FUNCTION__<<": "<<insidev<<",  ").flush();
		dbl3 xmin = vImg.X0();
		dbl3 dx = vImg.dx();
		forAllkji_(vImg)
		{  if(this->isAfter (dbl3(i+0.5,j+0.5,k+0.5)*dx+xmin))  vImg(i,j,k)  += insidev; }
	}
 };


 class cylinder : public shape
 {
	dbl3 p1, p2;  double rr;
	double mag_p12;
 public:
	cylinder(stringstream & ins)
	{//http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
		ins>>p1 >> p2 >>rr >>insidev>>outsidev;
		mag_p12=mag(p2-p1);
		cout <<"cylinder: p1:"<<p1<<"    p2="<<p2<<"   r ="<<rr<<"   inside "<<insidev<<"   outside: "<<outsidev<<endl;
	};

	int value(dbl3 ij) { return ( mag((ij-p1)^(ij-p2)) <= rr*mag_p12 )  ?  insidev : outsidev; }
 };



 class layer : public shape
 {/// paraPlates, input mX dy mZ y1, mX and mZ are slopes, dy is separation, y1 is elevation of bottom plane
	double p1,p2; dbl3 nrm; ///< p1: location on bottom plane. p2: location on top plane. nrm: normal vector to planes, along which locations are recorded
  public:
	layer(stringstream & ins)
	: nrm(0,1,0)
	{	dbl3 Pl; //point through bottom plate, 
		ins>>Pl>>nrm>>insidev>>outsidev; p2=mag(nrm); nrm/=p2;  p1=Pl&nrm; p2+=p1;
		cout <<"  layer  p1="<<Pl<<",  height="<<p2-p1<<",  normal="<<nrm<<",  insidev="<<insidev<<endl;
	};
	int value(dbl3 ij)    {  return  (ij&nrm)>=p1 && (ij&nrm)<p2  ? insidev : outsidev;  }
	int isBefore(dbl3 ij) {  return  (ij&nrm) < p1; }// used for network cut
	int isAfter (dbl3 ij) {  return  (ij&nrm) >= p2;  }
 };


 // outdated
 class paraPlates : public shape
 {/// paraPlates, input mX dy mZ y1, mX and mZ are slopes, dy is separation, y1 is elevation of bottom plane
	dbl3 p1,p2;
  public:
	paraPlates(stringstream & ins)
	{
		ins>>p1;   p2 = p1;   p1.y=0.; ins>>p1.y;
		cout <<"\nparaPlates: slope1,mX="<<p1.x<<"    separation,dy="<<p2.y<<"   slope2,mZ="<<p1.z<<"   shift, y1:"<<p1.y<<endl;
	};
	int value(dbl3 ij)    {  return ( ij.y > p1.x*ij.x+p1.z*ij.z+p1.y 
		                           && ij.y < p2.x*ij.x+p2.z*ij.z+p2.y )  ? insidev : outsidev;  }
	int isBefore(dbl3 ij) {  return ( ij.y < p2.x*ij.x+p2.z*ij.z+p2.y ); }// used for network cut
	int isAfter (dbl3 ij) {  return ( ij.y > p1.x*ij.x+p1.z*ij.z+p1.y );  }
 };
 class kube : public shape
 {/// TODO  generalize to hexahedron, TODO optimize
	dbl3 p1,p2;
 public:
	kube(stringstream & ins)
	{
		ins>>p1>>p2>>insidev>>outsidev;
		cout <<"\nkube: "<<p1<<" + "<<p2<<" in: "<<insidev<<"  out: "<<outsidev<<endl;
		p2+=p1;
	};
	int value(dbl3 ij)  {  return ij.x >= p1.x && ij.y >= p1.y && ij.z >= p1.z &&   ij.x < p2.x && ij.y < p2.y && ij.z < p2.z  
		                             ? insidev : outsidev;  }
 };

 class plate : public shape
 {//! plane capped by sphere
	dbl3 p1,p2, po_;
  public:

	plate(stringstream & ins) 
	{cout<<"Error: fix me saqdakjoigfgfgfg "<<endl;
		ins>>p1>>po_; p2=p1;
		cout <<"\n plate: slope_x="<<p1.x<<"  r_cap="<<p1.y<<"  slope_z="<<p1.z<<"\n"<<endl;
		p1.y=0.;    po_.z = po_.y;  po_.y = (p1.x*po_.x+p1.z*po_.z+p1.y);
		cout <<"plate: a1="<<p1.x<<"   a2="<<p2.x<<"    r="<<p1.y<<"    b2="<<p2.y<<endl;
	};
	int value(dbl3 ij)  {  return  ( ij.y > p1.x*ij.x + p1.z*ij.z+p1.y && mag(ij-po_) <= p2.y )  ?  insidev : outsidev;  }
 };

 class sphere : public shape
 {
	dbl3 p1; double r2;
  public:
	sphere(stringstream & ins)
	{
		double rr;
		ins>>p1 >>rr >>insidev>>outsidev;		r2=rr*rr;
		cout <<"sphere: p1="<<p1<<"    r^2="<<std::sqrt(r2)<<endl;
	}
	int value(dbl3 ij)  {  return   ( magSqr(ij-p1)<r2 )  ?  insidev : outsidev;  }


	template<typename T> void setIn(voxelImageT<T> & vImg)
	{
		dbl3 xmin = vImg.X0();
		dbl3 dx   = vImg.dx();
		//int3 nnn  = vImg.size3()
		//dbl3 cn = p1/vImg.dx();
		//dbl3 rn(rr/vImg.dx()[0],rr/vImg.dx()[1],rr/vImg.dx()[2]);
		//dbl3 rnSqr=rn*rn;

		//if(outsidev==invalidv)
		//{
			//OMPFor()
			//for (int k = max(0,int(cn[2]-rn[2])); k <=  min(int(cn[2]+rn[2]+0.5),nnn[2]); ++k)
			//{	double rz=c-cn[2], ry=sqrt(rnSqr[2]-rz*rz);
				//const int nn1=std::min(lround(cn[1]-ry),nnn[2]);
				//for (int j = std::max(lround(cn[1]+ry),0); j<=nn1; ++j)
				//{	double rx=sqrt(rnSqr[1]-rz*rz-ry*ry);
					//std::fill( vImg(i,j,k), ForwardIt last, const T& value );
				//}
			//}
		//}


		forAllkji_(vImg)
		{  int vv = this->value(dbl3(i+0.5,j+0.5,k+0.5)*dx+xmin);  if(vv!=shape::invalidv) vImg(i,j,k)  = vv; }
	}
	template<typename T> void addTo(voxelImageT<T> & vImg)
	{
		dbl3 xmin = vImg.X0();
		dbl3 dx = vImg.dx();
		forAllkji_(vImg)
		{  int vv = this->value(dbl3(i+0.5,j+0.5,k+0.5)*dx+xmin);  if(vv!=shape::invalidv) vImg(i,j,k) += vv; }
	}
 };



#define _SHAPERATE(_operate_)  \
	if(ins.peek()=='?') { ins.str("\"s(phere), p(late, capped), f(flat-plates), c(ylinder) or k(ube)\" position..."); return true; } \
	std::string tmpc;  ins>>tmpc; \
	(cout <<__FUNCTION__<<" "<<tmpc<<",  ").flush(); \
	if      (tmpc[0]=='s')  sphere(ins)._operate_(vImg); \
	else if (tmpc[0]=='p')  plate(ins)._operate_(vImg); \
	else if (tmpc[0]=='f')  paraPlates(ins)._operate_(vImg); \
	else if (tmpc[0]=='k')  kube(ins)._operate_(vImg); \
	else if (tmpc[0]=='l')  layer(ins)._operate_(vImg); \
	else if (tmpc[0]=='c')  cylinder(ins)._operate_(vImg); \
	else cout <<"unsupported shape type: "<<tmpc<<"\n"<<endl 

 template<typename T>  bool Paint( stringstream& ins, voxelImageT<T> & vImg)  {
	_SHAPERATE(setIn);
	return true;
 }

 template<typename T>  bool PaintAdd( stringstream& ins, voxelImageT<T> & vImg) {
	_SHAPERATE(addTo);
	return true;
 }

 template<typename T>  bool PaintBefore( stringstream& ins, voxelImageT<T> & vImg) {
	_SHAPERATE(setBefor);
	return true;
 }

 template<typename T>  bool PaintAfter( stringstream& ins, voxelImageT<T> & vImg) {
	_SHAPERATE(setAfter);
	return true;
 }

 template<typename T>  bool PaintAddBefore( stringstream& ins, voxelImageT<T> & vImg) {
	_SHAPERATE(addBefor);
	return true;
 }

 template<typename T>  bool PaintAddAfter( stringstream& ins, voxelImageT<T> & vImg) {
	_SHAPERATE(addAfter);
	return true;
 }

}

