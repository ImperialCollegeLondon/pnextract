
/*-------------------------------------------------------------------------*\

This file is part of voxelImage library, a C++ template library  
developed by Ali Qaseminejad Raeini for handelling 3D raw images.


Please see our website for relavant literature making use of this code:
http://www3.imperial.ac.uk/earthscienceandengineering/research/perm/porescalemodelling

For further information please contact us by email:
Ali Q Raeini: a.q.raeini@imperial.ac.uk

\*-------------------------------------------------------------------------*/

#include <sys/stat.h>
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
	template<typename T> void setIn(voxelImageT<T> & vxlImage)
	{	(cout <<__FUNCTION__<<": "<<insidev<<",  ").flush();
		dbl3 xmin = vxlImage.X0();
		dbl3 dx = vxlImage.dx();
		forAllkji_(vxlImage)
		{  int vv = this->value(dbl3(i+0.5,j+0.5,k+0.5)*dx+xmin);  if(vv!=shape::invalidv) vxlImage(i,j,k)  = vv; }
	}
	template<typename T> void addTo(voxelImageT<T> & vxlImage)
	{	(cout <<__FUNCTION__<<": "<<insidev<<",  ").flush();
		dbl3 xmin = vxlImage.X0();
		dbl3 dx = vxlImage.dx();
		forAllkji_(vxlImage)
		{  int vv = this->value(dbl3(i+0.5,j+0.5,k+0.5)*dx+xmin);  if(vv!=shape::invalidv) vxlImage(i,j,k) += vv; }
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

 class paraPlates : public shape
 {/// paraPlates, input mX dy mZ y1, mX and mZ are slopes, dy is separation, y1 is elevation of bottom plane
	dbl3 p1,p2;
 public:
	paraPlates(stringstream & ins)
	{
		ins>>p1;   p2 = p1;   p1.y=0.0; ins>>p1.y;
		cout <<"\nparaPlates: slope1,mX="<<p1.x<<"    separation,dy="<<p2.y<<"   slope2,mZ="<<p1.z<<"   shift, y1:"<<p1.y<<endl;
	};
	int value(dbl3 ij)    {  return ( ij.y > p1.x*ij.x+p1.z*ij.z+p1.y 
		                           && ij.y < p2.x*ij.x+p2.z*ij.z+p2.y )  ? insidev : outsidev;  }
 };

 class plate : public shape
 {//! plane capped by sphere
	dbl3 p1,p2, po_;
  public:

	plate(stringstream & ins) 
	{cout<<"Error: fix me saqdakjoigfgfgfg "<<endl;
		ins>>p1>>po_; p2=p1;
		cout <<"\n plate: slope_x="<<p1.x<<"  r_cap="<<p1.y<<"  slope_z="<<p1.z<<"\n"<<endl;
		p1.y=0.0;    po_.z = po_.y;  po_.y = (p1.x*po_.x+p1.z*po_.z+p1.y);
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
 };


 template<typename T>  bool shapeToVoxel( stringstream & ins, voxelImageT<T> & vxlImage)
 {
	if(ins.peek()=='?') { ins.str("s for sphere, p for capped plate, f for parallel plates, c for cylinder"); return true;}

	std::string tmpc;  ins>>tmpc;
	(cout <<__FUNCTION__<<" "<<tmpc<<",  ").flush();

	if      (tmpc[0]=='s')  sphere(ins).setIn(vxlImage); 
	else if (tmpc[0]=='p')  plate(ins).setIn(vxlImage); ///. plate with sphere cap
	else if (tmpc[0]=='f')  paraPlates(ins).setIn(vxlImage);
	else if (tmpc[0]=='c')  cylinder(ins).setIn(vxlImage);
	else cout <<"unsupported shape type: "<<tmpc<<"\n"<<endl; 
	return true;
 }

 template<typename T>  bool shapeToVoxelAdd( stringstream & ins, voxelImageT<T> & vxlImage)
 {
	if(ins.peek()=='?') { ins.str("s for sphere, p for capped plate, f for parallel plates, c for cylinder"); return true;}

	std::string tmpc;  ins>>tmpc;
	(cout <<__FUNCTION__<<" "<<tmpc<<",  ").flush();

	if      (tmpc[0]=='s')  sphere(ins).addTo(vxlImage); 
	else if (tmpc[0]=='p')  plate(ins).addTo(vxlImage); ///. plate with sphere cap
	else if (tmpc[0]=='f')  paraPlates(ins).addTo(vxlImage);
	else if (tmpc[0]=='c')  cylinder(ins).addTo(vxlImage);
	else cout <<"unsupported shape type: "<<tmpc<<"\n"<<endl; 
	return true;
 }

}

