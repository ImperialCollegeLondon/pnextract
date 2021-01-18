#ifndef TYPSES_H
#define TYPSES_H


// Convinience vector classes used by network extraction, flow simulation 
// and other codes developed by Ali Qaseminejad Raeini.
// Main template classes defined here: var3, var2, piece and lazyvec 





#include <iomanip>
#include <sstream>
#include <fstream>
#include <array>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <map>
#include <iostream>
#include <assert.h>
#include <regex>

#ifndef verySmall
	#define verySmall  1.0e-31
	#define   maxT(Typ)  (std::numeric_limits<Typ>::max())
	#define   minT(Typ)  (std::numeric_limits<Typ>::min())
	#define  iminT(Typ) (int(std::numeric_limits<Typ>::min()))
	#define  imaxT(Typ) (int(std::numeric_limits<Typ>::max()))
	#define  fmaxT(Typ) (float(std::numeric_limits<Typ>::max()))
	#define  dmaxT(Typ) (double(std::numeric_limits<Typ>::max()))
	#define  epsT(Typ)  std::numeric_limits<Typ>::epsilon()
#endif

#if (__cplusplus < 201403L )
	template< bool B, class T = void >
	using enable_if_t = typename std::enable_if<B,T>::type;
#else
	using std::enable_if_t;
#endif




#include "globals.h"


//! 3D vector class template
template<typename T>
class var3
{
 public:
	T	x, y, z;
	var3() = default;// Warning: not zero initialized in opt mode
	template< typename std::enable_if<std::is_arithmetic<T>::value,int>::type = 0> 
	var3(T r)  { x = r;  y = r;  z = r; } // use this to zero-initialize
	var3(T r, T s, T t)  { x = r;  y = s;  z = t; }
	var3(const T* d)  { x = d[0];  y = d[1];  z = d[2]; }
	template<typename T2>
	var3(var3<T2> n)  { x = n.x;  y = n.y;  z = n.z; }
	//var3(std::array<int,3> n)  { x = n[0];  y = n[1];  z = n[2]; }
#ifdef VMMLIB__VECTOR__HPP
	var3(const Vctr<3,T>& v3) 	{ x = v3[ 0 ];	y = v3[ 1 ];	z = v3[ 2 ]; 	}
#endif

	T&       operator [](long k) {	return ((&x)[k]); }
	const T& operator [](long k) const  { return ((&x)[k]); }
	T _0() const  { return x; }
	T _1() const  { return y; }
	T _2() const  { return z; }


    //operator  int() const { return x; }
    //operator  double() const { return x; }

	var3&  operator += (const var3& v)  { x += v.x;  y += v.y;  z += v.z;  return (*this); }
	var3&  operator -= (const var3& v)  { x -= v.x;  y -= v.y;  z -= v.z;  return (*this); }
	var3&  operator += (const T& t)  { x += t;  y += t;  z += t;  return (*this); } // clumsy
	var3&  operator -= (const T& t)  { x -= t;  y -= t;  z -= t;  return (*this); } // clumsy
	var3&  operator *= (const double t) { x *= t;  y *= t;  z *= t;  return (*this); }
	var3&  operator /= (const double t) { double f = 1.0 / t;  x *= f;  y *= f;  z *= f;  return (*this); }
	var3&  operator ^= (const var3& v)  { double r, s;  r=y*v.z-z*v.y;  s=z*v.x-x*v.z;  z=x*v.y-y*v.x;  x=r; y=s; 	return (*this); }
	var3&  operator *= (const var3& v)  { x *= v.x;  y *= v.y;  z *= v.z;  return (*this); }
	var3   operator -  (void) const          { return (var3(-x, -y, -z)); }
	var3   operator +  (const var3& v) const { return (var3(x+v.x, y+v.y, z+v.z)); }
	var3   operator -  (const var3& v) const { return (var3(x-v.x, y-v.y, z-v.z)); }
	var3   operator +  (const T& t) const { return (var3(x+t, y+t, z+t)); } // clumsy
	var3   operator -  (const T& t) const { return (var3(x-t, y-t, z-t)); } // clumsy
	var3   operator *  (const double t)const { return (var3(x*t, y*t, z*t)); }
	var3   operator /  (const double t)const { double f = 1.0 / t;  return (var3(x*f, y*f, z*f)); }
	double operator &  (const var3& v) const { return (x*v.x+y*v.y+z*v.z); }
	var3   operator ^  (const var3& v) const { return (var3(y*v.z-z*v.y,  z*v.x-x*v.z,  x*v.y-y*v.x)); }
	var3   operator *  (const var3& v) const { return (var3(x*v.x, y*v.y, z*v.z)); }
	var3   operator /  (const var3& v) const { return (var3(x/v.x, y/v.y, z/v.z)); }
	bool   operator == (const var3& v) const { return ((x-v.x)*(x-v.x) < verySmall) && ((y-v.y)*(y-v.y) < verySmall) && ((z-v.z)*(z-v.z) < verySmall); }
	bool   operator != (const var3& v) const { return ((x-v.x)*(x-v.x) >= verySmall) || ((y-v.y)*(y-v.y) >= verySmall) || ((z-v.z)*(z-v.z) >= verySmall); }
};
typedef  var3<int> int3;
typedef  var3<var3<int> > int3x3;

typedef   var3<float> float3;
typedef   var3<double> dbl3;

template<class T>  var3<T> rotateAroundLine(var3<T> y, double gamma,  var3<T> n, var3<T> x)
{///. rotate y around line passing through x, in the direction of n, http://inside.mines.edu/~gmurray/ArbitraryAxisRotation
	double s = sinf(gamma),   c = cosf(gamma);
	double k = 1.0 - c;
	return var3<T>(
	 	( x.x*(n.y*n.y+n.z*n.z) - n.x*( x.y*n.y+x.z*n.z-n.x*y.x- n.y*y.y-n.z*y.z ) )*k + y.x*c + (-x.z*n.y+x.y*n.z-n.z*y.y+n.y*y.z )*s,
		( x.y*(n.x*n.x+n.z*n.z) - n.y*( x.x*n.x+x.z*n.z-n.x*y.x- n.y*y.y-n.z*y.z ) )*k + y.y*c + ( x.z*n.x-x.x*n.z+n.z*y.x-n.x*y.z )*s,
		( x.z*(n.x*n.x+n.y*n.y) - n.z*( x.x*n.x+x.y*n.y-n.x*y.x- n.y*y.y-n.z*y.z ) )*k + y.z*c + (-x.y*n.x+x.x*n.y-n.y*y.x+n.x*y.y )*s );
}
template<class T>  var3<T> rotateAroundVec(const var3<T> y, double gamma, var3<T> n)
{///. rotate y around n (line passing through centre, in the direction of n) http://inside.mines.edu/~gmurray/ArbitraryAxisRotation
	double s = sinf(gamma),   c = cosf(gamma);
	double k = 1.0 - c;
	return var3<T>(
		(  - n.x*( -n.x*y.x- n.y*y.y-n.z*y.z ) )*k + y.x*c + (n.y*y.z-n.z*y.y)*s,
		(  - n.y*( -n.x*y.x- n.y*y.y-n.z*y.z ) )*k + y.y*c + (n.z*y.x-n.x*y.z)*s,
		(  - n.z*( -n.x*y.x- n.y*y.y-n.z*y.z ) )*k + y.z*c + (n.x*y.y-n.y*y.x)*s );
}

//! 2D vector class template
template<typename T>
class var2
{
 public:
	T a;//union {T	first;  T x;};
	T b;//union {T	second; T y;};
	var2() = default;//not zero initialized 
	var2(T r, T s)  { a = r;  b = s;}
	var2(int r)  { a = r;  b = r; }  // use this to zero-initialize //ERROR REMOVE ME
	var2(const T* d)  { a = d[0];  b = d[1]; }
	var2(std::pair<T,T> d)  { a = d.first;  b = d.second; }

	T&       operator [](long k) {	return ((&a)[k]); }
	const T& operator [](long k) const  { return ((&a)[k]); }

	explicit operator int() const { return a; }
	explicit operator double() const { return a; }
	explicit operator std::pair<T,T>&() const { return *this; }

	T x() const  { return a; }
	T y() const  { return b; }

	T first() const  { return a; }
	T second() const  { return b; }

	var2&  operator += (const var2& v)  { a += v.b;  b += v.b;  return (*this); }
	var2&  operator -= (const var2& v)  { a -= v.b;  b -= v.b;   return (*this); }
	var2&  operator *= (const double t) { a *= t;  b *= t;   return (*this); }
	var2&  operator /= (const double t) { double f = 1.0 / t;  a *= f;  b *= f;  return (*this); }
	var2&  operator *= (const var2& v)  { a *= v.b;  b *= v.b;   return (*this); }
	var2   operator -  (void) const          { return (var2(-a, -b)); }
	var2   operator +  (const var2& v) const { return (var2(a+v.b, b+v.b)); }
	var2   operator -  (const var2& v) const { return (var2(a-v.b, b-v.b)); }
	var2   operator *  (const double t)const { return (var2(a*t, b*t)); }
	var2   operator /  (const double t)const { double f = 1.0 / t;  return (var2(a*f, b*f)); }
	double operator &  (const var2& v) const { return (a*v.b+b*v.b); }
	var2   operator *  (const var2& v) const { return (var2(a*v.b, b*v.b)); }
	var2   operator /  (const var2& v) const { return (var2(a/v.b, b/v.b)); }
	bool   operator == (const var2& v) const { return ((a-v.b)*(a-v.b) < verySmall) && ((b-v.b)*(b-v.b) < verySmall) ; }
	bool   operator != (const var2& v) const { return ((a-v.b)*(a-v.b) >= verySmall) || ((b-v.b)*(b-v.b) >= verySmall) ; }
};

inline  dbl3  operator *(const int3& v, const dbl3& d) { return dbl3(v.x*d.x, v.y*d.y, v.z*d.z); } // boost int to T
template<typename T>  var3<T>  operator *(double t, const var3<T>& v) { return var3<T>(t*v.x, t*v.y, t*v.z); }
template<typename T>  T  mag(const var3<T>& v)             { return std::sqrt(v.x*v.x+v.y*v.y+v.z*v.z); }
template<typename T>  T  sumvars(const var3<T>& v)             { return (v.x+v.y+v.z); }
template<typename T>  double      magSqr(const var3<T>& v) { return (v.x*v.x+v.y*v.y+v.z*v.z); }
template<typename T>  var3<T>  norm(const var3<T>& v)      { return  v/(mag(v)+1.0e-300); }


//. component-wise max and min
template<typename T> var3<T> max(const var3<T>& a, const var3<T>& b)  { return var3<T>(std::max(a.x,b.x),std::max(a.y,b.y),std::max(a.z,b.z)); }
template<typename T> var3<T> min(const var3<T>& a, const var3<T>& b)  { return var3<T>(std::min(a.x,b.x),std::min(a.y,b.y),std::min(a.z,b.z)); }

template<class T>  T  toScalar(const T& v)       { return v; } // NOW IMPLICITly CONVERTED
template<class T>  T  toScalar(const var3<T>& v) { return mag(v); }





//! a piece of a contiguous array, used for efficient and prettier pass of array contents than using iterators
template<class T>
class piece  
{
  public:

	piece(): d(0), dn(0) {};
	piece(T* dd, int n): d(dd), dn(d+n)     {};
	piece(T* dd, T* de): d(dd), dn(de)      {};
	piece(const piece& p): d(p.d), dn(p.dn) {};//! Note this is different from = operator note data hold by piece are not const unless piece is const itself
	piece(std::vector<T>& vs): d(&vs[0]), dn(d+vs.size()) {};
	void reset(T* dd, int n)     { d=dd; dn=d+n; };
	void reset(const piece& vs)     { d=&vs[0]; dn=d+vs.size(); };//! note data hold by piece are not const unless piece is const itself
	void reset(std::vector<T>& vs)  { d=&vs[0]; dn=d+vs.size(); };

	T* begin() const {return d;};
	T* end()   const {return dn;};
	const T& back()   const {return *(dn-1);};
	const T* cbegin() const {return d;};
	const T* cend()   const {return dn;};
	T* operator ()()   const {return d;};
	//const T& operator [](int i) const {return d[i];};
	//const T* operator ()() const {return d;};
	T& operator [](int i) const {return d[i];};
	size_t size() const        {return dn-d;};
	size_t capacity()          {return dn-d;}
	T* data() {return d;}
	const T* data() const {return d;}

	piece& fill(const T& v)  { std::fill(d, dn, v);  return (*this); }
	piece& operator =(const piece& v)  { ensure(size()==v.size(), "in operator =, piece sizes should be the same"); std::copy(v.d, v.dn, d);  return (*this); }
	piece& operator +=(const piece& v)  { for(auto& a:*this){ a += v[&a-d];};  return (*this); }
	piece& operator -=(const piece& v)  { for(auto& a:*this){ a -= v[&a-d];};  return (*this); }
	piece& operator *=(const piece& v)  { for(auto& a:*this){ a *= v[&a-d];};  return (*this); }
	piece& operator /=(const piece& v)  { for(auto& a:*this){ a /= v[&a-d];};  return (*this); }
	piece& operator +=(T v)             { for(auto& a:*this){ a += v;};  return (*this); }
	piece& operator -=(T v)             { for(auto& a:*this){ a -= v;};  return (*this); }
	piece& operator *=(double t)        { for(auto& a:*this){ a *= t;};  return (*this); }
	piece& operator *=(int t)           { for(auto& a:*this){ a *= t;};  return (*this); }
	piece& operator /=(double t)        { return (*this)*=(1.0/t); }
	T sum() const                       { T sm=(T()*0.0); for(auto a:*this){ sm += a;}  return sm; }
	T avg() const                       { return sum()/size(); }

  //protected:
	T*  d;
	T*  dn;
};

//! lazyvec can be used to hold a contiguous array, similar to std::vallarray class but not that ugly constructor of vallarray
template<class T>
class lazyvec: public piece<T> 
{

 public:
	using piece<T>::d;
	using piece<T>::dn;

	lazyvec(): piece<T>(0,0) {};
	lazyvec(int siz): piece<T>(new T[siz],siz) {};
	lazyvec(size_t siz, const T& val): piece<T>(new T[siz],siz) {  std::fill(d, dn, val);}
	lazyvec(const piece<T>& v): piece<T>(new T[v.size()],v.size())  {  std::copy(v.d, v.dn, d); }
	lazyvec(const lazyvec& v): piece<T>(new T[v.size()],v.size())  {  std::copy(v.d, v.dn, d); }
	lazyvec(const std::vector<T>& v): piece<T>(new T[v.size()],v.size())  {  std::copy(&v[0], &*v.end(), d); }
	lazyvec(lazyvec&& v): piece<T>(v.d,v.size())  {  v.d=0; v.dn=0; }
	lazyvec(const T* dd, int nn): piece<T>(new T[nn],nn)            { std::copy(dd, dd+nn, d); }
	lazyvec(const T* dd, const T* de): piece<T>(new T[de-dd],de-dd)    { std::copy(dd, de, d);}
	lazyvec(const T* dd, const T* de, T(*func)(const T&)): piece<T>(new T[de-dd],de-dd)    { std::transform(dd, de, d, func);}
	~lazyvec(){ if(d) delete[] d; /*(std::cout<<'~').flush();*/ };

	lazyvec(lazyvec& v, bool move): piece<T>(v.d,v.size())  {  v.d=0; v.dn=0; assert(move); } //same as above but enforced by move

	void operator =(lazyvec&& v){ eat(v); };
	void eat(lazyvec& v)       { { if(d) delete[] d; }  d = v.d;  dn=v.dn; v.d=0; v.dn=0; /*(std::cout<<'$').flush();*/ };
	void operator =(const piece<T>& v){ { if(d) delete[] d; }  d = new T[v.size()]; std::copy(v.d, v.dn, d);  dn=d+v.size(); };
	void operator =(const lazyvec& v){ { if(d) delete[] d; }  d = new T[v.size()]; std::copy(v.d, v.dn, d);  dn=d+v.size();  /*(std::cout<<'&').flush();*/ };
	void operator =(const std::vector<T>& v){ { if(d) delete[] d; }  d = new T[v.size()]; std::copy(&v[0], &*v.end(), d);  dn=d+v.size(); };
	void operator =(const std::initializer_list<T>& v){ { if(d) delete[] d; }  d = new T[v.size()]; std::copy(v.begin(), v.end(), d);  dn=d+v.size(); };


	lazyvec& operator +=(const piece<T>& v)  { this->piece<T>::operator+=(v);  return (*this); }
	lazyvec& operator -=(const piece<T>& v)  { this->piece<T>::operator-=(v);   return (*this); }
	lazyvec& operator *=(const piece<T>& v)  { this->piece<T>::operator*=(v);   return (*this); }
	template<typename T2> lazyvec& operator *=(const piece<T2>& v)  { for(auto& a:*this){ a *= v[&a-d];};  return (*this); }
	template<typename T2> lazyvec& operator /=(const piece<T2>& v)  { for(auto& a:*this){ a /= v[&a-d];};  return (*this); }
	lazyvec& operator +=(T v)  { for(auto& a:*this){ a += v;};  return (*this); }
	lazyvec& operator -=(T v)  { for(auto& a:*this){ a -= v;};  return (*this); }
	lazyvec& operator *=(double t)             { for(auto& a:*this){ a *= t;};  return (*this); }
	lazyvec& operator *=(int t)             { for(auto& a:*this){ a *= t;};  return (*this); }
	lazyvec& operator /=(double t)  { return (*this)*=(1.0/t); }


	void resize(int nn)
	{ { if(d) delete[] d; }   if(nn) {d=new T[nn]; dn=d+nn;} else {d=0; dn=0;} }
	void resize(int nn,const T& val)
	{ { if(d) delete[] d; }   if(nn) {d=new T[nn]; dn=d+nn;  std::fill(d,dn, val); } else {d=0; dn=0;} }
	void pbak(T& vj)
	{	if(d){ T* od=d;  d=new T[dn+1-d]; std::copy(od,dn,d);
				 dn=d+dn-od+1;    *(dn-1)=vj;      delete[] od; }
		else { d=new T[1];      *d=vj;           dn=d+1;      }
	}
	void pbak(const piece<T> vs)
	{	if(d){ T* od=d;               d=new T[dn+vs.size()-d];       std::copy(od,dn, d);
				 dn=d+dn-od+vs.size();  std::copy(vs.d, vs.dn, dn-vs.size());  delete[] od; }
		else { d=new T[vs.size()];    std::copy(vs.d, vs.dn, d);          dn=d+vs.size(); }
	}
};

template<class T>	lazyvec<T>             operator -(const piece<T>& dis)   { lazyvec<T> tmp(dis);  for(auto& a:tmp){ a = -a;}  return tmp;  }
template<class T>	lazyvec<T>             operator +(const piece<T>& dis, const piece<T>& v)  { lazyvec<T> tmp(dis); return tmp+=v;  }
template<class T>	lazyvec<T>             operator -(const piece<T>& dis, const piece<T>& v)  { lazyvec<T> tmp(dis); return tmp-=v;  }
//template<class T>	lazyvec<T>             operator *(const piece<T>& dis, const piece<T>& v)  { lazyvec<T> tmp(dis); tmp*=v;  return tmp;  }
template<class T, typename T2> lazyvec<T2>  operator *(const T dis, const piece<T2>& v) { lazyvec<T2> tmp(v); return tmp*=dis;  }
template<class T> lazyvec<var3<T>>  operator *(const piece<T>& dis, const piece<var3<T>>& v) { lazyvec<var3<T>> tmp(v); return tmp*=dis;  }
//template<class T, typename T2> lazyvec<T>  operator *(const piece<T>& dis, const piece<T2>& v) { lazyvec<T> tmp(dis); return tmp*=v;  }
template<class T, typename T2> lazyvec<T>  operator /(const piece<T>& dis, const piece<T2>& v) { lazyvec<T> tmp(dis); return tmp/=v;  }
template<class T>	lazyvec<T>             operator +(const piece<T>& dis,T t)        { lazyvec<T> tmp(dis); tmp+=t;     return tmp; }
template<class T>	lazyvec<T>             operator -(const piece<T>& dis,T t)        { lazyvec<T> tmp(dis); tmp-=t;     return tmp; }
template<class T>	lazyvec<T>             operator *(const piece<T>& dis, double t)  { lazyvec<T> tmp(dis); tmp*=t;     return tmp; }
template<class T>	lazyvec<T>             operator *(const piece<T>& dis, int t)     { lazyvec<T> tmp(dis); tmp*=t;     return tmp; }
template<class T>	lazyvec<T>             operator /(const piece<T>& dis, double t)  { lazyvec<T> tmp(dis); tmp*=1.0/t; return tmp; }
template<class T>	lazyvec<T>             operator /(double t, const piece<T>& dis)  { lazyvec<T> tmp(dis); for(auto& a:tmp){ a = t/a;} return tmp; }
template<class T>	lazyvec<T>   mag(const piece<var3<T>>& dis)  { lazyvec<T> tmp(dis.size()); for(int i=0;i<dis.size();++i){ tmp[i] = mag(dis[i]); } return tmp; }



template<class T> lazyvec<T>  operator &(const piece<var3<T>>& dis, const piece<var3<T>>& v) { lazyvec<T> tmp(dis.size()); for(int i=0;i<dis.size();++i){ tmp[i] = dis[i]&v[i]; } return tmp; }



//! lazyvec but with a defult value to hold a uniform vector
//! a name and transformation information, used in #svplot only for now
template<class T>
class varsORv: public lazyvec<T>
{
  public:
	T dfult;
	std::string name;
	T Xa,  Xp; //! scale, shift
	using piece<T>::d;
	using piece<T>::dn;
	T(*transf)(T);
	T ax,  px;

	varsORv()                        :                      dfult(0)         , Xa(1)     , Xp(0)     , transf([](double a){return a;}), ax(1)     , px(0)     {}
	varsORv(const T& val)            :                      dfult(val)       , Xa(1)     , Xp(0)     , transf([](double a){return a;}), ax(1)     , px(0)     {}
	varsORv(size_t siz, const T& val): lazyvec<T>(siz,val), dfult(val)       , Xa(1)     , Xp(0)     , transf([](double a){return a;}), ax(1)     , px(0)     {}
	varsORv(const lazyvec<T>& vs)    : lazyvec<T>(vs),      dfult(vs.back()) , Xa(1)     , Xp(0)     , transf([](double a){return a;}), ax(1)     , px(0)     {}
	varsORv(const varsORv& vs)       : lazyvec<T>(vs),      dfult(vs.dfult)  , Xa(vs.Xa) , Xp(vs.Xp) , transf(vs.transf)              , ax(vs.ax) , px(vs.px) {}
	varsORv(const std::vector<T>& vs): lazyvec<T>(vs),      dfult(vs.dfult)  , Xa(vs.Xa) , Xp(vs.Xp) , transf(vs.transf)              , ax(vs.ax) , px(vs.px) {}
	varsORv(varsORv&& vs)            : lazyvec<T>(vs),      dfult(vs.dfult)  , Xa(vs.Xa) , Xp(vs.Xp) , transf(vs.transf)              , ax(vs.ax) , px(vs.px) {}
	varsORv(const T* dd, int nn)     : lazyvec<T>(dd,nn),   dfult(dd[nn-1])  , Xa(1)     , Xp(0)     , transf([](double a){return a;}), ax(1)     , px(0)     {}
	varsORv(const T* dd, const T* de): lazyvec<T>(dd,de),   dfult(*(de-1))   , Xa(1)     , Xp(0)     , transf([](double a){return a;}), ax(1)     , px(0)     {}
	varsORv(const T* dd, const T* de, T(*fn)(T)): lazyvec<T>(dd,de,fn)       , Xa(1)     , Xp(0)     , transf([](double a){return a;}), ax(1)     , px(0)     {}
	~varsORv(){};

	//varsORv(lazyvec<T>& vs, bool move): lazyvec<T>(vs,move)  {} //same as above but enforced by move
	void operator =(const piece<T>& v)  { { if(d) delete[] d; }  d = new T[v.size()]; std::copy(v.d, v.dn, d);  dn=d+v.size(); };
	void operator =(const lazyvec<T>& v){ { if(d) delete[] d; }  d = new T[v.size()]; std::copy(v.d, v.dn, d);  dn=d+v.size(); };
	void operator =(const varsORv& v)   { { if(d) delete[] d; }  d = new T[v.size()]; std::copy(v.d, v.dn, d);  dn=d+v.size(); dfult=v.dfult; Xa=v.Xa; Xp=v.Xp; transf=v.transf; ax=v.ax; px=v.px; };
	void operator =(const std::vector<T>& v){ { if(d) delete[] d; }  d = new T[v.size()]; std::copy(&v[0], &*v.end(), d);  dn=d+v.size(); };
	void operator =(const std::initializer_list<T>& v){ { if(d) delete[] d; }  d = new T[v.size()]; std::copy(v.begin(), v.end(), d);  dn=d+v.size(); };

	void operator =(lazyvec<T>&& v){ lazyvec<T>::eat(v); };

	T& operator [](int i)             { if(d+i<dn) return d[i]; else return dfult;};
	const T& operator [](int i) const { if(d+i<dn) return d[i]; else return dfult;};
	//T& operator ()(int i)             {return   (d+i<dn ? Xp+Xa*transf(d[i]+px) : dfult);};
	//const T& operator ()(int i) const {return   (d+i<dn ? Xp+Xa*transf(d[i]+px) : dfult);};
	T scalefrom01(T val) const {return   Xp+Xa*transf(val*ax+px);};

	void rescaledata(T d0, T del)	                       { Xp=d0; Xa=del; T* dd=d-1;  while(++dd<dn)  *dd = (*dd-d0)/del;  }
	void rescaledata(T d0, T del, T(*fn)(T))             { Xp=d0; Xa=del; T* dd=d-1;  while(++dd<dn)  *dd = fn((*dd-Xp)/Xa);  }
	void rescaledata(T d0, T del, T(*fn)(T), T p0, T a0) { Xp=d0; Xa=del; px=fn((p0-Xp)/Xa); ax=fn((p0+a0-Xp)/Xa)-px;   T* dd=d-1;
	                                                       while(++dd<dn)  *dd = (fn((*dd-Xp)/Xa)-px)/ax;  }
	//void transformdata(T d0, T del, T(*fn)(T), T p0)	{
		//Xp=d0; Xa=del; px=fn((p0-d0)/del); T* dd=d-1;  while(++dd<dn)	
			//*dd = fn((*dd-d0)/del)-px;}

};
//template<class T>
	//T rescaleval(T val, T d0, T del, T(*fn)(T), T p0, T a0) { T px=fn((p0-d0)/del); T ax=fn((p0+a0-d0)/del)-px;  return (fn((val-d0)/del)-px)/ax;}
inline double rescaleval(double val, double d0, double del, double(*fn)(double), double p0, double a0) { double px=fn((p0-d0)/del); double ax=fn((p0+a0-d0)/del)-px;  return (fn((val-d0)/del)-px)/ax;}

template<class T>  lazyvec<T> operator -(T v, lazyvec<T> vs)  { lazyvec<T> tmp(vs); for(auto& a:tmp){ a =v-a;};  return tmp; }

template<class T>  var3<T>    round(const var3<T>& vec) 	{ 	var3<T> res(vec); 	res.x=round(vec.x); res.y=round(vec.y); res.z=round(vec.z); 	return res; 	}
template<class T>  lazyvec<T> round(const lazyvec<T>& vec) 	{ 	lazyvec<T> res(vec); 	for(auto& v:res) v=round(v); 	return res; 	}

//! moves iss into oss then deletes iss, but efficiently
//template<class T> void transferto(lazyvec<lazyvec<T> >& oss, size_t i, lazyvec<lazyvec<T> >&& iss, size_t j=0)
//{	auto* os = oss.d+i-1;	auto* is = iss.d+j-1;	while ((++is)<iss.dn) { (++os)->d=is->d; os->dn=is->dn; is->d=0;	}  /*iss.d=0;*/	}




#define for_i(_vector_m_)  for(size_t i=0;i<_vector_m_.size();++i)
//#define for_e(_vector_m_)  for(const auto&e:_vector_m_)
#ifndef for_
 #define for_(_vector_m_, _i_m_)  for(size_t _i_m_=0;_i_m_<_vector_m_.size();++_i_m_)
#endif
#define fori0to(_n_m_)  for(int i=0;i<_n_m_;++i)
#define forin(_i_bgn_m,_i_end_m)  for(int i=_i_bgn_m;i<_i_end_m;++i)
#define forint(_ind_, _i_bgn_m,_i_end_m)  for(int _ind_=_i_bgn_m; _ind_<_i_end_m; ++_ind_)

using float2 = var2<float>;
using dbl2 = var2<double>;
using int2 = var2<int>;
using int2x2 = var2<var2<int> > ;

 template<size_t N> using  Dbl_ = std::array<double,N>;
 template<class T> using  vars = lazyvec<T>; // or std::vector<T>;
 using  dbls = vars<double>;
 using  floats = vars<float>;

 typedef   vars<dbl3>   dbl3s;
 using  ints = vars<int>;
 using  uchars = vars<unsigned char>;


//! lambda help
typedef double(*dblFunc)(const double&);  //[](const double& a){return a;}
typedef float(*floatFunc)(const float&);



#endif // TYPSES_H





#ifndef TYPSES_OERATIONS_H
#define TYPSES_OERATIONS_H



template<class C> int len(C vec) { return vec.size(); }

template<class T>
double sumdbl(const piece<T>& ps)  { double sm=1.0e-300; for(auto a:ps){ sm += a;}  return sm; }
template<class T>// error add dbg ize checks
double sumdbl(const piece<T>& ps, const piece<T>& ws)  { double sm=1.0e-300; const T* p = ps.d-1, *w=ws.d-1; while(++p<ps.dn){ sm += *w * *(++p);}  return sm; }
inline double sumdbl(const dbl3& ps, const dbl3& ws)  { return ps.x+ws.x + ps.y+ws.y + ps.z+ws.z; }
inline double sumdbl(const dbl3& ps)  { return ps.x+ ps.y+ ps.z; }

//! Note: these can be made more efficient if assuming non-empty vec
template<class T>                 T sumvars(const piece<T>& ps)  {
	T sm = (T()*0.0); const T* p = ps.d-1; while(++p<ps.dn){ sm += *p;}  return sm; }
template<typename T, typename T2> T sumvars(const piece<T>& ps, const piece<T2>& ws)  {
	T sm = (T()*0.0); const T* p = ps.d-1; const T2 *w=ws.d-1; while(++p<ps.dn){ sm += *(++w) * *p;}  return sm; }
	
template<class T>                 T sumSqrs(const piece<T>& ps)  {
	T sm = (T()*0.0); const T* p = ps.d-1; while(++p<ps.dn){ sm += *p * *p;}  return sm; }
template<typename T, typename T2> T sumSqrs(const piece<T>& ps, const piece<T2>& ws)  {
	T sm = (T()*0.0); const T* p = ps.d-1; const T2 *w=ws.d-1; while(++p<ps.dn){ sm += *(++w) * *p * *p;}  return sm; }


template<class T> T sum(const piece<T>& pis)  { return sumvars(pis); }
template<class T>	lazyvec<T>   sum(piece<piece<T>> ps)  { lazyvec<T> sm(ps[0]); while(++ps.d<ps.dn){ sm += *ps.d;}  return sm;  } // does not work for size zero
template<class T>	T avg(const piece<T>& dis)       { return sum(dis)/dis.size(); }

template<class T> T min(const piece<T>& pis)  { return *std::min_element(pis.d,pis.dn); }
template<class T> T max(const piece<T>& pis)  { return *std::max_element(pis.d,pis.dn); }
template<class T, typename TFunc> void transform(piece<T>& pis, TFunc func)  { for(T& dr:pis) dr=func(dr); }
template<typename T> void transform(T* d, const T* dn, T(*func)(const T&))  { --d; while (++d<dn) *d=func(*d); }
//template<class T> void transform(piece<T>& pis, float(*func)(float&))  { for(float& dr:pis) dr=func(dr); }

//double (*func)(const double&) = ([](const double& d) {return d;}) // this is just  a note/reminder on old  function / lambda...


template<typename Type>
std::vector<Type>& operator *= (std::vector<Type>& vs, double scale)
{ 	for(Type& v : vs) v*=scale; 		return vs; 	}


template<class T> std::ostream&  operator << (std::ostream& out, const var3<T>& pos)
{
	//std::ios_base::fmtflags flgs=out.flags();
	//out.flags(std::ios::showpoint | std::ios::scientific);out<< std::setprecision(5);
	out   << pos.x   <<" " << pos.y   <<" " << pos.z;
	//out.flags(flgs);
	return out;
}
template<class T> std::ofstream& operator << (std::ofstream& out, const var3<T>& pos)
{
	std::ios_base::fmtflags flgs=out.flags();
	//out.flags(std::ios::showpoint | std::ios::scientific);	 << std::setprecision(8)
	out  << pos.x   <<' ' << pos.y   <<' ' << pos.z;
	out.flags(flgs);
	return out;
}


template<class T> std::istream&                    operator >> (std::istream& in, var3<T>& pos)
{	in >> pos.x >> pos.y >> pos.z;	return in;	}

//inline std::istream&                               operator >> (std::istream& in, int3& i3)
//{    in >> i3[0] >> i3[1] >> i3[2];  return in;	}


template<typename T1, typename T2> std::istream&   operator >> (std::istream& in, std::pair<T1,T2>& ab)
{ 	in >> ab.first >> ab.second;      return in; 	}

template<typename T1, typename T2>  std::ostream&  operator << (std::ostream& out, const std::pair<T1,T2>& ab)
{ 	out << ab.first<<" "<< ab.second; return out; 	}

template<typename T>    std::istream&              operator >> (std::istream& in, var2<T>& ab)
{ 	in >> ab.a >> ab.b;     return in; 	}

template<typename T>  std::ostream&                operator << (std::ostream& out, const var2<T>& ab)
{ 	out << ab.a<<" "<< ab.b; return out; 	}

template<typename T>  std::ostream &               operator << (std::ostream & out, const std::vector<T>& vec)
{
	if(vec.size() && vec.size()<10)  for (auto v : vec) out << v << '\t';
	else                             for (auto v : vec) out << v << '\n';
	return out;
}

template<typename T1, typename T2> std::ostream & operator << (std::ostream & out, const std::map<T1,T2>& vec)
{ 	for (auto v : vec) out << v << '\n'; 	return out; 	}


template<typename T>  std::ostream &              operator << (std::ostream & out, const piece<T>& vec)
{
	if(vec.size() && vec.size()<10 ) for (auto v : vec) out << v << '\t';
	else                             for (auto v : vec) out << v << '\n';
	return out;
}


template<typename T> std::istream & operator >> (std::istream & inp, piece<T> tabl) // Warning relying on copy constructor not to deep copy
{
	for (; tabl.d<tabl.end();++tabl.d)   inp >> *tabl.d;
	//if (inp.fail()) std::cout<<"Error while reading array"<<std::endl;
	return inp;
}

//std::ostream & operator << (std::ostream & out, const std::valarray<std::valarray<double> >& vecvec)
//{
	//for (size_t i=0; i<vecvec[0].size();++i)
	//{ for (size_t j=0; j<vecvec.size();++j) out << vecvec[j][i] << ' ';   out << '\n'; }
	//out << '\n';
	//return out;
//}


inline std::istream & operator >> (std::istream & inp, std::vector< std::pair<double, double> >& tabl)
{
	for (size_t i=0; i<tabl.size();++i)
		inp >> tabl[i].first >> tabl[i].second ;
	//if (inp.fail()) std::cout<<"Error while reading array"<<std::endl;
	return inp;
}



inline void  replaceInFromTo(std::string& str, const std::string& frm, const std::string & to)
{	//return str = std::regex_replace( str, std::regex(frm), to );
	size_t po = 0;
	while((po = str.find(frm, po)) != std::string::npos) {
		str.replace(po, frm.length(), to);  po += to.length(); }
}
inline void  replaceInBetweenTo(std::string& str, const std::string& bgn, const std::string& end, const std::string & to)
{	//no greedy macth
	size_t po = 0, po0=0;
	while((po = str.find(bgn, po)) != std::string::npos) { po0=po; po=po+bgn.length();
		while((po = str.find(end, po)) != std::string::npos) {
			str.replace(po0, po-po0, to); po = po0+to.length();  break; } }
}

inline std::string replaceFromTo(std::string  str, const std::string& frm, const std::string & to)
{
	return std::regex_replace(str, std::regex(frm), to);
	//size_t po = 0;
	//while((po = str.find(frm, po)) != std::string::npos) { str.replace(po, frm.length(), to);  po += to.length();   }
	//return str;
}

inline std::string baseName(const std::string& fName)
{
    if(fName.find_last_of('.') != std::string::npos)
        return fName.substr(0,fName.find_last_of('.'));
    return fName;
}



template<typename T> bool _1At(T n, int ibgn)  {  return n&(1<<ibgn);  }
inline               bool _1At(unsigned int n, int ibgn)  {  return n&(1<<ibgn);  }// for debugger
template<typename T> bool _1In(T n, int ibgn, int lnt)  {  return (n>>ibgn)&((1<<lnt)-1);  }





template<typename T, template<typename ...> class C1, template<typename ...> class C2>
vars<vars<T> > transpose(const C1<C2<T> >& vecvec)
{
	if(!vecvec.size()) return vars<vars<T> >();
	vars<vars<T> > trans(vecvec[0].size(),vars<T>(vecvec.size()));
	for (size_t i=0; i<vecvec[0].size();++i)
		for (size_t j=0; j<vecvec.size();++j) trans[i][j] = vecvec[j][i] ;
	return trans;
}

template<typename T, template<typename ...> class C>
piece<T> allof(C<T>& vs) {return piece<T>(&*vs.begin(), &*vs.end());}



//! class table IO, usage: cout<< tableIO(stdVecVec)  (not TableIO)
template<typename T, template<typename ...> class C1, template<typename ...> class C2>
class TableIO
{ public:
	TableIO(const C1<C2<T> >& vecvec, std::vector<std::string> hdrs, bool transpose, char sepr) :  vss_(vecvec),  hds_(hdrs),  transpose_(transpose), sep_(sepr) {};
	const C1<C2<T> >&   vss_;
	std::vector<std::string> hds_;
	const char        transpose_;
	const char        sep_;
};
template<typename T, template<typename ...> class C1, template<typename ...> class C2>
TableIO<T,C1,C2> tableIO(const C1<C2<T> >& vecvec, std::vector<std::string> hdrs=std::vector<std::string>(), bool transpose=true, char sepr='\t') {return TableIO<T,C1,C2>(vecvec,hdrs,transpose,sepr);}

template<typename T, template<typename ...> class C1, template<typename ...> class C2>
std::ostream & operator << (std::ostream & out, const TableIO<T,C1,C2>& tbl)
{
	if(tbl.hds_.size()==tbl.vss_.size()) { for (size_t j=0; j<tbl.vss_.size();++j) out << tbl.hds_[j]<<tbl.sep_;  }   	out<<'\n';
	for (size_t i=0; i<tbl.vss_[0].size();++i)
	{ for (size_t j=0; j<tbl.vss_.size();++j) out<<std::setw(12)<<std::left<< tbl.vss_[j][i]<<tbl.sep_;   out<<'\n'; }
	out << '\n';
	return out;
}
template<typename T, template<typename ...> class C1, template<typename ...> class C2>
inline std::istream & operator >> (std::istream & inp, TableIO<T,C1,C2>& tbl)
{
	inp >> std::ws;
	std::string row, item;
	if( std::is_arithmetic<T>::value && isalpha(inp.peek()))
	{	std::stringstream ss( row );
		while ( getline ( ss, item, tbl.sep_ ) ) tbl.hds_.push_back(item);
	}
	std::vector<std::vector<T> > res;
	while( getline( inp, row ) )
	{	std::vector<T> R;	std::stringstream ss( row );
		while ( getline ( ss, item, tbl.sep_ ) ) R.push_back( strTo<T>(item) );
		res.push_back( R );
	}
	if(res.size()) tbl.vss_.resize(res[0].size());   else return inp;
	for (size_t i=0; i<tbl.vss_.size();++i)
	{	tbl.vss_[i].resize(res.size());
		if(res[i].size()!=tbl.vss_.size()) std::cout<<"Error table sizes don't match"<<std::endl;
		for (size_t j=0; j<tbl.vss_[i].size();++j) tbl.vss_[i][j]=res[j][i];
	}
	return inp;
}




// CLEAN UP ME PLEASE:



//template<class T> piece<T>   insidvars(piece<T>& vs) {return piece<T>(vs.d+1, vs.dn-1);}
template<class T> lazyvec<T>  diffVars(piece<T>& vs) { lazyvec<T> res(vs);  piece<T>(res.d+1,res.dn)-=vs; res[0]=res[1]; return res;}





#ifdef VMMLIB__VECTOR__HPP
 template< size_t M, typename T >
 Vctr<M,T>::Vctr(const var3<T>& v3) 	{ array[ 0 ] = v3.x;	array[ 1 ] = v3.y;	array[ 2 ] = v3.z; 	}
#endif







template<class T>
vars<dbls> distribution(const piece<T>  xs, const piece<T> ws, int nBins=64)
{
	vars<dbls> distrib(3,dbls(nBins, 0.0));
	double minU=min(xs);
	double maxU=max(xs);
	double deltaU=(maxU-minU)/nBins+1.0e-72;

	for (int i=0; i<nBins; ++i)	distrib[0][i] = minU+deltaU/2+i*deltaU;

	for (int i=0; i<ws.size(); ++i)
	{
		int distInd=std::min(int((xs[i]-minU)/deltaU+0.5),nBins-1);
		++distrib[1][distInd];
		distrib[2][distInd]+=ws[i];
	}

	distrib[1]/=distrib[1].sum()*deltaU;
	distrib[2]/=distrib[2].sum()*deltaU;

	return distrib;
}




inline double linearInterpolate(double x, double x1, double x2, double y1, double y2) 	 { 	 return y1+(x-x1)*(y2-y1)/(x2-x1); 	 }

inline double averageCDF(double xMin, double xMax, const std::vector< std::pair<double, double> >& tabl)
{/// used to get average of a part of a distribution
	std::vector< std::pair<double, double> >::const_iterator itr = tabl.begin();
	while((++itr)->first <=  xMin && itr != tabl.end()) ;

	if( !(itr->first >= xMin && (itr-1)->first <=  xMin) )
	{std::cout<<"\n LookUp Error:  "	<< xMin-1.0<<"     " << tabl.begin()->first-1.0<<"  "  << tabl.begin()->second<<"    " << tabl.rbegin()->first-1.0<<"  " << tabl.rbegin()->second
			 <<"\n\n "<<(itr == tabl.end())<<std::endl; exit(-1);}

	assert(itr->first > xMin && (itr-1)->first <=  xMin);
	
	double wSum=itr->first-xMin;
	double wxy1=(itr->first-xMin)*0.5*(itr->second+linearInterpolate(xMin,(itr-1)->first, itr->first, (itr-1)->second, itr->second));
	
	while((++itr)->first <=  xMax && itr != tabl.end())
	{
		wSum+=itr->first-(itr-1)->first;
		wxy1+=(itr->first-(itr-1)->first)*0.5*(itr->second+(itr-1)->second);
		//cout<<itr->first<<": "<<itr->second<<"     "<<0.5*(itr->second+(itr-1)->second)<<" s/ "<<(itr->first-(itr-1)->first)<<"   =  "<<wxy1/wSum<<std::endl;
	}
	if(itr != tabl.end())
	{
		wSum+=itr->first-xMax;
		wxy1+=(itr->first-xMax)*0.5*(itr->second+linearInterpolate(xMin,(itr-1)->first, itr->first, (itr-1)->second, itr->second));
	}


	return wxy1/(wSum+1.0e-64);

}


#endif //TYPSES_OERATIONS_H




/*
//- Debugging:
 Edit file: /usr/share/gcc-5/python/libstdcxx/v6/printers.py, for gdb pretty printing
class varsPrinter:
    "Print a vars"

    class _iterator(Iterator):
        def __init__ (self, start, finish):
            self.item = start
            self.finish = finish
            self.count = 0

        def __iter__(self):
            return self

        def __next__(self):
            count = self.count
            self.count = self.count + 1

            if self.item == self.finish:
               raise StopIteration
            elt = self.item.dereference()
            self.item = self.item + 1
            return ('[%d]' % count, elt)

    def __init__(self, typename, val):
        self.typename = typename
        self.val = val

    def children(self):
        return self._iterator(self.val['d'], self.val['dn'])

    def to_string(self):
        start = self.val['d']
        finish = self.val['dn']
        return ('%s of length %d' % (self.typename, int (finish - start)))

    def display_hint(self):
        return 'array'


in build_libstdcxx_dictionary, add:
    libstdcxx_printer.add_container('', 'piece', varsPrinter)
    libstdcxx_printer.add_container('', 'lazyvec', varsPrinter)
*/


