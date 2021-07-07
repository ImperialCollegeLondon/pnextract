#ifndef TYPSES_H
#define TYPSES_H


// Convinience vector classes used by network extraction, flow simulation 
// and other codes developed by Ali Qaseminejad Raeini.
// Main template classes defined here: var3, var2, piece and Vars 


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
	#define verySmall  1e-31
	#define   maxT(Typ)  (std::numeric_limits<Typ>::max())
	#define   minT(Typ)  (std::numeric_limits<Typ>::min())
	#define  iminT(Typ) (int(std::numeric_limits<Typ>::min()))
	#define  dminT(Typ) (double(std::numeric_limits<Typ>::min()))
	#define  imaxT(Typ) (int(std::numeric_limits<Typ>::max()))
	#define  fmaxT(Typ) (float(std::numeric_limits<Typ>::max()))
	#define  dmaxT(Typ) (double(std::numeric_limits<Typ>::max()))
	#define  epsT(Typ)  std::numeric_limits<Typ>::epsilon()
#endif

#if (__cplusplus < 201403L)
	template< bool B, class T = void >
	using enable_if_t = typename std::enable_if<B,T>::type;
#else
	using std::enable_if_t;
#endif


#ifdef OpenMP
 #ifdef _debugCompile_
  #define OMPragma(_args_) 
 #else
  #define OMPragma(_args_) _Pragma(_args_)
 #endif
#else
 #define OMPragma(_args_)
#endif
#define OMPFor(_args_)  OMPragma(TOSTRING(omp parallel for _args_))


#define for_i(_vector_m_)  for(size_t i=0; i<_vector_m_.size(); ++i)
#ifndef for_
 #define for_(_vector_m_, _i_m_)  for(size_t _i_m_=0;_i_m_<_vector_m_.size(); ++_i_m_)
#endif
#define fori0to(_n_m_)  for(int i=0; i<_n_m_; ++i)
#define forin(_i_bgn_m,_i_end_m)  for(int i=_i_bgn_m; i<_i_end_m; ++i)
#define forint(_ind_, _i_bgn_m,_i_end_m)  for(int _ind_=_i_bgn_m; _ind_<_i_end_m; ++_ind_)



constexpr double PI = 3.14159265358979;


#include "globals.h"


//! 3D vector class template
template<class T>
class var3  {
 public:
	T	x, y, z;
	var3() = default;// Warning: not zero initialized in opt mode
	template< typename std::enable_if<std::is_arithmetic<T>::value,int>::type = 0> 
	explicit var3(T r)        { x = r;     y = r;     z = r; } // use this to zero-initialize 
	var3(T r, T s, T t)       { x = r;     y = s;     z = t; }
	explicit var3(const T* d)          { x = d[0];  y = d[1];  z = d[2]; }
	template<class U>
	var3(var3<U> n)           { x = n.x;   y = n.y;   z = n.z; }
	#ifdef VMMLIB__VECTOR__HPP
	var3(const Vctr<3,T>& v3) { x = v3[0];	y = v3[1]; z = v3[2]; }
	#endif

	T&       operator [](long k)        { return ((&x)[k]); }
	const T& operator [](long k) const  { return ((&x)[k]); }
	T _0()                       const  { return x; }
	T _1()                       const  { return y; }
	T _2()                       const  { return z; }


    //operator  int() const { return x; }
    //operator  double() const { return x; }

	var3&  operator += (const var3& v)  { x += v.x;  y += v.y;  z += v.z;  return  *this; }
	var3&  operator -= (const var3& v)  { x -= v.x;  y -= v.y;  z -= v.z;  return  *this; }
	var3&  operator += (const T& t)     { x += t;    y += t;    z += t;    return  *this; } // clumsy
	var3&  operator -= (const T& t)     { x -= t;    y -= t;    z -= t;    return  *this; } // clumsy
	var3&  operator *= (const double t) { x *= t;    y *= t;    z *= t;    return  *this; }
	var3&  operator /= (const double t) { return (*this *= 1./t);  }
	var3&  operator ^= (const var3& v)  { double r=y*v.z-z*v.y, s=z*v.x-x*v.z;  z=x*v.y-y*v.x;  x=r; y=s; 	return *this; }
	var3&  operator *= (const var3& v)  { x *= v.x;  y *= v.y;  z *= v.z;  return *this; }
	var3   operator -  (void)          const { return  var3(-x, -y, -z); }
	var3   operator +  (const var3& v) const { return  var3(x+v.x, y+v.y, z+v.z); }
	var3   operator -  (const var3& v) const { return  var3(x-v.x, y-v.y, z-v.z); }
	var3   operator +  (const T& t)    const { return  var3(x+t, y+t, z+t); } // clumsy
	var3   operator -  (const T& t)    const { return  var3(x-t, y-t, z-t); } // clumsy
	var3   operator *  (const double t)const { return  var3(x*t, y*t, z*t); }
	var3   operator /  (const double t)const { double  f = 1./t;  return var3(x*f, y*f, z*f); }
	double operator &  (const var3& v) const { return  x*v.x+y*v.y+z*v.z; }
	var3   operator ^  (const var3& v) const { return  var3(y*v.z-z*v.y,  z*v.x-x*v.z,  x*v.y-y*v.x); }
	var3   operator *  (const var3& v) const { return  var3(x*v.x, y*v.y, z*v.z); }
	var3   operator /  (const var3& v) const { return  var3(x/v.x, y/v.y, z/v.z); }
	bool   operator == (const var3& v) const { return  (x-v.x)*(x-v.x) + (y-v.y)*(y-v.y) + (z-v.z)*(z-v.z) < verySmall; }
	bool   operator != (const var3& v) const { return  !(v==*this); }
};

typedef  var3<int>        int3;
typedef  var3<var3<int> > int3x3;
typedef  var3<float>      float3;
typedef  var3<double>     dbl3;

template<class T>  var3<T> rotateAroundLine(var3<T> b, double gamma,  var3<T> n, var3<T> a)  {
	//! rotate y around line passing through x, in the direction of n, http://inside.mines.edu/~gmurray/ArbitraryAxisRotation
	double s = sinf(gamma),   c = cosf(gamma),  nb = n.x*b.x + n.y*b.y +n.z*b.z, lc = 1.-c;
	return var3<T>(
	 	( a.x*(n.y*n.y+n.z*n.z) - n.x*( a.y*n.y+a.z*n.z - nb ) )*lc + b.x*c + (-a.z*n.y+a.y*n.z-n.z*b.y+n.y*b.z )*s,
		( a.y*(n.x*n.x+n.z*n.z) - n.y*( a.x*n.x+a.z*n.z - nb ) )*lc + b.y*c + ( a.z*n.x-a.x*n.z+n.z*b.x-n.x*b.z )*s,
		( a.z*(n.x*n.x+n.y*n.y) - n.z*( a.x*n.x+a.y*n.y - nb ) )*lc + b.z*c + (-a.y*n.x+a.x*n.y-n.y*b.x+n.x*b.y )*s );
}
template<class T>  var3<T> rotateAroundVec(const var3<T> b, double gamma, var3<T> n)  {
	//! Rotate b around line in the direction of n passing through centre, http://inside.mines.edu/~gmurray/ArbitraryAxisRotation
	double sg = sinf(gamma),   cg = cosf(gamma),  nb = (n.x*b.x + n.y*b.y + n.z*b.z)*(1.-cg);
	return var3<T>( n.x*nb + b.x*cg + (n.y*b.z-n.z*b.y)*sg,
	                n.y*nb + b.y*cg + (n.z*b.x-n.x*b.z)*sg,
	                n.z*nb + b.z*cg + (n.x*b.y-n.y*b.x)*sg );
}




//! 2D vector class template
template<class T>
class var2  {
 public:
	T a;//union {T	first;  T x; };
	T b;//union {T	second; T y; };
	var2() = default;//not zero initialized 
	var2(T r, T s)          { a = r;  b = s; }
	var2(int r)             { a = r;  b = r; }  // use this to zero-initialize //ERROR REMOVE ME
	var2(const T* d)        { a = d[0];  b = d[1]; }
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

	var2&  operator += (const var2& v)  { a += v.b;  b += v.b;  return *this; }
	var2&  operator -= (const var2& v)  { a -= v.b;  b -= v.b;   return *this; }
	var2&  operator *= (const double t) { a *= t;  b *= t;   return *this; }
	var2&  operator /= (const double t) { double f = 1./t;  a *= f;  b *= f;  return *this; }
	var2&  operator *= (const var2& v)  { a *= v.b;  b *= v.b;   return *this; }
	var2   operator -  (void) const          { return (var2(-a, -b)); }
	var2   operator +  (const var2& v) const { return (var2(a+v.b, b+v.b)); }
	var2   operator -  (const var2& v) const { return (var2(a-v.b, b-v.b)); }
	var2   operator *  (const double t)const { return (var2(a*t, b*t)); }
	var2   operator /  (const double t)const { double f = 1./t;  return (var2(a*f, b*f)); }
	double operator &  (const var2& v) const { return (a*v.b+b*v.b); }
	var2   operator *  (const var2& v) const { return (var2(a*v.b, b*v.b)); }
	var2   operator /  (const var2& v) const { return (var2(a/v.b, b/v.b)); }
	bool   operator == (const var2& v) const { return ((a-v.b)*(a-v.b) < verySmall) && ((b-v.b)*(b-v.b) < verySmall) ; }
	bool   operator != (const var2& v) const { return ((a-v.b)*(a-v.b) >= verySmall) || ((b-v.b)*(b-v.b) >= verySmall) ; }
};

typedef   var2<float> float2;
typedef   var2<double> dbl2;
typedef   var2<int> int2;

template<class T>  var3<T> operator *(double t, const var3<T>& v) { return var3<T>(t*v.x, t*v.y, t*v.z); }
inline  dbl3   operator *(const int3& v, const dbl3& d) { return dbl3(v.x*d.x, v.y*d.y, v.z*d.z); } // boost int to T
template<class T>  T       mag(const var3<T>& v)        { return std::sqrt(v.x*v.x+v.y*v.y+v.z*v.z); }
template<class T>  T       sum(const var3<T>& v)        { return (v.x+v.y+v.z); }
template<class T>  double  magSqr(const var3<T>& v)     { return (v.x*v.x+v.y*v.y+v.z*v.z); }
template<class T>  var3<T> norm(const var3<T>& v)       { return  v/(mag(v)+1e-300); }


//. component-wise max and min
template<class T> var3<T> max(const var3<T>& a, const var3<T>& b)  { return var3<T>(std::max(a.x,b.x),std::max(a.y,b.y),std::max(a.z,b.z)); }
template<class T> var3<T> min(const var3<T>& a, const var3<T>& b)  { return var3<T>(std::min(a.x,b.x),std::min(a.y,b.y),std::min(a.z,b.z)); }

template<class T>  T  toScalar(const T& v)       { return v; } // NOW IMPLICITly CONVERTED
template<class T>  T  toScalar(const var3<T>& v) { return mag(v); }





//! a piece of a contiguous array, used for efficient and prettier pass of array contents than using iterators, similar to std::string_view
template<class T>
class piece  {
 public:

	piece() noexcept: d(0), dn(0) {};
	piece(T* dd, int n) noexcept: d(dd), dn(d+n)     {};
	piece(T* dd, T* de) noexcept: d(dd), dn(de)      {};
	piece(const piece& p) noexcept: d(p.d), dn(p.dn) {};//! Note this is different from = operator note data hold by piece are not const unless piece is const itself
	piece(std::vector<T>& vs) noexcept: d(&vs[0]), dn(d+vs.size()) {};
	void reset(T* dd, int n)     { d=dd; dn=d+n; };
	void reset(const piece& vs)     { d=&vs[0]; dn=d+vs.size(); };//! note data hold by piece are not const unless piece is const itself
	void reset(std::vector<T>& vs)  { d=&vs[0]; dn=d+vs.size(); };

	T* begin() const { return d; };
	T* end()   const { return dn; };
	const T& back()   const { return *(dn-1); };
	const T* cbegin() const { return d; };
	const T* cend()   const { return dn; };
	T* operator ()()   const { return d; };
	//const T& operator [](int i) const { return d[i]; };
	//const T* operator ()() const { return d; };
	T& operator [](int i) const { return d[i]; };
	size_t size() const        { return dn-d; };
	size_t capacity()          { return dn-d; }
	T* data() { return d; }
	const T* data() const { return d; }

	piece& fill(const T& v)  { std::fill(d, dn, v);  return *this; }
	piece& operator =(const piece& v)  { ensure(size()==v.size(), "in operator =, piece sizes should be the same"); std::copy(v.d, v.dn, d);  return *this; }
	piece& operator +=(const piece& v)  { for(auto& a:*this){ a += v[&a-d]; };  return *this; }
	piece& operator -=(const piece& v)  { for(auto& a:*this){ a -= v[&a-d]; };  return *this; }
	piece& operator *=(const piece& v)  { for(auto& a:*this){ a *= v[&a-d]; };  return *this; }
	piece& operator /=(const piece& v)  { for(auto& a:*this){ a /= v[&a-d]; };  return *this; }
	piece& operator +=(T v)             { for(auto& a:*this){ a += v; };  return *this; }
	piece& operator -=(T v)             { for(auto& a:*this){ a -= v; };  return *this; }
	piece& operator *=(double t)        { for(auto& a:*this){ a *= t; };  return *this; }
	piece& operator *=(int t)           { for(auto& a:*this){ a *= t; };  return *this; }
	piece& operator /=(double t)        { return (*this)*=(1./t); }
	T sum() const                       { T sm=(T()*0.); for(auto a:*this){ sm += a; }  return sm; }
	T avg() const                       { return sum()/size(); }
	//T avgdbl() const                       { return sumdbl()/size(); }

  //protected:
	T*  d;
	T*  dn;
};

//! Vars can be used to hold a contiguous array, similar to std::vallarray class but without the ugly constructor of vallarray
template<class T>
class Vars: public piece<T>  {

 public:
	using piece<T>::d;
	using piece<T>::dn;

	Vars(): piece<T>(0,0) {};
	Vars(int siz): piece<T>(new T[siz],siz) {};
	Vars(size_t siz, const T& val): piece<T>(new T[siz],siz) {  std::fill(d, dn, val); }
	Vars(const piece<T>& v): piece<T>(new T[v.size()],v.size())  {  std::copy(v.d, v.dn, d); }
	Vars(const Vars& v): piece<T>(new T[v.size()],v.size())   {  std::copy(v.d, v.dn, d); }
	Vars(const std::vector<T>& v): piece<T>(new T[v.size()],v.size())  {  std::copy(&v[0], &*v.end(), d); }
	Vars(Vars&& v)  noexcept: piece<T>(v.d,v.size()) {  v.d=0; v.dn=0; }
	Vars(const T* dd, int nn): piece<T>(new T[nn],nn)            { std::copy(dd, dd+nn, d); }
	Vars(const T* dd, const T* de): piece<T>(new T[de-dd],de-dd)    { std::copy(dd, de, d); }
	Vars(const T* dd, const T* de, T(*func)(const T&)): piece<T>(new T[de-dd],de-dd)    { std::transform(dd, de, d, func); }
	~Vars(){ if(d) delete[] d; /*(std::cout<<'~').flush();*/ };

	Vars(Vars& v, bool move): piece<T>(v.d,v.size())  {  v.d=0; v.dn=0; assert(move); } //same as above but enforced by move

	void operator =(Vars&& v){ eat(v); };
	void eat(Vars& v)       { { if(d) delete[] d; }  d = v.d;  dn=v.dn; v.d=0; v.dn=0; /*(std::cout<<'$').flush();*/ };
	void operator =(const piece<T>& v){ { if(d) delete[] d; }  d = new T[v.size()]; std::copy(v.d, v.dn, d);  dn=d+v.size(); };
	void operator =(const Vars& v){ { if(d) delete[] d; }  d = new T[v.size()]; std::copy(v.d, v.dn, d);  dn=d+v.size();  /*(std::cout<<'&').flush();*/ };
	void operator =(const std::vector<T>& v){ { if(d) delete[] d; }  d = new T[v.size()]; std::copy(&v[0], &*v.end(), d);  dn=d+v.size(); };
	void operator =(const std::initializer_list<T>& v){ { if(d) delete[] d; }  d = new T[v.size()]; std::copy(v.begin(), v.end(), d);  dn=d+v.size(); };


	Vars& operator +=(const piece<T>& v)  { this->piece<T>::operator+=(v);  return *this; }
	Vars& operator -=(const piece<T>& v)  { this->piece<T>::operator-=(v);   return *this; }
	Vars& operator *=(const piece<T>& v)  { this->piece<T>::operator*=(v);   return *this; }
	template<class U> Vars& operator *=(const piece<U>& v)  { for(auto& a:*this){ a *= v[&a-d]; };  return *this; }
	template<class U> Vars& operator /=(const piece<U>& v)  { for(auto& a:*this){ a /= v[&a-d]; };  return *this; }
	Vars& operator +=(T v)       { for(auto& a:*this){ a += v; };  return *this; }
	Vars& operator -=(T v)       { for(auto& a:*this){ a -= v; };  return *this; }
	Vars& operator *=(double t)  { for(auto& a:*this){ a *= t; };  return *this; }
	Vars& operator *=(int t)     { for(auto& a:*this){ a *= t; };  return *this; }
	Vars& operator /=(double t)  { return (*this)*=(1./t); }


	void resize(int nn)  {
		{ if(d) delete[] d; }   if(nn) {d=new T[nn]; dn=d+nn; } else {d=0; dn=0; } }
	void resize(int nn,const T& val)  {
		{ if(d) delete[] d; }   if(nn) {d=new T[nn]; dn=d+nn;  std::fill(d,dn, val); } else {d=0; dn=0; } }
	void pbak(T& vj) {// inefficient, use std::vector instead for dyn_array
		if(d){ T* od=d;  d=new T[dn+1-d]; std::copy(od,dn,d);
				 dn=d+dn-od+1;    *(dn-1)=vj;      delete[] od; }
		else { d=new T[1];      *d=vj;           dn=d+1;      }
	}
	void pbak(const piece<T> vs){
		if(d){ T* od=d;               d=new T[dn+vs.size()-d];       std::copy(od,dn, d);
				 dn=d+dn-od+vs.size();  std::copy(vs.d, vs.dn, dn-vs.size());  delete[] od; }
		else { d=new T[vs.size()];    std::copy(vs.d, vs.dn, d);          dn=d+vs.size(); }
	}
};

typedef   Vars<double>         dbls;
typedef   Vars<float>          floats;
typedef   Vars<dbl3>           dbl3s;
typedef   Vars<int>            ints;
typedef   Vars<unsigned char>  uchars;
template<class T>  using  vars = Vars<T>; // or std::vector<T>;

//! lambda help
typedef double(*dblFunc)(const double&);  //[](const double& a){ return a; }
typedef float(*floatFunc)(const float&);

template<size_t N> using  Dbl_ = std::array<double,N>;


template<class T>	Vars<T>        operator -(const piece<T>& dis)   { Vars<T> rt(dis);  for(auto& a:rt){ a = -a; }  return rt;  }
template<class T>	Vars<T>        operator +(const piece<T>& dis, const piece<T>& v)  { Vars<T> rt(dis); return rt+=v;  }
template<class T>	Vars<T>        operator -(const piece<T>& dis, const piece<T>& v)  { Vars<T> rt(dis); return rt-=v;  }
//template<class T>	Vars<T>        operator *(const piece<T>& dis, const piece<T>& v)  { Vars<T> rt(dis); rt*=v;  return rt;  }
template<class T, class U> Vars<U> operator *(const T dis, const piece<U>& v) { Vars<U> rt(v); return rt*=dis;  }
template<class T> Vars<var3<T>>    operator *(const piece<T>& dis, const piece<var3<T>>& v) { Vars<var3<T>> rt(v); return rt*=dis;  }
//template<class T, typename U> Vars<T>  operator *(const piece<T>& dis, const piece<U>& v) { Vars<T> rt(dis); return rt*=v;  }
template<class T, class U> Vars<T> operator /(const piece<T>& dis, const piece<U>& v) { Vars<T> rt(dis);  return rt/=v;  }
template<class T>	Vars<T>        operator +(const piece<T>& dis,T t)        { Vars<T> rt(dis); rt+=t;    return rt; }
template<class T>	Vars<T>        operator -(const piece<T>& dis,T t)        { Vars<T> rt(dis); rt-=t;    return rt; }
template<class T>	Vars<T>        operator *(const piece<T>& dis, double t)  { Vars<T> rt(dis); rt*=t;    return rt; }
template<class T>	Vars<T>        operator *(const piece<T>& dis, int t)     { Vars<T> rt(dis); rt*=t;    return rt; }
template<class T>	Vars<T>        operator /(const piece<T>& dis, double t)  { Vars<T> rt(dis); rt*=1./t; return rt; }
template<class T>	Vars<T>        operator /(double t, const piece<T>& dis)  { Vars<T> rt(dis); for(auto& a:rt){ a = t/a; } return rt; }
template<class T>	Vars<T>   mag(const piece<var3<T>>& dis)  { Vars<T> rt(dis.size()); for_(dis,i){ rt[i] = mag(dis[i]); } return rt; }



template<class T> Vars<T>  operator &(const piece<var3<T>>& dis, const piece<var3<T>>& v) { Vars<T> rt(dis.size()); for_i(rt){ rt[i] = dis[i]&v[i]; } return rt; }



//! same as Vars but with a default value to hold a uniform vector,
//! a name and transformation information, used in #svplot only for now
template<class T>
class varsORv: public Vars<T> {
 public:
	T dfult;
	std::string name;
	T Xa,  Xp; //! scale, shift
	using piece<T>::d;
	using piece<T>::dn;
	T (*transf)(T);
	T ax,  px; ///  scale, shift before transform

	varsORv()                        :                      dfult(0)         , Xa(1)     , Xp(0)     , transf([](double a){ return a; }), ax(1)     , px(0)     {}
	varsORv(const T& val)            :                      dfult(val)       , Xa(1)     , Xp(0)     , transf([](double a){ return a; }), ax(1)     , px(0)     {}
	varsORv(size_t siz, const T& val): Vars<T>(siz,val), dfult(val)       , Xa(1)     , Xp(0)     , transf([](double a){ return a; }), ax(1)     , px(0)     {}
	varsORv(const Vars<T>& vs)    : Vars<T>(vs),      dfult(vs.back()) , Xa(1)     , Xp(0)     , transf([](double a){ return a; }), ax(1)     , px(0)     {}
	varsORv(const varsORv& vs)       : Vars<T>(vs),      dfult(vs.dfult)  , Xa(vs.Xa) , Xp(vs.Xp) , transf(vs.transf)              , ax(vs.ax) , px(vs.px) {}
	varsORv(const std::vector<T>& vs): Vars<T>(vs),      dfult(vs.dfult)  , Xa(vs.Xa) , Xp(vs.Xp) , transf(vs.transf)              , ax(vs.ax) , px(vs.px) {}
	varsORv(varsORv&& vs)  noexcept  : Vars<T>(vs),      dfult(vs.dfult)  , Xa(vs.Xa) , Xp(vs.Xp) , transf(vs.transf)              , ax(vs.ax) , px(vs.px) {} ;
	varsORv(const T* dd, int nn)     : Vars<T>(dd,nn),   dfult(dd[nn-1])  , Xa(1)     , Xp(0)     , transf([](double a){ return a; }), ax(1)     , px(0)     {}
	varsORv(const T* dd, const T* de): Vars<T>(dd,de),   dfult(*(de-1))   , Xa(1)     , Xp(0)     , transf([](double a){ return a; }), ax(1)     , px(0)     {}
	varsORv(const T* dd, const T* de, T(*fn)(T)): Vars<T>(dd,de,fn)       , Xa(1)     , Xp(0)     , transf([](double a){ return a; }), ax(1)     , px(0)     {}
	~varsORv(){};

	//varsORv(Vars<T>& vs, bool move): Vars<T>(vs,move)  {} //same as above but enforced by move
	void operator =(const piece<T>& v)  { { if(d) delete[] d; }  d = new T[v.size()]; std::copy(v.d, v.dn, d);  dn=d+v.size(); };
	void operator =(const Vars<T>& v){ { if(d) delete[] d; }  d = new T[v.size()]; std::copy(v.d, v.dn, d);  dn=d+v.size(); };
	void operator =(const varsORv& v)   { { if(d) delete[] d; }  d = new T[v.size()]; std::copy(v.d, v.dn, d);  dn=d+v.size(); dfult=v.dfult; Xa=v.Xa; Xp=v.Xp; transf=v.transf; ax=v.ax; px=v.px; };
	void operator =(const std::vector<T>& v){ { if(d) delete[] d; }  d = new T[v.size()]; std::copy(&v[0], &*v.end(), d);  dn=d+v.size(); };
	void operator =(const std::initializer_list<T>& v){ { if(d) delete[] d; }  d = new T[v.size()]; std::copy(v.begin(), v.end(), d);  dn=d+v.size(); };

	void operator =(Vars<T>&& v){ Vars<T>::eat(v); };

	T& operator [](int i)             { if(d+i<dn) return d[i]; else return dfult; };
	const T& operator [](int i) const { if(d+i<dn) return d[i]; else return dfult; };
	//T& operator ()(int i)             { return   (d+i<dn ? Xp+Xa*transf(d[i]+px) : dfult); };
	//const T& operator ()(int i) const { return   (d+i<dn ? Xp+Xa*transf(d[i]+px) : dfult); };
	T scalefrom01(T val) const { return   Xp+Xa*transf(val*ax+px); };

	void rescaledata(T d0, T del)	                     { Xp=d0; Xa=del; T* dd=d-1;  while(++dd<dn)  *dd = (*dd-d0)/del;  }
	void rescaledata(T d0, T del, T(*fn)(T))             { Xp=d0; Xa=del; T* dd=d-1;  while(++dd<dn)  *dd = fn((*dd-Xp)/Xa);  }
	void rescaledata(T d0, T del, T(*fn)(T), T p0, T a0) { Xp=d0; Xa=del; px=fn((p0-Xp)/Xa); ax=fn((p0+a0-Xp)/Xa)-px;   T* dd=d-1;
	                                                       while(++dd<dn)  *dd = (fn((*dd-Xp)/Xa)-px)/ax;  }
	//void transformdata(T d0, T del, T(*fn)(T), T p0)	{
		//Xp=d0; Xa=del; px=fn((p0-d0)/del); T* dd=d-1;  while(++dd<dn)	
			//*dd = fn((*dd-d0)/del)-px; }
};

//template<class T>
	//T rescaleval(T val, T d0, T del, T(*fn)(T), T p0, T a0) { T px=fn((p0-d0)/del); T ax=fn((p0+a0-d0)/del)-px;  return (fn((val-d0)/del)-px)/ax; }
inline double rescaleval(double val, double d0, double del, double(*fn)(double), double p0, double a0) { double px=fn((p0-d0)/del); double ax=fn((p0+a0-d0)/del)-px;  return (fn((val-d0)/del)-px)/ax; }

template<class T>  Vars<T> operator -(T v, Vars<T> vs)  { Vars<T> rt(vs); for(auto& a:rt){ a =v-a; };  return rt; }

template<class T>  var3<T> round(const var3<T>& vec) 	{ 	var3<T> rt(vec); 	rt.x=round(vec.x); rt.y=round(vec.y); rt.z=round(vec.z); 	return rt; 	}
template<class T>  Vars<T> round(const Vars<T>& vec) 	{ 	Vars<T> rt(vec); 	for(auto& v:rt) v=round(v); 	return rt; 	}

//! moves iss into oss then deletes iss, but efficiently
//template<class T> void transferto(Vars<Vars<T> >& oss, size_t i, Vars<Vars<T> >&& iss, size_t j=0)
//{	auto* os = oss.d+i-1;	auto* is = iss.d+j-1;	while ((++is)<iss.dn) { (++os)->d=is->d; os->dn=is->dn; is->d=0;	}  /*iss.d=0;*/	}



#endif // TYPSES_H





#ifndef TYPSES_OERATIONS_H
#define TYPSES_OERATIONS_H



template<class C> int len(C vec) { return vec.size(); } // casts size to int

template<class T> double sumdbl(const piece<T>& ps)  { double sm=1e-300; for(auto a:ps){ sm += a; }  return sm; }
template<class T> double sumdbl(const piece<T>& ps, const piece<T>& ws)  { double sm=1e-300; const T* p = ps.d-1, *w=ws.d-1; while(++p<ps.dn){ sm += *w * *(++p); }  return sm; }
inline            double sumdbl(const dbl3& ps, const dbl3& ws)  { return ps.x+ws.x + ps.y+ws.y + ps.z+ws.z; }
inline            double sumdbl(const dbl3& ps)                  { return ps.x+ ps.y+ ps.z; }

template<class T> T sumq(piece<T> ps)  {	T sm = *ps.d; while(++ps.d<ps.dn){ sm += *ps.d; }  return sm; }//! quick sum, assumes non-empty vec

template<class T>          T sum(const piece<T>& ps)  {	T sm = (T()*0.); const T* p= ps.d-1; while(++p<ps.dn){ sm+= *p; }  return sm; }
template<class T, class U> T sum(const piece<T>& ps, const piece<U>& ws)  {  T sm= (T()*0.);  T* p= ps.d-1;  U *w=ws.d-1; while(++p<ps.dn){ sm+= *(++w) * *p; }  return sm; }
template<class T> vars<T> sum(piece<piece<T>> ps)  { Vars<T> sm(ps[0]); while(++ps.d<ps.dn){ sm+= *ps.d; }  return sm;  } // does not work for size zero
	
template<class T>           T sumSqrs(const piece<T>& ps)  {  T sm= (T()*0.); const T* p= ps.d-1; while(++p<ps.dn){ sm+= *p * *p; }  return sm; }
template<class T, class U>  T sumSqrs(const piece<T>& ps, const piece<U>& ws)  { //TODO remove refs
	T sm = (T()*0.); const T* p = ps.d-1; const U *w=ws.d-1; while(++p<ps.dn){ sm += *(++w) * *p * *p; }  return sm; }
template<class T>   double sumdblSqrs(const piece<T>& ps)  {  double sm=0.; const T* p= ps.d-1; while(++p<ps.dn){ sm+= double(*p) * *p; }  return sm; }


template<class T>	T avg(const piece<T>& dis)  { return sum(dis)/dis.size(); }
template<class T>   T min(const piece<T>& pis)  { return *std::min_element(pis.d,pis.dn); }
template<class T>   T max(const piece<T>& pis)  { return *std::max_element(pis.d,pis.dn); }
template<class T, class TFunc> void transform(piece<T>& pis, TFunc func)  { for(T& dr:pis) dr=func(dr); }
template<class T>              void transform(T* d, const T* dn, T(*func)(const T&))  { --d; while (++d<dn) *d=func(*d); }
//template<class T> void transform(piece<T>& pis, float(*func)(float&))  { for(float& dr:pis) dr=func(dr); }

//double (*func)(const double&) = ([](const double& d) { return d; }) // this is just  a note/reminder on old  function / lambda...


template<class T> std::vector<T>& operator *= (std::vector<T>& vs, double scale)  {  for(T& v: vs) v*=scale;   return vs;  }



template<class T, template<class ...> class C1, template<class ...> class C2>
vars<vars<T> > transpose(const C1<C2<T> >& vecvec)  {
	if(!vecvec.size()) return vars<vars<T> >();
	vars<vars<T> > trans(vecvec[0].size(),vars<T>(vecvec.size()));
	for (size_t i=0; i<vecvec[0].size(); ++i)
		for (size_t j=0; j<vecvec.size(); ++j) trans[i][j] = vecvec[j][i] ;
	return trans;
}

template<class T, size_t N, template<class ...> class C1>
vars<vars<T> > transpose(const C1<std::array<T,N> >& vecvec)  {
	if(!vecvec.size()) return vars<vars<T> >();
	vars<vars<T> > trans(vecvec[0].size(),vars<T>(vecvec.size()));
	for (size_t i=0; i<vecvec[0].size(); ++i)
		for (size_t j=0; j<vecvec.size(); ++j) trans[i][j] = vecvec[j][i] ;
	return trans;
}

template<class T, template<class ...> class C>  piece<T> allof(C<T>& vs) { return piece<T>(&*vs.begin(), &*vs.end()); }




//! class table IO, usage: cout<< tableIO(stdVecVec)  (not TableIO, for auto detecting template paarams)
template<class T, template<class ...> class C1, template<class ...> class C2>
class TableIO  { public:

	TableIO(const C1<C2<T> >& tbl, std::vector<std::string> hdrs=std::vector<std::string>(), bool transpose=true, char sepr='\t') 
	:	vss_(tbl),  hds_(hdrs),  transpose_(transpose), sep_(sepr) {}; // careful: tbl lifetime: don't pass rvalue tbl

	std::string hdr(size_t i) const { return hds_.size()>i ? hds_[i] : ""; }

	const C1<C2<T> >&        vss_; // ref is bad idea
	std::vector<std::string> hds_;
	const char               transpose_;
	const char               sep_;
};


//! class table IO, usage:    cout<< pieces().apnd(vec,"nam") ;
//template<class T>
//class piecesIO : public TableIO<T, std::vector, piece > { public:
	//using TableIO<T, std::vector, piece >::hds_;
	//piecesIO(bool transpose=true, char sepr='\t') :  TableIO<T,std::vector, piece>(vss_,std::vector<std::string>(),transpose,sepr) {};
	//piecesIO& apnd(piece<T> vs, std::string nam="")  { vss_.push_back(vs); if(hds_.size()||nam.size()) {hds_.push_back(nam); }  return *this; }
	//std::vector<piece<T> >   vss_;
//};

//! class Table for IO, usage: cout<< Table().apnd(vec(...),"nam") ;
template<class T, template<class ...> class C2>
class Table : public TableIO<T, std::vector, C2 > { public:
	using TableIO<T, std::vector, C2>::hds_;
	Table(bool transpose=true, char sepr='\t') :  TableIO<T,std::vector, C2>(vss_,std::vector<std::string>(),transpose,sepr) {};
	Table& apnd(const C2<T>& vs, std::string nam="")  { vss_.push_back(vs); if(hds_.size()||nam.size()) {hds_.push_back(nam); }  return *this; }
	Table& apnd(const TableIO<T, std::vector, C2>& tbl) { 
		for_i(tbl.vss_) { vss_.push_back(tbl.vss_[i]); if(hds_.size()&&tbl.hds_.size()) hds_.push_back(tbl.hds_[i]); }  return *this; }
	std::vector<C2<T> >   vss_;
};

template<class T, template<class ...> class C1, template<class ...> class C2>
TableIO<T,C1,C2> tableIO(const C1<C2<T> >& vecvec, std::vector<std::string> hdrs=std::vector<std::string>(), bool transpose=true, char sepr='\t') {
	return TableIO<T,C1,C2>(vecvec,hdrs,transpose,sepr);  }

template<class T, template<class ...> class C1, template<class ...> class C2>
                          std::ostream& operator << (std::ostream& out, const TableIO<T,C1,C2>& tbl) {
	if(tbl.hds_.size()==tbl.vss_.size()) {
		for_(tbl.hds_,j) { out<<std::setw(12)<<std::left<< tbl.hds_[j]   <<tbl.sep_; }  out<<'\n'; }
	for_(tbl.vss_[0],i) {
		for_(tbl.vss_,j) { out<<std::setw(12)<<std::left<< tbl.vss_[j][i]<<tbl.sep_; }  out<<'\n'; }
	out << '\n';
	return out;
}

template<class T, template<class ...> class C2>
                          std::istream& operator >> (std::istream& ins, Table<T,C2>& tbl) { // inefficient!!!
	ins >> std::ws;
	std::string row, item;
	if( std::is_arithmetic<T>::value && isalpha(ins.peek()))  {
		getline( ins, row );
		std::stringstream ss( row );
		while ( getline( ss, item, tbl.sep_ ) )   tbl.hds_.push_back(item);
	}
	std::vector<std::vector<T> > res;
	while( getline( ins, row ))  {
		std::vector<T> Ro;	std::stringstream ss( row );
		while ( getline( ss, item, tbl.sep_ ) )   Ro.push_back( strTo<T>(item) );
		res.push_back( Ro );
	}
	if(res.size()) tbl.vss_.resize(res[0].size());   else return ins;
	for (size_t i=0; i<tbl.vss_.size(); ++i)  {
		tbl.vss_[i].resize(res.size());
		if(res[i].size()!=tbl.vss_.size()) std::cout<<"Error table sizes don't match"<<std::endl;
		for (size_t j=0; j<tbl.vss_[i].size(); ++j) tbl.vss_[i][j]=res[j][i];
	}
	return ins;
}




template<class T>          std::ostream& operator << (std::ostream& out, const var3<T>& pos) {
	out << pos.x<<' '<< pos.y<<' '<<pos.z;   return out;  }

template<class T>          std::istream& operator >> (std::istream& ins, var3<T>& pos) {
	ins >> pos.x >> pos.y >> pos.z;	return ins;	}

template<class T, class U> std::istream& operator >> (std::istream& ins, std::pair<T,U>& ab) {
 	ins >> ab.first >> ab.second;   return ins; 	}

template<class T, class U> std::ostream& operator << (std::ostream& out, const std::pair<T,U>& ab) {
 	out << ab.first<<' '<<ab.second; return out; 	}

template<class T>          std::istream& operator >> (std::istream& ins, var2<T>& ab) {
 	ins >> ab.a >> ab.b;     return ins; 	}

template<class T>          std::ostream& operator << (std::ostream& out, const var2<T>& ab) {
 	out << ab.a<<" "<< ab.b; return out; 	}

template<class T>          std::ostream& operator << (std::ostream& out, const std::vector<T>& vec) {
	if(vec.size()<16)  for (auto v : vec) out << v << '\t';
	else               for (auto v : vec) out << v << '\n';
	return out;
}

template<class T>          std::istream& operator >> (std::istream& ins, std::vector<T>& vec) {
	if(vec.size()) 
		for (size_t i=0; i<vec.size(); ++i)  ins >> vec[i];
	else { T t; while (ins>>t) vec.push_back(t); }
	return ins;
}

template<class T, class U> std::ostream& operator << (std::ostream& out, const std::map<T,U>& vec) {
 	for (auto v : vec) { out << v << '\n'; }  return out; 	}


template<class T>          std::ostream& operator << (std::ostream& out, const piece<T>& vec)  {
	if(vec.size()<16 ) for (auto v : vec) out << v << '\t';
	else               for (auto v : vec) out << v << '\n';
	return out;
}


template<class T>          std::istream& operator >> (std::istream& ins, piece<T> vec) {// Warning relying on copy constructor not to deep copy
	for (; vec.d<vec.end(); ++vec.d)   ins >> *vec.d;
	//if (ins.fail()) std::cout<<"Error while reading array"<<std::endl;
	return ins;
}

template<class T, size_t N> std::istream& operator >> (std::istream& ins, std::array<T,N>& vec) {
	for(auto& v:vec)      { ins >> v; }   return ins;   }
template<class T, size_t N> std::ostream& operator << (std::ostream& out, const std::array<T,N>& vec)  {
	for (const auto& v: vec) { out << v << '\t'; } return out;  }




inline void        replaceInFromTo(std::string& str, const std::string& frm, const std::string& to) {
	//return str = std::regex_replace( str, std::regex(frm), to );
	size_t po = 0;
	while((po = str.find(frm, po)) != std::string::npos) {
		str.replace(po, frm.length(), to);  po += to.length(); }
}

inline void        replaceInBetweenTo(std::string& str, const std::string& bgn, const std::string& end, const std::string & to) {
	//no greedy macth
	size_t po = 0, po0=0;
	while((po = str.find(bgn, po)) != std::string::npos) { po0=po; po=po+bgn.length();
		while((po = str.find(end, po)) != std::string::npos) {
			str.replace(po0, po-po0, to); po = po0+to.length();  break; } }
}

inline std::string replaceFromTo(const std::string& str, const std::string& frm, const std::string& to) {
	return std::regex_replace(str, std::regex(frm), to);
	//size_t po = 0;
	//while((po = str.find(frm, po)) != std::string::npos) { str.replace(po, frm.length(), to);  po += to.length();   }
	//return str;
}

inline std::string baseName(const std::string& fName) { // removes file suffix
	auto dloc=fName.find_last_of('.'); 	
    if( dloc != std::string::npos)    {
    	auto sl=fName.find_last_of('/');  
    	if( sl!=std::string::npos && sl>dloc) return fName; 
		return fName.substr(0,dloc);
	}
	return fName;
}




template<class T> bool _1At(T n, int ibgn)  {  return n&(1<<ibgn);  }
inline            bool _1At(unsigned int n, int ibgn)  {  return n&(1<<ibgn);  }// for debugger otherwise redundant
template<class T> bool _1In(T n, int ibgn, int lnt)  {  return (n>>ibgn)&((1<<lnt)-1);  }





// CLEAN UP PLEASE:

inline double roundec(double x,int d)  { double scale=std::pow(10,int(std::log10(x)))/std::pow(10,d); return std::round(x/scale)*scale; } // round x to d significant digits

//template<class T> piece<T>   insidvars(piece<T>& vs) { return piece<T>(vs.d+1, vs.dn-1); }
template<class T> Vars<T>  diffVars(piece<T> vs) { Vars<T> rt(vs);  piece<T>(rt.d+1,rt.dn)-=vs; rt[0]=rt[1]; return rt; }
template<class T> Vars<T>  movingAvg(piece<T> vs) { // apply rolling 3 moving average X <- ( X-- + X + X++)/3
	Vars<T> rt(vs.size());  T* d=rt.d; *d=0.5*(vs[0]+vs[1]); 
	for(T* o=vs.d+2; o<vs.dn;++o) { *(++d)=(1./3.)*(*(o-2)+*(o-1)+*o); } 
	*(++d)=0.5*(*(vs.dn-2)+*(vs.dn-1)); return rt; } // TODO: use vs storage instead of o

template<class T> T closer(T cv, T v1, T v2) { return (abs(cv-v1)<abs(cv-v2)) ?  v1 : v2; }
//template<class T> Vars<T>  biMovingAvg(piece<T> vs) { // apply rolling 3 moving average X <- ( X-- + X + X++)/3
	//if(len(vs)<7) return vs; constexpr double p3=(1./3.);	Vars<T> res(vs.size());  T* rt=res.d; 

	//*rt=0.5*(vs[0]+vs[1]); *++rt=p3*(vs[0]+vs[1]+vs[2]);  *++rt=p3*(vs[1]+vs[2]+vs[3]); 
	//vs.dn-=3; T* o=vs.d+3;
	//for(; o<vs.dn;++o) { *++rt = closer(*o,0.25*(*(o-3)+*(o-2)+*(o-1)+*o), 0.25*(*(o)+*(o+1)+*(o+2)+*(o+3)));  }  --o;
	//*++rt=p3*(o[0]+o[1]+o[2]); *++rt=p3*(o[1]+o[2]+o[3]); *++rt=0.5*(o[2]+o[3]);  
	//return res; 
//}


// TODO test BC handling
template<class T> Vars<T>  biMovingAvg(piece<T> vs,int krnl, int bKrnl=1) { // apply rolling 3 moving average X <- ( X-- + X + X++)/3
	if(len(vs)<7) return vs; 
	Vars<T> res(vs.size());  T* rt=res.d-1;  const double pl=1./(1+krnl);	  bKrnl=std::min(bKrnl,krnl);
	for(int i=0; i<krnl; ++i)  { int krp = std::max(i,bKrnl)+1; *++rt=sumq(piece<T>(vs.d,krp)/krp); }
	vs.dn-=krnl; T* o=vs.d+krnl;
	for(; o<vs.dn;++o) { *++rt=closer(*o, sumq(piece<T>(o-krnl,1+krnl))*pl, sumq(piece<T>(o,1+krnl))*pl); }  --o; vs.dn+=krnl;
	for(int i=krnl; i>0; --i) { int krp = std::max(i,bKrnl)+1;  *++rt=sumq(piece<T>(vs.dn-krp, krp))/krp; }
	return res; 
}

// TODO test BC handling
template<class T> void  NaNsToMean(piece<T> vs) {
	double sum=0.; size_t count=0;
	for(auto v:vs)  if(v==v) { sum+=v; ++count; } //not NaNs
	T mean = sum/count;
	for(auto& v:vs)  if(!(v==v)) { v=mean; } //NaNs
}


template<class T> T med_closer(vars<T> vs, T orig)  {
	if (vs.size()&1==0) {
		 const auto nth = vs.begin() + vs.size()/2 - 1;
		 std::nth_element(vs.begin(),   nth , vs.end());   const auto e1 = *nth;
		 std::nth_element(vs.begin(), ++nth , vs.end());   const auto e2 = *nth;
		 return closer(orig, e1, e2);
	} else { const auto nth = vs.begin()+vs.size()/2;   std::nth_element(vs.begin(), nth, vs.end());    return *nth;  }
}
template<class T> T med_odd(vars<T> vs)  { const auto nth = vs.begin()+vs.size()/2;  std::nth_element(vs.begin(), nth, vs.end());  return *nth;  }
template<class T> Vars<T>  median(piece<T> vs,int krnl, int bKrnl=1) { // apply rolling 3 moving average X <- ( X-- + X + X++)/3
	if(len(vs)<7) return vs;	
	Vars<T> res(vs.size());  T* rt=res.d-1;	    bKrnl=std::min(bKrnl,krnl);
	for(int i=0; i<krnl; ++i)   *++rt=med_odd(vars<T>(vs.d,2*std::max(i,bKrnl)+1));
	vs.dn-=krnl; T* o=vs.d+krnl;
	for(; o<vs.dn;++o) { *++rt=med_odd(vars<T>(o-krnl,1+2*krnl)); }  --o; vs.dn+=krnl;  // vars -> piece ?
	for(int i=krnl; i>0; --i) { int krn = std::max(i,bKrnl);  *++rt=med_odd(vars<T>(vs.dn-2*krn-1, 2*krn+1)); }
	return res; 
}





#ifdef VMMLIB__VECTOR__HPP
 template< size_t M, typename T >
 Vctr<M,T>::Vctr(const var3<T>& v3) 	{ array[ 0 ] = v3.x;	array[ 1 ] = v3.y;	array[ 2 ] = v3.z; 	}
#endif







template<class T>
vars<dbls> distribution(const piece<T>  xs, const piece<T> ws, int nBins=64)  {
	vars<dbls>  distrib(3, dbls(nBins,0.));
	double minU=min(xs),   maxU=max(xs),   deltaU=(maxU-minU)/nBins+1e-72;

	for (int i=0; i<nBins; ++i)	distrib[0][i] = minU+deltaU/2+i*deltaU;

	for (size_t i=0; i<ws.size(); ++i)  {
		int distInd=std::min(int((xs[i]-minU)/deltaU+0.5),nBins-1);
		++distrib[1][distInd];
		distrib[2][distInd]+=ws[i];
	}

	distrib[1]/=distrib[1].sum()*deltaU;
	distrib[2]/=distrib[2].sum()*deltaU;

	return distrib;
}




inline double linearInterpolate(double x, double x1, double x2, double y1, double y2) 	 { 	 return y1+(x-x1)*(y2-y1)/(x2-x1); 	 }

inline double averageCDF(double xMin, double xMax, const std::vector< std::pair<double, double> >& tabl)  {
	/// used to get average of a part of a distribution
	auto itr = tabl.begin();
	while((++itr)->first <=  xMin && itr != tabl.end()) ;

	if( !(itr->first >= xMin && (itr-1)->first <=  xMin))  {
		std::cout<<"\n LookUp Error:  "	<< xMin-1.<<"     " << tabl.begin()->first-1.<<"  "  << tabl.begin()->second<<"    " << tabl.rbegin()->first-1.<<"  " << tabl.rbegin()->second
			 <<"\n\n "<<(itr == tabl.end())<<std::endl; exit(-1);  }

	assert(itr->first > xMin && (itr-1)->first <=  xMin);
	
	double wSum=itr->first-xMin;
	double wxy1=(itr->first-xMin)*0.5*(itr->second+linearInterpolate(xMin,(itr-1)->first, itr->first, (itr-1)->second, itr->second));
	
	while((++itr)->first <=  xMax && itr != tabl.end())  {
		wSum+=itr->first-(itr-1)->first;
		wxy1+=(itr->first-(itr-1)->first)*0.5*(itr->second+(itr-1)->second);
		//cout<<itr->first<<": "<<itr->second<<"     "<<0.5*(itr->second+(itr-1)->second)<<" s/ "<<(itr->first-(itr-1)->first)<<"   =  "<<wxy1/wSum<<std::endl;
	}
	if(itr != tabl.end())  {
		wSum+=itr->first-xMax;
		wxy1+=(itr->first-xMax)*0.5*(itr->second+linearInterpolate(xMin,(itr-1)->first, itr->first, (itr-1)->second, itr->second));
	}


	return wxy1/(wSum+1e-64);

}


//double variance(const piece<double>& vs, piece<double> ws)  {  return sumSqrs(vs-vs.avg(), ws)/sum(ws);  } // biased
template<class T, class T2>
double corelCoeff(piece<T> X, piece<T> Y, piece<T2> ws) { return sum((X-X.avg())*(Y-Y.avg())*ws) / sqrt(sumSqrs(X-X.avg(), ws) * sumSqrs(Y-Y.avg(), ws)); }
//template<class T>
//double covarianceDbl(piece<T> X, piece<T> Y) { return sumdbl((X-X.avgdbl())*(Y-Y.avgdbl())) / (X.size()-1); }
//template<class T>
//double   varianceDbl(piece<T> X) { return sumdblSqrs(X-X.avgdbl()) / (X.size()-1); }


inline double sqr(double a) { return a*a; }

#endif //TYPSES_OERATIONS_H




/* //- Debugging:
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
    libstdcxx_printer.add_container('', 'Vars', varsPrinter)
*/


