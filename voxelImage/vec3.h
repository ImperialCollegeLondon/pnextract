#ifndef VEC3_H
#define VEC3_H

#include <iomanip>
#include <sstream>
#include <fstream>
#include <array>
#include <cmath>
#include <cstring>
#include <algorithm>

#ifndef verySmall
	#define verySmall  1.0e-31
#endif


//typedef  int Int;
typedef  long long Long;
typedef  std::array<int,3> int3;
//typedef  std::array<Int,3> int3;
typedef  std::array<int3,3> int3x3;


class vec3
{
public:

	double	x;
	double	y;
	double	z;


	vec3() {}
	vec3(double r, double s, double t)  { x = r;  y = s;  z = t; }
	//vec3(int3 n)  { x = n[0];  y = n[1];  z = n[2]; }
	vec3(int3 n)  { x = n[0];  y = n[1];  z = n[2]; }
	vec3(double* d)  { x = d[0];  y = d[1];  z = d[2]; }

	vec3& Set(double r, double s, double t)  { x = r;  y = s;  z = t;  return (*this); }

	double&       operator [](long k) {	return ((&x)[k]); }

	const double& operator [](long k) const  { return ((&x)[k]); }

	vec3& operator +=(const vec3& v)  { x += v.x;  y += v.y;  z += v.z;  return (*this); }
	vec3& operator -=(const vec3& v)  { x -= v.x;  y -= v.y;  z -= v.z;  return (*this); }
	vec3& operator *=(double t)  { x *= t;  y *= t;  z *= t;  return (*this); }
	vec3& operator /=(double t)  {  double f = 1.0 / t;  x *= f;  y *= f;  z *= f;  return (*this); }
	vec3& operator ^=(const vec3& v)  { double r, s;  r=y*v.z-z*v.y;  s=z*v.x-x*v.z;  z=x*v.y-y*v.x;  x=r; y=s; 	return (*this); }
	vec3& operator *=(const vec3& v)  { x *= v.x;  y *= v.y;  z *= v.z;  return (*this); }
	vec3  operator -(void) const  { return (vec3(-x, -y, -z)); }
	vec3  operator +(const vec3& v) const  { return (vec3(x+v.x, y+v.y, z+v.z)); }
	vec3  operator -(const vec3& v) const  { return (vec3(x-v.x, y-v.y, z-v.z)); }
	vec3  operator *(double t) const  { return (vec3(x*t, y*t, z*t)); }
	vec3  operator /(double t) const  { double f = 1.0 / t;  return (vec3(x*f, y*f, z*f)); }
	double operator &(const vec3& v) const  { return (x*v.x+y*v.y+z*v.z); }
	vec3  operator ^(const vec3& v) const  { return (vec3(y*v.z-z*v.y,  z*v.x-x*v.z,  x*v.y-y*v.x)); }
	vec3  operator *(const vec3& v) const  { return (vec3(x*v.x, y*v.y, z*v.z)); }
	vec3  operator /(const vec3& v) const  { return (vec3(x/v.x, y/v.y, z/v.z)); }
	//~ bool operator ==(const vec3& v) const  { return ((x == v.x) && (y == v.y) && (z == v.z)); }
	//~ bool operator !=(const vec3& v) const  { return ((x != v.x) || (y != v.y) || (z != v.z)); }
	bool operator ==(const vec3& v) const  { return ((x-v.x)*(x-v.x) < verySmall) && ((y-v.y)*(y-v.y) < verySmall) && ((z-v.z)*(z-v.z) < verySmall); }
	bool operator !=(const vec3& v) const  { return ((x-v.x)*(x-v.x) >= verySmall) || ((y-v.y)*(y-v.y) >= verySmall) || ((z-v.z)*(z-v.z) >= verySmall); }

};

inline vec3 rotateAroundLine(vec3 y, double gamma,  vec3 n, vec3 x)
{///. rotate y around line passing through x, in the direction of n, http://inside.mines.edu/~gmurray/ArbitraryAxisRotation
	double s = sinf(gamma),   c = cosf(gamma);
	double k = 1.0 - c;
	return vec3(
	 	( x.x*(n.y*n.y+n.z*n.z) - n.x*( x.y*n.y+x.z*n.z-n.x*y.x- n.y*y.y-n.z*y.z ) )*k + y.x*c + (-x.z*n.y+x.y*n.z-n.z*y.y+n.y*y.z )*s,
		( x.y*(n.x*n.x+n.z*n.z) - n.y*( x.x*n.x+x.z*n.z-n.x*y.x- n.y*y.y-n.z*y.z ) )*k + y.y*c + ( x.z*n.x-x.x*n.z+n.z*y.x-n.x*y.z )*s,
		( x.z*(n.x*n.x+n.y*n.y) - n.z*( x.x*n.x+x.y*n.y-n.x*y.x- n.y*y.y-n.z*y.z ) )*k + y.z*c + (-x.y*n.x+x.x*n.y-n.y*y.x+n.x*y.y )*s );
}
inline vec3 rotateAroundVec(const vec3 y, double gamma, vec3 n)
{///. rotate y around n (line passing through centre, in the direction of n) http://inside.mines.edu/~gmurray/ArbitraryAxisRotation
	double s = sinf(gamma),   c = cosf(gamma);
	double k = 1.0 - c;
	return vec3(
		(  - n.x*( -n.x*y.x- n.y*y.y-n.z*y.z ) )*k + y.x*c + (n.y*y.z-n.z*y.y)*s,
		(  - n.y*( -n.x*y.x- n.y*y.y-n.z*y.z ) )*k + y.y*c + (n.z*y.x-n.x*y.z)*s,
		(  - n.z*( -n.x*y.x- n.y*y.y-n.z*y.z ) )*k + y.z*c + (n.x*y.y-n.y*y.x)*s );
}



inline vec3 operator*(double t, const vec3& v) { return vec3(t*v.x, t*v.y, t*v.z); }

inline double mag(const vec3& v) { return std::sqrt(v.x*v.x+v.y*v.y+v.z*v.z); }

inline double magSqr(const vec3& v) { return (v.x*v.x+v.y*v.y+v.z*v.z); }
inline vec3 norm(const vec3& v)  { return  v/(mag(v)+1.0e-300); }



inline vec3  operator*( int3 n, vec3 x) { return vec3(n[0]*x[0], n[1]*x[1], n[2]*x[2]);}
inline int3 operator-( int3 n, int3 x) { return int3{{n[0]-x[0], n[1]-x[1], n[2]-x[2]}};}
inline int3 operator*( int s, int3 n) { return int3{{s*n[0], s*n[1], s*n[2]}};}
inline int3 operator*( double s, int3 n) { return int3{{int(s*n[0]), int(s*n[1]), int(s*n[2])}};}
inline int3 operator/( int3 n, int s) { return int3{{n[0]/s, n[1]/s, n[2]/s}};}
//inline int3& operator+=( int3& n, int3 x) { n[0]+=x[0]; n[1]+=x[1]; n[2]+=x[2]; return n;}
inline int3& operator+=( int3& n, int3 x) { n[0]+=x[0]; n[1]+=x[1]; n[2]+=x[2]; return n;}






inline std::ostream& operator<< (std::ostream& out, const vec3& node)
{
	std::ios_base::fmtflags flgs=out.flags();
	out.flags(std::ios::showpoint | std::ios::scientific);	
	out << std::setprecision(5)  << node.x   <<" " << node.y   <<" " << node.z;
	out.flags(flgs);
	return out;
}
inline std::ofstream& operator<< (std::ofstream& out, const vec3& node)
{
	std::ios_base::fmtflags flgs=out.flags();
	out.flags(std::ios::showpoint | std::ios::scientific);	
	out << std::setprecision(8)  << node.x   <<" " << node.y   <<" " << node.z;
	out.flags(flgs);
	return out;
}

inline std::ostream& operator<< (std::ostream& out, const int3& ijk)
{
    out << ijk[0] <<" "<< ijk[1]<<" " << ijk[2];
	return out;
}
inline std::istream& operator>> (std::istream& in, vec3& node)
{
	in >> node.x >> node.y >> node.z;
	return in;
}
inline std::istream& operator>> (std::istream& in, int3& ijk)
{
    in >> ijk[0] >> ijk[1] >> ijk[2];
    return in;
}


#if cplusplus > 201103L
using dblpair = std::pair<double,double>;
#else
#define dblpair  std::pair<double,double>
#endif

template<typename T1, typename T2>
inline std::istream& operator>> (std::istream& in, std::pair<T1,T2>& interval)
{
	in >> interval.first >> interval.second;
	return in;
}

template<typename T1, typename T2>
inline std::ostream& operator<< (std::ostream& out, std::pair<T1,T2>& interval)
{
	out  << interval.first<<" "<< interval.second;
	return out;
}

template<typename Type>
inline std::vector<Type>& operator*= (std::vector<Type>& vs, double scale)
{
	for(Type& v : vs) v*=scale;
	return vs;
}

template<typename T>
std::ostream & operator << (std::ostream & outstream, const std::vector<T>& vec)
{
	if(vec.size() && vec.size()<10)  for (auto v : vec) outstream << v << '\t';
	else                             for (auto v : vec) outstream << v << '\n';
	return outstream;
}



#ifdef VMMLIB__VECTOR__HPP
template< size_t M, typename T >
Vctr<M,T>::Vctr(const vec3& v3)
{
	array[ 0 ] = v3.x;
	array[ 1 ] = v3.y;
	array[ 2 ] = v3.z;
}
#endif


#ifndef TOSTR
#define TOSTR
template<typename T> std::string toStr(const T& n){  std::ostringstream stm;  stm<<n;  return stm.str();  }
#endif



#if __cplusplus >= 201103L
 template<class T> using  svec = std::vector<T>;
#else
 #define svec  std::vector
#endif



template<class T>
class piece
{
  protected:
	piece(): d(0), dn(0) {};
  public:

	piece(T* dd, int n): d(dd), dn(d+n) {};
	piece(const piece& p): d(p.d), dn(p.dn) {};//! note data hold by piece are not const

	T* begin() const {return d;};
	T* end()   const {return dn;};
	const T* cbegin() const {return d;};
	const T* cend()   const {return dn;};
	T* operator()()   const {return d;};
	//const T& operator[](int i) const {return d[i];};
	T& operator[](int i) const {return d[i];};
	//const T* operator()() const {return d;};
	size_t size() const {return dn-d;};




	piece& operator +=(const piece& v)  { for(auto& a:*this){ a += v[&a-d];};  return (*this); }
	piece& operator -=(const piece& v)  { for(auto& a:*this){ a -= v[&a-d];};  return (*this); }
	piece& operator *=(const piece& v)  { for(auto& a:*this){ a *= v[&a-d];};  return (*this); }
	piece& operator /=(const piece& v)  { for(auto& a:*this){ a /= v[&a-d];};  return (*this); }
	piece& operator *=(T t)             { for(auto& a:*this){ a *= t;};  return (*this); }
	piece& operator /=(T t)  { return (*this)*=(1.0/t); }
	T sum() const  { T sm=0; for(auto a:*this){ sm += a;}  return sm; }


  //protected:
	T*  d;
	T*  dn;
};
#define forvitr(vec,b,n) for(auto vi=piece<std::remove_reference<decltype(vec[0])>::type>(&vec[0]+b, n);vi()<vi.cend();++vi.d)

template<class T>
double sumdbl(piece<T> ps)  { double sm=1.0e-300; for(auto a:ps){ sm += a;}  return sm; }

template<class T>
class lazyvec: public piece<T>
{
	using piece<T>::d;
	using piece<T>::dn;

 public:

	lazyvec(): piece<T>(0,0) {};
	lazyvec(int siz): piece<T>(new T[siz],siz) {};  
	lazyvec(size_t siz, T val): piece<T>(new T[siz],siz) {  std::fill(d, dn, val);}
	lazyvec(const lazyvec& v): piece<T>(new T[v.size()],v.size())  {  std::memcpy(d, v.d, v.size()*sizeof(T)); }
	lazyvec(lazyvec&& v): piece<T>(v.d,v.size())  {  v.d=0; }
	lazyvec(const T* dd, int nn): piece<T>(new T[nn],nn)            { std::memcpy(d, dd, nn*sizeof(T)); }
	lazyvec(const T* dd, const T* de): piece<T>(new T[de-dd],de-dd)    { std::memcpy(d, dd, (de-dd)*sizeof(T));}
	~lazyvec(){ if(d) delete[] d; };

	void operator=(const lazyvec& v){ { if(d) delete[] d; }  d = new T[v.size()]; std::memcpy(d, v.d, v.size()*sizeof(T));  dn=d+v.size(); };

	lazyvec  operator -(void) const  { lazyvec res(*this);  for(auto& a:*this){ a = -a;};  return (*this);  }
	lazyvec  operator +(const piece<T>& v) const { return lazyvec(*this)+=v;  }
	lazyvec  operator -(const piece<T>& v) const { return lazyvec(*this)-=v;  }
	lazyvec  operator *(const piece<T>& v) const { return lazyvec(*this)*=v;  }
	lazyvec  operator /(const piece<T>& v) const { return lazyvec(*this)/=v;  }
	lazyvec  operator *(T t) const { return lazyvec(*this)*=t;  }
	lazyvec  operator /(T t) const { return lazyvec(*this)*=1.0/t;  }


   void resize(int nn)  
   { { if(d) delete[] d; }   if(nn) {d=new T[nn]; dn=d+nn;} else {d=0; dn=0;} }
   void pbak(T& vj)
   {	if(d) 
		{ T* od=d;  d=new T[dn+1-od];  
			std::memcpy(d, od, (dn-od)*sizeof(T));
			dn=dn+1-od; *(dn-1)=vj;  delete[] od;  }
		else    { d=new T[1];   *d=vj;   dn=d+1; }
	}
	

};
//template<typename T>
//std::string toStr(const piece<T>& vec, char spr='\t')
//{
	//std::ostringstream outstream; 
	//if(sizeof(T)<=8)  for (auto v : vec) outstream << v << spr;
	//else              for (auto v : vec) outstream << v << '\n';
	//return outstream.str();
//}
template<typename T>
std::ostream & operator << (std::ostream & outstream, const piece<T>& vec)
{
	if(vec.size() && vec.size()<10 )  for (auto v : vec) outstream << v << '\t';
	else                             for (auto v : vec) outstream << v << '\n';
	return outstream;
}





#endif


