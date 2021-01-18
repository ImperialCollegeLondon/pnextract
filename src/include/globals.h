#ifndef GLOBALS_SkipH
#define GLOBALS_SkipH
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <map>
// definition of global class and macros used to debugging and testing  

const static double PI = 3.14159265358979;



//using _s = std::to_string  is bad in decimal notation
#ifndef Str_SkipH
#define Str_SkipH
	template<typename T> std::string toStr(const T& n){  std::ostringstream ss;  ss<<n;  return ss.str();  } //outdated
	template<typename T> std::string _s(const T& n){  std::ostringstream ss;  ss<<n;  return ss.str();  }
template<typename T> T strTo(const std::string &s){  std::istringstream ss(s);  T t;  ss>>t;  return t; }
	//template<typename T> T fileExt(const T& path) { auto const pos = path.rfind(".");  return (pos==T::npos) ?  T{} : path.substr(pos); }
	template<typename T> bool hasExt(const T& path, size_t lnt, const char* ext) { return path.size()>lnt && path.compare(path.size()-lnt,lnt,ext) == 0; }

	//stringify macro args after (sub)macro expansion
	#define STRINGIFY(x) #x
	#define TOSTRING(x) STRINGIFY(x)
#endif // Str_SkipH



#ifndef Tst_SkipH
#define Tst_SkipH
	//- Testing macros

	inline bool _cerr_(std::string hdr="",std::string msg="", bool xit=false) // for debugger breakpoints: don't optimize out please !!!
		{	 if(xit) throw std::runtime_error(msg);  else std::cerr<< hdr+msg <<std::endl; return true; }

// Variable argument macro trick
	 #define ERR_HDR(isOk, _xit_) "\n  "+_xit_?"Error":"Warning" \
				" in "+ std::string(__FUNCTION__)+", "+std::string(__FILE__)+":"+_s(__LINE__) \
				, std::string(": { ")+std::string(#isOk)
	 #define ensure1(isOk)           (!((isOk)|| _cerr_(ERR_HDR(isOk,0)+" }")))
	 #define ensure2(isOk, msg)      (!((isOk)|| _cerr_(ERR_HDR(isOk,0)+"   '"+msg+"'  }")))
	 #define ensure3(isOk, msg, xit) (!((isOk)|| _cerr_(ERR_HDR(isOk,xit)+"   '"+msg+"'  }", xit)))
 #define GET_MACRO3(_1,_2,_3,NAME,...) NAME

//! Validation/production phase ensure/assert. Usage:
//!   \code{.cpp} ensure(condition, "message", throw_on_error=false); \endcode
 #define ensure(...)  GET_MACRO3(__VA_ARGS__, ensure3, ensure2,ensure1, "Only 1 to 3 args please")(__VA_ARGS__)
#endif // Tst_SkipH


#ifndef Dbg_SkipH
#define Dbg_SkipH
	class globals
	{
		static int    debugLevel_;
		static size_t typeCounter_;//dynamic type index, not used yet, see https://codereview.stackexchange.com/questions/44936/unique-type-id-in-c
	public:
		template<typename T>
		static size_t getTypeCounter()
		{
			static size_t id = typeCounter_++;
			return id;
		}
		static int getDebugLevel()  {  return debugLevel_;  }
		static void setDebugLevel(int dbgLvl)  {  debugLevel_=dbgLvl;  }

	};
	#if (__cplusplus >= 201402L)
	template<typename T> const auto  TypeID = globals::getTypeCounter<T>;
	#else
	#define  TypeID  globals::getTypeCounter<T>
	#endif

	#define debugLevel    globals::getDebugLevel()


//#define _debugCompile_
//! Debug error message, similar to ensure, but printed only if _debugCompile_ is defined and debugLevel is not zero
	#ifdef _debugCompile_
  //#define d_assert(...) ensure((!debugLevel)||__VA_ARGS__)
  #define d_assert(...) if(debugLevel) ensure(__VA_ARGS__)
	#else
	  #define d_assert(...)
#endif


#endif // Dbg_SkipH



#endif
