#ifndef GLOBALS_SkipH
#define GLOBALS_SkipH

#include <sstream>
#include <iostream>
#include <stdexcept>
#include <memory>
#include <map>
#include <filesystem>  // Include the standard filesystem header

namespace fs = std::filesystem;  // Create an alias for convenience

inline std::string getpwd() { return fs::current_path().string(); }
inline int mkdirs(const std::string& dir) { std::error_code ec; fs::create_directories(dir, ec); return ec.value(); }
inline int chdir(const std::string& dir) { std::error_code ec; fs::current_path(dir, ec); return ec.value(); }

#define _TRY_(_syscmnd_) std::string((_syscmnd_==0) ? " succeed " : " failed ")

// Non-folding brackets for namespace '{' and '}', to be used in early stages of code development
#define _begins_       {
#define _end_of_(sec)  }

// used by SiRun and VxlPro
#define KeyHint(_args_desc_)  if(ins.peek()=='?') { ins.str(_args_desc_); return 0; }

#if defined __has_cpp_attribute 
    #if __has_cpp_attribute(fallthrough)
        #define fallThrough [[fallthrough]]
    #else
        #define fallThrough
    #endif
#else
    #define fallThrough
#endif

#ifdef FOOL_GEANY_TO_COLOUR_HIGHLIGHT
class vector {};  class for_ {};  class for_i {};  class fori0to {};
class iterator {}; class fluidf {}; class Elem{}; class dbl {}; class string {};
class map {}; class T {}; class endl {}; class cout {}; class cerr {};
namespace std {}
#endif 

#ifndef Str_SkipH  // string utilities
#define Str_SkipH
    //using _s = std::to_string  is bad in decimal notation
    template<typename T> std::string toStr(const T& n){  std::ostringstream ss;  ss<<n;  return ss.str();  } //outdated
    template<typename T> std::string _s(const T& n){  std::ostringstream ss;  ss<<n;  return ss.str();  }
    template<typename T> T strTo(const std::string &s){  std::istringstream ss(s);  T t;  ss>>t;  return t; }
    template<typename T> bool hasExt(const T& path, size_t siz, const char* ext) { return path.size()>siz && path.compare(path.size()-siz, siz, ext) == 0; }
    template<typename T> bool hasExt(const T& path, const std::string& ext) { return path.size()>ext.size() && path.compare(path.size()-ext.size(), ext.size(), ext) == 0; }

    //stringify macro args after (sub)macro expansion
    #define STRINGIFY(x) #x
    #define TOSTRING(x) STRINGIFY(x)

    #define __FILEBASENAME__ std::string(strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__, strcspn(strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__, "."))
#endif // Str_SkipH

#ifndef Tst_SkipH  // ensure, alert (and ifnot hack)  runtime checks 
#define Tst_SkipH
    // Run-time check and testing macros
    inline bool _cerr_(std::string hdr="", std::string msg="", int xit=0) // for debugger breakpoints: don't optimize out please !!!
        {   if(xit) throw std::runtime_error(hdr+msg);  else std::cerr<< hdr+msg <<std::endl; return true; }

    #define ERR_HDR(_xit_)  std::string(_xit_?"\n\n  Error":"\n  Warning") \
                +" in "+ std::string(__FUNCTION__)+", "+std::string(__FILE__)+":"+_s(__LINE__)  \
                ,std::string(":  ")
    #define ensure1(isOk)           (!((isOk)|| _cerr_(ERR_HDR(0)+std::string(#isOk))))
    #define ensure2(isOk, msg)      (!((isOk)|| _cerr_(ERR_HDR(0)+  "   \""+msg+"\"")))
    #define ensure3(isOk, msg, xit) (!((isOk)|| _cerr_(ERR_HDR(xit)+"\n   "+msg+"\n", xit)))
    #define GET_MACRO3(_1,_2,_3,NAME,...) NAME

    //! Validation/production phase ensure/assert. Usage:
    //!   \code{.cpp} ensure(condition, "message", int throw_on_error=false); \endcode
    #define ensure(...)  GET_MACRO3(__VA_ARGS__, ensure3, ensure2, ensure1, "Only 1 to 3 args please")(__VA_ARGS__)
    #define alert(...)  GET_MACRO3(false,__VA_ARGS__, ensure3, ensure2, "Only 1 to 2 args please")(false,__VA_ARGS__)

    #define ifnot(isOk, msg)  if(!(isOk)&& std::cout<<msg<<endl)

#endif // Tst_SkipH

#ifdef _InitGlobals  // Allow declaration and definition of static variables in the same file
    #define _Extern 
    #define _Eq(...)  = __VA_ARGS__
#else
    #define _Extern extern 
    #define _Eq(...)  
#endif

//////////////////////////// OUTDATED \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\.
#ifdef _debugCompile_ // by default do not use the obsolete debugLevel/dbgAsrt
    template<class T> int debuglevel_(T level) {  static int l_=0;  if (level>=0) { l_=level; }  return l_;  }
    #define debugLevel    debuglevel_(-1)
    #define IN_DBG(...) __VA_ARGS__
    #define dbgAsrt(...) (debugLevel<=0 || ensure(__VA_ARGS__))  // debug assert with message
#else
    #define debugLevel  0
    #define dbgAsrt(...)
    #define IN_DBG(...)  
#endif // _debugCompile_
// Note use addr2line more for debugging

#endif // GLOBALS_SkipH

