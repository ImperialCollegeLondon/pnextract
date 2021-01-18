#ifndef INPUTFILE_H
#define INPUTFILE_H

// Input data file used by 3D image processing, network extraction, 
// flow simulation  and other codes
// Developed by Ali Q. Raeini. See the documentation of the 
// relevant codes for user guids and contact details, 



#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <set>
#include <vector>
#include <array>
#include <algorithm>
#include <functional>
#include <map>
#include <string>
#include "globals.h"


#if defined _MSC_VER
 #include <direct.h>
#elif defined __GNUC__
 #include <sys/types.h>
 #include <sys/stat.h>
 #include <unistd.h>
#endif

#define  _TRY_(_syscmnd_) ((std::array<std::string,2>{"failed.", "succeed."})[(_syscmnd_)==0])

template<class T> 
using  stvec = std::vector<T>;
using  ststr = std::string;
using  isstr = std::istringstream;



inline int readInt(isstr& in) {int n = 0; in>>n; return n;}
inline bool readBoolOr(std::string st, std::istream& in) { in>>st; return st[0]=='T' || st[0]=='t' || st[0]=='Y' || st[0]=='y' || st[0]=='1';}

#ifndef for_i
#define for_i(_vector_m_)  for(size_t i=0;i<_vector_m_.size();++i)
#endif

class InputFile //! InputFile is a general input file reader, with some flexibility to chose the keyword endings etc
{
public:
	using string = std::string;

	InputFile(bool multiline=true)
	:	informative(true), multiline_(multiline)
	{
		data_.push_back({string("end"), string()});//. add empty data at the end
	};

   InputFile(const InputFile& input, const string& title)
	:	data_(input.data_), fileName_(input.fileName_), folder_(input.folder_), name_(title)
	  , informative(input.informative), multiline_(input.multiline_)
	{
		setKeyword("name",title);
		#ifdef Dbg_SkipH
		globals::setDebugLevel(input.getOr(debugLevel,"debugLevel"));
		#endif
	}

	InputFile(const string& fnam, bool multiline=true, bool inform=true, bool init=true)
	:	informative(inform), multiline_(multiline)
	{
		read(fnam);

		for_i(data_)
			if( (data_[i].first=="include" || data_[i].first=="append" ) && !data_[i].second.empty() )
			{
				data_[i].first = "included";
				read(data_[i].second);
			}
		if(init)
		{
			string wdir=getOr(string(),"workingDir");
			if(!wdir.empty() && wdir!="PWD" &&  wdir!="pwd")
			{
				if(wdir=="inputDir")  {
					size_t eslshp=fnam.find_last_of("\\/")+1;
					if(eslshp<fnam.size())  wdir=fnam.substr(0,eslshp);  }
				std::cout<<"Changing working directory: "<<wdir<<",  "<<
				_TRY_(chdir(wdir.c_str()))            << std::endl;
			}

			setTitle(fnam);
		}
		if(informative) std::cout<<std::endl;
	};

	bool safeGetline(std::istream& is, string& sr, bool noeq=false)
	{
		sr.clear();
		auto begl = is.tellg();
		std::istream::sentry se(is, true);
		std::streambuf& sf = *is.rdbuf();
		for(;;) {
			int cr = sf.sbumpc(), cn=1;
			switch (cr) {
				case '\\': sr+=char(sf.sbumpc()); break; //. read next
				#ifdef HASH_ENDS_LINE // backward compatibility with porenflow
				case '#': return false;
				case '/':	if(sf.sgetc()!='/') {  sr += '/';  break;  }
				#else
				case '/':	if(sf.sgetc()!='/') {  sr += '/';  break;  }
				case '#':
				#endif
				case '%':
					while (sf.sbumpc()!='\n');
					sr += '\t';
					return multiline_;

				case '=':
				case ':':
					if(noeq)  {  sr.clear(); is.seekg(begl);  return false;  }
					sr += '\t';         return true;;
				case ',':  sr += '\t';  break;

				case EOF:  if(sr.empty())  is.setstate(std::ios::eofbit);  return false;

				case '{': case '}': cr='\t';  break;
				//case '{':  cr='\t'; do{ sr+=cr; cr=sf.sbumpc(); }while(cr!='}' && cr!=EOF); cr='\t';  break;
				case '\'': cr='\t'; do{ sr+=cr; cr=sf.sbumpc(); }while( cr!='\''&& cr!=EOF );   cr='\t';  break;
				case '"':  cr='\t'; do{ sr+=cr; cr=sf.sbumpc(); }while( cr!='"' && cr!=EOF );   cr='\t';  break;
				case '[':  cr='\t'; do{ sr+=cr; cr=sf.sbumpc(); cn+=int(cr=='[')-(cr==']'); }while(cn && cr!=EOF);  cr='\t';  break;
				case '(':  cr='\t'; do{ sr+=cr; cr=sf.sbumpc(); cn+=int(cr=='(')-(cr==')'); }while(cn && cr!=EOF);  cr='\t';  break;

				case '\r':
					if(sf.sgetc()=='\n') sf.sbumpc();
				case '\n':
					cr=sf.sgetc();
					return (cr=='\n' || cr=='\r') ? false  : multiline_; //! double new lines are always treated as end of keyword


				case ';':  return  false;
				default:   sr += char(cr);
			}
		}
	}


	int read(const string& fnam, int importance=2)
	{
		fileName_ = fnam;

		if(informative) (std::cout<< "/""/ -*- C -*- input: "+fnam+" ").flush();
		std::ifstream in(fnam);
		if (in)  return read(in);

		ensure(importance==0, "can not open " + fnam, importance);
		return 0;
	}

	int read(std::istream& in)
	{
		if(data_.size() && data_.back().first=="end") data_.pop_back();

		string prev("NO_KEYWORD_READ");
		while(in.good())
		{
			string key, kydata, bufr;
			bool readnext=safeGetline(in,bufr,false);
			if(bufr.empty())            continue;

			{
				size_t bgn=bufr.find_first_not_of(" \t");
				if (bgn == string::npos)   continue ;
				size_t lst= bufr.find_last_not_of(" \t")+1;

				size_t endKy=
				 #ifndef ALLOW_SPACE_IN_KEY // allowing space in key is a source of user error, deactivated by default
				  std::min(bufr.find_first_of(":= \t", bgn+1), lst);
				 #else
				  bufr.find_first_of(":=", bgn+1);  if(endKy==string::npos) endKy = std::min(bufr.find_first_of(" \t",bgn+1), lst);
				 #endif

				key = bufr.substr(bgn, endKy-bgn);
				if(key=="end")				break;
				ensure(key.size()>0 && key.size()<100, " after '"+prev+"'\n'"+bufr+"'\n@["+_s(endKy)+" "+_s(bgn)+"]",-1);

				bgn=bufr.find_first_not_of(" \t",endKy);
				if (bgn != string::npos) kydata = bufr.substr(bgn, lst-bgn);
			}
			while(readnext)
			{
				readnext=safeGetline(in,bufr, true);
				size_t bgn=bufr.find_first_not_of(" \t");
				if (bgn==string::npos) continue;
				size_t lst=bufr.find_last_not_of(" \t");
				if(kydata.size())
				{
					if(kydata[0]!='\n') kydata = "\n"+kydata;
					kydata += "\n";
				}
				kydata += bufr.substr(bgn, lst-bgn+1) ;
			}

			prev = key;

			#ifdef Dbg_SkipH
			if (key=="debugLevel")  globals::setDebugLevel(atoi(kydata.c_str()));
			#endif

			data_.push_back({key, kydata});
		}

		data_.push_back({string("end"), string()});//. add empty data at the end
		return 1;
	}


	void setTitle(string fnam="")
	{	// call this after setting name and/or prefix
		string prf=getOr(string(),"prefix");
		if(prf.size())
		{
			folder_.resize(0);
			size_t slashloc=prf.find_first_of("\\/");
			if (slashloc<prf.size()-1)
			{
				folder_=prf.substr(0,slashloc+1);
				std::cout<<"Creating folder: "<<folder_<<", "
				#if defined(_WIN32)
					<< _TRY_(mkdir(folder_.c_str())) //. check also _mkdir
				#else
					<< _TRY_(mkdir(folder_.c_str(), 0733))  // note  0777 is octal
				#endif
					<<std::endl;
				prf=prf.substr(slashloc+1);
			}

			if (prf.size()>1 && (prf.back()=='/' || prf.back()=='\\') )
			{
				folder_+=prf;
				std::cout<<"Creating folder: "<<folder_<<"  "
					#if defined(_WIN32)
						<< mkdir(folder_.c_str()) //. check also _mkdir
					#else
						<< mkdir(folder_.c_str(), 0733) // note  0777 is octal
					#endif
					<<std::endl;
				prf="";
			}
		}

		if( lookup("name",name_) || lookup("TITLE",name_) || lookup("title",name_) ) name_ = prf+name_;
		else if(prf.size()) name_ = prf;
		else
		{	prf = getOr(getOr(string(""),"network"),"networkFile");
			if (prf.empty()) { prf = fnam; }
			if (prf.empty()) { if(!lookup("ElementDataFile",prf)) prf = fileName_; }
			if (prf.size()>7 && prf.substr(prf.size()-3,3)==".gz") prf = prf.substr(0,prf.size()-3);
			size_t dl=prf.find_last_of(".");   if (dl<prf.size()) prf.erase(dl);
			size_t sl=prf.find_last_of("\\/"); if (sl<prf.size()) prf=prf.substr(sl+1);
			name_ = prf;
		}
		std::cout<< " output prefix: " << name_<<" ";
	}


	#ifdef Dbg_SkipH
	 #define _debugInfo_(pref_) if(debugLevel) std::cout<<pref_+key+":"+data_[i].second<<std::endl
	#else
	 #define _debugInfo_(pref_) 
	#endif


	void echoKeywords(std::ostream& out=std::cout) const
	{	int po=out.tellp();
		if(po==0)  out<<"{/""/ -*- C -*- "<<outputName()<<" input keywords:\n";
		out<<"\n";
		for_i(data_) {
			out <<" "<< data_[i].first  << ":\t";
			if(data_[i].second.find_first_of("[:;(\"\'") != string::npos)
				out<<'['<<data_[i].second  <<']'<< "\n\n";
			else out<< data_[i].second  << "\n\n";
		}
		if(po==0)  out<<"}";
		out<< std::endl;
	}
	void addKeyword(string key, string val)
	{
		ensure(data_.back().first=="end","input file handelling");
		data_.back() = {key,val};
		data_.push_back({string("end"), string()});//. add empty data at the end
	}
	void setKeyword(string key, string val, bool overwrite=true)
	{
		for_i(data_)  if(data_[i].first == key)  { if(overwrite) { data_[i].second = val; _debugInfo_("Setting "); } return; }
		addKeyword(key,val);
	}
	string& operator[](string key)
	{
		for_i(data_)  if(data_[i].first == key)  { _debugInfo_("[] "); return data_[i].second; }
		addKeyword(key,"");
		ensure(data_.size()>=2 && data_[data_.size()-2].first==key);
		return data_[data_.size()-2].second;
	}
	void renameKeys(string key, string newkey)
	{
		for_i(data_)  if(data_[i].first == key)  { _debugInfo_("Renaming "); data_[i].first = newkey; }
		return;
	}

	const string& keyvals(const string& key, int importance=0) const
	{
		for_i(data_) if(data_[i].first == key) { _debugInfo_("Reading ");     return (data_[i].second);  }
		Assert(importance<1, key, "missing keyword", importance>1);
		return data_.back().second;//. empty string
	}

	bool getData(isstr& iss, const string& key, int importance=0) const
	{
		iss.clear();
		for_i(data_)  if(data_[i].first == key) { _debugInfo_("Reading ");  iss.str(data_[i].second);  return true; }
		Assert(importance<1, key, "missing keyword", importance>1);
		return false;
	}
	bool giv(const string& key, isstr& iss, int importance=0) const { return getData(iss, key, importance); }  //! get me

	template<class T> bool lookup(const string& key, T& var) const
	{
		isstr iss;
		if(getData(iss, key))  {  iss>>var; return true; }
		else  return false;
	}
	bool lookup(const string& key, bool& var) const
	{
		isstr iss;
		if(getData(iss, key)){  char c; iss>>c;  var=(c=='T'||c=='t'||c=='Y'||c=='y'||c=='1'); return true;  }
		else  return false;
	}
	template<class T> T  lookupOr(const string& key, T var)  const {  lookup(key, var);  return var;  } // similar to openfoam
	template<class T> bool getVar(T& var, const string& key) const {  return lookup(key,var);  } // obselete
	template<class T> T     getOr(T var, const string& key)  const {  return lookupOr(key, var);  } // obselete



	void Assert(bool isOK, const string& key, const string message="", bool severe = true) const
	{	if(!isOK)  { std::cerr<<"\n\n"<<(severe?"Error":"Warning")<<" in file "+fileName_+", keyword:\n  "+key<<": "+keyvals(key,0)+";\n  "+message+"\n"<<std::endl;
			if (severe) exit(EXIT_FAILURE); }
	}

	void checkEndOfData(isstr& iss, const string& key, const string message="", bool severe = true) const
	{  Assert(!iss.fail(), key,"Incomplete/wrong data", severe);  char c;  Assert(!iss.get(c), key,"Too much data", severe);  
	}

	string outputName() const { return folder_+name_; }
	string prefix() const { return folder_; }
	string name() const { return name_; }
	string fileName() const { return fileName_; }


	const std::vector< std::pair<string,string> >&  data() const {return data_;};

protected:

	std::vector< std::pair<string,string> >   data_;
	string                                    fileName_;
	string                                    folder_;
	string                                    name_;

public:	//. extra:
	bool                           informative;
	bool                           multiline_;
};



#endif





//TODO: move to separate file
#ifndef ONDEMANDSTREAM_H
#define ONDEMANDSTREAM_H

#define defaultFloat  std::resetiosflags(std::ios_base::floatfield)<<std::setprecision(6)

/// output stream to write both to std::cout, and .prt file
class mstream
{
 public:
	enum: unsigned char {PRTF = 1, STDO = 2};
	std::ofstream prt_;
	unsigned char outps; /// output .prt std::out
	mstream(std::string fnam, unsigned char po=PRTF|STDO)
	: outps(po) { if (!fnam.empty()) prt_.open(fnam.c_str()); if(!prt_) outps&=~PRTF;};
	~mstream(void){};
	mstream& operator<<(std::ostream& (*fn)(std::ostream&)) {if(outps&PRTF) fn(prt_); if(outps&STDO) fn(std::cout); return *this;}
	std::ofstream& fileStream() {return prt_;};
};

template <class T>
mstream& operator<< (mstream& st, T val)
{
  if(st.outps & mstream::PRTF) st.prt_ << val;
  if(st.outps & mstream::STDO) std::cout<< val;
  return st;
}

/// output stream to write .dbg  file for debugging
class OnDemandStream
{
	public:
	std::ofstream prt_;
	bool opened;
	OnDemandStream(): opened(false){};
	OnDemandStream(std::string fnam,int dbgMode) 
	: opened(dbgMode) { if (dbgMode && !fnam.empty()) { prt_.open(fnam);  opened=true;} };
	~OnDemandStream(void){};
	OnDemandStream& operator<<(std::ostream& (*fn)(std::ostream&)) { fn(prt_); return *this;}
	void open(std::string fnam) {if (!fnam.empty()) {prt_.open(fnam); opened=true;}};
	void close() {if (opened) {prt_.close(); opened=false;}};
	void flush() {if (opened) {prt_.flush();} std::cout.flush(); };

	//#
	thread_local 
	//#
		static OnDemandStream    dbgFile;

};

#define outD OnDemandStream::dbgFile


template <class T>
OnDemandStream& operator<< (OnDemandStream& st, T val)
{
	#ifdef _debugCompile_
	if (st.opened)	{  (st.prt_) << val;  st.prt_.flush(); }
	#endif
  return st;
}
template <>
inline OnDemandStream& operator<< <double> (OnDemandStream& st, double val)
{
	#ifdef _debugCompile_
	if (st.opened)  {  if(val>1.e10)  (st.prt_)<<"+~";  else if(val<-1.e10)  (st.prt_)<<"-~";  else  (st.prt_)<<val;  st.flush();  }
	#endif
	return st;
}


#endif
