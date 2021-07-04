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



template<class T>  using  stvec    = std::vector<T>;
template<class T>  using  stvecs = std::vector<std::vector<T>>;
using  ststr    = std::string;
using  isstr    = std::istringstream;



inline int readInt(isstr& in) {int n=0; in>>n; return n; }
inline bool readBoolOr(std::string st, std::istream& in) { in>>st; return st[0]=='T' || st[0]=='t' || st[0]=='Y' || st[0]=='y' || st[0]=='1'; }

#ifndef for_i
#define for_i(_vector_m_)  for(size_t i=0; i<_vector_m_.size(); ++i)
#endif

class InputFile {//! InputFile is a general input file reader, with some flexibility to chose the keyword endings etc
 public:
	using string = std::string;

	InputFile(bool multiline=true)
	:	informative(true), multiline_(multiline)  {
		data_.push_back({string("end"), string()});//. add empty data at the end
	}


	InputFile(const string& fnam, bool multiline=true, bool inform=true, bool init=true)
	:	informative(inform), multiline_(multiline)  {

		if(fnam.empty()) {  data_.push_back({string("end"), string()});  return; }//  initialize empty InputFile on demand

		read(fnam);

		for_i(data_)
		 if( (data_[i].first=="include" || data_[i].first=="append" ) && !data_[i].second.empty())  {
			data_[i].first = "included";
			read(data_[i].second);
		 }

		if(init)  {
			string wdir=getOr("workingDir", string());
			if(!wdir.empty() && wdir!="PWD" &&  wdir!="pwd")  {
				if(wdir=="inputDir")  {
					size_t sl=fnam.find_last_of("\\/")+1;
					if(sl<fnam.size())  wdir=fnam.substr(0,sl);  }
				std::cout<<"Changing working directory: "<<wdir<<",  "<<
				_TRY_(chdir(wdir))            << std::endl;
			}

			setTitle(fnam);
		}
		if(informative) std::cout<<std::endl;
	}

	InputFile(const string& kwrds, const string& nam, bool multiline=false)
	:	informative(true), multiline_(multiline)  {  std::istringstream ss(kwrds);  read(ss, nam);  }

	InputFile(std::istream& in, const string& nam, bool multiline=false)
	:	informative(true), multiline_(multiline) 	{ 	read(in, nam);	}

   InputFile(const InputFile& input, const string& nam)
	:	data_(input.data_), fileName_(input.fileName_), folder_(input.folder_), name_(nam),
		informative(input.informative), multiline_(input.multiline_)  {	set("name",nam);	}

	bool safeGetline(std::istream& is, string& sr, bool noeq=false)  { /// \return readnext,  true: keep reading same key, if returned sr not empty
		sr.clear();
		auto begl = is.tellg();
		std::istream::sentry se(is, true);
		std::streambuf& sf = *is.rdbuf();
		for(;;) {
			int cr = sf.sbumpc(), cn=1;
			switch (cr) {
				case '\\': sr+=char(sf.sbumpc()); break; //. read next
				#ifdef HASH_ENDS_LINE // backward compatibility with poreflow
				case '#': return false;
				case '/':	if(sf.sgetc()!='/') {  sr += '/';  break;  } fallThrough;
				#else
				case '/':	if(sf.sgetc()!='/') {  sr += '/';  break;  } fallThrough;
				case '#': fallThrough;
				#endif
				case '%':
					while ((cr=sf.sbumpc())!='\n' && cr!=EOF);
					sr += '\t';
					return multiline_;

				case '=': case ':':
					if(noeq)  {  sr.clear(); is.seekg(begl);  return false;  }
					sr += '\t';         break;
				case ',':  sr += '\t';  break;  //! comma is treated as white space 

				case EOF:  if(sr.empty())  is.setstate(std::ios::eofbit);  return false;

				case '{':
					if(sr.size() && !noeq) { sf.sungetc(); return true; } // try again with empty noeq=true (or empty sr to skip {})
					else if(noeq){ cr='\t';  do{ sr+=cr; cr=sf.sbumpc(); cn+=int(cr=='{')-(cr=='}'); }while(cn && cr!=EOF);  }
					fallThrough;
				case '}':  sr += '\t';  break; //return  false;
				//case '{':  cr='\t'; do{ sr+=cr; cr=sf.sbumpc(); }while(cr!='}' && cr!=EOF); cr='\t';  break;
				case '\'': cr='\t'; do{ sr+=cr; cr=sf.sbumpc(); }while( cr!='\''&& cr!=EOF );     break;
				case '"':  cr='\t'; do{ sr+=cr; cr=sf.sbumpc(); }while( cr!='"' && cr!=EOF );     break;
				case '[':  cr='\t'; do{ sr+=cr; cr=sf.sbumpc(); cn+=int(cr=='[')-(cr==']'); }while(cn && cr!=EOF);    break;
				case '(':  cr='\t'; do{ sr+=cr; cr=sf.sbumpc(); cn+=int(cr=='(')-(cr==')'); }while(cn && cr!=EOF);    break;
				case '<':  cr='\t'; do{         cr=sf.sbumpc(); cn+=int(cr=='<')-(cr=='>'); }while(cn && cr!=EOF);    break; // primitive: filter out  xml tags

				case '\r':
					if(sf.sgetc()=='\n') sf.sbumpc(); fallThrough;
				case '\n':
					cr=sf.sgetc();
					return (cr=='\n' || cr=='\r') ? false  : multiline_; //! double new lines are always treated as end of keyword


				case ';':  return  false;
				default :  sr += char(cr); 
			}
		}
	}


	int read(const string& fnam, int importance=2)  {

		if(informative) (std::cout<< "/""/-*-C++-*-input: "+fnam+";  ").flush();
		std::ifstream in(fnam);
		if (in)  return read(in, fnam);

		ensure(importance==0, "can not open " + fnam, importance);
		return 0;
	}

	int read(std::istream& in, string fnam="")  {

		if(fnam.size()) fileName_ = fnam;
		if(data_.size() && data_.back().first=="end") data_.pop_back();

		string prev("NO_KEYWORD_READ");
		while(in.good())  {
			string key, kydata, bufr;
			bool readnext=safeGetline(in,bufr,false);
			if(bufr.empty())                           continue;

			{
				size_t bgn=bufr.find_first_not_of(" \t");
				if (bgn == string::npos)                continue;

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
			while(readnext)  {
				readnext=safeGetline(in,bufr, true);
				size_t bgn=bufr.find_first_not_of(" \t");
				if (bgn==string::npos) continue;
				size_t lst=bufr.find_last_not_of(" \t");
				if(kydata.size())  {
					if(kydata[0]!='\n') kydata = "\n"+kydata;
					kydata += "\n";
				}
				kydata += bufr.substr(bgn, lst-bgn+1) ;
			}

			prev = key;

			IN_DBG ( if (key=="debugLevel")  debuglevel_(atoi(kydata.c_str())); )

			data_.push_back({key, kydata});
		}

		data_.push_back({string("end"), string()});//. add empty data at the end
		return 0;
	}


	void setTitle(string fnam="")   {	// call this after setting name and/or prefix
		string prf=getOr("prefix", string());
		if(prf.size())  {
			folder_.resize(0);
			size_t slashloc=prf.find_last_of("\\/");
			if (slashloc<prf.size()-1)  {
				folder_=prf.substr(0,slashloc+1);
				std::cout<<"Creating folder: "<<folder_<<", "
					<< _TRY_(mkdirs(folder_))
					<<std::endl;
				prf=prf.substr(slashloc+1);
			}

			if (prf.size()>1 && (prf.back()=='/' || prf.back()=='\\'))  {
				folder_+=prf;
				std::cout<<"Creating folder: "<<folder_<<"  "
					<< _TRY_(mkdirs(folder_))
					<<std::endl;
				prf="";
			}
		}

		if( lookup("name",name_) || lookup("TITLE",name_) || lookup("title",name_) ) name_ = prf+name_;
		else if(prf.size()) name_ = prf;
		else
		{	prf = getOr("network", getOr("networkFile", string("")));
			if (prf.empty()) { prf = fnam; }
			if (prf.empty()) { if(!lookup("ElementDataFile",prf)) prf = fileName(); }
			if (prf.size()>7 && prf.substr(prf.size()-3,3)==".gz") prf = prf.substr(0,prf.size()-3); // sync bname
			size_t dl=prf.find_last_of(".");   if (dl<prf.size()) prf.erase(dl);
			size_t sl=prf.find_last_of("\\/"); if (sl<prf.size()) prf=prf.substr(sl+1);
			name_ = prf;
		}
		std::cout<< " output prefix: " << name_<<" ";
	}


	#ifdef _debugCompile_
	 #define _debugInfo_(pref_) if(debugLevel) std::cout<<pref_+key+":"+data_[i].second<<std::endl
	#else
	 #define _debugInfo_(pref_) 
	#endif


	void echoKeywords(std::ostream& out=std::cout) const  {
		char el=';';
		if(out.tellp()==0)  {el='\n'; out<<"{/""/ -*- C -*- "<<outputName()<<" input keywords:\n\n"; }
		for_i(data_) {
			out <<" "<< data_[i].first  << ":\t";
			if(data_[i].second.find_first_of("[:;(\"\'{")!=string::npos)  out<<"{ "<<data_[i].second<<" }"<<el<<"\n";
			else out<< data_[i].second  <<el<< "\n";
		}
		if(el=='\n')  out<<"}"<<std::endl;
	}

	void add(string key, string val)  {
		ensure(data_.back().first=="end","input file handelling, overwriting "+data_.back().first);
		data_.back() = {key,val};
		data_.push_back({string("end"), string()});//. add empty data at the end
	}

	void set(string key, string val, bool overwrite=true)  {
		for_i(data_)  if(data_[i].first == key)  { if(overwrite) { _debugInfo_("Setting "); data_[i].second = val; } return; }
		add(key,val);
	}
	void setDefault(string key, string val)  { set(key,val,false); }

	string& operator[](string key)  { /// not recommended: cannot be refactored
		for_i(data_)  if(data_[i].first == key)  { _debugInfo_("[] "); return data_[i].second; }
		add(key,"");
		ensure(data_.size()>=2 && data_[data_.size()-2].first==key);
		return data_[data_.size()-2].second;
	}

	void renameKeys(string key, string newkey)  {
		for_i(data_)  if(data_[i].first == key)  { _debugInfo_("Renaming "); data_[i].first = newkey; }
	}



	const string& kwrd(const string& key, int importance=0) const  { //! get
		for_i(data_) if(data_[i].first == key) { _debugInfo_("Reading ");     return (data_[i].second);  }
		Assert(importance<1, key, "missing keyword", importance>1);
		return data_.back().second;//. empty string
	}

	bool giv(const string& key, isstr& iss, int importance=0) const  { //! give me
		iss.clear();
		for_i(data_)  if(data_[i].first == key) { _debugInfo_("Reading ");  iss.str(data_[i].second);  return true; }
		Assert(importance<1, key, "missing keyword", importance>1);
		return false;
	}
	template<class T> bool giv(const string& key, T& var, int importance=0) const  {//! give me
		isstr iss;
		if(giv(key, iss,importance))  {  iss>>var; return true; }
		else  return false;
	}
	bool giv(const string& key, bool& var, int importance=0) const {
		isstr iss;
		if(giv(key, iss, importance)){  char c; iss>>c;  var=(c=='T'||c=='t'||c=='Y'||c=='y'||c=='1'); return true;  }
		else  return false;
	}
	template<class T> bool lookup(const string& key, T& var, int importance=0) const {  return giv(key,var);  } /// similar to openfoam
	template<class T> T    getOr(const string& key, T var)  const {  giv(key, var);  return var;  }



	void Assert(bool isOK, const string& key, const string message="", bool severe = true) const  {
		if(!isOK)  { std::cerr<<"\n\n"<<(severe?"Error":"Warning")<<" in file "+fileName()+", keyword:\n  "+key<<": "+kwrd(key,0)+";\n  "+message+"\n"<<std::endl;
			if (severe) exit(EXIT_FAILURE); }
	}

	void checkEndOfData(isstr& iss, const string& key, const string message="", bool severe = true) const  {
		Assert(!iss.fail(), key,"Incomplete/wrong data", severe);  char c;  Assert(!iss.get(c), key,"Too much data", severe);  
	}

	string outputName() const { return folder_+name_; }
	string prefix()     const { return folder_; }
	string name()       const { return name_; }
	string fileName()   const { return fileName_; } // optional


	const std::vector< std::pair<string,string>>&  data() const { return data_; };

protected:

	std::vector<std::pair<string,string>>   data_; //TODO use together with unordered_multimap<string_view,string*> to speed up search (?)
	string                                  fileName_; // optional
	string                                  folder_;
	string                                  name_;

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
class mstream  {
 public:
	enum: short { PRTF = 1, STDO = 2 }; // .prt and stdout

	thread_local 
		static std::ofstream  dbgFile; /// output stream to write .dbg  file for debugging

	mstream(std::string fnam, unsigned char po=PRTF|STDO)
	:	outps(po) { if (!fnam.empty() && outps&PRTF) prt_.open(fnam); if(!prt_) outps&=~PRTF; };

	mstream& operator<<(std::ostream& (*fn)(std::ostream&))  {
		if(outps&PRTF) { fn(prt_);  }  if(outps&STDO) { fn(std::cout); std::cout.flush(); }  return *this;  }

	template <class T>  mstream& operator<< (T val)  {
		if(outps&PRTF) { prt_<<val; }  if(outps&STDO) { std::cout<< val; }   return *this;  }

	std::ofstream& fileStream() { return prt_; }

 private:
	std::ofstream prt_;
	unsigned char outps; /// output .prt std::out
};


#define outD mstream::dbgFile


#endif
