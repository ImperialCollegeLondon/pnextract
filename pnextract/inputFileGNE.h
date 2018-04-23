#ifndef INPUTFILEGNE_H
#define INPUTFILEGNE_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <set>
#include <vector>
#include <algorithm>
#include <functional>
#include <map>
#include <string>
#include <vec3.h>

using namespace std;


extern int debugLevel;



///. for debugger breakpoints: don't mess up !
inline void ErrorIfGNE(bool cond, string message="") {	if(cond)  {
	if (debugLevel)  cout<<" Error  "<< message <<endl;  }  }

//#define _debugCompile_
#ifndef _debugCompile_
	#define d_assert(cond) {}
#else
 #ifndef d_assert
	#define d_assert(cond)   ErrorIfGNE(!(cond),  \
					std::string("in ")+ std::string(__FUNCTION__)+  \
					", "+std::string(__FILE__)+":"+ toStr(__LINE__) \
					+std::string(":  ")+std::string(#cond) );
 #endif
#endif








class InputFileNE
{
	bool safeGetline(std::istream& is, std::string& t)
	{
		t.clear();
		std::istream::sentry se(is, true);
		std::streambuf* sb = is.rdbuf();
		for(;;) {
			int c = sb->sbumpc();
			switch (c) {
				case '\n':
					return false; ///. stopReading after new line
				case '%': /// allow Matlab comments
				case '/': /// allow C++ single-line comments
				case '#': /// allow Python and bash comments
					while (sb->sbumpc()!='\n');
					return false;
				case '\r':
					if(sb->sgetc() == '\n')    sb->sbumpc();
					return false;
				case EOF: // Also handle the case when the last line has no line ending
					if(t.empty())   is.setstate(std::ios::eofbit);
					return false;
				case ';':  return false;
				case '=':  t += '\t';   break;
				case ':':  t += '\t';   break;
				default:   t += (char)c;
			}
		}
	}

public:


	InputFileNE(const InputFileNE& input)
		:	_parsedData(input._parsedData), _fileName(input._fileName), _title(input._title), verbose(input.verbose)
	{
		if(verbose) std::cout<< " copied " << _title <<std::endl;
	}

   InputFileNE(const InputFileNE& input, const std::string& title)
	:	_parsedData(input._parsedData), _fileName(input._fileName), _title(title), verbose(input.verbose)
	{
		setKeyword("title",title);
		renameKeys("stage1","processed stage1 , ...");
		renameKeys("stage1_1","processed stage1_1 , ...");
		debugLevel=input.getOr(0,"debugLevel");

		if(verbose) std::cout<< " created " << _title <<std::endl;
	}

	InputFileNE(const std::string& fileName, int readfile=true)
	{
		_fileName = fileName;

		if(readfile) read(fileName);
		_parsedData.push_back(std::pair< std::string, std::string >("end", ""));///. add empty data

		if (!(getVar(_title, "title") || getVar(_title, "TITLE") || getVar(_title, "baseName")))
		{  _title = getOr(getOr(fileName,"ElementDataFile"),"imageFile");
			size_t dotloc=_title.find_last_of("."); if (dotloc<_title.size()) _title.erase(_title.find_last_of("."), string::npos);
			size_t slashloc=_title.find_last_of("\\/"); if (slashloc<_title.size()) _title=_title.substr(slashloc+1);
		}
		std::cout<< " output base name: " << _title <<std::endl;
		debugLevel=getOr(0,"debugLevel");

		if(verbose) std::cout<< " read " << _title <<std::endl;
	};

	void addKeyword(std::string key, std::string data)
	{
		if (key!="end")
		{
			if (_parsedData.rbegin()->first!="end") std::cout<<"Error in input data handelling"<<std::endl;
			_parsedData.rbegin()->first=key;   _parsedData.rbegin()->second=data;
		}

		_parsedData.push_back(std::pair< std::string, std::string >("end", ""));///. add empty data
	}
	void setKeyword(std::string key, std::string data)
	{
		if(verbose) cout<<"setting keyword "<<key<<" "<<data<<endl;
		for(unsigned i = 0; i < _parsedData.size(); ++i)
			if(_parsedData[i].first == key)
			{	_parsedData[i].second = data; 		return;	}

		addKeyword(key,data);///. empty string
	}
	void renameKeys(std::string key, std::string newkey)
	{
		if(verbose) cout<<"replacing keys "<<key<<" with "<<newkey<<endl;
		for(unsigned i = 0; i < _parsedData.size(); ++i)
			if(_parsedData[i].first == key)
				_parsedData[i].first = newkey; 
		return;

	}


   void read(std::string fileName, int importance=2)
	{
		verbose = false;

		std::string prevKeyword("NO_KEYWORD_READ");
		cout<< "Loading file: " << fileName; cout.flush();
		std::ifstream in(fileName.c_str());
		if (!in) {	std::cerr << "\n\n"<<(importance>1 ? "Error" : "Warning")<<": Unable to open input file, " << fileName << endl;	if (importance>1) exit(-1); }

		while(in.good())
		{
			std::string keyword, dataString;
			std::string bufferStr;
			bool oktocontinue=safeGetline(in,bufferStr);
			if(bufferStr.empty()) continue;


			{
				size_t beginKey=bufferStr.find_first_not_of(" \t");
				if (beginKey == std::string::npos) continue ;
				size_t keyEnding= bufferStr.find_last_not_of(" \t")+1;
				size_t endKey=bufferStr.find_first_of(":=");
				if (endKey == std::string::npos) endKey = bufferStr.find_first_of(" \t", beginKey+1);
				if (endKey == std::string::npos) continue ;

				keyword = bufferStr.substr(beginKey, endKey-beginKey);

				beginKey=bufferStr.find_first_not_of(" \t",endKey);
		
				if (beginKey != std::string::npos) dataString = bufferStr.substr(beginKey, keyEnding-beginKey)+'\n';

				if(keyword.size()<1 || keyword.size()>100)	{std::cerr << "\n\nError: Data file contains errors after keyword: " << prevKeyword << endl << bufferStr;	exit(-1);}
			}

			if (keyword=="end")					break;

			{
				while(oktocontinue)
				{
					oktocontinue=safeGetline(in,bufferStr);
					size_t beginKey=bufferStr.find_first_not_of(" \t");
					if (beginKey==std::string::npos) continue;
					size_t keyEnding=bufferStr.find_last_not_of(" \t");
					dataString += bufferStr.substr(beginKey, keyEnding-beginKey+1) +'\n';
				}



				prevKeyword = keyword;
				if (keyword=="debugLevel")
				{
					debugLevel=atoi(dataString.c_str());
					verbose = debugLevel%2;
				}

				_parsedData.push_back(std::pair<std::string, std::string>(keyword, dataString));
				if(verbose) cout<< keyword << ":\n" << _parsedData.rbegin()->second << endl ;
			}
		}

		in.close();
		cout<< endl;

	}





	inline void echoKeywords(std::ostream& out) const
	{
		for(unsigned i = 0; i < _parsedData.size(); ++i)
			out <<" "<< _parsedData[i].first  << ":\t"
				 << _parsedData[i].second  << ";"<< endl << endl;
	}

	inline const std::string& keywordData(const std::string& keyword, int importance=0) const
	{
		if(verbose) cout<<"Reading "<<keyword<<endl;
		for(unsigned i = 0; i < _parsedData.size(); ++i)
			if(_parsedData[i].first == keyword)   return (_parsedData[i].second);

		Assert(importance<1, keyword, "missing keyword", importance>1);
		return _parsedData[_parsedData.size()-1].second;///. empty string
	}
	bool getData(std::istringstream& data, const std::string& keyword, int importance=0) const
	{
		data.clear();
		for(unsigned i = 0; i < _parsedData.size(); ++i)
			if(_parsedData[i].first == keyword) {
				data.str(_parsedData[i].second);	return true; }

		Assert(importance<1, keyword, "missing keyword", importance>1);
		return false;
	}

	bool getVar(bool& var, const string& keyword) const
	{
		istringstream data;
		if(getData(data, keyword))
		{
			char varC;  data>>varC;  var = (varC=='T' || varC=='t' || varC=='Y' || varC=='y' || varC=='1');
			if(verbose) cout<<keyword<<" = "<<var<<";"<<endl;
			return true;
		}
		else
			return false;
	}

	template<typename Type>
    bool getVar(Type& var, const string& keyword) const
    {
		istringstream data;
		if(getData(data, keyword))
		{	data>>var;	if(verbose) std::cout<<" "<<keyword<<" = "<<var<<";"<<endl; return true; }
		else						return false;
	}

    template<typename Type>
	Type getOr(Type var, const string& keyword)	 const {  getVar(var, keyword);  return var;  }


	inline bool Assert(bool isOK, const std::string& keyword, const std::string extraMessage="", bool severe = true) const
	{
	  if (!isOK)
	  { std::cout << endl
			  << "========================================================"   << endl
			  << "Error: " << extraMessage                                    << endl
			  << "  File: " << _fileName                                     << endl
			  << "  Keyword: "<< keyword <<"  "<< keywordData(keyword,0)<<";" << endl
			  << "========================================================"   << endl;
		 if (severe) exit(-1);
	  }
	  return  isOK;
	}

	inline void checkEndOfData(std::istringstream& data, const std::string& keyword, bool severe = true) const
	{
		 Assert(!data.fail(), keyword,"Incomplete/bad data", severe);
		 char leftOvers; 	 data >> leftOvers;
		 Assert(!data, keyword,"Too much data", severe);
	}
	std::string baseName() const { return _title; }
	std::string fileName() const { return _fileName; }

protected:

	std::vector< std::pair< std::string, std::string > >        _parsedData;
	std::string                                                 _fileName;
	std::string                                                 _title;
	bool verbose;

};





#endif

