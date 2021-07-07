#ifndef TIMER_H
#define TIMER_H

#include <ctime>
#include <iostream>
#include <unordered_map>

// see performance.cpp for high resolution clock

struct Timer
{
	Timer() : elapsed(0.) {}
	void  tic() { clock_=clock(); }
	void  toc() { elapsed += double(clock()-clock_)/double(CLOCKS_PER_SEC); }
	clock_t     clock_;
	double      elapsed;
};

class Watch
{
	double&      elapsed_;
	clock_t     clock_;
 public:
	Watch(double& elapsed) : elapsed_(elapsed), clock_(clock()) {}
	~Watch() { elapsed_ += double(clock()-clock_)/double(CLOCKS_PER_SEC); }
};

//struct TimeScope
//{
	//TimeScope(double& elapsed) : passed_(elapsed), clock_(clock()) {}
	//~TimeScope() { elapsed += double(clock()-clock_)/double(CLOCKS_PER_SEC); }
	//clock_t  clock_;
	//double&  elapsed;
//};

#ifdef _PROFILE_
#define _tic_  timer.tic();
#define _toc_  timer.toc();
#else
#define _tic_  
#define _toc_  
#endif

#ifdef _PROFILE_
#define _tic_t(timeri)  timeri.tic();
#define _toc_t(timeri)  timeri.toc();
#else
#define _tic_t(timeri) 
#define _toc_t(timeri)  
#endif

inline double elapsedTime(clock_t start){  return (double(clock()-start)/CLOCKS_PER_SEC);  }

// used in network extraction code

class Timing  {
 public:
	Timing() : clock_(clock()) {}
	~Timing() {
		(*this)("");// insert last task
		std::cout<<"cpu times:\n";
		for(auto& itr:times) std::cout<<"  "<<itr.first<<": "<<itr.second<<std::endl;
	}

	void operator()(const std::string& cs)
	{
		if(task_.size())
		  times.insert({task_,0.}).first->second += elapsedTime(clock_); 
		clock_ = clock();
		task_ = cs;
	}

	clock_t               clock_;
	std::string           task_;
	std::unordered_map<std::string,double> times;
};




class Stats  {
	std::map<std::string,int> counts_;
 public:
	~Stats()  {
		std::cout<<"Stats:\n";
		for(auto& it:counts_) std::cout<<"  "<<it.first<<": "<<it.second<<std::endl;  }
	void count(const std::string& cs) {  counts_.insert({cs,0}).first->second += 1;  }
};	

_Extern	thread_local Stats _Stats;


#endif
