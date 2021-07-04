#ifndef ELEMENT_H
#define ELEMENT_H


#include "typses.h"
#include <map>
#include <iostream>
#include <array>


const static int LEVEL_MAX = 32767;
#define nAprox 1





class medialBall;
class voxel  {
 public:
	voxel(short ii, short jj, short kk, float rR) : ball(nullptr), i(ii), j(jj), k(kk), R(rR) {};
	voxel() : ball(nullptr),i(-1), j(-1), k(-1),R(0) {};
 public:
	medialBall* ball;
	short i, j, k;
	float R;
};

//class vxlface
//{
//public:
	//vxlface(short i, short j, short k, char dirr ) :  _ijk({{i,j,k}}),dir(dirr) {};//,hasBall(false)
	//vxlface() : _ijk({{-1,-1,-1}}),dir(0) {};
	//const array<short,3>& ijk() {return _ijk;};
	//array<short,3> ijkp() {array<short,3> ijkp1 = _ijk; ijkp1[dir]+=1; return ijkp1;};
//private:
	//array<short,3> _ijk;
	//char dir;
//};



class segment {
 public:
	voxel* v(int i){return segV+(i-start);};
	int           start;
	unsigned char value;
	voxel*         segV;
};


class segments  {
 public:
	segments(): s(nullptr), cnt(0) {};
	void reSize(int size){ 
		if (!s) s = new segment[size];
		else {std::cout<<" \nError in segments "<<size_t(s)<<std::endl; s = new segment[size];} }

	~segments() { if (s && cnt) {delete[] s; s = nullptr;} }

	voxel* vxl(int i) const  {
		for (int p=0; p<cnt; ++p)  if (i>=s[p].start && i<s[p+1].start)	return s[p].v(i);
		return nullptr; }

	segment* s;
	int cnt;
};





#define _mp5 -0.5

/// maximal-sphere class, describing spheres, which are generated for voxels on the medial surface
class medialBall  {
 public:
	medialBall() = delete;
	medialBall(short t) : vxl(nullptr), fi(-10000), fj(-0.5), fk(-10000),type(t), R(-10000), nKids(0), nNeis(0), corId(0), boss(this), kids(nullptr), neis(nullptr){};
	medialBall(voxel* v, short t) : vxl(v), fi(v->i-_mp5), fj(v->j-_mp5), fk(v->k-_mp5), type(t), R(v->R), nKids(0), nNeis(0), corId(0), boss(this), kids(nullptr), neis(nullptr){};
	~medialBall()  { if(kids) delete[] kids; if(neis) delete[] neis; }


	short level() const { if (this==boss)  return 1;  else  return boss->level()+1; }
	bool inParents(medialBall* vj) const {
		if (vj==boss)         return true;
		else  if(boss==this)  return false;
		else                  return boss->inParents(vj);
	}


	medialBall* mastrSphere() const  {  if (this==boss) return boss; else return boss->mastrSphere();  }

	//inline const medialBall* smallestP() const
	//{	const medialBall* pMin = this;
		//const medialBall* v = this;
		//while (v != v->boss)
		//{	pMin = v->R < pMin->R ? v : pMin;   v=v->boss;  }
		//return pMin;
	//}

	bool isNei(const medialBall* bal) const {
		for (int ii=0; ii<nNeis; ++ii)  if(neis[ii]==bal)  return true;
		return false;
	}
	//bool neisNotNeiWith(const medialBall* bal) const
	//{	for (int ii=0;ii<nNeis;++ii)  if (neis[ii]->isNei(bal))  return false;
		//return true;
	//}
	//medialBall*& neiRef(const medialBall* bal)
	//{	for (int ii=0;ii<nNeis;++ii)  if (neis[ii]==bal)  return neis[ii];
		//std::cout<<"no nei found! "<<int(nNeis)<<std::endl;
		//return neis[1000];
	//}
	void removeKidBoss(medialBall* kid)  {
		if (nKids>1) {
			int ii=-1;
			for (int kk=0;kk<nKids;++kk)    if(kids[kk]!=kid)   kids[++ii]=kids[kk];
			nKids=ii+1;
		}
		else if (nKids==1 && kids[0]==kid) { delete[] kids; kids=nullptr; nKids=0; }
		kid->boss = boss;
	}
	//bool isKid(const medialBall* kid) const
	//{	for (int ii=0;ii<nKids;++ii)  if (kids[ii]==kid)  return  true;
		//return false;
	//}
	//medialBall*& getKid(medialBall* kid)
	//{	for (int ii=0;ii<nKids;++ii)  if (kids[ii]==kid)  return kids[ii];
		//std::cout<<"no kid found! "<<" i"<<fi<<" j"<<fj<<" k"<<fk<<" "<<nKids<<std::endl;
		//return kids[1000];
	//}




	void addNei(medialBall* vj)   {
		if (nNeis>0)  {
			medialBall** ownNeis = neis;
			neis = new medialBall*[nNeis+1];
			int kk=-1;
			while(++kk<nNeis)	 neis[kk] = ownNeis[kk];
			delete[] ownNeis;
			neis[kk] = vj;
			++nNeis;
		}
		else  {
			neis    = new medialBall*[1];
			neis[0] = vj;
			nNeis   = 1;
		}
	}
public:

	float _0() const  { return fi; }
	float _1() const  { return fj; }
	float _2() const  { return fk; }
	dbl3 node() const { return dbl3(fi,fj,fk); }

	 voxel*          vxl;
	 float          fi, fj, fk;
	 short          type;
	 float          R;          ///< radius
	 unsigned short nKids;
	 unsigned short ikid;
	 unsigned short nNeis;
	 unsigned short mark;
	 unsigned short corId;
	 medialBall*    boss;
	 medialBall**   kids;
	 medialBall**   neis;


};

inline bool operator != ( const medialBall & a, const voxel& b) {

	return (short(a.fi) != b.i || short(a.fj) != b.j ||  short(a.fk) != b.k);  }

inline dbl3 operator-(const medialBall& a, const medialBall& b) {
	return dbl3(a.fi-b.fi,a.fj-b.fj,a.fk-b.fk);  }




class poreNE {
	poreNE(const poreNE&){};
 public:
	poreNE(): surfaceArea(0), volumn(0)
	{}
	dbl3 node() const     { return mb->node(); }
	int contact(int pid)  {
		auto it = contacts.find(pid);
		if (it!=contacts.end())			 return it->second;
		else { std::cout<<" E:noCntct ";  return -1;}
	}
	double radius() const { return  mb->R; }

 public:

	int surfaceArea;
	int volumn;

	std::map<int,int> contacts;
	const medialBall* mb;
};




class throatNE {
 public:
	throatNE(int trotid, int elm1, int elm2):  tid(trotid), e1(elm1),e2(elm2), surfaceArea(0), volumn(0), cachInd(0), nCrnrs(0),CrosArea(0.,0.,0.)
	{}
	int nToxel2Balls() const  {
		int nBs=0;  for(auto vi:toxels2){ if (vi->ball) ++nBs; }  return nBs; }

	double radius() const { return toxels1.empty() ? mb22()->R : 0.5*(mb22()->R+mb11()->R); }
	int neip(int i) const {return (&e1)[i];}

	dbl3 node() const { return mb22()->node(); }


 public:
	int tid;
	int e1;
	int e2;

	int surfaceArea;
	int volumn;

	unsigned short cachInd;

	short nCrnrs;

	//std::vector<vxlface>  vxfaces;
	std::vector<voxel*>  toxels2;
	std::vector<voxel*>  toxels1;
	const medialBall* mb22() const { return toxels2.front()->ball;}
	const medialBall* mb11() const { return toxels1.empty() ? nullptr :  toxels1.front()->ball; }
	medialBall* mb2Ch() {return toxels2.front()->ball;}
	medialBall* mb1Ch() {return toxels1.front()->ball;}


	dbl3 CrosArea;
	//int3x3 C;
};


inline double distSqr(const medialBall* i, const medialBall* j) {
	return (i->fi-j->fi)*(i->fi-j->fi)+(i->fj-j->fj)*(i->fj-j->fj)+(i->fk-j->fk)*(i->fk-j->fk); }

inline double distSqr(const medialBall* i, const voxel* j) {
	return (i->fi-j->i-_mp5)*(i->fi-j->i-_mp5)+(i->fj-j->j-_mp5)*(i->fj-j->j-_mp5)+(i->fk-j->k-_mp5)*(i->fk-j->k-_mp5);  }

inline double distSqr(const voxel& i, const voxel& j) {
 	return (i.i-j.i)*(i.i-j.i)+(i.j-j.j)*(i.j-j.j)+(i.k-j.k)*(i.k-j.k);  }


class metaballcomparer {
 public:
	bool operator()(const voxel* a, const voxel* b) const { return (a->R) > (b->R); }
};
inline double dist(const medialBall* i, const medialBall* j) {
	return sqrt(double((i->fi-j->fi)*(i->fi-j->fi)+(i->fj-j->fj)*(i->fj-j->fj)+(i->fk-j->fk)*(i->fk-j->fk)));
}
#endif

