#ifndef ELEMENT_H
#define ELEMENT_H


#include "vec3.h"
#include <map>
#include <iostream>
#include <array>

const int LEVEL_MAX = 32767;
#define nAprox 1





class medialBall;
class voxel
{/// voxel 
public:
	voxel(const voxel& v) :  ball(v.ball),i(v.i), j(v.j), k(v.k),R(v.R) {};//,hasBall(false)
	voxel() : ball(NULL),i(-1), j(-1), k(-1),R(0) {};
public:
	medialBall* ball;
	short i, j, k;
	float R; /// placed here for memory and performance reasons 
};

class vxlface
{
public:
	vxlface(short i, short j, short k, char dirr ) :  _ijk({{i,j,k}}),dir(dirr) {};//,hasBall(false)
	vxlface() : _ijk({{-1,-1,-1}}),dir(0) {};
	const array<short,3>& ijk() {return _ijk;};
	array<short,3> ijkp() {array<short,3> ijkp1 = _ijk; ijkp1[dir]+=1; return ijkp1;};
private:
	array<short,3> _ijk;
	char dir;
};

///  segment is a list of voxels along a x axis of the image, which have the same value,
///  it used to save memory and computational time.
class segment
{
public:
	voxel* v(int i){return segV+(i-start);};
	int start;
	unsigned char value;///. don't know why the order is important when including c headers !
	voxel* segV;
};

/// #segments is the list of all "#segment"s in a x axis in the image
class segments
{
public:
 segments(): s(NULL),cnt(0){};
 void reSize(int size){ if (s == NULL) {s = new segment[size];} else {std::cout<<" \nError in segments "<<size_t(s)<<std::endl; s = new segment[size];} }

 ~segments(){ if (s != NULL && cnt) {delete[] s;s = NULL;} }

 voxel* vxl(int i) const
 {  for (int p=0; p<cnt; ++p) { if ( (i>=s[p].start) && (i<s[p+1].start) )	return s[p].v(i); } 	return NULL; }

 segment* s;
 int cnt;
};




/// used to correct for half-voxel distacne
#define _pp5 0.0
#define _mp5 -0.5

/// maximal-sphere class
class medialBall
{
  public:

	medialBall(short t) : vxl(NULL), fi(-10000), fj(-0.5), fk(-10000),type(t), R(-10000), nKids(0), nNeis(0), corId(0), boss(this), kids(NULL), neis(NULL){};
	medialBall(voxel* v, short t) : vxl(v), fi(v->i-_mp5), fj(v->j-_mp5), fk(v->k-_mp5), type(t), R(v->R), nKids(0), nNeis(0), corId(0), boss(this), kids(NULL), neis(NULL){};
	~medialBall(){ if(kids) delete[] kids; if(neis) delete[] neis;};


	short level() const
	{	if (this == boss)	    {return 1;}//counter=0;
		else  if (boss != NULL)	return boss->level()+1;
		else {std::cout<<"F";std::cout.flush(); return -10000;}
	}
	bool inParents(medialBall* vj) const
	{	if (vj == boss)			{return true;}//counter=0;
		else  if (boss == this)	return false;
		else  if (boss != NULL)	return boss->inParents(vj);
		else {std::cout<<"F";std::cout.flush(); return false;}
	}

	medialBall* mastrVxl() const
	{	if (this == boss) 	return boss;
		else 				return boss->mastrVxl();
	};

	inline const medialBall* smallestP() const
	{	register const medialBall* pvMin = this;
		register const medialBall* v = this;
		while (v != v->boss)
		{	pvMin = v->R < pvMin->R ? v : pvMin;   v=v->boss;  }
		return pvMin;
	};

	bool isNei(const medialBall* neib) const
	{	for (int ii=0;ii<nNeis;++ii)  if (neis[ii]==neib)	return true;
		return false;
	};
	bool neisNotNeiWith(const medialBall* neib) const
	{	for (int ii=0;ii<nNeis;++ii)  if (neis[ii]->isNei(neib))  return false;
		return true;
	};
	medialBall*& neiRef(const medialBall* neib)
	{	for (int ii=0;ii<nNeis;++ii)  if (neis[ii]==neib)  return neis[ii];
		std::cout<<"no nei found! "<<int(nNeis)<<std::endl;
		return neis[1000];
	};
	void removeKidBoss(medialBall* kid)
	{
        if (nKids>1)
        {
            int ii=-1;
            for (int kk=0;kk<nKids;++kk)    if(kids[kk]!=kid)   kids[++ii]=kids[kk];
            nKids=ii+1;
        }
        else if (nKids==1 && kids[0]==kid) {delete[] kids; kids=NULL; nKids=0;}
        kid->boss = boss;
	};
	bool isKid(const medialBall* kid) const
	{	for (int ii=0;ii<nKids;++ii)  if (kids[ii]==kid)  return  true;
		return false;
	};
	medialBall*& getKid(medialBall* kid)
	{	for (int ii=0;ii<nKids;++ii)  if (kids[ii]==kid)  return kids[ii];
		std::cout<<"no kid found! "<<" i"<<fi<<" j"<<fj<<" k"<<fk<<" "<<nKids<<std::endl;
		return kids[1000];
	};




   void addNei(medialBall* vj)
   {
		if (nNeis>0)
		{
			medialBall** ownNeis = neis;
			neis=new medialBall*[nNeis+1];
			int kk=-1;
			while(++kk<nNeis)		neis[kk]=ownNeis[kk];
			delete[] ownNeis;
			neis[kk]=vj;
			++nNeis;
		}
		else
		{
			neis=new medialBall*[1];
			neis[0]=vj;
			nNeis=1;
		}
	}
public:

	 voxel* vxl; 
	 float fi,fj,fk;
	 short type;
	 float R; /// const but gcc-4.7 won't like it
	 unsigned short nKids, ikid; ///. ikid is improving computational time only
	 unsigned short nNeis;
	 unsigned short mark, corId;
	 medialBall* boss;
	 medialBall** kids;
	 medialBall** neis;


};

inline bool operator != ( const medialBall & a, const voxel& b)
{
	return (short(a.fi+_pp5) != b.i || short(a.fj+_pp5) != b.j ||  short(a.fk+_pp5) != b.k);
}

inline vec3 operator-(const medialBall& a, const medialBall& b)
{	return vec3(a.fi-b.fi,a.fj-b.fj,a.fk-b.fk);}











/// base element class for pores and throats, not necassary but anyway.
class ElementO
{
public:
 ElementO() : surfaceArea(0), volumn(0), nSec(3), cachInd(0) {}
 virtual ~ElementO() {}
 virtual bool ispore() = 0;
 int surfaceArea;
 int volumn;

 medialBall* balls;

	int subid;                  /// corner id
	static int nTroatSubids;    /// to be added to subid to get subid2. each throat is assumed to have 2 halves

	const unsigned short nSec;  /// number of corners
	unsigned short cachInd;     /// used for counting

};




class poreElementI: public ElementO
{
 poreElementI(const poreElementI&){};
public:
 poreElementI() :   vMin(255), vMax(0) {};
 virtual ~poreElementI() {};
 bool ispore() {return true;};
 int contact(int pid)
 { 	  std::map<int,int>::iterator it = contacts.find(pid);
	  if (it!=contacts.end())			 return it->second;
	  else { std::cout<<" E:noCntct ";  return -1;}
}
 double radius() const {return  mb->R;};

public:

 std::map<int,int> contacts;
 const medialBall* mb;
 unsigned char vMin,vMax;
 double vDist[nAprox],vAvg[nAprox];
};








class throatElementI : public ElementO
{
private:
 throatElementI(){};
public:
 throatElementI(int trotid, int elm1, int elm2): tid(trotid), e1(elm1),e2(elm2),nCrnrs(0),CrosArea(0.0,0.0,0.0),C{{ {{0,0,0}},{{0,0,0}},{{0,0,0}}}}, vMin(255), vMax(0) {}//mb2(NULL),mb1(NULL),
 virtual ~throatElementI() {}
 bool ispore() {return false;}
 int nToxel2Balls() const
 {	int nBs=0;
	for (std::vector<voxel*>::const_iterator  vitr = toxels2.begin(); vitr<toxels2.end();++vitr) if ((*vitr)->ball) ++nBs;
	return nBs;
 }
 double radius() const { return toxels1.empty() ? mb22()->R : 0.5*(mb22()->R+mb11()->R); }
 int neip(int i) const {return (&e1)[i];}
public:
 int tid; //. tid needed when sorting
 int e1; /// first pore index
 int e2; /// second pore index

 short nCrnrs;

 std::vector<vxlface>  faces;
 std::vector<voxel*>  toxels2;
 std::vector<voxel*>  toxels1;
 const medialBall* mb22() const {return (*toxels2.begin())->ball;}
 const medialBall* mb11() const  { return toxels1.empty() ? NULL :  (*toxels1.begin())->ball; }
 medialBall* mb2Ch() {return (*toxels2.begin())->ball;}
 medialBall* mb1Ch() {return (*toxels1.begin())->ball;}


 vec3 CrosArea;
 int3x3 C;
 unsigned char vMin,vMax;
 double vDist[nAprox],vAvg[nAprox];
};


inline double distSqr(const medialBall* i, const medialBall* j)
{	return (i->fi-j->fi)*(i->fi-j->fi)+(i->fj-j->fj)*(i->fj-j->fj)+(i->fk-j->fk)*(i->fk-j->fk); }

inline double distSqr(const medialBall* i, const voxel* j)
{	return (i->fi-j->i-_mp5)*(i->fi-j->i-_mp5)+(i->fj-j->j-_mp5)*(i->fj-j->j-_mp5)+(i->fk-j->k-_mp5)*(i->fk-j->k-_mp5);  }

inline double distSqr(const voxel& i, const voxel& j)
{ 	return (i.i-j.i)*(i.i-j.i)+(i.j-j.j)*(i.j-j.j)+(i.k-j.k)*(i.k-j.k);  }


class metaballcomparer
{ public:
	bool operator()(const voxel* a, const voxel* b) const { return (a->R) > (b->R); }
};
inline double dist(const medialBall* i, const medialBall* j) // distance between two voxels
{
	return sqrt(double((i->fi-j->fi)*(i->fi-j->fi)+(i->fj-j->fj)*(i->fj-j->fj)+(i->fk-j->fk)*(i->fk-j->fk)));
}
#endif

