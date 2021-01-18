#ifndef voxelImageManip_H
#define voxelImageManip_H

#include "blockNet.h"


//class mapComparer  {  public:
//bool operator() (pair<const int,short>& i1, pair<const int,short> i2) {return i1.second<i2.second;}  };

using namespace std; // std::pair, vector map

unsigned int growPores_X2(voxelField<int>&  VElems, int bgn, int lst, int porValue)
{
	voxelField<int> voxls = VElems;
	long long nChanges(0);

	OMPragma("omp parallel for reduction(+:nChanges)")
	for ( int k = 1; k<int(voxls.nz())-1 ; ++k )
	{
	  for ( int j = 1; j<int(voxls.ny())-1 ; ++j )
		for ( int i = 1; i<int(voxls.nx())-1 ; ++i )
		{
			const int* pijk = &voxls(i,j,k);
		  if (*pijk == porValue)
		  {
				   if ( bgn <= voxls.v_i(1,pijk) && voxls.v_i(1,pijk) <= lst )
				  { VElems(i,j,k) = voxls.v_i(1,pijk); ++nChanges; }
			  else if ( bgn <= voxls.v_i(-1,pijk) && voxls.v_i(-1,pijk) <= lst )
				  { VElems(i,j,k) = voxls.v_i(-1,pijk); ++nChanges; }
			  else if ( bgn <= voxls.v_j(1,pijk) && voxls.v_j(1,pijk) <= lst )
				  { VElems(i,j,k) = voxls.v_j(1,pijk); ++nChanges; }
			  else if ( bgn <= voxls.v_j(-1,pijk) && voxls.v_j(-1,pijk) <= lst )
				  { VElems(i,j,k) = voxls.v_j(-1,pijk); ++nChanges; }
			  else if ( bgn <= voxls.v_k(1,pijk) && voxls.v_k(1,pijk) <= lst )
				  { VElems(i,j,k) = voxls.v_k(1,pijk); ++nChanges; }
			  else if ( bgn <= voxls.v_k(-1,pijk) && voxls.v_k(-1,pijk) <= lst )
				  { VElems(i,j,k) = voxls.v_k(-1,pijk); ++nChanges; }
		  }
	  }
	}
	cout<<"  ngrowX3:"<<nChanges<<","; nChanges=0;

	OMPragma("omp parallel for reduction(+:nChanges)")
	for ( int k = 1; k<int(VElems.nz())-1 ; ++k )
	{
	  for ( int j = 1; j<int(VElems.ny())-1 ; ++j )
		for ( int i = 1; i<int(VElems.nx())-1 ; ++i )
		{
			const int* pijk = &VElems(i,j,k);
		  if (*pijk == porValue)
		  {
				   if ( bgn <= VElems.v_i(1,pijk) && VElems.v_i(1,pijk) <= lst )
				  { VElems(i,j,k) = VElems.v_i(1,pijk); ++nChanges; }
			  else if ( bgn <= VElems.v_i(-1,pijk) && VElems.v_i(-1,pijk) <= lst )
				  { VElems(i,j,k) = VElems.v_i(-1,pijk); ++nChanges; }
			  else if ( bgn <= VElems.v_j(1,pijk) && VElems.v_j(1,pijk) <= lst )
				  { VElems(i,j,k) = VElems.v_j(1,pijk); ++nChanges; }
			  else if ( bgn <= VElems.v_j(-1,pijk) && VElems.v_j(-1,pijk) <= lst )
				  { VElems(i,j,k) = VElems.v_j(-1,pijk); ++nChanges; }
			  else if ( bgn <= VElems.v_k(1,pijk) && VElems.v_k(1,pijk) <= lst )
				  { VElems(i,j,k) = VElems.v_k(1,pijk); ++nChanges; }
			  else if ( bgn <= VElems.v_k(-1,pijk) && VElems.v_k(-1,pijk) <= lst )
				  { VElems(i,j,k) = VElems.v_k(-1,pijk); ++nChanges; }
		  }
	  }
	}
	cout<<nChanges<<","; nChanges=0;

	OMPragma("omp parallel for reduction(+:nChanges)")
	for ( int k = int(VElems.nz())-2; k>=1 ; --k )
	{
	  for ( int j = int(VElems.ny())-2; j>=1 ; --j )
		for ( int i = int(VElems.nx())-2; i>=1 ; --i )
		{
			const int* pijk = &VElems(i,j,k);
		  if (*pijk == porValue)
		  {
				   if ( bgn <= VElems.v_i(1,pijk) && VElems.v_i(1,pijk) <= lst )
				  { VElems(i,j,k) = VElems.v_i(1,pijk); ++nChanges; }
			  else if ( bgn <= VElems.v_i(-1,pijk) && VElems.v_i(-1,pijk) <= lst )
				  { VElems(i,j,k) = VElems.v_i(-1,pijk); ++nChanges; }
			  else if ( bgn <= VElems.v_j(1,pijk) && VElems.v_j(1,pijk) <= lst )
				  { VElems(i,j,k) = VElems.v_j(1,pijk); ++nChanges; }
			  else if ( bgn <= VElems.v_j(-1,pijk) && VElems.v_j(-1,pijk) <= lst )
				  { VElems(i,j,k) = VElems.v_j(-1,pijk); ++nChanges; }
			  else if ( bgn <= VElems.v_k(1,pijk) && VElems.v_k(1,pijk) <= lst )
				  { VElems(i,j,k) = VElems.v_k(1,pijk); ++nChanges; }
			  else if ( bgn <= VElems.v_k(-1,pijk) && VElems.v_k(-1,pijk) <= lst )
				  { VElems(i,j,k) = VElems.v_k(-1,pijk); ++nChanges; }
		  }
	  }

	}
	cout<<"  ngrowX2:"<<nChanges<<"  ";
	return nChanges;
}




void growPores(voxelField<int>&  VElems, int bgn, int lst, int porValue)
{



	const voxelField<int> voxls = VElems;
	long long nChanges(0);

	OMPragma("omp parallel for reduction(+:nChanges)")
	for ( int k = 1; k<int(voxls.nz())-1 ; ++k )
	{
	  for ( int j = 1; j<int(voxls.ny())-1 ; ++j )
		for ( int i = 1; i<int(voxls.nx())-1 ; ++i )
		{
			 if (VElems(i,j,k) == porValue)
			 {
				const int* pijk = &voxls(i,j,k);
			  if ( bgn <= voxls.v_i(1,pijk) && voxls.v_i(1,pijk) <= lst )
				  { VElems(i,j,k) = voxls.v_i(1,pijk); ++nChanges; }
			  else if ( bgn <= voxls.v_i(-1,pijk) && voxls.v_i(-1,pijk) <= lst )
				  { VElems(i,j,k) = voxls.v_i(-1,pijk); ++nChanges; }
			  else if ( bgn <= voxls.v_j(1,pijk) && voxls.v_j(1,pijk) <= lst )
				  { VElems(i,j,k) = voxls.v_j(1,pijk); ++nChanges; }
			  else if ( bgn <= voxls.v_j(-1,pijk) && voxls.v_j(-1,pijk) <= lst )
				  { VElems(i,j,k) = voxls.v_j(-1,pijk); ++nChanges; }
			  else if ( bgn <= voxls.v_k(1,pijk) && voxls.v_k(1,pijk) <= lst )
				  { VElems(i,j,k) = voxls.v_k(1,pijk); ++nChanges; }
			  else if ( bgn <= voxls.v_k(-1,pijk) && voxls.v_k(-1,pijk) <= lst )
				  { VElems(i,j,k) = voxls.v_k(-1,pijk); ++nChanges; }

		 }
	  }
	}
	cout<<"  ngrowPors:"<<nChanges<<"  ";
}




void retreatPoresMedian(const inputDataNE & cg, voxelField<int>&  VElems, long bgn,  long lst,
		const vector<poreNE*>& poreIs, long unassigned)
{

	voxelField<int> voxls = VElems;
	long long nChanges(0);
	OMPragma("omp parallel for reduction(+:nChanges)")
	for (short k = 1; k <= cg.nz; ++k)
	{for (short j = 1; j <= cg.ny; ++j)
	 {const segments& s = cg.segs_[k-1][j-1];
	  for (short ix = 0; ix<s.cnt; ++ix)
	  {	for (short i = s.s[ix].start+1; i <= s.s[ix+1].start; ++i)
		{
			const int* pijk = &voxls(i,j,k);
		  long pID = *pijk;
		  short nSameID = 0;
		  short nDifferentID = 0;

		  if (pID>= bgn && pID <= lst)
		  {
			 if (voxls.v_i(-1,pijk) == pID) nSameID++;
			 else if (bgn <= voxls.v_i(-1,pijk) && lst>= voxls.v_i(-1,pijk))
				nDifferentID++;
			 if (voxls.v_i(1,pijk) == pID) nSameID++;
			 else if (bgn <= voxls.v_i(1,pijk) && lst>= voxls.v_i(1,pijk))
				nDifferentID++;
			 if (voxls.v_j(-1,pijk) == pID) nSameID++;
			 else if (bgn <= voxls.v_j(-1,pijk) && lst>= voxls.v_j(-1,pijk))
				nDifferentID++;
			 if (voxls.v_j(1,pijk) == pID) nSameID++;
			 else if (bgn <= voxls.v_j(1,pijk) && lst>= voxls.v_j(1,pijk))
				nDifferentID++;
			 if (voxls.v_k(-1,pijk) == pID) nSameID++;
			 else if (bgn <= voxls.v_k(-1,pijk) && lst>= voxls.v_k(-1,pijk))
				nDifferentID++;
			 if (voxls.v_k(1,pijk) == pID) nSameID++;
			 else if (bgn <= voxls.v_k(1,pijk) && lst>= voxls.v_k(1,pijk))
				nDifferentID++;

			if (nDifferentID > 0 && nSameID>0)
			{
				VElems(i,j,k) = unassigned;
				++nChanges;
			}
		  }

		}
	  }
	 }
	}

	(cout<<"  nRetreat:"<<nChanges).flush();

}


void growPoresMedStrict(const inputDataNE & cg, voxelField<int>&  VElems, long bgn,  long lst,
		const vector<poreNE*>& poreIs, long rawValue)
{

	voxelField<int> voxls = VElems;
	long long nChanges(0);
	OMPragma("omp parallel for reduction(+:nChanges)")
	for (short k = 1; k <= cg.nz; ++k)
	{for (short j = 1; j <= cg.ny; ++j)
	 {const segments& s = cg.segs_[k-1][j-1];
	  for (short ix = 0; ix<s.cnt; ++ix)
	  {	for (short i = s.s[ix].start+1; i <= s.s[ix+1].start; ++i)
		{
		  long pID = voxls(i,j,k);
		  const int* pijk = &voxls(i,j,k);

		  if (pID == rawValue)
		  {
			  float R = cg.segs_[k-1][j-1].vxl(i-1)->R;
			  short nDifferentID = 0;
			 if (bgn <= voxls.v_i(-1,pijk) && lst>= voxls.v_i(-1,pijk) && cg.segs_[k-1][j-1].vxl(i-2)->R >=R)
				nDifferentID++;
			 if (bgn <= voxls.v_i(1,pijk) && lst>= voxls.v_i(1,pijk) && cg.segs_[k-1][j-1].vxl(i)->R >= R)
				nDifferentID++;
			 if (bgn <= voxls.v_j(-1,pijk) && lst>= voxls.v_j(-1,pijk) && cg.segs_[k-1][j-2].vxl(i-1)->R >= R)
				nDifferentID++;
			 if (bgn <= voxls.v_j(1,pijk) && lst>= voxls.v_j(1,pijk) && cg.segs_[k-1][j].vxl(i-1)->R >= R)
				nDifferentID++;
			 if (bgn <= voxls.v_k(-1,pijk) && lst>= voxls.v_k(-1,pijk) && cg.segs_[k-2][j-1].vxl(i-1)->R >= R)
				nDifferentID++;
			 if (bgn <= voxls.v_k(1,pijk) && lst>= voxls.v_k(1,pijk) && cg.segs_[k][j-1].vxl(i-1)->R >= R)
				nDifferentID++;

			if (nDifferentID >= 3)
			{
			 map<int,short> neis;

			 long
			 neI = voxls.v_i(-1,pijk);
			 if (neI != pID && bgn <= neI && lst>= neI && cg.segs_[k-1][j-1].vxl(i-2)->R > R) 	 ++(neis.insert(pair<int,short>(neI,0)).first->second);
			 neI = voxls.v_i(1,pijk);
			 if (neI != pID && bgn <= neI && lst>= neI && cg.segs_[k-1][j-1].vxl(i)->R > R) 	 ++(neis.insert(pair<int,short>(neI,0)).first->second);
			 neI = voxls.v_j(-1,pijk);
			 if (neI != pID && bgn <= neI && lst>= neI && cg.segs_[k-1][j-2].vxl(i-1)->R > R) 	 ++(neis.insert(pair<int,short>(neI,0)).first->second);
			 neI = voxls.v_j(1,pijk);
			 if (neI != pID && bgn <= neI && lst>= neI && cg.segs_[k-1][j].vxl(i-1)->R > R) 	 ++(neis.insert(pair<int,short>(neI,0)).first->second);
			 neI = voxls.v_k(-1,pijk);
			 if (neI != pID && bgn <= neI && lst>= neI && cg.segs_[k-2][j-1].vxl(i-1)->R > R) 	 ++(neis.insert(pair<int,short>(neI,0)).first->second);
			 neI = voxls.v_k(1,pijk);
			 if (neI != pID && bgn <= neI && lst>= neI && cg.segs_[k][j-1].vxl(i-1)->R > R) 	 ++(neis.insert(pair<int,short>(neI,0)).first->second);

			 map<int,short>::iterator neitr = max_element(neis.begin(), neis.end(), mapComparer<int>());
			 if (neitr->second >= 3)
			 {
				++nChanges;
				VElems(i,j,k) = neitr->first;
			 }
		   }
		  }

		}
	  }
	 }
	}

	cout<<"  ngMedStrict:"<<nChanges<<" ";

}

void growPoresMedian(const inputDataNE & cg, voxelField<int>&  VElems, long bgn,  long lst,
		const vector<poreNE*>& poreIs, long rawValue)
{

	const voxelField<int> voxls = VElems;
	long long nChanges(0);
	OMPragma("omp parallel for reduction(+:nChanges)")
	for (short k = 1; k <= cg.nz; ++k)
	{for (short j = 1; j <= cg.ny; ++j)
	 {const segments& s = cg.segs_[k-1][j-1];
	  for (short ix = 0; ix<s.cnt; ++ix)
	  {	for (short i = s.s[ix].start+1; i <= s.s[ix+1].start; ++i)
		{
			//voxel* v=s.s[ix].v(i-1);
			//if (v)

			 const int* pijk = &voxls(i,j,k);
		  long pID = *pijk;

		  if (pID == rawValue)
		  {
			 float R = cg.segs_[k-1][j-1].vxl(i-1)->R;

			 short nDifferentID = 0;
			 if (bgn <= voxls.v_i(-1,pijk) && lst>= voxls.v_i(-1,pijk) && cg.segs_[k-1][j-1].vxl(i-2)->R >R)
				nDifferentID++;
			 if (bgn <= voxls.v_i(1,pijk) && lst>= voxls.v_i(1,pijk) && cg.segs_[k-1][j-1].vxl(i)->R > R)
				nDifferentID++;
			 if (bgn <= voxls.v_j(-1,pijk) && lst>= voxls.v_j(-1,pijk) && cg.segs_[k-1][j-2].vxl(i-1)->R > R)
				nDifferentID++;
			 if (bgn <= voxls.v_j(1,pijk) && lst>= voxls.v_j(1,pijk) && cg.segs_[k-1][j].vxl(i-1)->R > R)
				nDifferentID++;
			 if (bgn <= voxls.v_k(-1,pijk) && lst>= voxls.v_k(-1,pijk) && cg.segs_[k-2][j-1].vxl(i-1)->R > R)
				nDifferentID++;
			 if (bgn <= voxls.v_k(1,pijk) && lst>= voxls.v_k(1,pijk) && cg.segs_[k][j-1].vxl(i-1)->R > R)
				nDifferentID++;

			if (nDifferentID >= 2)
			{
			 map<int,short> neis;

			 long
			 neI = voxls.v_i(-1,pijk);
			 if (neI != pID && bgn <= neI && lst>= neI && cg.segs_[k-1][j-1].vxl(i-2)->R > R) 	 ++(neis.insert(pair<int,short>(neI,0)).first->second);
			 neI = voxls.v_i(1,pijk);
			 if (neI != pID && bgn <= neI && lst>= neI && cg.segs_[k-1][j-1].vxl(i)->R > R) 	 ++(neis.insert(pair<int,short>(neI,0)).first->second);
			 neI = voxls.v_j(-1,pijk);
			 if (neI != pID && bgn <= neI && lst>= neI && cg.segs_[k-1][j-2].vxl(i-1)->R > R) 	 ++(neis.insert(pair<int,short>(neI,0)).first->second);
			 neI = voxls.v_j(1,pijk);
			 if (neI != pID && bgn <= neI && lst>= neI && cg.segs_[k-1][j].vxl(i-1)->R > R) 	 ++(neis.insert(pair<int,short>(neI,0)).first->second);
			 neI = voxls.v_k(-1,pijk);
			 if (neI != pID && bgn <= neI && lst>= neI && cg.segs_[k-2][j-1].vxl(i-1)->R > R) 	 ++(neis.insert(pair<int,short>(neI,0)).first->second);
			 neI = voxls.v_k(1,pijk);
			 if (neI != pID && bgn <= neI && lst>= neI && cg.segs_[k][j-1].vxl(i-1)->R > R) 	 ++(neis.insert(pair<int,short>(neI,0)).first->second);

			 map<int,short>::iterator neitr = max_element(neis.begin(), neis.end(), mapComparer<int>());
			 if (neitr->second >= 2)
			 {
				++nChanges;
				VElems(i,j,k) = neitr->first;
			 }
		   }
		  }

		}
	  }
	 }
	}

	cout<<"  ngMedian:"<<nChanges<<" ";

}




void growPoresMedEqs(const inputDataNE & cg, voxelField<int>&  VElems, long bgn,  long lst,
		const vector<poreNE*>& poreIs, long rawValue)
{

	const voxelField<int> voxls = VElems;
	long long nChanges(0);
	OMPragma("omp parallel for reduction(+:nChanges)")
	for (short k = 1; k <= cg.nz; ++k)
	{for (short j = 1; j <= cg.ny; ++j)
	 {const segments& s = cg.segs_[k-1][j-1];
	  for (short ix = 0; ix<s.cnt; ++ix)
	  {	for (short i = s.s[ix].start+1; i <= s.s[ix+1].start; ++i)
		{
			//voxel* v=s.s[ix].v(i-1);
			//if (v)

			const int* pijk = &voxls(i,j,k);
		  long pID = *pijk;

		  if (pID == rawValue)
		  {
			 float R = cg.segs_[k-1][j-1].vxl(i-1)->R;

			 short nDifferentID = 0;
			 if (bgn <= voxls.v_i(-1,pijk) && lst>= voxls.v_i(-1,pijk) && cg.segs_[k-1][j-1].vxl(i-2)->R >= R)
				nDifferentID++;
			 if (bgn <= voxls.v_i(1,pijk) && lst>= voxls.v_i(1,pijk) && cg.segs_[k-1][j-1].vxl(i)->R >= R)
				nDifferentID++;
			 if (bgn <= voxls.v_j(-1,pijk) && lst>= voxls.v_j(-1,pijk) && cg.segs_[k-1][j-2].vxl(i-1)->R >= R)
				nDifferentID++;
			 if (bgn <= voxls.v_j(1,pijk) && lst>= voxls.v_j(1,pijk) && cg.segs_[k-1][j].vxl(i-1)->R >= R)
				nDifferentID++;
			 if (bgn <= voxls.v_k(-1,pijk) && lst>= voxls.v_k(-1,pijk) && cg.segs_[k-2][j-1].vxl(i-1)->R >= R)
				nDifferentID++;
			 if (bgn <= voxls.v_k(1,pijk) && lst>= voxls.v_k(1,pijk) && cg.segs_[k][j-1].vxl(i-1)->R >= R)
				nDifferentID++;

			if (nDifferentID >= 2)
			{
			 map<int,short> neis;

			 long
			 neI = voxls.v_i(-1,pijk);
			 if (neI != pID && bgn <= neI && lst>= neI && cg.segs_[k-1][j-1].vxl(i-2)->R >= R) 	 ++(neis.insert(pair<int,short>(neI,0)).first->second);
			 neI = voxls.v_i(1,pijk);
			 if (neI != pID && bgn <= neI && lst>= neI && cg.segs_[k-1][j-1].vxl(i)->R >= R) 	 ++(neis.insert(pair<int,short>(neI,0)).first->second);
			 neI = voxls.v_j(-1,pijk);
			 if (neI != pID && bgn <= neI && lst>= neI && cg.segs_[k-1][j-2].vxl(i-1)->R >= R) 	 ++(neis.insert(pair<int,short>(neI,0)).first->second);
			 neI = voxls.v_j(1,pijk);
			 if (neI != pID && bgn <= neI && lst>= neI && cg.segs_[k-1][j].vxl(i-1)->R >= R) 	 ++(neis.insert(pair<int,short>(neI,0)).first->second);
			 neI = voxls.v_k(-1,pijk);
			 if (neI != pID && bgn <= neI && lst>= neI && cg.segs_[k-2][j-1].vxl(i-1)->R >= R) 	 ++(neis.insert(pair<int,short>(neI,0)).first->second);
			 neI = voxls.v_k(1,pijk);
			 if (neI != pID && bgn <= neI && lst>= neI && cg.segs_[k][j-1].vxl(i-1)->R >= R) 	 ++(neis.insert(pair<int,short>(neI,0)).first->second);

			 map<int,short>::iterator neitr = max_element(neis.begin(), neis.end(), mapComparer<int>());
			 if (neitr->second >= 2)
			 {
				++nChanges;
				VElems(i,j,k) = neitr->first;
			 }
		   }
		  }

		}
	  }
	 }
	}

	cout<<"  ngMedEqs:"<<nChanges<<"  ";

}



void growPoresMedEqsLoose(const inputDataNE & cg, voxelField<int>&  VElems, long bgn,  long lst,
		const vector<poreNE*>& poreIs, long rawValue)
{

	voxelField<int> voxls = VElems;
	long long nChanges(0);
	OMPragma("omp parallel for reduction(+:nChanges)")
	for (short k = 1; k <= cg.nz; ++k)
	{for (short j = 1; j <= cg.ny; ++j)
	 {const segments& s = cg.segs_[k-1][j-1];
	  for (short ix = 0; ix<s.cnt; ++ix)
	  {	for (short i = s.s[ix].start+1; i <= s.s[ix+1].start; ++i)
		{
			//voxel* v=s.s[ix].v(i-1);
			//if (v)

			const int* pijk = &voxls(i,j,k);
		  long pID = *pijk;

		  if (pID == rawValue)
		  {
//			 float R = cg.segs_[k-1][j-1].vxl(i-1)->R;

			 short nDifferentID = 0;
			 if (bgn <= voxls.v_i(-1,pijk) && lst>= voxls.v_i(-1,pijk))
				nDifferentID++;
			 if (bgn <= voxls.v_i(1,pijk) && lst>= voxls.v_i(1,pijk))
				nDifferentID++;
			 if (bgn <= voxls.v_j(-1,pijk) && lst>= voxls.v_j(-1,pijk))
				nDifferentID++;
			 if (bgn <= voxls.v_j(1,pijk) && lst>= voxls.v_j(1,pijk))
				nDifferentID++;
			 if (bgn <= voxls.v_k(-1,pijk) && lst>= voxls.v_k(-1,pijk))
				nDifferentID++;
			 if (bgn <= voxls.v_k(1,pijk) && lst>= voxls.v_k(1,pijk))
				nDifferentID++;

			if (nDifferentID >= 2)
			{
			 map<int,short> neis;

			 long
			 neI = voxls.v_i(-1,pijk);	 if (neI != pID && bgn <= neI && lst>= neI) 	 ++(neis.insert(pair<int,short>(neI,0)).first->second);
			 neI = voxls.v_i(1,pijk);	 if (neI != pID && bgn <= neI && lst>= neI) 	 ++(neis.insert(pair<int,short>(neI,0)).first->second);
			 neI = voxls.v_j(-1,pijk);	 if (neI != pID && bgn <= neI && lst>= neI) 	 ++(neis.insert(pair<int,short>(neI,0)).first->second);
			 neI = voxls.v_j(1,pijk);	 if (neI != pID && bgn <= neI && lst>= neI) 	 ++(neis.insert(pair<int,short>(neI,0)).first->second);
			 neI = voxls.v_k(-1,pijk);	 if (neI != pID && bgn <= neI && lst>= neI) 	 ++(neis.insert(pair<int,short>(neI,0)).first->second);
			 neI = voxls.v_k(1,pijk);	 if (neI != pID && bgn <= neI && lst>= neI) 	 ++(neis.insert(pair<int,short>(neI,0)).first->second);

			 map<int,short>::iterator neitr = max_element(neis.begin(), neis.end(), mapComparer<int>());
			 if (neitr->second >= 2)
			 {
				++nChanges;
				VElems(i,j,k) = neitr->first;
			 }
		   }
		  }

		}
	  }
	 }
	}

	cout<<"  ngMedLoose:"<<nChanges<<"  ";

}











void medianElem(const inputDataNE & cg, voxelField<int>&  VElems, long bgn,  long lst,
		const vector<poreNE*>& poreIs)
{
	voxelField<int> voxls = VElems;
	long long nChanges(0);
	OMPragma("omp parallel for reduction(+:nChanges)")
	for (short k = 1; k <= cg.nz; ++k)
	{for (short j = 1; j <= cg.ny; ++j)
	 {const segments& s = cg.segs_[k-1][j-1];
	  for (short ix = 0; ix<s.cnt; ++ix)
	  {	for (short i = s.s[ix].start+1; i <= s.s[ix+1].start; ++i)
		{
			const int* pijk = &voxls(i,j,k);
		  long pID = *pijk;
		  short nSameID = 0;
		  short nDifferentID = 0;

		  if (pID>= bgn && pID <= lst)
		  {
			 if (voxls.v_i(-1,pijk) == pID) nSameID++;
			 else if (bgn <= voxls.v_i(-1,pijk) && lst>= voxls.v_i(-1,pijk))
				nDifferentID++;
			 if (voxls.v_i(1,pijk) == pID) nSameID++;
			 else if (bgn <= voxls.v_i(1,pijk) && lst>= voxls.v_i(1,pijk))
				nDifferentID++;
			 if (voxls.v_j(-1,pijk) == pID) nSameID++;
			 else if (bgn <= voxls.v_j(-1,pijk) && lst>= voxls.v_j(-1,pijk))
				nDifferentID++;
			 if (voxls.v_j(1,pijk) == pID) nSameID++;
			 else if (bgn <= voxls.v_j(1,pijk) && lst>= voxls.v_j(1,pijk))
				nDifferentID++;
			 if (voxls.v_k(-1,pijk) == pID) nSameID++;
			 else if (bgn <= voxls.v_k(-1,pijk) && lst>= voxls.v_k(-1,pijk))
				nDifferentID++;
			 if (voxls.v_k(1,pijk) == pID) nSameID++;
			 else if (bgn <= voxls.v_k(1,pijk) && lst>= voxls.v_k(1,pijk))
				nDifferentID++;

			if (nDifferentID > nSameID)
			{
			 map<int,short> neis;


			 long
			 neI = voxls.v_i(-1,pijk);	 if (neI != pID && bgn <= neI && lst>= neI) 	 ++(neis.insert(pair<int,short>(neI,0)).first->second);
			 neI = voxls.v_i(1,pijk);	 if (neI != pID && bgn <= neI && lst>= neI) 	 ++(neis.insert(pair<int,short>(neI,0)).first->second);
			 neI = voxls.v_j(-1,pijk);	 if (neI != pID && bgn <= neI && lst>= neI) 	 ++(neis.insert(pair<int,short>(neI,0)).first->second);
			 neI = voxls.v_j(1,pijk);	 if (neI != pID && bgn <= neI && lst>= neI) 	 ++(neis.insert(pair<int,short>(neI,0)).first->second);
			 neI = voxls.v_k(-1,pijk);	 if (neI != pID && bgn <= neI && lst>= neI) 	 ++(neis.insert(pair<int,short>(neI,0)).first->second);
			 neI = voxls.v_k(1,pijk);	 if (neI != pID && bgn <= neI && lst>= neI) 	 ++(neis.insert(pair<int,short>(neI,0)).first->second);
			 for (map<int,short>::iterator neitr = neis.begin();neitr != neis.end();++neitr)
			 {
				if (neitr->second > nSameID)
				{
					++nChanges;
					VElems(i,j,k) = neitr->first;
					nSameID = neitr->second;
				}
			 }
			}
		  }

		}
	  }
	 }
	}

	cout<<"  nMedian:"<<nChanges<<" ";

}






#endif
