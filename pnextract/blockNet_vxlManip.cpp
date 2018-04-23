#ifndef voxelImageManip_H
#define voxelImageManip_H



//class mapComparer  {  public:
//bool operator() (std::pair<const int,short>& i1, std::pair<const int,short> i2) {return i1.second<i2.second;}  };


unsigned int growPores_X2(voxelField<int>&  VElems, int min, int max, int porValue)
	{
	voxelField<int> voxls = VElems;
	register long long nChanges(0);

	for ( int k = 1; k<int(voxls.size3()[2])-1 ; k++ )
	{
	  for ( int j = 1; j<int(voxls.size3()[1])-1 ; j++ )
	  {
			for ( int i = 1; i<int(voxls.size3()[0])-1 ; i++ )
			{
					const int* pijk = &voxls(i,j,k);
			  if (*pijk == porValue)
			  {
						 if ( min <= voxls.v_i(1,pijk) && voxls.v_i(1,pijk) <= max )
					  { VElems(i,j,k) = voxls.v_i(1,pijk); ++nChanges; }
				  else if ( min <= voxls.v_i(-1,pijk) && voxls.v_i(-1,pijk) <= max )
					  { VElems(i,j,k) = voxls.v_i(-1,pijk); ++nChanges; }
				  else if ( min <= voxls.v_j(1,pijk) && voxls.v_j(1,pijk) <= max )
					  { VElems(i,j,k) = voxls.v_j(1,pijk); ++nChanges; }
				  else if ( min <= voxls.v_j(-1,pijk) && voxls.v_j(-1,pijk) <= max )
					  { VElems(i,j,k) = voxls.v_j(-1,pijk); ++nChanges; }
				  else if ( min <= voxls.v_k(1,pijk) && voxls.v_k(1,pijk) <= max )
					  { VElems(i,j,k) = voxls.v_k(1,pijk); ++nChanges; }
				  else if ( min <= voxls.v_k(-1,pijk) && voxls.v_k(-1,pijk) <= max )
					  { VElems(i,j,k) = voxls.v_k(-1,pijk); ++nChanges; }

				  if (VElems(i,j,k) != porValue)
				  {
						if ( VElems(i+1,j,k) == porValue )
					{ VElems(i+1,j,k) = VElems(i,j,k); ++nChanges; }
				 else if ( VElems(i-1,j,k) == porValue )
					{ VElems(i-1,j,k) = VElems(i,j,k); ++nChanges; }
				 else if ( VElems(i,j+1,k) == porValue )
					{ VElems(i,j+1,k) = VElems(i,j,k); ++nChanges; }
				 else if ( VElems(i,j-1,k) == porValue )
					{ VElems(i,j-1,k) = VElems(i,j,k); ++nChanges; }
				 else if ( VElems(i,j,k+1) == porValue )
					{ VElems(i,j,k+1) = VElems(i,j,k); ++nChanges; }
				 else if ( VElems(i,j,k-1) == porValue )
					{ VElems(i,j,k-1) = VElems(i,j,k); ++nChanges; }
			 }
		  }
		  }
	  }
	}
	cout<<"  ngrowX2:"<<nChanges<<"  ";
	return nChanges;
}




void growPores(voxelField<int>&  VElems, int min, int max, int porValue)
{

	const voxelField<int> voxls = VElems;
	register long long nChanges(0);

	for ( int k = 1; k<int(voxls.size3()[2])-1 ; k++ )
	{
	  for ( int j = 1; j<int(voxls.size3()[1])-1 ; j++ )
	  {
			for ( int i = 1; i<int(voxls.size3()[0])-1 ; i++ )
			{
				 if (VElems(i,j,k) == porValue)
				 {
					const int* pijk = &voxls(i,j,k);
				  if ( min <= voxls.v_i(1,pijk) && voxls.v_i(1,pijk) <= max )
					  { VElems(i,j,k) = voxls.v_i(1,pijk); ++nChanges; }
				  else if ( min <= voxls.v_i(-1,pijk) && voxls.v_i(-1,pijk) <= max )
					  { VElems(i,j,k) = voxls.v_i(-1,pijk); ++nChanges; }
				  else if ( min <= voxls.v_j(1,pijk) && voxls.v_j(1,pijk) <= max )
					  { VElems(i,j,k) = voxls.v_j(1,pijk); ++nChanges; }
				  else if ( min <= voxls.v_j(-1,pijk) && voxls.v_j(-1,pijk) <= max )
					  { VElems(i,j,k) = voxls.v_j(-1,pijk); ++nChanges; }
				  else if ( min <= voxls.v_k(1,pijk) && voxls.v_k(1,pijk) <= max )
					  { VElems(i,j,k) = voxls.v_k(1,pijk); ++nChanges; }
				  else if ( min <= voxls.v_k(-1,pijk) && voxls.v_k(-1,pijk) <= max )
					  { VElems(i,j,k) = voxls.v_k(-1,pijk); ++nChanges; }

			 }
		  }
	  }
	}
	cout<<"  ngrowPors:"<<nChanges<<"  ";
}




void retreatPoresMedian(const inputDataNE & cg, voxelField<int>&  VElems, long min,  long max,
		const std::vector<poreElementI*>& poreIs, long rawValue)
{
	voxelField<int> voxls = VElems;
	register long long nChanges(0);
	for (short k = 1; k <= cg.nz; ++k)
	{for (short j = 1; j <= cg.ny; ++j)
	 {const segments& s = cg.segs_[k-1][j-1];
	  for (short ix = 0; ix<s.cnt; ++ix)
	  {	for (short i = s.s[ix].start+1; i <= s.s[ix+1].start; ++i)
		{
			const int* pijk = &voxls(i,j,k);
		  register long pID = *pijk;
		  register short nSameID = 0;
		  register short nDifferentID = 0;

		  if (pID>= min && pID <= max)
		  {
			 if (voxls.v_i(-1,pijk) == pID) nSameID++;
			 else if (min <= voxls.v_i(-1,pijk) && max>= voxls.v_i(-1,pijk))
				nDifferentID++;
			 if (voxls.v_i(1,pijk) == pID) nSameID++;
			 else if (min <= voxls.v_i(1,pijk) && max>= voxls.v_i(1,pijk))
				nDifferentID++;
			 if (voxls.v_j(-1,pijk) == pID) nSameID++;
			 else if (min <= voxls.v_j(-1,pijk) && max>= voxls.v_j(-1,pijk))
				nDifferentID++;
			 if (voxls.v_j(1,pijk) == pID) nSameID++;
			 else if (min <= voxls.v_j(1,pijk) && max>= voxls.v_j(1,pijk))
				nDifferentID++;
			 if (voxls.v_k(-1,pijk) == pID) nSameID++;
			 else if (min <= voxls.v_k(-1,pijk) && max>= voxls.v_k(-1,pijk))
				nDifferentID++;
			 if (voxls.v_k(1,pijk) == pID) nSameID++;
			 else if (min <= voxls.v_k(1,pijk) && max>= voxls.v_k(1,pijk))
				nDifferentID++;

			if (nDifferentID > 0 && nSameID>0)
			{
				VElems(i,j,k) = rawValue;
				++nChanges;
			}
		  }

		}
	  }
	 }
	}

	(cout<<"  nRetreat:"<<nChanges).flush();

}


void growPoresMedStrict(const inputDataNE & cg, voxelField<int>&  VElems, long min,  long max,
		const std::vector<poreElementI*>& poreIs, long rawValue)
{
	voxelField<int> voxls = VElems;
	register long long nChanges(0);
	for (short k = 1; k <= cg.nz; ++k)
	{for (short j = 1; j <= cg.ny; ++j)
	 {const segments& s = cg.segs_[k-1][j-1];
	  for (short ix = 0; ix<s.cnt; ++ix)
	  {	for (short i = s.s[ix].start+1; i <= s.s[ix+1].start; ++i)
		{
		  register long pID = voxls(i,j,k);
		  const int* pijk = &voxls(i,j,k);

		  if (pID == rawValue)
		  {
			  float R = cg.segs_[k-1][j-1].vxl(i-1)->R;
			  register short nDifferentID = 0;
			 if (min <= voxls.v_i(-1,pijk) && max>= voxls.v_i(-1,pijk) && cg.segs_[k-1][j-1].vxl(i-2)->R >=R)
				nDifferentID++;
			 if (min <= voxls.v_i(1,pijk) && max>= voxls.v_i(1,pijk) && cg.segs_[k-1][j-1].vxl(i)->R >= R)
				nDifferentID++;
			 if (min <= voxls.v_j(-1,pijk) && max>= voxls.v_j(-1,pijk) && cg.segs_[k-1][j-2].vxl(i-1)->R >= R)
				nDifferentID++;
			 if (min <= voxls.v_j(1,pijk) && max>= voxls.v_j(1,pijk) && cg.segs_[k-1][j].vxl(i-1)->R >= R)
				nDifferentID++;
			 if (min <= voxls.v_k(-1,pijk) && max>= voxls.v_k(-1,pijk) && cg.segs_[k-2][j-1].vxl(i-1)->R >= R)
				nDifferentID++;
			 if (min <= voxls.v_k(1,pijk) && max>= voxls.v_k(1,pijk) && cg.segs_[k][j-1].vxl(i-1)->R >= R)
				nDifferentID++;

			if (nDifferentID >= 3)
			{
			 std::map<int,short> neis;///.  ID-counter

			 register long
			 neiPID = voxls.v_i(-1,pijk);
			 if (neiPID != pID && min <= neiPID && max>= neiPID && cg.segs_[k-1][j-1].vxl(i-2)->R > R) 	 ++(neis.insert(std::pair<int,short>(neiPID,0)).first->second);
			 neiPID = voxls.v_i(1,pijk);
			 if (neiPID != pID && min <= neiPID && max>= neiPID && cg.segs_[k-1][j-1].vxl(i)->R > R) 	 ++(neis.insert(std::pair<int,short>(neiPID,0)).first->second);
			 neiPID = voxls.v_j(-1,pijk);
			 if (neiPID != pID && min <= neiPID && max>= neiPID && cg.segs_[k-1][j-2].vxl(i-1)->R > R) 	 ++(neis.insert(std::pair<int,short>(neiPID,0)).first->second);
			 neiPID = voxls.v_j(1,pijk);
			 if (neiPID != pID && min <= neiPID && max>= neiPID && cg.segs_[k-1][j].vxl(i-1)->R > R) 	 ++(neis.insert(std::pair<int,short>(neiPID,0)).first->second);
			 neiPID = voxls.v_k(-1,pijk);
			 if (neiPID != pID && min <= neiPID && max>= neiPID && cg.segs_[k-2][j-1].vxl(i-1)->R > R) 	 ++(neis.insert(std::pair<int,short>(neiPID,0)).first->second);
			 neiPID = voxls.v_k(1,pijk);
			 if (neiPID != pID && min <= neiPID && max>= neiPID && cg.segs_[k][j-1].vxl(i-1)->R > R) 	 ++(neis.insert(std::pair<int,short>(neiPID,0)).first->second);

			 std::map<int,short>::iterator neitr = std::max_element(neis.begin(), neis.end(), mapComparer<int>());
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

void growPoresMedian(const inputDataNE & cg, voxelField<int>&  VElems, long min,  long max,
		const std::vector<poreElementI*>& poreIs, long rawValue)
{
	const voxelField<int> voxls = VElems;
	register long long nChanges(0);
	for (short k = 1; k <= cg.nz; ++k)
	{for (short j = 1; j <= cg.ny; ++j)
	 {const segments& s = cg.segs_[k-1][j-1];
	  for (short ix = 0; ix<s.cnt; ++ix)
	  {	for (short i = s.s[ix].start+1; i <= s.s[ix+1].start; ++i)
		{
			//~ voxel* v=s.s[ix].v(i-1);
			//~ if (v)

			 const int* pijk = &voxls(i,j,k);
		  register long pID = *pijk;

		  if (pID == rawValue)
		  {
			 float R = cg.segs_[k-1][j-1].vxl(i-1)->R;

			 register short nDifferentID = 0;
			 if (min <= voxls.v_i(-1,pijk) && max>= voxls.v_i(-1,pijk) && cg.segs_[k-1][j-1].vxl(i-2)->R >R)
				nDifferentID++;
			 if (min <= voxls.v_i(1,pijk) && max>= voxls.v_i(1,pijk) && cg.segs_[k-1][j-1].vxl(i)->R > R)
				nDifferentID++;
			 if (min <= voxls.v_j(-1,pijk) && max>= voxls.v_j(-1,pijk) && cg.segs_[k-1][j-2].vxl(i-1)->R > R)
				nDifferentID++;
			 if (min <= voxls.v_j(1,pijk) && max>= voxls.v_j(1,pijk) && cg.segs_[k-1][j].vxl(i-1)->R > R)
				nDifferentID++;
			 if (min <= voxls.v_k(-1,pijk) && max>= voxls.v_k(-1,pijk) && cg.segs_[k-2][j-1].vxl(i-1)->R > R)
				nDifferentID++;
			 if (min <= voxls.v_k(1,pijk) && max>= voxls.v_k(1,pijk) && cg.segs_[k][j-1].vxl(i-1)->R > R)
				nDifferentID++;

			if (nDifferentID >= 2)
			{
			 std::map<int,short> neis;///.  ID-counter

			 register long
			 neiPID = voxls.v_i(-1,pijk);
			 if (neiPID != pID && min <= neiPID && max>= neiPID && cg.segs_[k-1][j-1].vxl(i-2)->R > R) 	 ++(neis.insert(std::pair<int,short>(neiPID,0)).first->second);
			 neiPID = voxls.v_i(1,pijk);
			 if (neiPID != pID && min <= neiPID && max>= neiPID && cg.segs_[k-1][j-1].vxl(i)->R > R) 	 ++(neis.insert(std::pair<int,short>(neiPID,0)).first->second);
			 neiPID = voxls.v_j(-1,pijk);
			 if (neiPID != pID && min <= neiPID && max>= neiPID && cg.segs_[k-1][j-2].vxl(i-1)->R > R) 	 ++(neis.insert(std::pair<int,short>(neiPID,0)).first->second);
			 neiPID = voxls.v_j(1,pijk);
			 if (neiPID != pID && min <= neiPID && max>= neiPID && cg.segs_[k-1][j].vxl(i-1)->R > R) 	 ++(neis.insert(std::pair<int,short>(neiPID,0)).first->second);
			 neiPID = voxls.v_k(-1,pijk);
			 if (neiPID != pID && min <= neiPID && max>= neiPID && cg.segs_[k-2][j-1].vxl(i-1)->R > R) 	 ++(neis.insert(std::pair<int,short>(neiPID,0)).first->second);
			 neiPID = voxls.v_k(1,pijk);
			 if (neiPID != pID && min <= neiPID && max>= neiPID && cg.segs_[k][j-1].vxl(i-1)->R > R) 	 ++(neis.insert(std::pair<int,short>(neiPID,0)).first->second);

			 std::map<int,short>::iterator neitr = std::max_element(neis.begin(), neis.end(), mapComparer<int>());
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




void growPoresMedEqs(const inputDataNE & cg, voxelField<int>&  VElems, long min,  long max,
		const std::vector<poreElementI*>& poreIs, long rawValue)
{
	const voxelField<int> voxls = VElems;
	register long long nChanges(0);
	for (short k = 1; k <= cg.nz; ++k)
	{for (short j = 1; j <= cg.ny; ++j)
	 {const segments& s = cg.segs_[k-1][j-1];
	  for (short ix = 0; ix<s.cnt; ++ix)
	  {	for (short i = s.s[ix].start+1; i <= s.s[ix+1].start; ++i)
		{
			//~ voxel* v=s.s[ix].v(i-1);
			//~ if (v)

			const int* pijk = &voxls(i,j,k);
		  register long pID = *pijk;

		  if (pID == rawValue)
		  {
			 float R = cg.segs_[k-1][j-1].vxl(i-1)->R;

			 register short nDifferentID = 0;
			 if (min <= voxls.v_i(-1,pijk) && max>= voxls.v_i(-1,pijk) && cg.segs_[k-1][j-1].vxl(i-2)->R >= R)
				nDifferentID++;
			 if (min <= voxls.v_i(1,pijk) && max>= voxls.v_i(1,pijk) && cg.segs_[k-1][j-1].vxl(i)->R >= R)
				nDifferentID++;
			 if (min <= voxls.v_j(-1,pijk) && max>= voxls.v_j(-1,pijk) && cg.segs_[k-1][j-2].vxl(i-1)->R >= R)
				nDifferentID++;
			 if (min <= voxls.v_j(1,pijk) && max>= voxls.v_j(1,pijk) && cg.segs_[k-1][j].vxl(i-1)->R >= R)
				nDifferentID++;
			 if (min <= voxls.v_k(-1,pijk) && max>= voxls.v_k(-1,pijk) && cg.segs_[k-2][j-1].vxl(i-1)->R >= R)
				nDifferentID++;
			 if (min <= voxls.v_k(1,pijk) && max>= voxls.v_k(1,pijk) && cg.segs_[k][j-1].vxl(i-1)->R >= R)
				nDifferentID++;

			if (nDifferentID >= 2)
			{
			 std::map<int,short> neis;///.  ID-counter

			 register long
			 neiPID = voxls.v_i(-1,pijk);
			 if (neiPID != pID && min <= neiPID && max>= neiPID && cg.segs_[k-1][j-1].vxl(i-2)->R >= R) 	 ++(neis.insert(std::pair<int,short>(neiPID,0)).first->second);
			 neiPID = voxls.v_i(1,pijk);
			 if (neiPID != pID && min <= neiPID && max>= neiPID && cg.segs_[k-1][j-1].vxl(i)->R >= R) 	 ++(neis.insert(std::pair<int,short>(neiPID,0)).first->second);
			 neiPID = voxls.v_j(-1,pijk);
			 if (neiPID != pID && min <= neiPID && max>= neiPID && cg.segs_[k-1][j-2].vxl(i-1)->R >= R) 	 ++(neis.insert(std::pair<int,short>(neiPID,0)).first->second);
			 neiPID = voxls.v_j(1,pijk);
			 if (neiPID != pID && min <= neiPID && max>= neiPID && cg.segs_[k-1][j].vxl(i-1)->R >= R) 	 ++(neis.insert(std::pair<int,short>(neiPID,0)).first->second);
			 neiPID = voxls.v_k(-1,pijk);
			 if (neiPID != pID && min <= neiPID && max>= neiPID && cg.segs_[k-2][j-1].vxl(i-1)->R >= R) 	 ++(neis.insert(std::pair<int,short>(neiPID,0)).first->second);
			 neiPID = voxls.v_k(1,pijk);
			 if (neiPID != pID && min <= neiPID && max>= neiPID && cg.segs_[k][j-1].vxl(i-1)->R >= R) 	 ++(neis.insert(std::pair<int,short>(neiPID,0)).first->second);

			 std::map<int,short>::iterator neitr = std::max_element(neis.begin(), neis.end(), mapComparer<int>());
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



void growPoresMedEqsLoose(const inputDataNE & cg, voxelField<int>&  VElems, long min,  long max,
		const std::vector<poreElementI*>& poreIs, long rawValue)
{
	voxelField<int> voxls = VElems;
	register long long nChanges(0);
	for (short k = 1; k <= cg.nz; ++k)
	{for (short j = 1; j <= cg.ny; ++j)
	 {const segments& s = cg.segs_[k-1][j-1];
	  for (short ix = 0; ix<s.cnt; ++ix)
	  {	for (short i = s.s[ix].start+1; i <= s.s[ix+1].start; ++i)
		{
			//~ voxel* v=s.s[ix].v(i-1);
			//~ if (v)

			const int* pijk = &voxls(i,j,k);
		  register long pID = *pijk;

		  if (pID == rawValue)
		  {
//			 float R = cg.segs_[k-1][j-1].vxl(i-1)->R;

			 register short nDifferentID = 0;
			 if (min <= voxls.v_i(-1,pijk) && max>= voxls.v_i(-1,pijk))
				nDifferentID++;
			 if (min <= voxls.v_i(1,pijk) && max>= voxls.v_i(1,pijk))
				nDifferentID++;
			 if (min <= voxls.v_j(-1,pijk) && max>= voxls.v_j(-1,pijk))
				nDifferentID++;
			 if (min <= voxls.v_j(1,pijk) && max>= voxls.v_j(1,pijk))
				nDifferentID++;
			 if (min <= voxls.v_k(-1,pijk) && max>= voxls.v_k(-1,pijk))
				nDifferentID++;
			 if (min <= voxls.v_k(1,pijk) && max>= voxls.v_k(1,pijk))
				nDifferentID++;

			if (nDifferentID >= 2)
			{
			 std::map<int,short> neis;///.  ID-counter

			 register long
			 neiPID = voxls.v_i(-1,pijk);
			 if (neiPID != pID && min <= neiPID && max>= neiPID) 	 ++(neis.insert(std::pair<int,short>(neiPID,0)).first->second);
			 neiPID = voxls.v_i(1,pijk);
			 if (neiPID != pID && min <= neiPID && max>= neiPID) 	 ++(neis.insert(std::pair<int,short>(neiPID,0)).first->second);
			 neiPID = voxls.v_j(-1,pijk);
			 if (neiPID != pID && min <= neiPID && max>= neiPID) 	 ++(neis.insert(std::pair<int,short>(neiPID,0)).first->second);
			 neiPID = voxls.v_j(1,pijk);
			 if (neiPID != pID && min <= neiPID && max>= neiPID) 	 ++(neis.insert(std::pair<int,short>(neiPID,0)).first->second);
			 neiPID = voxls.v_k(-1,pijk);
			 if (neiPID != pID && min <= neiPID && max>= neiPID) 	 ++(neis.insert(std::pair<int,short>(neiPID,0)).first->second);
			 neiPID = voxls.v_k(1,pijk);
			 if (neiPID != pID && min <= neiPID && max>= neiPID) 	 ++(neis.insert(std::pair<int,short>(neiPID,0)).first->second);

			 std::map<int,short>::iterator neitr = std::max_element(neis.begin(), neis.end(), mapComparer<int>());
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











void medianElem(const inputDataNE & cg, voxelField<int>&  VElems, long min,  long max,
		const std::vector<poreElementI*>& poreIs)
{
	voxelField<int> voxls = VElems;
	register long long nChanges(0);
	for (short k = 1; k <= cg.nz; ++k)
	{for (short j = 1; j <= cg.ny; ++j)
	 {const segments& s = cg.segs_[k-1][j-1];
	  for (short ix = 0; ix<s.cnt; ++ix)
	  {	for (short i = s.s[ix].start+1; i <= s.s[ix+1].start; ++i)
		{
			const int* pijk = &voxls(i,j,k);
		  register long pID = *pijk;
		  register short nSameID = 0;
		  register short nDifferentID = 0;

		  if (pID>= min && pID <= max)
		  {
			 if (voxls.v_i(-1,pijk) == pID) nSameID++;
			 else if (min <= voxls.v_i(-1,pijk) && max>= voxls.v_i(-1,pijk))
				nDifferentID++;
			 if (voxls.v_i(1,pijk) == pID) nSameID++;
			 else if (min <= voxls.v_i(1,pijk) && max>= voxls.v_i(1,pijk))
				nDifferentID++;
			 if (voxls.v_j(-1,pijk) == pID) nSameID++;
			 else if (min <= voxls.v_j(-1,pijk) && max>= voxls.v_j(-1,pijk))
				nDifferentID++;
			 if (voxls.v_j(1,pijk) == pID) nSameID++;
			 else if (min <= voxls.v_j(1,pijk) && max>= voxls.v_j(1,pijk))
				nDifferentID++;
			 if (voxls.v_k(-1,pijk) == pID) nSameID++;
			 else if (min <= voxls.v_k(-1,pijk) && max>= voxls.v_k(-1,pijk))
				nDifferentID++;
			 if (voxls.v_k(1,pijk) == pID) nSameID++;
			 else if (min <= voxls.v_k(1,pijk) && max>= voxls.v_k(1,pijk))
				nDifferentID++;

			if (nDifferentID > nSameID)
			{
			 std::map<int,short> neis;///.  ID-counter

			 register long
			 neiPID = voxls.v_i(-1,pijk);
			 if (neiPID != pID && min <= neiPID && max>= neiPID) 	 ++(neis.insert(std::pair<int,short>(neiPID,0)).first->second);
			 neiPID = voxls.v_i(1,pijk);
			 if (neiPID != pID && min <= neiPID && max>= neiPID) 	 ++(neis.insert(std::pair<int,short>(neiPID,0)).first->second);
			 neiPID = voxls.v_j(-1,pijk);
			 if (neiPID != pID && min <= neiPID && max>= neiPID) 	 ++(neis.insert(std::pair<int,short>(neiPID,0)).first->second);
			 neiPID = voxls.v_j(1,pijk);
			 if (neiPID != pID && min <= neiPID && max>= neiPID) 	 ++(neis.insert(std::pair<int,short>(neiPID,0)).first->second);
			 neiPID = voxls.v_k(-1,pijk);
			 if (neiPID != pID && min <= neiPID && max>= neiPID) 	 ++(neis.insert(std::pair<int,short>(neiPID,0)).first->second);
			 neiPID = voxls.v_k(1,pijk);
			 if (neiPID != pID && min <= neiPID && max>= neiPID) 	 ++(neis.insert(std::pair<int,short>(neiPID,0)).first->second);
			 for (std::map<int,short>::iterator neitr = neis.begin();neitr != neis.end();++neitr)
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
