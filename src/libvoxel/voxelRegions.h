/*-------------------------------------------------------------------------*\

This file is part of voxelImage library, a C++ template library  
developed by Ali Qaseminejad Raeini for handelling 3D raw images.


Please see our website for relavant literature making use of this code:
https://www.imperial.ac.uk/earth-science/research/research-groups/pore-scale-modelling/

For further information please contact us by email:
Ali Q Raeini: a.q.raeini@imperial.ac.uk

\*-------------------------------------------------------------------------*/


#ifndef voxelRegions_H
#define voxelRegions_H

using namespace std;

#define min7Nei(_vxls, vp_) \
		min(*vp_,min(min( \
		min(_vxls.v_i(-1,vp_),_vxls.v_i( 1,vp_)),\
		min(_vxls.v_j(-1,vp_),_vxls.v_j( 1,vp_))),\
		min(_vxls.v_k(-1,vp_),_vxls.v_k( 1,vp_))))
#define min5Nei(_vxls, vp_) \
		min(*vp_, min( \
		min(_vxls.v_i(-1,vp_),_vxls.v_i( 1,vp_)),\
		min(_vxls.v_j(-1,vp_),_vxls.v_j( 1,vp_))))


inline int recurseMin(const vector<int>& lblMap, int lbv)
{
	while(lblMap[lbv]<lbv) lbv=lblMap[lbv];
	return lbv;
}


template<typename T>
voxelImageT<int> labelImage(const voxelImageT<T>& vImage, const T minvv=0, const T maxvv=std::numeric_limits<T>::max())  {//  TODO to be tested

	//voxels.write("dump1.mhd");
	int3 n_2 = vImage.size3();
	static const int bigN=255*255*255*126;

	cout<<" labeling Image between "<<int(minvv)<<"  "<<int(maxvv)<<endl;

	voxelImageT<int> lbls(n_2[0]+2,n_2[1]+2,n_2[2]+2, bigN); //+2 is needed for min5Nei
	lbls.X0Ch()=vImage.X0();
	lbls.dxCh()=vImage.dx();
	forAllkji_(vImage) if(minvv<=vImage(i,j,k) && vImage(i,j,k)<=maxvv) lbls(i+1,j+1,k+1)= i+j*n_2[0]+2; ///. int won't work for large images, go slice by slice

	///. reduce size in each ***slice**
	forAllkvp_(lbls)  if(*vp<bigN) *vp=min5Nei(lbls,vp);
	forAllk_vp_(lbls)  if(*vp<bigN) *vp=min5Nei(lbls,vp);
	forAllkvp_(lbls)  if(*vp<bigN) *vp=min5Nei(lbls,vp);
	forAllk_vp_(lbls)  if(*vp<bigN) *vp=min5Nei(lbls,vp);
	forAllkvp_(lbls)  if(*vp<bigN) *vp=min5Nei(lbls,vp);
	forAllk_vp_(lbls)  if(*vp<bigN) *vp=min5Nei(lbls,vp);
	forAllkvp_(lbls)  if(*vp<bigN) *vp=min5Nei(lbls,vp);
	forAllk_vp_(lbls)  if(*vp<bigN) *vp=min5Nei(lbls,vp);
	forAllkvp_(lbls)  if(*vp<bigN) *vp=min5Nei(lbls,vp);
	forAllk_vp_(lbls)  if(*vp<bigN) *vp=min5Nei(lbls,vp);
	forAllkvp_(lbls)  if(*vp<bigN) *vp=min5Nei(lbls,vp);


	{ // accumulate slice labels 2020-06-20
		forAllkvp_(lbls)  if(*vp<bigN) *vp=min(*vp,*(vp-1));
		int oldvv=-1, il=0;
		forAllkji(lbls)   //RLE compress
		{	int vv = lbls(i,j,k);
			if(vv<bigN) { if(vv!=oldvv) { oldvv=lbls(i,j,k); lbls(i,j,k) = ++il; } else lbls(i,j,k) = il; }
			else oldvv=-1;
		} 
		cout<<"max_il: "<<il<<endl;
	}

	//OMPragma("omp parallel for")
	//for (int k=1; k<(lbls).nz()-1; ++k)
	//{	int mxl=lbls.nij_+1;
		//vector<int> lblMap(lbls.nij_,mxl);
		//vector<int> lblMapp(lbls.nij_,mxl);
		//for(auto* vp=&lbls(0,0,k), *_ve_=&lbls(0,0,k+1); vp<_ve_; ++vp)
			 //if(*vp<bigN) { int lmin=recurseMin(lblMap,min5Nei(lbls,vp)); lblMap[*vp] = lmin;  *vp=lmin; }

		//int nReg=1;
		//for(auto& lbl: lblMap) if(lbl<mxl) { lbl=recurseMin(lblMap,lbl); if(lblMapp[lbl]>nReg) lblMapp[lbl]=++nReg; }

		//for(auto* vp=&lbls(0,0,k), *_ve_=&lbls(0,0,k+1); vp<_ve_; ++vp)
			 //if(*vp<mxl) *vp=lblMapp[lblMap[*vp]];
	//}

	///////. restore global lbl ind
	//vector<int> max_kLbls(n_2[2]+2,0);
	//forAllkvp_(lbls) if(*vp<bigN) { *vp=min5Nei(lbls,vp); max_kLbls[k]=max(max_kLbls[k],*vp); }
	//{int sumLbl(0);	for(auto& lblb: max_kLbls) {sumLbl+=lblb+2; lblb=sumLbl; }}
	//forAllkvp_(lbls) if(*vp<bigN) (*vp)+=max_kLbls[k];
	//cout<<"max_kLbls: "<<max_kLbls<<endl;

	//lbls.write("lll3.raw");

	long long sumlbl = 0, sumlblOld = -1; int itr=0;
	while(sumlbl!=sumlblOld)
	{	sumlblOld=sumlbl;

		if(++itr==1)
		{
		lbls.rotate('z');
		forAllvp_(lbls)  if(*vp<bigN) *vp=min7Nei(lbls,vp);
		forAll_vp_(lbls)  if(*vp<bigN) *vp=min7Nei(lbls,vp);
		forAllvp_(lbls)  if(*vp<bigN) *vp=min7Nei(lbls,vp);
		forAll_vp_(lbls)  if(*vp<bigN) *vp=min7Nei(lbls,vp);
		forAllvp_(lbls)  if(*vp<bigN) *vp=min7Nei(lbls,vp);
		forAll_vp_(lbls)  if(*vp<bigN) *vp=min7Nei(lbls,vp);
		forAllvp_(lbls)  if(*vp<bigN) *vp=min7Nei(lbls,vp);
		forAll_vp_(lbls)  if(*vp<bigN) *vp=min7Nei(lbls,vp);
		forAllvp_(lbls)  if(*vp<bigN) *vp=min7Nei(lbls,vp);
		forAll_vp_(lbls)  if(*vp<bigN) *vp=min7Nei(lbls,vp);
		forAllvp_(lbls)  if(*vp<bigN) *vp=min7Nei(lbls,vp);
		forAll_vp_(lbls)  if(*vp<bigN) *vp=min7Nei(lbls,vp);
		forAllvp_(lbls)  if(*vp<bigN) *vp=min7Nei(lbls,vp);
		forAll_vp_(lbls)  if(*vp<bigN) *vp=min7Nei(lbls,vp);
		forAllvp_(lbls)  if(*vp<bigN) *vp=min7Nei(lbls,vp);
		forAll_vp_(lbls)  if(*vp<bigN) *vp=min7Nei(lbls,vp);
		forAllvp_(lbls)  if(*vp<bigN) *vp=min7Nei(lbls,vp);
		forAll_vp_(lbls)  if(*vp<bigN) *vp=min7Nei(lbls,vp);
		forAllvp_(lbls)  if(*vp<bigN) *vp=min7Nei(lbls,vp);
		forAll_vp_(lbls)  if(*vp<bigN) *vp=min7Nei(lbls,vp);
		forAllvp_(lbls)  if(*vp<bigN) *vp=min7Nei(lbls,vp);
		forAll_vp_(lbls)  if(*vp<bigN) *vp=min7Nei(lbls,vp);
		forAllvp_(lbls)  if(*vp<bigN) *vp=min7Nei(lbls,vp);
		forAll_vp_(lbls)  if(*vp<bigN) *vp=min7Nei(lbls,vp);
		forAllvp_(lbls)  if(*vp<bigN) *vp=min7Nei(lbls,vp);
		forAllvp_(lbls)  if(*vp<bigN) *vp=min7Nei(lbls,vp); //! increase these if code doesn't work!!!!
		lbls.rotate('z');
		//lbls.write("lll4.raw");
		}
		forAll_vp_(lbls)  if(*vp<bigN) *vp=min7Nei(lbls,vp);
		forAllvp_(lbls)  if(*vp<bigN) *vp=min7Nei(lbls,vp);
		forAll_vp_(lbls)  if(*vp<bigN) *vp=min7Nei(lbls,vp);
		forAllvp_(lbls)  if(*vp<bigN) *vp=min7Nei(lbls,vp);


		///. compress
		int mxl=0; sumlbl=0;
		OMPragma("omp parallel for reduction(max:mxl) reduction(+:sumlbl)")
		forAllcp(lbls) if((*cp)<bigN) { mxl = max(mxl,*cp);  sumlbl+=*cp; }
		++mxl;
		vector<int> lblMap(mxl+3,mxl);// size should be larger than value?
		vector<int> lblMapp(mxl+3,mxl);
		cout<<" mxl:"<<mxl<<" sumlbl:"<<sumlbl<<endl;
		forAllvp_seq(lbls)
		 if(*vp<bigN)  { int lmin=recurseMin(lblMap,min7Nei(lbls,vp));   lblMap[*vp]=lmin; *vp=lmin; }
		for(auto& lbl: lblMap) {if(lbl<mxl) lbl=recurseMin(lblMap,lbl); }
		for(auto& lbl: lblMap) {if(lbl<mxl) lbl=recurseMin(lblMap,lbl); }
		for(auto& lbl: lblMap) {if(lbl<mxl) lbl=recurseMin(lblMap,lbl); }

		int nReg=1;
		for(auto& lbl: lblMap) if(lbl<mxl) { lbl=recurseMin(lblMap,lbl); if(lblMapp[lbl]>nReg) lblMapp[lbl]=++nReg; }

		forAllvp_(lbls)  if(*vp<mxl) *vp=lblMapp[lblMap[*vp]];
	}




	forAllvp_(lbls)  if(*vp>=bigN) *vp=-1;

	lbls.cropD(int3(1,1,1),int3(n_2[0]+1,n_2[1]+1,n_2[2]+1),0,0,false);
	return lbls;

}


template<typename T>
void keepLargest0(voxelImageT<T>& vImage, const T minvv=0, const T maxvv=0)  {
	const voxelImageT<int> lbls = labelImage(vImage,minvv,maxvv);
	const auto* lbl0=&lbls(0);
	int mxl = 0;
	OMPragma("omp parallel for reduction(max:mxl)")
	forAllcp(lbls) if(*cp>=0) mxl = max(mxl,*cp);
	vector<long long> lblsN(mxl+1,0);
	forAllcp_seq_(lbls)  if(*cp>=0)  { const T vv=vImage(cp-lbl0);   if(minvv<=vv&& vv<=maxvv) ++lblsN[*cp]; }

	int lrglbl = distance(lblsN.begin(), max_element(lblsN.begin(), lblsN.end()));
	//lbls.write("dumplbls.raw");
	int3 n = vImage.size3();
	cout<<" connected fraction: "<<double(lblsN[lrglbl])/(n[0]*n[1]*n[2])<<endl;
	forAllvp_(lbls) if(*vp!=lrglbl)	{ const T vv=vImage(vp-lbl0);   if(minvv<=vv&& vv<=maxvv) vImage(vp-lbl0) = maxT(T)-1;	}//! isolated=254
}


#endif
