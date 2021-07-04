
#include "inputData.h"
#include "medialSurf.h"
#include "typses.h"




#include "medialRadius.cpp"

medialSurface::medialSurface(inputDataNE& cfg)//, double vmvLimRelF, double crossAreaf
:	cg_(cfg), segs_(cfg.segs_), ToBeAssigned(0)  {
	setDefaults(-5.); 	//set _minRp to negative value if not provided by user, to be re-assigned in calc_distmaps(),  to be synced with setDefaults()

	nx = cfg.nx;
	ny = cfg.ny;
	nz = cfg.nz;

	size_t nvoxls=0; // local copy for omp
	OMPragma("omp parallel for reduction(+:nvoxls)")
	for (short k = 0; k < nz; ++k)
	 for (short j = 0; j < ny; ++j)  {	const segments& s = cg_.segs_[k][j];
		for (short ix = 0; ix<s.cnt; ++ix)
			if (s.s[ix].value == 0) nvoxls+=s.s[ix+1].start-s.s[ix].start;
	 }
	nVxls = nvoxls;
	invalidSeg.start=-10000;
	invalidSeg.value=255;

}

void medialSurface::setDefaults(double avgR)  {
	/// minRPore/Rnoise is a different keyword, of its own, not part of medialSurfaceSettings.    It is supposed to be the only adjustable parameter for the regular users.       The default value of minRPore is 1.75.

	/// medialSurfaceSettings are for advanced users and its defaults are chosen (if not provided by the user) based on the minRPore value.


	/// To extract a network that has a lower network coordination number, you can decrease the values of lenNf (e.g. 0.4), and vmvRadRelNf (e.g. 1.05) and increase the values of nRSmoothing (e.g. 9), RCorsf (e.g. 0.2) and RCors (e.g. 2.5).  You can also consider increasing minRPore (also named Rnoise) keyword to let say 2..     Every change you make you need to check that the network produces reasonable results as these are sensitive parameters and do not behave linearly. We should not be woried about the high coordination number as long as the network predicts the physical properties correctly, but not everybody agrees with me here!


	_minRp=min(1.25, avgR*0.25)+0.5;
	if (cg_.giv("Rnoise"+_s(0), _minRp) || cg_.giv("minRPore", _minRp) || cg_.giv("Rnoise", _minRp)) cout<< " minimum pore radius: " << _minRp <<endl;
	else  cout<<" keyword \"minRPore\" not found, default value ("<<abs(_minRp)<<") will be used"<<endl;

	_clipROutx=0.05;
	_clipROutyz=0.98;
	_midRf=0.7;
	_MSNoise=1.*abs(_minRp)+1.;
	_lenNf=0.6;
	_vmvRadRelNf=1.1;
	_nRSmoothing=3;
	_RCorsnf=0.15;
	_RCorsn=abs(_minRp);



	if(cg_.nBP6==6)	 _clipROutyz=_clipROutx;

    std::istringstream keywrdData;
    if (cg_.giv("medialSurfaceSettings"+_s(0), keywrdData) || cg_.giv("medialSurfaceSettings", keywrdData))  {
	 	keywrdData  >>_clipROutx >>_clipROutyz  >>_midRf >>_MSNoise  >>_lenNf >>_vmvRadRelNf >>_nRSmoothing >>_RCorsnf >>_RCorsn;
	}

	if(_minRp<0.) cout<<" Default setting, will be updated after distance map computation:\n";
	cout<<"  minRPore     : "<< abs(_minRp)<<";\n";
	cout<<"  medialSurfaceSettings: " << _clipROutx<<"  "<<_clipROutyz<<"  "<<_midRf<<"  "<<_MSNoise<<"  "<< _lenNf<<"  "<<_vmvRadRelNf<<"  "<<_nRSmoothing<<"  "<<_RCorsnf<<"  " << _RCorsn<<endl;

	cout<<"  medialSurfaceSettings:\n"
		<<"   clipROutx     : "<< _clipROutx<<"\n"
		<<"   clipROutyz    : "<< _clipROutyz<<"\n"
		<<"   midRFrac      : "<< _midRf<<"\n"
		<<"   RMedSurfNoise : "<< _MSNoise<<"\n"
		<<"   lenNf         : "<< _lenNf<<"\n"
		<<"   vmvRadRelNf   : "<< _vmvRadRelNf<<"\n"
		<<"   nRSmoothing   : "<< _nRSmoothing<<"\n"
		<<"   RCorsf  : "<< _RCorsnf<<"\n"
		<<"   RCors   : "<< _RCorsn<<"\n"
		<<endl;
}


void medialSurface::buildvoxelspace()  { ///  Build voxelspace -- memory for void/active voxels
	cout<<"\nProcessing "<<cg_._rockTypes[0].name <<" voxels:"<<endl;
	cout<<" Creating "<<nVxls<<" voxels with index: "<< int(0);cout.flush();


	vxlSpace.resize(nVxls);

	std::vector<voxel>::iterator p = vxlSpace.begin();
	for (int iz = 0; iz<nz; ++iz)  {
		for (int iy = 0; iy<ny; ++iy)  {
			const segments& s = segs_[iz][iy];
			for (int ix = 0; ix<s.cnt; ++ix)  {  if (s.s[ix].value == 0)
				 for (int i = s.s[ix].start; i<s.s[ix+1].start; ++i)  {
					p->i = i;
					p->j = iy;
					p->k = iz;
					++p;
				 }
			}

		}
	}

	if ( nVxls != size_t(p-vxlSpace.begin()) ) cout<<"\n Error created "<<size_t(p-vxlSpace.begin())<<" voxels "<<endl;

	cout<<endl;

	/// Link voxels to segments
	p = vxlSpace.begin();
	for (int iz = 0; iz<nz; ++iz)
		for (int iy = 0; iy<ny; ++iy)  {	segments& s = segs_[iz][iy];
			for (int ix = 0; ix<s.cnt; ++ix)
				if (s.s[ix].value == 0)  {
					s.s[ix].segV = &*p;
					p += s.s[ix+1].start - s.s[ix].start;
				}
		}
}


void medialSurface::paradox_pre_removeincludedballI() //to remove the included maximal bals
{ /// Remove maximal-balls, leave one in each adjacent voxesl. This saves time when sorting in paradoxremoveincludedballI()
	if (!nVxls) { return; }

	cout<< " pre-remove included balls: out of " <<vxlSpace.size(); cout.flush();

	int ndel = 0;

	for (int kk = 0; kk<nz; kk += 2)  { for (int jj = 0; jj<ny; jj += 2)
	  { const segments& s = segs_[kk][jj];
		for (int ix = 0; ix<s.cnt; ++ix)  { if (s.s[ix].value == 0)
		  { for (int ii = s.s[ix].start; ii<s.s[ix+1].start; ii += 2)  {
				voxel* smallers[8] = {nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr};
				int counter = -1;
				float maxRRR = 0;
				voxel* maxRPV = nullptr;


				for (int c = 0; c<2 ; ++c)
				 for (int b = 0; b<2 ; ++b)
					for (int a = 0; a<2; ++a)  {
						voxel* vi = vxl(ii+a, jj+b, kk+c);
						if (vi != nullptr && vi->ball==&ToBeAssigned )  {
							if (vi->R>maxRRR)  {
								if (maxRPV)	smallers[++counter] = maxRPV;
								maxRRR = vi->R;
								maxRPV = vi;
							}
							else
							{
								smallers[++counter] = vi;
							}
						}
					}
				++counter;
				ndel += counter;
				while(counter> 0)	if (smallers[--counter]) {  smallers[counter]->ball=nullptr; }

			}
		  }
		}
	  }
	}

	nBalls -= ndel;
	cout<< ",   removed = " << ndel << " remained = " << nBalls <<endl;

}


void medialSurface::paradoxremoveincludedballI()
{/// Remove included balls.  What remains are called maximal spheres
	if (!nVxls) { return; }


	cout<< " sorting... ";
	std::vector<voxel* >tvs; tvs.reserve(nBalls);//allocate memory for void voxels
	{
		std::vector<voxel>::iterator ti = vxlSpace.begin()-1;
		std::vector<voxel>::iterator tend = vxlSpace.end();
		while (++ti<tend) if (ti->ball) tvs.push_back(&(*(ti)));
	}
	cout<< tvs.size() <<" balls"<<endl;
	sort(tvs.begin(), tvs.end(), metaballcomparer());

	cout<< " remove included balls:"; cout.flush();

	int ndel = 0;
	auto vpp = tvs.begin(),   end = tvs.end();
	while (vpp < end)	{

		voxel* vi = *vpp;
		if (!vi->ball)	{	++vpp;	continue;	}

		const int   x = vi->i;
		const int   y = vi->j;
		const int   z = vi->k;
		const float ri=vi->R;
		const float ripinc = ri+0.55;//.+RPreDelete
		const float mbmbDist = _RCorsnf*ri+_RCorsn;

		int ex, ey, ez;
		ex = ripinc;
		for (int a = -ex; a <=  ex; ++a)  {
			ey = std::sqrt(ripinc*ripinc-a*a);
			for (int b = -ey; b <=  ey; ++b)  {
				ez = sqrt(ripinc*ripinc-a*a-b*b);//sqrts(r2i)+1-a-b;
				for (int c = -ez; c <=  ez; ++c)  {
					voxel* vj = vxl(x+a, y+b, z+c);
					if ((vj != nullptr) && (vj->ball) && (vj != vi))
					{
						const float rj = vj->R;
						if(rj <= ri)  {
							float D = sqrtf(a*a+b*b+c*c);

							if ( D < mbmbDist || (D+rj<ripinc+_MSNoise) )
							{
								vj->ball = nullptr;
								++ndel;
							}
						}
					}
				}
			}
		}

		++vpp;
		if ((vpp-tvs.begin())%10000 == 0)	cout<< "\r  remove = " << ndel;
	}
	cout<< "\r  removed = " << ndel << " remained = " << tvs.size()-ndel << " balls"<<endl;
	nBalls -= ndel;
 	return ;
}



void medialSurface::moveUphill(medialBall* b_i) // const
{	/// Refines the maximal-sphere location and radius,





	
	const voxel* vi = vxl(b_i->fi, b_i->fj, b_i->fk);
	dbl3 disp(0.,0.,0.);
	{
		const voxel* vjm = vxl(vi->i-1, vi->j, vi->k);
		const voxel* vjp = vxl(vi->i+1, vi->j, vi->k);
		if (vjm && vjp)  {
			float gp = vjp->R-vi->R;
			float gm = vi->R-vjm->R;
			if(abs(gp-gm)>0.01)  disp.x = max(-0.49, min(0.49, -0.5*(gp+gm)/(gp-gm)));
		}
	}
	{
		const voxel* vjm = vxl(vi->i, vi->j-1, vi->k);
		const voxel* vjp = vxl(vi->i, vi->j+1, vi->k);
		if (vjm && vjp)  {
			float gp = vjp->R-vi->R;
			float gm = vi->R-vjm->R;
			if(abs(gp-gm)>0.01)  disp.y = max(-0.49, min(0.49, -0.5*(gp+gm)/(gp-gm)));
		}
	}
	{
		const voxel* vjm = vxl(vi->i, vi->j, vi->k-1);
		const voxel* vjp = vxl(vi->i, vi->j, vi->k+1);
		if (vjm && vjp)  {
			float gp = vjp->R-vi->R;
			float gm = vi->R-vjm->R;
			if(abs(gp-gm)>0.01)  disp.z = max(-0.49, min(0.49, -0.5*(gp+gm)/(gp-gm)));
		}
	}
	if(b_i!=b_i->boss)  {
		dbl3 BosKidVec=*b_i - *(b_i->boss);	
		disp -= 0.95*((BosKidVec&disp)/(magSqr(BosKidVec)+1e-12))*BosKidVec;
	}
	b_i->fi=vi->i-_mp5+disp.x;   b_i->fj=vi->j-_mp5+disp.y;   b_i->fk=vi->k-_mp5+disp.z;
	b_i->R = vi->R+0.95*mag(disp);
}


void medialSurface::moveUphillp1(medialBall* bi) // const
{/// Refines the maximal-sphere location, by moving it uphil the gradient of the distance-map, potentially relocates to new voxels




	const voxel* vi = vxl(bi->fi, bi->fj, bi->fk);
	dbl3 disp(0.,0.,0.), grad(0.,0.,0.);

	{
		const voxel* vjm = vxl(vi->i-1, vi->j, vi->k);
		const voxel* vjp = vxl(vi->i+1, vi->j, vi->k);
		if (vjm && vjp)  {
			float gp = vjp->R-vi->R;
			float gm = vi->R-vjm->R;
			grad.x = 0.5*(gp+gm);
			if(abs(gp-gm)>0.01)  disp.x = max(-0.59, min(0.59, -0.5*(gp+gm)/(gp-gm)));
		}
	}
	{
		const voxel* vjm = vxl(vi->i, vi->j-1, vi->k);
		const voxel* vjp = vxl(vi->i, vi->j+1, vi->k);
		if (vjm && vjp)  {
			float gp = vjp->R-vi->R;
			float gm = vi->R-vjm->R;
			grad.y = 0.5*(gp+gm);
			if(abs(gp-gm)>0.01)  disp.y = max(-0.59, min(0.59, -0.5*(gp+gm)/(gp-gm)));
		}
	}
	{
		const voxel* vjm = vxl(vi->i, vi->j, vi->k-1);
		const voxel* vjp = vxl(vi->i, vi->j, vi->k+1);
		if (vjm && vjp)  {
			float gp = vjp->R-vi->R;
			float gm = vi->R-vjm->R;
			grad.z = 0.5*(gp+gm);
			if(abs(gp-gm)>0.01)  disp.z = max(-0.59, min(0.59, -0.5*(gp+gm)/(gp-gm)));
		}
	}
	disp+=1.4*grad;
	
	if(bi!=bi->boss)  {
		dbl3 BosKidVec=*bi - *(bi->boss);
		disp -= 0.5*((BosKidVec&disp)/(magSqr(BosKidVec)+1e-12))*BosKidVec; 
	}
	disp/=(0.55*mag(disp)+0.05);

	voxel* vxlj=vxl(bi->fi+disp[0],bi->fj+disp[1],bi->fk+disp[2]);
	if(vxlj && vi!=vxlj && vxlj->R>vi->R && vxlj->ball==nullptr )  {
			bi->fi=vxlj->i-_mp5;
			bi->fj=vxlj->j-_mp5;
			bi->fk=vxlj->k-_mp5;
			bi->R=vxlj->R;
			bi->vxl->ball=nullptr;		bi->vxl=vxlj;
			vxlj->ball=bi;
			//++nrelocations;
	}

		//cout<<nrelocations<<" relocations "<<endl;
}



void makeFriend(medialBall* vi, medialBall* vj)  {
	if (vj->R > vi->R) {medialBall* tmp=vi; vi=vj; vj=tmp;}
	if ((vi->R < 1.5*vj->R)  && (!vi->isNei(vj)) && (!vi->inParents(vj)) && !vj->inParents(vi))  {
		vi->addNei(vj);
		vj->addNei(vi);
	}
}




/*inline double cosAngleWithBossPerD(const medialBall* a, const medialBall* b)  {
	dbl3 v1=*a-*b;
	if (b==b->boss) return 1./(mag(v1));
	dbl3 v2=*b-*(b->boss);
	double dotProd=v2&v1;
	return sqrt(dotProd*dotProd/(magSqr(v1)*magSqr(v1)*magSqr(v2)));
}*/




void medialSurface::competeForParent(medialBall* vi, medialBall* vj)  {///## Select parent sphere between each nearby maximal-spheres to construct their hirachy


	const double noise=_MSNoise;

	double ri = vi->R;
	double rj = vj->R;
	double riSqr = ri*ri;
	double rjSqr = rj*rj;
	double dSqr = distSqr(vi,vj);

	double wsinv=1./(riSqr+rjSqr);
	//  double wi(dSqr+rjSqrlim-rlim), wj(dSqr+rlim-rjSqrlim);



	const voxel* middlevxl=vxl( wsinv*(vi->fi*rjSqr+vj->fi*riSqr), wsinv*(vi->fj*rjSqr+vj->fj*riSqr), wsinv*(vi->fk*rjSqr+vj->fk*riSqr) );
	if ( middlevxl && middlevxl->R>min(ri,rj)*_midRf-0.5 && 1.01*sqrt(dSqr)<ri+rj+1.+1.*noise )  {



		if (vj->boss == vj)
		{	if (vi->mastrSphere() != vj)  {
				if (ri >= rj)       				vj->boss = vi;
				else if (vi->boss->R <= rj)		vi->boss = vj;
				else if (ri >= rj-noise && ri*_vmvRadRelNf+1.*noise>=rj) vj->boss = vi;
		}	}
		else if (vi->boss == vi)
		{   if (vj->mastrSphere() != vi)  {
				if (rj>=ri)         				vi->boss = vj;
				else if (vj->boss->R <= ri)		vj->boss = vi;
				else if (rj>=ri-noise && rj*_vmvRadRelNf+1.*noise>=ri) vi->boss = vj;
		}	}


		medialBall* mvi=vi->mastrSphere();
		medialBall* mvj=vj->mastrSphere();
		if (mvi != vj && mvj != vi)  {
	

		 if (mvi==mvj)  {

			short leveli=vi->level();
			short levelj=vj->level();






			if ( leveli+1 < levelj &&
					(vj->boss->R-vj->R+2.*noise)/(dist(vj->boss,vj)+0.25) <  (vi->R-vj->R+2.*noise+0.01)/(dist(vi,vj)+0.2)
				) vj->boss = vi;
			else if ( leveli > levelj+1 &&
							(vi->boss->R-vi->R+2.*noise)/(dist(vi->boss,vi)+0.25) <  (vj->R-vi->R+2.*noise+0.01)/(dist(vj,vi)+0.2)
						) vi->boss = vj;
			else  { //if (leveli == levelj   or  leveli==levelj+1 ...)  

				     if ( leveli > levelj && (vi->boss->R-vi->R+2.*noise)/(dist(vi->boss,vi)+1.2) <  (vj->R-vi->R+2.*noise)/(dist(vj,vi)+1.3)
						&& !vj->inParents(vi) ) vi->boss = vj;
				else if ( leveli < levelj && (vj->boss->R-vj->R+2.*noise)/(dist(vj->boss,vj)+1.2) <  (vi->R-vj->R+2.*noise)/(dist(vi,vj)+1.3)
						&& !vi->inParents(vj) ) vj->boss = vi;
				else if(middlevxl && middlevxl->R>=0.45*(ri+rj)-1. && sqrt(dSqr)<(ri+rj)*0.5+2. )
					makeFriend(vi,vj);
			}
			if(vi->mastrSphere() != vj->mastrSphere()) {cout<<"sdsdsds"<<endl; exit(-1);}
		 }
		 else  {//if (mvi!=mvj)  




			if (distSqr(mvi,mvj) <= _lenNf*(0.5*(mvi->R+mvj->R)+2.*noise)*(0.5*(mvi->R+mvj->R)+2.*noise))  {

				if(mvi->R < mvj->R)  {
					medialBall* tmp=vi; vi=vj; vj=tmp;     tmp=mvi; mvi=mvj; mvj=tmp;
				}

			   if (mvj->R < _vmvRadRelNf * vj->R +noise &&  mvj->R < _vmvRadRelNf * vi->R +noise && mvj->R < _vmvRadRelNf * vi->boss->R +noise)  {
				  while (vj != vj->boss  &&  mvj->R < _vmvRadRelNf*vj->boss->R+noise )  {
					 medialBall* pvj = vj->boss;		 vj->boss = vi;			vi=vj;			vj=pvj;
					 cout<<'.';
				  }
				  if(vj->boss==vj && vi->mastrSphere()!=vj)   vj->boss = vi;
			   }

			}

			if (vi != vj->boss)
			{

				mvi=vi->mastrSphere();
				mvj=vj->mastrSphere();

				short leveli=vi->level();
				short levelj=vj->level();

				float distAvg = dist(mvj,mvi)*1.+0.5*noise;
				while ( leveli>=levelj && (vi->boss->R-vi->R+0.55*noise)/(dist(mvi,vi)+distAvg) <  (vj->R-vi->R+0.5*noise)/(dist(mvj,vi)+distAvg))  {
				  medialBall* pvi = vi->boss;		vi->boss = vj;	vj=vi;	vi=pvi; ++levelj; --leveli;
				}
				while ( levelj>=leveli && (vj->boss->R-vj->R+0.55*noise)/(dist(mvj,vj)+distAvg) <  (vi->R-vj->R+0.5*noise)/(dist(mvi,vj)+distAvg) )  {
				  medialBall* pvj = vj->boss;		vj->boss = vi;	vi=vj;	vj=pvj; ++leveli; --levelj;
				}

				makeFriend(vi,vj);

			}
		 }
		}
	}

}


void medialSurface::findBoss(medialBall* vi)  {




	const float  x = vi->fi,   y = vi->fj,   z = vi->fk;
	const float  ripp = vi->R*0.6+2.*_MSNoise+2.;
	const float  ex = x+ripp;
	for (float xpa=2.*x-ex; xpa<=ex; xpa+=1.0f)  {	float  ey = y+sqrt(ripp*ripp-(xpa-x)*(xpa-x));
			for (float ypb=2.*y-ey; ypb<=ey; ypb+=1.0f)  {	float ez = z+sqrt(ripp*ripp-(xpa-x)*(xpa-x)-(y-ypb)*(y-ypb));
				for (float zpc=2.*z-ez; zpc<=ez; zpc+=1.0f)  {voxel* vj = this->vxl(xpa, ypb, zpc);
					if ((vj != nullptr) && vj->ball && (*vi != *vj) )
					{ //--------------------------------------------------------
							competeForParent(&*vi,vj->ball);
					} //--------------------------------------------------------
				}
			}
	}

}



void medialSurface::createBallsAndHierarchy()  {/// Create distance map, maximal-spheres, and their hirarchy (medial-surface connectivity)





//	mediaAxes medAxis(*this, minRP, clipOutSideBallFraction, clipOutSideBallFraction*0.5+0.);


	buildvoxelspace();


	calc_distmaps();




   for (int i=0;i<_nRSmoothing;++i)        smoothRadius();






	nBalls = 0;   double rBalls = 0.;
	std::vector<voxel>::iterator vit = vxlSpace.begin()-1;
	const std::vector<voxel>::iterator vend = vxlSpace.end();
	while (++vit < vend)  {
		if (vit->R >= _minRp) {vit->ball = &(ToBeAssigned);	++nBalls; rBalls+=vit->R;}
		else 	              {vit->ball = nullptr;}
	}
	cout<< "\n  number of potential maximal spheres: " << nBalls << ",  average radius = "<<rBalls/nBalls<<endl;




	paradox_pre_removeincludedballI();

	paradoxremoveincludedballI();




	cout<< " collecting maximal balls out of "<<nBalls<<endl;





	std::vector<voxel*>tvs;tvs.reserve(nBalls);
	{
		std::vector<voxel>::iterator ti = vxlSpace.begin()-1;
		std::vector<voxel>::iterator tend = vxlSpace.end();
		while (++ti<tend) if (ti->ball)  {
			if ((ti->R)>=_minRp)		tvs.push_back(&(*(ti)));     else 	cout<<"  sdsd ";
		}

		cout<< " sorting "<<int(tvs.size())<<" maximal balls"<<endl;
		sort(tvs.begin(), tvs.end(), metaballcomparer());
	}

	ballSpace.reserve(nBalls);
	{
		std::vector<voxel*>::iterator ti = tvs.begin()-1;
		std::vector<voxel*>::iterator tend = tvs.end();
		while (++ti<tend)  {
			ballSpace.emplace_back(*ti,0);
			(*ti)->ball = &*(ballSpace.rbegin());
		}
	}



	const std::vector<medialBall>::iterator voxend = ballSpace.end();

	{
		std::vector<medialBall>::iterator vi = ballSpace.begin()-1;
		while (++vi != voxend)    moveUphill(&*vi);
	}
	{
		std::vector<medialBall>::iterator vi = ballSpace.begin()-1;
		while (++vi != voxend)    moveUphillp1(&*vi);
	}

	{
		std::vector<medialBall>::iterator vi = ballSpace.begin()-1;
		while (++vi != voxend)    moveUphill(&*vi);
	}

	cout<< " creating ball hierarchy:";  cout.flush();
	const std::vector<medialBall>::iterator voxp = ballSpace.begin();
	{
		std::vector<medialBall>::iterator vi = ballSpace.begin();
		while (vi != voxend)  {
			findBoss(&*vi);
			if ((vi-voxp)%100000 == 0)   {cout<< "\r   ball: " << int(vi-voxp); cout.flush();}
			++vi;
		}
		cout<< "\r   ball: " << int(vi-voxp)<<endl;
	}




}



