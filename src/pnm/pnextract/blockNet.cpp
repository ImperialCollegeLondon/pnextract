
#include "blockNet.h"

//#include "vxlImage_manip.h"

//clock_t myTime::start = clock();




int blockNetwork::nBP6 = 2; //   will be set in inputDataNE::init...
int blockNetwork::SideImax = 1; //=nBP6-1




void blockNetwork::createMedialSurface(medialSurface*& srf, inputDataNE& cfg, size_t startValue)  { /// calls medialSurface::createBallsAndHierarchy(), ...


	{

	//  clipROutx clipROutyz  midRFrac  RMedSurfNoise  lenNf   vmvRadRelNf  nRSmoothing   RCorsf   RCors

	medialSurface* medSurf = new medialSurface(cfg);
	medSurf->createBallsAndHierarchy();
	srf=(medSurf);

  }

}



void blockNetwork::CreateVElem(size_t startValue)  { ///### map pore labels from maximal spheres to the image blockNetwork::VElems, standing for VoxelElements.





	cout<< "\n\nCreating pore elements:"<<endl;



	if(poreIs.empty())  {
		sides[0] = voxel(-1, cg.ny/2, cg.nz/2, _SideRadius); //minX
		poreNE* tmp=new poreNE();  tmp->volumn=cg.ny*cg.nz;  tmp->mb=new medialBall(&sides[0],-1);
		poreIs.push_back(tmp);

		sides[1] = voxel(cg.nx+1, cg.ny/2, cg.nz/2, _SideRadius); //maxX
		tmp = new poreNE();  tmp->volumn = cg.ny*cg.nz; tmp->mb=new medialBall(&sides[1],-1);
		poreIs.push_back(tmp);


		if(nBP6==6) {
			sides[2] = voxel(cg.nx/2, -1, cg.nz/2, _SideRadius); //minY
			tmp = new poreNE();  tmp->volumn = cg.nx*cg.nz; tmp->mb=new medialBall(&sides[2],-1);
			poreIs.push_back(tmp);

			sides[3] = voxel(cg.nx/2, cg.ny+1, cg.nz/2, _SideRadius); //maxY
			tmp = new poreNE();  tmp->volumn = cg.nx*cg.nz; tmp->mb=new medialBall(&sides[3],-1);
			poreIs.push_back(tmp);

			sides[4] = voxel(cg.nx/2, cg.ny/2, -1, _SideRadius); //minZ
			tmp = new poreNE();  tmp->volumn = cg.ny*cg.nx; tmp->mb=new medialBall(&sides[4],-1);
			poreIs.push_back(tmp);

			sides[5] = voxel(cg.nx/2, cg.ny/2, cg.nz+1, _SideRadius);//maxZ
			tmp = new poreNE();  tmp->volumn = cg.ny*cg.nx; tmp->mb=new medialBall(&sides[5],-1);
			poreIs.push_back(tmp);
		}
	}


	VElems.reset(cg.nx+2,cg.ny+2,cg.nz+2,-257);
	VElems.X0Ch()=cg.VImage.X0()-cg.VImage.dx();  VElems.dxCh()=cg.VImage.dx();
	int nVVs=0;
	for (int iz = 0; iz<cg.nz; ++iz)  {for (int iy = 0; iy<cg.ny; ++iy)
	 {const segments& s = cg.segs_[iz][iy];
	  for (int ix = 0; ix<s.cnt; ++ix)  {int value=-1-int(s.s[ix].value);
		nVVs=max(nVVs,value);
	  	for (int i = s.s[ix].start; i<s.s[ix+1].start; ++i)  {
			 VElems(i+1,iy+1,iz+1) = value;
		}
	  }
	 }
	} ++nVVs;



	VElems.setSlice('i',0,       0);
	VElems.setSlice('i',cg.nx+1, 1);
	if(nBP6==6){
		VElems.setSlice('j',0,       2); //TODO: check
		VElems.setSlice('j',cg.ny+1, 3);
		VElems.setSlice('k',0,       4);
		VElems.setSlice('k',cg.nz+1, 5);
	} else {
		VElems.setSlice('j',0,      -1-2-nVVs); //TODO: check, this means negatives can not be sued to uniquely identify pores
		VElems.setSlice('j',cg.ny+1,-1-3-nVVs);
		VElems.setSlice('k',0,      -1-4-nVVs);
		VElems.setSlice('k',cg.nz+1,-1-5-nVVs);
	}


	firstPore = 2;




	{

      int uasyned = -1;

		const std::vector<medialBall>& balspc = srf->ballSpace;



		firstPores=(poreIs.size());
		for(const auto& bi:balspc)  {
		  if (bi.boss == &bi)  {
				VElems(bi.fi+1, bi.fj+1, bi.fk+1) = poreIs.size();
				poreNE* tmp = new poreNE();
				tmp->mb = &bi;
				poreIs.push_back(tmp);
		  }
		}
		cout<<"\n created "<<poreIs.size()<<" pores (+boundaries)  for up to index "<<0<<endl;
		lastPores=(poreIs.size()-1);



	  const int firstPoreInd=firstPores;
	  const int lastPoreInd=lastPores;
		cout<<" mapping pores indices to image, for index "<<0<< ":  "<<firstPores<<" to "<<lastPores<<",  unasigned:"<<uasyned<<endl;



		for(const auto& bi:balspc)  if (bi.boss)  {
			medialBall* mastrSphere = bi.mastrSphere();
			float apmi(mastrSphere->fi), bpmi(mastrSphere->fj), cpmi(mastrSphere->fk);
			int VElemV = VElems(apmi+1,bpmi+1,cpmi+1);
			ensure(VElemV>0 && VElemV<len(poreIs));


			const float   x = bi.fi,   y = bi.fj,   z = bi.fk;
			apmi =         x - mastrSphere->fi; bpmi = y - mastrSphere->fj;   cpmi = z - mastrSphere->fk;
			float R = bi.R;
			int r2 = std::max(R*0.25-1.,1.001)*std::max(R*0.25-1.,1.001);

			float ex = sqrt(r2);
			for (float xpa = max((x-ex),0.5f); xpa <=  min((x+ex),cg.nx-0.5f); xpa+=1.0f)  {
				 float ey = sqrt(r2-(xpa-x)*(xpa-x));
				 for (float ypb = max((y-ey),0.5f); ypb <=  min((y+ey),cg.ny-0.5f); ypb+=1.0f)  { 
					float ez = sqrt(r2-(xpa-x)*(xpa-x)-(ypb-y)*(ypb-y));
				   for (float zpc = max((z-ez),0.5f); zpc <=  min((z+ez),cg.nz-0.5f); zpc+=1.0f)  {
					 int idj = VElems(xpa+1,ypb+1,zpc+1);
					   if      (idj == (-1-int(0)) )   VElems(xpa+1,ypb+1,zpc+1) = VElemV;
					   else if (VElemV != idj && (firstPoreInd<=idj && idj<=lastPoreInd))  {
						 voxel* vj=srf->vxl(xpa,ypb,zpc);
						 if (!vj->ball && vj->R<R)  {
						 	 const medialBall* mvj = poreIs[idj]->mb;
						 	 float  amj = xpa-mvj->fi,   bmj = ypb-mvj->fj,   cmj = zpc-mvj->fk;
						 	 float  ami = xpa - mastrSphere->fi,   bmi = ypb - mastrSphere->fj,   cmi = zpc - mastrSphere->fk;
							 if ( ami*ami+bmi*bmi+cmi*cmi < amj*amj+bmj*bmj+cmj*cmj)
						 		VElems(xpa+1,ypb+1,zpc+1) = VElemV;//elem->id;
						 }

					   }
				   }
				 }
			}
		}








		growPoresMedStrict(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPoresMedStrict(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPoresMedStrict(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPoresMedian(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPoresMedStrict(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPoresMedian(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPoresMedStrict(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPoresMedian(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPoresMedStrict(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		cout<<endl;
		growPoresMedian(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPoresMedStrict(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPoresMedian(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPoresMedStrict(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPoresMedian(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPoresMedStrict(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		cout<<endl;
		growPoresMedEqs(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPoresMedian(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPoresMedStrict(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPoresMedian(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPoresMedian(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPoresMedStrict(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPoresMedEqs(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPoresMedian(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPoresMedian(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPoresMedStrict(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		cout<<endl;
		growPoresMedEqs(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPoresMedEqs(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPoresMedEqs(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPoresMedEqs(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPoresMedEqs(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPoresMedEqsLoose(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPoresMedEqs(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		cout<<endl;
		growPoresMedEqs(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPoresMedEqsLoose(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPoresMedEqs(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPoresMedEqsLoose(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPoresMedEqs(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPoresMedEqsLoose(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPoresMedEqs(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPoresMedEqsLoose(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPoresMedEqs(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPoresMedEqsLoose(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPoresMedEqs(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPoresMedEqsLoose(cg, VElems,  firstPores, lastPores, poreIs, uasyned);

		cout<<endl;

		growPores(VElems, firstPores, lastPores, uasyned);
		growPoresMedian(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPoresMedEqs(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPores(VElems, firstPores, lastPores, uasyned);
		growPoresMedian(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPores(VElems, firstPores, lastPores, uasyned);
		growPores(VElems, firstPores, lastPores, uasyned);
		growPores(VElems, firstPores, lastPores, uasyned);
		cout<<endl;
		growPores(VElems, firstPores, lastPores, uasyned);
		growPores(VElems, firstPores, lastPores, uasyned);
		growPores(VElems, firstPores, lastPores, uasyned);
		growPoresMedian(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPores(VElems, firstPores, lastPores, uasyned);
		growPores(VElems, firstPores, lastPores, uasyned);
		growPores(VElems, firstPores, lastPores, uasyned);
		cout<<endl;
		growPores(VElems, firstPores, lastPores, uasyned);
		growPores(VElems, firstPores, lastPores, uasyned);
		growPores(VElems, firstPores, lastPores, uasyned);
		growPores(VElems, firstPores, lastPores, uasyned);
		growPores(VElems, firstPores, lastPores, uasyned);
		growPores(VElems, firstPores, lastPores, uasyned);
		cout<<endl;
		growPores(VElems, firstPores, lastPores, uasyned);
		growPores_X2(VElems, firstPores, lastPores, uasyned);
		growPores_X2(VElems, firstPores, lastPores, uasyned);
		growPores_X2(VElems, firstPores, lastPores, uasyned);
		growPores_X2(VElems, firstPores, lastPores, uasyned);
		cout<<endl;



		medianElem(cg, VElems,  firstPores, lastPores, poreIs);
		medianElem(cg, VElems,  firstPores, lastPores, poreIs);
		medianElem(cg, VElems,  firstPores, lastPores, poreIs);

		medianElem(cg, VElems,  firstPores, lastPores, poreIs);
		medianElem(cg, VElems,  firstPores, lastPores, poreIs);

		cout<<endl;

		growPores(VElems, firstPores, lastPores, uasyned);
		while(growPores_X2(VElems, firstPores, lastPores, uasyned));
		growPores(VElems, firstPores, lastPores, uasyned);

		cout<<endl;


		retreatPoresMedian(cg, VElems,  firstPores, lastPores, poreIs, uasyned);



		for(const auto& bi:balspc) {
		      medialBall* mastrSphere = bi.mastrSphere();
			  VElems(bi.fi+1, bi.fj+1, bi.fk+1) = VElems(mastrSphere->fi+1, mastrSphere->fj+1, mastrSphere->fk+1);
		}
		growPoresMedian(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPoresMedEqs(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPoresMedEqs(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPoresMedEqs(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPoresMedian(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPoresMedEqs(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPoresMedEqs(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPoresMedEqsLoose(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPoresMedEqs(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPoresMedEqsLoose(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPoresMedEqs(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPoresMedEqsLoose(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPoresMedEqs(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPores(VElems, firstPores, lastPores, uasyned);
		growPores_X2(VElems, firstPores, lastPores, uasyned);
		growPoresMedEqs(cg, VElems,  firstPores, lastPores, poreIs, uasyned);
		growPores(VElems, firstPores, lastPores, uasyned);
		growPores_X2(VElems, firstPores, lastPores, uasyned);

		cout<<endl;

		medianElem(cg, VElems,  firstPores, lastPores, poreIs);
		for(const auto& bi:balspc) {
		      medialBall* mastrSphere = bi.mastrSphere();
			  VElems(bi.fi+1, bi.fj+1, bi.fk+1) = VElems(mastrSphere->fi+1, mastrSphere->fj+1, mastrSphere->fk+1);
        }
		growPoresMedEqsLoose(cg, VElems,  firstPores, lastPores, poreIs, uasyned);

		cout<<endl;
		cout<<endl;


	}



}




void blockNetwork::createNewThroats(medialSurface*& srf)  {

 //vector<int> iThroatFaces; iThroatFaces.reserve(throatIs.size()+poreIs.size()*10);
 { cout<<"\nlooking for connections, iFirst = "<< throatIs.size()<< " ..  "; cout.flush();


     for (int i=0; i<VElems.nx()-1 ; i++ )
   for (int k=1; k<VElems.nz()-1; ++k)
    for (int j=1; j<VElems.ny()-1; ++j)  {
		 int p1ID = VElems(i,j,k);
		 if(p1ID>= 0)  {
			int p2ID = VElems(i+1,j,k);
			if (p1ID != p2ID)
			{ if (p2ID>= 0)
			  {
				bool dir = p2ID>p1ID;
				if (!dir) {int tmp=p1ID; p1ID=p2ID; p2ID=tmp;}

				poreNE* p1 = poreIs[p1ID];
				poreNE* p2 = poreIs[p2ID];

				if (!p1) cout<<"  ERROR p1ID"<<p1ID<<endl;
				int tIDNext = throatIs.size();
				std::pair<std::map<int,int>::iterator,bool> ret = p2->contacts.insert(std::pair<int,int>( p1ID,tIDNext));
				if(ret.second)
				{

					if (ret.second != p1->contacts.insert(std::pair<int,int>( p2ID,tIDNext)).second ) cout<<"Errorpnm821-2 "<<p1ID<<" "<<p2ID<<endl;
					throatIs.push_back(new throatNE(tIDNext,p1ID,p2ID));
					//if (p1ID<2 &&p2ID<2) cout <<"Errskziw1: p1ID<2 &&p2ID<2: tid: "<< tIDNext<<"  "<<i<<" "<<j<<" "<<k<<" " <<"  "<<p1->mb->fi<<" "<<p1->mb->fj<<" "<<p1->mb->fk<<" "<<"   "<<p2->mb->fi<<" "<<p2->mb->fj<<" "<<p2->mb->fk<<" "<<endl;
				}
				throatNE* trot = throatIs[ret.first->second];
				trot->CrosArea.x += (2*dir-1);
				//trot->C[0] += int3{{i-1,j-1,k-1}};
			  }

			  if (p1ID>= firstPore)	++(poreIs[p1ID]->surfaceArea);
			  if (p2ID>= firstPore)	++(poreIs[p2ID]->surfaceArea);
			}
			
			if (p1ID>= firstPore) { ++(poreIs[p1ID]->volumn);}
		 }
     }



	(cout<<throatIs.size()<<" .. ").flush();


    for (int j=0; j<VElems.ny()-1 ; j++ ) //     y >-<
   for (int k=1; k<VElems.nz()-1; ++k)
     for (int i=1; i<VElems.nx()-1; ++i)  {
		 int p1ID = VElems(i,j,k);
		 if(p1ID>= 0)  {
			int p2ID = VElems(i,j+1,k);
			if (p1ID != p2ID)  { if (p2ID>= 0)
			  {
				bool dir = p2ID>p1ID;
				if (!dir) {int tmp=p1ID; p1ID=p2ID; p2ID=tmp;}

				poreNE* p1 = poreIs[p1ID];
				poreNE* p2 = poreIs[p2ID];

				if (!p1) cout<<"  ERROR p1ID"<<p1ID<<endl;
				int tIDNext = throatIs.size();
				std::pair<std::map<int,int>::iterator,bool> ret = p2->contacts.insert(std::pair<int,int>( p1ID,tIDNext));
				if(ret.second)  {
					if (ret.second != p1->contacts.insert(std::pair<int,int>( p2ID,tIDNext)).second ) cout<<"Errorpnm821-2 "<<p1ID<<" "<<p2ID<<endl;
					throatIs.push_back(new throatNE(tIDNext,p1ID,p2ID));
					//if (p1ID<2 &&p2ID<2) cout <<"Errskziw2: p1ID<2 &&p2ID<2: tid: "<< tIDNext<<endl;
				}
				throatNE* trot = throatIs[ret.first->second];
				trot->CrosArea.y += (2*dir-1);
				//trot->C[1] += int3{{i-1,j-1,k-1}};
			  }

			  if (p1ID>= firstPore)	++(poreIs[p1ID]->surfaceArea);
			  if (p2ID>= firstPore)	++(poreIs[p2ID]->surfaceArea);
			}
		 }
     }



	(cout<<throatIs.size()<<" .. ").flush();


   for (int k=0; k<VElems.nz()-1; ++k) //     z >-<
    for (int j=1; j<VElems.ny()-1; ++j)
     for (int i=1; i<VElems.nx()-1; ++i)  {
		 int p1ID = VElems(i,j,k);
		 if(p1ID>= 0)  {
			int p2ID = VElems(i,j,k+1);
			if (p1ID != p2ID)  { if (p2ID>= 0)
			  {
				bool dir = p2ID>p1ID;
				if (!dir) {int tmp=p1ID; p1ID=p2ID; p2ID=tmp;}
				poreNE* p1 = poreIs[p1ID];
				poreNE* p2 = poreIs[p2ID];

				if (!p1) cout<<"  ERROR p1ID"<<p1ID<<endl;
				int tIDNext = throatIs.size();
				std::pair< std::map<int,int>::iterator, bool > ret = p2->contacts.insert(std::pair<int,int>( p1ID,tIDNext));
				if(ret.second)  {
					if (ret.second != p1->contacts.insert(std::pair<int,int>( p2ID,tIDNext)).second ) cout<<"Errorpnm821-2 "<<p1ID<<" "<<p2ID<<endl;
					throatIs.push_back(new throatNE(tIDNext,p1ID,p2ID));
					//if (p1ID<2 &&p2ID<2) cout <<"Errskziw3: p1ID<2 &&p2ID<2: tid: "<< tIDNext<<endl;
				}
				throatNE* trot = throatIs[ret.first->second];
				trot->CrosArea.z += (2*dir-1);
				//trot->C[2] += int3{{i-1,j-1,k-1}};
			  }

			  if (p1ID>= firstPore)	++(poreIs[p1ID]->surfaceArea);
			  if (p2ID>= firstPore)	++(poreIs[p2ID]->surfaceArea);
			}
		 }
     }




	(cout<<throatIs.size()<<", ").flush();


 }

	cout<<" nThroats: "<< throatIs.size()<<endl;




	nNodes=poreIs.size();
	nTrots=throatIs.size();
	nElems=nNodes+nTrots;

	cout<<"nElems:  "<<nElems<<" = "<<poreIs.size()<<" + "<<throatIs.size()<<endl;





	throadAdditBalls.reserve(throatIs.size()*5);// to improve the efficiency when later generating and storing new maximal balls for throat surfaces


	cout<<"\ncalc throat properties: "; cout.flush();

  { cout<<" collecting face voxels,  ";cout.flush();

	for (auto tr: throatIs)  {
		tr->toxels2.reserve(abs(tr->CrosArea[0])+abs(tr->CrosArea[1])+abs(tr->CrosArea[2])+1);
		tr->toxels1.reserve(abs(tr->CrosArea[0])+abs(tr->CrosArea[1])+abs(tr->CrosArea[2])+1);
	}



	int nMultiTouchErrors = 0;
	forAllkji_1(VElems)  {
		int p1ID = VElems(i,j,k);
		if (p1ID>= firstPore)  {
			std::set<int> neis;
			int neiPID;

			neiPID = VElems(i-1,j,k); 	if ((p1ID != neiPID) && (neiPID>= 0) ) neis.insert(neiPID);
			neiPID = VElems(i+1,j,k); 	if ((p1ID != neiPID) && (neiPID>= 0) ) neis.insert(neiPID);
			neiPID = VElems(i,j-1,k); 	if ((p1ID != neiPID) && (neiPID>= 0) ) neis.insert(neiPID);
			neiPID = VElems(i,j+1,k); 	if ((p1ID != neiPID) && (neiPID>= 0) ) neis.insert(neiPID);
			neiPID = VElems(i,j,k-1); 	if ((p1ID != neiPID) && (neiPID>= 0) ) neis.insert(neiPID);
			neiPID = VElems(i,j,k+1); 	if ((p1ID != neiPID) && (neiPID>= 0) ) neis.insert(neiPID);
			for (int nei:neis)  {
				poreNE* p1 = poreIs[p1ID];
				throatNE* trot = throatIs[p1->contacts[nei]];
				if(p1ID>nei)  {
					dbgAsrt(srf->vxl(i-1,j-1,k-1));
					trot->toxels2.push_back( (srf->vxl(i-1,j-1,k-1)) ); 	
				}
				else {
					dbgAsrt(srf->vxl(i-1,j-1,k-1));
					trot->toxels1.push_back( (srf->vxl(i-1,j-1,k-1)) );
				}
			}
			if (neis.size()>1) ++nMultiTouchErrors;
		}
   }
	if (nMultiTouchErrors>0)  cout<<"\n   Warning: "<< nMultiTouchErrors <<" voxels being in touch to more than two pores"<<endl;
  }



	for (auto tr: throatIs)  {
        if (tr->toxels2.empty() || tr->toxels2.empty() ) (cout<<"  ERROR1017, toxls size:"<<tr->toxels2.size()<<" "<<tr->toxels1.size()<<"  ").flush();

		sort(tr->toxels2.begin(), tr->toxels2.end(), metaballcomparer());
		sort(tr->toxels1.begin(), tr->toxels1.end(), metaballcomparer());
	}




	cout<<" calculating throat radii";cout.flush();
	for (throatNE* tr : throatIs)  { 
		if (tr->toxels2.size()>0)  {
			voxel* tvox2=*(tr->toxels2.begin());// get the largest distance map throat voxel
			if (tvox2->ball!=nullptr)  {
				medialBall* vbi = (*tr->toxels2.begin())->ball;
				vbi->type = 5;

						medialBall* mvi=vbi->mastrSphere();
						if (mvi!=nullptr && mvi!=vbi && tr->e2 != VElems(mvi->fi+1, mvi->fj+1, mvi->fk+1))	cout<<" Dmb2rrr  "<< VElems(mvi->fi+1, mvi->fj+1, mvi->fk+1)<<"   ";
			}
			else  {
				tvox2->ball = new medialBall(tvox2, 15); throadAdditBalls.push_back(tvox2->ball);
				srf->moveUphill(tvox2->ball);
			}


		}



      if (!tr->toxels1.empty())  {
			sort(tr->toxels1.begin(), tr->toxels1.end(), metaballcomparer());
			voxel* tvox1=*(tr->toxels1.begin());
			if (tvox1->ball)  {
				medialBall* vbi = tvox1->ball;
				vbi->type = 6;

				medialBall* mvi=vbi->mastrSphere();
				if (mvi && mvi!=vbi && tr->e1 != VElems(mvi->fi+1, mvi->fj+1, mvi->fk+1))		cout<<" Dmb1rrr  "<< VElems(mvi->fi+1, mvi->fj+1, mvi->fk+1) <<"   ";
			}
		  else  {
            tvox1->ball = new medialBall(tvox1,16); throadAdditBalls.push_back(tvox1->ball);
            srf->moveUphill(tvox1->ball);
		  }
	  }


	}

	cout<<"."<<endl;

}








#include "blockNet_write_cnm.cpp"
#include "blockNet_vxlManip.cpp"
