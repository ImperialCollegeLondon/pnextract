
#include "blockNet.h"

//~ #include "vxlImage_manip.h"

//~ clock_t myTime::start = clock();







void blockNetwork::createBallsAndHierarchy(medialSurface*& refs, inputDataNE& cfg, size_t startValue)
{ /// calls medialSurface::createBallsAndHierarchy(), ...



  {


//  clipROutx clipROutyz  midRFrac  RMedSurfNoise  lenNf   vmvRadRelNf  nRSmoothing   RCorsf   RCors
    refs=(new medialSurface(cfg));



    refs->createBallsAndHierarchy();


  }

}



void blockNetwork::CreateVElem(size_t startValue)
{ /// map pore labels to the image

	cout<< "\n\nCreating pore elements:"<<endl;
	outletV.i = cg.nx+1;
	outletV.j = cg.ny/2;
	outletV.k = cg.nz/2;
	outletV.R = cg.nx/4;

	inletV.i = -1;
	inletV.j = cg.ny/2;
	inletV.k = cg.nz/2;
	inletV.R = cg.nx/4;

	VElems.reset(cg.nx+2,cg.ny+2,cg.nz+2,-257); 
	VElems.X0Ch()=cg.VImage.X0()-cg.VImage.dx();  VElems.dxCh()=cg.VImage.dx();
	int nVVs=0;
	for (int iz = 0; iz<cg.nz; ++iz)
	{for (int iy = 0; iy<cg.ny; ++iy)
	 {const segments& s = cg.segs_[iz][iy];
	  for (int ix = 0; ix<s.cnt; ++ix)
	  {register int value=-1-int(s.s[ix].value);
		nVVs=max(nVVs,value);
	  	for (int i = s.s[ix].start; i<s.s[ix+1].start; ++i)
		{
			 VElems(i+1,iy+1,iz+1) = value;
		}
	  }
	 }
	} ++nVVs;

const int	Left  =0,		Right =1,	
				Bottom=nVVs+2,	Top   =nVVs+3,	
				Back  =nVVs+4,	Front =nVVs+5;

	VElems.setSlice('i',0,         Left  );
	VElems.setSlice('i',cg.nx+1,   Right );
	VElems.setSlice('j',0,      -1-Bottom);
	VElems.setSlice('j',cg.ny+1,-1-Top   );
	VElems.setSlice('k',0,      -1-Back  );
	VElems.setSlice('k',cg.nz+1,-1-Front );




	//for (int iz = 0; iz<cg.nz; ++iz)
		//for (int iy = 0; iy<cg.ny; ++iy)		{
			//VElems(0,iy+1,iz+1) = 0;
			//VElems(cg.nx+1,iy+1,iz+1) = 1;	}


	int nPores = 0;
	if(poreIs.empty())
	{ /// add inlet
		 poreElementI* tmp = new poreElementI();
		tmp->volumn = cg.ny*cg.nz;
		poreIs.push_back(tmp);
		tmp->mb=new medialBall(&inletV,-1);

	  /// add outlet
		tmp = new poreElementI();
		tmp->volumn = cg.ny*cg.nz;
		tmp->mb=new medialBall(&outletV,-1);
		poreIs.push_back(tmp);
	}




	firstPore = 2;



	///.  copy pores to V Elems

	{

      int uasyned = -1;

		const std::vector<medialBall>& balspc = rfs->ballSpace;


	  ///.  (create and) copy micro pores to poreIs
	  {
		firstPores=(poreIs.size());

		std::vector<medialBall>::const_iterator bi = balspc.begin();
		std::vector<medialBall>::const_iterator voxend = balspc.end();
		while (bi<voxend)
		{
		  if (bi->boss == &*bi)
		  { ++nPores;
			VElems(bi->fi+1+_pp5,bi->fj+1+_pp5,bi->fk+1+_pp5) = poreIs.size();
			poreElementI* tmp = new poreElementI();
			tmp->mb = &*bi;
			poreIs.push_back(tmp);
		  }
		  ++bi;
		}
		cout<<"\n created "<<nPores<<" pores  for up to index "<<0<<endl;

		lastPores=(poreIs.size()-1);

	  }


	  register const int firstPoreInd=firstPores;
	  register const int lastPoreInd=lastPores;
		cout<<" mapping pores indices to image for index "<<0<< ":  "<<firstPores<<" to "<<lastPores<<",  uasyned:"<<uasyned<<endl;

		int nWarnings = 0;
		std::vector<medialBall>::const_iterator bi = balspc.begin();
		const std::vector<medialBall>::const_iterator voxend = balspc.end();
		while (bi<voxend)
		{
		  if (bi->boss != NULL)
		  {
			  medialBall* mastrVxl = bi->mastrVxl();
			  register float apmi(mastrVxl->fi+_pp5), bpmi(mastrVxl->fj+_pp5), cpmi(mastrVxl->fk+_pp5);
			  int VElemV = VElems(apmi+1,bpmi+1,cpmi+1);


			  poreElementI *elem = poreIs[VElemV];

			  if (elem->ispore())
			  { ///.  transfer to EVElems
				register const float   x = bi->fi+_pp5,   y = bi->fj+_pp5,   z = bi->fk+_pp5;
				apmi =         x - mastrVxl->fi; bpmi = y - mastrVxl->fj;   cpmi = z - mastrVxl->fk;
				register float R = bi->R;
				register int r2 = std::max(R*0.25-1.0,1.001)*std::max(R*0.25-1.0,1.001);

				float ex = sqrt(r2);
				for (float xpa = max((x-ex),0.5f); xpa <=  min((x+ex),cg.nx-0.5f); xpa+=1.0f)
				{ float ey = sqrt(r2-(xpa-x)*(xpa-x));
				 for (float ypb = max((y-ey),0.5f); ypb <=  min((y+ey),cg.ny-0.5f); ypb+=1.0f)
				 { float ez = sqrt(r2-(xpa-x)*(xpa-x)-(ypb-y)*(ypb-y));
				   for (float zpc = max((z-ez),0.5f); zpc <=  min((z+ez),cg.nz-0.5f); zpc+=1.0f)
				   {
					 register int idj = VElems(xpa+1,ypb+1,zpc+1);
					   if (idj == (-1-int(0))  )
					   {
						  VElems(xpa+1,ypb+1,zpc+1) = VElemV;//elem->id;
					   }
					   else if (VElemV != idj && (firstPoreInd<=idj && idj<=lastPoreInd))
					   {
						voxel* vj=rfs->vxl(xpa,ypb,zpc);
						if (!vj->ball && vj->R<R)
						{
							 const medialBall* mvj = poreIs[idj]->mb;
							 register float  amj = xpa-mvj->fi,   bmj = ypb-mvj->fj,   cmj = zpc-mvj->fk;
							 register float  ami = xpa - mastrVxl->fi,   bmi = ypb - mastrVxl->fj,   cmi = zpc - mastrVxl->fk;
							if ( ami*ami+bmi*bmi+cmi*cmi < amj*amj+bmj*bmj+cmj*cmj)
								VElems(xpa+1,ypb+1,zpc+1) = VElemV;//elem->id;
						}

					   }
				   }
				 }
				}
			 }
			 else
				cout<<"error: @ ball "<< bi->fi<<" "<< bi->fj<<" "<<  bi->fk<<", expecting pore, found throat "	<<endl;

		  }
		  bi++;
		}
		if (nWarnings>0) cout<<" expecting pore for "<<nWarnings<<" balls, found throats"<<endl;



//        continue;
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


		for (bi = balspc.begin();bi<voxend; ++bi)
		{
		      medialBall* mastrVxl = bi->mastrVxl();
			  VElems(bi->fi+1+_pp5, bi->fj+1+_pp5, bi->fk+1+_pp5) = VElems(mastrVxl->fi+1+_pp5, mastrVxl->fj+1+_pp5, mastrVxl->fk+1+_pp5);
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
		for (bi = balspc.begin();bi<voxend; ++bi)
		{
		      medialBall* mastrVxl = bi->mastrVxl();
			  VElems(bi->fi+1+_pp5, bi->fj+1+_pp5, bi->fk+1+_pp5) = VElems(mastrVxl->fi+1+_pp5, mastrVxl->fj+1+_pp5, mastrVxl->fk+1+_pp5);
        }
		growPoresMedEqsLoose(cg, VElems,  firstPores, lastPores, poreIs, uasyned);

		cout<<endl;
		cout<<endl;


	}



}




void blockNetwork::createNewThroats(medialSurface*& refs)
{  /// looks for connections, calculates areas


 //vector<int> iThroatFaces; iThroatFaces.reserve(throatIs.size()+poreIs.size()*10);
 { cout<<"\nlooking for connections, iFirst = "<< throatIs.size()<< " ..  "; cout.flush();


     for ( int i = 0; i<int(VElems.size3()[0])-1 ; i++ ) ///. the order affects the throat orders
   for (int k = 1; k<int(VElems.size3()[2])-1 ; k++ ) ///.  mark throats x >-<
    for ( int j = 1; j<int(VElems.size3()[1])-1 ; j++ )
     {
		 int p1ID = VElems(i,j,k);
		 if(p1ID>= 0)
		 {
			int p2ID = VElems(i+1,j,k);
			if (p1ID != p2ID)
			{ if (p2ID>= 0)
			  {
				bool dir = p2ID>p1ID;
				if (!dir) {int tmp=p1ID; p1ID=p2ID; p2ID=tmp;}

				poreElementI* p1 = poreIs[p1ID];
				poreElementI* p2 = poreIs[p2ID];

				if (!p1) cout<<"  ERROR p1ID"<<p1ID<<endl;
				int tIDNext = throatIs.size();
				std::pair<std::map<int,int>::iterator,bool> ret = p2->contacts.insert(std::pair<int,int>( p1ID,tIDNext));
				if(ret.second)
				{
					if (ret.second != p1->contacts.insert(std::pair<int,int>( p2ID,tIDNext)).second ) cout<<"Errorpnm821-2 "<<p1ID<<" "<<p2ID<<endl;
					throatIs.push_back(new throatElementI(tIDNext,p1ID,p2ID));
					if (p1ID<2 &&p2ID<2) cout <<"Errskziw1: p1ID<2 &&p2ID<2: tid: "<< tIDNext<<"  "<<i<<" "<<j<<" "<<k<<" "
							  <<"  "<<p1->mb->fi<<" "<<p1->mb->fj<<" "<<p1->mb->fk<<" "<<"   "<<p2->mb->fi<<" "<<p2->mb->fj<<" "<<p2->mb->fk<<" "<<endl;
				}
				throatElementI* trot = throatIs[ret.first->second];
				trot->CrosArea.x += (2*dir-1);
				trot->C[0] += int3{{i-1,j-1,k-1}};
			  }

			  if (p1ID>= firstPore)	++(poreIs[p1ID]->surfaceArea);
			  if (p2ID>= firstPore)	++(poreIs[p2ID]->surfaceArea);
			}
			if (p1ID>= firstPore) { ++(poreIs[p1ID]->volumn);}
		 }
     }



	(cout<<throatIs.size()<<" .. ").flush();

    for ( int j = 0; j<int(VElems.size3()[1])-1 ; j++ )
   for ( int k = 1; k<int(VElems.size3()[2])-1 ; k++ ) ///.  mark throats y >-<
     for ( int i = 1; i<int(VElems.size3()[0])-1 ; i++ )
     {
		 int p1ID = VElems(i,j,k);
		 if(p1ID>= 0)
		 {
			int p2ID = VElems(i,j+1,k);
			if (p1ID != p2ID)
			{ if (p2ID>= 0)
			  {
				bool dir = p2ID>p1ID;
				if (!dir) {int tmp=p1ID; p1ID=p2ID; p2ID=tmp;}

				poreElementI* p1 = poreIs[p1ID];
				poreElementI* p2 = poreIs[p2ID];

				if (!p1) cout<<"  ERROR p1ID"<<p1ID<<endl;
				int tIDNext = throatIs.size();
				std::pair<std::map<int,int>::iterator,bool> ret = p2->contacts.insert(std::pair<int,int>( p1ID,tIDNext));
				if(ret.second)
				{
					if (ret.second != p1->contacts.insert(std::pair<int,int>( p2ID,tIDNext)).second ) cout<<"Errorpnm821-2 "<<p1ID<<" "<<p2ID<<endl;
					throatIs.push_back(new throatElementI(tIDNext,p1ID,p2ID));
					if (p1ID<2 &&p2ID<2) cout <<"Errskziw2: p1ID<2 &&p2ID<2: tid: "<< tIDNext<<endl;
				}
				throatElementI* trot = throatIs[ret.first->second];
				trot->CrosArea.y += (2*dir-1);
				trot->C[1] += int3{{i-1,j-1,k-1}};
			  }

			  if (p1ID>= firstPore)	++(poreIs[p1ID]->surfaceArea);
			  if (p2ID>= firstPore)	++(poreIs[p2ID]->surfaceArea);
			}
		 }
     }



	(cout<<throatIs.size()<<" .. ").flush();

   for ( int k = 0; k<int(VElems.size3()[2])-1 ; k++ ) ///.  mark throats z >-<
    for ( int j = 1; j<int(VElems.size3()[1])-1 ; j++ )
     for ( int i = 1; i<int(VElems.size3()[0])-1 ; i++ )
     {
		 int p1ID = VElems(i,j,k);
		 if(p1ID>= 0)
		 {
			int p2ID = VElems(i,j,k+1);
			if (p1ID != p2ID)
			{ if (p2ID>= 0)
			  {
				bool dir = p2ID>p1ID;
				if (!dir) {int tmp=p1ID; p1ID=p2ID; p2ID=tmp;}
				poreElementI* p1 = poreIs[p1ID];
				poreElementI* p2 = poreIs[p2ID];

				if (!p1) cout<<"  ERROR p1ID"<<p1ID<<endl;
				int tIDNext = throatIs.size();
				std::pair< std::map<int,int>::iterator, bool > ret = p2->contacts.insert(std::pair<int,int>( p1ID,tIDNext));
				if(ret.second)
				{
					if (ret.second != p1->contacts.insert(std::pair<int,int>( p2ID,tIDNext)).second ) cout<<"Errorpnm821-2 "<<p1ID<<" "<<p2ID<<endl;
					throatIs.push_back(new throatElementI(tIDNext,p1ID,p2ID));
					if (p1ID<2 &&p2ID<2) cout <<"Errskziw3: p1ID<2 &&p2ID<2: tid: "<< tIDNext<<endl;
				}
				throatElementI* trot = throatIs[ret.first->second];
				trot->CrosArea.z += (2*dir-1);
				trot->C[2] += int3{{i-1,j-1,k-1}};
			  }

			  if (p1ID>= firstPore)	++(poreIs[p1ID]->surfaceArea);
			  if (p2ID>= firstPore)	++(poreIs[p2ID]->surfaceArea);
			}
		 }
     }




	(cout<<throatIs.size()<<", ").flush();


 }

	cout<<" nThroats: "<< throatIs.size()<<endl;






	throadAdditBalls.reserve(throatIs.size()*5);

///. collect throat voxels, calculating inscribed radius
  cout<<"\ncalc throat properties: "; cout.flush();

  {cout<<" collecting face voxels,  ";cout.flush();

	for (std::vector<throatElementI*>::iterator titr = throatIs.begin();titr<throatIs.end();++titr)
	{
		(*titr)->toxels2.reserve(abs((*titr)->CrosArea[0])+abs((*titr)->CrosArea[1])+abs((*titr)->CrosArea[2])+1);//+1 is safety not needed
		(*titr)->toxels1.reserve(abs((*titr)->CrosArea[0])+abs((*titr)->CrosArea[1])+abs((*titr)->CrosArea[2])+1);//+1 is safety not needed
	}

	int nMultiTouchErrors = 0;
   for (int k = 1; k<int(VElems.size3()[2])-1 ; k++ ) ///.  mark throats x >-<
   {for ( int j = 1; j<int(VElems.size3()[1])-1 ; j++ )
    {for ( int i = 1; i<int(VElems.size3()[0])-1 ; i++ )
     {
		register int p1ID = VElems(i,j,k);
		if (p1ID>= firstPore)
		{
			 std::set<int> neis;
			 register int neiPID;
			 neiPID = VElems(i-1,j,k); 	if ((p1ID != neiPID) && (neiPID>= 0) ) neis.insert(neiPID);
			 neiPID = VElems(i+1,j,k); 	if ((p1ID != neiPID) && (neiPID>= 0) ) neis.insert(neiPID);
			 neiPID = VElems(i,j-1,k); 	if ((p1ID != neiPID) && (neiPID>= 0) ) neis.insert(neiPID);
			 neiPID = VElems(i,j+1,k); 	if ((p1ID != neiPID) && (neiPID>= 0) ) neis.insert(neiPID);
			 neiPID = VElems(i,j,k-1); 	if ((p1ID != neiPID) && (neiPID>= 0) ) neis.insert(neiPID);
			 neiPID = VElems(i,j,k+1); 	if ((p1ID != neiPID) && (neiPID>= 0) ) neis.insert(neiPID);
			for (std::set<int>::iterator neitr = neis.begin();neitr != neis.end();++neitr)
			{
				 poreElementI* p1 = poreIs[p1ID];
				 throatElementI* trot = throatIs[p1->contacts[*neitr]];
				 if(p1ID>*neitr) 
				 {
					d_assert(refs->vxl(i-1,j-1,k-1));
					trot->toxels2.push_back( (refs->vxl(i-1,j-1,k-1)) );
				 }
				 else 
				 {
					d_assert(refs->vxl(i-1,j-1,k-1));
					trot->toxels1.push_back( (refs->vxl(i-1,j-1,k-1)) );
				 }
			}
			if (neis.size()>1) ++nMultiTouchErrors;
		}
     }
    }
   }
	if (nMultiTouchErrors>0)  (cout<<"\n   Warning: "<< nMultiTouchErrors <<" voxels being in touch to more than two pores)\n      ").flush();
  }

	for (std::vector<throatElementI*>::iterator titr = throatIs.begin();titr<throatIs.end();++titr)
	{
        if ((*titr)->toxels2.empty() || (*titr)->toxels2.empty() )
            cout<<"  ERROR1017 "<<titr-throatIs.begin()<<", toxls size:"<<(*titr)->toxels2.size()<<" "<<(*titr)->toxels1.size()<<"  "
		    ;cout.flush();

		sort((*titr)->toxels2.begin(), (*titr)->toxels2.end(), metaballcomparer());
		sort((*titr)->toxels1.begin(), (*titr)->toxels1.end(), metaballcomparer());
	//~ cout<<titr - throatIs.begin()<< " n toxels:  "<<(*titr)->toxels1.size()<<"  "<<(*titr)->toxels2.size()<<endl;

	}


	nPores=poreIs.size();
	nTrots=throatIs.size();
	nElems=nPores+nTrots;

	cout<<"nElems:  "<<poreIs.size()<<" + "<<throatIs.size()<<" = "<<nElems<<endl;




  {cout<<" calculating throat radii";cout.flush();
	for (throatElementI* tr : throatIs )
	{ if (tr->toxels2.size()>0)
	  {
		voxel* tvox2=*(tr->toxels2.begin());
		if (tvox2->ball)
		{
			medialBall* vbi = (*tr->toxels2.begin())->ball;
			vbi->type = 5;

					medialBall* mvi=vbi->mastrVxl();
					if (mvi!=NULL && mvi!=vbi && tr->e2 != VElems(mvi->fi+1+_pp5, mvi->fj+1+_pp5, mvi->fk+1+_pp5))	cout<<" Dmb2rrr  "<< VElems(mvi->fi+1+_pp5, mvi->fj+1+_pp5, mvi->fk+1+_pp5)<<"   ";
		}
		else
		{
			tvox2->ball = new medialBall(tvox2, 15); throadAdditBalls.push_back(tvox2->ball);
			refs->moveUphill(tvox2->ball);
		}


	  }



      if (!tr->toxels1.empty())
      {
			sort(tr->toxels1.begin(), tr->toxels1.end(), metaballcomparer());
			voxel* tvox1=*(tr->toxels1.begin());
			if (tvox1->ball)
			{
				medialBall* vbi = tvox1->ball;
				vbi->type = 6;

				medialBall* mvi=vbi->mastrVxl();
				if (mvi!=NULL && mvi!=vbi && tr->e1 != VElems(mvi->fi+1+_pp5, mvi->fj+1+_pp5, mvi->fk+1+_pp5))		cout<<" Dmb1rrr  "<< VElems(mvi->fi+1+_pp5, mvi->fj+1+_pp5, mvi->fk+1+_pp5) <<"   ";
			}
		  else
		  {
            tvox1->ball = new medialBall(tvox1,16); throadAdditBalls.push_back(tvox1->ball);
            refs->moveUphill(tvox1->ball);
		  }
	  }


	}

	cout<<"."<<endl;
  }







}








#include "blockNet_write_cnm.cpp"
#include "blockNet_vxlManip.cpp"
