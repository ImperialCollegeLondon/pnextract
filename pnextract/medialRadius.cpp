#ifndef MEDIALAXIS_H
#define MEDIALAXIS_H





inline bool operator == (const node& a, const node & b)
{
	return (a.i == b.i && a.j == b.j &&  a.k == b.k);
}



 void medialSurface::calc_distmap(voxel * vit, unsigned char vValue, const voxelImage& vxls, std::vector<std::vector<node> >& oldAliens) const
 {

	const int i = vit->i, j = vit->j, k = vit->k;

	node nalien(i,j,-nz);

	int
	epxMax = 2*nx;//. XXXXXX WARNING

	double frz2 = epxMax*epxMax;
	int frz1 = 2, fry1 = 0;

	if(k>0)  ///. z search bounds
	{
		if(vValue!=vxls(i, j, k-1))
		{
			nalien.i = i;  nalien.j = j;  nalien.k = k-1;
			frz2 = 1.0;
			frz1 = 3;
			//fry1 = 0;
		}
		else
		{

			if(j>0)  ///. y search bounds
			{
				if(vValue!=vxls(i, j-1, k))
				{
					nalien.i = i;  nalien.j = j-1;  nalien.k = k;
					frz2 = 1.0;
					frz1 = 3;
					//fry1 = 0;
				}
				else
				{
					const node & nalienOldi = oldAliens[j-1][i];

					int neilienDistSqr = (nalienOldi.i-i)*(nalienOldi.i-i)+(nalienOldi.j-j)*(nalienOldi.j-j)+(nalienOldi.k-k)*(nalienOldi.k-k);
					frz2 = neilienDistSqr+0.0;
					frz1 = -sqrt(neilienDistSqr)-1;
					fry1 = nalienOldi.j-j-1;
					nalien.i = nalienOldi.i;  nalien.j = nalienOldi.j;  nalien.k = nalienOldi.k;

				}

			}

			const node & nalienOldi = oldAliens[j][i];
			int neilienDistSqr = (nalienOldi.i-i)*(nalienOldi.i-i)+(nalienOldi.j-j)*(nalienOldi.j-j)+(nalienOldi.k-k)*(nalienOldi.k-k);
			if (neilienDistSqr<frz2)
			{
				nalien.i = nalienOldi.i;  nalien.j = nalienOldi.j;  nalien.k = nalienOldi.k;

				frz2 = neilienDistSqr+0.0;
				frz1 = nalienOldi.k-k-1;
				if(j == 0) fry1 = -sqrt(frz2)-1;
			}

		}

	}
	else if(j>0)  ///. y search bounds
	{
		if(vValue!=vxls(i, j-1, k))
		{
			nalien.i = i;  nalien.j = j-1;  nalien.k = k;
			frz2 = 1.0;
			frz1 = 3;
			//fry1 = 0;
		}
		else
		{
			const node & nalienOldi = oldAliens[j-1][i];

			int neilienDistSqr = (nalienOldi.i-i)*(nalienOldi.i-i)+(nalienOldi.j-j)*(nalienOldi.j-j)+(nalienOldi.k-k)*(nalienOldi.k-k);
			frz2 = neilienDistSqr+0.0;
			frz1 = 0;
			fry1 = nalienOldi.j-j-1;
			nalien.i = nalienOldi.i;  nalien.j = nalienOldi.j;  nalien.k = nalienOldi.k;

		}

	}
	else //if(	epxMax == 2*nx)//. XXXXXX WARNING
	{
		if( isInside(nextSegg(i, j, k).start) )
		{
			epxMax = nextSegg(i, j, k).start-i;
			if( isInside(segg(i, j, k).start-1) && i-(segg(i, j, k).start-1)<epxMax )
			{
				epxMax = i-(segg(i, j, k).start-1);
				nalien.i = (segg(i, j, k).start-1);  nalien.j = j;  nalien.k = k;
			}
			else
			{
				nalien.i = nextSegg(i, j, k).start;  nalien.j = j;  nalien.k = k;
			}
		}
		else if( isInside(segg(i, j, k).start-1) )
		{
			epxMax = std::min(i-segg(i, j, k).start+1,epxMax);
			nalien.i = (segg(i, j, k).start-1);  nalien.j = j;  nalien.k = k;
		}
		else if( isInside(i) )
		{
			epxMax = std::min(nx,std::min(ny,nz))+1;
		}
		else
		{
			cout<<"\n\n Err \n\n";

		}
		frz2 = epxMax*epxMax+0.0;
		frz1 = -epxMax;
		fry1 = -epxMax;
	}


	if(epxMax <= 0)  	cout<<i<<" "<<j<<" "<<k<<"    "<<segg(i, j, k).start	<<" "<<(nextSegg(i, j, k)).start<<" "<<endl;






	for (int c = std::max(frz1-1,-k); c <=  min(int(sqrt(frz2))+1,nz-k-1); ++c)
	{
	  //if( isInside(i, j, k+c))
	  {const int blim=min(int(sqrt(frz2-c*c)+1.001),ny-j-1);
		for (int b = std::max(std::max(int(-sqrt(frz2-c*c)),fry1)-1,-j); b<=blim   ; ++b)
	   {

		  //if(isJInside(j+b))
		  {
			if (vValue!=vxls(i, j+b, k+c)) // find a failure voxel
			{
				if ((b*b+c*c) < frz2)
				{
					frz2 = b*b+c*c;
					nalien.i = i;  nalien.j = j+b;  nalien.k = k+c;
				}

			}
			else
			{
				const segment& s = segg(i, j+b, k+c);
				if(s.start>0)
				{
					int a = (s.start-1-i);
					if((a*a+b*b+c*c) < frz2 )
					{
						frz2 = a*a+b*b+c*c;
						nalien.i = i+a;  nalien.j = j+b;  nalien.k = k+c;
					}
				}


				if((&s+1)->start<nx)
				{
					int a = ((&s+1)->start-i);
					if((a*a+b*b+c*c) < frz2 )
					{
						frz2 = a*a+b*b+c*c;
						nalien.i = i+a;  nalien.j = j+b;  nalien.k = k+c;
					}
				}

			}
		  }

		}}
	}


	if ( ! isInside(nalien.i,nalien.j,nalien.k))
	{

			nalien.i = (i<nx/2) ? -nx/4-1  : nx*5/4+1;
			nalien.j = (j<ny/2) ? -ny/4-1  : ny*5/4+1;
			nalien.k = (k<nz/2) ? -nz/4-1  : nz*5/4+1;
			vit->R = sqrt((nalien.i-i)*(nalien.i-i) + (nalien.j-j)*(nalien.j-j) + (nalien.k-k)*(nalien.k-k)) - 0.5;

	}
	else
	{

		short dx=abs(nalien.i-i),  dy=abs(nalien.j-j),  dz=abs(nalien.k-k);

		double limit = sqrt(dx*dx + dy*dy + dz*dz) - 0.5;
        double iSqr=min((j+2),(ny-j+1)); ///. allow allways 2 extra void voxels around the image
		if (iSqr<limit)
            limit=max((1.0-_clipROutyz)*limit+_clipROutyz*iSqr,0.01);
        iSqr=min((k+2),(nz-k+1));
		if (iSqr<limit)
            limit=max((1.0-_clipROutyz)*limit+_clipROutyz*iSqr,0.01);
        iSqr=min((i+2),(nx-i+1));
		if (iSqr<limit)
            limit=max((1.0-_clipROutx )*limit+_clipROutx*iSqr,0.1);
		vit->R = limit ;


//		double limit = dx*dx + dy*dy + dz*dz; limit-=sqrt(limit)-0.25;
//        double iSqr=min((j+1),(ny-j)); iSqr*=iSqr; ///. allow allways 0.5 extra void voxels around the image
//		if (iSqr<limit)
//            limit=max((1.0-_clipROutyz)*limit+_clipROutyz*iSqr,0.1);
//        iSqr=min((k+1),(nz-k)); iSqr*=iSqr;
//		if (iSqr<limit)
//            limit=max((1.0-_clipROutyz)*limit+_clipROutyz*iSqr,0.1);
//        iSqr=min((i+1),(nx-i)); iSqr*=iSqr;
//		if (iSqr<limit)
//            limit=max((1.0-_clipROutx )*limit+_clipROutx*iSqr,0.1);
//		vit->R = sqrt(limit);

//		if (iSqr<limit)
//            limit=max((1.0-_clipROutx )*(limit-1.0)+_clipROutx*iSqr,0.1);
//		else if (dx*dx<limit) limit+=_clipROutx*(dx*dx-limit+1);
//		if (jSqr<limit)        limit+=_clipROutyz*(jSqr-limit+1);
//		else if (dy*dy<limit) limit+=_clipROutyz*(dy*dy-limit+1);
//		if (kSqr<limit)        limit+=_clipROutyz*(kSqr-limit+1);
//		else if (dz*dz<limit) limit+=_clipROutyz*(dz*dz-limit+1);


		if (frz2 <= 0) cout<<"WTF frz2 = "<<frz2<<endl;
		if (nalien.i<-2000 || limit > 16000000)
		{cout<<"Error i = "<<nalien.i<<endl;
			cout<<"frz2 "<<frz2<<endl;
			cout<<"frz1 "<<frz1<<endl;
			cout<<"i "<<i<<"  j "<<j<<"  k "<<k<<endl;
			cout<<"oldAliens[j][i]. i "<<oldAliens[j][i].i<<"  j "<<oldAliens[j][i].j<<"  k "<<oldAliens[j][i].k<<endl;
			exit(0);
		}
  }




		oldAliens[j][i] = nalien;

 }


voxelImage segToVxlMesh(const medialSurface & ref)
{/// converts segments back to voxelImage
	voxelImage vxls(ref.nx,ref.ny,ref.nz,255);
 	for (int iz = 0; iz<ref.nz; ++iz)
 	{
 		for (int iy = 0; iy<ref.ny; ++iy)
 		{
 			const segments& s = ref.segs_[iz][iy];
 			for (int ix = 0; ix<s.cnt; ++ix)
 			{
				std::fill (&vxls(s.s[ix].start,iy,iz),&vxls(s.s[ix+1].start,iy,iz),s.s[ix].value);
 			}
 		}
 	}
 	return vxls;
}



void medialSurface::calc_distmaps() //search  MBs at each voxel
{
	cout<< " computing distance map for index "<<int(0); cout.flush();


	nBalls = 0;
	if (!nVxls) { cout<<" no voxels no balls,\n"<<endl; return; }

	voxelImage vxls = segToVxlMesh(*this);
	double rBalls = 0.0;

	std::vector<std::vector<node> > oldAliens(ny+1,std::vector<node>(nx));
	for (int j = 0; j<ny+1; ++j)
	for (int i = 0; i<nx; ++i)
	{
		oldAliens[j][i].i = i;
		oldAliens[j][i].j = j;
		oldAliens[j][i].k = -nz/2-1;
	}

	size_t nvxls10th=pow(10,int(log10(10000+vxlSpace.size()/100)));
	std::vector<voxel>::iterator vit = vxlSpace.begin();
	for (size_t i = 0; i<nVxls; ++i)
	{
		calc_distmap(&*vit, 0, vxls, oldAliens);
		if (i%nvxls10th == 0)	{cout<< "\r  found " << i << " balls  radius = " <<vit->R ;cout.flush();}
		if (vit->R >= _minRp) {vit->ball = &(ToBeAssigned);	++nBalls; rBalls+=vit->R;}
		else 	    {vit->ball = NULL;}
		++vit;
	}
	cout<< "\n  created " << nBalls << " balls,  average radius = "<<rBalls/nBalls<<endl;
 	return ;
 }




void medialSurface::smoothRadius()
{/// smooth distance map, optional 
	cout<<" smoothing R, "; cout.flush();
	
	/// find Gaussian smoothed dR
	std::vector<double> delRrr(vxlSpace.size(),0.0);
 	for (short k = 0; k < nz; ++k)
	 for (short j = 0; j < ny; ++j)
 	 {
		const segments& s = cg_.segs_[k][j];
		for (short ix = 0; ix<s.cnt; ++ix)
		if (s.s[ix].value == 0)
		{/////////////
		  segment& seg=s.s[ix];
		  for (short i = seg.start; i < s.s[ix+1].start; ++i)
		  {
			 double sumR=0.0; int counter=0;
		    for ( int kk = max(k-1,0); kk<min(k+2,nz) ; ++kk )
		    for ( int jj = max(j-1,0); jj<min(j+2,ny) ; ++jj )
		    { int ii = max(i-1,0);
				const segment* segbc=cg_.segptr(ii,jj,kk);
				if (segbc->value != 0 && (segbc+1)->value == 0)		{++segbc; ii=segbc->start;}
				if (segbc->value == 0)
				{
					int ii2=min((segbc+1)->start,i+2);
					voxel* vxlj=segbc->segV+(ii-segbc->start);
					for ( ; ii<ii2 ; ++ii )
					{
						sumR += vxlj->R;	++vxlj;
						counter+=1;
					}
				}
			  }
			  delRrr[seg.segV+(i-seg.start)-(&vxlSpace[0])] = 4.0*sumR/(3*counter+27)-seg.segV[i-seg.start].R;
		   }
		}///////////////
	 }

	/// relax and assign R
 	for (short k = 0; k < nz; ++k)
	 for (short j = 0; j < ny; ++j)
 	 {
		const segments& s = cg_.segs_[k][j];
		for (short ix = 0; ix<s.cnt; ++ix)
		if (s.s[ix].value == 0)
		{/////////////
		  segment& seg=s.s[ix];
		  for (short i = seg.start; i < s.s[ix+1].start; ++i)
		  {
			 double sumDelR=0.0; int counter=0;
		    for ( int kk = max(k-1,0); kk<min(k+2,nz) ; ++kk )
		    for ( int jj = max(j-1,0); jj<min(j+2,ny) ; ++jj )
		    { int ii = max(i-1,0);
				const segment* segbc=cg_.segptr(ii,jj,kk);
				if (segbc->value != 0 && (segbc+1)->value == 0)		{++segbc; ii=segbc->start;}
				if (segbc->value == 0)
				{
					int ii2=min((segbc+1)->start,i+2);
					voxel* vxlj=segbc->segV+(ii-segbc->start);
					for ( ; ii<ii2 ; ++ii )
					{
						sumDelR += delRrr[vxlj-(&vxlSpace[0])];	++vxlj;++counter;
					}
				}
			  }
			  
			  //delRrr[seg.segV+(i-seg.start)-(&vxlSpace[0])] -= 1.2*sumR/counter;
			  seg.segV[i-seg.start].R  += min(max(0.02* (delRrr[seg.segV+(i-seg.start)-(&vxlSpace[0])] - 0.99*2.0*sumDelR/(1*counter+27)),-0.005),0.01);
		   }
		}///////////////
	 }


	{/// remove maximal-spheres with small radii
		std::vector<voxel>::iterator vi = vxlSpace.begin()-1;
		std::vector<voxel>::iterator vend = vxlSpace.end();
		while(++vi<vend)  if(vi->ball && vi->R<_minRp) { vi->ball=NULL;  --nBalls; }
	}

	{/// report max distance map
		float maxrrr=0;
		std::vector<voxel>::iterator ti = vxlSpace.begin()-1;
		std::vector<voxel>::iterator tend = vxlSpace.end();
		while (++ti<tend) maxrrr = max(maxrrr,ti->R);
		cout<< " maxrrr " << maxrrr << endl;
	}

}


#endif
