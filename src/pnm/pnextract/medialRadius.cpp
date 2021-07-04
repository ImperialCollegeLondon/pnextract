#ifndef MEDIALAXIS_H
#define MEDIALAXIS_H








 void medialSurface::calc_distmap(voxel * vit, unsigned char vValue, const voxelImage& vxls, std::vector<std::vector<node> >& oldAliens) const
 {






	const int i = vit->i, j = vit->j, k = vit->k;

	node nalien(i,j,-nz);

	int
	epxMax = 2*nx;


	double frz2 = epxMax*epxMax;
	int frz1 = 2, fry1 = 0;



	if(k>0)  {
		if(vValue!=vxls(i, j, k-1))  {
			nalien.i = i;  nalien.j = j;  nalien.k = k-1;
			frz2 = 1.;
			frz1 = 3;
			//fry1 = 0;
		}
		else
		{

			if(j>0)
			{
				if(vValue!=vxls(i, j-1, k))  {
					nalien.i = i;  nalien.j = j-1;  nalien.k = k;
					frz2 = 1.;
					frz1 = 3;
					//fry1 = 0;
				}
				else
				{
					const node & nalienOldi = oldAliens[j-1][i];

					int neilienDistSqr = (nalienOldi.i-i)*(nalienOldi.i-i)+(nalienOldi.j-j)*(nalienOldi.j-j)+(nalienOldi.k-k)*(nalienOldi.k-k);
					frz2 = neilienDistSqr+0.;
					frz1 = -sqrt(neilienDistSqr)-1;
					fry1 = nalienOldi.j-j-1;
					nalien.i = nalienOldi.i;  nalien.j = nalienOldi.j;  nalien.k = nalienOldi.k;

				}

			}

			const node & nalienOldi = oldAliens[j][i];
			int neilienDistSqr = (nalienOldi.i-i)*(nalienOldi.i-i)+(nalienOldi.j-j)*(nalienOldi.j-j)+(nalienOldi.k-k)*(nalienOldi.k-k);
			if (neilienDistSqr<frz2)  {
				nalien.i = nalienOldi.i;  nalien.j = nalienOldi.j;  nalien.k = nalienOldi.k;

				frz2 = neilienDistSqr+0.;
				frz1 = nalienOldi.k-k-1;
				if(j == 0) fry1 = -sqrt(frz2)-1;
			}

		}

	}
	else if(j>0)  {
		if(vValue!=vxls(i, j-1, k))  {
			nalien.i = i;  nalien.j = j-1;  nalien.k = k;
			frz2 = 1.;
			frz1 = 3;
			//fry1 = 0;
		}
		else
		{
			const node & nalienOldi = oldAliens[j-1][i];

			int neilienDistSqr = (nalienOldi.i-i)*(nalienOldi.i-i)+(nalienOldi.j-j)*(nalienOldi.j-j)+(nalienOldi.k-k)*(nalienOldi.k-k);
			frz2 = neilienDistSqr+0.;
			frz1 = 0;
			fry1 = nalienOldi.j-j-1;
			nalien.i = nalienOldi.i;  nalien.j = nalienOldi.j;  nalien.k = nalienOldi.k;

		}

	}
	else //if(	epxMax == 2*nx)//. XXXXXX WARNING
	{
		if( isInside(nextSegg(i, j, k).start))  {
			epxMax = nextSegg(i, j, k).start-i;
			if( isInside(segg(i, j, k).start-1) && i-(segg(i, j, k).start-1)<epxMax)  {
				epxMax = i-(segg(i, j, k).start-1);
				nalien.i = (segg(i, j, k).start-1);  nalien.j = j;  nalien.k = k;
			}
			else
			{
				nalien.i = nextSegg(i, j, k).start;  nalien.j = j;  nalien.k = k;
			}
		}
		else if( isInside(segg(i, j, k).start-1))  {
			epxMax = std::min(i-segg(i, j, k).start+1,epxMax);
			nalien.i = (segg(i, j, k).start-1);  nalien.j = j;  nalien.k = k;
		}
		else if( isInside(i))  {
			epxMax = std::min(nx,std::min(ny,nz))+1;
		}
		else
		{
			cout<<"\n\n Error: outside voxel \n\n";
		}
		frz2 = epxMax*epxMax+0.;
		frz1 = -epxMax;
		fry1 = -epxMax;
	}


	if(epxMax <= 0)  	cout<<i<<" "<<j<<" "<<k<<"    "<<segg(i, j, k).start	<<" "<<(nextSegg(i, j, k)).start<<" "<<endl;








	for (int c = std::max(frz1-1,-k); c <=  min(int(sqrt(frz2))+1,nz-k-1); ++c)	{
	  //if( isInside(i, j, k+c))  
	  {const int blim=min(int(sqrt(frz2-c*c)+1.001),ny-j-1);
		for (int b = std::max(std::max(int(-sqrt(frz2-c*c)),fry1)-1,-j); b<=blim  ; ++b)  {

		  //if(isJInside(j+b)) //
		  {
			if (vValue!=vxls(i, j+b, k+c))  {
				if ((b*b+c*c) < frz2)  {
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
					if((a*a+b*b+c*c) < frz2)  {
						frz2 = a*a+b*b+c*c;
						nalien.i = i+a;  nalien.j = j+b;  nalien.k = k+c;
					}
				}


				if((&s+1)->start<nx)
				{
					int a = ((&s+1)->start-i);
					if((a*a+b*b+c*c) < frz2)  {
						frz2 = a*a+b*b+c*c;
						nalien.i = i+a;  nalien.j = j+b;  nalien.k = k+c;
					}
				}

			}
		  }

		}}
	}



	if ( ! isInside(nalien.i,nalien.j,nalien.k))  {
			nalien.i = (i<nx/2) ? -nx/4-1  : nx*5/4+1;
			nalien.j = (j<ny/2) ? -ny/4-1  : ny*5/4+1;
			nalien.k = (k<nz/2) ? -nz/4-1  : nz*5/4+1;
			vit->R = sqrt((nalien.i-i)*(nalien.i-i) + (nalien.j-j)*(nalien.j-j) + (nalien.k-k)*(nalien.k-k)) - 0.5;
	}
	else {
		int dx=abs(nalien.i-i),  dy=abs(nalien.j-j),  dz=abs(nalien.k-k);

		double limit = sqrt(dx*dx + dy*dy + dz*dz) - 0.5;
        double iSqr=min((j+2),(ny-j+1));
		if (iSqr<limit)
            limit=max((1.-_clipROutyz)*limit+_clipROutyz*iSqr,0.01);
        iSqr=min((k+2),(nz-k+1));
		if (iSqr<limit)
            limit=max((1.-_clipROutyz)*limit+_clipROutyz*iSqr,0.01);
        iSqr=min((i+2),(nx-i+1));
		if (iSqr<limit)
            limit=max((1.-_clipROutx )*limit+_clipROutx*iSqr,0.1);
		vit->R = limit ;



		if (frz2 <= 0) cout<<"WTF frz2 = "<<frz2<<endl;
		if (nalien.i<-2000 || limit > 16000000)  {cout<<"Error i = "<<nalien.i<<endl;
			cout<<"frz2 "<<frz2<<endl;
			cout<<"frz1 "<<frz1<<endl;
			cout<<"i "<<i<<"  j "<<j<<"  k "<<k<<endl;
			cout<<"oldAliens[j][i]. i "<<oldAliens[j][i].i<<"  j "<<oldAliens[j][i].j<<"  k "<<oldAliens[j][i].k<<endl;
			exit(0);
		}
  }




		oldAliens[j][i] = nalien;

 }


voxelImage segToVxlMesh(const medialSurface & ref)  {/// converts segments back to voxelImage
	voxelImage vxls(ref.nx,ref.ny,ref.nz,255);
	for (int iz = 0; iz<ref.nz; ++iz)  {
		for (int iy = 0; iy<ref.ny; ++iy)  {
			const segments& s = ref.segs_[iz][iy];
			for (int ix = 0; ix<s.cnt; ++ix)  {
				std::fill (&vxls(s.s[ix].start,iy,iz),&vxls(s.s[ix+1].start,iy,iz),s.s[ix].value);
			}
		}
	}
	return vxls;
}



void medialSurface::calc_distmaps() //search  MBs at each voxel
{
	cout<< " computing distance map for index "<<int(0); cout.flush();

	if (!nVxls) { cout<<" no voxels no balls,\n"<<endl; return; }

	voxelImage vxls = segToVxlMesh(*this);
	double rBalls = 0.;




	OMPragma("omp parallel reduction(+:rBalls)")  {
		std::vector<std::vector<node> > oldAliens(ny+1,std::vector<node>(nx));
		for (int j=0; j<ny+1; ++j)
		for (int i=0; i<nx; ++i)  {
			oldAliens[j][i].i = i;
			oldAliens[j][i].j = j;
			oldAliens[j][i].k = -nz/2-1;
		}

		size_t nvxls10th=max(10*int(vxlSpace.size()/200),1);
		const voxel* const vnd = &*vxlSpace.end();;
		OMPragma("omp for")
		for (voxel* vit = &vxlSpace[0]; vit<vnd; ++vit)  {
			calc_distmap(&*vit, 0, vxls, oldAliens);
			if (size_t(vit)%nvxls10th == 0)  {  ( cout<< "\r  distance map / sphere radius = " <<vit->R ).flush();  }
			rBalls+=vit->R;
		}
	}
	cout<< "\n  average distance map = "<<rBalls/nVxls<<endl;

	if(_minRp<0.)  {
		setDefaults(rBalls/nVxls);
	}


	return ;
 }




void medialSurface::smoothRadius()  {

	(cout<<" smoothing R  ").flush();
	






	std::vector<float> delRrr(vxlSpace.size(),0.0f);  (cout<<"*").flush();
	OMPFor()
	for (int k = 0; k < nz; ++k) {
	 for (int j = 0; j<ny; ++j)  {
		const segments& s = cg_.segs_[k][j];
		for (int ix = 0; ix<s.cnt; ++ix)
		if (s.s[ix].value == 0)  {
		  segment& seg=s.s[ix];
		  for (int i = seg.start; i < s.s[ix+1].start; ++i)  {
			 double sumR=0.; int counter=0;
		    for (int kk = max(k-1,0); kk<min(k+2,nz); ++kk)
		      for (int jj = max(j-1,0); jj<min(j+2,ny); ++jj)  {
					int ii = max(i-1,0);
					const segment* segbc=cg_.segptr(ii,jj,kk);
					if (segbc->value != 0 && (segbc+1)->value == 0)  { ++segbc; ii=segbc->start; }
					if (segbc->value == 0)  {
						int ii2=min((segbc+1)->start,i+2);
						voxel* vxlj=segbc->segV+(ii-segbc->start);
						for ( ; ii<ii2; ++ii)  {
							sumR += vxlj->R;	++vxlj;
							counter+=1;
						}
					}
				}




				delRrr[seg.segV+(i-seg.start)-(&vxlSpace[0])] = 4.*sumR/(3*counter+27)-seg.segV[i-seg.start].R;
		  }
		}
	  }
	 }  (cout<<"*").flush();



	OMPFor()
	for (int k = 0; k < nz; ++k) {
	 for (int j = 0; j<ny; ++j)  {
		const segments& s = cg_.segs_[k][j];
		for (int ix = 0; ix<s.cnt; ++ix)
		if (s.s[ix].value == 0)  {
		  segment& seg=s.s[ix];
		  for (int i = seg.start; i < s.s[ix+1].start; ++i)  {
			 double sumDelR=0.; int counter=0;
		    for (int kk = max(k-1,0); kk<min(k+2,nz); ++kk)
		    for (int jj = max(j-1,0); jj<min(j+2,ny); ++jj)  {
				int ii = max(i-1, 0);
				const segment* segbc=cg_.segptr(ii,jj,kk);
				if (segbc->value != 0 && (segbc+1)->value == 0)		{++segbc; ii=segbc->start;}
				if (segbc->value == 0)  {
					int ii2=min((segbc+1)->start,i+2);
					voxel* vxlj=segbc->segV+(ii-segbc->start);
					for ( ; ii<ii2; ++ii)  {
						sumDelR += delRrr[vxlj-(&vxlSpace[0])];	++vxlj; ++counter;
					}
				}
			  }



			  seg.segV[i-seg.start].R  += min(max(0.02* (delRrr[seg.segV+(i-seg.start)-(&vxlSpace[0])] - 0.99*2.*sumDelR/(1*counter+27)),-0.005),0.01);
		   }
		} }
	 }  (cout<<"*").flush();

	{/// Finally, report max distance map, to confirm that distance map is not changed too much
		float maxrrr=0;
		OMPragma("omp parallel for reduction(max:maxrrr)")
		for(auto ti = vxlSpace.begin(); ti<vxlSpace.end(); ++ti) maxrrr = max(maxrrr,ti->R);
		cout<< " maxrrr " << maxrrr << endl;
	}

}


#endif
