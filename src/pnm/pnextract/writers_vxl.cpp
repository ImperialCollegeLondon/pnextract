
#include "blockNet.h"

///@cond INTERNAL


voxelField<float> ballRadiiToVoxel(const blockNetwork& mpn)  {
	cout<< " write_radius   "<<endl;
	voxelField<float> vfild(mpn.cg.nx,mpn.cg.ny,mpn.cg.nz,0);
	{
		const medialSurface & rf = *mpn.srf;

		float maxr=0.,  minr=0.;
		for(const auto& vi: rf.vxlSpace)  {
			vfild(vi.i,vi.j,vi.k) = vi.R;
			maxr=max(maxr,vi.R);
			minr=min(minr,vi.R);
		}
		cout<< "  radius  [" << minr <<"  " << maxr <<"]   nVs:" << rf.vxlSpace.size() <<endl;
	}
	return vfild;
}



/*//voxelImageT<int> VElemsPlusThroats(const blockNetwork& mpn)  {
	//cout<< " VElemsPlusThroats   "<<endl;
	//voxelImageT<int> vfild(mpn.VElems);
	//int throatValue=mpn.poreIs.size()+1000000;
	//for(const auto tr: mpn.throatIs)
		//for (const auto  vi: tr->toxels2)
			//vfild(vi->i+1,vi->j+1,vi->k+1) = throatValue;
	//return vfild;
//}*/





///.  Originally written by Tom Bultreys {

voxelField<int> VThroats(const blockNetwork& mpn)  {
	cout<< " write_throats   "<<endl;
	voxelField<int> vfild (mpn.VElems.size3() , 0);
	for(const auto tr: mpn.throatIs)  {
		int throatValue = tr->tid + 1;
		for (const auto  vi: tr->toxels2)
			vfild(vi->i+1,vi->j+1,vi->k+1) = throatValue;
	}
	return vfild;
}

voxelField<int> VThroats(const blockNetwork& mpn, int beginSlice, int endSlice)  {
	cout<< " write_throats   "<<endl;
   int x,y,zdummy;
	mpn.VElems.getSize(x,y,zdummy);
	voxelField<int> vfild (int3(x,y,endSlice - beginSlice) , 0);
	for(const auto tr: mpn.throatIs)  {
		int throatValue = tr->tid + 1;
		for (const auto  vi: tr->toxels2)
			if( vi->k+1 >= beginSlice && vi->k+1 < endSlice)
				vfild(vi->i+1, vi->j+1, vi->k+1-beginSlice) = throatValue;
	}
	return vfild;
}



voxelField<int> poreMaxBalls(const blockNetwork& mpn)  {
	cout<< " write_poreMaxBalls   "<<endl;
	voxelField<int> vfild(mpn.cg.nx, mpn.cg.ny, mpn.cg.nz,0);

	const medialSurface & rf = *mpn.srf;

	for(const medialBall& vi: rf.ballSpace)  if(vi.level()==1)  {
		int x= vi.fi,   y= vi.fj,   z= vi.fk;
		float rlim = vi.R*vi.R;

		/// absorb 2R range balls
		int ex, ey, ez;
		ex = 1*sqrt(rlim);
		for (int a=-ex; a<=ex; ++a)  {
		  ey = sqrt(1*rlim-a*a);
		  for (int b=-ey; b<=ey; ++b)  {
			 ez = sqrt(1*rlim-a*a-b*b);
			 for (int c=-ez; c<=ez; ++c)
				  if (rf.isInside(x+a, y+b, z+c) && vfild(x+a,y+b,z+c) == 0)
						vfild(x+a,y+b,z+c) = mpn.VElems(x,y,z);
		  }
		}
	}
	return vfild;
}

voxelField<int> poreMaxBalls(const blockNetwork& mpn, int firstSlice, int lastSlice)  {
	cout<< " write_poreMaxBalls   "<<endl;

	voxelField<int> vfild(mpn.cg.nx, mpn.cg.ny, lastSlice - firstSlice, 0);

	const medialSurface & rf = *mpn.srf;

	for(const medialBall& vi: rf.ballSpace)  if (vi.level()==1)  {
		int x= vi.fi,   y= vi.fj,   z= vi.fk;
		float rlim = vi.R*vi.R;

		/// absorb 2R range balls
		int ex, ey, ez;
		ex = 1*sqrt(rlim);
		// if( (x + (1 + ex) < firstSlice) || (x - ( 1 + ex) > lastSlice) )  continue;
		for (int a=-ex; a<=ex; ++a)  {
		  ey = sqrt(1*rlim-a*a);
		  for (int b=-ey; b<=ey; ++b)  {
			 ez = sqrt(1*rlim-a*a-b*b);
			 for (int c=-ez; c<=ez; ++c)
				if (rf.isInside(x+a, y+b, z+c) && (z+c >= firstSlice && z+c < lastSlice) 
						&& vfild(x+a,y+b,(z-firstSlice)+c) == 0)
					vfild(x+a,y+b,z-firstSlice+c) = mpn.VElems(x,y,z);

		  }
		}
	}
	return vfild;
}


voxelField<int> throatMaxBalls(const blockNetwork& mpn)  {
	cout<< " write_throatMaxBalls   "<<endl;
	voxelField<int> vfild(mpn.cg.nx,mpn.cg.ny,mpn.cg.nz,0);

	const medialSurface & rf = *mpn.srf;

	for(const auto tr: mpn.throatIs)  {
		const medialBall* vi1 = tr->mb11(),  *vi2 = tr->mb22(), *vi;

		if(vi1 && vi2)  {	// select largest inscribed sphere
			if(vi1->R > vi2->R)  vi = tr->mb11();
			else                 vi = tr->mb22();  }
		else if (vi1)          vi = tr->mb11();
		else                   vi = tr->mb22();

		int x= vi->fi,   y= vi->fj,   z= vi->fk;
		float rlim = vi->R*vi->R;

		/// absorb 2R range balls
		int ex, ey, ez;
		ex = 1*sqrt(rlim);
		for (int a=-ex; a<=ex; ++a)  {
		 ey = sqrt(1*rlim-a*a);
		 for (int b=-ey; b<=ey; ++b)  {
			ez = sqrt(1*rlim-a*a-b*b);
			for (int c=-ez; c<=ez; ++c)
				if (rf.isInside(x+a, y+b, z+c) && vfild(x+a,y+b,z+c) == 0)
					vfild(x+a,y+b,z+c) = tr->tid + 1;
		 }
		}
	}
	return vfild;
}


 voxelField<int> throatMaxBalls(const blockNetwork& mpn, int firstSlice, int lastSlice)  {
	cout<< " write_throatMaxBalls   "<<endl;
	voxelField<int> vfild(mpn.cg.nx,mpn.cg.ny,lastSlice - firstSlice,0);
	const medialSurface & rf = *mpn.srf;

	for(const auto tr: mpn.throatIs)  {
		const medialBall* vi1 = tr->mb11(),  *vi2 = tr->mb22(),  *vi;

		if(vi1 && vi2)  {		//select largest inscribed sphere
			 if(vi1->R > vi2->R)   vi = tr->mb11();
			 else                  vi = tr->mb22();  }
		else if (vi1)     vi = tr->mb11();
		else              vi = tr->mb22();

		int x= vi->fi,   y= vi->fj,   z= vi->fk;
		float rlim = vi->R*vi->R;

		/// absorb 2R range balls
		int ex, ey, ez;
		ex = 1*sqrt(rlim);
		for (int a=-ex; a<=ex; ++a)  {
			ey = sqrt(1*rlim-a*a);
			for (int b=-ey; b<=ey; ++b)  {
				ez = sqrt(1*rlim-a*a-b*b);
				for (int c=-ez; c<=ez; ++c)
					if (rf.isInside(x+a, y+b, z+c) && (z+c >= firstSlice && 
							z+c < lastSlice) && vfild(x+a,y+b,z-firstSlice+c) == 0)
						vfild(x+a,y+b,z-firstSlice+c) = tr->tid + 1;
			 }
		}
	}
	return vfild;
}
///. } Tom


///@endcond

