#ifndef ImageIO_H
#define ImageIO_H


#include "blockNet.h"








voxelField<float> ballRadiiToVoxel(const blockNetwork& mpn)
{
	voxelField<float> vfild(mpn.cg.nx,mpn.cg.ny,mpn.cg.nz,0);
	{
		const medialSurface & rf = *mpn.rfs;

		{///. remove voxels with NULL balls
			float maxr=0.0;
			float minr=0.0;
			std::vector<voxel>::const_iterator ti = rf.vxlSpace.begin()-1;
			std::vector<voxel>::const_iterator tend = rf.vxlSpace.end();
			while (++ti<tend)
			{
				{
					vfild(ti->i,ti->j,ti->k) = ti->R;
						 maxr=max(maxr,ti->R);
						 minr=min(minr,ti->R);
				}
			}
			cout<< "  radius  [" << minr <<"  " << maxr <<"]   nVs:" << rf.vxlSpace.size() <<endl;
		}
	}

 	return vfild;
}



voxelField<int> VElemsPlusThroats(const blockNetwork& mpn)
{
	voxelField<int> vfild = mpn.VElems;
 	cout<< " VElemsPlusThroats   "<<endl   ;
	register int throatValue=mpn.poreIs.size()+1000000;
	std::vector<throatElementI*>::const_iterator ti = mpn.throatIs.begin();
	std::vector<throatElementI*>::const_iterator tend = mpn.throatIs.end();
	while (ti<tend)
	{
		register std::vector<voxel*>::iterator  vitr = (*ti)->toxels2.begin();
		for ( ; vitr<(*ti)->toxels2.end();++vitr)
			vfild((*vitr)->i+1,(*vitr)->j+1,(*vitr)->k+1) = throatValue;
		++ti;
	}


 	return vfild;
}




///. Written by Tom Bultreys {

voxelField<int> VThroats(const blockNetwork& mpn)
{
    voxelField<int> vfild (mpn.VElems.size3() , 0);
    cout<< " VThroats   "<<endl   ;
    //register int throatValue=mpn.poreIs.size()+1000000;
    std::vector<throatElementI*>::const_iterator ti = mpn.throatIs.begin();
    std::vector<throatElementI*>::const_iterator tend = mpn.throatIs.end();
    while (ti<tend)
    {
        int throatValue = (*ti)->tid + 1;
        register std::vector<voxel*>::iterator  vitr = (*ti)->toxels2.begin();
        for ( ; vitr<(*ti)->toxels2.end();++vitr)
            vfild((*vitr)->i+1,(*vitr)->j+1,(*vitr)->k+1) = throatValue;
        ++ti;
    }


    return vfild;
}


voxelField<int> VThroats(const blockNetwork& mpn, int beginSlice, int endSlice)
{
    int x,y,zdummy;
    mpn.VElems.getSize(x,y,zdummy);
    voxelField<int> vfild ({{x,y,endSlice - beginSlice}} , 0);
    cout<< " VThroats   "<<endl   ;
    //register int throatValue=mpn.poreIs.size()+1000000;
    std::vector<throatElementI*>::const_iterator ti = mpn.throatIs.begin();
    std::vector<throatElementI*>::const_iterator tend = mpn.throatIs.end();
    while (ti<tend)
    {
        int throatValue = (*ti)->tid + 1;
        register std::vector<voxel*>::iterator  vitr = (*ti)->toxels2.begin();
        for ( ; vitr<(*ti)->toxels2.end();++vitr)
            if( (*vitr)->k+1 >= beginSlice && (*vitr)->k+1 < endSlice)
                vfild((*vitr)->i+1,(*vitr)->j+1,(*vitr)->k + 1 - beginSlice) = throatValue;
        ++ti;
    }


    return vfild;
}
 voxelField<int> poreMaxBalls(const blockNetwork& mpn)
 {
    voxelField<int> vfild(mpn.cg.nx, mpn.cg.ny, mpn.cg.nz,0);
    {
        const medialSurface & rf = *mpn.rfs;

        std::vector<medialBall>::const_iterator vi = rf.ballSpace.begin();
        std::vector<medialBall>::const_iterator voxend = rf.ballSpace.end();
        while (vi<voxend)
        {
           int x,y,z;
           x = vi->fi+_pp5;
           y = vi->fj+_pp5;
           z = vi->fk+_pp5;
             {

               float rlim = vi->R*vi->R;//

               /// absorb 2R range balls
               int ex, ey, ez;
               ex = 1*sqrt(rlim);
               for (int a = -ex; a <=  ex; ++a)
               { ey = sqrt(1*rlim-a*a);
                 for (int b = -ey; b <=  ey; ++b)
                 { ez = sqrt(1*rlim-a*a-b*b);
                   for (int c = -ez; c <=  ez; ++c)
                   {
                       if (rf.isInside(x+a, y+b, z+c) && vfild(x+a,y+b,z+c) == 0 && vi->level()==1)
                      // vfild(x+a,y+b,z+c) = vi->level();
                           vfild(x+a,y+b,z+c) = mpn.VElems(x,y,z);

                   }
                 }
               }
             }
             vi++;
        }

  }
    return vfild;
}


 voxelField<int> poreMaxBalls(const blockNetwork& mpn, int firstSlice, int lastSlice)
 {
    voxelField<int> vfild(mpn.cg.nx, mpn.cg.ny, lastSlice - firstSlice, 0);
    {
        const medialSurface & rf = *mpn.rfs;

        std::vector<medialBall>::const_iterator vi = rf.ballSpace.begin();
        std::vector<medialBall>::const_iterator voxend = rf.ballSpace.end();
        while (vi<voxend)
        {
           int x,y,z;
           x = vi->fi+_pp5;
           y = vi->fj+_pp5;
           z = vi->fk+_pp5;
             {

               float rlim = vi->R*vi->R;//

               /// absorb 2R range balls
               int ex, ey, ez;
               ex = 1*sqrt(rlim);
              // if( (x + (1 + ex) < firstSlice) || (x - ( 1 + ex) > lastSlice) )
              //     continue;
               for (int a = -ex; a <=  ex; ++a)
               { ey = sqrt(1*rlim-a*a);
                 for (int b = -ey; b <=  ey; ++b)
                 { ez = sqrt(1*rlim-a*a-b*b);
                   for (int c = -ez; c <=  ez; ++c)
                   {
                       if (rf.isInside(x+a, y+b, z+c) && (z+c >= firstSlice && z+c < lastSlice) && vfild(x+a,y+b,(z-firstSlice)+c) == 0 && vi->level()==1)
                      // vfild(x+a,y+b,z+c) = vi->level();
                           vfild(x+a,y+b,z-firstSlice+c) = mpn.VElems(x,y,z);

                   }
                 }
               }
             }
             vi++;
        }

  }
    return vfild;
}

 voxelField<int> throatMaxBalls(const blockNetwork& mpn)
 {
    voxelField<int> vfild(mpn.cg.nx,mpn.cg.ny,mpn.cg.nz,0);
    {
        const medialSurface & rf = *mpn.rfs;

        std::vector<throatElementI *>::const_iterator thr = mpn.throatIs.begin();
        std::vector<throatElementI *>::const_iterator thrEnd = mpn.throatIs.end();
        while (thr<thrEnd)
        {
           int x,y,z;
          const medialBall* vi1 = (*thr)->mb11();
          const medialBall* vi2 = (*thr)->mb22();
          const medialBall* vi;

          //select largest inscribed sphere
          if(vi1 && vi2){
              if(vi1->R > vi2->R)
                  vi = (*thr)->mb11();
              else
                  vi = (*thr)->mb22();
            }
          else if (vi1)
              vi = (*thr)->mb11();
          else
              vi = (*thr)->mb22();

           x = vi->fi +_pp5;
           y = vi->fj +_pp5;
           z = vi->fk +_pp5;

             {

               float rlim = vi->R*vi->R;//

               /// absorb 2R range balls
               int ex, ey, ez;
               ex = 1*sqrt(rlim);
               for (int a = -ex; a <=  ex; ++a)
               { ey = sqrt(1*rlim-a*a);
                 for (int b = -ey; b <=  ey; ++b)
                 { ez = sqrt(1*rlim-a*a-b*b);
                   for (int c = -ez; c <=  ez; ++c)
                   {
                       if (rf.isInside(x+a, y+b, z+c) && vfild(x+a,y+b,z+c) == 0){
                           int val = (*thr)->tid + 1;
                           vfild(x+a,y+b,z+c) = val;
                       }

                   }
                 }
               }
             }
             thr++;
        }

  }
    return vfild;
}


 voxelField<int> throatMaxBalls(const blockNetwork& mpn, int firstSlice, int lastSlice)
 {
    voxelField<int> vfild(mpn.cg.nx,mpn.cg.ny,lastSlice - firstSlice,0);
    {
        const medialSurface & rf = *mpn.rfs;

        std::vector<throatElementI *>::const_iterator thr = mpn.throatIs.begin();
        std::vector<throatElementI *>::const_iterator thrEnd = mpn.throatIs.end();
        while (thr<thrEnd)
        {
           int x,y,z;
          const medialBall* vi1 = (*thr)->mb11();
          const medialBall* vi2 = (*thr)->mb22();
          const medialBall* vi;

          //select largest inscribed sphere
          if(vi1 && vi2){
              if(vi1->R > vi2->R)
                  vi = (*thr)->mb11();
              else
                  vi = (*thr)->mb22();
            }
          else if (vi1)
              vi = (*thr)->mb11();
          else
              vi = (*thr)->mb22();

           x = vi->fi +_pp5;
           y = vi->fj +_pp5;
           z = vi->fk +_pp5;

             {

               float rlim = vi->R*vi->R;//

               /// absorb 2R range balls
               int ex, ey, ez;
               ex = 1*sqrt(rlim);
               for (int a = -ex; a <=  ex; ++a)
               { ey = sqrt(1*rlim-a*a);
                 for (int b = -ey; b <=  ey; ++b)
                 { ez = sqrt(1*rlim-a*a-b*b);
                   for (int c = -ez; c <=  ez; ++c)
                   {
                       if (rf.isInside(x+a, y+b, z+c) && (z+c >= firstSlice && z+c < lastSlice) && vfild(x+a,y+b,z-firstSlice+c) == 0){
                           int val = (*thr)->tid + 1;
                           vfild(x+a,y+b,z-firstSlice+c) = val;
                       }

                   }
                 }
               }
             }
             thr++;
        }

  }
    return vfild;
}
///. } Tom 







#endif
