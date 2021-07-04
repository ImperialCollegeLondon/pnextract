#ifndef VTUWRITER_H
#define VTUWRITER_H

/*---------------------------------------------------------------------------*\
written by:
	Ali Q Raeini  email: a.q.raeini@imperial.ac.uk  and
	Tom Bultreys
\*---------------------------------------------------------------------------*/



#include <sstream>
#include "voxelImage.h"
#include "ElementGNE.h"








voxelImage segToVxlMesh(const medialSurface & ref);
voxelField<float> ballRadiiToVoxel(const blockNetwork& mpn);
voxelField<int> VElemsPlusThroats(const blockNetwork& mpn);

///- written by Tom Bultreys:
 voxelField<int> VThroats(const blockNetwork& mpn);
 voxelField<int> VThroats(const blockNetwork& mpn, int firstSlice, int lastSlice);
 voxelField<int> poreMaxBalls(const blockNetwork& mpn);
 voxelField<int> poreMaxBalls(const blockNetwork& mpn, int firstSlice, int lastSlice);
 voxelField<int> throatMaxBalls(const blockNetwork& mpn);
 voxelField<int> throatMaxBalls(const blockNetwork& mpn, int firstSlice, int lastSlice);



void vtuWriteMbMbs(std::string bNam, const stvec<medialBall>& ballSpace,  const stvec<poreNE*> poreIs, const voxelField<int>&  VElems, double dx, dbl3 X0);
void vtuWriteThroatMbMbs(std::string baseNam,   const stvec<throatNE*>& throatIs,   const stvec<poreNE*> poreIs,  const voxelField<int>&  VElems, double dx, dbl3 X0);


void vtuWritePores(std::string bNam,  const stvec<poreNE*>& poreIs, const stvec<throatNE*>& throatIs, double dx, dbl3 X0);
void vtuWriteTHroatSpheres(std::string bNam,  const stvec<poreNE*>& poreIs, const stvec<throatNE*>& throatIs, double dx, dbl3 X0);
void vtuWriteThroats(std::string bNam,  const stvec<poreNE*>& poreIs, const stvec<throatNE*>& throatIs, double dx, dbl3 X0);






#endif

