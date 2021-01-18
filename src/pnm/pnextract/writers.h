#ifndef VTUWRITER_H
#define VTUWRITER_H

/*---------------------------------------------------------------------------*\
Developed by:
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



void vtuWriteMbMbs(std::string suffix, const std::vector<medialBall*>& ballSpace,  const  std::vector<poreNE*> poreIs, const voxelField<int>&  VElems, double dx, dbl3 X0);
void vtuWriteThroatMbMbs(std::string baseName,   const std::vector<throatNE*>& throatIs,   const  std::vector<poreNE*> poreIs,  const voxelField<int>&  VElems, double dx, dbl3 X0);

void vtuWriteMedialSurface(std::string suffix, const std::vector<medialBall*>& ballSpace,  const  std::vector<poreNE*> poreIs, const voxelField<int>&  VElems, double dx, dbl3 X0);

void vtuWritePores(std::string suffix,  const std::vector<poreNE*>& poreIs, const std::vector<throatNE*>& throatIs, double dx, dbl3 X0);
void vtuWriteTHroatSpheres(std::string suffix,  const std::vector<poreNE*>& poreIs, const std::vector<throatNE*>& throatIs, double dx, dbl3 X0);
void vtuWriteThroats(std::string suffix,  const std::vector<poreNE*>& poreIs, const std::vector<throatNE*>& throatIs, double dx, dbl3 X0);






#endif

