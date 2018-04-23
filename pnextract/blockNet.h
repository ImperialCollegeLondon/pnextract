#ifndef POROUSBLOCK_H
#define POROUSBLOCK_H

///. nAprox is the approximation order for discretising continuous functions


#include <map>
#include <fstream>
#include "voxelImage.h"
#include "inputData.h"
#include "medialSurf.h"
//~ #include "util.h"


#define RLEVEL3 0.7
#define SubElem0 nElems

int nextract(inputDataNE& cfg);

class blockNetwork
{
public:
 blockNetwork(medialSurface*& rs, const inputDataNE& c) : rfs(rs), cg(c), maxNCors(5) {};




	void createBallsAndHierarchy(medialSurface*& refs, inputDataNE& cfg, size_t startValue);
	void collectAllballs(medialSurface*& refs, inputDataNE& cfg, size_t startValue);
	void createchildHierarchy();
	void CreateVElem(size_t startValue);
	void createNewThroats(medialSurface*& refs);
	void findParentInElem(medialBall* vi, int e1s2, int& EbadElem, int& EbadBoss, int& nAllElem, int& nAllBoss);
	void createThroatBallConnectivity(medialSurface*& refs);
	void createCornerHirarchy();









	void writeStatoilFormat() const;



public:

	medialSurface*& rfs;
	// std::vector<mediaAxes*> mas;
	int firstPore;
	int firstPores;
	int lastPores;
	const inputDataNE& cg;
	voxel inletV, outletV;

	voxelImageT<int>  VElems; ///.  -vxlV-1,0,1

	std::vector<poreElementI*> poreIs;
	std::vector<throatElementI*> throatIs;

	std::vector<medialBall*> allSpace;
	vector<medialBall*> throadAdditBalls;



	int nPores;
	int nTrots;
	int nElems;



	const int maxNCors;

};





//void growPores_XX(voxelField<int>&  VElems, int min, int max, int porValue);
//void shrinkPores(voxelField<int>&  VElems, int min, int max, int porValue);

unsigned int growPores_X2(voxelField<int>&  VElems, int min, int max, int porValue);


void growPores(voxelField<int>&  VElems, int min, int max, int porValue);

void retreatPoresMedian(const inputDataNE & cg, voxelField<int>&  VElems, long min,  long max,
		const std::vector<poreElementI*>& poreIs, long rawValue);

void growPoresMedStrict(const inputDataNE & cg, voxelField<int>&  VElems, long min,  long max,
		const std::vector<poreElementI*>& poreIs, long rawValue);

void growPoresMedian(const inputDataNE & cg, voxelField<int>&  VElems, long min,  long max,
		const std::vector<poreElementI*>& poreIs, long rawValue);

void growPoresMedEqs(const inputDataNE & cg, voxelField<int>&  VElems, long min,  long max,
		const std::vector<poreElementI*>& poreIs, long rawValue);

void growPoresMedEqsLoose(const inputDataNE & cg, voxelField<int>&  VElems, long min,  long max,
		const std::vector<poreElementI*>& poreIs, long rawValue);






void medianElem(const inputDataNE & cg, voxelField<int>&  VElems, long min,  long max,
		const std::vector<poreElementI*>& poreIs);





#endif
