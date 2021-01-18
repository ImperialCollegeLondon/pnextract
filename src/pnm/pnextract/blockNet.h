#ifndef POROUSBLOCK_H
#define POROUSBLOCK_H




#include <map>
#include <fstream>


#include "inputData.h"
#include "medialSurf.h"



int nextract(inputDataNE& cfg, bool verbos=true);



class blockNetwork
{
 public:

	blockNetwork(medialSurface*& rs, const inputDataNE& c)
	: srf(rs), cg(c), maxNCors(5) {};


	void createMedialSurface(medialSurface*& srf, inputDataNE& cfg, size_t startValue);
	void collectAllballs(medialSurface*& srf, inputDataNE& cfg, size_t startValue);
	void createchildHierarchy();
	void CreateVElem(size_t startValue);
	void createNewThroats(medialSurface*& srf);
	void findParentInElem(medialBall* vi, int e1s2, int& EbadElem, int& EbadBoss, int& nAllElem, int& nAllBoss);
	void createThroatBallConnectivity(medialSurface*& srf);
	void createCornerHirarchy();






	void writePNM() const;





public:

	medialSurface*&   srf;
	// std::vector<mediaAxes*> mas;
	int firstPore;
	int firstPores;
	int lastPores;
	const inputDataNE& cg;
	std::array<voxel,8> sides;

	voxelImageT<int>  VElems;

	std::vector<poreNE*> poreIs;
	std::vector<throatNE*> throatIs;

	std::vector<medialBall*> allSpace;
	vector<medialBall*> throadAdditBalls;



	int nNodes;
	int nTrots;
	int nElems;



	const int maxNCors;

};





//void growPores_XX(voxelField<int>&  VElems, int min, int max, int porValue);
//void shrinkPores(voxelField<int>&  VElems, int min, int max, int porValue);

unsigned int growPores_X2(voxelField<int>&  VElems, int min, int max, int porValue);


void growPores(voxelField<int>&  VElems, int min, int max, int porValue);

void retreatPoresMedian(const inputDataNE & cg, voxelField<int>&  VElems, long min,  long max,
		const std::vector<poreNE*>& poreIs, long rawValue);

void growPoresMedStrict(const inputDataNE & cg, voxelField<int>&  VElems, long min,  long max,
		const std::vector<poreNE*>& poreIs, long rawValue);

void growPoresMedian(const inputDataNE & cg, voxelField<int>&  VElems, long min,  long max,
		const std::vector<poreNE*>& poreIs, long rawValue);

void growPoresMedEqs(const inputDataNE & cg, voxelField<int>&  VElems, long min,  long max,
		const std::vector<poreNE*>& poreIs, long rawValue);

void growPoresMedEqsLoose(const inputDataNE & cg, voxelField<int>&  VElems, long min,  long max,
		const std::vector<poreNE*>& poreIs, long rawValue);






void medianElem(const inputDataNE & cg, voxelField<int>&  VElems, long min,  long max,
		const std::vector<poreNE*>& poreIs);





#endif
