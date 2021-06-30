#ifndef BLOCKNET_H
#define BLOCKNET_H




#include <map>
#include <fstream>


#include "inputData.h"
#include "medialSurf.h"

#define _SideRadius 1
#define RLEVEL3 0.7

#define NINE 7
#define _r10(_id)  ((_id)>>3)
#define _pc10(_id) ((_id)-(((_id) >> 3)<<3))
#define _x10(_id)  ((_id)<<3)

#define SubElem0 nElems



#define _frmtd_(_d_m,_p_m,_w_m) setw(_w_m)<<(1./_p_m)*round((_d_m)*_p_m)
#define _frmtd3_(_d_m,_p_m,_w_m)  _frmtd_((_d_m)._0(),_p_m,_w_m) <<" "<< _frmtd_((_d_m)._1(),_p_m,_w_m) <<" "<< _frmtd_((_d_m)._2(),_p_m,_w_m)




int nextract(inputDataNE& cfg, bool verbos=true);



class blockNetwork
{
 public:

	blockNetwork(medialSurface*& rs, const inputDataNE& cfg)
	: srf(rs), cg(cfg), maxNCors(5) { nBP6=cfg.nBP6; SideImax=nBP6-1; };


	void createMedialSurface(medialSurface*& srf, inputDataNE& cfg, size_t startValue);
	void collectAllballs(medialSurface*& srf, inputDataNE& cfg, size_t startValue);
	void createchildHierarchy();
	void CreateVElem(size_t startValue);
	void createNewThroats(medialSurface*& srf);
	void findParentInElem(medialBall* vi, int e1s2, int& EbadElem, int& EbadBoss, int& nAllElem, int& nAllBoss);
	void createThroatBallConnectivity(medialSurface*& srf);
	void createCornerHirarchy();






	void writePNM() const;




	static int SideImax;
	static int nBP6; // 1+SideImax

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
