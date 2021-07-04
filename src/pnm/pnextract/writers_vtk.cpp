
/*---------------------------------------------------------------------------*\
written by: Ali Q Raeini  email: a.q.raeini@imperial.ac.uk
Imperial College, Total project on network extraction
https://www.imperial.ac.uk/earth-science/research/research-groups/pore-scale-modelling/
\*---------------------------------------------------------------------------*/




/// @cond INTERNAL
#define _nl_  ((i&31)==31 ? '\n' : ' ')


#include "typses.h"
#include "blockNet.h"
#include "writers.h"

#define _pi  3.14159265359
#define RadVTKFaCT  1.



using namespace std;


int findOrInsert(vector<dbl3>&  pts, dbl3& pt)  {
	for (auto p=pts.crbegin(); p!=pts.crbegin()+64 && p!=pts.crend(); ++p)
		if (*p==pt) return int(pts.rend()-p)-1;
	pts.push_back(pt);
	return pts.size()-1;
}



void insertHalfCorneroints2(vector<dbl3>&  points, vector<int>& cellPoints, 
		dbl3 c1, dbl3 c2, dbl3 nE1, dbl3 nE2, double apxDist_inrR, double apxDist_outR,
		double conAng1, double conAng2, double rOuter, double hafAng, double CARelax = 1.)  {
	dbl3 dd = c2-c1;
	dbl3 normal = dd/(mag(dd)+1e-32);

	vector<dbl3> hcPoints(8);

	dbl3 e1 = c1+nE1*rOuter;//. edgePoint
	dbl3 nE11 = rotateAroundVec(nE1,hafAng,normal);

	double gama = abs(hafAng)*CARelax;
	double appex_iInnerR = apxDist_inrR*(cos(gama)+(sin(gama)/cos(gama+conAng1))*(sin(gama+conAng1)-1));
	double appex_iOuterR = apxDist_outR*(cos(gama)+(sin(gama)/cos(gama+conAng2))*(sin(gama+conAng2)-1));

	if(hafAng>0)  {
		hcPoints[0] = e1-appex_iInnerR*nE1;
		hcPoints[1] = e1-appex_iOuterR*nE1;
		hcPoints[2] = e1-apxDist_outR*nE11;
		hcPoints[3] = e1-apxDist_inrR*nE11;  }
	else  {
		hcPoints[3] = e1-appex_iInnerR*nE1;
		hcPoints[2] = e1-appex_iOuterR*nE1;
		hcPoints[1] = e1-apxDist_outR*nE11;
		hcPoints[0] = e1-apxDist_inrR*nE11;  }



	e1 = c2+nE2*rOuter;///. edgePoint
	nE11 = rotateAroundVec(nE2,hafAng,normal);
	if(hafAng>0)  {
		hcPoints[4] = e1-appex_iInnerR*nE2;
		hcPoints[5] = e1-appex_iOuterR*nE2;
		hcPoints[6] = e1-apxDist_outR*nE11;
		hcPoints[7] = e1-apxDist_inrR*nE11;  }
	else  {
		hcPoints[7] = e1-appex_iInnerR*nE2;
		hcPoints[6] = e1-appex_iOuterR*nE2;
		hcPoints[5] = e1-apxDist_outR*nE11;
		hcPoints[4] = e1-apxDist_inrR*nE11;  }

	for (int i=0;i<8;++i)  ///. 8 points each elem
		cellPoints.push_back(findOrInsert(points, hcPoints[i]));
}





string vtkWriter_start(size_t nPoints, size_t nCells)  {
	return "<?xml version = \"1.\"?>\n"
	       "<VTKFile type = \"UnstructuredGrid\" version = \"0.1\" byte_order = \"LittleEndian\">\n"
	       " <UnstructuredGrid>"
	       "  <Piece NumberOfPoints = \""+_s(nPoints)+"\" NumberOfCells = \""+_s(nCells)+"\" >\n";
}
string  vtkWriter_finish()  {
	return "  </Piece>\n"  
	       " </UnstructuredGrid>\n"
	       "</VTKFile>\n";
 }




template<typename T>
void writeVtuArray(ofstream& outp, string name, const vector<T> & data, string typeStr="Float32")  {
	outp<<"\t\t<DataArray type = \""<<typeStr<<"\" Name = \""<<name<<"\" format = \"ascii\">\n";
	for_i(data)  {  outp << data[i] << _nl_;  }
	outp<<"\t\t</DataArray>"<<endl;
}

#define writeVtu_i(outp, name, data, Expr, typeStr)  \
	outp<<"\t\t<DataArray type = \"" typeStr "\" Name = \""<<name<<"\" format = \"ascii\">\n";  \
	for_i(data)  {  outp << Expr << _nl_;  }  \
	outp<<"\n\t\t</DataArray>\n";




void addThroatMesh(
	int ind,
	const vector<throatNE*>& throatIs,
	size_t trIndx,
	const vector<poreNE*>& poreIs,
	vector<dbl3>& points,
	vector<int>& subTypes,
	vector<size_t>& cellTrots,
	vector<int>& cellPors,
	vector<int>& cellPoints,
	double scaleFactor,
	bool visualizeCorners,
	unsigned int thetaResulution
 )  {
	const throatNE * elem = throatIs[trIndx];
	int porInd=(ind==1) ? elem->e1 : elem->e2;
	const poreNE * elemp =  poreIs[porInd];


	dbl3 c1(elem->mb22()->fi,elem->mb22()->fj,elem->mb22()->fk);
	dbl3 c2(elemp->mb->fi,elemp->mb->fj,elemp->mb->fk);
	double r = elem->radius()*scaleFactor;


	if (elem->e1<2 && ind==1)  {
		double throatIncLength=1.1*elem->radius()+2e-9;
		c1.y = c2.y;//y is wrong
		c1.z = c2.z;
		c2.y += 1e-9;	c2.z += 1e-9;
		if (c1.x<c2.x)    c1.x = c2.x - throatIncLength;
		else              c1.x = c2.x + throatIncLength;
	}
	if (elem->e2<2 && ind==2)  {
		double throatIncLength=1.1*elem->radius()+2e-9;
		c2.y = c1.y;//y is wrong
		c2.z = c1.z;
		c2.y += 1e-9;	c2.z += 1e-9;
		if (c2.x<c1.x)  c2.x = c1.x- throatIncLength;
		else            c2.x = c1.x+ throatIncLength;
	}

	dbl3 c1c2 = c2-c1;
	dbl3 ncc = c1c2/(mag(c1c2)+1e-33);
	dbl3 nCE(0.,0.,0.);	///. pick a corner point for the first subElem
	for(size_t i = 0;i<3;i++){ if (abs(ncc[i])<0.6){nCE[i] = 1.;break; } ; }

	nCE = ncc^nCE;
	nCE=rotateAroundVec(nCE, _pi/4., ncc);
	nCE = nCE/(mag(nCE)+1e-33);

	int thetaResulutionp2 = (thetaResulution+1)/2;

	for(int i = 0; i < 2*thetaResulutionp2; ++i)  {
		double hafAng = 0.5*_pi/thetaResulutionp2;
		nCE = rotateAroundVec(nCE, 2*hafAng, ncc);
		insertHalfCorneroints2( points,cellPoints,c1,c2,  -nCE, -nCE,  0,  r,  0., 0.6*(_pi-hafAng),   0,  hafAng);
		subTypes.push_back(5);
		cellTrots.push_back(trIndx);
		cellPors.push_back(porInd);

		insertHalfCorneroints2( points,cellPoints,c1,c2,  -nCE, -nCE,  0,  r,  0., 0.6*(_pi-hafAng),  0,  -hafAng);
		subTypes.push_back(5);
		cellTrots.push_back(trIndx);
		cellPors.push_back(porInd);
	}
}




void vtuWriteThroats(string basNam, const vector<poreNE*>& poreIs, const vector<throatNE*>& throatIs, double dx, dbl3 X0)  {

	vector<dbl3> points;
	vector<int> subTypes;
	vector<size_t> cellTrots;

	vector<int> cellPoints;
	points.reserve(throatIs.size()*150);
	cellPoints.reserve(throatIs.size()*300);
	subTypes.reserve(throatIs.size()*50);
	cellTrots.reserve(throatIs.size()*50);
	vector<int> cellPors;
	cellPors.reserve(throatIs.size()*50);



	for(size_t i = 0; i < throatIs.size(); ++i)  {
		addThroatMesh(1, throatIs,i,poreIs,points,subTypes,cellTrots,cellPors,cellPoints,0.2,false,4);
		addThroatMesh(2, throatIs,i,poreIs,points,subTypes,cellTrots,cellPors,cellPoints,0.2,false,4);
	}


	ofstream outp(basNam+".vtu");
	outp<<vtkWriter_start(points.size(),subTypes.size());

	outp<<"\t<Points>\n";
	outp<<"\t\t<DataArray type = \"Float32\" NumberOfComponents = \"3\" format = \"ascii\">\n";
	for_i(points)  {	outp << points[i]*dx+X0<< _nl_;  }
	outp<<"\n\t\t</DataArray>\n";
	outp<<"\t</Points>\n";
	
	outp<<"\t<Cells>\n";
	writeVtuArray(outp,"connectivity",cellPoints,"Int32");
	writeVtu_i(outp, "offsets",subTypes, 8*i+8,"Int32")
	writeVtu_i(outp, "types"  ,subTypes, 12,"UInt8")
	outp<<"\t</Cells>\n";

	outp<<"\t<CellData Scalars = \"index\">\n";
	writeVtuArray( outp, "trotIndex",  cellTrots, "Int32");
	writeVtuArray( outp, "index",  cellPors, "Int32");
	outp<<"\t</CellData>\n";

	outp<<vtkWriter_finish();
}








void AddCylinder(	dbl3& c1, dbl3& c2, double r,
	size_t& poreIndx,
 	vector<dbl3>& points,
	vector<int>& subTypes,
	vector<int>& cellPores,
	vector<float>& alpha,
	vector<int>& cellPoints,
	double& scaleFactor,
	unsigned int& thetaResulution
) {
	dbl3 c1c2 = c2-c1;
	dbl3 ncc = c1c2/(mag(c1c2)+1e-32);
	dbl3 nCE(0.,0.,0.);	///. pick a corner point for the first subElem
	for(size_t i = 0;i<3;i++){ if (ncc[i]<0.6){nCE[i] = 1.;break; } }

	nCE = ncc^nCE;
	nCE = nCE/(mag(nCE)+1e-32);

	int thetaResulutionp2 = (thetaResulution+1)/2;

	for(int i = 0; i < thetaResulutionp2; ++i)  {
		double hafAng = _pi/thetaResulutionp2;
		double hafAngleAzim = 0.5*_pi/thetaResulutionp2;
		nCE = rotateAroundVec(nCE, 2*hafAng, ncc);
		dbl3 lAzimuth1 = ncc^nCE; ///. normal to ncc and nE1
		dbl3 nCE2 = rotateAroundVec(nCE, thetaResulutionp2*hafAngleAzim, lAzimuth1);///. edge-centre ncc vector
		for(int j = 0; j < thetaResulutionp2; ++j)  {

			dbl3 nCE1 = nCE2;///. edge-centre ncc vector
			nCE2 = rotateAroundVec(nCE2, hafAngleAzim*2., lAzimuth1);///. edge-centre ncc vector

			insertHalfCorneroints2( points,cellPoints,c1,c2,  -nCE1,  -nCE2,  1e-18,  r, 0.,0.,   0,  hafAng,0.); ///. Warning: CA is not implemented for spheres
				subTypes.push_back(0);
			cellPores.push_back(poreIndx);

			insertHalfCorneroints2( points,cellPoints,c1,c2,  -nCE1,  -nCE2,  1e-18,  r,  0.,0.,   0,  -hafAng,0.);///. Warning: CA is not implemented for spheres
				subTypes.push_back(0);
			cellPores.push_back(poreIndx);
		}
	}
}


void addSpherePoreMesh(
	const vector<poreNE*>& poreIs,
	size_t poreIndx,
	const vector<throatNE*>& throatIs,
	vector<dbl3>& points,
	vector<int>& subTypes,
	vector<int>& cellPores,
	vector<float>& alpha,
	vector<int>& cellPoints,
	double scaleFactor,
	unsigned int thetaResulution
) {
	const poreNE * elem = poreIs[poreIndx];

	dbl3 c1(elem->mb->fi,elem->mb->fj,elem->mb->fk);
	dbl3 c2(elem->mb->fi,elem->mb->fj,elem->mb->fk);
	c2.y += 1e-12;
	double r = (elem->mb->R)*scaleFactor;
	AddCylinder(c1, c2, r, poreIndx, points, subTypes, cellPores, alpha, cellPoints, scaleFactor, thetaResulution );
}


void addCylinderThroatMesh(
	const vector<throatNE*>& throatIs,
	size_t trotIndx,
	const vector<poreNE*>& poreIs,
	vector<dbl3>& points,
	vector<int>& subTypes,
	vector<int>& cellPores,
	vector<float>& alpha,
	vector<int>& cellPoints,
	double scaleFactor,
	unsigned int thetaResulution
 ) {
	const throatNE * elem = throatIs[trotIndx];

	dbl3 c2(elem->mb22()->fi,elem->mb22()->fj,elem->mb22()->fk);
	dbl3 c1 = c2;//elem->e1>SideImax ? dbl3(elem->mb1->i,elem->mb1->j,elem->mb1->k) : c2;
	c2.y += 1e-11;
	double r = elem->mb22()->R*scaleFactor;
	AddCylinder(c1, c2, r, trotIndx, points, subTypes, cellPores, alpha, cellPoints, scaleFactor, thetaResulution );
}



void vtuWriteTHroatSpheres(string basNam, const vector<poreNE*>& poreIs, const vector<throatNE*>& throatIs, double dx, dbl3 X0)  {

	vector<dbl3> points;
	vector<int> subTypes;
	vector<int> cellPores;

	vector<int> cellPoints;
	points.reserve(poreIs.size()*150);
	cellPoints.reserve(poreIs.size()*300);
	subTypes.reserve(poreIs.size()*50);
	cellPores.reserve(poreIs.size()*50);
	vector<float> alpha;
	alpha.reserve(poreIs.size()*50);


	for(size_t i = 0; i <  throatIs.size(); ++i)
		addCylinderThroatMesh(throatIs,i,poreIs,points,subTypes,cellPores,alpha,cellPoints,RadVTKFaCT,8);


	ofstream outp(basNam+".vtu");
	outp<<vtkWriter_start(points.size(),subTypes.size());

	outp<<"\t<Points>\n";
	outp<<"\t\t<DataArray type = \"Float32\" NumberOfComponents = \"3\" format = \"ascii\">\n";
	for_i(points)  {	outp << points[i]*dx+X0 << _nl_;  }
	outp<<"\n\t\t</DataArray>\n";
	outp<<"\t</Points>\n";

	outp<<"\t<Cells>\n";	
	writeVtuArray(outp,"connectivity",cellPoints,"Int32");
	writeVtu_i(outp, "offsets",subTypes, 8*i+8,"Int32");
	writeVtu_i(outp, "types"  ,subTypes, 12,"UInt8");
	outp<<"\t</Cells>\n";

	outp<<"\t<CellData Scalars = \"subType\">\n";
	writeVtuArray( outp, "subType",  subTypes, "Int32");
	writeVtuArray( outp, "index",  cellPores, "Int32");
	outp<<"\t</CellData>\n";

	outp<<vtkWriter_finish();
}




void vtuWritePores(string basNam, const vector<poreNE*>& poreIs, const vector<throatNE*>& throatIs, double dx, dbl3 X0)  {

	vector<dbl3> points;
	vector<int> subTypes;
	vector<int> cellPores;

	vector<int> cellPoints;
	//vector<FacePoints> facePoints;
	//vector<CellFaces> cellFaces;
	//vector<FaceCells> faceCells;
	points.reserve(poreIs.size()*200);
	cellPoints.reserve(poreIs.size()*400);
	subTypes.reserve(poreIs.size()*60);
	cellPores.reserve(poreIs.size()*60);
	vector<float> alpha;
	alpha.reserve(poreIs.size()*60);


	for(size_t i = 2; i <  poreIs.size(); ++i)
	  addSpherePoreMesh(poreIs,i,throatIs,points,subTypes,cellPores,alpha,cellPoints,RadVTKFaCT,8);


	ofstream outp(basNam+".vtu");
	outp<<vtkWriter_start(points.size(),subTypes.size());

	outp<<"\t<Points>\n";
	outp<<"\t\t<DataArray type = \"Float32\" NumberOfComponents = \"3\" format = \"ascii\">\n";
	for_i(points)  {	outp << points[i]*dx+X0<< " " << _nl_;  }
	outp<<"\n\t\t</DataArray>\n";
	outp<<"\t</Points>\n";

	outp<<"\t<Cells>\n";
	writeVtuArray(outp,"connectivity",cellPoints,"Int32");
	writeVtu_i(outp, "offsets",subTypes, 8*i+8,"Int32");
	writeVtu_i(outp, "types"  ,subTypes, 12,"UInt8");
	outp<<"\t</Cells>\n";

	outp<<"\t<CellData Scalars = \"alpha\">\n";
	writeVtuArray( outp, "subType",  subTypes, "Int32");
	writeVtuArray( outp, "index",  cellPores, "Int32");
	outp<<"\t</CellData>\n";

	outp<<vtkWriter_finish();
}










void addMbMbMesh(
	const medialBall& vi,    const medialBall* vj,
	size_t poreIndx,         vector<dbl3>& points,
	vector<int>& subTypes,   vector<int>& pointTyps,
	vector<int>& cellPores,  vector<float>& radius,
	vector<int>& cellPoints
 )  {

	dbl3 c1(vi.fi,vi.fj,vi.fk);
	dbl3 c2(vj->fi,vj->fj,vj->fk);
	if (c1==c2) {c1.z-=0.2;c2.z+=0.2; }

	cellPoints.push_back(findOrInsert(points, c1));
	if(radius.size()<points.size()) {radius.push_back(vi.R);  pointTyps.push_back(vi.type); }
	cellPoints.push_back(findOrInsert(points, c2));
	if(radius.size()<points.size()) {radius.push_back(vj->R); pointTyps.push_back(vj->type); }

	cellPores.push_back(poreIndx); //. 2 points each elem
	subTypes.push_back(vi.corId); //. 2 points each elem
}


void vtuWriteMbMbs(string basNam, const vector<medialBall>& ballSpace, const  vector<poreNE*> poreIs, const voxelField<int>&  VElems, double dx, dbl3 X0 )  {
	cout<<"\nwrite_hierarchy >> "+basNam+".vtu  "<<len(ballSpace)<<endl;

	vector<dbl3> points;
	vector<int> subTypes;
	vector<int> pointTyps;
	vector<int> cellPores;
	vector<float> radius;

	vector<int> cellPoints;
	//vector<FacePoints> facePoints;
	//vector<CellFaces> cellFaces;
	//vector<FaceCells> faceCells;
	points.reserve(ballSpace.size()*2);
	cellPoints.reserve(ballSpace.size()*2);
	subTypes.reserve(ballSpace.size()*2);
	pointTyps.reserve(ballSpace.size()*2);
	cellPores.reserve(ballSpace.size()*2);
	radius.reserve(ballSpace.size()*2);


	for (const auto& vi:ballSpace)  if (vi.boss)  {
		medialBall& mvi=*vi.mastrSphere();
		int id=VElems(mvi.fi+1, mvi.fj+1, mvi.fk+1);
		addMbMbMesh(vi, vi.boss, id,points,subTypes,pointTyps,cellPores,radius,cellPoints);
	}


	ofstream outp(basNam+".vtu");
	outp<<vtkWriter_start(points.size(),subTypes.size());

	outp<<"\t<Points>\n";
	outp<<"\t\t<DataArray type = \"Float32\" NumberOfComponents = \"3\" format = \"ascii\">\n";
	for_i(points)  { outp << points[i]*dx+X0 << _nl_;  }
	outp<<"\n\t\t</DataArray>\n";
	outp<<"\t</Points>\n";

	outp<<"\t<Cells>\n";
	writeVtuArray(outp,"connectivity",cellPoints,"Int32");
	writeVtu_i(outp, "offsets",subTypes, 2*i+2,"Int32");
	writeVtu_i(outp, "types"  ,subTypes, 3,"UInt8");
	outp<<"\t</Cells>\n";


	outp<<"\t<CellData Scalars = \"alpha\">\n";
	writeVtuArray( outp, "index",  cellPores, "Int32");
	writeVtuArray( outp, "subType",  subTypes, "Int32");
	outp<<"\t</CellData>\n";

	outp<<"\t<PointData>\n";
	writeVtuArray( outp, "radius",  radius);
	writeVtuArray( outp, "cellType",  pointTyps, "Int32");
	outp<<"\t</PointData>\n";

	outp<<vtkWriter_finish();
}



void vtuWriteThroatMbMbs(string basNam,      /// ********* throat *********
 const vector<throatNE*>& throatIs,
 const  vector<poreNE*> poreIs,
 const voxelField<int>&  VElems, double dx, dbl3 X0
 )  {
	cout<<"\nwrite_throatHierarchy >> "+basNam+".vtu"<<endl;

	vector<dbl3> points;
	vector<int> subTypes;
	vector<int> pointTyps;
	vector<int> cellPores;
	vector<float> radius;

	vector<int> cellPoints;
	//vector<FacePoints> facePoints;	//vector<CellFaces> cellFaces;	//vector<FaceCells> faceCells;
	points.reserve(throatIs.size()*20);
	cellPoints.reserve(throatIs.size()*20);
	subTypes.reserve(throatIs.size()*20);
	pointTyps.reserve(throatIs.size()*20);
	cellPores.reserve(throatIs.size()*20);
	radius.reserve(throatIs.size()*20);

	int E2msfc=0, E_msfd=0, E8msgl=0, Dkfhj2=0, Dkfhj1=0;

	for_i(throatIs)  {  const auto& tr = *throatIs[i];
		const medialBall *mb2 = tr.mb22(),   *mb1 = tr.mb11();

		if (mb1)
			addMbMbMesh(*mb1, mb2 ,i+poreIs.size()*10, points,subTypes,pointTyps,cellPores,radius,cellPoints);
		else 										if (tr.e1>2) cout<<" EioMbk ";

		medialBall* mvi=mb2->mastrSphere();
													if(tr.e2 != VElems(mvi->fi+1, mvi->fj+1, mvi->fk+1)) ++Dkfhj1;
													if(mb2->boss==mb2)      ++E_msfd;

		while (mb2->boss != mb2)  {
			addMbMbMesh(*mb2, mb2->boss, tr.e2, points,subTypes,pointTyps,cellPores,radius,cellPoints);
			mb2 = mb2->boss;
		}

													if(mb1 && mb1->boss==mb1) ++E2msfc;
		while ( mb1 && mb1->boss != mb1)  {
			mvi=mb1->mastrSphere();
													if(tr.e1 != VElems(mvi->fi+1, mvi->fj+1, mvi->fk+1)) ++Dkfhj2;
			addMbMbMesh(*mb1, mb1->boss, tr.e1, points,subTypes,pointTyps,cellPores,radius,cellPoints);
			mb1 = mb1->boss;
		}

													if (mb1==mb2) ++E8msgl;
	}

	if (E2msfc||E_msfd||Dkfhj2||E8msgl||Dkfhj1)
		cout<<"\n *********** bads: mb1->boss==mb1 -> "<<E2msfc<<",   mb2->boss==mb2 -> "<<E_msfd<<",   e1!=VElem -> "<<Dkfhj2<<",   mb1==mb2 -> "<<E8msgl<<",   e2!=VElem -> "<<Dkfhj1<<endl;

	ofstream outp(basNam+".vtu");

	outp<<vtkWriter_start(points.size(),subTypes.size());


	outp<<"\t<Points>\n";
	outp<<"\t\t<DataArray type = \"Float32\" NumberOfComponents = \"3\" format = \"ascii\">\n";
	for_i(points) outp << points[i]*dx+X0<< _nl_;
	outp<<"\n\t\t</DataArray>\n";
	outp<<"\t</Points>\n";

	outp<<"\t<Cells>\n";
	writeVtuArray(outp,"connectivity",cellPoints,"Int32");
	writeVtu_i(outp, "offsets",subTypes, 2*i+2,"Int32");
	writeVtu_i(outp, "types"  ,subTypes, 3,"UInt8");
	outp<<"\t</Cells>\n";


	outp<<"\t<PointData>\n";
	writeVtuArray(outp, "radius", radius);
	writeVtuArray(outp, "cellType", pointTyps);
	outp<<"\t</PointData>\n";

	outp<<"\t<CellData Scalars = \"subType\">\n";
	writeVtuArray( outp, "index",  cellPores, "Int32");
	writeVtuArray( outp, "subType",  subTypes, "Int32");
	outp<<"\t</CellData>\n";

	outp<<vtkWriter_finish();
}





/// @endcond

