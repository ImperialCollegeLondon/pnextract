
/*---------------------------------------------------------------------------*\
written by: Ali Q Raeini  email: a.q.raeini@imperial.ac.uk
Imperial College, Total project on network extraction
https://www.imperial.ac.uk/earth-science/research/research-groups/pore-scale-modelling/
\*---------------------------------------------------------------------------*/








#include <vector>
#include <cmath>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <set>
#include <cassert>

using namespace std;


#include "typses.h"
#include "blockNet.h"
#include "writers.h"

#define _pi  3.14159265359
#define RadVTKFaCT  1.





dbl3 rotate( dbl3 n, dbl3 x, dbl3 y, double gamma)  {
	return dbl3
	  (	( x[0]*(n[1]*n[1]+n[2]*n[2]) - n[0]*( x[1]*n[1]+x[2]*n[2]-n[0]*y[0]- n[1]*y[1]-n[2]*y[2] ) )*(1-cos(gamma)) + y[0]*cos(gamma) + (-x[2]*n[1]+x[1]*n[2]-n[2]*y[1]+n[1]*y[2] )*sin(gamma),
		( x[1]*(n[0]*n[0]+n[2]*n[2]) - n[1]*( x[0]*n[0]+x[2]*n[2]-n[0]*y[0]- n[1]*y[1]-n[2]*y[2] ) )*(1-cos(gamma)) + y[1]*cos(gamma) + ( x[2]*n[0]-x[0]*n[2]+n[2]*y[0]-n[0]*y[2] )*sin(gamma),
		( x[2]*(n[0]*n[0]+n[1]*n[1]) - n[2]*( x[0]*n[0]+x[1]*n[1]-n[0]*y[0]- n[1]*y[1]-n[2]*y[2] ) )*(1-cos(gamma)) + y[2]*cos(gamma) + (-x[1]*n[0]+x[0]*n[1]-n[1]*y[0]+n[0]*y[1] )*sin(gamma)
	  );
}
dbl3 rotateAroundVec( dbl3 n, dbl3 y, double gamma)  {
	return dbl3
	  (	(  - n[0]*( -n[0]*y[0]- n[1]*y[1]-n[2]*y[2] ) )*(1-cos(gamma)) + y[0]*cos(gamma) + (n[1]*y[2]-n[2]*y[1])*sin(gamma),
		(  - n[1]*( -n[0]*y[0]- n[1]*y[1]-n[2]*y[2] ) )*(1-cos(gamma)) + y[1]*cos(gamma) + (n[2]*y[0]-n[0]*y[2])*sin(gamma),
		(  - n[2]*( -n[0]*y[0]- n[1]*y[1]-n[2]*y[2] ) )*(1-cos(gamma)) + y[2]*cos(gamma) + (n[0]*y[1]-n[1]*y[0])*sin(gamma)
	  );
	//-n*(n&y)*(1-cos(gamma)) + y*cos(gamma) + n^y*sin(gamma);
}




int findOrInsertPoint(vector<dbl3>&  points, dbl3& point)  {

	for (vector<dbl3>::const_reverse_iterator  rip = points.rbegin(); rip != points.rbegin()+64 && rip != points.rend(); ++rip)
		if (*rip == point) return int(points.rend()-rip)-1;

	points.push_back(point);
	return points.size()-1;
}



void insertHalfCorneroints2(vector<dbl3>&  points, vector<int>& cellPoints, dbl3 c1, dbl3 c2, dbl3 nE1, dbl3 nE2, double appexDist_inrR, double appexDist_outR,
							double conAng1, double conAng2, double rOuter, double hafAng, double CARelax = 1.)  {
	//const double convertToRad = _pi/180.;
	dbl3 dd = c2-c1;
	dbl3 normal = dd/(mag(dd)+1e-32);

	vector<dbl3> hcPoints(8);

	dbl3 e1 = c1+nE1*rOuter;//. edgePoint
	dbl3 nE11 = rotateAroundVec(normal,nE1,hafAng);

	///. radii of curvature

	double gama = abs(hafAng)*CARelax;
	double appex_iInnerR = appexDist_inrR*(cos(gama)+(sin(gama)/cos(gama+conAng1))*(sin(gama+conAng1)-1));
	double appex_iOuterR = appexDist_outR*(cos(gama)+(sin(gama)/cos(gama+conAng2))*(sin(gama+conAng2)-1));

	if(hafAng>0)  {
		hcPoints[0] = e1-appex_iInnerR*nE1;
		hcPoints[1] = e1-appex_iOuterR*nE1;
		hcPoints[2] = e1-appexDist_outR*nE11;
		hcPoints[3] = e1-appexDist_inrR*nE11;
	}
	else
	{
		hcPoints[3] = e1-appex_iInnerR*nE1;
		hcPoints[2] = e1-appex_iOuterR*nE1;
		hcPoints[1] = e1-appexDist_outR*nE11;
		hcPoints[0] = e1-appexDist_inrR*nE11;
	}



	 e1 = c2+nE2*rOuter;///. edgePoint
	 nE11 = rotateAroundVec(normal,nE2,hafAng);
	if(hafAng>0)  {
		hcPoints[4] = e1-appex_iInnerR*nE2;
		hcPoints[5] = e1-appex_iOuterR*nE2;
		hcPoints[6] = e1-appexDist_outR*nE11;
		hcPoints[7] = e1-appexDist_inrR*nE11;
	}
	else
	{
		hcPoints[7] = e1-appex_iInnerR*nE2;
		hcPoints[6] = e1-appex_iOuterR*nE2;
		hcPoints[5] = e1-appexDist_outR*nE11;
		hcPoints[4] = e1-appexDist_inrR*nE11;
	}



	for (int i=0;i<8;++i)  {///. 8 points each elem
		int id=(findOrInsertPoint(points, hcPoints[i]));
		cellPoints.push_back(id);
	}

}






//////////////////////////////////////////////////////////
 string vtkWriter_start(size_t nPoints, size_t nCells)  {
		stringstream  str;
		str<<"<?xml version = \"1.\"?>\n"
		   <<"<VTKFile type = \"UnstructuredGrid\" version = \"0.1\" byte_order = \"LittleEndian\">\n"
		   <<"  <UnstructuredGrid>"
		   <<"    <Piece NumberOfPoints = \""<<nPoints<<"\" NumberOfCells = \""<<nCells<<"\" >\n";
		return str.str();
	}
	 string  vtkWriter_finish()  {
		stringstream  str;
		str<<"    </Piece>\n"
		   <<"  </UnstructuredGrid>\n"
		   <<"</VTKFile>\n";
		return str.str();
	}




template<typename Type>
void writeCellData(ofstream& outp, string name, const vector<Type> & data, string typeStr="Float32")  {
	outp<<"        <DataArray type = \""<<typeStr<<"\" Name = \""<<name<<"\" format = \"ascii\">\n";
    for(size_t i = 0; i < data.size(); ++i)  {
        outp << data[i] << " ";
        if ((i+1)%20 == 0)     outp << "\n";
    }
	outp<<"        </DataArray>"<<endl;
}






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
		c1[1] = c2[1];//y is wrong
		c1[2] = c2[2];


		c2[1] += 1e-9;	c2[2] += 1e-9;

		if (c1[0]<c2[0])  {
			c1[0] = c2[0] - throatIncLength;
		}
		else
		{
			c1[0] = c2[0] + throatIncLength;
		}
	}
	if (elem->e2<2 && ind==2)  {
		double throatIncLength=1.1*elem->radius()+2e-9;

		c2[1] = c1[1];//y is wrong
		c2[2] = c1[2];

		c2[1] += 1e-9;	c2[2] += 1e-9;

		if (c2[0]<c1[0])  {
			c2[0] = c1[0]-  throatIncLength;
		}
		else
		{
			c2[0] = c1[0]+ throatIncLength;
		}
	}

	dbl3 c1c2 = c2-c1;
	dbl3 ncc = c1c2/(mag(c1c2)+1e-33);
	dbl3 nCE(0.,0.,0.);	///. pick a corner point for the first subElem
	for(size_t i = 0;i<3;i++){ if (abs(ncc[i])<0.6){nCE[i] = 1.;break;} ;}

	nCE = ncc^nCE;
	nCE=rotateAroundVec( ncc, nCE, _pi/4.);
	nCE = nCE/(mag(nCE)+1e-33);

	{
		int thetaResulutionp2 = (thetaResulution+1)/2;

		for(int i = 0; i < 2*thetaResulutionp2; ++i)  {
			double hafAng = 0.5*_pi/thetaResulutionp2;
			nCE = rotateAroundVec(ncc,nCE,2*hafAng);
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



	outp<<"      <Points>\n";
	outp<<"        <DataArray type = \"Float32\" NumberOfComponents = \"3\" format = \"ascii\">\n";
    for(size_t i = 0; i < points.size(); ++i)  {
        outp << (points[i][0]*dx+X0[0])<< " " << (points[i][1]*dx+X0[1])<< " " << (points[i][2]*dx+X0[2])<< " ";
        if ((i+1)%20 == 0)     outp << "\n";
    }
	outp<<"\n        </DataArray>\n";
	outp<<"      </Points>\n";

	outp<<"      <Cells>\n";///// \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/
	outp<<"        <DataArray type = \"Int32\" Name = \"connectivity\" format = \"ascii\">\n";
    for(size_t i = 0; i < cellPoints.size(); ++i)  {
        outp << cellPoints[i] << " ";
        if ((i+1)%20 == 0)     outp << "\n";
    }

	outp<<"\n        </DataArray>\n";

	outp<<"	<DataArray type = \"Int32\" Name = \"offsets\" format = \"ascii\">\n";
    for(size_t i = 0; i < subTypes.size(); ++i)  {
        outp << 8*i+8 << " ";
        if ((i+1)%20 == 0)     outp << "\n";
    }
	outp<<"\n        </DataArray>\n";

	outp<<"	<DataArray type = \"UInt8\" Name = \"types\" format = \"ascii\">\n";
    for(size_t i = 0; i < subTypes.size(); ++i)  {
        outp << 12 << " ";
        if ((i+1)%20 == 0)     outp << "\n";
    }
	outp<<"\n        </DataArray>\n";
	outp<<"      </Cells>\n";// //  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//




	outp<<"      <CellData Scalars = \"index\">\n"; //////////////////////////////////////


	writeCellData( outp, "trotIndex",  cellTrots, "Int32");
	writeCellData( outp, "index",  cellPors, "Int32");


	outp<<"      </CellData>\n"; /////////////////////////////////////////////////////////




    outp<<vtkWriter_finish();
    outp.close();
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
)  {
	dbl3 c1c2 = c2-c1;
	dbl3 ncc = c1c2/(mag(c1c2)+1e-32);
	dbl3 nCE(0.,0.,0.);	///. pick a corner point for the first subElem
	for(size_t i = 0;i<3;i++){ if (ncc[i]<0.6){nCE[i] = 1.;break;} }

	nCE = ncc^nCE;
	nCE = nCE/(mag(nCE)+1e-32);

    //const Polygon* polyShape = dynamic_cast< const Polygon* >(elem->shape());
    //if(0 && polyShape) ///. for now ignore corners    { }	else
	{
		int thetaResulutionp2 = (thetaResulution+1)/2;

		for(int i = 0; i < thetaResulutionp2; ++i)  {
			double hafAng = _pi/thetaResulutionp2;
			double hafAngleAzim = 0.5*_pi/thetaResulutionp2;
			nCE = rotateAroundVec(ncc,nCE,2*hafAng);
			dbl3 lAzimuth1 = ncc^nCE; ///. normal to ncc and nE1
			dbl3 nCE2 = rotateAroundVec(lAzimuth1,nCE,thetaResulutionp2*hafAngleAzim);///. edge-centre ncc vector
			for(int j = 0; j < thetaResulutionp2; ++j)  {

				dbl3 nCE1 = nCE2;///. edge-centre ncc vector
				nCE2 = rotateAroundVec(lAzimuth1,nCE2,hafAngleAzim*2.);///. edge-centre ncc vector

				insertHalfCorneroints2( points,cellPoints,c1,c2,  -nCE1,  -nCE2,  1e-18,  r, 0.,0.,   0,  hafAng,0.); ///. Warning: CA is not implemented for spheres
					subTypes.push_back(0);
				cellPores.push_back(poreIndx);
				//alpha.push_back(elem->shape()->containCOil());

				insertHalfCorneroints2( points,cellPoints,c1,c2,  -nCE1,  -nCE2,  1e-18,  r,  0.,0.,   0,  -hafAng,0.);///. Warning: CA is not implemented for spheres
					subTypes.push_back(0);
				cellPores.push_back(poreIndx);
				//alpha.push_back(elem->shape()->containCOil());
			}
		}
	}

}


void addSpherePoreMesh
(
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
)  {
	const poreNE * elem = poreIs[poreIndx];

	dbl3 c1(elem->mb->fi,elem->mb->fj,elem->mb->fk);
	dbl3 c2(elem->mb->fi,elem->mb->fj,elem->mb->fk);
	c2[1] += 1e-12;
	double r = (elem->mb->R)*scaleFactor;
	AddCylinder(c1, c2, r, poreIndx, points, subTypes, cellPores, alpha, cellPoints, scaleFactor, thetaResulution );
}


void addCylinderThroatMesh
(
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
)  {
	const throatNE * elem = throatIs[trotIndx];

	dbl3 c2(elem->mb22()->fi,elem->mb22()->fj,elem->mb22()->fk);
	dbl3 c1 = c2;//elem->e1>SideImax ? dbl3(elem->mb1->i,elem->mb1->j,elem->mb1->k) : c2;
	c2[1] += 1e-11;
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


    for(size_t i = 0; i <  throatIs.size(); ++i)  {
        addCylinderThroatMesh(throatIs,i,poreIs,points,subTypes,cellPores,alpha,cellPoints,RadVTKFaCT,8);
    }


    ofstream outp(basNam+".vtu");

	outp<<vtkWriter_start(points.size(),subTypes.size());




	outp<<"      <Points>\n";
	outp<<"        <DataArray type = \"Float32\" NumberOfComponents = \"3\" format = \"ascii\">\n";
    for(size_t i = 0; i < points.size(); ++i)  {
        outp << (points[i][0]*dx+X0[0])<< " " << (points[i][1]*dx+X0[1])<< " " << (points[i][2]*dx+X0[2])<< " ";
        if ((i+1)%20 == 0)     outp << "\n";
    }
	outp<<"\n        </DataArray>\n";
	outp<<"      </Points>\n";

	outp<<"      <Cells>\n";///// \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/
	outp<<"        <DataArray type = \"Int32\" Name = \"connectivity\" format = \"ascii\">\n";
    for(size_t i = 0; i < cellPoints.size(); ++i)  {
        outp << cellPoints[i] << " ";
        if ((i+1)%20 == 0)     outp << "\n";
    }
	outp<<"\n        </DataArray>\n";

	outp<<"	<DataArray type = \"Int32\" Name = \"offsets\" format = \"ascii\">\n";
    for(size_t i = 0; i < subTypes.size(); ++i)  {
        outp << 8*i+8 << " ";
        if ((i+1)%20 == 0)     outp << "\n";
    }
	outp<<"\n        </DataArray>\n";

	outp<<"	<DataArray type = \"UInt8\" Name = \"types\" format = \"ascii\">\n";
    for(size_t i = 0; i < subTypes.size(); ++i)  {
        outp << 12 << " ";
        if ((i+1)%20 == 0)     outp << "\n";
    }
	outp<<"\n        </DataArray>\n";
	outp<<"      </Cells>\n";// //  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//



	outp<<"      <CellData Scalars = \"alpha\">\n"; //////////////////////////////////////

	vector<float> saturation(cellPores.size(),0.);
	vector<float> pis_c(cellPores.size(),0.);
	vector<float> p(cellPores.size(),0.);
	vector<float> p_o(cellPores.size(),0.);
	vector<float> p_w(cellPores.size(),0.);
	vector<float> elem_type(cellPores.size(),0.);
	writeCellData( outp, "subType",  subTypes);
	writeCellData( outp, "type",  elem_type);
	writeCellData( outp, "index",  cellPores, "Int32");

	outp<<"      </CellData>\n"; /////////////////////////////////////////////////////////




    outp<<vtkWriter_finish();
    outp.close();
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


	for(size_t i = 2; i <  poreIs.size(); ++i)  {
	  addSpherePoreMesh(poreIs,i,throatIs,points,subTypes,cellPores,alpha,cellPoints,RadVTKFaCT,8);
	}


	ofstream outp(basNam+".vtu");

	outp<<vtkWriter_start(points.size(),subTypes.size());




	outp<<"      <Points>\n";
	outp<<"        <DataArray type = \"Float32\" NumberOfComponents = \"3\" format = \"ascii\">\n";
    for(size_t i = 0; i < points.size(); ++i)  {
        outp << (points[i][0]*dx+X0[0])<< " " << (points[i][1]*dx+X0[1])<< " " << (points[i][2]*dx+X0[2])<< " ";
        if ((i+1)%20 == 0)     outp << "\n";
    }
	outp<<"\n        </DataArray>\n";
	outp<<"      </Points>\n";

	outp<<"      <Cells>\n";///// \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/
	outp<<"        <DataArray type = \"Int32\" Name = \"connectivity\" format = \"ascii\">\n";
    for(size_t i = 0; i < cellPoints.size(); ++i)  {
        outp << cellPoints[i] << " ";
        if ((i+1)%20 == 0)     outp << "\n";
    }
	outp<<"\n        </DataArray>\n";

	outp<<"	<DataArray type = \"Int32\" Name = \"offsets\" format = \"ascii\">\n";
    for(size_t i = 0; i < subTypes.size(); ++i)  {
        outp << 8*i+8 << " ";
        if ((i+1)%20 == 0)     outp << "\n";
    }
	outp<<"\n        </DataArray>\n";

	outp<<"	<DataArray type = \"UInt8\" Name = \"types\" format = \"ascii\">\n";
    for(size_t i = 0; i < subTypes.size(); ++i)  {
        outp << 12 << " ";
        if ((i+1)%20 == 0)     outp << "\n";
    }
	outp<<"\n        </DataArray>\n";
	outp<<"      </Cells>\n";// //  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//



	outp<<"      <CellData Scalars = \"alpha\">\n"; //////////////////////////////////////

	vector<float> saturation(cellPores.size(),0.);
	vector<float> pis_c(cellPores.size(),0.);
	vector<float> p(cellPores.size(),0.);
	vector<float> p_o(cellPores.size(),0.);
	vector<float> p_w(cellPores.size(),0.);
	vector<float> elem_type(cellPores.size(),0.);
	writeCellData( outp, "subType",  subTypes);
	writeCellData( outp, "type",  elem_type);
	writeCellData( outp, "index",  cellPores, "Int32");

	outp<<"      </CellData>\n"; /////////////////////////////////////////////////////////




    outp<<vtkWriter_finish();
    outp.close();
}















// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // /


void addMbMbMesh(
	 const medialBall* vi,
	 const medialBall* vj,
	 size_t poreIndx,
	 vector<dbl3>& points,
	 vector<int>& subTypes,
	 vector<int>& cellType,
	 vector<int>& cellPores,
	 vector<float>& radius,
	 vector<int>& cellPoints
	 )  {

	dbl3 c1(vi->fi,vi->fj,vi->fk);
	dbl3 c2(vj->fi,vj->fj,vj->fk);
	if (c1==c2) {c1.z-=0.2;c2.z+=0.2;}

	cellPoints.push_back(findOrInsertPoint(points, c1));
	if(radius.size()<points.size()) {radius.push_back(vi->R);  cellType.push_back(vi->type); }
	cellPoints.push_back(findOrInsertPoint(points, c2));
	if(radius.size()<points.size()) {radius.push_back(vj->R); cellType.push_back(vj->type);}

	cellPores.push_back(poreIndx); ///. 2 points each elem
	subTypes.push_back(vi->corId); ///. 2 points each elem
//	cellType.push_back(vi->type); ///. 2 points each elem
	//if(distSqr(vi,vi->boss)< sqrt(vi->boss->limit)*1.2+2-sqrt(vi->limit))  *subTypes.rbegin()=-1;

}

void addMbMbMbMesh(
	 const medialBall* vi,
	 const medialBall* vj,
	 const medialBall* iKid,
	 const medialBall* jKid,
	 size_t poreIndx,
	 vector<dbl3>& points,
	 vector<int>& subTypes,
	 vector<int>& cellType,
	 vector<int>& cellPores,
	 vector<float>& radius,
	 vector<int>& cellPoints
	 )  {

	dbl3 c1(vi->fi,vi->fj,vi->fk);
	dbl3 c2(iKid->fi,iKid->fj,iKid->fk);
	if (c1==c2) {c1.z-=0.2;c2.z+=0.2;}
	dbl3 c3(jKid->fi,jKid->fj,jKid->fk);
	if (c2==c3) {c2.y-=0.2;c3.y+=0.2;}
	dbl3 c4(vj->fi,vj->fj,vj->fk);
	//if (c2==c3) {c2.y-=0.2;c3.y+=0.2;}

	cellPoints.push_back(findOrInsertPoint(points, c1));
	if(radius.size()<points.size()) {radius.push_back(vi->R);  cellType.push_back(vi->type); }
	cellPoints.push_back(findOrInsertPoint(points, c2));
	if(radius.size()<points.size()) {radius.push_back(iKid->R); cellType.push_back(iKid->type);}
	cellPoints.push_back(findOrInsertPoint(points, c3));
	if(radius.size()<points.size()) {radius.push_back(jKid->R); cellType.push_back(jKid->type);}
	cellPoints.push_back(findOrInsertPoint(points, c4));
	if(radius.size()<points.size()) {radius.push_back(vj->R); cellType.push_back(vj->type);}

	cellPores.push_back(poreIndx); ///. 3 points each elem
	subTypes.push_back(max(iKid->corId,jKid->corId)); ///. 3 points each elem
//	cellType.push_back(max(iKid->type,jKid->type)); ///. 3 points each elem
	//if(distSqr(vi,vi->boss)< sqrt(vi->boss->limit)*1.2+2-sqrt(vi->limit))  *subTypes.rbegin()=-1;

}



void vtuWriteMbMbs(string basNam, const vector<medialBall*>& ballSpace, const  vector<poreNE*> poreIs, const voxelField<int>&  VElems, double dx, dbl3 X0 )  {

	vector<dbl3> points;
	vector<int> subTypes;
	vector<int> cellType;
	vector<int> cellPores;
	vector<float> radius;

	vector<int> cellPoints;
	//vector<FacePoints> facePoints;
	//vector<CellFaces> cellFaces;
	//vector<FaceCells> faceCells;
    points.reserve(ballSpace.size()*2);
    cellPoints.reserve(ballSpace.size()*2);
    subTypes.reserve(ballSpace.size()*2);
    cellType.reserve(ballSpace.size()*2);
    cellPores.reserve(ballSpace.size()*2);
    radius.reserve(ballSpace.size()*2);


	vector<medialBall*>::const_iterator vi = ballSpace.begin();
	vector<medialBall*>::const_iterator vend = ballSpace.end();
	while (vi<vend)  {
		if ((*vi)->boss != nullptr)  {
			medialBall* mvi=(*vi)->mastrSphere();
			int id=VElems(mvi->fi+1, mvi->fj+1, mvi->fk+1);
			addMbMbMesh(*vi, (*vi)->boss , id,points,subTypes,cellType,cellPores,radius,cellPoints);

			//for (int i=0;i<vi->nKids;++i)
			//{
				//if(sqrt(distSqr(vi->kids[i],&*vi))< sqrt(vi->limit)*1.1+1-sqrt(vi->kids[i]->limit))
					//{addMbMbMesh(&*vi, vi->boss , id,points,subTypes,cellPores,radius,cellPoints);*radius.rbegin()=-1;}
				//for (int j=i+1;j<vi->nKids;++j)
				//{
						 //if(sqrt(distSqr(vi->kids[i], vi->kids[j])) < sqrt(vi->kids[i]->limit)*1.1+1-sqrt(vi->kids[j]->limit))
						//{addMbMbMesh(vi->kids[i], vi->kids[j] , id,points,subTypes,cellPores,radius,cellPoints);*radius.rbegin()=-1;}
					//else if(sqrt(distSqr(vi->kids[i], vi->kids[j])) < sqrt(vi->kids[j]->limit)*1.1+1-sqrt(vi->kids[i]->limit))
						//{addMbMbMesh(vi->kids[i], vi->kids[j] , id,points,subTypes,cellPores,radius,cellPoints);*radius.rbegin()=-1;}
				//}
//~
//~
			//}

		}
		++vi;
    }


    ofstream outp(basNam+".vtu");

	outp<<vtkWriter_start(points.size(),subTypes.size());




	outp<<"      <Points>\n";
	outp<<"        <DataArray type = \"Float32\" NumberOfComponents = \"3\" format = \"ascii\">\n";
    for(size_t i = 0; i < points.size(); ++i)  {
        outp << (points[i][0]*dx+X0[0])<< " " << (points[i][1]*dx+X0[1])<< " " << (points[i][2]*dx+X0[2])<< " ";
        if ((i+1)%20 == 0)     outp << "\n";
    }
	outp<<"\n        </DataArray>\n";
	outp<<"      </Points>\n";

	outp<<"      <Cells>\n";///// \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/
	outp<<"        <DataArray type = \"Int32\" Name = \"connectivity\" format = \"ascii\">\n";
    for(size_t i = 0; i < cellPoints.size(); ++i)  {
        outp << cellPoints[i] << " ";
        if ((i+1)%20 == 0)     outp << "\n";
    }
	outp<<"\n        </DataArray>\n";

	outp<<"	<DataArray type = \"Int32\" Name = \"offsets\" format = \"ascii\">\n";
    for(size_t i = 0; i < subTypes.size(); ++i)  {
        outp << 2*i+2 << " ";
        if ((i+1)%20 == 0)     outp << "\n";
    }
	outp<<"\n        </DataArray>\n";


	outp<<"	<DataArray type = \"UInt8\" Name = \"types\" format = \"ascii\">\n";
    for(size_t i = 0; i < subTypes.size(); ++i)  {
        outp << 3 << " ";
        if ((i+1)%20 == 0)     outp << "\n";
    }
	outp<<"\n        </DataArray>\n";
	outp<<"      </Cells>\n";// //  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//


	outp<<"      <CellData Scalars = \"alpha\">\n"; //////////////////////////////////////


	writeCellData( outp, "index",  cellPores, "Int32");
//	writeCellData( outp, "radius",  radius);
	writeCellData( outp, "subType",  subTypes);
//	writeCellData( outp, "cellType",  cellType);
	outp<<"      </CellData>\n"; /////////////////////////////////////////////////////////

	outp<<"      <PointData>\n"; //////////////////////////////////////
        outp<<"        <DataArray type = \""<<"Float32"<<"\" Name = \""<<"radius"<<"\" format = \"ascii\">\n";
        for(size_t i = 0; i < radius.size(); ++i) { outp<<radius[i]<<" ";  if((i+1)%20 == 0) outp << "\n"; }
        outp<<"        </DataArray>"<<endl;

        outp<<"        <DataArray type = \""<<"Float32"<<"\" Name = \""<<"cellType"<<"\" format = \"ascii\">\n";
        for(size_t i = 0; i < cellType.size(); ++i) { outp<<cellType[i]<<" ";  if((i+1)%20 == 0) outp << "\n"; }
        outp<<"        </DataArray>"<<endl;
	outp<<"      </PointData>\n"; /////////////////////////////////////////////////////////



    outp<<vtkWriter_finish();
    outp.close();
}




void vtuWriteMedialSurface(string basNam,
 const vector<medialBall*>& ballSpace,
 const  vector<poreNE*> poreIs,
 const voxelField<int>&  VElems, double dx, dbl3 X0
 )  {

	vector<dbl3> points;
	vector<int> subTypes;
	vector<int> cellPores;
	vector<int> cellType;
	vector<float> radius;

	vector<int> cellPoints;
	//vector<FacePoints> facePoints;
	//vector<CellFaces> cellFaces;
	//vector<FaceCells> faceCells;
    points.reserve(ballSpace.size()*2);
    cellPoints.reserve(ballSpace.size()*2);
    subTypes.reserve(ballSpace.size()*2);
    cellType.reserve(ballSpace.size()*2);
    cellPores.reserve(ballSpace.size()*2);
    radius.reserve(ballSpace.size()*2);


	//vector<medialBall*>::const_iterator vi = ballSpace.begin();
	//vector<medialBall*>::const_iterator vend = ballSpace.end();
	//while (vi<vend)  {
		//if ((*vi)->boss != nullptr)  {
			//medialBall* mvi=(*vi)->mastrSphere();
			//int id=VElems(mvi->fi+1, mvi->fj+1, mvi->fk+1);

			//for (int i=0;i<(*vi)->nNeis;++i)  {
					//{addMbMbMbMesh((*vi)->boss, *vi, (*vi)->neis[i], id,points,subTypes,cellType,cellPores,radius,cellPoints);}
					//if((*vi)->boss!=(*vi)->neis[i]->boss)  {addMbMbMbMesh((*vi)->boss, (*vi)->neis[i]->boss, (*vi), (*vi)->neis[i], id,points,subTypes,cellType,cellPores,radius,cellPoints);}
			//}
		//}
		//++vi;
    //}


    ofstream outp(basNam+".vtu");

	outp<<vtkWriter_start(points.size(),subTypes.size());




	outp<<"      <Points>\n";
	outp<<"        <DataArray type = \"Float32\" NumberOfComponents = \"3\" format = \"ascii\">\n";
    for(size_t i = 0; i < points.size(); ++i)  {
        outp << (points[i][0]*dx+X0[0])<< " " << (points[i][1]*dx+X0[1])<< " " << (points[i][2]*dx+X0[2])<< " ";
        if ((i+1)%20 == 0)     outp << "\n";
    }
	outp<<"\n        </DataArray>\n";
	outp<<"      </Points>\n";

	outp<<"      <Cells>\n";///// \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\.
	outp<<"        <DataArray type = \"Int32\" Name = \"connectivity\" format = \"ascii\">\n";
    for(size_t i = 0; i < cellPoints.size(); ++i)  {
        outp << cellPoints[i] << " ";
        if ((i+1)%20 == 0)     outp << "\n";
    }
	outp<<"\n        </DataArray>\n";

	outp<<"	<DataArray type = \"Int32\" Name = \"offsets\" format = \"ascii\">\n";
    for(size_t i = 0; i < subTypes.size(); ++i)  {
        outp << 4*i+4 << " ";
        if ((i+1)%20 == 0)     outp << "\n";
    }
	outp<<"\n        </DataArray>\n";


	outp<<"	<DataArray type = \"UInt8\" Name = \"types\" format = \"ascii\">\n";
    for(size_t i = 0; i < subTypes.size(); ++i)  {
        outp << 5 << " ";
        if ((i+1)%20 == 0)     outp << "\n";
    }
	outp<<"\n        </DataArray>\n";
	outp<<"      </Cells>\n";// //  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//




	outp<<"      <CellData Scalars = \"alpha\">\n"; //////////////////////////////////////



	writeCellData( outp, "index",  cellPores, "Int32");
//	writeCellData( outp, "radius",  radius);
	writeCellData( outp, "subType",  subTypes);
//	writeCellData( outp, "cellType",  cellType);
	outp<<"      </CellData>\n"; /////////////////////////////////////////////////////////

	outp<<"      <PointData>\n"; //////////////////////////////////////
        outp<<"        <DataArray type = \""<<"Float32"<<"\" Name = \""<<"radius"<<"\" format = \"ascii\">\n";
        for(size_t i = 0; i < radius.size(); ++i) { outp<<radius[i]<<" ";  if((i+1)%20 == 0) outp << "\n"; }
        outp<<"        </DataArray>"<<endl;

        outp<<"        <DataArray type = \""<<"Float32"<<"\" Name = \""<<"cellType"<<"\" format = \"ascii\">\n";
        for(size_t i = 0; i < cellType.size(); ++i) { outp<<cellType[i]<<" ";  if((i+1)%20 == 0) outp << "\n"; }
        outp<<"        </DataArray>"<<endl;
	outp<<"      </PointData>\n"; /////////////////////////////////////////////////////////



    outp<<vtkWriter_finish();
    outp.close();
}





void vtuWriteThroatMbMbs(string basNam,      /// ********* throat *********
 const vector<throatNE*>& throatIs,
 const  vector<poreNE*> poreIs,
 const voxelField<int>&  VElems, double dx, dbl3 X0
 )  {

	vector<dbl3> points;
	vector<int> subTypes;
	vector<int> cellType;
	vector<int> cellPores;
	vector<float> radius;

	vector<int> cellPoints;
	//vector<FacePoints> facePoints;
	//vector<CellFaces> cellFaces;
	//vector<FaceCells> faceCells;
    points.reserve(throatIs.size()*20);
    cellPoints.reserve(throatIs.size()*20);
    subTypes.reserve(throatIs.size()*20);
    cellType.reserve(throatIs.size()*20);
    cellPores.reserve(throatIs.size()*20);
    radius.reserve(throatIs.size()*20);

	int E2msfc=0, E_msfd=0, E8msgl=0, Dkfhj2=0, Dkfhj1=0;

	vector<throatNE*>::const_iterator tbgn = throatIs.begin();
	vector<throatNE*>::const_iterator ti = throatIs.begin();
	vector<throatNE*>::const_iterator tend = throatIs.end();
	//int poreIndex(0);
	while (ti<tend)  {




		const medialBall * mb2 = (*ti)->mb22();
		const medialBall * mb1 = (*ti)->mb11();
		if (mb1 != nullptr)
			addMbMbMesh(mb1, mb2 ,(ti-tbgn)+poreIs.size()*10,points,subTypes,cellType,cellPores,radius,cellPoints);
		else if ((*ti)->e1>2) cout<<" EioMbk ";

		medialBall* mvi=mb2->mastrSphere();
		if ((*ti)->e2 != VElems(mvi->fi+1, mvi->fj+1, mvi->fk+1)) ++Dkfhj1;//cout<<" Dkfhj1  "<< VElems[mvi->k+1][mvi->j+1][mvi->i+1]<<"   ";
		if (mb2->boss == mb2) ++E_msfd;

		while (mb2->boss != mb2)  {

			addMbMbMesh( mb2 , mb2->boss , (*ti)->e2,points,subTypes,cellType,cellPores,radius,cellPoints);
			mb2 = mb2->boss;

		}

		if (mb1 != nullptr && mb1->boss == mb1) ++E2msfc;
		while ( mb1 != nullptr && mb1->boss != mb1)  {

			mvi=mb1->mastrSphere();
			if ((*ti)->e1 != VElems(mvi->fi+1, mvi->fj+1, mvi->fk+1)) ++Dkfhj2;//cout<<" Dkfhj2  "<< VElems[mvi->k+1][mvi->j+1][mvi->i+1]<<"   ";

			addMbMbMesh(mb1 , mb1->boss , (*ti)->e1,points,subTypes,cellType,cellPores,radius,cellPoints);
			mb1 = mb1->boss;

		}


		if (mb1==mb2) ++E8msgl;//cout<<"\n E8msgl "<<mb1->limit<<" "<<((mb1->i == (*ti)->mb1->i)&&(mb1->j == (*ti)->mb1->j))<<"     ";


		++ti;
    }

	if (E2msfc+E_msfd+Dkfhj2+E8msgl+Dkfhj1)
	   cout<<"\n ******************** bads mb1->boss==mb1: "<<E2msfc<<"   mb2->boss==mb2 = "<<E_msfd<<"   e1!=VElem = "<<Dkfhj2<<"   mb1==mb2 = "<<E8msgl<<"   e2!=VElem = "<<Dkfhj1<<endl;

	ofstream outp(basNam+".vtu");

	outp<<vtkWriter_start(points.size(),subTypes.size());




	outp<<"      <Points>\n";
	outp<<"        <DataArray type = \"Float32\" NumberOfComponents = \"3\" format = \"ascii\">\n";
    for(size_t i = 0; i < points.size(); ++i)  {
        outp << (points[i][0]*dx+X0[0])<< " " << (points[i][1]*dx+X0[1])<< " " << (points[i][2]*dx+X0[2])<< " ";
        if ((i+1)%20 == 0)     outp << "\n";
    }
	outp<<"\n        </DataArray>\n";
	outp<<"      </Points>\n";

	outp<<"      <Cells>\n";///// \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/
	outp<<"        <DataArray type = \"Int32\" Name = \"connectivity\" format = \"ascii\">\n";
    for(size_t i = 0; i < cellPoints.size(); ++i)  {
        outp << cellPoints[i] << " ";
        if ((i+1)%20 == 0)     outp << "\n";
    }
	outp<<"\n        </DataArray>\n";

	outp<<"	<DataArray type = \"Int32\" Name = \"offsets\" format = \"ascii\">\n";
    for(size_t i = 0; i < subTypes.size(); ++i)  {
        outp << 2*i+2 << " ";
        if ((i+1)%20 == 0)     outp << "\n";
    }
	outp<<"\n        </DataArray>\n";

	outp<<"	<DataArray type = \"UInt8\" Name = \"types\" format = \"ascii\">\n";
    for(size_t i = 0; i < subTypes.size(); ++i)  {
        outp << 3 << " ";
        if ((i+1)%20 == 0)     outp << "\n";
    }
	outp<<"\n        </DataArray>\n";
	outp<<"      </Cells>\n";// //  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//


	outp<<"      <PointData>\n"; //////////////////////////////////////
        outp<<"        <DataArray type = \""<<"Float32"<<"\" Name = \""<<"radius"<<"\" format = \"ascii\">\n";
        for(size_t i = 0; i < radius.size(); ++i) { outp<<radius[i]<<" ";  if((i+1)%20 == 0) outp << "\n"; }
        outp<<"        </DataArray>"<<endl;

        outp<<"        <DataArray type = \""<<"Float32"<<"\" Name = \""<<"cellType"<<"\" format = \"ascii\">\n";
        for(size_t i = 0; i < cellType.size(); ++i) { outp<<cellType[i]<<" ";  if((i+1)%20 == 0) outp << "\n"; }
        outp<<"        </DataArray>"<<endl;
	outp<<"      </PointData>\n"; /////////////////////////////////////////////////////////



	outp<<"      <CellData Scalars = \"alpha\">\n"; //////////////////////////////////////

	writeCellData( outp, "index",  cellPores, "Int32");
//	writeCellData( outp, "radius",  radius);
	writeCellData( outp, "subType",  subTypes);
	outp<<"      </CellData>\n"; /////////////////////////////////////////////////////////




    outp<<vtkWriter_finish();
    outp.close();
}










//---------------------------------------------------------------------'
