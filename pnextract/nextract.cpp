/**************************************************************************************


gnextract code version 0.3, March. 2017

For more information, please visit Imperial College pore-scale modelling website:
http://www.imperial.ac.uk/earth-science/research/research-groups/perm/research/pore-scale-modelling/
or contact Ali Q. Raeini by email: a.qaseminejad-raeini@imperial.ac.uk
***************************************************************************************/

#include <array>



#include "blockNet.h"
#include "writers.h"



#ifndef MAIN
#ifndef EXE_VAR
 #define EXE_VAR  5
#endif

inline void usage(int detailed=1)
{

	std::string  netExeVers[] = {"generalized, beta","generalized","x","x","conventional, beta","conventional","xx","xx"};
	std::cout<<"Pore Network Extraction Code, "<<netExeVers[EXE_VAR]<<" version, compile date: " << __DATE__ << std::endl;
	if(detailed)
	{
		std::cout<<" For more information, please visit Imperial College pore-scale modelling website:"<<std::endl
			 <<"http://www.imperial.ac.uk/earth-science/research/research-groups/perm/research/pore-scale-modelling"<<std::endl
			 <<"or contact Ali Q. Raeini by email: a.qaseminejad-raeini@imperial.ac.uk"<<endl;
		std::cout<<"\nUsage:"<< std::endl;
		std::cout<<"  gnextract vxlImage.mhd\n"<< std::endl;
	}
}

int debugLevel=0;
int main(int argc, char* argv[])
{

	usage(0);

	string inputFileName;
	if(argc>1)  inputFileName = argv[1];
	else  {     cout << "Please input data file: ";    getline(cin, inputFileName);    }
	if (inputFileName.empty()) { usage(); inputFileName = "vxlImage.mhd"; }
	if (inputFileName=="-h") { usage(); exit(0);}
	srand(1001);

	inputDataNE cfg(inputFileName);
	return nextract(cfg);
}

#endif

int nextract(inputDataNE& cfg)
{




		if (cfg.getOr(false,"write_all"))
		{
			cfg.addKeyword("write_radius","true");
			cfg.addKeyword("write_statistics","true");
			cfg.addKeyword("write_elements","true");
				 cfg.addKeyword("write_poreMaxBalls","true");
				 cfg.addKeyword("write_throatMaxBalls","true");
				 cfg.addKeyword("write_throats","true");
			cfg.addKeyword("write_poroats","true");
			cfg.addKeyword("write_subElements","true");
			cfg.addKeyword("write_hierarchy","true");
			cfg.addKeyword("write_medialSurface","true");
			cfg.addKeyword("write_throatHierarchy","true");
			cfg.addKeyword("write_corners","true");
			cfg.addKeyword("write_vtkNetwork","true");
			cfg.addKeyword("write_StatoilFormat","true");
		}





	medialSurface* refs;
	blockNetwork mpn(refs,cfg);
	mpn.createBallsAndHierarchy(refs,cfg,0);

			if (cfg.getOr(false,"write_radius"))	ballRadiiToVoxel(mpn).writeBin(cfg.baseName()+"_radius"+cfg.imgfrmt);

	mpn.CreateVElem(0);

			if (cfg.getOr(true,"write_elements"))	mpn.VElems.write(cfg.baseName()+"_VElems.mhd");

	mpn.createNewThroats(refs);



		mpn.writeStatoilFormat();

			if (cfg.getOr(false,"write_hierarchy"))
																vtuWriteMbMbs(cfg.baseName()+"_mbHierarchy"+toStr(0), mpn.allSpace,  mpn.poreIs,  mpn.VElems, cfg.precision, mpn.VElems.X0()+mpn.VElems.dx());
			if (cfg.getOr(false,"write_medialSurface"))
																 vtuWriteMedialSurface(cfg.baseName()+"_medialSurface"+toStr(0), mpn.allSpace,  mpn.poreIs,  mpn.VElems, cfg.precision,mpn.VElems.X0()+mpn.VElems.dx());
			if (cfg.getOr(false,"write_poroats"))        VElemsPlusThroats(mpn).writeBin(cfg.baseName()+"_poroats"+cfg.imgfrmt,1,cfg.nx,1,cfg.ny,1,cfg.nz);
			if (cfg.getOr(false,"write_throatHierarchy")) vtuWriteThroatMbMbs(cfg.baseName()+"_throatHierarchy", mpn.throatIs,  mpn.poreIs,  mpn.VElems, cfg.precision,mpn.VElems.X0()+mpn.VElems.dx());
			if (cfg.getOr(false,"write_vtkNetwork"))	vtuWritePores(cfg.baseName()+"_pores",  mpn.poreIs, mpn.throatIs, cfg.precision, mpn.VElems.X0()+mpn.VElems.dx());
			if (cfg.getOr(false,"write_vtkNetwork"))	vtuWriteTHroatSpheres(cfg.baseName()+"_throatsBalls",  mpn.poreIs, mpn.throatIs, cfg.precision, mpn.VElems.X0()+mpn.VElems.dx());



	int outputBlockSize = 0; //. developed by Tom Bultreys
	if(cfg.getVar(outputBlockSize,"outputBlockSize"))
		cout << "OutputBlockSize:" << outputBlockSize << endl;

        if (!outputBlockSize) {
            if (cfg.getOr(false,"write_throats"))                  VThroats(mpn).writeBin(cfg.baseName()+"_throats"+cfg.imgfrmt,0,cfg.nx,0,cfg.ny,0,cfg.nz);
            if (cfg.getOr(false,"write_poreMaxBalls"))             poreMaxBalls(mpn).writeBin(cfg.baseName()+"_poreMBs"+cfg.imgfrmt,0,cfg.nx,0,cfg.ny,0,cfg.nz);
            if (cfg.getOr(false,"write_throatMaxBalls"))           throatMaxBalls(mpn).writeBin(cfg.baseName()+"_throatMBs"+cfg.imgfrmt,0,cfg.nx,0,cfg.ny,0,cfg.nz);
        } else {

            int blockNumber = 1;
            int beginSlice = 0;
            int endSlice = 0;

            while (endSlice < cfg.nz-1){
                cout << " WRITING BLOCK \n";
                beginSlice = (blockNumber-1) * outputBlockSize;
                endSlice = min(blockNumber * outputBlockSize, cfg.nx-1);
                if (cfg.getOr(false,"write_throats"))                  VThroats(mpn, beginSlice, endSlice).writeBin(cfg.baseName()+"_throats" + toStr(blockNumber) + cfg.imgfrmt,0 , cfg.nx, 0,cfg.ny,0, endSlice-beginSlice);
                if (cfg.getOr(false,"write_poreMaxBalls"))               poreMaxBalls(mpn, beginSlice, endSlice).writeBin(cfg.baseName()+"_poreMBs" + toStr(blockNumber) + cfg.imgfrmt, 0 , cfg.nx, 0,cfg.ny,0, endSlice-beginSlice);
                if (cfg.getOr(false,"write_throatMaxBalls"))           throatMaxBalls(mpn, beginSlice, endSlice).writeBin(cfg.baseName()+"_throatMBs" + toStr(blockNumber) + cfg.imgfrmt,0 , cfg.nx, 0,cfg.ny,0, endSlice-beginSlice);
                blockNumber ++;
            }
		  }

																if (cfg.getOr(false,"write_vtkNetwork"))	vtuWriteThroats(cfg.baseName()+"_throats",  mpn.poreIs, mpn.throatIs, cfg.precision, mpn.VElems.X0()+mpn.VElems.dx());




	cout<<endl<<cfg.baseName()<<endl;
	cout<<"***  " <<mpn.poreIs.size()-2<<" pores, "<<mpn.throatIs.size()<<" throats,   ratio: "<<double(mpn.throatIs.size())/(mpn.poreIs.size()-2.0+1e-6)<<"  ***"<<endl;
	cout<<"end"<<endl;

	return 0;
}

