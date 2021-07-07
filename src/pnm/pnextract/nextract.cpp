/**************************************************************************************


pnextract algorithm version 0.3, March. 2017

For more information, please visit Imperial College pore-scale modelling website:
https://www.imperial.ac.uk/earth-science/research/research-groups/pore-scale-modelling/
or contact Ali Q. Raeini by email: a.qaseminejad-raeini@imperial.ac.uk
***************************************************************************************/

#include <array>



#ifndef RELEASE_DATE
 #define RELEASE_DATE  __DATE__
#endif


#ifndef MAIN
#define _InitGlobals
#endif

#include "blockNet.h"
#include "writers.h"
#include "profilers.h"

using namespace std;

void GNMAddCorrelations(string netName, dbl2 minRtForDNS, dbl2 scaleCrlByDNS);




#ifndef MAIN

inline void usage(int detailed=1)  {

	cout<<"Pore Network Extraction: pnextract version " << RELEASE_DATE << endl;
	if(detailed)  {
		cout<<"\nUsage:"
		      "\n  pnextract vxlImage.mhd    #  extract network"
		      "\n  pnextract -g vxlImage.mhd # -generate vxlImage.mhd \n"<< endl;

		cout<<" For more information, please visit Imperial College pore-scale modelling website:\n"
		      "https://www.imperial.ac.uk/earth-science/research/research-groups/pore-scale-modelling\n"
			   "or contact Ali Q. Raeini: a.q.raeini@imperial.ac.uk"<<endl;
	}
}

int main(int argc, char* argv[])  {

	usage(0);

	string arg1;
	if(argc>1)  arg1 = argv[1];
	else  {     cout << "Please input data file: ";    getline(cin, arg1);    }
	if (arg1.empty()) { usage(); arg1 = "vxlImage.mhd"; }
	if (arg1=="-h") { usage(); exit(0);}
	if (arg1=="-g") return createSample_input_nextract(string(argc>2? argv[2] : ""), arg1);
	srand(1001);
	try
	{
			auto fnam=arg1;
			if (hasExt(arg1,".am") || hasExt(arg1,".tif") || hasExt(arg1,".gz") || hasExt(arg1,".raw")) arg1="";
			inputDataNE cfg(arg1);
			if(arg1.empty()) cfg.set("ElementDataFile",fnam);
			
				nextract(cfg,true);

	}
	catch (exception &exc) {  cerr<<"\n\n Error in pnextract: \n   "<<exc.what()<<"\n Aborting! \n" <<endl;	return 1; }
	catch (...)            {  cerr<<"\n\n Unknown Error in pnextract! \n Aborting! \n"<<endl;	return 1; }

	return 0;
}

#endif

int nextract(inputDataNE& cfg, bool verbose)  {
	
	if (!cfg.getOr("overwrite", true) && ifstream(cfg.netName()+".gnm")) { // overwrite by default,  sync with WRITENETWORKG(), sync with tests
		cout<<"\n File "+cfg.netName()+".gnm"+" exists, delete it to regenerate, or set keyword `overwrite T;`\n"<<endl;  return 2; }

	(cout<< "InputData: ").flush();
	cfg.echoKeywords(std::cout);

		if (cfg.getOr("write_all", false))  { // use `write_all` as a rememberable alternative for all other visualization keywords
			cfg.add("write_radius","true");
			cfg.add("write_statistics","true");
			cfg.add("write_elements","true");
				 cfg.add("write_poreMaxBalls","true");
				 cfg.add("write_throatMaxBalls","true");
				 cfg.add("write_throats","true");
			//cfg.add("write_poroats","true");   leads to seg fault
			cfg.add("write_hierarchy","true");
			cfg.add("write_throatHierarchy","true");
			cfg.add("write_vtkNetwork","true");
			cfg.add("write_cnm","true");
			cfg.add("write_fullThroats","true");
		}
		if (cfg.getOr("write_Statoil", false) || cfg.getOr("write_StatoilFormat", false)) cfg.add("write_cnm","true");


	Timing tim;
															tim("Init");
	cfg.init(verbose); // process input file
	cfg.readImage(); // read image 
	cfg.createSegments();  // RLE compress image 

															tim("createMedialSurface");
	medialSurface* srf;
	blockNetwork mpn(srf, cfg);
	mpn.createMedialSurface(srf, cfg,0);

															tim("write");
			if (cfg.getOr("write_radius", false))	ballRadiiToVoxel(mpn).writeBin(cfg.name()+"_radius"+cfg.imgfrmt);

															tim("CreateVElem");
	mpn.CreateVElem(0);

															tim("write");
			if (cfg.getOr("write_elements", true))	mpn.VElems.write(cfg.name()+"_VElems.mhd");

															tim("write");
	mpn.createNewThroats(srf);



		mpn.writePNM();

																tim("write");
			if (cfg.getOr("write_hierarchy", false))
																vtuWriteMbMbs(cfg.name()+"_mbHierarchy"+_s(0), srf->ballSpace,  mpn.poreIs,  mpn.VElems, cfg.vxlSize, mpn.VElems.X0()+mpn.VElems.dx());

			//if (cfg.getOr("write_poroats", false))         VElemsPlusThroats(mpn).writeBin(cfg.name()+"_poroats"+cfg.imgfrmt,1,cfg.nx,1,cfg.ny,1,cfg.nz);
			if (cfg.getOr("write_throatHierarchy", false)) vtuWriteThroatMbMbs(cfg.name()+"_throatHierarchy", mpn.throatIs,  mpn.poreIs,  mpn.VElems, cfg.vxlSize,mpn.VElems.X0()+mpn.VElems.dx());
			if (cfg.getOr("write_vtkNetwork", false))      vtuWritePores(cfg.name()+"_pores",  mpn.poreIs, mpn.throatIs, cfg.vxlSize, mpn.VElems.X0()+mpn.VElems.dx());
			if (cfg.getOr("write_vtkNetwork", false))      vtuWriteTHroatSpheres(cfg.name()+"_throatsBalls",  mpn.poreIs, mpn.throatIs, cfg.vxlSize, mpn.VElems.X0()+mpn.VElems.dx());



		int outputBlockSize = 0; /// keywords write_throats, write_poreMaxBalls and write_throatMaxBalls developed by Tom Bultreys
		if(cfg.giv("outputBlockSize", outputBlockSize)) cout << "OutputBlockSize:" << outputBlockSize << endl;
        if (!outputBlockSize) {
            if (cfg.getOr("write_throats", false))             VThroats(mpn).writeBin(cfg.name()+"_throats"+cfg.imgfrmt,0,cfg.nx,0,cfg.ny,0,cfg.nz);
            if (cfg.getOr("write_poreMaxBalls", false))        poreMaxBalls(mpn).writeBin(cfg.name()+"_poreMBs"+cfg.imgfrmt,0,cfg.nx,0,cfg.ny,0,cfg.nz);
            if (cfg.getOr("write_throatMaxBalls", false))      throatMaxBalls(mpn).writeBin(cfg.name()+"_throatMBs"+cfg.imgfrmt,0,cfg.nx,0,cfg.ny,0,cfg.nz);
         } else {
            int blockNumber = 1,  beginSlice = 0,  endSlice = 0;
            while (endSlice < cfg.nz-1){
                cout << " WRITING BLOCK \n";
                beginSlice = (blockNumber-1) * outputBlockSize;
                endSlice = min(blockNumber * outputBlockSize, cfg.nx-1);
                if (cfg.getOr("write_throats", false))         VThroats(mpn, beginSlice, endSlice).writeBin(cfg.name()+"_throats" + _s(blockNumber) + cfg.imgfrmt,0 , cfg.nx, 0,cfg.ny,0, endSlice-beginSlice);
                if (cfg.getOr("write_poreMaxBalls", false))    poreMaxBalls(mpn, beginSlice, endSlice).writeBin(cfg.name()+"_poreMBs" + _s(blockNumber) + cfg.imgfrmt, 0 , cfg.nx, 0,cfg.ny,0, endSlice-beginSlice);
                if (cfg.getOr("write_throatMaxBalls", false))  throatMaxBalls(mpn, beginSlice, endSlice).writeBin(cfg.name()+"_throatMBs" + _s(blockNumber) + cfg.imgfrmt,0 , cfg.nx, 0,cfg.ny,0, endSlice-beginSlice);
                blockNumber ++;
            }
		}

																tim("write");
																if (cfg.getOr("write_vtkNetwork", false))	vtuWriteThroats(cfg.name()+"_throats",  mpn.poreIs, mpn.throatIs, cfg.vxlSize, mpn.VElems.X0()+mpn.VElems.dx());

	cout<<endl<<cfg.name()<<endl;
	cout<<"***  " <<mpn.poreIs.size()<<"-"<<mpn.nBP6<<" pores, "<<mpn.throatIs.size()<<" throats,   ratio: "<<double(mpn.throatIs.size())/(mpn.poreIs.size()-6.+1e-6)<<"  ***"<<endl;
	cout<<"end"<<endl;


 return 0;
}

