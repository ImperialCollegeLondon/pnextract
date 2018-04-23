
#include "blockNet.h"






inline double randomG () // to randomly distribute the shape factors
{
	double x1, x2, w, y;
	do{
		x1 = 2.0 * rand()/RAND_MAX - 1.0;
		x2 = 2.0 * rand()/RAND_MAX - 1.0;
		w = x1 * x1 + x2 * x2;
	} while (w >= 1.0);
	w = sqrt ( (-2.0 * log(w)) / w);
	y = 0.00625 * (x1 * w + 5.0);
	if (y > 0.049)	y = 0.0625;
	return y;
}


void  blockNetwork::writeStatoilFormat() const
{

	vector<double> t_radiuss(nTrots,0.0);//
	vector<double> t_shapeFacts(nTrots,0.0);//
	vector<double> t_lengthP1toP2s(nTrots,0.0);
	vector<double> t_lp1s(nTrots,0.0);//
	vector<double> t_lp2s(nTrots,0.0);//
	vector<double> t_physlength(nTrots,0.0);
	
	vector<double> p_radiuss(nPores,0.0);
	vector<double> p_shape1s(nPores,0.0);
	vector<double> p_physlengths(nPores,0.0);


	cout<<"\ncalcThroats:"<<endl;
	int lengthP1toP2Warnings = 0;
	double nBelowAllowedG(0.0), nAboveAllowedG(0.0), totalArea(0.0);


	for (int tid=0; tid<nTrots; ++tid)
	{
	    throatElementI* trot = throatIs[tid];

		double lthroat = 0;
		double lp1 = 0;
		double lp2 = 0;


		if (trot->surfaceArea == 0)		trot->surfaceArea = 6;
		if (mag(trot->CrosArea) <0.01)		trot->CrosArea[0] = 0.1;


		double lpt1 = (   (trot->e1 <2)  ?  (trot->e1 == 0 ? trot->mb22()->fi:cg.nx-trot->mb22()->fi) : dist (poreIs[trot->e1]->mb, trot->mb22())   );
		double rp1 = std::max(   (trot->e1 <2 ) ? trot->mb22()->R : poreIs[trot->e1]->mb->R , 1.0f );
		double lpt2 = (   (trot->e2 <2)  ?  (trot->e2 == 0?trot->mb22()->fi:cg.nx-trot->mb22()->fi) : dist (poreIs[trot->e2]->mb, trot->mb22())   );
		double rp2 = std::max(   (trot->e2 <2 ) ? trot->mb22()->R : poreIs[trot->e2]->mb->R , 1.0f );


		double rr = std::max(trot->mb22()->R,0.5f);
		rr = std::min(std::min(rr,rp1),rp2);
		t_radiuss[tid] = rr+0.5*(0.5-double(rand())/RAND_MAX);

		double lengthP1toP2 = lpt1+lpt2;
		if (lengthP1toP2 < 3.0) 	lengthP1toP2 = 3.01, ++lengthP1toP2Warnings;

		lp1 = lpt1*0.67;
		lp2 = lpt2*0.67;
		if (trot->e1 < 2 )	lp1 = 1;
		if (trot->e2 < 2 )	lp2 = 1;
		lthroat = lengthP1toP2-lp1-lp2;


		if (lthroat < 0.0000001) 	lthroat = 1;

		t_shapeFacts[tid] = rr*rr/4.0/mag(trot->CrosArea);  ///. +3 is to cancel the fact that rr is bigger than rt by 0.5 voxels


		if (t_shapeFacts[tid]>= 0.09 )
			{t_shapeFacts[tid] = std::min(0.079,t_shapeFacts[tid]/2.0); nAboveAllowedG += mag(trot->CrosArea);}///. Probably shared throat, temporary fix to handle in the flow code

		if (t_shapeFacts[tid]<0.01)
			{t_shapeFacts[tid] = std::max(randomG(),0.01); nBelowAllowedG += mag(trot->CrosArea);}

		totalArea += mag(trot->CrosArea);

		t_lengthP1toP2s[tid] = lengthP1toP2*1.0; ///. wrong for now
		t_lp1s[tid] = lp1*1.0;
		t_lp2s[tid] = lp2*1.0;
		t_physlength[tid] = lthroat*1;
	}
	cout<<  " P1-to-P2 length < 3    for " <<lengthP1toP2Warnings<<" throats"<<endl;
	cout<<" shapefactor: belowAllowedG "<<nBelowAllowedG/totalArea*100<<"%   aboveAllowedG "<<nAboveAllowedG/totalArea*100<<"%"<<endl;


	 cout<<"calc Pores"<<endl;


	for (int pid=2; pid<nPores; ++pid)
	{
		 poreElementI* por = poreIs[pid];
		double radius = por->mb->R;
		p_radiuss[pid] = radius;
		if (por->surfaceArea<1) por->surfaceArea = 6;
		if (por->volumn<1) por->volumn = 1;



		//double poreLengthMax = radius;
		double shapeFactor(5.0e-38), SumTArea(1.0e-36);
		//double SumTRad(1.0e-36),poreLengthAvg(0.0);
		for (std::map<int,int>::const_iterator bi = por->contacts.begin(); bi != por->contacts.end(); ++bi)
		{
			throatElementI* trot = throatIs[bi->second];
			//double prtrtLength = dist(por->mb, trot->mb22());
			//poreLengthMax = prtrtLength>poreLengthMax ? prtrtLength : poreLengthMax;

			//SumTRad += abs(trot->radius);
			//poreLengthAvg += prtrtLength*abs(trot->radius);
			shapeFactor += t_shapeFacts[bi->second]*mag(trot->CrosArea);
			SumTArea += mag(trot->CrosArea);
		}
		shapeFactor /= SumTArea;
		//poreLengthAvg /= SumTRad;
		//poreLengthAvg = std::max(poreLengthAvg*2.0-radius/3.0,2.0*sqrt(3.0)/3.0*radius); ///. Warning: should have been 2.0
		//double shapeFactor = por->radius*por->radius*poreLengthAvg/4.0/(por->volumn); ///. To compare with the simple G above

		double porA(radius*radius/4.0/shapeFactor);
		SumTArea += porA;
		double pVol = por->volumn;
		por->volumn = pVol * porA/SumTArea;

		for (std::map<int,int>::const_iterator bi = por->contacts.begin(); bi != por->contacts.end(); ++bi)
		{
			throatElementI* trot = throatIs[bi->second];
			trot->volumn += pVol * mag(trot->CrosArea)/SumTArea;
		}

		//if ( shapeFactor>= 0.09 )
			//{shapeFactor = std::max(0.0485,(por->volumn*poreLengthMax*2/por->surfaceArea/por->surfaceArea)); nAboveAllowedG += por->volumn;}///. Probably shared throat, temporary fix to handle in the flow code

		//if (shapeFactor<0.01)
			//{shapeFactor = std::max(por->volumn*poreLengthMax*2/por->surfaceArea/por->surfaceArea,0.01); nBelowAllowedG += por->volumn;}

		p_shape1s[pid] = shapeFactor;//(por->volumn*poreLengthMax*2/por->surfaceArea/por->surfaceArea);
		//~ por->crossSecArea = (por->volumn)/std::max(poreLengthAvg,por->radius);
	}



	const double vxllength = cg.precision;
	cout<<"Writing throats";cout.flush();

	{ ///.  _link1.dat
		FILE* fil = fopen((cg.baseName() + "_link1.dat").c_str(), "w");
		fprintf(fil, "%6d\n", int(throatIs.size()));
		for (int tid = 0; tid < int(throatIs.size()); ++tid)
		{	const throatElementI* trot = throatIs[tid];

			fprintf(fil, "%6d %6d %6d %E %E %E\n", tid+1, int(trot->e1-1), int(trot->e2-1),
						  trot->radius()*vxllength, t_shapeFacts[tid], t_lengthP1toP2s[tid]*vxllength);
		}
		fclose(fil);
	}

	{///.  _link2.dat
		FILE* fil = fopen((cg.baseName() + "_link2.dat").c_str(), "w");
		for (int tid = 0; tid < int(throatIs.size()); ++tid)
		{	const throatElementI* trot = throatIs[tid];

			fprintf(fil, "%6d %6d %6d %E %E %E %E %E\n", tid+1,  int(trot->e1-1), int(trot->e2-1),
			  t_lp1s[tid]*vxllength, t_lp2s[tid]*vxllength, t_physlength[tid]*vxllength, trot->volumn*vxllength*vxllength*vxllength, 0.0);
		}
		fclose(fil);
	}
	cout<<".\n";cout.flush();


	cout<<"Writing pores";cout.flush();
	{///.  _node1.dat
		FILE* fil = fopen((cg.baseName() + "_node1.dat").c_str(), "w");
		fprintf(fil, "%6d  %E  %E  %E\n", int(poreIs.size())-2, cg.nx*vxllength, cg.ny*vxllength, cg.nz*vxllength);

		for (int pid = 2; pid < int(poreIs.size()); ++pid)///. 0th and first are inlet and outlet elements
		{	const poreElementI* por = poreIs[pid];///. 0th and first are inlet and outlet elements

			fprintf(fil, "%6d %E %E %E %3d", int(pid-1), por->mb->fi*vxllength+cg.X0[0], por->mb->fj*vxllength+cg.X0[1], por->mb->fk*vxllength+cg.X0[2], int(por->contacts.size()));

			int inlet = 0,outlet = 0;
			for (std::map<int,int>::const_iterator bi = por->contacts.begin(); bi != por->contacts.end(); ++bi)
			{
				const throatElementI* tb = throatIs[bi->second];
				if (tb->e1 == pid)
				{
					if (tb->e2 == 0)		 inlet = 1;
					else if (tb->e2 == 1)  outlet = 1;

					fprintf(fil, "\t%5d", int(tb->e2-1));
				}
				else
				{
					if (tb->e1 == 0)		 inlet = 1;
					else if (tb->e1 == 1)  outlet = 1;
					fprintf(fil, "\t%5d", int(tb->e1-1));
				}
			}

			fprintf(fil, "\t%5d", (inlet) ); ///. .dummy
			fprintf(fil, "\t%5d", (outlet) ); ///. .dummy

			for (std::map<int,int>::const_iterator bi = por->contacts.begin(); bi != por->contacts.end(); ++bi)
			{
				fprintf(fil, "\t%5d", int(bi->second)+1);
			}
			fputs("\n", fil);
		}
		fclose(fil);
	}

	{///.  _node2.dat
		FILE* fil = fopen((cg.baseName() + "_node2.dat").c_str(), "w");
		for (int pid = 2; pid < int(poreIs.size()); ++pid)///. 0th and first are inlet and outlet elements
		{	const poreElementI* por = poreIs[pid];///. 0th and first are inlet and outlet elements

			fprintf(fil, "%6d %E %E %E %E\n", int(pid-1),  por->volumn*vxllength*vxllength*vxllength,  por->radius()*vxllength, p_shape1s[pid], 0.0);

		}
		fclose(fil);
	}

	cout<<".\n";cout.flush();

}








