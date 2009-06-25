//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    CppExternalEffects.cc
//
// CREATED
//    06/27/2008
//
// DESCRIPTION
//    The base class for C++ implementation of a external effects 
//    during the transport of particles through the external field. 
//    It should be sub-classed on Python level and implement
//    setupEffects(Bunch* bunch)
//    finalizeEffects(Bunch* bunch)
//    applyEffects(Bunch* bunch, int index, 
//	                            double* y_in_vct, double* y_out_vct, 
//														  double t, double t_step, 
//														  OrbitUtils::BaseFieldSource* fieldSource)
//    methods.
//    The results of these methods will be available from the c++ level.
//    This is an example of embedding Python in C++ Orbit level.
//
//     
//         
//       
//     
//
///////////////////////////////////////////////////////////////////////////
#include "orbit_mpi.hh"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

#include "RecordEvolution.hh"
#include "RungeKuttaTracker.hh"



#define x0(i) 	 Coords->attArr(i)[0]
#define y0(i)  	 Coords->attArr(i)[2]
#define z0(i)  	 Coords->attArr(i)[4]
#define px0(i)   Coords->attArr(i)[1]
#define py0(i)   Coords->attArr(i)[3]
#define pz0(i)   Coords->attArr(i)[5]

#define t_step_ev(i)   	Evol->attArr(i)[num_plot+1]
#define x0_ev(i)   		Evol->attArr(i)[num_plot+2]
#define y0_ev(i)   		Evol->attArr(i)[num_plot+3]
#define z0_ev(i)   		Evol->attArr(i)[num_plot+4]

using namespace LaserStripping;
using namespace OrbitUtils;



RecordEvolution::RecordEvolution(char* effect,int ind_effect, int num)
{
	

	setRankSetup(1);
	setRankMemorize(1);
	setRankFinalize(1);
	setRankApply(1);

	num_plot=num;
	index_effect=ind_effect;
	effect_name=effect;


}




RecordEvolution::~RecordEvolution() 
{
}






void RecordEvolution::setupEffects(Bunch* bunch){
		

	setup_par = true;

	
	
	if(bunch->hasParticleAttributes("Evolution")==1)
		bunch->removeParticleAttributes("Evolution");
	
	std::map<std::string,double> part_attr_dict;
	part_attr_dict["size"] = num_plot+5;
	bunch->addParticleAttributes("Evolution",part_attr_dict);
	
	Evol = bunch->getParticleAttributes("Evolution");
	
	
	if (bunch->hasParticleAttributes("pq_coords")==0)	{
	std::map<std::string,double> part_attr_dict;
	part_attr_dict["size"] = 6;
	bunch->addParticleAttributes("pq_coords",part_attr_dict);
	}

	Coords = bunch->getParticleAttributes("pq_coords");
	RecEff = bunch->getParticleAttributes(effect_name);
			

	for (int i=0; i<bunch->getSize();i++)                        
	Evol->attArr(i)[0] = RecEff->attArr(i)[index_effect];

	

	
}



void RecordEvolution::memorizeInitParams(Bunch* bunch){
	
	for (int i=0; i<bunch->getSize();i++)		{
		x0(i)=bunch->coordArr()[i][0];
		y0(i)=bunch->coordArr()[i][2];
		z0(i)=bunch->coordArr()[i][4];
		
		px0(i)=bunch->coordArr()[i][1];
		py0(i)=bunch->coordArr()[i][3];
		pz0(i)=bunch->coordArr()[i][5];

	}
	
	
}
	
		

	
void RecordEvolution::finalizeEffects(Bunch* bunch) {

	if(bunch->hasParticleAttributes("pq_coords")==1)
		bunch->removeParticleAttributes("pq_coords");

}







void RecordEvolution::applyEffects(Bunch* bunch, int index, 
	                            double* y_in_vct, double* y_out_vct, 
														  double t, double t_step, 
														  BaseFieldSource* fieldSource,
															RungeKuttaTracker* tracker)			{

	
	
	

	if (setup_par == true)	{

	Num = tracker->getStepsNumber();
	
	if(Num<num_plot)	{cout<<"The number of evolution points must be equal or less then the number of steps \n"; abort();}
	if(Num%num_plot!=0)	{cout<<"The number of steps must be divisible by the number of evolution point \n"; abort();}
	
	
	for (int i=0; i<bunch->getSize();i++)	{
	x0_ev(i) = x0(i);
	y0_ev(i) = y0(i);
	z0_ev(i) = z0(i);
	t_step_ev(i) = Num*t_step/num_plot;
	}
	
	setup_par = false;
	t_in = t;

	}

	
	


	
if(int((t-t_in+t_step)/t_step+0.5)%(Num/num_plot) == 0)		
for (int i=0; i<bunch->getSize();i++)	{
	Evol->attArr(i)[int((t-t_in+t_step)/t_step+0.5)/(Num/num_plot)] = RecEff->attArr(i)[index_effect];

}




}
