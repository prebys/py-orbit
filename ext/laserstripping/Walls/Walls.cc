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



#include "Walls.hh"
#include "RungeKuttaTracker.hh"







#define x0(i) 	 Coords->attArr(i)[0]
#define y0(i)  	 Coords->attArr(i)[2]
#define z0(i)  	 Coords->attArr(i)[4]

#define x(i) 	 bunch->coordArr()[i][0]
#define y(i)  	 bunch->coordArr()[i][2]
#define z(i)  	 bunch->coordArr()[i][4]





using namespace LaserStripping;
using namespace OrbitUtils;


Walls::Walls()
{

}




Walls::~Walls() 
{

}
















void Walls::setupEffects(Bunch* bunch){
	

	
	if (bunch->hasParticleAttributes("pq_coords")==0)	{
	std::map<std::string,double> part_attr_dict;
	part_attr_dict["size"] = 6;
	bunch->addParticleAttributes("pq_coords",part_attr_dict);
	}

	Coords = bunch->getParticleAttributes("pq_coords");


}



void Walls::memorizeInitParams(Bunch* bunch){
	
	for (int i=0; i<bunch->getSize();i++)		{
		x0(i)=bunch->coordArr()[i][0];
		y0(i)=bunch->coordArr()[i][2];
		z0(i)=bunch->coordArr()[i][4];

	}
	
	
}
		
	
void Walls::finalizeEffects(Bunch* bunch) {
	
	if(bunch->hasParticleAttributes("pq_coords")==1)
		bunch->removeParticleAttributes("pq_coords");
	
}





/*

bool Walls::crossPlane(int i,Bunch* bunch,double xs, double ys,double zs,double vx, double vy,double vz){
	
	double vr = (x0(i)-x(i))*vx+(y0(i)-y(i))*vy+(z0(i)-z(i))*vz;
	
	
	if(vr==0)
		return false;
	
	else
	{
		double coeff = (vx*(xs-x0(i))+vy*(ys-y0(i))+vz*(zs-z0(i)))/vr;
		double x_cross = x0(i)+(x0(i)-x(i))*coeff;
		double y_cross = y0(i)+(y0(i)-y(i))*coeff;
		double z_cross = z0(i)+(z0(i)-z(i))*coeff;
		
		if ((x0(i)<=x_cross<=x(i))&&(y0(i)<=y_cross<=y(i))&&(z0(i)<=z_cross<=z(i)))
			return true;
		else
			return false;
		
		
	}
	

}

*/



bool Walls::crossSurface(int i,Bunch* bunch)	{

	
if(
/**************************		
		((z0(i)<=2.936809)&&(2.936809<=z(i))&&(0.118921>=x(i))&&(x(i)>=0.026719))		//position of thick foil
		
		||
		
		((z0(i)<=4.442887)&&(4.442887<=z(i))&&((x(i)+0.09679)*(x(i)+0.09679)+(y(i)-0)*(y(i)-0)<=0.0103476))	//position of ring flange
***************************/	
		
		((z0(i)<=0.25)&&(0.25<=z(i)))
		

)
	
	return true;
else
	return false;

}







void Walls::applyEffects(Bunch* bunch, int index, 
	                            double* y_in_vct, double* y_out_vct, 
														  double t, double t_step, 
														  BaseFieldSource* fieldSource,
															RungeKuttaTracker* tracker)			{


	
		for (int i=0; i<bunch->getSize();i++)	{

			if(crossSurface(i,bunch))
				bunch->deleteParticleFast(i);
		}	

	
}






