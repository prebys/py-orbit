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

#include "PrintExtEffects.hh"
#include "RungeKuttaTracker.hh"






using namespace LaserStripping;
using namespace OrbitUtils;



PrintExtEffects::PrintExtEffects(int i,char* addr)
{

	
	ORBIT_MPI_Comm_size(MPI_COMM_WORLD, &size_MPI);
	ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank_MPI);
	
	max_print_par=i;
	print_par=max_print_par;
	addr_print=addr;

}




PrintExtEffects::~PrintExtEffects() 
{
}








void PrintExtEffects::setupEffects(Bunch* bunch){		

}
		

	
void PrintExtEffects::finalizeEffects(Bunch* bunch) {
	
}







void PrintExtEffects::applyEffects(Bunch* bunch, int index, 
	                            double* y_in_vct, double* y_out_vct, 
														  double t, double t_step, 
														  BaseFieldSource* fieldSource,
															RungeKuttaTracker* tracker)			{


	PopAttr = (AtomPopulations*) bunch->getParticleAttributes("Populations");

	
		for (int i=0; i<bunch->getSize();i++)	{

			
			
			if (print_par==max_print_par)	{
			if (i==bunch->getSize()-1)	print_par=0;
			sprintf(addr_name,"%s%i.dat",addr_print,i*size_MPI+rank_MPI);		
			ofstream file(addr_name,ios::app);
			
			file<<t<<"\t";
			for (int j=1; j<PopAttr->getAttSize();j++)	
				file<<PopAttr->attArr(i)[j]<<"\t";
			
			file<<PopAttr->attArr(i)[0]<<"\n";
			file.close();
				
			}
			if (i==bunch->getSize()-1)	print_par++;

		
			
/*			
			//This function prints evolution of populations in file
			if (print_par==max_print_par)	{
			if (i==bunch->getSize()-1) print_par=0;
			sprintf(addr_name,"%s%i.dat",addr_print,i*size_MPI+rank_MPI);
			ofstream file(addr_name,ios::app);
			file<<t<<"\t";
			for(int n=1;n<levels+1;n++)	file<<Re(i,n)*Re(i,n)+Im(i,n)*Im(i,n)<<"\t";
			double sum=0;	for(int n=1;n<levels+1;n++)	sum+=Re(i,n)*Re(i,n)+Im(i,n)*Im(i,n);
			file<<sum<<"\n";
			file.close();	
			}
			if ((max_print_par!=-2)&&(i==bunch->getSize()-1))	print_par++;
*/			
			
			
/*			
			if (print_par==max_print_par)	{
			if (i==bunch->getSize()-1) print_par=0;
			sprintf(addr_name,"%s%i.dat",addr_print,i*size_MPI+rank_MPI);
			ofstream file(addr_name,ios::app);
			file<<t<<"\t";
			for(int n=1;n<levels+1;n++)	file<<Re(i,n,n)<<"\t";
			double sum=0;	for(int n=1;n<levels+1;n++)	sum+=Re(i,n,n);
			file<<sum<<"\n";
			file.close();			
			}
			if ((max_print_par!=-2)&&(i==bunch->getSize()-1))	print_par++;
			
			
			
 */
			
			
			
		
		}	


}



























