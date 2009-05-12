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



#include "DM_noLaserField.hh"
#include "RungeKuttaTracker.hh"
#include "OrbitConst.hh"
#include "LorentzTransformationEM.hh"






#define pop(i,n) PopAttr->attArr(i)[n]		//i-part index, n,m-attr index


#define x0(i) 	 Coords->attArr(i)[0]
#define y0(i)  	 Coords->attArr(i)[2]
#define z0(i)  	 Coords->attArr(i)[4]
#define px0(i)   Coords->attArr(i)[1]
#define py0(i)   Coords->attArr(i)[3]
#define pz0(i)   Coords->attArr(i)[5]





using namespace LaserStripping;
using namespace OrbitUtils;


DM_noLaserField::DM_noLaserField(HydrogenStarkParam* Stark)
{
	setName("unnamed");
	

	StarkEffect=Stark;

	levels=	StarkEffect->getStates()*(1+StarkEffect->getStates())*(1+2*StarkEffect->getStates())/6;
	


//allocating memory for koefficients of 4-th order Runge-Kutta method and other koeeficients of the master equation
k_RungeKutt=new double*[levels+1]; for (int i=0;i<levels+1;i++) k_RungeKutt[i]=new double[5];
gamma_ij=new double*[levels+1];	for (int i=0;i<levels+1;i++)	gamma_ij[i]=new double[levels+1];
Gamma_i=new double[levels+1];
E_i=new double[levels+1];


if(StarkEffect->getPyWrapper() != NULL){
		Py_INCREF(StarkEffect->getPyWrapper());
	}


}

DM_noLaserField::~DM_noLaserField()
{

	for (int i=0;i<levels+1;i++) delete [] k_RungeKutt[i]; 	delete	[]	k_RungeKutt;
	delete [] E_i;
	delete [] Gamma_i;
	for (int i=0;i<levels+1;i++)	delete	[]	gamma_ij[i];	delete	[]	gamma_ij;

	
	

	
	if(StarkEffect->getPyWrapper() == NULL){
		delete StarkEffect;
	} else {
		Py_XDECREF(StarkEffect->getPyWrapper());
	}

	

}

void DM_noLaserField::CalcPopulations(int i, Bunch* bunch)	{
		
	PopAttr->attArr(i)[0] = 0;
	for (int j=1; j<levels + 1;j++)	
	PopAttr->attArr(i)[0] += PopAttr->attArr(i)[j];


	
}




void DM_noLaserField::setupEffects(Bunch* bunch){	
	
	
	if(bunch->hasParticleAttributes("Amplitudes")==1)	{
		bunch->removeParticleAttributes("Amplitudes");
	}

		if(bunch->hasParticleAttributes("Populations")==0)	{
			std::map<std::string,double> part_attr_dict;
			part_attr_dict["size"] = levels+1;
			bunch->addParticleAttributes("Populations",part_attr_dict);
			
			for (int i=0; i<bunch->getSize();i++)
			bunch->getParticleAttributes("Populations")->attValue(i,1) = 1;
		}
	
	
	
	
	if (bunch->hasParticleAttributes("pq_coords")==0)	{
	std::map<std::string,double> part_attr_dict;
	part_attr_dict["size"] = 6;
	bunch->addParticleAttributes("pq_coords",part_attr_dict);
	}
	
	Coords = bunch->getParticleAttributes("pq_coords");
	PopAttr = bunch->getParticleAttributes("Populations");
	
	
	for (int i=0; i<bunch->getSize();i++)
	CalcPopulations(i, bunch);


}
		

void DM_noLaserField::memorizeInitParams(Bunch* bunch){
	
	for (int i=0; i<bunch->getSize();i++)		{
		x0(i)=bunch->coordArr()[i][0];
		y0(i)=bunch->coordArr()[i][2];
		z0(i)=bunch->coordArr()[i][4];
		
		px0(i)=bunch->coordArr()[i][1];
		py0(i)=bunch->coordArr()[i][3];
		pz0(i)=bunch->coordArr()[i][5];

	}
	
	
}
	
void DM_noLaserField::finalizeEffects(Bunch* bunch){
	
	if(bunch->hasParticleAttributes("pq_coords")==1)
		bunch->removeParticleAttributes("pq_coords");
	
}







void DM_noLaserField::applyEffects(Bunch* bunch, int index, 
	                            double* y_in_vct, double* y_out_vct, 
														  double t, double t_step, 
														  BaseFieldSource* fieldSource,
															RungeKuttaTracker* tracker)			{




	
		for (int i=0; i<bunch->getSize();i++)	{

			
			//	This function gives parameters Ez_stat	Ex_las[1...3]	Ey_las[1...3]	Ez_las[1...3]	
			//in natural unts (Volt per meter)	in the frame of particle				
			GetParticleFrameFields(i, t,t_step, bunch,fieldSource);
	
			//	This function gives parameters Ez_stat	Ex_las[1...3]	Ey_las[1...3]	Ez_las[1...3]	t_part	omega_part	part_t_step 
			//in atomic units in frame of particle		
			GetParticleFrameParameters(i,t,t_step,bunch);	

			//	This function provides step solution for density matrix using rk4	method
			AmplSolver4step(i,bunch);	
		
			CalcPopulations(i, bunch);
		
		}	

//	cout<<scientific<<setprecision(20)<<bunch->x(0)<<"\t"<<bunch->y(0)<<"\t"<<bunch->z(0)<<"\n";

	
}
















void DM_noLaserField::GetParticleFrameFields(int i,double t,double t_step,  Bunch* bunch,  BaseFieldSource* fieldSource)	{
	
	double** xyz = bunch->coordArr();
	double Ez;
	
		fieldSource->getElectricMagneticField(x0(i),y0(i),z0(i),t,Ex_stat,Ey_stat,Ez_stat,Bx_stat,By_stat,Bz_stat);		
		LorentzTransformationEM::transform(bunch->getMass(),px0(i),py0(i),pz0(i),Ex_stat,Ey_stat,Ez_stat,Bx_stat,By_stat,Bz_stat);

		
		
	
Ez=sqrt(Ex_stat*Ex_stat+Ey_stat*Ey_stat+Ez_stat*Ez_stat);
	
	Ex_stat=0;	
	Ey_stat=0;		
	Ez_stat=Ez;
	
			
	
}








void	DM_noLaserField::GetParticleFrameParameters(int i, double t,double t_step, Bunch* bunch)	{
	
		
	double ta=2.418884326505e-17;			//atomic unit of time
	double Ea=5.14220642e011;				//Atomic unit of electric field
	double m=bunch->getMass();

//This line calculates relyativistic factor-Gamma
double gamma=sqrt(m*m+px0(i)*px0(i)+py0(i)*py0(i)+pz0(i)*pz0(i))/m;
double coeff=1./(gamma*ta);

part_t_step=t_step*coeff;	//time step in frame of particle (in atomic units)
t_part=t*coeff;	


Ex_stat/=Ea;	
Ey_stat/=Ea;		
Ez_stat/=Ea;



StarkEffect->SetE(Ez_stat);	
	
}














void DM_noLaserField::AmplSolver4step(int i, Bunch* bunch)	{
	
	
	

	tcomplex mu_x,mu_y,mu_z;
	double dt,sum;
		

	
	
	
	
			for(int n=1; n<levels+1;n++)	
				StarkEffect->GetEnergyAutoionization(n,E_i[n], Gamma_i[n]);
			

					
			for(int n=2; n<levels+1;n++)
			for(int m=1; m<n;m++)	{
				
				StarkEffect->GetDipoleTransition(n,m,mu_x,mu_y,mu_z);
												
				StarkEffect->GetRelax(n,m,gamma_ij[n][m]);

			}
			



		

		for(int j=1; j<5; j++)
		for(int m=1; m<levels+1;m++)	{	
			
			if (j==4)	dt=part_t_step;	else dt=part_t_step/2.;
			
			k_RungeKutt[m][j]=0;
			
			for(int k=m+1;k<levels+1;k++)	k_RungeKutt[m][j]+=gamma_ij[k][m]*(pop(i,k)+k_RungeKutt[k][j-1]*dt); 
			sum=0; 
			for(int k=1;k<m;k++)			sum+=gamma_ij[m][k];
			
			sum+=Gamma_i[m];
			k_RungeKutt[m][j]-=(pop(i,m)+k_RungeKutt[m][j-1]*dt)*sum;


					
		}


		
		
		
		
	for(int m=1;m<levels+1;m++)	
		pop(i,m) += part_t_step*(k_RungeKutt[m][1]+2.*k_RungeKutt[m][2]+2.*k_RungeKutt[m][3]+k_RungeKutt[m][4])/6.;	


		

	
}




















