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






#define Re(i,n,m) AmplAttr->attArr(i)[(n-1)*levels+m]		//i-part index, n,m-attr index
#define Im(i,n,m) AmplAttr->attArr(i)[(n-1)*levels+m+levels*levels]

#define x0(i) 	 Coords->attArr(i)[0]
#define y0(i)  	 Coords->attArr(i)[2]
#define z0(i)  	 Coords->attArr(i)[4]
#define px0(i)   Coords->attArr(i)[1]
#define py0(i)   Coords->attArr(i)[3]
#define pz0(i)   Coords->attArr(i)[5]

#define dm(i,n,m) tcomplex(AmplAttr->attArr(i)[(n-1)*levels+m],AmplAttr->attArr(i)[(n-1)*levels+m+levels*levels])
#define k_rk(j,n,m) k_RungeKutt[j][(n-1)*levels+m]




using namespace LaserStripping;
using namespace OrbitUtils;


DM_noLaserField::DM_noLaserField(HydrogenStarkParam* Stark)
{
	setName("unnamed");
	

	StarkEffect=Stark;

	levels=	StarkEffect->getStates()*(1+StarkEffect->getStates())*(1+2*StarkEffect->getStates())/6;
	


//allocating memory for koefficients of 4-th order Runge-Kutta method and other koeeficients of the master equation
k_RungeKutt=new tcomplex**[levels+1];	for (int i=0;i<levels+1;i++)	k_RungeKutt[i]=new tcomplex*[levels+1]; for (int i=0;i<levels+1;i++) for (int j=0;j<levels+1;j++)	k_RungeKutt[i][j]=new tcomplex[5];
gamma_ij=new double*[levels+1];	for (int i=0;i<levels+1;i++)	gamma_ij[i]=new double[levels+1];
Gamma_i=new double[levels+1];
E_i=new double[levels+1];


if(StarkEffect->getPyWrapper() != NULL){
		Py_INCREF(StarkEffect->getPyWrapper());
	}


}

DM_noLaserField::~DM_noLaserField()
{

	for (int i=0;i<levels+1;i++) for (int j=0;j<levels+1;j++)	delete [] k_RungeKutt[i][j]; for (int i=0;i<levels+1;i++)	delete [] k_RungeKutt[i];	delete	[]	k_RungeKutt;
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
	for (int j=1; j<levels + 1;j++)	{
	PopAttr->attArr(i)[j] =  Re(i,j,j);
	PopAttr->attArr(i)[0] += PopAttr->attArr(i)[j];
	}
	
}



void DM_noLaserField::setupEffects(Bunch* bunch){		
	
	AmplAttr = (WaveFunctionAmplitudes*) bunch->getParticleAttributes("Amplitudes");
	PopAttr = (AtomPopulations*) bunch->getParticleAttributes("Populations");
	PopAttr = (AtomPopulations*) bunch->getParticleAttributes("Populations");

	t_part=0;
	
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




Ex_stat/=Ea;	
Ey_stat/=Ea;		
Ez_stat/=Ea;


part_t_step=t_step/gamma/ta;	//time step in frame of particle (in atomic units)


StarkEffect->SetE(Ez_stat);	
	
}














void DM_noLaserField::AmplSolver4step(int i, Bunch* bunch)	{
	
	
	

	tcomplex z,z1,z2,mu_x,mu_y,mu_z;
	double dt,sum;
		

	
	
	
	
			for(int n=1; n<levels+1;n++)	
				StarkEffect->GetEnergyAutoionization(n,E_i[n], Gamma_i[n]);
			

					
			for(int n=2; n<levels+1;n++)
			for(int m=1; m<n;m++)	{
				
				StarkEffect->GetDipoleTransition(n,m,mu_x,mu_y,mu_z);
												
				StarkEffect->GetRelax(n,m,gamma_ij[n][m]);
//				cout<<"delta_res= "<<Ez_las[1]<<"\n";
			}
			

			//cout<<"delta_res= "<<omega_part-(E_i[7]-E_i[1])<<"\n";

					
		

		for(int j=1; j<5; j++)
		for(int n=1; n<levels+1;n++)
		for(int m=1; m<n+1;m++)	{	
			
			if (j==4)	dt=part_t_step;	else dt=part_t_step/2.;
			



			if(n==m){
			for(int k=m+1;k<levels+1;k++)	k_RungeKutt[n][m][j]+=gamma_ij[k][m]*(dm(i,k,k)+k_RungeKutt[k][k][j-1]*dt);
			sum=0; 
			for(int k=1;k<m;k++)			sum+=gamma_ij[m][k];
			k_RungeKutt[n][m][j]-=(dm(i,n,m)+k_RungeKutt[n][m][j-1]*dt)*sum;
			}
			
			if(n!=m){
			sum=0;
			for(int k=1;k<m;k++)			sum+=gamma_ij[m][k];
			for(int k=1;k<n;k++)			sum+=gamma_ij[n][k];
			k_RungeKutt[n][m][j]-=(dm(i,n,m)+k_RungeKutt[n][m][j-1]*dt)*sum/2.;
			}
	
			
			k_RungeKutt[n][m][j]-=(Gamma_i[n]+Gamma_i[m])*(dm(i,n,m)+k_RungeKutt[n][m][j-1]*dt)/2.;

		
			k_RungeKutt[m][n][j]=conj(k_RungeKutt[n][m][j]);

			
			
		}


		
		
		
		
		
	for(int n=1;n<levels+1;n++)
	for(int m=1;m<n+1;m++)	{
		
		
		z=(k_RungeKutt[n][m][1]+2.*k_RungeKutt[n][m][2]+2.*k_RungeKutt[n][m][3]+k_RungeKutt[n][m][4])/6.;	
		Re(i,n,m)+=z.real()*part_t_step;
		Im(i,n,m)+=z.imag()*part_t_step;
		
		Re(i,m,n)=Re(i,n,m);
		Im(i,m,n)=-Im(i,n,m);

		
	}

	
	t_part+=part_t_step;

	
}




















