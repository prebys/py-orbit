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
#include "MathPolinomial.hh"


#include "LasStripExternalEffects.hh"
#include "RungeKuttaTracker.hh"
#include "OrbitConst.hh"
#include "LorentzTransformationEM.hh"
#include "HydrogenStarkParam.hh"





#define Re(i,n,m) bunch->arrAttr[i][(n-1)*levels+m]		//i-part index, n,m-attr index
#define Im(i,n,m) bunch->arrAttr[i][(n-1)*levels+m+levels*levels]
#define time(i)  bunch->arrAttr[i][2*levels*levels+1]
#define x0(i)  bunch->arrAttr[i][2*levels*levels+2]
#define y0(i)  bunch->arrAttr[i][2*levels*levels+3]
#define z0(i)  bunch->arrAttr[i][2*levels*levels+4]
#define px0(i)  bunch->arrAttr[i][2*levels*levels+5]
#define py0(i)  bunch->arrAttr[i][2*levels*levels+6]
#define pz0(i)  bunch->arrAttr[i][2*levels*levels+7]
#define dm(i,n,m) tcomplex(bunch->arrAttr[i][(n-1)*levels+m],bunch->arrAttr[i][(n-1)*levels+m+levels*levels])
#define k_rk(j,n,m) k_RungeKutt[j][(n-1)*levels+m]




using namespace LaserStripping;
using namespace OrbitUtils;


LasStripExternalEffects::LasStripExternalEffects(BaseLaserFieldSource*	BaseLaserField,char* addressEG,int states,double par_res)
{
	setName("unnamed");
	
//	for(int i=0;i<11;i++)
//	for(int j=0;j<11;j++)
//	cout<<MathPolinomial::Factorial(i)<<"\n";

	LaserField=BaseLaserField;
	Parameter_resonance=par_res;
	levels=states*(1+states)*(1+2*states)/6;

	HydrogenStarkParam::ReadData(addressEG,states);

	
	
	
//allocating memory for koefficients of 4-th order Runge-Kutta method and other koeeficients of the master equation
k_RungeKutt=new tcomplex**[levels+1];	for (int i=0;i<levels+1;i++)	k_RungeKutt[i]=new tcomplex*[levels+1]; for (int i=0;i<levels+1;i++) for (int j=0;j<levels+1;j++)	k_RungeKutt[i][j]=new tcomplex[5];
exp_mu_El=new tcomplex**[levels+1];	for (int i=0;i<levels+1;i++)	exp_mu_El[i]=new tcomplex*[levels+1]; for (int i=0;i<levels+1;i++) for (int j=0;j<levels+1;j++)	exp_mu_El[i][j]=new tcomplex[5];
mu_Elas=new tcomplex**[levels+1];	for (int i=0;i<levels+1;i++)	mu_Elas[i]=new tcomplex*[levels+1];		for (int i=0;i<levels+1;i++) for (int j=0;j<levels+1;j++)	mu_Elas[i][j]=new tcomplex[3];
gamma_ij=new double*[levels+1];	for (int i=0;i<levels+1;i++)	gamma_ij[i]=new double[levels+1];
cond=new bool*[levels+1];	for (int i=0;i<levels+1;i++)	cond[i]=new bool[levels+1];
Gamma_i=new double[levels+1];
E_i=new double[levels+1];


}


LasStripExternalEffects::~LasStripExternalEffects()
{

	for (int i=0;i<levels+1;i++) for (int j=0;j<levels+1;j++)	delete [] k_RungeKutt[i][j]; for (int i=0;i<levels+1;i++)	delete [] k_RungeKutt[i];	delete	[]	k_RungeKutt;
	for (int i=0;i<levels+1;i++) for (int j=0;j<levels+1;j++)	delete [] exp_mu_El[i][j]; for (int i=0;i<levels+1;i++)	delete [] exp_mu_El[i];	delete	[]	exp_mu_El;		
	for (int i=0;i<levels+1;i++) for (int j=0;j<levels+1;j++)	delete [] mu_Elas[i][j]; for (int i=0;i<levels+1;i++)	delete [] mu_Elas[i];	delete	[]	mu_Elas;
	delete [] E_i;
	delete [] Gamma_i;
	for (int i=0;i<levels+1;i++)	delete	[]	gamma_ij[i];	delete	[]	gamma_ij;
	for (int i=0;i<levels+1;i++)	delete	[]	cond[i];		delete	[]	cond;

}




void LasStripExternalEffects::setupEffects(Bunch* bunch){
	
	for (int i=0; i<bunch->getSizeGlobal();i++)		{
		x0(i)=bunch->coordArr()[i][0];
		y0(i)=bunch->coordArr()[i][2];
		z0(i)=bunch->coordArr()[i][4];
		
		px0(i)=bunch->coordArr()[i][1];
		py0(i)=bunch->coordArr()[i][3];
		pz0(i)=bunch->coordArr()[i][5];
		
		time(i)=0.;
		//All other attributes of particle are initiated in Python script

	}
	
}
		

	
void LasStripExternalEffects::finalizeEffects(Bunch* bunch){
	

}



void LasStripExternalEffects::applyEffects(Bunch* bunch, int index, 
	                            double* y_in_vct, double* y_out_vct, 
														  double t, double t_step, 
														  BaseFieldSource* fieldSource,
															RungeKuttaTracker* tracker)			{




	
		for (int i=0; i<bunch->getSizeGlobal();i++)	{

			
			//	This function gives parameters Ez_stat	Ex_las[1...3]	Ey_las[1...3]	Ez_las[1...3]	
			//in natural unts (Volt per meter)	in the frame of particle				
			GetFrameParticleFields(i, t, bunch,fieldSource);
	
			//	This function gives parameters Ez_stat	Ex_las[1...3]	Ey_las[1...3]	Ez_las[1...3]	t_part	omega_part	part_t_step 
			//in atomic units in frame of particle		
			GetFrameParticleParameters(i,t,t_step,bunch);	
	
			
			
			ofstream file("/home/tg4/workspace/PyOrbit/ext/laserstripping/working_dir/data_ampl.txt",ios::app);
			file<<t<<"\t";
			for(int n=1;n<levels+1;n++)	file<<Re(i,n,n)<<"\t";
			double sum=0;	for(int n=1;n<levels+1;n++)	sum+=Re(i,n,n);
			file<<sum<<"\n";
			file.close();			
		
				
			//	This function provides step solution for density matrix using rk4	method
			AmplSolver4step(i,bunch);	
		
	
		
		}	

//	cout<<scientific<<setprecision(20)<<bunch->x(0)<<"\t"<<bunch->y(0)<<"\t"<<bunch->z(0)<<"\n";

	
}









double	LasStripExternalEffects::RotateElectricFields(double Exs, double Eys, double Ezs,tcomplex& Exl,tcomplex& Eyl,tcomplex& Ezl){
	
	double Exy2=Exs*Exs+Eys*Eys;
	double E2=Exs*Exs+Eys*Eys+Ezs*Ezs;
	double E=sqrt(E2);

	
	if(Exy2>1e-40*E2)	{
		
		tcomplex Exll=Exl;
		tcomplex Eyll=Eyl;
		tcomplex Ezll=Ezl;
		
		Exl=(-Exs*Exy2*Ezll + Eyll*Exs*Eys*(-E + Ezs) + Exll*(E*Eys*Eys + Exs*Exs*Ezs))/(E*Exy2);
		Eyl=(-Eys*Exy2*Ezll + Exll*Exs*Eys*(-E + Ezs) + Eyll*(E*Exs*Exs + Eys*Eys*Ezs))/(E*Exy2);
		Ezl=(Exll*Exs + Eyll*Eys + Ezll*Ezs)/E;
		
	}
	
	else	{
		
		Exl=Exl;
		Eyl=Eyl;
		Ezl=Ezl;

	}
	
return E;

}












void LasStripExternalEffects::GetFrameParticleFields(int i,double t,  Bunch* bunch,  BaseFieldSource* fieldSource)	{
	
	double** xyz = bunch->coordArr();
	double Ez;
	
		fieldSource->getElectricMagneticField(x0(i),y0(i),z0(i),t,Ex_stat,Ey_stat,Ez_stat,Bx_stat,By_stat,Bz_stat);		
		LorentzTransformationEM::transform(bunch->getMass(),px0(i),py0(i),pz0(i),Ex_stat,Ey_stat,Ez_stat,Bx_stat,By_stat,Bz_stat);
					
	for (int j=0; j<3;j++)	{
							
		LaserField->getLaserElectricMagneticField(x0(i)+j*(xyz[i][0]-x0(i))/2,y0(i)+j*(xyz[i][2]-y0(i))/2,z0(i)+j*(xyz[i][4]-z0(i))/2,
				t,Ex_las[j],Ey_las[j],Ez_las[j],Bx_las[j],By_las[j],Bz_las[j]);
			
	LorentzTransformationEM::complex_transform(bunch->getMass(),
											px0(i),py0(i),pz0(i),
																		 Ex_las[j],Ey_las[j],Ez_las[j],
																		 Bx_las[j],By_las[j],Bz_las[j]);	
	
	
	Ez=RotateElectricFields(Ex_stat,Ey_stat,Ez_stat,Ex_las[j],Ey_las[j],Ez_las[j]);
	
	}
	
	Ex_stat=0;	
	Ey_stat=0;		
	Ez_stat=Ez;
	
			
	
}








void	LasStripExternalEffects::GetFrameParticleParameters(int i, double t,double t_step, Bunch* bunch)	{
	
		
	double ta=2.418884326505e-17;			//atomic unit of time
	double Ea=5.14220642e11;				//Atomic unit of electric field

//This line calculates relyativistic factor-Gamma
double gamma=sqrt(bunch->getMass()*bunch->getMass()+px0(i)*px0(i)+py0(i)*py0(i)+pz0(i)*pz0(i))/bunch->getMass();




for(int j=0;j<3;j++)	
{		
	Ex_las[j]/=Ea; 
	Ey_las[j]/=Ea; 
	Ez_las[j]/=Ea; }

Ex_stat/=Ea;	
Ey_stat/=Ea;		
Ez_stat/=Ea;


part_t_step=t_step/gamma/ta;	//time step in frame of particle (in atomic units)
t_part=time(i);					//time  in frame of particle (in atomic units). 
omega_part=ta*LaserField->getFrequencyOmega(x0(i),y0(i),z0(i),px0(i),py0(i),pz0(i),t);		// frequensy of laser in particle frame (in atomic units)

HydrogenStarkParam::SetE(Ez_stat);	
	
}














void LasStripExternalEffects::AmplSolver4step(int i, Bunch* bunch)	{
	
	
	

	tcomplex z,z1,z2,mu_x,mu_y,mu_z;
	double dt;
		

	
	
	
	
			for(int n=1; n<levels+1;n++)	
				HydrogenStarkParam::GetEnergyAutoionization(n,E_i[n], Gamma_i[n]);
			

					
			for(int n=2; n<levels+1;n++)
			for(int m=1; m<n;m++)	{
				
				HydrogenStarkParam::GetDipoleTransition(n,m,mu_x,mu_y,mu_z);
				
				for(int j=0; j<3;j++)
					mu_Elas[n][m][j]=mu_x*conj(Ex_las[j])+mu_y*conj(Ey_las[j])+mu_z*conj(Ez_las[j]);	
				
				cond[n][m]=fabs(fabs(E_i[n]-E_i[m])-omega_part)<Parameter_resonance*abs(mu_Elas[n][m][1]);	///THIS CRITERII SHOULD BE CHANGED
				HydrogenStarkParam::GetRelax(n,m,gamma_ij[n][m]);
			
			}
			
				




		
		for(int n=2; n<levels+1;n++)
		for(int m=1; m<n;m++)		
			if (cond[n][m])		{
				z1=exp(J*t_part*fabs(E_i[n]-E_i[m]));		
				z2=exp(J*part_t_step*fabs(E_i[n]-E_i[m])/2.);
			
			exp_mu_El[n][m][1]=z1*mu_Elas[n][m][0];			
			exp_mu_El[n][m][2]=z1*z2*mu_Elas[n][m][1];
			exp_mu_El[n][m][3]=exp_mu_El[n][m][2];
			exp_mu_El[n][m][4]=z1*z2*z2*mu_Elas[n][m][2];	
			
		}


					
		

		for(int j=1; j<5; j++)
		for(int n=1; n<levels+1;n++)
		for(int m=1; m<n+1;m++)	{	
			
			if (j==4)	dt=part_t_step;	else dt=part_t_step/2.;

			k_RungeKutt[n][m][j]*=0.;
			for(int k=1;k<n;k++)			if (cond[n][k]) 	k_RungeKutt[n][m][j]+=exp_mu_El[n][k][j]*(dm(i,k,m)+k_RungeKutt[k][m][j-1]*dt);	
			for(int k=n+1;k<levels+1;k++)	if (cond[k][n]) 	k_RungeKutt[n][m][j]+=conj(exp_mu_El[k][n][j])*(dm(i,k,m)+k_RungeKutt[k][m][j-1]*dt);
			for(int k=1;k<m;k++)			if (cond[m][k]) 	k_RungeKutt[n][m][j]-=conj(exp_mu_El[m][k][j])*(dm(i,n,k)+k_RungeKutt[n][k][j-1]*dt);
			for(int k=m+1;k<levels+1;k++)	if (cond[k][m]) 	k_RungeKutt[n][m][j]-=exp_mu_El[k][m][j]*(dm(i,n,k)+k_RungeKutt[n][k][j-1]*dt);
			k_RungeKutt[n][m][j]*=J/2.;
			



			if(n==m){
			for(int k=m+1;k<levels+1;k++)	k_RungeKutt[n][m][j]+=gamma_ij[k][m]*(dm(i,k,k)+k_RungeKutt[k][k][j-1]*dt);
			for(int k=1;k<m;k++)			k_RungeKutt[n][m][j]-=gamma_ij[m][k]*(dm(i,n,m)+k_RungeKutt[n][m][j-1]*dt);
			}
			
			if(n!=m){
			for(int k=1;k<m;k++)			k_RungeKutt[n][m][j]-=gamma_ij[m][k]*(dm(i,n,m)+k_RungeKutt[n][m][j-1]*dt)/2.;
			for(int k=1;k<n;k++)			k_RungeKutt[n][m][j]-=gamma_ij[n][k]*(dm(i,n,m)+k_RungeKutt[n][m][j-1]*dt)/2.;
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

	time(i)+=part_t_step;
	
	x0(i)=bunch->coordArr()[i][0];
	y0(i)=bunch->coordArr()[i][2];
	z0(i)=bunch->coordArr()[i][4];
	
	px0(i)=bunch->coordArr()[i][1];
	py0(i)=bunch->coordArr()[i][3];
	pz0(i)=bunch->coordArr()[i][5];

	

}




