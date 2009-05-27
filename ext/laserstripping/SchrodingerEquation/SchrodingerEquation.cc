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



#include "SchrodingerEquation.hh"
#include "RungeKuttaTracker.hh"
#include "OrbitConst.hh"
#include "LorentzTransformationEM.hh"




#define Re(i,m) AmplAttr->attArr(i)[m]		//i-part index, n,m-attr index
#define Im(i,m) AmplAttr->attArr(i)[m+levels]

#define x0(i) 	 Coords->attArr(i)[0]
#define y0(i)  	 Coords->attArr(i)[2]
#define z0(i)  	 Coords->attArr(i)[4]
#define px0(i)   Coords->attArr(i)[1]
#define py0(i)   Coords->attArr(i)[3]
#define pz0(i)   Coords->attArr(i)[5]

#define x(i) 	 bunch->coordArr()[i][0]
#define y(i)  	 bunch->coordArr()[i][2]
#define z(i)  	 bunch->coordArr()[i][4]
#define px(i)    bunch->coordArr()[i][1]
#define py(i)    bunch->coordArr()[i][3]
#define pz(i)    bunch->coordArr()[i][5]

#define a(i,m) tcomplex(AmplAttr->attArr(i)[m],AmplAttr->attArr(i)[m+levels])





using namespace LaserStripping;
using namespace OrbitUtils;


SchrodingerEquation::SchrodingerEquation(BaseLaserFieldSource*	BaseLaserField, HydrogenStarkParam* Stark,double par_res)
{
	setName("unnamed");
	


	StarkEffect=Stark;
	LaserField=BaseLaserField;
	Parameter_resonance=par_res;
	levels=	StarkEffect->getStates()*(1+StarkEffect->getStates())*(1+2*StarkEffect->getStates())/6;
	
	zero_cross=true;


	
	
//allocating memory for koefficients of 4-th order Runge-Kutta method and other koeeficients of the master equation
k_RungeKutt=new tcomplex*[levels+1];	for (int i=0;i<levels+1;i++)	k_RungeKutt[i]=new tcomplex[5];
exp_mu_El=new tcomplex**[levels+1];	for (int i=0;i<levels+1;i++)	exp_mu_El[i]=new tcomplex*[levels+1]; for (int i=0;i<levels+1;i++) for (int j=0;j<levels+1;j++)	exp_mu_El[i][j]=new tcomplex[5];
mu_Elas=new tcomplex**[levels+1];	for (int i=0;i<levels+1;i++)	mu_Elas[i]=new tcomplex*[levels+1];		for (int i=0;i<levels+1;i++) for (int j=0;j<levels+1;j++)	mu_Elas[i][j]=new tcomplex[3];
cond=new bool*[levels+1];	for (int i=0;i<levels+1;i++)	cond[i]=new bool[levels+1];
Gamma_i=new double[levels+1];
E_i=new double[levels+1];

if(LaserField->getPyWrapper() != NULL){
		Py_INCREF(LaserField->getPyWrapper());
	}


if(StarkEffect->getPyWrapper() != NULL){
		Py_INCREF(StarkEffect->getPyWrapper());
	}

}














SchrodingerEquation::~SchrodingerEquation()
{

	for (int i=0;i<levels+1;i++) delete	[]	k_RungeKutt[i];		delete [] k_RungeKutt;
	for (int i=0;i<levels+1;i++) for (int j=0;j<levels+1;j++)	delete [] exp_mu_El[i][j]; for (int i=0;i<levels+1;i++)	delete [] exp_mu_El[i];	delete	[]	exp_mu_El;		
	for (int i=0;i<levels+1;i++) for (int j=0;j<levels+1;j++)	delete [] mu_Elas[i][j]; for (int i=0;i<levels+1;i++)	delete [] mu_Elas[i];	delete	[]	mu_Elas;
	delete [] E_i;
	delete [] Gamma_i;
	for (int i=0;i<levels+1;i++)	delete	[]	cond[i];		delete	[]	cond;
	
	delete [] nx;
	delete [] ny;
	delete [] nz;
	

	
	if(LaserField->getPyWrapper() == NULL){
		delete LaserField;
	} else { 
		Py_XDECREF(LaserField->getPyWrapper());
	}
		
	
	if(StarkEffect->getPyWrapper() == NULL){
		delete StarkEffect; 
	} else {
		Py_XDECREF(StarkEffect->getPyWrapper());
	}



}


void SchrodingerEquation::CalcPopulations(int i, Bunch* bunch)	{
		
	PopAttr->attArr(i)[0] = 1;
	for (int j=1; j<levels + 1;j++)	{
	PopAttr->attArr(i)[j] =  Re(i,j)*Re(i,j)+Im(i,j)*Im(i,j);
	PopAttr->attArr(i)[0] -= PopAttr->attArr(i)[j];
	}
	
}




void SchrodingerEquation::setupEffects(Bunch* bunch){	
	
	install_field_dir=true;
	

	if(bunch->hasParticleAttributes("Amplitudes")==0)	{
		std::map<std::string,double> part_attr_dict;
		part_attr_dict["size"] = 2*levels+1;
		bunch->addParticleAttributes("Amplitudes",part_attr_dict);
		
		for (int i=0; i<bunch->getSize();i++)
		bunch->getParticleAttributes("Amplitudes")->attValue(i,1) = 1;
	}

	
		if(bunch->hasParticleAttributes("Populations")==0)	{
			std::map<std::string,double> part_attr_dict;
			part_attr_dict["size"] = levels+1;
			bunch->addParticleAttributes("Populations",part_attr_dict);
		}
	
	
	
	
	
	if (bunch->hasParticleAttributes("pq_coords")==0)	{
	std::map<std::string,double> part_attr_dict;
	part_attr_dict["size"] = 6;
	bunch->addParticleAttributes("pq_coords",part_attr_dict);
	}
	
	Coords = bunch->getParticleAttributes("pq_coords");
	AmplAttr = bunch->getParticleAttributes("Amplitudes");
	PopAttr = bunch->getParticleAttributes("Populations");
	

	nx=new double[bunch->getSize()];
	ny=new double[bunch->getSize()];
	nz=new double[bunch->getSize()];
	
	for (int i=0; i<bunch->getSize();i++)
	CalcPopulations(i, bunch);
	


}
		


void SchrodingerEquation::memorizeInitParams(Bunch* bunch){
	
	for (int i=0; i<bunch->getSize();i++)		{
		x0(i)=bunch->coordArr()[i][0];
		y0(i)=bunch->coordArr()[i][2];
		z0(i)=bunch->coordArr()[i][4];
		
		px0(i)=bunch->coordArr()[i][1];
		py0(i)=bunch->coordArr()[i][3];
		pz0(i)=bunch->coordArr()[i][5];

	}
	
	
}

	
void SchrodingerEquation::finalizeEffects(Bunch* bunch){
	
	if(bunch->hasParticleAttributes("pq_coords")==1)
		bunch->removeParticleAttributes("pq_coords");
	
}







void SchrodingerEquation::applyEffects(Bunch* bunch, int index, 
	                            double* y_in_vct, double* y_out_vct, 
														  double t, double t_step, 
														  BaseFieldSource* fieldSource,
															RungeKuttaTracker* tracker)			{




	
			for (int i=0; i<bunch->getSize();i++)	{


			//	This function gives parameters Ez_stat	Ex_las[1...3]	Ey_las[1...3]	Ez_las[1...3]	
			//in natural unts (Volt per meter)	in the frame of particle

			GetParticleFrameFields(i, t, t_step,bunch,fieldSource);
	
			//	This function gives parameters Ez_stat	Ex_las[1...3]	Ey_las[1...3]	Ez_las[1...3]	t_part	omega_part	part_t_step 
			//in atomic units in frame of particle		
			GetParticleFrameParameters(i,t,t_step,bunch);	
				
			//	This function provides step solution for density matrix using rk4	method
			AmplSolver4step(i,bunch);	
			
			CalcPopulations(i, bunch);
		
		}	

}








void	SchrodingerEquation::RotateElectricFields(double nx,double ny,double nz, tcomplex& Exl,tcomplex& Eyl,tcomplex& Ezl){
	
	double nxy2=nx*nx+ny*ny;
	double n2=nxy2+nz*nz;
	double n=sqrt(n2);



		
		tcomplex Exll=Exl;
		tcomplex Eyll=Eyl;
		tcomplex Ezll=Ezl;
		
		
		if(nxy2>1e-40*n2)	{
			
		Exl=(-nx*nxy2*Ezll + Eyll*nx*ny*(-n + nz) + Exll*(n*ny*ny + nx*nx*nz))/(n*nxy2);
		Eyl=(-ny*nxy2*Ezll + Exll*nx*ny*(-n + nz) + Eyll*(n*nx*nx + ny*ny*nz))/(n*nxy2);
		Ezl=(Exll*nx + Eyll*ny + Ezll*nz)/n;

		}
		

		else	{
			
			Exl=Exl;
			Eyl=Eyl;
			Ezl=Ezl;

		}
	
return ;

}












void SchrodingerEquation::GetParticleFrameFields(int i,double t, double t_step,  Bunch* bunch,  BaseFieldSource* fieldSource)	{
	


		fieldSource->getElectricMagneticField(x0(i),y0(i),z0(i),t,Ex_stat,Ey_stat,Ez_stat,Bx_stat,By_stat,Bz_stat);		
		LorentzTransformationEM::transform(bunch->getMass(),px0(i),py0(i),pz0(i),Ex_stat,Ey_stat,Ez_stat,Bx_stat,By_stat,Bz_stat);	

		if (zero_cross&&install_field_dir)	{
		nx[i] = Ex_stat;	
		ny[i] = Ey_stat;	
		nz[i] = Ez_stat;
		install_field_dir=false;
		}	
		
		
		if (!zero_cross)	{
		nx[i] = Ex_stat;	
		ny[i] = Ey_stat;	
		nz[i] = Ez_stat;
		}
		
		
	for (int j=0; j<3;j++)	{
							
		LaserField->getLaserElectricMagneticField(x0(i)+j*(x(i)-x0(i))/2,y0(i)+j*(y(i)-y0(i))/2,z0(i)+j*(z(i)-z0(i))/2,
				t+j*t_step/2,Ex_las[j],Ey_las[j],Ez_las[j],Bx_las[j],By_las[j],Bz_las[j]);
		

			
	LorentzTransformationEM::complex_transform(bunch->getMass(),
											px0(i),py0(i),pz0(i),
																		 Ex_las[j],Ey_las[j],Ez_las[j],
																		 Bx_las[j],By_las[j],Bz_las[j]);	


	RotateElectricFields(nx[i],ny[i],nz[i],Ex_las[j],Ey_las[j],Ez_las[j]);
	
	}

	
	Ez_stat=sqrt(Ex_stat*Ex_stat+Ey_stat*Ey_stat+Ez_stat*Ez_stat);
	Ex_stat=0;	
	Ey_stat=0;		

	
			
}








void	SchrodingerEquation::GetParticleFrameParameters(int i, double t,double t_step, Bunch* bunch)	{
	
		
	double ta=2.418884326505e-17;			//atomic unit of time
	double Ea=5.14220642e011;				//Atomic unit of electric field
	double m=bunch->getMass();

//This line calculates relyativistic factor-Gamma
double gamma=sqrt(m*m+px0(i)*px0(i)+py0(i)*py0(i)+pz0(i)*pz0(i))/m;
double coeff=1./(gamma*ta);

part_t_step=t_step*coeff;	//time step in frame of particle (in atomic units)
t_part=t*coeff;	



for(int j=0;j<3;j++)	
{		
	Ex_las[j]/=Ea; 
	Ey_las[j]/=Ea; 
	Ez_las[j]/=Ea; }

Ex_stat/=Ea;	
Ey_stat/=Ea;		
Ez_stat/=Ea;



omega_part=gamma*ta*LaserField->getFrequencyOmega(m,x0(i),y0(i),z0(i),px0(i),py0(i),pz0(i),t);		// frequensy of laser in particle frame (in atomic units)

StarkEffect->SetE(Ez_stat);	


	
}














void SchrodingerEquation::AmplSolver4step(int i, Bunch* bunch)	{
	

	

	tcomplex z,z1,z2,mu_x,mu_y,mu_z;
	double dt;
		

	
	
	
	
			for(int n=1; n<levels+1;n++)	
				StarkEffect->GetEnergyAutoionization(n,E_i[n], Gamma_i[n]);
			

					
			for(int n=2; n<levels+1;n++)
			for(int m=1; m<n;m++)	{
				
				StarkEffect->GetDipoleTransition(n,m,mu_x,mu_y,mu_z);
				
				for(int j=0; j<3;j++)
					mu_Elas[n][m][j]=(mu_x*conj(Ex_las[j])+mu_y*conj(Ey_las[j])+mu_z*conj(Ez_las[j]))*(J/2.);	
				
				cond[n][m]=fabs(fabs(E_i[n]-E_i[m])-omega_part)<Parameter_resonance*abs(2.*mu_Elas[n][m][1]);	
								

			}

//			cout<<"delta_res= "<<omega_part-(E_i[7]-E_i[1])<<"\n";



		
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
		for(int m=1; m<levels+1; m++)	{

			
			if (j==4)	dt=part_t_step;	else dt=part_t_step/2.;

			k_RungeKutt[m][j]=tcomplex(0.);
			for(int k=1;k<m;k++)			if (cond[m][k])		 k_RungeKutt[m][j]-=exp_mu_El[m][k][j]*(a(i,k)+k_RungeKutt[k][j-1]*dt);	
			for(int k=m+1;k<levels+1;k++)	if (cond[k][m]) 	 k_RungeKutt[m][j]+=conj(exp_mu_El[k][m][j])*(a(i,k)+k_RungeKutt[k][j-1]*dt); 
			
			
			k_RungeKutt[m][j]-=Gamma_i[m]*(a(i,m)+k_RungeKutt[m][j-1]*dt)/2.;
			

	
		}


		
		
		
		
		
	for(int m=1;m<levels+1;m++)	{
				
		z=(k_RungeKutt[m][1]+2.*k_RungeKutt[m][2]+2.*k_RungeKutt[m][3]+k_RungeKutt[m][4])/6.;	
		Re(i,m)+=z.real()*part_t_step;
		Im(i,m)+=z.imag()*part_t_step;

		
	}
	
	



	
	

}




















