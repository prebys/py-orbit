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






#include "TwoLevelAtom.hh"
#include "RungeKuttaTracker.hh"
#include "OrbitConst.hh"
#include "LorentzTransformationEM.hh"






#define Re1(i) AmplAttr->attArr(i)[1]		//i-part index, n,m-attr index
#define Im1(i) AmplAttr->attArr(i)[2]
#define Re2(i) AmplAttr->attArr(i)[3]		//i-part index, n,m-attr index
#define Im2(i) AmplAttr->attArr(i)[4]

#define time(i)  AmplAttr->attArr(i)[5]
#define x0(i) 	 AmplAttr->attArr(i)[6]
#define y0(i)  	 AmplAttr->attArr(i)[7]
#define z0(i)  	 AmplAttr->attArr(i)[8]
#define px0(i)   AmplAttr->attArr(i)[9]
#define py0(i)   AmplAttr->attArr(i)[10]
#define pz0(i)   AmplAttr->attArr(i)[11]
#define Ampl_1(i) tcomplex(AmplAttr->attArr(i)[1],AmplAttr->attArr(i)[2])
#define Ampl_2(i) tcomplex(AmplAttr->attArr(i)[3],AmplAttr->attArr(i)[4])

//#define k_rk(j,n,m) k_RungeKutt[j][(n-1)*levels+m]




using namespace LaserStripping;
using namespace OrbitUtils;


TwoLevelAtom::TwoLevelAtom(BaseLaserFieldSource*	BaseLaserField, double delta_E, double dipole_tr)
{

	setName("unnamed");
	
	LaserField=BaseLaserField;
	d_Energy=delta_E;
	dip_transition=dipole_tr;
	
	if(LaserField->getPyWrapper() != NULL){
			Py_INCREF(LaserField->getPyWrapper());
		}

}




TwoLevelAtom::~TwoLevelAtom() 
{
	
			if(LaserField->getPyWrapper() == NULL){
				delete LaserField;
			} else {
				Py_XDECREF(LaserField->getPyWrapper());
			}

}













void TwoLevelAtom::CalcPopulations(int i, Bunch* bunch)	{
		
	PopAttr->attArr(i)[1] =  Re1(i)*Re1(i)+Im1(i)*Im1(i);
	PopAttr->attArr(i)[2] =  Re2(i)*Re2(i)+Im2(i)*Im2(i);
	PopAttr->attArr(i)[0] = PopAttr->attArr(i)[1]+PopAttr->attArr(i)[2];
	
}




void TwoLevelAtom::setupEffects(Bunch* bunch){		

	AmplAttr = (WaveFunctionAmplitudes*) bunch->getParticleAttributes("Amplitudes");
	PopAttr = (AtomPopulations*) bunch->getParticleAttributes("Populations");

	
	for (int i=0; i<bunch->getSize();i++)		{
		x0(i)=bunch->coordArr()[i][0];
		y0(i)=bunch->coordArr()[i][2];
		z0(i)=bunch->coordArr()[i][4];
		
		px0(i)=bunch->coordArr()[i][1];
		py0(i)=bunch->coordArr()[i][3];
		pz0(i)=bunch->coordArr()[i][5];
		
		time(i)=0.;
		//All other attributes of particle are initiated in Python script

		CalcPopulations(i, bunch);
	}
	
}
		

	
void TwoLevelAtom::finalizeEffects(Bunch* bunch) {}







void TwoLevelAtom::applyEffects(Bunch* bunch, int index, 
	                            double* y_in_vct, double* y_out_vct, 
														  double t, double t_step, 
														  BaseFieldSource* fieldSource,
															RungeKuttaTracker* tracker)			{




	
		for (int i=0; i<bunch->getSize();i++)	{

			//	This function gives parameters Ez_stat	Ex_las[1...3]	Ey_las[1...3]	Ez_las[1...3]	
			//in natural unts (Volt per meter)	in the frame of particle				
			GetParticleFrameFields(i, t,t_step, bunch);
	
			//	This function gives parameters Ez_stat	Ex_las[1...3]	Ey_las[1...3]	Ez_las[1...3]	t_part	omega_part	part_t_step 
			//in atomic units in frame of particle		
			GetParticleFrameParameters(i,t,t_step,bunch);	
				
			
			//	This function provides step solution for density matrix using rk4	method
			AmplSolver4step(i,bunch);	
			
			CalcPopulations(i, bunch);
			

		}	
//	cout<<scientific<<setprecision(20)<<bunch->x(0)<<"\t"<<bunch->y(0)<<"\t"<<bunch->z(0)<<"\n";


	
}















void TwoLevelAtom::GetParticleFrameFields(int i,double t,double t_step,  Bunch* bunch)	{
	
	double** xyz = bunch->coordArr();
	double E_x,E_y,E_z,nx,ny,nz,Eabs;

	
	for (int j=0; j<3;j++)	{
							
		LaserField->getLaserElectricMagneticField(x0(i)+j*(xyz[i][0]-x0(i))/2,y0(i)+j*(xyz[i][2]-y0(i))/2,z0(i)+j*(xyz[i][4]-z0(i))/2,
				t+j*t_step/2,Ex_las[j],Ey_las[j],Ez_las[j],Bx_las[j],By_las[j],Bz_las[j]);

	LorentzTransformationEM::complex_electric_transform(bunch->getMass(),
											px0(i),py0(i),pz0(i),
																		 Ex_las[j],Ey_las[j],Ez_las[j],
																		 Bx_las[j],By_las[j],Bz_las[j]);	

	if(j==0)		{
		E_x=abs(Ex_las[j]);
		E_y=abs(Ey_las[j]);
		E_z=abs(Ez_las[j]);
		Eabs=sqrt(pow(E_x,2)+pow(E_y,2)+pow(E_z,2));
	}

	
	Ez_las[j]=(Ex_las[j]*E_x+Ey_las[j]*E_y+Ez_las[j]*E_z)/Eabs;

	
	}
			
}








void	TwoLevelAtom::GetParticleFrameParameters(int i, double t,double t_step, Bunch* bunch)	{
	
		
	double ta=2.418884326505e-17;			//atomic unit of time
	double Ea=5.14220642e+011;				//Atomic unit of electric field
	double m=bunch->getMass();

//This line calculates relyativistic factor-Gamma
double gamma=sqrt(m*m+px0(i)*px0(i)+py0(i)*py0(i)+pz0(i)*pz0(i))/m;



//Convertion form SI units to atomic units
for(int j=0;j<3;j++)		
	Ez_las[j]/=Ea; 


part_t_step=t_step/gamma/ta;	//time step in frame of particle (in atomic units)
t_part=time(i);					//time  in frame of particle (in atomic units) 




	
}














void TwoLevelAtom::AmplSolver4step(int i, Bunch* bunch)	{
	
	
	

	tcomplex z1,z2;
	double dt;
		

	
	
					
	for(int j=0; j<3;j++)	mu_Elas[j]=dip_transition*conj(Ez_las[j])*J/2.;	
	

			
			z1=exp(J*t_part*fabs(d_Energy));		
			z2=exp(J*part_t_step*fabs(d_Energy)/2.);

			
			exp_mu_El[1]=z1*mu_Elas[0];			
			exp_mu_El[2]=z1*z2*mu_Elas[1];
			exp_mu_El[3]=exp_mu_El[2];
			exp_mu_El[4]=z1*z2*z2*mu_Elas[2];
			
						

			for(int j=1; j<5; j++)	{	
				
				if (j==4)	dt=part_t_step;	else dt=part_t_step/2.;

				k_RungeKutt_1[j]=conj(exp_mu_El[j])*(Ampl_2(i)+k_RungeKutt_2[j-1]*dt);	
				k_RungeKutt_2[j]=-exp_mu_El[j]*(Ampl_1(i)+k_RungeKutt_1[j-1]*dt);	
			
			}
					
			
				
		z1=(k_RungeKutt_1[1]+2.*k_RungeKutt_1[2]+2.*k_RungeKutt_1[3]+k_RungeKutt_1[4])/6.;	
		z2=(k_RungeKutt_2[1]+2.*k_RungeKutt_2[2]+2.*k_RungeKutt_2[3]+k_RungeKutt_2[4])/6.;
		
		Re1(i)+=z1.real()*part_t_step;
		Im1(i)+=z1.imag()*part_t_step;
		
		Re2(i)+=z2.real()*part_t_step;
		Im2(i)+=z2.imag()*part_t_step;

		
			
			
	
	
	time(i)+=part_t_step;
	
	x0(i)=bunch->coordArr()[i][0];
	y0(i)=bunch->coordArr()[i][2];
	z0(i)=bunch->coordArr()[i][4];
	
	px0(i)=bunch->coordArr()[i][1];
	py0(i)=bunch->coordArr()[i][3];
	pz0(i)=bunch->coordArr()[i][5];
	
	
	
	

}




















