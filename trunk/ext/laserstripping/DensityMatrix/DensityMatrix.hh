#ifndef DENSITYMATRIX_HH_
#define DENSITYMATRIX_HH_



#include "Python.h"

#include "ExternalEffects.hh"
#include "BaseLaserFieldSource.hh"
#include "HydrogenStarkParam.hh"


using namespace TrackerRK4;

namespace LaserStripping{
	
	class  DensityMatrix: public ExternalEffects
	{
		public:
		
			/** Constructor */
			DensityMatrix(OrbitUtils::BaseLaserFieldSource*	BaseLaserField, HydrogenStarkParam* Stark,double par_res);

			
			/** Destructor. */
			~DensityMatrix();
			
		/** Method that initialise and defines parameters of printing */
		void SetupPrint(int i,char* addr_print);
		
		/** It initializes effects. */
		void setupEffects(Bunch* bunch);
		
		/** It finalizes effects. */
		void finalizeEffects(Bunch* bunch);

		/** It applies the external effects to a particle with certain index. */
		void applyEffects(Bunch* bunch, int index, 
	                            double* y_in_vct, double* y_out_vct, 
														  double t, double t_step, 
														  OrbitUtils::BaseFieldSource* fieldSource,
															RungeKuttaTracker* tracker);	
		

		

		  private:

		  
			  
			  OrbitUtils::BaseLaserFieldSource*	LaserField;
			  HydrogenStarkParam* StarkEffect;
			  
			  //this array is used on each step of solution of density matrix equation at definite field  
			  tcomplex*** exp_mu_El;
			  tcomplex*** k_RungeKutt;
			  tcomplex*** mu_Elas;
			  double* Gamma_i;
			  double* E_i;
			  double** gamma_ij;
			  bool** cond;

			  
			  int levels;
			  double Parameter_resonance;
			  
		  

			  			 
			  //time and frequensy of laser in frame of particle
			  double omega_part;
			  double part_t_step;
			  double t_part;
			  double phasa_part;
			  					  
			  //Fields in the particle frame  
			  tcomplex Ex_las[3];
			  tcomplex Ey_las[3];
			  tcomplex Ez_las[3];
			  
			  tcomplex Bx_las[3];
			  tcomplex By_las[3];
			  tcomplex Bz_las[3];
			  
			  double Ex_stat;
			  double Ey_stat;
			  double Ez_stat;
			  
			  double Bx_stat;
			  double By_stat;
			  double Bz_stat;
			  
			  int print_par;
			  int max_print_par;
			  char* addr_print;
			  			 

			  	
				/**Solver for Amplitudes**/
				void AmplSolver4step(int i,Bunch* bunch);
							
				/*this all parameters in frame of particle in atomic units necesarry for applying external effect  */
				void GetParticleFrameParameters(int i, double t, double t_step, Bunch* bunch);
				
				/*this method gives laser and static fields transformed by rotation relatively z axes  */
				void GetParticleFrameFields(int i, double t,double t_step,   Bunch* bunch,  OrbitUtils::BaseFieldSource* fieldSource);
				
				/*This method provides rotational transformation of statis and laser field in frame of particle to z axes */
				double	RotateElectricFields(double Ex_s, double Ey_s, double Ez_s,tcomplex& Ex_l,tcomplex& Ey_l,tcomplex& Ez_l);
				

				
	};
};



#endif /*DENSITYMATRIX_HH_*/


