#ifndef SCHRODINGEREQUATION_HH_
#define SCHRODINGEREQUATION_HH_



#include "Python.h"

#include "ExternalEffects.hh"
#include "BaseLaserFieldSource.hh"
#include "HydrogenStarkParam.hh"
#include "ParticleAttributes.hh"


using namespace TrackerRK4;

namespace LaserStripping{
	
	class  SchrodingerEquation: public ExternalEffects
	{
		public:
		
			/** Constructor */
			SchrodingerEquation(OrbitUtils::BaseLaserFieldSource*	BaseLaserField, HydrogenStarkParam* Stark,double par_res);

			
			/** Destructor. */
			~SchrodingerEquation();
			

		
		/** It initializes effects. */
		void setupEffects(Bunch* bunch);
		
		/*it memorizes initial coordinates and impulses before rk step*/
		void memorizeInitParams(Bunch* bunch);
		
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
			  tcomplex** k_RungeKutt;
			  tcomplex*** mu_Elas;
			  double* Gamma_i;
			  double* E_i;
			  bool** cond;
			  

			  
			  int levels;
			  double Parameter_resonance;
			  
			  ParticleAttributes* AmplAttr;
			  ParticleAttributes* PopAttr;
			  ParticleAttributes* Coords;
			  
		  

			  			 
			  //time and frequensy of laser in frame of particle
			  double omega_part;
			  double part_t_step;
			  double t_part;

			  					  
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
			  
			  double* nx;
			  double* ny;
			  double* nz;
			  
			  bool zero_cross;
			  bool install_field_dir;
			  

			  	
				/**Solver for Amplitudes**/
				void AmplSolver4step(int i,Bunch* bunch);
				
				/**Calculates populations**/
				void CalcPopulations(int i,Bunch* bunch);
							
				/*this all parameters in frame of particle in atomic units necesarry for applying external effect  */
				void GetParticleFrameParameters(int i, double t, double t_step, Bunch* bunch);
				
				/*this method gives laser and static fields transformed by rotation relatively z axes  */
				void GetParticleFrameFields(int i, double t,double t_step,  Bunch* bunch,  OrbitUtils::BaseFieldSource* fieldSource);
				
				/*This method provides rotational transformation of statis and laser field in frame of particle to z axes */
				void	RotateElectricFields(double nx, double ny, double nz, tcomplex& Ex_l,tcomplex& Ey_l,tcomplex& Ez_l);
				

				
	};
};



#endif /*SCHRODINGEREQUATION_HH_*/


