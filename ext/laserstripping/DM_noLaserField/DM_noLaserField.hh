#ifndef DM_NOLASERFIELD_HH_
#define DM_NOLASERFIELD_HH_



#include "Python.h"

#include "ExternalEffects.hh"
#include "BaseLaserFieldSource.hh"
#include "HydrogenStarkParam.hh"
#include "WaveFunctionAmplitudes.hh"
#include "AtomPopulations.hh"
#include "pq_coordinates.hh"


using namespace TrackerRK4;

namespace LaserStripping{
	
	class  DM_noLaserField: public ExternalEffects
	{
		public:
		
			/** Constructor */
			DM_noLaserField(HydrogenStarkParam* Stark);

			
			/** Destructor. */
			~DM_noLaserField();
			
		
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

		  
			  

			  HydrogenStarkParam* StarkEffect;
			  
			  //this array is used on each step of solution of density matrix equation at definite field  

			  tcomplex*** k_RungeKutt;

			  double* Gamma_i;
			  double* E_i;
			  double** gamma_ij;


			  
			  int levels;

			  
			  WaveFunctionAmplitudes* AmplAttr;
			  AtomPopulations* PopAttr;
			  pq_coordinates* Coords;
			  

 			 
			  //time and frequensy of laser in frame of particle

			  double part_t_step;
			  double t_part;

			  					  
			  //Fields in the particle frame  
		  
			  double Ex_stat;
			  double Ey_stat;
			  double Ez_stat;
			  
			  double Bx_stat;
			  double By_stat;
			  double Bz_stat;
			  


			  	
				/**Solver for Amplitudes**/
				void AmplSolver4step(int i,Bunch* bunch);
				
				/**Calculates populations**/
				void CalcPopulations(int i,Bunch* bunch);
							
				/*this all parameters in frame of particle in atomic units necesarry for applying external effect  */
				void GetParticleFrameParameters(int i, double t, double t_step, Bunch* bunch);
				
				/*this method gives laser and static fields transformed by rotation relatively z axes  */
				void GetParticleFrameFields(int i, double t,double t_step,   Bunch* bunch,  OrbitUtils::BaseFieldSource* fieldSource);


				
	};
};



#endif /*DM_NOLASERFIELD_HH_*/


