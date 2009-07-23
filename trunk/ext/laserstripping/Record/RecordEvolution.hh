#ifndef RECORDEVOLUTION_HH_
#define RECORDEVOLUTION_HH_


#include "Python.h"

#include "ExternalEffects.hh"
#include "BaseLaserFieldSource.hh"
#include "ParticleAttributes.hh"
#include <string>



using namespace TrackerRK4;

namespace LaserStripping{
	
	class  RecordEvolution: public ExternalEffects
	{
		public:
		
			/** Constructor. */
			RecordEvolution(std::string effect,int ind_effect, int num);
			
			/** Destructor. */
			~RecordEvolution();
			
		
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

		  		  
			int num_plot;
			int index_effect;
			std::string effect_name;
			double t_per;
			double t_step;
			double t_in;
			int Num;
			bool setup_par;
			
			
			ParticleAttributes* RecEff;
			ParticleAttributes* Coords;
			ParticleAttributes* Evol;
			   
				
	};
};




#endif /*RECORDEVOLUTION_HH_*/
