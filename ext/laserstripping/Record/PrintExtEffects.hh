#ifndef PRINTEXTEFFECTS_HH_
#define PRINTEXTEFFECTS_HH_



#include "Python.h"

#include "ExternalEffects.hh"
#include "BaseLaserFieldSource.hh"


using namespace TrackerRK4;

namespace LaserStripping{
	
	class  PrintExtEffects: public ExternalEffects
	{
		public:
		
			/** Constructor. */
			PrintExtEffects(char* eff_name,int i,char* addr_print);
			
			/** Destructor. */
			~PrintExtEffects();
			
		
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

		  

			  int Num;
			  int num_print;
			  char* addr_print;
			  char* eff_name;
			  char addr_name[1024];
			  
			  int rank_MPI,size_MPI;
			  bool setup_par;
			  double t_in;
			  
			  			 


								

				
	};
};



#endif /*PRINTEXTEFFECTS_HH_*/


