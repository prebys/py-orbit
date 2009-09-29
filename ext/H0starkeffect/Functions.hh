#ifndef FUNCTIONS_HH_
#define FUNCTIONS_HH_


#include "tcomplex.hh"
#include "CppPyWrapper.hh"
#include <string>






		
		class  Functions: public OrbitUtils::CppPyWrapper
		{

			
		public:
		/*constructor*/	
		Functions(int n1,int n2, int m, int n_stepss, double stepp, int nsumm);
		
		/*destructor*/	
		~Functions();
		
		void setupE(double E);
			


		double getMfunctionModulus(double E, double Gamma, double reZ, double imZ);


		


		private:  
			
			int n1;
			int n2;
			int m;
			
			int n_steps;
			double step;
			int nsum;
			
			double E;
			


		
		};




#endif /*FUNCTIONS_HH_*/
