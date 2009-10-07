#ifndef FUNCTIONS_HH_
#define FUNCTIONS_HH_


#include "tcomplex.hh"
#include "CppPyWrapper.hh"
#include <string>
#include <mpc.h>






		
		class  Functions: public OrbitUtils::CppPyWrapper
		{

			
		public:
		/*constructor*/	
		Functions(int n11,int n22, int mm, int point11);
		
		/*destructor*/	
		~Functions();
		
		int setupPrecision(std::string field,std::string c_energy, std::string c_Z1);


			


		std::string getM(std::string c_F,std::string c_energy, std::string c_Z1);


		


		private:  
			

			
			int n1;
			int n2;
			int m;
			int n;
			
			int point1;
			int point2;
			int point4;
			int nsum1;
			
			int prec_bit;
			int precision;
			int out_len;
			

			

			
//			tcomplex* C;
			mpc_t *C;
			
			mpc_t Z1;
			mpc_t Z2;
			mpc_t E;
			mpc_t M;
			
			mpc_t temp;


			

			


		
		};




#endif /*FUNCTIONS_HH_*/
