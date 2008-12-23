
#include <complex>

typedef std::complex<double>	tcomplex;



#ifndef HYDROGENSTARKPARAM_HH_
#define HYDROGENSTARKPARAM_HH_




		
		class  HydrogenStarkParam
		{

			
		public:
			/*this method gives koefficient of autoionization as a function of field  */
		static	void GetEnergyAutoionization(int k,double& E, double& Gamma);
			
			/*this method gives frequendy transitions as a function of field */
		static	void GetDipoleTransition(int k,int ks,tcomplex& mu_x,tcomplex& mu_y,tcomplex& mu_z);	
			
			/*this method gives frequendy transitions as a function of field */
		static	void GetRelax(int k,int ks,double& relax);
		
			/*this method gives frequendy transitions as a function of field */
		static	void ReadData(char* addressEG,int states);
		
			/*this method setups electrostatic field */
		static	void SetE(double E);
		
		
		
		
		
		private:  
		  //parameters of dipole transition, energy, lifetime, spontaneous relaxation, 
		  
		 static double*** dipole_transition_x;
		 static double*** dipole_transition_y;
		 static double*** dipole_transition_z;
		 static double*** gamma_spontaneous_relax;
		 static double** energy;
		 static double** gamma_autoionization;
		  
		 static double delta_F;
		 static int n_data;
		 static double Ez_stat;
		
		};




#endif /*HYDROGENSTARKPARAM_HH_*/
