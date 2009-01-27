//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    PyBaseFieldSource.hh
//
// CREATED
//    04/21/2003
//
// DESCRIPTION
//    The base class for Python implementation of a field source. 
//    It should be sub-classed on Python level and implement 
//    getElectricField(x,y,z,t) and getMagneticField (x,y,z,t) methods.
//    The results of these methods will be available from the c++ level.
//    This is an example of embedding Python in C++ Orbit level.
//
//    ATTENTION: Using this class in real calculations is not wise! 
//               It is slow, because it delegates the field calculations
//               to the Python level. It is for prototyping and 
//               debugging only.
//
///////////////////////////////////////////////////////////////////////////
#ifndef HERMITE_GAUSSIAN_LF_MODE_HH_
#define HERMITE_GAUSSIAN_LF_MODE_HH_

#include "Python.h"

#include "BaseLaserFieldSource.hh"

#include <complex>
typedef std::complex<double>	tcomplex;

namespace OrbitUtils{
	
	class  HermiteGaussianLFmode: public BaseLaserFieldSource{
		public:
		
			HermiteGaussianLFmode(double Cnm,int n,int m,double wx,double wy,double f_x,double f_y,double la);
			~HermiteGaussianLFmode();
		
			void getLaserElectricMagneticField(double x, double y, double z, double t, 
					tcomplex& E_x, tcomplex& E_y, tcomplex& E_z,
					tcomplex& H_x, tcomplex& H_y, tcomplex& H_z);

					
			tcomplex getNonOrientedU(int n, int m, double x, double y, double z, double t);
		
			
			double getFrequencyOmega(double m, double x, double y, double z, double px, double py, double pz, double t);
			
			void setLaserFieldOrientation(double x_0, double y_0, double z_0,
										double k_x, double k_y, double k_z,
										double m_x, double m_y, double m_z,
										double n_Ex, double n_Ey, double n_Ez);
			
			

			
		private:
			
			
			double FroissartStoraTestOmega(double x, double y, double z, double px, double py, double pz, double t);
			double HermiteGaussianOmega(double m, double x, double y, double z, double px, double py, double pz, double t);
			tcomplex E_FroissartStora(double x, double y, double z, double t);
			
			

			
			  //Laser parameters
			  tcomplex Unm;
			  double Laser_lambda;
			  int	n_moda;
			  int	m_moda;
			  double wx;
			  double wy;
			  double fx;
			  double fy;
			
			  //parameters of orientation of laser field
			  double x0,y0,z0,kx,ky,kz,mx,my,mz,nEx,nEy,nEz;
			  double nHx, nHy, nHz;
			

	};
};



#endif /*HERMITE_GAUSSIAN_LF_MODE_HH_*/
