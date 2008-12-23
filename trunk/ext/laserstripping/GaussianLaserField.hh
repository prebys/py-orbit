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
#ifndef GAUSSIANLASERFIELD_HH_
#define GAUSSIANLASERFIELD_HH_

#include "Python.h"

#include "BaseLaserFieldSource.hh"

#include <complex>
typedef std::complex<double>	tcomplex;

namespace OrbitUtils{
	
	class  GaussianLaserField: public BaseLaserFieldSource{
		public:
		
			GaussianLaserField();
			~GaussianLaserField();
		
			void getLaserElectricField(double x, double y, double z, double t, tcomplex& f_x, tcomplex& f_y, tcomplex& f_z);
			void getLaserMagneticField(double x, double y, double z, double t, tcomplex& f_x, tcomplex& f_y, tcomplex& f_z);
			double getFrequencyOmega(double x, double y, double z, double px, double py, double pz, double t);
			
			
			/** It defines parameters of the laser beam**/
			void	setLaserHalfAngle(double a);
			double	getLaserHalfAngle();
			void	setLaserPower(double a);
			double	getLaserPower();
			void	setLaser_omega(double a);
			double	getLaser_omega();
			
		private:
			
			  //Laser parameters
			  double Laser_omega;
			  double LaserPower;
			  double Laser_half_angle;
			

	};
};



#endif /*GAUSSIANLASERFIELD_HH_*/
