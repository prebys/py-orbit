//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    BaseFieldSource.hh
//
// CREATED
//    04/21/2003
//
// DESCRIPTION
//    The base class for field source. It should be sub-classed. The units
//    for E and B are unknown at this point. They should be defined in 
//    subclasses.
//
///////////////////////////////////////////////////////////////////////////
#ifndef BASELASERFIELDSOURCE_HH_
#define BASELASERFIELDSOURCE_HH_

#include "CppPyWrapper.hh"
#include <complex>
typedef std::complex<double>	tcomplex;
#define J tcomplex(0.,1.)


namespace OrbitUtils{
	class  BaseLaserFieldSource: public CppPyWrapper
	{
	public:
		
		BaseLaserFieldSource();
		virtual ~BaseLaserFieldSource();
		
		virtual void getLaserElectricField(double x, double y, double z, double t, tcomplex& f_x, tcomplex& f_y, tcomplex& f_z);
		virtual void getLaserMagneticField(double x, double y, double z, double t, tcomplex& f_x, tcomplex& f_y, tcomplex& f_z);
		virtual double getFrequencyOmega(double x, double y, double z, double px, double py, double pz, double t);

	};
};

#endif /*BASELASERFIELDSOURCE_HH_*/
