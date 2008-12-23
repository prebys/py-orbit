//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    BaseFieldSource.cc
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


#include "BaseLaserFieldSource.hh"

using namespace OrbitUtils;

BaseLaserFieldSource::BaseLaserFieldSource(): CppPyWrapper(NULL)
{
}

BaseLaserFieldSource::~BaseLaserFieldSource()
{
}

void BaseLaserFieldSource::getLaserElectricField(double x, double y, double z, double t, tcomplex& f_x, tcomplex& f_y, tcomplex& f_z)
{
  f_x = 0.0; f_y = 0.0; f_z = 0.0; 
}

void BaseLaserFieldSource::getLaserMagneticField(double x, double y, double z, double t, tcomplex& f_x, tcomplex& f_y, tcomplex& f_z)
{
  f_x = 0.0; f_y = 0.0; f_z = 0.0; 
}

double BaseLaserFieldSource::getFrequencyOmega(double x, double y, double z,double px, double py, double pz, double t)
{
	return 0.;
}


