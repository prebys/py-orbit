//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    CppBaseFieldSource.cc
//
// CREATED
//    04/21/2003
//
// DESCRIPTION
//    The base class for Python implementation of a field source. 
//    It should be sub-classed on Python level and implements 
//    getElectricField(x,y,z,t) and getMagneticField (x,y,z,t) methods.
//    The results of these methods will be available from the c++ level.
//    This is an example of embedding Python in C++ Orbit level.
//
//
///////////////////////////////////////////////////////////////////////////
#include "LSFieldSource.hh"

#include "orbit_mpi.hh"
#include <iostream>

using namespace OrbitUtils;

LSFieldSource::LSFieldSource()
{
}

LSFieldSource::~LSFieldSource()
{
}

void LSFieldSource::getElectricMagneticField(double x, double y, double z, double t, 
		double& E_x, double& E_y, double& E_z,
		double& H_x, double& H_y, double& H_z)

{	  
	 E_x = 0.0; E_y = 0.0; E_z = 0.0; 
	 H_x = 300.0/10000; H_y = 0.0; H_z = 0.0; 
}


