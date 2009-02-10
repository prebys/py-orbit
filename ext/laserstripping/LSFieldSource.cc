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
flag_Const_field=false;
flag_LSfield=true;
}

LSFieldSource::LSFieldSource(double E_x,double E_y,double E_z,double B_x,double B_y,double B_z)
{
	flag_Const_field=true;
	flag_LSfield=false;
	
	ExConst=E_x;
	EyConst=E_y;
	EzConst=E_z;
	
	BxConst=B_x;
	ByConst=B_y;
	BzConst=B_z;
	

}


LSFieldSource::~LSFieldSource()
{
}

void LSFieldSource::getElectricMagneticField(double x, double y, double z, double t, 
		double& E_x, double& E_y, double& E_z,
		double& H_x, double& H_y, double& H_z)

{
	
if	(flag_Const_field) getConstField(E_x, E_y, E_z, H_x, H_y, H_z);
if	(flag_LSfield) getLSfield(x, y, z, t, E_x, E_y, E_z, H_x, H_y, H_z);

}




void LSFieldSource::getConstField(
		double& E_x, double& E_y, double& E_z,
		double& H_x, double& H_y, double& H_z)

{	  
	 E_x = ExConst; E_y = EyConst; E_z = EzConst; 
	 H_x = BxConst; H_y = ByConst; H_z = BzConst; 

}

void LSFieldSource::getLSfield(double x, double y, double z, double t, 
		double& E_x, double& E_y, double& E_z,
		double& H_x, double& H_y, double& H_z)

{	  
	 E_x = 0.0; E_y = 0.0; E_z = 0.0; 
	 H_x = 0.0; H_y = 0.0; H_z = 0.0; 

}



