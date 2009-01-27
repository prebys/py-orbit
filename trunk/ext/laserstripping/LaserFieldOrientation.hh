//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    LorentzTransformationEM.hh
//
// AUTHOR
//    T. Gorlov
//
// CREATED
//    06/28/2008
//
// DESCRIPTION
//    This class provides Lorentz transformations for the electromagnetic field
//    from the laboratory frame to the particle rest frame.
//    mass - mass of the particle in GeV
//    px,py,pz - momentum of the particle the lab frame in GeV/c
//    E_x,E_y,E_z - components of the electric field V/m (parameters are replaced in place) 
//    B_x,B_y,B_z - components of the magnetic field [T] (parameters are replaced in place)   
//
//    OrbitConst::c in [m/sec]
///////////////////////////////////////////////////////////////////////////
#ifndef LASERFIELDORIENTATION_HH_
#define LASERFIELDORIENTATION_HH_


#include <cmath>

class  LaserFieldOrientation
{
public:
	
  static  void	OrientCoordinates(double& x, double& y,double& z,
		  double x0, double y0,double z0,
		  double kx, double ky,double kz,
		  double mx, double my,double mz);
  

  
};






#endif /*LASERFIELDORIENTATION_HH_*/
