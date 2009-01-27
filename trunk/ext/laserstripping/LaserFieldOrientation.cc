//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    LorentzTransformationEM.cc
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
//
///////////////////////////////////////////////////////////////////////////
#include "LaserFieldOrientation.hh"


void LaserFieldOrientation::OrientCoordinates(double& x, double& y,double& z,
		  double x0, double y0,double z0,
		  double kx, double ky,double kz,
		  double mx, double my,double mz)	{
	
double xt=x-x0;
double yt=y-y0;
double zt=z-z0;

double a=sqrt(pow(kx*my-ky*mx,2)+pow(kx*mz-kz*mx,2)+pow(ky*mz-kz*my,2));
double b=sqrt(kx*kx+ky*ky+kz*kz);

x=	(xt*((ky*ky+kz*kz)*mx-kx*(ky*my+kz*mz))
	+yt*((kz*kz+kx*kx)*my-ky*(kz*mz+kx*mx))
	+zt*((kx*kx+ky*ky)*mz-kz*(kx*mx+ky*my)))/(a*b);

y=	(xt*(ky*mz-kz*my)
	+yt*(kz*mx-kx*mz)
	+zt*(kx*my-ky*mx))/a;

z=	(kx*xt+ky*yt+kz*zt)/b;


}




