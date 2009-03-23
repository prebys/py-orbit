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
#include "HermiteGaussianLFmode.hh"
#include "MathPolynomial.hh"

#include "orbit_mpi.hh"
#include "OrbitConst.hh"
#include "LaserFieldOrientation.hh"
#include <iostream>

#define J tcomplex(0.,1.)

using namespace OrbitUtils;

HermiteGaussianLFmode::HermiteGaussianLFmode(double Cnm,int n,int m,double w_x,double w_y,double f_x,double f_y,double la)
{
	Unm=Cnm*sqrt(2*OrbitConst::c*OrbitConst::permeability)*sqrt(pow(w_x,2*n+1)*pow(w_y,2*m+1)*MathPolynomial::Factorial(n)*MathPolynomial::Factorial(m)*pow(2,n+m+1)/(MathPolynomial::Factorial(2*n)*MathPolynomial::Factorial(2*m)*OrbitConst::PI));

	Laser_lambda=la;
	fx=f_x;
	fy=f_y;
	wx=w_x;
	wy=w_y;
	n_moda=n;
	m_moda=m;
	

	//default orientation of laser field
	nEx=1;nEy=0;nEz=0;
	nHx=0;nHy=1;nHz=0;
	kx=0;ky=0;kz=1;
	mx=1;my=0;mz=0;
	
	
	
}

HermiteGaussianLFmode::~HermiteGaussianLFmode()
{
}









void HermiteGaussianLFmode::setLaserFieldOrientation(double x_0, double y_0, double z_0,
													double k_x, double k_y, double k_z,
													double m_x, double m_y, double m_z,
													double n_Ex, double n_Ey, double n_Ez)	
{
	double a=sqrt(n_Ex*n_Ex+n_Ey*n_Ey+n_Ez*n_Ez);

	
	x0=x_0;	y0=y_0;	z0=z_0;	//Point in space 
	kx=k_x;	ky=k_y;	kz=k_z;	//Wave vector
	mx=m_x;	my=m_y;	mz=m_z; //Vector of moda orientation
	
	nEx=n_Ex/a;	
	nEy=n_Ey/a;	
	nEz=n_Ez/a;	//Vector of polarization of electric field
	
double	n_Hx=nEz*ky-nEy*kz;
double	n_Hy=nEx*kz-nEz*kx;
double	n_Hz=nEy*kx-nEx*ky;
	
double b=sqrt(n_Hx*n_Hx+n_Hy*n_Hy+n_Hz*n_Hz);
	
	nHx=n_Hx/b;
	nHy=n_Hy/b;
	nHz=n_Hz/b;
	

}









double HermiteGaussianLFmode::HermiteGaussianOmega(double m,double x, double y, double z, double px, double py, double pz, double t){
	
	
	double a=OrbitConst::c/sqrt(m*m+px*px+py*py+pz*pz);
	
	double vx=px*a;
	double vy=py*a;
	double vz=pz*a;
	
	double k=2*OrbitConst::PI/Laser_lambda;

double ax=k*k*pow(wx,4)+4*pow(fx-z,2);
double ay=k*k*pow(wy,4)+4*pow(fy-z,2);


	return 	k*(OrbitConst::c-vz-4*k*k*vz*(pow(wx*wx*x/ax,2)+pow(wy*wy*y/ay,2)) + (vz*(wx*wx+ 2*x*x) + 4*vx*x*(fx - z))/ax + (vz*(wy*wy + 2*y*y) + 4*vy*y*(fy - z))/ay);

}












void HermiteGaussianLFmode::getLaserElectricMagneticField(double x, double y, double z, double t, 
		tcomplex& E_x, tcomplex& E_y, tcomplex& E_z,
		tcomplex& H_x, tcomplex& H_y, tcomplex& H_z){
	


LaserFieldOrientation::OrientCoordinates(x,y,z,x0,y0,z0,kx,ky,kz,mx,my,mz);
tcomplex	E=Unm*getNonOrientedU(n_moda,m_moda,x,y,z,t);


tcomplex	H=E/OrbitConst::c;



	E_x=E*nEx;	E_y=E*nEy;	E_z=E*nEz;
	H_x=H*nHx;	H_y=H*nHy;	H_z=H*nHz;
	
}







double HermiteGaussianLFmode::getFrequencyOmega(double m, double x, double y, double z, double px, double py, double pz, double t){

	
	LaserFieldOrientation::OrientCoordinates(x,y,z,x0,y0,z0,kx,ky,kz,mx,my,mz);
	LaserFieldOrientation::OrientCoordinates(px,py,pz,0.,0.,0.,kx,ky,kz,mx,my,mz);
	

	return HermiteGaussianOmega(m,x,y,z,px,py,pz,t);

}











tcomplex HermiteGaussianLFmode::getNonOrientedU(int n, int m, double x, double y, double z, double t){

	
	if (n==0&&m==0)	{
		
		double a=Laser_lambda/OrbitConst::PI;
		tcomplex k=tcomplex(0.,2*(z-t*OrbitConst::c)/a);
		tcomplex funx=pow(tcomplex(wx*wx,-(z-fx)*a),-1);
		tcomplex funy=pow(tcomplex(wy*wy,-(z-fy)*a),-1);
		
		return	sqrt(funx*funy)*exp(-x*x*funx-y*y*funy-k);
	
	}	else	{
	
	
	
	tcomplex k=2*OrbitConst::PI/Laser_lambda;
	tcomplex funx=sqrt(tcomplex(wx*wx-2.0*J*(z-fx)/k));
	tcomplex funy=sqrt(tcomplex(wy*wy-2.0*J*(z-fy)/k));
	tcomplex xf=x/funx;
	tcomplex yf=y/funy;
	
	return	pow(funx,-n-1)*pow(funy,-m-1)*exp(-xf*xf-yf*yf-J*k*(z-t*OrbitConst::c))*MathPolynomial::ComplexHermite(n,xf)*MathPolynomial::ComplexHermite(m,yf);
	
	}
	

	
}










