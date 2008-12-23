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
#include "GaussianLaserField.hh"

#include "orbit_mpi.hh"
#include "OrbitConst.hh"
#include <iostream>

#define J tcomplex(0.,1.)

using namespace OrbitUtils;

GaussianLaserField::GaussianLaserField()
{
}

GaussianLaserField::~GaussianLaserField()
{
}



void	GaussianLaserField::setLaserHalfAngle(double a)	{	Laser_half_angle=a;}

double	GaussianLaserField::getLaserHalfAngle()	{	return Laser_half_angle;}

void	GaussianLaserField::setLaserPower(double a)	{	LaserPower=a;}

double	GaussianLaserField::getLaserPower()	{	return LaserPower;}

void	GaussianLaserField::setLaser_omega(double a)	{	Laser_omega=a;}

double	GaussianLaserField::getLaser_omega()	{	return Laser_omega;}












void GaussianLaserField::getLaserElectricField(double x, double y, double z, double t, tcomplex& E_x, tcomplex& E_y, tcomplex& E_z){
	
	
	double ta=2.418884326505e-17;			//atomic unit of time
	double Ea=5.14220642e11;				//Atomic unit of electric field
	t/=ta;
	
		double Rabi=1e+12*ta;					//Rabi frequensy (in atomic units)
		double Elas=64*Rabi*sqrt(2.)/27;		//Amplitude of laser field (in atonic units)
//		Ex_las=0.5773502691896258*Elas;			//components of laser field (in atomic units)
//		Ey_las=0.5773502691896258*Elas;	
//		Ez_las=0.5773502691896258*Elas;

		double Ex_las=Elas*Ea;
		double Ey_las=0.5773502691896258*Elas*Ea;
		double Ez_las=0.5773502691896258*Elas*Ea;
	
	

		double gamma_sweep=-OrbitConst::PI*Rabi*Rabi/2/log(1-0.87);
//		double omega_part=0.44417915081800191102;				// frequensy of laser in particle frame (in atomic units) 
//		double phasa=omega_part*t/ta;
		
		double phasa=4/9.*t+gamma_sweep*t*t/2;
	
	E_x=Ex_las*exp(J*phasa); E_y=0*Ey_las*exp(J*phasa); E_z=0*Ez_las*exp(J*phasa);		

}


void GaussianLaserField::getLaserMagneticField(double x, double y, double z, double t, tcomplex& B_x, tcomplex& B_y, tcomplex& B_z){
	
	B_x=tcomplex(0.0); B_y=tcomplex(0.0); B_z=tcomplex(0.0);
}


double GaussianLaserField::getFrequencyOmega(double x, double y, double z, double px, double py, double pz, double t){
	
	double ta=2.418884326505e-17;			//atomic unit of time
	t/=ta;
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////froissar stora test      froissar stora test      froissar stora test      froissar stora test    ///////////////   
	//next eight lines are testing parameters of laser electric field (in atomic units)
/*
	double Rabi=1e+12*ta;					//Rabi frequensy (in atomic units)
	double Elas=64*Rabi*sqrt(2.)/27;		//Amplitude of laser field (in atonic units)
	Ex_las=0.5773502691896258*Elas;			//components of laser field (in atomic units)
	Ey_las=0.5773502691896258*Elas;	
	Ez_las=0.5773502691896258*Elas;

	Ex_las=Elas;
	Ey_las=0;
	Ez_las=0;


	Ez_stat=0.00005;


	double gamma_sweep=-OrbitConst::PI*Rabi*Rabi/2/log(1-0.9);

	omega_part=0*0.44421741823240407099+4/9.+gamma_sweep*(t_part-10000*(2*OrbitConst::PI/Rabi)/2);				// frequensy of laser in particle frame (in atomic units) 
*/
	/////////////froissar stora test      froissar stora test      froissar stora test      froissar stora test ///////////          
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	

	double Rabi=1e+12*ta;
	double gamma_sweep=-OrbitConst::PI*Rabi*Rabi/2/log(1-0.87);
	return (4/9.+gamma_sweep*t)/ta;
	
}

//void 	LaserField::GetLabLaserField(double x, double y, double z, double t, tcomplex& E_x, tcomplex& E_y, tcomplex& E_z,tcomplex& B_x, tcomplex& B_y, tcomplex& B_z){
	
	/**Here must be the code of conversion of fields and coordinates from laser frame to lab frame**/
	
//	ElliptGaussian(x,y,z,t,E_x,E_y,E_z,B_x,B_y,B_z);
	
	/**Here must be the code of conversion of fields and coordinates from laser frame to lab frame**/
	
//}


