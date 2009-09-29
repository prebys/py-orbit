//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    HydrogenStarkParam.cc
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
#include "orbit_mpi.hh"
#include "Functions.hh"
#include <cmath>
#include <fstream>
#include <iostream>



using namespace OrbitUtils;
using namespace std;




Functions::Functions(int n11,int n22, int mm, int n_stepss, double stepp, int nsumm)	{		

	n1 = n11;
	n2 = n22;
	m = abs(mm);
	
	n_steps = n_stepss;
	step = stepp;
	nsum = nsumm;
	
	
}





Functions::~Functions()	{

}




void Functions::setupE(double EE)	{
	
	E = EE;
	
	return;
}

double Functions::getMfunctionModulus(double Energy, double Gamma, double reZ, double imZ){
	

		
	return E;
}


