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
#ifndef LSFIELDSOURCE_HH_
#define LSFIELDSOURCE_HH_

#include "Python.h"

#include "BaseFieldSource.hh"

namespace OrbitUtils{
	
	class  LSFieldSource: public BaseFieldSource
	{
		public:
		
			LSFieldSource();
			LSFieldSource(double E_x,double E_y,double E_z,double B_x,double B_y,double B_z);
			~LSFieldSource();
		
			void getElectricMagneticField(double x, double y, double z, double t, double& E_x, double& E_y, double& E_z, double& H_x, double& H_y, double& H_z);

		private:
			
			void getConstField(double& E_x, double& E_y, double& E_z, double& H_x, double& H_y, double& H_z);
			void getLSfield(double x, double y, double z, double t,double& E_x, double& E_y, double& E_z, double& H_x, double& H_y, double& H_z);

			
			bool flag_Const_field;
			bool flag_LSfield;
			
		
			double ExConst;
			double EyConst;
			double EzConst;
			double BxConst;
			double ByConst;
			double BzConst;
			
			
	};
};




#endif /*LSFIELDSOURCE_HH_*/



