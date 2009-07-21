#ifndef FIELDROTATION_HH_
#define FIELDROTATION_HH_


#include "tcomplex.hh"


class  FieldRotation
{
public:

	 FieldRotation();
	 ~FieldRotation();

static double RotateElectricFields(double vx,double vy,double vz, tcomplex& Exl,tcomplex& Eyl,tcomplex& Ezl);


};

///////////////////////////////////////////////////////////////////////////
//
// END OF FILE
//
///////////////////////////////////////////////////////////////////////////


#endif /*FIELDROTATION_HH_*/
