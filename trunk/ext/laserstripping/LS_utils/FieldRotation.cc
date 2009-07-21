#include "FieldRotation.hh"




double	FieldRotation::RotateElectricFields(double vx,double vy,double vz, tcomplex& Exl,tcomplex& Eyl,tcomplex& Ezl){
	
	double vxy2=vx*vx+vy*vy;
	double v2=vxy2+vz*vz;
	double v=sqrt(v2);



		
		tcomplex Exll=Exl;
		tcomplex Eyll=Eyl;
		tcomplex Ezll=Ezl;
		
		
		if(vxy2>1e-40*v2)	{
			
		Exl=(-vx*vxy2*Ezll + Eyll*vx*vy*(-v + vz) + Exll*(v*vy*vy + vx*vx*vz))/(v*vxy2);
		Eyl=(-vy*vxy2*Ezll + Exll*vx*vy*(-v + vz) + Eyll*(v*vx*vx + vy*vy*vz))/(v*vxy2);
		Ezl=(Exll*vx + Eyll*vy + Ezll*vz)/v;

		}
		

		else	{
			
			Exl=Exl;
			Eyl=Eyl;
			Ezl=Ezl;

		}
	
return v;

}
