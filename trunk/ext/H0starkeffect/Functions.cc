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
#include <iomanip>




using namespace OrbitUtils;
using namespace std;




Functions::Functions(int n11,int n22, int mm, int point11)	{		
	

	n1 = n11;
	n2 = n22;
	m = abs(mm);
	n = n1 + n2 + m + 1; 
	
	point1 = point11;
	point2 = point1*point1;
	point4 = point1*point1*point1*point1;
	
	nsum1 = 10;
	precision = 10;
	prec_bit = (int)(precision*3.3219280948873626);
	
	C = new mpc_t[nsum1 + 1];
	for (int i=0;i<nsum1+1;i++)	
		mpc_init2(C[i], prec_bit);

	


}





Functions::~Functions()	{
	
	for (int i=0;i<nsum1+1;i++)	
		mpc_clear(C[i]);
	delete [] C;



}




int Functions::setupPrecision(std::string c_field,std::string c_energy, std::string c_Z1)	{

	for (int i=0;i<nsum1+1;i++)	
		mpc_clear(C[i]);
	delete [] C;


	std::string line = "empty line", Mstr;
	
	int sum = 0;
	int prec = 1;
	int _temp;
	int count = 0;
	int exponent;
	
	mpfr_t temp_mpfr;
	mpc_t temp_mpc;
	
	for (int j=0;j<2;j++) 	{

		
		
	precision = 10;
	nsum1 = 10;


		
	while (true)	{
		
		
		nsum1 += 10*sum;
		precision += 10*prec;
		
		prec_bit = (int)(precision*3.3219280948873626);
		
		C = new mpc_t[nsum1 + 1];
		for (int i=0;i<nsum1+1;i++)	
			mpc_init2(C[i], prec_bit);

		if(j==0) out_len = 10;
		if(j==1) out_len = exponent;


		Mstr = getM(c_field,c_energy, c_Z1);
		
		if (line != Mstr)	{
			line.clear();
			line = Mstr;		
//			cout<<nsum1<<"  "<<precision<<"  "<<line<<"\n";
			count = 0;
			
		} else {
			
			_temp = sum;
			sum = prec;
			prec = _temp;
			count++;
		}
			
		if (count == 2)
		{break;	}
		else	{
	
		for (int i=0;i<nsum1+1;i++)	
			mpc_clear(C[i]);
		delete [] C;
		}

	
	}

		mpfr_init2(temp_mpfr, prec_bit);
		mpc_init2(temp_mpc, prec_bit);
	
		mpc_set_str(temp_mpc, &Mstr[0], 10, MPC_RNDNN);
		mpc_abs(temp_mpfr,temp_mpc,GMP_RNDN);
		mpfr_log10(temp_mpfr, temp_mpfr,GMP_RNDN);
		
		exponent = (int)mpfr_get_d(temp_mpfr,GMP_RNDN);
		if (j==1) cout<<"exponent= "<<exponent<<"\n";
		
		mpfr_clear(temp_mpfr);
		mpc_clear(temp_mpc);
		
	}

	
	out_len = precision;
	
	return prec_bit + 10;
}














std::string Functions::getM(std::string c_F, std::string c_energy, std::string c_Z1){
	
	mpfr_t F;
	mpc_t E;
	mpc_t Z1;
	mpc_t M;
	mpc_t temp;
	std::string str;
	
	mpfr_init2(F, prec_bit);
	mpc_init2(E, prec_bit);
	mpc_init2(Z1, prec_bit);
	mpc_init2(M, prec_bit);
	mpc_init2(temp, prec_bit);
	
	
	
	mpfr_set_str(F, &c_F[0], 10, GMP_RNDN);
	mpc_set_str(E, &c_energy[0], 10, MPC_RNDNN);
	mpc_set_str(Z1, &c_Z1[0], 10, MPC_RNDNN);


	//Calculation of M - function and its derivative at the point mu = point1
	
	//Calculation of C[0]
	mpc_set_ui(C[0], 1, MPC_RNDNN);

	
	//Calculation of C[1]*mu*mu
	mpc_mul_si(C[1],Z1,-point2, MPC_RNDNN);	
	mpc_div_ui(C[1],C[1],4*(1 + m), MPC_RNDNN);

	//Calculation of C[3]*mu*mu*mu*mu*mu*mu
	mpc_mul(C[2],Z1,Z1,MPC_RNDNN); 
	mpc_mul_ui(temp,E,8*(1 + m),MPC_RNDNN); 
	mpc_sub(C[2],C[2],temp,MPC_RNDNN); 
	mpc_mul_ui(C[2],C[2],point4,MPC_RNDNN); 
	mpc_div_ui(C[2],C[2],32*(m + 1)*(m + 2),MPC_RNDNN);	
	

	//Calculation of C[i]*pow(mu,2*i)
	for (int i=3;i<nsum1+1;i++)	{
		
		mpc_mul_fr(C[i],C[i-3],F,MPC_RNDNN); 
		mpc_mul_ui(C[i],C[i],point4,MPC_RNDNN); 
		mpc_mul(temp,C[i-1],Z1,MPC_RNDNN);
		mpc_sub(C[i],C[i],temp,MPC_RNDNN);
		mpc_mul(temp,C[i-2],E,MPC_RNDNN); 
		mpc_mul_ui(temp,temp,2*point2,MPC_RNDNN);
		mpc_sub(C[i],C[i],temp,MPC_RNDNN);
		mpc_mul_ui(C[i],C[i],point2,MPC_RNDNN);
		mpc_div_ui(C[i],C[i],4*i*(i + m),MPC_RNDNN);
		
	}

	//Calculation of M(mu) 
	mpc_set_ui(M, 0, MPC_RNDNN);
	for (int i=nsum1;i>-1;i--)	
		mpc_add(M,M,C[i],MPC_RNDNN);
	
	mpc_set_ui(temp,1 + 2*m, MPC_RNDNN);
	mpc_div_ui(temp,temp,2, MPC_RNDNN);
	mpc_set_ui(C[0], point1, MPC_RNDNN);
	mpc_pow(temp, C[0],temp, MPC_RNDNN);
	
	mpc_mul(M,M,temp,MPC_RNDNN);
	
	str = (std::string)mpc_get_str(10, out_len, M, MPC_RNDNN);
	
	mpfr_clear(F);
	mpc_clear(E);
	mpc_clear(Z1);
	mpc_clear(M);
	mpc_clear(temp);
	
	
	return str;

/*	
	tcomplex Z1 = tcomplex(reZ1, imZ1);
	tcomplex E = tcomplex(Energy, -Gamma/2.0);
	tcomplex sum = 0,sum1 = 0; 
	tcomplex M, M_der;


	//Calculation of M - function and its derivative at the point mu = "step"
	C[0] = 1.0;	
	C[1] = -point1*point1*Z1/(4.*(1 + m));	
	C[2] = point1*point1*point1*point1*(Z1*Z1 - 8.*(1 + m)*E)/(32.*(m + 1)*(m + 2));		

	
	for (int i=3;i<nsum1+1;i++)	
		C[i]=point1*point1*(point1*point1*point1*point1*F*C[i-3]-Z1*C[i-1]-2.*point1*point1*E*C[i-2])/(4.*i*(i + m));
	

	for (int i=nsum1;i>-1;i--)	{
		sum += C[i];
		sum1 += C[i]*(2*i + m + 0.5);

	}
	
	
	M = pow(point1,m)*sqrt(point1)*sum;
	M_der = pow(point1,m-1)*sqrt(point1)*sum1;
	

	
	
	cout<<setprecision(20)<<M.real()<<"\n";
	cout<<setprecision(20)<<M_der.real()<<"\n\n";
	*/
	
	

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
	
	
	
	
/*	
	
	for (int j=1;j<10;j++)	{
	
	

	double a1 = j*step;
	double a2 = a1*a1;
	double a3 = a2*a1;
	double a4 = a3*a1;
	double a5 = a4*a1;
	double a6 = a5*a1;

	
	tcomplex a = 2*a4*E - a6*F + (0.25 - m*m) + a2*Z1;
	tcomplex b = 8*a3*E - 6*a5*F +2*a1*Z1;
	tcomplex c = 12*a2*E - 15*a4*F + Z1;
	tcomplex d = 8*a1*E - 20*a5*F;
	tcomplex e = 2.*E - 15*a2*F;
	tcomplex f = -6*a1*F;
	
	return 0;
}


double Functions::getMfunctionModulus(double Energy, double Gamma, double reZ1, double imZ1){
	
	mpc_set_d_d(E, Energy,-Gamma/2, MPC_RNDNN);
	mpc_set_d_d(Z1, reZ1, imZ1, MPC_RNDNN);
	
	//Calculation of M - function and its derivative at the point mu = point1
	
	//Calculation of C[0]
	mpc_set_ui(C[0], 1, MPC_RNDNN);
	
	//Calculation of C[1]*mu*mu
	mpc_mul_si(C[1],Z1,-point2, MPC_RNDNN);	
	mpc_div_ui(C[1],C[1],4*(1 + m), MPC_RNDNN);
	
	//Calculation of C[3]*mu*mu*mu*mu*mu*mu
	mpc_mul(C[2],Z1,Z1,MPC_RNDNN); 
	mpc_mul_ui(temp,E,8*(1 + m),MPC_RNDNN); 
	mpc_sub(C[2],C[2],temp,MPC_RNDNN); 
	mpc_mul_ui(C[2],C[2],point4,MPC_RNDNN); 
	mpc_div_ui(C[2],C[2],32*(m + 1)*(m + 2),MPC_RNDNN);	
	
	
	//Calculation of C[i]*pow(mu,2*i)
	for (int i=3;i<nsum1+1;i++)	{
		
		mpc_mul(C[i],C[i-3],F,MPC_RNDNN); 
		mpc_mul_ui(C[i],C[i],point4,MPC_RNDNN); 
		mpc_mul(temp,C[i-1],Z1,MPC_RNDNN);
		mpc_sub(C[i],C[i],temp,MPC_RNDNN);
		mpc_mul(temp,C[i-2],E,MPC_RNDNN); 
		mpc_mul_ui(temp,temp,2*point2,MPC_RNDNN);
		mpc_sub(C[i],C[i],temp,MPC_RNDNN);
		mpc_mul_ui(C[i],C[i],point2,MPC_RNDNN);
		mpc_div_ui(C[i],C[i],4*i*(i + m),MPC_RNDNN);
		
	}
	//		C[i]=point1*point1*(point1*point1*point1*point1*F*C[i-3]-Z1*C[i-1]-2.*point1*point1*E*C[i-2])/(4.*i*(i + m));

	//Calculation of M(mu) 
	mpc_set_ui(M, 0, MPC_RNDNN);
	for (int i=nsum1;i>-1;i--)	
		mpc_add(M,M,C[i],MPC_RNDNN);
	
	mpc_set_d(temp, pow(point1*1.,m)*sqrt(point1*1.), MPC_RNDNN);
	mpc_mul(M,M,temp,MPC_RNDNN);


	
	cout<<mpc_get_str(10, 20, M, MPC_RNDNN)<<"\n";


	
/*	
	tcomplex Z1 = tcomplex(reZ1, imZ1);
	tcomplex E = tcomplex(Energy, -Gamma/2.0);
	tcomplex sum = 0,sum1 = 0; 
	tcomplex M, M_der;


	//Calculation of M - function and its derivative at the point mu = "step"
	C[0] = 1.0;	
	C[1] = -point1*point1*Z1/(4.*(1 + m));	
	C[2] = point1*point1*point1*point1*(Z1*Z1 - 8.*(1 + m)*E)/(32.*(m + 1)*(m + 2));		

	
	for (int i=3;i<nsum1+1;i++)	
		C[i]=point1*point1*(point1*point1*point1*point1*F*C[i-3]-Z1*C[i-1]-2.*point1*point1*E*C[i-2])/(4.*i*(i + m));
	

	for (int i=nsum1;i>-1;i--)	{
		sum += C[i];
		sum1 += C[i]*(2*i + m + 0.5);

	}
	
	
	M = pow(point1,m)*sqrt(point1)*sum;
	M_der = pow(point1,m-1)*sqrt(point1)*sum1;
	

	
	
	cout<<setprecision(20)<<M.real()<<"\n";
	cout<<setprecision(20)<<M_der.real()<<"\n\n";
	*/
	
	

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
	
	
	
	
/*	
	
	for (int j=1;j<10;j++)	{
	
	

	double a1 = j*step;
	double a2 = a1*a1;
	double a3 = a2*a1;
	double a4 = a3*a1;
	double a5 = a4*a1;
	double a6 = a5*a1;

	
	tcomplex a = 2*a4*E - a6*F + (0.25 - m*m) + a2*Z1;
	tcomplex b = 8*a3*E - 6*a5*F +2*a1*Z1;
	tcomplex c = 12*a2*E - 15*a4*F + Z1;
	tcomplex d = 8*a1*E - 20*a5*F;
	tcomplex e = 2.*E - 15*a2*F;
	tcomplex f = -6*a1*F;
	
	C[0] = M;
	C[1] = M_der;
	C[2] = -C[0]*a/(2*a2);
	C[3] = -(4*a1*C[2] + C[0]*b + C[1]*a)/(6*a2);
	C[4] = -(2.*C[2] + 12*a1*C[3] + C[0]*c + C[1]*b + C[2]*a)/(12*a2);
	C[5] = -(d*C[0] + 6.*C[3] + 24*a1*C[4] + c*C[1] + b*C[2] + a*C[3])/(20*a2);
	C[6] = -(e*C[0] + d*C[1] + 12.*C[4] + 40*a1*C[5] + c*C[2] + b*C[3] + a*C[4])/(30*a2);
	C[7] = -(f*C[0] + e*C[1] + d*C[2] + 20.*C[5] + 60*a1*C[6] + c*C[3] + b*C[4] + a*C[5])/(42*a2);
	
	for (int i=8;i<nsum+1;i++)
		C[i] = -(-F*C[i-8] + f*C[i-7] + e*C[i-6] + d*C[i-5] + ((i-2.)*(i-3.) + a)*C[i-2] + 2*(i-1)*(i-2)*a1*C[i-1] + c*C[i-4] + b*C[i-3] + a*C[i-2])/(i*(i-1)*a2);
	
	M = 0;
	M_der = 0;
	for (int i=nsum;i>-1;i--)	{
//		cout<<setprecision(20)<<C[i].real()<<"\n";
		M += C[i];
		M_der += i*1.*C[i];

	}
	cout<<setprecision(20)<<M.real()<<"\n";
	cout<<setprecision(20)<<M_der.real()<<"\n\n";
	
	}
	
*/

	


}





















