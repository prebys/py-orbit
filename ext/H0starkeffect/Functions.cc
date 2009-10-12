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




Functions::Functions(int n11,int n22, int mm, int point11, int _err_exp)	{		
	

	n1 = n11;
	n2 = n22;
	m = abs(mm);
	n = n1 + n2 + m + 1; 
	
	point1 = point11;
	point2 = point1*point1;
	
	out_lenN = 1000;
	pointN1 = 200;
	pointN2 = pointN1*pointN1;
	
	
	precisionN = 3000;
	prec_bitN = (int)(precisionN*3.3219280948873626);

	err_exp = _err_exp;
	
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




long int Functions::calcPrecisionForM(std::string c_field,std::string c_energy, std::string c_Z1)	{

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
		if(j==1) out_len = exponent - err_exp;


		Mstr = getM(c_field,c_energy, c_Z1);
		
		if (line != Mstr)	{
			line.clear();
			line = Mstr;		
			cout<<out_len<<" "<<nsum1<<"  "<<precision<<"  "<<line<<"\n";
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
		if (exponent<5) exponent = 5;
		if (j==1) {cout<<"exponent= "<<exponent<<"\n";cout<<"nsum1= "<<nsum1<<"\n";}
		
		mpfr_clear(temp_mpfr);
		mpc_clear(temp_mpc);
		
	}	
	
	
//	out_len = precision;
	
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
	mpc_mul_ui(C[2],C[2],point2,MPC_RNDNN);
	mpc_mul_ui(C[2],C[2],point2,MPC_RNDNN);
	mpc_div_ui(C[2],C[2],32*(m + 1)*(m + 2),MPC_RNDNN);	
	

	//Calculation of C[i]*pow(mu,2*i)
	for (int i=3;i<nsum1+1;i++)	{
		
		mpc_mul_fr(C[i],C[i-3],F,MPC_RNDNN); 
		mpc_mul_ui(C[i],C[i],point2,MPC_RNDNN); 
		mpc_mul_ui(C[i],C[i],point2,MPC_RNDNN); 
		mpc_mul(temp,C[i-1],Z1,MPC_RNDNN);
		mpc_sub(C[i],C[i],temp,MPC_RNDNN);
		mpc_mul(temp,C[i-2],E,MPC_RNDNN); 
		mpc_mul_ui(temp,temp,2*point2,MPC_RNDNN);
		mpc_sub(C[i],C[i],temp,MPC_RNDNN);
		mpc_mul_ui(C[i],C[i],point2,MPC_RNDNN);
		mpc_div_ui(C[i],C[i],4,MPC_RNDNN);
		mpc_div_ui(C[i],C[i],i,MPC_RNDNN);
		mpc_div_ui(C[i],C[i],i + m,MPC_RNDNN);
		
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

/*
int Functions::getN(std::string& str_N, std::string& str_derN,std::string c_F, std::string c_energy, std::string c_Z2)	{
	
	
	for (int i=0;i<nsum1+1;i++)	
		mpc_clear(C[i]);
	delete [] C;


	std::string line = "empty line";
	std::string line1 = "empty line";

	
	int sum = 0;
	int prec = 1;
	int _temp;
	int count = 0;

	
	mpfr_t temp_mpfr;
	mpc_t temp_mpc;
	
	point1 = 100;
	point2 = point1*point1;
	point3 = point1*point1*point1;
	point4 = point1*point1*point1*point1;

		
		
	precision = 10;
	nsum1 = 10;


		
	while (true)	{
		
		
		nsum1 += 10*sum;
		precision += 10*prec;
		
		prec_bit = (int)(precision*3.3219280948873626);
		
		C = new mpc_t[nsum1 + 1];
		for (int i=0;i<nsum1+1;i++)	
			mpc_init2(C[i], prec_bit);

//		out_len = 200;
		getfastN(str_N, str_derN, c_F,  c_energy,c_Z2);

		
		if ((line != str_N)||(line1 != str_derN))	{
			line.clear();
			line1.clear();
			line = str_N;
			line1 = str_derN;	
			cout<<nsum1<<"  "<<precision<<"\n"<<str_N<<"\n"<<str_derN<<"\n\n";
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

	
	
	
//	out_len = precision;
	
	return prec_bit + 10;
}
*/


long int Functions::getN(std::string& str_N, std::string& str_derN,std::string c_F, std::string c_energy, std::string c_Z2)	{
	
	long int crit;
	


	
	while(true)	{
		crit =  out_lenN + precisionN - getfastN(str_N,  str_derN, c_F,  c_energy,  c_Z2);
		if (crit>=700)	{pointN1++; pointN2 = pointN1*pointN1;}
		if (crit<=300)	{pointN1--; pointN2 = pointN1*pointN1;}
		if ((crit>300)&&(crit<700))	break;
	}
	cout<<"crit = "<<crit<<"\n";
	return (int)(out_lenN*3.3219280948873626);
}



long int Functions::getfastN(std::string& str_N, std::string& str_der_N,std::string c_F, std::string c_energy, std::string c_Z2)	{


	std::string temp_str_N = "empty";
	std::string temp_str_der_N = "empty";
	
	
	mpfr_t F;
	mpfr_t Ci_max;
	mpfr_t absCi;
	mpfr_t abs_sum_N;
	mpc_t E;
	mpc_t Z2;
	mpc_t N;
	mpc_t der_N;
	mpc_t sum_der_N;
	mpc_t sum_N;
	mpc_t temp;
	mpc_t power;
	mpc_t Ci_3;
	mpc_t Ci_2;
	mpc_t Ci_1;
	mpc_t Ci;
	



	
	mpfr_init2(F, prec_bitN);
	mpfr_init2(abs_sum_N,prec_bitN);
	mpfr_init2(Ci_max, prec_bitN);
	mpfr_init2(absCi, prec_bitN);
	mpc_init2(E, prec_bitN);
	mpc_init2(Z2, prec_bitN);
	mpc_init2(N, prec_bitN);
	mpc_init2(der_N, prec_bitN);
	mpc_init2(sum_der_N, prec_bitN);
	mpc_init2(sum_N, prec_bitN);
	mpc_init2(temp, prec_bitN);
	mpc_init2(power, prec_bitN);
	mpc_init2(Ci_3, prec_bitN);
	mpc_init2(Ci_2, prec_bitN);
	mpc_init2(Ci_1, prec_bitN);
	mpc_init2(Ci, prec_bitN);
	
	

	mpfr_set_str(F, &c_F[0], 10, GMP_RNDN);
	mpc_set_str(E, &c_energy[0], 10, MPC_RNDNN);
	mpc_set_str(Z2, &c_Z2[0], 10, MPC_RNDNN);

	mpfr_neg(F,F,GMP_RNDN);	
	

	mpc_set_ui(Ci_3, 0, MPC_RNDNN);
	mpc_set_ui(Ci_2, 0, MPC_RNDNN);
	mpc_set_ui(Ci_1, 1, MPC_RNDNN);

		


	mpfr_set_ui(Ci_max,1,GMP_RNDN);
	mpc_set_ui(sum_N,1,MPC_RNDNN);	
	mpc_set_ui(sum_der_N,2*m + 1,MPC_RNDNN);
	mpc_div_ui(sum_der_N,sum_der_N,2,MPC_RNDNN);
		for (long int i=1;i<1000000000;i++)	{
			
			mpc_mul_fr(Ci,Ci_3,F,MPC_RNDNN); 
			mpc_mul_ui(Ci,Ci,pointN2,MPC_RNDNN); 
			mpc_mul_ui(Ci,Ci,pointN2,MPC_RNDNN);
			mpc_mul(temp,Ci_1,Z2,MPC_RNDNN);
			mpc_sub(Ci,Ci,temp,MPC_RNDNN);
			mpc_mul(temp,Ci_2,E,MPC_RNDNN); 
			mpc_mul_ui(temp,temp,2*pointN2,MPC_RNDNN);
			mpc_sub(Ci,Ci,temp,MPC_RNDNN);
			mpc_mul_ui(Ci,Ci,pointN2,MPC_RNDNN);
			mpc_div_ui(Ci,Ci,4,MPC_RNDNN);
			mpc_div_ui(Ci,Ci,i,MPC_RNDNN);
			mpc_div_ui(Ci,Ci,i + m,MPC_RNDNN);
			
			mpc_abs(absCi,Ci,GMP_RNDN);
			if(mpfr_greater_p(absCi,Ci_max) != 0)
				mpfr_set(Ci_max,absCi,GMP_RNDN);
			
			mpc_set_ui(temp,4,MPC_RNDNN);
			mpc_mul_ui(temp,temp,i,MPC_RNDNN);
			mpc_add_ui(temp,temp,2*m + 1,MPC_RNDNN);
			mpc_div_ui(temp,temp,2,MPC_RNDNN);
			mpc_mul(temp,temp,Ci,MPC_RNDNN);
			
			
			mpc_add(sum_der_N,sum_der_N,temp,MPC_RNDNN);
			mpc_add(sum_N,sum_N,Ci,MPC_RNDNN);
			
		
			
			if (i%100==0)	{
				std::cout<<"pointN1 = "<<pointN1<<" i= "<<i<<" "<<"sum_der_N = "<<mpc_get_str(10, out_lenN, sum_der_N, MPC_RNDNN)<<"  sum = "<<mpc_get_str(10, out_lenN, sum_N, MPC_RNDNN)<<"\n";
			
			str_N = mpc_get_str(10, out_lenN, sum_N, MPC_RNDNN);
			str_der_N = mpc_get_str(10, out_lenN, sum_der_N, MPC_RNDNN);
			
			if((str_N != temp_str_N)&&(temp_str_der_N != str_der_N))
			{
				temp_str_N = str_N;
				temp_str_der_N = temp_str_der_N;
			}
			else
				break;
			}
			
			
			
			
			mpc_set(Ci_3,Ci_2,MPC_RNDNN);
			mpc_set(Ci_2,Ci_1,MPC_RNDNN);
			mpc_set(Ci_1,Ci,MPC_RNDNN);
			
		}
		
		mpc_set_ui(power,1 + 2*m, MPC_RNDNN);
		mpc_div_ui(power,power,2, MPC_RNDNN);
		mpc_set_ui(temp,pointN1, MPC_RNDNN);
		mpc_pow(temp, temp, power, MPC_RNDNN);
		mpc_mul(N, sum_N, temp, MPC_RNDNN);
		
		mpc_mul(der_N,sum_der_N,temp, MPC_RNDNN);
		mpc_div_ui(der_N,der_N,pointN1, MPC_RNDNN);
		

	
//		mpc_set_ui(der_N,1 + 2*m, MPC_RNDNN);
		
		str_N = mpc_get_str(10, out_lenN, N, MPC_RNDNN);
		str_der_N = mpc_get_str(10, out_lenN, der_N, MPC_RNDNN);
		
		
		mpfr_log10(Ci_max,Ci_max,GMP_RNDN);
		mpc_abs(abs_sum_N,sum_N,GMP_RNDN);
		mpfr_log10(abs_sum_N,absCi,GMP_RNDN);
		long int poww = (long int)mpfr_get_d(Ci_max,GMP_RNDN) - (long int)mpfr_get_d(abs_sum_N,GMP_RNDN);
		
		mpfr_clear(F);
		mpfr_clear(Ci_max);
		mpfr_clear(abs_sum_N);
		mpfr_clear(absCi);
		mpc_clear(E);
		mpc_clear(Z2);
		mpc_clear(N);
		mpc_clear(der_N);
		mpc_clear(sum_der_N);
		mpc_clear(sum_N);
		mpc_clear(temp);
		mpc_clear(power);
		mpc_clear(Ci_3);
		mpc_clear(Ci_2);
		mpc_clear(Ci_1);
		mpc_clear(Ci);
	


	
	return poww;
}

/*
void Functions::getfastN(std::string& str_N, std::string& str_derN,std::string c_F, std::string c_energy, std::string c_Z2){
	
	mpfr_t F;
	mpc_t E;
	mpc_t Z2;
	mpc_t N;
	mpc_t derN;
	mpc_t temp;
//	char nameEG[MAX_LENGTH_ADDRESS];

	
	mpfr_init2(F, prec_bit);
	mpc_init2(E, prec_bit);
	mpc_init2(Z2, prec_bit);
	mpc_init2(N, prec_bit);
	mpc_init2(derN, prec_bit);
	mpc_init2(temp, prec_bit);
	
	
	
	mpfr_set_str(F, &c_F[0], 10, GMP_RNDN);
	mpc_set_str(E, &c_energy[0], 10, MPC_RNDNN);
	mpc_set_str(Z2, &c_Z2[0], 10, MPC_RNDNN);

	mpfr_neg(F,F,GMP_RNDN);	
	//Calculation of M - function and its derivative at the point mu = point1
	
	//Calculation of C[0]
	mpc_set_ui(C[0], 1, MPC_RNDNN);

	
	//Calculation of C[1]*mu*mu
	mpc_mul_si(C[1],Z2,-point2, MPC_RNDNN);	
	mpc_div_ui(C[1],C[1],4*(1 + m), MPC_RNDNN);

	//Calculation of C[3]*mu*mu*mu*mu*mu*mu
	mpc_mul(C[2],Z2,Z2,MPC_RNDNN); 
	mpc_mul_ui(temp,E,8*(1 + m),MPC_RNDNN); 
	mpc_sub(C[2],C[2],temp,MPC_RNDNN); 
	mpc_mul_ui(C[2],C[2],point4,MPC_RNDNN); 
	mpc_div_ui(C[2],C[2],32*(m + 1)*(m + 2),MPC_RNDNN);	
	

	//Calculation of C[i]*pow(mu,2*i)
	for (int i=3;i<nsum1+1;i++)	{
		
		mpc_mul_fr(C[i],C[i-3],F,MPC_RNDNN); 
		mpc_mul_ui(C[i],C[i],point4,MPC_RNDNN); 
		mpc_mul(temp,C[i-1],Z2,MPC_RNDNN);
		mpc_sub(C[i],C[i],temp,MPC_RNDNN);
		mpc_mul(temp,C[i-2],E,MPC_RNDNN); 
		mpc_mul_ui(temp,temp,2*point2,MPC_RNDNN);
		mpc_sub(C[i],C[i],temp,MPC_RNDNN);
		mpc_mul_ui(C[i],C[i],point2,MPC_RNDNN);
		mpc_div_ui(C[i],C[i],4*i*(i + m),MPC_RNDNN);
		
	}

	//Calculation of M(mu) 
	mpc_set_ui(N, 0, MPC_RNDNN);
	mpc_set_ui(derN, 0, MPC_RNDNN);
	for (int i=nsum1;i>-1;i--)		{
		mpc_add(N,N,C[i],MPC_RNDNN);
		mpc_mul_ui(temp,C[i],4*i+2*m+1, MPC_RNDNN);
		mpc_div_ui(temp,temp,2, MPC_RNDNN);
		mpc_add(derN,derN,temp,MPC_RNDNN);
	}
	
	mpc_set_ui(temp,1 + 2*m, MPC_RNDNN);
	mpc_div_ui(temp,temp,2, MPC_RNDNN);
	mpc_set_ui(C[0], point1, MPC_RNDNN);
	mpc_pow(temp, C[0],temp, MPC_RNDNN);
	
	mpc_mul(N,N,temp,MPC_RNDNN);
	
	mpc_div_ui(temp, temp, point1, MPC_RNDNN);
	
	mpc_mul(derN,derN,temp,MPC_RNDNN);
	
	str_N = mpc_get_str(10, out_len, N, MPC_RNDNN);
	str_derN = mpc_get_str(10, out_len, derN, MPC_RNDNN);
	
	
	mpfr_clear(F);
	mpc_clear(E);
	mpc_clear(Z2);
	mpc_clear(N);
	mpc_clear(derN);
	mpc_clear(temp);
	
	
	return;
}

*/


long int Functions::get_a(std::string& str_a, std::string& str_der_a,std::string c_F, std::string c_energy, std::string c_Z2)	{
	
	
	for (int i=0;i<nsum1+1;i++)	
		mpc_clear(C[i]);
	delete [] C;


	std::string line = "empty line";
	std::string line1 = "empty line";

	
	int sum = 0;
	int prec = 1;
	int _temp;
	int count = 0;

	
	mpfr_t temp_mpfr;
	mpc_t temp_mpc;
	


		
		
	precision = 10;
	nsum1 = 10;


		
	while (true)	{
		
		
		nsum1 += 10*sum;
		precision += 10*prec;
		
		prec_bit = (int)(precision*3.3219280948873626);
		
		C = new mpc_t[nsum1 + 1];
		for (int i=0;i<nsum1+1;i++)	
			mpc_init2(C[i], prec_bit);

//		out_len = 200;
		getfast_a(str_a, str_der_a, c_F,  c_energy,c_Z2);


		if ((line != str_a)||(line1 != str_der_a))	{
			line.clear();
			line1.clear();
			line = str_a;
			line1 = str_der_a;	
				cout<<out_len<<" "<<nsum1<<"  "<<precision<<"\n"<<str_a<<"\n"<<str_der_a<<"\n\n";
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

	
	
	

	
	return prec_bit + 10;
}




void Functions::getfast_a(std::string& str_a, std::string& str_der_a,std::string c_F, std::string c_energy, std::string c_Z2){
	
	
	mpfr_t F;
	mpfr_t sqrt_F;
	mpc_t E;
	mpc_t Z2;
	mpc_t a;
	mpc_t der_a;
	mpc_t temp;
	mpfr_t temp1;
	mpc_t temp_exp;
	mpc_t sum_a;
	mpc_t sum_an;
	
	mpc_t k1;
	mpc_t k2;
	mpc_t k3;

	
	mpc_init2(k1, prec_bit);
	mpc_init2(k2, prec_bit);
	mpc_init2(k3, prec_bit);
	mpfr_init2(F, prec_bit);
	mpfr_init2(sqrt_F, prec_bit);
	mpc_init2(E, prec_bit);
	mpc_init2(Z2, prec_bit);
	mpc_init2(a, prec_bit);
	mpc_init2(der_a, prec_bit);
	mpc_init2(temp, prec_bit);
	mpfr_init2(temp1, prec_bit);
	mpc_init2(temp_exp, prec_bit);
	mpc_init2(sum_a, prec_bit);
	mpc_init2(sum_an, prec_bit);
	
	mpfr_set_str(F, &c_F[0], 10, GMP_RNDN);
	mpc_set_str(E, &c_energy[0], 10, MPC_RNDNN);
	mpc_set_str(Z2, &c_Z2[0], 10, MPC_RNDNN);
	
	
	
	
	//calculation of sqrt(F)
	mpfr_sqrt(sqrt_F,F,GMP_RNDN);
	
	//calculation of constant k1
	mpc_mul(k1,E,E,MPC_RNDNN);
	mpc_div_fr(k1,k1,F,MPC_RNDNN);
	mpc_sub(k1,k1,Z2,MPC_RNDNN);
	mpc_mul_ui(k1,k1,4*point2,MPC_RNDNN);
	
	//calculation of constant k2
	mpc_set_ui_ui(k2,0,8,MPC_RNDNN);
	mpc_mul_ui(k2,k2,point1,MPC_RNDNN);
	mpc_mul_ui(k2,k2,point2,MPC_RNDNN);
	mpc_mul_fr(k2,k2,sqrt_F,MPC_RNDNN);
	mpc_neg(k2,k2,MPC_RNDNN);
	
	//calculation of constant k3
	mpc_set_ui_ui(k3,0,8*point1,MPC_RNDNN);
	mpc_mul(k3,k3,E,MPC_RNDNN);
	mpc_div_fr(k3,k3,sqrt_F,MPC_RNDNN);
	
	


		
	//Calculation of C[1]/point1
	mpc_set_ui(C[1],1, MPC_RNDNN);
	mpc_div_ui(C[1], C[1], point1, MPC_RNDNN);

	
	//Calculation of C[2]/point1^2
	mpc_mul(C[2],C[1],k1, MPC_RNDNN);	
	mpc_div(C[2],C[2],k2,MPC_RNDNN);

	
	//Calculation of C[3]/point1^3
	mpc_mul(temp,C[2],k1,MPC_RNDNN);
	mpc_mul(C[3],k3,C[1],MPC_RNDNN);
	mpc_add(C[3],C[3],temp,MPC_RNDNN);
	mpc_div(C[3],C[3],k2,MPC_RNDNN);
	
	
	//Calculation of C[i]/point1^i
	for (int i=4;i<nsum1+1;i++)	{
		
		mpc_mul_si(C[i],C[i-3],4*m*m -(5-2*i)*(5-2*i),MPC_RNDNN);
		mpc_mul_ui(temp,k3,i - 2,MPC_RNDNN);
		mpc_mul(temp,temp,C[i-2],MPC_RNDNN);
		mpc_add(C[i],C[i], temp,MPC_RNDNN);
		mpc_mul(temp, C[i-1],k1,MPC_RNDNN);
		mpc_add(C[i],C[i],temp,MPC_RNDNN);
		mpc_div(C[i],C[i],k2,MPC_RNDNN);
		mpc_div_ui(C[i],C[i],i - 1,MPC_RNDNN);

	}
	
	mpc_set_ui(sum_a, 0, MPC_RNDNN);
	mpc_set_ui(sum_an, 0, MPC_RNDNN);
	for (int i=nsum1;i>0;i--)	{	
		mpc_add(sum_a,sum_a,C[i],MPC_RNDNN);
		mpc_mul_ui(temp,C[i],i,MPC_RNDNN);
		mpc_add(sum_an,sum_an,temp,MPC_RNDNN);
	}
	
	
	mpfr_mul_ui(temp1,sqrt_F,point1,GMP_RNDN);
	mpfr_mul_ui(temp1,temp1,point2,GMP_RNDN);
	mpfr_div_ui(temp1,temp1,3,GMP_RNDN);
	mpc_mul_ui(temp,E,point1,MPC_RNDNN);
	mpc_div_fr(temp,temp,sqrt_F,MPC_RNDNN);
	mpc_add_fr(temp, temp,temp1,MPC_RNDNN);
	mpc_mul_i(temp,temp,1,MPC_RNDNN);
	mpc_exp(temp_exp,temp,MPC_RNDNN);
	

	mpc_mul(a,sum_a,temp_exp,MPC_RNDNN);
	
	
	
	mpfr_mul_ui(temp1,sqrt_F,point1,GMP_RNDN);
	mpfr_mul_ui(temp1,temp1,point2,GMP_RNDN);
	mpc_mul_ui(temp,E,point1,MPC_RNDNN);
	mpc_div_fr(temp,temp,sqrt_F,MPC_RNDNN);
	mpc_add_fr(temp, temp,temp1,MPC_RNDNN);
	mpc_mul_i(temp,temp,1,MPC_RNDNN);
	
	mpc_mul(der_a,sum_a,temp,MPC_RNDNN);
	mpc_sub(der_a,der_a,sum_an,MPC_RNDNN);
	mpc_mul(der_a,der_a,temp_exp,MPC_RNDNN);
	mpc_div_ui(der_a,der_a,point1,MPC_RNDNN);

	
	str_a = mpc_get_str(10, out_len, a, MPC_RNDNN);
	str_der_a = mpc_get_str(10, out_len, der_a, MPC_RNDNN);
	
//	cout<<"str_a = "<<str_a<<" "<<"str_der_a = "<<str_der_a<<"\n\n";
	
	mpc_clear(k1);
	mpc_clear(k2);
	mpc_clear(k3);
	mpfr_clear(F);
	mpfr_clear(sqrt_F);
	mpc_clear(E);
	mpc_clear(Z2);
	mpc_clear(a);
	mpc_clear(der_a);
	mpc_clear(temp);
	mpc_clear(temp_exp);
	mpfr_clear(temp1);
	mpc_clear(sum_a);
	mpc_clear(sum_an);
	
	return;
	
}


long int Functions::get_b(std::string& str_b, std::string& str_der_b,std::string c_F, std::string c_energy, std::string c_Z2)	{
	
	
	for (int i=0;i<nsum1+1;i++)	
		mpc_clear(C[i]);
	delete [] C;


	std::string line = "empty line";
	std::string line1 = "empty line";

	
	int sum = 0;
	int prec = 1;
	int _temp;
	int count = 0;

	
	mpfr_t temp_mpfr;
	mpc_t temp_mpc;
	


		
		
	precision = 10;
	nsum1 = 10;


		
	while (true)	{
		
		
		nsum1 += 10*sum;
		precision += 10*prec;
		
		prec_bit = (int)(precision*3.3219280948873626);
		
		C = new mpc_t[nsum1 + 1];
		for (int i=0;i<nsum1+1;i++)	
			mpc_init2(C[i], prec_bit);

//		out_len = 200;
		getfast_b(str_b, str_der_b, c_F,  c_energy,c_Z2);

		
		if ((line != str_b)||(line1 != str_der_b))	{
			line.clear();
			line1.clear();
			line = str_b;
			line1 = str_der_b;	
//			cout<<nsum1<<"  "<<precision<<"\n"<<str_N<<"\n"<<str_derN<<"\n\n";
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

	

	
	return prec_bit + 10;
}




void Functions::getfast_b(std::string& str_b, std::string& str_der_b,std::string c_F, std::string c_energy, std::string c_Z2){
	
	
	mpfr_t F;
	mpfr_t sqrt_F;
	mpc_t E;
	mpc_t Z2;
	mpc_t b;
	mpc_t der_b;
	mpc_t temp;
	mpfr_t temp1;
	mpc_t temp_exp;
	mpc_t sum_b;
	mpc_t sum_bn;
	
	mpc_t k1;
	mpc_t k2;
	mpc_t k3;

	
	mpc_init2(k1, prec_bit);
	mpc_init2(k2, prec_bit);
	mpc_init2(k3, prec_bit);
	mpfr_init2(F, prec_bit);
	mpfr_init2(sqrt_F, prec_bit);
	mpc_init2(E, prec_bit);
	mpc_init2(Z2, prec_bit);
	mpc_init2(b, prec_bit);
	mpc_init2(der_b, prec_bit);
	mpc_init2(temp, prec_bit);
	mpfr_init2(temp1, prec_bit);
	mpc_init2(temp_exp, prec_bit);
	mpc_init2(sum_b, prec_bit);
	mpc_init2(sum_bn, prec_bit);
	
	mpfr_set_str(F, &c_F[0], 10, GMP_RNDN);
	mpc_set_str(E, &c_energy[0], 10, MPC_RNDNN);
	mpc_set_str(Z2, &c_Z2[0], 10, MPC_RNDNN);
	
	
	
	
	//calculation of sqrt(F)
	mpfr_sqrt(sqrt_F,F,GMP_RNDN);
	
	//calculation of constant k1
	mpc_mul(k1,E,E,MPC_RNDNN);
	mpc_div_fr(k1,k1,F,MPC_RNDNN);
	mpc_sub(k1,k1,Z2,MPC_RNDNN);
	mpc_mul_ui(k1,k1,4*point2,MPC_RNDNN);
	
	//calculation of constant k2
	mpc_set_ui_ui(k2,0,8,MPC_RNDNN);
	mpc_mul_ui(k2,k2,point1,MPC_RNDNN);
	mpc_mul_ui(k2,k2,point2,MPC_RNDNN);
	mpc_mul_fr(k2,k2,sqrt_F,MPC_RNDNN);


	
	//calculation of constant k3
	mpc_set_ui_ui(k3,0,8*point1,MPC_RNDNN);
	mpc_mul(k3,k3,E,MPC_RNDNN);
	mpc_div_fr(k3,k3,sqrt_F,MPC_RNDNN);
	
	


		
	//Calculation of C[1]/point1
	mpc_set_ui(C[1],1, MPC_RNDNN);
	mpc_div_ui(C[1], C[1], point1, MPC_RNDNN);

	
	//Calculation of C[2]/point1^2
	mpc_mul(C[2],C[1],k1, MPC_RNDNN);	
	mpc_div(C[2],C[2],k2,MPC_RNDNN);

	
	//Calculation of C[3]/point1^3
	mpc_mul(C[3],C[2],k1,MPC_RNDNN);
	mpc_mul(temp,k3,C[1],MPC_RNDNN);
	mpc_sub(C[3],C[3],temp,MPC_RNDNN);
	mpc_div(C[3],C[3],k2,MPC_RNDNN);
	
	
	//Calculation of C[i]/point1^i
	for (int i=4;i<nsum1+1;i++)	{
		
		mpc_mul_si(C[i],C[i-3],4*m*m -(5-2*i)*(5-2*i),MPC_RNDNN);
		mpc_mul_ui(temp,k3,i - 2,MPC_RNDNN);
		mpc_mul(temp,temp,C[i-2],MPC_RNDNN);
		mpc_sub(C[i],C[i], temp,MPC_RNDNN);
		mpc_mul(temp, C[i-1],k1,MPC_RNDNN);
		mpc_add(C[i],C[i],temp,MPC_RNDNN);
		mpc_div(C[i],C[i],k2,MPC_RNDNN);
		mpc_div_ui(C[i],C[i],i - 1,MPC_RNDNN);

	}
	
	mpc_set_ui(sum_b, 0, MPC_RNDNN);
	mpc_set_ui(sum_bn, 0, MPC_RNDNN);
	for (int i=nsum1;i>0;i--)	{	
		mpc_add(sum_b,sum_b,C[i],MPC_RNDNN);
		mpc_mul_ui(temp,C[i],i,MPC_RNDNN);
		mpc_add(sum_bn,sum_bn,temp,MPC_RNDNN);
	}
	
	
	mpfr_mul_ui(temp1,sqrt_F,point1,GMP_RNDN);
	mpfr_mul_ui(temp1,temp1,point2,GMP_RNDN);
	mpfr_div_ui(temp1,temp1,3,GMP_RNDN);
	mpc_mul_ui(temp,E,point1,MPC_RNDNN);
	mpc_div_fr(temp,temp,sqrt_F,MPC_RNDNN);
	mpc_add_fr(temp, temp,temp1,MPC_RNDNN);
	mpc_mul_i(temp,temp,-1,MPC_RNDNN);
	mpc_exp(temp_exp,temp,MPC_RNDNN);
	

	mpc_mul(b,sum_b,temp_exp,MPC_RNDNN);

	
	
	
	mpfr_mul_ui(temp1,sqrt_F,point1,GMP_RNDN);
	mpfr_mul_ui(temp1,temp1,point2,GMP_RNDN);
	mpc_mul_ui(temp,E,point1,MPC_RNDNN);
	mpc_div_fr(temp,temp,sqrt_F,MPC_RNDNN);
	mpc_add_fr(temp, temp,temp1,MPC_RNDNN);
	mpc_mul_i(temp,temp,-1,MPC_RNDNN);
	
	mpc_mul(der_b,sum_b,temp,MPC_RNDNN);
	mpc_sub(der_b,der_b,sum_bn,MPC_RNDNN);
	mpc_mul(der_b,der_b,temp_exp,MPC_RNDNN);
	mpc_div_ui(der_b,der_b,point1,MPC_RNDNN);

	
	str_b = mpc_get_str(10, out_len, b, MPC_RNDNN);
	str_der_b = mpc_get_str(10, out_len, der_b, MPC_RNDNN);
	
	mpc_clear(k1);
	mpc_clear(k2);
	mpc_clear(k3);
	mpfr_clear(F);
	mpfr_clear(sqrt_F);
	mpc_clear(E);
	mpc_clear(Z2);
	mpc_clear(b);
	mpc_clear(der_b);
	mpc_clear(temp);
	mpc_clear(temp_exp);
	mpfr_clear(temp1);
	mpc_clear(sum_b);
	mpc_clear(sum_bn);
	
	return;
	
}





