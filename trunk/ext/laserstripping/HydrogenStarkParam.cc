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
#include "HydrogenStarkParam.hh"
#include <cmath>
#include <fstream>



double*** HydrogenStarkParam::dipole_transition_x = NULL;
double*** HydrogenStarkParam::dipole_transition_y = NULL;
double*** HydrogenStarkParam::dipole_transition_z = NULL;
double*** HydrogenStarkParam::gamma_spontaneous_relax = NULL;
double** HydrogenStarkParam::energy = NULL;
double** HydrogenStarkParam::gamma_autoionization = NULL;
 
double HydrogenStarkParam::delta_F = 0.0;
int HydrogenStarkParam::n_data = 0;
double HydrogenStarkParam::Ez_stat = 0.0;



void HydrogenStarkParam::GetDipoleTransition(int k,int ks,tcomplex& mu_x,tcomplex& mu_y,tcomplex& mu_z){

	int i=(int)(Ez_stat/delta_F);
	double c=(Ez_stat-i*delta_F)/delta_F;
	
mu_x=dipole_transition_x[k][ks][i]+c*(dipole_transition_x[k][ks][i+1]-dipole_transition_x[k][ks][i]);
mu_y=dipole_transition_y[k][ks][i]+c*(dipole_transition_y[k][ks][i+1]-dipole_transition_y[k][ks][i]);	mu_y*=tcomplex(0.,1.);
mu_z=dipole_transition_z[k][ks][i]+c*(dipole_transition_z[k][ks][i+1]-dipole_transition_z[k][ks][i]);



}



void HydrogenStarkParam::GetRelax(int k,int ks,double& relax){
	
	int i=(int)(Ez_stat/delta_F);
	double c=(Ez_stat-i*delta_F)/delta_F;
	
relax=fabs(gamma_spontaneous_relax[k][ks][i]+c*(gamma_spontaneous_relax[k][ks][i+1]-gamma_spontaneous_relax[k][ks][i]));

}

void	HydrogenStarkParam::SetE(double E){	Ez_stat=E;	}




void HydrogenStarkParam::GetEnergyAutoionization(int k,double& E, double& Gamma){
	
	int i=(int)(Ez_stat/delta_F);
	double c=(Ez_stat-i*delta_F)/delta_F;

	
	if ((gamma_autoionization[k][i]>1e-20)&&(gamma_autoionization[k][i+1]>1e-20))
	Gamma=gamma_autoionization[k][i]*pow(gamma_autoionization[k][i+1]/gamma_autoionization[k][i],c);
	else	Gamma=0;
	
	E=energy[k][i]+c*(energy[k][i+1]-energy[k][i]);
	

}


	void HydrogenStarkParam::ReadData(char* addressEG,int states)	{
	
	

	std::ifstream file;
	double F,alpha=7.297352570e-3;
	int k,ks,fi;
	char nameEG[1024];

	int levels=states*(1+states)*(1+2*states)/6;

	
//this function measures parameter of input files (length, delta_F) using groung level file 000.txt
sprintf(nameEG,"%s000.txt",addressEG);	
file.open(nameEG);file>>F>>F>>F>>delta_F; file.close();
file.open(nameEG);fi=0;	while(!file.eof())	{file>>F>>F>>F;fi++;} file.close();n_data=fi-1;


//allocating memory for dynamic massive of data that will be read fron files
energy=new double*[levels+1];	for (int i=0;i<levels+1;i++)	energy[i]=new double[n_data+10];
gamma_autoionization=new double*[levels+1];	for (int i=0;i<levels+1;i++)	gamma_autoionization[i]=new double[n_data+10];
gamma_spontaneous_relax=new double**[levels+1];	for (int i=0;i<levels+1;i++)	gamma_spontaneous_relax[i]=new double*[levels+1]; for (int i=0;i<levels+1;i++) for (int j=0;j<levels+1;j++)	gamma_spontaneous_relax[i][j]=new double[n_data+10];
dipole_transition_x=new double**[levels+1];	for (int i=0;i<levels+1;i++)	dipole_transition_x[i]=new double*[levels+1]; for (int i=0;i<levels+1;i++) for (int j=0;j<levels+1;j++)	dipole_transition_x[i][j]=new double[n_data+10];
dipole_transition_y=new double**[levels+1];	for (int i=0;i<levels+1;i++)	dipole_transition_y[i]=new double*[levels+1]; for (int i=0;i<levels+1;i++) for (int j=0;j<levels+1;j++)	dipole_transition_y[i][j]=new double[n_data+10];
dipole_transition_z=new double**[levels+1];	for (int i=0;i<levels+1;i++)	dipole_transition_z[i]=new double*[levels+1]; for (int i=0;i<levels+1;i++) for (int j=0;j<levels+1;j++)	dipole_transition_z[i][j]=new double[n_data+10];	





	//this loop reads energies and autoionization coefficients
	for(int n=1;n<states+1;n++)
	for(int m=-(n-1);m<(n-1)+1;m++)
	for(int n1=0;n1<n-abs(m)-1+1;n1++)	{
		
		k=1+m*n+n1+(n*n*n-n)/3-m*abs(m-1)/2;
		
sprintf(nameEG,"%s%i%i%i.txt",addressEG,n1,n-n1-abs(m)-1,abs(m));		
file.open(nameEG);for (fi=0; fi<n_data; fi++)	{file>>F>>energy[k][fi]>>gamma_autoionization[k][fi];}	file.close();

	}
	

	
	
	//this double loop reads dipole transitions anf fills spontaneous relaxation coefficients
	for(int n=1;n<states+1;n++)
	for(int m=-(n-1);m<(n-1)+1;m++)
	for(int n1=0;n1<n-abs(m)-1+1;n1++)

	for(int ns=1;ns<states+1;ns++)
	for(int ms=-(ns-1);ms<(ns-1)+1;ms++)
	for(int n1s=0;n1s<ns-abs(ms)-1+1;n1s++)		{
			
			k=1+m*n+n1+(n*n*n-n)/3-m*abs(m-1)/2;
			ks=1+ms*ns+n1s+(ns*ns*ns-ns)/3-ms*abs(ms-1)/2;
			
	sprintf(nameEG,"%s%i%i%i---%i%i%i.txt",addressEG,n1,n-n1-abs(m)-1,m,n1s,ns-n1s-abs(ms)-1,ms);
	file.open(nameEG);
	for (fi=0; fi<n_data; fi++)	{file>>F>>dipole_transition_x[k][ks][fi]>>dipole_transition_y[k][ks][fi]>>dipole_transition_z[k][ks][fi];
	
	// this condition assumes that probability of spontaneous (see next line) and indused tansition  between levels with the same principal quantum number n is sero
	if(ns==n)	{dipole_transition_x[k][ks][fi]=0;dipole_transition_y[k][ks][fi]=0;dipole_transition_z[k][ks][fi]=0;}
	
	//this loop fills spontaneous relaxation (transition) of atom in relative atomic units
	gamma_spontaneous_relax[k][ks][fi]=fabs((4*alpha*alpha*alpha/3)*pow(energy[k][fi]-energy[ks][fi],3)*(pow(dipole_transition_x[k][ks][fi],2)+pow(dipole_transition_y[k][ks][fi],2)+pow(dipole_transition_z[k][ks][fi],2)));
	
	}
	file.close();
	

//		cout<<setprecision(20)<<n1<<n-n1-abs(m)-1<<m<<"---"<<n1s<<ns-n1s-abs(ms)-1<<ms<<"  dx="<<dipole_transition_x[k][ks][1]<<"  dy="<<dipole_transition_y[k][ks][1]<<"  dz="<<dipole_transition_z[k][ks][1]<<"\n";
		}


	
	
	
	
}



