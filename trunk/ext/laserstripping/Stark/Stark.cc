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
#include "BufferStore.hh"

#include "Stark.hh"
#include "LorentzTransformationEM.hh"
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>


using namespace OrbitUtils;
using namespace std;

inline int convert3to1level(int n,int n1, int m){
	return 1+m*n+n1+(n*n*n-n)/3-m*abs(m-1)/2;
}




Stark::Stark(char* addressEG,int states)
{

	int rank_MPI,size_MPI;
	ORBIT_MPI_Comm_size(MPI_COMM_WORLD, &size_MPI);
	ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank_MPI);
	
	std::ifstream file;
	double F;
	double alpha=7.297352570e-3;
	int k,ks,fi,nn1,nn2,mm,nn1s,nn2s,mms;
	char nameEG[1024];
	string dump,dump2;
	int ff;
	
	st=states;
	const_relax = 4*alpha*alpha*alpha/3;
	levels=st*(1+st)*(1+2*st)/6;

	
	//this function measures parameter of input files (length, delta_F) using groung level file 000.txt

	if(rank_MPI == 0) {
		sprintf(nameEG,"%s000.txt",addressEG);	
		file.open(nameEG);	file>>F>>F>>F>>delta_F; file.close();
		file.open(nameEG);fi=0;	while(!file.eof())	{file>>F>>F>>F;fi++;} file.close();n_data=fi-1;
	}
	
  ORBIT_MPI_Bcast(&delta_F, 1, MPI_DOUBLE,0,MPI_COMM_WORLD);	
  ORBIT_MPI_Bcast(&n_data, 1, MPI_INT,0,MPI_COMM_WORLD);
	
  
  
	//allocating memory for dynamic massive of data that will be read fron files
  	field_thresh = new double[levels+1];
  	En = new double[levels+1];
  	Gamman = new double[levels+1];
  
	energy=new double*[levels+1];	
	for (int i=0;i<levels+1;i++)	energy[i]=new double[n_data+10];
	
	gamma_autoionization=new double*[levels+1];	
	for (int i=0;i<levels+1;i++)	gamma_autoionization[i]=new double[n_data+10];
	
	
	
	
	
	
	if(rank_MPI == 0) {
		sprintf(nameEG,"%stransitions.txt",addressEG);
		file.open(nameEG);	file>>ff>>ff>>ff>>dump>>ff>>ff>>ff; fi=0;	while(dump2!=dump)	{file>>dump2;fi++;} file.close(); order_trans = (fi-4)/3;
	}
	
	ORBIT_MPI_Bcast(&order_trans, 1, MPI_INT,0,MPI_COMM_WORLD);
	  
	
	//allocating memory for dynamic massive of data that will be read fron files	
	E_pow_n = new double[order_trans];
	
	dipole_transition_x=new double**[levels+1];	
	for (int i=0;i<levels+1;i++)	dipole_transition_x[i]=new double*[levels+1]; 
	for (int i=0;i<levels+1;i++) for (int j=0;j<levels+1;j++)	dipole_transition_x[i][j]=new double[order_trans];
	
	dipole_transition_y=new double**[levels+1];	
	for (int i=0;i<levels+1;i++)	dipole_transition_y[i]=new double*[levels+1]; 
	for (int i=0;i<levels+1;i++) for (int j=0;j<levels+1;j++)	dipole_transition_y[i][j]=new double[order_trans];
	
	dipole_transition_z=new double**[levels+1];	
	for (int i=0;i<levels+1;i++)	dipole_transition_z[i]=new double*[levels+1]; 
	for (int i=0;i<levels+1;i++) for (int j=0;j<levels+1;j++)	dipole_transition_z[i][j]=new double[order_trans];	
	
	int buff_index = 0;
	
  double* dump_arr = BufferStore::getBufferStore()->getFreeDoubleArr(buff_index,n_data+10);
  
  

	
	//this double loop reads dipole transitions anf fills spontaneous relaxation coefficients
	if(rank_MPI == 0) {
	
		std::ifstream file_in;
		sprintf(nameEG,"%stransitions.txt",addressEG);
		file_in.open(nameEG); 

		for(int n=1;n<states+1;n++){
			for(int m=-(n-1);m<(n-1)+1;m++){
				for(int n1=0;n1<n-abs(m)-1+1;n1++){
					
					for(int ns=1;ns<states+1;ns++){
						for(int ms=-(ns-1);ms<(ns-1)+1;ms++){
							for(int n1s=0;n1s<ns-abs(ms)-1+1;n1s++)		{
								
				
								k=convert3to1level(n,n1,m);
								ks=convert3to1level(ns,n1s,ms);
								
								while(true)	{
								file_in>>nn1>>nn2>>mm>>dump>>nn1s>>nn2s>>mms;
						
								if(nn1==n1&&nn2==n-n1-fabs(m)-1&&mm==m&&nn1s==n1s&&nn2s==ns-n1s-fabs(ms)-1&&mms==ms)	{
									
										for (fi=0; fi<order_trans; fi++) file_in>>dipole_transition_x[k][ks][fi];
										for (fi=0; fi<order_trans; fi++) file_in>>dipole_transition_y[k][ks][fi];
										for (fi=0; fi<order_trans; fi++) file_in>>dipole_transition_z[k][ks][fi];
																				
									break;
									
								}
								else for (fi=0; fi<order_trans; fi++) file_in>>F>>F>>F;
								}
									

								}

	
							}
						}
					}
				}
			}
		file_in.close();
		}
	

	
	
	
	
	for(int n=1;n<states+1;n++){
		for(int m=-(n-1);m<(n-1)+1;m++){
			for(int n1=0;n1<n-abs(m)-1+1;n1++){
				
				for(int ns=1;ns<states+1;ns++){
					for(int ms=-(ns-1);ms<(ns-1)+1;ms++){
						for(int n1s=0;n1s<ns-abs(ms)-1+1;n1s++)	{
							
							k=convert3to1level(n,n1,m);
							ks=convert3to1level(ns,n1s,ms);
													
						
							for (fi=0; fi<order_trans; fi++) ORBIT_MPI_Bcast(&dipole_transition_x[k][ks][fi], 1, MPI_DOUBLE,0,MPI_COMM_WORLD);
							for (fi=0; fi<order_trans; fi++) ORBIT_MPI_Bcast(&dipole_transition_y[k][ks][fi], 1, MPI_DOUBLE,0,MPI_COMM_WORLD);
							for (fi=0; fi<order_trans; fi++) ORBIT_MPI_Bcast(&dipole_transition_z[k][ks][fi], 1, MPI_DOUBLE,0,MPI_COMM_WORLD);

								
								
							}		
						}
					}
				}
			}
		}

	
	
	
	
	
	
	
	
	//this loop reads energies and autoionization coefficients
	if(rank_MPI == 0) {
		for(int n=1;n<states+1;n++){
			for(int m=-(n-1);m<(n-1)+1;m++){
				for(int n1=0;n1<n-abs(m)-1+1;n1++)	{
					std::ifstream file_in;
					k=convert3to1level(n,n1,m);		
					sprintf(nameEG,"%s%i%i%i.txt",addressEG,n1,n-n1-abs(m)-1,abs(m));		
					file_in.open(nameEG);for (fi=0; fi<n_data; fi++)	{file_in>>F>>energy[k][fi]>>gamma_autoionization[k][fi];}	file_in.close();
					
					file.open(nameEG);	while(!file.eof())	file>>field_thresh[k]>>F>>F; file.close();
					
				}
			}
		}
	}
	
	
	ORBIT_MPI_Bcast(field_thresh,levels+1 ,MPI_DOUBLE,0,MPI_COMM_WORLD);
	

	
	for(int n=1;n<states+1;n++){
		for(int m=-(n-1);m<(n-1)+1;m++){	
			for(int n1=0;n1<n-abs(m)-1+1;n1++)	{
				k=convert3to1level(n,n1,m);	
				
				
				for(fi=0; fi<n_data; fi++){
					dump_arr[fi] = energy[k][fi];
				}
				ORBIT_MPI_Bcast(dump_arr,n_data , MPI_DOUBLE,0,MPI_COMM_WORLD);
				for(fi=0; fi<n_data; fi++){
					energy[k][fi] = dump_arr[fi];
				}		
				
				for(fi=0; fi<n_data; fi++){
					dump_arr[fi] = gamma_autoionization[k][fi];
				}
				ORBIT_MPI_Bcast(dump_arr,n_data , MPI_DOUBLE,0,MPI_COMM_WORLD);
				for(fi=0; fi<n_data; fi++){
					gamma_autoionization[k][fi] = dump_arr[fi];
				}		
				
			}
		}
	}	
	


	
	BufferStore::getBufferStore()->setUnusedDoubleArr(buff_index);				
}


Stark::~Stark()	{

	
	for (int i=0;i<levels+1;i++)	delete	[]	energy[i];	delete	[]	energy;
	for (int i=0;i<levels+1;i++)	delete	[]	gamma_autoionization[i];	delete	[]	gamma_autoionization;
	for (int i=0;i<levels+1;i++) for (int j=0;j<levels+1;j++)	delete [] dipole_transition_x[i][j]; for (int i=0;i<levels+1;i++)	delete [] dipole_transition_x[i];	delete	[]	dipole_transition_x;
	for (int i=0;i<levels+1;i++) for (int j=0;j<levels+1;j++)	delete [] dipole_transition_y[i][j]; for (int i=0;i<levels+1;i++)	delete [] dipole_transition_y[i];	delete	[]	dipole_transition_y;
	for (int i=0;i<levels+1;i++) for (int j=0;j<levels+1;j++)	delete [] dipole_transition_z[i][j]; for (int i=0;i<levels+1;i++)	delete [] dipole_transition_z[i];	delete	[]	dipole_transition_z;
//	for (int i=0;i<levels+1;i++) for (int j=0;j<levels+1;j++)	delete [] gamma_spontaneous_relax[i][j]; for (int i=0;i<levels+1;i++)	delete [] gamma_spontaneous_relax[i];	delete	[]	gamma_spontaneous_relax;
	delete [] field_thresh;
	delete [] E_pow_n;

}





void Stark::getTransition(int k,int ks,tcomplex& mu_x,tcomplex& mu_y,tcomplex& mu_z){
	
	
	double mu_xd=0, mu_yd=0, mu_zd=0;
	
	E_pow_n[0] = 1;
	for (int fi=1; fi<order_trans; fi++)  E_pow_n[fi] *= Ez_stat;
	
	
	for (int fi=0; fi<order_trans; fi++) mu_xd += dipole_transition_x[k][ks][fi]*E_pow_n[fi];
	for (int fi=0; fi<order_trans; fi++) mu_yd += dipole_transition_y[k][ks][fi]*E_pow_n[fi];
	for (int fi=0; fi<order_trans; fi++) mu_zd += dipole_transition_z[k][ks][fi]*E_pow_n[fi];
			
	mu_x = tcomplex(mu_xd,0);
	mu_y = tcomplex(0,mu_yd);
	mu_z = tcomplex(mu_zd,0);
	
}



double Stark::getRelax(int k,int ks){
	
	double mu_xd, mu_yd, mu_zd;
	
	double absdE = fabs(En[k]-En[ks]);
	
	E_pow_n[0] = 1;
	for (int fi=1; fi<order_trans; fi++)  E_pow_n[fi] *= Ez_stat;
	
	for (int fi=0; fi<order_trans; fi++) mu_xd += dipole_transition_x[k][ks][fi]*E_pow_n[fi];
	for (int fi=0; fi<order_trans; fi++) mu_yd += dipole_transition_y[k][ks][fi]*E_pow_n[fi];
	for (int fi=0; fi<order_trans; fi++) mu_zd += dipole_transition_z[k][ks][fi]*E_pow_n[fi];
		
	return const_relax*absdE*absdE*absdE*(mu_xd*mu_xd+mu_yd*mu_yd+mu_zd*mu_zd);
}



double Stark::getRelaxTransition(int k,int ks,tcomplex& mu_x,tcomplex& mu_y,tcomplex& mu_z){
	
	double mu_xd, mu_yd, mu_zd;
	double absdE = fabs(En[k]-En[ks]);
	
	E_pow_n[0] = 1;
	for (int fi=1; fi<order_trans; fi++)  E_pow_n[fi] *= Ez_stat;
	
	
	for (int fi=0; fi<order_trans; fi++) mu_xd += dipole_transition_x[k][ks][fi]*E_pow_n[fi];
	for (int fi=0; fi<order_trans; fi++) mu_yd += dipole_transition_y[k][ks][fi]*E_pow_n[fi];
	for (int fi=0; fi<order_trans; fi++) mu_zd += dipole_transition_z[k][ks][fi]*E_pow_n[fi];
	
	
	mu_x = tcomplex(mu_xd,0);
	mu_y = tcomplex(0,mu_yd);
	mu_z = tcomplex(mu_zd,0);
		
	return const_relax*absdE*absdE*absdE*(mu_xd*mu_xd+mu_yd*mu_yd+mu_zd*mu_zd);
}




void	Stark::SetE(double E){	
	
	Ez_stat=E;	
	
	double c = Ez_stat/delta_F;
	iEz=(int)c;
	cEz=c-iEz;
	
	for(int k=1; k<levels+1;k++)	{
	if(Ez_stat>field_thresh[k])	{En[k] = 0; Gamman[k] = 0;}	
	else	{
	En[k] = energy[k][iEz]+cEz*(energy[k][iEz+1]-energy[k][iEz]);
	
	if(gamma_autoionization[k][iEz]<1e-30)
	Gamman[k] = 0;
	else	
	Gamman[k] = gamma_autoionization[k][iEz]*pow(gamma_autoionization[k][iEz+1]/gamma_autoionization[k][iEz],cEz);
	}
	}


}

int Stark::getStates() {return st;}



double	Stark::getStarkEnergy(
	double mass, int n1, int n2, int m, 
	double Ex,double Ey,double Ez,
	double Bx,double By,double Bz,
	double px,double py,double pz){



LorentzTransformationEM::transform(mass,px,py,pz,Ex,Ey,Ez,Bx,By,Bz);

SetE(sqrt(Ex*Ex+Ey*Ey+Ez*Ez)/5.14220642e011);	

return En[convert3to1level(n1+n2+abs(m)+1,n1,m)];

	}
	
	
	
	

