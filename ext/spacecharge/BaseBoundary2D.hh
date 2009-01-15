#ifndef SC_BASE_BOUNDARY_2D_H
#define SC_BASE_BOUNDARY_2D_H

//pyORBIT utils
#include "CppPyWrapper.hh"

#include "Grid2D.hh"

using namespace std;

/** 
The BaseBoundary2D class defines a boundary geometry
and calculates the potential created by charges on the boundary 
surface.
*/

class BaseBoundary2D: public OrbitUtils::CppPyWrapper
{
	public:
	
		/** Constructor */
		BaseBoundary2D(int nPoints, int nModes);
		
		/** Destructor */
		virtual ~BaseBoundary2D();
		
		/** Returns the maximal value of the grid in x-axis */   
		double getMaxX();
		
		/** Returns the minimal value of the grid in x-axis */   	
		double getMinX();
		
		/** Returns the maximal value of the grid in y-axis */ 	
		double getMaxY();
		
		/** Returns the minimal value of the grid in y-axis */ 	
		double getMinY();
		
		/** Adds potential from the boundary to the external grid */
		virtual void addBoundaryPotential(Grid2D* rhoGrid,Grid2D*  phiGrid);
		
		/** Returns the number of boundary points */
		int getNumberOfPoints();
		
		/** Returns the number of modes in the free-space field */
		int getNumberOfModes();
		
		/** Returns 0 if the boundary is not initialized */
		int isInitialized();
		
		/** Returns the x-coordinate of the boundary point with an index i */
		double getBoundaryPointX(int i);
		
		/** Returns the y-coordinate of the boundary point with an index i */
		double getBoundaryPointY(int i);
		
		/** Sets the boundary point with index to (x,y) */
	  void setBoundaryPoint(int index, double x, double y);
		
		/** Initializes all arrays related to boundary points  */
		void initializeBPs();
		
	protected:
	
		/** Calculates an inverse matrix  */
		void _gaussjinv(double **a, int n);
		
		/** Calculates all of the LSQM functions at one point */
		double* lsq_fuctions(double x, double y);	
		
	protected:
		
		//initialized - 1 if not it is 0
		int initialized_;
		
		//The limits of the boundary in [m]
		double xMin_,xMax_;
		double yMin_,yMax_; 
		
		//number of points on the boundary
		int nPoints_;
		
		
		//The arrays with boundary points 
		double* bArrX_;
		double* bArrY_;
		
		//-------------------------------------------
		//data relaited to the free-space solutions
		//-------------------------------------------
		//LSQM array and matrix 
		//2*nModes_+1 - number of functions in LSQM 
		
		//number of modes in the free-space field
		int nModes_;	
		
		double BPrnorm_;
		
		double* bnd_phi_arr_;
		
		double* lsq_func_vctr_;
		double* lsq_coef_vctr_;
		
		int* ipiv_tmp;
		int* indxr_tmp;
		int* indxc_tmp;	
		
		double** LSQ_matrix_;
		double** tmp_matrix_;
};

#endif
//endif for SC_BASE_BOUNDARY_2D_H
