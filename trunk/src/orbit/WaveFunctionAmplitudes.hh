//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//   MyAttr.hh
//
// AUTHOR
//    T. Gorlov
//
// CREATED
//    07/14/2005
//
// DESCRIPTION
//    A subclass of the particle attributes class. This is container
//    for a macrosize of macro-particles in the bunch.
//
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
//
// INCLUDE FILES
//
///////////////////////////////////////////////////////////////////////////
#ifndef MYATTR_H
#define MYATTR_H

///////////////////////////////////////////////////////////////////////////
//
// INCLUDE FILES
//
///////////////////////////////////////////////////////////////////////////
#include <string>

///////////////////////////////////////////////////////////////////////////
//
// CLASS NAME
//     WaveFunctionAmplitudes
//
///////////////////////////////////////////////////////////////////////////
#include "ParticleAttributes.hh"

class WaveFunctionAmplitudes : public ParticleAttributes
{
public:
  //--------------------------------------
  //the public methods of the ParticleMacroSize class
  //--------------------------------------
	
	/** Constructor. This Attribute describe complex coefficients of Wave functions.
	  * The defailt size is 400. 
		*/
	WaveFunctionAmplitudes(Bunch* bunch);
	
	/** This Attribute describe complex coefficients of Wave functions.
	  * User can specify the number of variables that he wants to reserve.
		*/
	WaveFunctionAmplitudes(Bunch* bunch, int size_in);
	
  ~WaveFunctionAmplitudes();
  
  double& Re0(int particle_index);
  double& Im0(int particle_index);
	
  int getAttSize();
	
private:
	int size;
	
	
};

///////////////////////////////////////////////////////////////////////////
//
// END OF FILE
//
///////////////////////////////////////////////////////////////////////////

#endif
