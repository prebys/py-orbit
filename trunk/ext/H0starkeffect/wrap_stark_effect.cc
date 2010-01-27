#include "orbit_mpi.hh"

#include "wrap_functions.hh"
#include "wrap_wave_function.hh"
#include "wrap_funct.hh"

static PyMethodDef StarkeffectMethods[] = { {NULL,NULL} };

#ifdef __cplusplus
extern "C" {
#endif

  void initstarkeffect(){
    //create new module
    PyObject* module = Py_InitModule("starkeffect",StarkeffectMethods);
		wrap_functions::initFunctions(module);
		wrap_wave_function::initWaveFunction(module);
		wrap_funct::initFunct(module);

		
  }
	
#ifdef __cplusplus
}
#endif
