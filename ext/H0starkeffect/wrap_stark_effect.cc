#include "orbit_mpi.hh"

#include "wrap_functions.hh"


static PyMethodDef StarkeffectMethods[] = { {NULL,NULL} };

#ifdef __cplusplus
extern "C" {
#endif

  void initstarkeffect(){
    //create new module
    PyObject* module = Py_InitModule("starkeffect",StarkeffectMethods);
		wrap_functions::initFunctions(module);

		
  }
	
#ifdef __cplusplus
}
#endif
