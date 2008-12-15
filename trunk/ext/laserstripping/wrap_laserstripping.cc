#include "orbit_mpi.hh"

#include "wrap_las_strip_external_effects.hh"
#include "wrap_cpp_base_field_source.hh"

static PyMethodDef laserStrippingMethods[] = { {NULL,NULL} };

#ifdef __cplusplus
extern "C" {
#endif

  void initlaserstripping(){
    //create new module
    PyObject* module = Py_InitModule("laserstripping",laserStrippingMethods);
		wrap_laserstripping_las_strip_external_effects::initLasStripExternalEffects(module);
		wrap_laserstripping_cpp_base_field_source::initCppBaseFieldSource(module);

  }
	
#ifdef __cplusplus
}
#endif

