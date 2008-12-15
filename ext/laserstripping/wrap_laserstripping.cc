#include "orbit_mpi.hh"

#include "wrap_las_strip_external_effects.hh"
#include "wrap_cpp_base_field_source.hh"

static PyMethodDef trackerrk4Methods[] = { {NULL,NULL} };

#ifdef __cplusplus
extern "C" {
#endif

  void initlaserstripping(){
    //create new module
    PyObject* module = Py_InitModule("laserstripping",trackerrk4Methods);
		wrap_trackerrk4_las_strip_external_effects::initLasStripExternalEffects(module);
		wrap_trackerrk4_cpp_base_field_source::initCppBaseFieldSource(module);

  }
	
#ifdef __cplusplus
}
#endif

