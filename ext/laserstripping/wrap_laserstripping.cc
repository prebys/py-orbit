#include "orbit_mpi.hh"

#include "wrap_las_strip_external_effects.hh"
#include "wrap_ls_field_source.hh"
#include "wrap_hermite_gaussian_lf_mode.hh"
#include "wrap_las_field_container.hh"
#include "wrap_hydrogen_stark_param.hh"

static PyMethodDef laserStrippingMethods[] = { {NULL,NULL} };

#ifdef __cplusplus
extern "C" {
#endif

  void initlaserstripping(){
    //create new module
    PyObject* module = Py_InitModule("laserstripping",laserStrippingMethods);
		wrap_laserstripping_las_strip_external_effects::initLasStripExternalEffects(module);
		wrap_laserstripping_ls_field_source::initLSFieldSource(module);
		wrap_hermite_gaussian_lf_mode::initHermiteGaussianLFmode(module);
		wrap_las_field_container::initLaserFieldContainer(module);
		wrap_hydrogen_stark_param::initHydrogenStarkParam(module);
		
  }
	
#ifdef __cplusplus
}
#endif
