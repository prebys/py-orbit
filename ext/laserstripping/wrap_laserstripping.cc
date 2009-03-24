#include "orbit_mpi.hh"

#include "wrap_density_matrix.hh"
#include "wrap_schrodinger_equation.hh"
#include "wrap_ls_field_source.hh"
#include "wrap_hermite_gaussian_lf_mode.hh"
#include "wrap_las_field_container.hh"
#include "wrap_hydrogen_stark_param.hh"
#include "wrap_two_level_atom.hh"
#include "wrap_froissart_stora_lf.hh"

static PyMethodDef laserStrippingMethods[] = { {NULL,NULL} };

#ifdef __cplusplus
extern "C" {
#endif

  void initlaserstripping(){
    //create new module
    PyObject* module = Py_InitModule("laserstripping",laserStrippingMethods);
		wrap_density_matrix::initDensityMatrix(module);
		wrap_schrodinger_equation::initSchrodingerEquation(module);
		wrap_ls_field_source::initLSFieldSource(module);
		wrap_hermite_gaussian_lf_mode::initHermiteGaussianLFmode(module);
		wrap_las_field_container::initLaserFieldContainer(module);
		wrap_hydrogen_stark_param::initHydrogenStarkParam(module);
		wrap_two_level_atom::initTwoLevelAtom(module);
		wrap_froissart_stora_lf::initFroissartStoraLF(module);
		
  }
	
#ifdef __cplusplus
}
#endif