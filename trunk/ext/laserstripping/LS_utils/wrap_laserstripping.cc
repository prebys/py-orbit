#include "orbit_mpi.hh"

#include "wrap_density_matrix.hh"
#include "wrap_DM_noLaserField.hh"
#include "wrap_schrodinger_equation.hh"
#include "wrap_regular_grid_fs.hh"
#include "wrap_hermite_gaussian_lf_mode.hh"
#include "wrap_las_field_container.hh"
#include "wrap_hydrogen_stark_param.hh"
#include "wrap_two_level_atom.hh"
#include "wrap_froissart_stora_lf.hh"
#include "wrap_quad_em_field.hh"
#include "wrap_const_em_field.hh"
#include "wrap_print_ext_effects.hh"
#include "wrap_record_evolution.hh"
#include "wrap_walls.hh"
#include "wrap_fringe_field.hh"

static PyMethodDef laserStrippingMethods[] = { {NULL,NULL} };

#ifdef __cplusplus
extern "C" {
#endif

  void initlaserstripping(){
    //create new module
    PyObject* module = Py_InitModule("laserstripping",laserStrippingMethods);
		wrap_density_matrix::initDensityMatrix(module);
		wrap_DM_noLaserField::initDM_noLaserField(module);
		wrap_schrodinger_equation::initSchrodingerEquation(module);
		wrap_regular_grid_fs::initRegularGridFS(module);
		wrap_hermite_gaussian_lf_mode::initHermiteGaussianLFmode(module);
		wrap_las_field_container::initLaserFieldContainer(module);
		wrap_hydrogen_stark_param::initHydrogenStarkParam(module);
		wrap_two_level_atom::initTwoLevelAtom(module);
		wrap_froissart_stora_lf::initFroissartStoraLF(module);
		wrap_quad_em_field::initQuadEMfield(module);
		wrap_const_em_field::initConstEMfield(module);
		wrap_print_ext_effects::initPrintExtEffects(module);
		wrap_record_evolution::initRecordEvolution(module);
		wrap_walls::initWalls(module);
		wrap_fringe_field::initFringeField(module);
		
  }
	
#ifdef __cplusplus
}
#endif
