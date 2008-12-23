#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_utils.hh"
#include "wrap_gaussian_laser_field.hh"

#include <iostream>
#include <string>

#include "GaussianLaserField.hh"

using namespace OrbitUtils;
using namespace wrap_orbit_utils;

namespace wrap_gasussian_laser_field{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif

	//---------------------------------------------------------
	//Python CppBaseFieldSource class definition
	//---------------------------------------------------------

	//constructor for python class wrapping CppBaseFieldSource instance
	//It never will be called directly
	static PyObject* GaussianLaserField_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		return (PyObject *) self;
	}

  //initializator for python  CppBaseFieldSource class
  //this is implementation of the __init__ method
  static int GaussianLaserField_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
		self->cpp_obj = new GaussianLaserField();	//constructor with parameters is nedeed
    return 0;
  }
  
  
  
  /*
 	void	LasStripExternalEffects::setLaserHalfAngle(double a);
 	double	LasStripExternalEffects::getLaserHalfAngle();
 	void	LasStripExternalEffects::setLaserPower(double a);
 	double	LasStripExternalEffects::getLaserPower();
 	void	LasStripExternalEffects::setLaser_lambda(double a);
 	double	LasStripExternalEffects::getLaser_lambda();
   */
   
 static PyObject* GaussianLaserField_setLaserHalfAngle(PyObject *self, PyObject *args){
  		GaussianLaserField* LaserField = (GaussianLaserField*)((pyORBIT_Object*) self)->cpp_obj;
  		  double half_Angle;	       
  		        if(!PyArg_ParseTuple(	args,"d:",&half_Angle)){
  		          error("LaserExternalEfects - setLaserHalfAngle(half_Angle) - param. is needed");
  		        }
  		        else LaserField->setLaserHalfAngle(half_Angle);
  		      
  	   
  		    Py_INCREF(Py_None);
  		    return Py_None;	  
   }		

 static PyObject* GaussianLaserField_setLaserPower(PyObject *self, PyObject *args){
	 GaussianLaserField* LaserField = (GaussianLaserField*)((pyORBIT_Object*) self)->cpp_obj;
 		  double LaserPower;	       
 		        if(!PyArg_ParseTuple(	args,"d:",&LaserPower)){
 		          error("LaserExternalEfects - setLaserPower(LaserPower) - param. is needed");
 		        }
 		        else LaserField->setLaserPower(LaserPower);
 		      
 	   
 		    Py_INCREF(Py_None);
 		    return Py_None;	  
 }

 static PyObject* GaussianLaserField_setLaser_lambda(PyObject *self, PyObject *args){
	 GaussianLaserField* LaserField = (GaussianLaserField*)((pyORBIT_Object*) self)->cpp_obj;
 		  double Laser_lambda;	       
 		        if(!PyArg_ParseTuple(	args,"d:",&Laser_lambda)){
 		          error("LaserExternalEfects - setLaser_lambda(Laser_lambda) - param. is needed");
 		        }
 		        else LaserField->setLaserPower(Laser_lambda);
 		      
 	   
 		    Py_INCREF(Py_None);
 		    return Py_None;	  
 }




 static PyObject* GaussianLaserField_getLaserHalfAngle(PyObject *self, PyObject *args){
	 GaussianLaserField* LaserField = (GaussianLaserField*)((pyORBIT_Object*) self)->cpp_obj;
 				       
 		  return Py_BuildValue("d",LaserField->getLaserHalfAngle()); 
 }		
     

 static PyObject* GaussianLaserField_getLaserPower(PyObject *self, PyObject *args){
	 GaussianLaserField* LaserField = (GaussianLaserField*)((pyORBIT_Object*) self)->cpp_obj;
 				       
 		  return Py_BuildValue("d",LaserField->getLaserPower()); 
 }		    

 static PyObject* GaussianLaserField_getLaser_lambda(PyObject *self, PyObject *args){
	 GaussianLaserField* LaserField = (GaussianLaserField*)((pyORBIT_Object*) self)->cpp_obj;
 				       
 		  return Py_BuildValue("d",LaserField->getLaser_omega()); 
 }		    


  
  
  
  
  
  

  //-----------------------------------------------------
  //destructor for python PyBaseFieldSource class (__del__ method).
  //-----------------------------------------------------
  static void GaussianLaserField_del(pyORBIT_Object* self){
		//std::cerr<<"The CppBaseFieldSource __del__ has been called!"<<std::endl;
		delete ((GaussianLaserField*)self->cpp_obj);
		self->ob_type->tp_free((PyObject*)self);
  }

	// defenition of the methods of the python PyBaseFieldSource wrapper class
	// they will be vailable from python level
  static PyMethodDef GaussianLaserFieldClassMethods[] = {
			{ "setLaserHalfAngle",         GaussianLaserField_setLaserHalfAngle,         METH_VARARGS,"Sets or returns the name of effects."},
			{ "setLaserPower",         GaussianLaserField_setLaserPower,         METH_VARARGS,"Sets or returns the name of effects."},
			{ "setLaser_lambda",         GaussianLaserField_setLaser_lambda,         METH_VARARGS,"Sets or returns the name of effects."},
			{ "getLaserHalfAngle",         GaussianLaserField_getLaserHalfAngle,         METH_VARARGS,"Sets or returns the name of effects."},
			{ "getLaserPower",         GaussianLaserField_getLaserPower,         METH_VARARGS,"Sets or returns the name of effects."},
			{ "getLaser_lambda",         GaussianLaserField_getLaser_lambda,         METH_VARARGS,"Sets or returns the name of effects."},
    {NULL}
  };

	// defenition of the memebers of the python PyBaseFieldSource wrapper class
	// they will be vailable from python level
	static PyMemberDef GaussianLaserFieldClassMembers [] = {
		{NULL}
	};

	//new python PyBaseFieldSource wrapper type definition
	static PyTypeObject pyORBIT_GaussianLaserField_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"CppBaseFieldSource", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) GaussianLaserField_del , /*tp_dealloc*/
		0, /*tp_print*/
		0, /*tp_getattr*/
		0, /*tp_setattr*/
		0, /*tp_compare*/
		0, /*tp_repr*/
		0, /*tp_as_number*/
		0, /*tp_as_sequence*/
		0, /*tp_as_mapping*/
		0, /*tp_hash */
		0, /*tp_call*/
		0, /*tp_str*/
		0, /*tp_getattro*/
		0, /*tp_setattro*/
		0, /*tp_as_buffer*/
		Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
		"The GaussianLaserField python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		GaussianLaserFieldClassMethods, /* tp_methods */
		GaussianLaserFieldClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) GaussianLaserField_init, /* tp_init */
		0, /* tp_alloc */
		GaussianLaserField_new, /* tp_new */
	};	

	//--------------------------------------------------
	//Initialization function of the pyPyBaseFieldSource class
	//It will be called from Bunch wrapper initialization
	//--------------------------------------------------
  void initGaussianLaserField(PyObject* module){
		if (PyType_Ready(&pyORBIT_GaussianLaserField_Type) < 0) return;
		Py_INCREF(&pyORBIT_GaussianLaserField_Type);
		PyModule_AddObject(module, "GaussianLaserField", (PyObject *)&pyORBIT_GaussianLaserField_Type);
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_utils_cpp_base_field_source
}
