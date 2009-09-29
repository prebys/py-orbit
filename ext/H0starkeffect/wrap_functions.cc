#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_utils.hh"
#include "wrap_functions.hh"

#include <iostream>
#include <string>
#include <cmath>


#include "Functions.hh"

//using namespace OrbitUtils;
using namespace wrap_orbit_utils;

namespace wrap_Functions{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif

	//---------------------------------------------------------
	//Python CppBaseFieldSource class definition
	//---------------------------------------------------------

	//constructor for python class wrapping CppBaseFieldSource instance
	//It never will be called directly
	static PyObject* Functions_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		return (PyObject *) self;
	}

  //initializator for python  CppBaseFieldSource class
  //this is implementation of the __init__ method
  static int Functions_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){

	  int n1,n2,m;
	  
	  int n_steps;
      double step;
      int nsum;


		 if(!PyArg_ParseTuple(	args,"iiiidi:",&n1, &n2, &m, &n_steps, &step, &nsum)){
			 		          error("Functions(n1, n2, m, address) - params. are needed");
			 			 		        }  
		 else	{
			 

		 self->cpp_obj =  new  Functions(n1,n2,m, n_steps, step, nsum);
		 ((Functions*) self->cpp_obj)->setPyWrapper((PyObject*) self);
		 }
	

    return 0;
    
    
  }
  
  
  
  
  static PyObject* Functions_getMfunctionModulus(PyObject *self, PyObject *args){
	  Functions* cpp_Functions = (Functions*)((pyORBIT_Object*) self)->cpp_obj;
  				       
 
	   double E;
	   double Gamma;
	   double reZ;
	   double imZ;
       
       double val;
       
       

           //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
           if(!PyArg_ParseTuple(	args,"dddd:",&E,&Gamma,&reZ, &imZ))
             error(" getMfunctionModulus(k, step, nsum) - parameters are needed");
           else
           val = cpp_Functions->getMfunctionModulus(E, Gamma, reZ, imZ);
           return Py_BuildValue("d",val);
  }
  
  
  static PyObject* Functions_setupE(PyObject *self, PyObject *args){
	  Functions* cpp_Functions = (Functions*)((pyORBIT_Object*) self)->cpp_obj;			       
       double E;
           //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
           if(!PyArg_ParseTuple(	args,"d:",&E))
             error(" SetupsetupE(E - parameter is needed");
           else 	  
        	   cpp_Functions->setupE(E);
         
  		    Py_INCREF(Py_None);
  		    return Py_None;
  }
  
  
  	
  

  

  //-----------------------------------------------------
  //destructor for python PyBaseFieldSource class (__del__ method).
  //-----------------------------------------------------
  static void Functions_del(pyORBIT_Object* self){
		//std::cerr<<"The CppBaseFieldSource __del__ has been called!"<<std::endl;
		delete ((Functions*)self->cpp_obj);
		self->ob_type->tp_free((PyObject*)self);
  }

	// defenition of the methods of the python PyBaseFieldSource wrapper class
	// they will be vailable from python level
  static PyMethodDef FunctionsClassMethods[] = {
		    { "getMfunctionModulus",  Functions_getMfunctionModulus,         		METH_VARARGS,"gets M-function Modulus"},
			{ "setupE",         Functions_setupE,				         			METH_VARARGS,"Sets the parameter of electric field."},
/*			{ "setLaserPower",         HermiteGaussianLFmode_setLaserPower,         METH_VARARGS,"Sets or returns the name of effects."},
			{ "setLaser_lambda",         HermiteGaussianLFmode_setLaser_lambda,         METH_VARARGS,"Sets or returns the name of effects."},
			{ "getLaserHalfAngle",         HermiteGaussianLFmode_getLaserHalfAngle,         METH_VARARGS,"Sets or returns the name of effects."},
			{ "getLaserPower",         HermiteGaussianLFmode_getLaserPower,         METH_VARARGS,"Sets or returns the name of effects."},
			{ "getLaser_lambda",         HermiteGaussianLFmode_getLaser_lambda,         METH_VARARGS,"Sets or returns the name of effects."},
			{ "setLaserFieldOrientation",  HermiteGaussianLFmode_setLaserFieldOrientation,         METH_VARARGS,"Sets or returns the name of effects."},
*/
    {NULL}
  };

	// defenition of the memebers of the python PyBaseFieldSource wrapper class
	// they will be vailable from python level
	static PyMemberDef FunctionsClassMembers [] = {
		{NULL}
	};

	//new python PyBaseFieldSource wrapper type definition
	static PyTypeObject pyORBIT_Functions_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"Functions", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) Functions_del , /*tp_dealloc*/
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
		"The Functions python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		FunctionsClassMethods, /* tp_methods */
		FunctionsClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) Functions_init, /* tp_init */
		0, /* tp_alloc */
		Functions_new, /* tp_new */
	};	

	//--------------------------------------------------
	//Initialization function of the pyPyBaseFieldSource class
	//It will be called from Bunch wrapper initialization
	//--------------------------------------------------
  void initFunctions(PyObject* module){
		if (PyType_Ready(&pyORBIT_Functions_Type) < 0) return;
		Py_INCREF(&pyORBIT_Functions_Type);
		PyModule_AddObject(module, "Functions", (PyObject *)&pyORBIT_Functions_Type);
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_utils_cpp_base_field_source
}
