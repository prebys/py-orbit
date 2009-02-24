#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_utils.hh"
#include "wrap_ls_field_source.hh"

#include <iostream>
#include <string>

#include "LSFieldSource.hh"

using namespace OrbitUtils;
using namespace wrap_orbit_utils;

namespace wrap_ls_field_source{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif

	//---------------------------------------------------------
	//Python CppBaseFieldSource class definition
	//---------------------------------------------------------

	//constructor for python class wrapping CppBaseFieldSource instance
	//It never will be called directly
	static PyObject* LSFieldSource_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		return (PyObject *) self;
	}

  //initializator for python  CppBaseFieldSource class
  //this is implementation of the __init__ method
  static int LSFieldSource_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
	  
//	  self->cpp_obj = new LSFieldSource();
	 
	  double E_x;
	  double E_y;
	  double E_z;
	  double B_x;
	  double B_y;
	  double B_z;
	  int nVars = PyTuple_Size(args);
	     

//	 if (PyArg_ParseTuple(	args,"")) {self->cpp_obj = new LSFieldSource();}
//	 if(PyArg_ParseTuple(args,"dddddd:",&E_x,&E_y,&E_z,&B_x,&B_y,&B_z))	{ self->cpp_obj = new LSFieldSource(E_x,E_y,E_z,B_x,B_y,B_z);}
	  
	  if(nVars==0)	 if (!PyArg_ParseTuple(	args,"")) {} else {self->cpp_obj = new LSFieldSource();}
	  if(nVars==6)	 if (!PyArg_ParseTuple(args,"dddddd:",&E_x,&E_y,&E_z,&B_x,&B_y,&B_z)) 
	  {error("Parameters  (E_x,E_y,E_z,B_x,B_y,B_z) -are needed");} else {self->cpp_obj = new LSFieldSource(E_x,E_y,E_z,B_x,B_y,B_z);}

		  
		
    return 0;
  }
  


  //-----------------------------------------------------
  //destructor for python PyBaseFieldSource class (__del__ method).
  //-----------------------------------------------------
  static void LSFieldSource_del(pyORBIT_Object* self){
		//std::cerr<<"The LSFieldSource __del__ has been called!"<<std::endl;
		delete ((LSFieldSource*)self->cpp_obj);
		self->ob_type->tp_free((PyObject*)self);
  }

	// defenition of the methods of the python PyBaseFieldSource wrapper class
	// they will be vailable from python level
  static PyMethodDef LSFieldSourceClassMethods[] = {
    {NULL}
  };

	// defenition of the memebers of the python PyBaseFieldSource wrapper class
	// they will be vailable from python level
	static PyMemberDef LSFieldSourceClassMembers [] = {
		{NULL}
	};

	//new python PyBaseFieldSource wrapper type definition
	static PyTypeObject pyORBIT_LSFieldSource_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"LSFieldSource", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) LSFieldSource_del , /*tp_dealloc*/
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
		"The LSFieldSource python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		LSFieldSourceClassMethods, /* tp_methods */
		LSFieldSourceClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) LSFieldSource_init, /* tp_init */
		0, /* tp_alloc */
		LSFieldSource_new, /* tp_new */
	};	

	//--------------------------------------------------
	//Initialization function of the pyPyBaseFieldSource class
	//It will be called from Bunch wrapper initialization
	//--------------------------------------------------
  void initLSFieldSource(PyObject* module){
		if (PyType_Ready(&pyORBIT_LSFieldSource_Type) < 0) return;
		Py_INCREF(&pyORBIT_LSFieldSource_Type);
		PyModule_AddObject(module, "LSFieldSource", (PyObject *)&pyORBIT_LSFieldSource_Type);
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_utils_cpp_base_field_source
}
