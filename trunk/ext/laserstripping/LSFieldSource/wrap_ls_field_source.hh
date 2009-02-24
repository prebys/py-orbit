#ifndef WRAP_LS_FIELD_SOURCE_HH_
#define WRAP_LS_FIELD_SOURCE_HH_

#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_ls_field_source{
    void initLSFieldSource(PyObject* module);
  }

#ifdef __cplusplus
}
#endif


#endif /*WRAP_LS_FIELD_SOURCE_HH_*/

