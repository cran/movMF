#include <R.h>
#include <R_ext/Rdynload.h>

#include "movMF.h"

#define CDEF(name, n)  {#name, (DL_FUNC) &name, n}

static const R_CMethodDef cMethods[] = {
    CDEF(mycfG, 4),
    CDEF(mycfP, 4),
    CDEF(my0F1, 5),
    CDEF(rW, 4),
    {NULL, NULL, 0}
};

void R_init_movMF(DllInfo *dll)
{
    R_useDynamicSymbols(dll, FALSE);
    R_registerRoutines(dll, cMethods, NULL, NULL, NULL);
}
