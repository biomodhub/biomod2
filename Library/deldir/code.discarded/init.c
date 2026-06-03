#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(binsrt)(double *, double *, double *, int *, int *,
                             int *, double *, double *, int *, int *);
extern void F77_NAME(intri)(double *, double *, double *, double *, int *,
                            int *);
extern void F77_NAME(master)(double *, double *, double *, int *, int *,
                             int *, int *, double *, double *, int *,
                             double *, double *, int *, double *, int *);
extern void F77_NAME(mnnd)(double *, double *, int *, double *, double *);
extern void F77_NAME(cross)(double *, double *, int *, double *);

static const R_FortranMethodDef FortranEntries[] = {
    {"binsrt", (DL_FUNC) &F77_NAME(binsrt), 10},
    {"intri",  (DL_FUNC) &F77_NAME(intri),   6},
    {"master", (DL_FUNC) &F77_NAME(master), 15},
    {"mnnd",   (DL_FUNC) &F77_NAME(mnnd),    5},
    {"cross",  (DL_FUNC) &F77_NAME(cross),   4},
    {NULL, NULL, 0}
};

void R_init_deldir(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

