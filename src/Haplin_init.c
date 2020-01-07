#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

// #define DISTS_API extern "C"
/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void JohnsonFitR(double *xnp,double *xmp,double *x0p,double *xkp,double *xpp,
	double *gammap,double *deltap,double *xip,double *lambdap,int *typep);
extern void pJohnsonR(double *xp,double *gammap,double *deltap,double *xip,
	double *lambdap,int *typep,int *Np,double *valuep);
extern void uJohnsonR(double *xp,double *gammap,double *deltap,double *xip,
	double *lambdap,int *typep,int *Np,double *valuep);

static const R_CMethodDef CEntries[] = {
    {"JohnsonFitR", (DL_FUNC) &JohnsonFitR, 10},
    {"pJohnsonR",   (DL_FUNC) &pJohnsonR,    8},
    {"uJohnsonR",   (DL_FUNC) &uJohnsonR,    8},
    {NULL, NULL, 0}
};

void R_init_Haplin(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

