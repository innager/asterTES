#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP A0A1(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP A0A1locR(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP PuieR(SEXP, SEXP, SEXP, SEXP);
extern SEXP PumsR(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP PuUnR(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP Pxy1R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"A0A1",     (DL_FUNC) &A0A1,     13},
    {"A0A1locR", (DL_FUNC) &A0A1locR, 13},
    {"PuieR",    (DL_FUNC) &PuieR,     4},
    {"PumsR",    (DL_FUNC) &PumsR,     6},
    {"PuUnR",    (DL_FUNC) &PuUnR,     6},
    {"Pxy1R",    (DL_FUNC) &Pxy1R,     8},
    {NULL, NULL, 0}
};

void R_init_asterTES(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
