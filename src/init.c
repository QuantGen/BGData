#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <stdlib.h> // for NULL

extern SEXP summarize(SEXP);
extern SEXP rayOLS(SEXP, SEXP);
extern SEXP preprocess(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef callMethods[] = {
    {"summarize", (DL_FUNC) &summarize, 1},
    {"rayOLS", (DL_FUNC) &rayOLS, 2},
    {"preprocess", (DL_FUNC) &preprocess, 5},
    {NULL, NULL, 0}
};

void R_init_BGData(DllInfo *dll) {
    R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
