#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <stdlib.h> // for NULL

extern SEXP summarize(SEXP);

static const R_CallMethodDef callMethods[] = {
    {"C_summarize", (DL_FUNC) &summarize, 1},
    {NULL, NULL, 0}
};

void R_init_BGData(DllInfo *dll) {
    R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
