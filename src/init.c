#include "summarize.h"
#include "rayOLS.h"
#include "preprocess.h"
#include "fitLSYS.h"

#include <R_ext/Rdynload.h>

static const R_CallMethodDef callMethods[] = {
    {"summarize", (DL_FUNC) &summarize, 1},
    {"rayOLS", (DL_FUNC) &rayOLS, 2},
    {"preprocess", (DL_FUNC) &preprocess, 4},
    {"fitLSYS", (DL_FUNC) &fitLSYS, 9},
    {NULL, NULL, 0}
};

void R_init_BGData(DllInfo *dll) {
    R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
