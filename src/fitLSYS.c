#include "fitLSYS.h"

SEXP fitLSYS(SEXP C, SEXP rhs, SEXP b, SEXP active, SEXP RSS, SEXP maxIter, SEXP tolerance) {
    int p = Rf_ncols(C);
    R_xlen_t nActive = Rf_xlength(active);
    int nIter = Rf_asInteger(maxIter);
    double tol = Rf_asReal(tolerance);
    double *pC = REAL(C);
    double *prhs = REAL(rhs);
    b = PROTECT(Rf_duplicate(b));
    double *pb = REAL(b);
    int *pactive = INTEGER(active);
    double oldRSS = Rf_asReal(RSS);
    double newRSS = oldRSS;
    for (int iter = 0; iter < nIter; iter++) {
        oldRSS = newRSS;
        for (int j = 0; j < nActive; j++) { // loop over active predictors
            int k = pactive[j];
            double Ckk = pC[k * (p + 1)];
            double offset = 0;
            for (int m = 0; m < nActive; m++) {
                int n = pactive[m];
                offset += pC[p * k + n] * pb[n];
            }
            offset -= Ckk * pb[k];
            double rhs_offset = prhs[k] - offset;
            double sol = rhs_offset / Ckk;
            newRSS += (pow(sol, 2) - pow(pb[k], 2)) * Ckk - 2 * (sol - pb[k]) * rhs_offset;
            pb[k] = sol;
        }
        if (((oldRSS - newRSS) / oldRSS) < tol) {
            break;
        }
    }
    // Creating a list to return results
    SEXP list = PROTECT(Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(list, 0, b);
    SET_VECTOR_ELT(list, 1, Rf_ScalarReal(newRSS));
    UNPROTECT(2); // b, list
    return list;
}
