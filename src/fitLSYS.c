#include "fitLSYS.h"

SEXP fitLSYS(SEXP nCol, SEXP nActive, SEXP C, SEXP rhs, SEXP b, SEXP active, SEXP RSS, SEXP maxIter, SEXP tolerance) {
    int p;
    int m;
    int n;
    int q;
    int nIter = 0;
    int j = 0;
    int iter = 0;
    int k = 0;
    int *pactive;
    double Ckk;
    double offset;
    double rhs_offset;
    double tol;
    double sol;
    double *pRSS;
    double *pC;
    double *prhs;
    double *pb;
    double RSS0;
    SEXP list;
    p = Rf_asInteger(nCol);
    q = Rf_asInteger(nActive);
    nIter = Rf_asInteger(maxIter);
    tol = Rf_asReal(tolerance);
    PROTECT(C = Rf_coerceVector(C, REALSXP));
    pC = REAL(C);
    PROTECT(rhs = Rf_coerceVector(rhs, REALSXP));
    prhs = REAL(rhs);
    PROTECT(b = Rf_coerceVector(b, REALSXP));
    pb = REAL(b);
    PROTECT(active = Rf_coerceVector(active, INTSXP));
    pactive = INTEGER(active);
    PROTECT(RSS = Rf_coerceVector(RSS, REALSXP));
    pRSS = REAL(RSS);
    RSS0 = pRSS[0] + 0.0;
    while (iter < nIter) {
        iter += 1;
        RSS0 = pRSS[0] + 0.0;
        for (j = 0; j < q; j++) { // loop over active predictors
            k = pactive[j];
            Ckk = pC[k * (p + 1)];
            offset = 0.0;
            for (m = 0; m < q; m++) {
                n = pactive[m];
                offset += pC[p * k + n] * pb[n];
            }
            offset -= Ckk * pb[k];
            rhs_offset = prhs[k] - offset;
            sol = rhs_offset / Ckk;
            pRSS[0] += (pow(sol, 2) - pow(pb[k], 2)) * Ckk -2 * (sol - pb[k]) * rhs_offset;
            pb[k] = sol;
        }
        if (((RSS0 - pRSS[0]) / RSS0) < tol) {
            break;
        }
    }
    // Creating a list to return results
    PROTECT(list = Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(list, 0, b);
    SET_VECTOR_ELT(list, 1, RSS);
    UNPROTECT(6);
    return list;
 }
