#define R_NO_REMAP

#ifdef _OPENMP
#include <omp.h>
#endif
#include <R.h>
#include <Rinternals.h>
#include <stddef.h>

void preprocess_int(int *in, int nrows, int ncols, double *out, int center, double *centers, int computeCenters, int scale, double *scales, int computeScales, int impute) {
    #pragma omp parallel for schedule(static) default(none) shared(NA_INTEGER, NA_REAL, in, nrows, ncols, out, center, centers, computeCenters, scale, scales, computeScales, impute)
    for (ptrdiff_t j = 0; j < ncols; j++) {
        double mean;
        if (computeCenters || computeScales || impute) {
            double sum = 0;
            double sumsq = 0;
            ptrdiff_t n = 0;
            for (ptrdiff_t i = 0; i < nrows; i++) {
                int *cin = in + j * nrows + i;
                if (*cin != NA_INTEGER) {
                    sum += *cin;
                    sumsq += *cin * *cin;
                    n++;
                }
            }
            mean = sum / n;
            if (computeCenters) {
                centers[j] = mean;
            }
            if (computeScales) {
                scales[j] = sqrt((sumsq - (sum * sum) / n) / (n - 1));
            }
        }
        for (ptrdiff_t i = 0; i < nrows; i++) {
            int *cin = in + j * nrows + i;
            double *cout = out + j * nrows + i;
            if (*cin == NA_INTEGER) {
                if (impute) {
                    if (center) {
                        *cout = 0;
                    } else {
                        *cout = mean;
                    }
                } else {
                    *cout = NA_REAL;
                }
            } else {
                *cout = *cin;
                if (center) {
                    *cout -= centers[j];
                }
                if (scale) {
                    *cout /= scales[j];
                }
            }
        }
    }
}

void preprocess_real(double *in, int nrows, int ncols, double *out, int center, double *centers, int computeCenters, int scale, double *scales, int computeScales, int impute) {
    #pragma omp parallel for schedule(static) default(none) shared(NA_REAL, in, nrows, ncols, out, center, centers, computeCenters, scale, scales, computeScales, impute)
    for (ptrdiff_t j = 0; j < ncols; j++) {
        double mean;
        if (computeCenters || computeScales || impute) {
            double sum = 0;
            double sumsq = 0;
            ptrdiff_t n = 0;
            for (ptrdiff_t i = 0; i < nrows; i++) {
                double *cin = in + j * nrows + i;
                if (!ISNAN(*cin)) {
                    sum += *cin;
                    sumsq += *cin * *cin;
                    n++;
                }
            }
            mean = sum / n;
            if (computeCenters) {
                centers[j] = mean;
            }
            if (computeScales) {
                scales[j] = sqrt((sumsq - (sum * sum) / n) / (n - 1));
            }
        }
        for (ptrdiff_t i = 0; i < nrows; i++) {
            double *cin = in + j * nrows + i;
            double *cout = out + j * nrows + i;
            *cout = *cin;
            if (ISNA(*cin)) {
                if (impute) {
                    if (center) {
                        *cout = 0;
                    } else {
                        *cout = mean;
                    }
                }
            } else {
                if (center) {
                    *cout -= centers[j];
                }
                if (scale) {
                    *cout /= scales[j];
                }
            }
        }
    }
}

SEXP lazy_coerce_to_real(SEXP sIn, SEXP sInPlace) {
    R_xlen_t length = Rf_xlength(sIn);
    int inplace = Rf_asLogical(sInPlace);
    SEXP sOut;
    switch(TYPEOF(sIn)) {
    case REALSXP:
        if (inplace) {
            sOut = PROTECT(sIn); // protection only needed for UNPROTECT count
        } else {
            sOut = PROTECT(Rf_allocVector(REALSXP, length));
            DUPLICATE_ATTRIB(sOut, sIn);
        }
        break;
    case INTSXP:
        if (inplace) {
            Rf_warning("inplace is not supported for integer vectors");
        }
        sOut = PROTECT(Rf_allocVector(REALSXP, length));
        DUPLICATE_ATTRIB(sOut, sIn);
        break;
    }
    UNPROTECT(1); // sOut
    return sOut;
}

SEXP preprocess(SEXP sIn, SEXP sCenter, SEXP sScale, SEXP sImpute, SEXP sInPlace) {
    int nprotect = 0;
    int nrows = Rf_nrows(sIn);
    int ncols = Rf_ncols(sIn);
    int center;
    SEXP sCenters;
    double *centers;
    int computeCenters;
    switch(TYPEOF(sCenter)) {
    case LGLSXP:
        center = Rf_asLogical(sCenter);
        if (center) {
            sCenters = PROTECT(Rf_allocVector(REALSXP, ncols));
            nprotect++;
            centers = REAL(sCenters);
            computeCenters = 1;
        } else {
            centers = NULL;
            computeCenters = 0;
        }
        break;
    case REALSXP:
        center = 1;
        sCenters = PROTECT(Rf_duplicate(sCenter));
        nprotect++;
        centers = REAL(sCenters);
        computeCenters = 0;
        break;
    }
    int scale;
    SEXP sScales;
    double *scales;
    int computeScales;
    switch(TYPEOF(sScale)) {
    case LGLSXP:
        scale = Rf_asLogical(sScale);
        if (scale) {
            sScales = PROTECT(Rf_allocVector(REALSXP, ncols));
            nprotect++;
            scales = REAL(sScales);
            computeScales = 1;
        } else {
            scales = NULL;
            computeScales = 0;
        }
        break;
    case REALSXP:
        scale = 1;
        sScales = PROTECT(Rf_duplicate(sScale));
        nprotect++;
        scales = REAL(sScales);
        computeScales = 0;
        break;
    }
    int impute = Rf_asLogical(sImpute);
    SEXP sOut = PROTECT(lazy_coerce_to_real(sIn, sInPlace));
    nprotect++;
    switch(TYPEOF(sIn)) {
    case REALSXP:
        preprocess_real(
            REAL(sIn),
            nrows,
            ncols,
            REAL(sOut),
            center,
            centers,
            computeCenters,
            scale,
            scales,
            computeScales,
            impute
        );
        break;
    case INTSXP:
        preprocess_int(
            INTEGER(sIn),
            nrows,
            ncols,
            REAL(sOut),
            center,
            centers,
            computeCenters,
            scale,
            scales,
            computeScales,
            impute
        );
        break;
    }
    if (center) {
        Rf_setAttrib(sOut, Rf_install("scaled:center"), sCenters);
    }
    if (scale) {
        Rf_setAttrib(sOut, Rf_install("scaled:scale"), sScales);
    }
    UNPROTECT(nprotect);
    return sOut;
}
