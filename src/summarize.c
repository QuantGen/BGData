#define R_NO_REMAP

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

SEXP summarize_real(SEXP X) {
    // Get dimensions of X
    int nrow = Rf_nrows(X);
    int ncol = Rf_ncols(X);
    // Get data pointer
    double *X_data = REAL(X);
    // Allocate output matrix
    SEXP out = PROTECT(Rf_allocMatrix(REALSXP, ncol, 3));
    // Iterate over columns of X
    int col_idx = 0;
    int row_idx = 0;
    for (col_idx = 0; col_idx < ncol; col_idx++) {
        // Compute number of non-missing values (n), and
        // Compute column sum (xt1), and
        // Compute column sum of squares (xtx)
        R_xlen_t n = 0;
        double xt1 = 0;
        double xtx = 0;
        for (row_idx = 0; row_idx < nrow; row_idx++) {
            double x_val = X_data[row_idx + (col_idx * nrow)];
            if (!ISNA(x_val)) {
                n++;
                xt1 += x_val;
                xtx += x_val * x_val;
            }
        }
        double freq_na;
        double allele_freq;
        double sd;
        if (n) {
            // Center xtx
            xtx -= (xt1 * xt1) / n;
            // Compute summary statistics
            freq_na = (nrow - n) / (double) nrow;
            allele_freq = xt1 / n / 2;
            sd = sqrt(xtx / (n - 1));
        } else {
            freq_na = 1;
            allele_freq = NA_REAL;
            sd = NA_REAL;
        }
        // Write results into output matrix
        REAL(out)[col_idx] = freq_na;
        REAL(out)[col_idx + ncol] = allele_freq;
        REAL(out)[col_idx + (2 * ncol)] = sd;
    }
    UNPROTECT(1);
    return out;
}

SEXP summarize_integer(SEXP X) {
    // Get dimensions of X
    int nrow = Rf_nrows(X);
    int ncol = Rf_ncols(X);
    // Get data pointer
    int *X_data = INTEGER(X);
    // Allocate output matrix
    SEXP out = PROTECT(Rf_allocMatrix(REALSXP, ncol, 3));
    // Iterate over columns of X
    int col_idx = 0;
    int row_idx = 0;
    for (col_idx = 0; col_idx < ncol; col_idx++) {
        // Compute number of non-missing values (n), and
        // Compute column sum (xt1), and
        // Compute column sum of squares (xtx)
        R_xlen_t n = 0;
        double xt1 = 0;
        double xtx = 0;
        for (row_idx = 0; row_idx < nrow; row_idx++) {
            int x_val = X_data[row_idx + (col_idx * nrow)];
            if (x_val != NA_INTEGER) {
                n++;
                xt1 += x_val;
                xtx += x_val * x_val;
            }
        }
        double freq_na;
        double allele_freq;
        double sd;
        if (n) {
            // Center xtx
            xtx -= (xt1 * xt1) / n;
            // Compute summary statistics
            freq_na = (nrow - n) / (double) nrow;
            allele_freq = xt1 / n / 2;
            sd = sqrt(xtx / (n - 1));
        } else {
            freq_na = 1;
            allele_freq = NA_REAL;
            sd = NA_REAL;
        }
        // Write results into output matrix
        REAL(out)[col_idx] = freq_na;
        REAL(out)[col_idx + ncol] = allele_freq;
        REAL(out)[col_idx + (2 * ncol)] = sd;
    }
    UNPROTECT(1);
    return out;
}

SEXP summarize(SEXP X) {
    // Dispatch to real or integer function
    // TODO: Macro-based generics
    switch (TYPEOF(X)) {
        case REALSXP:
            return summarize_real(X);
            break;
        case INTSXP:
            return summarize_integer(X);
            break;
        default:
            Rf_error("X needs to be a numeric matrix");
            break;
    }
}
