#define R_NO_REMAP

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

SEXP rayOLS_real(SEXP X, SEXP y) {
    // Get dimensions of X
    int X_nrow = Rf_nrows(X);
    int X_ncol = Rf_ncols(X);
    // Check if dimensions match
    R_xlen_t y_length = Rf_xlength(y);
    if (X_nrow != y_length) {
        Rf_error("The number of rows in X and the length of y need to match\n");
    }
    // Allocate output matrix
    SEXP out = PROTECT(Rf_allocMatrix(REALSXP, X_ncol, 4));
    // Get data pointers
    double *X_data = REAL(X);
    double *y_data = REAL(y);
    // Iterate over columns of X
    for (int col_idx = 0; col_idx < X_ncol; col_idx++) {
        // Compute number of non-missing values in both x and y (n), and
        // Compute sum of x (xt1) for centering x, and
        // Compute sum of y (yt1) for centering y, and
        // Compute sum of products of x and y (xty) for Cov(x, y), and
        // Compute sum of squares of x (xtx) for Var(x), and
        // Compute sum of squares of y (yty) for RSS
        int n = 0;
        double xt1 = 0;
        double yt1 = 0;
        double xty = 0;
        double xtx = 0;
        double yty = 0;
        for (int row_idx = 0; row_idx < X_nrow; row_idx++) {
            double x_val = X_data[row_idx + (col_idx * X_nrow)];
            if (!(ISNA(x_val) || ISNA(y_data[row_idx]))) {
                n++;
                xt1 += x_val;
                yt1 += y_data[row_idx];
                xty += x_val * y_data[row_idx];
                xtx += x_val * x_val;
                yty += y_data[row_idx] * y_data[row_idx];
            }
        }
        // Center xty, xtx, and yty
        xty -= (xt1 * yt1) / n;
        xtx -= (xt1 * xt1) / n;
        yty -= (yt1 * yt1) / n;
        // Compute beta_1 as Cov(x, y) / Var(x)
        // For centered data, beta_0 will be 0: mean(y) - beta_1 * mean(x)
        double beta_1 = xty / xtx;
        // Compute remaining statistics
        double rss = yty - (xtx * pow(beta_1, 2));
        double se = sqrt((rss / (n - 2)) / xtx);
        double z_stat = beta_1 / se;
        double p_value = Rf_pt(fabs(z_stat), n - 2, 0, 0) * 2;
        // Write results
        REAL(out)[col_idx] = beta_1;
        REAL(out)[col_idx + X_ncol] = se;
        REAL(out)[col_idx + (2 * X_ncol)] = z_stat;
        REAL(out)[col_idx + (3 * X_ncol)] = p_value;
    }
    UNPROTECT(1);
    return out;
}

SEXP rayOLS_integer(SEXP X, SEXP y) {
    // Get dimensions of X
    int X_nrow = Rf_nrows(X);
    int X_ncol = Rf_ncols(X);
    // Check if dimensions match
    R_xlen_t y_length = Rf_xlength(y);
    if (X_nrow != y_length) {
        Rf_error("The number of rows in X and the length of y need to match\n");
    }
    // Allocate output matrix
    SEXP out = PROTECT(Rf_allocMatrix(REALSXP, X_ncol, 4));
    // Get data pointers
    int *X_data = INTEGER(X);
    double *y_data = REAL(y);
    // Iterate over columns of X
    for (int col_idx = 0; col_idx < X_ncol; col_idx++) {
        // Compute number of non-missing values in both x and y (n), and
        // Compute sum of x (xt1) for centering x, and
        // Compute sum of y (yt1) for centering y, and
        // Compute sum of products of x and y (xty) for Cov(x, y), and
        // Compute sum of squares of x (xtx) for Var(x), and
        // Compute sum of squares of y (yty) for RSS
        int n = 0;
        double xt1 = 0;
        double yt1 = 0;
        double xty = 0;
        double xtx = 0;
        double yty = 0;
        for (R_xlen_t row_idx = 0; row_idx < X_nrow; row_idx++) {
            int x_val = X_data[row_idx + (col_idx * X_nrow)];
            if (!(x_val == NA_INTEGER || ISNA(y_data[row_idx]))) {
                n++;
                xt1 += x_val;
                yt1 += y_data[row_idx];
                xty += x_val * y_data[row_idx];
                xtx += x_val * x_val;
                yty += y_data[row_idx] * y_data[row_idx];
            }
        }
        // Center xty, xtx, and yty
        xty -= (xt1 * yt1) / n;
        xtx -= (xt1 * xt1) / n;
        yty -= (yt1 * yt1) / n;
        // Compute beta_1 as Cov(x, y) / Var(x)
        // For centered data, beta_0 will be 0: mean(y) - beta_1 * mean(x)
        double beta_1 = xty / xtx;
        // Compute remaining statistics
        double rss = yty - (xtx * pow(beta_1, 2));
        double se = sqrt((rss / (n - 2)) / xtx);
        double z_stat = beta_1 / se;
        double p_value = Rf_pt(fabs(z_stat), n - 2, 0, 0) * 2;
        // Write results
        REAL(out)[col_idx] = beta_1;
        REAL(out)[col_idx + X_ncol] = se;
        REAL(out)[col_idx + (2 * X_ncol)] = z_stat;
        REAL(out)[col_idx + (3 * X_ncol)] = p_value;
    }
    UNPROTECT(1);
    return out;
}

SEXP rayOLS(SEXP X, SEXP y) {
    // Dispatch to real or integer function
    // TODO: Macro-based generics
    switch (TYPEOF(X)) {
        case REALSXP:
            return rayOLS_real(X, y);
            break;
        case INTSXP:
            return rayOLS_integer(X, y);
            break;
        default:
            Rf_error("x needs to be a numeric vector");
            break;
    }
}
