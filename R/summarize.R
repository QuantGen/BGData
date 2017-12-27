#' Generates Various Summary Statistics.
#'
#' Computes the frequency of missing values, the (minor) allele frequency, and
#' standard deviation of each column of `X`.
#'
#' @inheritSection BGData-package Memory-mapping
#' @inheritSection BGData-package Multi-level parallelism
#' @param X A matrix-like object, typically `@@geno` of a [BGData-class]
#' object.
#' @param i Indicates which rows of `X` should be used. Can be integer,
#' boolean, or character. By default, all rows are used.
#' @param j Indicates which columns of `X` should be used. Can be integer,
#' boolean, or character. By default, all columns are used.
#' @param bufferSize The number of columns of `X` that are brought into
#' physical memory for processing per core. If `NULL`, all elements in `j` are
#' used. Defaults to 5000.
#' @param nCores The number of cores (passed to [parallel::mclapply()]).
#' Defaults to the number of cores as detected by [parallel::detectCores()].
#' @param verbose Whether progress updates will be posted. Defaults to `FALSE`.
#' @return A `data.frame` with three columns: `freq_na` for frequencies of
#' missing values, `allele_freq` for (minor) allele frequencies, and `sd` for
#' standard deviations.
#' @example man/examples/summarize.R
#' @export
summarize <- function(X, i = seq_len(nrow(X)), j = seq_len(ncol(X)), bufferSize = 5000L, nCores = getOption("mc.cores", 2L), verbose = FALSE) {
    i <- convertIndexTypes(i, rownames(X))
    j <- convertIndexTypes(j, colnames(X))
    m <- chunkedApply(X = X, MARGIN = 2L, FUN = function(col) {
        freqNA <- mean(is.na(col))
        alleleFreq <- mean(col, na.rm = TRUE) / 2L
        sd <- stats::sd(col, na.rm = TRUE)
        cbind(freqNA, alleleFreq, sd)
    }, i = i, j = j, bufferSize = bufferSize, nCores = nCores, verbose = verbose)
    df <- data.frame(
        freq_na = m[1L, ],
        allele_freq = m[2L, ],
        sd = m[3L, ],
        stringsAsFactors = FALSE
    )
    rownames(df) <- colnames(X)[j]
    return(df)
}
