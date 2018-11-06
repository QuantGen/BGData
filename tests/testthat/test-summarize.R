context("summarize")

set.seed(1)

nRows <- 5
nCols <- 10
percentNA <- 0.15

summarize_R <- function(X) {
    res <- data.frame(
        freq_na = vector(mode = "double", length = ncol(X)),
        allele_freq = vector(mode = "double", length = ncol(X)),
        sd = vector(mode = "double", length = ncol(X))
    )
    for (col in seq_len(ncol(X))) {
        x <- X[, col]
        nMissing <- sum(is.na(x))
        res$freq_na[col] <- nMissing / length(x)
        res$allele_freq[col] <- sum(x, na.rm = TRUE) / ((length(x) - nMissing) * 2)
        res$sd[col] <- sd(x, na.rm = TRUE)
    }
    return(res)
}

test_that("summarize", {

    for (mode in c("integer", "double")) {

        X <- matrix(data = rnorm(nRows * nCols, sd = 100), nrow = nRows, ncol = nCols)
        X[sample(seq_along(X), size = as.integer(length(X) * percentNA))] <- NA
        storage.mode(X) <- mode

        expect_equal(
            summarize(X),
            summarize_R(X)
        )

    }

})
