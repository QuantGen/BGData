context("summarize")

set.seed(1)

n <- 5
p <- 10

summarize_R <- function(X, i = seq_len(nrow(X)), j = seq_len(ncol(X))) {
    res <- data.frame(
        freq_na = vector(mode = "double", length = length(j)),
        allele_freq = vector(mode = "double", length = length(j)),
        sd = vector(mode = "double", length = length(j))
    )
    for (col in seq_along(j)) {
        x <- X[i, j[col]]
        nMissing <- sum(is.na(x))
        res$freq_na[col] <- nMissing / length(x)
        res$allele_freq[col] <- sum(x, na.rm = TRUE) / ((length(x) - nMissing) * 2)
        res$sd[col] <- sd(x, na.rm = TRUE)
    }
    return(res)
}

test_that("summarize", {

    for (mode in c("integer", "double")) {

        X <- matrix(data = rnorm(n * p, sd = 100), nrow = n, ncol = p)
        X[sample(1:length(X), size = 5)] <- NA
        storage.mode(X) <- mode

        for (nCores in seq_len(2)) {

            hasCores(nCores)

            for (chunkSize in c(3, 6)) {

                expect_equal(summarize(X = X, chunkSize = chunkSize, nCores = nCores), summarize_R(X))
                expect_equal(summarize(X = X, i = c(1, 3), chunkSize = chunkSize, nCores = nCores), summarize_R(X, i = c(1, 3)))
                expect_equal(summarize(X = X, j = c(1, 3, 5), chunkSize = chunkSize, nCores = nCores), summarize_R(X, j =  c(1, 3, 5)))
                expect_equal(summarize(X = X, i = c(1, 3), j = c(1, 3, 5), chunkSize = chunkSize, nCores = nCores), summarize_R(X, i = c(1, 3), j = c(1, 3, 5)))
                expect_equal(summarize(X = X, nCores = nCores), summarize_R(X))

            }

        }

    }

})
