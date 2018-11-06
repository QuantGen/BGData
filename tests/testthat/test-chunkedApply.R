context("chunkedApply")

for (nCores in seq_len(2)) {

    test_that(paste("chunkedMap", "on", nCores, "cores"), {

        hasCores(nCores)

        X <- matrix(data = rnorm(50), nrow = 5, ncol = 10)

        for (chunkSize in c(5, 10)) {
            expect_equal(unlist(chunkedMap(X = X, FUN = rowSums, chunkBy = 1, chunkSize = chunkSize, nCores = nCores)), rowSums(X))
            expect_equal(unlist(chunkedMap(X = X, FUN = colSums, chunkSize = chunkSize, nCores = nCores)), colSums(X))
        }

    })

    test_that(paste("chunkedApply", "on", nCores, "cores"), {

        hasCores(nCores)

        X <- matrix(data = rnorm(50), nrow = 5, ncol = 10)

        for (chunkSize in c(5, 10)) {
            expect_equal(chunkedApply(X = X, MARGIN = 1, FUN = sum, chunkSize = chunkSize, nCores = nCores), rowSums(X))
            expect_equal(chunkedApply(X = X, MARGIN = 2, FUN = sum, chunkSize = chunkSize, nCores = nCores), colSums(X))
            expect_equal(chunkedApply(X = X, MARGIN = 1, FUN = sum, chunkSize = chunkSize, nCores = nCores), apply(X, 1, sum))
            expect_equal(chunkedApply(X = X, MARGIN = 2, FUN = sum, chunkSize = chunkSize, nCores = nCores), apply(X, 2, sum))
        }

    })

}
