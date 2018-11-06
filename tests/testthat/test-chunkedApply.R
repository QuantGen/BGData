context("chunkedApply")

set.seed(1)

nRows <- 5
nCols <- 10
nNAs <- 5

X <- matrix(data = rnorm(nRows * nCols, sd = 100), nrow = nRows, ncol = nCols)
X[sample(1:length(X), size = nNAs)] <- NA

test_that("chunkedMap", {

    for (nCores in seq_len(2)) {

        hasCores(nCores)

        for (chunkSize in c(5, 10)) {

            expect_equal(unlist(chunkedMap(X = X, FUN = rowSums, chunkBy = 1, chunkSize = chunkSize, nCores = nCores)), rowSums(X))
            expect_equal(unlist(chunkedMap(X = X, FUN = colSums, chunkSize = chunkSize, nCores = nCores)), colSums(X))

            expect_equal(unlist(chunkedMap(X = X, FUN = rowSums, chunkBy = 1, i = c(1, 3), chunkSize = chunkSize, nCores = nCores)), rowSums(X[c(1, 3), ]))
            expect_equal(unlist(chunkedMap(X = X, FUN = colSums, i = c(1, 3), chunkSize = chunkSize, nCores = nCores)), colSums(X[c(1, 3), ]))

            expect_equal(unlist(chunkedMap(X = X, FUN = rowSums, chunkBy = 1, j = c(1, 3, 5), chunkSize = chunkSize, nCores = nCores)), rowSums(X[, c(1, 3, 5)]))
            expect_equal(unlist(chunkedMap(X = X, FUN = colSums, j = c(1, 3, 5), chunkSize = chunkSize, nCores = nCores)), colSums(X[, c(1, 3, 5)]))

            expect_equal(unlist(chunkedMap(X = X, FUN = rowSums, chunkBy = 1, i = c(1, 3), j = c(1, 3, 5), chunkSize = chunkSize, nCores = nCores)), rowSums(X[c(1, 3), c(1, 3, 5)]))
            expect_equal(unlist(chunkedMap(X = X, FUN = colSums, i = c(1, 3), j = c(1, 3, 5), chunkSize = chunkSize, nCores = nCores)), colSums(X[c(1, 3), c(1, 3, 5)]))

        }
    }

})

test_that("chunkedApply", {

    for (nCores in seq_len(2)) {

        hasCores(nCores)

        for (chunkSize in c(5, 10)) {

            expect_equal(chunkedApply(X = X, MARGIN = 1, FUN = sum, chunkSize = chunkSize, nCores = nCores), apply(X, 1, sum))
            expect_equal(chunkedApply(X = X, MARGIN = 2, FUN = sum, chunkSize = chunkSize, nCores = nCores), apply(X, 2, sum))

            expect_equal(chunkedApply(X = X, MARGIN = 1, FUN = sum, i = c(1, 3), chunkSize = chunkSize, nCores = nCores), apply(X[c(1, 3), ], 1, sum))
            expect_equal(chunkedApply(X = X, MARGIN = 2, FUN = sum, i = c(1, 3), chunkSize = chunkSize, nCores = nCores), apply(X[c(1, 3), ], 2, sum))

            expect_equal(chunkedApply(X = X, MARGIN = 1, FUN = sum, j = c(1, 3, 5), chunkSize = chunkSize, nCores = nCores), apply(X[, c(1, 3, 5)], 1, sum))
            expect_equal(chunkedApply(X = X, MARGIN = 2, FUN = sum, j = c(1, 3, 5), chunkSize = chunkSize, nCores = nCores), apply(X[, c(1, 3, 5)], 2, sum))

            expect_equal(chunkedApply(X = X, MARGIN = 1, FUN = sum, i = c(1, 3), j = c(1, 3, 5), chunkSize = chunkSize, nCores = nCores), apply(X[c(1, 3), c(1, 3, 5)], 1, sum))
            expect_equal(chunkedApply(X = X, MARGIN = 2, FUN = sum, i = c(1, 3), j = c(1, 3, 5), chunkSize = chunkSize, nCores = nCores), apply(X[c(1, 3), c(1, 3, 5)], 2, sum))

        }

    }

})
