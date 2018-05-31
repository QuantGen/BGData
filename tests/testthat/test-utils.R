library(parallel)

context("utils")

testDir <- function() {
    paste0(tempdir(), "/BGData-", BGData:::randomString(), "/")
}

hasCores <- function(numCores) {
    # For CRAN
    if (Sys.getenv("_R_CHECK_LIMIT_CORES_") == TRUE || numCores > parallel::detectCores()) {
        skip("Not enough cores or number of cores capped for CRAN submission checks.")
    }
    # For WinBuilder
    if (.Platform$OS.type == "windows" && numCores > 1) {
        skip("mc.cores > 1 is not supported on Windows.")
    }
}

for (nCores in seq_len(2)) {

    test_that(paste("chunkedApply", "on", nCores, "cores"), {

        hasCores(nCores)

        X <- matrix(data = rnorm(50), nrow = 5, ncol = 10)

        for (bufferSize in c(5, 10)) {
            expect_equal(chunkedApply(X = X, MARGIN = 1, FUN = sum, bufferSize = bufferSize, nCores = nCores), rowSums(X))
            expect_equal(chunkedApply(X = X, MARGIN = 2, FUN = sum, bufferSize = bufferSize, nCores = nCores), colSums(X))
            expect_equal(chunkedApply(X = X, MARGIN = 1, FUN = sum, bufferSize = bufferSize, nCores = nCores), apply(X, 1, sum))
            expect_equal(chunkedApply(X = X, MARGIN = 2, FUN = sum, bufferSize = bufferSize, nCores = nCores), apply(X, 2, sum))
        }

    })

    test_that(paste("getG_symDMatrix", "on", nCores, "cores"), {

        hasCores(nCores)

        W <- matrix(data = rnorm(200), nrow = 10, ncol = 20)
        G1 <- tcrossprod(scale(W))
        G1 <- G1 / mean(diag(G1))

        G2 <- getG_symDMatrix(X = W, blockSize = ceiling(nrow(W) / 3), folderOut = testDir(), nCores = nCores)
        expect_equivalent(G2[], G1) # use equivalent to correct slight difference in NULL dimnames handling

    })

    test_that(paste("getGi", "on", nCores, "cores"), {

        hasCores(nCores)

        n <- 10
        p <- 100
        X <- matrix(data = rnorm(n * p), nrow = n, ncol = p)

        for (bufferSize in c(NULL, p, ceiling(p / 3))) {

            # both scalings
            G <- tcrossprod(scale(X))
            G <- G / mean(diag(G))
            G2 <- getG(X = X, scale = TRUE, scaleG = TRUE, bufferSize = bufferSize, nCores = nCores)
            expect_equal(G, G2, check.attributes = FALSE)

            # without scaling to average diagonal = 1 (scaleG)
            G <- tcrossprod(scale(X))
            G2 <- getG(X = X, scale = TRUE, scaleG = FALSE, bufferSize = bufferSize, nCores = nCores)
            expect_equal(G, G2, check.attributes = FALSE)

            # without scaling columns, but scaling average diagonal = 1 (scaleG)
            G <- tcrossprod(scale(X, center = TRUE, scale = FALSE))
            G <- G / mean(diag(G))
            G2 <- getG(X = X, scale = FALSE, scaleG = TRUE, bufferSize = bufferSize, nCores = nCores)

            expect_equal(G, G2, check.attributes = FALSE)

            # no scaling at all
            G <- tcrossprod(scale(X, center = TRUE, scale = FALSE))
            G2 <- getG(X = X, scale = FALSE, scaleG = FALSE, bufferSize = bufferSize, nCores = nCores)
            expect_equal(G, G2, check.attributes = FALSE)

            # neither scaling nor centering
            G <- tcrossprod(X)
            G2 <- getG(X = X, center = FALSE, scale = FALSE, scaleG = FALSE, bufferSize = bufferSize, nCores = nCores)
            expect_equal(G, G2, check.attributes = FALSE)

        }

        X[sample(1:length(X), size = 20)] <- NA
        G <- getG(X, nCores = nCores)
        expect_true(!any(is.na(G)))

    })


    test_that(paste("getGij", "on", nCores, "cores"), {

        hasCores(nCores)

        n <- 10
        p <- 100
        X <- matrix(data = rnorm(n * p), nrow = n, ncol = p)

        for (bufferSize in c(NULL, p, ceiling(p / 3))) {

            i <- sample(1:nrow(X), size = 3)
            i2 <- sample(1:nrow(X), size = 4)

            centers <- colMeans(X)
            scales <- apply(X, 2, sd) * sqrt((n - 1)/n)

            # all scalings
            G <- tcrossprod(scale(X))
            G <- G / mean(diag(G))
            G_12 <- getG(X = X, center = centers, scale = scales, scaleG = TRUE, i = i, i2 = i2, bufferSize = bufferSize, nCores = nCores)
            expect_equal(G[i, i2], G_12, check.attributes = FALSE)

            G_12 <- getG(X = X, center = centers, scale = scales, scaleG = TRUE, i = i, i2 = i, bufferSize = bufferSize, nCores = nCores)
            expect_equal(G[i, i], G_12, check.attributes = FALSE)

            # without scaling to average diagonal = 1
            G <- tcrossprod(scale(X) * sqrt(n/(n - 1)))
            G_12 <- getG(X = X, center = centers, scale = scales, scaleG = FALSE, i = i, i2 = i2, bufferSize = bufferSize, nCores = nCores)
            expect_equal(G[i, i2], G_12, check.attributes = FALSE)

            G_12 <- getG(X = X, center = centers, scale = scales, scaleG = FALSE, i = i, i2 = i, bufferSize = bufferSize, nCores = nCores)
            expect_equal(G[i, i], G_12, check.attributes = FALSE)

            # without scaling columns, but scaling average diagonal = 1
            scales <- rep(1, ncol(X))

            G <- tcrossprod(scale(X, center = TRUE, scale = FALSE))
            G <- G / ncol(X)
            G_12 <- getG(X = X, center = centers, scale = scales, scaleG = TRUE, i = i, i2 = i2, bufferSize = bufferSize, nCores = nCores)
            expect_equal(G[i, i2], G_12, check.attributes = FALSE)

            G_12 <- getG(X = X, center = centers, scale = scales, scaleG = TRUE, i = i, i2 = i, bufferSize = bufferSize, nCores = nCores)
            expect_equal(G[i, i], G_12, check.attributes = FALSE)

            # no scaling at all
            G <- tcrossprod(scale(X, center = TRUE, scale = FALSE))
            G_12 <- getG(X = X, center = centers, scale = scales, scaleG = FALSE, i = i, i2 = i2, bufferSize = bufferSize, nCores = nCores)
            expect_equal(G[i, i2], G_12, check.attributes = FALSE)

            G_12 <- getG(X = X, center = centers, scale = scales, scaleG = FALSE, i = i, i2 = i, bufferSize = bufferSize, nCores = nCores)
            expect_equal(G[i, i], G_12, check.attributes = FALSE)

        }
    })


    test_that(paste("summarize", "on", nCores, "cores"), {

        hasCores(nCores)

        genotypes <- matrix(data = c(0, 0, 1, 0, 2, 2, 1, 2, 0, 1, 2, 0, 0, 1, 2, 0, NA, 0), nrow = 3, ncol = 6)

        computeDummy <- function(i = seq_len(nrow(genotypes)), j = seq_len(ncol(genotypes))) {
            dummy <- data.frame(
                freq_na = vector(mode = "double", length = length(j)),
                allele_freq = vector(mode = "double", length = length(j)),
                sd = vector(mode = "double", length = length(j))
            )
            for (col in seq_len(length(j))) {
                Z <- genotypes[i, j[col]]
                NAs <- sum(is.na(Z))
                dummy$freq_na[col] <- NAs / length(Z)
                dummy$allele_freq[col] <- sum(Z, na.rm = TRUE) / ((length(Z) - NAs) * 2)
                dummy$sd[col] <- sd(Z, na.rm = TRUE)
            }
            return(dummy)
        }

        for (bufferSize in c(3, 6)) {
            expect_equal(summarize(X = genotypes, bufferSize = bufferSize, nCores = nCores), computeDummy())
            expect_equal(summarize(X = genotypes, i = c(1, 3), bufferSize = bufferSize, nCores = nCores), computeDummy(i = c(1, 3)))
            expect_equal(summarize(X = genotypes, j = c(1, 3, 5), bufferSize = bufferSize, nCores = nCores), computeDummy(j =  c(1, 3, 5)))
            expect_equal(summarize(X = genotypes, i = c(1, 3), j = c(1, 3, 5), bufferSize = bufferSize, nCores = nCores), computeDummy(i = c(1, 3), j = c(1, 3, 5)))
            expect_equal(summarize(X = genotypes, nCores = nCores), computeDummy())
        }

    })


    test_that(paste("GWAS", "on", nCores, "cores"), {

        hasCores(nCores)

        X <- matrix(data = rnorm(50), nrow = 5, ncol = 10)
        y <- data.frame(y = rnorm(5))

        DATA <- BGData(geno = X, pheno = y)

        # Without i
        comp <- t(apply(X, 2, function(z) {
            fm <- lsfit(x = cbind(z, 1), y = y, intercept = FALSE)
            ls.print(fm, print.it = FALSE)$coef.table[[1]][1, ]
        }))
        rownames(comp) <- colnames(DATA@geno)

        for (bufferSize in c(3, 6)) {
            fm <- GWAS(formula = y ~ 1, data = DATA, method = "lsfit", bufferSize = bufferSize, nCores = nCores)
            expect_equal(comp, fm)
        }

        # With i
        i <- seq_len(nrow(X))[-3]
        comp <- t(apply(X[i, ], 2, function(z) {
            fm <- lsfit(x = cbind(z, 1), y = y[i, ], intercept = FALSE)
            ls.print(fm, print.it = FALSE)$coef.table[[1]][1, ]
        }))
        rownames(comp) <- colnames(DATA@geno)

        for (bufferSize in c(3, 6)) {
            fm <- GWAS(formula = y ~ 1, data = DATA, method = "lsfit", i = i, bufferSize = bufferSize, nCores = nCores)
            expect_equal(comp, fm)
        }

        # With j
        j <- seq_len(ncol(X))[-3]
        comp <- t(apply(X[, j], 2, function(z) {
            fm <- lsfit(x = cbind(z, 1), y = y, intercept = FALSE)
            ls.print(fm, print.it = FALSE)$coef.table[[1]][1, ]
        }))
        rownames(comp) <- colnames(DATA@geno)[j]

        for (bufferSize in c(3, 6)) {
            fm <- GWAS(formula = y ~ 1, data = DATA, method = "lsfit", j = j, bufferSize = bufferSize, nCores = nCores)
            expect_equal(comp, fm)
        }

    })

}


test_that("normalizeType", {

    expect_equal(typeof(normalizeType("double")), "double")
    expect_equal(typeof(normalizeType(double())), "double")
    expect_equal(typeof(normalizeType("integer")), "integer")
    expect_equal(typeof(normalizeType(integer())), "integer")
    expect_equal(typeof(normalizeType("character")), "character")
    expect_equal(typeof(normalizeType(character())), "character")
    expect_equal(typeof(normalizeType("complex")), "complex")
    expect_equal(typeof(normalizeType(complex())), "complex")
    expect_warning(normalizeType("test"))
    expect_equal(suppressWarnings(typeof(normalizeType("test"))), "character")
    expect_equal(typeof(normalizeType(1)), "double")
    expect_equal(typeof(normalizeType(1L)), "integer")

})
