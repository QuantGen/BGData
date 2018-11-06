context("GWAS")

for (nCores in seq_len(2)) {

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

        for (chunkSize in c(3, 6)) {
            fm <- GWAS(formula = y ~ 1, data = DATA, method = "lsfit", chunkSize = chunkSize, nCores = nCores)
            expect_equal(comp, fm)
        }

        # With i
        i <- seq_len(nrow(X))[-3]
        comp <- t(apply(X[i, ], 2, function(z) {
            fm <- lsfit(x = cbind(z, 1), y = y[i, ], intercept = FALSE)
            ls.print(fm, print.it = FALSE)$coef.table[[1]][1, ]
        }))
        rownames(comp) <- colnames(DATA@geno)

        for (chunkSize in c(3, 6)) {
            fm <- GWAS(formula = y ~ 1, data = DATA, method = "lsfit", i = i, chunkSize = chunkSize, nCores = nCores)
            expect_equal(comp, fm)
        }

        # With j
        j <- seq_len(ncol(X))[-3]
        comp <- t(apply(X[, j], 2, function(z) {
            fm <- lsfit(x = cbind(z, 1), y = y, intercept = FALSE)
            ls.print(fm, print.it = FALSE)$coef.table[[1]][1, ]
        }))
        rownames(comp) <- colnames(DATA@geno)[j]

        for (chunkSize in c(3, 6)) {
            fm <- GWAS(formula = y ~ 1, data = DATA, method = "lsfit", j = j, chunkSize = chunkSize, nCores = nCores)
            expect_equal(comp, fm)
        }

    })

}
