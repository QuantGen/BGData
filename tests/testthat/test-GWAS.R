context("GWAS")

set.seed(1)

nRows <- 5
nCols <- 10
nNAs <- 5

X <- matrix(data = rnorm(nRows * nCols), nrow = nRows, ncol = nCols)
y <- data.frame(y = rnorm(nRows))

DATA <- BGData(geno = X, pheno = y)

lsfit_R <- function(X, y) {
    res <- apply(X, 2, function(x) {
        fm <- lsfit(x = cbind(x, 1), y = y, intercept = FALSE)
        ls.print(fm, print.it = FALSE)$coef.table[[1]][1, ]
    })
    res <- t(res)
    rownames(res) <- colnames(X)
    return(res)
}

test_that("lsfit", {

    expect_equal(GWAS(formula = y ~ 1, data = DATA, method = "lsfit"), lsfit_R(DATA@geno, DATA@pheno))

})
