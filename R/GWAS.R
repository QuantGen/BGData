#' Performs Single Marker Regressions Using BGData Objects.
#'
#' Implements single marker regressions. The regression model includes all the
#' covariates specified in the right-hand-side of the `formula` plus one column
#' of `@@geno` at a time. The data from the association tests is obtained from
#' a [BGData-class] object.
#'
#' @inheritSection BGData-package File-backed matrices
#' @inheritSection BGData-package Multi-level parallelism
#' @param formula The formula for the GWAS model without including the marker,
#' e.g. `y ~ 1` or `y ~ factor(sex) + age`. The variables included in the
#' formula must be in the `@@pheno` object of the [BGData-class].
#' @param data A [BGData-class] object.
#' @param method The regression method to be used. Currently, the following
#' methods are implemented: [stats::lm()], [stats::lm.fit()], [stats::lsfit()],
#' [stats::glm()], [lme4::lmer()], and [SKAT::SKAT()]. Defaults to `lsfit`.
#' @param i Indicates which rows of `@@geno` should be used. Can be integer,
#' boolean, or character. By default, all rows are used.
#' @param j Indicates which columns of `@@geno` should be used. Can be integer,
#' boolean, or character. By default, all columns are used.
#' @param chunkSize The number of columns of `@@geno` that are brought into
#' physical memory for processing per core. If `NULL`, all elements in `j` are
#' used. Defaults to 5000.
#' @param nCores The number of cores (passed to [parallel::mclapply()]).
#' Defaults to the number of cores as detected by [parallel::detectCores()].
#' @param verbose Whether progress updates will be posted. Defaults to `FALSE`.
#' @param ... Additional arguments for chunkedApply and regression method.
#' @return The same matrix that would be returned by `coef(summary(model))`.
#' @example man/examples/GWAS.R
#' @export
GWAS <- function(formula, data, method = "lsfit", i = seq_len(nrow(data@geno)), j = seq_len(ncol(data@geno)), chunkSize = 5000L, nCores = getOption("mc.cores", 2L), verbose = FALSE, ...) {

    if (class(data) != "BGData") {
        stop("data must BGData")
    }

    if (!method %in% c("lm", "lm.fit", "lsfit", "glm", "lmer", "SKAT", "rayOLS")) {
        stop("Only lm, lm.fit, lsfit, glm, lmer, SKAT, and rayOLS have been implemented so far.")
    }

    i <- crochet::convertIndex(data@geno, i, "i")
    j <- crochet::convertIndex(data@geno, j, "j")

    if (method == "lsfit") {
        OUT <- GWAS.lsfit(formula = formula, data = data, i = i, j = j, chunkSize = chunkSize, nCores = nCores, verbose = verbose, ...)
    } else if (method == "rayOLS") {
        if (length(attr(stats::terms(formula), "term.labels")) > 0L) {
            stop("method rayOLS can only be used with y~1 formula, if you want to add covariates pre-adjust your phenotype.")
        }
        OUT <- GWAS.rayOLS(formula = formula, data = data, i = i, j = j, chunkSize = chunkSize, nCores = nCores, verbose = verbose, ...)
    } else if (method == "SKAT") {
        if (!requireNamespace("SKAT", quietly = TRUE)) {
            stop("SKAT needed for this function to work. Please install it.", call. = FALSE)
        }
        OUT <- GWAS.SKAT(formula = formula, data = data, i = i, j = j, verbose = verbose, ...)
    } else {
        if (method == "lmer") {
            if (!requireNamespace("lme4", quietly = TRUE)) {
                stop("lme4 needed for this function to work. Please install it.", call. = FALSE)
            }
            FUN <- lme4::lmer
        } else {
            FUN <- match.fun(method)
        }
        pheno <- data@pheno
        GWAS.model <- stats::update(stats::as.formula(formula), ".~z+.")
        OUT <- chunkedApply(X = data@geno, MARGIN = 2L, FUN = function(col, ...) {
            pheno$z <- col
            fm <- FUN(GWAS.model, data = pheno, ...)
            getCoefficients(fm)
        }, i = i, j = j, chunkSize = chunkSize, nCores = nCores, verbose = verbose, ...)
        colnames(OUT) <- colnames(data@geno)[j]
        OUT <- t(OUT)
    }

    return(OUT)
}


# the GWAS method for rayOLS
GWAS.rayOLS <- function(formula, data, i = seq_len(nrow(data@geno)), j = seq_len(ncol(data@geno)), chunkSize = 5000L, nCores = getOption("mc.cores", 2L), verbose = FALSE, ...) {
    y <- data@pheno[i, as.character(stats::terms(formula)[[2L]]), drop = TRUE]
    y <- y - mean(y, na.rm = TRUE)
    n <- length(y)
    Int <- rep(1, n)
    SSy <- sum(y^2, na.rm = TRUE)
    isNAY <- which(is.na(y))
    res <- chunkedApply(X = data@geno, MARGIN = 2L, FUN = rayOLS, i = i, j = j, chunkSize = chunkSize, nCores = nCores, verbose = verbose, y = y, Int = Int, SSy = SSy, n = n, isNAY = isNAY, ...)
    return(t(res))
}


GWAS.lsfit <- function(formula, data, i = seq_len(nrow(data@geno)), j = seq_len(ncol(data@geno)), chunkSize = 5000L, nCores = getOption("mc.cores", 2L), verbose = FALSE, ...) {

    # subset of model.frame has bizarre scoping issues
    frame <- stats::model.frame(formula = formula, data = data@pheno)[i, , drop = FALSE]
    model <- stats::model.matrix(formula, frame)
    model <- cbind(1L, model) # Reserve space for marker column

    y <- data@pheno[i, as.character(stats::terms(formula)[[2L]]), drop = TRUE]

    res <- chunkedApply(X = data@geno, MARGIN = 2L, FUN = function(col, ...) {
        model[, 1L] <- col
        fm <- stats::lsfit(x = model, y = y, intercept = FALSE)
        stats::ls.print(fm, print.it = FALSE)$coef.table[[1L]][1L, ]
    }, i = i, j = j, chunkSize = chunkSize, nCores = nCores, verbose = verbose, ...)
    colnames(res) <- colnames(data@geno)[j]
    res <- t(res)

    return(res)
}


# formula: the formula for the GWAS model without including the markers, e.g.
# y~1 or y~factor(sex)+age
# all the variables in the formula must be in data@pheno data (BGData)
# containing slots @pheno and @geno
# groups: a vector mapping markers into groups (can be integer, character or
# factor)
GWAS.SKAT <- function(formula, data, groups, i = seq_len(nrow(data@geno)), j = seq_len(ncol(data@geno)), verbose = FALSE, ...) {

    uniqueGroups <- unique(groups)

    OUT <- matrix(data = double(), nrow = length(uniqueGroups), ncol = 2L)
    colnames(OUT) <- c("nMrk", "p-value")
    rownames(OUT) <- uniqueGroups

    H0 <- SKAT::SKAT_Null_Model(formula, data = data@pheno[i, , drop = FALSE], ...)

    for (group in seq_along(uniqueGroups)) {
        Z <- data@geno[i, groups == uniqueGroups[group], drop = FALSE]
        fm <- SKAT::SKAT(Z = Z, obj = H0, ...)
        OUT[group, ] <- c(ncol(Z), fm$p.value)
        if (verbose) {
            message("Group ", group, " of ", length(uniqueGroups), " ...")
        }
    }

    return(OUT)
}


# OLS for the regression y=xb+e (data is assumed to be pre-adjusted by non-genetic effects
rayOLS <- function(y, x, SSy, Int, n, isNAY) {
    isNAX <- which(is.na(x))
    isNAXY <- unique(c(isNAX, isNAY))
    SSy <- SSy - sum(y[isNAX]^2, na.rm = TRUE)
    if (length(isNAXY) > 0) {
        y[isNAXY] <- 0
        x[isNAXY] <- 0
    }
    n <- n - length(isNAXY)
    # crossproducts
    sX <- crossprod(Int, x)
    sY <- crossprod(Int, y)
    XX <- crossprod(x) - sX * sX / n
    Xy <- crossprod(x, y) - sX * sY / n
    # solution and SE
    sol <- Xy / XX
    RSS <- SSy - XX * sol^2
    SE <- sqrt(RSS / (n - 2L) / XX)
    z_stat <- sol / SE
    p_val <- stats::pt(q = abs(z_stat), df = n - 2L, lower.tail = FALSE) * 2L
    return(c(sol, SE, z_stat, p_val))
}


getCoefficients <- function(x) {
    UseMethod("getCoefficients")
}


getCoefficients.lm <- function(x) {
    summary(x)$coef[2L, ]
}


getCoefficients.glm <- function(x) {
    summary(x)$coef[2L, ]
}


getCoefficients.lmerMod <- function(x) {
    ans <- summary(x)$coef[2L, ]
    ans <- c(ans, c(1L - stats::pnorm(ans[3L])))
    return(ans)
}
