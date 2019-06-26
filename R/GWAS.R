GWAS <- function(formula, data, method = "lsfit", i = seq_len(nrow(geno(data))), j = seq_len(ncol(geno(data))), chunkSize = 5000L, nCores = getOption("mc.cores", 2L), verbose = FALSE, ...) {

    if (class(data) != "BGData") {
        stop("data must BGData")
    }

    if (!method %in% c("rayOLS", "lsfit", "lm", "lm.fit", "glm", "lmer", "SKAT")) {
        stop("Only rayOLS, lsfit, lm, lm.fit, glm, lmer, and SKAT have been implemented so far.")
    }

    i <- crochet::convertIndex(geno(data), i, "i")
    j <- crochet::convertIndex(geno(data), j, "j")

    if (method == "rayOLS") {
        if (length(labels(stats::terms(formula))) > 0L) {
            stop("method rayOLS can only be used with y~1 formula, if you want to add covariates pre-adjust your phenotype.")
        }
        OUT <- GWAS.rayOLS(formula = formula, data = data, i = i, j = j, chunkSize = chunkSize, nCores = nCores, verbose = verbose, ...)
    } else if (method == "lsfit") {
        OUT <- GWAS.lsfit(formula = formula, data = data, i = i, j = j, chunkSize = chunkSize, nCores = nCores, verbose = verbose, ...)
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
        GWAS.model <- stats::update(formula, ".~z+.")
        OUT <- chunkedApply(X = geno(data), MARGIN = 2L, FUN = function(col, ...) {
            df <- pheno(data)[i, , drop = FALSE]
            df$z <- col
            fm <- FUN(GWAS.model, data = df, ...)
            getCoefficients(fm)
        }, i = i, j = j, chunkSize = chunkSize, nCores = nCores, verbose = verbose, ...)
        OUT <- t(OUT)
        rownames(OUT) <- colnames(geno(data))[j]
    }

    return(OUT)
}

GWAS.rayOLS <- function(formula, data, i = seq_len(nrow(geno(data))), j = seq_len(ncol(geno(data))), chunkSize = 5000L, nCores = getOption("mc.cores", 2L), verbose = FALSE, ...) {
    y <- pheno(data)[i, getResponse(formula)]
    y <- as.numeric(y)
    res <- chunkedMap(X = geno(data), FUN = rayOLS, i = i, j = j, chunkSize = chunkSize, nCores = nCores, verbose = verbose, y = y, ...)
    res <- do.call(rbind, res)
    colnames(res) <- c("Estimate", "Std.Err", "t-value", "Pr(>|t|)", "n")
    rownames(res) <- colnames(geno(data))[j]
    return(res)
}

GWAS.lsfit <- function(formula, data, i = seq_len(nrow(geno(data))), j = seq_len(ncol(geno(data))), chunkSize = 5000L, nCores = getOption("mc.cores", 2L), verbose = FALSE, ...) {

    # The subset argument of model.frame is evaluated in the environment of the
    # formula, therefore subset after building the frame.
    frame <- stats::model.frame(formula = formula, data = pheno(data), na.action = stats::na.pass)[i, , drop = FALSE]
    model <- stats::model.matrix(formula, frame)

    y <- pheno(data)[i, getResponse(formula)]

    res <- chunkedApply(X = geno(data), MARGIN = 2L, FUN = function(col, ...) {
        fm <- stats::lsfit(x = cbind(col, model), y = y, intercept = FALSE)
        stats::ls.print(fm, print.it = FALSE)$coef.table[[1L]][1L, ]
    }, i = i, j = j, chunkSize = chunkSize, nCores = nCores, verbose = verbose, ...)
    res <- t(res)
    rownames(res) <- colnames(geno(data))[j]

    return(res)
}

# formula: the formula for the GWAS model without including the markers, e.g.
# y~1 or y~factor(sex)+age
# all the variables in the formula must be in data@pheno (BGData)
# containing slots @pheno and @geno
# groups: a vector mapping markers into groups (can be integer, character or
# factor)
GWAS.SKAT <- function(formula, data, groups, i = seq_len(nrow(geno(data))), j = seq_len(ncol(geno(data))), verbose = FALSE, ...) {

    uniqueGroups <- unique(groups)

    OUT <- matrix(data = double(), nrow = length(uniqueGroups), ncol = 2L)
    colnames(OUT) <- c("nMrk", "p-value")
    rownames(OUT) <- uniqueGroups

    H0 <- SKAT::SKAT_Null_Model(formula, data = pheno(data)[i, , drop = FALSE], ...)

    for (group in seq_along(uniqueGroups)) {
        Z <- geno(data)[i, groups == uniqueGroups[group], drop = FALSE]
        fm <- SKAT::SKAT(Z = Z, obj = H0, ...)
        OUT[group, ] <- c(ncol(Z), fm$p.value)
        if (verbose) {
            message("Group ", group, " of ", length(uniqueGroups), " ...")
        }
    }

    return(OUT)
}

rayOLS <- function(x, y) {
    .Call(C_rayOLS, x, y)
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

getResponse <- function(formula) {
    # Extract component from parse tree (see https://cran.r-project.org/doc/manuals/r-release/R-lang.html#Language-objects)
    sym <- formula[[2L]]
    # Convert symbol to character
    as.character(sym)
}
