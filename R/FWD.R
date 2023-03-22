FWD <- function(y, X, df = 20, tol = 1e-7, maxIter = 1000, centerImpute = TRUE, verbose = TRUE) {
    y <- y - mean(y)
    if (centerImpute) {
        X <- BGData::preprocess(X, center = TRUE, impute = TRUE)
    }
    if (is.null(colnames(X))) {
        colNames <- paste0("X", 1:ncol(X))
    } else {
        colNames <- colnames(X)
    }
    X <- cbind(1, X)
    df <- df + 1
    colNames <- c("Int", colNames)
    C <- crossprod(X)
    rhs <- crossprod(X, y)
    n <- length(y)
    p <- ncol(X)
    active <- rep(FALSE, p)
    names(active) <- colNames
    B <- matrix(data = 0, nrow = p, ncol = df)
    rownames(B) <- colNames
    RSS <- rep(NA_real_, df)
    DF <- rep(NA_real_, df)
    VARE <- rep(NA_real_, df)
    LogLik <- rep(NA_real_, df)
    AIC <- rep(NA_real_, df)
    BIC <- rep(NA_real_, df)
    path <- rep(NA_character_, df)
    active[1] <- TRUE
    B[1, 1] <- mean(y)
    RSS[1] <- sum((y - B[1, 1])^2)
    DF[1] <- 1
    VARE[1] <- RSS[1] / (n - DF[1])
    LogLik[1] <- -(n / 2) * log(2 * pi * VARE[1]) - RSS[1] / (2 * VARE[1])
    AIC[1] <- -2 * LogLik[1] + 2 * DF[1]
    BIC[1] <- -2 * LogLik[1] + log(n) * (DF[1] + 1)
    path[1] <- colNames[1]
    tol <- tol * RSS[1]
    for (i in 2:df) {
        tmp <- addOne(
            C = C,
            rhs = rhs,
            active = active,
            b = B[, i - 1],
            RSS = RSS[i - 1],
            maxIter = maxIter,
            tol = tol
        )
        B[, i] <- tmp[["b"]]
        if (length(tmp[["newPred"]]) > 0) {
            active[tmp[["newPred"]]] <- TRUE
            path[i] <- colNames[tmp[["newPred"]]]
        } else {
            path[i] <- NA
        }
        RSS[i] <- tmp[["RSS"]]
        DF[i] <- sum(active)
        VARE[i] <- RSS[i] / (n - DF[i])
        LogLik[i] <- -(n / 2) * log(2 * pi * VARE[i]) - RSS[i] / VARE[i] / 2
        AIC[i] <- -2 * LogLik[i] + 2 * (DF[i] + 1)
        BIC[i] <- -2 * LogLik[i] + log(n) * (DF[i] + 1)
        if (verbose) {
            message("  ", DF[i] - 1, " predictors, AIC=", round(AIC[i], 2))
        }
    }
    OUT <- list(
        B = B,
        path = data.frame(
            variable = path,
            RSS = RSS,
            LogLik = LogLik,
            VARE = VARE,
            DF = DF,
            AIC = AIC,
            BIC = BIC
        )
    )
    return(OUT)
}

addOne <- function(C, rhs, active, b, RSS, maxIter, tol) {
    activeSet <- which(active)
    inactiveSet <- which(!active)
    nActive <- length(activeSet)
    nInactive <- length(inactiveSet)
    # if model is not null
    if (nActive > 1) {
        RSSNew <- rep(NA_real_, nInactive)
        for (i in 1:nInactive) {
            fm <- fitSYS(
                C = C,
                rhs = rhs,
                b = b,
                active = c(inactiveSet[i], activeSet),
                RSS = RSS,
                maxIter = maxIter,
                tol = tol
            )
            RSSNew[i] <- fm[["RSS"]]
        }
        k <- which.min(RSSNew)
        fm <- fitSYS(
            C = C,
            rhs = rhs,
            b = b,
            active = c(inactiveSet[k], activeSet),
            RSS = RSS,
            maxIter = maxIter,
            tol = tol
        )
        ans <- list(b = fm[["b"]], newPred = inactiveSet[k], RSS = fm[["RSS"]])
    # if model is null
    } else {
        bOLS <- rhs / diag(C)
        dRSS <- diag(C) * bOLS^2
        k <- which.max(dRSS)
        b[k] <- bOLS[k]
        RSS <- RSS - bOLS[k]^2 * C[k, k]
        ans <- list(b = b, newPred = k, RSS = RSS)
    }
    return(ans)
}

fitSYS <- function(C, rhs, b, active, RSS, maxIter, tol) {
    active <- active - 1L # for the 0-based index
    ans <- .Call(C_fitLSYS, C, rhs, b, active, RSS, maxIter, tol)
    return(list(b = ans[[1]], RSS = ans[[2]]))
}
