FWD <- function(y, X, df = 20, tol = 1e-7, verbose = TRUE, maxIter = 1000, centerImpute = TRUE) {
    if (is.null(colnames(X))) {
        colnames(X) <- paste0("X", 1:ncol(X))
    }
    colNames <- colnames(X)
    y <- y - mean(y)
    if (centerImpute) {
        for (i in 1:ncol(X)) {
            X[, i] <- X[, i] - mean(X[, i], na.rm = TRUE)
            X[, i] <- ifelse(is.na(X[, i]), 0, X[, i])
        }
    }
    X <- cbind(1, X)
    df <- df + 1
    colnames(X) <- c("Int", colNames)
    colNames <- colnames(X)
    C <- crossprod(X)
    rhs <- crossprod(X,y)
    n <- length(y)
    p <- ncol(X)
    B <- matrix(nrow = p, ncol = df, 0)
    rownames(B) <- colNames
    B[1, 1] <- mean(y)
    RSS <- rep(NA, df)
    LogLik <- RSS
    VARE <- RSS
    AIC <- RSS
    DF <- RSS
    BIC <- RSS
    path <- rep(NA, df)
    RSS[1] <- sum((y - B[1, 1])^2)
    tol <- tol * RSS[1]
    DF[1] <- 1
    VARE[1] <- RSS[1] / (n - DF[1])
    LogLik[1] <- -(n / 2) * log(2 * pi * VARE[1]) - RSS[1] / (2 * VARE[1])
    AIC[1] <- -2 * LogLik[1] + 2 * DF[1]
    path[1] <- colNames[1]
    for (i in 2:df) {
        tmpB <- B[, i - 1]
        tmpRSS <- RSS[i - 1]
        tmp <- addOne(C, rhs, b = tmpB, RSS = tmpRSS, tol = tol, maxIter = maxIter)
        B[, i] <- tmp$b
        if (length(tmp$newPred) > 0) {
            path[i] <- colNames[tmp$newPred]
        } else {
            path[i] <- NA
        }
        RSS[i] <- tmp$RSS
        DF[i] <- sum(tmp$b != 0)
        VARE[i] <- RSS[i] / (n - DF[i])
        LogLik[i] <- -(n / 2) * log(2 * pi * VARE[i]) - RSS[i] / VARE[i] / 2
        AIC[i] <- -2 * LogLik[i] + 2 * (DF[i] + 1)
        BIC[i] <- -2 * LogLik[i] + log(n) * (DF[i] + 1)
        if (verbose) {
            message("  ", DF[i] - 1, " predictors, AIC=", round(AIC[i], 2))
        }
    }
    OUT <- list(B = B, path = path, RSS = RSS, LogLik = LogLik, VARE = VARE, DF = DF, AIC = AIC, BIC = BIC)
    return(OUT)
}

addOne <- function(C, rhs, RSS, b, tol = 1e-5, maxIter = 100) {
    notActive <- which(b == 0)
    active <- which(b != 0)
    q <- length(notActive)
    nActive <- length(active)
    b0 <- b + 0.0
    RSS0 <- RSS + 0.0
    # if model is null
    if (nActive == 0) {
        bOLS <- rhs / diag(C)
        dRSS <- diag(C) * bOLS^2
        k <- which.max(dRSS)
        b[k] <- bOLS[k]
        RSS <- RSS - bOLS^2 * C[k, k]
        ans <- list(b = b, newPred = notActive[k], RSS = RSS)
    # when model is not null
    } else {
        RSSNew <- rep(NA, q)
        for (i in 1:q) {
            RSS <- RSS0 + 0.0
            b <- b0 + 0.0
            fm <- fitSYS(C = C, rhs = rhs, RSS = RSS, b = b, tol = tol, maxIter = maxIter, active = c(notActive[i], active))
            RSSNew[i] <- fm$RSS
        }
        k <- which.min(RSSNew)
        b <- b0 + 0.0
        RSS <- RSS0 + 0.0
        fm <- fitSYS(C = C, rhs = rhs, RSS = RSS, b = b, tol = tol, maxIter = maxIter, active = c(notActive[k], active))
        ans <- list(b = fm$b, newPred = notActive[k], RSS = fm$RSS)
    }
    return(ans)
}

fitSYS <- function(C, rhs, b, active, RSS, tol, maxIter) {
    active <- active - 1 # for the 0-based index
    tmp_b <- b + 0.0
    tmp_RSS <- RSS + 0.0
    ans <- .Call("fitLSYS", ncol(C), length(active), C, rhs, tmp_b, active, tmp_RSS, maxIter, tol)
    return(list(b = ans[[1]], RSS = ans[[2]]))
}
