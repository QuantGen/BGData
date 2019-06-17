\name{GWAS}
\alias{GWAS}
\title{Performs Single Marker Regressions Using BGData Objects}
\description{
    Implements single marker regressions. The regression model includes all the
    covariates specified in the right-hand-side of the \code{formula} plus one
    column of \code{@geno} at a time. The data from the association tests is
    obtained from a \code{BGData} object.
}
\usage{
GWAS(formula, data, method = "lsfit", i = seq_len(nrow(data@geno)),
  j = seq_len(ncol(data@geno)), chunkSize = 5000L,
  nCores = getOption("mc.cores", 2L), verbose = FALSE, ...)
}
\arguments{
    \item{formula}{
        The formula for the GWAS model without including the marker, e.g.
        \code{y ~ 1} or \code{y ~ factor(sex) + age}. The variables included in
        the formula must be in the \code{@pheno} object of the \code{BGData}.
    }
    \item{data}{
        A \code{BGData} object.
    }
    \item{method}{
        The regression method to be used. Currently, the following methods are
        implemented: \code{rayOLS}, \code{lsfit}, \code{lm}, \code{lm.fit},
        \code{glm}, \code{lmer}, and \code{SKAT}. Defaults to \code{lsfit}.
    }
    \item{i}{
        Indicates which rows of \code{@geno} should be used. Can be integer,
        boolean, or character. By default, all rows are used.
    }
    \item{j}{
        Indicates which columns of \code{@geno} should be used. Can be integer,
        boolean, or character. By default, all columns are used.
    }
    \item{chunkSize}{
        The number of columns of \code{@geno} that are brought into physical
        memory for processing per core. If \code{NULL}, all elements in
        \code{j} are used. Defaults to 5000.
    }
    \item{nCores}{
        The number of cores (passed to \code{mclapply}). Defaults to the number
        of cores as detected by \code{detectCores}.
    }
    \item{verbose}{
        Whether progress updates will be posted. Defaults to \code{FALSE}.
    }
    \item{...}{
        Additional arguments for chunkedApply and regression method.
    }
}
\value{
    The same matrix that would be returned by \code{coef(summary(model))}.
}
\seealso{
    \code{\link{file-backed-matrices}} for more information on file-backed
    matrices. \code{\link{multi-level-parallelism}} for more information on
    multi-level parallelism. \code{\link{BGData-class}} for more information on
    the \code{BGData} class. \code{\link[stats]{lsfit}},
    \code{\link[stats]{lm}}, \code{\link[stats]{lm.fit}},
    \code{\link[stats]{glm}}, \code{\link[lme4]{lmer}}, and
    \code{\link[SKAT]{SKAT}} for more information on regression methods.
}
\examples{
# Restrict number of cores to 1 on Windows
if (.Platform$OS.type == "windows") {
    options(mc.cores = 1)
}

# Load example data
bg <- BGData:::loadExample()

# Perform a single marker regression
res1 <- GWAS(formula = FT10 ~ 1, data = bg)

# Draw a Manhattan plot
plot(-log10(res1[, 4]))

# Use lm instead of lsfit (the default)
res2 <- GWAS(formula = FT10 ~ 1, data = bg, method = "lm")

# Use glm instead of lsfit (the default)
y <- bg@pheno$FT10
bg@pheno$FT10.01 <- y > quantile(y, 0.8, na.rm = TRUE)
res3 <- GWAS(formula = FT10.01 ~ 1, data = bg, method = "glm")

# Perform a single marker regression on the first 50 markers (useful for
# distributed computing)
res4 <- GWAS(formula = FT10 ~ 1, data = bg, j = 1:50)
}
