\name{GWAS}
\alias{GWAS}
\title{Performs Single Marker Regressions Using BGData Objects}
\description{
    Implements single marker regressions. The regression model includes all the
    covariates specified in the right-hand-side of the \code{formula} plus one
    column of \code{@geno} at a time. The data from the association tests is
    obtained from a \linkS4class{BGData} object.
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
        the formula must be in the \code{@pheno} object of the
        \linkS4class{BGData}.
    }
    \item{data}{
        A \linkS4class{BGData} object.
    }
    \item{method}{
        The regression method to be used. Currently, the following methods are
        implemented: \code{rayOLS}, \code{\link[stats:lsfit]{stats::lsfit()}},
        \code{\link[stats:lm]{stats::lm()}},
        \code{\link[stats:lm.fit]{stats::lm.fit()}},
        \code{\link[stats:glm]{stats::glm()}},
        \code{\link[lme4:lmer]{lme4::lmer()}}, and
        \code{\link[SKAT:SKAT]{SKAT::SKAT()}}. Defaults to \code{lsfit}.
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
        The number of cores (passed to
        \code{\link[parallel:mclapply]{parallel::mclapply()}}).  Defaults to
        the number of cores as detected by
        \code{\link[parallel:detectCores]{parallel::detectCores()}}.
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
\section{File-backed matrices}{
    Functions with the \code{chunkSize} parameter work best with file-backed
    matrices such as \link[BEDMatrix:BEDMatrix-class]{BEDMatrix::BEDMatrix}
    objects. To avoid loading the whole, potentially very large matrix into
    memory, these functions will load chunks of the file-backed matrix into
    memory and perform the operations on one chunk at a time. The size of the
    chunks is determined by the \code{chunkSize} parameter. Care must be taken
    to not set \code{chunkSize} too high to avoid memory shortage, particularly
    when combined with parallel computing.
}
\section{Multi-level parallelism}{
    Functions with the \code{nCores}, \code{i}, and \code{j} parameters provide
    capabilities for both parallel and distributed computing.

    For parallel computing, \code{nCores} determines the number of cores the
    code is run on. Memory usage can be an issue for higher values of
    \code{nCores} as R is not particularly memory-efficient. As a rule of
    thumb, at least around \code{(nCores * object_size(chunk)) +
    object_size(result)} MB of total memory will be needed for operations
    on file-backed matrices, not including potential copies of your data that
    might be created (for example \code{\link[stats:lsfit]{stats::lsfit()}}
    runs \code{cbind(1, X)}). \code{i} and \code{j} can be used to include or
    exclude certain rows or columns. Internally, the
    \code{\link[parallel:mclapply]{parallel::mclapply()}} function is used and
    therefore parallel computing will not work on Windows machines.

    For distributed computing, \code{i} and \code{j} determine the subset of
    the input matrix that the code runs on. In an HPC environment, this can be
    used not just to include or exclude certain rows or columns, but also to
    partition the task among many nodes rather than cores. Scheduler-specific
    code and code to aggregate the results need to be written by the user. It
    is recommended to set \code{nCores} to \code{1} as nodes are often cheaper
    than cores.
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
