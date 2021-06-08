\name{GWAS}
\alias{GWAS}
\title{Performs Single Marker Regressions Using BGData Objects}
\description{
    Implements single marker regressions. The regression model includes all the
    covariates specified in the right-hand-side of the \code{formula} plus one
    column of the genotypes at a time. The data from the association tests is
    obtained from a \code{BGData} object.
}
\usage{
GWAS(formula, data, method = "lsfit", i = seq_len(nrow(geno(data))),
  j = seq_len(ncol(geno(data))), chunkSize = 5000L,
  nCores = getOption("mc.cores", 2L), verbose = FALSE, ...)
}
\arguments{
    \item{formula}{
        The formula for the GWAS model without the variant, e.g. \code{y ~ 1}
        or \code{y ~ factor(sex) + age}. The variables included in the formula
        must be column names in the sample information of the \code{BGData}
        object.
    }
    \item{data}{
        A \code{BGData} object.
    }
    \item{method}{
        The regression method to be used. Currently, the following methods are
        implemented: \code{rayOLS} (see below), \code{lsfit}, \code{lm},
        \code{lm.fit}, \code{glm}, \code{lmer}, and \code{SKAT}. Defaults to
        \code{lsfit}.
    }
    \item{i}{
        Indicates which rows of the genotypes should be used. Can be integer,
        boolean, or character. By default, all rows are used.
    }
    \item{j}{
        Indicates which columns of the genotypes should be used. Can be
        integer, boolean, or character. By default, all columns are used.
    }
    \item{chunkSize}{
        The number of columns of the genotypes that are brought into physical
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
\details{
    The \code{rayOLS} method is a regression through the origin that can only
    be used with a \code{y ~ 1} formula, i.e. it only allows for one
    quantitative response variable \code{y} and one variant at a time as an
    explanatory variable (the variant is not included in the formula, hence
    \code{1} is used as a dummy). If covariates are needed, consider
    preadjustment of \code{y}. Among the provided methods, it is by far the
    fastest.

    Some regression methods may require the data to not contain columns with
    variance 0 or too many missing values. We suggest running \code{summarize}
    to detect variants that do not clear the desired minor-allele frequency and
    rate of missing genotype calls, and filtering these variants out using the
    \code{j} parameter of the \code{GWAS} function (see example below).
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

# Detect variants that do not pass MAF and missingness thresholds
summaries <- summarize(geno(bg))
maf <- ifelse(summaries$allele_freq > 0.5, 1 - summaries$allele_freq,
    summaries$allele_freq)
exclusions <- maf < 0.01 | summaries$freq_na > 0.05

# Perform a single marker regression
res1 <- GWAS(formula = FT10 ~ 1, data = bg, j = !exclusions)

# Draw a Manhattan plot
plot(-log10(res1[, 4]))

# Use lm instead of lsfit (the default)
res2 <- GWAS(formula = FT10 ~ 1, data = bg, method = "lm", j = !exclusions)

# Use glm instead of lsfit (the default)
y <- pheno(bg)$FT10
pheno(bg)$FT10.01 <- y > quantile(y, 0.8, na.rm = TRUE)
res3 <- GWAS(formula = FT10.01 ~ 1, data = bg, method = "glm", j = !exclusions)

# Perform a single marker regression on the first 50 markers (useful for
# distributed computing)
res4 <- GWAS(formula = FT10 ~ 1, data = bg, j = 1:50)
}
