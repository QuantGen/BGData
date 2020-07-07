\name{summarize}
\alias{summarize}
\title{Generates Various Summary Statistics}
\description{
    Computes the frequency of missing values, the (minor) allele frequency, and
    standard deviation of each column of \code{X}.
}
\usage{
summarize(X, i = seq_len(nrow(X)), j = seq_len(ncol(X)),
  chunkSize = 5000L, nCores = getOption("mc.cores", 2L),
  verbose = FALSE)
}
\arguments{
    \item{X}{
        A matrix-like object, typically the genotypes of a \code{BGData}
        object.
    }
    \item{i}{
        Indicates which rows of \code{X} should be used. Can be integer,
        boolean, or character. By default, all rows are used.
    }
    \item{j}{
        Indicates which columns of \code{X} should be used. Can be integer,
        boolean, or character. By default, all columns are used.
    }
    \item{chunkSize}{
        The number of columns of \code{X} that are brought into physical memory
        for processing per core. If \code{NULL}, all elements in \code{j} are
        used. Defaults to 5000.
    }
    \item{nCores}{
        The number of cores (passed to \code{mclapply}). Defaults to the
        number of cores as detected by \code{detectCores}.
    }
    \item{verbose}{
        Whether progress updates will be posted. Defaults to \code{FALSE}.
    }
}
\value{
    A \code{data.frame} with three columns: \code{freq_na} for frequencies of
    missing values, \code{allele_freq} for allele frequencies of the counted
    allele, and \code{sd} for standard deviations.
}
\seealso{
    \code{\link{file-backed-matrices}} for more information on file-backed
    matrices. \code{\link{multi-level-parallelism}} for more information on
    multi-level parallelism. \code{\link{BGData-class}} for more information on
    the \code{BGData} class.
}
\examples{
# Restrict number of cores to 1 on Windows
if (.Platform$OS.type == "windows") {
    options(mc.cores = 1)
}

# Load example data
bg <- BGData:::loadExample()

# Summarize the whole dataset
sum1 <- summarize(X = geno(bg))

# Summarize the first 50 individuals
sum2 <- summarize(X = geno(bg), i = 1:50)

# Summarize the first 1000 markers (useful for distributed computing)
sum3 <- summarize(X = geno(bg), j = 1:100)

# Summarize the first 50 individuals on the first 1000 markers
sum4 <- summarize(X = geno(bg), i = 1:50, j = 1:100)

# Summarize by names
sum5 <- summarize(X = geno(bg), j = c("snp81233_C", "snp81234_C", "snp81235_T"))
}
