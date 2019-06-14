\name{readRAW}
\alias{readRAW}
\alias{readRAW_matrix}
\alias{readRAW_big.matrix}
\title{Creates a BGData Object From a .raw File or a .ped-Like File}
\usage{
readRAW(fileIn, header = TRUE, dataType = integer(), n = NULL,
  p = NULL, sep = "", na.strings = "NA", nColSkip = 6L,
  idCol = c(1L, 2L), nNodes = NULL, linked.by = "rows",
  folderOut = paste0("BGData_", sub("\\\\.[[:alnum:]]+$", "",
  basename(fileIn))), outputType = "byte", dimorder = if (linked.by ==
  "rows") 2L:1L else 1L:2L, verbose = FALSE)

readRAW_matrix(fileIn, header = TRUE, dataType = integer(), n = NULL,
  p = NULL, sep = "", na.strings = "NA", nColSkip = 6L,
  idCol = c(1L, 2L), verbose = FALSE)

readRAW_big.matrix(fileIn, header = TRUE, dataType = integer(),
  n = NULL, p = NULL, sep = "", na.strings = "NA", nColSkip = 6L,
  idCol = c(1L, 2L), folderOut = paste0("BGData_",
  sub("\\\\.[[:alnum:]]+$", "", basename(fileIn))), outputType = "char",
  verbose = FALSE)
}
\arguments{
\item{fileIn}{The path to the plaintext file.}

\item{header}{Whether \code{fileIn} contains a header. Defaults to \code{TRUE}.}

\item{dataType}{The coding type of genotypes in \code{fileIn}. Use \code{integer()} or
\code{double()} for numeric coding. Alpha-numeric coding is currently not
supported for \code{\link[=readRAW]{readRAW()}} and \code{\link[=readRAW_big.matrix]{readRAW_big.matrix()}}: use the \code{--recodeA}
option of PLINK to convert the .ped file into a .raw file. Defaults to
\code{integer()}.}

\item{n}{The number of individuals. Auto-detect if \code{NULL}. Defaults to
\code{NULL}.}

\item{p}{The number of markers. Auto-detect if \code{NULL}. Defaults to \code{NULL}.}

\item{sep}{The field separator character. Values on each line of the file
are separated by this character. If \code{sep = ""} (the default for \code{\link[=readRAW]{readRAW()}}
the separator is "white space", that is one or more spaces, tabs, newlines
or carriage returns.}

\item{na.strings}{The character string used in the plaintext file to denote
missing value. Defaults to \code{NA}.}

\item{nColSkip}{The number of columns to be skipped to reach the genotype
information in the file. Defaults to \code{6}.}

\item{idCol}{The index of the ID column. If more than one index is given,
both columns will be concatenated with "_". Defaults to \code{c(1, 2)}, i.e. a
concatenation of the first two columns.}

\item{nNodes}{The number of nodes to create. Auto-detect if \code{NULL}. Defaults
to \code{NULL}.}

\item{linked.by}{If \code{columns} a column-linked matrix
(\link[LinkedMatrix:ColumnLinkedMatrix-class]{LinkedMatrix::ColumnLinkedMatrix}) is created, if \code{rows} a
row-linked matrix (\link[LinkedMatrix:RowLinkedMatrix-class]{LinkedMatrix::RowLinkedMatrix}). Defaults to
\code{rows}.}

\item{folderOut}{The path to the folder where to save the binary files.
Defaults to the name of the input file (\code{fileIn}) without extension prefixed
with "BGData_".}

\item{outputType}{The \code{vmode} for \code{ff} and \code{type} for
\link[bigmemory:big.matrix-class]{bigmemory::big.matrix}) objects. Default to \code{byte} for \code{ff} and
\code{char} for \link[bigmemory:big.matrix-class]{bigmemory::big.matrix} objects.}

\item{dimorder}{The physical layout of the underlying \code{ff} object of each
node.}

\item{verbose}{Whether progress updates will be posted. Defaults to \code{FALSE}.}
}
\description{
Creates a \linkS4class{BGData} object from a .raw file (generated with \code{--recodeA}
in \href{https://www.cog-genomics.org/plink2}{PLINK}). Other text-based file
formats are supported as well by tweaking some of the parameters as long as
the records of individuals are in rows, and phenotypes, covariates and
markers are in columns.
}
\details{
The data included in the first couple of columns (up to \code{nColSkip}) is used
to populate the \code{@pheno} slot of a \linkS4class{BGData} object, and the remaining
columns are used to fill the \code{@geno} slot. If the first row contains a
header (\code{header = TRUE}), data in this row is used to determine the column
names for \code{@pheno} and \code{@geno}.

\code{@geno} can take several forms, depending on the function that is called
(\code{readRAW}, \code{readRAW_matrix}, or \code{readRAW_big.matrix}). The following
sections illustrate each function in detail.
}
\section{readRAW}{

Genotypes are stored in a \link[LinkedMatrix:LinkedMatrix-class]{LinkedMatrix::LinkedMatrix} object where
each node is an \code{ff} instance. Multiple \code{ff} files are used because the
array size in \code{ff} is limited to the largest integer which can be
represented on the system (\code{.Machine$integer.max}) and for genetic data this
limitation is often exceeded. The \link[LinkedMatrix:LinkedMatrix-class]{LinkedMatrix::LinkedMatrix} package
makes it possible to link several \code{ff} files together by columns or by rows
and treat them similarly to a single matrix. By default a
\link[LinkedMatrix:ColumnLinkedMatrix-class]{LinkedMatrix::ColumnLinkedMatrix} is used for \code{@geno}, but the user
can modify this using the \code{linked.by} argument. The number of nodes to
generate is either specified by the user using the \code{nNodes} argument or
determined internally so that each \code{ff} object has a number of cells that is
smaller than \code{.Machine$integer.max / 1.2}. A folder (see \code{folderOut}) that
contains the binary flat files (named \code{geno_*.bin}) and an external
representation of the \linkS4class{BGData} object in \code{BGData.RData} is created.
}

\section{readRAW_matrix}{

Genotypes are stored in a regular \code{matrix} object. Therefore, this function
will only work if the .raw file is small enough to fit into memory.
}

\section{readRAW_big.matrix}{

Genotypes are stored in a filebacked \link[bigmemory:big.matrix-class]{bigmemory::big.matrix} object.
A folder (see \code{folderOut}) that contains the binary flat file (named
\code{BGData.bin}), a descriptor file (named \code{BGData.desc}), and an external
representation of the \linkS4class{BGData} object in \code{BGData.RData} are created.
}

\section{Reloading a BGData object}{

To reload a \linkS4class{BGData} object, it is recommended to use the
\code{\link[=load.BGData]{load.BGData()}} function instead of the \code{\link[base:load]{base::load()}} function as
\code{\link[base:load]{base::load()}} does not initialize \code{ff} objects or attach
\link[bigmemory:big.matrix-class]{bigmemory::big.matrix} objects.
}

\examples{
# Path to example data
path <- system.file("extdata", package = "BGData")

# Convert RAW files of chromosome 1 to a BGData object
bg <- readRAW(fileIn = paste0(path, "/chr1.raw"))
}
\seealso{
\code{\link[=load.BGData]{load.BGData()}} to load a previously saved \linkS4class{BGData} object,
\code{\link[=as.BGData]{as.BGData()}} to create \linkS4class{BGData} objects from non-text files (e.g. BED
files).
}
