# Convert ff_matrix into an S4 class
setOldClass("ff_matrix")


# Convert BEDMatrix into an S4 class
setOldClass("BEDMatrix")


setClassUnion("geno", c("LinkedMatrix", "BEDMatrix", "big.matrix", "ff_matrix", "matrix"))


#' An S4 class to represent phenotype and genotype data.
#'
#' This class is inspired by the phenotype/genotype file format PED and its
#' binary companion (also known as BED) of
#' \href{http://pngu.mgh.harvard.edu/~purcell/plink/}{PLINK}. It is used by
#' several functions of this package such as \code{\link{GWAS}} for performing a
#' Genome Wide Association Study or \code{\link{getG}} for calculating a genomic
#' relationship matrix.
#'
#' There are several ways to create an instance of this class:  \enumerate{
#' \item from a raw PED file (or a PED-like file) using \code{\link{readPED}},
#' \code{\link{readPED.matrix}}, or \code{\link{readPED.big.matrix}}. \item from
#' a BED file using \code{\link[=as.BGData.BEDMatrix]{as.BGData}}. \item from
#' arbitrary phenotype/genotype data using one of the constructors
#' \code{BGData(...)} or \code{new("BGData",...)}. \item from a previously saved
#' \code{BGData} object using \code{\link{load.BGData}}. \item from multiple
#' files (even a mixture of different file types) using
#' \code{\link[=LinkedMatrix-class]{LinkedMatrix}} }
#'
#' A PED file can be recoded to a raw PED in
#' \href{http://pngu.mgh.harvard.edu/~purcell/plink/}{PLINK} using \code{plink
#' --file myfile --recodeA}, or converted to a BED file using \code{plink --file
#' myfile --make-bed}. Conversely, a BED file can be transformed back to a PED
#' file using \code{plink --bfile myfile --recode} or to a raw PED file using
#' \code{plink --bfile myfile --recodeA} without losing information.
#'
#' @slot geno A \code{geno} object that contains genotypes. \code{geno} is a
#'   class union of several matrix-like objects, many of them suitable for very
#'   large datasets. Currently supported are
#'   \code{\link[=LinkedMatrix-class]{LinkedMatrix}}, \code{BEDMatrix},
#'   \code{\link[=big.matrix-class]{big.matrix}}, \code{ff_matrix}, and
#'   \code{matrix}.
#' @slot pheno A \code{\link{data.frame}} that contains phenotypes.
#' @slot map A \code{\link{data.frame}} that contains a genetic map.
#' @export BGData
#' @exportClass BGData
BGData <- setClass("BGData", slots = c(geno = "geno", pheno = "data.frame", map = "data.frame"))


#' Creates a new \code{\link[=BGData-class]{BGData}} instance.
#'
#' This method is run when a \code{\link[=BGData-class]{BGData}} object is
#' created using \code{BGData(...)} or \code{new("BGData",...)}.
#'
#' @param .Object The \code{\link[=BGData-class]{BGData}} instance to be
#'   initialized. This argument is passed in by R and can be ignored, but still
#'   needs to be documented.
#' @param geno A \code{geno} object that contains genotypes. \code{geno} is a
#'   class union of several matrix-like objects, many of them suitable for very
#'   large datasets. Currently supported are
#'   \code{\link[=LinkedMatrix-class]{LinkedMatrix}}, \code{BEDMatrix},
#'   \code{\link[=big.matrix-class]{big.matrix}}, \code{ff_matrix}, and
#'   \code{matrix}.
#' @param pheno A \code{\link{data.frame}} that contains phenotypes. A stub that
#'   only contains an \code{IID} column populated with the rownames of
#'   \code{@@geno} will be generated if missing.
#' @param map A \code{\link{data.frame}} that contains a genetic map. A stub
#'   that only contains a \code{mrk} column populated with the colnames of
#'   \code{@@geno} will be generated if missing.
#' @export
setMethod("initialize", "BGData", function(.Object, geno, pheno, map) {
    if (!is(geno, "geno")) {
        stop("Only LinkedMatrix, BEDMatrix, big.matrix, ff_matrix, or regular matrix objects are allowed for geno.")
    }
    if (is.null(colnames(geno))) {
        colnames(geno) <- paste0("mrk_", seq_len(ncol(geno)))
    }
    if (is.null(rownames(geno))) {
        rownames(geno) <- paste0("id_", seq_len(nrow(geno)))
    }
    if (missing(pheno)) {
        pheno <- data.frame(IID = rownames(geno), stringsAsFactors = FALSE)
    }
    if (missing(map)) {
        map <- data.frame(mrk = colnames(geno), stringsAsFactors = FALSE)
    }
    .Object@geno <- geno
    .Object@pheno <- pheno
    .Object@map <- map
    return(.Object)
})


pedDims <- function(fileIn, header, n, p, sep = "", nColSkip = 6) {
    if (is.null(n)) {
        n <- getLineCount(fileIn, header)
    }
    if (header) {
        headerLine <- getFileHeader(fileIn, sep)
        p <- length(headerLine) - nColSkip
    } else {
        if (is.null(p)) {
            p <- getColumnCount(fileIn, sep) - nColSkip
        }
    }
    return(list(n = n, p = p))
}


parsePED <- function(BGData, fileIn, header, dataType, nColSkip = 6, idCol = c(1, 2), sep = "", na.strings = "NA", verbose = FALSE, ...) {

    p <- ncol(BGData@geno)
    pedFile <- gzfile(fileIn, open = "r")

    # Update colnames
    if (header) {
        headerLine <- scan(pedFile, nlines = 1, what = character(), sep = sep, quiet = TRUE)
        colnames(BGData@pheno) <- headerLine[seq_len(nColSkip)]
        colnames(BGData@geno) <- headerLine[-(seq_len(nColSkip))]
    }

    # Parse file
    j <- seq_len(p)
    for (i in seq_len(nrow(BGData@geno))) {
        xSkip <- scan(pedFile, n = nColSkip, what = character(), sep = sep, quiet = TRUE)
        x <- scan(pedFile, n = p, what = dataType, sep = sep, na.strings = na.strings, quiet = TRUE)
        BGData@pheno[i, ] <- xSkip
        BGData@geno <- `[<-`(BGData@geno, i, j, ..., value = x)
        if (verbose) {
            message("Subject ", i, " / ", nrow(BGData@geno))
        }
    }
    close(pedFile)

    # Update rownames
    IDs <- apply(BGData@pheno[, idCol, drop = FALSE], 1, paste, collapse = "_")
    rownames(BGData@pheno) <- IDs
    rownames(BGData@geno) <- IDs

    # Convert types in pheno
    BGData@pheno[] <- lapply(BGData@pheno, utils::type.convert, as.is = TRUE)

    return(BGData)
}


#' Creates a \code{\link[=BGData-class]{BGData}} object from a raw PED file
#' (generated with \code{--recodeA} in PLINK) or a PED-like file and stores the
#' genotypes in a \code{\link[=LinkedMatrix-class]{LinkedMatrix}} that contains
#' memory-mapped \code{ff} nodes.
#'
#' \code{readPED} assumes that the plaintext file (\code{fileIn}) contains
#' records of individuals in rows, and phenotypes, covariates and markers in
#' columns. The columns included in the first couple of columns
#' (\code{seq_len(nColSkip)}) are used to populate the \code{@@pheno} slot of a
#' \code{\link[=BGData-class]{BGData}} object, and the remaining columns are
#' used to fill the \code{@@geno} slot. If the first row contains a header
#' (\code{header=TRUE}), data in this row is used to determine variables names
#' for \code{@@pheno} and marker names for \code{@@map} and \code{@@geno}.
#'
#' Genotypes are stored in a \code{\link[=LinkedMatrix-class]{LinkedMatrix}}
#' object where each node is an \code{ff} instance rather than a single
#' \code{ff} object. This is because the array size in \code{ff} is limited to
#' the largest integer which can be represented on the system
#' (\code{.Machine$integer.max}) and for genetic data this limitation is often
#' exceeded. The \code{\link[=LinkedMatrix-class]{LinkedMatrix}} package makes
#' it possible to link several \code{ff} files together by columns or by rows
#' and treat them similarly to a single matrix. By default a
#' \code{\link[=ColumnLinkedMatrix-class]{ColumnLinkedMatrix}} is used for
#' \code{@@geno}, but the user can modify this using the \code{linked.by}
#' argument. The number of nodes to generate is either specified by the user
#' using the \code{nNodes} argument or determined internally so that each
#' \code{ff} object has a number of cells that is smaller than
#' \code{.Machine$integer.max / 1.2}.
#'
#' \code{readPED} creates a folder (see \code{folderOut}) that contains the
#' binary flat files (named \code{geno_*.bin}) and an external representation of
#' the \code{\link[=BGData-class]{BGData}} object in \code{BGData.RData}. A
#' \code{\link[=BGData-class]{BGData}} object can be reloaded using
#' \code{load.BGData} (the regular \code{load} function will only work if the
#' working directory is set to the path that contains the binary flat files).
#'
#' @param fileIn The path to the plaintext file.
#' @param header If TRUE, the file contains a header.
#' @param dataType The coding of genotypes. Use \code{integer()} or
#'   \code{double()} for numeric coding. Character coding is currently not
#'   supported: use the \code{--recodeA} option of PLINK to convert the PED file
#'   into a raw file.
#' @param n The number of individuals.
#' @param p The number of markers.
#' @param sep The field separator character. Values on each line of the file
#'   are separated by this character. If \code{sep = ""} (the default for
#'   \code{readPED}) the separator is "white space", that is one or more spaces,
#'   tabs, newlines or carriage returns.
#' @param na.strings The character string used in the plaintext file to denote
#'   missing value.
#' @param nColSkip The number of columns to be skipped to reach the genotype
#'   information in the file.
#' @param idCol The index of the ID column. If more than one index is given,
#'   both columns will be concatenated with "_".
#' @param nNodes The number of nodes to create.
#' @param linked.by If \code{columns} a column-linked matrix
#'   (\code{\link[=ColumnLinkedMatrix-class]{ColumnLinkedMatrix}}) is created,
#'   if \code{rows} a row-linked matrix
#'   (\code{\link[=RowLinkedMatrix-class]{RowLinkedMatrix}}).
#' @param folderOut The path to the folder where to save the binary files.
#' @param dimorder The physical layout of the underlying \code{ff} object of
#'   each node.
#' @param verbose Whether progress updates will be posted. Defaults to
#'   \code{FALSE}.
#' @seealso \code{\link[=BGData-class]{BGData}},
#'   \code{\link[=LinkedMatrix-class]{LinkedMatrix}},
#'   \code{\link[=ColumnLinkedMatrix-class]{ColumnLinkedMatrix}},
#'   \code{\link[=RowLinkedMatrix-class]{RowLinkedMatrix}}, \code{\link[ff]{ff}}
#' @export
readPED <- function(fileIn, header, dataType, n = NULL, p = NULL, sep = "", na.strings = "NA", nColSkip = 6, idCol = c(1, 2), nNodes = NULL, linked.by = "rows", folderOut = paste0("BGData_", sub("\\.[[:alnum:]]+$", "", basename(fileIn))), dimorder = if (linked.by == "rows") 2:1 else 1:2, verbose = FALSE) {

    # Create output directory
    if (file.exists(folderOut)) {
        stop(paste("Output folder", folderOut, "already exists. Please move it or pick a different one."))
    }
    dir.create(folderOut)

    dims <- pedDims(fileIn = fileIn, header = header, n = n, p = p, sep = sep, nColSkip = nColSkip)

    # Determine number of nodes
    if (is.null(nNodes)) {
        if (linked.by == "columns") {
            chunkSize <- min(dims$p, floor(.Machine$integer.max / dims$n / 1.2))
            nNodes <- ceiling(dims$p / chunkSize)
        } else {
            chunkSize <- min(dims$n, floor(.Machine$integer.max / dims$p / 1.2))
            nNodes <- ceiling(dims$n / chunkSize)
        }
    } else {
        if (linked.by == "columns") {
            chunkSize <- ceiling(dims$p / nNodes)
            if (chunkSize * dims$n >= .Machine$integer.max / 1.2) {
              stop("More nodes are needed")
            }
        } else {
            chunkSize <- ceiling(dims$n / nNodes)
            if (chunkSize * dims$p >= .Machine$integer.max / 1.2) {
              stop("More nodes are needed")
            }
        }
    }

    dataType <- normalizeType(dataType)
    if (!typeof(dataType) %in% c("integer", "double")) {
        stop("dataType must be either integer() or double()")
    }
    if (!linked.by %in% c("columns", "rows")) {
        stop("linked.by must be either columns or rows")
    }

    vmode <- ifelse(typeof(dataType) == "integer", "byte", "double")

    # Prepare geno
    geno <- LinkedMatrix::LinkedMatrix(nrow = dims$n, ncol = dims$p, nNodes = nNodes, linkedBy = linked.by, nodeInitializer = ffNodeInitializer, vmode = vmode, folderOut = folderOut, dimorder = dimorder)

    # Generate nodes
    nodes <- LinkedMatrix::nodes(geno)

    # Generate index
    index <- LinkedMatrix::index(geno)

    # Prepare pheno
    pheno <- as.data.frame(matrix(nrow = dims$n, ncol = nColSkip), stringsAsFactors = FALSE)

    # Construct BGData object
    BGData <- new("BGData", geno = geno, pheno = pheno)

    # Parse PED file
    BGData <- parsePED(BGData = BGData, fileIn = fileIn, header = header, dataType = dataType, nColSkip = nColSkip, idCol = idCol, sep = sep, na.strings = na.strings, nodes = nodes, index = index, verbose = verbose)

    # Save BGData object
    attr(BGData, "origFile") <- list(path = fileIn, dataType = typeof(dataType))
    attr(BGData, "dateCreated") <- date()
    save(BGData, file = paste0(folderOut, "/BGData.RData"))

    return(BGData)
}


#' Creates an \code{\link[=BGData-class]{BGData}} object from a raw PED file
#' (generated with \code{--recodeA} in PLINK) or a PED-like file and stores the
#' genotypes in an in-memory \code{matrix}.
#'
#' \code{readPED.matrix} assumes that the plaintext file (\code{fileIn})
#' contains records of individuals in rows, and phenotypes, covariates and
#' markers in columns. The columns included in columns \code{seq_len(nColSkip)}
#' are used to populate the slot \code{@@pheno} of a
#' \code{\link[=BGData-class]{BGData}} object, and the remaining columns are
#' used to fill the slot \code{@@geno}. If the first row contains a header
#' (\code{header=TRUE}), data in this row is used to determine variables names
#' for \code{@@pheno} and marker names for \code{@@map} and \code{@@geno}.
#'
#' Genotypes are stored in a regular \code{matrix} object. Therefore, this
#' function will only work if the raw PED file is small enough to fit into
#' memory.
#'
#' @param fileIn The path to the plaintext file.
#' @param header If TRUE, the file contains a header.
#' @param dataType The coding of genotypes. Use \code{character()} for A/C/G/T
#'   or \code{integer()} for numeric coding.
#' @param n The number of individuals.
#' @param p The number of markers.
#' @param sep The field separator character. Values on each line of the file
#'   are separated by this character. If \code{sep = ""} (the default for
#'   \code{readPED.matrix}) the separator is "white space", that is one or more
#'   spaces, tabs, newlines or carriage returns.
#' @param na.strings The character string used in the plaintext file to denote
#'   missing value.
#' @param nColSkip The number of columns to be skipped to reach the genotype
#'   information in the file.
#' @param idCol The index of the ID column. If more than one index is given,
#'   both columns will be concatenated with "_".
#' @param verbose Whether progress updates will be posted. Defaults to
#'   \code{FALSE}.
#' @return Returns a \code{\link[=BGData-class]{BGData}} object.
#' @seealso \code{\link[=BGData-class]{BGData}}
#' @export
readPED.matrix <- function(fileIn, header, dataType, n = NULL, p = NULL, sep = "", na.strings = "NA", nColSkip = 6, idCol = c(1, 2), verbose = FALSE) {

    dims <- pedDims(fileIn = fileIn, header = header, n = n, p = p, sep = sep, nColSkip = nColSkip)

    dataType <- normalizeType(dataType)

    # Prepare geno
    geno <- matrix(nrow = dims$n, ncol = dims$p)

    # Prepare pheno
    pheno <- as.data.frame(matrix(nrow = dims$n, ncol = nColSkip), stringsAsFactors = FALSE)

    # Construct BGData object
    BGData <- new("BGData", geno = geno, pheno = pheno)

    # Parse PED file
    BGData <- parsePED(BGData = BGData, fileIn = fileIn, header = header, dataType = dataType, nColSkip = nColSkip, idCol = idCol, sep = sep, na.strings = na.strings, verbose = verbose)

    return(BGData)
}


#' Creates a \code{\link[=BGData-class]{BGData}} object from a raw PED file
#' (generated with \code{--recodeA} in PLINK) or a PED-like file and stores the
#' genotypes in a \code{\link[=big.matrix-class]{big.matrix}}.
#'
#' \code{readPED.matrix} assumes that the plaintext file (\code{fileIn})
#' contains records of individuals in rows, and phenotypes, covariates and
#' markers in columns. The columns included in columns \code{seq_len(nColSkip)} are
#' used to populate the slot \code{@@pheno} of a
#' \code{\link[=BGData-class]{BGData}} object, and the remaining columns are
#' used to fill the slot \code{@@geno}. If the first row contains a header
#' (\code{header=TRUE}), data in this row is used to determine variables names
#' for \code{@@pheno} and marker names for \code{@@map} and \code{@@geno}.
#'
#' Genotypes are stored in a filebacked
#' \code{\link[=big.matrix-class]{big.matrix}} object. \code{readPED.big.matrix}
#' creates a folder (see \code{folderOut}) that contains the binary flat file
#' (named \code{BGData.bin}), a descriptor file (named \code{BGData.desc}), and
#' an external representation of the \code{\link[=BGData-class]{BGData}} object
#' in \code{BGData.RData}. A \code{\link[=BGData-class]{BGData}} object can be
#' reloaded using \code{load.BGData} (the regular \code{load} function will not
#' work unless the \code{\link[=big.matrix-class]{big.matrix}} instance is
#' manually attached using \code{\link[bigmemory]{attach.big.matrix}}).
#'
#' @param fileIn The path to the plaintext file.
#' @param header If TRUE, the file contains a header.
#' @param dataType The coding of genotypes. Use \code{character()} for A/C/G/T
#'   or \code{integer()} for numeric coding.
#' @param n The number of individuals.
#' @param p The number of markers.
#' @param sep The field separator character. Values on each line of the file
#'   are separated by this character. If \code{sep = ""} (the default for
#'   \code{readPED.big.matrix}) the separator is "white space", that is one or
#'   more spaces, tabs, newlines or carriage returns.
#' @param na.strings The character string used in the plaintext file to denote
#'   missing value.
#' @param nColSkip The number of columns to be skipped to reach the genotype
#'   information in the file.
#' @param idCol The index of the ID column. If more than one index is given,
#'   both columns will be concatenated with "_".
#' @param folderOut The path to the folder where to save the binary files.
#' @param verbose Whether progress updates will be posted. Defaults to
#'   \code{FALSE}.
#' @return Returns a \code{\link[=BGData-class]{BGData}} object.
#' @seealso \code{\link[=BGData-class]{BGData}}
#' @export
readPED.big.matrix <- function(fileIn, header, dataType, n = NULL, p = NULL, sep = "", na.strings = "NA", nColSkip = 6, idCol = c(1, 2), folderOut = paste0("BGData_", sub("\\.[[:alnum:]]+$", "", basename(fileIn))), verbose = FALSE) {

    if (file.exists(folderOut)) {
        stop(paste("Output folder", folderOut, "already exists. Please move it or pick a different one."))
    }

    dataType <- normalizeType(dataType)
    if (typeof(dataType) == "double") {
        type <- "double"
    } else if (typeof(dataType) == "integer") {
        type <- "char"
    } else {
        stop("dataType must be either integer() or double()")
    }

    dims <- pedDims(fileIn = fileIn, header = header, n = n, p = p, sep = sep, nColSkip = nColSkip)

    options(bigmemory.typecast.warning = FALSE)
    options(bigmemory.allow.dimnames = TRUE)

    # Create output directory
    dir.create(folderOut)

    # Prepare geno
    geno <- bigmemory::filebacked.big.matrix(nrow = dims$n, ncol = dims$p, type = type, backingpath = folderOut, backingfile = "BGData.bin", descriptorfile = "BGData.desc")

    # Prepare pheno
    pheno <- as.data.frame(matrix(nrow = dims$n, ncol = nColSkip), stringsAsFactors = FALSE)

    # Construct BGData object
    BGData <- new("BGData", geno = geno, pheno = pheno)

    # Parse PED file
    BGData <- parsePED(BGData = BGData, fileIn = fileIn, header = header, dataType = dataType, nColSkip = nColSkip, idCol = idCol, sep = sep, na.strings = na.strings, verbose = verbose)

    # Save BGData object
    attr(BGData, "origFile") <- list(path = fileIn, dataType = typeof(dataType))
    attr(BGData, "dateCreated") <- date()
    save(BGData, file = paste0(folderOut, "/BGData.RData"))

    return(BGData)
}


#' S3 generic to convert other objects into \code{\link[=BGData-class]{BGData}}
#' objects.
#'
#' @param x An object.
#' @param ... Additional arguments (see method).
#' @return A BGData object.
#' @export
as.BGData <- function(x, ...) {
    UseMethod("as.BGData")
}


#' Converts a \code{BEDMatrix} object to a \code{\link[=BGData-class]{BGData}}
#' object.
#'
#' If a FAM file (which corresponds to the first six columns of a PED file) of
#' the same name and in the same directory as the BED file exists, the
#' \code{@@pheno} slot will be populated with the data stored in that file.
#' Otherwise a stub that only contains an \code{IID} column populated with the
#' rownames of \code{@@geno} will be generated.
#'
#' PED only allows for a single phenotype. If more phenotypes are required it is
#' possible to store them in an
#' \href{http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#pheno}{alternate
#' phenotype file}. The path to such a file can be provided with
#' \code{alternatePhenotypeFile} and will be merged with the data in the
#' \code{@@pheno} slot.
#'
#' @param x A \code{BEDMatrix} object.
#' @param alternatePhenotypeFile Path to an
#'   \href{http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#pheno}{alternate
#'    phenotype file}.
#' @param ... Additional arguments to the \code{read.table} or \code{fread}
#'   call (if data.table package is installed) call to parse the alternate pheno
#'   file.
#' @return A \code{\link[=BGData-class]{BGData}} object.
#' @export
as.BGData.BEDMatrix <- function(x, alternatePhenotypeFile = NULL, ...) {
    # Path to BED file
    bedPath <- attr(x, "path")
    # Path to FAM file
    famPath <- sub(".bed", ".fam", bedPath)
    # Path to BIM file
    bimPath <- sub(".bed", ".bim", bedPath)
    # Read in pheno file
    if (file.exists(famPath)) {
        message("Extracting phenotypes from FAM file...")
        if (requireNamespace("data.table", quietly = TRUE)) {
            pheno <- data.table::fread(famPath, col.names = c(
                "FID",
                "IID",
                "PAT",
                "MAT",
                "SEX",
                "PHENOTYPE"
            ), data.table = FALSE, showProgress = FALSE)
        } else {
            pheno <- utils::read.table(famPath, col.names = c(
                "FID",
                "IID",
                "PAT",
                "MAT",
                "SEX",
                "PHENOTYPE"
            ), stringsAsFactors = FALSE)
        }
    } else {
        splits <- strsplit(rownames(x), "_")
        pheno <- data.frame(FID = sapply(splits, "[", 1), IID = sapply(splits, "[", 2), stringsAsFactors = FALSE)
    }
    # Read in map file
    if (file.exists(bimPath)) {
        message("Extracting map from BIM file...")
        if (requireNamespace("data.table", quietly = TRUE)) {
            map <- data.table::fread(bimPath, col.names = c(
                "chromosome",
                "snp_id",
                "genetic_distance",
                "base_pair_position",
                "allele_1",
                "allele_2"
            ), data.table = FALSE, showProgress = FALSE)
        } else {
            map <- utils::read.table(bimPath, col.names = c(
                "chromosome",
                "snp_id",
                "genetic_distance",
                "base_pair_position",
                "allele_1",
                "allele_2"
            ), stringsAsFactors = FALSE)
        }
    } else {
        splits <- strsplit(colnames(x), "_")
        map <- data.frame(bimPath, snp_id = sapply(splits, function(x) {
            paste0(x[seq_len(length(x) - 1)], collapse = "_")
        }), allele_1 = sapply(splits, function(x) {
            x[length(x)]
        }), stringsAsFactors = FALSE)
    }
    if (!is.null(alternatePhenotypeFile)) {
        if (!file.exists(alternatePhenotypeFile)) {
            stop("Alternate phenotype file does not exist.")
        } else {
            message("Merging alternate phenotype file...")
            if (requireNamespace("data.table", quietly = TRUE)) {
                alternatePhenotypes <- data.table::fread(alternatePhenotypeFile, data.table = FALSE, showProgress = FALSE, ...)
            } else {
                # Check if the file has a header, i.e. if the first row starts with
                # an FID and an IID entry (unfortunately, treating the first row as
                # colnames after read.table screws with the types)
                hasHeader = FALSE
                if (grepl("FID\\s+IID", readLines(alternatePhenotypeFile, n = 1))) {
                    hasHeader = TRUE
                }
                alternatePhenotypes <- utils::read.table(alternatePhenotypeFile, header = hasHeader, stringsAsFactors = FALSE, ...)
            }
            # Add artificial sort column to preserve order after merging
            # (merge's `sort = FALSE` order is unspecified)
            pheno$.sortColumn <- seq_len(nrow(pheno))
            # Merge phenotypes and alternate phenotypes
            pheno <- merge(pheno, alternatePhenotypes, by = c("FID", "IID"), all.x = TRUE)
            # Reorder phenotypes to match original order and delete artificial
            # column
            pheno <- pheno[order(pheno$.sortColumn), ]
            pheno <- pheno[, names(pheno) != ".sortColumn"]
        }
    }
    BGData(geno = x, pheno = pheno, map = map)
}


#' Loads \code{\link[=BGData-class]{BGData}} (and other) objects from .RData
#' files.
#'
#' This function is similar to \code{load}, but also initializes the different
#' types of objects that the \code{@@geno} slot of a
#' \code{\link[=BGData-class]{BGData}} object can take. Currently supported are
#' \code{ff_matrix}, \code{\link[=big.matrix-class]{big.matrix}}, and
#' \code{BEDMatrix} objects. If the object is of type
#' \code{\link[=LinkedMatrix-class]{LinkedMatrix}}, all nodes will be
#' initialized with their appropriate method.
#'
#' @param file The name of the .RData file to be loaded.
#' @param envir The environment where to load the data.
#' @export
load.BGData <- function(file, envir = parent.frame()) {
    # Load data into new environment
    loadingEnv <- new.env()
    load(file = file, envir = loadingEnv)
    names <- ls(envir = loadingEnv)
    for (name in names) {
        object <- get(name, envir = loadingEnv)
        # Load genotypes of BGData objects
        if (class(object) == "BGData") {
            object@geno <- loadGeno(object@geno, path = dirname(file))
        }
        # Assign object to envir
        assign(name, object, envir = envir)
    }
    message("Loaded objects: ", paste0(names, collapse = ", "))
}


loadGeno <- function(x, ...) {
    UseMethod("loadGeno")
}


loadGeno.LinkedMatrix <- function(x, path, ...) {
    for (i in seq_len(LinkedMatrix::nNodes(x))) {
        x[[i]] <- loadGeno(x[[i]], path = path)
    }
    return(x)
}


# Absolute paths to ff files are not stored, so the ff objects have to be
# loaded from the same directory as the RData file.
loadGeno.ff_matrix <- function(x, path, ...) {
    # Store current working directory and set working directory to path
    cwd <- getwd()
    setwd(path)
    # Open ff object
    ff::open.ff(x)
    # Restore the working directory
    setwd(cwd)
    return(x)
}


loadGeno.big.matrix <- function(x, path, ...) {
    return(bigmemory::attach.big.matrix(paste0(path, .Platform$file.sep, "BGData.desc")))
}


loadGeno.BEDMatrix <- function(x, ...) {
    dnames <- attr(x, "dnames")
    dims <- attr(x, "dims")
    path <- attr(x, "path")
    x <- BEDMatrix::BEDMatrix(path = path, n = dims[1], p = dims[2])
    dimnames(x) <- dnames
    return(x)
}


loadGeno.default <- function(x, ...) {
    return(x)
}


ffNodeInitializer <- function(nodeIndex, nrow, ncol, vmode, folderOut, ...) {
    filename <- paste0("geno_", nodeIndex, ".bin")
    node <- ff::ff(dim = c(nrow, ncol), vmode = vmode, filename = paste0(folderOut, .Platform$file.sep, filename), ...)
    # Change ff path to a relative one
    bit::physical(node)$pattern <- "ff"
    bit::physical(node)$filename <- filename
    return(node)
}
