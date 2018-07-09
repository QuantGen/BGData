# Convert ff_matrix into an S4 class
setOldClass("ff_matrix")


#' An Abstract S4 Class Union of Matrix-Like Types.
#'
#' [geno-class] is a class union of several matrix-like types, many of them
#' suitable for very large datasets. Currently supported are
#' [LinkedMatrix::LinkedMatrix-class], [BEDMatrix::BEDMatrix-class],
#' [bigmemory::big.matrix-class], `ff_matrix`, and `matrix`.
#'
#' @seealso The `@@geno` slot of [BGData-class] that accepts [geno-class]
#' objects.
setClassUnion("geno", c("LinkedMatrix", "BEDMatrix", "big.matrix", "ff_matrix", "matrix"))


#' An S4 Class to Represent Phenotype and Genotype Data.
#'
#' This class is inspired by the phenotype/genotype file format .raw and its
#' binary companion (also known as .bed) of
#' [PLINK](https://www.cog-genomics.org/plink2). It is used by several
#' functions of this package such as [GWAS()] for performing a Genome Wide
#' Association Study or [getG()] for calculating a genomic relationship matrix.
#'
#' There are several ways to create an instance of this class:
#' * from arbitrary phenotype/genotype data using one of the constructors
#' `[BGData(...)][initialize,BGData-method]` or `[new("BGData",
#' ...)][initialize,BGData-method]`.
#' * from a BED file using [as.BGData()].
#' * from a previously saved [BGData-class] object using [load.BGData()].
#' * from multiple files (even a mixture of different file types) using
#' [LinkedMatrix::LinkedMatrix-class].
#' * from a .raw file (or a .ped-like file) using [readRAW()],
#' [readRAW_matrix()], or [readRAW_big.matrix()].
#'
#' A .ped file can be recoded to a .raw file in
#' [PLINK](https://www.cog-genomics.org/plink2) using `plink --file myfile
#' --recodeA`, or converted to a BED file using `plink --file myfile
#' --make-bed`. Conversely, a BED file can be transformed back to a .ped file
#' using `plink --bfile myfile --recode` or to a .raw file using `plink --bfile
#' myfile --recodeA` without losing information.
#'
#' @slot geno A [geno-class] object that contains genotypes. [geno-class] is a
#' class union of several matrix-like types, many of them suitable for very
#' large datasets. Currently supported are [LinkedMatrix::LinkedMatrix-class],
#' [BEDMatrix::BEDMatrix-class], [bigmemory::big.matrix-class], `ff_matrix`,
#' and `matrix`.
#' @slot pheno A `data.frame` that contains phenotypes.
#' @slot map A `data.frame` that contains a genetic map.
#' @example man/examples/BGData.R
#' @export BGData
#' @exportClass BGData
BGData <- setClass("BGData", slots = c(geno = "geno", pheno = "data.frame", map = "data.frame"))


#' Creates a New BGData Instance.
#'
#' This method is run when a [BGData-class] object is created using
#' `BGData(...)` or `new("BGData", ...)`.
#'
#' @param .Object The [BGData-class] instance to be initialized. This argument
#' is passed in by R and can be ignored, but still needs to be documented.
#' @param geno A [geno-class] object that contains genotypes. [geno-class] is a
#' class union of several matrix-like types, many of them suitable for very
#' large datasets. Currently supported are [LinkedMatrix::LinkedMatrix-class],
#' [BEDMatrix::BEDMatrix-class], [bigmemory::big.matrix-class], `ff_matrix`,
#' and `matrix`.
#' @param pheno A `data.frame` that contains phenotypes. A stub that only
#' contains an `IID` column populated with the rownames of `@@geno` will be
#' generated if missing.
#' @param map A `data.frame` that contains a genetic map. A stub that only
#' contains a `mrk` column populated with the colnames of `@@geno` will be
#' generated if missing.
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


pedDims <- function(fileIn, header, n, p, sep = "", nColSkip = 6L) {
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


parseRAW <- function(BGData, fileIn, header, dataType, nColSkip = 6L, idCol = c(1L, 2L), sep = "", na.strings = "NA", verbose = FALSE) {

    p <- ncol(BGData@geno)
    pedFile <- file(fileIn, open = "r")

    # Update colnames
    if (header) {
        headerLine <- scan(pedFile, nlines = 1L, what = character(), sep = sep, quiet = TRUE)
        colnames(BGData@pheno) <- headerLine[seq_len(nColSkip)]
        colnames(BGData@geno) <- headerLine[-(seq_len(nColSkip))]
    }

    # Parse file
    for (i in seq_len(nrow(BGData@geno))) {
        xSkip <- scan(pedFile, n = nColSkip, what = character(), sep = sep, quiet = TRUE)
        x <- scan(pedFile, n = p, what = dataType, sep = sep, na.strings = na.strings, quiet = TRUE)
        BGData@pheno[i, ] <- xSkip
        BGData@geno[i, ] <- x
        if (verbose) {
            message("Subject ", i, " / ", nrow(BGData@geno))
        }
    }
    close(pedFile)

    # Update rownames
    IDs <- apply(BGData@pheno[, idCol, drop = FALSE], 1L, paste, collapse = "_")
    rownames(BGData@pheno) <- IDs
    rownames(BGData@geno) <- IDs

    # Convert types in pheno
    BGData@pheno[] <- lapply(BGData@pheno, utils::type.convert, as.is = TRUE)

    return(BGData)
}


#' Creates a BGData Object From a .raw File or a .ped-Like File.
#'
#' Creates a [BGData-class] object from a .raw file (generated with `--recodeA`
#' in [PLINK](https://www.cog-genomics.org/plink2)). Other text-based file
#' formats are supported as well by tweaking some of the parameters as long as
#' the records of individuals are in rows, and phenotypes, covariates and
#' markers are in columns.
#'
#' The data included in the first couple of columns (up to `nColSkip`) is used
#' to populate the `@@pheno` slot of a [BGData-class] object, and the remaining
#' columns are used to fill the `@@geno` slot. If the first row contains a
#' header (`header = TRUE`), data in this row is used to determine the column
#' names for `@@pheno` and `@@geno`.
#'
#' `@@geno` can take several forms, depending on the function that is called
#' (`readRAW`, `readRAW_matrix`, or `readRAW_big.matrix`). The following
#' sections illustrate each function in detail.
#'
#' @section readRAW:
#' Genotypes are stored in a [LinkedMatrix::LinkedMatrix-class] object where
#' each node is an `ff` instance. Multiple `ff` files are used because the
#' array size in `ff` is limited to the largest integer which can be
#' represented on the system (`.Machine$integer.max`) and for genetic data this
#' limitation is often exceeded. The [LinkedMatrix::LinkedMatrix-class] package
#' makes it possible to link several `ff` files together by columns or by rows
#' and treat them similarly to a single matrix. By default a
#' [LinkedMatrix::ColumnLinkedMatrix-class] is used for `@@geno`, but the user
#' can modify this using the `linked.by` argument. The number of nodes to
#' generate is either specified by the user using the `nNodes` argument or
#' determined internally so that each `ff` object has a number of cells that is
#' smaller than `.Machine$integer.max / 1.2`. A folder (see `folderOut`) that
#' contains the binary flat files (named `geno_*.bin`) and an external
#' representation of the [BGData-class] object in `BGData.RData` is created.
#'
#' @section readRAW_matrix:
#' Genotypes are stored in a regular `matrix` object. Therefore, this function
#' will only work if the .raw file is small enough to fit into memory.
#'
#' @section readRAW_big.matrix:
#' Genotypes are stored in a filebacked [bigmemory::big.matrix-class] object.
#' A folder (see `folderOut`) that contains the binary flat file (named
#' `BGData.bin`), a descriptor file (named `BGData.desc`), and an external
#' representation of the [BGData-class] object in `BGData.RData` are created.
#'
#' @section Reloading a BGData object:
#' To reload a [BGData-class] object, it is recommended to use the
#' [load.BGData()] function instead of the [base::load()] function as
#' [base::load()] does not initialize `ff` objects or attach
#' [bigmemory::big.matrix-class] objects.
#'
#' @param fileIn The path to the plaintext file.
#' @param header Whether `fileIn` contains a header. Defaults to `TRUE`.
#' @param dataType The coding type of genotypes in `fileIn`. Use `integer()` or
#' `double()` for numeric coding. Alpha-numeric coding is currently not
#' supported for [readRAW()] and [readRAW_big.matrix()]: use the `--recodeA`
#' option of PLINK to convert the .ped file into a .raw file. Defaults to
#' `integer()`.
#' @param n The number of individuals. Auto-detect if `NULL`. Defaults to
#' `NULL`.
#' @param p The number of markers. Auto-detect if `NULL`. Defaults to `NULL`.
#' @param sep The field separator character. Values on each line of the file
#' are separated by this character. If `sep = ""` (the default for [readRAW()]
#' the separator is "white space", that is one or more spaces, tabs, newlines
#' or carriage returns.
#' @param na.strings The character string used in the plaintext file to denote
#' missing value. Defaults to `NA`.
#' @param nColSkip The number of columns to be skipped to reach the genotype
#' information in the file. Defaults to `6`.
#' @param idCol The index of the ID column. If more than one index is given,
#' both columns will be concatenated with "_". Defaults to `c(1, 2)`, i.e. a
#' concatenation of the first two columns.
#' @param nNodes The number of nodes to create. Auto-detect if `NULL`. Defaults
#' to `NULL`.
#' @param linked.by If `columns` a column-linked matrix
#' ([LinkedMatrix::ColumnLinkedMatrix-class]) is created, if `rows` a
#' row-linked matrix ([LinkedMatrix::RowLinkedMatrix-class]). Defaults to
#' `rows`.
#' @param folderOut The path to the folder where to save the binary files.
#' Defaults to the name of the input file (`fileIn`) without extension prefixed
#' with "BGData_".
#' @param outputType The `vmode` for `ff` and `type` for
#' [bigmemory::big.matrix-class]) objects. Default to `byte` for `ff` and
#' `char` for [bigmemory::big.matrix-class] objects.
#' @param dimorder The physical layout of the underlying `ff` object of each
#' node.
#' @param verbose Whether progress updates will be posted. Defaults to `FALSE`.
#' @seealso [load.BGData()] to load a previously saved [BGData-class] object,
#' [as.BGData()] to create [BGData-class] objects from non-text files (e.g. BED
#' files).
#' @example man/examples/readRAW.R
#' @export
readRAW <- function(fileIn, header = TRUE, dataType = integer(), n = NULL, p = NULL, sep = "", na.strings = "NA", nColSkip = 6L, idCol = c(1L, 2L), nNodes = NULL, linked.by = "rows", folderOut = paste0("BGData_", sub("\\.[[:alnum:]]+$", "", basename(fileIn))), outputType = "byte", dimorder = if (linked.by == "rows") 2L:1L else 1L:2L, verbose = FALSE) {

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

    # Prepare geno
    geno <- LinkedMatrix::LinkedMatrix(nrow = dims$n, ncol = dims$p, nNodes = nNodes, linkedBy = linked.by, nodeInitializer = ffNodeInitializer, vmode = outputType, folderOut = folderOut, dimorder = dimorder)

    # Prepare pheno
    pheno <- as.data.frame(matrix(nrow = dims$n, ncol = nColSkip), stringsAsFactors = FALSE)

    # Construct BGData object
    BGData <- new("BGData", geno = geno, pheno = pheno)

    # Parse .raw file
    BGData <- parseRAW(BGData = BGData, fileIn = fileIn, header = header, dataType = dataType, nColSkip = nColSkip, idCol = idCol, sep = sep, na.strings = na.strings, verbose = verbose)

    # Save BGData object
    attr(BGData, "origFile") <- list(path = fileIn, dataType = typeof(dataType))
    attr(BGData, "dateCreated") <- date()
    save(BGData, file = paste0(folderOut, "/BGData.RData"))

    return(BGData)
}


#' @rdname readRAW
#' @export
readRAW_matrix <- function(fileIn, header = TRUE, dataType = integer(), n = NULL, p = NULL, sep = "", na.strings = "NA", nColSkip = 6L, idCol = c(1L, 2L), verbose = FALSE) {

    dims <- pedDims(fileIn = fileIn, header = header, n = n, p = p, sep = sep, nColSkip = nColSkip)

    dataType <- normalizeType(dataType)

    # Prepare geno
    geno <- matrix(nrow = dims$n, ncol = dims$p)

    # Prepare pheno
    pheno <- as.data.frame(matrix(nrow = dims$n, ncol = nColSkip), stringsAsFactors = FALSE)

    # Construct BGData object
    BGData <- new("BGData", geno = geno, pheno = pheno)

    # Parse .raw file
    BGData <- parseRAW(BGData = BGData, fileIn = fileIn, header = header, dataType = dataType, nColSkip = nColSkip, idCol = idCol, sep = sep, na.strings = na.strings, verbose = verbose)

    return(BGData)
}


#' @rdname readRAW
#' @export
readRAW_big.matrix <- function(fileIn, header = TRUE, dataType = integer(), n = NULL, p = NULL, sep = "", na.strings = "NA", nColSkip = 6L, idCol = c(1L, 2L), folderOut = paste0("BGData_", sub("\\.[[:alnum:]]+$", "", basename(fileIn))), outputType = "char", verbose = FALSE) {

    if (file.exists(folderOut)) {
        stop(paste("Output folder", folderOut, "already exists. Please move it or pick a different one."))
    }

    dataType <- normalizeType(dataType)
    if (!typeof(dataType) %in% c("integer", "double")) {
        stop("dataType must be either integer() or double()")
    }

    dims <- pedDims(fileIn = fileIn, header = header, n = n, p = p, sep = sep, nColSkip = nColSkip)

    options(bigmemory.typecast.warning = FALSE)
    options(bigmemory.allow.dimnames = TRUE)

    # Create output directory
    dir.create(folderOut)

    # Prepare geno
    geno <- bigmemory::filebacked.big.matrix(nrow = dims$n, ncol = dims$p, type = outputType, backingpath = folderOut, backingfile = "BGData.bin", descriptorfile = "BGData.desc")

    # Prepare pheno
    pheno <- as.data.frame(matrix(nrow = dims$n, ncol = nColSkip), stringsAsFactors = FALSE)

    # Construct BGData object
    BGData <- new("BGData", geno = geno, pheno = pheno)

    # Parse .raw file
    BGData <- parseRAW(BGData = BGData, fileIn = fileIn, header = header, dataType = dataType, nColSkip = nColSkip, idCol = idCol, sep = sep, na.strings = na.strings, verbose = verbose)

    # Save BGData object
    attr(BGData, "origFile") <- list(path = fileIn, dataType = typeof(dataType))
    attr(BGData, "dateCreated") <- date()
    save(BGData, file = paste0(folderOut, "/BGData.RData"))

    return(BGData)
}


loadFamFile <- function(path) {
    if (!file.exists(path)) {
        stop(path, " not found")
    }
    message("Extracting phenotypes from .fam file...")
    if (requireNamespace("data.table", quietly = TRUE)) {
        pheno <- data.table::fread(path, col.names = c(
            "FID",
            "IID",
            "PAT",
            "MAT",
            "SEX",
            "PHENOTYPE"
        ), data.table = FALSE, showProgress = FALSE)
    } else {
        pheno <- utils::read.table(path, col.names = c(
            "FID",
            "IID",
            "PAT",
            "MAT",
            "SEX",
            "PHENOTYPE"
        ), stringsAsFactors = FALSE)
    }
    return(pheno)
}


generatePheno <- function(x) {
    # Extract path to BED file
    bedPath <- attr(x, "path")
    # Try to load .fam file, generate pheno otherwise
    ex <- try({
        pheno <- loadFamFile(sub(".bed", ".fam", bedPath))
    }, silent = TRUE)
    if (class(ex) == "try-error") {
        splits <- strsplit(rownames(x), "_")
        pheno <- data.frame(FID = sapply(splits, "[", 1L), IID = sapply(splits, "[", 2L), stringsAsFactors = FALSE)
    }
    rownames(pheno) <- paste0(pheno$FID, "_", pheno$IID)
    return(pheno)
}


loadBimFile <- function(path) {
    if (!file.exists(path)) {
        stop(path, " not found")
    }
    message("Extracting map from .bim file...")
    if (requireNamespace("data.table", quietly = TRUE)) {
        map <- data.table::fread(path, col.names = c(
            "chromosome",
            "snp_id",
            "genetic_distance",
            "base_pair_position",
            "allele_1",
            "allele_2"
        ), data.table = FALSE, showProgress = FALSE)
    } else {
        map <- utils::read.table(path, col.names = c(
            "chromosome",
            "snp_id",
            "genetic_distance",
            "base_pair_position",
            "allele_1",
            "allele_2"
        ), stringsAsFactors = FALSE)
    }
    return(map)
}


generateMap <- function(x) {
    # Extract path to BED file
    bedPath <- attr(x, "path")
    # Try to load .fam file, generate pheno otherwise
    ex <- try({
        map <- loadBimFile(sub(".bed", ".bim", bedPath))
    }, silent = TRUE)
    if (class(ex) == "try-error") {
        splits <- strsplit(colnames(x), "_")
        map <- data.frame(
            snp_id = sapply(splits, function(x) {
                paste0(x[seq_len(length(x) - 1L)], collapse = "_")
            }),
            allele_1 = sapply(splits, function(x) {
                x[length(x)]
            }),
            stringsAsFactors = FALSE
        )
    }
    rownames(map) <- paste0(map$snp_id, "_", map$allele_1)
    return(map)
}


loadAlternatePhenotypeFile <- function(path, ...) {
    if (!file.exists(path)) {
        stop("Alternate phenotype file does not exist.")
    } else {
        message("Merging alternate phenotype file...")
        if (requireNamespace("data.table", quietly = TRUE)) {
            alternatePhenotypes <- data.table::fread(path, data.table = FALSE, showProgress = FALSE, ...)
        } else {
            # Check if the file has a header, i.e. if the first row starts with
            # an FID and an IID entry
            hasHeader = FALSE
            if (grepl("FID\\s+IID", readLines(path, n = 1L))) {
                hasHeader = TRUE
            }
            alternatePhenotypes <- utils::read.table(path, header = hasHeader, stringsAsFactors = FALSE, ...)
        }
    }
    return(alternatePhenotypes)
}


#' Merge Two Data Frames Keeping the Order of the First
#'
#' This is a simplified version of [base::merge()] useful for merging
#' additional data into a [BGData-class] object while keeping the order of the
#' data in the [BGData-class] object.
#'
#' @param x Data frame
#' @param y Data frame
#' @param by Specifications of the columns used for merging. Defaults to the
#' first two columns of the data frame, which traditionally has the family ID
#' and the individual ID.
#' @return Merged data frame
#' @export
orderedMerge <- function(x, y, by = c(1L, 2L)) {
    # Add artificial sort column to preserve order after merging
    # (merge's `sort = FALSE` order is unspecified)
    x$.sortColumn <- seq_len(nrow(x))
    # Merge phenotypes and alternate phenotypes
    merged <- merge(x, y, by = by, all.x = TRUE)
    # Reorder phenotypes to match original order and delete artificial
    # column
    merged <- merged[order(merged$.sortColumn), ]
    merged <- merged[, names(merged) != ".sortColumn"]
    # Restore rownames (assuming order is retained and no rows disappear...)
    rownames(merged) <- rownames(x)
    return(merged)
}


#' Convert Other Objects to BGData Objects.
#'
#' Converts other objects to [BGData-class] objects by loading supplementary
#' phenotypes and map files referenced by the object to be used for the
#' `@@pheno` and `@@map` slot, respectively. Currently supported are
#' [BEDMatrix::BEDMatrix-class] objects, plain or nested in
#' [LinkedMatrix::ColumnLinkedMatrix-class] objects.
#'
#' The .ped and .raw formats only allows for a single phenotype. If more
#' phenotypes are required it is possible to store them in an [alternate
#' phenotype file](https://www.cog-genomics.org/plink2/input#pheno). The path
#' to such a file can be provided with `alternatePhenotypeFile` and will be
#' merged with the data in the `@@pheno` slot. The first and second columns of
#' that file must contain family and within-family IDs, respectively.
#'
#' For [BEDMatrix::BEDMatrix-class] objects: If a .fam file (which corresponds
#' to the first six columns of a .ped or .raw file) of the same name and in the
#' same directory as the BED file exists, the `@@pheno` slot will be populated
#' with the data stored in that file. Otherwise a stub that only contains an
#' `IID` column populated with the rownames of `@@geno` will be generated. The
#' same will happen for a .bim file for the `@@map` slot.
#'
#' For [LinkedMatrix::ColumnLinkedMatrix-class] objects: See the case for
#' [BEDMatrix::BEDMatrix-class] objects, but only the .fam file of the first
#' node of the [LinkedMatrix::LinkedMatrix-class] will be read and used for the
#' `@@pheno` slot, and the .bim files of all nodes will be combined and used
#' for the `@@map` slot.
#'
#' @param x An object. Currently supported are [BEDMatrix::BEDMatrix-class]
#' objects, plain or nested in [LinkedMatrix::ColumnLinkedMatrix-class]
#' objects.
#' @param alternatePhenotypeFile Path to an [alternate phenotype
#' file](https://www.cog-genomics.org/plink2/input#pheno).
#' @param ... Additional arguments to the [utils::read.table()] or
#' [data.table::fread()] call (if data.table package is installed) call to
#' parse the alternate pheno file.
#' @return A [BGData-class] object.
#' @seealso [readRAW()] to convert text files to [BGData-class] objects.
#' @example man/examples/as.BGData.R
#' @export
as.BGData <- function(x, alternatePhenotypeFile = NULL, ...) {
    UseMethod("as.BGData")
}


#' @rdname as.BGData
#' @export
as.BGData.BEDMatrix <- function(x, alternatePhenotypeFile = NULL, ...) {
    # Read in pheno file
    fam <- generatePheno(x)
    # Read in map file
    map <- generateMap(x)
    # Load and merge alternate phenotype file
    if (!is.null(alternatePhenotypeFile)) {
        alternatePhenotypes <- loadAlternatePhenotypeFile(alternatePhenotypeFile, ...)
        fam <- orderedMerge(fam, alternatePhenotypes)
    }
    BGData(geno = x, pheno = fam, map = map)
}


#' @rdname as.BGData
#' @export
as.BGData.ColumnLinkedMatrix <- function(x, alternatePhenotypeFile = NULL, ...) {
    n <- LinkedMatrix::nNodes(x)
    # For now, all elements have to be of type BEDMatrix
    if (!all(sapply(x, function(node) class(node)) == "BEDMatrix")) {
        stop("Only BEDMatrix instances are supported as elements of the LinkedMatrix right now.")
    }
    # Read in the fam file of the first node
    message("Extracting phenotypes from .fam file, assuming that the .fam file of the first BEDMatrix instance is representative of all the other nodes...")
    fam <- suppressMessages(generatePheno(x[[1L]]))
    # Read in map files
    message("Extracting map from .bim files...")
    map <- do.call(base::rbind, lapply(x, function(node) {
        suppressMessages(generateMap(node))
    }))
    # Load and merge alternate phenotype file
    if (!is.null(alternatePhenotypeFile)) {
        alternatePhenotypes <- loadAlternatePhenotypeFile(alternatePhenotypeFile, ...)
        fam <- orderedMerge(fam, alternatePhenotypes)
    }
    BGData(geno = x, pheno = fam, map = map)
}


#' @rdname as.BGData
#' @export
as.BGData.RowLinkedMatrix <- function(x, alternatePhenotypeFile = NULL, ...) {
    n <- LinkedMatrix::nNodes(x)
    # For now, all elements have to be of type BEDMatrix
    if (!all(sapply(x, function(node) class(node)) == "BEDMatrix")) {
        stop("Only BEDMatrix instances are supported as elements of the LinkedMatrix right now.")
    }
    # Read in the fam files
    message("Extracting phenotypes from .fam files...")
    fam <- do.call(base::rbind, lapply(x, function(node) {
        suppressMessages(generatePheno(node))
    }))
    # Read in the map file of the first node
    message("Extracting map from .bim file, assuming that the .bim file of the first BEDMatrix instance is representative of all the other nodes...")
    map <- suppressMessages(generateMap(x[[1L]]))
    # Load and merge alternate phenotype file
    if (!is.null(alternatePhenotypeFile)) {
        alternatePhenotypes <- loadAlternatePhenotypeFile(alternatePhenotypeFile, ...)
        fam <- orderedMerge(fam, alternatePhenotypes)
    }
    BGData(geno = x, pheno = fam, map = map)
}


#' Loads BGData (and Other) Objects from .RData Files.
#'
#' This function is similar to [base::load()], but also initializes the
#' different types of objects that the `@@geno` slot of a [BGData-class] object
#' can take. Currently supported are `ff_matrix`,
#' [bigmemory::big.matrix-class], and [BEDMatrix::BEDMatrix-class] objects. If
#' the object is of type [LinkedMatrix::LinkedMatrix-class], all nodes will be
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
        # Initialize genotypes of BGData objects
        if (class(object) == "BGData") {
            object@geno <- initializeGeno(object@geno, path = dirname(file))
        }
        # Assign object to envir
        assign(name, object, envir = envir)
    }
    message("Loaded objects: ", paste0(names, collapse = ", "))
}


initializeGeno <- function(x, ...) {
    UseMethod("initializeGeno")
}


initializeGeno.LinkedMatrix <- function(x, path, ...) {
    for (i in seq_len(LinkedMatrix::nNodes(x))) {
        x[[i]] <- initializeGeno(x[[i]], path = path)
    }
    return(x)
}


# Absolute paths to ff files are not stored, so the ff objects have to be
# loaded from the same directory as the RData file.
initializeGeno.ff_matrix <- function(x, path, ...) {
    # Store current working directory and set working directory to path
    cwd <- getwd()
    setwd(path)
    # Open ff object
    ff::open.ff(x)
    # Restore the working directory
    setwd(cwd)
    return(x)
}


initializeGeno.big.matrix <- function(x, path, ...) {
    return(bigmemory::attach.big.matrix(paste0(path, "/BGData.desc")))
}


initializeGeno.BEDMatrix <- function(x, ...) {
    dnames <- attr(x, "dnames")
    dims <- attr(x, "dims")
    path <- attr(x, "path")
    x <- BEDMatrix::BEDMatrix(path = path, n = dims[1L], p = dims[2L])
    dimnames(x) <- dnames
    return(x)
}


initializeGeno.default <- function(x, ...) {
    return(x)
}


ffNodeInitializer <- function(nodeIndex, nrow, ncol, vmode, folderOut, ...) {
    filename <- paste0("geno_", nodeIndex, ".bin")
    node <- ff::ff(dim = c(nrow, ncol), vmode = vmode, filename = paste0(folderOut, "/", filename), ...)
    # Change ff path to a relative one
    bit::physical(node)$filename <- filename
    return(node)
}
