# Convert ff_matrix into an S4 class
setOldClass("ff_matrix")


# Convert BEDMatrix into an S4 class
setOldClass("BEDMatrix")


#' An Abstract S4 Class Union of Matrix-Like Types.
#'
#' [geno-class] is a class union of several matrix-like types, many of them
#' suitable for very large datasets. Currently supported are
#' [LinkedMatrix::LinkedMatrix-class], [BEDMatrix::BEDMatrix],
#' [bigmemory::big.matrix-class], `ff_matrix`, and `matrix`.
#'
#' @seealso The `@@geno` slot of [BGData-class] that accepts [geno-class]
#' objects.
setClassUnion("geno", c("LinkedMatrix", "BEDMatrix", "big.matrix", "ff_matrix", "matrix"))


#' An S4 Class to Represent Phenotype and Genotype Data.
#'
#' This class is inspired by the phenotype/genotype file format PED and its
#' binary companion (also known as BED) of
#' [PLINK](http://pngu.mgh.harvard.edu/~purcell/plink/). It is used by several
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
#' * from a raw PED file (or a PED-like file) using [readPED()],
#' [readPED.matrix()], or [readPED.big.matrix()].
#'
#' A PED file can be recoded to a raw PED in
#' [PLINK](http://pngu.mgh.harvard.edu/~purcell/plink/) using `plink --file
#' myfile --recodeA`, or converted to a BED file using `plink --file myfile
#' --make-bed`. Conversely, a BED file can be transformed back to a PED file
#' using `plink --bfile myfile --recode` or to a raw PED file using `plink
#' --bfile myfile --recodeA` without losing information.
#'
#' @slot geno A [geno-class] object that contains genotypes. [geno-class] is a
#' class union of several matrix-like types, many of them suitable for very
#' large datasets. Currently supported are [LinkedMatrix::LinkedMatrix-class],
#' [BEDMatrix::BEDMatrix], [bigmemory::big.matrix-class], `ff_matrix`, and
#' `matrix`.
#' @slot pheno A `data.frame` that contains phenotypes.
#' @slot map A `data.frame` that contains a genetic map.
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
#' [BEDMatrix::BEDMatrix], [bigmemory::big.matrix-class], `ff_matrix`, and
#' `matrix`.
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
    pedFile <- file(fileIn, open = "r")

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


#' Creates a BGData Object From a Raw PED File or a PED-Like File.
#'
#' Creates a [BGData-class] object from a raw PED file (generated with
#' `--recodeA` in [PLINK](http://pngu.mgh.harvard.edu/~purcell/plink/)). Other
#' text-based file formats are supported as well by tweaking some of the
#' parameters as long as the records of individuals are in rows, and
#' phenotypes, covariates and markers are in columns.
#'
#' The data included in the first couple of columns (up to `nColSkip`) is used
#' to populate the `@@pheno` slot of a [BGData-class] object, and the remaining
#' columns are used to fill the `@@geno` slot. If the first row contains a
#' header (`header = TRUE`), data in this row is used to determine the column
#' names for `@@pheno` and `@@geno`.
#'
#' `@@geno` can take several forms, depending on the function that is called
#' (`readPED`, `readPED.matrix`, or `readPED.big.matrix`). The following
#' sections illustrate each function in detail.
#'
#' @section readPED:
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
#' @section readPED.matrix:
#' Genotypes are stored in a regular `matrix` object. Therefore, this function
#' will only work if the raw PED file is small enough to fit into memory.
#'
#' @section readPED.big.matrix:
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
#' supported for [readPED()] and [readPED.big.matrix()]: use the `--recodeA`
#' option of PLINK to convert the PED file into a raw file. Defaults to
#' `integer()`.
#' @param n The number of individuals. Auto-detect if `NULL`. Defaults to
#' `NULL`.
#' @param p The number of markers. Auto-detect if `NULL`. Defaults to `NULL`.
#' @param sep The field separator character. Values on each line of the file
#' are separated by this character. If `sep = ""` (the default for [readPED()]
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
#' @export
readPED <- function(fileIn, header = TRUE, dataType = integer(), n = NULL, p = NULL, sep = "", na.strings = "NA", nColSkip = 6, idCol = c(1, 2), nNodes = NULL, linked.by = "rows", folderOut = paste0("BGData_", sub("\\.[[:alnum:]]+$", "", basename(fileIn))), outputType = "byte", dimorder = if (linked.by == "rows") 2:1 else 1:2, verbose = FALSE) {

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


#' @rdname readPED
#' @export
readPED.matrix <- function(fileIn, header = TRUE, dataType = integer(), n = NULL, p = NULL, sep = "", na.strings = "NA", nColSkip = 6, idCol = c(1, 2), verbose = FALSE) {

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


#' @rdname readPED
#' @export
readPED.big.matrix <- function(fileIn, header = TRUE, dataType = integer(), n = NULL, p = NULL, sep = "", na.strings = "NA", nColSkip = 6, idCol = c(1, 2), folderOut = paste0("BGData_", sub("\\.[[:alnum:]]+$", "", basename(fileIn))), outputType = "char", verbose = FALSE) {

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

    # Parse PED file
    BGData <- parsePED(BGData = BGData, fileIn = fileIn, header = header, dataType = dataType, nColSkip = nColSkip, idCol = idCol, sep = sep, na.strings = na.strings, verbose = verbose)

    # Save BGData object
    attr(BGData, "origFile") <- list(path = fileIn, dataType = typeof(dataType))
    attr(BGData, "dateCreated") <- date()
    save(BGData, file = paste0(folderOut, "/BGData.RData"))

    return(BGData)
}


loadFamFile <- function(path) {
    if (file.exists(path)) {
        message("Extracting phenotypes from FAM file...")
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
    } else {
        splits <- strsplit(rownames(x), "_")
        pheno <- data.frame(FID = sapply(splits, "[", 1), IID = sapply(splits, "[", 2), stringsAsFactors = FALSE)
    }
    return(pheno)
}


loadBimFile <- function(path) {
    if (file.exists(path)) {
        message("Extracting map from BIM file...")
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
    } else {
        splits <- strsplit(colnames(x), "_")
        map <- data.frame(path, snp_id = sapply(splits, function(x) {
            paste0(x[seq_len(length(x) - 1)], collapse = "_")
        }), allele_1 = sapply(splits, function(x) {
            x[length(x)]
        }), stringsAsFactors = FALSE)
    }
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
            if (grepl("FID\\s+IID", readLines(path, n = 1))) {
                hasHeader = TRUE
            }
            alternatePhenotypes <- utils::read.table(path, header = hasHeader, stringsAsFactors = FALSE, ...)
        }
    }
    return(alternatePhenotypes)
}


mergeAlternatePhenotypes <- function(pheno, alternatePhenotypes) {
    # Add artificial sort column to preserve order after merging
    # (merge's `sort = FALSE` order is unspecified)
    pheno$.sortColumn <- seq_len(nrow(pheno))
    # Merge phenotypes and alternate phenotypes
    pheno <- merge(pheno, alternatePhenotypes, by = c(1L, 2L), all.x = TRUE)
    # Reorder phenotypes to match original order and delete artificial
    # column
    pheno <- pheno[order(pheno$.sortColumn), ]
    pheno <- pheno[, names(pheno) != ".sortColumn"]
    return(pheno)
}


#' Convert Other Objects to BGData Objects.
#'
#' Converts other objects to [BGData-class] objects by loading supplementary
#' phenotypes and map files referenced by the object to be used for the
#' `@@pheno` and `@@map` slot, respectively. Currently supported are
#' [BEDMatrix::BEDMatrix] objects, plain or nested in
#' [LinkedMatrix::ColumnLinkedMatrix-class] objects.
#'
#' The PED format only allows for a single phenotype. If more phenotypes are
#' required it is possible to store them in an [alternate phenotype
#' file](http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#pheno). The path
#' to such a file can be provided with `alternatePhenotypeFile` and will be
#' merged with the data in the `@@pheno` slot.
#'
#' For [BEDMatrix::BEDMatrix] objects: If a FAM file (which corresponds to the
#' first six columns of a PED file) of the same name and in the same directory
#' as the BED file exists, the `@@pheno` slot will be populated with the data
#' stored in that file.  Otherwise a stub that only contains an `IID` column
#' populated with the rownames of `@@geno` will be generated. The same will
#' happen for a BIM file for the `@@map` slot.
#'
#' For [LinkedMatrix::ColumnLinkedMatrix-class] objects: See the case for
#' [BEDMatrix::BEDMatrix] objects, but only the FAM file of the first node of
#' the [LinkedMatrix::LinkedMatrix-class] will be read and used for the
#' `@@pheno` slot, and the BIM files of all nodes will be combined and used for
#' the `@@map` slot.
#'
#' @param x An object. Currently supported are [BEDMatrix::BEDMatrix] objects,
#' plain or nested in [LinkedMatrix::ColumnLinkedMatrix-class] objects.
#' @param alternatePhenotypeFile Path to an [alternate phenotype
#' file](http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#pheno).
#' @param ... Additional arguments to the [utils::read.table()] or
#' [data.table::fread()] call (if data.table package is installed) call to
#' parse the alternate pheno file.
#' @return A [BGData-class] object.
#' @seealso [readPED()] to convert text files to [BGData-class] objects.
#' @export
as.BGData <- function(x, alternatePhenotypeFile = NULL, ...) {
    UseMethod("as.BGData")
}


#' @rdname as.BGData
#' @export
as.BGData.BEDMatrix <- function(x, alternatePhenotypeFile = NULL, ...) {
    # Extract path to BED file
    bedPath <- attr(x, "path")
    # Read in pheno file
    fam <- loadFamFile(sub(".bed", ".fam", bedPath))
    # Read in map file
    map <- loadBimFile(sub(".bed", ".bim", bedPath))
    # Load and merge alternate phenotype file
    if (!is.null(alternatePhenotypeFile)) {
        alternatePhenotypes <- loadAlternatePhenotypeFile(alternatePhenotypeFile, ...)
        fam <- mergeAlternatePhenotypes(fam, alternatePhenotypes)
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
    message("Extracting phenotypes from FAM file, assuming that the FAM file of the first BEDMatrix instance is representative of all the other nodes...")
    fam <- suppressMessages(loadFamFile(sub(".bed", ".fam", attr(x[[1]], "path"))))
    # Read in map files
    message("Extracting map from BIM files...")
    map <- do.call("rbind", lapply(x, function(node) {
        suppressMessages(loadBimFile(sub(".bed", ".bim", attr(node, "path"))))
    }))
    # Load and merge alternate phenotype file
    if (!is.null(alternatePhenotypeFile)) {
        alternatePhenotypes <- loadAlternatePhenotypeFile(alternatePhenotypeFile, ...)
        fam <- mergeAlternatePhenotypes(fam, alternatePhenotypes)
    }
    BGData(geno = x, pheno = fam, map = map)
}


#' Loads BGData (and Other) Objects from .RData Files.
#'
#' This function is similar to [base::load()], but also initializes the
#' different types of objects that the `@@geno` slot of a [BGData-class] object
#' can take. Currently supported are `ff_matrix`,
#' [bigmemory::big.matrix-class], and [BEDMatrix::BEDMatrix] objects. If the
#' object is of type [LinkedMatrix::LinkedMatrix-class], all nodes will be
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
