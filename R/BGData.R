# Convert ff_matrix into an S4 class
setOldClass("ff_matrix")

setClassUnion("geno", c("LinkedMatrix", "BEDMatrix", "big.matrix", "ff_matrix", "matrix"))

setClass("BGData", slots = c(geno = "geno", pheno = "data.frame", map = "data.frame"))

BGData <- function(geno, pheno = NULL, map = NULL) {
    if (!is(geno, "geno")) {
        stop("Only LinkedMatrix, BEDMatrix, big.matrix, ff_matrix, or regular matrix objects are allowed for geno.")
    }
    if (is.null(pheno)) {
        if (is.null(rownames(geno))) {
            sampleIDs <- paste0("sample_", seq_len(nrow(geno)))
        } else {
            sampleIDs <- rownames(geno)
        }
        pheno <- data.frame(sample_id = sampleIDs, row.names = sampleIDs, stringsAsFactors = FALSE)
    }
    if (is.null(map)) {
        if (is.null(colnames(geno))) {
            variantIDs <- paste0("variant_", seq_len(ncol(geno)))
        } else {
            variantIDs <- colnames(geno)
        }
        map <- data.frame(variant_id = variantIDs, row.names = variantIDs, stringsAsFactors = FALSE)
    }
    obj <- new("BGData", geno = geno, pheno = pheno, map = map)
    return(obj)
}

setValidity("BGData", function(object) {
    if (nrow(slot(object, "geno")) != nrow(slot(object, "pheno"))) {
        return("Number of rows of geno and number of rows of pheno do not match.")
    }
    # Do not assume that geno has row names, but if it does, it should match
    # the row names of pheno
    if (!is.null(rownames(slot(object, "geno"))) && any(rownames(slot(object, "geno")) != rownames(slot(object, "pheno")))) {
        warning("Row names of geno and row names of pheno do not match.")
    }
    if (ncol(slot(object, "geno")) != nrow(slot(object, "map"))) {
        return("Number of columns of geno and number of rows of map do not match.")
    }
    # Do not assume that geno has column names, but if it does, it should match
    # the row names of map
    if (!is.null(colnames(slot(object, "geno"))) && any(colnames(slot(object, "geno")) != rownames(slot(object, "map")))) {
        warning("Column names of geno and row names of map do not match.")
    }
    return(TRUE)
})

setGeneric("geno", function(x) standardGeneric("geno"))
setMethod("geno", "BGData", function(x) slot(x, "geno"))

setGeneric("geno<-", function(x, value) standardGeneric("geno<-"))
setMethod("geno<-", "BGData", function(x, value) {
    slot(x, "geno") <- value
    validObject(x)
    x
})

setGeneric("pheno", function(x) standardGeneric("pheno"))
setMethod("pheno", "BGData", function(x) slot(x, "pheno"))

setGeneric("pheno<-", function(x, value) standardGeneric("pheno<-"))
setMethod("pheno<-", "BGData", function(x, value) {
    slot(x, "pheno") <- value
    validObject(x)
    x
})

setGeneric("map", function(x) standardGeneric("map"))
setMethod("map", "BGData", function(x) slot(x, "map"))

setGeneric("map<-", function(x, value) standardGeneric("map<-"))
setMethod("map<-", "BGData", function(x, value) {
    slot(x, "map") <- value
    validObject(x)
    x
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

    p <- ncol(geno(BGData))
    pedFile <- file(fileIn, open = "r")

    # Update colnames
    if (header) {
        headerLine <- scan(pedFile, nlines = 1L, what = character(), sep = sep, quiet = TRUE)
        # Suppress warnings here to not get in trouble with validity method
        suppressWarnings(colnames(pheno(BGData)) <- headerLine[seq_len(nColSkip)])
        suppressWarnings(colnames(geno(BGData)) <- headerLine[-(seq_len(nColSkip))])
        suppressWarnings(rownames(map(BGData)) <- colnames(geno(BGData)))
    }

    # Parse file
    for (i in seq_len(nrow(geno(BGData)))) {
        xSkip <- scan(pedFile, n = nColSkip, what = character(), sep = sep, quiet = TRUE)
        x <- scan(pedFile, n = p, what = dataType, sep = sep, na.strings = na.strings, quiet = TRUE)
        pheno(BGData)[i, ] <- xSkip
        geno(BGData)[i, ] <- x
        if (verbose) {
            message("Subject ", i, " / ", nrow(geno(BGData)))
        }
    }
    close(pedFile)

    # Update rownames
    IDs <- apply(pheno(BGData)[, idCol, drop = FALSE], 1L, paste, collapse = "_")
    rownames(pheno(BGData)) <- IDs
    rownames(geno(BGData)) <- IDs

    # Convert types in pheno
    pheno(BGData)[] <- lapply(pheno(BGData), utils::type.convert, as.is = TRUE)

    return(BGData)
}

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
    BGData <- BGData(geno = geno, pheno = pheno)

    # Parse .raw file
    BGData <- parseRAW(BGData = BGData, fileIn = fileIn, header = header, dataType = dataType, nColSkip = nColSkip, idCol = idCol, sep = sep, na.strings = na.strings, verbose = verbose)

    # Save BGData object
    attr(BGData, "origFile") <- list(path = fileIn, dataType = typeof(dataType))
    attr(BGData, "dateCreated") <- date()
    save(BGData, file = paste0(folderOut, "/BGData.RData"))

    return(BGData)
}

readRAW_matrix <- function(fileIn, header = TRUE, dataType = integer(), n = NULL, p = NULL, sep = "", na.strings = "NA", nColSkip = 6L, idCol = c(1L, 2L), verbose = FALSE) {

    dims <- pedDims(fileIn = fileIn, header = header, n = n, p = p, sep = sep, nColSkip = nColSkip)

    dataType <- normalizeType(dataType)

    # Prepare geno
    geno <- matrix(nrow = dims$n, ncol = dims$p)

    # Prepare pheno
    pheno <- as.data.frame(matrix(nrow = dims$n, ncol = nColSkip), stringsAsFactors = FALSE)

    # Construct BGData object
    BGData <- BGData(geno = geno, pheno = pheno)

    # Parse .raw file
    BGData <- parseRAW(BGData = BGData, fileIn = fileIn, header = header, dataType = dataType, nColSkip = nColSkip, idCol = idCol, sep = sep, na.strings = na.strings, verbose = verbose)

    return(BGData)
}

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
    BGData <- BGData(geno = geno, pheno = pheno)

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
    # It was considered to read the PHENOTYPE column as double, but the PLINK
    # documentation mentions non-numeric case/control values.
    if (requireNamespace("data.table", quietly = TRUE)) {
        pheno <- data.table::fread(path, col.names = c(
            "FID",
            "IID",
            "PAT",
            "MAT",
            "SEX",
            "PHENOTYPE"
        ), colClasses = "character", data.table = FALSE, showProgress = FALSE)
    } else {
        pheno <- utils::read.table(path, col.names = c(
            "FID",
            "IID",
            "PAT",
            "MAT",
            "SEX",
            "PHENOTYPE"
        ), colClasses = "character", stringsAsFactors = FALSE)
    }
    return(pheno)
}

generatePheno <- function(x) {
    # Extract path to BED file
    bedPath <- attr(x, "path")
    # Try to load .fam file, generate pheno otherwise
    ex <- try({
        pheno <- loadFamFile(sub("\\.bed", "\\.fam", bedPath))
    }, silent = TRUE)
    if (inherits(ex, "try-error")) {
        # x may not have rownames (e.g., when a BEDMatrix is created using the
        # n parameter)
        if (is.null(rownames(x))) {
            pheno <- data.frame(FID = "0", IID = as.character(1:nrow(x)), stringsAsFactors = FALSE)
        } else {
            # Make no assumptions about the structure of the rownames of x
            # here, i.e., do not try to extract FID and IID.
            pheno <- data.frame(FID = "0", IID = rownames(x), stringsAsFactors = FALSE)
        }
    }
    # Preserve rownames of x (if not NULL)
    rownames(pheno) <- rownames(x)
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
        ), colClasses = c(
            "character",
            "character",
            "integer",
            "integer",
            "character",
            "character"
        ), data.table = FALSE, showProgress = FALSE)
    } else {
        map <- utils::read.table(path, col.names = c(
            "chromosome",
            "snp_id",
            "genetic_distance",
            "base_pair_position",
            "allele_1",
            "allele_2"
        ), colClasses = c(
            "character",
            "character",
            "integer",
            "integer",
            "character",
            "character"
        ), stringsAsFactors = FALSE)
    }
    return(map)
}

generateMap <- function(x) {
    # Extract path to BED file
    bedPath <- attr(x, "path")
    # Try to load .bim file, generate map otherwise
    ex <- try({
        map <- loadBimFile(sub("\\.bed", "\\.bim", bedPath))
    }, silent = TRUE)
    if (inherits(ex, "try-error")) {
        # x may not have colnames (e.g., when a BEDMatrix is created using the
        # p parameter)
        if (is.null(colnames(x))) {
            map <- data.frame(snp_id = as.character(1:ncol(x)), stringsAsFactors = FALSE)
        } else {
            # Make no assumptions about the structure of the colnames of x
            # here, i.e., do not try to extract minor allele.
            map <- data.frame(snp_id = colnames(x), stringsAsFactors = FALSE)
        }
    }
    # Preserve colnames of x (if not NULL)
    rownames(map) <- colnames(x)
    return(map)
}

loadAlternatePhenotypeFile <- function(path, ...) {
    if (!file.exists(path)) {
        stop("Alternate phenotype file does not exist.")
    } else {
        message("Merging alternate phenotype file...")
        if (requireNamespace("data.table", quietly = TRUE)) {
            alternatePhenotypes <- data.table::fread(path, colClasses = list(
                character = 1:2
            ), data.table = FALSE, showProgress = FALSE, ...)
        } else {
            # Check if the file has a header, i.e. if the first row starts with
            # an FID and an IID entry
            hasHeader = FALSE
            if (grepl("FID\\s+IID", readLines(path, n = 1L))) {
                hasHeader = TRUE
            }
            alternatePhenotypes <- utils::read.table(path, header = hasHeader, stringsAsFactors = FALSE, ...)
            alternatePhenotypes[[1]] <- as.character(alternatePhenotypes[[1]]) # FID
            alternatePhenotypes[[2]] <- as.character(alternatePhenotypes[[2]]) # IID
        }
    }
    return(alternatePhenotypes)
}

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

as.BGData <- function(x, alternatePhenotypeFile = NULL, ...) {
    UseMethod("as.BGData")
}

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

as.BGData.ColumnLinkedMatrix <- function(x, alternatePhenotypeFile = NULL, ...) {
    n <- LinkedMatrix::nNodes(x)
    # For now, all elements have to be of type BEDMatrix
    if (!all(vapply(x, inherits, TRUE, "BEDMatrix"))) {
        stop("Only BEDMatrix instances are supported as elements of the LinkedMatrix right now.")
    }
    # Read in the fam file of the first node
    message("Extracting phenotypes from .fam file, assuming that the .fam file of the first BEDMatrix instance is representative of all the other nodes...")
    fam <- suppressMessages(generatePheno(x[[1L]]))
    # Read in map files
    message("Extracting map from .bim files...")
    map <- do.call(rbind, lapply(x, function(node) {
        suppressMessages(generateMap(node))
    }))
    # Load and merge alternate phenotype file
    if (!is.null(alternatePhenotypeFile)) {
        alternatePhenotypes <- loadAlternatePhenotypeFile(alternatePhenotypeFile, ...)
        fam <- orderedMerge(fam, alternatePhenotypes)
    }
    BGData(geno = x, pheno = fam, map = map)
}

as.BGData.RowLinkedMatrix <- function(x, alternatePhenotypeFile = NULL, ...) {
    n <- LinkedMatrix::nNodes(x)
    # For now, all elements have to be of type BEDMatrix
    if (!all(vapply(x, inherits, TRUE, "BEDMatrix"))) {
        stop("Only BEDMatrix instances are supported as elements of the LinkedMatrix right now.")
    }
    # Read in the fam files
    message("Extracting phenotypes from .fam files...")
    fam <- do.call(rbind, lapply(x, function(node) {
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

load.BGData <- function(file, envir = parent.frame()) {
    # Load data into new environment
    loadingEnv <- new.env()
    load(file = file, envir = loadingEnv)
    names <- ls(envir = loadingEnv)
    for (name in names) {
        object <- get(name, envir = loadingEnv)
        # Initialize genotypes of BGData objects
        if (inherits(object, "BGData")) {
            geno(object) <- initializeGeno(geno(object), path = dirname(file))
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
