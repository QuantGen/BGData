# Convert ff_matrix into an S4 class
setOldClass("ff_matrix")


# Convert BEDMatrix into an S4 class
setOldClass("BEDMatrix")


setClassUnion("geno", c("LinkedMatrix", "BEDMatrix", "big.matrix", "ff_matrix", "matrix"))


#' An S4 class to represent GWAS data.
#' 
#' @slot geno A \code{geno} object (\code{\linkS4class{LinkedMatrix}},
#'   \code{BEDMatrix}, \code{\linkS4class{big.matrix}} \code{ff_matrix}, or
#'   \code{\link{matrix}}) that contains genotypes.
#' @slot pheno A \code{\link{data.frame}} that contains phenotypes.
#' @slot map A \code{\link{data.frame}} that contains a genetic map.
#' @export BGData
#' @exportClass BGData
BGData <- setClass("BGData", slots = c(geno = "geno", pheno = "data.frame", map = "data.frame"))


#' Creates a new \code{BGData} instance.
#' 
#' @param .Object The \code{ColumnLinkedMatrix} instance to be initialized.
#' @param geno A \code{geno} object (\code{\linkS4class{LinkedMatrix}},
#'   \code{BEDMatrix}, \code{\linkS4class{big.matrix}}, \code{ff_matrix}, or
#' \code{\link{matrix}}) that contains genotypes.
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
        colnames(geno) <- paste0("mrk_", 1:ncol(geno))
    }
    if (is.null(rownames(geno))) {
        rownames(geno) <- paste0("id_", 1:nrow(geno))
    }
    if (missing(pheno)) {
        pheno <- data.frame(IID = rownames(geno))
        rownames(pheno) <- rownames(geno)
    }
    if (missing(map)) {
        map <- data.frame(mrk = colnames(geno))
        rownames(map) <- colnames(geno)
    }
    .Object@geno <- geno
    .Object@pheno <- pheno
    .Object@map <- map
    return(.Object)
})


pedDims <- function(fileIn, header, n, p, nColSkip = 6) {
    if (is.null(n)) {
        n <- getLineCount(fileIn, header)
    }
    if (header) {
        headerLine <- getFileHeader(fileIn)
        p <- length(headerLine) - nColSkip
    } else {
        if (is.null(p)) {
            p <- getColumnCount(fileIn) - nColSkip
        }
    }
    return(list(n = n, p = p))
}


parsePED <- function(BGData, fileIn, header, dataType, nColSkip = 6, idCol = c(1, 2), na.strings = "NA", verbose = FALSE, ...) {

    p <- ncol(BGData@geno)
    pedFile <- gzfile(fileIn, open = "r")

    # Update colnames
    if (header) {
        headerLine <- scan(pedFile, nlines = 1, what = character(), quiet = TRUE)
        colnames(BGData@pheno) <- headerLine[1:nColSkip]
        colnames(BGData@geno) <- headerLine[-(1:nColSkip)]
    }

    # Find replacement method to avoid repeated lookups
    lookupEx <- try({
        replace <- getMethod("[<-", class(BGData@geno))
    }, silent = TRUE)
    if (class(lookupEx) == "try-error") replace <- `[<-`

    # Parse file
    j <- 1:p
    for (i in 1:nrow(BGData@geno)) {
        time <- proc.time()
        xSkip <- scan(pedFile, n = nColSkip, what = character(), quiet = TRUE)
        x <- scan(pedFile, n = p, what = dataType, na.strings = na.strings, quiet = TRUE)
        BGData@pheno[i, ] <- xSkip
        BGData@geno <- replace(BGData@geno, i, j, ..., value = x)
        if (verbose) {
            cat("Subject", i, " ", round(proc.time()[3] - time[3], 3), "sec / subject.", "\n")
        }
    }
    close(pedFile)

    # Update rownames
    IDs <- apply(BGData@pheno[, idCol, drop = FALSE], 1, paste, collapse = "_")
    rownames(BGData@pheno) <- IDs
    rownames(BGData@geno) <- IDs

    # Convert types in pheno
    BGData@pheno[] <- lapply(BGData@pheno, type.convert, as.is = TRUE)

    return(BGData)
}


#' Creates a memory-mapped \code{\linkS4class{BGData}} object from a plaintext 
#' raw PED file (generated with \code{--recodeA} in PLINK) or a PED-like file.
#' 
#' \code{readPED} assumes that the plaintext file (\code{fileIn}) contains 
#' records of individuals in rows, and phenotypes, covariates and markers in 
#' columns. The columns included in the first couple of columns 
#' (\code{1:nColSkip}) are used to populate the \code{@@pheno} slot of a 
#' \code{\linkS4class{BGData}} object, and the remaining columns are used to 
#' fill the \code{@@geno} slot. If the first row contains a header 
#' (\code{header=TRUE}), data in this row is used to determine variables names 
#' for \code{@@pheno} and marker names for \code{@@map} and \code{@@geno}.
#' 
#' Genotypes are stored in a \code{\linkS4class{LinkedMatrix}} object, where 
#' each node is an \code{ff} instance. By default a column-linked 
#' (\code{\linkS4class{ColumnLinkedMatrix}}) is used for \code{@@geno}, but the 
#' user can modify this using the \code{linked.by} argument. The number of nodes
#' is either specified by the user using the \code{nNodes} argument or 
#' determined internally so that each \code{ff} object has a number of cells 
#' that is smaller than \code{.Machine$integer.max / 1.2} (the array limit of 
#' \code{ff}).
#' 
#' \code{readPED} creates a folder (see \code{folderOut}) that contains the 
#' binary flat files (named \code{geno_*.bin}) and an external representation of
#' the \code{\linkS4class{BGData}} object in \code{BGData.RData}. A 
#' \code{\linkS4class{BGData}} object can be reloaded using \code{load.BGData} 
#' (the regular \code{load} function will only work if the working directory is 
#' set to the path that contains the binary flat files).
#' 
#' @param fileIn The path to the plaintext file.
#' @param header If TRUE, the file contains a header.
#' @param dataType The coding of genotypes. Use \code{integer()} or 
#'   \code{double()} for numeric coding. Character coding is currently not 
#'   supported: use the \code{--recodeA} option of PLINK to convert the PED file
#'   into a raw file.
#' @param n The number of individuals.
#' @param p The number of markers.
#' @param na.strings The character string used in the plaintext file to denote 
#'   missing value.
#' @param nColSkip The number of columns to be skipped to reach the genotype 
#'   information in the file.
#' @param idCol The index of the ID column. If more than one index is given, 
#'   both columns will be concatenated with "_".
#' @param verbose If TRUE, progress updates will be posted.
#' @param nNodes The number of nodes to create.
#' @param linked.by If \code{columns} a column-linked matrix 
#'   (\code{\linkS4class{ColumnLinkedMatrix}}) is created, if \code{rows} a 
#'   row-linked matrix (\code{\linkS4class{RowLinkedMatrix}}).
#' @param folderOut The path to the folder where to save the binary files.
#' @param dimorder The physical layout of the underlying \code{ff} object of 
#'   each node.
#' @seealso \code{\linkS4class{BGData}}, \code{\linkS4class{LinkedMatrix}}, 
#'   \code{\linkS4class{ColumnLinkedMatrix}}, 
#'   \code{\linkS4class{RowLinkedMatrix}}, \code{\link[ff]{ff}}
#' @export
readPED <- function(fileIn, header, dataType, n = NULL, p = NULL, na.strings = "NA", nColSkip = 6, idCol = c(1, 2), verbose = FALSE, nNodes = NULL, linked.by = "rows", folderOut = paste("BGData_", sub("\\.[[:alnum:]]+$", "", basename(fileIn)), sep = ""), dimorder = if (linked.by == "rows") 2:1 else 1:2) {

    if (file.exists(folderOut)) {
        stop(paste("Output folder", folderOut, "already exists. Please move it or pick a different one."))
    }

    dataType <- normalizeType(dataType)
    if (!typeof(dataType) %in% c("integer", "double")) {
        stop("dataType must be either integer() or double()")
    }

    if (!linked.by %in% c("columns", "rows")) {
        stop("linked.by must be either columns or rows")
    }

    class <- ifelse(linked.by == "columns", "ColumnLinkedMatrix", "RowLinkedMatrix")
    vmode <- ifelse(typeof(dataType) == "integer", "byte", "double")

    dims <- pedDims(fileIn = fileIn, header = header, n = n, p = p, nColSkip = nColSkip)

    # Determine chunk size and number of nodes
    if (is.null(nNodes)) {
        if (class == "RowLinkedMatrix") {
            chunkSize <- min(dims$n, floor(.Machine$integer.max / dims$p / 1.2))
            nNodes <- ceiling(dims$n / chunkSize)
        } else {
            chunkSize <- min(dims$p, floor(.Machine$integer.max / dims$n / 1.2))
            nNodes <- ceiling(dims$p / chunkSize)
        }
    } else {
        if (class == "RowLinkedMatrix") {
            chunkSize <- ceiling(dims$n / nNodes)
            if (chunkSize * dims$p >= .Machine$integer.max / 1.2) {
              stop("More nodes are needed")
            }
        } else {
            chunkSize <- ceiling(dims$p / nNodes)
            if (chunkSize * dims$n >= .Machine$integer.max / 1.2) {
              stop("More nodes are needed")
            }
        }
    }

    # Determine dimorder for ff
    if (is.null(dimorder)) {
        if (class == "RowLinkedMatrix") {
            dimorder <- 2:1
        } else {
            dimorder <- 1:2
        }
    }

    # Create output directory
    dir.create(folderOut)

    # Initialize list
    geno <- new(class)
    end <- 0
    if (class == "RowLinkedMatrix") {
        for (i in 1:nNodes) {
            ini <- end + 1
            end <- min(dims$n, ini + chunkSize - 1)
            filename <- paste0("geno_", i, ".bin")
            geno[[i]] <- ff(vmode = vmode, dim = c((end - ini + 1), dims$p), dimorder = dimorder, filename = paste0(folderOut, .Platform$file.sep, filename))
            # Change ff path to a relative one
            physical(geno[[i]])$pattern <- "ff"
            physical(geno[[i]])$filename <- filename
        }
    } else {
        for (i in 1:nNodes) {
            ini <- end + 1
            end <- min(dims$p, ini + chunkSize - 1)
            filename <- paste0("geno_", i, ".bin")
            geno[[i]] <- ff(vmode = vmode, dim = c(dims$n, (end - ini + 1)), dimorder = dimorder, filename = paste0(folderOut, .Platform$file.sep, filename))
            # Change ff path to a relative one
            physical(geno[[i]])$pattern <- "ff"
            physical(geno[[i]])$filename <- filename
        }
    }

    # Generate nodes
    nodes <- LinkedMatrix::nodes(geno)

    # Generate index
    index <- LinkedMatrix::index(geno)

    # Prepare pheno
    pheno <- as.data.frame(matrix(nrow = dims$n, ncol = nColSkip))

    # Construct BGData object
    BGData <- new("BGData", geno = geno, pheno = pheno)

    # Parse PED file
    BGData <- parsePED(BGData = BGData, fileIn = fileIn, header = header, dataType = dataType, nColSkip = nColSkip, na.strings = na.strings, verbose = verbose, nodes = nodes, index = index)

    # Save BGData object
    attr(BGData, "origFile") <- list(path = fileIn, dataType = typeof(dataType))
    attr(BGData, "dateCreated") <- date()
    save(BGData, file = paste(folderOut, "/BGData.RData", sep = ""))

    return(BGData)
}


#' Creates a \code{\linkS4class{BGData}} object from a plaintext PED-like file.
#' 
#' \code{readPED.matrix} assumes that the plaintext file (\code{fileIn}) 
#' contains records of individuals in rows, and phenotypes, covariates and 
#' markers in columns. The columns included in columns \code{1:nColSkip} are 
#' used to populate the slot \code{@@pheno} of a \code{\linkS4class{BGData}}
#' object, and the remaining columns are used to fill the slot \code{@@geno}. If
#' the first row contains a header (\code{header=TRUE}), data in this row is
#' used to determine variables names for \code{@@pheno} and marker names for
#' \code{@@map} and \code{@@geno}.
#' 
#' @param fileIn The path to the plaintext file.
#' @param header If TRUE, the file contains a header.
#' @param dataType The coding of genotypes. Use \code{character()} for A/C/G/T 
#'   or \code{integer()} for numeric coding.
#' @param n The number of individuals.
#' @param p The number of markers.
#' @param na.strings The character string used in the plaintext file to denote 
#'   missing value.
#' @param nColSkip The number of columns to be skipped to reach the genotype 
#'   information in the file.
#' @param idCol The index of the ID column. If more than one index is given, 
#'   both columns will be concatenated with "_".
#' @param verbose If TRUE, progress updates will be posted.
#' @return Returns a \code{\linkS4class{BGData}} object.
#' @seealso \code{\linkS4class{BGData}}
#' @export
readPED.matrix <- function(fileIn, header, dataType, n = NULL, p = NULL, na.strings = "NA", nColSkip = 6, idCol = c(1, 2), verbose = FALSE) {

    dims <- pedDims(fileIn = fileIn, header = header, n = n, p = p, nColSkip = nColSkip)

    dataType <- normalizeType(dataType)

    # Prepare geno
    geno <- matrix(nrow = dims$n, ncol = dims$p)

    # Prepare pheno
    pheno <- as.data.frame(matrix(nrow = dims$n, ncol = nColSkip))

    # Construct BGData object
    BGData <- new("BGData", geno = geno, pheno = pheno)

    # Parse PED file
    BGData <- parsePED(BGData = BGData, fileIn = fileIn, header = header, dataType = dataType, nColSkip = nColSkip, na.strings = na.strings, verbose = verbose)

    return(BGData)
}


#' Creates a \code{\linkS4class{BGData}} object from a plaintext PED-like file.
#' 
#' @param fileIn The path to the plaintext file.
#' @param header If TRUE, the file contains a header.
#' @param dataType The coding of genotypes. Use \code{character()} for A/C/G/T 
#'   or \code{integer()} for numeric coding.
#' @param n The number of individuals.
#' @param p The number of markers.
#' @param na.strings The character string used in the plaintext file to denote 
#'   missing value.
#' @param nColSkip The number of columns to be skipped to reach the genotype 
#'   information in the file.
#' @param idCol The index of the ID column. If more than one index is given, 
#'   both columns will be concatenated with "_".
#' @param verbose If TRUE, progress updates will be posted.
#' @param folderOut The path to the folder where to save the binary files.
#' @return Returns a \code{\linkS4class{BGData}} object.
#' @seealso \code{\linkS4class{BGData}}
#' @export
readPED.big.matrix <- function(fileIn, header, dataType, n = NULL, p = NULL, na.strings = "NA", nColSkip = 6, idCol = c(1, 2), verbose = FALSE, folderOut = paste("BGData_", sub("\\.[[:alnum:]]+$", "", basename(fileIn)), sep = "")) {

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

    dims <- pedDims(fileIn = fileIn, header = header, n = n, p = p, nColSkip = nColSkip)

    options(bigmemory.typecast.warning = FALSE)
    options(bigmemory.allow.dimnames = TRUE)

    # Create output directory
    dir.create(folderOut)

    # Prepare geno
    geno <- bigmemory::filebacked.big.matrix(nrow = dims$n, ncol = dims$p, type = type, backingpath = folderOut, backingfile = "BGData.bin", descriptorfile = "BGData.desc")

    # Prepare pheno
    pheno <- as.data.frame(matrix(nrow = dims$n, ncol = nColSkip))

    # Construct BGData object
    BGData <- new("BGData", geno = geno, pheno = pheno)

    # Parse PED file
    BGData <- parsePED(BGData = BGData, fileIn = fileIn, header = header, dataType = dataType, nColSkip = nColSkip, na.strings = na.strings, verbose = verbose)

    # Save BGData object
    attr(BGData, "origFile") <- list(path = fileIn, dataType = typeof(dataType))
    attr(BGData, "dateCreated") <- date()
    save(BGData, file = paste(folderOut, "/BGData.RData", sep = ""))

    return(BGData)
}


#' Loads BGData objects.
#' 
#' If the value in the \code{@@geno} slot is of class 
#' \code{\linkS4class{LinkedMatrix}} and the nodes contain \code{ff} objects 
#' \code{load.BGData} will attempt to open them. If the value is of class 
#' \code{\linkS4class{big.matrix}}, it will attempt to attach the matrix.
#' 
#' Genotypes are stored in a filebacked \code{\linkS4class{big.matrix}} object. 
#' \code{readPED.big.matrix} creates a folder (see \code{folderOut}) that
#' contains the binary flat file (named \code{BGData.bin}), a descriptor file
#' (named \code{BGData.desc}), and an external representation of the
#' \code{\linkS4class{BGData}} object in \code{BGData.RData}. A 
#' \code{\linkS4class{BGData}} object can be reloaded using \code{load.BGData} 
#' (the regular \code{load} function will not work unless the
#' \code{\linkS4class{big.matrix}} instance is manually attached using
#' \code[bigmemory]{attach.big.matrix}).
#' 
#' @param file The name of the .RData file to be loaded.
#' @param envir The environment where to load the data.
#' @export
load.BGData <- function(file, envir = parent.frame()) {

    # Determine object name and class
    lsOLD <- ls()
    load(file = file)
    lsNEW <- ls()
    objectName <- lsNEW[(!lsNEW %in% lsOLD) & (lsNEW != "lsOLD")]
    object <- get(objectName)
    objectClass <- class(object)

    if (objectClass != "BGData") {
        stop("Object class must be BGData")
    }

    message(paste0("Loaded object ", objectName, " of class ", objectClass))

    if (inherits(class(object@geno), "LinkedMatrix")) {

        # Store current working directory and set working directory to directory of file
        cwd <- getwd()
        setwd(dirname(file))

        # Open all nodes for reading (we do not store absolute paths to ff files, so this
        # has to happen in the same working directory)
        nNodes <- LinkedMatrix::nNodes(object@geno)
        for (i in seq_len(nNodes)) {
            node <- object@geno[[i]]
            if (inherits(node, "ff_matrix")) {
                message(paste0("Opening flat file ", i, "..."))
                open(node)
            }
        }

        # Restore the working directory
        setwd(cwd)

    } else if (class(object@geno) == "big.matrix") {
        object@geno <- bigmemory::attach.big.matrix(paste0(dirname(file), "/", "BGData.desc"))
    }

    # Send the object to envir
    assign(objectName, object, envir = envir)
}
