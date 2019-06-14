getLineCount <- function(path, header) {
    file <- file(path, open = "r")
    n <- 0L
    while (length(readLines(file, n = 1L)) > 0L) {
        n <- n + 1L
    }
    if (header) {
        n <- n - 1L
    }
    close(file)
    return(n)
}

getFileHeader <- function(path, sep = "") {
    file <- file(path, open = "r")
    header <- scan(file, nlines = 1L, what = character(), sep = sep, quiet = TRUE)
    close(file)
    return(header)
}

getColumnCount <- function(path, sep = "") {
    header <- getFileHeader(path, sep)
    p <- length(header)
    return(p)
}

randomString <- function() {
    paste(sample(c(0L:9L, letters, LETTERS), size = 5L, replace = TRUE), collapse = "")
}

normalizeType <- function(val) {
    type <- typeof(val)
    # detect strings
    if (type == "character" && length(val) > 0L) {
        # convert to type if type and value match
        convert <- try(vector(mode = val), silent = TRUE)
        if (class(convert) == "try-error") {
            # return a character type if conversion failed
            warning("could no convert type, using character instead")
            character()
        } else {
            # return conversion result otherwise
            convert
        }
        # value doesn't contain type information and can be handled by typeof
    } else {
        val
    }
}

loadExample <- function() {
    path <- system.file("extdata", package = "BGData")
    message("Loading chromosomes as BED files...")
    m <- do.call(LinkedMatrix::ColumnLinkedMatrix, lapply(c("chr1", "chr2", "chr3"), function(chr) {
        suppressMessages(BEDMatrix::BEDMatrix(paste0(path, "/", chr)))
    }))
    as.BGData(m, alternatePhenotypeFile = paste0(path, "/pheno.txt"))
}
