#' S3 generic to convert other objects into BGData objects.
#' 
#' @param x An object.
#' @param ... Additional arguments (see method).
#' @return A BGData object.
#' @export
as.BGData <- function(x, ...) {
    UseMethod("as.BGData")
}


#' Convert a BEDMatrix object into a BGData object.
#' 
#' @param x A BEDMatrix object.
#' @param alternatePhenotypeFile Path to an alternate phenotype file (see PLINK documentation).
#' @param ... Additional arguments to the BEDMatrix constructor (currently unused).
#' @return A BGData object.
#' @export
as.BGData.BEDMatrix <- function(x, alternatePhenotypeFile = NULL, ...) {
    # Path to BED file
    bedPath <- attr(x, "path")
    # Path to FAM file
    famPath <- sub(".bed", ".fam", bedPath)
    # Read in pheno file
    if (file.exists(famPath)) {
        message("Extracting phenotypes from FAM file...")
        pheno <- read.table(famPath, col.names = c(
            "Family_ID",
            "Individual_ID",
            "Paternal_ID",
            "Maternal_ID",
            "Sex",
            "Phenotype"
        ), stringsAsFactors = FALSE)
    } else {
        pheno <- data.frame(IID = rownames(x))
        rownames(pheno) <- rownames(x)
    }
    if (!is.null(alternatePhenotypeFile)) {
        if (!file.exists(alternatePhenotypeFile)) {
            stop("Alternate phenotype file does not exist.")
        } else {
            message("Merging alternate phenotype file...")
            # Check if the file has a header, i.e. if the first row starts with
            # an FID and an IID entry (unfortunately, treating the first row as
            # colnames after read.table screws with the types)
            hasHeader = FALSE
            if (grepl("FID\\s+IID", readLines(alternatePhenotypeFile, n = 1))) {
                hasHeader = TRUE
            }
            alternatePhenotypes <- read.table(alternatePhenotypeFile, header = hasHeader, stringsAsFactors = FALSE)
            # Merge phenotypes and alternate phenotypes
            pheno <- merge(pheno, alternatePhenotypes, by.x = c("Family_ID", "Individual_ID"), by.y = c("FID", "IID"), all.x = TRUE)
        }
    }
    BGData(geno = x, pheno = pheno)
}
