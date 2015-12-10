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
#' @param ... Additional arguments to the BEDMatrix constructor (currently unused).
#' @return A BGData object.
#' @export
as.BGData.BEDMatrix <- function(x, ...) {
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
    BGData(geno = x, pheno = pheno)
}
