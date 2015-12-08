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
    BGData(geno = x)
}
