preprocess <- function(X, center = FALSE, scale = FALSE, impute = FALSE) {
    .Call(C_preprocess, X, center, scale, impute, FALSE)
}
