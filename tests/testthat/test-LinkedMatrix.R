# Prepare dummy data
genotypes <- matrix(c(4, 4, 4, 3, 2, 3, 1, 2, 1), nrow = 3, ncol = 3)
rownames(genotypes) <- paste0("id_", 1:nrow(genotypes))
colnames(genotypes) <- paste0("mrk_", 1:ncol(genotypes))

nodes <- function(dim, nNodes) {
    chunkSizes <- vector(mode = "integer", length = nNodes)
    for (node in 1:nNodes) {
        chunkSizes[node] <- floor(dim / nNodes)
    }
    for (node in seq_len(dim %% nNodes)) {
        chunkSizes[node] <- chunkSizes[node] + 1
    }
    nodes <- matrix(nrow = nNodes, ncol = 2, NA)
    idx <- 1
    for (node in 1:nNodes) {
        nodes[node, 1] <- idx
        idx <- idx + chunkSizes[node]
        nodes[node, 2] <- idx - 1
    }
    return(nodes)
}

createLinkedMatrix <- function(class, nNodes) {
    linkedMatrix <- new(class)
    nodes <- nodes(ifelse(class == "ColumnLinkedMatrix", ncol(genotypes), nrow(genotypes)), nNodes)
    for (node in 1:nNodes) {
        if (class == "ColumnLinkedMatrix") {
            linkedMatrix[[node]] <- genotypes[, nodes[node, 1]:nodes[node, 2], drop = FALSE]
        } else {
            linkedMatrix[[node]] <- genotypes[nodes[node, 1]:nodes[node, 2], , drop = FALSE]
        }
    }
    linkedMatrix[] <- genotypes
    rownames(linkedMatrix) <- paste0("id_", 1:nrow(genotypes))
    colnames(linkedMatrix) <- paste0("mrk_", 1:ncol(genotypes))
    return(linkedMatrix)
}

for (class in c("ColumnLinkedMatrix", "RowLinkedMatrix")) {

    context(class)

    for (nNodes in 1:ifelse(class == "ColumnLinkedMatrix", ncol(genotypes), nrow(genotypes))) {

        context(paste0(class, " with ", nNodes, " nodes"))

        # Prepare LinkedMatrix object
        linkedMatrix <- createLinkedMatrix(class, nNodes)

        test_that("apply", {
            expect_equal(apply(linkedMatrix, 1, sum), base::apply(genotypes, 1, sum))
            expect_equal(apply(linkedMatrix, 2, sum), base::apply(genotypes, 2, sum))
            expect_equal(apply(linkedMatrix, 1, summary), base::apply(genotypes, 1, summary))
            expect_equal(apply(linkedMatrix, 2, summary), base::apply(genotypes, 2, summary))
        })

        test_that("colMeans", {
            expect_equal(colMeans(linkedMatrix), base::colMeans(genotypes))

            # Introduce NA
            genotypes_na <- genotypes
            genotypes_na[1, 1] <- NA
            linkedMatrix[1, 1] <- NA

            expect_warning(colMeans(linkedMatrix))
            expect_equal(colMeans(linkedMatrix, na.rm = FALSE), base::colMeans(genotypes_na, na.rm = FALSE))
            expect_equal(colMeans(linkedMatrix, na.rm = TRUE), base::colMeans(genotypes_na, na.rm = TRUE))

            # Revert NA
            linkedMatrix[] <- genotypes
        })

        test_that("colSums", {
            expect_equal(colSums(linkedMatrix), base::colSums(genotypes))

            # Introduce NA
            genotypes_na <- genotypes
            genotypes_na[1, 1] <- NA
            linkedMatrix[1, 1] <- NA

            expect_warning(colSums(linkedMatrix))
            expect_equal(colSums(linkedMatrix, na.rm = FALSE), base::colSums(genotypes_na, na.rm = FALSE))
            expect_equal(colSums(linkedMatrix, na.rm = TRUE), base::colSums(genotypes_na, na.rm = TRUE))

            # Revert NA
            linkedMatrix[] <- genotypes
        })

        test_that("rowMeans", {
            expect_equal(rowMeans(linkedMatrix), base::rowMeans(genotypes))

            # Introduce NA
            genotypes_na <- genotypes
            genotypes_na[1, 1] <- NA
            linkedMatrix[1, 1] <- NA
            expect_warning(rowMeans(linkedMatrix))
            expect_equal(rowMeans(linkedMatrix, na.rm = FALSE), base::rowMeans(genotypes_na, na.rm = FALSE))
            expect_equal(rowMeans(linkedMatrix, na.rm = TRUE), base::rowMeans(genotypes_na, na.rm = TRUE))

            # Revert NA
            linkedMatrix[] <- genotypes
        })

        test_that("rowSums", {
            expect_equal(rowSums(linkedMatrix), base::rowSums(genotypes))

            # Introduce NA
            genotypes_na <- genotypes
            genotypes_na[1, 1] <- NA
            linkedMatrix[1, 1] <- NA

            expect_warning(rowSums(linkedMatrix))
            expect_equal(rowSums(linkedMatrix, na.rm = FALSE), base::rowSums(genotypes_na, na.rm = FALSE))
            expect_equal(rowSums(linkedMatrix, na.rm = TRUE), base::rowSums(genotypes_na, na.rm = TRUE))

            # Revert NA
            linkedMatrix[] <- genotypes
        })

    }

} 
