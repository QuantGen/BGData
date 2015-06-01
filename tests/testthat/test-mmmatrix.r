library(BGData)

# Prepare dummy data
genotypes <- matrix(c(4, 4, 4, 3, 2, 3, 1, 2, 1), nrow = 3, ncol = 3)
colnames(genotypes) <- paste0('mrk_', 1:3)
rownames(genotypes) <- paste0('id_', 1:3)

for (class in c('cmmMatrix', 'rmmMatrix')) {

    context(class)

    # Prepare mmMatrix object
    geno <- new(class, nrow = 3, ncol = 3)
    geno[] <- genotypes
    colnames(geno) <- paste0('mrk_', 1:3)
    rownames(geno) <- paste0('id_', 1:3)

    test_that("subsetting", {

        expect_true(all.equal(geno[], genotypes))

        # expect_true(all.equal(geno[1], genotypes[1])) Not implemented yet
        expect_true(all.equal(geno[1, ], genotypes[1, ]))
        expect_true(all.equal(geno[, 1], genotypes[, 1]))
        expect_true(all.equal(geno[1, 1], genotypes[1, 1]))
        expect_true(all.equal(geno[1, , drop = FALSE], genotypes[1, , drop = FALSE]))
        expect_true(all.equal(geno[, 1, drop = FALSE], genotypes[, 1, drop = FALSE]))
        expect_true(all.equal(geno[1, 1, drop = FALSE], genotypes[1, 1, drop = FALSE]))

        # expect_true(all.equal(geno[1:2], genotypes[1:2]))  Not implemented yet
        expect_true(all.equal(geno[1:2, ], genotypes[1:2, ]))
        expect_true(all.equal(geno[, 1:2], genotypes[, 1:2]))
        expect_true(all.equal(geno[1:2, 1:2], genotypes[1:2, 1:2]))
        expect_true(all.equal(geno[1:2, , drop = FALSE], genotypes[1:2, , drop = FALSE]))
        expect_true(all.equal(geno[, 1:2, drop = FALSE], genotypes[, 1:2, drop = FALSE]))
        expect_true(all.equal(geno[1:2, 1:2, drop = FALSE], genotypes[1:2, 1:2, drop = FALSE]))

        # expect_true(all.equal(geno[2:1], genotypes[2:1])) Not implemented yet
        expect_true(all.equal(geno[2:1, ], genotypes[2:1, ]))
        expect_true(all.equal(geno[, 2:1], genotypes[, 2:1]))
        expect_true(all.equal(geno[2:1, 2:1], genotypes[2:1, 2:1]))
        expect_true(all.equal(geno[2:1, , drop = FALSE], genotypes[2:1, , drop = FALSE]))
        expect_true(all.equal(geno[, 2:1, drop = FALSE], genotypes[, 2:1, drop = FALSE]))
        expect_true(all.equal(geno[2:1, 2:1, drop = FALSE], genotypes[2:1, 2:1, drop = FALSE]))

        # expect_true(all.equal(geno[c(3, 1)], genotypes[c(3, 1)])) Not implemented yet
        expect_true(all.equal(geno[c(3, 1), ], genotypes[c(3, 1), ]))
        expect_true(all.equal(geno[, c(3, 1)], genotypes[, c(3, 1)]))
        expect_true(all.equal(geno[c(3, 1), c(3, 1)], genotypes[c(3, 1), c(3, 1)]))
        expect_true(all.equal(geno[c(3, 1), , drop = FALSE], genotypes[c(3, 1), , drop = FALSE]))
        expect_true(all.equal(geno[, c(3, 1), drop = FALSE], genotypes[, c(3, 1), drop = FALSE]))
        expect_true(all.equal(geno[c(3, 1), c(3, 1), drop = FALSE], genotypes[c(3, 1), c(3, 1), drop = FALSE]))

        # expect_true(all.equal(geno[genotypes > 1], genotypes[genotypes > 1])) Not implemented yet
        expect_true(all.equal(geno[c(TRUE, FALSE, TRUE), ], genotypes[c(TRUE, FALSE, TRUE), ]))
        expect_true(all.equal(geno[, c(TRUE, FALSE, TRUE)], genotypes[, c(TRUE, FALSE, TRUE)]))
        expect_true(all.equal(geno[c(TRUE, FALSE, TRUE), c(TRUE, FALSE, TRUE)], genotypes[c(TRUE, FALSE, TRUE), c(TRUE, FALSE, TRUE)]))
        expect_true(all.equal(geno[c(TRUE, FALSE, TRUE), , drop = FALSE], genotypes[c(TRUE, FALSE, TRUE), , drop = FALSE]))
        expect_true(all.equal(geno[, c(TRUE, FALSE, TRUE), drop = FALSE], genotypes[, c(TRUE, FALSE, TRUE), drop = FALSE]))
        expect_true(all.equal(geno[c(TRUE, FALSE, TRUE), c(TRUE, FALSE, TRUE), drop = FALSE], genotypes[c(TRUE, FALSE, TRUE), c(TRUE, FALSE, TRUE), drop = FALSE]))

        expect_true(all.equal(geno['id_1', ], genotypes['id_1', ]))
        expect_true(all.equal(geno[, 'mrk_1'], genotypes[, 'mrk_1']))
        expect_true(all.equal(geno['id_1', 'mrk_1'], genotypes['id_1', 'mrk_1']))
        expect_true(all.equal(geno['id_1', , drop = FALSE], genotypes['id_1', , drop = FALSE]))
        expect_true(all.equal(geno[, 'mrk_1', drop = FALSE], genotypes[, 'mrk_1', drop = FALSE]))
        expect_true(all.equal(geno['id_1', 'mrk_1', drop = FALSE], genotypes['id_1', 'mrk_1', drop = FALSE]))

        expect_true(all.equal(geno[c('id_1', 'id_2'), ], genotypes[c('id_1', 'id_2'), ]))
        expect_true(all.equal(geno[, c('mrk_1', 'mrk_2')], genotypes[, c('mrk_1', 'mrk_2')]))
        expect_true(all.equal(geno[c('id_1', 'id_2'), c('mrk_1', 'mrk_2')], genotypes[c('id_1', 'id_2'), c('mrk_1', 'mrk_2')]))
        expect_true(all.equal(geno[c('id_1', 'id_2'), , drop = FALSE], genotypes[c('id_1', 'id_2'), , drop = FALSE]))
        expect_true(all.equal(geno[, c('mrk_1', 'mrk_2'), drop = FALSE], genotypes[, c('mrk_1', 'mrk_2'), drop = FALSE]))
        expect_true(all.equal(geno[c('id_1', 'id_2'), c('mrk_1', 'mrk_2'), drop = FALSE], genotypes[c('id_1', 'id_2'), c('mrk_1', 'mrk_2'), drop = FALSE]))

        expect_true(all.equal(geno[c('id_2', 'id_1'), ], genotypes[c('id_2', 'id_1'), ]))
        expect_true(all.equal(geno[, c('mrk_2', 'mrk_1')], genotypes[, c('mrk_2', 'mrk_1')]))
        expect_true(all.equal(geno[c('id_2', 'id_1'), c('mrk_2', 'mrk_1')], genotypes[c('id_2', 'id_1'), c('mrk_2', 'mrk_1')]))
        expect_true(all.equal(geno[c('id_2', 'id_1'), , drop = FALSE], genotypes[c('id_2', 'id_1'), , drop = FALSE]))
        expect_true(all.equal(geno[, c('mrk_2', 'mrk_1'), drop = FALSE], genotypes[, c('mrk_2', 'mrk_1'), drop = FALSE]))
        expect_true(all.equal(geno[c('id_2', 'id_1'), c('mrk_2', 'mrk_1'), drop = FALSE], genotypes[c('id_2', 'id_1'), c('mrk_2', 'mrk_1'), drop = FALSE]))

        expect_true(all.equal(geno[c('id_3', 'id_1'), ], genotypes[c('id_3', 'id_1'), ]))
        expect_true(all.equal(geno[, c('mrk_3', 'mrk_1')], genotypes[, c('mrk_3', 'mrk_1')]))
        expect_true(all.equal(geno[c('id_3', 'id_1'), c('mrk_3', 'mrk_1')], genotypes[c('id_3', 'id_1'), c('mrk_3', 'mrk_1')]))
        expect_true(all.equal(geno[c('id_3', 'id_1'), , drop = FALSE], genotypes[c('id_3', 'id_1'), , drop = FALSE]))
        expect_true(all.equal(geno[, c('mrk_3', 'mrk_1'), drop = FALSE], genotypes[, c('mrk_3', 'mrk_1'), drop = FALSE]))
        expect_true(all.equal(geno[c('id_3', 'id_1'), c('mrk_3', 'mrk_1'), drop = FALSE], genotypes[c('id_3', 'id_1'), c('mrk_3', 'mrk_1'), drop = FALSE]))

    })

    test_that("replacement", {

        # Generate new genotypes for replacement
        replacement <- matrix(c(3, 1, 3, 2, 4, 3, 1, 1, 2), nrow = 3, ncol = 3)
        colnames(replacement) <- paste0('mrk_', 1:3)
        rownames(replacement) <- paste0('id_', 1:3)
        comparison <- genotypes

        testAndRestore <- function (label) {
            expect_equal(all.equal(geno[], comparison), TRUE, label=label)
            geno[] <- genotypes # no environment change necessary
            assign('comparison', genotypes, parent.frame())
        }

        geno[] <- replacement
        comparison[] <- replacement
        testAndRestore('[]')

        geno[1, ] <- replacement[1, ]
        comparison[1, ] <- replacement[1, ]
        testAndRestore('[1, ]')
        geno[, 1] <- replacement[, 1]
        comparison[, 1] <- replacement[, 1]
        testAndRestore('[, 1]')
        geno[1, 1] <- replacement[1, 1]
        comparison[1, 1] <- replacement[1, 1]
        testAndRestore('[1, 1]')

        geno[1:2, ] <- replacement[1:2, ]
        comparison[1:2, ] <- replacement[1:2, ]
        testAndRestore('[1:2, ]')
        geno[, 1:2] <- replacement[, 1:2]
        comparison[, 1:2] <- replacement[, 1:2]
        testAndRestore('[, 1:2]')
        geno[1:2, 1:2] <- replacement[1:2, 1:2]
        comparison[1:2, 1:2] <- replacement[1:2, 1:2]
        testAndRestore('[1:2, 1:2]')

        geno[2:1, ] <- replacement[2:1, ]
        comparison[2:1, ] <- replacement[2:1, ]
        testAndRestore('[2:1, ]')
        geno[, 2:1] <- replacement[, 2:1]
        comparison[, 2:1] <- replacement[, 2:1]
        testAndRestore('[, 2:1]')
        geno[2:1, 2:1] <- replacement[2:1, 2:1]
        comparison[2:1, 2:1] <- replacement[2:1, 2:1]
        testAndRestore('[2:1, 2:1]')

        geno[c(3, 1), ] <- replacement[c(3, 1), ]
        comparison[c(3, 1), ] <- replacement[c(3, 1), ]
        testAndRestore('[c(3, 1), ]')
        geno[, c(3, 1)] <- replacement[, c(3, 1)]
        comparison[, c(3, 1)] <- replacement[, c(3, 1)]
        testAndRestore('[, c(3, 1)]')
        geno[c(3, 1), c(3, 1)] <- replacement[c(3, 1), c(3, 1)]
        comparison[c(3, 1), c(3, 1)] <- replacement[c(3, 1), c(3, 1)]
        testAndRestore('[c(3, 1), c(3, 1)]')

    })

    test_that("dim", {
        expect_equal(dim(geno), dim(genotypes))
    })

    test_that("apply", {

        expect_true(all.equal(colMeans(geno), colMeans(genotypes)))
        expect_true(all.equal(colSums(geno), colSums(genotypes)))
        expect_true(all.equal(rowMeans(geno), rowMeans(genotypes)))
        expect_true(all.equal(rowSums(geno), rowSums(genotypes)))

        # Introduce NA
        genotypes_na <- genotypes
        genotypes_na[1,1] <- NA
        geno[1,1] <- NA

        expect_warning(colMeans(geno))
        expect_warning(colSums(geno))
        expect_warning(rowMeans(geno))
        expect_warning(rowSums(geno))

        expect_true(all.equal(colMeans(geno, na.rm=FALSE), colMeans(genotypes_na, na.rm=FALSE)))
        expect_true(all.equal(colSums(geno, na.rm=FALSE), colSums(genotypes_na, na.rm=FALSE)))
        expect_true(all.equal(rowMeans(geno, na.rm=FALSE), rowMeans(genotypes_na, na.rm=FALSE)))
        expect_true(all.equal(rowSums(geno, na.rm=FALSE), rowSums(genotypes_na, na.rm=FALSE)))

        expect_true(all.equal(colMeans(geno, na.rm=TRUE), colMeans(genotypes_na, na.rm=TRUE)))
        expect_true(all.equal(colSums(geno, na.rm=TRUE), colSums(genotypes_na, na.rm=TRUE)))
        expect_true(all.equal(rowMeans(geno, na.rm=TRUE), rowMeans(genotypes_na, na.rm=TRUE)))
        expect_true(all.equal(rowSums(geno, na.rm=TRUE), rowSums(genotypes_na, na.rm=TRUE)))

        # Revert NA
        geno[] <- genotypes

    })

}
