# Load example data
bg <- BGData:::loadExample()

# Compute a scaled genomic relationship matrix from centered and scaled
# genotypes
g1 <- getG(X = bg@geno)

# Disable scaling of G
g2 <- getG(X = bg@geno, scaleG = FALSE)

# Disable centering of genotypes
g3 <- getG(X = bg@geno, center = FALSE)

# Disable scaling of genotypes
g4 <- getG(X = bg@geno, scale = FALSE)

# Provide own scales
scales <- chunkedApply(X = bg@geno, MARGIN = 2, FUN = sd)
g4 <- getG(X = bg@geno, scale = scales)

# Provide own centers
centers <- chunkedApply(X = bg@geno, MARGIN = 2, FUN = mean)
g5 <- getG(X = bg@geno, center = centers)

# Only use the first 50 individuals (useful to account for population structure)
g6 <- getG(X = bg@geno, i = 1:50)

# Only use the first 100 markers (useful to ignore some markers)
g7 <- getG(X = bg@geno, j = 1:100)

# Compute unscaled G matrix by combining blocks of $XX_{i2}'$ where $X_{i2}$ is
# a horizontal partition of X. This is useful for distributed computing as each
# block can be computed in parallel. Centers and scales need to be precomputed.
block1 <- getG(X = bg@geno, i2 = 1:100, center = centers, scale = scales)
block2 <- getG(X = bg@geno, i2 = 101:199, center = centers, scale = scales)
g8 <- cbind(block1, block2)

# Compute unscaled G matrix by combining blocks of $X_{i}X_{i2}'$ where both
# $X_{i}$ and $X_{i2}$ are horizontal partitions of X. Similarly to the example
# above, this is useful for distributed computing, in particular to compute
# very large G matrices. Centers and scales need to be precomputed. This
# approach is similar to the one taken by the symDMatrix package, but the
# symDMatrix package adds memory-mapped blocks, only stores the upper side of
# the triangular matrix, and provides a type that allows for indexing as if the
# full G matrix is in memory.
block11 <- getG(X = bg@geno, i = 1:100, i2 = 1:100, center = centers, scale = scales)
block12 <- getG(X = bg@geno, i = 1:100, i2 = 101:199, center = centers, scale = scales)
block21 <- getG(X = bg@geno, i = 101:199, i2 = 1:100, center = centers, scale = scales)
block22 <- getG(X = bg@geno, i = 101:199, i2 = 101:199, center = centers, scale = scales)
g9 <- rbind(
    cbind(block11, block12),
    cbind(block21, block22)
)
