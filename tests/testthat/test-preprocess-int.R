context("preprocess for integers")

# Parameters
n <- 250
p <- 50
length <- n * p
nmiss <- 100

# Data
set.seed(4711)
X <- sample(0:9, size = length, replace = TRUE)
dim(X) <- c(n, p)
missing <- sample(seq_len(length), size = nmiss)
X[missing] <- NA

centers <- colMeans(X, na.rm = TRUE)
scales <- apply(X, 2, sd, na.rm = TRUE)

# No operation
expect_equal(
    scale(X, center = FALSE, scale = FALSE),
    preprocess(X, center = FALSE, scale = FALSE, impute = FALSE)
)

# Tests without imputation

# Compute centers and scales
expect_equal(
    scale(X, center = TRUE, scale = TRUE),
    preprocess(X, center = TRUE, scale = TRUE, impute = FALSE)
)
expect_equal(
    scale(X, center = TRUE, scale = FALSE),
    preprocess(X, center = TRUE, scale = FALSE, impute = FALSE)
)
expect_equal(
    scale(X, center = FALSE, scale = scales), # scale() uses root mean squares if 'center = FALSE'
    preprocess(X, center = FALSE, scale = TRUE, impute = FALSE)
)

# Provide own centers and scales
expect_equal(
    scale(X, center = centers, scale = scales),
    preprocess(X, center = centers, scale = scales, impute = FALSE)
)
expect_equal(
    scale(X, center = centers, scale = FALSE),
    preprocess(X, center = centers, scale = FALSE, impute = FALSE)
)
expect_equal(
    scale(X, center = FALSE, scale = scales),
    preprocess(X, center = FALSE, scale = scales, impute = FALSE)
)

# Provide own centers, compute scales
expect_equal(
    scale(X, center = centers, scale = TRUE),
    preprocess(X, center = centers, scale = TRUE, impute = FALSE)
)

# Provide own scales, compute centers
expect_equal(
    scale(X, center = TRUE, scale = scales),
    preprocess(X, center = TRUE, scale = scales, impute = FALSE)
)


# Tests with imputation

# center = TRUE and impute = TRUE means impute by 0
expect_equal(
    {
        W <- scale(X, center = TRUE, scale = FALSE)
        W[missing] <- 0
        W
    },
    preprocess(X, center = TRUE, scale = FALSE, impute = TRUE)
)

# Given centers and impute = TRUE means impute by 0
expect_equal(
    {
        W <- scale(X, center = centers, scale = FALSE)
        W[missing] <- 0
        W
    },
    preprocess(X, center = centers, scale = FALSE, impute = TRUE)
)

# center = FALSE and impute = TRUE means impute by mean
expect_equal(
    {
        means <- rep(colMeans(X, na.rm = TRUE), each = n)
        W <- X
        W[missing] <- means[missing]
        W
    },
    preprocess(X, center = FALSE, scale = FALSE, impute = TRUE)
)
