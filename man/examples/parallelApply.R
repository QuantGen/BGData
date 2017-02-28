# Load example data
bg <- BGData:::loadExample()

# Compute standard deviation of columns
parallelApply(X = bg@geno, MARGIN = 2, FUN = sd)
