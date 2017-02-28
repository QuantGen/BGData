# Load example data
bg <- BGData:::loadExample()

# Compute xx' in parallel
tcrossprod_parallel(x = bg@geno)

# Compute xy' in parallel (see getG)
tcrossprod_parallel(x = bg@geno, y = bg@geno[1:50, ])

# Compute x'x in parallel
crossprod_parallel(x = bg@geno)
