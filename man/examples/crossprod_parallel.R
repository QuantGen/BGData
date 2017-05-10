# Restrict number of cores to 1 on Windows
if (.Platform$OS.type == "windows") {
    options(mc.cores = 1)
}

# Load example data
bg <- BGData:::loadExample()

# Compute xx' in parallel
tcrossprod_parallel(x = bg@geno)

# Compute xy' in parallel (see getG)
tcrossprod_parallel(x = bg@geno, y = bg@geno[1:50, ])

# Compute x'x in parallel
crossprod_parallel(x = bg@geno)
