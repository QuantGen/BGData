# Restrict number of cores to 1 on Windows
if (.Platform$OS.type == "windows") {
    options(mc.cores = 1)
}

# Load example data
bg <- BGData:::loadExample()

# Summarize the whole dataset
sum1 <- summarize(X = bg@geno)

# Summarize the first 50 individuals
sum2 <- summarize(X = bg@geno, i = 1:50)

# Summarize the first 1000 markers (useful for distributed computing)
sum3 <- summarize(X = bg@geno, j = 1:100)

# Summarize the first 50 individuals on the first 1000 markers
sum4 <- summarize(X = bg@geno, i = 1:50, j = 1:100)

# Summarize by names
sum5 <- summarize(X = bg@geno, j = c("snp81233_C", "snp81234_C", "snp81235_T"))
