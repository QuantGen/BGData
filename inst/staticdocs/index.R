sd_section(
    "Classes",
    "",
    c(
        "BGData-class"
    )
)

sd_section(
    "Creating a BGData object from a PED (or PED-like) file",
    "",
    c(
        "readPED",
        "readPED.big.matrix",
        "readPED.matrix"
    )
)

sd_section(
    "Creating a BGData object from a BED file",
    "The BED file needs to be wrapped by a BEDMatrix first.",
    c(
        "as.BGData.BEDMatrix"
    )
)

sd_section(
    "Creating a BGData object from a previously saved BGData object",
    "",
    c(
        "load.BGData"
    )
)

sd_section(
    "Creating a BGData object from arbitrary genotype/phenotype data",
    "",
    c(
        "initialize,BGData-method"
    )
)

sd_section(
    "Functions",
    "",
    c(
        "GWAS",
        "getG",
        "getG.symDMatrix",
        "summarize",
        "simPED",
        "chunkedApply",
        "parallelApply",
        "crossprod.parallel",
        "tcrossprod.parallel"
    )
)
