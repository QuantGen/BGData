BGData: A Suite of Packages for Analysis of Big Genomic Data
============================================================

[![Travis-CI Build Status](https://travis-ci.org/QuantGen/BGData.svg?branch=master)](https://travis-ci.org/QuantGen/BGData)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/BGData)](https://cran.r-project.org/package=BGData)
[![Rdoc](http://www.rdocumentation.org/badges/version/BGData)](http://www.rdocumentation.org/packages/BGData)

**Contact**: gruenebe@msu.edu gustavoc@msu.edu

Modern genomic datasets are big (large *n*), high-dimensional (large *p*), and multi-layered. The challenges that need to be addressed are memory requirements and computational demands. Our goal is to develop software that will enable researchers to carry out analyses with big genomic data within the R environment.

We have identified several approaches to tackle those challanges within R:

- Memory mapping: The data is stored in on the hard drive and users can read in smaller chunks when they are needed.
- Linked arrays: For very large datasets a single memory-mapped array may not be enough or convenient. A linked array is an array whose content is distributed over multiple memory-mapped nodes.
- Multiple dispatch: Methods are presented to users so that they can treat these arrays pretty much as if they were RAM arrays.
- Multi-level parallelism: Exploit multi-core and multi-node computing.
- Inputs: Users can create these arrays from standard formats (e.g., PED, BED).

The BGData package is an umbrella package that comprises several packages: [BEDMatrix](https://github.com/QuantGen/BEDMatrix), [LinkedMatrix](https://github.com/QuantGen/LinkedMatrix), and [symDMatrix](https://github.com/QuantGen/symDMatrix). It features scalable and efficient computational methods for large genomic datasets such as genome-wide association studies (GWAS) or genomic relationship matrices (G matrix). It also contains a data structure called `BGData` that holds genotypes in the `@geno` slot, phenotypes in the `@pheno` slot, and additional information in the `@map` slot.


Installation
------------

The BGData package is not available on [CRAN](http://cran.r-project.org/) yet, but it can be installed directly from GitHub using the [devtools](https://github.com/hadley/devtools) package:

    # install.packages("devtools")
    devtools::install_github("QuantGen/BGData")


Examples
--------

### Creating a BGData object from a PLINK BED file

This example uses a [BED version](https://www.cog-genomics.org/plink2/formats#bed) of the `mice` dataset that is included in the BGLR package. See the [mice.bed gist](https://gist.github.com/agrueneberg/812564cbe860db4ee864d019be940aaf) for instructions on how it was generated.

Load the BGData package:

    > library(BGData)

Load the `mice` dataset as BEDMatrix:

    > bed <- BEDMatrix(system.file("extdata", "mice.bed", package = "BEDMatrix"))
    Extracting number of individuals and rownames from FAM file...
    Extracting number of markers and colnames from BIM file...

Display structure of BEDMatrix object:

    > str(bed)
    BEDMatrix: 1814 x 10346 [BEDMatrix/extdata/mice.bed]

Convert BEDMatrix object to BGData object and read in phenotypes:

    > bgd <- as.BGData(bed, alternatePhenotypeFile = system.file("extdata", "mice.pheno", package = "BEDMatrix"))
    Extracting phenotypes from FAM file...
    Merging alternate phenotype file...

Display structure of BGData object:

    > str(bgd)
    Formal class "BGData" [package "BGData"] with 3 slots
    ..@ geno :BEDMatrix: 1814 x 10346 [/home/agrueneberg/.pkgs/R/BEDMatrix/extdata/mice.bed]
    ..@ pheno:"data.frame":       1814 obs. of  22 variables:
    .. ..$ FID                   : chr [1:1814] "A048005080" "A048006063" "A048006555" "A048007096" ...
    .. ..$ IID                   : chr [1:1814] "A048005080" "A048006063" "A048006555" "A048007096" ...
    .. ..$ PAT                   : int [1:1814] 0 0 0 0 0 0 0 0 0 0 ...
    .. ..$ MAT                   : int [1:1814] 0 0 0 0 0 0 0 0 0 0 ...
    .. ..$ SEX                   : int [1:1814] 2 1 1 1 2 1 1 1 1 2 ...
    .. ..$ PHENOTYPE             : int [1:1814] -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 ...
    .. ..$ PROJECT.NAME          : chr [1:1814] "HS_mouse" "HS_mouse" "HS_mouse" "HS_mouse" ...
    .. ..$ PHENOTYPE.NAME        : chr [1:1814] "Obesity" "Obesity" "Obesity" "Obesity" ...
    .. ..$ Obesity.BMI           : num [1:1814] -0.52 -0.401 -0.527 -0.415 -0.546 ...
    .. ..$ Obesity.BodyLength    : num [1:1814] 8.2 8.2 8.1 7.6 7.8 6.7 6 8.8 8.2 7.5 ...
    .. ..$ Date.Month            : int [1:1814] 5 3 4 3 1 3 3 4 3 2 ...
    .. ..$ Date.Year             : int [1:1814] 2003 2003 2003 2003 2003 2003 2003 2003 2003 2003 ...
    .. ..$ Date.Season           : chr [1:1814] "spring" "spring" "spring" "spring" ...
    .. ..$ Date.StudyStartSeconds: int [1:1814] 9072000 5443200 6652800 3628800 -518400 4838400 6048000 8467200 6048000 3024000 ...
    .. ..$ Date.Hour             : int [1:1814] 0 0 0 0 0 0 0 0 0 0 ...
    .. ..$ Date.StudyDay         : int [1:1814] 106 64 78 43 -5 57 71 99 71 36 ...
    .. ..$ GENDER                : chr [1:1814] "F" "M" "M" "M" ...
    .. ..$ EndNormalBW           : num [1:1814] 25.3 31.6 28.2 27.1 20.9 22.8 14.4 27.7 26.1 20.9 ...
    .. ..$ CoatColour            : chr [1:1814] "black" "dark.brown" "silverado" "chocolate" ...
    .. ..$ CageDensity           : int [1:1814] 5 5 4 4 6 7 3 4 4 3 ...
    .. ..$ Litter                : int [1:1814] 2 4 1 4 1 3 1 1 3 1 ...
    .. ..$ cage                  : chr [1:1814] "19F" "13C" "15C" "10B" ...
    ..@ map  :"data.frame":       10346 obs. of  1 variable:
    .. ..$ mrk: chr [1:10346] "rs3683945_G" "rs3707673_G" "rs6269442_G" "rs6336442_G" ...


### Creating a BGData object from a regular matrix

For small datasets, memory-mapping may not be necessary. In those cases, a BGData object can be created by manually providing `@geno`, `@pheno`, and `@map`.

`@geno` accepts matrix-like objects of various types (`LinkedMatrix`, `BEDMatrix`, `ff`, `big.matrix`, `matrix`), and optionally `@pheno` and `@map` accept data frames.

    bgd <- BGData(geno = genotypes, pheno = phenotypes, map = map)


### Creating a BGData object from a PLINK RAW file

This example uses a RAW version of the `mice` dataset that is included in the BGLR package. The dataset can be downloaded from https://github.com/QuantGen/BGData/raw/data/mice.raw.gz or generated using the script in the [mice.raw gist](https://gist.github.com/agrueneberg/f00e52c7d2cc2d666609e983df8ec3d7).

To load the BGData package:

    > library(BGData)

A `BGData` object can be generated from any plaintext file that stores individuals in rows, phenotypes in the first couple of columns (including an identifier for each individual), and genotypes as single allele dosage numbers in the remaining columns. This structure is intentionally similar to a [PED file](https://www.cog-genomics.org/plink2/formats#ped) that further restricts the structure of the phenotype section. The `mice` dataset that we will use as an example is only PED-like: there are more than six initial columns and the order of the columns is not according to the specification, but the BGData package is flexible enough to read it thanks to the `nColSkip` and `idCol` parameters.

    > bgd <- readPED(fileIn = "mice.raw.gz", header = TRUE, dataType = integer(), nColSkip = 17, idCol = 1)

Display structure of BGData object:

    > str(bgd)
    Formal class "BGData" [package "BGData"] with 3 slots
    ..@ geno :Formal class "RowLinkedMatrix" [package "LinkedMatrix"] with 1 slot
    .. .. ..@ .Data:List of 1
    .. .. .. ..$ : list()
    .. .. .. .. ..- attr(*, "physical")=Class "ff_pointer" <externalptr> 
    .. .. .. .. .. ..- attr(*, "vmode")= chr "byte"
    .. .. .. .. .. ..- attr(*, "maxlength")= int 18767644
    .. .. .. .. .. ..- attr(*, "pattern")= chr "ff"
    .. .. .. .. .. ..- attr(*, "filename")= chr "geno_1.bin"
    .. .. .. .. .. ..- attr(*, "pagesize")= int 65536
    .. .. .. .. .. ..- attr(*, "finalizer")= chr "close"
    .. .. .. .. .. ..- attr(*, "finonexit")= logi TRUE
    .. .. .. .. .. ..- attr(*, "readonly")= logi FALSE
    .. .. .. .. .. ..- attr(*, "caching")= chr "mmnoflush"
    .. .. .. .. ..- attr(*, "virtual")= list()
    .. .. .. .. .. ..- attr(*, "Length")= int 18767644
    .. .. .. .. .. ..- attr(*, "Dim")= int [1:2] 1814 10346
    .. .. .. .. .. ..- attr(*, "Dimorder")= int [1:2] 2 1
    .. .. .. .. .. ..- attr(*, "Symmetric")= logi FALSE
    .. .. .. .. .. ..- attr(*, "Dimnames")=List of 2
    .. .. .. .. .. .. ..$ : chr [1:1814] "A048005080_HS_mouse" "A048006063_HS_mouse" "A048006555_HS_mouse" "A048007096_HS_mouse" ...
    .. .. .. .. .. .. ..$ : chr [1:10346] "rs3683945_G" "rs3707673_G" "rs6269442_G" "rs6336442_G" ...
    .. .. .. .. .. - attr(*, "class") =  chr [1:3] "ff_matrix" "ff_array" "ff"
    ..@ pheno:"data.frame":       1814 obs. of  17 variables:
    .. ..$ SUBJECT.NAME          : chr [1:1814] "A048005080" "A048006063" "A048006555" "A048007096" ...
    .. ..$ PROJECT.NAME          : chr [1:1814] "HS_mouse" "HS_mouse" "HS_mouse" "HS_mouse" ...
    .. ..$ PHENOTYPE.NAME        : chr [1:1814] "Obesity" "Obesity" "Obesity" "Obesity" ...
    .. ..$ Obesity.BMI           : num [1:1814] -0.52 -0.401 -0.527 -0.415 -0.546 ...
    .. ..$ Obesity.BodyLength    : num [1:1814] 8.2 8.2 8.1 7.6 7.8 6.7 6 8.8 8.2 7.5 ...
    .. ..$ Date.Month            : int [1:1814] 5 3 4 3 1 3 3 4 3 2 ...
    .. ..$ Date.Year             : int [1:1814] 2003 2003 2003 2003 2003 2003 2003 2003 2003 2003 ...
    .. ..$ Date.Season           : chr [1:1814] "spring" "spring" "spring" "spring" ...
    .. ..$ Date.StudyStartSeconds: int [1:1814] 9072000 5443200 6652800 3628800 -518400 4838400 6048000 8467200 6048000 3024000 ...
    .. ..$ Date.Hour             : int [1:1814] 0 0 0 0 0 0 0 0 0 0 ...
    .. ..$ Date.StudyDay         : int [1:1814] 106 64 78 43 -5 57 71 99 71 36 ...
    .. ..$ GENDER                : chr [1:1814] "F" "M" "M" "M" ...
    .. ..$ EndNormalBW           : num [1:1814] 25.3 31.6 28.2 27.1 20.9 22.8 14.4 27.7 26.1 20.9 ...
    .. ..$ CoatColour            : chr [1:1814] "black" "dark.brown" "silverado" "chocolate" ...
    .. ..$ CageDensity           : int [1:1814] 5 5 4 4 6 7 3 4 4 3 ...
    .. ..$ Litter                : int [1:1814] 2 4 1 4 1 3 1 1 3 1 ...
    .. ..$ cage                  : chr [1:1814] "19F" "13C" "15C" "10B" ...
    ..@ map  :"data.frame":       10346 obs. of  1 variable:
    .. ..$ mrk: chr [1:10346] "mrk_1" "mrk_2" "mrk_3" "mrk_4" ...


The genotypes are internally stored as `ff` objects, part of the [ff package](http://cran.r-project.org/web/packages/ff/index.html). This package implements memory mapped arrays and provides a very fast implementation of indexing operations, which allows accessing cells of the array almost at the same speed as accessing those cells in a regular matrix object that is held in RAM. However, with `ff` the array size is limited to the size of an integer; with genomic data we often exceed this. We therefore developed a new package [`LinkedMatrix`](https://github.com/QuantGen/LinkedMatrix) and two new classes `RowLinkedMatrix` and `ColumnLinkedMatrix` that combine several matrix-like objects into a data structure that acts like a regular matrix. This package is used to overcome the limitations of `ff` by linking multiple `ff` objects together, either by columns (`ColumnLinkedMatrix`) or by rows (`RowLinkedMatrix`).
    
The files that back `ff` objects can be opened in other scientific environments such as Julia as well:

    fileIn = open("geno_1.bin", "r")
    X = mmap_array(Int8, (5, 10), fileIn)


### Saving a BGData object

A BGData object can be saved like any other R object using the `save` function:

    > save(bgd, file = "BGData.RData")


### Loading a BGData object

The genotypes in a `BGData` object can be of various types, some of which need to be initialized in a particular way. The `load.BGData` takes care of reloading a saved BGData object properly:

    > load.BGData("BGData.RData")
    Loaded objects: bgd


### Exploring operators

The `@geno` slot of a BGData object supports several matrix-like objects. A matrix-like object is an object that implements key methods such as subsetting, replacement, and others so that it looks and feels like a regular matrix in R even though its data may be stored on the filesystem and may be never read into memory in its entirety at a given time. This allows for convenient analysis of large datasets, with seamless integration into the rest of R's capabilities.

    # Subsetting
    bgd@geno[1, ]
    bgd@geno[, 1]
    bgd@geno[1:10, ]
    bgd@geno[, 1:10]
    bgd@geno[1, 1]

    # Replacement
    bgd@geno[1, 1] <- NA

    # Other methods
    dim(bgd@geno)
    nrow(bgd@geno)
    ncol(bgd@geno)
    rownames(bgd@geno)
    colnames(bgd@geno)
    dimnames(bgd@geno)


### Summarizing data

Use `chunkedApply` to count missing values (among others):

    countNAs <- chunkedApply(X = bgd, MARGIN = 2, FUN = function(x) sum(is.na(x)), bufferSize = 500)

Use the `summarize` function to calculate minor allele frequencies and frequency of missing values:

    summarize(bgd@geno)


### Running GWASes with different regression methods

A data structure for genomic data is useful when defining methods that act on both phenotype and genotype information. We have implemented a `GWAS` function that supports various regression methods. The formula takes phenotypes from `@pheno` and inserts one marker at a time.

    # lsfit (the default method)
    fmLM <- GWAS(formula = Obesity.BMI ~ GENDER + Litter, data = bgd)

    # lm
    fmLM <- GWAS(formula = Obesity.BMI ~ GENDER + Litter, data = bgd, method = "lm")

    # glm
    bgd@pheno$GENDER01 <- ifelse(bgd@pheno$GENDER == "M", 1, 0)
    fmGLM <- GWAS(formula = GENDER01 ~ Obesity.BMI, data = bgd, method = "glm", family = "binomial")

    # lmer
    fmLMER <- GWAS(formula = Obesity.BMI ~ GENDER + Litter + (1|cage), data = bgd, method = "lmer")

    # SKAT
    groups <- ceiling(1:ncol(bgd@geno) / 5)
    fmSKAT <- GWAS(formula = Obesity.BMI ~ GENDER + Litter, data = bgd, method = "SKAT", groups = groups)


### Generating the G Matrix

    G <- getG(bgd@geno)


Example Dataset
---------------

The example dataset in the `inst/extdata` folder was generated from the 250k SNP and phenotype data in [Atwell et al. (2010)](http://www.nature.com/nature/journal/v465/n7298/full/nature08800.html). Only the first 300 SNPs of chromosome 1, 2, and 3 were included to keep the size of the example dataset small. [PLINK](https://www.cog-genomics.org/plink2) was used to convert the data to BED and RAW files. `FT10` has been chosen as a phenotype and is provided as an [alternate phenotype file](https://www.cog-genomics.org/plink2/input#pheno). The file is intentionally shuffled to demonstrate that the additional phenotypes are put in the same order as the rest of the phenotypes.
