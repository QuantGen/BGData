BGData: A Suite of Packages for Analysis of Big Genomic Data
============================================================

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/BGData)](https://CRAN.R-project.org/package=BGData)
[![Rdoc](http://www.rdocumentation.org/badges/version/BGData)](http://www.rdocumentation.org/packages/BGData)
[![Travis-CI Build Status](https://travis-ci.org/QuantGen/BGData.svg?branch=master)](https://travis-ci.org/QuantGen/BGData)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/8xac0jfmwhrj0gc3?svg=true)](https://ci.appveyor.com/project/agrueneberg/bgdata)
[![Coverage status](https://codecov.io/gh/QuantGen/BGData/branch/master/graph/badge.svg)](https://codecov.io/github/QuantGen/BGData?branch=master)

BGData is an R package that provides scalable and efficient computational methods for large genomic datasets, e.g., genome-wide association studies (GWAS) or genomic relationship matrices (G matrices). It also contains a data structure called `BGData` that holds genotypes in the `@geno` slot, phenotypes in the `@pheno` slot, and additional information in the `@map` slot.

Modern genomic datasets are big (large *n*), high-dimensional (large *p*), and multi-layered. The challenges that need to be addressed are memory requirements and computational demands. Our goal is to develop software that will enable researchers to carry out analyses with big genomic data within the R environment.

We have identified several approaches to tackle those challenges within R:

- File-backed matrices: The data is stored in on the hard drive and users can read in smaller chunks when they are needed.
- Linked arrays: For very large datasets a single file-backed array may not be enough or convenient. A linked array is an array whose content is distributed over multiple file-backed nodes.
- Multiple dispatch: Methods are presented to users so that they can treat these arrays pretty much as if they were RAM arrays.
- Multi-level parallelism: Exploit multi-core and multi-node computing.
- Inputs: Users can create these arrays from standard formats (e.g., PLINK .bed).

The BGData package is an umbrella package that comprises several packages: [BEDMatrix](https://CRAN.R-project.org/package=BEDMatrix), [LinkedMatrix](https://CRAN.R-project.org/package=LinkedMatrix), and [symDMatrix](https://CRAN.R-project.org/package=symDMatrix).


Examples
--------

### Loading the package

Load the BGData package:

```R
library(BGData)
```

### Inspecting the example dataset

The `inst/extdata` folder contains example files that were generated from the 250k SNP and phenotype data in [Atwell et al. (2010)](http://www.nature.com/nature/journal/v465/n7298/full/nature08800.html). Only the first 300 SNPs of chromosome 1, 2, and 3 were included to keep the size of the example dataset small enough for CRAN. [PLINK](https://www.cog-genomics.org/plink2) was used to convert the data to [.bed](https://www.cog-genomics.org/plink2/input#bed) and [.raw](https://www.cog-genomics.org/plink2/input#raw) files. `FT10` has been chosen as a phenotype and is provided as an [alternate phenotype file](https://www.cog-genomics.org/plink2/input#pheno). The file is intentionally shuffled to demonstrate that the additional phenotypes are put in the same order as the rest of the phenotypes.

```R
path <- system.file("extdata", package = "BGData")
list.files(path)
#>  [1] "chr1.bed"  "chr1.bim"  "chr1.fam"  "chr1.raw"  "chr2.bed"  "chr2.bim"
#>  [7] "chr2.fam"  "chr2.raw"  "chr3.bed"  "chr3.bim"  "chr3.fam"  "chr3.raw"
#> [13] "pheno.txt"
```

### Loading example dataset

#### Loading individual PLINK .bed files

Load the .bed file for chromosome 1 (chr1.bed) using the [BEDMatrix](https://CRAN.R-project.org/package=BEDMatrix) package:

```R
chr1 <- BEDMatrix(paste0(path, "/chr1.bed"))
#> Extracting number of individuals and rownames from .fam file...
#> Extracting number of markers and colnames from .bim file...
```

`BEDMatrix` objects behave similarly to regular matrices:

```R
dim(chr1)
#> [1] 199 300
rownames(chr1)[1:10]
#> [1] "5837_5837" "6008_6008" "6009_6009" "6016_6016" "6040_6040" "6042_6042"
#> [7] "6043_6043" "6046_6046" "6064_6064" "6074_6074"
colnames(chr1)[1:10]
#> [1] "snp1_T"  "snp2_G"  "snp3_A"  "snp4_T"  "snp5_G"  "snp6_T"  "snp7_C"
#> [8] "snp8_C"  "snp9_C"  "snp10_G"
chr1["6008_6008", "snp5_G"]
#> [1] 0
```

#### Linking multiple BEDMatrix objects together

Load the other two .bed files:

```R
chr2 <- BEDMatrix(paste0(path, "/chr2.bed"))
#> Extracting number of individuals and rownames from .fam file...
#> Extracting number of markers and colnames from .bim file...
chr3 <- BEDMatrix(paste0(path, "/chr3.bed"))
#> Extracting number of individuals and rownames from .fam file...
#> Extracting number of markers and colnames from .bim file...
```

Combine the BEDMatrix objects by columns using the [LinkedMatrix](https://CRAN.R-project.org/package=LinkedMatrix) to avoid the inconvenience of having three separate matrices:

```R
wg <- ColumnLinkedMatrix(chr1, chr2, chr3)
```

Just like `BEDMatrix` objects, `LinkedMatrix` objects also behave similarly to regular matrices:

```R
dim(wg)
#> [1] 199 900
rownames(wg)[1:10]
#> [1] "5837_5837" "6008_6008" "6009_6009" "6016_6016" "6040_6040" "6042_6042"
#> [7] "6043_6043" "6046_6046" "6064_6064" "6074_6074"
colnames(wg)[1:10]
#> [1] "snp1_T"  "snp2_G"  "snp3_A"  "snp4_T"  "snp5_G"  "snp6_T"  "snp7_C"
#> [8] "snp8_C"  "snp9_C"  "snp10_G"
wg["6008_6008", "snp5_G"]
#> [1] 0
```

### Creating a BGData object

`BGData` objects can be created from individual `BEDMatrix` objects or a collection of `BEDMatrix` objects as a `LinkedMatrix` object using the `as.BGData()` function. This will read the .fam and .bim file that comes with the .bed files. The `alternatePhenotypeFile` parameter points to the file that contains the `FT10` phenotype:

```R
bg <- as.BGData(wg, alternatePhenotypeFile = paste0(path, "/pheno.txt"))
#> Extracting phenotypes from .fam file, assuming that the .fam file of the first BEDMatrix instance is representative of all the other nodes...
#> Extracting map from .bim files...
#> Merging alternate phenotype file...
```

The `bg` object will have the `LinkedMatrix` object in the `@geno` slot, the .fam file augmented by the `FT10` phenotype in the `@pheno` slot, and the .bim file in the `@map` slot.

```R
str(bg)
#> Formal class 'BGData' [package "BGData"] with 3 slots
#>   ..@ geno :Formal class 'ColumnLinkedMatrix' [package "LinkedMatrix"] with 1 slot
#>   .. .. ..@ .Data:List of 3
#>   .. .. .. ..$ :BEDMatrix: 199 x 300 [/home/agrueneberg/.pkgs/R/BGData/extdata/chr1.bed]
#>   .. .. .. ..$ :BEDMatrix: 199 x 300 [/home/agrueneberg/.pkgs/R/BGData/extdata/chr2.bed]
#>   .. .. .. ..$ :BEDMatrix: 199 x 300 [/home/agrueneberg/.pkgs/R/BGData/extdata/chr3.bed]
#>   ..@ pheno:'data.frame':       199 obs. of  7 variables:
#>   .. ..$ FID      : int [1:199] 5837 6008 6009 6016 6040 6042 6043 6046 6064 6074 ...
#>   .. ..$ IID      : int [1:199] 5837 6008 6009 6016 6040 6042 6043 6046 6064 6074 ...
#>   .. ..$ PAT      : int [1:199] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. ..$ MAT      : int [1:199] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. ..$ SEX      : int [1:199] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. ..$ PHENOTYPE: int [1:199] -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 ...
#>   .. ..$ FT10     : num [1:199] 57 60 98 75 71 56 90 93 96 91 ...
#>   ..@ map  :'data.frame':       900 obs. of  6 variables:
#>   .. ..$ chromosome        : int [1:900] 1 1 1 1 1 1 1 1 1 1 ...
#>   .. ..$ snp_id            : chr [1:900] "snp1" "snp2" "snp3" "snp4" ...
#>   .. ..$ genetic_distance  : int [1:900] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. ..$ base_pair_position: int [1:900] 657 3102 4648 4880 5975 6063 6449 6514 6603 6768 ...
#>   .. ..$ allele_1          : chr [1:900] "T" "G" "A" "T" ...
#>   .. ..$ allele_2          : chr [1:900] "C" "A" "C" "C" ...
```

### Saving a BGData object

A BGData object can be saved like any other R object using the `save` function:

```R
save(bg, file = "BGData.RData")
```

### Loading a BGData object

The genotypes in a `BGData` object can be of various types, some of which need to be initialized in a particular way. The `load.BGData` takes care of reloading a saved BGData object properly:

```R
load.BGData("BGData.RData")
#> Loaded objects: bg
```

### Summarizing data

Use `chunkedApply` to count missing values (among others):

```R
countNAs <- chunkedApply(X = bg@geno, MARGIN = 2, FUN = function(x) sum(is.na(x)))
```

Use the `summarize` function to calculate minor allele frequencies and frequency of missing values:

```R
summarize(bg@geno)
```

### Running GWASes with different regression methods

A data structure for genomic data is useful when defining methods that act on both phenotype and genotype information. We have implemented a `GWAS` function that supports various regression methods. The formula takes phenotypes from `@pheno` and inserts one marker at a time.

```R
gwas <- GWAS(formula = FT10 ~ 1, data = bg)
```

### Generating the G Matrix

```R
G <- getG(bg@geno)
```


Installation
------------

Install the stable version from CRAN:

```R
install.packages("BGData")
```

Alternatively, install the development version from GitHub:

```R
# install.packages("remotes")
remotes::install_github("QuantGen/BGData")
```


Contribute
----------

- Issue Tracker: https://github.com/QuantGen/BGData/issues
- Source Code: https://github.com/QuantGen/BGData


Documentation
-------------

Further documentation can be found on [RDocumentation](http://www.rdocumentation.org/packages/BGData).
