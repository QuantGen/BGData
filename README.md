BGData: Memory Mapped Matrices and Data-Structures for Genomic Data for R
=========================================================================

[![Travis-CI Build Status](https://travis-ci.org/QuantGen/BGData.svg?branch=master)](https://travis-ci.org/QuantGen/BGData)

**Contact**:     gruenebe@msu.edu   gustavoc@msu.edu

Genetic data can be very large and holding data in RAM is often not feasible. One approach to overcome this restriction are [memory mapped files](http://en.wikipedia.org/wiki/Memory-mapped_file) that store the actual data on the hard drive and only read in smaller chunks when they are needed.

The [ff package for R](http://cran.r-project.org/web/packages/ff/index.html) implements memory mapped arrays and provides a very fast implementation of indexing operations, which allows accessing cells of the array almost at the same speed as accessing those cells in a regular matrix object that is held in RAM. However, with `ff` the array size is limited to the size of an integer; with genomic data we often exceed this.

We therefore developed a new package [`LinkedMatrix`](https://github.com/QuantGen/LinkedMatrix) and two new classes `RowLinkedMatrix` and `ColumnLinkedMatrix` that combine several matrix-like objects into a data structure that acts like a regular matrix. This package is used to overcome the limitations of `ff` by linking multiple `ff` objects together, either by columns (`ColumnLinkedMatrix`) or by rows (`RowLinkedMatrix`).

The `BGData` package contains a data structure `BGData` that holds genotypes in the `@geno` slot, phenotypes in the `@pheno` slot, and additional information in the `@map` slot. In addition, we have developed several methods that operate on this data structure, for example to compute genomic relationship matrices, etc.

![Conceptual Diagram](https://docs.google.com/drawings/d/1m2bV3-woWrO9F9_RXxw30FlzUORXzcnaIPJUlxc-MMk/pub?w=739&h=559)


Installation
------------
The BGData package is not available on [CRAN](http://cran.r-project.org/) yet, but it can be installed directly from GitHub using the [devtools](https://github.com/hadley/devtools) package:

```R
# install.packages("devtools")
devtools::install_github("QuantGen/BGData")
```

Alternatively, you can download the most recent version as a bundle for [Windows](https://github.com/QuantGen/BGData/archive/master.zip) or [Mac OS / Linux](https://github.com/QuantGen/BGData/archive/master.tar.gz).


Introduction
------------

### Prerequisites

The `BGData` package has to be installed to follow along. The example data used in this introduction was generated from the `mice` dataset of the `BGLR` package and serves as an example on how to create a `BGData` object from a plaintext file. We have uploaded it for convenience and it can be downloaded at https://github.com/QuantGen/BGData/raw/data/mice.raw.gz, but you can also generate it in the following way:

```R
# Writes genotypes in an ASCII file
library(BGLR)
data(mice)
gzmice <- gzfile('mice.raw.gz', 'w')
write.table(cbind(mice.pheno, mice.X), gzmice, quote=F, row.names=F)
close(gzmice)
```

### Converting a plaintext file (e.g. in PED format) to a `BGData` object
A `BGData` object can be generated from any plaintext file that stores individuals in rows, phenotypes in the first couple of columns (including an identifier for each individual), and genotypes---either coded as characters (e.g. `A`, `C`, `G`, `T`) or numbers (e.g. `1` for `A`, `2` for `C`, etc.), but preferably as single allele dosage numbers---in the remaining columns. This structure is intentionally similar to a [PED file](http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped) that further restricts the structure of the phenotype section. The `mice` dataset that we will use as an example is only PED-like: there are more than six initial columns and the order of the columns is not according to the specification, but the BGData package is flexible enough to read it thanks to the `nColSkip` and `idCol` parameters.

A `BGData` object has three slots: `@pheno`, `@geno`, and `@map`. The phenotypes go into `@pheno`, the genotypes into `@geno`, and `@map` is filled with a placeholder that can be replaced later on.

```R
library(BGData)
BGData <- readPED(fileIn='mice.raw.gz', header=TRUE, dataType=integer(), nColSkip=17, idCol=1)
str(BGData)
```

### Reloading a `BGData` object from the filesystem
The genotypes in a `BGData` object are backed by an efficient binary representation of the original dataset on the filesystem, a `LinkedMatrix`. By default `readPED` stores this representation as files called `geno_*.bin` in the current working directory in a folder that starts with `BGData_` followed by the filename without extension. To reload a `BGData` object from the filesystem, load the accompanying `BGData.RData` file in that directory using the `load.BGData` function.

```R
rm(BGData)
load.BGData('BGData_mice.ped/BGData.RData')
dim(BGData@geno)
```

`BGData` objects can also be loaded using the regular `load` function, but you have to change your current working directory to the one that contains the file for it to work.

```R
rm(BGData)
# The working directory must be the one where the binary files (geno_*.bin) are saved.
setwd('BGData_mice.ped/')
load('BGData.RData')
dim(BGData@geno)
```

### Exploring operators
A `LinkedMatrix` object (the datatype of the contents of the `@geno` slot of a `BGData` object) behaves like any other matrix, even though its data is stored on the filesystem and is never read into memory in its entirety at a given time. This allows for convenient analysis of large datasets, with seamless integration into the rest of R's capabilities. Subsetting, replacement, and other functions have been implemented.

```R
# Subsetting
BGData@geno[1, ]
BGData@geno[, 1]
BGData@geno[1:10, ]
BGData@geno[, 1:10]
BGData@geno[1, 1]

# Replacement
BGData@geno[1, 1] <- NA

# Other methods
dim(BGData@geno)
nrow(BGData@geno)
ncol(BGData@geno)
rownames(BGData@geno)
colnames(BGData@geno)
dimnames(BGData@geno)

# Summaries
summarize(BGData@geno)

# apply
countNAs <- chunkedApply(X=BGData@geno, MARGIN=2, FUN=function(x) sum(is.na(x)), bufferSize = 500)
```

### Running GWASes with different regression methods
A data structure for genomic data is useful when defining methods that act on both phenotype and genotype information. We have implemented a `GWAS` function that supports various regression methods and plotting. The formula takes phenotypes from `@pheno` and inserts one marker at a time.

```R
# lm (set plot=TRUE to get a Manhattan plot of the p values)
fmLM <- GWAS(formula=Obesity.BMI~GENDER+Litter, data=BGData, method='lm', plot=T)

# glm (set plot=TRUE to get a Manhattan plot of the p values)
BGData@pheno$GENDER01 <- ifelse(BGData@pheno[, 'GENDER'] == 'M', 1, 0)
fmGLM <- GWAS(formula=GENDER01~Obesity.BMI, data=BGData, method='glm', family='binomial')

# lmer (set plot=TRUE to get a Manhattan plot of the p values)
fmLMER <- GWAS(formula=Obesity.BMI~GENDER+Litter+(1|cage), data=BGData, method='lmer')

# SKAT (set plot=TRUE to get a Manhattan plot of the p values)
groups <- ceiling(1:ncol(BGData@geno) / 5)
fmSKAT <- GWAS(formula=Obesity.BMI~GENDER+Litter, data=BGData, method='SKAT', groups=groups)
```

### Adding columns to phenotype data
`@pheno` is stored as a regular data frame and can be modified. These changes will also be available in functions such as `GWAS`.

```R
newPheno <- data.frame(SUBJECT.NAME=BGData@pheno$SUBJECT.NAME, NEWCOLUMN=1:nrow(BGData@pheno))
BGData@pheno <- merge(BGData@pheno, newPheno, by='SUBJECT.NAME', sort=FALSE)
```

### Generating the G Matrix
```R
G <- getG(BGData@geno)
```

### BGData can be used with genotypes in RAM (e.g., using matrix)
```R
n <- 50
p <- 100
X <- matrix(nrow=n, ncol=p, sample(0:2, size=n*p, replace=TRUE))
colnames(X) <- paste0('mrk_', 1:p)
rownames(X) <- paste0('id_', 1:n)
Y <- data.frame(id=rownames(X), age=rnorm(n), y=rnorm(n))
MAP <- data.frame(mrk_id=colnames(X))

# Creating a BGData object from objects in RAM
BGData <- BGData(geno=X, pheno=Y, map=MAP)

# Computing a genomic relationship matrix
G <- getG(BGData@geno)

# Single-marker regression
SMR <- GWAS(y~age, data=BGData, method='lm', plot=TRUE)
```

### Reading the `ff` objects in [Julia](http://julialang.org/) using `mmap_array`
```julia
fileIn = open("geno_1.bin", "r")
X = mmap_array(Int8, (5, 10), fileIn)
```
