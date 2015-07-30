BGData
======

# Memory Mapped Matrices and Data-Structures for Genomic Data for R



## Description
Genetic data can be very large and holding data in RAM is often not feasible. One approach to overcome this restriction are [memory mapped files](http://en.wikipedia.org/wiki/Memory-mapped_file) that store the actual data on the hard drive and only read in smaller chunks when they are needed.

The [ff package for R](http://cran.r-project.org/web/packages/ff/index.html) implements memory mapped arrays and provides a very fast implementation of indexing operations, which allows accessing cells of the array almost at the same speed as accessing those cells in a regular matrix object that is held in RAM. However, with `ff` the array size is limited to the size of an integer; with genomic data we often exceed this.



We are therefore developing new classes (`rmmMatrix` and `cmmMatrix`) which are essentially collections of `ff` objects. In these classes we distribute a matrix either by rows (`rmmMatrix`) or columns (`cmmMatrix`) into multiple `ff` objects. We have developed indexing and many other methods that allow the user to deal with these objects as if they were regular matrices. In addition we have developed methods that can take `rmmMatrix` or `cmmMatrix` as input to compute genomic relationship matrices, etc.

The classes `cmmMatrix` and `rmmMatrix` were designed to hold genotype data. The class `BGData` contains three slots `@geno`, `@pheno` and `@map`, and can be used to hold GWAS data.

![Conceptual Diagram](https://docs.google.com/drawings/d/1m2bV3-woWrO9F9_RXxw30FlzUORXzcnaIPJUlxc-MMk/pub?w=739&h=559)

## Developers
- Gustavo de los Campos (gdeloscampos@gmail.com)
- Alexander Grüneberg (alexander.grueneberg@googlemail.com)
- Paulino Perez (perpdgo@gmail.com)

## Classes & Methods
- `rmmMatrix`: row-distributed matrix (collection of `ff` objects)
- `cmmMatrix`: column-distributed matrix (collection of `ff` objects)
- `mmMatrix` : distributed matrix (a class-union of `rmmMatrix`  and `cmmMatrix`)
- `BGData`: a structure to hold genotype and phenotype data. Objects of this class contain three slots: `@geno` (`rmmMatrix`, `cmmMatrix`, `ff_matrix` or `matrix`), `@pheno` (data.frame), and `@map` (data.frame).

### Methods Implemented for `rmmMatrix` and `cmmMatrix`
- `[` and `[<-` for subsetting and replacement, respectively
- `dim(x)`, `nrow(x)`, `ncol(x) `  
- `rownames(x)`, `colnames(x)`, `dimnames(x)`
- `colSums(x)`, `colMeans(x)`, `rowSums(x)`, `rowMeans(x)`
- `summary(x)`
- `apply(x)` with the same parameters and similar behavior than the generic function apply()
- `chunks(x)` returns information about how the `rmmMatrix` or `cmmMatrix` is split into chunks (each chunk is an `ff` object)
- `colindexes(x, columns)` returns the global (in the `cmmMatrix` object) and local (in each of the `ff` objects that constitute the chunks of the `cmmMatrix`) indexes for a set of columns
- `rowindexes(x, columns)` returns the global (in the `rmmMatrix` object) and local (in each of the `ff` objects that constitute the chunks of the `rmmMatrix`) indexes for a set of rows
- `getG(x, ...)` computes a genomic relationship matrix (XX') with options for centering and scaling
- `as.matrix(x)` converts an `mmMatrix` to a matrix (if small enough)

### Methods Implemented for `BGData`
- Both `readPED` and `readPED.matrix` create a `BGData` object from a plaintext file containing the phenotypes and genotypes (individuals in rows, phenotypes in the first few columns, markers in the remaining columns, e.g. the raw format in [PLINK](http://pngu.mgh.harvard.edu/~purcell/plink/dataman.shtml)). `readPED` stores genotype information in an `rmmMatrix` or `cmmMatrix` (dependending on the value of the `distributed.by` parameter) while `readPED.matrix` uses a regular matrix.
- `GWAS` uses a `BGData` object to conduct single marker association tests using regression methods such as `lm()`, `glm()` or `lmer()`


## Installation

BGData is not available on [CRAN](http://cran.r-project.org/) yet. However, it can be installed directly from GitHub using the [devtools](https://github.com/hadley/devtools) package.

1. Install `devtools` package: `install.packages('devtools')`
2. Load `devtools` package: `library(devtools)`
3. Install `BGData` package from GitHub: `install_github('QuantGen/BGData')`
4. Load `BGData` package: `library(BGData)`

Alternatively, you can download the most recent version as a bundle for [Windows](https://github.com/QuantGen/BGData/archive/master.zip) or [Mac OS / Linux](https://github.com/QuantGen/BGData/archive/master.tar.gz).


## Introduction

### Prerequisites

The `BGData` package has to be installed to follow along. The example data used in this introduction was generated from the `mice` dataset of the `BGLR` package and serves as an example on how to create a `BGData` object from a plaintext file. We have uploaded it for convenience and it can be downloaded at https://github.com/QuantGen/BGData/raw/data/mice.ped.gz, but you can also generate it in the following way:

```R
library(BGLR)
data(mice)
write.table(cbind(mice.pheno, mice.X), 'mice.ped', quote=FALSE,
            row.names=FALSE)
```

### Converting a plaintext file (e.g. in PED format) to a `BGData` object
A `BGData` object can be generated from any plaintext file that stores individuals in rows, phenotypes in the first couple of columns (including an identifier for each individual), and genotypes either coded as characters (`A`, `C`, `G`, `T`) or numbers (`1` for `A`, `2` for `C`, etc.) in the remaining columns. Genotypes will be scanned into `BGData` object exactly as it is coded in the text file as long as there is only one column for each marker. This structure is similar to a [PED file](http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped) that further restricts the structure of the phenotype section. The `mice` dataset that we will use as an example is only PED-like: there are more than six initial columns and the order of the columns is not according to the specification, but the BGData package is flexible enough to read it thanks to the `nColSkip` and `idCol` parameters.

A `BGData` object has three slots: `@pheno`, `@geno`, and `@map`. The phenotypes go into `@pheno`, the genotypes into `@geno`, and `@map` is filled with a placeholder that can be replaced later on.

```R
BGData <- readPED(fileIn='mice.ped.gz', header=TRUE,
                  dataType=integer(), nColSkip=17, idCol=1)
head(BGData@pheno)
dim(BGData@geno)
dim(BGData@map)
```

### Reloading a `BGData` object from the filesystem
The genotypes in a `BGData` object are backed by an efficient binary representation of the original dataset on the filesystem, an `mmMatrix`. By default `readPED` stores this representation as files called `geno_*.bin` in the current working directory in a folder that starts with `BGData_` followed by the filename without extension. To reload a `BGData` object from the filesystem, load the accompanying `BGData.RData` file in that directory. There is one caveat, though: you have to change your current working directory to the one that contains the file for it to work.

```R
rm(BGData)
# Note: The working directory must be the one where the binary files 
#       (geno_*.bin) are saved.
setwd('BGData_mice.ped/')
load('BGData.RData')
head(BGData@pheno)
dim(BGData@geno)
dim(BGData@map)
```

### Exploring operators
An `mmMatrix` object (the datatype of the contents of the `@geno` slot of a `BGData` object) behaves like any other matrix, even though its data is stored on the filesystem and is never read into memory in its entirety at a given time. This allows for convenient analysis of large datasets, with seamless integration into the rest of R's capabilities. Subsetting, replacement, and other functions have been implemented.

```R
# Subsetting
BGData@geno[1,]
BGData@geno[,1]
BGData@geno[1:10,]
BGData@geno[,1:10]
BGData@geno[1,1]

# Replacement
BGData@geno[1,1] <- NA

# Other generic
dim(BGData@geno)
nrow(BGData@geno)
ncol(BGData@geno)
rownames(BGData@geno)
colnames(BGData@geno)
dimnames(BGData@geno)

# Summaries
rowSums(BGData@geno)
colSums(BGData@geno)
rowMeans(BGData@geno)
colMeans(BGData@geno)

# apply
countNAs <- apply(X=BGData@geno, MARGIN=2, FUN=function(x) sum(is.na(x)))
```

### Running GWASes with different regression methods
A data structure for genomic data is useful when defining methods that act on both phenotype and genotype information. We have implemented a `GWAS` function that supports various regression methods and plotting. The formula takes phenotypes from `@pheno` and inserts one marker at a time.

```R
# lm (set plot=TRUE to get a Manhattan plot of the p values)
fmLM <- GWAS(formula=Obesity.BMI~GENDER+Litter, data=BGData, method='lm')

# glm (set plot=TRUE to get a Manhattan plot of the p values)
BGData@pheno$GENDER01 <- ifelse(BGData@pheno[,'GENDER'] == 'M', 1, 0)
fmGLM <- GWAS(formula=GENDER01~Obesity.BMI, data=BGData, method='glm',
              family='binomial')

# lmer (set plot=TRUE to get a Manhattan plot of the p values)
fmLMER <- GWAS(formula=Obesity.BMI~GENDER+Litter+(1|cage), 
               data=BGData, method='lmer')

# SKAT (set plot=TRUE to get a Manhattan plot of the p values) 
groups <- ceiling(1:ncol(genData@geno) / 5)
fmSKAT <- GWAS(formula=Obesity.BMI~GENDER+Litter, 
               data=genData, method='SKAT', groups=groups)
```

### Adding columns to phenotype data
`@pheno` is stored as a regular data frame and can be modified. These changes will also be available in functions such as `GWAS`.

```R
newPheno <- data.frame(SUBJECT.NAME=BGData@pheno$SUBJECT.NAME,
                       NEWCOLUMN=1:nrow(BGData@pheno))
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
