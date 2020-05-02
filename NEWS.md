# BGData 2.1.0.9000

- Follow [Bioconductor S4 practices][1].
  - If you have used `new()` to create `BGData` instances, please use the
    `BGData()` constructor function instead.
  - If you have used `@` to access the slots of `BGData` instances, please use
    the `geno()`, `pheno()`, and `map()` accessors instead.
- `BGData()`:
  - Do not create dimnames for `geno` as this object is likely shared.
  - Check if `geno` has row names before creating `pheno` stub.
  - Check if `geno` has column names before creating `map` stub.
  - Rename `IID` in `pheno` stub to `sample_id`.
  - Rename `mrk` in `map` stub to `variant_id`.
  - Change format of rownames for `pheno` stub to a sequence starting with
    `sample_` and rownames for `map` stub to a sequence starting with
    `variant_` if `geno` does not have dimnames.
- `as.BGData()`:
  - Force column classes when loading .fam and .bim files.
  - Force `FID` and `IID` columns to be of type `character` when loading
    alternate phenotype files.
  - Do not make assumptions about the structure of dimnames of a BEDMatrix
    object if it is passed without .fam and .bim file unless they are `NULL`.
- Add validity tests for `BGData` objects:
  - Check if number of rows of `geno` matches number of rows of `pheno`.
  - Check if number of columns of `geno` matches number of rows of `map`.
  - Warn if the row names of `pheno` do not match the row names of `geno`.
  - Warn if the row names of `map` do not match the column names of `geno`.
- Add `preprocess()` function for fast centering, scaling, and imputation.
- `GWAS()`: Return number of records used for each variant and allele
  frequencies in `rayOLS`.
- Update citation instructions.


# BGData 2.1.0

- Add `chunkedMap()` function.
- Improve error handling in `chunkedMap()` and `chunkedApply()`.
- `summarize()`: Improve performance.
- `GWAS()`: Improve performance of `rayOLS` method.
- `GWAS()`: Fix bug when computing p-values for methods other than rayOLS,
  lsfit, or SKAT when `i` is used to subset samples.
- `GWAS()`: Fix wrong results in `lsfit` method when covariates with missing
  values are used.
- `as.BGData()`: Fix bug loading .fam and .bim files when path contains the
  word `bed`.


# BGData 2.0.0

## Breaking Changes

- Rename `bufferSize` to `chunkSize`.
- Remove `nTasks` parameter from `chunkedApply()` and methods based on it.
- Remove `crossprods` function.

## Other Changes

- Change chunking strategy to improve parallelism: instead of loading a subset
  of `chunkSize` in the main process, load a subset of `chunkSize` in the each
  fork. That way `nTasks` is not necessary anymore and the same code can be
  used for one core and multiple cores.
- Add `findRelated()` function for use with matrices and symDMatrix objects.
- Add `orderedMerge()` function that allows for phenotypes to be easily merged
  into a BGData object.
- Performance improvements in `getG()` function: use single shared memory
  matrix to collect results.
- Performance improvements in `rayOLS` method in `GWAS()` function.
- `getG_symDMatrix()`: Support version 2 of symDMatrix package.
- `getG_symDMatrix()`: Add `chunkSize` parameter.
- `getG_symDMatrix()`: Add `minVar` parameter.
- `as.BGData()`: Use rownames of BEDMatrix object as rownames for pheno, and
  colnames of BEDMatrix object as rownames for map.
- Include process ID in verbose output if `nCores` > 1.

## Bug Fixes

- `getG_symDMatrix()`: Fix scaling error when `scale = FALSE`.
- `getG_symDMatrix()`: Compute block indices correctly for out-of-order,
  non-sequential indices.
- `getG_symDMatrix()`: Do not include centers and scales in attributes anymore
  because the influence of `j` and `minVar` is difficult to retain.


# BGData 1.0.0

Initial release.

[1]: https://bioconductor.org/help/course-materials/2017/Zurich/S4-classes-and-methods.html
