# BGData dev

- Check if `pheno` and `map` are data.frames when creating a BGData object.


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
