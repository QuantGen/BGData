# BGData 1.0.0.9000

## Breaking Changes

- Rename `bufferSize` to `chunkSize`.
- Remove `nTasks` parameter from `chunkedApply()` and methods based on it (see
  above).

## Other Changes

- Change chunking strategy to improve parallelism: instead of loading a subset
  of `chunkSize` in the main process, load a subset of `chunkSize` in the each
  fork. That way `nTasks` is not necessary anymore and the same code can be
  used for one core and multiple cores.
- Add `findRelated()` function for use with matrices and symDMatrix objects.
- getG_symDMatrix: Add chunkSize parameter.
- getG_symDMatrix: Compute block indices correctly for out-of-order,
  non-sequential indices.
- as.BGData: Use rownames of BEDMatrix object as rownames for pheno, and
  colnames of BEDMatrix object as rownames for map.
- orderedMerge: Retain rownames.
- Include process ID in verbose output if nCores > 1.


# BGData 1.0.0

Initial release.
