# BGData 1.0.0.9000

- Change buffering strategy to improve parallelism: instead of buffering
  `bufferSize` in the main process, buffer `bufferSize` in the forks. That way
  `nTasks` is not necessary anymore and the buffering code is the same for one
  core and for multiple cores.
- Remove `nTasks` parameter from `chunkedApply()` and methods based on it.
- getG_symDMatrix: Compute block indices correctly for out-of-order,
  non-sequential indices.
- as.BGData: Use rownames of BEDMatrix object as rownames for pheno, and
  colnames of BEDMatrix object as rownames for map.
- orderedMerge: Retain rownames.


# BGData 1.0.0

Initial release.
