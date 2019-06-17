\name{multi-level-parallelism}
\alias{multi-level-parallelism}
\title{Multi-Level Parallelism}
\description{
    Functions with the \code{nCores}, \code{i}, and \code{j} parameters provide
    capabilities for both parallel and distributed computing.

    For parallel computing, \code{nCores} determines the number of cores the
    code is run on. Memory usage can be an issue for higher values of
    \code{nCores} as R is not particularly memory-efficient. As a rule of
    thumb, at least around \code{(nCores * object_size(chunk)) +
    object_size(result)} MB of total memory will be needed for operations
    on file-backed matrices, not including potential copies of your data that
    might be created (for example \code{lsfit} runs \code{cbind(1, X)}).
    \code{i} and \code{j} can be used to include or exclude certain rows or
    columns. Internally, the \code{mclapply} function is used and therefore
    parallel computing will not work on Windows machines.

    For distributed computing, \code{i} and \code{j} determine the subset of
    the input matrix that the code runs on. In an HPC environment, this can be
    used not just to include or exclude certain rows or columns, but also to
    partition the task among many nodes rather than cores. Scheduler-specific
    code and code to aggregate the results need to be written by the user. It
    is recommended to set \code{nCores} to \code{1} as nodes are often cheaper
    than cores.
}
\seealso{
    \code{\link[parallel]{mclapply}} to learn more about the function used to
    implement parallel computing. \code{\link[parallel]{detectCores}} to detect
    the number of available cores.
}
