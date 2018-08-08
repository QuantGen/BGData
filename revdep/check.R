#!/usr/bin/env Rscript

library(devtools)

revdep_check(bioconductor = TRUE)
revdep_check_save_summary()
revdep_check_print_problems()
