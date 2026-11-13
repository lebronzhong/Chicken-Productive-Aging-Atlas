#!/usr/bin/env Rscript
# Test expression-defined transcript sets using the recorded enrichment statistic.
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file)), "run_gwas_enrichment.R"))
run_gwas_enrichment(commandArgs(trailingOnly = TRUE))
