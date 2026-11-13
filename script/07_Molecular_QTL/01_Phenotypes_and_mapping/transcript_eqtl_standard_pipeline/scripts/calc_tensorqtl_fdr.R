#!/usr/bin/env Rscript
# Chicken Productive Aging Atlas | 07_Molecular_QTL
# Purpose: Calc tensorqtl fdr.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript calc_tensorqtl_fdr.R <cis_qtl.txt.gz> <output_txt> <fdr>")
}

file_perm <- args[1]
file_output <- args[2]
fdr_threshold <- as.numeric(args[3])

if (!file.exists(file_perm)) {
  stop(paste("Input file does not exist:", file_perm))
}

eqtl <- read.table(gzfile(file_perm), sep = "\t", header = TRUE, check.names = FALSE)
cat("* Number of molecular phenotypes =", nrow(eqtl), "\n")
cat("* Correlation between beta-approximation and empirical p-values =",
    round(cor(eqtl$pval_beta, eqtl$pval_perm), 4), "\n")

qval <- p.adjust(eqtl$pval_beta, method = "BH")
set0 <- eqtl[qval <= fdr_threshold, ]
set1 <- eqtl[qval > fdr_threshold, ]

if (nrow(set0) == 0) {
  stop("No significant phenotypes were detected at the requested FDR threshold.")
}
if (nrow(set1) == 0) {
  stop("All phenotypes are significant; cannot derive a separating empirical threshold robustly.")
}

pthreshold <- (sort(set1$pval_perm)[1] - sort(-1.0 * set0$pval_perm)[1]) / 2
cat("* Corrected empirical p-value threshold =", pthreshold, "\n")

nthresholds <- qbeta(
  pthreshold,
  eqtl$beta_shape1,
  eqtl$beta_shape2,
  ncp = 0,
  lower.tail = TRUE,
  log.p = FALSE
)

eqtl$qval <- qval
eqtl$pval_nominal_threshold <- nthresholds
eqtl$is_eGene <- (eqtl$pval_nominal < eqtl$pval_nominal_threshold) & (eqtl$qval <= fdr_threshold)

write.table(eqtl, file_output, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Done\n")
