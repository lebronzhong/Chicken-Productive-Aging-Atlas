#!/usr/bin/env Rscript
# Chicken Productive Aging Atlas | 02_Expression_dynamics
# Purpose: Within tissue zscore.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("Usage: Rscript SCRIPT MATRIX_TSV SAMPLE_MANIFEST_TSV TISSUE OUTPUT_TSV")
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file)), "read_age_samples.R"))
aligned <- read_age_samples(args[1], args[2], args[3])
output_file <- args[4]
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
count_matrix <- aligned$matrix
database <- data.frame(condition = factor(paste0("wk", aligned$samples$age)), row.names = aligned$samples$sample)
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = database,
                              design = ~ condition)
dds <- DESeq(dds)
normalized_counts <- counts(dds, normalized=TRUE)
gene_means <- rowMeans(normalized_counts)
gene_sd <- apply(normalized_counts, 1, sd)
z_scores <- (normalized_counts - gene_means) / gene_sd

write.table(z_scores, output_file,quote = F,col.names = T,row.names = T,sep = "\t")
