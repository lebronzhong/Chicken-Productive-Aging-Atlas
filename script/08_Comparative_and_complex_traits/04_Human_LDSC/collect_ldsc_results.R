#!/usr/bin/env Rscript
# Chicken Productive Aging Atlas | 08_Comparative_and_complex_traits
# Purpose: Collect ldsc results.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

library(dplyr)
library(readr)
library(stringr)
library(purrr)

# Paths.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: Rscript collect_ldsc_results.R TRAIT_LIST RESULTS_DIR OUTPUT_PREFIX")
merge_file <- args[1]
results_dir <- args[2]
output_prefix <- args[3]
dir.create(dirname(output_prefix), recursive=TRUE, showWarnings=FALSE)

# Read GWAS IDs.
gwas_ids <- read_lines(merge_file)

# Read one GWAS result file.
read_ldsc <- function(gwas_id) {
  file <- file.path(results_dir, paste0(gwas_id, ".results"))
  if (!file.exists(file)) {
    warning(paste("File not found:", file))
    return(NULL)
  }
  df <- read_tsv(file, col_types = cols())
  # The Category column identifies the tissue annotation.
  df <- df %>% select(Category, Enrichment_p)
  colnames(df)[2] <- gwas_id   # Name the Enrichment_p column after the GWAS ID.
  return(df)
}

# Read and merge all results.
all_data <- map(gwas_ids, read_ldsc) %>%
  compact() %>%                # Remove NULL results.
  reduce(full_join, by = "Category")

# Save results.
out_file <- paste0(output_prefix, ".tsv")
write.table(all_data,out_file,quote=F,row.names=F,col.names=T)

out_file <- paste0(output_prefix, ".xlsx")
library(openxlsx)
write.xlsx(all_data, paste0(output_prefix, ".xlsx"))

cat("Merged results written to:", out_file, "\n")
