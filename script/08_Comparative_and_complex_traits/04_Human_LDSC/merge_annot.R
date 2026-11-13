#!/usr/bin/env Rscript
# Chicken Productive Aging Atlas | 08_Comparative_and_complex_traits
# Purpose: Merge annot.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

library(data.table)

# Path configuration.
baseline_dir <- "/path/to/chicken-aging/LDSC/ready_data/baseline_by_chicken_aging"
annot_dir <- "/path/to/chicken-aging/LDSC/24.3.annot"
merged_dir <- "/path/to/chicken-aging/LDSC/24.35.merged_annot" 
tissue_list <- scan("/path/to/chicken-aging/LDSC/shell/tissue_list", what = "", sep = "\n")

for(chr in 1:22) {
  # Read the baseline annotation.
  baseline_file <- file.path(baseline_dir, paste0("baseline.", chr, ".annot.gz"))
  baseline <- fread(cmd = paste("zcat", baseline_file), header = TRUE)
  colnames(baseline) <- "BASELINE"

  # Read and merge tissue annotations.
  for(tissue in tissue_list) {
    tissue_file <- file.path(annot_dir, paste0(tissue, ".", chr, ".annot.gz"))  
    annot <- fread(cmd = paste("zcat", tissue_file), header = TRUE)
    colnames(annot) <- toupper(tissue)  # Use uppercase column names.
    baseline <- cbind(baseline, annot)
  }

  # Write the merged compressed annotation.
  out_file <- file.path(merged_dir, paste0("merged.", chr, ".annot.gz"))        
  fwrite(baseline, out_file, sep = "\t", quote = FALSE, col.names = TRUE)       

  cat("chr", chr, "merged:", ncol(baseline), "columns, ", nrow(baseline), "rows\n")
}
