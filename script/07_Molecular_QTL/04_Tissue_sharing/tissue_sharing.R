#!/usr/bin/env Rscript
# Chicken Productive Aging Atlas | 07_Molecular_QTL
# Purpose: Tissue sharing.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(svglite)
  library(ragg)
})

# Fig. 7e-f: cross-tissue clustering of cis-QTL effects.
# Effect correlations are calculated only for identical phenotype_id + variant_id
# pairs observed in both tissues. This avoids correlating effects of different SNPs.

root_dir <- "/path/to/chicken-aging/analysis/QTL_summary"
out_dir <- file.path(root_dir, "Fig7_e_f")
source_dir <- file.path(out_dir, "source_data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript tissue_sharing.R TISSUE_METADATA_TSV")
metadata <- read.delim(args[[1]], stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(all(c("tissue", "system", "color_hex", "source_label", "system_color") %in% names(metadata)))
stopifnot(nrow(metadata) >= 2L, !anyDuplicated(metadata$tissue), !anyNA(metadata))
tissues <- metadata$tissue
systems <- setNames(metadata$system, tissues)

system_metadata <- unique(metadata[, c("system", "system_color")])
if (anyDuplicated(system_metadata$system)) stop("Each system requires one consistent color")
system_cols <- setNames(system_metadata$system_color, system_metadata$system)

tissue_cols <- setNames(metadata$color_hex, tissues)
fwrite(metadata[, c("tissue", "color_hex", "source_label")],
       file.path(source_dir, "Fig7e_tissue_colour_mapping.tsv"), sep = "\t")

file_map <- list(
  eQTL = setNames(file.path("/path/to/chicken-aging/qtl/inputs", paste0("ChickenGTEx.", tissues, ".cis_qtl.txt")), tissues),
  teQTL = setNames(file.path(root_dir, "teQTL", paste0("ChickenAging.", tissues, ".cis_qtl.txt")), tissues),
  tuQTL = setNames(file.path(root_dir, "tuQTL", paste0("ChickenAging.", tissues, "_PSI.cis_qtl.txt")), tissues)
)

missing_files <- unlist(file_map)[!file.exists(unlist(file_map))]
if (length(missing_files) > 0) stop("Missing input files: ", paste(missing_files, collapse = "; "))

read_effects <- function(paths) {
  out <- lapply(names(paths), function(tissue) {
    x <- fread(paths[[tissue]], select = c("phenotype_id", "variant_id", "slope"))
    x <- x[is.finite(slope) & !is.na(phenotype_id) & !is.na(variant_id)]
    x[, pair_id := paste(phenotype_id, variant_id, sep = "||")]
    x <- x[, .(slope = mean(slope)), by = pair_id]
    setkey(x, pair_id)
    x
  })
  names(out) <- names(paths)
  out
}

pairwise_correlation <- function(effect_list, type) {
  n <- length(effect_list)
  rho <- matrix(NA_real_, n, n, dimnames = list(tissues, tissues))
  n_pairs <- matrix(NA_integer_, n, n, dimnames = list(tissues, tissues))
  diag(rho) <- 1
  for (i in seq_len(n)) {
    n_pairs[i, i] <- nrow(effect_list[[i]])
    if (i < n) for (j in (i + 1):n) {
      z <- merge(
        effect_list[[i]], effect_list[[j]], by = "pair_id",
        suffixes = c("_i", "_j"), all = FALSE
      )
      r <- if (nrow(z) >= 30) {
        suppressWarnings(cor(z$slope_i, z$slope_j, method = "spearman", use = "complete.obs"))
      } else {
        NA_real_
      }
      rho[i, j] <- rho[j, i] <- r
      n_pairs[i, j] <- n_pairs[j, i] <- nrow(z)
    }
  }
  if (anyNA(rho)) stop(type, ": insufficient identical phenotype-variant pairs for one or more tissue pairs")
  list(rho = rho, n = n_pairs)
}

effects <- lapply(file_map, read_effects)
stats <- lapply(names(effects), function(type) pairwise_correlation(effects[[type]], type))
names(stats) <- names(effects)

write_matrix <- function(mat, filename, row_name = "tissue") {
  x <- data.frame(setNames(list(rownames(mat)), row_name), mat, check.names = FALSE)
  fwrite(x, file.path(source_dir, filename), sep = "\t", na = "NA")
}

for (type in names(stats)) {
  write_matrix(stats[[type]]$rho, paste0(type, "_tissue_Spearman_rho.tsv"))
  write_matrix(stats[[type]]$n, paste0(type, "_shared_phenotype_variant_pair_counts.tsv"))
}

# Complete-linkage clustering of 1 - Spearman rho.
hclusts <- lapply(stats, function(x) hclust(as.dist(1 - x$rho), method = "complete"))
