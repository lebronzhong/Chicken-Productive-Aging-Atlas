#!/usr/bin/env Rscript
# Chicken Productive Aging Atlas | 07_Molecular_QTL
# Purpose: Usage to qtl bed.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript usage_to_qtl_bed.R <psi_file> <gtf_file> <out_file>")
}

psi_file <- args[1]
gtf_file <- args[2]
out_file <- args[3]

# -----------------------------
# Input and output files.
# -----------------------------
#psi_file <- "/path/to/chicken-aging/qtl/inputs/Tissue_transcript_TPM/salmon_pheno.ChickenGTEx.Duodenum.TPM_isoform.filtered.psi"
#gtf_file <- "/path/to/chicken-aging/qtl/inputs/Chicken_GTEx_input/union.topo.gtf"
#out_file <- "/path/to/chicken-aging/qtl/inputs/Tissue_transcript_TPM/salmon_pheno.ChickenGTEx.Duodenum.TPM_isoform.filtered.tensorqtl.bed.gz"

# -----------------------------
# 1. Read the PSI matrix.
# The header contains sample names but no label for the first column.
# -----------------------------
cat("Reading PSI matrix...\n")

psi <- fread(
  psi_file,
  header = TRUE,
  sep = "\t",
  na.strings = c("nan", "NA", "")
)

# Extract transcript_id.
psi[, transcript_id := sub("^.*;", "", gene_tx)]

cat("PSI rows:", nrow(psi), "\n")

# -----------------------------
# 2. Read transcript features from the GTF.
# Extract transcript and gene IDs, chromosome, start and end.
# -----------------------------
cat("Reading transcript coordinates from GTF...\n")

gtf <- fread(
  gtf_file,
  sep = "\t",
  header = FALSE,
  quote = "",
  fill = TRUE
)

if (ncol(gtf) < 9) {
  stop("GTF file does not have at least 9 columns.")
}

setnames(gtf, c("chr","source","feature","start","end","score","strand","frame","attribute"))

gtf_tx <- gtf[feature == "transcript", .(chr, start, end, attribute)]

# Extract transcript_id.
gtf_tx[, transcript_id := sub('.*transcript_id "([^"]+)".*', '\\1', attribute)]

# Remove duplicates.
gtf_tx <- unique(gtf_tx[, .(transcript_id, chr, start, end)])

cat("Transcript entries in GTF:", nrow(gtf_tx), "\n")

# -----------------------------
# 3. Join PSI values to transcript coordinates.
# -----------------------------
merged <- merge(
  gtf_tx,
  psi,
  by = "transcript_id",
  all = FALSE
)

cat("Merged PSI rows:", nrow(merged), "\n")

# -----------------------------
# 4. Replace missing abundance values with zero.
# Operate on sample columns only.
# -----------------------------
meta_cols <- c("transcript_id", "chr", "start", "end", "gene_tx")
sample_cols <- setdiff(colnames(merged), meta_cols)

for (cc in sample_cols) {
  set(merged, i = which(is.na(merged[[cc]])), j = cc, value = 0)
}

# -----------------------------
# 5. Arrange columns in TensorQTL phenotype BED format.
# #chr start end phenotype_id sample1 sample2 ...
# -----------------------------
bed <- merged[, c("chr", "start", "end", "transcript_id", sample_cols), with = FALSE]
setnames(bed, c("chr", "start", "end", "transcript_id"), c("#chr", "start", "end", "phenotype_id"))

# Sort records.
bed <- bed[order(`#chr`, start, end, phenotype_id)]

# -----------------------------
# 6. Write output.
# -----------------------------
tmp_out <- sub("\\.gz$", "", out_file)
fwrite(bed, file = tmp_out, sep = "\t", quote = FALSE, na = "0")
system(paste("gzip -f", shQuote(tmp_out)))

cat("Written TensorQTL BED:", out_file, "\n")
cat("Done.\n")

