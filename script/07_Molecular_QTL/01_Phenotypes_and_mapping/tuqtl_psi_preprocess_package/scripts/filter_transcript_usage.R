#!/usr/bin/env Rscript
# Chicken Productive Aging Atlas | 07_Molecular_QTL
# Purpose: Filter transcript usage.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7) {
  stop("Usage: Rscript filter_transcript_usage.R <gene_tpm_file> <tx2gene_file> <tx_bed_file> <psi_file> <gene_tpm_out> <tx_bed_out> <psi_out>")
}

gene_tpm_file <- args[1]
tx2gene_file <- args[2]
tx_bed_file <- args[3]
psi_file <- args[4]
gene_tpm_out <- args[5]
tx_bed_out <- args[6]
psi_out <- args[7]

#gene_tpm_file <- "/path/to/chicken-aging/qtl/inputs/Gene_expression/ChickenGTEx.Duodenum.gene_tpm.tsv.gz"
#tx2gene_file  <- "/path/to/chicken-aging/qtl/inputs/shell/tx2gene_union.csv"
#tx_bed_file   <- "/path/to/chicken-aging/qtl/inputs/Tissue_transcript_TPM/salmon_pheno.ChickenGTEx.Duodenum.TPM.bed.gz"
#psi_file      <- "/path/to/chicken-aging/qtl/inputs/Tissue_transcript_TPM/salmon_pheno.ChickenGTEx.Duodenum.TPM_isoform.psi"

#gene_tpm_out  <- "/path/to/chicken-aging/qtl/inputs/Gene_expression/ChickenGTEx.Duodenum.gene_tpm.filtered.tsv.gz"
#tx_bed_out    <- "/path/to/chicken-aging/qtl/inputs/Tissue_transcript_TPM/salmon_pheno.ChickenGTEx.Duodenum.TPM.filtered.bed.gz"
#psi_out       <- "/path/to/chicken-aging/qtl/inputs/Tissue_transcript_TPM/salmon_pheno.ChickenGTEx.Duodenum.TPM_isoform.filtered.psi"

# -----------------------------
# 2. Read gene TPM and filter by median abundance.
# -----------------------------
cat("Reading gene TPM matrix...\n")
gene_tpm <- fread(cmd = paste("zcat", shQuote(gene_tpm_file)))

if (!"GENEID" %in% colnames(gene_tpm)) {
  stop("GENEID column not found in gene TPM matrix.")
}

sample_cols_gene <- setdiff(colnames(gene_tpm), "GENEID")

gene_tpm[, median_tpm := apply(.SD, 1, median, na.rm = TRUE), .SDcols = sample_cols_gene]
gene_keep <- gene_tpm[median_tpm > 0.1, GENEID]

cat("Total genes           :", nrow(gene_tpm), "\n")
cat("Genes kept (median>0.1):", length(gene_keep), "\n")

gene_tpm_filtered <- gene_tpm[GENEID %in% gene_keep, c("GENEID", sample_cols_gene), with = FALSE]

fwrite(gene_tpm_filtered, file = gene_tpm_out, sep = "\t", quote = FALSE)
cat("Written filtered gene TPM:", gene_tpm_out, "\n")

# -----------------------------
# 3. Map retained genes to transcripts using tx2gene.
# -----------------------------
cat("Reading tx2gene...\n")
tx2gene <- fread(tx2gene_file)

if (!all(c("TXNAME", "GENEID") %in% colnames(tx2gene))) {
  stop("tx2gene file must contain TXNAME and GENEID columns.")
}

tx_keep <- unique(tx2gene[GENEID %in% gene_keep, TXNAME])

cat("Total tx2gene rows    :", nrow(tx2gene), "\n")
cat("Transcripts kept      :", length(tx_keep), "\n")

# -----------------------------
# 4. Filter the transcript TPM BED.
# -----------------------------
cat("Reading transcript TPM BED...\n")
tx_bed <- fread(cmd = paste("zcat", shQuote(tx_bed_file)))

if (!"phenotype_id" %in% colnames(tx_bed)) {
  stop("phenotype_id column not found in transcript BED.")
}

tx_bed_filtered <- tx_bed[phenotype_id %in% tx_keep]

fwrite(tx_bed_filtered, file = tx_bed_out, sep = "\t", quote = FALSE)
cat("Written filtered transcript TPM BED:", tx_bed_out, "\n")
cat("Transcript BED kept rows:", nrow(tx_bed_filtered), "\n")

# -----------------------------
# 5. Filter the PSI matrix.
# The first PSI column has no header and contains GENEID;TRANSCRIPT identifiers.
# -----------------------------
cat("Reading PSI matrix...\n")

psi_header <- scan(psi_file, what = character(), sep = "\t", nlines = 1, quiet = TRUE)

psi_colnames <- c("gene_tx", psi_header)

psi <- fread(
  psi_file,
  header = FALSE,
  skip = 1,
  sep = "\t",
  na.strings = c("nan", "NA", "")
)

if (ncol(psi) != length(psi_colnames)) {
  stop(sprintf(
    "PSI file column mismatch: data has %d columns, header expects %d columns",
    ncol(psi), length(psi_colnames)
  ))
}

# 5. Assign column names.
setnames(psi, psi_colnames)

# 6. Extract the transcript ID after the semicolon.
psi[, transcript_id := sub("^.*;", "", gene_tx)]

psi_filtered <- psi[transcript_id %in% tx_keep]
psi_filtered[, transcript_id := NULL]

fwrite(psi_filtered, file = psi_out, sep = "\t", quote = FALSE, na = "nan")
cat("Written filtered PSI matrix:", psi_out, "\n")
cat("PSI kept rows:", nrow(psi_filtered), "\n")
cat("Done.\n")

