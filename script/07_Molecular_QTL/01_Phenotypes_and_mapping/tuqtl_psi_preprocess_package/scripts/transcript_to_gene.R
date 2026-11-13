#!/usr/bin/env Rscript
# Chicken Productive Aging Atlas | 07_Molecular_QTL
# Purpose: Transcript to gene.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

suppressPackageStartupMessages({
  library(data.table)
  library(tximport)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript transcript_to_gene.R <fam> <salmon_dir> <tx2gene_csv> <out_prefix> <countsFromAbundance>")
}

fam_file <- args[1]
salmon_dir <- args[2]
tx2gene_file <- args[3]
out_prefix <- args[4]
counts_from_abundance <- args[5]   # no / scaledTPM / lengthScaledTPM

cat("fam_file: ", fam_file, "\n")
cat("salmon_dir: ", salmon_dir, "\n")
cat("tx2gene_file: ", tx2gene_file, "\n")
cat("out_prefix: ", out_prefix, "\n")
cat("countsFromAbundance: ", counts_from_abundance, "\n")

# 1. Read sample IDs from the second FAM column (IID).
fam <- fread(fam_file, header = FALSE)
if (ncol(fam) < 2) stop("FAM file must have at least 2 columns.")
samples <- fam[[2]]

cat("Samples in FAM:", length(samples), "\n")
cat("Unique samples   :", length(unique(samples)), "\n")

# 2. Construct quant.sf paths.
files <- file.path(salmon_dir, samples, "quant.sf")
names(files) <- samples

# 3. Check for missing input files.
exists_flag <- file.exists(files)
if (!all(exists_flag)) {
  missing_df <- data.frame(
    sample = names(files)[!exists_flag],
    file = files[!exists_flag],
    stringsAsFactors = FALSE
  )
  fwrite(missing_df, paste0(out_prefix, ".missing_quantsf.tsv"), sep = "\t")
  stop(sprintf("Missing %d quant.sf files. See %s.missing_quantsf.tsv",
               sum(!exists_flag), out_prefix))
}

cat("All quant.sf files found.\n")

# 4. Read the transcript-to-gene mapping.
tx2gene <- fread(tx2gene_file)
req_cols <- c("TXNAME", "GENEID")
if (!all(req_cols %in% colnames(tx2gene))) {
  stop("tx2gene file must contain columns: TXNAME, GENEID")
}
tx2gene <- tx2gene[, .(TXNAME, GENEID)]

cat("tx2gene rows:", nrow(tx2gene), "\n")
cat("unique transcripts:", uniqueN(tx2gene$TXNAME), "\n")
cat("unique genes      :", uniqueN(tx2gene$GENEID), "\n")

# 5. tximport
txi <- tximport(
  files = files,
  type = "salmon",
  tx2gene = tx2gene,
  countsFromAbundance = counts_from_abundance,
  ignoreTxVersion = FALSE,
  ignoreAfterBar = FALSE
)

# 6. Export results.
counts_dt <- as.data.table(txi$counts, keep.rownames = "GENEID")
abundance_dt <- as.data.table(txi$abundance, keep.rownames = "GENEID")
length_dt <- as.data.table(txi$length, keep.rownames = "GENEID")

fwrite(counts_dt, paste0(out_prefix, ".gene_counts.tsv.gz"), sep = "\t")
fwrite(abundance_dt, paste0(out_prefix, ".gene_tpm.tsv.gz"), sep = "\t")
fwrite(length_dt, paste0(out_prefix, ".gene_length.tsv.gz"), sep = "\t")

cat("Output written:\n")
cat("  ", paste0(out_prefix, ".gene_counts.tsv.gz"), "\n")
cat("  ", paste0(out_prefix, ".gene_tpm.tsv.gz"), "\n")
cat("  ", paste0(out_prefix, ".gene_length.tsv.gz"), "\n")
cat("Done.\n")

