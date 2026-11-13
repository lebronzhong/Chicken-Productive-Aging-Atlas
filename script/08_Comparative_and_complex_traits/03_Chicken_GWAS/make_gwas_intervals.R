#!/usr/bin/env Rscript
# Chicken Productive Aging Atlas | 08_Comparative_and_complex_traits
# Purpose: Make gwas intervals.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

# Set input and output files.

infile <- "/path/to/chicken-aging/ont/shell/10_20240903.gffcompare/9.6_final_gtf/9.9_finalgtf_annotated_rename.gtf"

outfile <- "/path/to/chicken-aging/ont/shell/23.LD_score_regression/Transcript_LDSC.20k.txt"

# Read the GTF.

gtf <- read.table(infile, sep="\t", header=FALSE, stringsAsFactors=FALSE, quote="")

# Standard GTF columns.

colnames(gtf) <- c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

transcripts <- subset(gtf, feature == "transcript")

library(stringr)

transcripts$transcript_id <- str_match(transcripts$attribute, 'transcript_id "([^"]+)"')[,2]

transcripts$new_start <- pmax(1, transcripts$start - 20000)

transcripts$new_end <- transcripts$end + 20000

result <- transcripts[, c("chr", "new_start", "new_end", "transcript_id")]

write.table(result, file=outfile, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

