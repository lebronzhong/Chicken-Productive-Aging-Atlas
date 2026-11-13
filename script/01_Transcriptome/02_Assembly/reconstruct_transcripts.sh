#!/usr/bin/env bash
# Chicken Productive Aging Atlas | 01_Transcriptome
# Purpose: Transcript reconstruction and annotation merge.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

set -euo pipefail
bam="${1:?Usage: bash reconstruct_transcripts.sh BAM GENOME_FASTA REFERENCE_GTF OUTPUT_DIR}"
genome="${2:?Provide the reference FASTA}"
annotation="${3:?Provide the reference annotation}"
output_dir="${4:?Provide the output directory}"
isoquant.py --reference "${genome}" --genedb "${annotation}" --complete_genedb \
    --bam "${bam}" --data_type nanopore -o "${output_dir}" --threads 20
