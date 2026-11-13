#!/usr/bin/env bash
# Chicken Productive Aging Atlas | 01_Transcriptome
# Purpose: Read processing.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

set -euo pipefail
sample_dir="${1:?Usage: bash process_reads.sh SAMPLE_DIRECTORY GENOME_FASTA}"
genome="${2:?Provide the reference FASTA}"
mkdir -p "${sample_dir}/2_pychopper" "${sample_dir}/3_minimap2" "${sample_dir}/4_samtools"
pychopper -t 20 -Q 7 -z 50 \
    -r "${sample_dir}/2_pychopper/report.pdf" \
    -u "${sample_dir}/2_pychopper/unclassified.fq" \
    -w "${sample_dir}/2_pychopper/rescued.fq" \
    "${sample_dir}/pass.fq.gz" "${sample_dir}/2_pychopper/pass_full_length.fq"
minimap2 -ax splice -uf -k14 "${genome}" \
    "${sample_dir}/2_pychopper/pass_full_length.fq" > "${sample_dir}/3_minimap2/aln.sam"
samtools sort -@ 20 -l 4 -O bam -T "${sample_dir}/4_samtools/tem_prefix" \
    -o "${sample_dir}/4_samtools/aln_sort.bam" "${sample_dir}/3_minimap2/aln.sam"
samtools index "${sample_dir}/4_samtools/aln_sort.bam"
samtools flagstat "${sample_dir}/4_samtools/aln_sort.bam" > "${sample_dir}/4_samtools/flagstat.txt"
