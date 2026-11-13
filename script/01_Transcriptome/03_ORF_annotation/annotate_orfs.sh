#!/usr/bin/env bash
# Chicken Productive Aging Atlas | 01_Transcriptome
# Purpose: ORF annotation.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

set -euo pipefail
stage="${1:?Usage: bash annotate_orfs.sh STAGE ARGUMENTS; see module README}"
shift
td_bin="${TRANSDECODER_BIN:?Set TRANSDECODER_BIN to the TransDecoder installation}"
case "${stage}" in
    sequence)
        gtf="${1:?Provide the study GTF}"; genome="${2:?Provide the reference FASTA}"
        prefix="${3:?Provide an output prefix}"
        "${td_bin}/util/gtf_genome_to_cdna_fasta.pl" "${gtf}" "${genome}" > "${prefix}.fasta"
        "${td_bin}/util/gtf_to_alignment_gff3.pl" "${gtf}" > "${prefix}.gff3"
        ;;
    long_orfs)
        cdna="${1:?Provide cDNA FASTA}"; output_dir="${2:?Provide output directory}"
        "${td_bin}/TransDecoder.LongOrfs" -t "${cdna}" --output_dir "${output_dir}"
        ;;
    homology)
        peptides="${1:?Provide longest_orfs.pep}"; database="${2:?Provide Swiss-Prot DIAMOND database}"
        output="${3:?Provide the output table}"
        diamond blastp --query "${peptides}" --db "${database}" --outfmt 6 --out "${output}" \
            --max-target-seqs 1 --evalue 1e-5 --threads 20 --more-sensitive
        ;;
    domains)
        peptides="${1:?Provide longest_orfs.pep}"; database="${2:?Provide Pfam HMM database}"
        output="${3:?Provide the output table}"
        hmmsearch --cpu 24 -E 1e-10 --domtblout "${output}" "${database}" "${peptides}"
        ;;
    predict)
        cdna="${1:?Provide cDNA FASTA}"; pfam="${2:?Provide Pfam hits}"
        blastp="${3:?Provide Swiss-Prot hits}"; output_dir="${4:?Provide the prediction directory}"
        "${td_bin}/TransDecoder.LongOrfs" -t "${cdna}" --output_dir "${output_dir}"
        "${td_bin}/TransDecoder.Predict" -t "${cdna}" --retain_pfam_hits "${pfam}" \
            --retain_blastp_hits "${blastp}" --output_dir "${output_dir}"
        ;;
    project)
        orfs="${1:?Provide predicted ORF GFF3}"; alignments="${2:?Provide transcript alignment GFF3}"
        cdna="${3:?Provide cDNA FASTA}"; prefix="${4:?Provide an output prefix}"
        "${td_bin}/util/cdna_alignment_orf_to_genome_orf.pl" "${orfs}" "${alignments}" "${cdna}" > "${prefix}.gff3"
        gffread "${prefix}.gff3" -T -o "${prefix}.gtf"
        ;;
    *) printf 'Unknown stage: %s\n' "${stage}" >&2; exit 2 ;;
esac
