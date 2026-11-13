#!/usr/bin/env bash
# Chicken Productive Aging Atlas | 06_ASE_ASTS
# Purpose: WGS variant calling.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

set -euo pipefail
stage="${1:?Usage: bash wgs_variant_calling.sh STAGE ARGUMENTS; see module README}"
shift
case "${stage}" in
    trim)
        sample="${1:?Provide sample ID}"; r1="${2:?Provide R1 FASTQ}"; r2="${3:?Provide R2 FASTQ}"
        output_dir="${4:?Provide the output directory}"; mkdir -p "${output_dir}"
        fastp --thread=20 -i "${r1}" -I "${r2}" \
            -o "${output_dir}/${sample}_clean_1.fq.gz" -O "${output_dir}/${sample}_clean_2.fq.gz" \
            -j "${output_dir}/${sample}.json" -h "${output_dir}/${sample}.html"
        ;;
    align)
        sample="${1:?Provide sample ID}"; genome="${2:?Provide reference FASTA}"
        r1="${3:?Provide clean R1 FASTQ}"; r2="${4:?Provide clean R2 FASTQ}"; output="${5:?Provide output BAM}"
        bwa mem -t 20 -M -R "@RG\tID:${sample}\tLB:LIB_${sample}\tPL:ILLUMINA\tSM:${sample}" \
            "${genome}" "${r1}" "${r2}" | samtools sort -@ 20 -O bam -o "${output}"
        samtools index "${output}"
        ;;
    deduplicate)
        bam="${1:?Provide sorted BAM}"; output="${2:?Provide deduplicated BAM}"; metrics="${3:?Provide metrics file}"
        gatk MarkDuplicates -I "${bam}" -O "${output}" -M "${metrics}" --CREATE_INDEX true
        samtools flagstat "${output}" > "${output}.flagstat.txt"
        samtools depth -q 20 -Q 20 "${output}" \
            | perl -alne '{$pos{$F[0]}++;$depth{$F[0]}+=$F[2]} END{print "$_\t$pos{$_}\t$depth{$_}" foreach sort keys %pos}' \
            > "${output}.depth_summary.txt"
        ;;
    call)
        genome="${1:?Provide reference FASTA}"; bam="${2:?Provide deduplicated BAM}"; output="${3:?Provide output gVCF}"
        gatk HaplotypeCaller -R "${genome}" -I "${bam}" -O "${output}" \
            --min-base-quality-score 20 --minimum-mapping-quality 20 -ERC GVCF -native-pair-hmm-threads 20
        ;;
    combine)
        genome="${1:?Provide reference FASTA}"; manifest="${2:?Provide one gVCF path per line}"; output="${3:?Provide combined gVCF}"
        variants=()
        while IFS= read -r path || [[ -n "${path}" ]]; do
            [[ -z "${path}" ]] && continue
            variants+=(--variant "${path}")
        done < "${manifest}"
        ((${#variants[@]} > 0)) || { echo 'Empty gVCF manifest' >&2; exit 2; }
        gatk CombineGVCFs -R "${genome}" "${variants[@]}" -O "${output}"
        ;;
    genotype)
        genome="${1:?Provide reference FASTA}"; gvcf="${2:?Provide combined gVCF}"; output="${3:?Provide output VCF}"
        gatk GenotypeGVCFs -R "${genome}" --variant "${gvcf}" --include-non-variant-sites --output "${output}"
        ;;
    *) printf 'Unknown stage: %s\n' "${stage}" >&2; exit 2 ;;
esac
