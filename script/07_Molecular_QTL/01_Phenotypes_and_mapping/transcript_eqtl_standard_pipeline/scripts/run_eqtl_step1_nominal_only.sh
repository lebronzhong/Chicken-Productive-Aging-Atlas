#!/usr/bin/env bash
# Chicken Productive Aging Atlas | 07_Molecular_QTL
# Purpose: Run eqtl step1 nominal only.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

set -euo pipefail

usage() {
  echo "Usage: bash run_eqtl_step1_nominal_only.sh -g <genotype_prefix> -p <phenotype_bed_gz> -c <covariates_txt> -t <name> -o <outdir>"
  exit 1
}

while getopts ":g:p:c:t:o:" opt; do
  case "$opt" in
    g) genotype="$OPTARG" ;;
    p) phenotype="$OPTARG" ;;
    c) covariates="$OPTARG" ;;
    t) tissue="$OPTARG" ;;
    o) outdir="$OPTARG" ;;
    *) usage ;;
  esac
done

[[ -z "${genotype:-}" || -z "${phenotype:-}" || -z "${covariates:-}" || -z "${tissue:-}" || -z "${outdir:-}" ]] && usage

for ext in bed bim fam; do
  [[ -f "${genotype}.${ext}" ]] || { echo "Missing genotype file: ${genotype}.${ext}" >&2; exit 1; }
done
[[ -f "$phenotype" ]] || { echo "Missing phenotype file: $phenotype" >&2; exit 1; }
[[ -f "$covariates" ]] || { echo "Missing covariates file: $covariates" >&2; exit 1; }

mkdir -p "$outdir"
prefix="${outdir}/ChickenAging.${tissue}"

python3 -m tensorqtl \
  "$genotype" "$phenotype" "$prefix" \
  --mode cis_nominal \
  --covariates "$covariates" \
  --maf_threshold 0.05


