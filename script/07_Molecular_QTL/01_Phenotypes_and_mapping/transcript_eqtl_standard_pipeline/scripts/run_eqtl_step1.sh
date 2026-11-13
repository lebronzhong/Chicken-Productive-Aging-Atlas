#!/usr/bin/env bash
# Chicken Productive Aging Atlas | 07_Molecular_QTL
# Purpose: Run eqtl step1.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

set -euo pipefail

usage() {
  cat <<USAGE
Usage:
  bash run_eqtl_step1.sh \
    -g <genotype_prefix> \
    -p <phenotype_bed_gz> \
    -c <covariates_txt> \
    -t <tissue_name> \
    -o <outdir>

Required inputs:
  -g  PLINK genotype prefix. Example: /path/to/genotype/ChickenGTEx.{tissue}
      The following files must exist: .bed .bim .fam
  -p  Molecular phenotype BED.gz file used by TensorQTL
  -c  Covariates file (samples in columns, covariates in rows, first column ID)
  -t  Tissue or analysis name used in output prefix
  -o  Output directory
USAGE
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

log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
log "Running cis nominal mapping"
python3 -m tensorqtl "$genotype" "$phenotype" "$prefix" \
  --mode cis_nominal --covariates "$covariates" --maf_threshold 0.05
log "Running cis permutation mapping"
python3 -m tensorqtl "$genotype" "$phenotype" "$prefix" \
  --mode cis --covariates "$covariates" --maf_threshold 0.05
log "Calculating FDR and phenotype-specific nominal thresholds"
Rscript "$SCRIPT_DIR/calc_tensorqtl_fdr.R" "${prefix}.cis_qtl.txt.gz" "${prefix}.cis_qtl.txt" 0.05
awk 'BEGIN{FS=OFS="\t"} NR==1 || $0 ~ /\tTRUE$/ {print}' \
  "${prefix}.cis_qtl.txt" > "${prefix}.cis_qtl.fdr.txt"
log "Step 1 completed"
