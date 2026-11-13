#!/usr/bin/env bash
# Chicken Productive Aging Atlas | 07_Molecular_QTL
# Purpose: iso-eQTL mapping.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

set -euo pipefail

tissue="${1:?Usage: bash run_transcript_eqtl_standard.sh TISSUE [WORK_DIR]}"
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
work_dir="${2:-.}"
cd -- "${work_dir}"

GENO="ChickenGTEx.${tissue}"
RAW_BED="salmon_pheno.ChickenGTEx.${tissue}.TPM.bed.gz"
COV="${tissue}.combined_covariates.txt"

OUTROOT="${tissue}_standard_qtl"
PREP_DIR="${OUTROOT}/phenotype"
CHR_DIR="${OUTROOT}/phenotype_by_chr"
RES_DIR="${OUTROOT}/results"
PREP_BED="${PREP_DIR}/salmon_pheno.ChickenGTEx.${tissue}.TPM.clean.norm.bed.gz"

mkdir -p "$PREP_DIR" "$CHR_DIR" "$RES_DIR"

echo "[1] Prepare clean normalized phenotype BED"
python3 "${script_dir}/prepare_transcript_qtl_bed.py" \
  --input-bed "$RAW_BED" \
  --bim "${GENO}.bim" \
  --output-bed-gz "$PREP_BED" \
  --tpm-threshold 0.1 \
  --min-sample-frac 0.2

echo "[2] Check samples"
python3 "${script_dir}/check_inputs.py" \
  --genotype "$GENO" \
  --phenotype "$PREP_BED" \
  --covariates "$COV"

echo "[3] Split phenotype by chromosome"
python3 "${script_dir}/split_bed_by_chr.py" \
  --input-bed-gz "$PREP_BED" \
  --outdir "$CHR_DIR" \
  --prefix "salmon_pheno.ChickenGTEx.${tissue}.TPM.clean.norm"

echo "[4] Run cis_nominal by chromosome"
for bed in "${CHR_DIR}"/*.bed.gz; do
  chr=$(basename "$bed" | sed -E 's/.*\.chr([^.]*)\.bed\.gz/\1/')
  echo "  Running chromosome ${chr}"

  bash "${script_dir}/run_eqtl_step1_nominal_only.sh" \
    -g "$GENO" \
    -p "$bed" \
    -c "$COV" \
    -t "${tissue}_chr${chr}" \
    -o "$RES_DIR/chr${chr}"
done

echo "[5] Run cis permutation and FDR on full clean phenotype"
bash "${script_dir}/run_eqtl_step1.sh" \
  -g "$GENO" \
  -p "$PREP_BED" \
  -c "$COV" \
  -t "$tissue" \
  -o "$RES_DIR/full"

echo "[done] Results are in ${OUTROOT}"


