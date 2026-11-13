#!/usr/bin/env bash
# Chicken Productive Aging Atlas | 07_Molecular_QTL
# Purpose: Qc qtl inputs.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

set -euo pipefail

GENO=${1:?genotype prefix required}
BED=${2:?phenotype bed.gz required}
COV=${3:?covariates file required}

echo "==== genotype chromosomes ===="
cut -f1 "${GENO}.bim" | sort -V | uniq -c

echo "==== phenotype chromosomes ===="
zcat "$BED" | tail -n +2 | cut -f1 | sort -V | uniq -c

echo "==== shared chromosomes ===="
comm -12 \
  <(cut -f1 "${GENO}.bim" | sort -V | uniq) \
  <(zcat "$BED" | tail -n +2 | cut -f1 | sort -V | uniq)

echo "==== phenotype header ===="
zcat "$BED" | head -1 | cut -f1-10

echo "==== sample check ===="
python3 scripts/check_inputs.py \
  --genotype "$GENO" \
  --phenotype "$BED" \
  --covariates "$COV"


