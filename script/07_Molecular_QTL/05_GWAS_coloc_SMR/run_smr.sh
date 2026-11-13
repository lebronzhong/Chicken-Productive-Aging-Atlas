#!/usr/bin/env bash
# Chicken Productive Aging Atlas | 07_Molecular_QTL
# Run SMR/HEIDI with harmonized GWAS summary statistics and prepared BESD input.
set -euo pipefail
if [[ $# -lt 5 || $# -gt 6 ]]; then
    echo "Usage: $0 SMR_EXECUTABLE PLINK_PREFIX GWAS_SUMMARY BESD_PREFIX OUTPUT_PREFIX [THREADS]" >&2
    exit 2
fi
smr="$1"; bfile="$2"; gwas="$3"; besd="$4"; out="$5"; threads="${6:-1}"
[[ -x "$smr" ]] || { echo "SMR executable not found: $smr" >&2; exit 2; }
[[ "$threads" =~ ^[1-9][0-9]*$ ]] || { echo "THREADS must be positive" >&2; exit 2; }
for file in "$bfile.bed" "$bfile.bim" "$bfile.fam" "$gwas" "$besd.besd" "$besd.esi" "$besd.epi"; do
    [[ -s "$file" ]] || { echo "Missing or empty input: $file" >&2; exit 2; }
done
# The study GWAS input uses the eight-column SMR summary schema.
awk 'NR==1 {if(NF!=8 || $1!="SNP" || $2!="A1" || $3!="A2" || $4!="freq" || $5!="b" || $6!="se" || $7!="p" || $8!="n") exit 1; next}
     NF!=8 {exit 1} END {if(NR<2) exit 1}' "$gwas" || { echo "Expected GWAS columns: SNP A1 A2 freq b se p n" >&2; exit 2; }
mkdir -p "$(dirname "$out")"
"$smr" --bfile "$bfile" --gwas-summary "$gwas" --beqtl-summary "$besd" \
    --out "${out}.SMR" --thread-num "$threads" --diff-freq 0.8 \
    --diff-freq-prop 0.01 --heidi-mtd 1 --peqtl-smr 1e-6
