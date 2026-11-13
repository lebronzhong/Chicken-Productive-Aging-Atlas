#!/usr/bin/env bash
# Chicken Productive Aging Atlas | 08_Comparative_and_complex_traits
# Prepare LDSC summary statistics using the recovered study command.
# Inputs: LDSC directory, GWAS summary file, HapMap3 allele list, output prefix.
set -euo pipefail
if [[ $# -ne 4 ]]; then
    echo "Usage: $0 LDSC_DIR GWAS_SUMMARY HM3_ALLELE_LIST OUTPUT_PREFIX" >&2
    exit 2
fi
ldsc="$1"; sumstats="$2"; alleles="$3"; out="$4"
for file in "$ldsc/munge_sumstats.py" "$sumstats" "$alleles"; do
    [[ -s "$file" ]] || { echo "Missing input: $file" >&2; exit 2; }
done
mkdir -p "$(dirname "$out")"
python "$ldsc/munge_sumstats.py" --sumstats "$sumstats" --frq AF1 \
    --chunksize 500000 --merge-alleles "$alleles" --out "$out"
