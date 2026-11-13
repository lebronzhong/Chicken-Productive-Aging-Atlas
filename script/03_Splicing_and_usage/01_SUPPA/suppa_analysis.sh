#!/usr/bin/env bash
# Chicken Productive Aging Atlas | 03_Splicing_and_usage
# Purpose: SUPPA events and isoform fractions.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

set -euo pipefail
stage="${1:?Usage: bash suppa_analysis.sh STAGE GTF OUTPUT_PREFIX [TPM_MATRIX]}"
annotation="${2:?Provide the study or reference GTF}"
prefix="${3:?Provide the output prefix}"
suppa="${SUPPA_SCRIPT:?Set SUPPA_SCRIPT to the installed suppa.py}"
case "${stage}" in
    events) python "${suppa}" generateEvents -i "${annotation}" -o "${prefix}" -f ioe -e SE SS MX RI FL ;;
    isoforms) python "${suppa}" generateEvents -i "${annotation}" -o "${prefix}" -f ioi ;;
    psi) python "${suppa}" psiPerIsoform -g "${annotation}" -e "${4:?Provide the TPM matrix}" -o "${prefix}" ;;
    *) printf 'Unknown stage: %s\n' "${stage}" >&2; exit 2 ;;
esac
