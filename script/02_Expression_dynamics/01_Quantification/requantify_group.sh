#!/usr/bin/env bash
# Chicken Productive Aging Atlas | 02_Expression_dynamics
# Purpose: Requantification with the study annotation.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

set -euo pipefail
tissue="${1:?Usage: bash requantify_group.sh TISSUE INPUT_YAML GENOME_FASTA STUDY_GTF_OR_DB OUTPUT_ROOT}"
input_yaml="${2:?Provide the grouped-library YAML}"
genome="${3:?Provide the reference FASTA}"
annotation="${4:?Provide the newly constructed study GTF or its IsoQuant database}"
output_root="${5:?Provide the output root}"
isoquant.py -d nanopore --yaml "${input_yaml}" --genedb "${annotation}" \
    --reference "${genome}" --output "${output_root}/${tissue}" --threads 24
