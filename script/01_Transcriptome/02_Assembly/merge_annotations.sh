#!/usr/bin/env bash
# Chicken Productive Aging Atlas | 01_Transcriptome
# Purpose: Transcript reconstruction and annotation merge.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

set -euo pipefail
annotation="${1:?Usage: bash merge_annotations.sh REFERENCE_GTF GTF_LIST OUTPUT_PREFIX}"
gtf_list="${2:?Provide one input GTF path per line}"
output_prefix="${3:?Provide the output prefix}"
gffcompare -r "${annotation}" -o "${output_prefix}" -i "${gtf_list}"
