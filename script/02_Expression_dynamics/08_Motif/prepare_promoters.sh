#!/usr/bin/env bash
# Chicken Productive Aging Atlas | 02_Expression_dynamics
# Purpose: Prepare promoters.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

TSI_path=/path/to/chicken-aging/analysis/motif
JASPAR_DB=/path/to/chicken-aging/software/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.meme
HOCOMOCO_DB=/path/to/chicken-aging/software/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme
REF_path=/path/to/chicken-aging/reference/galGal6/Ensemble/Gallus_gallus.GRCg6a.dna_sm.toplevel.fa


for i in {1..10}
do

bedtools getfasta -fi ${REF_path} -bed ${TSI_path}/Cluster_${i}_genes.promoters.bed -fo ${TSI_path}/Cluster_${i}_genes.promoters.fa

done
