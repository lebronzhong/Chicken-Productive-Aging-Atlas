#!/usr/bin/env bash
# Chicken Productive Aging Atlas | 02_Expression_dynamics
# Purpose: Motif enrichment.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

export PATH=$HOME/meme/bin:$HOME/meme/libexec/meme-5.5.7:$PATH

TSI_path=/path/to/chicken-aging/analysis/motif
TSI_motif_path=/path/to/chicken-aging/ont/shell/14.5.LOESS_Cluster_ALL/2.Cluster/Cluster_Results_V3/MEME_motif
JASPAR_DB=/path/to/chicken-aging/software/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.meme
HOCOMOCO_DB=/path/to/chicken-aging/software/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme
REF_path=/path/to/chicken-aging/reference/galGal6/Ensemble/Gallus_gallus.GRCg6a.dna_sm.toplevel.fa
out_path=/path/to/chicken-aging/ont/shell/14.5.LOESS_Cluster_ALL/2.Cluster/Cluster_Results_V3/MEME_motif/MEME

for i in {1..10}
do

mkdir ${out_path}/Cluster_${i}_genes

ame --verbose 1 --oc ${out_path}/Cluster_${i}_genes \
${TSI_path}/Cluster_${i}_genes.promoters.fa \
${JASPAR_DB}

done
