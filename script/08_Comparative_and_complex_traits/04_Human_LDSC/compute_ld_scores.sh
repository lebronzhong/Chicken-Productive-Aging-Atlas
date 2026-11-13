#!/usr/bin/env bash
# Chicken Productive Aging Atlas | 08_Comparative_and_complex_traits
# Purpose: Compute ld scores.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

LDSC=/path/to/chicken-aging/software/ldsc
gene_coord_file=/path/to/chicken-aging/LDSC/ready_data/ENSG_coord.txt
plink_file=/path/to/chicken-aging/LDSC/ready_data/1000G_EUR_Phase3_plink
hapmap3snp=/path/to/chicken-aging/LDSC/ready_data/hapmap3_snps 
annot_file=/path/to/chicken-aging/LDSC/24.3.annot
in_dir=/path/to/chicken-aging/LDSC/24.2.hg_liftover
LDscore=/path/to/chicken-aging/LDSC/24.4.LDscore_new
baseline=/path/to/chicken-aging/LDSC/ready_data/baseline_by_chicken_aging
merged_annot=/path/to/chicken-aging/LDSC/24.35.merged_annot    

for chr in {1..22}
do
python ${LDSC}/ldsc.py \
        --l2 \
        --bfile ${plink_file}/1000G.EUR.QC.${chr} \
        --ld-wind-cm 1 \
        --annot ${merged_annot}/merged.${chr}.annot.gz \
        --thin-annot \
        --out ${LDscore}/merged.${chr} \
        --print-snps ${hapmap3snp}/hm.${chr}.snp
done
