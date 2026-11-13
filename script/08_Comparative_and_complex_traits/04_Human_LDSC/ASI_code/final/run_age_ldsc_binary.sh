#!/usr/bin/env bash
# Chicken Productive Aging Atlas | 08_Comparative_and_complex_traits
# Purpose: Run age ldsc binary.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

LDSC=/path/to/chicken-aging/software/ldsc
gene_coord_file=/path/to/chicken-aging/LDSC/ready_data/ENSG_coord.txt
plink_file=/path/to/chicken-aging/LDSC/ready_data/1000G_EUR_Phase3_plink
hapmap3snp=/path/to/chicken-aging/LDSC/ready_data/hapmap3_snps 
in_dir=/path/to/chicken-aging/LDSC/24.7.chicken_tissue_ASI/4.top5_to_hg19genelist_ensg
annot_file=/path/to/chicken-aging/LDSC/24.8.ASI_annot
merged_annot=/path/to/chicken-aging/LDSC/24.8.5.ASI_merged_annot
LDscore=/path/to/chicken-aging/LDSC/24.9.ASI_LDscore
result=/path/to/chicken-aging/LDSC/24.10.LDSC_for_tissue_ASI_final

sumstat=/path/to/chicken-aging/LDSC/GWAS_sumstats_binary       
frq_file=/path/to/chicken-aging/LDSC/ready_data/1000G_Phase3_frq

for i in `cat /path/to/chicken-aging/LDSC/shell/ASI_code/final/final_TSI2.txt`
do
date

python ${LDSC}/ldsc.py \
  --h2 ${sumstat}/${i}.sumstats.gz \
  --ref-ld-chr ${LDscore}/merged_ASI. \
  --out ${result}/${i} \
  --overlap-annot \
  --frqfile-chr /path/to/chicken-aging/LDSC/ready_data/1000G_Phase3_frq/1000G.EUR.QC. \
  --w-ld-chr /path/to/chicken-aging/LDSC/ready_data/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
  --print-coefficients

date
done
