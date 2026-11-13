#!/usr/bin/env bash
# Chicken Productive Aging Atlas | 08_Comparative_and_complex_traits
# Purpose: Human trait heritability enrichment.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

LDSC=/path/to/chicken-aging/software/ldsc
gene_coord_file=/path/to/chicken-aging/LDSC/ready_data/ENSG_coord.txt
plink_file=/path/to/chicken-aging/LDSC/ready_data/1000G_EUR_Phase3_plink
hapmap3snp=/path/to/chicken-aging/LDSC/ready_data/hapmap3_snps 
annot_file=/path/to/chicken-aging/LDSC/24.3.annot
in_dir=/path/to/chicken-aging/LDSC/24.2.hg_liftover
LDscore=/path/to/chicken-aging/LDSC/24.4.LDscore
result=/path/to/chicken-aging/LDSC/24.5.LDSC_for_tissue        

sumstat=/path/to/chicken-aging/LDSC/GWAS_sumstats
frq_file=/path/to/chicken-aging/LDSC/ready_data/1000G_Phase3_frq

for i in `cat /path/to/chicken-aging/LDSC/GWAS_filename`       
do
date

python ${LDSC}/ldsc.py \
  --h2 ${sumstat}/${i}.sumstats.gz \
  --ref-ld-chr /path/to/chicken-aging/LDSC/24.4.LDscore_new/merged. \
  --out ${result}/${i} \
  --overlap-annot \
  --frqfile-chr /path/to/chicken-aging/LDSC/ready_data/1000G_Phase3_frq/1000G.EUR.QC. \
  --w-ld-chr /path/to/chicken-aging/LDSC/ready_data/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
  --print-coefficients

date
done
