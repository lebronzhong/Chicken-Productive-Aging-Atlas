#!/usr/bin/env bash
# Chicken Productive Aging Atlas | 08_Comparative_and_complex_traits
# Purpose: Make base.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

LDSC=/path/to/chicken-aging/software/ldsc
gene_coord_file=/path/to/chicken-aging/LDSC/ready_data/ENSG_coord.txt
plink_file=/path/to/chicken-aging/LDSC/ready_data/1000G_EUR_Phase3_plink
annot_file=/path/to/chicken-aging/LDSC/24.3.annot
baseline=/path/to/chicken-aging/LDSC/ready_data/baseline_by_chicken_aging

for chr in {1..22}
do
    bimfile=${plink_file}/1000G.EUR.QC.${chr}.bim
    {
        echo "BASELINE"
        awk '{print 1}' "$bimfile"
    } > ${baseline}/baseline.${chr}.annot

    gzip -f ${baseline}/baseline.${chr}.annot
done
