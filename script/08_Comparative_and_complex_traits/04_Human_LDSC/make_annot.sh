#!/usr/bin/env bash
# Chicken Productive Aging Atlas | 08_Comparative_and_complex_traits
# Purpose: Make annot.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

LDSC=/path/to/chicken-aging/software/ldsc
gene_coord_file=/path/to/chicken-aging/LDSC/ready_data/ENSG_coord.txt
plink_file=/path/to/chicken-aging/LDSC/ready_data/1000G_EUR_Phase3_plink
annot_file=/path/to/chicken-aging/LDSC/24.3.annot
in_dir=/path/to/chicken-aging/LDSC/24.2.hg_liftover

for tissue in `cat /path/to/chicken-aging/LDSC/shell/tissue_list`
do
for i in {1..22}
do
python ${LDSC}/make_annot.py \
                --gene-set-file ${in_dir}/${tissue}_genelist.txt \
                --gene-coord-file ${gene_coord_file} \
                --windowsize 100000 \
                --bimfile ${plink_file}/1000G.EUR.QC.${i}.bim \
                --annot-file ${annot_file}/${tissue}.${i}.annot.gz
done
done
