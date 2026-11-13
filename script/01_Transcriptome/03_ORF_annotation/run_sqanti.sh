#!/usr/bin/env bash
# Chicken Productive Aging Atlas | 01_Transcriptome
# Purpose: Run sqanti.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

export PYTHONPATH="/path/to/chicken-aging/software/cDNA_Cupcake-master/sequence/:$PYTHONPATH"

soft_path=/path/to/chicken-aging/software/SQANTI3-5.1.2
ref_path=/path/to/chicken-aging/reference/galGal6/Ensemble/Gallus_gallus.GRCg6a.dna_sm.toplevel.fa
gtf_path=/path/to/chicken-aging/reference/galGal6/Ensemble/Gallus_gallus.GRCg6a.106.gtf

### Using gffcompare to make ONT downstream analysis ###

InFilePath=/path/to/chicken-aging/ont/shell/10.gffcompare/GTF_gather_extended
OutFilePath=/path/to/chicken-aging/ont/shell/10.5.SQANTI/ALL

python ${soft_path}/sqanti3_qc.py /path/to/chicken-aging/ont/shell/10.5.SQANTI/822_merge_ALL_corrected.combined.gtf \
${gtf_path} \
${ref_path} \
-o 822_merge_ALL.combined \
-d ${OutFilePath} \
--cpus 20 --report both
