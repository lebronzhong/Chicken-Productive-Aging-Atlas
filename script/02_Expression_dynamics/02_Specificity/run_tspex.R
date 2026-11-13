#!/usr/bin/env Rscript
# Chicken Productive Aging Atlas | 02_Expression_dynamics
# Purpose: Specificity indices.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

### parse parameters
ARGS <- commandArgs(trailingOnly = TRUE)
tspex= ARGS[1] # where tspex is install (tspex source: https://github.com/apcamargo/tspex)
tpm_output = ARGS[2] # where TPM matrix is needed
prefix = ARGS[3] # output prefix

### run tspex
idx=c("counts", "tau", "gini", "simpson", "shannon_specificity", "roku_specificity", "tsi", "zscore", "spm", "spm_dpm", "js_specificity", "js_specificity_dpm")
for (x in idx){
 cmd=paste0(tspex, " --log ", tpm_output," ", prefix,".", x,".txt ", x)
 system(cmd)
}
