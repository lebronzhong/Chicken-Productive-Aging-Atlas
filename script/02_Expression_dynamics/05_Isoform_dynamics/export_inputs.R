#!/usr/bin/env Rscript
# Chicken Productive Aging Atlas | 02_Expression_dynamics
# Purpose: Export inputs.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

.libPaths(c(file.path(R.home(),'library'),.libPaths()))
library(data.table)
args <- commandArgs(TRUE)
old <- args[1]
destination <- args[2]
setDTthreads(1)
packages <- c('ggplot2','patchwork','svglite','systemfonts','jsonlite')
available <- vapply(packages,requireNamespace,logical(1),quietly=FALSE)
print(available)
stopifnot(all(available))
stopifnot(capabilities('cairo'))
tissues <- readLines(file.path(old,'01_background/tissues.txt'))
manifest <- fread(file.path(old,'01_background/sample_manifest.tsv'))
sample_order <- unique(manifest[,.(sample,age,individual)])
sample_order[,bird_number:=as.integer(sub('.*_','',individual))]
setorder(sample_order,age,bird_number)
samples <- sample_order$sample
stopifnot(length(samples)==20)
fwrite(sample_order,file.path(destination,'01_inputs/sample_order.tsv'),sep='\t')
for(tissue in tissues) {
  folder <- file.path(destination,'01_inputs',tissue)
  dir.create(folder,recursive=TRUE,showWarnings=FALSE)
  for(feature in c('tx','gene')) {
    source <- readRDS(file.path(old,'02_age_expression',tissue,paste0(feature,'_four_age_data.rds')))
    stopifnot(setequal(colnames(source$normalized),samples),nrow(source$metadata)==20)
    expression <- as.data.table(log2(source$normalized[,samples]+.5),keep.rownames='feature_id')
    fwrite(expression,file.path(folder,paste0(feature,'_log_expression.tsv.gz')),sep='\t',compress='gzip')
    file.copy(file.path(old,'02_age_expression',tissue,paste0(feature,'_four_age.tsv.gz')),folder,overwrite=FALSE)
  }
  cat('EXPORTED',tissue,'\n')
}
writeLines(capture.output(sessionInfo()),file.path(destination,'08_report/R_sessionInfo.txt'))
