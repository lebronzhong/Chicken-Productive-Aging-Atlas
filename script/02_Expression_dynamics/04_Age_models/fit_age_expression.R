#!/usr/bin/env Rscript
# Chicken Productive Aging Atlas | 02_Expression_dynamics
# Purpose: Age-related expression model.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

.libPaths(c(file.path(R.home(),'library'),.libPaths()))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(DESeq2))
args <- commandArgs(TRUE)
setwd(args[1])
tissue <- args[2]
setDTthreads(1)
directory <- file.path('02_age_expression',tissue)
dir.create(directory,recursive=TRUE,showWarnings=FALSE)
manifest <- fread('01_background/sample_manifest.tsv')[get('tissue')==args[2]]
stopifnot(nrow(manifest)==20)
metadata <- as.data.frame(manifest)
rownames(metadata) <- metadata$sample
metadata$age <- factor(metadata$age,levels=c(20,50,80,100))
contrasts <- list(c(50,20),c(80,20),c(100,20),c(80,50),c(100,50),c(100,80))
summaries <- list()
for(feature in c('gene','tx','legacy_gene')) {
  cat('START',tissue,feature,format(Sys.time()),'\n')
  input <- fread(file.path('01_background',tissue,paste0(feature,'_counts.tsv.gz')))
  identifiers <- input[[1]]
  counts <- as.matrix(input[,metadata$sample,with=FALSE])
  rownames(counts) <- identifiers
  storage.mode(counts) <- 'integer'
  cpm <- sweep(counts,2,colSums(counts),'/')*1e6
  support <- sapply(levels(metadata$age),function(age) rowSums(counts[,metadata$age==age,drop=FALSE]>=3 & cpm[,metadata$age==age,drop=FALSE]>=1))
  keep <- rowSums(support>=3)>=1
  counts <- counts[keep,,drop=FALSE]
  support <- support[keep,,drop=FALSE]
  fit <- DESeqDataSetFromMatrix(counts,metadata,~age)
  fit <- estimateSizeFactors(fit,type='poscounts')
  fit <- DESeq(fit,quiet=TRUE,parallel=FALSE,minReplicatesForReplace=Inf)
  output <- data.table(feature_id=rownames(fit),tissue=tissue)
  for(age in c('20','50','80','100')) output[[paste0('support_',age)]] <- support[,age]
  normalized <- counts(fit,normalized=TRUE)
  for(age in c('20','50','80','100')) output[[paste0('median_',age)]] <- apply(normalized[,metadata$age==age,drop=FALSE],1,median)
  for(contrast in contrasts) {
    tag <- paste0(contrast[1],'vs',contrast[2])
    result <- results(fit,contrast=c('age',as.character(contrast[1]),as.character(contrast[2])),alpha=.05)
    for(column in c('log2FoldChange','lfcSE','pvalue','padj')) output[[paste0(tag,'_',column)]] <- result[[column]]
  }
  lrt <- nbinomLRT(fit,reduced=~1,quiet=TRUE)
  overall <- results(lrt,alpha=.05)
  output[, `:=`(overall_pvalue=overall$pvalue,overall_padj=overall$padj)]
  fwrite(output,file.path(directory,paste0(feature,'_four_age.tsv.gz')),sep='\t',compress='gzip')
  saveRDS(list(counts=counts,normalized=normalized,support=support,metadata=metadata,size_factors=sizeFactors(fit)),file.path(directory,paste0(feature,'_four_age_data.rds')),compress=FALSE)
  primary_columns <- metadata$age %in% c('50','100')
  primary_keep <- rowSums(support[,c('50','100'),drop=FALSE]>=3)>=1
  primary_meta <- droplevels(metadata[primary_columns,,drop=FALSE])
  primary_counts <- counts[primary_keep,primary_columns,drop=FALSE]
  primary <- DESeqDataSetFromMatrix(primary_counts,primary_meta,~age)
  primary <- estimateSizeFactors(primary,type='poscounts')
  primary <- DESeq(primary,quiet=TRUE,parallel=FALSE,minReplicatesForReplace=Inf)
  primary_result <- results(primary,contrast=c('age','100','50'),alpha=.05)
  primary_table <- as.data.table(as.data.frame(primary_result),keep.rownames='feature_id')
  primary_table[, tissue:=tissue]
  primary_table[, `:=`(support_50=support[primary_keep,'50'],support_100=support[primary_keep,'100'])]
  primary_normalized <- counts(primary,normalized=TRUE)
  primary_log <- as.data.table(log2(primary_normalized+.5),keep.rownames='feature_id')
  fwrite(primary_log,file.path(directory,paste0(feature,'_50_100_log_expression.tsv.gz')),sep='\t',compress='gzip')
  primary_table[, `:=`(median_50=apply(primary_normalized[,primary_meta$age=='50',drop=FALSE],1,median),median_100=apply(primary_normalized[,primary_meta$age=='100',drop=FALSE],1,median))]
  fwrite(primary_table,file.path(directory,paste0(feature,'_50_100.tsv.gz')),sep='\t',compress='gzip')
  saveRDS(list(counts=primary_counts,normalized=primary_normalized,support=support[primary_keep,c('50','100'),drop=FALSE],metadata=primary_meta,size_factors=sizeFactors(primary)),file.path(directory,paste0(feature,'_50_100_data.rds')),compress=FALSE)
  summaries[[feature]] <- data.table(tissue=tissue,feature=feature,total_rows=length(identifiers),eligible_four_age=sum(keep),eligible_50_100=sum(primary_keep),testable_50_100=sum(is.finite(primary_result$padj)),significant_50_100=sum(primary_result$padj<.05 & abs(primary_result$log2FoldChange)>=1,na.rm=TRUE),four_age_significant=sum(overall$padj<.05,na.rm=TRUE))
  rm(fit,lrt,primary,counts,normalized,primary_normalized)
  gc()
}
fwrite(rbindlist(summaries),file.path(directory,'summary.tsv'),sep='\t')
writeLines(capture.output(sessionInfo()),file.path(directory,'sessionInfo.txt'))
cat('DONE',tissue,format(Sys.time()),'\n')
