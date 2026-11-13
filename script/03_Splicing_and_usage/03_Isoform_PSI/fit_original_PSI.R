#!/usr/bin/env Rscript
# Chicken Productive Aging Atlas | 03_Splicing_and_usage
# Purpose: PSI statistics.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

suppressPackageStartupMessages({library(data.table);library(limma)})
args<-commandArgs(TRUE);root<-normalizePath(args[1]);tissue<-args[2];setDTthreads(1)
inp<-file.path(root,'inputs',tissue);out<-file.path(root,'results',tissue);dir.create(out,recursive=TRUE,showWarnings=FALSE)
ann<-fread(file.path(inp,'features.tsv.gz'));d<-fread(file.path(inp,'original_PSI.tsv.gz'));sm<-fread(file.path(inp,'samples.tsv'))
stopifnot(identical(d$psi_id,ann$psi_id),identical(names(d)[-1],sm$sample),all(table(sm$age)==5))
y<-as.matrix(d[ann$testable==TRUE,-1]);rownames(y)<-ann[testable==TRUE]$psi_id
age<-factor(sm$age,levels=c(20,50,80,100));design<-model.matrix(~0+age);colnames(design)<-paste0('a',levels(age))
contrasts<-makeContrasts(`50vs20`=a50-a20,`80vs20`=a80-a20,`100vs20`=a100-a20,`80vs50`=a80-a50,`100vs50`=a100-a50,`100vs80`=a100-a80,levels=design)
allres<-list()
for(scale in c('arcsine_sqrt','raw_PSI')){
 z<-if(scale=='arcsine_sqrt')asin(sqrt(pmin(pmax(y,0),1)))else y
 base<-lmFit(z,design)
 fit<-eBayes(contrasts.fit(base,contrasts),robust=TRUE,trend=TRUE)
 for(j in seq_len(ncol(contrasts))){
  r<-data.table(psi_id=rownames(y),scale=scale,contrast=colnames(contrasts)[j],coefficient=fit$coefficients[,j],SE=sqrt(fit$s2.post)*fit$stdev.unscaled[,j],t=fit$t[,j],PValue=fit$p.value[,j],df=fit$df.total)
  allres[[length(allres)+1]]<-r
 }
 # Three independent spanning contrasts give the 3-df four-age test.
 om<-eBayes(contrasts.fit(base,contrasts[,1:3,drop=FALSE]),robust=TRUE,trend=TRUE)
 allres[[length(allres)+1]]<-data.table(psi_id=rownames(y),scale=scale,contrast='four_age',coefficient=NA_real_,SE=NA_real_,t=NA_real_,PValue=om$F.p.value,df=om$df.total,F=om$F)
}
fwrite(rbindlist(allres,fill=TRUE),file.path(out,'PSI_tests.tsv.gz'),sep='\t',na='NA')
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'));writeLines('SUCCESS',file.path(out,'SUCCESS'));cat('PSI MODEL COMPLETE',tissue,nrow(y),'\n')
