#!/usr/bin/env Rscript
# Chicken Productive Aging Atlas | 04_FIM
# Purpose: Usage recovery.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

suppressPackageStartupMessages({library(data.table);library(limma);library(edgeR)})
args<-commandArgs(TRUE);root<-normalizePath(args[1]);tissue<-args[2];setDTthreads(1)
inp<-file.path(root,'inputs',tissue);out<-file.path(root,'results',tissue)
readmat<-function(p){d<-fread(p);m<-as.matrix(d[,-1]);rownames(m)<-d[[1]];storage.mode(m)<-'double';m}
write<-function(d,n)fwrite(d,file.path(out,n),sep='\t',na='NA')
for(layer in c('original_isoform','event')){
 y<-readmat(file.path(inp,layer,'FIM_test_matrix.tsv.gz'));design<-model.matrix(~factor(rep(c('control','FIM'),each=5),levels=c('control','FIM')))
 if(nrow(y)){
  fit<-eBayes(lmFit(asin(sqrt(pmin(pmax(y,0),1))),design),robust=TRUE,trend=TRUE)
  write(data.table(feature_id=rownames(y),coefficient=fit$coefficients[,2],SE=sqrt(fit$s2.post)*fit$stdev.unscaled[,2],t=fit$t[,2],PValue=fit$p.value[,2],df=fit$df.total),paste0(layer,'_FIM_tests.tsv.gz'))
 }
 if(layer=='event'){
  y<-readmat(file.path(inp,layer,'natural_test_matrix.tsv.gz'));age<-factor(rep(c(20,50,80,100),each=5));dd<-model.matrix(~0+age);colnames(dd)<-paste0('a',levels(age))
  cc<-makeContrasts(`50vs20`=a50-a20,`80vs20`=a80-a20,`100vs20`=a100-a20,`80vs50`=a80-a50,`100vs50`=a100-a50,`100vs80`=a100-a80,levels=dd)
  base<-lmFit(asin(sqrt(pmin(pmax(y,0),1))),dd);fit<-eBayes(contrasts.fit(base,cc),robust=TRUE,trend=TRUE)
  rr<-lapply(seq_len(6),function(j)data.table(feature_id=rownames(y),contrast=colnames(cc)[j],coefficient=fit$coefficients[,j],SE=sqrt(fit$s2.post)*fit$stdev.unscaled[,j],t=fit$t[,j],PValue=fit$p.value[,j],df=fit$df.total))
  om<-eBayes(contrasts.fit(base,cc[,1:3]),robust=TRUE,trend=TRUE);rr[[7]]<-data.table(feature_id=rownames(y),contrast='four_age',PValue=om$F.p.value,F=om$F,df=om$df.total)
  write(rbindlist(rr,fill=TRUE),'event_natural_tests.tsv.gz')
 }
 cat('LIMMA_COMPLETE',tissue,layer,'\n');flush.console()
}
x<-readmat(file.path(inp,'count_DTU','FIM_count_matrix.tsv.gz'));ann<-fread(file.path(inp,'count_DTU','model_annotation.tsv'));lib<-fread(file.path(inp,'count_DTU','library_sizes.tsv'));stopifnot(identical(rownames(x),ann$transcript_id))
design<-model.matrix(~factor(rep(c('control','FIM'),each=5),levels=c('control','FIM')))
y<-DGEList(x,genes=as.data.frame(ann),lib.size=lib$library_size);y<-calcNormFactors(y);y<-estimateDisp(y,design,robust=TRUE);fit<-glmQLFit(y,design,robust=TRUE,legacy=FALSE)
ds<-diffSplice(fit,coef=2,geneid='gene_id',exonid='transcript_id',verbose=FALSE)
gr<-as.data.table(ds$gene.genes);gr[,`:=`(F=as.numeric(ds$gene.F),PValue=as.numeric(ds$gene.F.p.value))]
tr<-as.data.table(ds$genes);tr[,`:=`(relative_logFC=as.numeric(ds$coefficients),PValue=as.numeric(ds$p.value))]
write(gr,'count_DTU_FIM_gene_tests.tsv.gz');write(tr,'count_DTU_FIM_feature_tests.tsv.gz')
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'));writeLines('SUCCESS',file.path(out,'SUCCESS'));cat('COMPLETE',tissue,'\n')
