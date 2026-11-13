#!/usr/bin/env Rscript
# Chicken Productive Aging Atlas | 03_Splicing_and_usage
# Purpose: Count-based DTU.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

suppressPackageStartupMessages({library(edgeR);library(data.table)})
args<-commandArgs(TRUE);root<-normalizePath(args[1]);tissue<-args[2]
setDTthreads(2);set.seed(20261025)
out<-file.path(root,'results',tissue);dir.create(out,recursive=TRUE,showWarnings=FALSE)
readmat<-function(path){d<-fread(path);ids<-d[[1]];m<-as.matrix(d[,-1]);rownames(m)<-ids;storage.mode(m)<-'double';m}
input_dir<-file.path(root,'inputs','validated_fig3',tissue)
stopifnot(file.exists(file.path(input_dir,'VALIDATED')))
tx<-readmat(file.path(input_dir,'tx_counts.tsv.gz'))
ge<-readmat(file.path(input_dir,'gene_counts.tsv.gz'))
mp<-fread(file.path(root,'inputs','canonical_transcript_mapping.tsv.gz'));mp<-mp[match(rownames(tx),transcript_id)]
sm<-fread(file.path(input_dir,'samples.tsv'));stopifnot(identical(colnames(tx),sm$sample),identical(colnames(ge),sm$sample),!anyDuplicated(mp$transcript_id))
write_tsv<-function(d,name)fwrite(d,file.path(out,name),sep='\t',na='NA')
fit_one<-function(ages,label){
  cat('START',tissue,label,format(Sys.time()),'\n');flush.console()
  ss<-sm[age %in% ages];x<-tx[,ss$sample,drop=FALSE];g<-ge[,ss$sample,drop=FALSE]
  age<-factor(ss$age,levels=ages);design<-model.matrix(~age)
  stopifnot(qr(design)$rank==length(ages),all(table(age)==5))
  totals<-rowsum(x,mp$gene_id,reorder=FALSE);den<-totals[match(mp$gene_id,rownames(totals)),,drop=FALSE]
  fraction<-x/den;fraction[den==0]<-NA_real_
  group_mean<-sapply(ages,function(a)rowMeans(fraction[,ss$age==a,drop=FALSE],na.rm=TRUE));colnames(group_mean)<-paste0('mean_',ages)
  tx_support<-sapply(ages,function(a)rowSums(x[,ss$age==a,drop=FALSE]>=5 & fraction[,ss$age==a,drop=FALSE]>=.05,na.rm=TRUE)>=3)
  gene_support<-sapply(ages,function(a)rowSums(totals[,ss$age==a,drop=FALSE]>=20)>=3)
  keep_tx<-rowSums(tx_support)>0
  supported_gene<-rownames(totals)[rowSums(gene_support)==length(ages)]
  nret<-table(mp$gene_id[keep_tx]);eligible<-intersect(supported_gene,names(nret)[nret>=2])
  eff<-data.table(transcript_id=mp$transcript_id,gene_id=mp$gene_id,eligible_isoform=keep_tx,gene_DTU_eligible=mp$gene_id %in% eligible,group_mean)
  eff[,usage_range:=apply(.SD,1,function(z)if(all(is.finite(z)))max(z)-min(z) else NA_real_),.SDcols=paste0('mean_',ages)]
  write_tsv(eff,paste0(label,'_isoform_effects.tsv.gz'))
  # Denominator includes every observed canonical transcript, including low-use isoforms.
  write_tsv(data.table(transcript_id=mp$transcript_id,gene_id=mp$gene_id,fraction),paste0(label,'_sample_fractions.tsv.gz'))
  bk<-data.table(gene_id=union(rownames(g),rownames(totals)))
  bk[,total_tx:=as.integer(table(mp$gene_id)[gene_id])]
  bk[,retained_tx:=as.integer(nret[gene_id])];bk[is.na(retained_tx),retained_tx:=0L]
  bk[,DTU_testable:=gene_id %in% eligible]
  bk[,DTU_reason:=fifelse(DTU_testable,'tested',fifelse(!gene_id %in% rownames(totals),'no_transcript_count',fifelse(!gene_id %in% supported_gene,'insufficient_gene_counts_in_one_or_more_ages','fewer_than_two_supported_isoforms')))]
  gene_y<-DGEList(g);gene_keep<-filterByExpr(gene_y,group=age,min.count=10,min.total.count=15)
  bk[,DGE_testable:=gene_id %in% rownames(g)[gene_keep]]
  write_tsv(bk,paste0(label,'_testability.tsv.gz'))
  gene_y<-calcNormFactors(gene_y[gene_keep,,keep.lib.sizes=FALSE]);gene_y<-estimateDisp(gene_y,design,robust=TRUE)
  gene_fit<-glmQLFit(gene_y,design,robust=TRUE,legacy=FALSE)
  dge_omni<-glmQLFTest(gene_fit,coef=2:ncol(design))$table
  dge_omni<-data.table(gene_id=rownames(dge_omni),dge_omni)
  dge_omni[,q_tissue:=p.adjust(PValue,'BH')];write_tsv(dge_omni,paste0(label,'_DGE_omnibus.tsv.gz'))
  # Rare isoforms are pooled as a residual category, so the model conserves each gene's total.
  selected<-keep_tx & mp$gene_id %in% eligible
  xfit<-x[selected,,drop=FALSE];ann<-data.table(gene_id=mp$gene_id[selected],transcript_id=rownames(xfit),residual_category=FALSE)
  remaining<-rowsum(x[!selected & mp$gene_id %in% eligible,,drop=FALSE],mp$gene_id[!selected & mp$gene_id %in% eligible],reorder=FALSE)
  remaining<-remaining[rowSums(remaining)>0,,drop=FALSE]
  if(nrow(remaining)){
    ann<-rbind(ann,data.table(gene_id=rownames(remaining),transcript_id=paste0('OTHER_',rownames(remaining)),residual_category=TRUE))
    rownames(remaining)<-paste0('OTHER_',rownames(remaining));xfit<-rbind(xfit,remaining)
  }
  stopifnot(max(abs(rowsum(xfit,ann$gene_id,reorder=FALSE)-totals[unique(ann$gene_id),]))<1e-8)
  y<-DGEList(xfit,genes=as.data.frame(ann),lib.size=colSums(x));y<-calcNormFactors(y)
  y<-estimateDisp(y,design,robust=TRUE);fit<-glmQLFit(y,design,robust=TRUE,legacy=FALSE)
  pairs<-t(combn(ages,2));allgene<-list();alltx<-list();alldge<-list()
  for(j in seq_len(nrow(pairs))){
    younger<-pairs[j,1];older<-pairs[j,2];nm<-paste0(older,'vs',younger)
    contrast<-rep(0,ncol(design));if(older!=ages[1])contrast[match(older,ages)]<-1;if(younger!=ages[1])contrast[match(younger,ages)]<--1
    cat('DTU',tissue,label,nm,format(Sys.time()),'\n');flush.console()
    ds<-diffSplice(fit,contrast=contrast,geneid='gene_id',exonid='transcript_id',verbose=FALSE)
    gr<-as.data.table(ds$gene.genes);gr[,`:=`(contrast=nm,F=as.numeric(ds$gene.F),PValue=as.numeric(ds$gene.F.p.value))]
    tr<-as.data.table(ds$genes);tr[,`:=`(contrast=nm,relative_logFC=as.numeric(ds$coefficients),PValue=as.numeric(ds$p.value))]
    tr[,delta_fraction:=eff[[paste0('mean_',older)]][match(transcript_id,eff$transcript_id)]-eff[[paste0('mean_',younger)]][match(transcript_id,eff$transcript_id)]]
    ae<-tr[!grepl('^OTHER_',transcript_id),.(max_abs_delta=max(abs(delta_fraction),na.rm=TRUE)),by=gene_id]
    gr<-merge(gr,ae,by='gene_id',all.x=TRUE);gr[,q_tissue:=p.adjust(PValue,'BH')]
    allgene[[nm]]<-gr;alltx[[nm]]<-tr
    dg<-glmQLFTest(gene_fit,contrast=contrast)$table;dg<-data.table(gene_id=rownames(dg),contrast=nm,dg);alldge[[nm]]<-dg
  }
  gene_results<-rbindlist(allgene);tx_results<-rbindlist(alltx);dge_results<-rbindlist(alldge)
  write_tsv(gene_results,paste0(label,'_DTU_contrasts.tsv.gz'));write_tsv(tx_results,paste0(label,'_isoform_contrasts.tsv.gz'));write_tsv(dge_results,paste0(label,'_DGE_contrasts.tsv.gz'))
  # Bonferroni combination of three spanning contrasts tests the global equality null.
  # The contrasts themselves need not be statistically independent for this union bound.
  reference_contrasts<-paste0(ages[-1],'vs',ages[1])
  omni<-gene_results[contrast %in% reference_contrasts,.(PValue=if(all(is.finite(PValue)))min(1,length(reference_contrasts)*min(PValue)) else NA_real_,n_reference_tests=.N),by=gene_id]
  range_gene<-eff[eligible_isoform & gene_DTU_eligible,.(max_usage_range=max(usage_range,na.rm=TRUE)),by=gene_id]
  omni<-merge(omni,range_gene,by='gene_id');omni[,q_tissue:=p.adjust(PValue,'BH')]
  write_tsv(omni,paste0(label,'_DTU_omnibus.tsv.gz'))
  write_tsv(data.table(sample=ss$sample,age=ss$age,batch=ss$batch,tx_library=colSums(x),gene_library=colSums(g),tx_TMM=y$samples$norm.factors,gene_TMM=gene_y$samples$norm.factors),paste0(label,'_normalization.tsv'))
  normalized_gene<-cpm(gene_y,log=FALSE);normalized_tx<-sweep(x,2,y$samples$lib.size*y$samples$norm.factors,'/')*1e6
  write_tsv(data.table(gene_id=rownames(normalized_gene),normalized_gene),paste0(label,'_gene_CPM.tsv.gz'))
  write_tsv(data.table(transcript_id=rownames(normalized_tx),normalized_tx),paste0(label,'_transcript_CPM.tsv.gz'))
  saveRDS(list(gene_fit=gene_fit,tx_fit=fit,design=design,samples=ss,settings=list(min_tx_count=5,min_tx_fraction=.05,min_replicates=3,min_gene_count=20,ages=ages)),file.path(out,paste0(label,'_models.rds')),compress='gzip')
  cat('DONE',tissue,label,'eligible genes',length(eligible),'isoforms',sum(selected),'\n');flush.console()
}
fit_one(c(20,50,80,100),'four_age')
fit_one(c(50,100),'same_batch_50_100')
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
writeLines('SUCCESS',file.path(out,'SUCCESS'))
