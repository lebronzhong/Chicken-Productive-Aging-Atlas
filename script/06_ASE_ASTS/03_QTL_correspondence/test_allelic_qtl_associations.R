#!/usr/bin/env Rscript
# Chicken Productive Aging Atlas | 06_ASE_ASTS
# Purpose: Allelic-effect and QTL association.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

if(.Platform$OS.type=='windows') invisible(Sys.setlocale('LC_CTYPE','English_United States.utf8'))
suppressPackageStartupMessages(library(data.table))
read_table<-function(path){
 if(.Platform$OS.type!='windows' && grepl('\\.gz$',path)){
  tmp<-tempfile(pattern='verified_gunzip_',fileext='.tsv')
  on.exit(unlink(tmp),add=TRUE)
  rc<-system2('gzip',args=c('-dc','--',shQuote(path)),stdout=tmp)
  if(rc!=0L)stop('Decompression failed; refusing partial data: ',path)
  withCallingHandlers(fread(tmp),warning=function(w)stop('Input-read warning: ',conditionMessage(w)))
 } else fread(path)
}
args<-commandArgs(TRUE);root<-args[1];out<-if(length(args)>1)args[2]else file.path(root,'results','analysis');set.seed(20261026)
check_counts<-function(d,expected_col){
 auditfile<-file.path(out,'analysis_eligibility_audit.tsv')
 if(!file.exists(auditfile))return(invisible(NULL)) # Synthetic validation only.
 a<-fread(auditfile)[status=='ready',.(tissue,qtl_type,expected=get(expected_col))]
 actual<-d[version=='all450',.N,by=.(tissue,qtl_type)]
 z<-merge(a,actual,by=c('tissue','qtl_type'),all=TRUE)
 if(anyNA(z)||any(z$expected!=z$N))stop('Loaded row counts differ from independently generated eligibility audit: ',expected_col)
 cat('VERIFIED_INPUT_COUNTS',expected_col,nrow(z),sum(z$N),'\n')
}
cluster_vcov<-function(X,resid,bread,block){
 u<-rowsum(X*as.numeric(resid),block,reorder=FALSE);G<-nrow(u);N<-nrow(X);K<-ncol(X)
 if(G<20||N<=K)stop('fewer_than_20_genomic_blocks_or_insufficient_df')
 list(V=(G/(G-1))*((N-1)/(N-K))*bread%*%crossprod(u)%*%bread,G=G)
}
adjusted_or<-function(d){
 warnings<-character();ctrl<-c('log_coverage','log_donors','log_qtl_structures','MAF','log_expression')
 ctrl<-ctrl[vapply(d[,..ctrl],function(x)sd(x)>0,logical(1))]
 if(uniqueN(d$positive)<2||uniqueN(d$outcome)<2)return(list(status='no_outcome_or_exposure_variation'))
 # A zero exposure/outcome cell permits complete or quasi-complete separation.
 # glm can report convergence with a finite iteration-dependent coefficient;
 # sandwich Wald inference is invalid at that boundary.
 if(any(table(factor(d$positive,levels=0:1),factor(d$outcome,levels=0:1))==0))return(list(status='separation_zero_cell_not_estimable'))
 fit<-tryCatch(withCallingHandlers(glm(reformulate(c('positive',ctrl),response='outcome'),data=d,family=binomial()),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)e)
 if(inherits(fit,'error'))return(list(status=conditionMessage(fit)))
 if(!fit$converged||length(warnings))return(list(status=paste(c('unstable',warnings),collapse=';')))
 X<-model.matrix(fit);if(fit$rank<ncol(X))return(list(status='rank_deficient'))
 cv<-tryCatch(cluster_vcov(X,residuals(fit,type='response'),solve(crossprod(X,X*fit$weights)),d$block),error=function(e)e)
 if(inherits(cv,'error'))return(list(status=conditionMessage(cv)))
 b<-coef(fit)['positive'];se<-sqrt(cv$V['positive','positive']);crit<-qt(.975,cv$G-1)
 list(status='ok',adjusted_OR=exp(b),low=exp(b-crit*se),high=exp(b+crit*se),p=2*pt(-abs(b/se),cv$G-1),n_blocks=cv$G)
}
e<-read_table(file.path(out,'site_enrichment_input.tsv.gz'));expr<-read_table(file.path(root,'input/expression_covariate.tsv.gz'));e<-merge(e,expr,by=c('tissue','gene'),all.x=TRUE)
check_counts(e,'confirmed_gene_SNPs')
enrich<-list();nn<-0L
for(ver in unique(e$version))for(t in unique(e$tissue))for(qt0 in c('teQTL','tuQTL'))for(assay in c('ASE','ASTS'))for(label in c('QTL_supported','CS_supported','PIP_0.1','PIP_0.9'))for(selection in c('all_unique_gene_reporter_SNPs','one_SNP_per_gene')){
 d<-copy(e[version==ver&tissue==t&qtl_type==qt0&single_gene_reporter==TRUE]);prefix<-tolower(assay)
 d[,donors:=get(paste0(prefix,'_donors'))];d<-d[donors>=2];d[,positive:=as.integer(get(paste0(prefix,'_sig_donors'))>=2)]
 if(selection=='one_SNP_per_gene'){setorder(d,-median_coverage,variant_id);d<-d[!duplicated(gene)]}
 d[,outcome:=as.integer(if(label=='PIP_0.1')max_PIP>=.1 else if(label=='PIP_0.9')max_PIP>=.9 else get(label))]
 d[,`:=`(log_coverage=log1p(median_coverage),log_donors=log1p(donors),log_qtl_structures=log1p(n_QTL_tested_structures),MAF=pmin(af,1-af),log_expression=log1p(mean_natural_CPM))]
 Nbefore<-nrow(d);d<-d[is.finite(log_expression)&is.finite(MAF)]
 tb<-table(factor(d$positive,levels=0:1),factor(d$outcome,levels=0:1));ft<-tryCatch(fisher.test(tb),error=function(e)NULL)
 a<-if(nrow(d)>=20)adjusted_or(d)else list(status='too_few_sites')
 nn<-nn+1L;enrich[[nn]]<-c(list(version=ver,tissue=t,qtl_type=qt0,assay=assay,label=label,selection=selection,n=nrow(d),n_removed_missing_covariates=Nbefore-nrow(d),both=tb[2,2],allelic_only=tb[2,1],QTL_only=tb[1,2],neither=tb[1,1],Fisher_OR=if(is.null(ft))NA_real_ else unname(ft$estimate),Fisher_p_descriptive=if(is.null(ft))NA_real_ else ft$p.value),a)
}
er<-rbindlist(enrich,fill=TRUE);if(!'p'%in%names(er))er[,p:=NA_real_]
er[,q72:=p.adjust(p,'BH',n=72),by=.(version,label,selection)];fwrite(er,file.path(out,'site_enrichment_results.tsv'),sep='\t')
f<-read_table(file.path(out,'signed_effect_input.tsv.gz'));corr<-list();picked<-list();nn<-0L
check_counts(f,'exact_QTL_pairs')
for(ver in unique(f$version))for(t in unique(f$tissue))for(qt0 in c('teQTL','tuQTL'))for(endpoint in c('gene_ASE','transcript_usage','transcript_ASE'))for(subset0 in c('all_tested_pairs','QTL_supported_pairs')){
 d<-copy(f[version==ver&tissue==t&qtl_type==qt0&orientation_validated==TRUE&single_gene_reporter==TRUE&primary_structure_eligible==TRUE])
 if(subset0=='QTL_supported_pairs')d<-d[nominal_QTL_supported==TRUE]
 if(endpoint=='gene_ASE')d<-d[ase_donors>=2&is.finite(mean_ASE_log2)]
 if(endpoint=='transcript_usage')d<-d[usage_donors>=2&is.finite(mean_delta_usage)]
 if(endpoint=='transcript_ASE')d<-d[tx_ASE_donors>=2&is.finite(mean_tx_ASE_log2)]
 # Representative selection uses RNA coverage, with identifier tie-breaks, never effect agreement.
 setorder(d,-median_coverage,-mean_tx_coverage,variant_id,phenotype_id,na.last=TRUE);d<-d[!duplicated(gene)]
 d[,allelic_effect:=get(switch(endpoint,gene_ASE='mean_ASE_log2',transcript_usage='mean_delta_usage',transcript_ASE='mean_tx_ASE_log2'))]
 d<-d[is.finite(beta_ALT_provisional)&is.finite(allelic_effect)];d[,`:=`(endpoint=endpoint,subset=subset0)]
 if(nrow(d))picked[[length(picked)+1L]]<-d
 ans<-list(status='too_few_genes_or_constant_effect');G<-uniqueN(d$block)
 if(nrow(d)>=20&&G>=20&&sd(d$allelic_effect)>0&&sd(d$beta_ALT_provisional)>0){
  xr<-as.numeric(scale(rank(d$beta_ALT_provisional)));yr<-as.numeric(scale(rank(d$allelic_effect)));X<-cbind(Intercept=1,QTL_rank=xr);fit<-lm.fit(X,yr)
  cv<-cluster_vcov(X,fit$residuals,solve(crossprod(X)),d$block);rho<-fit$coefficients[2];se<-sqrt(cv$V[2,2]);p<-2*pt(-abs(rho/se),cv$G-1)
  blocks<-split(seq_len(nrow(d)),d$block);boot<-replicate(1000,{ix<-unlist(blocks[sample.int(length(blocks),length(blocks),replace=TRUE)],use.names=FALSE);suppressWarnings(cor(d$beta_ALT_provisional[ix],d$allelic_effect[ix],method='spearman'))})
  ci<-quantile(boot,c(.025,.975),na.rm=TRUE);ok<-d$beta_ALT_provisional*d$allelic_effect!=0
  ans<-list(status='ok',Spearman_rho=unname(rho),bootstrap_low=ci[1],bootstrap_high=ci[2],cluster_rank_regression_p=p,nonzero_pairs=sum(ok),direction_concordance=if(any(ok))mean(d$beta_ALT_provisional[ok]*d$allelic_effect[ok]>0)else NA_real_,bootstrap_replicates=sum(is.finite(boot)))
 }
 nn<-nn+1L;corr[[nn]]<-c(list(version=ver,tissue=t,qtl_type=qt0,endpoint=endpoint,subset=subset0,n_genes=nrow(d),n_blocks=G),ans)
}
cr<-rbindlist(corr,fill=TRUE);if(!'cluster_rank_regression_p'%in%names(cr))cr[,cluster_rank_regression_p:=NA_real_]
cr[,q108:=p.adjust(cluster_rank_regression_p,'BH',n=108),by=.(version,subset)];fwrite(cr,file.path(out,'signed_effect_correlations.tsv'),sep='\t')
if(length(picked))fwrite(rbindlist(picked,fill=TRUE),file.path(out,'correlation_selected_points.tsv.gz'),sep='\t')
capture.output(sessionInfo(),file=file.path(out,'R_sessionInfo.txt'))
writeLines('Interpretation gate: existing allele-assigned counts; no full read-mapping bias correction. Cluster units are fixed 1Mb genomic bins, not inferred LD blocks. Correlation P is the cluster-robust standardized-rank regression test; interval is a 1000-replicate genomic-block bootstrap. All effects use verified BIM-A1 to genomic ALT alignment.',file.path(out,'STATISTICAL_SCOPE.txt'))
