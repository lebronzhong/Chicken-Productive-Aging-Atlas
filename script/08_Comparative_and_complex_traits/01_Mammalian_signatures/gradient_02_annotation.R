#!/usr/bin/env Rscript
# Chicken Productive Aging Atlas | 08_Comparative_and_complex_traits
# Purpose: Gradient 02 annotation.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

root <- '/path/to/chicken-aging/analysis/conservation'
suppressPackageStartupMessages(library(AnnotationDbi))
suppressPackageStartupMessages(library(org.Gg.eg.db))
suppressPackageStartupMessages(library(GO.db))
m <- read.delim(file.path(root,'results/mapped_tested_genes.tsv'),check.names=FALSE)
valid <- intersect(m$gene_id,keys(org.Gg.eg.db,keytype='ENSEMBL'))
# Installed annotation has no direct matches to these GRCg6a IDs.
# Use only exact canonical chicken symbols with a one-to-one Entrez bridge.
valid_symbols <- intersect(unique(m$chicken_symbol),keys(org.Gg.eg.db,keytype='SYMBOL'))
stopifnot(length(valid_symbols)>0)
bridge <- AnnotationDbi::select(org.Gg.eg.db,keys=valid_symbols,keytype='SYMBOL',columns='ENTREZID')
bridge <- unique(bridge[!is.na(bridge$ENTREZID),])
bridge <- merge(m[,c('gene_id','chicken_symbol')],bridge,by.x='chicken_symbol',by.y='SYMBOL')
ambiguous <- unique(bridge$gene_id[duplicated(bridge$gene_id)|duplicated(bridge$gene_id,fromLast=TRUE)|duplicated(bridge$ENTREZID)|duplicated(bridge$ENTREZID,fromLast=TRUE)])
bridge <- bridge[!bridge$gene_id%in%ambiguous,]
bridge$mapping_method <- 'unique_exact_chicken_symbol_to_Entrez'
write.table(bridge,file.path(root,'reference/chicken_annotation_ID_bridge.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
x <- AnnotationDbi::select(org.Gg.eg.db,keys=bridge$ENTREZID,keytype='ENTREZID',columns=c('GOALL','ONTOLOGYALL'))
x <- x[!is.na(x$GOALL)&x$ONTOLOGYALL=='BP',]
x <- merge(bridge[,c('gene_id','ENTREZID')],x,by='ENTREZID')
x<-unique(x[,c('gene_id','ENTREZID','GOALL')]);names(x)<-c('gene_id','chicken_entrez_id','go_id')
terms<-AnnotationDbi::select(GO.db,keys=unique(x$go_id),keytype='GOID',columns=c('TERM','ONTOLOGY'))
x<-merge(x,terms,by.x='go_id',by.y='GOID',all.x=TRUE)
write.table(x,file.path(root,'reference/chicken_GO_BP_annotations.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
status<-data.frame(gene_id=m$gene_id,has_ensembl_annotation=m$gene_id%in%valid,has_exact_symbol_bridge=m$gene_id%in%bridge$gene_id,ambiguous_entrez=m$gene_id%in%ambiguous,has_BP_annotation=m$gene_id%in%x$gene_id)
write.table(status,file.path(root,'reference/annotation_coverage.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(metadata(org.Gg.eg.db),file.path(root,'reference/org_Gg_metadata.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(root,'reference/annotation_sessionInfo.txt'))
cat('ANNOTATION: mapped=',nrow(m),' Ensembl matched=',length(valid),' unique symbol bridge=',nrow(bridge),' BP annotated=',length(unique(x$gene_id)),' ambiguous=',length(ambiguous),'\n')
