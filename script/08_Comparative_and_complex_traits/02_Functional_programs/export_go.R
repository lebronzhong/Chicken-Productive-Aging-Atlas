#!/usr/bin/env Rscript
# Chicken Productive Aging Atlas | 08_Comparative_and_complex_traits
# Purpose: Export go.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

args <- commandArgs(TRUE)
root <- args[1]
suppressPackageStartupMessages({library(AnnotationDbi); library(org.Gg.eg.db); library(GO.db)})
reference <- read.delim(file.path(root, '01_reference/mapped_tested_genes.tsv'))
valid_symbols <- intersect(reference$chicken_symbol, keys(org.Gg.eg.db, keytype = 'SYMBOL'))
bridge <- unique(AnnotationDbi::select(org.Gg.eg.db, keys = valid_symbols, keytype = 'SYMBOL', columns = 'ENTREZID'))
bridge <- bridge[!is.na(bridge$ENTREZID), ]
bridge <- merge(reference[, c('gene_id', 'chicken_symbol')], bridge, by.x = 'chicken_symbol', by.y = 'SYMBOL')
ambiguous <- duplicated(bridge$gene_id) | duplicated(bridge$gene_id, fromLast = TRUE) | duplicated(bridge$ENTREZID) | duplicated(bridge$ENTREZID, fromLast = TRUE)
bridge <- bridge[!ambiguous, ]
annotation <- AnnotationDbi::select(org.Gg.eg.db, keys = bridge$ENTREZID, keytype = 'ENTREZID', columns = c('GOALL', 'ONTOLOGYALL'))
annotation <- unique(annotation[!is.na(annotation$GOALL) & annotation$ONTOLOGYALL == 'BP', c('ENTREZID', 'GOALL')])
annotation <- merge(bridge, annotation, by = 'ENTREZID')
terms <- AnnotationDbi::select(GO.db, keys = unique(annotation$GOALL), keytype = 'GOID', columns = 'TERM')
annotation <- merge(annotation, terms, by.x = 'GOALL', by.y = 'GOID')
names(annotation)[names(annotation) == 'GOALL'] <- 'go_id'
write.table(annotation, gzfile(file.path(root, '01_reference/full_GO_BP.tsv.gz')), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(bridge, file.path(root, '01_reference/GO_symbol_bridge.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(metadata(org.Gg.eg.db), file.path(root, '02_audit/GO_metadata.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(root, '02_audit/annotation_sessionInfo.txt'))
