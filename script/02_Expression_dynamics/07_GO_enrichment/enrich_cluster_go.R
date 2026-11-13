#!/usr/bin/env Rscript
# Chicken Productive Aging Atlas | 02_Expression_dynamics
# Purpose: GO enrichment.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(AnnotationDbi)
  library(org.Gg.eg.db)
  library(GO.db)
})
data.table::setDTthreads(4)
stopifnot(as.character(packageVersion("org.Gg.eg.db")) == "3.22.0", as.character(packageVersion("GO.db")) == "3.22.0")
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(2L, 3L)) stop("Usage: Rscript enrich_cluster_go.R INPUT_DIR OUTPUT_DIR [EXPECTED_RECORDS]")
inputs <- args[1]
root <- args[2]
expected_records <- if (length(args) == 3L) as.integer(args[3]) else NA_integer_
if (length(args) == 3L && (is.na(expected_records) || expected_records < 1L)) stop("EXPECTED_RECORDS must be a positive integer")
for (folder in c("01_inputs", "02_audit", "03_pooled", "04_tissue", "05_summary", "06_figures", "07_report")) {
  dir.create(file.path(root, folder), recursive = TRUE, showWarnings = FALSE)
}
save_table <- function(value, relative_path) fwrite(value, file.path(root, relative_path), sep = "\t", na = "NA")
classification <- fread(file.path(inputs, "transcript_classification.tsv.gz"))
old <- fread(file.path(inputs, "canonical_old_cluster_assignments.tsv.gz"))
aliases <- c(Medulla = "Medulla_oblongata", Olfactory = "Olfactory_bulb", Thigh = "Thigh_muscle")
old[tissue %in% names(aliases), tissue := unname(aliases[tissue])]
stopifnot(!anyDuplicated(classification, by = c("tissue", "transcript_id")), !anyDuplicated(old, by = c("tissue", "transcript_id")))
canonical <- fread(file.path(inputs, "canonical_transcript_mapping.tsv.gz"), select = c("transcript_id", "gene_id"))
stopifnot(!anyDuplicated(canonical$transcript_id))
mapped <- merge(classification[, .(transcript_id, gene_id)], canonical, by = "transcript_id", suffixes = c("", "_canonical"))
stopifnot(nrow(mapped) == nrow(classification), all(mapped$gene_id == mapped$gene_id_canonical))
eligible <- merge(classification[common_background == TRUE], old[distinct_old_clusters == 1 & !is.na(old_cluster)], by = c("tissue", "transcript_id"))
eligible[, cluster := as.integer(old_cluster)]
primary <- eligible[age_associated == TRUE]
plotted <- fread(file.path(inputs, "primary_cluster_z.tsv.gz"), select = c("tissue", "transcript_id", "gene_id", "old_cluster"))
stopifnot(nrow(primary) > 0L, all(primary$cluster %in% seq_len(10L)))
if (!is.na(expected_records)) stopifnot(nrow(primary) == expected_records)
comparison <- merge(primary[, .(tissue, transcript_id, gene_id, old_cluster)], plotted, by = c("tissue", "transcript_id", "gene_id", "old_cluster"))
stopifnot(nrow(comparison) == nrow(primary), nrow(plotted) == nrow(primary))
save_table(primary[, .(tissue, transcript_id, gene_id, symbol, cluster, class)], "01_inputs/current_cluster_members.tsv.gz")
primary_genes <- sort(unique(primary$gene_id))
all_genes <- sort(unique(eligible$gene_id))
symbols <- fread(file.path(inputs, "gene_symbols.tsv"))[gene_id %in% all_genes]
stopifnot(!anyDuplicated(symbols$gene_id))
message("Reconstructed members: ", nrow(primary), "; unique genes: ", length(primary_genes))

direct_keys <- intersect(all_genes, keys(org.Gg.eg.db, keytype = "ENSEMBL"))
direct <- as.data.table(AnnotationDbi::select(org.Gg.eg.db, keys = direct_keys, keytype = "ENSEMBL", columns = "ENTREZID"))
setnames(direct, "ENSEMBL", "gene_id")
direct <- unique(direct[!is.na(ENTREZID)])
symbol_keys <- intersect(symbols$symbol, keys(org.Gg.eg.db, keytype = "SYMBOL"))
symbol_links <- as.data.table(AnnotationDbi::select(org.Gg.eg.db, keys = symbol_keys, keytype = "SYMBOL", columns = "ENTREZID"))
symbol_links <- unique(merge(symbols, symbol_links[!is.na(ENTREZID)], by.x = "symbol", by.y = "SYMBOL"))
direct_summary <- direct[, .(direct_n = uniqueN(ENTREZID), direct_ids = paste(sort(unique(ENTREZID)), collapse = ";")), by = gene_id]
symbol_summary <- symbol_links[, .(symbol_n = uniqueN(ENTREZID), symbol_ids = paste(sort(unique(ENTREZID)), collapse = ";")), by = gene_id]
bridge <- merge(data.table(gene_id = all_genes), symbols, by = "gene_id", all.x = TRUE)
bridge <- merge(bridge, direct_summary, by = "gene_id", all.x = TRUE)
bridge <- merge(bridge, symbol_summary, by = "gene_id", all.x = TRUE)
bridge[is.na(direct_n), direct_n := 0L]
bridge[is.na(symbol_n), symbol_n := 0L]
bridge[, `:=`(ENTREZID = NA_character_, mapping_status = "unmapped")]
bridge[direct_n == 1, `:=`(ENTREZID = direct_ids, mapping_status = "direct_ensembl")]
bridge[direct_n == 0 & symbol_n == 1, `:=`(ENTREZID = symbol_ids, mapping_status = "symbol_fallback")]
bridge[direct_n > 1 | (direct_n == 0 & symbol_n > 1), mapping_status := "ambiguous_mapping"]
conflicts <- symbol_links[gene_id %in% bridge[direct_n == 1 & symbol_n > 0, gene_id]]
conflicts <- merge(conflicts, bridge[, .(gene_id, direct_ids)], by = "gene_id")
conflicting_genes <- conflicts[, .(agree = any(ENTREZID == direct_ids)), by = gene_id][agree == FALSE, gene_id]
bridge[gene_id %in% conflicting_genes, mapping_status := "direct_symbol_conflict"]
bridge[!mapping_status %in% c("direct_ensembl", "symbol_fallback"), ENTREZID := NA_character_]
duplicate_ids <- bridge[!is.na(ENTREZID), .N, by = ENTREZID][N > 1, ENTREZID]
bridge[ENTREZID %in% duplicate_ids, `:=`(mapping_status = "multiple_gene_ids_per_entrez", ENTREZID = NA_character_)]
save_table(bridge, "02_audit/gene_annotation_mapping_audit.tsv")
save_table(bridge[, .N, by = mapping_status], "02_audit/mapping_status_counts.tsv")
stopifnot(!anyDuplicated(bridge[!is.na(ENTREZID)]$ENTREZID))
go <- as.data.table(AnnotationDbi::select(org.Gg.eg.db, keys = bridge[!is.na(ENTREZID), ENTREZID], keytype = "ENTREZID", columns = c("GOALL", "ONTOLOGYALL")))
go <- unique(go[!is.na(GOALL), .(ENTREZID, go_id = GOALL, ontology = ONTOLOGYALL)])
terms <- as.data.table(AnnotationDbi::select(GO.db, keys = unique(go$go_id), keytype = "GOID", columns = c("TERM", "ONTOLOGY")))
setnames(terms, c("GOID", "TERM", "ONTOLOGY"), c("go_id", "term", "term_ontology"))
go <- merge(go, bridge[!is.na(ENTREZID), .(gene_id, ENTREZID)], by = "ENTREZID")
go <- merge(go, terms, by = "go_id")
stopifnot(all(go$ontology == go$term_ontology), !anyNA(go$term))
go <- unique(go[, .(gene_id, ENTREZID, go_id, term, ontology)])
save_table(go, "01_inputs/full_gene_GOALL_annotation.tsv.gz")
save_table(as.data.table(metadata(org.Gg.eg.db)), "02_audit/org_Gg_metadata.tsv")
save_table(as.data.table(metadata(GO.db)), "02_audit/GO_metadata.tsv")
save_table(data.table(package = c("org.Gg.eg.db", "GO.db"), version = c(as.character(packageVersion("org.Gg.eg.db")), as.character(packageVersion("GO.db")))), "02_audit/annotation_versions.tsv")

gene_covariates <- eligible[, .(eligible_isoforms = uniqueN(transcript_id), eligible_tissues = uniqueN(tissue)), by = gene_id]
gene_covariates[, isoform_bin := as.character(cut(eligible_isoforms, c(0, 1, 2, 4, 8, Inf), labels = c("1", "2", "3-4", "5-8", "9+")))]
gene_covariates[, tissue_bin := as.character(cut(eligible_tissues, c(0, 5, 15, 25, 33), labels = c("1-5", "6-15", "16-25", "26-33")))]
gene_covariates[, stratum := paste(isoform_bin, tissue_bin, sep = ":")]
save_table(gene_covariates, "02_audit/gene_detection_covariates.tsv")
gene_clusters <- unique(primary[, .(gene_id, cluster)])
save_table(gene_clusters[, .(n_clusters = .N, clusters = paste(sort(cluster), collapse = ";")), by = gene_id], "02_audit/gene_cluster_multiplicity.tsv")
coverage <- list()
pooled_results <- list()
tissue_results <- list()

test_context <- function(selected_members, background_ids, ontology_name, context, tissue_name = "Pooled", use_strata = FALSE) {
  annotated_ids <- sort(intersect(background_ids, go[ontology == ontology_name, gene_id]))
  annotation <- go[ontology == ontology_name & gene_id %in% annotated_ids]
  term_counts <- annotation[, .(background_hits = uniqueN(gene_id), term = term[1]), by = go_id]
  term_counts <- term_counts[background_hits >= 10 & background_hits <= 2000 & background_hits <= length(annotated_ids) / 2]
  annotation <- annotation[go_id %in% term_counts$go_id]
  term_ids <- term_counts$go_id
  gene_term <- sparseMatrix(i = match(annotation$gene_id, annotated_ids), j = match(annotation$go_id, term_ids), x = 1,
                            dims = c(length(annotated_ids), length(term_ids)))
  selected_pairs <- unique(selected_members[gene_id %in% annotated_ids, .(gene_id, cluster)])
  memberships <- sparseMatrix(i = match(selected_pairs$gene_id, annotated_ids), j = selected_pairs$cluster, x = 1,
                              dims = c(length(annotated_ids), 10L))
  overlap <- as.matrix(crossprod(gene_term, memberships))
  foreground_sizes <- as.numeric(colSums(memberships))
  population <- length(annotated_ids)
  hits <- as.numeric(colSums(gene_term))
  results <- rbindlist(lapply(seq_len(10), function(cluster_id) {
    foreground <- foreground_sizes[cluster_id]
    data.table(context = context, ontology = ontology_name, tissue = tissue_name, cluster = cluster_id,
               go_id = term_ids, term = term_counts$term, hits = overlap[, cluster_id], foreground_genes = foreground,
               background_hits = hits, background_genes = population,
               fold_enrichment = (overlap[, cluster_id] / foreground) / (hits / population),
               p_value = phyper(overlap[, cluster_id] - 1, hits, population - hits, foreground, lower.tail = FALSE))
  }))
  results[foreground_genes == 0, `:=`(p_value = 1, fold_enrichment = NA_real_)]
  if (use_strata) {
    covariates <- gene_covariates[match(annotated_ids, gene_id)]
    expected <- matrix(0, nrow = length(term_ids), ncol = 10)
    variance <- expected
    for (stratum_id in unique(covariates$stratum)) {
      positions <- which(covariates$stratum == stratum_id)
      stratum_size <- length(positions)
      positives <- as.numeric(colSums(gene_term[positions, , drop = FALSE]))
      draws <- as.numeric(colSums(memberships[positions, , drop = FALSE]))
      expected <- expected + outer(positives / stratum_size, draws)
      if (stratum_size > 1) {
        variance <- variance + outer((positives / stratum_size) * (1 - positives / stratum_size), draws * (stratum_size - draws) / (stratum_size - 1))
      }
    }
    results[, `:=`(stratified_expected = as.vector(expected), stratified_variance = as.vector(variance))]
    results[, stratified_z := (hits - stratified_expected - 0.5) / sqrt(stratified_variance)]
    results[, stratified_p := pnorm(stratified_z, lower.tail = FALSE)]
    results[stratified_variance <= 0, stratified_p := 1]
    results[, stratified_low_expected := stratified_expected < 5]
  }
  backgrounds <- list(context = context, ontology = ontology_name, tissue = tissue_name,
                      background_genes_before_annotation = length(unique(background_ids)), background_genes_annotated = population,
                      tested_terms = length(term_ids))
  denominators <- rbindlist(lapply(seq_len(10), function(cluster_id) {
    c(backgrounds, list(cluster = cluster_id, foreground_genes_before_annotation = uniqueN(selected_members[cluster == cluster_id, gene_id]),
                       foreground_genes_annotated = foreground_sizes[cluster_id], transcript_tissue_members = nrow(selected_members[cluster == cluster_id])))
  }))
  list(tests = results, denominators = denominators)
}

resume_tests <- "--resume-tests" %in% commandArgs(TRUE)
if (!resume_tests) {
for (ontology_name in c("BP", "MF", "CC")) {
  message("Pooled ", ontology_name)
  for (context in c("age_associated_primary", "common_detectable_sensitivity")) {
    background <- if (context == "age_associated_primary") primary_genes else all_genes
    result <- test_context(primary, background, ontology_name, context, use_strata = context == "age_associated_primary")
    result$tests[, q_value := p.adjust(p_value, "BH")]
    if ("stratified_p" %in% names(result$tests)) result$tests[, stratified_q := p.adjust(stratified_p, "BH")]
    pooled_results[[paste(context, ontology_name)]] <- result$tests
    coverage[[paste(context, ontology_name)]] <- result$denominators
  }
  local_results <- list()
  for (tissue_name in sort(unique(primary$tissue))) {
    selected <- primary[tissue == tissue_name]
    result <- test_context(selected, unique(selected$gene_id), ontology_name, "within_tissue_primary", tissue_name)
    local_results[[tissue_name]] <- result$tests
    coverage[[paste("tissue", ontology_name, tissue_name)]] <- result$denominators
  }
  local_results <- rbindlist(local_results)
  local_results[, q_value := p.adjust(p_value, "BH")]
  tissue_results[[ontology_name]] <- local_results
  message("Tissue tests ", ontology_name, ": ", nrow(local_results))
}
pooled <- rbindlist(pooled_results, fill = TRUE)
tissue_tests <- rbindlist(tissue_results)
save_table(pooled, "03_pooled/all_pooled_GO_tests.tsv.gz")
save_table(tissue_tests, "04_tissue/all_tissue_GO_tests.tsv.gz")
save_table(pooled[q_value < 0.05 & hits >= 3], "03_pooled/significant_pooled_GO.tsv")
save_table(tissue_tests[q_value < 0.05 & hits >= 3], "04_tissue/significant_tissue_GO.tsv.gz")
save_table(rbindlist(coverage), "02_audit/foreground_background_denominators.tsv")
} else {
  pooled <- fread(file.path(root, "03_pooled/all_pooled_GO_tests.tsv.gz"))
  tissue_tests <- fread(file.path(root, "04_tissue/all_tissue_GO_tests.tsv.gz"))
}

bp_primary <- pooled[context == "age_associated_primary" & ontology == "BP"]
candidate_pool <- bp_primary[q_value < 0.05 & stratified_q < 0.05 & !stratified_low_expected & hits >= 5 & fold_enrichment >= 1.2]
annotation_sets <- split(go[ontology == "BP" & gene_id %in% primary_genes, gene_id], go[ontology == "BP" & gene_id %in% primary_genes, go_id])
representatives <- list()
for (cluster_id in seq_len(10)) {
  candidates <- candidate_pool[cluster == cluster_id][order(stratified_q, q_value, -fold_enrichment, -hits)]
  chosen <- character()
  for (term_id in candidates$go_id) {
    redundant <- any(vapply(chosen, function(previous) {
      length(intersect(annotation_sets[[term_id]], annotation_sets[[previous]])) / length(union(annotation_sets[[term_id]], annotation_sets[[previous]])) > 0.5
    }, logical(1)))
    if (!redundant) chosen <- c(chosen, term_id)
    if (length(chosen) >= 3) break
  }
  selected <- candidates[match(chosen, go_id)]
  if (nrow(selected)) {
    selected[, representative_rank := seq_len(.N)]
    representatives[[as.character(cluster_id)]] <- selected
  }
}
representatives <- if (length(representatives)) rbindlist(representatives, fill = TRUE) else copy(bp_primary[0])
if (!"representative_rank" %in% names(representatives)) representatives[, representative_rank := integer()]
representatives[, supporting_tissues := vapply(seq_len(.N), function(row_index) {
  term_id <- go_id[row_index]
  cluster_id <- cluster[row_index]
  paste(tissue_tests[ontology == "BP" & go_id == term_id & cluster == cluster_id & q_value < 0.05 & hits >= 3, tissue], collapse = ";")
}, character(1))]
representatives[, supporting_tissue_count := lengths(strsplit(supporting_tissues, ";", fixed = TRUE))]
representatives[supporting_tissues == "", supporting_tissue_count := 0L]
broad <- pooled[context == "common_detectable_sensitivity" & ontology == "BP", .(cluster, go_id, broad_q = q_value, broad_fold = fold_enrichment)]
representatives <- merge(representatives, broad, by = c("cluster", "go_id"), all.x = TRUE)
setorder(representatives, cluster, representative_rank)
save_table(representatives, "05_summary/representative_BP_candidates.tsv")
save_table(candidate_pool, "05_summary/all_qualifying_BP_candidates.tsv")

tissue_representatives <- list()
for (cluster_id in seq_len(10)) {
  candidates <- tissue_tests[ontology == "BP" & cluster == cluster_id & q_value < 0.05 & hits >= 5 & fold_enrichment >= 1.2][order(q_value, -fold_enrichment)]
  chosen <- character()
  selected <- integer()
  for (row_index in seq_len(nrow(candidates))) {
    term_id <- candidates$go_id[row_index]
    redundant <- any(vapply(chosen, function(previous) {
      length(intersect(annotation_sets[[term_id]], annotation_sets[[previous]])) / length(union(annotation_sets[[term_id]], annotation_sets[[previous]])) > 0.5
    }, logical(1)))
    if (!redundant) {
      chosen <- c(chosen, term_id)
      selected <- c(selected, row_index)
    }
    if (length(chosen) >= 3) break
  }
  if (length(selected)) {
    retained <- candidates[selected]
    retained[, representative_rank := seq_len(.N)]
    tissue_representatives[[as.character(cluster_id)]] <- retained
  }
}
tissue_representatives <- rbindlist(tissue_representatives)
save_table(tissue_representatives, "05_summary/tissue_qualified_BP_candidates.tsv")

old_labels <- data.table(cluster = 1:10,
  old_label = c("protein metabolism", "RNA splicing", "RNA binding", "cytoskeletal binding", "RNA processing", "adenyl nucleotide binding", "RNP complex", "Cell adhesion", "ATP binding", "RNP complex"),
  mapped_term = c("protein metabolic process", "RNA splicing", "RNA binding", "cytoskeletal protein binding", "RNA processing", "adenyl nucleotide binding", "ribonucleoprotein complex", "cell adhesion", "ATP binding", "ribonucleoprotein complex"))
old_term_map <- unique(go[, .(go_id, term, ontology)])
old_audit <- merge(old_labels, old_term_map, by.x = "mapped_term", by.y = "term", all.x = TRUE)
old_audit <- merge(old_audit, pooled[context == "age_associated_primary"], by = c("cluster", "go_id", "ontology"), all.x = TRUE)
setorder(old_audit, cluster)
save_table(old_audit, "05_summary/original_function_label_audit.tsv")

summary <- primary[, .(transcript_tissue_members = .N, unique_transcripts = uniqueN(transcript_id), unique_genes = uniqueN(gene_id)), by = cluster]
summary <- merge(summary, bp_primary[, .(BP_terms_tested = .N, BP_terms_q_lt_05 = sum(q_value < 0.05 & hits >= 3),
                                         BP_terms_both_q_lt_05 = sum(q_value < 0.05 & stratified_q < 0.05 & hits >= 3 & !stratified_low_expected)), by = cluster], by = "cluster")
summary <- merge(summary, representatives[representative_rank == 1, .(cluster, candidate_function = term, go_id, hits, foreground_genes, fold_enrichment, q_value, stratified_q, broad_q, supporting_tissue_count)], by = "cluster", all.x = TRUE)
summary[is.na(candidate_function), candidate_function := "No robust pooled BP label; retain trajectory description"]
save_table(summary, "05_summary/ten_cluster_function_summary.tsv")
selected_terms <- unique(c(representatives$go_id, tissue_representatives$go_id))
save_table(tissue_tests[ontology == "BP" & go_id %in% selected_terms], "05_summary/representative_terms_all_tissues.tsv.gz")
hit_pairs <- unique(pooled[q_value < 0.05 & hits >= 3, .(cluster, go_id, ontology)])
hits <- merge(gene_clusters, go, by = "gene_id", allow.cartesian = TRUE)
hits <- merge(hits, hit_pairs, by = c("cluster", "go_id", "ontology"))
hits <- merge(hits, symbols, by = "gene_id", all.x = TRUE)
save_table(hits, "03_pooled/significant_term_hit_genes.tsv.gz")
save_table(pooled[, .(tests = .N, significant = sum(q_value < 0.05 & hits >= 3)), by = .(context, ontology)], "02_audit/test_family_sizes.tsv")
save_table(tissue_tests[, .(tests = .N, significant = sum(q_value < 0.05 & hits >= 3)), by = ontology], "02_audit/tissue_test_family_sizes.tsv")
fwrite(data.table(file = list.files(inputs, full.names = TRUE), md5 = unname(tools::md5sum(list.files(inputs, full.names = TRUE)))), file.path(root, "02_audit/input_md5.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(root, "02_audit/R_sessionInfo.txt"))
print(summary)
