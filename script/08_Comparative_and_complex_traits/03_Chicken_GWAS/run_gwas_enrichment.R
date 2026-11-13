# Shared sum-of-squared-effects enrichment with an explicit trait/column manifest.
run_gwas_enrichment <- function(args) {
  if (!length(args) %in% c(4L, 5L)) stop("Usage: Rscript ENTRY TRAIT_MANIFEST INTERVAL_BED OUTPUT_DIR GENESET_HELPER [N_PERM]")
  suppressPackageStartupMessages(library(data.table))
  manifest_file <- normalizePath(args[1], mustWork = TRUE)
  traits <- fread(manifest_file, colClasses = "character")
  required <- c("trait_id", "file", "chr_col", "pos_col", "p_col", "effect_col")
  if (!all(required %in% names(traits)) || nrow(traits) < 1L || anyDuplicated(traits$trait_id) || anyNA(traits[, ..required])) stop("Invalid trait manifest")
  if (any(!grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", traits$trait_id))) stop("Trait identifiers must be safe output names")
  nperm <- if (length(args) == 5L) as.integer(args[5]) else 1000L
  if (is.na(nperm) || nperm < 1L) stop("N_PERM must be positive")
  infile <- normalizePath(args[2], mustWork = TRUE)
  outdir <- args[3]; dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  source(args[4])
  for (i in seq_len(nrow(traits))) {
    spec <- traits[i]
    path <- spec$file
    if (!grepl("^(/|[A-Za-z]:[/\\])", path)) path <- file.path(dirname(manifest_file), path)
    raw <- fread(path)
    fields <- unlist(spec[, .(chr_col, pos_col, p_col, effect_col)], use.names = FALSE)
    if (!all(fields %in% names(raw))) stop("GWAS fields absent for ", spec$trait_id)
    # Preserve the recorded SNP coordinates and squared-effect statistic.
    GWAS <- data.frame(CHR = raw[[spec$chr_col]], BP = raw[[spec$pos_col]],
                       END = raw[[spec$pos_col]], P = raw[[spec$p_col]],
                       EFFECT = raw[[spec$effect_col]])
    work <- tempfile(pattern = "gwas_", tmpdir = outdir); dir.create(work)
    temp_bed <- file.path(work, "gwas.bed"); overlap <- file.path(work, "overlap.tsv")
    write.table(GWAS, temp_bed, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    status <- system2("bedtools", c("intersect", "-a", shQuote(temp_bed), "-b", shQuote(infile), "-wa", "-wb", "-loj", "-bed"), stdout = overlap)
    if (status != 0L) stop("bedtools failed for ", spec$trait_id)
    joined <- readLines(overlap)
    writeLines(gsub("\t\t", "\t", joined), overlap)
    GWAS_MDR <- as.data.frame(fread(overlap, header = FALSE, na.strings = ".", sep = "\t", colClasses = c("character","integer","integer","numeric","numeric","character","integer","integer","character","character","character")))
    if (ncol(GWAS_MDR) != 11L) stop("Interval BED must contain six columns")
    GWAS_MDR <- GWAS_MDR[, -6]
    GWAS_MDR <- GWAS_MDR[, -3]
    colnames(GWAS_MDR) <- c("CHR","BP","P","Zscore","Start","End","State","transcript_id","Tissue")
    GWAS_MDR$SNP_ID <- paste0(GWAS_MDR$CHR,":",GWAS_MDR$BP,":",GWAS_MDR$transcript_id,":",GWAS_MDR$Tissue)
    rsIDs <- as.character(GWAS_MDR$SNP_ID)
    GENEIDs <- as.character(GWAS_MDR$State)
    tSets <- unstack(data.frame(rsIDs = rsIDs, GENEIDs = GENEIDs), rsIDs ~ GENEIDs)
    GWAS_MDR <- GWAS_MDR[!duplicated(GWAS_MDR$SNP_ID), ]
    row.names(GWAS_MDR) <- GWAS_MDR$SNP_ID
    b2 <- (GWAS_MDR$Zscore)**2
    names(b2) <- rownames(GWAS_MDR)
    setTt2 <- setTest(stat = b2, sets = tSets, nperm = nperm, method = "sum")
    setTt2$FDR <- p.adjust(setTt2$p, method = "fdr")
    write.table(setTt2, file.path(outdir, paste0(spec$trait_id, "_GWAS_enrichment.txt")), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
    unlink(work, recursive = TRUE)
  }
}
