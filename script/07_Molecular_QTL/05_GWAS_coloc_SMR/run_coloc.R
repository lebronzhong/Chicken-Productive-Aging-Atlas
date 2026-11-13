#!/usr/bin/env Rscript
# Chicken Productive Aging Atlas | 07_Molecular_QTL
# Purpose: GWAS colocalization and SMR.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

# ============================================================
# Candidate-level colocalization
# Dual colocalization for candidate trait loci
#
# Use the paths supplied on the command line.
#
# Comparison 1:
#   GWAS vs tuQTL
#
# Comparison 2:
#   GWAS vs conventional gene-eQTL
#
# Both comparisons use exactly the same SNP set.
#
# coloc version expected: 5.2.3
# ============================================================

# ============================================================
# R library
# ============================================================

.libPaths(
    c(
        "/path/to/chicken-aging/software/Rlibs",
        .libPaths()
    )
)

suppressPackageStartupMessages(
    library(coloc)
)

cat(
    "coloc version:",
    as.character(
        packageVersion("coloc")
    ),
    "\n"
)

# ============================================================
# General settings
# ============================================================

options(
    stringsAsFactors = FALSE,
    warn = 1
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: Rscript run_coloc.R MANIFEST_TSV MATCHED_DIR OUTPUT_DIR")
MANIFEST <- args[1]
MATCH_DIR <- args[2]
OUT_DIR <- args[3]
if (!file.exists(MANIFEST)) stop("Candidate manifest not found: ", MANIFEST)
if (!dir.exists(MATCH_DIR)) stop("Matched input directory not found: ", MATCH_DIR)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

manifest <- read.table(
    MANIFEST,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

required_manifest <- c(
    "candidate_id",
    "tissue",
    "transcript_id",
    "parent_gene_id",
    "trait_num",
    "trait_name",
    "N_qtl"
)

missing_manifest <- setdiff(
    required_manifest,
    colnames(manifest)
)

if (length(missing_manifest) > 0) {

    stop(
        paste(
            "Manifest missing columns:",
            paste(
                missing_manifest,
                collapse = ", "
            )
        )
    )
}

cat(
    "Candidates:",
    nrow(manifest),
    "\n"
)

cat(
    "Manifest:",
    MANIFEST,
    "\n"
)

# ============================================================
# coloc priors
#
# Same priors for GWAS-tuQTL and GWAS-gene-eQTL
# ============================================================

P1  <- 1e-4
P2  <- 1e-4
P12 <- 1e-5

# ============================================================
# Required columns in matched input
# ============================================================

required_cols <- c(

    "snp",

    "gwas_beta",
    "gwas_se",
    "gwas_p",
    "gwas_maf",
    "gwas_N",

    "tu_beta",
    "tu_se",
    "tu_p",
    "tu_maf",

    "e_beta",
    "e_se",
    "e_p",
    "e_maf",

    "qtl_N"
)

# ============================================================
# Helper functions
# ============================================================

get_pp <- function(
    summary_vector,
    name
) {

    x <- summary_vector[name]

    if (length(x) == 0) {
        return(NA_real_)
    }

    as.numeric(x)
}

conditional_h4 <- function(
    h3,
    h4
) {

    denom <- h3 + h4

    if (
        is.na(denom) ||
        denom <= 0
    ) {
        return(NA_real_)
    }

    h4 / denom
}

# ============================================================
# Main
# ============================================================

summary_list <- list()

k <- 1

for (i in seq_len(nrow(manifest))) {

    meta <- manifest[i, ]

    CID <- meta$candidate_id

    infile <- file.path(
        MATCH_DIR,
        paste0(
            CID,
            ".matched.tsv.gz"
        )
    )

    cat(
        "\n",
        paste(
            rep("=", 100),
            collapse = ""
        ),
        "\n",
        sep = ""
    )

    cat(
        CID,
        "|",
        meta$tissue,
        "|",
        meta$trait_name,
        "|",
        meta$transcript_id,
        "|",
        meta$parent_gene_id,
        "\n"
    )

    cat(
        paste(
            rep("=", 100),
            collapse = ""
        ),
        "\n"
    )

    if (!file.exists(infile)) {

        stop(
            paste(
                "Matched file not found:",
                infile
            )
        )
    }

    # ========================================================
    # Read matched data
    # ========================================================

    dat <- read.table(
        gzfile(infile),
        header = TRUE,
        sep = "\t",
        quote = "",
        comment.char = "",
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

    missing_cols <- setdiff(
        required_cols,
        colnames(dat)
    )

    if (length(missing_cols) > 0) {

        stop(
            paste(
                CID,
                "missing columns:",
                paste(
                    missing_cols,
                    collapse = ", "
                )
            )
        )
    }

    cat(
        "Raw matched SNPs:",
        nrow(dat),
        "\n"
    )

    # ========================================================
    # Convert numeric columns
    # ========================================================

    numeric_cols <- setdiff(
        required_cols,
        "snp"
    )

    for (x in numeric_cols) {

        dat[[x]] <- suppressWarnings(
            as.numeric(
                dat[[x]]
            )
        )
    }

    # ========================================================
    # Complete-case QC
    #
    # IMPORTANT:
    # Keep exactly the same SNPs for BOTH coloc analyses.
    # ========================================================

    good <- (

        !is.na(dat$snp) &

        is.finite(dat$gwas_beta) &
        is.finite(dat$gwas_se) &

        is.finite(dat$tu_beta) &
        is.finite(dat$tu_se) &

        is.finite(dat$e_beta) &
        is.finite(dat$e_se) &

        is.finite(dat$gwas_maf) &
        is.finite(dat$tu_maf) &
        is.finite(dat$e_maf) &

        dat$gwas_se > 0 &
        dat$tu_se > 0 &
        dat$e_se > 0 &

        dat$gwas_maf > 0 &
        dat$gwas_maf < 0.5 &

        dat$tu_maf > 0 &
        dat$tu_maf < 0.5 &

        dat$e_maf > 0 &
        dat$e_maf < 0.5
    )

    dat <- dat[
        good,
        ,
        drop = FALSE
    ]

    # Remove duplicated SNP IDs if any
    dat <- dat[
        !duplicated(dat$snp),
        ,
        drop = FALSE
    ]

    cat(
        "QC-passed SNPs:",
        nrow(dat),
        "\n"
    )

    if (nrow(dat) < 50) {

        warning(
            paste(
                CID,
                "has only",
                nrow(dat),
                "SNPs after QC."
            )
        )
    }

    if (nrow(dat) < 10) {

        stop(
            paste(
                CID,
                "has too few SNPs for coloc."
            )
        )
    }

    # ========================================================
    # Sample size
    # ========================================================

    GWAS_N <- median(
        dat$gwas_N[
            is.finite(dat$gwas_N)
        ],
        na.rm = TRUE
    )

    QTL_N <- as.numeric(
        meta$N_qtl
    )

    if (!is.finite(GWAS_N)) {

        stop(
            paste(
                CID,
                "GWAS N is invalid."
            )
        )
    }

    if (!is.finite(QTL_N)) {

        stop(
            paste(
                CID,
                "QTL N is invalid."
            )
        )
    }

    cat(
        "GWAS N:",
        GWAS_N,
        "\n"
    )

    cat(
        "QTL N :",
        QTL_N,
        "\n"
    )

    # ========================================================
    # Dataset 1: GWAS
    # ========================================================

    D_GWAS <- list(

        beta =
            dat$gwas_beta,

        varbeta =
            dat$gwas_se^2,

        snp =
            dat$snp,

        MAF =
            dat$gwas_maf,

        N =
            GWAS_N,

        type =
            "quant"
    )

    # ========================================================
    # Dataset 2A: tuQTL
    # ========================================================

    D_TU <- list(

        beta =
            dat$tu_beta,

        varbeta =
            dat$tu_se^2,

        snp =
            dat$snp,

        MAF =
            dat$tu_maf,

        N =
            QTL_N,

        type =
            "quant"
    )

    # ========================================================
    # Dataset 2B: conventional gene-eQTL
    # ========================================================

    D_E <- list(

        beta =
            dat$e_beta,

        varbeta =
            dat$e_se^2,

        snp =
            dat$snp,

        MAF =
            dat$e_maf,

        N =
            QTL_N,

        type =
            "quant"
    )

    # ========================================================
    # COLOC 1
    # GWAS vs tuQTL
    # ========================================================

    cat(
        "\nRunning GWAS-tuQTL coloc...\n"
    )

    RES_TU <- coloc.abf(

        dataset1 = D_GWAS,
        dataset2 = D_TU,

        p1 = P1,
        p2 = P2,
        p12 = P12
    )

    # ========================================================
    # COLOC 2
    # GWAS vs gene-eQTL
    # ========================================================

    cat(
        "Running GWAS-gene-eQTL coloc...\n"
    )

    RES_E <- coloc.abf(

        dataset1 = D_GWAS,
        dataset2 = D_E,

        p1 = P1,
        p2 = P2,
        p12 = P12
    )

    # ========================================================
    # Extract posterior probabilities
    # ========================================================

    TU_H0 <- get_pp(
        RES_TU$summary,
        "PP.H0.abf"
    )

    TU_H1 <- get_pp(
        RES_TU$summary,
        "PP.H1.abf"
    )

    TU_H2 <- get_pp(
        RES_TU$summary,
        "PP.H2.abf"
    )

    TU_H3 <- get_pp(
        RES_TU$summary,
        "PP.H3.abf"
    )

    TU_H4 <- get_pp(
        RES_TU$summary,
        "PP.H4.abf"
    )

    E_H0 <- get_pp(
        RES_E$summary,
        "PP.H0.abf"
    )

    E_H1 <- get_pp(
        RES_E$summary,
        "PP.H1.abf"
    )

    E_H2 <- get_pp(
        RES_E$summary,
        "PP.H2.abf"
    )

    E_H3 <- get_pp(
        RES_E$summary,
        "PP.H3.abf"
    )

    E_H4 <- get_pp(
        RES_E$summary,
        "PP.H4.abf"
    )

    TU_COND_H4 <- conditional_h4(
        TU_H3,
        TU_H4
    )

    E_COND_H4 <- conditional_h4(
        E_H3,
        E_H4
    )

    # ========================================================
    # Top shared SNP posterior
    # ========================================================

    TU_TOP_INDEX <- which.max(
        RES_TU$results$SNP.PP.H4
    )

    TU_TOP_SNP <- (
        RES_TU$results$snp[
            TU_TOP_INDEX
        ]
    )

    TU_TOP_PP <- (
        RES_TU$results$SNP.PP.H4[
            TU_TOP_INDEX
        ]
    )

    E_TOP_INDEX <- which.max(
        RES_E$results$SNP.PP.H4
    )

    E_TOP_SNP <- (
        RES_E$results$snp[
            E_TOP_INDEX
        ]
    )

    E_TOP_PP <- (
        RES_E$results$SNP.PP.H4[
            E_TOP_INDEX
        ]
    )

    # ========================================================
    # Save SNP-level posterior
    # ========================================================

    TU_SNP_FILE <- file.path(
        OUT_DIR,
        paste0(
            CID,
            ".GWAS_tuQTL.snp_posteriors.tsv.gz"
        )
    )

    E_SNP_FILE <- file.path(
        OUT_DIR,
        paste0(
            CID,
            ".GWAS_gene_eQTL.snp_posteriors.tsv.gz"
        )
    )

    write.table(

        RES_TU$results,

        gzfile(
            TU_SNP_FILE
        ),

        quote = FALSE,
        sep = "\t",
        row.names = FALSE
    )

    write.table(

        RES_E$results,

        gzfile(
            E_SNP_FILE
        ),

        quote = FALSE,
        sep = "\t",
        row.names = FALSE
    )

    # ========================================================
    # Candidate-level summary
    # ========================================================

    summary_list[[k]] <- data.frame(

        candidate_id =
            CID,

        tissue =
            meta$tissue,

        trait_num =
            meta$trait_num,

        trait_name =
            meta$trait_name,

        transcript_id =
            meta$transcript_id,

        parent_gene_id =
            meta$parent_gene_id,

        nsnps =
            nrow(dat),

        GWAS_N =
            GWAS_N,

        QTL_N =
            QTL_N,

        # GWAS-tuQTL
        GWAS_tu_H0 =
            TU_H0,

        GWAS_tu_H1 =
            TU_H1,

        GWAS_tu_H2 =
            TU_H2,

        GWAS_tu_H3 =
            TU_H3,

        GWAS_tu_H4 =
            TU_H4,

        GWAS_tu_H4_conditional =
            TU_COND_H4,

        GWAS_tu_top_coloc_SNP =
            TU_TOP_SNP,

        GWAS_tu_top_SNP_PP_H4 =
            TU_TOP_PP,

        # GWAS-gene-eQTL
        GWAS_e_H0 =
            E_H0,

        GWAS_e_H1 =
            E_H1,

        GWAS_e_H2 =
            E_H2,

        GWAS_e_H3 =
            E_H3,

        GWAS_e_H4 =
            E_H4,

        GWAS_e_H4_conditional =
            E_COND_H4,

        GWAS_e_top_coloc_SNP =
            E_TOP_SNP,

        GWAS_e_top_SNP_PP_H4 =
            E_TOP_PP,

        # Comparison
        delta_H4 =
            TU_H4 - E_H4,

        delta_conditional_H4 =
            TU_COND_H4 -
            E_COND_H4,

        stringsAsFactors = FALSE
    )

    cat(
        "\n"
    )

    cat(
        "GWAS-tuQTL:\n"
    )

    cat(
        "  H3 =",
        format(
            TU_H3,
            digits = 6
        ),
        "\n"
    )

    cat(
        "  H4 =",
        format(
            TU_H4,
            digits = 6
        ),
        "\n"
    )

    cat(
        "  H4/(H3+H4) =",
        format(
            TU_COND_H4,
            digits = 6
        ),
        "\n"
    )

    cat(
        "  top SNP =",
        TU_TOP_SNP,
        "PP=",
        format(
            TU_TOP_PP,
            digits = 6
        ),
        "\n"
    )

    cat(
        "\nGWAS-gene-eQTL:\n"
    )

    cat(
        "  H3 =",
        format(
            E_H3,
            digits = 6
        ),
        "\n"
    )

    cat(
        "  H4 =",
        format(
            E_H4,
            digits = 6
        ),
        "\n"
    )

    cat(
        "  H4/(H3+H4) =",
        format(
            E_COND_H4,
            digits = 6
        ),
        "\n"
    )

    cat(
        "  top SNP =",
        E_TOP_SNP,
        "PP=",
        format(
            E_TOP_PP,
            digits = 6
        ),
        "\n"
    )

    cat(
        "\nDelta H4 =",
        format(
            TU_H4 - E_H4,
            digits = 6
        ),
        "\n"
    )

    k <- k + 1
}

# ============================================================
# Combine summary
# ============================================================

summary_df <- do.call(
    rbind,
    summary_list
)

# ============================================================
# Ranking
#
# Main objective:
# high GWAS-tuQTL H4
# low GWAS-gene-eQTL H4
#
# First rank by tu H4,
# second by delta H4.
# ============================================================

summary_df <- summary_df[
    order(
        -summary_df$GWAS_tu_H4,
        -summary_df$delta_H4
    ),
    ,
    drop = FALSE
]

SUMMARY_FILE <- file.path(
    OUT_DIR,
    "42C_coloc_summary.tsv"
)

write.table(

    summary_df,

    file = SUMMARY_FILE,

    quote = FALSE,
    sep = "\t",
    row.names = FALSE
)

# ============================================================
# Compact summary for inspection
# ============================================================

compact <- summary_df[
    ,
    c(
        "candidate_id",
        "tissue",
        "trait_name",
        "transcript_id",
        "parent_gene_id",
        "nsnps",

        "GWAS_tu_H3",
        "GWAS_tu_H4",
        "GWAS_tu_H4_conditional",

        "GWAS_e_H3",
        "GWAS_e_H4",
        "GWAS_e_H4_conditional",

        "delta_H4",

        "GWAS_tu_top_coloc_SNP",
        "GWAS_tu_top_SNP_PP_H4"
    )
]

COMPACT_FILE <- file.path(
    OUT_DIR,
    "42C_coloc_summary_compact.tsv"
)

write.table(

    compact,

    file = COMPACT_FILE,

    quote = FALSE,
    sep = "\t",
    row.names = FALSE
)

# ============================================================
# Final output
# ============================================================

cat(
    "\n",
    paste(
        rep("=", 110),
        collapse = ""
    ),
    "\n",
    sep = ""
)

cat(
    "FINAL COLOC SUMMARY\n"
)

cat(
    paste(
        rep("=", 110),
        collapse = ""
    ),
    "\n"
)

print(
    compact,
    row.names = FALSE
)

cat(
    "\nFull output:\n",
    SUMMARY_FILE,
    "\n"
)

cat(
    "\nCompact output:\n",
    COMPACT_FILE,
    "\n"
)

cat(
    "\nColocalization completed successfully.\n"
)
