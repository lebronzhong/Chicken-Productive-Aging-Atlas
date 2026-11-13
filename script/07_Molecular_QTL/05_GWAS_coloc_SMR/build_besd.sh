#!/usr/bin/env bash
# Chicken Productive Aging Atlas | 07_Molecular_QTL
# Purpose: Build besd.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

set -euo pipefail

# ============================================================
# BESD preparation
# Build tuQTL BESD files from declared inputs
#
# Input:
#   Stage39 significant-tuQTL nominal associations
#   PSI phenotype BED
#   original QTL BIM
#
# Output:
#   one BESD dataset per tissue:
#     *.besd
#     *.esi
#     *.epi
#
# Memory-safe:
#   - streaming zcat/awk
#   - GNU sort with limited RAM
#   - tissues processed serially
# ============================================================


if [[ $# -ne 4 ]]; then
    echo "Usage: $0 MANIFEST_TSV INPUT_DIR OUTPUT_DIR SMR_EXECUTABLE" >&2
    echo "Manifest columns: tissue, bim_file (tab-separated; no duplicate tissues)." >&2
    exit 2
fi
MANIFEST="$1"
INPUT="$2"
OUT="$3"
SMR="$4"
TMPROOT="${OUT}/tmp"
[[ -f "$MANIFEST" ]] || { echo "Missing manifest: $MANIFEST" >&2; exit 2; }
mapfile -t TISSUES < <(awk -F '\t' 'NR>1 {sub(/\r$/, ""); print $1}' "$MANIFEST")
[[ ${#TISSUES[@]} -gt 0 ]] || { echo "Empty manifest" >&2; exit 2; }
[[ "$(head -n1 "$MANIFEST" | tr -d '\r')" == $'tissue\tbim_file' ]] || { echo "Expected tissue and bim_file columns" >&2; exit 2; }
[[ $(printf '%s\n' "${TISSUES[@]}" | sort -u | wc -l) -eq ${#TISSUES[@]} ]] || { echo "Duplicate tissue" >&2; exit 2; }
mkdir -p "$OUT" "$TMPROOT"

# Memory control for GNU sort
#
# sort can spill to disk rather than occupying excessive RAM.
# ============================================================

SORT_MEM="512M"

export LC_ALL=C


# ============================================================
# Basic checks
# ============================================================

if [[ ! -x "${SMR}" ]]
then
    echo "ERROR: SMR executable not found or not executable:"
    echo "${SMR}"
    exit 1
fi


echo "============================================================"
echo "BESD preparation: Build tuQTL BESD"
echo "============================================================"
echo "SMR    : ${SMR}"
echo "INPUT  : ${INPUT}"
echo "OUTPUT : ${OUT}"
echo "TMP    : ${TMPROOT}"
echo "SORT RAM limit: ${SORT_MEM}"
echo


# ============================================================
# Resolve the original QTL BIM from the explicit manifest.
get_bim () {
    local T="$1"
    local VALUE
    VALUE=$(awk -F '\t' -v t="$T" 'NR>1 && $1==t {sub(/\r$/, "", $2); print $2}' "$MANIFEST")
    [[ -n "$VALUE" ]] || { echo "Missing BIM for $T" >&2; return 1; }
    if [[ "$VALUE" != /* ]]; then VALUE="$(dirname "$MANIFEST")/$VALUE"; fi
    printf '%s\n' "$VALUE"
}

# Main
# ============================================================

SUMMARY="${OUT}/00_BESD_build_summary.tsv"

echo -e \
"tissue\tphenotypes\tvariants\tmatched_variants\tbesd_bytes\tesi_bytes\tepi_bytes" \
> "${SUMMARY}"


for T in "${TISSUES[@]}"
do

    echo
    echo "======================================================================"
    echo "TISSUE: ${T}"
    echo "======================================================================"

    [[ -n "$T" && "$T" != */* && "$T" != *\\* && "$T" != . && "$T" != .. ]] || { echo "Invalid tissue label" >&2; exit 2; }
    BIM=$(get_bim "${T}")

    NOMINAL="${INPUT}/${T}/${T}.tuQTL.sig.nominal.tsv.gz"

    PHENO="${INPUT}/${T}/${T}.tuQTL.sig.PSI.bed"

    PREFIX="${OUT}/${T}.tuQTL"

    TMP="${TMPROOT}/${T}"


    # ========================================================
    # Input QC
    # ========================================================

    for F in \
        "${NOMINAL}" \
        "${PHENO}" \
        "${BIM}"
    do
        if [[ ! -s "${F}" ]]
        then
            echo "ERROR: missing or empty input:"
            echo "${F}"
            exit 1
        fi
    done


    echo "Manifest    : ${MANIFEST}"
    echo "Nominal     : ${NOMINAL}"
    echo "PSI BED     : ${PHENO}"
    echo "BIM         : ${BIM}"


    # ========================================================
    # Resume
    # ========================================================

    if [[ -s "${PREFIX}.besd" \
       && -s "${PREFIX}.esi" \
       && -s "${PREFIX}.epi" ]]
    then

        echo
        echo "[SKIP] Complete BESD already exists:"
        echo "${PREFIX}"

        continue
    fi


    # Remove incomplete old output if present
    rm -f \
        "${PREFIX}.besd" \
        "${PREFIX}.esi" \
        "${PREFIX}.epi"


    rm -rf "${TMP}"

    mkdir -p "${TMP}"


    # ========================================================
    # Verify gzip
    # ========================================================

    echo
    echo "[1/8] Check nominal gzip"

    gzip -t "${NOMINAL}"

    HEADER=$(zcat "${NOMINAL}" | head -1 || true)

    echo "Header:"
    echo "${HEADER}"
    [[ "${HEADER}" == $'phenotype_id\tvariant_id\tdistance\taf\tpval_nominal\tslope\tslope_se' ]] || {
        echo "Unexpected nominal schema; run prepare_smr_input.py first." >&2
        exit 2
    }


    # ========================================================
    # Convert Stage39 nominal format
    #
    # Stage39:
    # 1 phenotype_id
    # 2 variant_id
    # 3 distance
    # 4 af
    # 5 pval_nominal
    # 6 slope
    # 7 slope_se
    #
    # SMR fastQTL nominal:
    # phenotype variant distance p beta
    #
    # No header.
    # ========================================================

    echo
    echo "[2/8] Convert nominal -> fastQTL nominal format"

    zcat "${NOMINAL}" \
    | awk '
        BEGIN {
            OFS="\t"
        }

        NR > 1 {

            if (NF < 6) {
                print "ERROR: malformed nominal row at line " NR > "/dev/stderr"
                exit 42
            }

            print $1,$2,$3,$5,$6
        }
      ' \
    | gzip -c \
    > "${TMP}/fastqtl.nominal.gz"


    gzip -t \
        "${TMP}/fastqtl.nominal.gz"


    echo "Example:"
    zcat "${TMP}/fastqtl.nominal.gz" \
    | head -2 || true


    # ========================================================
    # Make initial BESD
    # ========================================================

    echo
    echo "[3/8] Make BESD"

    "${SMR}" \
        --eqtl-summary "${TMP}/fastqtl.nominal.gz" \
        --fastqtl-nominal-format \
        --make-besd \
        --out "${PREFIX}"


    if [[ ! -s "${PREFIX}.besd" ]]
    then
        echo "ERROR: BESD creation failed for ${T}"
        exit 1
    fi


    # ========================================================
    # Build unique SNP -> AF table
    #
    # This can contain millions of input rows, therefore:
    #   streaming extraction
    #   external sort
    #   only 512 MB sort RAM
    #
    # No pandas / Python memory usage.
    # ========================================================

    echo
    echo "[4/8] Extract unique SNP frequencies"

    zcat "${NOMINAL}" \
    | awk '
        BEGIN {
            OFS="\t"
        }

        NR > 1 {
            print $2,$4
        }
      ' \
    | sort \
        -T "${TMP}" \
        -S "${SORT_MEM}" \
        -k1,1 \
        -u \
    > "${TMP}/snp_af.tsv"


    N_VARIANTS=$(
        wc -l \
        < "${TMP}/snp_af.tsv"
    )

    echo "Unique nominal variants: ${N_VARIANTS}"


    # ========================================================
    # Normalize and sort BIM by SNP
    #
    # BIM output:
    # SNP CHR CM BP A1 A2
    # ========================================================

    echo
    echo "[5/8] Prepare BIM SNP map"

    awk '
        BEGIN {
            OFS="\t"
        }

        {
            print $2,$1,$3,$4,$5,$6
        }
      ' "${BIM}" \
    | sort \
        -T "${TMP}" \
        -S "${SORT_MEM}" \
        -k1,1 \
    > "${TMP}/bim_by_snp.tsv"


    # ========================================================
    # Check that all nominal variants occur in BIM
    # ========================================================

    cut -f1 "${TMP}/snp_af.tsv" \
    > "${TMP}/nominal_snps.txt"

    cut -f1 "${TMP}/bim_by_snp.tsv" \
    | sort \
        -T "${TMP}" \
        -S "${SORT_MEM}" \
        -u \
    > "${TMP}/bim_snps.txt"


    comm -23 \
        "${TMP}/nominal_snps.txt" \
        "${TMP}/bim_snps.txt" \
    > "${TMP}/variants_missing_from_bim.txt"


    N_MISSING=$(
        wc -l \
        < "${TMP}/variants_missing_from_bim.txt"
    )


    echo "Variants absent from BIM: ${N_MISSING}"


    if [[ "${N_MISSING}" -gt 0 ]]
    then

        echo "ERROR: nominal variants are missing from BIM."

        echo "First missing variants:"
        head "${TMP}/variants_missing_from_bim.txt"

        cp \
            "${TMP}/variants_missing_from_bim.txt" \
            "${OUT}/${T}.variants_missing_from_bim.txt"

        exit 1
    fi


    # ========================================================
    # Join AF with BIM
    #
    # Result:
    # 1 SNP
    # 2 AF
    # 3 CHR
    # 4 CM
    # 5 BP
    # 6 A1
    # 7 A2
    # ========================================================

    join \
        -t $'\t' \
        -1 1 \
        -2 1 \
        "${TMP}/snp_af.tsv" \
        "${TMP}/bim_by_snp.tsv" \
    > "${TMP}/snp_info.tsv"


    N_MATCHED=$(
        wc -l \
        < "${TMP}/snp_info.tsv"
    )


    echo "Matched variants: ${N_MATCHED}"


    if [[ "${N_MATCHED}" -ne "${N_VARIANTS}" ]]
    then

        echo "ERROR:"
        echo "Unique nominal SNP count != BIM-matched SNP count"

        exit 1
    fi


    # ========================================================
    # Build SMR ESI update file
    #
    # Format:
    # Chr SNP GD BP A1 A2 Freq
    # ========================================================

    awk '
        BEGIN {
            OFS="\t"
        }

        {
            print $3,$1,$4,$5,$6,$7,$2
        }
      ' "${TMP}/snp_info.tsv" \
    > "${TMP}/update.esi"


    # ========================================================
    # Build frequency update file
    #
    # Preserve format used by previous eQTL-SMR pipeline:
    #
    # SNP A1 A2 Freq
    # ========================================================

    awk '
        BEGIN {
            OFS="\t"
        }

        {
            print $1,$6,$7,$2
        }
      ' "${TMP}/snp_info.tsv" \
    > "${TMP}/update.freq"


    # ========================================================
    # Build transcript EPI directly from Stage39 PSI BED
    #
    # Format:
    # ProbeChr ProbeID GD ProbeBP Gene Orientation
    #
    # Strand is not available in the four-column Stage39 BED.
    # '+' is used as metadata only; SMR/HEIDI association
    # statistics do not depend on this field.
    # ========================================================

    echo
    echo "[6/8] Build transcript EPI"

    awk '
        BEGIN {
            OFS="\t"
        }

        NR > 1 {

            chr=$1

            sub(/^chr/, "", chr)

            print chr,$4,"0",$2,$4,"+"
        }
      ' "${PHENO}" \
    | sort \
        -T "${TMP}" \
        -S "${SORT_MEM}" \
        -k2,2 \
        -u \
    > "${TMP}/update.epi"


    N_PHENO=$(
        wc -l \
        < "${TMP}/update.epi"
    )

    echo "Phenotypes in EPI: ${N_PHENO}"


    # ========================================================
    # Update BESD metadata
    # ========================================================

    echo
    echo "[7/8] Update ESI / EPI / allele frequency"

    "${SMR}" \
        --beqtl-summary "${PREFIX}" \
        --update-esi "${TMP}/update.esi"

    "${SMR}" \
        --beqtl-summary "${PREFIX}" \
        --update-epi "${TMP}/update.epi"

    "${SMR}" \
        --beqtl-summary "${PREFIX}" \
        --update-freq "${TMP}/update.freq"


    # ========================================================
    # Final QC
    # ========================================================

    echo
    echo "[8/8] Final QC"

    for F in \
        "${PREFIX}.besd" \
        "${PREFIX}.esi" \
        "${PREFIX}.epi"
    do

        if [[ ! -s "${F}" ]]
        then
            echo "ERROR: missing final BESD component:"
            echo "${F}"
            exit 1
        fi

    done


    N_ESI=$(
        wc -l \
        < "${PREFIX}.esi"
    )

    N_EPI=$(
        wc -l \
        < "${PREFIX}.epi"
    )


    echo
    echo "Expected variants : ${N_VARIANTS}"
    echo "Final ESI rows    : ${N_ESI}"

    echo "Expected probes   : ${N_PHENO}"
    echo "Final EPI rows    : ${N_EPI}"


    # ESI and EPI written by SMR may include header depending on
    # version, therefore only report rather than hard-fail here.


    BESD_BYTES=$(
        stat -c%s "${PREFIX}.besd"
    )

    ESI_BYTES=$(
        stat -c%s "${PREFIX}.esi"
    )

    EPI_BYTES=$(
        stat -c%s "${PREFIX}.epi"
    )


    echo -e \
"${T}\t${N_PHENO}\t${N_VARIANTS}\t${N_MATCHED}\t${BESD_BYTES}\t${ESI_BYTES}\t${EPI_BYTES}" \
    >> "${SUMMARY}"


    echo
    echo "Completed: ${T}"

    ls -lh \
        "${PREFIX}.besd" \
        "${PREFIX}.esi" \
        "${PREFIX}.epi"


    # Remove large temporary files after successful completion
    rm -rf "${TMP}"

done


echo
echo "======================================================================"
echo "BESD preparation completed successfully"
echo "======================================================================"

echo
echo "Summary:"
column -t "${SUMMARY}" || cat "${SUMMARY}"

echo
echo "Output directory:"
echo "${OUT}"

