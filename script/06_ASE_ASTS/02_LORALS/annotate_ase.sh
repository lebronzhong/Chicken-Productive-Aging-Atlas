#!/usr/bin/env bash
# Chicken Productive Aging Atlas | 06_ASE_ASTS
# Purpose: Annotate ase.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

set -Eeuo pipefail

# 1
work=/path/to/chicken-aging/allelic
ase=${work}/ase
scripts=${work}/scripts
software=${work}/software
logs=${work}/logs

VCFDIR=${ase}/vcf
REFDIR=${software}/lorals_ref

MANIFEST=${scripts}/lorals_manifest_matched.tsv

GENE_BED=${REFDIR}/gene_exons.bed
NO_FILTER=${REFDIR}/no_filter.bed

POSTPROCESS=${scripts}/postprocess_ase.py

ENV=/path/to/chicken-aging/software/lorals_env

PYTHON=${ENV}/bin/python
BCFTOOLS=${ENV}/bin/bcftools
ANNOTATE_ASE=${ENV}/bin/annotate_ase

export PATH=${ENV}/bin:${PATH}
export LD_LIBRARY_PATH=${ENV}/lib:${LD_LIBRARY_PATH:-}

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export PYTHONNOUSERSITE=1

BLOCK_SIZE=10
MAX_PARALLEL=10
MIN_TOTAL_COUNT=20

mkdir -p \
    "${ase}/samples" \
    "${software}/tmp" \
    "${logs}"

cd "${work}"

# 2
for TOOL in \
    "${PYTHON}" \
    "${BCFTOOLS}" \
    "${ANNOTATE_ASE}"; do

    if [ ! -x "${TOOL}" ]; then
        echo "ERROR: executable not found: ${TOOL}"
        exit 1
    fi
done

for FILE in \
    "${MANIFEST}" \
    "${GENE_BED}" \
    "${NO_FILTER}" \
    "${POSTPROCESS}"; do

    if [ ! -s "${FILE}" ]; then
        echo "ERROR: required file missing: ${FILE}"
        exit 1
    fi
done

# 3
SITE_PACKAGES=$(
    "${PYTHON}" -c \
        'import site; print(site.getsitepackages()[0])'
)

LORALS_SOURCE=${SITE_PACKAGES}/lorals

if [ ! -d "${LORALS_SOURCE}" ]; then
    echo "ERROR: installed LORALS package not found"
    echo "${LORALS_SOURCE}"
    exit 1
fi

LORALS_PATCH_ROOT=${software}/tmp/lorals_annot_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
LORALS_ASTS=${LORALS_PATCH_ROOT}/lorals/asts.py

rm -rf "${LORALS_PATCH_ROOT}"
mkdir -p "${LORALS_PATCH_ROOT}"

cp -a \
    "${LORALS_SOURCE}" \
    "${LORALS_PATCH_ROOT}/lorals"

if [ ! -s "${LORALS_ASTS}" ]; then
    echo "ERROR: copied LORALS asts.py not found"
    echo "${LORALS_ASTS}"
    exit 1
fi

if ! grep -q \
    'statistic_location=utils.nan' \
    "${LORALS_ASTS}"; then

    sed -i \
        's/stats\.stats\.KstestResult(statistic=utils\.nan, pvalue=utils\.nan)/stats.stats.KstestResult(statistic=utils.nan, pvalue=utils.nan, statistic_location=utils.nan, statistic_sign=utils.nan)/g' \
        "${LORALS_ASTS}"
fi

if ! grep -q \
    'statistic_location=utils.nan' \
    "${LORALS_ASTS}"; then

    echo "ERROR: LORALS compatibility patch failed"

    grep -n \
        "KstestResult" \
        "${LORALS_ASTS}" \
        || true

    exit 1
fi

export PYTHONPATH="${LORALS_PATCH_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
export LORALS_PATCH_ROOT

# 4
echo "============================================================"
echo "Test patched LORALS"
echo "============================================================"

"${PYTHON}" -c '
import os
import lorals
from lorals.scripts import annotate_ase

expected = os.path.realpath(
    os.environ["LORALS_PATCH_ROOT"]
)

actual = os.path.realpath(
    lorals.__file__
)

print(
    "lorals loaded from =",
    actual
)

if not actual.startswith(
    expected + os.sep
):
    raise SystemExit(
        "ERROR: wrong LORALS package loaded: " +
        actual
    )

print(
    "annotate_ase import OK"
)
'

# 5
N_TOTAL=$(
    awk '
        NR > 1 {
            n++
        }

        END {
            print n + 0
        }
    ' "${MANIFEST}"
)

if ! [[ "${N_TOTAL}" =~ ^[0-9]+$ ]]; then
    echo "ERROR: invalid manifest sample count"
    echo "N_TOTAL = ${N_TOTAL}"
    exit 1
fi

if [ "${N_TOTAL}" -ne 450 ]; then
    echo "ERROR: expected 450 matched ONT samples"
    echo "Found = ${N_TOTAL}"
    exit 1
fi

START=$(( (SLURM_ARRAY_TASK_ID - 1) * BLOCK_SIZE + 1 ))
END=$(( SLURM_ARRAY_TASK_ID * BLOCK_SIZE ))

if [ "${END}" -gt "${N_TOTAL}" ]; then
    END=${N_TOTAL}
fi

if [ "${START}" -gt "${N_TOTAL}" ]; then
    echo "No samples assigned to array task ${SLURM_ARRAY_TASK_ID}"

    rm -rf "${LORALS_PATCH_ROOT}"

    exit 0
fi

ROW_LIST=${software}/tmp/ASE_annot_rows_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.txt

seq \
    "${START}" \
    "${END}" \
    > "${ROW_LIST}"

N_TASK=$(
    wc -l < "${ROW_LIST}" \
    | tr -d '[:space:]'
)

echo "============================================================"
echo "08_2 annotate_ase"
echo "Pre-filter = totalCount > ${MIN_TOTAL_COUNT}"
echo "Total samples = ${N_TOTAL}"
echo "Array task = ${SLURM_ARRAY_TASK_ID}"
echo "Rows = ${START}-${END}"
echo "Samples in this array = ${N_TASK}"
echo "xargs parallel = ${MAX_PARALLEL}"
echo "TIME = $(date)"
echo "============================================================"

# 6
process_one()
{
    ROW_ID=$1

    LINE=$(
        sed -n \
            "$((ROW_ID + 1))p" \
            "${MANIFEST}"
    )

    if [ -z "${LINE}" ]; then
        echo "ERROR: manifest row ${ROW_ID} missing"
        return 1
    fi

    IFS=$'\t' read -r \
        array_id \
        dataset \
        age \
        tissue \
        original_tissue \
        sample \
        replicate \
        vcf_sample \
        bam \
        <<< "${LINE}"

    if [ -z "${sample}" ] || \
       [ -z "${dataset}" ] || \
       [ -z "${tissue}" ] || \
       [ -z "${vcf_sample}" ]; then

        echo "ERROR: invalid manifest row ${ROW_ID}"
        echo "${LINE}"
        return 1
    fi

    SAMPLE_VCF_GZ=${VCFDIR}/${vcf_sample}.het.biallelic.vcf.gz

    ASE_OUT=${ase}/samples/${tissue}/${dataset}/${sample}

    ASE_RAW=${ASE_OUT}/${sample}.ase.raw.tsv
    ASE_RAW_COUNT20=${ASE_OUT}/${sample}.ase.raw.totalCount_gt20.tsv
    ASE_ANNOTATED=${ASE_OUT}/${sample}.ase.annotated.tsv
    ASE_QC=${ASE_OUT}/${sample}.ase.qc.tsv
    ASE_SIG=${ASE_OUT}/${sample}.ase.FDR005.tsv
    ASE_FOR_ASTS=${ASE_OUT}/${sample}.ase.for_asts.tsv
    DONE=${ASE_OUT}/${sample}.ASE_ANNOT.SUCCESS

    mkdir -p "${ASE_OUT}"

    echo "============================================================"
    echo "START annotate sample = ${sample}"
    echo "ROW = ${ROW_ID}"
    echo "dataset = ${dataset}"
    echo "age = ${age}"
    echo "tissue = ${tissue}"
    echo "vcf_sample = ${vcf_sample}"
    echo "TIME = $(date)"
    echo "============================================================"

    # 7
    if [ ! -s "${ASE_RAW}" ]; then
        echo "ERROR [${sample}]: raw ASE missing"
        echo "Run 08_1 first"
        echo "${ASE_RAW}"
        return 1
    fi

    if [ ! -s "${SAMPLE_VCF_GZ}" ]; then
        echo "ERROR [${sample}]: VCF missing"
        echo "${SAMPLE_VCF_GZ}"
        return 1
    fi

    if [ ! -s "${SAMPLE_VCF_GZ}.tbi" ]; then
        echo "ERROR [${sample}]: VCF index missing"
        echo "${SAMPLE_VCF_GZ}.tbi"
        return 1
    fi

    # 8
    FILTER_MARK=$(
        awk -F '\t' \
            '$1 == "totalCount_filter" {print $2; exit}' \
            "${DONE}" \
            2>/dev/null \
        || true
    )

    if [ -s "${DONE}" ] && \
       [ "${FILTER_MARK}" = ">20" ] && \
       [ -s "${ASE_ANNOTATED}" ] && \
       [ -s "${ASE_QC}" ] && \
       [ -s "${ASE_FOR_ASTS}" ]; then

        echo "[${sample}] SKIP annotate_ase: completed with totalCount > 20"

        return 0
    fi

    # 9
    RAW_LINES=$(
        wc -l < "${ASE_RAW}" \
        | tr -d '[:space:]'
    )

    if ! [[ "${RAW_LINES}" =~ ^[0-9]+$ ]]; then
        echo "ERROR [${sample}]: invalid raw ASE line count"
        echo "RAW_LINES = ${RAW_LINES}"
        return 1
    fi

    if [ "${RAW_LINES}" -gt 0 ]; then
        RAW_N=$(( RAW_LINES - 1 ))
    else
        RAW_N=0
    fi

    # 10
    TOTAL_COL=$(
        awk -F '\t' '
            NR == 1 {
                for (i = 1; i <= NF; i++) {
                    gsub(/\r/, "", $i)

                    if ($i == "totalCount") {
                        print i
                        exit
                    }
                }
            }
        ' "${ASE_RAW}" \
        | tr -d '[:space:]'
    )

    if ! [[ "${TOTAL_COL}" =~ ^[0-9]+$ ]]; then
        echo "ERROR [${sample}]: totalCount column not found"
        echo "Header:"
        head -n 1 "${ASE_RAW}"
        return 1
    fi

    # 11
    TMP_FILTER=${ASE_OUT}/${sample}.ase.raw.totalCount_gt20.tmp.${SLURM_ARRAY_JOB_ID}.${SLURM_ARRAY_TASK_ID}.${ROW_ID}

    rm -f "${TMP_FILTER}"

    awk \
        -F '\t' \
        -v OFS='\t' \
        -v col="${TOTAL_COL}" \
        -v cutoff="${MIN_TOTAL_COUNT}" \
        '
            NR == 1 {
                print
                next
            }

            (($col + 0) > cutoff) {
                print
            }
        ' \
        "${ASE_RAW}" \
        > "${TMP_FILTER}"

    if [ ! -s "${TMP_FILTER}" ]; then
        echo "ERROR [${sample}]: totalCount filtering failed"

        rm -f "${TMP_FILTER}"

        return 1
    fi

    mv -f \
        "${TMP_FILTER}" \
        "${ASE_RAW_COUNT20}"

    # 12
    FILTERED_LINES=$(
        wc -l < "${ASE_RAW_COUNT20}" \
        | tr -d '[:space:]'
    )

    if ! [[ "${FILTERED_LINES}" =~ ^[0-9]+$ ]]; then
        echo "ERROR [${sample}]: invalid filtered ASE line count"
        echo "FILTERED_LINES = ${FILTERED_LINES}"
        return 1
    fi

    if [ "${FILTERED_LINES}" -gt 0 ]; then
        FILTERED_N=$(( FILTERED_LINES - 1 ))
    else
        FILTERED_N=0
    fi

    echo "[${sample}] raw ASE variants = ${RAW_N}"
    echo "[${sample}] totalCount > 20 variants = ${FILTERED_N}"

    if [ "${RAW_N}" -gt 0 ]; then

        RETAINED_FRACTION=$(
            awk \
                -v filtered="${FILTERED_N}" \
                -v raw="${RAW_N}" \
                'BEGIN {printf "%.6f", filtered / raw}'
        )

        echo "[${sample}] retained fraction = ${RETAINED_FRACTION}"
    fi

    if [ "${FILTERED_N}" -le 0 ]; then
        echo "ERROR [${sample}]: no variants remain after totalCount > 20"
        return 1
    fi

    # 13
    TMP=${software}/tmp/ASE_annot_${sample}_${ROW_ID}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

    rm -rf "${TMP}"

    mkdir -p "${TMP}"

    SAMPLE_VCF=${TMP}/${sample}.vcf

    "${BCFTOOLS}" view \
        -Ov \
        -o "${SAMPLE_VCF}" \
        "${SAMPLE_VCF_GZ}"

    if [ ! -s "${SAMPLE_VCF}" ]; then
        echo "ERROR [${sample}]: VCF conversion failed"

        rm -rf "${TMP}"

        return 1
    fi

    # 14
    rm -f \
        "${ASE_ANNOTATED}" \
        "${ASE_QC}" \
        "${ASE_SIG}" \
        "${ASE_FOR_ASTS}" \
        "${DONE}"

    "${ANNOTATE_ASE}" \
        -i "${ASE_RAW_COUNT20}" \
        -b "${GENE_BED}" \
        -f "${SAMPLE_VCF}" \
        --blacklist "${NO_FILTER}" \
        --mapping "${NO_FILTER}" \
        -c 20 \
        -n 0.5 \
        -o "${ASE_ANNOTATED}" \
        -v error

    if [ ! -s "${ASE_ANNOTATED}" ]; then
        echo "ERROR [${sample}]: annotated ASE missing"

        rm -rf "${TMP}"

        return 1
    fi

    # 15
    "${PYTHON}" \
        "${POSTPROCESS}" \
        "${ASE_ANNOTATED}" \
        "${ASE_QC}" \
        "${ASE_SIG}" \
        "${ASE_FOR_ASTS}"

    if [ ! -s "${ASE_QC}" ]; then
        echo "ERROR [${sample}]: ASE QC output missing"

        rm -rf "${TMP}"

        return 1
    fi

    if [ ! -s "${ASE_FOR_ASTS}" ]; then
        echo "ERROR [${sample}]: ASE ASTS input missing"

        rm -rf "${TMP}"

        return 1
    fi

    # 16
    {
        echo -e "sample\t${sample}"
        echo -e "dataset\t${dataset}"
        echo -e "age\t${age}"
        echo -e "tissue\t${tissue}"
        echo -e "vcf_sample\t${vcf_sample}"
        echo -e "totalCount_filter\t>20"
        echo -e "raw_variants\t${RAW_N}"
        echo -e "filtered_variants\t${FILTERED_N}"
        echo -e "status\tSUCCESS"
        echo -e "finished\t$(date)"
    } > "${DONE}"

    rm -rf "${TMP}"

    echo "============================================================"
    echo "FINISHED annotate sample = ${sample}"
    echo "raw variants = ${RAW_N}"
    echo "totalCount > 20 = ${FILTERED_N}"
    echo "TIME = $(date)"
    echo "============================================================"
}

# 17
export -f process_one

export MANIFEST
export VCFDIR

export ase
export software

export GENE_BED
export NO_FILTER
export POSTPROCESS

export PYTHON
export BCFTOOLS
export ANNOTATE_ASE

export MIN_TOTAL_COUNT

export PATH
export PYTHONPATH
export LD_LIBRARY_PATH

export SLURM_ARRAY_JOB_ID
export SLURM_ARRAY_TASK_ID

# 18
echo "============================================================"
echo "Start xargs parallel processing"
echo "Array task = ${SLURM_ARRAY_TASK_ID}"
echo "Rows = ${START}-${END}"
echo "Samples = ${N_TASK}"
echo "MAX_PARALLEL = ${MAX_PARALLEL}"
echo "============================================================"

set +e

cat "${ROW_LIST}" \
    | xargs \
        -r \
        -n 1 \
        -P "${MAX_PARALLEL}" \
        bash -c \
        'set -euo pipefail; process_one "$1"' _

XARGS_STATUS=$?

set -e

# 19
FAILED=0

echo "============================================================"
echo "Validate current array outputs"
echo "============================================================"

while read -r ROW_ID; do

    [ -z "${ROW_ID}" ] && continue

    LINE=$(
        sed -n \
            "$((ROW_ID + 1))p" \
            "${MANIFEST}"
    )

    if [ -z "${LINE}" ]; then
        echo "FAILED: manifest row ${ROW_ID} missing"

        FAILED=$((FAILED + 1))

        continue
    fi

    IFS=$'\t' read -r \
        array_id \
        dataset \
        age \
        tissue \
        original_tissue \
        sample \
        replicate \
        vcf_sample \
        bam \
        <<< "${LINE}"

    ASE_OUT=${ase}/samples/${tissue}/${dataset}/${sample}

    if [ ! -s "${ASE_OUT}/${sample}.ase.for_asts.tsv" ] || \
       [ ! -s "${ASE_OUT}/${sample}.ASE_ANNOT.SUCCESS" ]; then

        echo "FAILED: ${sample} ASE annotation outputs missing"

        FAILED=$((FAILED + 1))

    else

        echo "SUCCESS: ${sample} ASE annotation outputs exist"
    fi

done < "${ROW_LIST}"

# 20
rm -f "${ROW_LIST}"

rm -rf "${LORALS_PATCH_ROOT}"

echo "============================================================"
echo "08_2 array summary"
echo "============================================================"
echo "Array task = ${SLURM_ARRAY_TASK_ID}"
echo "Rows = ${START}-${END}"
echo "Samples = ${N_TASK}"
echo "xargs status = ${XARGS_STATUS}"
echo "Failed samples = ${FAILED}"
echo "TIME = $(date)"
echo "============================================================"

if [ "${XARGS_STATUS}" -ne 0 ] || \
   [ "${FAILED}" -gt 0 ]; then

    echo "ERROR: 08_2 array ${SLURM_ARRAY_TASK_ID} failed"

    exit 1
fi

echo "08_2 array ${SLURM_ARRAY_TASK_ID} completed successfully"
