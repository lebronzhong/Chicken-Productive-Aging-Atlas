#!/usr/bin/env bash
# Chicken Productive Aging Atlas | 06_ASE_ASTS
# Purpose: ASE and ASTS analysis.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

set -Eeuo pipefail

############################################################
# 08-3b
# genome BAM + ASE_FOR_ASTS + transcriptome BAM
# -> calc_asts -m quant
############################################################

work=/path/to/chicken-aging/allelic
ase=${work}/ase
asts=${work}/asts
scripts=${work}/scripts
software=${work}/software
logs=${work}/logs

MANIFEST=${scripts}/lorals_manifest_matched.tsv

ENV=/path/to/chicken-aging/software/lorals_env
PYTHON=${ENV}/bin/python
SAMTOOLS=${ENV}/bin/samtools
CALC_ASTS=${ENV}/bin/calc_asts

export PATH=${ENV}/bin:${PATH}
export LD_LIBRARY_PATH=${ENV}/lib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export PYTHONNOUSERSITE=1

SAMPLES_PER_TASK=2
MAX_PARALLEL=2
THREADS_PER_SAMPLE=8

: "${BATCH_START:?BATCH_START is required}"
: "${BATCH_END:?BATCH_END is required}"
: "${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"

mkdir -p \
    "${asts}/samples" \
    "${software}/tmp" \
    "${logs}"

cd "${work}"

for TOOL in \
    "${PYTHON}" \
    "${SAMTOOLS}" \
    "${CALC_ASTS}"; do

    if [ ! -x "${TOOL}" ]; then
        echo "ERROR: executable not found: ${TOOL}"
        exit 1
    fi
done

if [ ! -s "${MANIFEST}" ]; then
    echo "ERROR: manifest missing: ${MANIFEST}"
    exit 1
fi

MANIFEST_LINES=$(
    wc -l < "${MANIFEST}" \
    | tr -d '[:space:]'
)
N_TOTAL=$(( MANIFEST_LINES - 1 ))

if [ "${N_TOTAL}" -ne 450 ]; then
    echo "ERROR: expected 450 matched ONT samples; found ${N_TOTAL}"
    exit 1
fi

if [ "${BATCH_START}" -lt 1 ] || \
   [ "${BATCH_END}" -gt "${N_TOTAL}" ] || \
   [ "${BATCH_START}" -gt "${BATCH_END}" ]; then

    echo "ERROR: invalid batch range ${BATCH_START}-${BATCH_END}"
    exit 1
fi

TMPROOT=${SLURM_TMPDIR:-/tmp}
if [ ! -d "${TMPROOT}" ] || [ ! -w "${TMPROOT}" ]; then
    TMPROOT=${software}/tmp
fi

JOB_TMP=${TMPROOT}/ASTS_3b_${SLURM_ARRAY_JOB_ID:-$$}_${SLURM_ARRAY_TASK_ID}
rm -rf "${JOB_TMP}"
mkdir -p "${JOB_TMP}"

cleanup_job_tmp()
{
    rm -rf "${JOB_TMP}"
}
trap cleanup_job_tmp EXIT

############################################################
# Private LORALS copy + SciPy compatibility patch
############################################################

SITE_PACKAGES=$(
    "${PYTHON}" -c \
        'import site; print(site.getsitepackages()[0])'
)
LORALS_SOURCE=${SITE_PACKAGES}/lorals
LORALS_PATCH_ROOT=${JOB_TMP}/lorals_patch
LORALS_ASTS=${LORALS_PATCH_ROOT}/lorals/asts.py

if [ ! -d "${LORALS_SOURCE}" ]; then
    echo "ERROR: installed LORALS package not found: ${LORALS_SOURCE}"
    exit 1
fi

mkdir -p "${LORALS_PATCH_ROOT}"
cp -a \
    "${LORALS_SOURCE}" \
    "${LORALS_PATCH_ROOT}/lorals"

if [ ! -s "${LORALS_ASTS}" ]; then
    echo "ERROR: private asts.py missing"
    exit 1
fi

"${PYTHON}" - "${LORALS_ASTS}" <<'PY'
import re
import sys
from pathlib import Path

path = Path(sys.argv[1])
text = path.read_text()

if (
    "statistic_location=utils.nan" in text
    and "statistic_sign=utils.nan" in text
):
    sys.exit(0)

pattern = re.compile(
    r"stats\.stats\.KstestResult"
    r"\(\s*"
    r"statistic\s*=\s*utils\.nan"
    r"\s*,\s*"
    r"pvalue\s*=\s*utils\.nan"
    r"\s*\)"
)

replacement = (
    "stats.stats.KstestResult("
    "statistic=utils.nan, "
    "pvalue=utils.nan, "
    "statistic_location=utils.nan, "
    "statistic_sign=utils.nan"
    ")"
)

new_text, n = pattern.subn(replacement, text)
if n == 0:
    raise SystemExit("ERROR: expected KstestResult expression not found")

path.write_text(new_text)
PY

export PYTHONPATH="${LORALS_PATCH_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
export LORALS_PATCH_ROOT

"${PYTHON}" - <<'PY'
import os
import lorals
from lorals.scripts import calc_asts

expected = os.path.realpath(os.environ["LORALS_PATCH_ROOT"])
actual = os.path.realpath(lorals.__file__)
print("lorals loaded from =", actual)

if not actual.startswith(expected + os.sep):
    raise SystemExit("ERROR: wrong LORALS package loaded")

print("calc_asts import OK")
PY

ROW_START=$(( BATCH_START + (SLURM_ARRAY_TASK_ID - 1) * SAMPLES_PER_TASK ))
ROW_END=$(( ROW_START + SAMPLES_PER_TASK - 1 ))

if [ "${ROW_START}" -gt "${BATCH_END}" ]; then
    echo "No samples assigned to array task ${SLURM_ARRAY_TASK_ID}"
    exit 0
fi

if [ "${ROW_END}" -gt "${BATCH_END}" ]; then
    ROW_END=${BATCH_END}
fi

ROW_LIST=${JOB_TMP}/rows.txt
seq \
    "${ROW_START}" \
    "${ROW_END}" \
    > "${ROW_LIST}"

N_TASK=$(
    wc -l < "${ROW_LIST}" \
    | tr -d '[:space:]'
)

echo "============================================================"
echo "08-3b calc_asts quant"
echo "batch = ${BATCH_START}-${BATCH_END}"
echo "array task = ${SLURM_ARRAY_TASK_ID}"
echo "manifest rows = ${ROW_START}-${ROW_END}"
echo "samples = ${N_TASK}"
echo "xargs parallel = ${MAX_PARALLEL}"
echo "tmp root = ${TMPROOT}"
echo "start = $(date)"
echo "============================================================"

process_one()
{
    local ROW_ID=$1
    local LINE
    local array_id dataset age tissue original_tissue sample replicate vcf_sample bam
    local ASE_OUT ASE_FOR_ASTS ASTS_OUT TRANSCRIPT_BAM
    local ASTS_PREFIX ASTS_RAW DONE_A DONE_B DONE_C
    local SAMPLE_TMP GBAM ASTS_FOUND

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

    ASE_OUT=${ase}/samples/${tissue}/${dataset}/${sample}
    ASE_FOR_ASTS=${ASE_OUT}/${sample}.ase.for_asts.tsv

    ASTS_OUT=${asts}/samples/${tissue}/${dataset}/${sample}
    TRANSCRIPT_BAM=${ASTS_OUT}/${sample}.transcriptome.sorted.bam
    ASTS_PREFIX=${ASTS_OUT}/${sample}
    ASTS_RAW=${ASTS_PREFIX}_asts_quant.tsv

    DONE_A=${ASTS_OUT}/${sample}.08_3a_transcript_bam.SUCCESS
    DONE_B=${ASTS_OUT}/${sample}.08_3b_calc_asts.SUCCESS
    DONE_C=${ASTS_OUT}/${sample}.08_3c_summary.SUCCESS

    if [ ! -s "${bam}" ]; then
        echo "ERROR [${sample}]: source genome BAM missing"
        return 1
    fi

    if [ ! -s "${ASE_FOR_ASTS}" ]; then
        echo "ERROR [${sample}]: ASE_FOR_ASTS missing"
        echo "${ASE_FOR_ASTS}"
        return 1
    fi

    if [ ! -s "${DONE_A}" ] || \
       [ ! -s "${TRANSCRIPT_BAM}" ] || \
       [ ! -s "${TRANSCRIPT_BAM}.bai" ]; then

        echo "ERROR [${sample}]: 08-3a transcriptome BAM incomplete"
        return 1
    fi

    "${SAMTOOLS}" quickcheck \
        -q \
        "${TRANSCRIPT_BAM}"

    if [ -s "${DONE_B}" ] && \
       [ -s "${ASTS_RAW}" ] && \
       [ "${DONE_B}" -nt "${ASE_FOR_ASTS}" ] && \
       [ "${DONE_B}" -nt "${TRANSCRIPT_BAM}" ] && \
       [ "${DONE_B}" -nt "${bam}" ] && \
       [ "${ASTS_RAW}" -nt "${ASE_FOR_ASTS}" ] && \
       [ "${ASTS_RAW}" -nt "${TRANSCRIPT_BAM}" ]; then

        echo "[${sample}] SKIP: current calc_asts output already exists"
        return 0
    fi

    # Any new 08-3b output makes old 08-3c summaries stale.
    rm -f \
        "${DONE_B}" \
        "${DONE_C}" \
        "${ASTS_OUT}/sample_gene_ase.tsv" \
        "${ASTS_OUT}/sample_gene_asts.tsv" \
        "${ASTS_OUT}/sample_gene_ase_FDR005.tsv" \
        "${ASTS_OUT}/sample_gene_asts_FDR005.tsv"

    find "${ASTS_OUT}" \
        -maxdepth 1 \
        -type f \
        -name "${sample}*asts*quant*.tsv" \
        -delete

    SAMPLE_TMP=${JOB_TMP}/${sample}_${ROW_ID}
    rm -rf "${SAMPLE_TMP}"
    mkdir -p "${SAMPLE_TMP}"

    GBAM=${SAMPLE_TMP}/${sample}.genome.filtered.bam

    echo "============================================================"
    echo "START 08-3b sample = ${sample}"
    echo "row = ${ROW_ID}"
    echo "time = $(date)"
    echo "============================================================"

    echo "[${sample}] Filter genome BAM for calc_asts"
    "${SAMTOOLS}" view \
        -@ "${THREADS_PER_SAMPLE}" \
        -b \
        -q 10 \
        -F 2308 \
        -o "${GBAM}" \
        "${bam}"

    if [ ! -s "${GBAM}" ]; then
        echo "ERROR [${sample}]: filtered genome BAM missing"
        return 1
    fi

    "${SAMTOOLS}" index \
        -@ "${THREADS_PER_SAMPLE}" \
        "${GBAM}"

    if [ ! -s "${GBAM}.bai" ]; then
        echo "ERROR [${sample}]: filtered genome BAM index missing"
        return 1
    fi

    echo "[${sample}] calc_asts -m quant"
    "${CALC_ASTS}" \
        -m quant \
        -b "${GBAM}" \
        -i "${ASE_FOR_ASTS}" \
        -x "${TRANSCRIPT_BAM}" \
        -o "${ASTS_PREFIX}" \
        -v error

    ASTS_RAW=${ASTS_PREFIX}_asts_quant.tsv

    if [ ! -s "${ASTS_RAW}" ]; then
        ASTS_FOUND=$(
            find "${ASTS_OUT}" \
                -maxdepth 1 \
                -type f \
                -name "${sample}*asts*quant*.tsv" \
                -print \
                -quit
        )

        if [ -n "${ASTS_FOUND}" ] && [ -s "${ASTS_FOUND}" ]; then
            mv -f \
                "${ASTS_FOUND}" \
                "${ASTS_RAW}"
        else
            echo "ERROR [${sample}]: calc_asts output missing"
            return 1
        fi
    fi

    {
        echo -e "sample\t${sample}"
        echo -e "row_id\t${ROW_ID}"
        echo -e "ase_for_asts\t${ASE_FOR_ASTS}"
        echo -e "transcript_bam\t${TRANSCRIPT_BAM}"
        echo -e "asts_raw\t${ASTS_RAW}"
        echo -e "status\tSUCCESS"
        echo -e "finished\t$(date)"
    } > "${DONE_B}"

    rm -rf "${SAMPLE_TMP}"

    echo "[${sample}] 08-3b SUCCESS"
}

export -f process_one
export MANIFEST ase asts
export PYTHON SAMTOOLS CALC_ASTS
export THREADS_PER_SAMPLE JOB_TMP
export PATH PYTHONPATH LD_LIBRARY_PATH

set +e
xargs \
    -r \
    -n 1 \
    -P "${MAX_PARALLEL}" \
    bash -c 'set -Eeuo pipefail; process_one "$1"' _ \
    < "${ROW_LIST}"
XARGS_STATUS=$?
set -e

FAILED=0
while read -r ROW_ID; do
    [ -z "${ROW_ID}" ] && continue

    LINE=$(
        sed -n \
            "$((ROW_ID + 1))p" \
            "${MANIFEST}"
    )

    IFS=$'\t' read -r \
        array_id dataset age tissue original_tissue sample replicate vcf_sample bam \
        <<< "${LINE}"

    ASTS_OUT=${asts}/samples/${tissue}/${dataset}/${sample}
    ASTS_RAW=${ASTS_OUT}/${sample}_asts_quant.tsv
    DONE_B=${ASTS_OUT}/${sample}.08_3b_calc_asts.SUCCESS

    if [ ! -s "${ASTS_RAW}" ] || [ ! -s "${DONE_B}" ]; then
        echo "FAILED 08-3b: ${sample}"
        FAILED=$((FAILED + 1))
    else
        echo "SUCCESS 08-3b: ${sample}"
    fi
done < "${ROW_LIST}"

echo "============================================================"
echo "08-3b summary"
echo "batch = ${BATCH_START}-${BATCH_END}"
echo "array task = ${SLURM_ARRAY_TASK_ID}"
echo "xargs status = ${XARGS_STATUS}"
echo "failed samples = ${FAILED}"
echo "finish = $(date)"
echo "============================================================"

if [ "${XARGS_STATUS}" -ne 0 ] || \
   [ "${FAILED}" -gt 0 ]; then

    exit 1
fi
