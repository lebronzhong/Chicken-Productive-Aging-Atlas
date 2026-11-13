#!/usr/bin/env bash
# Chicken Productive Aging Atlas | 06_ASE_ASTS
# Purpose: Summarize genes.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

set -Eeuo pipefail

############################################################
# 08-3c
# calc_asts raw output + transcript->gene mapping
# -> gene-level ASE / ASTS + BH-FDR
############################################################

work=/path/to/chicken-aging/allelic
asts=${work}/asts
scripts=${work}/scripts
software=${work}/software
logs=${work}/logs

REFDIR=${software}/lorals_ref
MANIFEST=${scripts}/lorals_manifest_matched.tsv
TRANSCRIPT_GENE=${REFDIR}/transcript_gene.tsv
SUMMARIZE=${scripts}/summarize_sample.py

ENV=/path/to/chicken-aging/software/lorals_env
PYTHON=${ENV}/bin/python

export PATH=${ENV}/bin:${PATH}
export LD_LIBRARY_PATH=${ENV}/lib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

SAMPLES_PER_TASK=2
MAX_PARALLEL=2

: "${BATCH_START:?BATCH_START is required}"
: "${BATCH_END:?BATCH_END is required}"
: "${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"

mkdir -p "${logs}"
cd "${work}"

if [ ! -x "${PYTHON}" ]; then
    echo "ERROR: python not executable: ${PYTHON}"
    exit 1
fi

for FILE in \
    "${MANIFEST}" \
    "${TRANSCRIPT_GENE}" \
    "${SUMMARIZE}"; do

    if [ ! -s "${FILE}" ]; then
        echo "ERROR: required file missing: ${FILE}"
        exit 1
    fi
done

MANIFEST_LINES=$(
    wc -l < "${MANIFEST}" \
    | tr -d '[:space:]'
)
N_TOTAL=$(( MANIFEST_LINES - 1 ))

if [ "${N_TOTAL}" -ne 450 ]; then
    echo "ERROR: expected 450 matched ONT samples; found ${N_TOTAL}"
    exit 1
fi

TMPROOT=${SLURM_TMPDIR:-/tmp}
if [ ! -d "${TMPROOT}" ] || [ ! -w "${TMPROOT}" ]; then
    TMPROOT=${software}/tmp
fi

JOB_TMP=${TMPROOT}/ASTS_3c_${SLURM_ARRAY_JOB_ID:-$$}_${SLURM_ARRAY_TASK_ID}
rm -rf "${JOB_TMP}"
mkdir -p "${JOB_TMP}"
trap 'rm -rf "${JOB_TMP}"' EXIT

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

echo "============================================================"
echo "08-3c gene-level summary"
echo "batch = ${BATCH_START}-${BATCH_END}"
echo "array task = ${SLURM_ARRAY_TASK_ID}"
echo "manifest rows = ${ROW_START}-${ROW_END}"
echo "start = $(date)"
echo "============================================================"

process_one()
{
    local ROW_ID=$1
    local LINE
    local array_id dataset age tissue original_tissue sample replicate vcf_sample bam
    local ASTS_OUT ASTS_RAW DONE_B DONE_C
    local FINAL_GENE_ASE FINAL_GENE_ASTS

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

    ASTS_OUT=${asts}/samples/${tissue}/${dataset}/${sample}
    ASTS_RAW=${ASTS_OUT}/${sample}_asts_quant.tsv
    DONE_B=${ASTS_OUT}/${sample}.08_3b_calc_asts.SUCCESS
    DONE_C=${ASTS_OUT}/${sample}.08_3c_summary.SUCCESS

    FINAL_GENE_ASE=${ASTS_OUT}/sample_gene_ase.tsv
    FINAL_GENE_ASTS=${ASTS_OUT}/sample_gene_asts.tsv

    if [ ! -s "${DONE_B}" ] || [ ! -s "${ASTS_RAW}" ]; then
        echo "ERROR [${sample}]: 08-3b output incomplete"
        return 1
    fi

    if [ -s "${DONE_C}" ] && \
       [ -s "${FINAL_GENE_ASE}" ] && \
       [ -s "${FINAL_GENE_ASTS}" ] && \
       [ "${DONE_C}" -nt "${ASTS_RAW}" ] && \
       [ "${FINAL_GENE_ASE}" -nt "${ASTS_RAW}" ] && \
       [ "${FINAL_GENE_ASTS}" -nt "${ASTS_RAW}" ]; then

        echo "[${sample}] SKIP: current 08-3c summary already exists"
        return 0
    fi

    rm -f \
        "${FINAL_GENE_ASE}" \
        "${FINAL_GENE_ASTS}" \
        "${ASTS_OUT}/sample_gene_ase_FDR005.tsv" \
        "${ASTS_OUT}/sample_gene_asts_FDR005.tsv" \
        "${DONE_C}"

    echo "[${sample}] gene-level ASE / ASTS summary"
    "${PYTHON}" \
        "${SUMMARIZE}" \
        "${ASTS_RAW}" \
        "${TRANSCRIPT_GENE}" \
        "${ASTS_OUT}" \
        "${sample}" \
        "${dataset}" \
        "${tissue}" \
        "${original_tissue}"

    if [ ! -s "${FINAL_GENE_ASE}" ]; then
        echo "ERROR [${sample}]: sample_gene_ase.tsv missing"
        return 1
    fi

    if [ ! -s "${FINAL_GENE_ASTS}" ]; then
        echo "ERROR [${sample}]: sample_gene_asts.tsv missing"
        return 1
    fi

    {
        echo -e "sample\t${sample}"
        echo -e "row_id\t${ROW_ID}"
        echo -e "asts_raw\t${ASTS_RAW}"
        echo -e "gene_ase\t${FINAL_GENE_ASE}"
        echo -e "gene_asts\t${FINAL_GENE_ASTS}"
        echo -e "status\tSUCCESS"
        echo -e "finished\t$(date)"
    } > "${DONE_C}"

    echo "[${sample}] 08-3c SUCCESS"
}

export -f process_one
export MANIFEST asts TRANSCRIPT_GENE SUMMARIZE PYTHON
export PATH LD_LIBRARY_PATH

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
    DONE_C=${ASTS_OUT}/${sample}.08_3c_summary.SUCCESS

    if [ ! -s "${DONE_C}" ] || \
       [ ! -s "${ASTS_OUT}/sample_gene_ase.tsv" ] || \
       [ ! -s "${ASTS_OUT}/sample_gene_asts.tsv" ]; then

        echo "FAILED 08-3c: ${sample}"
        FAILED=$((FAILED + 1))
    else
        echo "SUCCESS 08-3c: ${sample}"
    fi
done < "${ROW_LIST}"

echo "============================================================"
echo "08-3c summary"
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
