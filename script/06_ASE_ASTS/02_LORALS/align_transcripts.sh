#!/usr/bin/env bash
# Chicken Productive Aging Atlas | 06_ASE_ASTS
# Purpose: Align transcripts.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

set -Eeuo pipefail

############################################################
# 08-3a
# Each Slurm array task processes at most 2 samples in parallel.
# Supply BATCH_START and BATCH_END as environment variables.
############################################################

work=/path/to/chicken-aging/allelic
ase=${work}/ase
asts=${work}/asts
scripts=${work}/scripts
software=${work}/software
logs=${work}/logs

REFDIR=${software}/lorals_ref
MANIFEST=${scripts}/lorals_manifest_matched.tsv
TRANSCRIPTOME=${REFDIR}/transcriptome.fa
TRANSCRIPTOME_MMI=${REFDIR}/transcriptome.mmi

ENV=/path/to/chicken-aging/software/lorals_env
SAMTOOLS=${ENV}/bin/samtools
MINIMAP2=${ENV}/bin/minimap2

export PATH=${ENV}/bin:${PATH}
export LD_LIBRARY_PATH=${ENV}/lib:${LD_LIBRARY_PATH:-}

SAMPLES_PER_TASK=2
MAX_PARALLEL=2
THREADS_PER_SAMPLE=14
VIEW_THREADS=2
SORT_THREADS=8
SORT_MEM=1G

: "${BATCH_START:?BATCH_START is required}"
: "${BATCH_END:?BATCH_END is required}"
: "${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"

mkdir -p \
    "${asts}/samples" \
    "${software}/tmp" \
    "${logs}"

cd "${work}"

for TOOL in \
    "${SAMTOOLS}" \
    "${MINIMAP2}"; do

    if [ ! -x "${TOOL}" ]; then
        echo "ERROR: executable not found: ${TOOL}"
        exit 1
    fi
done

for FILE in \
    "${MANIFEST}" \
    "${TRANSCRIPTOME}" \
    "${TRANSCRIPTOME_MMI}"; do

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

JOB_TMP=${TMPROOT}/ASTS_3a_${SLURM_ARRAY_JOB_ID:-$$}_${SLURM_ARRAY_TASK_ID}
rm -rf "${JOB_TMP}"
mkdir -p "${JOB_TMP}"

cleanup_job_tmp()
{
    rm -rf "${JOB_TMP}"
}
trap cleanup_job_tmp EXIT

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
echo "08-3a prepare transcriptome BAM"
echo "batch = ${BATCH_START}-${BATCH_END}"
echo "array task = ${SLURM_ARRAY_TASK_ID}"
echo "manifest rows = ${ROW_START}-${ROW_END}"
echo "samples = ${N_TASK}"
echo "xargs parallel = ${MAX_PARALLEL}"
echo "threads/sample = ${THREADS_PER_SAMPLE}"
echo "tmp root = ${TMPROOT}"
df -h "${TMPROOT}" || true
echo "start = $(date)"
echo "============================================================"

process_one()
{
    local ROW_ID=$1
    local LINE
    local array_id dataset age tissue original_tissue sample replicate vcf_sample bam
    local ASTS_OUT TRANSCRIPT_BAM TRANSCRIPT_BAI FLAGSTAT DONE
    local SAMPLE_TMP GBAM FASTQ TRANSCRIPT_UNSORTED TRANSCRIPT_SORTED

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
       [ -z "${tissue}" ] || \
       [ -z "${dataset}" ] || \
       [ -z "${bam}" ]; then

        echo "ERROR: invalid manifest row ${ROW_ID}"
        echo "${LINE}"
        return 1
    fi

    if [ ! -s "${bam}" ]; then
        echo "ERROR [${sample}]: source genome BAM missing"
        echo "${bam}"
        return 1
    fi

    ASTS_OUT=${asts}/samples/${tissue}/${dataset}/${sample}
    TRANSCRIPT_BAM=${ASTS_OUT}/${sample}.transcriptome.sorted.bam
    TRANSCRIPT_BAI=${TRANSCRIPT_BAM}.bai
    FLAGSTAT=${ASTS_OUT}/${sample}.transcriptome.flagstat.txt
    DONE=${ASTS_OUT}/${sample}.08_3a_transcript_bam.SUCCESS

    mkdir -p "${ASTS_OUT}"

    if [ -s "${TRANSCRIPT_BAM}" ] && \
       [ -s "${TRANSCRIPT_BAI}" ] && \
       [ "${TRANSCRIPT_BAM}" -nt "${bam}" ] && \
       [ "${TRANSCRIPT_BAM}" -nt "${TRANSCRIPTOME_MMI}" ] && \
       "${SAMTOOLS}" quickcheck -q "${TRANSCRIPT_BAM}"; then

        if [ ! -s "${DONE}" ]; then
            {
                echo -e "sample\t${sample}"
                echo -e "row_id\t${ROW_ID}"
                echo -e "transcript_bam\t${TRANSCRIPT_BAM}"
                echo -e "status\tSUCCESS"
                echo -e "finished\t$(date)"
            } > "${DONE}"
        fi

        echo "[${sample}] SKIP: valid transcriptome BAM already exists"
        return 0
    fi

    rm -f \
        "${TRANSCRIPT_BAM}" \
        "${TRANSCRIPT_BAI}" \
        "${FLAGSTAT}" \
        "${DONE}"

    SAMPLE_TMP=${JOB_TMP}/${sample}_${ROW_ID}
    rm -rf "${SAMPLE_TMP}"
    mkdir -p "${SAMPLE_TMP}"

    GBAM=${SAMPLE_TMP}/${sample}.genome.filtered.bam
    FASTQ=${SAMPLE_TMP}/${sample}.fastq.gz
    TRANSCRIPT_UNSORTED=${SAMPLE_TMP}/${sample}.transcriptome.unsorted.bam
    TRANSCRIPT_SORTED=${SAMPLE_TMP}/${sample}.transcriptome.sorted.bam

    echo "============================================================"
    echo "START 08-3a sample = ${sample}"
    echo "row = ${ROW_ID}"
    echo "source BAM = ${bam}"
    echo "time = $(date)"
    echo "============================================================"

    echo "[${sample}] Filter genome BAM"
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

    echo "[${sample}] filtered BAM -> FASTQ"
    "${SAMTOOLS}" fastq \
        -@ "${THREADS_PER_SAMPLE}" \
        "${GBAM}" \
        | gzip -c \
        > "${FASTQ}"

    if [ ! -s "${FASTQ}" ]; then
        echo "ERROR [${sample}]: FASTQ generation failed"
        return 1
    fi

    echo "[${sample}] minimap2 -> transcriptome BAM"
    "${MINIMAP2}" \
        -t "${THREADS_PER_SAMPLE}" \
        -ax map-ont \
        --secondary=no \
        "${TRANSCRIPTOME_MMI}" \
        "${FASTQ}" \
        | "${SAMTOOLS}" view \
            -@ "${VIEW_THREADS}" \
            -b \
            -F 2308 \
            -o "${TRANSCRIPT_UNSORTED}" \
            -

    if [ ! -s "${TRANSCRIPT_UNSORTED}" ]; then
        echo "ERROR [${sample}]: unsorted transcriptome BAM missing"
        return 1
    fi

    echo "[${sample}] sort transcriptome BAM"
    "${SAMTOOLS}" sort \
        -@ "${SORT_THREADS}" \
        -m "${SORT_MEM}" \
        -o "${TRANSCRIPT_SORTED}" \
        "${TRANSCRIPT_UNSORTED}"

    if [ ! -s "${TRANSCRIPT_SORTED}" ]; then
        echo "ERROR [${sample}]: sorted transcriptome BAM missing"
        return 1
    fi

    "${SAMTOOLS}" quickcheck \
        -q \
        "${TRANSCRIPT_SORTED}"

    "${SAMTOOLS}" index \
        -@ "${SORT_THREADS}" \
        "${TRANSCRIPT_SORTED}"

    if [ ! -s "${TRANSCRIPT_SORTED}.bai" ]; then
        echo "ERROR [${sample}]: transcriptome BAM index missing"
        return 1
    fi

    "${SAMTOOLS}" flagstat \
        -@ 1 \
        "${TRANSCRIPT_SORTED}" \
        > "${FLAGSTAT}.tmp"

    mv -f \
        "${TRANSCRIPT_SORTED}" \
        "${TRANSCRIPT_BAM}"

    mv -f \
        "${TRANSCRIPT_SORTED}.bai" \
        "${TRANSCRIPT_BAI}"

    mv -f \
        "${FLAGSTAT}.tmp" \
        "${FLAGSTAT}"

    "${SAMTOOLS}" quickcheck \
        -q \
        "${TRANSCRIPT_BAM}"

    if [ ! -s "${TRANSCRIPT_BAM}" ] || \
       [ ! -s "${TRANSCRIPT_BAI}" ]; then

        echo "ERROR [${sample}]: persistent transcriptome BAM validation failed"
        return 1
    fi

    {
        echo -e "sample\t${sample}"
        echo -e "row_id\t${ROW_ID}"
        echo -e "source_bam\t${bam}"
        echo -e "transcript_bam\t${TRANSCRIPT_BAM}"
        echo -e "status\tSUCCESS"
        echo -e "finished\t$(date)"
    } > "${DONE}"

    # FASTQ, filtered genome BAM and unsorted BAM are no longer needed.
    rm -rf "${SAMPLE_TMP}"

    echo "[${sample}] 08-3a SUCCESS"
}

export -f process_one
export MANIFEST asts TRANSCRIPTOME_MMI
export SAMTOOLS MINIMAP2
export THREADS_PER_SAMPLE VIEW_THREADS SORT_THREADS SORT_MEM
export JOB_TMP
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
    TRANSCRIPT_BAM=${ASTS_OUT}/${sample}.transcriptome.sorted.bam
    DONE=${ASTS_OUT}/${sample}.08_3a_transcript_bam.SUCCESS

    if [ ! -s "${TRANSCRIPT_BAM}" ] || \
       [ ! -s "${TRANSCRIPT_BAM}.bai" ] || \
       [ ! -s "${DONE}" ] || \
       ! "${SAMTOOLS}" quickcheck -q "${TRANSCRIPT_BAM}"; then

        echo "FAILED 08-3a: ${sample}"
        FAILED=$((FAILED + 1))
    else
        echo "SUCCESS 08-3a: ${sample}"
    fi
done < "${ROW_LIST}"

echo "============================================================"
echo "08-3a summary"
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
