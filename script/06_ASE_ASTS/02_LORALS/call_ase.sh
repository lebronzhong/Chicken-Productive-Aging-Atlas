#!/usr/bin/env bash
# Chicken Productive Aging Atlas | 06_ASE_ASTS
# Purpose: ASE and ASTS analysis.
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
MANIFEST=${scripts}/lorals_manifest_matched.tsv

ENV=/path/to/chicken-aging/software/lorals_env

PYTHON=${ENV}/bin/python
SAMTOOLS=${ENV}/bin/samtools
BCFTOOLS=${ENV}/bin/bcftools
CALC_ASE=${ENV}/bin/calc_ase

export PATH=${ENV}/bin:${PATH}
export LD_LIBRARY_PATH=${ENV}/lib:${LD_LIBRARY_PATH:-}

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export PYTHONNOUSERSITE=1

BLOCK_SIZE=10
MAX_PARALLEL=10
THREADS_PER_SAMPLE=5

mkdir -p \
    "${ase}/samples" \
    "${software}/tmp" \
    "${logs}"

cd "${work}"

# 2
for TOOL in \
    "${PYTHON}" \
    "${SAMTOOLS}" \
    "${BCFTOOLS}" \
    "${CALC_ASE}"; do

    if [ ! -x "${TOOL}" ]; then
        echo "ERROR: executable not found"
        echo "${TOOL}"
        exit 1
    fi
done

if [ ! -s "${MANIFEST}" ]; then
    echo "ERROR: manifest not found"
    echo "${MANIFEST}"
    exit 1
fi

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

LORALS_PATCH_ROOT=${software}/tmp/lorals_calc_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

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

from lorals.scripts import calc_ase

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
    "calc_ase import OK"
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
    echo "17 original tissues + Magnum = 18 tissues"
    echo "18 tissues x 25 samples = 450"
    echo "Found = ${N_TOTAL}"
    exit 1
fi

# 6
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

ROW_LIST=${software}/tmp/ASE_calc_rows_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.txt

seq \
    "${START}" \
    "${END}" \
    > "${ROW_LIST}"

N_TASK=$(
    wc -l < "${ROW_LIST}"
)

N_TASK=$(
    echo "${N_TASK}" \
    | tr -d '[:space:]'
)

echo "============================================================"
echo "08_1 calc_ase"
echo "============================================================"
echo "Total samples = ${N_TOTAL}"
echo "BLOCK_SIZE = ${BLOCK_SIZE}"
echo "MAX_PARALLEL = ${MAX_PARALLEL}"
echo "Array task = ${SLURM_ARRAY_TASK_ID}"
echo "Rows = ${START}-${END}"
echo "Samples in this array = ${N_TASK}"
echo "TIME = $(date)"
echo "============================================================"

# 7
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
       [ -z "${tissue}" ] || \
       [ -z "${dataset}" ] || \
       [ -z "${vcf_sample}" ] || \
       [ -z "${bam}" ]; then

        echo "ERROR: invalid manifest row ${ROW_ID}"
        echo "${LINE}"
        return 1
    fi

    SAMPLE_VCF_GZ=${VCFDIR}/${vcf_sample}.het.biallelic.vcf.gz

    ASE_OUT=${ase}/samples/${tissue}/${dataset}/${sample}

    ASE_RAW=${ASE_OUT}/${sample}.ase.raw.tsv

    DONE=${ASE_OUT}/${sample}.ASE_CALC.SUCCESS

    mkdir -p "${ASE_OUT}"

    echo "============================================================"
    echo "START sample = ${sample}"
    echo "ROW = ${ROW_ID}"
    echo "dataset = ${dataset}"
    echo "age = ${age}"
    echo "tissue = ${tissue}"
    echo "original_tissue = ${original_tissue}"
    echo "replicate = ${replicate}"
    echo "vcf_sample = ${vcf_sample}"
    echo "bam = ${bam}"
    echo "TIME = $(date)"
    echo "============================================================"

    # 8
    if [ -s "${ASE_RAW}" ]; then

        echo "[${sample}] SKIP calc_ase"
        echo "[${sample}] Existing raw ASE protected"
        echo "${ASE_RAW}"

        {
            echo -e "sample\t${sample}"
            echo -e "dataset\t${dataset}"
            echo -e "age\t${age}"
            echo -e "tissue\t${tissue}"
            echo -e "vcf_sample\t${vcf_sample}"
            echo -e "status\tEXISTING_RAW_ASE"
            echo -e "raw_ase\t${ASE_RAW}"
            echo -e "finished\t$(date)"
        } > "${DONE}"

        return 0
    fi

    # 9
    if [ ! -s "${SAMPLE_VCF_GZ}" ]; then
        echo "ERROR [${sample}]: individual heterozygous VCF missing"
        echo "${SAMPLE_VCF_GZ}"
        return 1
    fi

    if [ ! -s "${SAMPLE_VCF_GZ}.tbi" ]; then
        echo "ERROR [${sample}]: individual VCF index missing"
        echo "${SAMPLE_VCF_GZ}.tbi"
        return 1
    fi

    if [ ! -s "${bam}" ]; then
        echo "ERROR [${sample}]: BAM missing"
        echo "${bam}"
        return 1
    fi

    # 10
    TMP=${software}/tmp/ASE_calc_${sample}_${ROW_ID}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

    rm -rf "${TMP}"

    mkdir -p "${TMP}"

    SAMPLE_VCF=${TMP}/${sample}.vcf

    GBAM=${TMP}/${sample}.genome.filtered.bam

    # 11
    "${BCFTOOLS}" view \
        -Ov \
        -o "${SAMPLE_VCF}" \
        "${SAMPLE_VCF_GZ}"

    if [ ! -s "${SAMPLE_VCF}" ]; then
        echo "ERROR [${sample}]: VCF conversion failed"

        rm -rf "${TMP}"

        return 1
    fi

    # 12
    "${SAMTOOLS}" view \
        -@ "${THREADS_PER_SAMPLE}" \
        -b \
        -q 10 \
        -F 2308 \
        -o "${GBAM}" \
        "${bam}"

    if [ ! -s "${GBAM}" ]; then
        echo "ERROR [${sample}]: filtered genome BAM missing"

        rm -rf "${TMP}"

        return 1
    fi

    "${SAMTOOLS}" index \
        -@ "${THREADS_PER_SAMPLE}" \
        "${GBAM}"

    if [ ! -s "${GBAM}.bai" ]; then
        echo "ERROR [${sample}]: filtered genome BAM index missing"

        rm -rf "${TMP}"

        return 1
    fi

    "${SAMTOOLS}" flagstat \
        -@ 1 \
        "${GBAM}" \
        > "${ASE_OUT}/${sample}.genome.flagstat.txt"

    # 13
    "${CALC_ASE}" \
        -b "${GBAM}" \
        -f "${SAMPLE_VCF}" \
        -o "${ASE_RAW}" \
        -v error

    if [ ! -s "${ASE_RAW}" ]; then
        echo "ERROR [${sample}]: calc_ase output missing"

        rm -rf "${TMP}"

        return 1
    fi

    # 14
    {
        echo -e "sample\t${sample}"
        echo -e "dataset\t${dataset}"
        echo -e "age\t${age}"
        echo -e "tissue\t${tissue}"
        echo -e "vcf_sample\t${vcf_sample}"
        echo -e "status\tSUCCESS"
        echo -e "raw_ase\t${ASE_RAW}"
        echo -e "finished\t$(date)"
    } > "${DONE}"

    rm -rf "${TMP}"

    echo "============================================================"
    echo "FINISHED sample = ${sample}"
    echo "raw ASE = ${ASE_RAW}"
    echo "TIME = $(date)"
    echo "============================================================"
}

# 15
export -f process_one

export MANIFEST
export VCFDIR

export ase
export software

export SAMTOOLS
export BCFTOOLS
export CALC_ASE

export THREADS_PER_SAMPLE

export PATH
export PYTHONPATH
export LD_LIBRARY_PATH

export SLURM_ARRAY_JOB_ID
export SLURM_ARRAY_TASK_ID

# 16
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

# 17
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

    ASE_RAW=${ase}/samples/${tissue}/${dataset}/${sample}/${sample}.ase.raw.tsv

    if [ ! -s "${ASE_RAW}" ]; then

        echo "FAILED: ${sample} raw ASE missing"

        FAILED=$((FAILED + 1))

    else

        echo "SUCCESS: ${sample} raw ASE exists"
    fi

done < "${ROW_LIST}"

# 18
rm -f "${ROW_LIST}"

rm -rf "${LORALS_PATCH_ROOT}"

echo "============================================================"
echo "08_1 array summary"
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

    echo "ERROR: 08_1 array ${SLURM_ARRAY_TASK_ID} failed"

    exit 1
fi

echo "08_1 array ${SLURM_ARRAY_TASK_ID} completed successfully"
