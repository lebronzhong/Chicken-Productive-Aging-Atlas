#!/usr/bin/env bash
# Chicken Productive Aging Atlas | 06_ASE_ASTS
# Purpose: Prepare vcfs.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

set -euo pipefail


############################################################
# 01. Define directories
############################################################

work=/path/to/chicken-aging/allelic

ase=${work}/ase
scripts=${work}/scripts
logs=${work}/logs

VCFDIR=${ase}/vcf
SUMMARY_DIR=${VCFDIR}/summary
TMPDIR=${VCFDIR}/tmp

SAMPLE_LIST=${scripts}/needed_vcf_samples.txt


############################################################
# 02. Input files
############################################################

VCF=/path/to/chicken-aging/WGS/CombineGVCF/Chicken_Aging_WGS.vcf.gz

REF=/path/to/chicken-aging/reference/galGal6/Gallus_gallus.GRCg6a.dna.toplevel.fa


############################################################
# 03. Conda environment
############################################################

ENV=/path/to/chicken-aging/software/lorals_env

BCFTOOLS=${ENV}/bin/bcftools

export PATH=${ENV}/bin:${PATH}
export LD_LIBRARY_PATH=${ENV}/lib:${LD_LIBRARY_PATH:-}


############################################################
# 04. Parallel settings
#
# 25 WGS samples
#
# Array 1:
#     samples 1-13
#
# Array 2:
#     samples 14-25
#
# Each array task:
#     maximum 13 samples simultaneously
############################################################

BLOCK_SIZE=13
MAX_PARALLEL=13
THREADS_PER_SAMPLE=3


############################################################
# 05. Working directory
############################################################

cd "${work}"


############################################################
# 06. Create directories
############################################################

mkdir -p "${VCFDIR}"
mkdir -p "${SUMMARY_DIR}"
mkdir -p "${TMPDIR}"
mkdir -p "${logs}"


############################################################
# 07. Basic information
############################################################

echo "============================================================"
echo "06 Individual heterozygous biallelic SNP VCF preparation"
echo "============================================================"
echo "HOST = $(hostname)"
echo "SLURM_JOB_ID = ${SLURM_JOB_ID}"
echo "SLURM_ARRAY_JOB_ID = ${SLURM_ARRAY_JOB_ID}"
echo "SLURM_ARRAY_TASK_ID = ${SLURM_ARRAY_TASK_ID}"
echo "SLURM_CPUS_PER_TASK = ${SLURM_CPUS_PER_TASK:-NA}"
echo "BLOCK_SIZE = ${BLOCK_SIZE}"
echo "MAX_PARALLEL = ${MAX_PARALLEL}"
echo "THREADS_PER_SAMPLE = ${THREADS_PER_SAMPLE}"
echo "TIME = $(date)"
echo "============================================================"


############################################################
# 08. Check software
############################################################

if [ ! -x "${BCFTOOLS}" ]; then
    echo "ERROR: bcftools not found"
    echo "${BCFTOOLS}"
    exit 1
fi

"${BCFTOOLS}" --version


############################################################
# 09. Check input files
############################################################

if [ ! -s "${VCF}" ]; then
    echo "ERROR: source VCF not found"
    echo "${VCF}"
    exit 1
fi

if [ ! -s "${SAMPLE_LIST}" ]; then
    echo "ERROR: sample list not found"
    echo "${SAMPLE_LIST}"
    exit 1
fi

if [ ! -s "${REF}" ]; then
    echo "ERROR: reference genome not found"
    echo "${REF}"
    exit 1
fi


############################################################
# 10. Check required samples exist in source VCF
############################################################

SOURCE_SAMPLE_LIST=${TMPDIR}/source_samples_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.txt

"${BCFTOOLS}" query \
    -l "${VCF}" \
    > "${SOURCE_SAMPLE_LIST}"

MISSING_SAMPLE=0

while read -r sample; do
    [ -z "${sample}" ] && continue

    if ! grep -Fxq "${sample}" "${SOURCE_SAMPLE_LIST}"; then
        echo "ERROR: sample not found in source VCF: ${sample}"
        MISSING_SAMPLE=$((MISSING_SAMPLE + 1))
    fi
done < "${SAMPLE_LIST}"

if [ "${MISSING_SAMPLE}" -gt 0 ]; then
    echo "ERROR: ${MISSING_SAMPLE} required samples are absent from source VCF"
    exit 1
fi


############################################################
# 11. Count WGS individuals
############################################################

N_TOTAL=$(grep -cve '^[[:space:]]*$' "${SAMPLE_LIST}")

echo "============================================================"
echo "Total required WGS samples = ${N_TOTAL}"
echo "============================================================"


############################################################
# 12. Determine sample range for this array task
############################################################

START=$(( (SLURM_ARRAY_TASK_ID - 1) * BLOCK_SIZE + 1 ))
END=$(( SLURM_ARRAY_TASK_ID * BLOCK_SIZE ))

if [ "${END}" -gt "${N_TOTAL}" ]; then
    END=${N_TOTAL}
fi

if [ "${START}" -gt "${N_TOTAL}" ]; then
    echo "No samples assigned to array ${SLURM_ARRAY_TASK_ID}"
    exit 0
fi

TASK_SAMPLE_LIST=${TMPDIR}/array_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.samples.txt

sed -n \
    "${START},${END}p" \
    "${SAMPLE_LIST}" \
    > "${TASK_SAMPLE_LIST}"

N_TASK=$(grep -cve '^[[:space:]]*$' "${TASK_SAMPLE_LIST}")

echo "============================================================"
echo "Array task = ${SLURM_ARRAY_TASK_ID}"
echo "Sample rows = ${START}-${END}"
echo "Samples in this array = ${N_TASK}"
echo "============================================================"

cat "${TASK_SAMPLE_LIST}"


############################################################
# 13. Function: process one WGS individual
############################################################

process_sample() {
    sample=$1

    echo "============================================================"
    echo "START sample = ${sample}"
    echo "HOST = $(hostname)"
    echo "TIME = $(date)"
    echo "============================================================"


    ########################################################
    # Output files
    ########################################################

    OUT=${VCFDIR}/${sample}.het.biallelic.vcf.gz
    SUMMARY=${SUMMARY_DIR}/${sample}.tsv

    TMP_SUBSET=${TMPDIR}/${sample}.${SLURM_ARRAY_JOB_ID}.${SLURM_ARRAY_TASK_ID}.subset.bcf
    TMP_HET=${TMPDIR}/${sample}.${SLURM_ARRAY_JOB_ID}.${SLURM_ARRAY_TASK_ID}.het.bcf
    GT_QC=${TMPDIR}/${sample}.${SLURM_ARRAY_JOB_ID}.${SLURM_ARRAY_TASK_ID}.GT_QC.tsv


    ########################################################
    # Remove previous incorrect / incomplete outputs
    ########################################################

    rm -f "${OUT}"
    rm -f "${OUT}.tbi"
    rm -f "${OUT}.csi"
    rm -f "${SUMMARY}"

    rm -f "${TMP_SUBSET}"
    rm -f "${TMP_SUBSET}.csi"

    rm -f "${TMP_HET}"
    rm -f "${TMP_HET}.csi"

    rm -f "${GT_QC}"


    ########################################################
    # STEP 1
    #
    # First subset ONE individual from the joint-genotyped VCF.
    #
    # Important:
    # Do NOT apply -g het in the same bcftools view command.
    ########################################################

    echo "[${sample}] STEP 1: subset individual"

    "${BCFTOOLS}" view \
        --threads "${THREADS_PER_SAMPLE}" \
        -s "${sample}" \
        -Ob \
        -o "${TMP_SUBSET}" \
        "${VCF}"

    if [ ! -s "${TMP_SUBSET}" ]; then
        echo "ERROR [${sample}]: subset BCF not generated"
        return 1
    fi


    ########################################################
    # Validate subset contains exactly one sample
    ########################################################

    N_SUBSET_SAMPLE=$(
        "${BCFTOOLS}" query \
            -l "${TMP_SUBSET}" \
        | wc -l
    )

    SUBSET_SAMPLE_NAME=$(
        "${BCFTOOLS}" query \
            -l "${TMP_SUBSET}"
    )

    if [ "${N_SUBSET_SAMPLE}" -ne 1 ]; then
        echo "ERROR [${sample}]: subset contains ${N_SUBSET_SAMPLE} samples"
        return 1
    fi

    if [ "${SUBSET_SAMPLE_NAME}" != "${sample}" ]; then
        echo "ERROR [${sample}]: wrong sample retained"
        echo "Expected = ${sample}"
        echo "Found = ${SUBSET_SAMPLE_NAME}"
        return 1
    fi


    ########################################################
    # STEP 2
    #
    # Now the file contains only ONE individual.
    #
    # Select:
    #     SNP
    #     biallelic
    #     heterozygous genotype
    #     PASS or unfiltered sites
    ########################################################

    echo "[${sample}] STEP 2: select individual heterozygous biallelic SNPs"

    "${BCFTOOLS}" view \
        --threads "${THREADS_PER_SAMPLE}" \
        -v snps \
        -m2 \
        -M2 \
        -g het \
        -f PASS,. \
        -Ob \
        -o "${TMP_HET}" \
        "${TMP_SUBSET}"

    if [ ! -s "${TMP_HET}" ]; then
        echo "ERROR [${sample}]: heterozygous BCF not generated"
        return 1
    fi


    ########################################################
    # STEP 3
    #
    # Normalize against galGal6 reference
    ########################################################

    echo "[${sample}] STEP 3: normalize variants"

    "${BCFTOOLS}" norm \
        --threads "${THREADS_PER_SAMPLE}" \
        -f "${REF}" \
        -Oz \
        -o "${OUT}" \
        "${TMP_HET}"

    if [ ! -s "${OUT}" ]; then
        echo "ERROR [${sample}]: output VCF not generated"
        return 1
    fi


    ########################################################
    # STEP 4
    #
    # Index final VCF
    ########################################################

    echo "[${sample}] STEP 4: index final VCF"

    "${BCFTOOLS}" index \
        --threads "${THREADS_PER_SAMPLE}" \
        -f \
        -t \
        "${OUT}"

    if [ ! -s "${OUT}.tbi" ]; then
        echo "ERROR [${sample}]: TBI index not generated"
        return 1
    fi


    ########################################################
    # STEP 5
    #
    # Count final variants using index
    ########################################################

    N_HET=$(
        "${BCFTOOLS}" index \
            -n "${OUT}"
    )

    if [ "${N_HET}" -eq 0 ]; then
        echo "ERROR [${sample}]: zero heterozygous SNPs"
        return 1
    fi

    echo "[${sample}] Heterozygous SNP records = ${N_HET}"


    ########################################################
    # STEP 6
    #
    # Full-file genotype QC
    #
    # The final VCF must contain ONLY:
    #
    #     0/1
    #     1/0
    #     0|1
    #     1|0
    #
    # Any:
    #
    #     0/0
    #     1/1
    #     ./.
    #     other GT
    #
    # causes the sample to fail.
    ########################################################

    echo "[${sample}] STEP 6: genotype QC"

    "${BCFTOOLS}" query \
        -f '[%GT\n]' \
        "${OUT}" \
        | awk '
            BEGIN {
                total = 0
                het = 0
                nonhet = 0
                phased = 0
                unphased = 0
            }

            {
                total++

                if (
                    $1 == "0/1" ||
                    $1 == "1/0"
                ) {
                    het++
                    unphased++
                } else if (
                    $1 == "0|1" ||
                    $1 == "1|0"
                ) {
                    het++
                    phased++
                } else {
                    nonhet++
                }
            }

            END {
                print total "\t" het "\t" nonhet "\t" phased "\t" unphased
            }
        ' \
        > "${GT_QC}"

    IFS=$'\t' read -r \
        N_QUERY \
        N_GT_HET \
        N_NONHET \
        N_PHASED \
        N_UNPHASED \
        < "${GT_QC}"


    ########################################################
    # Validate index count vs genotype count
    ########################################################

    if [ "${N_QUERY}" -ne "${N_HET}" ]; then
        echo "ERROR [${sample}]: index count and GT count differ"
        echo "Index count = ${N_HET}"
        echo "GT count = ${N_QUERY}"
        return 1
    fi


    ########################################################
    # Validate all records are heterozygous
    ########################################################

    if [ "${N_GT_HET}" -ne "${N_HET}" ]; then
        echo "ERROR [${sample}]: not all variants are heterozygous"
        echo "Total records = ${N_HET}"
        echo "Heterozygous GT = ${N_GT_HET}"
        return 1
    fi

    if [ "${N_NONHET}" -ne 0 ]; then
        echo "ERROR [${sample}]: ${N_NONHET} non-heterozygous records remain"

        echo "Example non-heterozygous genotypes:"

        "${BCFTOOLS}" query \
            -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' \
            "${OUT}" \
            | awk '
                $5 != "0/1" &&
                $5 != "1/0" &&
                $5 != "0|1" &&
                $5 != "1|0"
            ' \
            | head -n 20

        return 1
    fi


    ########################################################
    # STEP 7
    #
    # Additional site-level validation
    #
    # Confirm all sites remain:
    #     SNP
    #     biallelic
    ########################################################

    N_NON_SNP=$(
        "${BCFTOOLS}" view \
            -H \
            -e 'TYPE="snp"' \
            "${OUT}" \
        | wc -l
    )

    if [ "${N_NON_SNP}" -ne 0 ]; then
        echo "ERROR [${sample}]: ${N_NON_SNP} non-SNP records remain"
        return 1
    fi

    N_MULTI=$(
        "${BCFTOOLS}" view \
            -H \
            -m3 \
            "${OUT}" \
        | wc -l
    )

    if [ "${N_MULTI}" -ne 0 ]; then
        echo "ERROR [${sample}]: ${N_MULTI} multiallelic records remain"
        return 1
    fi


    ########################################################
    # STEP 8
    #
    # Save per-individual summary
    ########################################################

    {
        echo -e "sample\tn_het_biallelic_snp\tn_phased_het\tn_unphased_het\tn_nonhet"
        echo -e "${sample}\t${N_HET}\t${N_PHASED}\t${N_UNPHASED}\t${N_NONHET}"
    } > "${SUMMARY}"


    ########################################################
    # STEP 9
    #
    # Print first 20 final variants for inspection
    ########################################################

    echo "============================================================"
    echo "[${sample}] First 20 final heterozygous SNPs"
    echo "============================================================"

    "${BCFTOOLS}" query \
        -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' \
        "${OUT}" \
        | head -n 20


    ########################################################
    # STEP 10
    #
    # Remove temporary files
    ########################################################

    rm -f "${TMP_SUBSET}"
    rm -f "${TMP_SUBSET}.csi"

    rm -f "${TMP_HET}"
    rm -f "${TMP_HET}.csi"

    rm -f "${GT_QC}"


    ########################################################
    # Finish
    ########################################################

    echo "============================================================"
    echo "FINISHED sample = ${sample}"
    echo "Heterozygous biallelic SNPs = ${N_HET}"
    echo "Phased heterozygous SNPs = ${N_PHASED}"
    echo "Unphased heterozygous SNPs = ${N_UNPHASED}"
    echo "Non-heterozygous SNPs = ${N_NONHET}"
    echo "TIME = $(date)"
    echo "============================================================"
}


############################################################
# 14. Export function and variables
############################################################

export -f process_sample

export VCF
export REF

export VCFDIR
export SUMMARY_DIR
export TMPDIR

export BCFTOOLS

export THREADS_PER_SAMPLE

export PATH
export LD_LIBRARY_PATH


############################################################
# 15. Run WGS individuals in parallel
############################################################

echo "============================================================"
echo "Start xargs parallel processing"
echo "Array = ${SLURM_ARRAY_TASK_ID}"
echo "Samples = ${START}-${END}"
echo "MAX_PARALLEL = ${MAX_PARALLEL}"
echo "============================================================"

set +e

bash -c '
        set -o pipefail

        cat "'"${TASK_SAMPLE_LIST}"'" \
            | xargs \
                -r \
                -n 1 \
                -P '"${MAX_PARALLEL}"' \
                bash -c '"'"'
                    set -euo pipefail

                    process_sample "$1"
                '"'"' _
    '

XARGS_STATUS=$?

set -e


############################################################
# 16. Validate current array task outputs
############################################################

echo "============================================================"
echo "Validate array ${SLURM_ARRAY_TASK_ID}"
echo "============================================================"

FAILED=0

while read -r sample; do
    [ -z "${sample}" ] && continue

    OUT=${VCFDIR}/${sample}.het.biallelic.vcf.gz
    SUMMARY=${SUMMARY_DIR}/${sample}.tsv

    if [ ! -s "${OUT}" ] || \
       [ ! -s "${OUT}.tbi" ] || \
       [ ! -s "${SUMMARY}" ]; then

        echo "FAILED: ${sample}"

        FAILED=$((FAILED + 1))

        continue
    fi

    N=$(
        "${BCFTOOLS}" index \
            -n "${OUT}"
    )

    N_BAD=$(
        "${BCFTOOLS}" query \
            -f '[%GT\n]' \
            "${OUT}" \
            | awk '
                $1 != "0/1" &&
                $1 != "1/0" &&
                $1 != "0|1" &&
                $1 != "1|0"
            ' \
            | wc -l
    )

    if [ "${N_BAD}" -ne 0 ]; then
        echo -e "FAILED\t${sample}\t${N}\tnonhet=${N_BAD}"

        FAILED=$((FAILED + 1))
    else
        echo -e "SUCCESS\t${sample}\t${N}\tnonhet=0"
    fi

done < "${TASK_SAMPLE_LIST}"


############################################################
# 17. Clean task-level temporary files
############################################################

rm -f "${TASK_SAMPLE_LIST}"
rm -f "${SOURCE_SAMPLE_LIST}"


############################################################
# 18. Final status
############################################################

echo "============================================================"
echo "Array task = ${SLURM_ARRAY_TASK_ID}"
echo "Samples = ${N_TASK}"
echo "xargs status = ${XARGS_STATUS}"
echo "Failed = ${FAILED}"
echo "TIME = $(date)"
echo "============================================================"

if [ "${XARGS_STATUS}" -ne 0 ] || \
   [ "${FAILED}" -gt 0 ]; then

    echo "ERROR: array task ${SLURM_ARRAY_TASK_ID} failed"

    exit 1
fi

echo "06 array task ${SLURM_ARRAY_TASK_ID} successfully completed"
