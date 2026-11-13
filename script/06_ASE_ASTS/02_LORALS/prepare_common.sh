#!/usr/bin/env bash
# Chicken Productive Aging Atlas | 06_ASE_ASTS
# Purpose: Prepare common.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

set -euo pipefail

work=/path/to/chicken-aging/allelic

ase=${work}/ase
asts=${work}/asts
logs=${work}/logs
scripts=${work}/scripts
software=${work}/software

ONT=/path/to/chicken-aging/ont

VCF=/path/to/chicken-aging/WGS/CombineGVCF/Chicken_Aging_WGS.vcf.gz

REF=/path/to/chicken-aging/reference/galGal6/Gallus_gallus.GRCg6a.dna.toplevel.fa

GTF=/path/to/chicken-aging/ont/shell/10_20240903.gffcompare/0313_GTF_update/union.topo.gtf

REFDIR=${software}/lorals_ref
VCFDIR=${ase}/vcf

MANIFEST=${scripts}/lorals_manifest.tsv
MATCHED_MANIFEST=${scripts}/lorals_manifest_matched.tsv

# Load software
# Create directories
mkdir -p ${ase}
mkdir -p ${asts}
mkdir -p ${logs}
mkdir -p ${scripts}
mkdir -p ${software}
mkdir -p ${REFDIR}
mkdir -p ${VCFDIR}
mkdir -p ${VCFDIR}/summary

# Check software
echo "===== software ====="

which python
python --version

which minimap2
minimap2 --version

which calc_ase
which annotate_ase
which calc_asts

which samtools
samtools --version

which bcftools
bcftools --version | head -n 2

which bedtools
bedtools --version

# Reference FASTA index
if [ ! -f ${REF}.fai ]
then

    samtools faidx ${REF}

fi

# WGS VCF index
if [ ! -f ${VCF}.csi ] && [ ! -f ${VCF}.tbi ]
then

    bcftools index \
    --threads 20 \
    -f \
    ${VCF}

fi

# Prepare transcriptome reference
python \
${scripts}/01_prepare_reference.py \
${GTF} \
${REF} \
${REFDIR}

# Create dummy chicken blacklist/mappability BED
FIRST_CHR=$(head -n 1 ${REF}.fai | cut -f1)

FIRST_LEN=$(head -n 1 ${REF}.fai | cut -f2)

DUMMY_START=$((FIRST_LEN + 100000))

DUMMY_END=$((FIRST_LEN + 100001))

echo -e \
"${FIRST_CHR}\t${DUMMY_START}\t${DUMMY_END}\tDUMMY" \
> ${REFDIR}/no_filter.bed

# Get WGS sample names
bcftools query \
-l ${VCF} \
> ${scripts}/wgs_samples.txt

echo "===== WGS samples ====="

cat ${scripts}/wgs_samples.txt

# Build manifest
python \
${scripts}/02_make_manifest.py \
${ONT} \
${scripts}/wgs_samples.txt \
${MANIFEST} \
${MATCHED_MANIFEST}

echo "===== Manifest first 30 rows ====="

head -n 31 ${MANIFEST}

echo "===== Unmatched samples ====="

awk -F '\t' '
NR==1 || $8=="NA"
' ${MANIFEST}

echo "===== Matched manifest first 30 rows ====="

head -n 31 ${MATCHED_MANIFEST}

echo "===== Number of samples per tissue ====="

awk -F '\t' '
NR>1 {
    count[$4]++
}
END {
    for(tissue in count) {
        print tissue "\t" count[tissue]
    }
}
' ${MATCHED_MANIFEST} \
| sort

# Unique WGS samples required
awk -F '\t' '
NR>1 {
    print $8
}
' ${MATCHED_MANIFEST} \
| sort -u \
> ${scripts}/needed_vcf_samples.txt

echo "===== Required WGS samples ====="

cat ${scripts}/needed_vcf_samples.txt

# Count array sizes
N_VCF=$(wc -l \
< ${scripts}/needed_vcf_samples.txt)

N_ONT=$(awk '
END {
    print NR-1
}
' ${MATCHED_MANIFEST})

echo "================================"

echo "VCF array task number = ${N_VCF}"

echo "ASE/ASTS array task number = ${N_ONT}"

echo "================================"

echo "Preparation finished"
