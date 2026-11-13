#!/usr/bin/env bash
set -euo pipefail

usage() {
cat <<EOF
Usage:
  bash build_tuQTL_psi_phenotype.sh \
    --tissues "{tissue}" \
    --transcript-tpm-path <read_only_input_dir> \
    --output-path <output_dir> \
    --gene-expr-path <gene_expr_dir> \
    --genotype-path <dir> \
    --salmon-out <dir> \
    --tx2gene <tx2gene_union.csv> \
    --gtf <union.topo.gtf> \
    --suppa <suppa.py> \
    [--script-dir <dir>] \
    [--genotype-prefix <prefix>] \
    [--transcript-tpm-file <file>] \
    [--conda-env <environment>]
EOF
exit 1
}

TISSUES=""
TRANSCRIPT_TPM_PATH=""
OUTPUT_PATH=""
GENE_EXPR_PATH=""
GENOTYPE_PATH=""
SALMON_OUT=""
TX2GENE=""
GTF=""
SUPPA=""
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONDA_ENV=""
GENOTYPE_PREFIX=""
TPM_FILE=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --tissues) TISSUES="$2"; shift 2 ;;
    --genotype-prefix) GENOTYPE_PREFIX="$2"; shift 2 ;;
    --transcript-tpm-file) TPM_FILE="$2"; shift 2 ;;
    --transcript-tpm-path) TRANSCRIPT_TPM_PATH="$2"; shift 2 ;;
    --output-path) OUTPUT_PATH="$2"; shift 2 ;;
    --gene-expr-path) GENE_EXPR_PATH="$2"; shift 2 ;;
    --genotype-path) GENOTYPE_PATH="$2"; shift 2 ;;
    --salmon-out) SALMON_OUT="$2"; shift 2 ;;
    --tx2gene) TX2GENE="$2"; shift 2 ;;
    --gtf) GTF="$2"; shift 2 ;;
    --suppa) SUPPA="$2"; shift 2 ;;
    --script-dir) SCRIPT_DIR="$2"; shift 2 ;;
    --conda-env) CONDA_ENV="$2"; shift 2 ;;
    -h|--help) usage ;;
    *) echo "[ERROR] Unknown argument: $1" >&2; usage ;;
  esac
done

[[ -z "$TISSUES" ]] && { echo "[ERROR] --tissues is required"; usage; }
[[ -z "$TRANSCRIPT_TPM_PATH" ]] && { echo "[ERROR] --transcript-tpm-path is required"; usage; }
[[ -z "$OUTPUT_PATH" ]] && { echo "[ERROR] --output-path is required"; usage; }
[[ -z "$GENE_EXPR_PATH" ]] && { echo "[ERROR] --gene-expr-path is required"; usage; }
[[ -z "$GENOTYPE_PATH" ]] && { echo "[ERROR] --genotype-path is required"; usage; }
[[ -z "$SALMON_OUT" ]] && { echo "[ERROR] --salmon-out is required"; usage; }
[[ -z "$TX2GENE" ]] && { echo "[ERROR] --tx2gene is required"; usage; }
[[ -z "$GTF" ]] && { echo "[ERROR] --gtf is required"; usage; }
[[ -z "$SUPPA" ]] && { echo "[ERROR] --suppa is required"; usage; }

mkdir -p "$OUTPUT_PATH" "$GENE_EXPR_PATH"
touch "${OUTPUT_PATH}/.write_test" && rm -f "${OUTPUT_PATH}/.write_test"
touch "${GENE_EXPR_PATH}/.write_test" && rm -f "${GENE_EXPR_PATH}/.write_test"

if [[ -n "$CONDA_ENV" ]]; then
  command -v conda >/dev/null 2>&1 || { echo "Requested conda environment but conda is unavailable" >&2; exit 1; }
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate "$CONDA_ENV"
else
  echo "[INFO] Using the active environment."
fi

read -r -a tissue_array <<< "$TISSUES"
if [[ ( -n "$GENOTYPE_PREFIX" || -n "$TPM_FILE" ) && ${#tissue_array[@]} -ne 1 ]]; then
  echo "Explicit genotype/TPM files require exactly one tissue per invocation" >&2
  exit 1
fi
for tissue in "${tissue_array[@]}"; do
  echo "============================================================"
  echo "[INFO] Processing tissue: ${tissue}"
  echo "============================================================"

  # Input files
  bed="${TPM_FILE:-${TRANSCRIPT_TPM_PATH}/salmon_pheno.ChickenGTEx.${tissue}.TPM.bed.gz}"
  fam="${GENOTYPE_PREFIX:-${GENOTYPE_PATH}/ChickenGTEx.${tissue}}.fam"
  bim="${GENOTYPE_PREFIX:-${GENOTYPE_PATH}/ChickenGTEx.${tissue}}.bim"

  # Outputs in the supplied directories
  suppa_infile="${OUTPUT_PATH}/salmon_pheno.ChickenGTEx.${tissue}.TPM.suppa.txt"
  suppa_outprefix="${OUTPUT_PATH}/salmon_pheno.ChickenGTEx.${tissue}.TPM"
  psi_raw="${OUTPUT_PATH}/salmon_pheno.ChickenGTEx.${tissue}.TPM_isoform.psi"

  gene_tpm_file="${GENE_EXPR_PATH}/ChickenGTEx.${tissue}.gene_tpm.tsv.gz"
  gene_tpm_out="${GENE_EXPR_PATH}/ChickenGTEx.${tissue}.gene_tpm.filtered.tsv.gz"

  tx_bed_out="${OUTPUT_PATH}/salmon_pheno.ChickenGTEx.${tissue}.TPM.filtered.bed.gz"
  psi_filtered="${OUTPUT_PATH}/salmon_pheno.ChickenGTEx.${tissue}.TPM_isoform.filtered.psi"

  tensorqtl_bed="${OUTPUT_PATH}/salmon_pheno.ChickenGTEx.${tissue}.TPM_isoform.filtered.tensorqtl.bed.gz"
  tensorqtl_bed_bimchr="${OUTPUT_PATH}/salmon_pheno.ChickenGTEx.${tissue}.TPM_isoform.filtered.bimchr.tensorqtl.bed.gz"

  [[ -f "$bed" ]] || { echo "[ERROR] Missing transcript TPM BED: $bed" >&2; exit 1; }
  [[ -f "$fam" ]] || { echo "[ERROR] Missing FAM file: $fam" >&2; exit 1; }
  [[ -f "$bim" ]] || { echo "[ERROR] Missing BIM file: $bim" >&2; exit 1; }
  [[ -f "$TX2GENE" ]] || { echo "[ERROR] Missing tx2gene file: $TX2GENE" >&2; exit 1; }
  [[ -f "$GTF" ]] || { echo "[ERROR] Missing GTF file: $GTF" >&2; exit 1; }
  [[ -f "$SUPPA" ]] || { echo "[ERROR] Missing SUPPA script: $SUPPA" >&2; exit 1; }
  [[ -f "${SCRIPT_DIR}/transcript_to_gene.R" ]] || { echo "[ERROR] Missing ${SCRIPT_DIR}/transcript_to_gene.R" >&2; exit 1; }
  [[ -f "${SCRIPT_DIR}/filter_transcript_usage.R" ]] || { echo "[ERROR] Missing ${SCRIPT_DIR}/filter_transcript_usage.R" >&2; exit 1; }
  [[ -f "${SCRIPT_DIR}/usage_to_qtl_bed.R" ]] || { echo "[ERROR] Missing ${SCRIPT_DIR}/usage_to_qtl_bed.R" >&2; exit 1; }

  echo "[1/6] Build SUPPA expression input"
  zcat "$bed" | awk 'BEGIN{FS=OFS="\t"}
  NR==1{
      for(i=5;i<=NF;i++){
          printf $i;
          if(i<NF) printf OFS; else printf "\n";
      }
      next
  }
  {
      printf $4;
      for(i=5;i<=NF;i++) printf OFS $i;
      printf "\n";
  }' > "$suppa_infile"

  echo "[2/6] Run SUPPA psiPerIsoform"
  python "$SUPPA" psiPerIsoform     -g "$GTF"     -e "$suppa_infile"     -o "$suppa_outprefix"
  [[ -f "$psi_raw" ]] || { echo "[ERROR] SUPPA PSI output not found: $psi_raw" >&2; exit 1; }

  echo "[3/6] Generate gene-level TPM"
  Rscript "${SCRIPT_DIR}/transcript_to_gene.R"     "$fam"     "$SALMON_OUT"     "$TX2GENE"     "${GENE_EXPR_PATH}/ChickenGTEx.${tissue}"     no
  [[ -f "$gene_tpm_file" ]] || { echo "[ERROR] Gene TPM output not found: $gene_tpm_file" >&2; exit 1; }

  echo "[4/6] Jointly filter gene TPM, transcript TPM BED, and PSI"
  Rscript "${SCRIPT_DIR}/filter_transcript_usage.R"     "$gene_tpm_file"     "$TX2GENE"     "$bed"     "$psi_raw"     "$gene_tpm_out"     "$tx_bed_out"     "$psi_filtered"
  [[ -f "$gene_tpm_out" ]] || { echo "[ERROR] Filtered gene TPM not found: $gene_tpm_out" >&2; exit 1; }
  [[ -f "$tx_bed_out" ]] || { echo "[ERROR] Filtered transcript BED not found: $tx_bed_out" >&2; exit 1; }
  [[ -f "$psi_filtered" ]] || { echo "[ERROR] Filtered PSI not found: $psi_filtered" >&2; exit 1; }

  echo "[5/6] Convert filtered PSI to TensorQTL BED"
  Rscript "${SCRIPT_DIR}/usage_to_qtl_bed.R"     "$psi_filtered"     "$GTF"     "$tensorqtl_bed"
  [[ -f "$tensorqtl_bed" ]] || { echo "[ERROR] TensorQTL BED not found: $tensorqtl_bed" >&2; exit 1; }

  echo "[6/6] Keep chromosomes present in PLINK .bim"
  zcat "$tensorqtl_bed" | awk -v BIM="$bim" 'BEGIN{
      FS=OFS="\t";
      while ((getline line < BIM) > 0) {
          split(line, a, /[ \t]+/);
          chr=a[1];
          sub(/^chr/, "", chr);
          keep[chr]=1;
      }
      close(BIM);
  }
  NR==1{
      if($1=="chr") $1="#chr";
      print;
      next
  }
  {
      chr=$1;
      sub(/^chr/, "", chr);
      if (chr in keep) {
          $1=chr;
          print
      }
  }' | gzip -c > "$tensorqtl_bed_bimchr"

  echo "[QC] Before chromosome filtering:"
  zcat "$tensorqtl_bed" | tail -n +2 | wc -l
  echo "[QC] After .bim chromosome filtering:"
  zcat "$tensorqtl_bed_bimchr" | tail -n +2 | wc -l
  echo "[QC] Chromosomes in final BED:"
  zcat "$tensorqtl_bed_bimchr" | tail -n +2 | cut -f1 | sort -V | uniq -c

  echo "[DONE] ${tissue}"
  echo "[OUTPUT] Raw PSI: ${psi_raw}"
  echo "[OUTPUT] Filtered PSI: ${psi_filtered}"
  echo "[OUTPUT] Final TensorQTL BED for tuQTL: ${tensorqtl_bed_bimchr}"
  echo
done

echo "All tissues completed."
