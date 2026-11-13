# Molecular QTLs and trait integration

The 6 principal analysis entries are listed below. Grouped entries link to their component scripts; other files provide supporting functions.

| Analysis | Main code |
|---|---|
| iso-eQTL mapping | [run_transcript_eqtl_standard.sh](01_Phenotypes_and_mapping/transcript_eqtl_standard_pipeline/scripts/run_transcript_eqtl_standard.sh) |
| iso-uQTL preprocessing and mapping | [run_iso_uqtl.py](01_Phenotypes_and_mapping/iso_uqtl_mapping/run_iso_uqtl.py) |
| SuSiE fine-mapping | [run_susie.py](02_Fine_mapping/run_susie.py) |
| Conditional QTLs | [map_conditional_qtl.py](03_Conditional/map_conditional_qtl.py) |
| Tissue sharing | [tissue_sharing.R](04_Tissue_sharing/tissue_sharing.R) |
| GWAS colocalization and SMR | [run_coloc.R](05_GWAS_coloc_SMR/run_coloc.R); [run_smr.sh](05_GWAS_coloc_SMR/run_smr.sh) |

## Parameterized invocation

Replace each `{...}` placeholder with your input value.

```bash
python script/07_Molecular_QTL/02_Fine_mapping/run_susie.py --manifest "{input_manifest}" --tissue "{tissue}" --qtl "{teQTL_or_tuQTL}" --output "{output_dir}"
```

For programs that locate files relative to their script, prepare the working directory with `tools/prepare_workspace.py` and supply the [required inputs](../../docs/INPUTS.md). Install the external programs in the environment for the selected analysis.

`run_conditional_repair.py` and `empirical_fallback.py` provide alternative handling of non-finite P values in conditional analysis. Generate the input manifest for the dataset and analysis branch being run.

## Colocalization and SMR inputs

Use header-only templates in `config/templates/` to supply your study records. Paths inside the SMR input and BESD manifests are absolute or relative to the manifest file. Braced values below are placeholders.

```bash
python script/07_Molecular_QTL/05_GWAS_coloc_SMR/prepare_smr_input.py --manifest "{smr_tissues.tsv}" --output-dir "{prepared_input_dir}"
bash script/07_Molecular_QTL/05_GWAS_coloc_SMR/build_besd.sh "{besd_tissues.tsv}" "{prepared_input_dir}" "{besd_output_dir}" "{smr_executable}"
python script/07_Molecular_QTL/05_GWAS_coloc_SMR/extract_candidate_qtl.py --manifest "{coloc_candidates.tsv}" --input-dir "{prepared_input_dir}" --output-dir "{candidate_output_dir}"
Rscript script/07_Molecular_QTL/05_GWAS_coloc_SMR/run_coloc.R "{candidate_manifest.tsv}" "{matched_snp_dir}" "{output_dir}"
bash script/07_Molecular_QTL/05_GWAS_coloc_SMR/run_smr.sh "{smr_executable}" "{plink_prefix}" "{gwas_summary}" "{besd_prefix}" "{output_prefix}" "{threads}"
```

The SMR runner consumes prepared BESD files and a harmonized, uncompressed eight-column GWAS summary (`SNP A1 A2 freq b se p n`). BESD preparation validates the seven-column nominal schema produced by `prepare_smr_input.py`. Use distinct output directories for distinct input versions. Candidate IDs and sample sizes are supplied by the manifest; no particular tissue or candidate locus is embedded in the extraction code.

Colocalization requires allele-harmonized SNP tables for the selected molecular QTL and GWAS datasets, together with the candidate manifest. Supply the gene-eQTL comparison inputs required by the script.

## Fine-mapping input selection

Supply one record per tissue/QTL type using `config/templates/fine_mapping_inputs.tsv`. Genotype, covariate and phenotype paths are explicit, so dataset-specific exclusions are represented by the supplied input versions. The original sample-order and file checks are retained; the recorded chromosome count is 32 unless explicitly changed.

```bash
python script/07_Molecular_QTL/02_Fine_mapping/prepare_manifest.py --manifest "{fine_mapping_inputs.tsv}" --output "{fine_mapping_dir}/input_manifest.json" --expected-chromosomes 32
```

`study_tissues.json` records the conditional-analysis tissue lists; `input_partitions.json` records the credible-set input partition. These are study metadata and can be changed for another dataset without editing the analysis code.

Tissue sharing reads `04_Tissue_sharing/tissue_metadata.tsv` as its command-line argument. The included table preserves the original tissue membership and aliases; use a separate verified table for another comparison.

```bash
Rscript script/07_Molecular_QTL/04_Tissue_sharing/tissue_sharing.R "{tissue_metadata.tsv}"
```

## iso-uQTL pipeline

`01_Phenotypes_and_mapping/iso_uqtl_mapping/run_iso_uqtl.py` runs phenotype preprocessing, sample checks, chromosome-wise nominal tests and permutation/FDR analysis. Supply input paths in `config/templates/iso_uqtl_inputs.tsv`, either absolute or relative to the manifest. Select `phenotype_filter=variance` or `none` for each dataset.

The variance filter requires variance ≥ 1e-4, IQR ≥ 0.01 and nonzero fraction ≥ 0.2. TensorQTL uses MAF 0.05, and `post.calculate_qvalues` applies FDR 0.05. Nominal mapping processes all input chromosomes unless an explicit chromosome selection is supplied.

```bash
python script/07_Molecular_QTL/01_Phenotypes_and_mapping/iso_uqtl_mapping/run_iso_uqtl.py --manifest "{iso_uqtl_inputs.tsv}" --dry-run
python script/07_Molecular_QTL/01_Phenotypes_and_mapping/iso_uqtl_mapping/run_iso_uqtl.py --manifest "{iso_uqtl_inputs.tsv}" --tissue "{tissue}"
```

The preprocessing command uses the active environment unless a conda environment is explicitly requested. It accepts explicit genotype prefixes and transcript-TPM files; inputs are read in place and outputs go to the manifest's output directory.
