# Genotypes and allelic effects

The 3 principal analysis entries are listed below. Grouped entries link to their component scripts; other files provide supporting functions.

| Analysis | Main code |
|---|---|
| WGS variant calling | [wgs_variant_calling.sh](01_WGS/wgs_variant_calling.sh) |
| ASE and ASTS analysis | [call_ase.sh](02_LORALS/call_ase.sh); [call_asts.sh](02_LORALS/call_asts.sh) |
| Allelic-effect and QTL association | [test_allelic_qtl_associations.R](03_QTL_correspondence/test_allelic_qtl_associations.R) |

For programs that locate files relative to their script, prepare the working directory with `tools/prepare_workspace.py` and supply the [required inputs](../../docs/INPUTS.md). Install the external programs in the environment for the selected analysis.

`wgs_variant_calling.sh` consolidates `trim`, `align`, `deduplicate`, `call`, `combine` and `genotype`. Supply sample IDs and paths as arguments; `combine` reads one gVCF path per line. LORALS array wrappers still require their declared batch/array variables; cluster submission policies are configured by the user.

## Reference and sample preparation

`01_prepare_reference.py` builds the reference, and `02_make_manifest.py` matches ONT and WGS samples. `prepare_common.sh` runs these preparation steps. `lorals_datasets.json` contains the study dataset subdirectories, filename modes and tissue aliases; change this configuration for a different input layout rather than editing the matching algorithm.

```bash
python script/06_ASE_ASTS/02_LORALS/01_prepare_reference.py "{study_gtf}" "{reference_fasta}" "{reference_output_dir}"
python script/06_ASE_ASTS/02_LORALS/02_make_manifest.py "{ont_root}" "{wgs_sample_list}" "{all_samples.tsv}" "{matched_samples.tsv}" --config "{dataset_config.json}"
```

The allelic pipeline starts from the joint-genotyped VCF and uses BCFtools to select biallelic heterozygous SNPs for allele-specific tests and QTL correspondence.

## Select tissues and samples

`extract_allelic_sites.py` accepts optional `--tissues` and `--samples` TSVs and `--shards`/`--shard-index`; omission selects all manifest samples in one shard. `aggregate_requests.py` accepts the same selections and requires `--output`. `build_analysis_tables.py` accepts `--tissues`, `--output` and an optional alias JSON. `extract_nominal.py` selects a task by `--tissue` and `--qtl`. Supplied selection labels must exist in the input metadata. Subset runs should use distinct output directories.
