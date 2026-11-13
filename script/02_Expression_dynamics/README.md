# Expression specificity and aging dynamics

The 8 principal analysis entries are listed below. Grouped entries link to their component scripts; other files provide supporting functions.

| Analysis | Main code |
|---|---|
| Requantification with the study annotation | [requantify_group.sh](01_Quantification/requantify_group.sh) |
| Specificity indices | [run_tspex.R](02_Specificity/run_tspex.R) |
| Transcriptomic timing | [timing_score.R](03_Timing/timing_score.R) |
| Age-related expression model | [fit_age_expression.R](04_Age_models/fit_age_expression.R) |
| Isoform dynamics | [test_isoform_dynamics.py](05_Isoform_dynamics/test_isoform_dynamics.py) |
| Within-tissue and pooled trajectory clustering | Within each tissue (15 clusters): [within_tissue_zscore.R](06_Clustering/within_tissue_zscore.R) → [within_tissue_loess.R](06_Clustering/within_tissue_loess.R); pooled across tissues (10 clusters): [pooled_kmeans.py](06_Clustering/pooled_kmeans.py) |
| GO enrichment | [enrich_cluster_go.R](07_GO_enrichment/enrich_cluster_go.R) |
| Motif enrichment | [run_ame.sh](08_Motif/run_ame.sh) |

## Parameterized invocation

Replace each `{...}` placeholder with your input value.

```bash
bash script/02_Expression_dynamics/01_Quantification/requantify_group.sh "{tissue}" "{input_yaml}" "{genome_fasta}" "{study_gtf_or_db}" "{output_root}"
Rscript /scratch/chicken-aging/analysis/expression_models/fit_age_expression.R "{analysis_dir}" "{tissue}"
```

For programs that locate files relative to their script, prepare the working directory with `tools/prepare_workspace.py` and supply the [required inputs](../../docs/INPUTS.md). Install the external programs in the environment for the selected analysis.

`05_Isoform_dynamics/stage_expression_inputs.py` prepares the four-age input tables from `config/templates/expression_inputs.tsv`. Input paths are absolute or relative to the manifest; destination paths are relative to the analysis directory. Existing destination files are preserved.

```bash
python script/02_Expression_dynamics/05_Isoform_dynamics/stage_expression_inputs.py --manifest "{expression_inputs.tsv}" --output-dir "{analysis_dir}"
```

`compact_background.py` exports expression-background classifications. Optional `--expected-common` and `--expected-parent-undetected` arguments check expected record counts.

## Explicit sample order and GO inputs

Both single-tissue clustering scripts take `MATRIX_TSV SAMPLE_MANIFEST_TSV TISSUE OUTPUT_TSV`; the manifest has `sample`, `tissue` and `age` columns. Expression columns are matched by sample ID, with five replicates at each of 20, 50, 80 and 100 weeks. Within-tissue analysis uses LOESS and complete-linkage hierarchical clustering with 15 clusters; pooled analysis uses K-means with 10 clusters.

`enrich_cluster_go.R INPUT_DIR OUTPUT_DIR [EXPECTED_RECORDS]` reads five prepared tables for GO enrichment. `EXPECTED_RECORDS` optionally checks the input record count.
