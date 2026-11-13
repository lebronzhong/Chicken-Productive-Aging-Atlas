# Comparative analyses and complex traits

The 4 principal analysis entries are listed below. Grouped entries link to their component scripts; other files provide supporting functions.

| Analysis | Main code |
|---|---|
| Mammalian aging signatures | [gradient_01_genes.py](01_Mammalian_signatures/gradient_01_genes.py) |
| Functional programmes | [compare_functional_programs.py](02_Functional_programs/compare_functional_programs.py) |
| Chicken trait enrichment of expression and usage sets | Expression-defined sets: [gwas_enrichment_expression.R](03_Chicken_GWAS/gwas_enrichment_expression.R); usage-defined sets: [gwas_enrichment_usage.R](03_Chicken_GWAS/gwas_enrichment_usage.R) |
| Human trait heritability enrichment | [run_ldsc.sh](04_Human_LDSC/run_ldsc.sh) |

## Parameterized invocation

Replace each `{...}` placeholder with your input value.

```bash
Rscript /scratch/chicken-aging/GWAS/gwas_enrichment_expression.R "{trait_manifest.tsv}" "{expression_intervals.bed}" "{output_dir}" "{geneset_helper.R}"
Rscript /scratch/chicken-aging/GWAS/gwas_enrichment_usage.R "{trait_manifest.tsv}" "{usage_intervals.bed}" "{output_dir}" "{geneset_helper.R}"
```

For programs that locate files relative to their script, prepare the working directory with `tools/prepare_workspace.py` and supply the [required inputs](../../docs/INPUTS.md). Install the external programs in the environment for the selected analysis.

`04_Human_LDSC/munge_sumstats.sh` prepares GWAS statistics using the AF1 frequency column, chunk size 500000 and a HapMap3 allele list. Result collection takes a trait list, result directory and output prefix.

```bash
bash script/08_Comparative_and_complex_traits/04_Human_LDSC/munge_sumstats.sh "{ldsc_dir}" "{gwas_summary}" "{hm3_allele_list}" "{output_prefix}"
Rscript script/08_Comparative_and_complex_traits/04_Human_LDSC/collect_ldsc_results.R "{trait_list}" "{results_dir}" "{output_prefix}"
```

## Trait inputs and LDSC annotations

The expression and usage GWAS-enrichment entries both take `TRAIT_MANIFEST INTERVAL_BED OUTPUT_DIR GENESET_HELPER [N_PERM]`. Use `config/templates/gwas_traits.tsv`; chromosome, position, P-value and effect columns are named explicitly for each trait. All manifest rows are processed, with original SNP coordinates, sum-of-squared-effects statistic, default 1000 permutations and within-result FDR retained. Results use trait IDs rather than file-list positions.

The LDSC annotation builder uses prepared tissue/group gene sets and human gene coordinates. The [LDSC input guide](../../docs/LDSC_UPSTREAM.md) describes the 165-record annotation manifest and the annotation command.
