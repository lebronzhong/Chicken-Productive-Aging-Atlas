# LDSC annotation inputs

The annotation manifest contains 165 tissue/group records: 33 tissues for each of the four natural ages (20, 50, 80 and 100 weeks) and the FIM group. Each record links prepared t-statistic tables, selected transcript lists, chicken-to-human mapping tables and human Ensembl gene lists. `CR` in input filenames denotes the FIM group sampled at 100 weeks.

`build_ldsc_annotations.py` reads the prepared gene sets and calls LDSC annotation generation with a 100,000-bp window. Supply human gene coordinates and BIM files using the same genome assembly and chromosome naming convention. Prepare the t-statistics, transcript selection and orthologue mappings before this annotation step.

```bash
python script/08_Comparative_and_complex_traits/04_Human_LDSC/build_ldsc_annotations.py --manifest "{annotation_manifest.tsv}" --ldsc-dir "{ldsc_dir}" --gene-coordinates "{human_gene_coordinates}" --bim-pattern "{reference_dir}/1000G.EUR.QC.{chromosome}.bim" --chromosomes "{chromosomes.tsv}" --output-dir "{annotation_dir}"
```

Use `tools/prepare_workspace.py` to configure the project-root placeholder. For another dataset, supply a corresponding annotation manifest and chromosome list.
