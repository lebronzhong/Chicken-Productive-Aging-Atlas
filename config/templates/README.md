# Input manifest schemas

These tab-separated templates define the required columns. Add one row per sample, tissue, trait or candidate, as appropriate.

- `smr_tissues.tsv`: one tissue per row; explicit FDR file and nominal-Parquet/phenotype-BED patterns.
- `besd_tissues.tsv`: one tissue and the original QTL BIM path per row; exactly these two columns.
- `coloc_candidates.tsv`: unique candidate ID, tissue, transcript, gene, trait identity and positive integer QTL sample size.
- `expression_inputs.tsv`: source table and its destination relative to the four-age analysis directory.

Use consistent tissue and transcript identifiers across manifests and provide the sample size for each analysis.

- `fine_mapping_inputs.tsv`: explicit FDR file, genotype prefix, covariates and phenotype directory for each tissue/QTL type. Use legacy `teQTL`/`tuQTL` file labels in this manifest.
