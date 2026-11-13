# Study input metadata

This directory contains reference identities, sample metadata and analysis settings:

- `reference_assets.tsv` and `reference_mapping_runs.json`: reference versions and mapping parameters.
- `short_read_samples.tsv`: 21 matched short-read libraries.
- `model_records.json` and `model_tissue_counts.tsv`: model settings and the 33-tissue sample counts; use the model module's `study_samples.json` for tissue order and identifiers.
- `ame_recorded_parameters.json`: motif database and AME settings.
- `ldsc_age_annotations.tsv`: 165 tissue/group annotation records.
- `human_autosomes.tsv`: chromosome list for human autosomal LD-score analyses.

Supply the corresponding data files and configure paths for the local analysis directory. Input manifests for additional datasets can be prepared from `config/templates/`.
