# Input layout

Embedded analysis paths use `/path/to/chicken-aging` before configuration. `prepare_workspace.py` replaces that prefix with an absolute working directory and stages the selected code in its expected analysis locations.

| Working directory | Required resources |
|---|---|
| `ont/` | Tissue/group count and TPM tables, aligned reads, study GTF/database and SUPPA files. |
| `reference/` | GRCg6a sequence, annotations, projection references and protein/domain databases. |
| `analysis/` | Module inputs/outputs such as `specificity`, `expression_models`, `four_age`, `count_DTU`, `isoform_PSI`, `FIM_usage`, `conservation_programs` and `allelic_QTL`. |
| `model/experiment33/` | Coverage tracks, reference/exon tables, TFRecords, channel order, masks and model checkpoints. |
| `model/software/` | Baskerville checkout at the recorded commit, plus the Python environment referenced by the TFRecord wrapper. |
| `allelic/` and `WGS/` | WGS/ONT alignments, VCFs, LORALS input manifest and allele-assigned counts. |
| `qtl/` | Expression/usage phenotype packages, genotype/covariate files, full nominal statistics, fine-mapping manifests and GWAS inputs. |
| `GWAS/` | Chicken GWAS summary statistics, transcript interval sets and circular-shift enrichment outputs. |
| `LDSC/` | Human gene lists, coordinates, summary statistics, SNP/trait lists, LD-score reference panels and annotation outputs. |
| `software/` | External tool installations or symbolic links to the installations. |

Supply sequencing data, genotypes, phenotype matrices and model checkpoints separately. Study metadata are provided in `config/recorded_inputs/`, and manifest templates in `config/templates/`. Pass the grouped IsoQuant YAML to the quantification entry point.

Use consistent animal identifiers across tissues and distinct individuals across natural ages. Match sample labels between tables and manifests. Set the array/batch variables declared by the LORALS wrappers before submitting jobs.

The module README files describe the prepared tables consumed by each analysis. Keep reference versions, sample identifiers and transcript identifiers consistent across stages.
