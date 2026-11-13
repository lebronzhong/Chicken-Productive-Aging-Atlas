# Alternative splicing and transcript usage

The 4 principal analysis entries are listed below. Grouped entries link to their component scripts; other files provide supporting functions.

| Analysis | Main code |
|---|---|
| SUPPA events and isoform fractions | [suppa_analysis.sh](01_SUPPA/suppa_analysis.sh) |
| PSI statistics | [fit_original_PSI.R](03_Isoform_PSI/fit_original_PSI.R) |
| Count-based DTU | [fit_usage.R](02_Count_DTU/fit_usage.R) |
| Usage trajectories | [cluster_psi.py](04_Usage_clusters/cluster_psi.py) |

## Parameterized invocation

Replace each `{...}` placeholder with your input value.

```bash
SUPPA_SCRIPT="{installed_suppa.py}" bash script/03_Splicing_and_usage/01_SUPPA/suppa_analysis.sh events "{study_gtf}" "{output_prefix}"
SUPPA_SCRIPT="{installed_suppa.py}" bash script/03_Splicing_and_usage/01_SUPPA/suppa_analysis.sh psi "{study_gtf}" "{output_prefix}" "{tpm_matrix}"
```

For programs that locate files relative to their script, prepare the working directory with `tools/prepare_workspace.py` and supply the [required inputs](../../docs/INPUTS.md). Install the external programs in the environment for the selected analysis.

SUPPA generates separate files for the selected splicing-event classes. Event PSI and within-gene transcript fractions are analyzed separately.
