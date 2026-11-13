# Induced forced molting

The 2 principal analysis entries are listed below. Grouped entries link to their component scripts; other files provide supporting functions.

| Analysis | Main code |
|---|---|
| FIM expression changes and recovery | [fit_fim_expression.R](01_Expression_recovery/fit_fim_expression.R); [expression_recovery_index.py](01_Expression_recovery/expression_recovery_index.py) |
| Usage recovery | [fit_fim_usage.R](02_Usage_recovery/fit_fim_usage.R) |

For programs that locate files relative to their script, prepare the working directory with `tools/prepare_workspace.py` and supply the [required inputs](../../docs/INPUTS.md). Install the external programs in the environment for the selected analysis.
