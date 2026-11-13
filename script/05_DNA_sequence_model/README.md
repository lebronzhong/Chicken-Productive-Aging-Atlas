# DNA sequence model

The 3 principal analysis entries are listed below. Grouped entries link to their component scripts; other files provide supporting functions.

| Analysis | Main code |
|---|---|
| Model training data preparation | [prepare_model_inputs.py](prepare_model_inputs.py); [prepare_tfrecords_33.py](prepare_tfrecords_33.py) |
| Train the 33-tissue model | [train_main_model.sh](train_main_model.sh) |
| Evaluate coverage and end usage | [evaluate_model_33.py](evaluate_model_33.py) |

For programs that locate files relative to their script, prepare the working directory with `tools/prepare_workspace.py` and supply the [required inputs](../../docs/INPUTS.md). Install the external programs in the environment for the selected analysis.

`prepare_model_inputs.py` checks inputs, prepares coverage tracks and assembles training data for the 33-tissue model. Supply the coverage inputs, sequence windows and chromosome splits specified by the preparation scripts.

`study_samples.json` records the original tissue/channel order, tissue aliases and animal labels used during input preparation. Preserve this order when reproducing the 33-tissue model.
