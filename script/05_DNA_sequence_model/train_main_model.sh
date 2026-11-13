#!/usr/bin/env bash
# Chicken Productive Aging Atlas | 05_DNA_sequence_model
# Purpose: Train the 33-tissue model.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

set -euo pipefail
export BORZOI_ROOT=/path/to/chicken-aging/model
export PYTHONHASHSEED=20261013
cd "$BORZOI_ROOT/experiment33"
test -f gpu_scripts/params_medium_33.json
test -f data/tfrecords_v2/statistics.json
test ! -e runs/medium_33_full_v1
python -u gpu_scripts/train_seeded.py \
  -o runs/medium_33_full_v1 -l runs/medium_33_full_v1_tensorboard \
  gpu_scripts/params_medium_33.json data/tfrecords_v2
