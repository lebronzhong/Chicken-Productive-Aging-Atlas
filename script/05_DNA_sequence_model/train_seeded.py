#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 05_DNA_sequence_model
# Purpose: Train seeded.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

import os, runpy, sys
import tensorflow as tf

seed = 20261013
tf.keras.utils.set_random_seed(seed)
root = os.environ["BORZOI_ROOT"]
script = root + "/software/baskerville-tf/src/baskerville/scripts/hound_train.py"
print("FULL TRAINING; random seed:", seed, "; precision: float32", flush=True)
sys.argv = [script] + sys.argv[1:]
runpy.run_path(script, run_name="__main__")
