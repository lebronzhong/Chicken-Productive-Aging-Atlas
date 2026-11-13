#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 02_Expression_dynamics
# Purpose: Compact background.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

import pathlib
import gzip
import pandas as pd

import argparse

parser = argparse.ArgumentParser(description="Export the compact transcript classification table.")
parser.add_argument("--root", type=pathlib.Path, default=pathlib.Path(__file__).resolve().parent)
parser.add_argument("--expected-common", type=int)
parser.add_argument("--expected-parent-undetected", type=int)
args = parser.parse_args()
out = args.root
columns = [
    "tissue",
    "transcript_id",
    "gene_id",
    "symbol",
    "class",
    "age_q",
    "max_abs_lfc",
    "parent_age_q",
    "parent_max_abs_lfc",
    "common_background",
    "age_associated",
    "parent_not_detected",
    "all_age_supported",
    "structural_quality",
    "strict_quality",
]
totals = [0, 0, 0]
with gzip.open(
    out / "02_overall/transcript_classification.tsv.gz", "wt", encoding="utf-8", compresslevel=1
) as handle:
    for chunk_index, frame in enumerate(
        pd.read_csv(
            out / "02_overall/transcript_parent_four_age_results.tsv.gz",
            sep="\t",
            usecols=columns,
            chunksize=50000,
        )
    ):
        frame[columns].to_csv(handle, sep="\t", index=False, header=chunk_index == 0)
        totals[0] += len(frame)
        totals[1] += int(frame.common_background.sum())
        totals[2] += int(frame.parent_not_detected.sum())
if args.expected_common is not None:
    assert (
        totals[1] == args.expected_common
    ), "Common-background count differs from the expected value."
if args.expected_parent_undetected is not None:
    assert (
        totals[2] == args.expected_parent_undetected
    ), "Parent-undetected count differs from the expected value."
print("Compact transcript classification", totals, flush=True)
