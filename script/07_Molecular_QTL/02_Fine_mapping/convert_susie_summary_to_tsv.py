#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 07_Molecular_QTL
# Purpose: Convert susie summary to tsv.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

from pathlib import Path
import sys
import pandas as pd

if len(sys.argv) != 2:
    raise SystemExit("Usage: python3 convert_susie_summary_to_tsv.py <result_root>")
result_root = Path(sys.argv[1])
files = sorted(result_root.rglob("*.SuSiE_summary.parquet"))
if not files:
    raise SystemExit(f"No SuSiE summary parquet files found under: {result_root}")
print("Files found:", len(files))
for parquet_file in files:
    output_file = parquet_file.with_suffix(".tsv")
    print("Reading :", parquet_file)
    print("Writing :", output_file)
    df = pd.read_parquet(parquet_file)
    df.to_csv(output_file, sep="\t", index=False)
    n_phenotypes = df["phenotype_id"].nunique()
    n_credible_sets = df[["phenotype_id", "cs_id"]].drop_duplicates().shape[0]
    print("  rows           :", len(df))
    print("  phenotypes     :", n_phenotypes)
    print("  credible sets  :", n_credible_sets)
print("Done.")
