#!/usr/bin/env python3
import argparse
import gzip
import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--min-var", type=float, default=1e-4)
    parser.add_argument("--min-iqr", type=float, default=0.01)
    parser.add_argument("--min-nonzero-frac", type=float, default=0.2)
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t", compression="gzip")
    if df.columns[0] == "chr":
        df.rename(columns={df.columns[0]: "#chr"}, inplace=True)

    meta = ["#chr", "start", "end", "phenotype_id"]
    sample_cols = [c for c in df.columns if c not in meta]

    x = df[sample_cols].apply(pd.to_numeric, errors="coerce")
    x = x.replace([np.inf, -np.inf], np.nan)

    var = x.var(axis=1, skipna=True)
    q75 = x.quantile(0.75, axis=1, numeric_only=True)
    q25 = x.quantile(0.25, axis=1, numeric_only=True)
    iqr = q75 - q25
    nonzero_frac = (x > 0).sum(axis=1) / len(sample_cols)

    keep = var.ge(args.min_var) & iqr.ge(args.min_iqr) & nonzero_frac.ge(args.min_nonzero_frac)

    out = df.loc[keep].copy()

    with gzip.open(args.output, "wt") as f:
        out.to_csv(f, sep="\t", index=False)

    print(f"[input] {df.shape[0]} phenotypes")
    print(f"[kept ] {out.shape[0]} phenotypes")
    print(f"[drop ] {df.shape[0] - out.shape[0]} phenotypes")


if __name__ == "__main__":
    main()
