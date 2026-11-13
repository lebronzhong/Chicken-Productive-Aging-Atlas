#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 07_Molecular_QTL
# Purpose: Prepare transcript qtl bed.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

import argparse
import gzip
import sys
import numpy as np
import pandas as pd
from scipy.stats import rankdata, norm


def read_table_auto(path):
    if path.endswith(".gz"):
        return pd.read_csv(path, sep="\t", compression="gzip")
    return pd.read_csv(path, sep="\t")


def rank_inverse_normal(x):
    x = np.asarray(x, dtype=float)
    mask = np.isfinite(x)
    y = np.full_like(x, np.nan, dtype=float)
    if mask.sum() <= 2:
        return y
    vals = x[mask]
    if np.nanstd(vals) == 0:
        return y
    ranks = rankdata(vals, method="average")
    transformed = norm.ppf((ranks - 0.5) / len(vals))
    y[mask] = transformed
    return y


def chrom_key(x):
    s = str(x)
    try:
        return (0, int(s))
    except ValueError:
        return (1, s)


def main():
    parser = argparse.ArgumentParser(
        description="Prepare clean normalized transcript-level TensorQTL BED."
    )
    parser.add_argument("--input-bed", required=True)
    parser.add_argument("--bim", required=True)
    parser.add_argument("--output-bed-gz", required=True)
    parser.add_argument("--tpm-threshold", type=float, default=0.1)
    parser.add_argument("--min-sample-frac", type=float, default=0.2)
    parser.add_argument("--keep-raw-chroms", action="store_true")
    parser.add_argument("--no-int", action="store_true")
    args = parser.parse_args()
    df = read_table_auto(args.input_bed)
    cols = list(df.columns)
    if cols[0] in ["chr", "#chr"]:
        cols[0] = "#chr"
    else:
        raise ValueError(f"First column should be chr or #chr, but got: {cols[0]}")
    df.columns = cols
    required = ["#chr", "start", "end", "phenotype_id"]
    for c in required:
        if c not in df.columns:
            raise ValueError(f"Missing required column: {c}")
    sample_cols = [c for c in df.columns if c not in required]
    if len(sample_cols) == 0:
        raise ValueError("No sample columns found.")
    df["#chr"] = df["#chr"].astype(str)
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df = df.dropna(subset=["start", "end"])
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df = df[(df["start"] >= 0) & (df["end"] > df["start"])].copy()
    if not args.keep_raw_chroms:
        bim = pd.read_csv(args.bim, sep=r"\s+", header=None, usecols=[0], dtype=str)
        keep_chroms = set(bim[0].astype(str).unique())
        before = df.shape[0]
        df = df[df["#chr"].isin(keep_chroms)].copy()
        after = df.shape[0]
        print(
            f"[chrom filter] kept {after}/{before} phenotypes on genotype-covered chromosomes",
            file=sys.stderr,
        )
    expr = df[sample_cols].apply(pd.to_numeric, errors="coerce")
    expr = expr.replace([np.inf, -np.inf], np.nan)
    min_n = int(np.ceil(len(sample_cols) * args.min_sample_frac))
    detected = (expr > args.tpm_threshold).sum(axis=1) >= min_n
    before = df.shape[0]
    df = df.loc[detected].copy()
    expr = expr.loc[detected].copy()
    print(
        f"[expression filter] kept {df.shape[0]}/{before} phenotypes; threshold={args.tpm_threshold}, min_samples={min_n}",
        file=sys.stderr,
    )
    expr = np.log2(expr + 1.0)
    var = expr.var(axis=1, skipna=True)
    keep_var = var > 0
    before = df.shape[0]
    df = df.loc[keep_var].copy()
    expr = expr.loc[keep_var].copy()
    print(
        f"[variance filter] kept {df.shape[0]}/{before} non-zero variance phenotypes",
        file=sys.stderr,
    )
    if not args.no_int:
        transformed = np.vstack([rank_inverse_normal(row) for row in expr.to_numpy()])
        expr = pd.DataFrame(transformed, index=expr.index, columns=sample_cols)
    row_means = expr.mean(axis=1, skipna=True)
    expr = expr.T.fillna(row_means).T.fillna(0)
    out = pd.concat([df[required].reset_index(drop=True), expr.reset_index(drop=True)], axis=1)
    out["_chr_sort"] = out["#chr"].map(chrom_key)
    out = out.sort_values(by=["_chr_sort", "start", "end", "phenotype_id"]).drop(
        columns=["_chr_sort"]
    )
    with gzip.open(args.output_bed_gz, "wt") as f:
        out.to_csv(f, sep="\t", index=False)
    print(
        f"[done] wrote {args.output_bed_gz}; phenotypes={out.shape[0]}, samples={len(sample_cols)}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
