#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 04_FIM
# Purpose: FIM expression changes and recovery.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

import os
import argparse
import pandas as pd


def read_count_file(filepath):
    df = pd.read_csv(filepath, sep="\t")
    df = df.rename(columns={df.columns[0]: "transcript_id"})
    return df


def find_cols(columns, tissue, group):
    if group == "20":
        return [c for c in columns if c.startswith(f"{tissue}20-")]
    elif group == "100":
        return [c for c in columns if c.startswith(f"{tissue}100-")]
    elif group == "M":
        return [c for c in columns if c.startswith(f"{tissue}M-")]
    else:
        return []


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--direction_file",
        default="/path/to/chicken-aging/ont/shell/15.Transcripts_DESeq/DETs/1.Fig6B_transcript_direction_detail.tsv",
    )
    parser.add_argument(
        "--count_dir",
        default="/path/to/chicken-aging/ont/shell/12_20240919_TPM.count/count/5sample_over3_Tissue",
    )
    parser.add_argument(
        "--out_dir", default="/path/to/chicken-aging/ont/shell/21.Caloric_Restriction/Fig6-CRpart"
    )
    parser.add_argument(
        "--denom_cutoff",
        type=float,
        default=1.0,
        help="Exclude unstable ratios caused by small Expr20-minus-Expr100 differences.",
    )
    args = parser.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    direction = pd.read_csv(args.direction_file, sep="\t")
    reversal = direction[direction["category"] == "Reversal"].copy()
    all_res = []
    missing_log = []
    for tissue in sorted(reversal["tissue"].unique()):
        count_file = os.path.join(
            args.count_dir, f"{tissue}_merged_transcript_grouped_counts_filtered.tsv"
        )
        if not os.path.exists(count_file):
            print(f"[Missing count] {tissue}: {count_file}")
            missing_log.append([tissue, "missing_count_file"])
            continue
        count_df = read_count_file(count_file)
        cols_20 = find_cols(count_df.columns, tissue, "20")
        cols_100 = find_cols(count_df.columns, tissue, "100")
        cols_M = find_cols(count_df.columns, tissue, "M")
        if len(cols_20) == 0 or len(cols_100) == 0 or len(cols_M) == 0:
            print(f"[Skip] {tissue}: missing group columns")
            print(f"20wk: {cols_20}")
            print(f"100wk: {cols_100}")
            print(f"M/CR: {cols_M}")
            missing_log.append([tissue, "missing_group_columns"])
            continue
        sub = reversal[reversal["tissue"] == tissue][
            ["tissue", "transcript_id", "aging_change", "CR_change", "category"]
        ].copy()
        expr = count_df[["transcript_id"] + cols_20 + cols_100 + cols_M].copy()
        for c in cols_20 + cols_100 + cols_M:
            expr[c] = pd.to_numeric(expr[c], errors="coerce")
        expr["Expr_20wk"] = expr[cols_20].mean(axis=1)
        expr["Expr_100wk"] = expr[cols_100].mean(axis=1)
        expr["Expr_CR"] = expr[cols_M].mean(axis=1)
        expr = expr[["transcript_id", "Expr_20wk", "Expr_100wk", "Expr_CR"]]
        merged = sub.merge(expr, on="transcript_id", how="left")
        merged = merged.dropna(subset=["Expr_20wk", "Expr_100wk", "Expr_CR"]).copy()
        merged["denominator"] = merged["Expr_20wk"] - merged["Expr_100wk"]
        merged = merged[merged["denominator"].abs() > args.denom_cutoff].copy()
        merged["RecoveryScore"] = (merged["Expr_CR"] - merged["Expr_100wk"]) / (
            merged["Expr_20wk"] - merged["Expr_100wk"]
        )
        all_res.append(merged)
        print(
            f"[OK] {tissue}: "
            f"20={len(cols_20)}, 100={len(cols_100)}, M={len(cols_M)}, "
            f"RecoveryScore={len(merged)}"
        )
    if len(all_res) == 0:
        raise RuntimeError(
            "No RecoveryScore values were calculated; check file names, tissue labels and column names."
        )
    detail = pd.concat(all_res, ignore_index=True)
    summary = (
        detail.groupby("tissue")
        .agg(
            n=("RecoveryScore", "count"),
            median_recovery=("RecoveryScore", "median"),
            mean_recovery=("RecoveryScore", "mean"),
        )
        .reset_index()
    )
    detail_out = os.path.join(args.out_dir, "Fig6C_recovery_score_detail.tsv")
    summary_out = os.path.join(args.out_dir, "Fig6C_recovery_score_summary.tsv")
    missing_out = os.path.join(args.out_dir, "Fig6C_missing_log.tsv")
    detail.to_csv(detail_out, sep="\t", index=False)
    summary.to_csv(summary_out, sep="\t", index=False)
    if len(missing_log) > 0:
        pd.DataFrame(missing_log, columns=["tissue", "reason"]).to_csv(
            missing_out, sep="\t", index=False
        )
    print("\nDone.")
    print(f"Detail : {detail_out}")
    print(f"Summary: {summary_out}")


if __name__ == "__main__":
    main()
