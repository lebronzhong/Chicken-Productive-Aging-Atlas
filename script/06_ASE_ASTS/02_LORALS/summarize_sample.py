#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 06_ASE_ASTS
# Purpose: Summarize sample.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

import sys
import os
import numpy as np
import pandas as pd
from scipy.stats import binom_test, chi2_contingency

# 1
asts_file = sys.argv[1]
txgene_file = sys.argv[2]
outdir = sys.argv[3]
sample = sys.argv[4]
dataset = sys.argv[5]
tissue = sys.argv[6]
subtissue = sys.argv[7]
os.makedirs(outdir, exist_ok=True)
# 2
df = pd.read_csv(asts_file, sep="\t")
txgene = pd.read_csv(txgene_file, sep="\t")
df["transcript"] = df["transcript"].astype(str)
txgene["transcript"] = txgene["transcript"].astype(str)
df = df.merge(txgene, on="transcript", how="inner")
for c in ["position", "refCount", "altCount"]:
    df[c] = pd.to_numeric(df[c], errors="coerce")


# 3
def bh(p):
    p = np.asarray(p, dtype=float)
    if len(p) == 0:
        return np.array([])
    order = np.argsort(p)
    ranked = p[order]
    n = len(p)
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.minimum(q, 1.0)
    result = np.empty(n)
    result[order] = q
    return result


# 4
ase_columns = [
    "sample",
    "dataset",
    "tissue",
    "subtissue",
    "gene_id",
    "contig",
    "position",
    "refAllele",
    "altAllele",
    "refCount",
    "altCount",
    "total_coverage",
    "log2_ALT_REF",
    "p_value",
]
asts_columns = [
    "sample",
    "dataset",
    "tissue",
    "subtissue",
    "gene_id",
    "contig",
    "position",
    "refAllele",
    "altAllele",
    "total_coverage",
    "number_of_isoforms",
    "cohen_effect",
    "chi2",
    "dof",
    "p_value",
]
ase_rows = []
asts_rows = []
# 5
group_cols = ["gene", "contig", "position", "refAllele", "altAllele"]
for key, sub in df.groupby(group_cols, sort=False):
    gene = key[0]
    chrom = key[1]
    pos = key[2]
    ref = key[3]
    alt = key[4]
    ref_total = int(sub["refCount"].sum())
    alt_total = int(sub["altCount"].sum())
    total = ref_total + alt_total
    # 6
    if total >= 10:
        p_ase = binom_test(ref_total, n=total, p=0.5)
        log2fc = np.log2((alt_total + 0.5) / (ref_total + 0.5))
        ase_rows.append(
            {
                "sample": sample,
                "dataset": dataset,
                "tissue": tissue,
                "subtissue": subtissue,
                "gene_id": gene,
                "contig": chrom,
                "position": int(pos),
                "refAllele": ref,
                "altAllele": alt,
                "refCount": ref_total,
                "altCount": alt_total,
                "total_coverage": total,
                "log2_ALT_REF": log2fc,
                "p_value": p_ase,
            }
        )
    # 7
    tx = sub[(sub["refCount"] >= 10) | (sub["altCount"] >= 10)].copy()
    if len(tx) < 2:
        continue
    if tx["refCount"].sum() == 0 or tx["altCount"].sum() == 0:
        continue
    table = np.vstack([tx["refCount"].values, tx["altCount"].values])
    try:
        chi2, p_asts, dof, expected = chi2_contingency(table)
    except ValueError:
        continue
    total_asts = int(table.sum())
    effect = np.sqrt(chi2 / total_asts)
    asts_rows.append(
        {
            "sample": sample,
            "dataset": dataset,
            "tissue": tissue,
            "subtissue": subtissue,
            "gene_id": gene,
            "contig": chrom,
            "position": int(pos),
            "refAllele": ref,
            "altAllele": alt,
            "total_coverage": total_asts,
            "number_of_isoforms": len(tx),
            "cohen_effect": effect,
            "chi2": chi2,
            "dof": dof,
            "p_value": p_asts,
        }
    )
# 8
ase = pd.DataFrame(ase_rows, columns=ase_columns)
asts = pd.DataFrame(asts_rows, columns=asts_columns)
# 9
if len(ase) > 0:
    ase = (
        ase.sort_values(["gene_id", "total_coverage"], ascending=[True, False])
        .drop_duplicates("gene_id", keep="first")
        .reset_index(drop=True)
    )
    ase["FDR_BH"] = bh(ase["p_value"].values)
else:
    ase["FDR_BH"] = []
# 10
if len(asts) > 0:
    asts = (
        asts.sort_values(["gene_id", "total_coverage"], ascending=[True, False])
        .drop_duplicates("gene_id", keep="first")
        .reset_index(drop=True)
    )
    asts["FDR_BH"] = bh(asts["p_value"].values)
else:
    asts["FDR_BH"] = []
# 11
ase.to_csv(os.path.join(outdir, "sample_gene_ase.tsv"), sep="\t", index=False)
asts.to_csv(os.path.join(outdir, "sample_gene_asts.tsv"), sep="\t", index=False)
# 12
ase[ase["FDR_BH"] < 0.05].to_csv(
    os.path.join(outdir, "sample_gene_ase_FDR005.tsv"), sep="\t", index=False
)
asts[asts["FDR_BH"] < 0.05].to_csv(
    os.path.join(outdir, "sample_gene_asts_FDR005.tsv"), sep="\t", index=False
)
# 13
print("Gene ASE:", len(ase))
print("Gene ASTS:", len(asts))
print("ASTS FDR < 0.05:", int((asts["FDR_BH"] < 0.05).sum()))
