#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 06_ASE_ASTS
# Purpose: Postprocess ase.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

import sys
import numpy as np
import pandas as pd

# 1
inp = sys.argv[1]
out_qc = sys.argv[2]
out_sig = sys.argv[3]
out_asts = sys.argv[4]
# 2
df = pd.read_csv(inp, sep="\t")
original_columns = list(df.columns)
numeric_columns = [
    "refCount",
    "altCount",
    "totalCount",
    "rawDepth",
    "GENOTYPE_WARNING",
    "BLACKLIST",
    "MULTI_MAPPING",
    "OTHER_ALLELE_WARNING",
    "HIGH_INDEL_WARNING",
    "BINOM_P",
]
for col in numeric_columns:
    df[col] = pd.to_numeric(df[col], errors="coerce")
# 3
gene = df["GENE_ID"].astype(str)
keep = (
    (df["totalCount"] > 20)
    & (df["GENOTYPE_WARNING"] == 0)
    & (df["BLACKLIST"] == 0)
    & (df["MULTI_MAPPING"] == 0)
    & (df["OTHER_ALLELE_WARNING"] == 0)
    & (df["HIGH_INDEL_WARNING"] == 0)
    & (~gene.isin(["", ".", "nan", "None"]))
)
qc = df.loc[keep].copy()


# 4
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


# 5
qc["ASE_FDR_BH"] = np.nan
good = qc["BINOM_P"].notna()
if good.sum() > 0:
    qc.loc[good, "ASE_FDR_BH"] = bh(qc.loc[good, "BINOM_P"].values)
qc["REF_RATIO"] = qc["refCount"] / qc["totalCount"]
qc["ALT_RATIO"] = qc["altCount"] / qc["totalCount"]
# 6
qc.to_csv(out_qc, sep="\t", index=False)
qc.loc[qc["ASE_FDR_BH"] < 0.05].to_csv(out_sig, sep="\t", index=False)
# 7
asts_input = qc[original_columns].copy()
asts_input["BINOM_P_ADJUSTED"] = qc["ASE_FDR_BH"].values
asts_input.to_csv(out_asts, sep="\t", index=False)
# 8
print("Raw variants:", len(df))
print("QC variants:", len(qc))
print("ASE FDR < 0.05:", int((qc["ASE_FDR_BH"] < 0.05).sum()))
