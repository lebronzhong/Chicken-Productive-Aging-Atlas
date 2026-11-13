#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 03_Splicing_and_usage
# Purpose: Summarize psi.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

from pathlib import Path
import pandas as pd, numpy as np, json, shutil, gc
from scipy.stats import false_discovery_control, hypergeom

r = Path(__file__).resolve().parents[1]
old = r.parent / "count_DTU"
out = r / "results/summary"
out.mkdir(exist_ok=True)


def read(p):
    return pd.read_csv(p, sep="\t", low_memory=False)


def save(d, n):
    d.to_csv(out / n, sep="\t", index=False)


def adjust(d, p="PValue", q="q_global", method="bh"):
    d[q] = np.nan
    ix = np.isfinite(d[p])
    d.loc[ix, q] = false_discovery_control(d.loc[ix, p], method=method)
    return d


meta = read(r / "inputs/sample_manifest.tsv")
ts = sorted(meta.tissue.unique())
features = []
tests = []
for t in ts:
    assert (r / f"results/{t}/SUCCESS").exists(), t
    a = read(r / f"inputs/{t}/features.tsv.gz")
    a.insert(0, "tissue", t)
    features.append(a)
    z = read(r / f"results/{t}/PSI_tests.tsv.gz")
    z.insert(0, "tissue", t)
    tests.append(z)
f = pd.concat(features, ignore_index=True)
z = pd.concat(tests, ignore_index=True)
del features, tests
gc.collect()
z = z.merge(
    f[
        ["tissue", "psi_id", "gene_id", "transcript_id", "cluster", "usage_range"]
        + [f"mean_{a}" for a in [20, 50, 80, 100]]
    ],
    on=["tissue", "psi_id"],
    validate="many_to_one",
)
z["delta_PSI"] = np.nan
for con in z.contrast.unique():
    if con != "four_age":
        a, b = map(int, con.split("vs"))
        ix = z.contrast.eq(con)
        z.loc[ix, "delta_PSI"] = z.loc[ix, f"mean_{a}"] - z.loc[ix, f"mean_{b}"]
z["family"] = np.where(z.contrast.eq("four_age"), "four_age", "six_contrasts")
z["q_global"] = np.nan
z["q_global_BY"] = np.nan
for (scale, fam), ix in z.groupby(["scale", "family"]).groups.items():
    z.loc[ix, "q_global"] = false_discovery_control(z.loc[ix, "PValue"])
    z.loc[ix, "q_global_BY"] = false_discovery_control(z.loc[ix, "PValue"], method="by")
z["effect"] = np.where(z.contrast.eq("four_age"), z.usage_range, z.delta_PSI.abs())
z["significant"] = z.q_global.lt(0.05) & z.effect.ge(0.20)
save(z, "all_PSI_feature_tests.tsv.gz")
primary = z[z.scale.eq("arcsine_sqrt")].copy()
sensitivity = z[z.scale.eq("raw_PSI")][
    ["tissue", "psi_id", "contrast", "PValue", "q_global", "significant"]
].rename(
    columns={"PValue": "raw_PSI_P", "q_global": "raw_PSI_q", "significant": "raw_PSI_significant"}
)
save(
    primary[["tissue", "psi_id", "contrast", "q_global", "significant"]].merge(
        sensitivity, on=["tissue", "psi_id", "contrast"]
    ),
    "transformation_sensitivity.tsv.gz",
)
g = (
    primary.groupby(["tissue", "gene_id", "contrast", "family"])
    .agg(
        min_P=("PValue", "min"),
        n_tested_isoforms=("psi_id", "size"),
        max_effect=("effect", "max"),
        significant_feature=("significant", "any"),
    )
    .reset_index()
)
g["PValue"] = np.minimum(1, g.min_P * g.n_tested_isoforms)
g["q_global"] = np.nan
g["q_global_BY"] = np.nan
for fam, ix in g.groupby("family").groups.items():
    g.loc[ix, "q_global"] = false_discovery_control(g.loc[ix, "PValue"])
    g.loc[ix, "q_global_BY"] = false_discovery_control(g.loc[ix, "PValue"], method="by")
g["PSI_significant"] = g.q_global.lt(0.05) & g.significant_feature
save(g, "all_PSI_gene_tests.tsv.gz")
# Omnibus membership: labels and original members are fixed; testing adds annotations.
om = primary[primary.contrast.eq("four_age")][
    ["tissue", "psi_id", "PValue", "q_global", "q_global_BY", "significant"]
]
members = f.merge(om, on=["tissue", "psi_id"], how="left", validate="one_to_one")
members["significant"] = members.significant.fillna(False)
save(members, "original_clusters_with_new_evidence.tsv.gz")
units = (
    members.groupby("cluster")
    .agg(
        original_records=("psi_id", "size"),
        unique_transcripts=("transcript_id", "nunique"),
        unique_genes=("gene_id", "nunique"),
        testable_records=("testable", "sum"),
        supported_records=("significant", "sum"),
    )
    .reset_index()
)
save(units, "cluster_units.tsv")
com = (
    members.groupby(["tissue", "cluster"])
    .agg(
        original_records=("psi_id", "size"),
        testable_records=("testable", "sum"),
        supported_records=("significant", "sum"),
        unique_genes=("gene_id", "nunique"),
    )
    .reset_index()
)
com["supported_fraction_of_testable"] = com.supported_records / com.testable_records.replace(
    0, np.nan
)
com["supported_pattern_composition"] = com.supported_records / com.groupby(
    "tissue"
).supported_records.transform("sum").replace(0, np.nan)
save(com, "cluster_tissue_composition.tsv")
prof = []
for c, d in members.groupby("cluster"):
    for subset, k in [
        ("Original", d),
        ("Testable", d[d.testable]),
        ("Supported", d[d.significant]),
    ]:
        for a in [20, 50, 80, 100]:
            x = k[f"mean_{a}"].dropna()
            legacy = k[f"legacy_{a}"].dropna()
            prof.append(
                dict(
                    cluster=c,
                    subset=subset,
                    age=a,
                    records=len(k),
                    finite_records=len(x),
                    mean=x.mean(),
                    median=x.median(),
                    q25=x.quantile(0.25),
                    q75=x.quantile(0.75),
                    legacy_mean=legacy.mean(),
                )
            )
save(pd.DataFrame(prof), "cluster_profiles.tsv")
# DTE gene-level min-P and globally corrected source-feature tests, preserving original padj too.
dt = read(r / "inputs/original_DTE_20vs100.tsv.gz")
dt = dt[dt.gene_id.notna()].copy()
dt["q_global"] = np.nan
ix = dt.DTE_testable
dt.loc[ix, "q_global"] = false_discovery_control(dt.loc[ix, "pvalue"])
dt["DTE_feature_significant"] = dt.q_global.lt(0.05) & dt.log2FoldChange.abs().ge(0.5)
dt["legacy_feature_significant"] = dt.padj.lt(0.05) & dt.log2FoldChange.abs().ge(0.5)
dg = (
    dt[dt.DTE_testable]
    .groupby(["tissue", "gene_id"])
    .agg(
        min_P=("pvalue", "min"),
        n_tested_sources=("source_id", "size"),
        significant_feature=("DTE_feature_significant", "any"),
        legacy_DTE=("legacy_feature_significant", "any"),
    )
    .reset_index()
)
dg["DTE_P"] = np.minimum(1, dg.min_P * dg.n_tested_sources)
adjust(dg, "DTE_P", "DTE_q")
dg["DTE_significant"] = dg.DTE_q.lt(0.05) & dg.significant_feature
save(dt, "DTE_source_tests_100vs20.tsv.gz")
save(dg, "DTE_gene_tests_100vs20.tsv.gz")
pg = g[g.contrast.eq("100vs20")].rename(columns={"q_global": "PSI_q", "PValue": "PSI_P"})
universe = pd.concat(
    [f[["tissue", "gene_id"]], dt[["tissue", "gene_id"]]], ignore_index=True
).drop_duplicates()
common = universe.merge(dg, on=["tissue", "gene_id"], how="left").merge(
    pg[["tissue", "gene_id", "PSI_P", "PSI_q", "PSI_significant", "max_effect"]],
    on=["tissue", "gene_id"],
    how="left",
)
common["DTE_testable"] = common.DTE_q.notna()
common["PSI_testable"] = common.PSI_q.notna()
common["jointly_testable"] = common.DTE_testable & common.PSI_testable
common["category"] = np.select(
    [
        ~common.jointly_testable,
        common.DTE_significant.eq(True) & common.PSI_significant.eq(True),
        common.DTE_significant.eq(True),
        common.PSI_significant.eq(True),
    ],
    ["Not jointly testable", "Both", "DTE only", "PSI only"],
    default="Neither",
)
save(common, "DTE_PSI_common_background_100vs20.tsv.gz")
counts = (
    common.groupby(["tissue", "category"])
    .size()
    .unstack(fill_value=0)
    .reindex(
        columns=["DTE only", "PSI only", "Both", "Neither", "Not jointly testable"], fill_value=0
    )
)
counts["jointly_testable"] = counts[["DTE only", "PSI only", "Both", "Neither"]].sum(axis=1)
counts["significant_union"] = counts[["DTE only", "PSI only", "Both"]].sum(axis=1)
save(counts.reset_index(), "DTE_PSI_category_counts.tsv")
de = read(old / "results/summary/four_age_DGE_contrasts_all.tsv.gz").rename(
    columns={"q_global": "DGE_q", "PValue": "DGE_P"}
)
de4 = read(old / "results/summary/four_age_DGE_all.tsv.gz").rename(
    columns={"q_global": "DGE_q", "PValue": "DGE_P"}
)
de4["contrast"] = "four_age"
de = pd.concat([de, de4], ignore_index=True)
de = de[["tissue", "gene_id", "contrast", "DGE_P", "DGE_q", "logFC"]]
cor = g.rename(columns={"q_global": "PSI_q", "PValue": "PSI_P"}).merge(
    de, on=["tissue", "gene_id", "contrast"], how="outer"
)
fullbg = (
    pd.concat([universe, de[["tissue", "gene_id"]]], ignore_index=True)
    .drop_duplicates()
    .merge(pd.DataFrame({"contrast": g.contrast.unique()}), how="cross")
)
cor = fullbg.merge(cor, on=["tissue", "gene_id", "contrast"], how="left")
cor["jointly_testable"] = cor.PSI_q.notna() & cor.DGE_q.notna()
cor["DGE_significant"] = cor.DGE_q.lt(0.05)
cor["category"] = np.select(
    [
        ~cor.jointly_testable,
        cor.PSI_significant.eq(True) & cor.DGE_significant,
        cor.PSI_significant.eq(True),
        cor.DGE_significant,
    ],
    ["Not jointly testable", "Both", "PSI only", "DGE only"],
    default="Neither",
)
save(cor, "PSI_DGE_all_correspondence.tsv.gz")
cs = (
    cor.groupby(["tissue", "contrast", "category"])
    .size()
    .unstack(fill_value=0)
    .reindex(
        columns=["Both", "PSI only", "DGE only", "Neither", "Not jointly testable"], fill_value=0
    )
)
cs["PSI_changed_jointly_testable"] = cs["Both"] + cs["PSI only"]
cs["PSI_without_DGE_fraction"] = cs["PSI only"] / cs.PSI_changed_jointly_testable.replace(0, np.nan)
save(cs.reset_index(), "PSI_DGE_category_counts.tsv")
# Old SUPPA is a separate comparator; its P column was already gene-wise BH corrected.
ls = read(r / "inputs/legacy_SUPPA_tests.tsv.gz")
lc = primary[primary.contrast.ne("four_age")][
    ["tissue", "psi_id", "contrast", "delta_PSI", "q_global", "significant"]
].merge(ls, on=["tissue", "psi_id", "contrast"], how="left")
lc["delta_agreement"] = np.isclose(lc.delta_PSI, lc.legacy_delta, atol=1e-6, equal_nan=False)
save(lc, "legacy_SUPPA_comparison.tsv.gz")
# Direct Fig3 bridge: retain original expression-divergence and cluster assignments.
bp = old / "inputs/fig3_bridge/gene_tissue_tests.tsv.gz"
if bp.exists():
    b = read(bp)
    save(
        b.merge(g[g.contrast.eq("four_age")], on=["tissue", "gene_id"], how="left"),
        "Fig3_expression_divergence_PSI_bridge.tsv.gz",
    )
bm = old / "inputs/fig3_bridge/primary_age_associated_members.tsv.gz"
if bm.exists():
    b = read(bm)
    save(
        b.merge(
            members[
                [
                    "tissue",
                    "transcript_id",
                    "cluster",
                    "testable",
                    "q_global",
                    "significant",
                    "usage_range",
                ]
            ],
            on=["tissue", "transcript_id"],
            how="left",
            suffixes=("_Fig3", "_PSI"),
        ),
        "Fig3_cluster_PSI_bridge.tsv.gz",
    )
sens = []
for threshold in [0.1, 0.15, 0.2]:
    for con, d in primary.groupby("contrast"):
        sens.append(
            dict(
                contrast=con,
                threshold=threshold,
                testable_features=len(d),
                BH_significant=int((d.q_global.lt(0.05) & d.effect.ge(threshold)).sum()),
                BY_significant=int((d.q_global_BY.lt(0.05) & d.effect.ge(threshold)).sum()),
            )
        )
save(pd.DataFrame(sens), "effect_and_dependence_sensitivity.tsv")
summary = {
    "original_cluster_records": len(members),
    "testable_features": int(f.testable.sum()),
    "four_age_significant_features": int(members.significant.sum()),
    "four_age_significant_gene_tissue_records": int(
        g[g.contrast.eq("four_age")].PSI_significant.sum()
    ),
    "DTE_PSI_common_gene_tissue_records": int(common.jointly_testable.sum()),
    "six_cluster_original_records": int(members.cluster.isin([2, 3, 5, 6, 7, 8]).sum()),
}
(out / "analysis_summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary), flush=True)
