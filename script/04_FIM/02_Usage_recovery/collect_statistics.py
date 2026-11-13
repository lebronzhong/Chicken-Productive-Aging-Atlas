#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 04_FIM
# Purpose: Collect statistics.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

"""Apply prespecified FDR families; reuse, never refit, original natural-age tests."""

from pathlib import Path
import pandas as pd, numpy as np, json

r = Path(__file__).resolve().parents[1]
s = r / "results/summary"
s.mkdir(exist_ok=True)
tissues = sorted(pd.read_csv(r / "inputs/sample_manifest.tsv", sep="\t").tissue.unique())
families = []


def bh(x):
    a = np.asarray(x, float)
    q = np.full(len(a), np.nan)
    ix = np.flatnonzero(np.isfinite(a))
    o = ix[np.argsort(a[ix], kind="stable")]
    q[o] = np.minimum.accumulate((a[o] * len(o) / np.arange(1, len(o) + 1))[::-1])[::-1].clip(0, 1)
    return q


def save_family(d, name, p="PValue", q="q_global"):
    d[q] = bh(d[p])
    families.append(
        dict(
            family=name,
            n_tests=int(d[p].notna().sum()),
            BH_scope="all 33 tissues; no age-hit prescreening",
        )
    )
    d.to_csv(s / (name + ".tsv.gz"), sep="\t", index=False)
    return d


for layer in ["original_isoform", "event", "count_DTU"]:
    fname = layer + "_FIM_" + ("feature_tests" if layer == "count_DTU" else "tests") + ".tsv.gz"
    all = []
    for t in tissues:
        d = pd.read_csv(r / "results" / t / fname, sep="\t")
        d["tissue"] = t
        if layer == "count_DTU":
            d = d[~d.residual_category].copy()
            d["feature_id"] = d.transcript_id
        all.append(d)
    d = save_family(pd.concat(all, ignore_index=True), layer + "_FIM_all_tests")
    for t in tissues:
        d[d.tissue.eq(t)].to_csv(
            r / f"results/{t}/{layer}_FIM_adjusted.tsv.gz", sep="\t", index=False
        )
    if layer == "original_isoform":
        d["gene_id"] = d.feature_id.str.split(";").str[0]
        g = d.groupby(["tissue", "gene_id"]).PValue.agg(["min", "count"]).reset_index()
        g["PValue"] = (g["min"] * g["count"]).clip(upper=1)
        save_family(g, "original_isoform_FIM_gene_tests")
    if layer == "count_DTU":
        gg = []
        for t in tissues:
            g = pd.read_csv(r / f"results/{t}/count_DTU_FIM_gene_tests.tsv.gz", sep="\t")
            g["tissue"] = t
            f = pd.read_csv(r / f"inputs/{t}/count_DTU/features.tsv.gz", sep="\t")
            f["delta"] = (f.mean_FIM - f.mean_100).abs()
            effects = f[f.FIM_testable].groupby("gene_id").delta.max()
            g["max_abs_delta"] = g.gene_id.map(effects)
            gg.append(g)
        g = save_family(pd.concat(gg, ignore_index=True), "count_DTU_FIM_gene_tests")
        for t in tissues:
            g[g.tissue.eq(t)].to_csv(
                r / f"results/{t}/count_DTU_FIM_gene_adjusted.tsv.gz", sep="\t", index=False
            )
# New event natural-age tests: all six contrasts are one family, omnibus separate.
ev = pd.concat(
    [
        pd.read_csv(r / f"results/{t}/event_natural_tests.tsv.gz", sep="\t").assign(tissue=t)
        for t in tissues
    ],
    ignore_index=True,
)
ev = pd.concat(
    [
        save_family(ev[ev.contrast.ne("four_age")].copy(), "event_natural_six_contrasts"),
        save_family(ev[ev.contrast.eq("four_age")].copy(), "event_natural_omnibus"),
    ],
    ignore_index=True,
)
for t in tissues:
    ev[ev.tissue.eq(t)].to_csv(
        r / f"results/{t}/event_natural_adjusted.tsv.gz", sep="\t", index=False
    )
del ev
# Chunk original PSI archive (raw-PSI sensitivity fits are deliberately excluded).
pieces = []
p = r.parent / "isoform_PSI/results/summary/all_PSI_feature_tests.tsv.gz"
for d in pd.read_csv(p, sep="\t", chunksize=300000):
    d = d[d.scale.eq("arcsine_sqrt") & d.contrast.isin(["100vs20", "100vs50", "four_age"])].copy()
    if len(d):
        pieces.append(d)
d = pd.concat(pieces, ignore_index=True).rename(columns={"psi_id": "feature_id"})
assert not d.duplicated(["tissue", "feature_id", "contrast"]).any()
for t in tissues:
    d[d.tissue.eq(t)].to_csv(
        r / f"results/{t}/original_isoform_natural_reused.tsv.gz", sep="\t", index=False
    )
families.append(
    dict(
        family="original_isoform_natural_REUSED",
        n_tests="original atlas six-contrast family / separate omnibus",
        BH_scope="unchanged q_global from Fig4, arcsine_sqrt",
    )
)
del d, pieces
b = r.parent / "count_DTU/results/summary"
for fn, tag, key in [
    ("four_age_isoform_contrasts_all.tsv.gz", "feature", "transcript_id"),
    ("four_age_DTU_contrasts_all.tsv.gz", "gene", "gene_id"),
]:
    pieces = []
    for d in pd.read_csv(b / fn, sep="\t", chunksize=300000):
        d = d[d.contrast.isin(["100vs20", "100vs50"])].copy()
        if tag == "feature":
            d = d[d.real_isoform].copy()
            d["feature_id"] = d.transcript_id
        if len(d):
            pieces.append(d)
    d = pd.concat(pieces, ignore_index=True)
    for t in tissues:
        d[d.tissue.eq(t)].to_csv(
            r / f"results/{t}/count_DTU_natural_{tag}_reused.tsv.gz", sep="\t", index=False
        )
    d.to_csv(s / f"count_DTU_natural_{tag}_endpoints_reused.tsv.gz", sep="\t", index=False)
om = pd.read_csv(b / "four_age_DTU_all.tsv.gz", sep="\t")
for t in tissues:
    om[om.tissue.eq(t)].to_csv(
        r / f"results/{t}/count_DTU_natural_omnibus_reused.tsv.gz", sep="\t", index=False
    )
families.append(
    dict(
        family="count_DTU_natural_REUSED",
        n_tests="original gene / real-feature six-contrast families",
        BH_scope="unchanged Fig4 q_global",
    )
)
pd.DataFrame(families).to_csv(s / "multiple_testing_families.tsv", sep="\t", index=False)
(r / "results/STATISTICS_READY").write_text("SUCCESS")
