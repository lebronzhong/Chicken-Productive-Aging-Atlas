#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 03_Splicing_and_usage
# Purpose: Summarize results.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

from pathlib import Path
import pandas as pd, numpy as np, json
from scipy.stats import false_discovery_control

root = Path(__file__).resolve().parents[1]
out = root / "results/summary"
out.mkdir(exist_ok=True)
meta = pd.read_csv(root / "inputs/sample_manifest.tsv", sep="\t")
tissues = sorted(meta.tissue.unique())
assert all((root / "results" / t / "SUCCESS").exists() for t in tissues)


def write(d, name):
    d.to_csv(out / name, sep="\t", index=False)


def read_all(name):
    arr = []
    for t in tissues:
        d = pd.read_csv(root / "results" / t / name, sep="\t")
        d.insert(0, "tissue", t)
        arr.append(d)
    return pd.concat(arr, ignore_index=True)


def adjust(d, col="PValue", target="q_global", method="bh"):
    good = np.isfinite(d[col])
    d[target] = np.nan
    d.loc[good, target] = false_discovery_control(d.loc[good, col], method=method)
    return d


scopes = {}
for scope in ["four_age", "same_batch_50_100"]:
    du = adjust(read_all(f"{scope}_DTU_omnibus.tsv.gz"))
    adjust(du, target="q_global_BY", method="by")
    de = adjust(read_all(f"{scope}_DGE_omnibus.tsv.gz"))
    adjust(de, target="q_global_BY", method="by")
    du["DTU_q_significant"] = du.q_global.lt(0.05)
    du["DTU_effect_significant"] = du.DTU_q_significant & du.max_usage_range.ge(0.10)
    write(du, f"{scope}_DTU_all.tsv.gz")
    write(de, f"{scope}_DGE_all.tsv.gz")
    bk = read_all(f"{scope}_testability.tsv.gz")
    bk = bk.merge(
        du[
            [
                "tissue",
                "gene_id",
                "PValue",
                "q_global",
                "q_global_BY",
                "max_usage_range",
                "DTU_q_significant",
                "DTU_effect_significant",
            ]
        ],
        on=["tissue", "gene_id"],
        how="left",
    ).rename(columns={"PValue": "DTU_P", "q_global": "DTU_q", "q_global_BY": "DTU_q_BY"})
    bk = bk.merge(
        de[["tissue", "gene_id", "PValue", "q_global", "q_global_BY"]],
        on=["tissue", "gene_id"],
        how="left",
    ).rename(columns={"PValue": "DGE_P", "q_global": "DGE_q", "q_global_BY": "DGE_q_BY"})
    common = np.isfinite(bk.DGE_q) & np.isfinite(bk.DTU_q)
    bk["common_testable"] = common
    ds = bk.DGE_q.lt(0.05)
    us = bk.DTU_q.lt(0.05) & bk.max_usage_range.ge(0.10)
    bk["category"] = np.select(
        [~common, ds & us, ~ds & us, ds & ~us],
        ["Not jointly testable", "Both", "DTU only", "DGE only"],
        default="Neither",
    )
    bk["strict_DTU_without_significant_DGE"] = common & us & bk.DGE_q.ge(0.05)
    write(bk, f"{scope}_DGE_DTU_correspondence.tsv.gz")
    counts = (
        bk.groupby(["tissue", "category"])
        .size()
        .unstack(fill_value=0)
        .reindex(
            columns=["Both", "DTU only", "DGE only", "Neither", "Not jointly testable"],
            fill_value=0,
        )
    )
    counts["jointly_testable"] = counts[["Both", "DTU only", "DGE only", "Neither"]].sum(axis=1)
    for cat in ["Both", "DTU only", "DGE only", "Neither"]:
        counts[cat + "_fraction"] = counts[cat] / counts.jointly_testable
    write(counts.reset_index(), f"{scope}_category_counts.tsv")
    dc = adjust(read_all(f"{scope}_DTU_contrasts.tsv.gz"))
    ec = adjust(read_all(f"{scope}_DGE_contrasts.tsv.gz"))
    tc = read_all(f"{scope}_isoform_contrasts.tsv.gz")
    tc["real_isoform"] = ~tc.transcript_id.str.startswith("OTHER_")
    tc["q_global"] = np.nan
    good = tc.real_isoform & np.isfinite(tc.PValue)
    tc.loc[good, "q_global"] = false_discovery_control(tc.loc[good, "PValue"])
    write(dc, f"{scope}_DTU_contrasts_all.tsv.gz")
    write(ec, f"{scope}_DGE_contrasts_all.tsv.gz")
    write(tc, f"{scope}_isoform_contrasts_all.tsv.gz")
    cc = dc.merge(
        ec[["tissue", "gene_id", "contrast", "PValue", "q_global", "logFC"]],
        on=["tissue", "gene_id", "contrast"],
        how="outer",
        suffixes=("_DTU", "_DGE"),
    )
    cc["common_testable"] = np.isfinite(cc.q_global_DTU) & np.isfinite(cc.q_global_DGE)
    ds = cc.q_global_DGE.lt(0.05)
    us = cc.q_global_DTU.lt(0.05) & cc.max_abs_delta.ge(0.10)
    cc["category"] = np.select(
        [~cc.common_testable, ds & us, ~ds & us, ds & ~us],
        ["Not jointly testable", "Both", "DTU only", "DGE only"],
        default="Neither",
    )
    write(cc, f"{scope}_contrast_DGE_DTU_correspondence.tsv.gz")
    scopes[scope] = dict(du=du, de=de, bk=bk, dc=dc, ec=ec, tc=tc)
    print(scope, "DTU", us.sum(), "common", common.sum(), flush=True)
du = scopes["four_age"]["du"]
tc = scopes["four_age"]["tc"]
dc = scopes["four_age"]["dc"]
mp = pd.read_csv(
    root / "inputs/canonical_transcript_mapping.tsv.gz", sep="\t", low_memory=False
).set_index("transcript_id")
ann = pd.read_csv(
    root / "inputs/canonical_annotations.tsv.gz", sep="\t", low_memory=False
).set_index("transcript_id")
quality = pd.read_csv(root / "inputs/quality.tsv.gz", sep="\t", low_memory=False).set_index(
    "transcript_id"
)
# One pair per gene-tissue-contrast, chosen by the two opposite maximum fraction changes.
hits = dc[dc.q_global.lt(0.05) & dc.max_abs_delta.ge(0.10)][
    ["tissue", "gene_id", "contrast", "q_global", "max_abs_delta"]
]
eligible = tc[tc.real_isoform & tc.q_global.lt(0.05) & tc.delta_fraction.abs().ge(0.05)].merge(
    hits, on=["tissue", "gene_id", "contrast"], suffixes=("_isoform", "_gene")
)
pairrows = []
for key, d in eligible.groupby(["tissue", "gene_id", "contrast"]):
    pos = d[d.delta_fraction.gt(0)].sort_values(
        ["delta_fraction", "transcript_id"], ascending=[False, True]
    )
    neg = d[d.delta_fraction.lt(0)].sort_values(["delta_fraction", "transcript_id"])
    if pos.empty or neg.empty:
        continue
    up = pos.iloc[0]
    dn = neg.iloc[0]
    u = mp.loc[up.transcript_id]
    v = mp.loc[dn.transcript_id]
    ue = [tuple(map(int, x.split("-"))) for x in u.exons.split(";")]
    ve = [tuple(map(int, x.split("-"))) for x in v.exons.split(";")]
    uj = {(ue[i][1], ue[i + 1][0]) for i in range(len(ue) - 1)}
    vj = {(ve[i][1], ve[i + 1][0]) for i in range(len(ve) - 1)}
    plus = u.strand == "+"
    uf, ul = (ue[0], ue[-1]) if plus else (ue[-1], ue[0])
    vf, vl = (ve[0], ve[-1]) if plus else (ve[-1], ve[0])
    rec = dict(
        tissue=key[0],
        gene_id=key[1],
        contrast=key[2],
        up_isoform=up.transcript_id,
        down_isoform=dn.transcript_id,
        up_delta=up.delta_fraction,
        down_delta=dn.delta_fraction,
        gene_q=up.q_global_gene,
        up_q=up.q_global_isoform,
        down_q=dn.q_global_isoform,
        first_exon_changed=uf != vf,
        last_exon_changed=ul != vl,
        junction_chain_changed=uj != vj,
        shared_junctions=len(uj & vj),
        up_unique_junctions=len(uj - vj),
        down_unique_junctions=len(vj - uj),
        same_terminal_exons=uf == vf and ul == vl,
    )
    au = ann.loc[up.transcript_id]
    av = ann.loc[dn.transcript_id]
    rec["symbol"] = au.symbol
    for side, a in [("up", au), ("down", av)]:
        for field in [
            "structural_category",
            "subcategory",
            "ORF_aa",
            "ORF_unambiguous",
            "ORF_types",
            "predicted_NMD",
            "UTR5_nt",
            "UTR3_nt",
            "domain_assessable",
            "domains",
            "all_canonical",
            "RTS_stage",
            "latest_protein_agrees",
        ]:
            rec[f"{side}_{field}"] = a[field]
    rec["ORF_assessable"] = bool(
        au.ORF_unambiguous
        and av.ORF_unambiguous
        and np.isfinite(au.ORF_aa)
        and np.isfinite(av.ORF_aa)
    )
    rec["delta_ORF_aa"] = au.ORF_aa - av.ORF_aa if rec["ORF_assessable"] else np.nan
    rec["NMD_assessable"] = str(au.predicted_NMD).upper() in ["TRUE", "FALSE"] and str(
        av.predicted_NMD
    ).upper() in ["TRUE", "FALSE"]
    rec["NMD_gain"] = (
        rec["NMD_assessable"]
        and str(au.predicted_NMD).upper() == "TRUE"
        and str(av.predicted_NMD).upper() == "FALSE"
    )
    rec["NMD_loss"] = (
        rec["NMD_assessable"]
        and str(au.predicted_NMD).upper() == "FALSE"
        and str(av.predicted_NMD).upper() == "TRUE"
    )
    rec["domain_pair_assessable"] = bool(au.domain_assessable and av.domain_assessable)
    rec["domain_changed"] = rec["domain_pair_assessable"] and str(au.domains) != str(av.domains)
    rec["structure_strict"] = bool(
        au.all_canonical == "canonical"
        and av.all_canonical == "canonical"
        and str(au.RTS_stage).upper() == "FALSE"
        and str(av.RTS_stage).upper() == "FALSE"
        and "incomplete" not in str(au.structural_category)
        and "incomplete" not in str(av.structural_category)
    )
    rec["structure_strict"] = bool(
        rec["structure_strict"]
        and au.structural_category
        in ["full-splice_match", "novel_in_catalog", "novel_not_in_catalog"]
        and av.structural_category
        in ["full-splice_match", "novel_in_catalog", "novel_not_in_catalog"]
    )
    rec["ORF_assessable"] = bool(
        rec["ORF_assessable"]
        and rec["structure_strict"]
        and au.latest_protein_agrees
        and av.latest_protein_agrees
        and au.ORF_types == "complete"
        and av.ORF_types == "complete"
    )
    if not rec["ORF_assessable"]:
        rec["delta_ORF_aa"] = np.nan
    rec["NMD_assessable"] = bool(rec["NMD_assessable"] and rec["ORF_assessable"])
    rec["domain_pair_assessable"] = bool(rec["domain_pair_assessable"] and rec["ORF_assessable"])
    rec["exact_read_chain_evaluated"] = all(
        quality.loc[s, "lr_exact_chain_evaluation"] != "not_evaluated"
        for s in [up.transcript_id, dn.transcript_id]
    )
    pairrows.append(rec)
pairs = pd.DataFrame(pairrows)
write(pairs, "functional_pairs_all_contrasts.tsv.gz")
functional = []
for con, d in pairs.groupby("contrast"):
    for label, assess, positive in [
        ("Longer ORF favored", d.ORF_assessable, d.delta_ORF_aa.gt(0)),
        ("Shorter ORF favored", d.ORF_assessable, d.delta_ORF_aa.lt(0)),
        ("Predicted NMD gain", d.NMD_assessable, d.NMD_gain),
        ("Predicted NMD loss", d.NMD_assessable, d.NMD_loss),
        ("Domain difference", d.domain_pair_assessable, d.domain_changed),
        ("First exon difference", pd.Series(True, index=d.index), d.first_exon_changed),
        ("Last exon difference", pd.Series(True, index=d.index), d.last_exon_changed),
        ("Junction-chain difference", pd.Series(True, index=d.index), d.junction_chain_changed),
    ]:
        functional.append(
            dict(
                contrast=con,
                consequence=label,
                all_pairs=len(d),
                assessable=int(assess.sum()),
                changed=int((assess & positive).sum()),
                fraction_changed=(
                    float((assess & positive).sum() / assess.sum()) if assess.sum() else np.nan
                ),
                unique_genes=d.gene_id.nunique(),
            )
        )
write(pd.DataFrame(functional), "functional_consequence_summary.tsv")
# New dynamics retain individual-derived fractions; no zero imputation.
dynamic = []
significant = set(
    zip(du.loc[du.DTU_effect_significant, "tissue"], du.loc[du.DTU_effect_significant, "gene_id"])
)
localised = tc[tc.real_isoform & tc.q_global.lt(0.05) & tc.delta_fraction.abs().ge(0.10)]
localised_ids = set(zip(localised.tissue, localised.transcript_id))
for t in tissues:
    eff = pd.read_csv(root / "results" / t / "four_age_isoform_effects.tsv.gz", sep="\t")
    eff = eff[
        eff.eligible_isoform
        & eff.usage_range.ge(0.10)
        & eff.gene_id.map(lambda g: (t, g) in significant)
        & eff.transcript_id.map(lambda tr: (t, tr) in localised_ids)
    ].copy()
    eff.insert(0, "tissue", t)

    def classify(row):
        v = row[["mean_20", "mean_50", "mean_80", "mean_100"]].to_numpy(dtype=float)
        delta = np.diff(v)
        if (delta >= -0.05).all() and v[-1] - v[0] >= 0.10:
            return "Increasing"
        if (delta <= 0.05).all() and v[0] - v[-1] >= 0.10:
            return "Decreasing"
        if v.argmax() == 1 and v[1] - max(v[0], v[3]) >= 0.10:
            return "50-week peak"
        if v.argmax() == 2 and v[2] - max(v[0], v[3]) >= 0.10:
            return "80-week peak"
        if v.argmin() == 1 and min(v[0], v[3]) - v[1] >= 0.10:
            return "50-week trough"
        if v.argmin() == 2 and min(v[0], v[3]) - v[2] >= 0.10:
            return "80-week trough"
        return "Complex"

    eff["pattern"] = eff.apply(classify, axis=1)
    dynamic.append(eff)
dynamic = pd.concat(dynamic, ignore_index=True)
write(dynamic, "usage_dynamic_members.tsv.gz")
long = dynamic.melt(
    id_vars=["tissue", "gene_id", "transcript_id", "pattern"],
    value_vars=["mean_20", "mean_50", "mean_80", "mean_100"],
    var_name="age",
    value_name="fraction",
)
long.age = long.age.str.replace("mean_", "").astype(int)
profile = (
    long.groupby(["pattern", "age"])
    .fraction.agg(
        n="size", median="median", q25=lambda x: x.quantile(0.25), q75=lambda x: x.quantile(0.75)
    )
    .reset_index()
)
write(profile, "usage_dynamic_profiles.tsv")
write(
    dynamic.groupby(["tissue", "pattern"])
    .agg(records=("transcript_id", "size"), unique_genes=("gene_id", "nunique"))
    .reset_index(),
    "usage_dynamic_tissue_counts.tsv",
)
write(
    dynamic.groupby("pattern")
    .agg(
        transcript_tissue_records=("transcript_id", "size"),
        unique_transcripts=("transcript_id", "nunique"),
        unique_genes=("gene_id", "nunique"),
    )
    .reset_index(),
    "usage_dynamic_units.tsv",
)
# Exact contrast sensitivity and complete mapping to Fig3's final supplied-on-disk membership.
sens = []
for scope, objs in scopes.items():
    bk = objs["bk"]
    common = bk.common_testable
    for threshold in [0.05, 0.10, 0.20]:
        flag = bk.DTU_q.lt(0.05) & bk.max_usage_range.ge(threshold)
        sens.append(
            dict(
                scope=scope,
                threshold=threshold,
                common_genes=int(common.sum()),
                DTU=int((flag & common).sum()),
                DTU_no_significant_DGE=int((flag & common & bk.DGE_q.ge(0.05)).sum()),
                DTU_BY=int(
                    (common & bk.DTU_q_BY.lt(0.05) & bk.max_usage_range.ge(threshold)).sum()
                ),
            )
        )
write(pd.DataFrame(sens), "threshold_and_dependence_sensitivity.tsv")
bridge = root / "inputs/fig3_bridge"
if (bridge / "gene_tissue_tests.tsv.gz").exists():
    old = pd.read_csv(bridge / "gene_tissue_tests.tsv.gz", sep="\t", low_memory=False)
    merged = old.merge(
        scopes["four_age"]["bk"][
            [
                "tissue",
                "gene_id",
                "DTU_q",
                "max_usage_range",
                "DGE_q",
                "category",
                "common_testable",
            ]
        ],
        on=["tissue", "gene_id"],
        how="left",
    )
    write(merged, "Fig3_expression_divergence_DTU_bridge.tsv.gz")
    counts = (
        merged.groupby(["divergent", "category"], dropna=False)
        .size()
        .reset_index(name="gene_tissue_records")
    )
    write(counts, "Fig3_bridge_summary.tsv")
    oldmembers = pd.read_csv(
        bridge / "primary_age_associated_members.tsv.gz",
        sep="\t",
        usecols=["tissue", "transcript_id", "gene_id", "old_cluster"],
    )
    write(
        oldmembers.merge(
            dynamic[["tissue", "transcript_id", "pattern", "usage_range"]],
            on=["tissue", "transcript_id"],
            how="left",
        ),
        "Fig3_cluster_members_DTU_bridge.tsv.gz",
    )
# Candidate ranking uses known structure quality and a fixed endpoint, without claiming read validation.
candidates = pairs[
    pairs.contrast.eq("100vs20") & pairs.structure_strict & pairs.ORF_assessable
].copy()
candidates["effect"] = np.minimum(candidates.up_delta.abs(), candidates.down_delta.abs())
candidates = candidates.sort_values(
    ["same_terminal_exons", "effect", "gene_q", "tissue", "gene_id"],
    ascending=[False, False, True, True, True],
)
write(candidates, "case_candidate_ranking.tsv")
print("SUMMARY COMPLETE", len(pairs), "pairs", len(dynamic), "dynamic records", flush=True)
