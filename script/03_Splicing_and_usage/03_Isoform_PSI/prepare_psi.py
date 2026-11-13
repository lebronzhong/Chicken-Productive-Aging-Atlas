#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 03_Splicing_and_usage
# Purpose: Prepare psi.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

"""Preserve original SUPPA isoform PSI and cluster membership; no re-clustering."""

from pathlib import Path
import pandas as pd, numpy as np, json, re, hashlib, warnings

warnings.filterwarnings("ignore", category=RuntimeWarning)
root = Path(__file__).resolve().parents[1]
base = root.parent
old = base / "count_DTU"
out = root / "inputs"
out.mkdir(exist_ok=True)
meta = pd.read_csv(old / "inputs/sample_manifest.tsv", sep="\t")
meta["psi_sample"] = (
    meta["tissue"] + meta["age"].astype(str) + "-" + meta["individual"].str.rsplit("_", n=1).str[-1]
)
assert not meta.duplicated(["tissue", "psi_sample"]).any()
meta.to_csv(out / "sample_manifest.tsv", sep="\t", index=False)
mp = pd.read_csv(
    old / "inputs/canonical_transcript_mapping.tsv.gz", sep="\t", low_memory=False
).set_index("transcript_id")
clpath = (
    base
    / "13.5.PSI_Cluster_ALL/2.PSI_cluster/Cluster_Results_PSI/Clustered_transcripts_with_labels_update.txt"
)
cl = pd.read_csv(
    clpath,
    sep=r"\s+",
    skiprows=1,
    names=["item", "legacy_20", "legacy_50", "legacy_80", "legacy_100", "cluster"],
)
cl["tissue"] = cl.item.str.split("_ENSGALG", n=1).str[0]
# Novel gene prefixes are parsed using the exact longest tissue prefix.
ts = sorted(meta.tissue.unique(), key=len, reverse=True)


def split_item(x):
    for t in ts:
        if x.startswith(t + "_"):
            return t, x[len(t) + 1 :]
    raise ValueError(x)


parsed = cl.item.map(split_item)
cl["tissue"] = [x[0] for x in parsed]
cl["psi_id"] = [x[1] for x in parsed]
assert not cl.duplicated(["tissue", "psi_id"]).any()
cl.to_csv(out / "original_all_cluster_members.tsv.gz", sep="\t", index=False)
inventory = []
audits = []
case = []
legacy = []
dtes = []
for t in sorted(meta.tissue.unique()):
    ss = meta[meta.tissue.eq(t)].copy()
    frames = []
    for age in [20, 50, 80, 100]:
        p = (
            base
            / f"13.SUPPA_PSI/{t}/filtered_{t}_{age}wk.transcript_grouped_tpm_newtitle_isoform.psi"
        )
        d = pd.read_csv(p, sep="\t", index_col=0)
        expected = ss.loc[ss.age.eq(age), "psi_sample"].tolist()
        assert set(d.columns) == set(expected), (t, age, list(d.columns), expected)
        frames.append(d[expected])
        inventory.append(
            dict(
                path=str(p),
                bytes=p.stat().st_size,
                sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
            )
        )
    x = pd.concat(frames, axis=1)
    assert not x.index.duplicated().any()
    ids = pd.DataFrame({"psi_id": x.index})
    ids["gene_id"] = ids.psi_id.str.split(";").str[0]
    ids["source_id"] = ids.psi_id.str.split(";").str[-1]
    ids["transcript_id"] = ids.source_id.str.replace(r"\.p\d+$", "", regex=True)
    assert not ids.transcript_id.duplicated().any(), t
    assert set(ids.transcript_id) <= set(mp.index)
    ids["current_gene_id"] = mp.loc[ids.transcript_id, "gene_id"].to_numpy()
    ids["gene_mapping_consistent"] = ids.gene_id.eq(ids.current_gene_id)
    conflicted_genes = set(ids.loc[~ids.gene_mapping_consistent, "gene_id"])
    arr = x.to_numpy(float)
    assert np.all((arr[np.isfinite(arr)] >= 0) & (arr[np.isfinite(arr)] <= 1 + 1e-7))
    c = pd.read_csv(old / f"inputs/validated_fig3/{t}/tx_counts.tsv.gz", sep="\t", index_col=0)
    g = c.groupby(mp.loc[c.index, "gene_id"]).sum()
    supported = np.ones(len(g), dtype=bool)
    for age in [20, 50, 80, 100]:
        supported &= (g[ss.loc[ss.age.eq(age), "sample"]].ge(20).sum(axis=1) >= 3).to_numpy()
    gene_good = set(g.index[supported]) - conflicted_genes
    cv = c.reindex(ids.transcript_id)[ss["sample"]].to_numpy()
    tx_good = np.zeros(len(ids), dtype=bool)
    finite_good = np.ones(len(ids), dtype=bool)
    for age in [20, 50, 80, 100]:
        ix = np.flatnonzero(ss.age.to_numpy() == age)
        z = arr[:, ix]
        ids[f"n_{age}"] = np.isfinite(z).sum(axis=1)
        ids[f"mean_{age}"] = np.nanmean(z, axis=1)
        ids[f"sd_{age}"] = np.nanstd(z, axis=1, ddof=1)
        finite_good &= ids[f"n_{age}"].ge(3).to_numpy()
        tx_good |= ((cv[:, ix] >= 5) & (z >= 0.05)).sum(axis=1) >= 3
    ids["finite_3_each_age"] = finite_good
    ids["supported_isoform"] = tx_good
    ids["supported_gene"] = ids.gene_id.isin(gene_good)
    counts = ids[finite_good & tx_good].groupby("gene_id").size()
    ids["multi_isoform_gene"] = ids.gene_id.map(counts).fillna(0).ge(2)
    ids["variable"] = np.nanstd(arr, axis=1) > 1e-8
    ids["testable"] = (
        ids.finite_3_each_age
        & ids.supported_isoform
        & ids.supported_gene
        & ids.multi_isoform_gene
        & ids.variable
        & ids.gene_mapping_consistent
    )
    ids["usage_range"] = ids[[f"mean_{a}" for a in [20, 50, 80, 100]]].max(axis=1) - ids[
        [f"mean_{a}" for a in [20, 50, 80, 100]]
    ].min(axis=1)
    cc = cl[cl.tissue.eq(t)].drop(columns=["tissue", "item"])
    ids = ids.merge(cc, on="psi_id", how="left", validate="one_to_one")
    assert ids.cluster.notna().all()
    for a in [20, 50, 80, 100]:
        ids[f"imputation_difference_{a}"] = ids[f"mean_{a}"] - ids[f"legacy_{a}"]
    folder = out / t
    folder.mkdir(exist_ok=True)
    ids.to_csv(folder / "features.tsv.gz", sep="\t", index=False)
    x.index.name = "psi_id"
    x.columns = ss["sample"].tolist()
    x.to_csv(folder / "original_PSI.tsv.gz", sep="\t")
    ss.to_csv(folder / "samples.tsv", sep="\t", index=False)
    # Preserve original gene-wise corrected SUPPA P columns verbatim.
    for suffix, cons in [("", ["50vs20", "80vs50", "100vs80"]), ("_20vs100", ["100vs20"])]:
        p = base / f"13.SUPPA_PSI/{t}/filtered_{t}_diffSplice{suffix}.dpsi"
        d = pd.read_csv(p, sep="\t", index_col=0)
        for i, con in enumerate(cons):
            legacy.append(
                pd.DataFrame(
                    {
                        "tissue": t,
                        "psi_id": d.index,
                        "contrast": con,
                        "legacy_delta": d.iloc[:, 2 * i].values,
                        "legacy_gene_corrected_P": d.iloc[:, 2 * i + 1].values,
                    }
                )
            )
    p = base / f"15.Transcripts_DESeq/DETs/{t}_20vs100.txt"
    d = pd.read_csv(p, sep="\t", index_col=0)
    d["source_id"] = d.index
    d["transcript_id"] = d.source_id.str.replace(r"\.p\d+$", "", regex=True)
    d["gene_id"] = d.transcript_id.map(mp.gene_id)
    d["tissue"] = t
    d["DTE_testable"] = np.isfinite(d.padj) & np.isfinite(d.pvalue)
    dtes.append(d.reset_index(drop=True))
    audits.append(
        dict(
            tissue=t,
            original_features=len(ids),
            testable_features=int(ids.testable.sum()),
            testable_genes=ids.loc[ids.testable, "gene_id"].nunique(),
            mapping_conflict_features=int((~ids.gene_mapping_consistent).sum()),
            mapping_conflict_parent_genes=len(conflicted_genes),
            missing_values=int((~np.isfinite(arr)).sum()),
            partly_missing_features=int(
                ((ids[[f"n_{a}" for a in [20, 50, 80, 100]]] < 5).any(axis=1)).sum()
            ),
        )
    )
    print("PREPARED", t, len(ids), int(ids.testable.sum()), flush=True)
pd.DataFrame(inventory).to_csv(out / "PSI_input_manifest.tsv", sep="\t", index=False)
pd.DataFrame(audits).to_csv(out / "filter_audit.tsv", sep="\t", index=False)
pd.concat(legacy, ignore_index=True).to_csv(
    out / "legacy_SUPPA_tests.tsv.gz", sep="\t", index=False
)
pd.concat(dtes, ignore_index=True).to_csv(
    out / "original_DTE_20vs100.tsv.gz", sep="\t", index=False
)
(out / "PREPARATION_SUCCESS").write_text("SUCCESS\n")
