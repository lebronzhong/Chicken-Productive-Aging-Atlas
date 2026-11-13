#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 04_FIM
# Purpose: Prepare tissue.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

from pathlib import Path
import sys, pandas as pd, numpy as np, warnings, json
from scipy import sparse

warnings.filterwarnings("ignore", category=RuntimeWarning)
r = Path(__file__).resolve().parents[1]
t = sys.argv[1]
old = r.parent / "isoform_PSI"
inp = r / "inputs" / t
out = r / "results" / t
out.mkdir(exist_ok=True)
mp = pd.read_csv(r / "inputs/canonical_transcript_mapping.tsv.gz", sep="\t")
sm = pd.read_csv(inp / "samples.tsv", sep="\t")
x = np.load(inp / "counts.npz")["counts"]
tpm = np.load(inp / "tpm.npz")["TPM"]
assert x.shape == tpm.shape == (len(mp), 25)
groups = ["20", "50", "80", "100", "FIM"]
grp = sm.group.astype(str).to_numpy()
ix = {g: np.flatnonzero(grp == g) for g in groups}
assert all(len(v) == 5 for v in ix.values())
geneidx, genes = pd.factorize(mp.gene_id, sort=False)
gm = sparse.csr_matrix(
    (np.ones(len(mp)), (geneidx, np.arange(len(mp)))), shape=(len(genes), len(mp))
)
gt = gm @ x
gene_tpm = gm @ tpm
gpos = {g: i for i, g in enumerate(genes)}
txpos = {v: i for i, v in enumerate(mp.transcript_id)}


def summary(a, vals):
    for g, j in ix.items():
        a[f"mean_{g}"] = np.nanmean(vals[:, j], axis=1)
        a[f"sd_{g}"] = np.nanstd(vals[:, j], axis=1, ddof=1)
        a[f"n_{g}"] = np.isfinite(vals[:, j]).sum(axis=1)
    return a


def save_layer(name, a, v):
    folder = inp / name
    folder.mkdir(exist_ok=True)
    summary(a, v).to_csv(folder / "features.tsv.gz", sep="\t", index=False)
    np.savez_compressed(folder / "individuals.npz", values=v)
    keep = a.FIM_testable.to_numpy()
    pd.DataFrame(
        v[keep][:, np.r_[ix["100"], ix["FIM"]]],
        index=a.loc[keep, "feature_id"],
        columns=sm.loc[np.r_[ix["100"], ix["FIM"]], "sample"],
    ).rename_axis("feature_id").to_csv(folder / "FIM_test_matrix.tsv.gz", sep="\t")
    if name == "event":
        keep = a.natural_testable.to_numpy()
        pd.DataFrame(
            v[keep, :20], index=a.loc[keep, "feature_id"], columns=sm["sample"].iloc[:20]
        ).rename_axis("feature_id").to_csv(folder / "natural_test_matrix.tsv.gz", sep="\t")
    print(
        "LAYER",
        t,
        name,
        len(a),
        int(a.natural_testable.sum()),
        int(a.FIM_testable.sum()),
        flush=True,
    )


# Original SUPPA isoform PSI — preserve exact original values, including missingness.
a = pd.read_csv(old / f"inputs/{t}/features.tsv.gz", sep="\t")
ps = pd.read_csv(old / f"inputs/{t}/original_PSI.tsv.gz", sep="\t", index_col=0)
fimpath = r.parent / f"13.SUPPA_PSI/{t}/filtered_{t}_M.transcript_grouped_tpm_newtitle_isoform.psi"
fm = pd.read_csv(fimpath, sep="\t", index_col=0)
assert fm.index.is_unique
assert set(a.psi_id) <= set(fm.index)
v = np.c_[
    ps.loc[a.psi_id, sm["sample"].iloc[:20]].to_numpy(),
    fm.loc[a.psi_id, sm.source_sample.iloc[20:]].to_numpy(),
]
a = a.rename(columns={"testable": "natural_testable"})
a["feature_id"] = a.psi_id
txix = np.array([txpos[z] for z in a.transcript_id])
gi = np.array([gpos.get(z, -1) for z in a.gene_id])
gc = gt[np.maximum(gi, 0)]
gc[gi < 0] = 0
cv = x[txix]
finite = np.ones(len(a), bool)
genegood = np.ones(len(a), bool)
supported = np.zeros(len(a), bool)
for g in ["100", "FIM"]:
    j = ix[g]
    finite &= np.isfinite(v[:, j]).sum(axis=1) >= 3
    genegood &= (gc[:, j] >= 20).sum(axis=1) >= 3
    supported |= ((cv[:, j] >= 5) & (v[:, j] >= 0.05)).sum(axis=1) >= 3
badparents = set(a.loc[~a.gene_mapping_consistent, "gene_id"])
nr = a.loc[finite & supported].groupby("gene_id").size()
multi = a.gene_id.map(nr).fillna(0).ge(2).to_numpy()
variable = np.nanstd(v[:, np.r_[ix["100"], ix["FIM"]]], axis=1) > 1e-8
a["FIM_finite_supported"] = finite
a["FIM_gene_supported"] = genegood
a["FIM_isoform_supported"] = supported
a["FIM_multi_isoform"] = multi
a["FIM_variable"] = variable
a["FIM_testable"] = (
    finite
    & genegood
    & supported
    & multi
    & variable
    & a.gene_mapping_consistent
    & ~a.gene_id.isin(badparents)
)
a["FIM_testability_reason"] = np.select(
    [
        ~a.gene_mapping_consistent | a.gene_id.isin(badparents),
        ~finite,
        ~genegood,
        ~supported,
        ~multi,
        ~variable,
    ],
    [
        "parent_mapping_conflict",
        "insufficient_finite_PSI",
        "insufficient_gene_count",
        "insufficient_isoform_support",
        "fewer_than_two_supported_isoforms",
        "constant_PSI",
    ],
    default="tested",
)
save_layer("original_isoform", a, v)
legacy = pd.read_csv(
    r.parent / f"13.SUPPA_PSI/{t}/filtered_{t}_diffSplice_100vsM.dpsi", sep="\t", index_col=0
)
legacy.columns = ["legacy_delta_FIM", "legacy_gene_corrected_P"]
legacy.index.name = "feature_id"
legacy.to_csv(inp / "legacy_FIM_SUPPA.tsv.gz", sep="\t")
# Count fractions are a separate measure. Zero gene denominators remain missing.
den = gt[geneidx]
frac = np.divide(x, den, out=np.full(x.shape, np.nan), where=den > 0)
ag = (
    pd.read_csv(r.parent / f"count_DTU/results/{t}/four_age_isoform_effects.tsv.gz", sep="\t")
    .set_index("transcript_id")
    .loc[mp.transcript_id]
)
supported = np.zeros(len(mp), bool)
goodgene = np.ones(len(genes), bool)
for g in ["100", "FIM"]:
    j = ix[g]
    supported |= ((x[:, j] >= 5) & (frac[:, j] >= 0.05)).sum(axis=1) >= 3
    goodgene &= (gt[:, j] >= 20).sum(axis=1) >= 3
retained = np.bincount(geneidx, weights=supported, minlength=len(genes))
eligible = goodgene & (retained >= 2)
ac = mp[["transcript_id", "gene_id"]].copy()
ac["feature_id"] = ac.transcript_id
ac["natural_testable"] = ag.eligible_isoform.to_numpy() & ag.gene_DTU_eligible.to_numpy()
ac["FIM_testable"] = supported & eligible[geneidx]
ac["FIM_gene_testable"] = eligible[geneidx]
ac["FIM_isoform_supported"] = supported
ac["FIM_testability_reason"] = np.select(
    [~goodgene[geneidx], retained[geneidx] < 2, ~supported],
    ["insufficient_gene_count", "fewer_than_two_supported_isoforms", "pooled_rare_isoform"],
    default="tested",
)
save_layer("count_DTU", ac, frac)
fitix = np.flatnonzero(ac.FIM_testable)
fitx = x[fitix][:, np.r_[ix["100"], ix["FIM"]]]
ann = ac.loc[fitix, ["gene_id", "transcript_id"]].copy()
ann["residual_category"] = False
otherix = np.flatnonzero(~supported & eligible[geneidx])
remaining = gm[:, otherix] @ x[otherix][:, np.r_[ix["100"], ix["FIM"]]]
rg = np.flatnonzero(remaining.sum(axis=1) > 0)
if len(rg):
    fitx = np.r_[fitx, remaining[rg]]
    ann = pd.concat(
        [
            ann,
            pd.DataFrame(
                {
                    "gene_id": genes[rg],
                    "transcript_id": ["OTHER_" + z for z in genes[rg]],
                    "residual_category": True,
                }
            ),
        ],
        ignore_index=True,
    )
check = pd.DataFrame(fitx, index=ann.gene_id).groupby(level=0).sum()
assert np.allclose(
    check.to_numpy(), gt[[gpos[z] for z in check.index]][:, 15:25]
), "Pooled DTU denominator mismatch"
pd.DataFrame(fitx, index=ann.transcript_id, columns=sm["sample"].iloc[15:25]).rename_axis(
    "transcript_id"
).to_csv(inp / "count_DTU/FIM_count_matrix.tsv.gz", sep="\t")
ann.to_csv(inp / "count_DTU/model_annotation.tsv", sep="\t", index=False)
pd.DataFrame({"sample": sm["sample"].iloc[15:25], "library_size": x[:, 15:25].sum(axis=0)}).to_csv(
    inp / "count_DTU/library_sizes.tsv", sep="\t", index=False
)
# Event usage from canonical TPM and existing IOE structures; not read-based validation.
ed = pd.read_csv(r / "inputs/event_definitions.tsv.gz", sep="\t")
inc = sparse.load_npz(r / "inputs/event_inclusion_matrix.npz")
total = sparse.load_npz(r / "inputs/event_total_matrix.npz")
it = inc @ tpm
tt = total @ tpm
ev = np.divide(it, tt, out=np.full(it.shape, np.nan), where=tt > 0)
ic = inc @ x
tc = total @ x
ec = tc - ic
assert np.nanmin(ev) >= 0 and np.nanmax(ev) <= 1 + 1e-10
ae = ed[["event_id", "gene_id", "event_type", "scope"]].copy()
ae["feature_id"] = ae.event_id
for g, j in ix.items():
    ae["mean_denominator_count_" + g] = tc[:, j].mean(axis=1)
    ae["n_denominator_ge20_" + g] = (tc[:, j] >= 20).sum(axis=1)
for label, gs in [("natural", ["20", "50", "80", "100"]), ("FIM", ["100", "FIM"])]:
    finite = np.ones(len(ed), bool)
    covered = np.ones(len(ed), bool)
    both = [np.zeros(len(ed), bool), np.zeros(len(ed), bool)]
    for g in gs:
        j = ix[g]
        finite &= np.isfinite(ev[:, j]).sum(axis=1) >= 3
        covered &= (tc[:, j] >= 20).sum(axis=1) >= 3
        both[0] |= ((ic[:, j] >= 5) & (ev[:, j] >= 0.05)).sum(axis=1) >= 3
        both[1] |= ((ec[:, j] >= 5) & ((1 - ev[:, j]) >= 0.05)).sum(axis=1) >= 3
    variable = np.nanstd(ev[:, np.concatenate([ix[g] for g in gs])], axis=1) > 1e-8
    ae[label + "_testable"] = finite & covered & both[0] & both[1] & variable
    ae[label + "_testability_reason"] = np.select(
        [~finite, ~covered, ~both[0] | ~both[1], ~variable],
        [
            "insufficient_finite_PSI",
            "insufficient_event_denominator_counts",
            "insufficient_alternative_support",
            "constant_PSI",
        ],
        default="tested",
    )
save_layer("event", ae, ev)
(out / "PREPARED").write_text("SUCCESS")
