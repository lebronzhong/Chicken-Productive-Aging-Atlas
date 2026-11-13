#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 08_Comparative_and_complex_traits
# Purpose: Mammalian aging signatures.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

from pathlib import Path
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--expected-reference-genes", type=int)
parser.add_argument("--expected-strict-genes", type=int)
args = parser.parse_args()
import json, itertools, hashlib
import numpy as np, pandas as pd

root = Path("/path/to/chicken-aging/analysis/conservation")
prev = root.parent / "conservation_20261008"
out = root / "results"
sp = ["Mouse", "Rat", "Macaque", "Human"]
ref = pd.read_csv(
    prev / "reference/four_species_intersection.tsv", sep="\t", dtype={"Entrez.ID": str}
).rename(columns={"Entrez.ID": "mouse_entrez_id"})
assert len(ref) > 0 and not ref.mouse_entrez_id.duplicated().any()
if args.expected_reference_genes is not None:
    assert len(ref) == args.expected_reference_genes
up = np.column_stack([(ref[f"{s}_P.Adjusted"] < 0.05) & (ref[f"{s}_Slope"] > 0) for s in sp])
dn = np.column_stack([(ref[f"{s}_P.Adjusted"] < 0.05) & (ref[f"{s}_Slope"] < 0) for s in sp])
ref["n_up"] = up.sum(axis=1)
ref["n_down"] = dn.sum(axis=1)
ref["significant_conflict"] = (ref.n_up > 0) & (ref.n_down > 0)
ref["support_n"] = np.maximum(ref.n_up, ref.n_down)
ref["direction"] = np.where(
    ref.significant_conflict,
    "conflict",
    np.where(ref.n_up > 0, "up", np.where(ref.n_down > 0, "down", "none")),
)
ref["support_species"] = [
    ";".join(s for j, s in enumerate(sp) if up[i, j] or dn[i, j]) for i in range(len(ref))
]
ref["spans_rodents_primates"] = (up[:, :2].any(axis=1) | dn[:, :2].any(axis=1)) & (
    up[:, 2:].any(axis=1) | dn[:, 2:].any(axis=1)
)
ref.to_csv(out / "mammalian_support_all7807.tsv", sep="\t", index=False)
ref[ref.significant_conflict].to_csv(
    out / "significant_direction_conflicts.tsv", sep="\t", index=False
)
mpath = (
    root.parent / "conservation_20261007/metadata/chicken_mouse_one2one_highconfidence_entrez.tsv"
)
mp = pd.read_csv(mpath, sep="\t", dtype={"mouse_entrez_id": str})
m = mp.merge(ref, on="mouse_entrez_id", validate="one_to_one")
m.to_csv(out / "mapped_tested_genes.tsv", sep="\t", index=False)
grade = []
combos = []
tissue_rows = []
pairs = []
eligible = []
for k in [2, 3, 4]:
    for mode in ["at_least", "exactly"]:
        mask = (
            ref.support_n >= k if mode == "at_least" else ref.support_n == k
        ) & ~ref.significant_conflict
        r = ref[mask]
        ids = set(r.mouse_entrez_id)
        mm = m[m.mouse_entrez_id.isin(ids)]
        grade.append(
            dict(
                mode=mode,
                n_mammals=k,
                n_mammalian=len(r),
                n_up=int((r.direction == "up").sum()),
                n_down=int((r.direction == "down").sum()),
                n_spanning_orders=int(r.spans_rodents_primates.sum()),
                n_mapped_to_chicken=len(mm),
            )
        )
        r.to_csv(out / f"mammalian_{mode}{k}.tsv", sep="\t", index=False)
        if mode == "at_least":
            mm.to_csv(out / f"mapped_at_least{k}.tsv", sep="\t", index=False)
    for combo in itertools.combinations(sp, k):
        target = ";".join(combo)
        r = ref[(ref.support_n == k) & ~ref.significant_conflict & (ref.support_species == target)]
        combos.append(
            dict(
                n_mammals=k,
                species=target,
                n_genes=len(r),
                n_up=int((r.direction == "up").sum()),
                n_down=int((r.direction == "down").sum()),
            )
        )
pd.DataFrame(grade).to_csv(out / "mammalian_gradient_counts.tsv", sep="\t", index=False)
pd.DataFrame(combos).to_csv(out / "exact_species_combination_counts.tsv", sep="\t", index=False)
assert (
    next(x["n_mammalian"] for x in grade if x["mode"] == "at_least" and x["n_mammals"] == 4) == 695
)
for folder in sorted((prev / "direction_audit").iterdir()):
    if not folder.is_dir():
        continue
    tissue = folder.name
    a = pd.read_csv(folder / "refit_20vs100.tsv", sep="\t").dropna(
        subset=["log2FoldChange", "stat", "pvalue", "padj"]
    )
    b = pd.read_csv(folder / "refit_50vs100.tsv", sep="\t").dropna(
        subset=["log2FoldChange", "stat", "pvalue", "padj"]
    )
    z = a.merge(b, on="gene_id", suffixes=("_20", "_50"), validate="one_to_one").merge(
        m, on="gene_id", validate="one_to_one"
    )
    z["age_log2FC_20"] = -z.log2FoldChange_20
    z["age_log2FC_50"] = -z.log2FoldChange_50
    z["tissue"] = tissue
    keepcols = [
        "tissue",
        "gene_id",
        "mouse_entrez_id",
        "chicken_symbol",
        "mouse_symbol",
        "support_n",
        "direction",
        "significant_conflict",
        "support_species",
        "spans_rodents_primates",
        "age_log2FC_20",
        "padj_20",
        "age_log2FC_50",
        "padj_50",
    ]
    eligible.append(z[keepcols])
    same = ((z.direction == "up") & (z.age_log2FC_20 > 1) & (z.age_log2FC_50 > 1)) | (
        (z.direction == "down") & (z.age_log2FC_20 < -1) & (z.age_log2FC_50 < -1)
    )
    good = z[
        (z.padj_20 < 0.05)
        & (z.padj_50 < 0.05)
        & same
        & ~z.significant_conflict
        & (z.support_n >= 2)
    ].copy()
    pairs.append(good[keepcols])
    for k in [2, 3, 4]:
        q = good[good.support_n >= k]
        tissue_rows.append(
            dict(
                tissue=tissue,
                n_mammals=k,
                n_eligible_background=len(z),
                n_eligible_tier=int(((z.support_n >= k) & ~z.significant_conflict).sum()),
                n_candidates=len(q),
                n_up=int((q.direction == "up").sum()),
                n_down=int((q.direction == "down").sum()),
            )
        )
    print("GRADIENT", tissue, flush=True)
elig = pd.concat(eligible, ignore_index=True)
pr = pd.concat(pairs, ignore_index=True)
elig.to_csv(out / "tissue_eligible_genes.tsv.gz", sep="\t", index=False)
pr.to_csv(out / "chicken_supported_gene_tissue_pairs.tsv", sep="\t", index=False)
pd.DataFrame(tissue_rows).to_csv(out / "tissue_gradient_counts.tsv", sep="\t", index=False)
group = [
    "gene_id",
    "mouse_entrez_id",
    "chicken_symbol",
    "mouse_symbol",
    "support_n",
    "direction",
    "support_species",
    "spans_rodents_primates",
]
gs = (
    pr.groupby(group, dropna=False)
    .agg(n_tissues=("tissue", "nunique"), tissues=("tissue", lambda x: ";".join(sorted(x))))
    .reset_index()
    .sort_values(["n_tissues", "gene_id"], ascending=[False, True])
)
gs.to_csv(out / "chicken_supported_gene_summary.tsv", sep="\t", index=False)
counts = []
for k in [2, 3, 4]:
    q = gs[gs.support_n >= k]
    q.to_csv(out / f"chicken_supported_at_least{k}.tsv", sep="\t", index=False)
    exact = gs[gs.support_n == k]
    exact.to_csv(out / f"chicken_supported_exactly{k}.tsv", sep="\t", index=False)
    counts.append(
        dict(
            n_mammals=k,
            n_unique_genes=len(q),
            n_up=int((q.direction == "up").sum()),
            n_down=int((q.direction == "down").sum()),
            n_gene_tissue_pairs=int((pr.support_n >= k).sum()),
            n_in_3plus_tissues=int((q.n_tissues >= 3).sum()),
            n_in_5plus_tissues=int((q.n_tissues >= 5).sum()),
            n_exactly_k=len(exact),
        )
    )
old = pd.read_csv(prev / "results/robust_conserved_gene_tissue_pairs.tsv", sep="\t")
assert set(zip(old.tissue, old.gene_id)) == set(
    zip(pr.loc[pr.support_n == 4, "tissue"], pr.loc[pr.support_n == 4, "gene_id"])
)
assert len(eligible) == 33
pd.DataFrame(counts).to_csv(out / "chicken_gradient_counts.tsv", sep="\t", index=False)
summary = dict(
    common_mammal_universe=len(ref),
    mapped_universe=len(m),
    n_significant_conflicts=int(ref.significant_conflict.sum()),
    mammal=grade,
    chicken=counts,
    previous_four_species_reproduced=True,
)
(out / "gradient_summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2), flush=True)
