#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 02_Expression_dynamics
# Purpose: Conservation endpoint.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

import pathlib
import numpy as np
import pandas as pd
from scipy import stats

out = pathlib.Path(__file__).resolve().parent
folder = out / "09_conservation"
folder.mkdir(exist_ok=True)
source = out.parent / "expression_models/07_conservation/all_gene_tissue_evidence.tsv.gz"
legacy = pd.read_csv(source, sep="\t", low_memory=False)
columns = ["tissue", "gene_id", "support_n", "direction", "significant_conflict", "support_species"]
signatures = legacy[columns].drop_duplicates()
assert not signatures.duplicated(["tissue", "gene_id"]).any()
signatures["significant_conflict"] = (
    signatures.significant_conflict.astype(str).str.lower().eq("true")
)
signatures["direction_sign"] = signatures.direction.map({"up": 1, "down": -1})
signatures.to_csv(folder / "mammalian_signature_input.tsv.gz", sep="\t", index=False)
tx = pd.read_csv(
    out / "02_overall/transcript_parent_four_age_results.tsv.gz", sep="\t", low_memory=False
)
gene = pd.read_csv(out / "02_overall/gene_four_age_results.tsv.gz", sep="\t")
tag = "100vs20"
tx = tx.loc[tx.common_background & tx[f"{tag}_q"].notna()].copy()
tx["up"] = tx[f"{tag}_q"].lt(0.05) & tx[f"{tag}_log2FoldChange"].ge(1)
tx["down"] = tx[f"{tag}_q"].lt(0.05) & tx[f"{tag}_log2FoldChange"].le(-1)
tx["mean_age_median"] = tx[[f"median_{age}" for age in [20, 50, 80, 100]]].mean(axis=1)
directions = (
    tx.groupby(["tissue", "gene_id"])
    .agg(
        any_isoform_up=("up", "any"),
        any_isoform_down=("down", "any"),
        eligible_isoforms=("transcript_id", "size"),
    )
    .reset_index()
)
representatives = (
    tx.sort_values(["mean_age_median", "transcript_id"], ascending=[False, True])
    .drop_duplicates(["tissue", "gene_id"])[["tissue", "gene_id", "transcript_id", "up", "down"]]
    .rename(
        columns={
            "transcript_id": "representative_transcript",
            "up": "representative_up",
            "down": "representative_down",
        }
    )
)
background = (
    gene.loc[gene[f"{tag}_q"].notna()]
    .merge(signatures, on=["tissue", "gene_id"], validate="one_to_one")
    .merge(directions, on=["tissue", "gene_id"], validate="one_to_one")
    .merge(representatives, on=["tissue", "gene_id"], validate="one_to_one")
)
background["gene_up"] = background[f"{tag}_q"].lt(0.05) & background[f"{tag}_log2FoldChange"].ge(1)
background["gene_down"] = background[f"{tag}_q"].lt(0.05) & background[f"{tag}_log2FoldChange"].le(
    -1
)
background["gene_agrees"] = np.where(
    background.direction_sign.eq(1), background.gene_up, background.gene_down
)
background["isoform_agrees"] = np.where(
    background.direction_sign.eq(1), background.any_isoform_up, background.any_isoform_down
)
background["category"] = np.select(
    [
        background.gene_agrees & background.isoform_agrees,
        background.gene_agrees,
        background.isoform_agrees,
    ],
    ["Both", "Gene only", "Isoform only"],
    default="Neither",
)
background.to_csv(folder / "endpoint_common_background.tsv.gz", sep="\t", index=False)
generator = np.random.default_rng(20261023)
records, relations = [], []
for level in [2, 3, 4]:
    subset = background.loc[background.support_n.ge(level) & ~background.significant_conflict]
    assert subset.groupby("gene_id").direction_sign.nunique().max() == 1
    grouped = subset.groupby("gene_id").agg(
        direction=("direction_sign", "first"),
        gene_up=("gene_up", "sum"),
        gene_down=("gene_down", "sum"),
        any_isoform_up=("any_isoform_up", "sum"),
        any_isoform_down=("any_isoform_down", "sum"),
        representative_up=("representative_up", "sum"),
        representative_down=("representative_down", "sum"),
    )
    mammal_up = grouped.direction.to_numpy() > 0
    methods = ["gene", "any_isoform", "representative"]
    up = grouped[[method + "_up" for method in methods]].to_numpy(dtype=float)
    down = grouped[[method + "_down" for method in methods]].to_numpy(dtype=float)
    observed = np.where(mammal_up[:, None], up, down).sum(axis=0)
    null = np.empty((10000, 3))
    for permutation in range(10000):
        null[permutation] = np.where(generator.permutation(mammal_up)[:, None], up, down).sum(
            axis=0
        )
    for method_index, method in enumerate(methods):
        values = null[:, method_index]
        records.append(
            dict(
                contrast="100vs20_from_four_age_model",
                mammal_support=level,
                method=method,
                denominator_gene_tissue_events=len(subset),
                unique_genes=len(grouped),
                observed_events=int(observed[method_index]),
                null_mean=float(values.mean()),
                null_lower=float(np.quantile(values, 0.025)),
                null_upper=float(np.quantile(values, 0.975)),
                fold_over_null=float(observed[method_index] / values.mean()),
                permutation_p=float((1 + (values >= observed[method_index]).sum()) / 10001),
            )
        )
    for category, frame in subset.groupby("category"):
        relations.append(
            dict(
                mammal_support=level,
                category=category,
                events=len(frame),
                denominator=len(subset),
                fraction=len(frame) / len(subset),
            )
        )
result = pd.DataFrame(records)
result["q_9_tests"] = stats.false_discovery_control(result.permutation_p, method="bh")
result.to_csv(folder / "endpoint_direction_permutation.tsv", sep="\t", index=False)
pd.DataFrame(relations).to_csv(folder / "endpoint_relationships.tsv", sep="\t", index=False)
print(result.to_string(index=False))
