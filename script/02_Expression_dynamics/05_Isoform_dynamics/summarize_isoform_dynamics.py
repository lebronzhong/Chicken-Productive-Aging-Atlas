#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 02_Expression_dynamics
# Purpose: Summarize expression-cluster profiles and within-gene response patterns.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

import pathlib
import json
import numpy as np
import pandas as pd
from scipy import stats

out = pathlib.Path(__file__).resolve().parent
ages = [20, 50, 80, 100]
tx = pd.read_csv(
    out / "02_overall/transcript_parent_four_age_results.tsv.gz", sep="\t", low_memory=False
)
old = pd.read_csv(out / "05_clusters/canonical_old_cluster_assignments.tsv.gz", sep="\t")
members = tx.loc[tx.common_background & tx.age_associated].merge(
    old, on=["tissue", "transcript_id"], how="left", validate="one_to_one"
)
mapped = members.loc[members.old_cluster.notna() & members.distinct_old_clusters.eq(1)].copy()
values = np.log2(mapped[[f"median_{age}" for age in ages]].to_numpy() + 0.5)
deviation = values.std(axis=1)
valid = deviation > 0
mapped = mapped.loc[valid].copy()
values = values[valid]
standardized = (values - values.mean(axis=1, keepdims=True)) / values.std(axis=1, keepdims=True)
for age_index, age in enumerate(ages):
    mapped[f"z_{age}"] = standardized[:, age_index]
columns = [
    "tissue",
    "transcript_id",
    "gene_id",
    "symbol",
    "class",
    "old_cluster",
    "age_q",
    "parent_age_q",
] + [f"z_{age}" for age in ages]
mapped[columns].to_csv(
    out / "05_clusters/primary_age_associated_members.tsv.gz", sep="\t", index=False
)
long = mapped.melt(
    id_vars=["old_cluster", "transcript_id", "tissue"],
    value_vars=[f"z_{age}" for age in ages],
    var_name="age",
    value_name="z",
)
long["age"] = long.age.str.replace("z_", "", regex=False).astype(int)
profile = (
    long.groupby(["old_cluster", "age"])
    .agg(
        median=("z", "median"),
        q25=("z", lambda values: values.quantile(0.25)),
        q75=("z", lambda values: values.quantile(0.75)),
        events=("z", "size"),
    )
    .reset_index()
)
profile.to_csv(out / "05_clusters/primary_ten_cluster_profiles.tsv", sep="\t", index=False)
breadth = (
    mapped.groupby("transcript_id")
    .tissue.nunique()
    .value_counts()
    .sort_index()
    .rename_axis("n_tissues")
    .rename("transcripts")
    .reset_index()
)
breadth["cumulative_fraction"] = breadth.transcripts.cumsum() / breadth.transcripts.sum()
breadth.to_csv(out / "05_clusters/primary_cluster_tissue_breadth.tsv", sep="\t", index=False)
mapped.groupby(["old_cluster", "tissue"]).size().rename("events").reset_index().to_csv(
    out / "05_clusters/primary_cluster_tissue_composition.tsv", sep="\t", index=False
)
record = dict(
    primary_age_associated_common_events=len(members),
    mapped_to_unambiguous_original_clusters=len(mapped),
    primary_events_without_original_membership=int(members.old_cluster.isna().sum()),
    conflicting_original_membership=int(members.distinct_old_clusters.gt(1).sum()),
    method="Retain original cluster membership; restrict to current canonical common-background age-associated events; recompute median/IQR of row-standardized log2(age median normalized count +0.5). Not de novo clustering.",
)
json.dump(record, open(out / "05_clusters/primary_cluster_audit.json", "w"), indent=2)
within = pd.read_csv(out / "03_within_gene/gene_tissue_tests.tsv.gz", sep="\t")
tx_annotation = tx.drop_duplicates("transcript_id").set_index("transcript_id")["class"]
within["class_A"] = within.transcript_A.map(tx_annotation)
within["class_B"] = within.transcript_B.map(tx_annotation)
within["pair_class"] = [
    " + ".join(sorted([str(first), str(second)]))
    for first, second in zip(within.class_A, within.class_B)
]
within.loc[within.divergent].groupby("pair_class").agg(
    events=("gene_id", "size"),
    unique_genes=("gene_id", "nunique"),
    parent_not_detected=("parent_not_detected", "sum"),
).reset_index().to_csv(
    out / "03_within_gene/representative_pair_class_summary.tsv", sep="\t", index=False
)
selected_within = pd.read_csv(out / "06_examples/reviewed_within_cases.tsv", sep="\t")
selected_cross = pd.read_csv(out / "06_examples/reviewed_cross_cases.tsv", sep="\t")
samples = pd.read_csv(out / "01_inputs/sample_order.tsv", sep="\t")["sample"].tolist()
records = []
parent_intervals = []
for group, frame in [("within", selected_within), ("cross", selected_cross)]:
    for case_index, row in frame.iterrows():
        first_tissue = row.tissue if group == "within" else row.tissue_A
        second_tissue = first_tissue if group == "within" else row.tissue_B
        first_id = row.transcript_A if group == "within" else row.transcript_id
        second_id = row.transcript_B if group == "within" else row.transcript_id
        first = (
            pd.read_csv(
                out / f"01_inputs/{first_tissue}/tx_log_expression.tsv.gz", sep="\t", index_col=0
            )
            .loc[first_id, samples]
            .to_numpy()
        )
        second = (
            pd.read_csv(
                out / f"01_inputs/{second_tissue}/tx_log_expression.tsv.gz", sep="\t", index_col=0
            )
            .loc[second_id, samples]
            .to_numpy()
        )
        difference = (first - second).reshape(4, 5)
        means = difference.mean(axis=1)
        minimum_age, maximum_age = np.argmin(means), np.argmax(means)
        leave_out = []
        for age_index in range(4):
            for bird in range(5):
                changed = means.copy()
                changed[age_index] = np.delete(difference[age_index], bird).mean()
                leave_out.append(changed[maximum_age] - changed[minimum_age])
        parent_rows = tx.loc[tx.transcript_id.eq(first_id) & tx.tissue.eq(first_tissue)]
        records.append(
            dict(
                case=f"{group}_{case_index+1}",
                symbol=row.symbol,
                tissue_A=first_tissue,
                tissue_B=second_tissue,
                transcript_A=first_id,
                transcript_B=second_id,
                minimum_difference_age=ages[minimum_age],
                maximum_difference_age=ages[maximum_age],
                difference_range=float(np.ptp(means)),
                minimum_leave_one_bird_out_effect=float(min(leave_out)),
                all_leave_one_bird_out_direction_consistent=bool(min(leave_out) > 0),
                first_isoform_minimum_age_support=int(
                    parent_rows[[f"support_{age}" for age in ages]].min(axis=1).iloc[0]
                ),
                first_isoform_minimum_age_median=float(
                    parent_rows[[f"median_{age}" for age in ages]].min(axis=1).iloc[0]
                ),
            )
        )
        parent_first = (
            pd.read_csv(
                out / f"01_inputs/{first_tissue}/gene_log_expression.tsv.gz", sep="\t", index_col=0
            )
            .loc[row.gene_id, samples]
            .to_numpy()
        )
        if group == "cross":
            parent_second = (
                pd.read_csv(
                    out / f"01_inputs/{second_tissue}/gene_log_expression.tsv.gz",
                    sep="\t",
                    index_col=0,
                )
                .loc[row.gene_id, samples]
                .to_numpy()
            )
            parent_values = (parent_first - parent_second).reshape(4, 5)
        else:
            parent_values = parent_first.reshape(4, 5)
        for older_index in range(1, 4):
            for younger_index in range(older_index):
                older, younger = parent_values[older_index], parent_values[younger_index]
                first_variance, second_variance = older.var(ddof=1) / 5, younger.var(ddof=1) / 5
                variance = first_variance + second_variance
                degrees = variance**2 / (first_variance**2 / 4 + second_variance**2 / 4)
                effect = older.mean() - younger.mean()
                halfwidth = stats.t.ppf(0.975, degrees) * np.sqrt(variance)
                parent_intervals.append(
                    dict(
                        case=f"{group}_{case_index+1}",
                        symbol=row.symbol,
                        contrast=f"{ages[older_index]}vs{ages[younger_index]}",
                        parent_effect=float(effect),
                        CI95_lower=float(effect - halfwidth),
                        CI95_upper=float(effect + halfwidth),
                        scope="Descriptive unadjusted 95% interval; not equivalence or simultaneous inference",
                        effect_type=(
                            "within-bird tissue-difference age contrast"
                            if group == "cross"
                            else "parent log-expression age contrast"
                        ),
                    )
                )
pd.DataFrame(records).to_csv(
    out / "06_examples/case_quality_and_leave_one_out.tsv", sep="\t", index=False
)
pd.DataFrame(parent_intervals).to_csv(
    out / "06_examples/parent_effect_intervals.tsv", sep="\t", index=False
)
print("ADDITIONAL COMPLETE", record, flush=True)
