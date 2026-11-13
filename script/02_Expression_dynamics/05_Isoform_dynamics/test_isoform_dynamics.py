#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 02_Expression_dynamics
# Purpose: Isoform dynamics.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

import itertools
import json
import pathlib
import numpy as np
import pandas as pd
from scipy import stats

out = pathlib.Path(__file__).resolve().parent
ages = [20, 50, 80, 100]
contrasts = [(50, 20), (80, 20), (100, 20), (80, 50), (100, 50), (100, 80)]
tags = [f"{older}vs{younger}" for older, younger in contrasts]
tissues = (out / "01_inputs/tissues.txt").read_text().split()
sample_order = pd.read_csv(out / "01_inputs/sample_order.tsv", sep="\t")
samples = sample_order["sample"].tolist()
checks = []


def check(name, condition):
    checks.append({"check": name, "passed": bool(condition)})
    assert condition, name


def bh(values):
    values = np.asarray(values, dtype=float)
    valid = np.isfinite(values)
    result = np.full(values.shape, np.nan)
    if valid.any():
        result[valid] = stats.false_discovery_control(values[valid], method="bh")
    return result


def write(frame, filename):
    frame.to_csv(out / filename, sep="\t", index=False)


def welch_age(values):
    cube = np.asarray(values, dtype=float).reshape(-1, 4, 5)
    means = cube.mean(axis=2)
    variances = np.maximum(cube.var(axis=2, ddof=1), 1e-8)
    weights = 5 / variances
    weight_sum = weights.sum(axis=1)
    center = (weights * means).sum(axis=1) / weight_sum
    correction = (((1 - weights / weight_sum[:, None]) ** 2) / 4).sum(axis=1)
    statistic = (
        (weights * (means - center[:, None]) ** 2).sum(axis=1) / 3 / (1 + 4 / 15 * correction)
    )
    pvalue = stats.f.sf(statistic, 3, 5 / correction)
    return pvalue, np.ptp(means, axis=1), means


mapping = pd.read_csv(out / "01_inputs/canonical_transcript_mapping.tsv.gz", sep="\t")[
    ["transcript_id", "gene_id"]
]
quality = pd.read_csv(out / "01_inputs/quality.tsv.gz", sep="\t", low_memory=False).drop(
    columns="gene_id"
)
symbols = pd.read_csv(out / "01_inputs/gene_symbols.tsv", sep="\t").drop_duplicates("gene_id")
check("canonical_transcript_unique", mapping.transcript_id.is_unique)
check("quality_annotation_unique", quality.transcript_id.is_unique)
manifest = pd.read_csv(out / "01_inputs/sample_manifest.tsv", sep="\t")
check(
    "33_tissues_20_samples", len(manifest) == 660 and manifest.groupby("tissue").size().eq(20).all()
)
check("four_natural_ages_only", set(manifest.age) == set(ages))
check(
    "five_birds_per_tissue_age",
    manifest.groupby(["tissue", "age"]).individual.nunique().eq(5).all(),
)
check(
    "same_bird_labels_across_tissues", manifest.groupby("individual").tissue.nunique().eq(33).all()
)
if (out / "02_overall/overall_summary.tsv").exists():
    gene = pd.read_csv(out / "02_overall/gene_four_age_results.tsv.gz", sep="\t", low_memory=False)
    tx = pd.read_csv(
        out / "02_overall/transcript_parent_four_age_results.tsv.gz", sep="\t", low_memory=False
    )
else:
    tables = {}
    for feature, identifier in [("gene", "gene_id"), ("tx", "transcript_id")]:
        frame = pd.concat(
            [
                pd.read_csv(out / f"01_inputs/{tissue}/{feature}_four_age.tsv.gz", sep="\t").assign(
                    tissue=tissue
                )
                for tissue in tissues
            ],
            ignore_index=True,
        ).rename(columns={"feature_id": identifier})
        frame["age_q"] = bh(frame.overall_pvalue)
        for tag in tags:
            frame[f"{tag}_q"] = bh(frame[f"{tag}_pvalue"])
        frame["max_abs_lfc"] = frame[[f"{tag}_log2FoldChange" for tag in tags]].abs().max(axis=1)
        frame["age_associated"] = frame.age_q.lt(0.05) & frame.max_abs_lfc.ge(1)
        check(
            feature + "_unique_tissue_features", not frame.duplicated(["tissue", identifier]).any()
        )
        tables[feature] = frame
    gene = tables["gene"].merge(symbols, on="gene_id", how="left", validate="many_to_one")
    tx = (
        tables["tx"]
        .merge(mapping, on="transcript_id", how="left", validate="many_to_one")
        .merge(quality, on="transcript_id", how="left", validate="many_to_one")
        .merge(symbols, on="gene_id", how="left", validate="many_to_one")
    )
    parent_columns = (
        ["age_q", "max_abs_lfc", "age_associated"]
        + [f"support_{age}" for age in ages]
        + [f"median_{age}" for age in ages]
    )
    tx = tx.merge(
        gene[["gene_id", "tissue"] + parent_columns].rename(
            columns={column: "parent_" + column for column in parent_columns}
        ),
        on=["gene_id", "tissue"],
        how="left",
        validate="many_to_one",
    )
    joint_support = np.column_stack(
        [tx[f"support_{age}"].ge(3) & tx[f"parent_support_{age}"].ge(3) for age in ages]
    )
    tx["common_background"] = tx.age_q.notna() & tx.parent_age_q.notna() & joint_support.any(axis=1)
    tx["parent_not_detected"] = tx.common_background & tx.age_associated & tx.parent_age_q.ge(0.05)
    tx["all_age_supported"] = np.column_stack(
        [tx[f"support_{age}"].ge(3) & tx[f"parent_support_{age}"].ge(3) for age in ages]
    ).all(axis=1)
    tx["structural_quality"] = tx.sqanti_all_canonical.eq("canonical") & tx.sqanti_RTS_stage.astype(
        str
    ).str.lower().eq("false")
    tx["strict_quality"] = tx.structural_quality & np.column_stack(
        [
            tx[f"support_{age}"].ge(4)
            & tx[f"median_{age}"].ge(5)
            & tx[f"parent_support_{age}"].ge(3)
            for age in ages
        ]
    ).all(axis=1)
    write(gene, "02_overall/gene_four_age_results.tsv.gz")
    write(tx, "02_overall/transcript_parent_four_age_results.tsv.gz")
tx["parent_age_associated"] = tx.parent_age_associated.fillna(False).astype(bool)


def summarize_tx(frame, label):
    common = frame.loc[frame.common_background]
    associated = common.loc[common.age_associated]
    hidden = common.loc[common.parent_not_detected]
    return dict(
        tissue=label,
        common_tx_tissue_events=len(common),
        age_associated_tx_events=len(associated),
        parent_not_detected_events=len(hidden),
        hidden_fraction_among_associated=(
            len(hidden) / len(associated) if len(associated) else np.nan
        ),
        unique_associated_transcripts=associated.transcript_id.nunique(),
        unique_hidden_transcripts=hidden.transcript_id.nunique(),
        unique_hidden_genes=hidden.gene_id.nunique(),
    )


overall = summarize_tx(tx, "ALL")
write(
    pd.DataFrame([summarize_tx(frame, tissue) for tissue, frame in tx.groupby("tissue")]),
    "02_overall/per_tissue_summary.tsv",
)
write(pd.DataFrame([overall]), "02_overall/overall_summary.tsv")
class_records = [
    dict(structural_class=label, **summarize_tx(frame, "ALL"))
    for label, frame in tx.groupby("class")
]
write(pd.DataFrame(class_records), "02_overall/structural_class_summary.tsv")
thresholds = []
for subset, mask in [
    ("common", tx.common_background),
    ("all_four_ages_supported", tx.common_background & tx.all_age_supported),
    ("strict_quality", tx.common_background & tx.strict_quality),
]:
    for qcut in [0.01, 0.05]:
        for effect in [0, 0.5, 1, 1.5]:
            detected = mask & tx.age_q.lt(qcut) & tx.max_abs_lfc.ge(effect)
            hidden = detected & tx.parent_age_q.ge(qcut)
            thresholds.append(
                dict(
                    subset=subset,
                    q_cutoff=qcut,
                    effect_cutoff=effect,
                    denominator=int(detected.sum()),
                    numerator=int(hidden.sum()),
                    fraction=float(hidden.sum() / detected.sum()) if detected.sum() else np.nan,
                )
            )
write(pd.DataFrame(thresholds), "02_overall/threshold_background_summary.tsv")
common = tx.loc[tx.common_background].copy()
gene_events = (
    common.groupby(["tissue", "gene_id"])
    .agg(
        n_isoforms=("transcript_id", "nunique"),
        n_age_associated_isoforms=("age_associated", "sum"),
        gene_q=("parent_age_q", "first"),
        gene_associated=("parent_age_associated", "first"),
        gene_max_lfc=("parent_max_abs_lfc", "first"),
        symbol=("symbol", "first"),
    )
    .reset_index()
)
gene_events["status"] = np.select(
    [
        gene_events.gene_associated & gene_events.n_age_associated_isoforms.gt(0),
        gene_events.gene_associated,
        gene_events.n_age_associated_isoforms.gt(0),
    ],
    ["Both", "Gene only", "Isoform only"],
    default="Neither",
)
write(gene_events, "02_overall/common_gene_tissue_events.tsv.gz")
relationship = gene_events.groupby(["tissue", "status"]).size().rename("events").reset_index()
relationship["denominator"] = relationship.groupby("tissue").events.transform("sum")
relationship["fraction"] = relationship.events / relationship.denominator
write(relationship, "02_overall/gene_isoform_relationship_per_tissue.tsv")
check("relationship_denominator", relationship.events.sum() == len(gene_events))
contrast_records = []
for tag in tags:
    parent_q = gene.set_index(["tissue", "gene_id"])[f"{tag}_q"]
    paired_q = common.set_index(["tissue", "gene_id"]).index.map(parent_q)
    contrast_common = common[f"{tag}_q"].notna().to_numpy() & np.isfinite(paired_q)
    differential = (
        contrast_common
        & common[f"{tag}_q"].lt(0.05).to_numpy()
        & common[f"{tag}_log2FoldChange"].abs().ge(1).to_numpy()
    )
    hidden = differential & (np.asarray(paired_q) >= 0.05)
    contrast_records.append(
        dict(
            contrast=tag,
            testable_events=int(contrast_common.sum()),
            DE_events=int(differential.sum()),
            parent_not_detected=int(hidden.sum()),
            fraction=float(hidden.sum() / differential.sum()),
        )
    )
write(pd.DataFrame(contrast_records), "02_overall/six_contrasts_from_four_age_model.tsv")
print("OVERALL", json.dumps(overall), flush=True)
tx_ids = np.sort(common.transcript_id.unique())
gene_ids = np.sort(gene.gene_id.unique())
tx_lookup = pd.Series(np.arange(len(tx_ids)), index=tx_ids)
gene_lookup = pd.Series(np.arange(len(gene_ids)), index=gene_ids)
parents = (
    mapping.set_index("transcript_id").reindex(tx_ids).gene_id.map(gene_lookup).to_numpy(dtype=int)
)
tx_expression = np.full((len(tx_ids), len(tissues), 20), np.nan, dtype=np.float32)
gene_expression = np.full((len(gene_ids), len(tissues), 20), np.nan, dtype=np.float32)
tx_eligible = np.zeros((len(tx_ids), len(tissues)), dtype=bool)
gene_eligible = np.zeros((len(gene_ids), len(tissues)), dtype=bool)
tx_support = np.zeros((len(tx_ids), len(tissues), 4), dtype=bool)
gene_support = np.zeros((len(gene_ids), len(tissues), 4), dtype=bool)
tx_associated = np.zeros((len(tx_ids), len(tissues)), dtype=bool)
tx_strict = np.zeros((len(tx_ids), len(tissues)), dtype=bool)
tx_lfc = np.full((len(tx_ids), len(tissues), 6), np.nan, dtype=np.float32)
tx_contrast_de = np.zeros((len(tx_ids), len(tissues), 6), dtype=bool)
within_rows = []
for tissue_index, tissue in enumerate(tissues):
    frame = common.loc[common.tissue == tissue].copy()
    frame["index"] = frame.transcript_id.map(tx_lookup)
    indices = frame["index"].to_numpy()
    expression = pd.read_csv(
        out / f"01_inputs/{tissue}/tx_log_expression.tsv.gz", sep="\t", index_col=0
    )
    values = expression.reindex(frame.transcript_id)[samples].to_numpy(dtype=np.float32)
    check(tissue + "_complete_tx_expression", np.isfinite(values).all())
    tx_expression[indices, tissue_index, :] = values
    tx_eligible[indices, tissue_index] = True
    tx_associated[indices, tissue_index] = frame.age_associated
    tx_strict[indices, tissue_index] = frame.strict_quality
    tx_support[indices, tissue_index, :] = np.column_stack(
        [frame[f"support_{age}"].ge(3) & frame[f"parent_support_{age}"].ge(3) for age in ages]
    )
    tx_lfc[indices, tissue_index, :] = frame[[f"{tag}_log2FoldChange" for tag in tags]]
    tx_contrast_de[indices, tissue_index, :] = np.column_stack(
        [frame[f"{tag}_q"].lt(0.05) & frame[f"{tag}_log2FoldChange"].abs().ge(1) for tag in tags]
    )
    gene_frame = gene.loc[gene.tissue == tissue]
    gene_indices = gene_frame.gene_id.map(gene_lookup).to_numpy()
    expression = pd.read_csv(
        out / f"01_inputs/{tissue}/gene_log_expression.tsv.gz", sep="\t", index_col=0
    )
    gene_expression[gene_indices, tissue_index, :] = expression.reindex(gene_frame.gene_id)[
        samples
    ].to_numpy(dtype=np.float32)
    gene_eligible[gene_indices, tissue_index] = gene_frame.age_q.notna()
    gene_support[gene_indices, tissue_index, :] = np.column_stack(
        [gene_frame[f"support_{age}"].ge(3) for age in ages]
    )
    group_ids, first_ids, second_ids = [], [], []
    groups = list(frame.groupby("gene_id", sort=False))
    possible = {}
    for group_index, (gene_id, group) in enumerate(groups):
        if len(group) < 2:
            continue
        group_indices = group["index"].to_numpy()
        pairs = np.asarray(list(itertools.combinations(group_indices, 2)))
        group_ids.extend([group_index] * len(pairs))
        first_ids.extend(pairs[:, 0])
        second_ids.extend(pairs[:, 1])
        possible[group_index] = len(pairs)
    first_ids, second_ids, group_ids = map(np.asarray, [first_ids, second_ids, group_ids])
    tested = np.any(
        tx_support[first_ids, tissue_index] & tx_support[second_ids, tissue_index], axis=1
    )
    difference = tx_expression[first_ids, tissue_index, :].astype(float) - tx_expression[
        second_ids, tissue_index, :
    ].astype(float)
    pvalue, effect, means = welch_age(difference)
    pvalue[~tested] = 1
    associated = tx_associated[first_ids, tissue_index] | tx_associated[second_ids, tissue_index]
    eligible_signal = associated & (effect >= 1) & tested
    selected_p = np.where(eligible_signal, pvalue, 1)
    pair_table = pd.DataFrame(
        dict(
            group=group_ids,
            first=first_ids,
            second=second_ids,
            pvalue=pvalue,
            selected_p=selected_p,
            effect=effect,
            tested=tested,
        )
    )
    best = pair_table.sort_values("selected_p").drop_duplicates("group").set_index("group")
    strict_pair = tx_strict[first_ids, tissue_index] & tx_strict[second_ids, tissue_index]
    pair_table["strict_p"] = np.where(strict_pair & eligible_signal, pvalue, 1)
    strict_best = pair_table.sort_values("strict_p").drop_duplicates("group").set_index("group")
    tested_counts = pair_table.groupby("group").tested.sum()
    for group_index in best.index:
        gene_id, group = groups[group_index]
        row = best.loc[group_index]
        first, second = int(row["first"]), int(row["second"])
        opposite = (
            tx_contrast_de[first, tissue_index]
            & tx_contrast_de[second, tissue_index]
            & (tx_lfc[first, tissue_index] * tx_lfc[second, tissue_index] < 0)
        )
        profile_first = tx_expression[first, tissue_index].reshape(4, 5).mean(axis=1)
        profile_second = tx_expression[second, tissue_index].reshape(4, 5).mean(axis=1)
        strict_row = strict_best.loc[group_index]
        within_rows.append(
            dict(
                tissue=tissue,
                gene_id=gene_id,
                symbol=group.symbol.iloc[0],
                n_isoforms=len(group),
                n_possible_pairs=possible[group_index],
                n_tested_pairs=int(tested_counts.loc[group_index]),
                transcript_A=tx_ids[first],
                transcript_B=tx_ids[second],
                pair_p=row.pvalue,
                difference_range=row.effect,
                event_p=min(1, row.selected_p * possible[group_index]),
                gene_q=group.parent_age_q.iloc[0],
                gene_associated=group.parent_age_associated.iloc[0],
                both_isoforms_age_associated=bool(
                    tx_associated[first, tissue_index] & tx_associated[second, tissue_index]
                ),
                opposing_contrasts=",".join(np.asarray(tags)[opposite]),
                peak_age_A=ages[np.argmax(profile_first)],
                peak_age_B=ages[np.argmax(profile_second)],
                strict_pair=bool(tx_strict[first, tissue_index] & tx_strict[second, tissue_index]),
                strict_transcript_A=tx_ids[int(strict_row["first"])],
                strict_transcript_B=tx_ids[int(strict_row["second"])],
                strict_event_p=min(1, strict_row.strict_p * possible[group_index]),
                strict_eligible_isoforms=int(group.strict_quality.sum()),
            )
        )
    print("WITHIN", tissue, len(pair_table), flush=True)
within = pd.DataFrame(within_rows)
within["event_q"] = bh(within.event_p)
within["strict_event_q"] = bh(within.strict_event_p)
within["divergent"] = within.event_q.lt(0.05)
within["parent_not_detected"] = within.divergent & within.gene_q.ge(0.05)
within["opposing"] = within.divergent & within.opposing_contrasts.ne("")
within["pattern"] = np.select(
    [
        ~within.divergent,
        within.opposing,
        within.both_isoforms_age_associated & within.peak_age_A.ne(within.peak_age_B),
        within.both_isoforms_age_associated,
    ],
    [
        "No supported divergence",
        "Opposing supported contrast",
        "Different peak ages",
        "Both age-associated; other divergence",
    ],
    default="One age-associated isoform",
)
write(within, "03_within_gene/gene_tissue_tests.tsv.gz")
within_tissue = (
    within.groupby("tissue")
    .agg(
        eligible_gene_tissue_events=("gene_id", "size"),
        divergent_events=("divergent", "sum"),
        parent_not_detected_events=("parent_not_detected", "sum"),
        opposing_events=("opposing", "sum"),
        strict_eligible_events=(
            "strict_eligible_isoforms",
            lambda values: int((values >= 2).sum()),
        ),
        strict_divergent_events=("strict_event_q", lambda values: int((values < 0.05).sum())),
    )
    .reset_index()
)
within_tissue["divergent_fraction"] = (
    within_tissue.divergent_events / within_tissue.eligible_gene_tissue_events
)
within_tissue["hidden_fraction_among_divergent"] = (
    within_tissue.parent_not_detected_events / within_tissue.divergent_events
)
write(within_tissue, "03_within_gene/per_tissue_summary.tsv")
write(
    within.groupby(["pattern", "gene_associated"]).size().rename("events").reset_index(),
    "03_within_gene/pattern_summary.tsv",
)
check(
    "within_family_matches_multi_isoform_background",
    len(within) == int(gene_events.n_isoforms.ge(2).sum()),
)
print("WITHIN COMPLETE", len(within), int(within.divergent.sum()), flush=True)
pair_directory = out / "04_cross_tissue/pair_tests"
pair_directory.mkdir(exist_ok=True)
pair_definitions = list(itertools.combinations(range(len(tissues)), 2))
all_tx_p, all_gene_p, pair_counts = [], [], []
unique_min_p = np.ones(len(tx_ids))
unique_n_tests = np.zeros(len(tx_ids), dtype=int)
for pair_index, (first, second) in enumerate(pair_definitions):
    tx_rows = np.flatnonzero(
        tx_eligible[:, first]
        & tx_eligible[:, second]
        & np.any(tx_support[:, first] & tx_support[:, second], axis=1)
    )
    gene_rows = np.flatnonzero(
        gene_eligible[:, first]
        & gene_eligible[:, second]
        & np.any(gene_support[:, first] & gene_support[:, second], axis=1)
    )
    tx_p, tx_effect, tx_means = welch_age(
        tx_expression[tx_rows, first, :].astype(float)
        - tx_expression[tx_rows, second, :].astype(float)
    )
    gene_p, gene_effect, gene_means = welch_age(
        gene_expression[gene_rows, first, :].astype(float)
        - gene_expression[gene_rows, second, :].astype(float)
    )
    np.savez_compressed(
        pair_directory / f"pair_{pair_index:03d}.npz",
        tx_index=tx_rows,
        gene_index=gene_rows,
        tx_p=tx_p,
        tx_effect=tx_effect,
        tx_means=tx_means,
        gene_p=gene_p,
        gene_effect=gene_effect,
    )
    all_tx_p.append(tx_p)
    all_gene_p.append(gene_p)
    unique_min_p[tx_rows] = np.minimum(unique_min_p[tx_rows], tx_p)
    unique_n_tests[tx_rows] += 1
    pair_counts.append(
        dict(
            pair_index=pair_index,
            tissue_A=tissues[first],
            tissue_B=tissues[second],
            tested_transcripts=len(tx_rows),
            tested_genes=len(gene_rows),
        )
    )
    if pair_index % 80 == 0:
        print("CROSS TEST", pair_index, flush=True)
unique_q = np.full(len(tx_ids), np.nan)
unique_tested = unique_n_tests > 0
unique_q[unique_tested] = bh(
    np.minimum(1, unique_min_p[unique_tested] * unique_n_tests[unique_tested])
)
all_tx_q = bh(np.concatenate(all_tx_p))
all_gene_q = bh(np.concatenate(all_gene_p))
tx_offset = gene_offset = 0
cross_hits = []
per_tx = np.zeros((len(tx_ids), 4), dtype=int)
per_tissue_tested = np.zeros((len(tx_ids), len(tissues)), dtype=bool)
per_tissue_hit = np.zeros_like(per_tissue_tested)
per_tissue_hidden = np.zeros_like(per_tissue_tested)
for pair_index, (first, second) in enumerate(pair_definitions):
    with np.load(pair_directory / f"pair_{pair_index:03d}.npz") as source:
        data = {key: source[key] for key in source.files}
    tx_rows, gene_rows = data["tx_index"], data["gene_index"]
    tx_q = all_tx_q[tx_offset : tx_offset + len(tx_rows)]
    gene_q = all_gene_q[gene_offset : gene_offset + len(gene_rows)]
    tx_offset += len(tx_rows)
    gene_offset += len(gene_rows)
    parent_positions = (
        pd.Series(np.arange(len(gene_rows)), index=gene_rows).reindex(parents[tx_rows]).to_numpy()
    )
    check(f"cross_pair_{pair_index}_parent_available", np.isfinite(parent_positions).all())
    parent_q = gene_q[parent_positions.astype(int)]
    parent_effect = data["gene_effect"][parent_positions.astype(int)]
    age_associated = tx_associated[tx_rows, first] | tx_associated[tx_rows, second]
    hit = (tx_q < 0.05) & (unique_q[tx_rows] < 0.05) & (data["tx_effect"] >= 1) & age_associated
    opposite_matrix = (
        tx_contrast_de[tx_rows, first]
        & tx_contrast_de[tx_rows, second]
        & (tx_lfc[tx_rows, first] * tx_lfc[tx_rows, second] < 0)
    )
    opposite = hit & opposite_matrix.any(axis=1)
    hidden = hit & (parent_q >= 0.05)
    strict_pair = tx_strict[tx_rows, first] & tx_strict[tx_rows, second]
    per_tx[tx_rows, 0] += hit
    per_tx[tx_rows, 1] += opposite
    per_tx[tx_rows, 2] += hidden
    per_tx[tx_rows, 3] += hit & strict_pair
    for tissue_index in [first, second]:
        per_tissue_tested[tx_rows, tissue_index] = True
        per_tissue_hit[tx_rows[hit], tissue_index] = True
        per_tissue_hidden[tx_rows[hidden], tissue_index] = True
    pair_counts[pair_index].update(
        heterogeneous_transcripts=int(hit.sum()),
        opposing_transcripts=int(opposite.sum()),
        parent_not_detected=int(hidden.sum()),
        strict_eligible=int(strict_pair.sum()),
        strict_hits=int((hit & strict_pair).sum()),
    )
    np.savez_compressed(
        pair_directory / f"pair_{pair_index:03d}.npz", **data, tx_q=tx_q, gene_q=gene_q
    )
    if hit.any():
        rows = tx_rows[hit]
        cross_hits.append(
            pd.DataFrame(
                dict(
                    transcript_id=tx_ids[rows],
                    gene_id=gene_ids[parents[rows]],
                    tissue_A=tissues[first],
                    tissue_B=tissues[second],
                    pair_p=data["tx_p"][hit],
                    pair_q=tx_q[hit],
                    unique_transcript_q=unique_q[rows],
                    difference_range=data["tx_effect"][hit],
                    gene_pair_q=parent_q[hit],
                    gene_difference_range=parent_effect[hit],
                    parent_not_detected=hidden[hit],
                    opposing=opposite[hit],
                    strict_pair=strict_pair[hit],
                    opposing_contrasts=[
                        ",".join(np.asarray(tags)[mask]) for mask in opposite_matrix[hit]
                    ],
                )
            )
        )
cross = pd.concat(cross_hits, ignore_index=True) if cross_hits else pd.DataFrame()
cross = cross.merge(symbols, on="gene_id", how="left", validate="many_to_one").merge(
    quality, on="transcript_id", how="left", validate="many_to_one"
)
write(cross, "04_cross_tissue/significant_tissue_pairs.tsv.gz")
write(pd.DataFrame(pair_counts), "04_cross_tissue/tissue_pair_summary.tsv")
unique = pd.DataFrame(
    dict(
        transcript_id=tx_ids,
        gene_id=gene_ids[parents],
        tested_pairs=unique_n_tests,
        unique_transcript_q=unique_q,
        heterogeneous_pairs=per_tx[:, 0],
        opposing_pairs=per_tx[:, 1],
        parent_not_detected_pairs=per_tx[:, 2],
        strict_heterogeneous_pairs=per_tx[:, 3],
    )
)
write(unique, "04_cross_tissue/unique_transcript_tests.tsv.gz")
cross_tissue_summary = pd.DataFrame(
    dict(
        tissue=tissues,
        tested_transcripts=per_tissue_tested.sum(axis=0),
        heterogeneous_transcripts=per_tissue_hit.sum(axis=0),
        parent_not_detected_transcripts=per_tissue_hidden.sum(axis=0),
    )
)
cross_tissue_summary["heterogeneous_fraction"] = (
    cross_tissue_summary.heterogeneous_transcripts / cross_tissue_summary.tested_transcripts
)
write(cross_tissue_summary, "04_cross_tissue/per_tissue_summary.tsv")
np.savez_compressed(
    out / "04_cross_tissue/test_index.npz",
    tx_ids=tx_ids,
    gene_ids=gene_ids,
    tissues=np.array(tissues),
    samples=np.array(samples),
)
within_candidates = within.loc[within.divergent].copy()
within_candidates["score"] = within_candidates.difference_range * np.minimum(
    -np.log10(within_candidates.event_q.clip(lower=1e-100)), 30
)
within_candidates = within_candidates.merge(
    quality.add_suffix("_A"),
    left_on="transcript_A",
    right_on="transcript_id_A",
    how="left",
    validate="many_to_one",
).merge(
    quality.add_suffix("_B"),
    left_on="transcript_B",
    right_on="transcript_id_B",
    how="left",
    validate="many_to_one",
)
within_candidates.sort_values(
    ["strict_pair", "parent_not_detected", "opposing", "score"], ascending=False, inplace=True
)
write(within_candidates, "06_examples/within_gene_candidates.tsv.gz")
cross["score"] = cross.difference_range * np.minimum(-np.log10(cross.pair_q.clip(lower=1e-100)), 30)
cross.sort_values(
    ["strict_pair", "parent_not_detected", "opposing", "score"], ascending=False, inplace=True
)
write(cross.drop_duplicates("gene_id").head(200), "06_examples/cross_tissue_top_candidates.tsv")
classic = [
    "NEFM",
    "CALD1",
    "SLC15A1",
    "VPS29",
    "NRIP1",
    "LPL",
    "ACTA2",
    "SLC5A1",
    "POU1F1",
    "PCK1",
    "SORBS2",
    "MYO3A",
    "RYR2",
    "APOB",
    "PPARG",
    "PDE6B",
]
write(
    within_candidates.loc[within_candidates.symbol.isin(classic)],
    "06_examples/classic_within_candidates.tsv",
)
write(cross.loc[cross.symbol.isin(classic)], "06_examples/classic_cross_candidates.tsv.gz")
chosen_within = []
for symbol, tissue in [
    ("NEFM", "Medulla_oblongata"),
    ("SLC15A1", "Ileum"),
    ("VPS29", "Gallbladder"),
]:
    selected = within_candidates.loc[
        within_candidates.symbol.eq(symbol) & within_candidates.tissue.eq(tissue)
    ]
    if len(selected):
        chosen_within.append(selected.iloc[0])
if not chosen_within:
    chosen_within = [within_candidates.iloc[0]]
chosen_cross = []
for symbol, tissue_a, tissue_b in [
    ("CALD1", "Esophagus", "Lung"),
    ("LPL", "Adipose", "Gizzard"),
    ("SLC5A1", "Crop", "Duodenum"),
]:
    selected = cross.loc[
        cross.symbol.eq(symbol)
        & (
            ((cross.tissue_A == tissue_a) & (cross.tissue_B == tissue_b))
            | ((cross.tissue_A == tissue_b) & (cross.tissue_B == tissue_a))
        )
    ]
    if len(selected):
        chosen_cross.append(selected.iloc[0])
if not chosen_cross:
    chosen_cross = [cross.iloc[0]]
write(pd.DataFrame(chosen_within), "06_examples/reviewed_within_cases.tsv")
write(pd.DataFrame(chosen_cross), "06_examples/reviewed_cross_cases.tsv")
example_rows = []
for case_index, row in enumerate(chosen_within):
    tissue_index = tissues.index(row.tissue)
    features = [
        ("Parent gene", row.gene_id, gene_expression[int(gene_lookup[row.gene_id]), tissue_index]),
        (
            "Isoform A",
            row.transcript_A,
            tx_expression[int(tx_lookup[row.transcript_A]), tissue_index],
        ),
        (
            "Isoform B",
            row.transcript_B,
            tx_expression[int(tx_lookup[row.transcript_B]), tissue_index],
        ),
    ]
    for feature, identifier, values in features:
        for sample_index, value in enumerate(values):
            example_rows.append(
                dict(
                    case=f"within_{case_index+1}",
                    symbol=row.symbol,
                    tissue=row.tissue,
                    feature=feature,
                    feature_id=identifier,
                    sample=samples[sample_index],
                    age=ages[sample_index // 5],
                    bird=int(sample_order.bird_number.iloc[sample_index]),
                    log_expression=float(value),
                )
            )
for case_index, row in enumerate(chosen_cross):
    for tissue in [row.tissue_A, row.tissue_B]:
        tissue_index = tissues.index(tissue)
        features = [
            (
                "Parent gene",
                row.gene_id,
                gene_expression[int(gene_lookup[row.gene_id]), tissue_index],
            ),
            (
                "Isoform",
                row.transcript_id,
                tx_expression[int(tx_lookup[row.transcript_id]), tissue_index],
            ),
        ]
        for feature, identifier, values in features:
            for sample_index, value in enumerate(values):
                example_rows.append(
                    dict(
                        case=f"cross_{case_index+1}",
                        symbol=row.symbol,
                        tissue=tissue,
                        feature=feature,
                        feature_id=identifier,
                        sample=samples[sample_index],
                        age=ages[sample_index // 5],
                        bird=int(sample_order.bird_number.iloc[sample_index]),
                        log_expression=float(value),
                    )
                )
write(pd.DataFrame(example_rows), "06_examples/example_sample_expression.tsv")
summary = dict(
    ages=ages,
    tissues=len(tissues),
    samples=len(manifest),
    primary_model="four-age factor; no batch covariate",
    overall=overall,
    common_gene_tissue_events=len(gene_events),
    within_eligible_events=len(within),
    within_divergent_events=int(within.divergent.sum()),
    within_divergent_unique_genes=int(within.loc[within.divergent, "gene_id"].nunique()),
    within_parent_not_detected=int(within.parent_not_detected.sum()),
    within_opposing_events=int(within.opposing.sum()),
    within_strict_eligible_events=int(within.strict_eligible_isoforms.ge(2).sum()),
    within_strict_divergent_events=int(within.strict_event_q.lt(0.05).sum()),
    cross_tested_unique_transcripts=int(unique_tested.sum()),
    cross_heterogeneous_unique_transcripts=int((per_tx[:, 0] > 0).sum()),
    cross_opposing_unique_transcripts=int((per_tx[:, 1] > 0).sum()),
    cross_parent_not_detected_unique_transcripts=int((per_tx[:, 2] > 0).sum()),
    cross_strict_heterogeneous_unique_transcripts=int((per_tx[:, 3] > 0).sum()),
    cross_tested_transcript_tissue_pairs=int(
        sum(item["tested_transcripts"] for item in pair_counts)
    ),
)
check(
    "hidden_subset_of_associated",
    overall["parent_not_detected_events"]
    <= overall["age_associated_tx_events"]
    <= overall["common_tx_tissue_events"],
)
check(
    "cross_hit_subset_of_testable",
    summary["cross_heterogeneous_unique_transcripts"] <= summary["cross_tested_unique_transcripts"],
)
check(
    "within_hidden_subset_of_divergent",
    summary["within_parent_not_detected"] <= summary["within_divergent_events"] <= len(within),
)
json.dump(summary, open(out / "08_report/summary.json", "w"), indent=2)
json.dump(checks, open(out / "08_report/validation.json", "w"), indent=2)
print("COMPLETE", json.dumps(summary), flush=True)
