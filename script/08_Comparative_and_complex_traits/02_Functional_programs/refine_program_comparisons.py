#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 08_Comparative_and_complex_traits
# Purpose: Refine cross-species functional-program comparisons.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

from pathlib import Path
import sys
import json
import numpy as np
import pandas as pd
from scipy import stats

root = Path(sys.argv[1])
ages = [20, 50, 80, 100]


def read(name):
    return pd.read_csv(root / name, sep="\t", low_memory=False)


def save(frame, name):
    frame.to_csv(root / name, sep="\t", index=False)


def bh(values):
    values = np.asarray(values, dtype=float)
    result = np.full(values.shape, np.nan)
    valid = np.isfinite(values)
    if valid.any():
        result[valid] = stats.false_discovery_control(values[valid])
    return result


manifest = read("07_source_expression/sample_order.tsv")
samples = manifest["sample"].tolist()
members = read("04_tissue_function/module_gene_members.tsv.gz")
mammal = read("01_reference/mapped_tested_genes.tsv")
genes = read("07_source_expression/mapped_gene_four_age_tests.tsv.gz")
module_tests = read("04_tissue_function/module_age_tests.tsv")
signature = mammal[mammal.support_n.ge(2) & ~mammal.significant_conflict]
annotations = read("01_reference/original_full_GO_BP.tsv.gz")
definitions = read("04_tissue_function/module_definitions.tsv")
annotation_modules = annotations.merge(definitions[["module", "go_id"]], on="go_id")[
    ["gene_id", "module"]
].drop_duplicates()
generator = np.random.default_rng(20261025)
weights = np.zeros((5000, 4, 20))
for permutation_index in range(5000):
    order = np.arange(20) if permutation_index == 0 else generator.permutation(20)
    for age_index in range(4):
        weights[permutation_index, age_index, order[age_index * 5 : (age_index + 1) * 5]] = 0.2
competitive, competitive_scores, rank_results, gene_coherence = [], [], [], []
for tissue, tissue_members in members.groupby("tissue"):
    expression = read(f"07_source_expression/{tissue}_mapped_gene_log_expression.tsv.gz").set_index(
        "feature_id"
    )[samples]
    tissue_genes = genes[genes.tissue.eq(tissue)].set_index("gene_id")
    supported = (
        tissue_genes[[f"support_{age}" for age in ages]].ge(3).all(axis=1)
        & tissue_genes.age_q.notna()
    )
    valid = sorted(set(tissue_genes.index[supported]) & set(expression.index))
    expression = expression.loc[valid]
    centered = expression.sub(expression.iloc[:, :5].mean(axis=1), axis=0)
    bins = pd.qcut(expression.mean(axis=1).rank(method="first"), 5, labels=False)
    gene_effects = centered.iloc[:, 15:].mean(axis=1)
    for (module, tier, direction), group in tissue_members.groupby(["module", "tier", "direction"]):
        foreground = sorted(set(group.gene_id) & set(valid))
        if len(foreground) < 5:
            continue
        conserved = set(
            signature.loc[
                signature.support_n.ge(tier) & signature.direction.eq(direction), "gene_id"
            ]
        ) & set(valid)
        background_controls = conserved - set(foreground)
        matched_weights = pd.Series(0.0, index=valid)
        usable = True
        for bin_index, count in bins.loc[foreground].value_counts().items():
            controls = sorted(background_controls & set(bins.index[bins.eq(bin_index)]))
            if len(controls) < 3:
                usable = False
                break
            matched_weights.loc[controls] = count / len(foreground) / len(controls)
        trajectory = centered.loc[foreground].mean(axis=0).to_numpy()
        changes = centered.loc[foreground].to_numpy().reshape(len(foreground), 4, 5).mean(axis=2)
        for age_index in [1, 2, 3]:
            effects = changes[:, age_index]
            aligned = effects * (1 if direction == "up" else -1)
            leave_one_out = (effects.sum() - effects) / (len(effects) - 1)
            gene_coherence.append(
                dict(
                    tissue=tissue,
                    module=module,
                    tier=tier,
                    direction=direction,
                    age=ages[age_index],
                    n_genes=len(foreground),
                    same_direction_genes=int((aligned > 0).sum()),
                    fraction_same_direction=float((aligned > 0).mean()),
                    mean_effect=effects.mean(),
                    median_gene_effect=np.median(effects),
                    leave_one_gene_out_min=leave_one_out.min(),
                    leave_one_gene_out_max=leave_one_out.max(),
                )
            )
        if usable:
            comparator = centered.mul(matched_weights, axis=0).sum(axis=0).to_numpy()
            differential = trajectory - comparator
            perm_means = weights @ differential
            statistics = ((perm_means - differential.mean()) ** 2).sum(axis=1)
            probability = (1 + np.sum(statistics[1:] >= statistics[0] - 1e-12)) / len(statistics)
            means = differential.reshape(4, 5).mean(axis=1)
            competitive.append(
                dict(
                    tissue=tissue,
                    module=module,
                    tier=tier,
                    direction=direction,
                    n_genes=len(foreground),
                    n_controls=int(matched_weights.gt(0).sum()),
                    p_competitive_age=probability,
                    difference20=means[0],
                    difference50=means[1],
                    difference80=means[2],
                    difference100=means[3],
                    status="tested",
                )
            )
            for sample_index, sample in manifest.iterrows():
                competitive_scores.append(
                    dict(
                        tissue=tissue,
                        module=module,
                        tier=tier,
                        direction=direction,
                        sample=sample["sample"],
                        age=sample.age,
                        module_score=trajectory[sample_index],
                        matched_reference_score=comparator[sample_index],
                        difference=differential[sample_index],
                    )
                )
        else:
            competitive.append(
                dict(
                    tissue=tissue,
                    module=module,
                    tier=tier,
                    direction=direction,
                    n_genes=len(foreground),
                    n_controls=len(background_controls),
                    p_competitive_age=np.nan,
                    status="insufficient_same_direction_matched_controls",
                )
            )
    for tier in [2, 3, 4]:
        for direction in ["up", "down"]:
            foreground = set(
                signature.loc[
                    signature.support_n.ge(tier) & signature.direction.eq(direction), "gene_id"
                ]
            ) & set(valid)
            control = set(valid) - foreground
            if len(foreground) < 5 or len(control) < 5:
                continue
            age_stats = tissue_genes.loc[valid]
            signed_statistic = age_stats["100vs20_log2FoldChange"] / age_stats["100vs20_lfcSE"]
            finite = set(signed_statistic.index[np.isfinite(signed_statistic)])
            first = signed_statistic.loc[sorted(foreground & finite)].to_numpy()
            second = signed_statistic.loc[sorted(control & finite)].to_numpy()
            result = stats.mannwhitneyu(
                first, second, alternative="greater" if direction == "up" else "less"
            )
            rank_results.append(
                dict(
                    tissue=tissue,
                    tier=tier,
                    direction=direction,
                    n_signature=len(first),
                    n_background=len(second),
                    median_signed_wald=np.median(first),
                    rank_biserial=2 * result.statistic / (len(first) * len(second)) - 1,
                    p_rank=result.pvalue,
                    interpretation="competitive gene-level endpoint ranking; genes correlated; supporting only",
                )
            )
    print("REFINED", tissue, flush=True)
competitive = pd.DataFrame(competitive)
competitive["q_competitive_age_joint"] = bh(competitive.p_competitive_age)
save(competitive, "04_tissue_function/module_expression_matched_competitive_tests.tsv")
save(pd.DataFrame(competitive_scores), "04_tissue_function/module_matched_individual_scores.tsv.gz")
save(
    pd.DataFrame(gene_coherence), "04_tissue_function/module_member_direction_and_leave_one_out.tsv"
)
rank_results = pd.DataFrame(rank_results)
rank_results["q_rank_joint"] = bh(rank_results.p_rank)
save(rank_results, "04_tissue_function/full_background_directional_rank_support.tsv")
joined = module_tests.merge(
    competitive, on=["tissue", "module", "tier", "direction", "n_genes"], validate="one_to_one"
)
joined["dual_supported"] = joined.q_permutation_joint.lt(0.05) & joined.q_competitive_age_joint.lt(
    0.05
)
save(joined, "04_tissue_function/module_joint_evidence.tsv")
events = read("05_ont_increment/eligible_gene_tissue_events.tsv.gz")
module_increment = events.merge(annotation_modules, on="gene_id", validate="many_to_many")
module_increment = module_increment[
    module_increment.support_n.ge(2) & ~module_increment.significant_conflict
]
increment_rows = []
for tier in [2, 3, 4]:
    for keys, group in module_increment[module_increment.support_n.ge(tier)].groupby(
        ["quality", "module", "direction", "tissue"]
    ):
        quality, module, direction, tissue = keys
        increment_rows.append(
            dict(
                quality=quality,
                module=module,
                tier=tier,
                direction=direction,
                tissue=tissue,
                denominator_genes=group.gene_id.nunique(),
                parent_not_significant=group.loc[
                    group.parent_category.eq("parent_q_ge_0.05"), "gene_id"
                ].nunique(),
            )
        )
increment_frame = pd.DataFrame(increment_rows)
increment_frame["fraction_parent_not_significant"] = (
    increment_frame.parent_not_significant / increment_frame.denominator_genes
)
save(increment_frame, "05_ont_increment/module_tissue_increment.tsv")
within = read("06_candidates/all_within_gene_candidates.tsv")
module_names = annotation_modules.groupby("gene_id").module.agg(
    lambda values: ";".join(sorted(set(values)))
)
within["functional_modules"] = within.gene_id.map(module_names)
within["both_strict_isoforms_age_significant"] = within.age_q_A.lt(0.05) & within.age_q_B.lt(0.05)
save(within, "06_candidates/all_within_gene_candidates.tsv")
cross = read("06_candidates/conserved_cross_tissue_candidates.tsv.gz")
cross["rank_score"] = (
    cross.support_n * 10
    + cross.spans_rodents_primates.astype(int) * 5
    + cross.parent_category.eq("neither_tissue_detected").astype(int) * 10
)
cross["functional_modules"] = cross.gene_id.map(module_names)
cross = cross.sort_values(["rank_score", "pair_q"], ascending=[False, True])
save(cross, "06_candidates/conserved_cross_tissue_candidates.tsv.gz")
save(cross.drop_duplicates("gene_id").head(20), "06_candidates/cross_tissue_shortlist.tsv")
new_requests = []
for row in cross.drop_duplicates("gene_id").head(20).itertuples():
    for tissue in [row.tissue_A, row.tissue_B]:
        new_requests.extend(
            [
                dict(tissue=tissue, feature="gene", feature_id=row.gene_id, symbol=row.symbol),
                dict(tissue=tissue, feature="tx", feature_id=row.transcript_id, symbol=row.symbol),
            ]
        )
source = root.parent / "four_age"
expression_rows = []
for (tissue, feature), group in (
    pd.DataFrame(new_requests).drop_duplicates().groupby(["tissue", "feature"])
):
    frame = pd.read_csv(source / f"01_inputs/{tissue}/{feature}_log_expression.tsv.gz", sep="\t")
    frame = frame.merge(group, on="feature_id", validate="one_to_one")
    frame = frame.melt(
        id_vars=["tissue", "feature", "feature_id", "symbol"],
        value_vars=samples,
        var_name="sample",
        value_name="log_expression",
    ).merge(manifest, on="sample", validate="many_to_one")
    baseline = (
        frame[frame.age.eq(20)]
        .groupby(["tissue", "feature_id"])
        .log_expression.mean()
        .rename("baseline")
    )
    frame = frame.merge(baseline, on=["tissue", "feature_id"], validate="many_to_one")
    frame["relative20"] = frame.log_expression - frame.baseline
    expression_rows.append(frame.drop(columns="baseline"))
old_expression = read("06_candidates/candidate_individual_expression.tsv.gz")
expression = pd.concat([old_expression] + expression_rows).drop_duplicates(
    ["tissue", "feature", "feature_id", "sample"]
)
save(expression, "06_candidates/candidate_individual_expression.tsv.gz")
summary = json.loads((root / "08_report/summary.json").read_text())
summary.update(
    competitive_modules_tested=int(competitive.p_competitive_age.notna().sum()),
    dual_supported_module_tissue_tests=int(joined.dual_supported.sum()),
)
(root / "08_report/summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2), flush=True)
