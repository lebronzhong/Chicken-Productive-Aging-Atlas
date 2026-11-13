#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 08_Comparative_and_complex_traits
# Purpose: Tissue pairs.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

from pathlib import Path
import itertools
import json
import sys
import numpy as np
import pandas as pd
from scipy import stats

root = Path(sys.argv[1])


def read(name):
    return pd.read_csv(root / name, sep="\t", low_memory=False)


manifest = read("07_source_expression/sample_order.tsv")
interactions = read("04_tissue_function/module_age_by_tissue_interactions.tsv")
tissues = sorted(read("04_tissue_function/module_age_tests.tsv").tissue.unique())
generator = np.random.default_rng(20261026)
weights = np.zeros((5000, 4, 20))
for permutation_index in range(5000):
    order = np.arange(20) if permutation_index == 0 else generator.permutation(20)
    for age_index in range(4):
        weights[permutation_index, age_index, order[age_index * 5 : (age_index + 1) * 5]] = 0.2
records, source_rows = [], []
for module, direction in [
    ("Immune system process", "up"),
    ("Cell cycle", "down"),
    ("RNA processing", "down"),
]:
    specification = interactions[
        interactions.module.eq(module)
        & interactions.tier.eq(3)
        & interactions.direction.eq(direction)
    ].iloc[0]
    if specification.status != "tested":
        continue
    gene_ids = specification.common_gene_ids.split(";")
    tissue_scores = []
    for tissue in tissues:
        frame = read(f"07_source_expression/{tissue}_mapped_gene_log_expression.tsv.gz").set_index(
            "feature_id"
        )
        values = frame.loc[gene_ids, manifest["sample"]].to_numpy()
        score = (values - values[:, :5].mean(axis=1, keepdims=True)).mean(axis=0)
        tissue_scores.append(score)
        source_rows.extend(
            dict(
                module=module,
                tier=3,
                direction=direction,
                common_genes=len(gene_ids),
                tissue=tissue,
                sample=sample["sample"],
                age=sample.age,
                score=score[index],
            )
            for index, sample in manifest.iterrows()
        )
    tissue_scores = np.asarray(tissue_scores)
    indices = list(itertools.combinations(range(len(tissues)), 2))
    differences = np.array(
        [tissue_scores[first] - tissue_scores[second] for first, second in indices]
    ).T
    centered = differences - differences.mean(axis=0, keepdims=True)
    permuted_means = weights @ centered
    statistics = (permuted_means**2).sum(axis=1)
    probabilities = (1 + np.sum(statistics[1:] >= statistics[0] - 1e-12, axis=0)) / len(statistics)
    for pair_index, (first, second) in enumerate(indices):
        values = differences[:, pair_index]
        early, late = values[:5], values[15:]
        variance_early, variance_late = early.var(ddof=1) / 5, late.var(ddof=1) / 5
        degrees = (variance_early + variance_late) ** 2 / (
            (variance_early**2 + variance_late**2) / 4
        )
        effect = late.mean() - early.mean()
        margin = stats.t.ppf(0.975, degrees) * np.sqrt(variance_early + variance_late)
        means = values.reshape(4, 5).mean(axis=1)
        records.append(
            dict(
                module=module,
                tier=3,
                direction=direction,
                tissue_A=tissues[first],
                tissue_B=tissues[second],
                common_genes=len(gene_ids),
                p_interaction=probabilities[pair_index],
                difference50=means[1],
                difference80=means[2],
                difference100=means[3],
                effect100vs20=effect,
                ci_low=effect - margin,
                ci_high=effect + margin,
                p_100vs20=stats.ttest_ind(late, early, equal_var=False).pvalue,
            )
        )
result = pd.DataFrame(records)
result["q_interaction_joint"] = stats.false_discovery_control(result.p_interaction.to_numpy())
result["q_100vs20_joint"] = stats.false_discovery_control(result.p_100vs20.to_numpy())
result.to_csv(
    root / "04_tissue_function/core_modules_all_tissue_pair_tests.tsv", sep="\t", index=False
)
pd.DataFrame(source_rows).to_csv(
    root / "04_tissue_function/core_modules_common_gene_individual_scores.tsv.gz",
    sep="\t",
    index=False,
)
assert len(result) == 3 * 33 * 32 // 2
summary = json.loads((root / "08_report/summary.json").read_text())
summary.update(
    core_module_pair_tests=len(result),
    core_module_pair_interactions_significant=int(result.q_interaction_joint.lt(0.05).sum()),
)
(root / "08_report/summary.json").write_text(json.dumps(summary, indent=2))
print("All 1584 pairwise tests completed; within-bird pairing and identical gene sets retained.")
