#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 08_Comparative_and_complex_traits
# Purpose: Functional programmes.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

from pathlib import Path
import hashlib
import json
import sys
import numpy as np
import pandas as pd
from scipy import stats

root = Path(sys.argv[1])
source = root.parent / "four_age"
ages = [20, 50, 80, 100]
reference = root / "01_reference"
checks = []


def read(path, **kwargs):
    return pd.read_csv(path, sep="\t", low_memory=False, **kwargs)


def save(frame, name):
    frame.to_csv(root / name, sep="\t", index=False)


def check(name, condition):
    checks.append(dict(check=name, passed=bool(condition)))
    assert condition, name


def bh(values):
    values = np.asarray(values, dtype=float)
    result = np.full(len(values), np.nan)
    valid = np.isfinite(values)
    result[valid] = stats.false_discovery_control(values[valid])
    return result


def signature(frame, tier, direction):
    return set(
        frame.loc[
            (frame.support_n >= tier) & frame.direction.eq(direction) & ~frame.significant_conflict,
            "gene_id",
        ]
    )


def stratified(frame, foreground, positives):
    distribution = np.array([1.0])
    expected = 0.0
    for _, subset in frame.groupby("stratum", observed=True):
        members = set(subset.index)
        total, success, draws = len(members), len(members & positives), len(members & foreground)
        expected += success * draws / total
        if success and draws:
            probabilities = stats.hypergeom.pmf(
                np.arange(min(success, draws) + 1), total, success, draws
            )
            distribution = np.convolve(distribution, probabilities)
    observed = len(foreground & positives)
    probability = float(distribution[observed:].sum()) if observed < len(distribution) else 0.0
    return min(probability, 1.0), expected


def ora(background, foreground, positives):
    foreground, positives = foreground & background, positives & background
    total, draws, successes, overlap = (
        len(background),
        len(foreground),
        len(positives),
        len(foreground & positives),
    )
    valid = 0 < draws < total and 0 < successes < total
    table = np.array(
        [[overlap, draws - overlap], [successes - overlap, total - draws - successes + overlap]]
    )
    probability = stats.fisher_exact(table, alternative="greater").pvalue if valid else np.nan
    adjusted = table.astype(float) + 0.5
    log_odds = (
        np.log(adjusted[0, 0])
        + np.log(adjusted[1, 1])
        - np.log(adjusted[0, 1])
        - np.log(adjusted[1, 0])
    )
    standard_error = np.sqrt((1 / adjusted).sum())
    return dict(
        background_genes=total,
        foreground_genes=draws,
        signature_genes=successes,
        overlap_genes=overlap,
        fold_enrichment=(overlap / draws) / (successes / total) if draws and successes else np.nan,
        p_fisher=probability,
        odds_ratio_ha=np.exp(log_odds),
        or_ci_low_ha=np.exp(log_odds - 1.96 * standard_error),
        or_ci_high_ha=np.exp(log_odds + 1.96 * standard_error),
    )


def contrast(values, older_index):
    baseline, later = values[:5], values[older_index * 5 : (older_index + 1) * 5]
    difference = later.mean() - baseline.mean()
    variance_a, variance_b = baseline.var(ddof=1) / 5, later.var(ddof=1) / 5
    standard_error = np.sqrt(variance_a + variance_b)
    degrees = (variance_a + variance_b) ** 2 / max((variance_a**2 + variance_b**2) / 4, 1e-30)
    margin = stats.t.ppf(0.975, degrees) * standard_error
    probability = stats.ttest_ind(later, baseline, equal_var=False).pvalue
    return dict(
        effect=difference,
        ci_low=difference - margin,
        ci_high=difference + margin,
        pvalue=probability,
    )


print("READ AND AUDIT", flush=True)
mammal = read(reference / "mapped_tested_genes.tsv")
mammal_all = read(reference / "mammalian_support_all7807.tsv")
check("5502_unique_mapped_genes", len(mammal) == 5502 and mammal.gene_id.is_unique)
check("7807_independent_mammal_background", len(mammal_all) == 7807)
species = ["Mouse", "Rat", "Macaque", "Human"]
for frame in [mammal, mammal_all]:
    positive = sum(
        (frame[f"{item}_P.Adjusted"] < 0.05) & (frame[f"{item}_Slope"] > 0) for item in species
    )
    negative = sum(
        (frame[f"{item}_P.Adjusted"] < 0.05) & (frame[f"{item}_Slope"] < 0) for item in species
    )
    check(
        "species_support_recalculated_" + str(len(frame)),
        np.array_equal(frame.support_n, np.maximum(positive, negative)),
    )
    check(
        "direction_conflict_recalculated_" + str(len(frame)),
        np.array_equal(frame.significant_conflict, (positive > 0) & (negative > 0)),
    )
mapped = set(mammal.gene_id)
signature_lookup = {
    (tier, direction): signature(mammal, tier, direction)
    for tier in [2, 3, 4]
    for direction in ["up", "down"]
}
manifest = read(source / "01_inputs/sample_order.tsv")
sample_names = manifest["sample"].tolist()
check(
    "four_ages_five_birds",
    manifest.age.tolist() == np.repeat(ages, 5).tolist() and manifest.individual.is_unique,
)
save(manifest, "07_source_expression/sample_order.tsv")
classification = read(source / "02_overall/transcript_classification.tsv.gz")
classification = classification[classification.gene_id.isin(mapped)].copy()
check(
    "unique_transcript_tissue_records",
    not classification.duplicated(["tissue", "transcript_id"]).any(),
)
genes = read(source / "02_overall/gene_four_age_results.tsv.gz")
genes = genes[genes.gene_id.isin(mapped)].copy()
check("unique_gene_tissue_tests", not genes.duplicated(["tissue", "gene_id"]).any())
old = read(source / "05_clusters/canonical_old_cluster_assignments.tsv.gz")
old["tissue"] = old.tissue.replace(
    {"Medulla": "Medulla_oblongata", "Olfactory": "Olfactory_bulb", "Thigh": "Thigh_muscle"}
)
check("unique_old_cluster_key", not old.duplicated(["tissue", "transcript_id"]).any())
eligible = classification[classification.common_background].merge(
    old, on=["tissue", "transcript_id"], validate="one_to_one"
)
eligible = eligible[eligible.distinct_old_clusters.eq(1) & eligible.old_cluster.notna()].copy()
eligible["old_cluster"] = eligible.old_cluster.astype(int)
tissues = sorted(classification.tissue.unique())
check("33_tissue_aliases", len(tissues) == 33)
expression = {}
gene_support = {}
expression_means = []
for tissue in tissues:
    frame = read(source / f"01_inputs/{tissue}/gene_log_expression.tsv.gz").set_index("feature_id")
    frame = frame.loc[frame.index.isin(mapped), sample_names]
    check("finite_gene_expression_" + tissue, np.isfinite(frame.to_numpy()).all())
    expression[tissue] = frame
    tissue_genes = genes[genes.tissue.eq(tissue)]
    support_columns = [f"support_{age}" for age in ages]
    valid = tissue_genes[support_columns].ge(3).all(axis=1) & tissue_genes.age_q.notna()
    gene_support[tissue] = set(tissue_genes.loc[valid, "gene_id"]) & set(frame.index)
    expression_means.extend(
        dict(tissue=tissue, gene_id=gene, mean_log_expression=value)
        for gene, value in frame.mean(axis=1).items()
    )
    save(frame.reset_index(), f"07_source_expression/{tissue}_mapped_gene_log_expression.tsv.gz")
means = pd.DataFrame(expression_means)
save(means, "02_audit/gene_expression_covariates.tsv.gz")
save(genes, "07_source_expression/mapped_gene_four_age_tests.tsv.gz")
save(classification, "07_source_expression/mapped_transcript_classification.tsv.gz")
save(
    eligible[
        ["tissue", "transcript_id", "gene_id", "old_cluster", "age_associated", "strict_quality"]
    ],
    "03_cluster/eligible_members.tsv.gz",
)
print("CLUSTER REANALYSIS", flush=True)
cluster_rows = []
stratum_records = []
for background_name, records in [
    ("original_detectable", eligible),
    ("age_associated_conditional", eligible[eligible.age_associated]),
]:
    scopes = [("all", records)] + list(records.groupby("tissue"))
    for tissue, subset in scopes:
        counts = subset.groupby("gene_id").agg(
            n_isoforms=("transcript_id", "nunique"), n_tissues=("tissue", "nunique")
        )
        mean_lookup = (
            means.groupby("gene_id").mean_log_expression.mean()
            if tissue == "all"
            else means[means.tissue.eq(tissue)].set_index("gene_id").mean_log_expression
        )
        counts["expression"] = counts.index.map(mean_lookup)
        check(
            "complete_expression_covariate_" + background_name + "_" + tissue,
            counts.expression.notna().all(),
        )
        counts["expression_bin"] = pd.qcut(counts.expression.rank(method="first"), 3, labels=False)
        counts["stratum"] = (
            pd.cut(
                counts.n_isoforms, [0, 1, 2, 4, 8, np.inf], labels=["1", "2", "3-4", "5-8", "9+"]
            ).astype(str)
            + "|"
            + counts.expression_bin.astype(str)
        )
        if tissue == "all":
            counts["stratum"] += "|" + pd.cut(
                counts.n_tissues, [0, 1, 4, 12, 33], labels=["1", "2-4", "5-12", "13-33"]
            ).astype(str)
        stratum_records.append(
            counts.reset_index().assign(background=background_name, tissue=tissue)
        )
        background = set(counts.index)
        for cluster in range(1, 11):
            foreground = set(subset.loc[subset.old_cluster.eq(cluster), "gene_id"])
            for (tier, direction), signature_genes in signature_lookup.items():
                positives = signature_genes & background
                record = ora(background, foreground, positives)
                conditional_p, expected = stratified(counts, foreground, positives)
                record.update(
                    background=background_name,
                    tissue=tissue,
                    cluster=cluster,
                    tier=tier,
                    direction=direction,
                    p_matched=conditional_p,
                    expected_matched=expected,
                    fold_matched=record["overlap_genes"] / expected if expected else np.nan,
                    hit_gene_ids=";".join(sorted(foreground & positives)),
                )
                cluster_rows.append(record)
        print("CLUSTER", background_name, tissue, flush=True)
cluster_tests = pd.DataFrame(cluster_rows)
for background_name, indices in cluster_tests.groupby("background").groups.items():
    for method in ["fisher", "matched"]:
        cluster_tests.loc[indices, "q_" + method] = bh(cluster_tests.loc[indices, "p_" + method])
cluster_tests["robust"] = (
    cluster_tests.q_fisher.lt(0.05)
    & cluster_tests.q_matched.lt(0.05)
    & cluster_tests.overlap_genes.ge(3)
)
save(cluster_tests, "03_cluster/all_cluster_signature_tests.tsv")
save(cluster_tests[cluster_tests.robust], "03_cluster/robust_cluster_signature_tests.tsv")
save(pd.concat(stratum_records), "03_cluster/matching_strata.tsv.gz")
print("GO AND INDEPENDENT MODULES", flush=True)
annotation = read(reference / "original_full_GO_BP.tsv.gz").drop_duplicates(["gene_id", "go_id"])
go_sets = annotation.groupby("go_id").gene_id.agg(set).to_dict()
go_terms = annotation.drop_duplicates("go_id").set_index("go_id").TERM.to_dict()
annotated = set(annotation.gene_id)
go_rows = []
for (tier, direction), foreground in signature_lookup.items():
    for go_id, members in go_sets.items():
        if 10 <= len(members & annotated) <= 500:
            record = ora(annotated, foreground, members)
            record.update(tier=tier, direction=direction, go_id=go_id, term=go_terms[go_id])
            go_rows.append(record)
go_results = pd.DataFrame(go_rows)
go_results["q_joint"] = bh(go_results.p_fisher)
save(
    go_results.rename(
        columns={
            "foreground_genes": "mammal_signature_genes_in_background",
            "signature_genes": "term_genes_in_background",
        }
    ),
    "04_tissue_function/independent_mammal_full_GO.tsv.gz",
)
anchors = {
    "Immune system process": "GO:0002376",
    "Cell cycle": "GO:0007049",
    "RNA processing": "GO:0006396",
    "Translation": "GO:0006412",
    "Protein transport": "GO:0015031",
    "RNA export from nucleus": "GO:0006405",
    "Proteolysis": "GO:0006508",
    "Cytoskeleton organization": "GO:0007010",
    "Oxidative phosphorylation": "GO:0006119",
    "Mitochondrion organization": "GO:0007005",
    "Muscle contraction": "GO:0006936",
    "DNA repair": "GO:0006281",
}
save(
    pd.DataFrame(
        [
            dict(module=name, go_id=go_id, annotation_genes=len(go_sets.get(go_id, set())))
            for name, go_id in anchors.items()
        ]
    ),
    "04_tissue_function/module_definitions.tsv",
)
score_rows, score_tests, score_contrasts, module_members, interaction_rows = [], [], [], [], []
generator = np.random.default_rng(20261024)
permutations = np.array([generator.permutation(20) for _ in range(4999)])
weights = np.zeros((len(permutations) + 1, 4, 20))
for permutation_index, order in enumerate(np.vstack([np.arange(20), permutations])):
    for age_index in range(4):
        weights[permutation_index, age_index, order[age_index * 5 : (age_index + 1) * 5]] = 0.2
for module, go_id in anchors.items():
    for (tier, direction), conserved in signature_lookup.items():
        members = conserved & go_sets.get(go_id, set())
        module_key = f"{module}|{tier}|{direction}"
        if not members:
            continue
        for tissue in tissues:
            usable = sorted(members & gene_support[tissue])
            module_members.extend(
                dict(
                    module=module,
                    go_id=go_id,
                    tier=tier,
                    direction=direction,
                    tissue=tissue,
                    gene_id=gene,
                )
                for gene in usable
            )
            if len(usable) < 5:
                continue
            values = expression[tissue].loc[usable].to_numpy()
            centered = values - values[:, :5].mean(axis=1, keepdims=True)
            scores = centered.mean(axis=0)
            group_means = scores.reshape(4, 5).mean(axis=1)
            for sample_index, sample in manifest.iterrows():
                score_rows.append(
                    dict(
                        module=module,
                        tier=tier,
                        direction=direction,
                        tissue=tissue,
                        n_genes=len(usable),
                        sample=sample["sample"],
                        individual=sample.individual,
                        age=sample.age,
                        score=scores[sample_index],
                    )
                )
            test = stats.f_oneway(*scores.reshape(4, 5))
            permuted_means = weights @ scores
            statistics = ((permuted_means - scores.mean()) ** 2).sum(axis=1)
            p_permutation = (1 + np.sum(statistics[1:] >= statistics[0] - 1e-12)) / len(statistics)
            score_tests.append(
                dict(
                    module=module,
                    tier=tier,
                    direction=direction,
                    tissue=tissue,
                    n_genes=len(usable),
                    p_age=test.pvalue,
                    p_permutation=p_permutation,
                    mean20=group_means[0],
                    mean50=group_means[1],
                    mean80=group_means[2],
                    mean100=group_means[3],
                    peak_age=ages[int(group_means.argmax())],
                    trough_age=ages[int(group_means.argmin())],
                )
            )
            for age_index in [1, 2, 3]:
                record = contrast(scores, age_index)
                record.update(
                    module=module,
                    tier=tier,
                    direction=direction,
                    tissue=tissue,
                    n_genes=len(usable),
                    contrast=f"{ages[age_index]}vs20",
                )
                record["aligned_with_mammal"] = (
                    record["effect"] * (1 if direction == "up" else -1) > 0
                )
                score_contrasts.append(record)
        common = sorted(set.intersection(members, *[gene_support[tissue] for tissue in tissues]))
        if len(common) >= 5:
            scores_by_tissue = []
            for tissue in tissues:
                values = expression[tissue].loc[common].to_numpy()
                scores_by_tissue.append(
                    (values - values[:, :5].mean(axis=1, keepdims=True)).mean(axis=0)
                )
            scores_by_tissue = np.asarray(scores_by_tissue).T
            residual = scores_by_tissue - scores_by_tissue.mean(axis=1, keepdims=True)
            residual -= residual.mean(axis=0, keepdims=True)
            permuted = weights @ residual
            statistics = (permuted**2).sum(axis=(1, 2))
            probability = (1 + np.sum(statistics[1:] >= statistics[0] - 1e-12)) / len(statistics)
            interaction_rows.append(
                dict(
                    module=module,
                    tier=tier,
                    direction=direction,
                    common_genes=len(common),
                    tissues=33,
                    p_interaction=probability,
                    common_gene_ids=";".join(common),
                    status="tested",
                )
            )
        else:
            interaction_rows.append(
                dict(
                    module=module,
                    tier=tier,
                    direction=direction,
                    common_genes=len(common),
                    tissues=33,
                    p_interaction=np.nan,
                    common_gene_ids=";".join(common),
                    status="fewer_than_5_common_genes",
                )
            )
        print("MODULE", module_key, len(members), len(common), flush=True)
scores = pd.DataFrame(score_rows)
module_tests = pd.DataFrame(score_tests)
module_tests["q_age_joint"] = bh(module_tests.p_age)
module_tests["q_permutation_joint"] = bh(module_tests.p_permutation)
contrasts = pd.DataFrame(score_contrasts)
contrasts["q_joint"] = bh(contrasts.pvalue)
interactions = pd.DataFrame(interaction_rows)
interactions["q_interaction_joint"] = bh(interactions.p_interaction)
members_frame = pd.DataFrame(module_members)
save(scores, "04_tissue_function/module_individual_scores.tsv.gz")
save(module_tests, "04_tissue_function/module_age_tests.tsv")
save(contrasts, "04_tissue_function/module_vs20_contrasts.tsv")
save(interactions, "04_tissue_function/module_age_by_tissue_interactions.tsv")
save(members_frame, "04_tissue_function/module_gene_members.tsv.gz")
gene_cluster = eligible[eligible.age_associated][
    ["tissue", "gene_id", "old_cluster"]
].drop_duplicates()
carriers = (
    members_frame.merge(gene_cluster, on=["tissue", "gene_id"], validate="many_to_many")
    .groupby(["module", "tier", "direction", "tissue", "old_cluster"])
    .gene_id.nunique()
    .rename("unique_genes")
    .reset_index()
)
save(carriers, "04_tissue_function/module_cluster_carriers.tsv")
print("ONT INCREMENT", flush=True)
increment_rows, increment_events = [], []
for quality_name, selected in [
    (
        "common_detectable",
        classification[classification.common_background & classification.age_associated],
    ),
    (
        "strict_quality",
        classification[
            classification.common_background
            & classification.age_associated
            & classification.strict_quality
        ],
    ),
]:
    events = (
        selected.groupby(["tissue", "gene_id"])
        .agg(n_age_isoforms=("transcript_id", "nunique"))
        .reset_index()
    )
    events = events.merge(
        genes[["tissue", "gene_id", "age_q", "max_abs_lfc"]],
        on=["tissue", "gene_id"],
        validate="one_to_one",
    )
    events = events[events.age_q.notna()].copy()
    events["parent_category"] = np.select(
        [events.age_q.ge(0.05), events.max_abs_lfc.lt(1)],
        ["parent_q_ge_0.05", "parent_significant_small_effect"],
        default="parent_age_associated",
    )
    events = events.merge(
        mammal[["gene_id", "support_n", "direction", "significant_conflict"]],
        on="gene_id",
        validate="many_to_one",
    )
    increment_events.append(events.assign(quality=quality_name))
    for (tier, direction), conserved in signature_lookup.items():
        foreground = events[events.gene_id.isin(conserved)]
        for tissue, subset in [("all_gene_tissue_events", foreground)] + list(
            foreground.groupby("tissue")
        ):
            counts = subset.parent_category.value_counts()
            increment_rows.append(
                dict(
                    quality=quality_name,
                    tier=tier,
                    direction=direction,
                    tissue=tissue,
                    denominator_events=len(subset),
                    unique_genes=subset.gene_id.nunique(),
                    parent_nonsignificant=int(counts.get("parent_q_ge_0.05", 0)),
                    parent_significant_small_effect=int(
                        counts.get("parent_significant_small_effect", 0)
                    ),
                    parent_age_associated=int(counts.get("parent_age_associated", 0)),
                    fraction_parent_nonsignificant=(
                        counts.get("parent_q_ge_0.05", 0) / len(subset) if len(subset) else np.nan
                    ),
                    unique_genes_with_any_parent_nonsignificant=subset.loc[
                        subset.parent_category.eq("parent_q_ge_0.05"), "gene_id"
                    ].nunique(),
                )
            )
save(pd.DataFrame(increment_rows), "05_ont_increment/proportions_by_tissue_quality_tier.tsv")
save(pd.concat(increment_events), "05_ont_increment/eligible_gene_tissue_events.tsv.gz")
print("CANDIDATES", flush=True)
within = read(source / "03_within_gene/gene_tissue_tests.tsv.gz")
within = within[
    within.gene_id.isin(signature(mammal, 2, "up") | signature(mammal, 2, "down"))
    & within.strict_event_q.lt(0.05)
].copy()
within = within.merge(
    mammal[["gene_id", "support_n", "direction", "support_species", "spans_rodents_primates"]],
    on="gene_id",
    validate="many_to_one",
)
quality = classification[
    ["tissue", "transcript_id", "class", "age_q", "strict_quality", "parent_age_q"]
].merge(
    old[["tissue", "transcript_id", "old_cluster"]],
    on=["tissue", "transcript_id"],
    how="left",
    validate="one_to_one",
)
for side in ["A", "B"]:
    details = quality.rename(
        columns={column: column + "_" + side for column in quality.columns if column != "tissue"}
    )
    within = within.merge(
        details,
        left_on=["tissue", f"strict_transcript_{side}"],
        right_on=["tissue", f"transcript_id_{side}"],
        how="left",
        validate="many_to_one",
    )
regular = {"FSM", "ISM", "NIC", "NNC"}
within["eligible_case"] = (
    within.strict_quality_A.eq(True)
    & within.strict_quality_B.eq(True)
    & within.class_A.isin(regular)
    & within.class_B.isin(regular)
)
within["different_original_clusters"] = (
    within.old_cluster_A.notna()
    & within.old_cluster_B.notna()
    & within.old_cluster_A.ne(within.old_cluster_B)
)
within["parent_age_not_significant"] = within.parent_age_q_A.ge(0.05)
network = read(reference / "C7_original_network_exact_connections.tsv.gz")
network_keys = set(zip(network.tissue, network.gene_id))
within["original_C7_network_gene_tissue"] = [
    (tissue, gene) in network_keys for tissue, gene in zip(within.tissue, within.gene_id)
]
within["rank_score"] = (
    within.support_n * 10
    + within.spans_rodents_primates.astype(int) * 5
    + within.parent_age_not_significant.astype(int) * 10
    + within.different_original_clusters.astype(int) * 3
    + within.original_C7_network_gene_tissue.astype(int) * 2
)
within = within.sort_values(
    ["eligible_case", "rank_score", "strict_event_q"], ascending=[False, False, True]
)
save(within, "06_candidates/all_within_gene_candidates.tsv")
shortlist = within[within.eligible_case].drop_duplicates("gene_id").head(20)
named = within[
    within.symbol.isin(["EIF2S1", "MTFR1", "THOC7", "PCID2", "OPTN"]) & within.eligible_case
]
shortlist = pd.concat([shortlist, named]).drop_duplicates(["tissue", "gene_id"])
save(shortlist, "06_candidates/within_gene_shortlist.tsv")
cross = read(reference / "all_strict_opposing_candidates.tsv.gz")
cross = cross.merge(
    mammal[
        [
            "gene_id",
            "support_n",
            "direction",
            "significant_conflict",
            "support_species",
            "spans_rodents_primates",
        ]
    ],
    on="gene_id",
    validate="many_to_one",
)
cross = cross[cross.support_n.ge(2) & ~cross.significant_conflict].copy()
cross["rank_score"] = (
    cross.support_n * 10
    + cross.spans_rodents_primates.astype(int) * 5
    + cross.parent_category.eq("neither_tissue_detected").astype(int) * 10
)
cross = cross.sort_values(["rank_score", "pair_q"], ascending=[False, True])
save(cross, "06_candidates/conserved_cross_tissue_candidates.tsv.gz")
save(cross.drop_duplicates("gene_id").head(20), "06_candidates/cross_tissue_shortlist.tsv")
requests = []
for row in shortlist.itertuples():
    requests.extend(
        [
            dict(tissue=row.tissue, feature="gene", feature_id=row.gene_id, symbol=row.symbol),
            dict(
                tissue=row.tissue,
                feature="tx",
                feature_id=row.strict_transcript_A,
                symbol=row.symbol,
            ),
            dict(
                tissue=row.tissue,
                feature="tx",
                feature_id=row.strict_transcript_B,
                symbol=row.symbol,
            ),
        ]
    )
for row in cross.drop_duplicates("gene_id").head(10).itertuples():
    for tissue in [row.tissue_A, row.tissue_B]:
        requests.extend(
            [
                dict(tissue=tissue, feature="gene", feature_id=row.gene_id, symbol=row.symbol),
                dict(tissue=tissue, feature="tx", feature_id=row.transcript_id, symbol=row.symbol),
            ]
        )
requests = pd.DataFrame(requests).drop_duplicates()
candidate_expression = []
for (tissue, feature), group in requests.groupby(["tissue", "feature"]):
    frame = read(source / f"01_inputs/{tissue}/{feature}_log_expression.tsv.gz")
    frame = frame.merge(group, on="feature_id", validate="one_to_one")
    frame = frame.melt(
        id_vars=["tissue", "feature", "feature_id", "symbol"],
        value_vars=sample_names,
        var_name="sample",
        value_name="log_expression",
    ).merge(manifest, on="sample", validate="many_to_one")
    frame["relative20"] = frame.log_expression - frame.groupby(
        ["feature_id", "tissue"]
    ).log_expression.transform(lambda values: values[frame.loc[values.index, "age"].eq(20)].mean())
    candidate_expression.append(frame)
save(pd.concat(candidate_expression), "06_candidates/candidate_individual_expression.tsv.gz")
legacy = read(reference / "all_gene_tissue_evidence.tsv.gz")
legacy_rows = []
for keys, members in members_frame.groupby(["module", "tier", "direction", "tissue"]):
    module, tier, direction, tissue = keys
    evidence = legacy[legacy.tissue.eq(tissue) & legacy.gene_id.isin(members.gene_id)]
    legacy_rows.append(
        dict(
            module=module,
            tier=tier,
            direction=direction,
            tissue=tissue,
            module_genes=len(members),
            legacy_testable=len(evidence),
            pass20=int(evidence.pass_20.sum()),
            pass50=int(evidence.pass_50.sum()),
            pass_both=int((evidence.pass_20 & evidence.pass_50).sum()),
            opposite20=int(evidence.opposite_20.sum()),
            opposite50=int(evidence.opposite_50.sum()),
        )
    )
save(pd.DataFrame(legacy_rows), "04_tissue_function/legacy_reference_flags.tsv")
save(pd.DataFrame(checks), "02_audit/input_validation.tsv")
summary = dict(
    mapped_genes=len(mapped),
    eligible_cluster_gene_tissue=len(eligible[["tissue", "gene_id"]].drop_duplicates()),
    cluster_tests=len(cluster_tests),
    robust_cluster_tests=int(cluster_tests.robust.sum()),
    module_age_tests=len(module_tests),
    module_age_joint_significant=int(module_tests.q_permutation_joint.lt(0.05).sum()),
    interaction_tests=int(interactions.p_interaction.notna().sum()),
    significant_interactions=int(interactions.q_interaction_joint.lt(0.05).sum()),
    candidate_gene_tissues=int(within.eligible_case.sum()),
    candidate_unique_genes=within.loc[within.eligible_case, "gene_id"].nunique(),
    cross_tissue_candidates=len(cross),
    mammal_tissue_specific_reference_available=False,
    no_figures=True,
    permutations=4999,
)
(root / "08_report/summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2), flush=True)
