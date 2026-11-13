#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 02_Expression_dynamics
# Purpose: Age union analysis.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd


def tsi(values):
    logged = np.log1p(values)
    total = logged.sum(axis=-1, keepdims=True)
    return np.divide(logged, total, out=np.zeros_like(logged), where=total > 0)


def ratio(numerator, denominator):
    return float(numerator / denominator) if denominator else None


def to_genes(transcripts, parents, number):
    result = np.zeros((number,) + transcripts.shape[1:], dtype=bool)
    np.logical_or.at(result, parents, transcripts)
    return result


def categories(events, eligible):
    active = events.any(axis=2)
    number_ages = active.sum(axis=1)
    number_tissues = events.any(axis=1).sum(axis=1)
    result = np.full(len(events), "no_event", dtype=object)
    result[~eligible.any(axis=1)] = "not_detectable"
    result[(number_ages == 1) & eligible.all(axis=1)] = "single_age_all_four_assessable"
    result[(number_ages == 1) & ~eligible.all(axis=1)] = "single_age_incomplete_background"
    result[(number_ages >= 2) & (number_ages < 4) & (number_tissues == 1)] = (
        "multiple_ages_same_tissue"
    )
    result[(number_ages == 4) & (number_tissues == 1)] = "all_four_same_tissue"
    result[number_tissues >= 2] = "specific_tissue_changes"
    return result


def summarize(detect, specific, core, gene_specific, parents, gene_count, label, ages, tissues):
    eligible = detect.any(axis=2)
    gene_eligible = to_genes(eligible, parents, gene_count)
    gene_event = to_genes(specific.any(axis=2), parents, gene_count)
    gene_core = to_genes(core.any(axis=2), parents, gene_count)
    summaries = []
    tissue_rows = []
    for age_index in list(range(4)) + [None]:
        age_label = ages[age_index] if age_index is not None else "union"
        tx_background = eligible[:, age_index] if age_index is not None else eligible.any(axis=1)
        gene_background = (
            gene_eligible[:, age_index] if age_index is not None else gene_eligible.any(axis=1)
        )
        tx_positive = (
            specific[:, age_index, :].any(axis=1)
            if age_index is not None
            else specific.any(axis=(1, 2))
        )
        tx_core = (
            core[:, age_index, :].any(axis=1) if age_index is not None else core.any(axis=(1, 2))
        )
        gene_positive = (
            gene_event[:, age_index] if age_index is not None else gene_event.any(axis=1)
        )
        gene_core_positive = (
            gene_core[:, age_index] if age_index is not None else gene_core.any(axis=1)
        )
        gene_direct = (
            (gene_specific & gene_eligible)[:, age_index]
            if age_index is not None
            else (gene_specific & gene_eligible).any(axis=1)
        )
        gene_nonspecific = (
            (~gene_specific & gene_eligible)[:, age_index]
            if age_index is not None
            else (~gene_specific & gene_eligible).any(axis=1)
        )
        never_specific = gene_background & ~gene_direct
        strict_core = gene_core_positive & never_specific
        summary = dict(
            scenario=label,
            age=age_label,
            eligible_transcripts=int(tx_background.sum()),
            eligible_genes=int(gene_background.sum()),
            specific_transcripts=int(tx_positive.sum()),
            core_transcripts=int(tx_core.sum()),
            transcript_specific_parent_genes=int(gene_positive.sum()),
            gene_specific_parent_genes=int(gene_direct.sum()),
            core_parent_genes=int(gene_core_positive.sum()),
            nonspecific_in_any_comparable_age_genes=int(gene_nonspecific.sum()),
            never_specific_in_comparable_ages_genes=int(never_specific.sum()),
            strict_never_specific_core_genes=int(strict_core.sum()),
            core_genes_also_gene_specific_other_age=int((gene_core_positive & gene_direct).sum()),
            gene_only_union=int((gene_direct & ~gene_positive).sum()),
            both_levels_union=int((gene_direct & gene_positive).sum()),
            transcript_only_union=int((gene_positive & ~gene_direct).sum()),
            either_level_union=int((gene_positive | gene_direct).sum()),
            core_tx_over_eligible_tx=ratio(tx_core.sum(), tx_background.sum()),
            core_tx_over_specific_tx=ratio(tx_core.sum(), tx_positive.sum()),
            core_gene_over_eligible_gene=ratio(gene_core_positive.sum(), gene_background.sum()),
            core_gene_over_any_nonspecific_gene=ratio(
                gene_core_positive.sum(), gene_nonspecific.sum()
            ),
            strict_core_gene_over_never_specific_gene=ratio(
                strict_core.sum(), never_specific.sum()
            ),
        )
        assert (
            summary["core_parent_genes"]
            == summary["strict_never_specific_core_genes"]
            + summary["core_genes_also_gene_specific_other_age"]
        )
        assert (
            summary["either_level_union"]
            == summary["gene_only_union"]
            + summary["both_levels_union"]
            + summary["transcript_only_union"]
        )
        assert (
            summary["core_transcripts"]
            <= summary["specific_transcripts"]
            <= summary["eligible_transcripts"]
        )
        summaries.append(summary)
        for tissue_index, tissue in enumerate(tissues):
            background = (
                detect[:, age_index, tissue_index]
                if age_index is not None
                else detect[:, :, tissue_index].any(axis=1)
            )
            positive = (
                specific[:, age_index, tissue_index]
                if age_index is not None
                else specific[:, :, tissue_index].any(axis=1)
            )
            hidden = (
                core[:, age_index, tissue_index]
                if age_index is not None
                else core[:, :, tissue_index].any(axis=1)
            )
            gene_bg = np.unique(parents[background])
            gene_hidden = np.unique(parents[hidden])
            gene_positive_target = np.unique(parents[positive])
            tissue_nonspecific = ~gene_specific[parents] & detect[:, :, tissue_index]
            nonspecific_target = (
                tissue_nonspecific[:, age_index]
                if age_index is not None
                else tissue_nonspecific.any(axis=1)
            )
            nonspecific_parents = np.unique(parents[nonspecific_target])
            assert hidden.sum() <= positive.sum() <= background.sum()
            tissue_rows.append(
                dict(
                    scenario=label,
                    age=age_label,
                    tissue=tissue,
                    eligible_transcripts=int(background.sum()),
                    eligible_genes=len(gene_bg),
                    specific_transcripts=int(positive.sum()),
                    specific_parent_genes=len(gene_positive_target),
                    core_transcripts=int(hidden.sum()),
                    core_parent_genes=len(gene_hidden),
                    gene_nonspecific_in_any_eligible_age=len(nonspecific_parents),
                    core_tx_over_specific_tx=ratio(hidden.sum(), positive.sum()),
                    core_tx_over_eligible_tx=ratio(hidden.sum(), background.sum()),
                    core_gene_over_eligible_gene=ratio(len(gene_hidden), len(gene_bg)),
                    core_gene_over_any_nonspecific_gene=ratio(
                        len(gene_hidden), len(nonspecific_parents)
                    ),
                    global_eligible_transcripts=int(tx_background.sum()),
                    global_eligible_genes=int(gene_background.sum()),
                )
            )
    return summaries, tissue_rows, gene_eligible, gene_event, gene_core


def main():
    root = Path("/path/to/chicken-aging/analysis")
    source = root / "specificity/results"
    output = root / "Fig2_age_union_20261016/results"
    output.mkdir(parents=True, exist_ok=True)
    cache = np.load(source / "expression_background.npz")
    mapping = pd.read_csv(source / "canonical_transcript_mapping.tsv.gz", sep="\t")
    old = pd.read_csv(source / "primary_transcript_classification.tsv.gz", sep="\t")
    assert list(mapping.transcript_id) == list(old.transcript_id)
    parents = cache["parent"]
    gene_ids = cache["gene_ids"]
    ages = cache["ages"].tolist()
    tissues = cache["tissues"].tolist()
    assert np.array_equal(gene_ids[parents], mapping.gene_id.to_numpy())
    gene_present = cache["gene_present"]
    tx_expression = cache["tx_median"]
    gene_expression = cache["gene_median"]
    tx_tsi = tsi(tx_expression)
    gene_tsi = tsi(gene_expression)
    tx_support = cache["tx_support"]
    gene_support = cache["gene_support"]
    all_summaries = []
    all_tissues = []
    specifications = [
        ("primary", 0.8, 1, 3, "none"),
        ("fixed_four_ages", 0.8, 1, 3, "any_tissue"),
        ("fixed_four_ages_same_tissue", 0.8, 1, 3, "same_tissue"),
        ("TSI_0.7", 0.7, 1, 3, "none"),
        ("TSI_0.9", 0.9, 1, 3, "none"),
        ("TPM_0.1", 0.8, 0, 3, "none"),
        ("TPM_5", 0.8, 2, 3, "none"),
        ("reads_1", 0.8, 3, 3, "none"),
        ("reads_5", 0.8, 4, 3, "none"),
        ("replicates_2", 0.8, 1, 2, "none"),
        ("replicates_4", 0.8, 1, 4, "none"),
    ]
    for label, threshold, cutoff, repeats, fixed in specifications:
        detect = (
            (tx_support[cutoff] >= repeats)
            & (gene_support[cutoff][parents] >= repeats)
            & gene_present[parents, None, None]
        )
        if fixed == "any_tissue":
            detect &= detect.any(axis=2).all(axis=1)[:, None, None]
        elif fixed == "same_tissue":
            detect &= detect.all(axis=1)[:, None, :]
        specific = (tx_tsi > threshold) & detect
        gene_specific = gene_tsi.max(axis=2) > threshold
        core = specific & ~gene_specific[parents, :, None]
        assert (specific.sum(axis=2) <= 1).all()
        rows, tissue_rows, gene_eligible, gene_event, gene_core = summarize(
            detect, specific, core, gene_specific, parents, len(gene_ids), label, ages, tissues
        )
        all_summaries.extend(rows)
        all_tissues.extend(tissue_rows)
        if label == "primary":
            primary = dict(
                detect=detect.copy(),
                specific=specific.copy(),
                core=core.copy(),
                gene_specific=gene_specific.copy(),
                gene_eligible=gene_eligible.copy(),
                gene_event=gene_event.copy(),
                gene_core=gene_core.copy(),
            )
        print("Completed " + label, flush=True)
    summaries = pd.DataFrame(all_summaries)
    tissue_table = pd.DataFrame(all_tissues)
    summaries.to_csv(output / "overall_all_scenarios.tsv", sep="\t", index=False)
    tissue_table.to_csv(output / "by_tissue_all_scenarios.tsv", sep="\t", index=False)
    summaries[summaries.scenario == "primary"].to_csv(
        output / "primary_overall.tsv", sep="\t", index=False
    )
    tissue_table[tissue_table.scenario == "primary"].to_csv(
        output / "primary_by_tissue.tsv", sep="\t", index=False
    )
    eligible = primary["detect"].any(axis=2)
    specific = primary["specific"]
    core = primary["core"]
    primary_categories = categories(specific, eligible)
    core_categories = categories(core, eligible)
    detail = mapping[["transcript_id", "gene_id"]].copy()
    detail["eligible_age_count"] = eligible.sum(axis=1)
    detail["specific_age_count"] = specific.any(axis=2).sum(axis=1)
    detail["specific_tissue_count"] = specific.any(axis=1).sum(axis=1)
    detail["core_age_count"] = core.any(axis=2).sum(axis=1)
    detail["core_tissue_count"] = core.any(axis=1).sum(axis=1)
    detail["specificity_pattern"] = primary_categories
    detail["core_pattern"] = core_categories
    for age_index, age in enumerate(ages):
        active = specific[:, age_index].any(axis=1)
        detail[age + "_eligible"] = eligible[:, age_index]
        detail[age + "_tx_max_tsi"] = tx_tsi[:, age_index].max(axis=1)
        detail[age + "_tx_peak_tissue"] = np.array(tissues)[tx_tsi[:, age_index].argmax(axis=1)]
        detail[age + "_specific_tissue"] = np.where(
            active, np.array(tissues)[specific[:, age_index].argmax(axis=1)], ""
        )
        detail[age + "_parent_max_tsi"] = gene_tsi[:, age_index].max(axis=1)[parents]
        detail[age + "_core"] = core[:, age_index].any(axis=1)
    detail.to_csv(output / "transcript_age_classification.tsv.gz", sep="\t", index=False)
    event_rows = []
    for transcript_index, age_index, tissue_index in np.argwhere(specific):
        event_rows.append(
            dict(
                transcript_id=mapping.transcript_id.iloc[transcript_index],
                gene_id=mapping.gene_id.iloc[transcript_index],
                age=ages[age_index],
                tissue=tissues[tissue_index],
                tx_tsi=float(tx_tsi[transcript_index, age_index, tissue_index]),
                parent_max_tsi=float(gene_tsi[parents[transcript_index], age_index].max()),
                same_age_core=bool(core[transcript_index, age_index, tissue_index]),
                tx_support_samples=int(tx_support[1, transcript_index, age_index, tissue_index]),
                gene_support_samples=int(
                    gene_support[1, parents[transcript_index], age_index, tissue_index]
                ),
            )
        )
    events = pd.DataFrame(event_rows)
    assert not events.duplicated(["transcript_id", "age", "tissue"]).any()
    events.to_csv(output / "specific_events.tsv.gz", sep="\t", index=False)
    events[events.same_age_core].to_csv(
        output / "same_age_core_events.tsv.gz", sep="\t", index=False
    )
    category_rows = []
    for label, labels, event in [
        ("specificity", primary_categories, specific),
        ("core", core_categories, core),
    ]:
        selected = event.any(axis=(1, 2))
        for category in sorted(set(labels[selected])):
            mask = selected & (labels == category)
            category_rows.append(
                dict(
                    kind=label,
                    category=category,
                    transcripts=int(mask.sum()),
                    denominator=int(selected.sum()),
                    proportion=ratio(mask.sum(), selected.sum()),
                    parent_genes=len(np.unique(parents[mask])),
                )
            )
    pd.DataFrame(category_rows).to_csv(output / "cross_age_patterns.tsv", sep="\t", index=False)
    gene_detail = pd.DataFrame(dict(gene_id=gene_ids))
    for age_index, age in enumerate(ages):
        gene_detail[age + "_eligible"] = primary["gene_eligible"][:, age_index]
        gene_detail[age + "_gene_specific"] = (
            primary["gene_specific"][:, age_index] & primary["gene_eligible"][:, age_index]
        )
        gene_detail[age + "_has_specific_isoform"] = primary["gene_event"][:, age_index]
        gene_detail[age + "_has_core_isoform"] = primary["gene_core"][:, age_index]
    gene_detail["any_age_core"] = primary["gene_core"].any(axis=1)
    gene_detail["gene_specific_in_any_comparable_age"] = (
        primary["gene_specific"] & primary["gene_eligible"]
    ).any(axis=1)
    gene_detail["strict_never_specific_core"] = (
        gene_detail.any_age_core & ~gene_detail.gene_specific_in_any_comparable_age
    )
    gene_detail.to_csv(output / "gene_age_classification.tsv.gz", sep="\t", index=False)
    comparisons = []
    membership = []
    old_gene = (
        pd.read_csv(source / "primary_gene_classification.tsv.gz", sep="\t")
        .set_index("gene_id")
        .reindex(gene_ids)
    )
    old_gene_background = old_gene.eligible_transcripts.notna().to_numpy()
    new_gene_background = primary["gene_eligible"].any(axis=1)
    for event_name, old_positive, new_positive, labels, old_gene_positive, new_gene_positive in [
        (
            "specific",
            old.tx_specific.to_numpy(),
            specific.any(axis=(1, 2)),
            primary_categories,
            old_gene.specific_transcripts.fillna(0).to_numpy() > 0,
            primary["gene_event"].any(axis=1),
        ),
        (
            "core",
            old.parent_nonspecific_tx_specific.to_numpy(),
            core.any(axis=(1, 2)),
            core_categories,
            old_gene.hidden_transcripts.fillna(0).to_numpy() > 0,
            primary["gene_core"].any(axis=1),
        ),
    ]:
        status = np.select(
            [
                old_positive & new_positive,
                old_positive & ~new_positive,
                ~old_positive & new_positive,
            ],
            ["shared", "old_only", "new_only"],
            default="neither",
        )
        table = mapping[["transcript_id", "gene_id"]].copy()
        table["event"] = event_name
        table["comparison"] = status
        table["old_eligible"] = old.eligible.to_numpy()
        table["new_eligible"] = eligible.any(axis=1)
        table["new_age_pattern"] = labels
        table["gain_background"] = np.where(
            old.eligible.to_numpy(), "old_background_present", "outside_old_background"
        )
        membership.append(table)
        for level, old_bg, new_bg, old_event, new_event in [
            (
                "transcript",
                old.eligible.to_numpy(),
                eligible.any(axis=1),
                old_positive,
                new_positive,
            ),
            (
                "gene",
                old_gene_background,
                new_gene_background,
                old_gene_positive,
                new_gene_positive,
            ),
        ]:
            for scope, selection in [
                ("native_backgrounds", old_bg | new_bg),
                ("matched_background", old_bg & new_bg),
            ]:
                old_selected = old_event & selection
                new_selected = new_event & selection
                comparisons.append(
                    dict(
                        event=event_name,
                        level=level,
                        scope=scope,
                        background=int(selection.sum()),
                        old_background=int((old_bg & selection).sum()),
                        new_background=int((new_bg & selection).sum()),
                        old_positive=int(old_selected.sum()),
                        new_positive=int(new_selected.sum()),
                        shared=int((old_selected & new_selected).sum()),
                        old_only=int((old_selected & ~new_selected).sum()),
                        new_only=int((~old_selected & new_selected).sum()),
                        old_positive_over_background=ratio(
                            old_selected.sum(), (old_bg & selection).sum()
                        ),
                        new_positive_over_background=ratio(
                            new_selected.sum(), (new_bg & selection).sum()
                        ),
                    )
                )
        gene_comp = pd.DataFrame(
            dict(
                gene_id=gene_ids,
                old_eligible=old_gene_background,
                new_eligible=new_gene_background,
                old_positive=old_gene_positive,
                new_positive=new_gene_positive,
            )
        )
        gene_comp["comparison"] = np.select(
            [
                old_gene_positive & new_gene_positive,
                old_gene_positive & ~new_gene_positive,
                ~old_gene_positive & new_gene_positive,
            ],
            ["shared", "old_only", "new_only"],
            default="neither",
        )
        gene_comp.to_csv(output / f"gene_comparison_{event_name}.tsv.gz", sep="\t", index=False)
    comparison_detail = pd.concat(membership, ignore_index=True)
    comparison_detail.to_csv(
        output / "transcript_comparison_membership.tsv.gz", sep="\t", index=False
    )
    gains = (
        comparison_detail[comparison_detail.comparison == "new_only"]
        .groupby(["event", "gain_background", "new_age_pattern"])
        .size()
        .reset_index(name="transcripts")
    )
    gains.to_csv(output / "new_signal_decomposition.tsv", sep="\t", index=False)
    pd.DataFrame(comparisons).to_csv(output / "old_vs_union_comparison.tsv", sep="\t", index=False)
    for event_name in ["specific", "core"]:
        selected = comparison_detail[comparison_detail.event == event_name]
        for status in ["shared", "old_only", "new_only"]:
            selected[selected.comparison == status].to_csv(
                output / f"{event_name}_{status}_transcripts.tsv.gz", sep="\t", index=False
            )
    record = dict(
        status="passed",
        canonical_transcripts=len(mapping),
        ages=ages,
        tissues=tissues,
        method="within-age median, log1p TSI>0.8; same-age matched gene/transcript detection; union after paired classification",
        primary_detection="TPM>=1 and count>=3 in >=3/5 samples in that age and tissue for both levels",
        fixed_background="fixed_four_ages requires each age detectable in some tissue; fixed_four_ages_same_tissue requires the same tissue detectable in all four ages",
        comparable_gene_age="at least one eligible isoform-parent pair in that age",
        legacy_age_TSI="Not directly merged because legacy gene/transcript aggregation and ORF units differ; recomputed from previously verified unified sample counts.",
        current_counts_limitation="Existing unique-only quantification; discarded ambiguous reads are not recovered.",
        gene_union_warning="A gene can be gene-specific in one age and have a same-age core isoform in another; these union sets overlap.",
        expression_cache=str(source / "expression_background.npz"),
        source_cache_sha256=hashlib.sha256(
            (source / "expression_background.npz").read_bytes()
        ).hexdigest(),
        primary_overall=summaries[summaries.scenario == "primary"].to_dict(orient="records"),
        cross_age_patterns=category_rows,
        comparison=comparisons,
    )
    (output / "run_summary.json").write_text(json.dumps(record, indent=2, allow_nan=False))
    print(json.dumps(record["primary_overall"][-1]), flush=True)


if __name__ == "__main__":
    main()
