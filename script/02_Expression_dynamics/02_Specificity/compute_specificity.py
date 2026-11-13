#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 02_Expression_dynamics
# Purpose: Compute specificity.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

import argparse
import gzip
import json
import re
from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd


def parse_gtf(path):
    models = {}
    with path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9 or fields[2] not in ("transcript", "exon"):
                continue
            attrs = dict(re.findall(r'(\w+) "([^"]*)"', fields[8]))
            identifier = attrs.get("transcript_id", "")
            if not identifier.startswith("TCONS_"):
                continue
            model = models.setdefault(
                identifier, dict(gene=attrs["gene_id"], chrom=fields[0], strand=fields[6], exons=[])
            )
            if fields[2] == "exon":
                model["exons"].append((int(fields[3]), int(fields[4])))
    return models


def matrix(path):
    frame = pd.read_csv(path, sep="\t", index_col=0)
    assert frame.index.is_unique and frame.columns.is_unique, str(path)
    values = frame.to_numpy(dtype=np.float64)
    assert np.isfinite(values).all() and (values >= 0).all(), str(path)
    return frame


def share(values, transform):
    transformed = np.log1p(values) if transform == "log1p" else values
    totals = transformed.sum(axis=1, keepdims=True)
    return np.divide(transformed, totals, out=np.zeros_like(transformed), where=totals > 0)


def fraction(numerator, denominator):
    return float(numerator / denominator) if denominator else None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    cache = output / "normalized_counts"
    cache.mkdir(exist_ok=True)
    root = Path("/path/to/chicken-aging/ont")
    annot = root / "shell/10_20240903.gffcompare/0313_GTF_update"
    ages = ["20wk", "50wk", "80wk", "100wk"]
    with (
        root / "shell/16.Transcripts_TSI/TSI/1.merged_transcript_tpm_filtered_20wk.tsv"
    ).open() as handle:
        tissues = handle.readline().strip().split("\t")[1:]
    models = parse_gtf(annot / "union.topo.gtf")
    mapping = pd.read_csv(annot / "tx2gene.tsv", sep="\t", dtype=str)
    assert not mapping.transcript_id.duplicated().any()
    mapping = mapping.set_index("transcript_id").gene_id.to_dict()
    assert all(
        mapping.get(identifier) == model["gene"] for identifier, model in models.items()
    ), "Authoritative mapping/GTF disagreement"
    signatures = defaultdict(list)
    for identifier, model in models.items():
        exons = tuple(sorted(set(model["exons"])))
        assert exons
        signatures[(model["gene"], model["chrom"], model["strand"], exons)].append(identifier)
    rows = []
    full_to_index = {}
    for signature, identifiers in sorted(signatures.items(), key=lambda item: min(item[1])):
        bases = sorted({re.sub(r"\.p\d+$", "", identifier) for identifier in identifiers})
        canonical = bases[0]
        for identifier in identifiers:
            full_to_index[identifier] = len(rows)
        rows.append(
            dict(
                transcript_id=canonical,
                gene_id=signature[0],
                source_ids=",".join(sorted(identifiers)),
                base_ids=",".join(bases),
                n_source_ids=len(identifiers),
                chromosome=signature[1],
                strand=signature[2],
                exon_count=len(signature[3]),
                exons=";".join(f"{start}-{end}" for start, end in signature[3]),
            )
        )
    info = pd.DataFrame(rows)
    assert info.transcript_id.is_unique, "Same base ID has multiple structures"
    info.to_csv(output / "canonical_transcript_mapping.tsv.gz", sep="\t", index=False)
    gene_ids = sorted(info.gene_id.unique())
    gene_index = {gene: index for index, gene in enumerate(gene_ids)}
    parent = info.gene_id.map(gene_index).to_numpy()
    shape_tx = (len(info), len(ages), len(tissues))
    shape_gene = (len(gene_ids), len(ages), len(tissues))
    tx_mean = np.zeros(shape_tx, dtype=np.float32)
    tx_median = np.zeros(shape_tx, dtype=np.float32)
    gene_mean = np.zeros(shape_gene, dtype=np.float32)
    gene_median = np.zeros(shape_gene, dtype=np.float32)
    gene_sum_mean = np.zeros(shape_gene, dtype=np.float32)
    gene_sum_median = np.zeros(shape_gene, dtype=np.float32)
    cutoffs = [(0.1, 3), (1, 3), (5, 3), (1, 1), (1, 5)]
    tx_support = np.zeros((len(cutoffs),) + shape_tx, dtype=np.uint8)
    gene_support = np.zeros((len(cutoffs),) + shape_gene, dtype=np.uint8)
    sum_support = np.zeros((len(cutoffs),) + shape_gene, dtype=np.uint8)
    gene_present = np.ones(len(gene_ids), dtype=bool)
    input_qc = []
    for age_index, age in enumerate(ages):
        for tissue_index, tissue in enumerate(tissues):
            folder = root / (age + "_data") / (age + "_ONTdata") / tissue
            tx_frame = matrix(folder / (tissue + ".transcript_grouped_counts.tsv"))
            gene_frame = matrix(folder / (tissue + ".gene_grouped_counts.tsv"))
            assert list(tx_frame.columns) == list(gene_frame.columns)
            assert len(tx_frame.columns) == 5, "Expected five samples per cell"
            assert set(tx_frame.index) == set(full_to_index), "Incomplete transcript input"
            collapse = np.array([full_to_index[identifier] for identifier in tx_frame.index])
            tx_counts = np.zeros((len(info), len(tx_frame.columns)), dtype=np.float64)
            np.add.at(tx_counts, collapse, tx_frame.to_numpy())
            assert np.allclose(tx_counts.sum(axis=0), tx_frame.to_numpy().sum(axis=0))
            tx_total = tx_counts.sum(axis=0)
            gene_total = gene_frame.to_numpy().sum(axis=0)
            assert (tx_total > 0).all() and (gene_total > 0).all()
            present = np.array([gene in gene_frame.index for gene in gene_ids])
            gene_present &= present
            gene_counts = gene_frame.reindex(gene_ids).fillna(0).to_numpy()
            tx_tpm = tx_counts / tx_total * 1e6
            gene_tpm = gene_counts / gene_total * 1e6
            summed_counts = np.zeros_like(gene_counts)
            np.add.at(summed_counts, parent, tx_counts)
            sum_tpm = summed_counts / tx_total * 1e6
            tx_mean[:, age_index, tissue_index] = tx_tpm.mean(axis=1)
            tx_median[:, age_index, tissue_index] = np.median(tx_tpm, axis=1)
            gene_mean[:, age_index, tissue_index] = gene_tpm.mean(axis=1)
            gene_median[:, age_index, tissue_index] = np.median(gene_tpm, axis=1)
            gene_sum_mean[:, age_index, tissue_index] = sum_tpm.mean(axis=1)
            gene_sum_median[:, age_index, tissue_index] = np.median(sum_tpm, axis=1)
            for cutoff_index, (tpm_cutoff, read_cutoff) in enumerate(cutoffs):
                tx_support[cutoff_index, :, age_index, tissue_index] = (
                    (tx_tpm >= tpm_cutoff) & (tx_counts >= read_cutoff)
                ).sum(axis=1)
                gene_support[cutoff_index, :, age_index, tissue_index] = (
                    (gene_tpm >= tpm_cutoff) & (gene_counts >= read_cutoff)
                ).sum(axis=1)
                sum_support[cutoff_index, :, age_index, tissue_index] = (
                    (sum_tpm >= tpm_cutoff) & (summed_counts >= read_cutoff)
                ).sum(axis=1)
            np.savez_compressed(
                cache / f"{age}_{tissue}.npz",
                tx_counts=tx_counts.astype(np.float32),
                gene_counts=gene_counts.astype(np.float32),
                tx_total=tx_total,
                gene_total=gene_total,
                samples=np.array(tx_frame.columns),
            )
            input_qc.append(
                dict(
                    age=age,
                    tissue=tissue,
                    samples=";".join(tx_frame.columns),
                    full_tx_rows=len(tx_frame),
                    canonical_rows=len(info),
                    gene_rows=len(gene_frame),
                    missing_mapping_genes=int((~present).sum()),
                    assigned_tx_counts=float(tx_counts.sum()),
                    assigned_gene_counts=float(gene_frame.to_numpy().sum()),
                )
            )
            print(f"Normalized {age} {tissue}", flush=True)
    pd.DataFrame(input_qc).to_csv(output / "input_qc.tsv", sep="\t", index=False)
    np.savez_compressed(
        output / "expression_background.npz",
        tx_mean=tx_mean,
        tx_median=tx_median,
        gene_mean=gene_mean,
        gene_median=gene_median,
        gene_sum_mean=gene_sum_mean,
        gene_sum_median=gene_sum_median,
        tx_support=tx_support,
        gene_support=gene_support,
        sum_support=sum_support,
        gene_present=gene_present,
        parent=parent,
        gene_ids=np.array(gene_ids),
        ages=np.array(ages),
        tissues=np.array(tissues),
    )
    specifications = [dict(name="primary")]
    specifications += [
        dict(name=f"TSI_{threshold}", threshold=threshold) for threshold in [0.7, 0.9]
    ]
    specifications += [
        dict(name=f"TPM_{cutoff}", cutoff=index) for index, cutoff in [(0, 0.1), (2, 5)]
    ]
    specifications += [
        dict(name=f"reads_{cutoff}", cutoff=index) for index, cutoff in [(3, 1), (4, 5)]
    ]
    specifications += [dict(name=f"replicates_{number}", replicates=number) for number in [2, 4]]
    specifications += [dict(name=f"age_support_{number}", min_ages=number) for number in [1, 3, 4]]
    specifications += [
        dict(name="mean_instead_of_median", aggregation="mean"),
        dict(name="linear_TSI", transform="linear"),
        dict(name="gene_sum_of_transcripts", gene_source="sum"),
        dict(name="single_source_id_only", single_source=True),
    ]
    specifications += [
        dict(
            name=f"leave_out_{age}", age_indices=[index for index in range(4) if index != excluded]
        )
        for excluded, age in enumerate(ages)
    ]
    overall = []
    by_tissue = []
    primary = None
    for spec in specifications:
        params = dict(
            threshold=0.8,
            cutoff=1,
            replicates=3,
            min_ages=2,
            aggregation="median",
            transform="log1p",
            gene_source="direct",
            single_source=False,
            age_indices=list(range(4)),
        )
        params.update(spec)
        selected_ages = params["age_indices"]
        tx_values = (tx_median if params["aggregation"] == "median" else tx_mean)[
            :, selected_ages, :
        ].mean(axis=1)
        if params["gene_source"] == "sum":
            gene_values = (gene_sum_median if params["aggregation"] == "median" else gene_sum_mean)[
                :, selected_ages, :
            ].mean(axis=1)
            gene_rep = sum_support[params["cutoff"]]
        else:
            gene_values = (gene_median if params["aggregation"] == "median" else gene_mean)[
                :, selected_ages, :
            ].mean(axis=1)
            gene_rep = gene_support[params["cutoff"]]
        tx_detect = (tx_support[params["cutoff"]][:, selected_ages, :] >= params["replicates"]).sum(
            axis=1
        ) >= params["min_ages"]
        gene_detect = (gene_rep[:, selected_ages, :] >= params["replicates"]).sum(axis=1) >= params[
            "min_ages"
        ]
        eligible = (tx_detect & gene_detect[parent]).any(axis=1) & gene_present[parent]
        if params["single_source"]:
            eligible &= info.n_source_ids.to_numpy() == 1
        tx_tsi = share(tx_values, params["transform"])
        gene_tsi = share(gene_values, params["transform"])
        tx_specific = (
            (tx_tsi > params["threshold"]) & tx_detect & gene_detect[parent] & eligible[:, None]
        )
        gene_specific_any = gene_tsi.max(axis=1) > params["threshold"]
        hidden_any = tx_specific.any(axis=1) & ~gene_specific_any[parent]
        background_genes = np.unique(parent[eligible])
        hidden_genes = np.unique(parent[hidden_any])
        non_specific_genes = background_genes[~gene_specific_any[background_genes]]
        row = dict(
            scenario=params["name"],
            eligible_transcripts=int(eligible.sum()),
            eligible_genes=len(background_genes),
            non_specific_parent_genes=len(non_specific_genes),
            specific_transcripts=int(tx_specific.any(axis=1).sum()),
            hidden_transcripts=int(hidden_any.sum()),
            hidden_genes=len(hidden_genes),
            hidden_tx_over_all_eligible_tx=fraction(hidden_any.sum(), eligible.sum()),
            hidden_tx_over_specific_tx=fraction(hidden_any.sum(), tx_specific.any(axis=1).sum()),
            hidden_gene_over_all_eligible_gene=fraction(len(hidden_genes), len(background_genes)),
            hidden_gene_over_non_specific_gene=fraction(len(hidden_genes), len(non_specific_genes)),
        )
        overall.append(row)
        for tissue_index, tissue in enumerate(tissues):
            target_background = (
                eligible & tx_detect[:, tissue_index] & gene_detect[parent, tissue_index]
            )
            target_genes = np.unique(parent[target_background])
            hidden = hidden_any & tx_specific[:, tissue_index] & target_background
            weak = (
                tx_specific[:, tissue_index]
                & (gene_tsi[parent, tissue_index] <= params["threshold"])
                & target_background
            )
            target_non_specific = target_genes[~gene_specific_any[target_genes]]
            genes_hidden = np.unique(parent[hidden])
            by_tissue.append(
                dict(
                    scenario=params["name"],
                    tissue=tissue,
                    eligible_transcripts_in_tissue=int(target_background.sum()),
                    eligible_genes_in_tissue=len(target_genes),
                    non_specific_genes_in_tissue=len(target_non_specific),
                    hidden_transcripts=int(hidden.sum()),
                    hidden_genes=len(genes_hidden),
                    tissue_only_gene_nonspecific_transcripts=int(weak.sum()),
                    hidden_tx_over_tissue_tx=fraction(hidden.sum(), target_background.sum()),
                    hidden_gene_over_tissue_gene=fraction(len(genes_hidden), len(target_genes)),
                    hidden_gene_over_tissue_nonspecific_gene=fraction(
                        len(genes_hidden), len(target_non_specific)
                    ),
                    all_tissue_eligible_tx=int(eligible.sum()),
                    all_tissue_eligible_gene=len(background_genes),
                    hidden_tx_over_global_tx=fraction(hidden.sum(), eligible.sum()),
                    hidden_gene_over_global_gene=fraction(len(genes_hidden), len(background_genes)),
                )
            )
        if params["name"] == "primary":
            primary = dict(
                eligible=eligible,
                tx_specific=tx_specific,
                hidden=hidden_any,
                tx_values=tx_values,
                gene_values=gene_values,
                tx_tsi=tx_tsi,
                gene_tsi=gene_tsi,
                gene_specific_any=gene_specific_any,
                tx_detect=tx_detect,
                gene_detect=gene_detect,
            )
    pd.DataFrame(overall).to_csv(output / "overall_statistics.tsv", sep="\t", index=False)
    pd.DataFrame(by_tissue).to_csv(output / "by_tissue_statistics.tsv", sep="\t", index=False)
    details = info[["transcript_id", "gene_id", "n_source_ids", "exon_count"]].copy()
    details["eligible"] = primary["eligible"]
    details["tx_max_tsi"] = primary["tx_tsi"].max(axis=1)
    details["tx_preferred_tissue"] = np.array(tissues)[primary["tx_tsi"].argmax(axis=1)]
    details["parent_max_tsi"] = primary["gene_tsi"].max(axis=1)[parent]
    details["tx_specific"] = primary["tx_specific"].any(axis=1)
    details["parent_nonspecific_tx_specific"] = primary["hidden"]
    details.to_csv(output / "primary_transcript_classification.tsv.gz", sep="\t", index=False)
    gene_rows = []
    for gene, group in details[details.eligible].groupby("gene_id"):
        gene_rows.append(
            dict(
                gene_id=gene,
                eligible_transcripts=len(group),
                specific_transcripts=int(group.tx_specific.sum()),
                hidden_transcripts=int(group.parent_nonspecific_tx_specific.sum()),
                parent_max_tsi=float(group.parent_max_tsi.iloc[0]),
                preferred_tissues=",".join(
                    sorted(group.loc[group.tx_specific, "tx_preferred_tissue"].unique())
                ),
            )
        )
    pd.DataFrame(gene_rows).to_csv(
        output / "primary_gene_classification.tsv.gz", sep="\t", index=False
    )
    candidates = []
    for gene, group in details[details.parent_nonspecific_tx_specific].groupby("gene_id"):
        if group.tx_preferred_tissue.nunique() < 2:
            continue
        for transcript_index, item in group.iterrows():
            preferred = tissues.index(item.tx_preferred_tissue)
            consistent = int(
                (
                    (tx_support[4, transcript_index, :, preferred] >= 3)
                    & (tx_median[transcript_index].argmax(axis=1) == preferred)
                ).sum()
            )
            candidates.append(
                dict(
                    transcript_id=item.transcript_id,
                    gene_id=gene,
                    tissue=item.tx_preferred_tissue,
                    tx_tsi=float(item.tx_max_tsi),
                    parent_tsi=float(item.parent_max_tsi),
                    supported_ages_reads5_reps3=consistent,
                    n_source_ids=int(item.n_source_ids),
                    preferred_tpm=float(primary["tx_values"][transcript_index, preferred]),
                    status="expression_candidate_structure_support_not_yet_audited",
                )
            )
    pd.DataFrame(
        candidates,
        columns=[
            "transcript_id",
            "gene_id",
            "tissue",
            "tx_tsi",
            "parent_tsi",
            "supported_ages_reads5_reps3",
            "n_source_ids",
            "preferred_tpm",
            "status",
        ],
    ).to_csv(output / "opposing_isoform_candidates.tsv", sep="\t", index=False)
    method = dict(
        ages=ages,
        tissues=tissues,
        canonical_transcripts=len(info),
        mapped_genes=len(gene_ids),
        primary=dict(
            TSI="log1p(tissue_expression)/sum_over_33_tissues(log1p(tissue_expression))",
            expression="within_age_tissue_sample_median_then_equal_mean_over_four_ages",
            normalization="IsoQuant-style assigned counts per million, no length correction",
            TSI_threshold="strictly greater than 0.8",
            detection="TPM>=1 and count>=3 in >=3/5 samples in >=2/4 ages in at least one tissue",
            hidden="transcript passes TSI in a detectable tissue and parent max TSI<=0.8 across all tissues",
        ),
        limitation="Existing unique-only counts are collapsed without recovering discarded ambiguous reads; conclusions concern the jointly detectable quantified background.",
        source_mapping=str(annot / "tx2gene.tsv"),
        source_gtf=str(annot / "union.topo.gtf"),
        scenarios=specifications,
    )
    (output / "methods.json").write_text(json.dumps(method, ensure_ascii=False, indent=2))
    (output / "SUCCESS.json").write_text(
        json.dumps(dict(primary=overall[0], candidate_transcripts=len(candidates)), indent=2)
    )
    print(json.dumps(overall[0]), flush=True)


if __name__ == "__main__":
    main()
