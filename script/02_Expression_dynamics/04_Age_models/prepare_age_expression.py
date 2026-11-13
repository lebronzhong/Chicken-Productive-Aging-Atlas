#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 02_Expression_dynamics
# Purpose: Prepare.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

import hashlib
import json
import pathlib
import re
import shutil
import numpy as np
import pandas as pd

out = pathlib.Path(__file__).resolve().parent
root = pathlib.Path("/path/to/chicken-aging/ont")
base = root / "9.29_buchong/specificity/results"
for name in [
    "01_background",
    "02_age_expression",
    "03_clusters",
    "04_within_gene",
    "05_cross_tissue",
    "06_examples",
    "07_conservation",
    "08_figures",
    "09_report",
    "logs",
]:
    (out / name).mkdir(exist_ok=True)
mapping = pd.read_csv(base / "canonical_transcript_mapping.tsv.gz", sep="\t")
cache = np.load(base / "expression_background.npz", allow_pickle=True)
ages = cache["ages"].tolist()
tissues = cache["tissues"].tolist()
genes = cache["gene_ids"]
assert mapping.transcript_id.is_unique
assert np.array_equal(mapping.gene_id.to_numpy(), genes[cache["parent"]])
shutil.copy2(
    base / "canonical_transcript_mapping.tsv.gz",
    out / "01_background/canonical_transcript_mapping.tsv.gz",
)
quality = pd.read_csv(
    root
    / "9.29_buchong/Fig2_quality_integration_20261020/results/transcript_quality_integrated.tsv.gz",
    sep="\t",
    low_memory=False,
)
columns = ["transcript_id", "gene_id", "class"] + [
    column
    for column in quality
    if column
    in [
        "sqanti_all_canonical",
        "sqanti_RTS_stage",
        "sqanti_coding",
        "sqanti_predicted_NMD",
        "lr_exact_chain_evaluation",
    ]
]
quality[columns].to_csv(out / "01_background/quality.tsv.gz", sep="\t", index=False)
(out / "01_background/tissues.txt").write_text("\n".join(tissues) + "\n")
records, audits = [], []
for tissue in tissues:
    destination = out / "01_background" / tissue
    destination.mkdir(exist_ok=True)
    tx_cells, gene_cells, samples = [], [], []
    for age in ages:
        cell = np.load(base / f"normalized_counts/{age}_{tissue}.npz", allow_pickle=True)
        tx_cells.append(cell["tx_counts"])
        gene_cells.append(cell["gene_counts"])
        for source_sample in cell["samples"]:
            individual = re.search(r"(\d+)$", str(source_sample)).group(1)
            if age == "80wk" and "-" not in str(source_sample):
                individual = str(source_sample)[len(tissue) :]
            age_number = age.replace("wk", "")
            canonical = "wk" + age_number + individual
            samples.append(canonical)
            records.append(
                dict(
                    tissue=tissue,
                    sample=canonical,
                    source_sample=str(source_sample),
                    age=int(age_number),
                    individual=age_number + "_" + individual,
                    batch="B_50_100" if age_number in ["50", "100"] else "B_" + age_number,
                )
            )
    assert len(samples) == 20 and len(set(samples)) == 20
    tx_values = np.concatenate(tx_cells, axis=1)
    cached_gene = np.concatenate(gene_cells, axis=1)
    gene_path = root / f"shell/8.combine_count/Count_gather/{tissue}_counts_final.txt"
    gene_frame = pd.read_csv(gene_path, sep=r"\s+", index_col=0)
    assert set(samples).issubset(gene_frame.columns)
    natural_gene = gene_frame.loc[:, samples]
    compared_gene = natural_gene.reindex(genes).fillna(0).to_numpy()
    gene_equal = np.allclose(cached_gene, compared_gene, atol=1e-5)
    tx_source = (
        root / f"shell/12_20240919_TPM.count/count/{tissue}_merged_transcript_grouped_counts.tsv"
    )
    if tissue in ["Adipose", "Heart", "Retina"]:
        raw = pd.read_csv(tx_source, sep="\t", index_col=0)
        source_names = [record["source_sample"] for record in records if record["tissue"] == tissue]
        source_to_canonical = {
            source: row.transcript_id
            for row in mapping.itertuples()
            for source in row.source_ids.split(",")
        }
        selected = raw.loc[:, source_names]
        selected.index = selected.index.map(source_to_canonical)
        assert selected.index.notna().all()
        grouped = selected.groupby(level=0, sort=False).sum()
        absent = ~mapping.transcript_id.isin(grouped.index).to_numpy()
        collapsed = grouped.reindex(mapping.transcript_id).fillna(0).to_numpy()
        assert np.allclose(collapsed[~absent], tx_values[~absent], atol=1e-5), tissue
        pd.DataFrame(
            [
                dict(
                    tissue=tissue,
                    supplied_rows=len(grouped),
                    omitted_canonical_rows=int(absent.sum()),
                    omitted_assigned_counts=float(tx_values[absent].sum()),
                    max_omitted_sample_count=float(tx_values[absent].max()),
                    overlapping_counts_identical=True,
                )
            ]
        ).to_csv(destination / "legacy_prefilter_audit.tsv", sep="\t", index=False)
    for feature, values in [
        ("tx", tx_values),
        ("gene", cached_gene),
        ("legacy_gene", natural_gene.to_numpy()),
    ]:
        assert np.isfinite(values).all() and (values >= 0).all()
        deviation = np.abs(values - np.round(values)).max()
        assert deviation < 1e-5, (
            tissue,
            feature,
            "fractional counts require modelling review",
            deviation,
        )
    tx_frame = pd.DataFrame(
        np.round(tx_values).astype(np.int64), index=mapping.transcript_id, columns=samples
    )
    tx_frame.index.name = "feature_id"
    natural_gene.index.name = "feature_id"
    paired_gene = pd.DataFrame(np.round(cached_gene).astype(np.int64), index=genes, columns=samples)
    paired_gene.index.name = "feature_id"
    tx_frame.to_csv(destination / "tx_counts.tsv.gz", sep="\t")
    paired_gene.to_csv(destination / "gene_counts.tsv.gz", sep="\t")
    natural_gene.astype(np.int64).to_csv(destination / "legacy_gene_counts.tsv.gz", sep="\t")
    audits.append(
        dict(
            tissue=tissue,
            canonical_transcripts=len(mapping),
            gene_rows=len(paired_gene),
            legacy_gene_rows=len(natural_gene),
            natural_samples=len(samples),
            gene_cache_matches_supplied=gene_equal,
            legacy_missing_canonical_genes=int((~pd.Index(genes).isin(natural_gene.index)).sum()),
            raw_tx_checked=tissue in ["Adipose", "Heart", "Retina"],
            tx_assigned_reads=int(tx_values.sum()),
            gene_assigned_reads=int(cached_gene.sum()),
        )
    )
    print("prepared", tissue, flush=True)
pd.DataFrame(records).to_csv(out / "01_background/sample_manifest.tsv", sep="\t", index=False)
pd.DataFrame(audits).to_csv(out / "01_background/input_validation.tsv", sep="\t", index=False)
cluster_source = (
    root
    / "shell/14.5.LOESS_Cluster_ALL/2.Cluster/Cluster_Results_V3/Clustered_transcripts_with_labels.txt"
)
cluster = pd.read_csv(cluster_source, sep="\t")
source_to_canonical = {
    source: row.transcript_id
    for row in mapping.itertuples()
    for source in row.source_ids.split(",")
}
original_ids = cluster.iloc[:, 0].str.extract(r"^(.*)_(TCONS_\d+\.p\d+)$")
cluster["tissue"] = original_ids[0]
cluster["source_transcript_id"] = original_ids[1]
cluster["transcript_id"] = cluster.source_transcript_id.map(source_to_canonical)
cluster.to_csv(out / "03_clusters/old_members_mapped.tsv.gz", sep="\t", index=False)
centres = (
    root
    / "shell/14.5.LOESS_Cluster_ALL/2.Cluster/Cluster_Results_V3/Cluster_Mean_Expression_with_Counts.txt"
)
shutil.copy2(centres, out / "03_clusters/old_centres.tsv")
json.dump(
    {
        "status": "prepared",
        "tissues": len(tissues),
        "canonical_transcripts": len(mapping),
        "age_batch_confounded": True,
        "primary_contrast": "100wk_vs_50wk",
        "raw_tx_comparison_tissues": ["Adipose", "Heart", "Retina"],
    },
    open(out / "01_background/preparation.json", "w"),
    indent=2,
)
