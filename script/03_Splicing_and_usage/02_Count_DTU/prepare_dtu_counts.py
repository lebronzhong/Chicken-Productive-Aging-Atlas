#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 03_Splicing_and_usage
# Purpose: Prepare tianjin.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

from pathlib import Path
import pandas as pd, numpy as np, json, re, sys, hashlib

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT.parent
mp = pd.read_csv(ROOT / "inputs/canonical_transcript_mapping.tsv.gz", sep="\t", low_memory=False)
meta = pd.read_csv(ROOT / "inputs/sample_manifest.tsv", sep="\t")
srcmap = {s: r.transcript_id for r in mp.itertuples() for s in r.source_ids.split(",")}
assert len(srcmap) == sum(mp.n_source_ids)
parents = mp.set_index("transcript_id").gene_id
tissues = sys.argv[1:] or sorted(meta.tissue.unique())
records = []
for tissue in tissues:
    out = ROOT / "inputs" / tissue
    out.mkdir(exist_ok=True)
    sm = meta[meta.tissue.eq(tissue)].copy()
    assert len(sm) == 20 and sm["sample"].is_unique and sm.groupby("age").size().eq(5).all()
    path = BASE / f"12.20240919_TPM.count/count/{tissue}_merged_transcript_grouped_counts.tsv"
    tx = pd.read_csv(path, sep="\t", index_col=0)
    tx.index = tx.index.astype(str)
    assert tx.index.is_unique
    missing = set(sm.source_sample) - set(tx.columns)
    assert not missing, missing
    original_rows = len(tx)
    tx = tx.loc[:, sm.source_sample]
    tx.columns = sm["sample"]
    assert np.isfinite(tx.to_numpy()).all() and tx.to_numpy().min() >= 0
    bad = tx.index[~tx.index.isin(srcmap)]
    assert len(bad) == 0, bad[:10]
    annotation = pd.read_csv(
        BASE / "29.animal_compare/Transcript_List_All_info.txt",
        sep="\t",
        header=None,
        names=["source", "gene", "reference", "symbol"],
        dtype=str,
    )
    annotated = annotation.groupby("source").gene.agg(lambda x: ";".join(sorted(set(x.dropna()))))
    conflicts = [
        s
        for s in tx.index
        if s in annotated.index and annotated[s] != parents[srcmap[s]].split("^")[0]
    ]
    assert not conflicts, conflicts[:10]
    profile_duplicates = 0
    for canon, group in pd.Series(tx.index, index=tx.index).groupby(tx.index.map(srcmap)):
        if len(group) > 1:
            vals = tx.loc[group.values]
            nonzero = vals[vals.sum(axis=1) > 0]
            profile_duplicates += int(nonzero.duplicated().sum())
    tx.index = tx.index.map(srcmap)
    tx = tx.groupby(level=0, sort=False).sum()
    assert np.abs(tx.to_numpy() - np.round(tx.to_numpy())).max() < 1e-5
    tx = np.round(tx).astype("int64")
    tx.index.name = "transcript_id"
    tx.to_csv(out / "tx_counts.tsv.gz", sep="\t")
    gene = pd.read_csv(
        BASE / f"8.combine_count/Count_gather/{tissue}_counts_final.txt", sep=r"\s+", index_col=0
    )
    assert gene.index.is_unique and set(sm["sample"]).issubset(gene.columns)
    gene = gene.loc[:, sm["sample"]]
    assert np.isfinite(gene.to_numpy()).all() and gene.to_numpy().min() >= 0
    assert np.abs(gene.to_numpy() - np.round(gene.to_numpy())).max() < 1e-5
    gene = np.round(gene).astype("int64")
    gene.index.name = "gene_id"
    gene.to_csv(out / "gene_counts.tsv.gz", sep="\t")
    sm.to_csv(out / "samples.tsv", sep="\t", index=False)
    mp[mp.transcript_id.isin(tx.index)].to_csv(out / "mapping.tsv.gz", sep="\t", index=False)
    record = dict(
        tissue=tissue,
        source_rows=original_rows,
        canonical_rows=len(tx),
        gene_rows=len(gene),
        natural_samples=20,
        duplicate_nonzero_ORF_profiles=profile_duplicates,
        absent_from_tissue_file=int((~mp.transcript_id.isin(tx.index)).sum()),
        tx_total=int(tx.to_numpy().sum()),
        gene_total=int(gene.to_numpy().sum()),
        source=str(path),
        source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
    )
    (out / "preparation.json").write_text(json.dumps(record, indent=2))
    records.append(record)
    print(json.dumps(record), flush=True)
annotation.to_csv(ROOT / "inputs/source_annotation.tsv.gz", sep="\t", index=False)
pd.DataFrame(records).to_csv(
    ROOT / "audit/input_preparation_" + ("_".join(tissues) if len(tissues) < 3 else "all") + ".tsv",
    sep="\t",
    index=False,
)
