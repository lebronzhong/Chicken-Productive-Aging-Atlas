#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 08_Comparative_and_complex_traits
# Purpose: Build orthologs.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

from pathlib import Path
import json
import pandas as pd

root = Path("/path/to/chicken-aging/analysis/conservation")


def read(name):
    return pd.read_csv(root / "metadata" / name, sep="\t", dtype=str, keep_default_na=False)


x = read("chicken_mouse_orthologs_ensembl106.tsv")
raw_n = len(x)
x = x[
    (x.mmusculus_homolog_orthology_type == "ortholog_one2one")
    & (x.mmusculus_homolog_orthology_confidence == "1")
].copy()
x = x.rename(
    columns={
        "ensembl_gene_id": "gene_id",
        "external_gene_name": "chicken_symbol",
        "mmusculus_homolog_ensembl_gene": "mouse_ensembl_gene_id",
        "mmusculus_homolog_associated_gene_name": "mouse_symbol",
    }
)
assert not x.gene_id.duplicated().any()
assert not x.mouse_ensembl_gene_id.duplicated().any()
high_conf_n = len(x)
ent = read("mouse_entrez_ensembl106.tsv")
ent = ent[ent.entrezgene_id != ""][["ensembl_gene_id", "entrezgene_id"]].drop_duplicates()
# Restrict the identifier bridge to unique mappings in both directions.
ent = ent[~ent.ensembl_gene_id.duplicated(keep=False) & ~ent.entrezgene_id.duplicated(keep=False)]
ent = ent.rename(
    columns={"ensembl_gene_id": "mouse_ensembl_gene_id", "entrezgene_id": "mouse_entrez_id"}
)
x = x.merge(ent, on="mouse_ensembl_gene_id", how="inner", validate="one_to_one")
gtf = read("gene_annotation_GRCg6a_106.tsv")
x = x.merge(gtf[["gene_id", "gene_biotype"]], on="gene_id", how="inner", validate="one_to_one")
x.to_csv(
    root / "metadata" / "chicken_mouse_one2one_highconfidence_entrez.tsv", sep="\t", index=False
)
manifest = read("tissue_manifest.tsv")
coverage = []
for _, r in manifest.iterrows():
    ids = pd.read_csv(r.count_path, sep=r"\s+", usecols=[0], dtype=str).iloc[:, 0]
    coverage.append(
        dict(
            tissue=r.tissue, n_count_genes=len(ids), n_mappable_genes=int(ids.isin(x.gene_id).sum())
        )
    )
pd.DataFrame(coverage).to_csv(root / "qc" / "ortholog_coverage.tsv", sep="\t", index=False)
info = dict(
    raw_biomart_rows=raw_n,
    high_confidence_one2one_pairs=high_conf_n,
    unambiguous_entrez_and_gtf_pairs=len(x),
)
(root / "qc" / "ortholog_summary.json").write_text(json.dumps(info, indent=2))
print(json.dumps(info, indent=2))
