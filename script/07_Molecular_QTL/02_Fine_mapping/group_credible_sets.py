#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 07_Molecular_QTL
# Purpose: Group credible sets.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

"""Pool the accepted Fig6b credible sets across tissues, separately by QTL type."""

from pathlib import Path
import json, hashlib
import pandas as pd

R = Path(__file__).resolve().parent
O = R / "pooled_CS_counts"
O.mkdir(exist_ok=True)
OLD = set(json.loads((R / "input_partitions.json").read_text(encoding="utf8"))["co_inputs"])
m = pd.read_csv(R / "b_gene_revision/b_CS_group_membership.tsv.gz", sep="\t")
summary = []


class DSU:
    def __init__(self, n):
        self.p = list(range(n))
        self.size = [1] * n

    def find(self, a):
        while self.p[a] != a:
            self.p[a] = self.p[self.p[a]]
            a = self.p[a]
        return a

    def union(self, a, b):
        a, b = self.find(a), self.find(b)
        if a == b:
            return
        if self.size[a] < self.size[b]:
            a, b = b, a
        self.p[b] = a
        self.size[a] += self.size[b]


for q in ["teQTL", "tuQTL"]:
    z = m[m.qtl_type == q].copy().reset_index(drop=True)
    assert not z[["tissue", "phenotype_id", "cs_id"]].duplicated().any()
    global_dsu = DSU(len(z))
    gene_dsu = DSU(len(z))
    seen = {}
    gene_seen = {}
    n_rows = 0
    lead_variants = set()
    input_hashes = []
    for t, sub in z.groupby("tissue", sort=True):
        p = R / "inputs" / ("co_inputs" if t in OLD else "tj_inputs") / f"{t}.{q}.finemap.tsv.gz"
        input_hashes.append(dict(file=str(p), sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
        f = pd.read_csv(p, sep="\t")
        f = f.merge(
            sub[["phenotype_id", "cs_id", "gene_id"]].assign(node=sub.index),
            on=["phenotype_id", "cs_id"],
            how="inner",
            validate="many_to_one",
        )
        f = f.drop_duplicates(["node", "variant_id"])
        assert set(f.node) == set(sub.index)
        counts = f.groupby("node").size()
        assert counts.eq(sub.set_index(sub.index).n_cs_variants).all()
        lead = f.sort_values(["pip", "variant_id"], ascending=[False, True]).drop_duplicates("node")
        lead_variants.update(lead.variant_id)
        for node, var, gene in f[["node", "variant_id", "gene_id"]].itertuples(
            index=False, name=None
        ):
            if var in seen:
                global_dsu.union(node, seen[var])
            else:
                seen[var] = node
            key = (gene, var)
            if key in gene_seen:
                gene_dsu.union(node, gene_seen[key])
            else:
                gene_seen[key] = node
        n_rows += len(f)
    roots = [global_dsu.find(i) for i in range(len(z))]
    identifiers = {v: f"{q}_pooledCS_{j:06d}" for j, v in enumerate(sorted(set(roots)), 1)}
    z["pooled_cs_group"] = [identifiers[v] for v in roots]
    z["gene_restricted_root"] = [gene_dsu.find(i) for i in range(len(z))]
    z.to_csv(O / f"{q}_CS_membership.tsv.gz", sep="\t", index=False)
    groups = (
        z.groupby("pooled_cs_group")
        .agg(
            n_original_cs=("cs_id", "size"),
            n_tissues=("tissue", "nunique"),
            n_genes=("gene_id", "nunique"),
            n_transcripts=("phenotype_id", "nunique"),
        )
        .reset_index()
    )
    variant_counts = {}
    for var, node in seen.items():
        g = identifiers[global_dsu.find(node)]
        variant_counts[g] = variant_counts.get(g, 0) + 1
    groups["n_unique_variants"] = groups.pooled_cs_group.map(variant_counts)
    assert groups.n_original_cs.sum() == len(z)
    assert groups.n_unique_variants.sum() == len(seen)
    groups.to_csv(O / f"{q}_pooled_groups.tsv", sep="\t", index=False)
    # Every original per-gene/tissue component must remain intact after pooling.
    assert z.groupby(["tissue", "gene_id", "cs_group"]).pooled_cs_group.nunique().eq(1).all()
    row = dict(
        qtl_type=q,
        n_tissues=z.tissue.nunique(),
        n_original_cs=len(z),
        n_gene_tissue_cs_groups=z.groupby(["tissue", "gene_id"]).cs_group.nunique().sum(),
        n_pooled_cs_groups=len(groups),
        n_gene_restricted_pooled_groups=z.gene_restricted_root.nunique(),
        n_unique_cs_lead_snps=len(lead_variants),
        n_unique_cs_member_snps=len(seen),
        n_cs_variant_rows=n_rows,
        n_multitissue_groups=int((groups.n_tissues > 1).sum()),
        n_multigene_groups=int((groups.n_genes > 1).sum()),
    )
    summary.append(row)
    (O / f"{q}_input_hashes.json").write_text(json.dumps(input_hashes, indent=2), encoding="utf-8")
    print(
        json.dumps({k: int(v) if not isinstance(v, str) else v for k, v in row.items()}), flush=True
    )
pd.DataFrame(summary).to_csv(O / "pooled_summary.tsv", sep="\t", index=False)
(O / "README.txt").write_text(
    "Sources and eligibility match the approved Fig6b gene-level revision. Within each QTL type, pool all 18 tissues and all genes, connect original transcript credible sets sharing >=1 identical variant ID, and count connected components, including transitive links. n_pooled_cs_groups is the genome-wide pooled overlap-group total; it is NOT a count of conditionally independent loci. n_gene_restricted_pooled_groups retains gene identity and therefore counts gene-signal units, not distinct genomic loci. n_unique_cs_lead_snps counts distinct maximum-PIP SNP IDs per original CS (ties: lexicographic variant ID). Neither SNP count is an independent-locus count. See membership and group tables for traceability.\n",
    encoding="utf-8",
)
