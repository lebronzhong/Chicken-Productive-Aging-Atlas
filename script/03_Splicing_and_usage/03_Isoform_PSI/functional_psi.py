#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 03_Splicing_and_usage
# Purpose: Functional psi.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

from pathlib import Path
import pandas as pd, numpy as np, json
from scipy.stats import hypergeom, false_discovery_control

r = Path(__file__).resolve().parents[1]
old = r.parent / "count_DTU"
out = r / "results/summary"


def read(p):
    return pd.read_csv(p, sep="\t", low_memory=False)


def save(d, n):
    d.to_csv(out / n, sep="\t", index=False)


m = read(out / "original_clusters_with_new_evidence.tsv.gz")
go = (
    read(old / "inputs/gene_GO.tsv.gz")
    .dropna(subset=["gene_id", "GO"])
    .drop_duplicates(["gene_id", "GO"])
)
info = go[["GO", "ontology", "TERM"]].drop_duplicates("GO").set_index("GO")
sets = go.groupby("GO").gene_id.agg(set).to_dict()
annotated = set(go.gene_id)
enrich = []
coverage = []
sets_rows = []


def test_set(scope, tissue, cluster, foreground, background):
    fg = set(foreground)
    bg = set(background)
    assert fg <= bg
    universe = bg & annotated
    target = fg & universe
    B = len(universe)
    N = len(target)
    coverage.append(
        dict(
            scope=scope,
            tissue=tissue,
            cluster=cluster,
            background_genes=len(bg),
            annotated_background=B,
            foreground_genes=len(fg),
            annotated_foreground=N,
        )
    )
    sets_rows.append(
        dict(
            scope=scope,
            tissue=tissue,
            cluster=cluster,
            foreground_genes=";".join(sorted(fg)),
            background_genes=";".join(sorted(bg)),
        )
    )
    if N < 5:
        return
    for term, s in sets.items():
        ss = s & universe
        K = len(ss)
        if not 10 <= K <= 500:
            continue
        k = len(ss & target)
        enrich.append(
            dict(
                scope=scope,
                tissue=tissue,
                cluster=cluster,
                GO=term,
                ontology=info.loc[term, "ontology"],
                term=info.loc[term, "TERM"],
                background=B,
                target=N,
                background_term=K,
                hits=k,
                PValue=float(hypergeom.sf(k - 1, B, K, N)),
                fold_enrichment=(k / N) / (K / B),
                hit_genes=";".join(sorted(ss & target)),
            )
        )


for t, d in m.groupby("tissue"):
    bg = set(d.loc[d.testable, "gene_id"])
    for c in range(1, 11):
        test_set(
            "tissue_supported", t, c, set(d.loc[d.cluster.eq(c) & d.significant, "gene_id"]), bg
        )
# A common eligibility rule for foreground and background of the legacy >=17-tissue analysis.
eligible = m[m.testable]
shared = eligible.groupby("transcript_id").tissue.nunique()
shared_ids = set(shared[shared >= 17].index)
shared_bg = set(eligible.loc[eligible.transcript_id.isin(shared_ids), "gene_id"])
for c in range(1, 11):
    dd = eligible[eligible.cluster.eq(c)]
    n = dd.groupby("transcript_id").tissue.nunique()
    fg = set(dd.loc[dd.transcript_id.isin(n[n >= 17].index), "gene_id"])
    test_set("shared17_testable_descriptive", "Across_tissues", c, fg, shared_bg)
e = pd.DataFrame(enrich)
e["q_global"] = np.nan
# Families are declared separately for the primary tissue-supported and legacy-sharing questions.
for scope, ix in e.groupby("scope").groups.items():
    e.loc[ix, "q_global"] = false_discovery_control(e.loc[ix, "PValue"])
save(e, "GO_all_tests.tsv.gz")
save(e[e.q_global.lt(0.05) & e.hits.ge(3)].sort_values("q_global"), "GO_significant.tsv")
save(pd.DataFrame(coverage), "GO_background_coverage.tsv")
save(pd.DataFrame(sets_rows), "GO_foreground_background_genes.tsv.gz")
# One gain/loss pair per gene-tissue-contrast, selected solely by signed original PSI changes.
z = read(out / "all_PSI_feature_tests.tsv.gz")
z = z[z.scale.eq("arcsine_sqrt") & z.contrast.ne("four_age") & z.significant].copy()
g = read(out / "all_PSI_gene_tests.tsv.gz")
hits = g[g.PSI_significant & g.contrast.ne("four_age")][
    ["tissue", "gene_id", "contrast", "q_global"]
].rename(columns={"q_global": "gene_q"})
z = z.merge(hits, on=["tissue", "gene_id", "contrast"])
mp = read(old / "inputs/canonical_transcript_mapping.tsv.gz").set_index("transcript_id")
ann = read(old / "inputs/canonical_annotations.tsv.gz").set_index("transcript_id")
rows = []


def true(x):
    return str(x).upper() == "TRUE"


def exons(x):
    return [tuple(map(int, e.split("-"))) for e in str(x).split(";")]


def junctions(e):
    return {(e[i][1], e[i + 1][0]) for i in range(len(e) - 1)}


for key, d in z.groupby(["tissue", "gene_id", "contrast"]):
    up = d[d.delta_PSI.gt(0)].sort_values(["delta_PSI", "transcript_id"], ascending=[False, True])
    dn = d[d.delta_PSI.lt(0)].sort_values(["delta_PSI", "transcript_id"])
    if up.empty or dn.empty:
        continue
    u, v = up.iloc[0], dn.iloc[0]
    um, vm = mp.loc[u.transcript_id], mp.loc[v.transcript_id]
    a, b = ann.loc[u.transcript_id], ann.loc[v.transcript_id]
    ue, ve = exons(um.exons), exons(vm.exons)
    uj, vj = junctions(ue), junctions(ve)
    uf, ul = (ue[0], ue[-1]) if um.strand == "+" else (ue[-1], ue[0])
    vf, vl = (ve[0], ve[-1]) if vm.strand == "+" else (ve[-1], ve[0])
    rec = dict(
        tissue=key[0],
        gene_id=key[1],
        contrast=key[2],
        up_isoform=u.transcript_id,
        down_isoform=v.transcript_id,
        up_delta=u.delta_PSI,
        down_delta=v.delta_PSI,
        gene_q=u.gene_q,
        up_q=u.q_global,
        down_q=v.q_global,
        symbol=a.symbol,
        first_exon_changed=uf != vf,
        last_exon_changed=ul != vl,
        junction_chain_changed=uj != vj,
        internal_exon_changed=ue[1:-1] != ve[1:-1],
        same_terminal_exons=uf == vf and ul == vl,
    )
    for side, x in [("up", a), ("down", b)]:
        for field in [
            "structural_category",
            "subcategory",
            "ORF_aa",
            "ORF_unambiguous",
            "ORF_types",
            "predicted_NMD",
            "UTR5_nt",
            "UTR3_nt",
            "domain_assessable",
            "domains",
            "all_canonical",
            "RTS_stage",
            "latest_protein_agrees",
        ]:
            rec[f"{side}_{field}"] = x[field]
    rec["structure_strict"] = all(
        x.all_canonical == "canonical"
        and str(x.RTS_stage).upper() == "FALSE"
        and x.structural_category
        in ["full-splice_match", "novel_in_catalog", "novel_not_in_catalog"]
        for x in [a, b]
    )
    rec["ORF_assessable"] = rec["structure_strict"] and all(
        true(x.ORF_unambiguous)
        and true(x.latest_protein_agrees)
        and x.ORF_types == "complete"
        and np.isfinite(x.ORF_aa)
        for x in [a, b]
    )
    rec["delta_ORF_aa"] = a.ORF_aa - b.ORF_aa if rec["ORF_assessable"] else np.nan
    rec["NMD_assessable"] = rec["ORF_assessable"] and all(
        str(x.predicted_NMD).upper() in ["TRUE", "FALSE"] for x in [a, b]
    )
    rec["NMD_gain"] = rec["NMD_assessable"] and true(a.predicted_NMD) and not true(b.predicted_NMD)
    rec["NMD_loss"] = rec["NMD_assessable"] and not true(a.predicted_NMD) and true(b.predicted_NMD)
    rec["domain_pair_assessable"] = (
        rec["ORF_assessable"] and true(a.domain_assessable) and true(b.domain_assessable)
    )
    rec["domain_changed"] = rec["domain_pair_assessable"] and str(a.domains) != str(b.domains)
    for tag, col in [("UTR5", "UTR5_nt"), ("UTR3", "UTR3_nt")]:
        rec[tag + "_assessable"] = (
            rec["ORF_assessable"] and np.isfinite(a[col]) and np.isfinite(b[col])
        )
        rec["delta_" + tag + "_nt"] = a[col] - b[col] if rec[tag + "_assessable"] else np.nan
    rows.append(rec)
p = pd.DataFrame(rows)
save(p, "PSI_directional_structure_pairs.tsv.gz")
fun = []
for con, d in p.groupby("contrast"):
    specs = [
        ("Longer predicted ORF", "ORF_assessable", d.delta_ORF_aa.gt(0)),
        ("Shorter predicted ORF", "ORF_assessable", d.delta_ORF_aa.lt(0)),
        ("Predicted NMD gain", "NMD_assessable", d.NMD_gain),
        ("Predicted NMD loss", "NMD_assessable", d.NMD_loss),
        ("Domain-set difference", "domain_pair_assessable", d.domain_changed),
        ("First-exon difference", "structure_strict", d.first_exon_changed),
        ("Last-exon difference", "structure_strict", d.last_exon_changed),
        ("Internal-exon difference", "structure_strict", d.internal_exon_changed),
        ("Junction-chain difference", "structure_strict", d.junction_chain_changed),
        ("5-prime UTR difference", "UTR5_assessable", d.delta_UTR5_nt.ne(0)),
        ("3-prime UTR difference", "UTR3_assessable", d.delta_UTR3_nt.ne(0)),
    ]
    for name, flag, positive in specs:
        n = int(d[flag].sum())
        k = int((d[flag] & positive).sum())
        fun.append(
            dict(
                contrast=con,
                consequence=name,
                all_pairs=len(d),
                assessable=n,
                changed=k,
                fraction_changed=k / n if n else np.nan,
                assessable_fraction=n / len(d),
            )
        )
save(pd.DataFrame(fun), "PSI_functional_consequence_summary.tsv")
static = m[m.cluster.isin([2, 3, 5, 6, 7, 8])][
    ["tissue", "transcript_id", "cluster", "significant"]
].merge(
    ann.reset_index()[["transcript_id", "ORF_aa", "predicted_NMD", "UTR3_nt"]], on="transcript_id"
)
static = static.merge(mp.reset_index()[["transcript_id", "exon_count"]], on="transcript_id")
save(static, "legacy_cluster_static_attributes.tsv.gz")
print("FUNCTIONAL COMPLETE", len(e), "GO tests", len(p), "directional pairs", flush=True)
