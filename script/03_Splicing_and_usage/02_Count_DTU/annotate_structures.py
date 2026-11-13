#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 03_Splicing_and_usage
# Purpose: Annotate structures.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

from pathlib import Path
import pandas as pd, numpy as np, re, json, collections

root = Path(__file__).resolve().parents[1]
base = root.parent
mp = pd.read_csv(root / "inputs/canonical_transcript_mapping.tsv.gz", sep="\t", low_memory=False)
srcmap = {s: r.transcript_id for r in mp.itertuples() for s in r.source_ids.split(",")}
gtf = base / "10_20240903.gffcompare/9.6_final_gtf/9.9_finalgtf_annotated_rename.gtf"
exons = collections.defaultdict(list)
parents = {}
coords = {}
with gtf.open() as f:
    for line in f:
        if line.startswith("#"):
            continue
        v = line.rstrip("\n").split("\t")
        if len(v) < 9 or v[2] != "exon":
            continue
        a = dict(re.findall(r'(\S+) "([^"]+)"', v[8]))
        s = a.get("transcript_id")
        if s not in srcmap:
            continue
        exons[s].append((int(v[3]), int(v[4])))
        parents[s] = a["gene_id"]
        coords[s] = (v[0], v[6])
conf = []
for r in mp.itertuples():
    expected = sorted(tuple(map(int, x.split("-"))) for x in r.exons.split(";"))
    for s in r.source_ids.split(","):
        if (
            sorted(exons.get(s, [])) != expected
            or parents.get(s) != r.gene_id
            or coords.get(s) != (str(r.chromosome), r.strand)
        ):
            conf.append(
                dict(
                    transcript_id=r.transcript_id,
                    source_id=s,
                    gtf_gene=parents.get(s),
                    canonical_gene=r.gene_id,
                )
            )
pd.DataFrame(conf).to_csv(root / "audit/gtf_mapping_conflicts.tsv", sep="\t", index=False)
assert len(conf) == 0, ("GTF conflicts", len(conf))
print("GTF structure and parent validation passed", len(mp), flush=True)
sqpath = (
    base
    / "10.5_20240903_SQANTI/Swiss/9.4_combined.gtf_rmpoint_gffcompare_Only_Swiss_CDS_classification.txt"
)
sq = pd.read_csv(sqpath, sep="\t", low_memory=False).set_index("isoform")
latest = (
    base
    / "10.5_20260314_SQANTI/9.4_combined.gtf_rmpoint_gffcompare_Only_Swiss_CDS_corrected.clean_classification.txt"
)
sql = pd.read_csv(latest, sep="\t", low_memory=False).set_index("isoform")
pep = (
    base
    / "9_20240903.TransDecoder/OUT_Swiss/9.4_combined.gtf_rmpoint_gffcompare.fasta.transdecoder.pep"
)
orfs = {}
seqs = {}
s = None
with pep.open() as f:
    for line in f:
        if line.startswith(">"):
            s = line[1:].split()[0]
            m = re.search(r"ORF type:(\S+)", line)
            orfs[s] = m.group(1) if m else "unknown"
            seqs[s] = ""
        elif s:
            seqs[s] += line.strip().rstrip("*")
domains = collections.defaultdict(set)
domain_rows = []
with (base / "9_20240903.TransDecoder/pfam.domtblout").open() as f:
    for line in f:
        if line.startswith("#"):
            continue
        v = line.split()
        if len(v) < 22 or float(v[12]) > 0.001:
            continue
        # Domains belong to the actual HMM target protein, which may precede ORF refinement.
        domain_rows.append(
            dict(
                source_id=v[0],
                target_length=int(v[2]),
                domain=v[3],
                pfam=v[4],
                i_evalue=float(v[12]),
                ali_start=int(v[17]),
                ali_end=int(v[18]),
            )
        )
domains_df = pd.DataFrame(domain_rows)
domains_df.to_csv(root / "inputs/pfam_raw_hits.tsv.gz", sep="\t", index=False)
domain_target_lengths = collections.defaultdict(set)
for dr in domains_df.itertuples():
    domain_target_lengths[dr.source_id].add(dr.target_length)
for r in domains_df.itertuples():
    if r.source_id in seqs and len(seqs[r.source_id]) == r.target_length:
        domains[r.source_id].add(r.pfam)
fields = [
    "structural_category",
    "subcategory",
    "all_canonical",
    "RTS_stage",
    "predicted_NMD",
    "length",
    "exons",
    "CDS_start",
    "CDS_end",
    "CDS_genomic_start",
    "CDS_genomic_end",
    "diff_to_TSS",
    "diff_to_TTS",
    "within_CAGE_peak",
    "within_polyA_site",
    "polyA_motif_found",
    "FL",
]
rows = []
gene_symbols = {}
with (
    Path("/path/to/chicken-aging/reference/galGal6/Ensemble/Gallus_gallus.GRCg6a.106.gtf")
).open() as f:
    for line in f:
        if line.startswith("#"):
            continue
        v = line.split("\t")
        if len(v) < 9 or v[2] != "gene":
            continue
        a = dict(re.findall(r'(\S+) "([^"]+)"', v[8]))
        gene_symbols[a["gene_id"]] = a.get("gene_name", a["gene_id"])
for r in mp.itertuples():
    ids = r.source_ids.split(",")
    present = [s for s in ids if s in sq.index]
    q = sq.loc[present]
    rec = dict(
        transcript_id=r.transcript_id,
        gene_id=r.gene_id,
        symbol=gene_symbols.get(r.gene_id, r.gene_id),
        annotation_sources=len(present),
        source_ids=r.source_ids,
    )
    for field in fields:
        values = sorted(set(str(x) for x in q[field].dropna())) if field in q else []
        rec[field] = ";".join(values) if values else None
    proteins = set(str(x).rstrip("*") for x in q.ORF_seq.dropna()) if "ORF_seq" in q else set()
    rec["n_predicted_proteins"] = len(proteins)
    rec["ORF_unambiguous"] = len(proteins) == 1
    rec["ORF_aa"] = len(next(iter(proteins))) if len(proteins) == 1 else None
    rec["ORF_types"] = ";".join(sorted(set(orfs[s] for s in ids if s in orfs)))
    rec["transdecoder_aa"] = ";".join(
        str(x) for x in sorted(set(len(seqs[s]) for s in ids if s in seqs))
    )
    rec["transdecoder_matches_SQANTI"] = len(proteins) == 1 and all(
        seqs.get(s) == next(iter(proteins)) for s in ids
    )
    # Only label domains assessable if the searched sequence matches final SQANTI sequence exactly.
    matching = [
        s
        for s in ids
        if s in seqs
        and len(proteins) == 1
        and seqs[s] == next(iter(proteins))
        and len(seqs[s]) in domain_target_lengths[s]
    ]
    rec["domain_assessable"] = bool(matching)
    rec["domains"] = (
        ";".join(sorted(set().union(*(domains[s] for s in matching)))) if matching else None
    )

    def num(k):
        try:
            return float(rec[k])
        except:
            return np.nan

    rec["UTR5_nt"] = num("CDS_start") - 1
    rec["UTR3_nt"] = num("length") - num("CDS_end")
    ql = sql.loc[[s for s in ids if s in sql.index]]
    lp = set(str(x).rstrip("*") for x in ql.ORF_seq.dropna()) if "ORF_seq" in ql else set()
    rec["latest_protein_agrees"] = len(proteins) == 1 and lp == proteins
    rows.append(rec)
result = pd.DataFrame(rows)
result.to_csv(root / "inputs/canonical_annotations.tsv.gz", sep="\t", index=False)
pd.DataFrame([{"gene_id": k, "symbol": v} for k, v in gene_symbols.items()]).to_csv(
    root / "inputs/gene_symbols.tsv", sep="\t", index=False
)
print("Annotated", len(result), "unambiguous proteins", result.ORF_unambiguous.sum(), flush=True)
