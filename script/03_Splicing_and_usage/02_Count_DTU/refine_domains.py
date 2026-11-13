#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 03_Splicing_and_usage
# Purpose: Refine domains.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

from pathlib import Path
import pandas as pd, collections, json

root = Path(__file__).resolve().parents[1]
base = root.parent
ann = pd.read_csv(root / "inputs/canonical_annotations.tsv.gz", sep="\t", low_memory=False)
sq = pd.read_csv(
    base
    / "10.5_20240903_SQANTI/Swiss/9.4_combined.gtf_rmpoint_gffcompare_Only_Swiss_CDS_classification.txt",
    sep="\t",
    usecols=["isoform", "ORF_seq"],
    low_memory=False,
).set_index("isoform")
targets = set(s for r in ann.itertuples() for s in r.source_ids.split(","))
seqs = {}
current = None
seq = []
path = (
    base
    / "9_20240903.TransDecoder/OUT/9.4_combined.gtf_rmpoint_gffcompare.fasta.transdecoder_dir/longest_orfs.pep"
)
with path.open() as f:
    for line in f:
        if line.startswith(">"):
            if current in targets:
                seqs[current] = "".join(seq)
            current = line[1:].split()[0]
            seq = []
        elif current in targets:
            seq.append(line.strip())
if current in targets:
    seqs[current] = "".join(seq)
raw = pd.read_csv(root / "inputs/pfam_raw_hits.tsv.gz", sep="\t")
hits = collections.defaultdict(list)
for r in raw.itertuples():
    hits[r.source_id].append(r)
provenance = []
assess = []
domainlist = []
partiallist = []
for r in ann.itertuples():
    ids = r.source_ids.split(",")
    proteins = set(
        str(sq.loc[s, "ORF_seq"]).rstrip("*")
        for s in ids
        if s in sq.index and pd.notna(sq.loc[s, "ORF_seq"])
    )
    mappings = []
    if len(proteins) == 1:
        protein = next(iter(proteins))
        for s in ids:
            old = seqs.get(s, "")
            if (
                protein
                and old.count(protein) == 1
                and all(h.target_length == len(old) for h in hits[s])
            ):
                offset = old.index(protein)
                retained = []
                partial = []
                for h in hits[s]:
                    if h.ali_start >= offset + 1 and h.ali_end <= offset + len(protein):
                        retained.append(h.pfam)
                    elif h.ali_end > offset and h.ali_start <= offset + len(protein):
                        partial.append(h.pfam)
                mappings.append(
                    (
                        s,
                        offset,
                        len(protein),
                        tuple(sorted(set(retained))),
                        tuple(sorted(set(partial))),
                    )
                )
    consensus = set((a[3], a[4]) for a in mappings)
    ok = bool(mappings and len(consensus) == 1)
    assess.append(ok)
    domainlist.append(";".join(mappings[0][3]) if ok else None)
    partiallist.append(";".join(mappings[0][4]) if ok else None)
    for s, offset, n, ds, ps in mappings:
        provenance.append(
            dict(
                transcript_id=r.transcript_id,
                source_id=s,
                searched_peptide_length=len(seqs[s]),
                final_ORF_length=n,
                offset_in_searched_peptide=offset,
                retained_domains=";".join(ds),
                partial_domains=";".join(ps),
                consensus_assessable=ok,
            )
        )
ann["domain_assessable"] = assess
ann["domains"] = domainlist
ann["partial_domains"] = partiallist
ann.to_csv(root / "inputs/canonical_annotations.tsv.gz", sep="\t", index=False)
pd.DataFrame(provenance).to_csv(
    root / "inputs/domain_coordinate_provenance.tsv.gz", sep="\t", index=False
)
print("domain assessable", sum(assess), "/", len(assess), flush=True)
