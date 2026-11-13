#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 05_DNA_sequence_model
# Purpose: Prepare test masks 33.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

import pathlib, json, hashlib, collections
import numpy as np, pandas as pd, pyBigWig

R = pathlib.Path("/path/to/chicken-aging/model")
D = R / "experiment33/data/tfrecords_v2"
O = R / "experiment33/test_masks_v1"
O.mkdir(exist_ok=True)
w = pd.read_csv(
    D / "sequences.bed",
    sep="\t",
    header=None,
    names=["chrom", "start", "end", "split"],
    dtype={"chrom": str},
)
w = w[w.split == "test"].reset_index(drop=True)
assert len(w) == 223 and set(w.chrom) == {"5", "11", "17"}
w.index.name = "window_index"
w.to_csv(O / "test_windows.tsv", sep="\t")
tissues = pd.read_csv(D / "targets.txt", sep="\t").identifier.tolist()


def merge(a):
    out = []
    for x, y in sorted(a):
        if out and x <= out[-1][1]:
            out[-1] = (out[-1][0], max(y, out[-1][1]))
        else:
            out.append((x, y))
    return out


def subtract(a, b):
    a = merge(a)
    for u, v in merge(b):
        z = []
        for x, y in a:
            if v <= x or u >= y:
                z.append((x, y))
                continue
            if x < u:
                z.append((x, u))
            if v < y:
                z.append((v, y))
        a = z
    return a


def bins(reg, start, full=True):
    return sorted(
        {
            k
            for a, b in reg
            for k in range(
                max(0, (a - start + (31 if full else 0)) // 32),
                min(12288, (b - start + (0 if full else 31)) // 32),
            )
        }
    )


ex = pd.read_csv(
    R / "experiment33/transfer_bundle_33t25_v2/reference/transcript_exons.tsv",
    sep="\t",
    dtype={"chrom": str},
)
ex = ex[ex.chrom.isin(set(w.chrom))]
bounds = (
    ex.groupby(["chrom", "gene_id"]).agg(lo=("start0", "min"), hi=("end0", "max")).reset_index()
)
genes = []
events = []
skip = collections.Counter()
for (chrom, gene), g in ex.groupby(["chrom", "gene_id"], sort=True):
    lo = int(g.start0.min())
    hi = int(g.end0.max())
    ww = w[(w.chrom == chrom) & (w.start <= lo) & (w.end >= hi)]
    if len(ww) != 1:
        skip["not_one_complete_window"] += 1
        continue
    overlap = bounds[
        (bounds.chrom == chrom) & (bounds.gene_id != gene) & (bounds.lo < hi) & (bounds.hi > lo)
    ]
    if len(overlap):
        skip["overlap_other_gene_span"] += 1
        continue
    if g.strand.nunique() != 1:
        skip["mixed_strand"] += 1
        continue
    wi = int(ww.index[0])
    start = int(ww.iloc[0].start)
    strand = g.strand.iloc[0]
    tx = {
        t: merge([(int(x.start0), int(x.end0)) for x in q.itertuples()])
        for t, q in g.groupby("transcript_id")
    }
    gb = bins([ab for v in tx.values() for ab in v], start, False)
    if len(gb) < 3:
        continue
    genes.append(dict(gene=gene, chrom=chrom, window_index=wi, bins=gb))
    unique = set(tuple(v) for v in tx.values())
    if len(unique) < 2:
        continue
    common = list(next(iter(unique)))
    for v in unique:
        common = merge(
            [(max(a, c), min(b, d)) for a, b in common for c, d in v if max(a, c) < min(b, d)]
        )
    cb = bins(common, start)
    if len(cb) < 3:
        continue
    for kind in ["TSS", "APA"]:
        low = (kind == "TSS") == (strand == "+")
        ends = {t: (v[0][0] if low else v[-1][1]) for t, v in tx.items()}
        extreme = min(ends.values()) if low else max(ends.values())
        group = [t for t in tx if ends[t] == extreme]
        other = [t for t in tx if t not in group]
        if not other:
            continue
        if kind == "APA":
            anchor = lambda t: tx[t][-1][0] if strand == "+" else tx[t][0][1]
            if not ({anchor(t) for t in group} & {anchor(t) for t in other}):
                continue
        edge = [tx[t][0] if low else tx[t][-1] for t in group]
        specific = subtract(edge, [ab for t in other for ab in tx[t]])
        sb = bins(specific, start)
        if len(sb) < 3:
            continue
        events.append(
            dict(
                event=gene + "__" + kind,
                gene=gene,
                type=kind,
                window_index=wi,
                specific=sb,
                common=cb,
                group_transcripts=group,
                other_transcripts=other,
            )
        )
assert genes and events
manifest = dict(
    genes=genes,
    events=events,
    tissues=tissues,
    n_windows=len(w),
    sequence_length=393216,
    bin_width=32,
    sequence_bed_sha256=hashlib.sha256((D / "sequences.bed").read_bytes()).hexdigest(),
    definition="Gene union exon-overlapping bins; TSS/APA extreme-end-specific/full common exon bins, identical to prior structural definitions. Expression eligibility fixed from test observations only; masks built from annotation without viewing predictions.",
)
(O / "test_masks.json").write_text(json.dumps(manifest))
print(
    json.dumps(
        dict(
            genes=len(genes),
            events=len(events),
            event_types=dict(collections.Counter(e["type"] for e in events)),
            n_windows=len(w),
            n_tissues=len(tissues),
        )
    ),
    flush=True,
)
