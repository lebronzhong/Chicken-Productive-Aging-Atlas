#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 05_DNA_sequence_model
# Purpose: Make unstranded v1.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

import argparse, collections, csv, json, pathlib, time
import numpy as np
import pysam, pyBigWig

ROOT = pathlib.Path(__file__).resolve().parent
rows = list(
    csv.DictReader(open(ROOT / "samples.matched.tsv", encoding="utf-8-sig"), delimiter="\t")
)
p = argparse.ArgumentParser()
p.add_argument("--index", type=int, required=True)
args = p.parse_args()
s = rows[args.index]
final = ROOT / "tracks_unstranded_v1" / s["ont_sample_id"]
assert not final.exists(), "Completed output exists"
tmp = final.with_name(final.name + ".partial")
tmp.mkdir(parents=True, exist_ok=False)
t = time.time()
counts = collections.Counter()
chromstats = []


def keep(r):
    return not r.flag & (4 | 256 | 512 | 2048) and r.mapping_quality >= 20


with pysam.AlignmentFile(s["bam_path"], "rb", threads=2) as bam:
    header = list(zip(bam.references, bam.lengths))
    bw = pyBigWig.open(str(tmp / "coverage.bw"), "w")
    bw.addHeader(header)
    try:
        for chrom, n in header:
            delta = np.zeros(n + 1, dtype=np.int32)
            expected = 0
            reads = 0
            for r in bam.fetch(chrom):
                counts["mapped_records"] += 1
                if not keep(r):
                    continue
                reads += 1
                counts["retained_primary_reads"] += 1
                if r.is_duplicate:
                    counts["duplicate_flag_retained"] += 1
                for a, b in r.get_blocks():
                    assert 0 <= a <= b <= n
                    delta[a] += 1
                    delta[b] -= 1
                    expected += b - a
            cov = np.cumsum(delta[:-1], dtype=np.int64)
            assert int(delta.sum()) == 0 and int(cov.sum()) == expected and cov.min(initial=0) >= 0
            for offset in range(0, n, 1000000):
                a = cov[offset : offset + 1000000]
                bounds = np.r_[0, np.flatnonzero(np.diff(a)) + 1, len(a)]
                starts = bounds[:-1]
                ends = bounds[1:]
                values = a[starts]
                nz = values > 0
                if np.any(nz):
                    ss = (starts[nz] + offset).tolist()
                    bw.addEntries(
                        [chrom] * len(ss),
                        ss,
                        ends=(ends[nz] + offset).tolist(),
                        values=values[nz].astype(float).tolist(),
                    )
            del cov, delta
            chromstats.append({"chrom": chrom, "reads": reads, "aligned_bases": expected})
            if reads:
                print(s["ont_sample_id"], chrom, reads, flush=True)
    finally:
        bw.close()
assert counts["retained_primary_reads"] > 0
# Independent region checks use CIGAR directly rather than get_blocks.
checked = 0
with pysam.AlignmentFile(s["bam_path"], "rb") as bam, pyBigWig.open(str(tmp / "coverage.bw")) as bw:
    assert list(bw.chroms().items()) == header
    for chrom, n in header[:10]:
        anchor = None
        for j, r in enumerate(bam.fetch(chrom)):
            if keep(r):
                anchor = r.reference_start
                break
            if j >= 9999:
                break
        if anchor is None:
            continue
        lo = anchor
        hi = min(n, lo + 2048)
        expected = np.zeros(hi - lo, dtype=np.int64)
        for r in bam.fetch(chrom, lo, hi):
            if not keep(r):
                continue
            cursor = r.reference_start
            for op, length in r.cigartuples or []:
                if op in (0, 7, 8):
                    a = max(lo, cursor)
                    b = min(hi, cursor + length)
                    if a < b:
                        expected[a - lo : b - lo] += 1
                if op in (0, 2, 3, 7, 8):
                    cursor += length
        observed = np.nan_to_num(bw.values(chrom, lo, hi, numpy=True), nan=0)
        assert np.array_equal(observed, expected), (chrom, lo, hi)
        checked += 1
assert checked > 0
meta = {
    "sample": s,
    "stats": dict(counts),
    "chromosomes": chromstats,
    "elapsed_seconds": time.time() - t,
    "bigwig_bytes": (tmp / "coverage.bw").stat().st_size,
    "independent_regions_passed": checked,
    "mode": "unstranded_raw_depth",
    "filter": "Primary, mapped, not QC-fail, MAPQ>=20; duplicate flag retained; CIGAR N/D excluded.",
    "normalization_denominator": counts["retained_primary_reads"],
    "biological_validation": "Read orientation is not used; source mapping sensitivity remains a separate QC.",
}
(tmp / "metadata.json").write_text(json.dumps(meta, indent=2))
tmp.rename(final)
print("COMPLETE", s["ont_sample_id"], checked, flush=True)
