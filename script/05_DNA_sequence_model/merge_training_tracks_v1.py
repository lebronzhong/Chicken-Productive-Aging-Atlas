#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 05_DNA_sequence_model
# Purpose: Merge training tracks v1.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

import json, pathlib, os, sys, concurrent.futures, time
import numpy as np
import pyBigWig

ROOT = pathlib.Path(__file__).resolve().parent
OUT = ROOT / "training_inputs_v1"
inputs = json.load(open(OUT / "merge_inputs.json"))
chroms = []
for line in list(open(OUT / "chromosome_splits.tsv"))[1:]:
    x = line.split()
    chroms.append((x[0], int(x[1])))


def merge(tissue):
    start = time.time()
    safe = tissue.replace(" ", "_")
    folder = OUT / "tissue_tracks" / safe
    folder.mkdir(parents=True, exist_ok=False)
    target = folder / "coverage.bw.partial"
    handles = [pyBigWig.open(x["bigwig"]) for x in inputs[tissue]]
    bw = pyBigWig.open(str(target), "w")
    bw.addHeader(chroms)
    checks = []
    try:
        for chrom, n in chroms:
            for lo in range(0, n, 500000):
                hi = min(n, lo + 500000)
                values = np.zeros(hi - lo, dtype=np.float64)
                for h, item in zip(handles, inputs[tissue]):
                    raw = h.values(chrom, lo, hi, numpy=True)
                    assert not np.isinf(raw).any()
                    raw = np.nan_to_num(raw, nan=0)
                    assert (raw >= 0).all()
                    values += raw.astype(np.float64) * item["scale"] / 25
                values = values.astype(np.float32)
                assert np.isfinite(values).all()
                bounds = np.r_[0, np.flatnonzero(np.diff(values)) + 1, len(values)]
                ss = bounds[:-1]
                ee = bounds[1:]
                v = values[ss]
                nz = v > 0
                if np.any(nz):
                    starts = (ss[nz] + lo).tolist()
                    bw.addEntries(
                        [chrom] * len(starts),
                        starts,
                        ends=(ee[nz] + lo).tolist(),
                        values=v[nz].astype(float).tolist(),
                    )
                if lo == 0:
                    # Check an actual nonzero segment if present; otherwise a zero segment.
                    nzpos = np.flatnonzero(values > 0)
                    p = int(nzpos[0]) if len(nzpos) else 0
                    checks.append((chrom, lo + p, lo + min(p + 128, len(values))))
        bw.close()
        bw = None
        with pyBigWig.open(str(target)) as check:
            assert list(check.chroms().items()) == chroms
            for chrom, a, b in checks:
                expected = np.zeros(b - a)
                for h, item in zip(handles, inputs[tissue]):
                    expected += (
                        np.nan_to_num(h.values(chrom, a, b, numpy=True), nan=0).astype(float)
                        * item["scale"]
                        / 25
                    )
                actual = np.nan_to_num(check.values(chrom, a, b, numpy=True), nan=0)
                assert np.allclose(expected, actual, rtol=2e-6, atol=1e-6)
    finally:
        if bw is not None:
            bw.close()
        for h in handles:
            h.close()
    target.rename(folder / "coverage.bw")
    meta = {
        "tissue": tissue,
        "training_animals": 25,
        "checks": len(checks),
        "bytes": (folder / "coverage.bw").stat().st_size,
        "elapsed_seconds": time.time() - start,
        "normalization": "Mean of per-animal depth scaled to10million autosomal primary reads",
        "source_inputs": inputs[tissue],
    }
    (folder / "metadata.json").write_text(json.dumps(meta, indent=2))
    print(tissue, "COMPLETE", meta["bytes"], round(meta["elapsed_seconds"]), flush=True)
    return meta


with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
    results = list(pool.map(merge, inputs))
(OUT / "merge_summary.json").write_text(json.dumps(results, indent=2))
print("COMPLETE:33 tissue labels; all25 animals; chromosome holdout", flush=True)
