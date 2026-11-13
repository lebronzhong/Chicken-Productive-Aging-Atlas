#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 07_Molecular_QTL
# Purpose: Extract candidate qtl.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

from pathlib import Path
import gzip
import csv
from collections import defaultdict

import argparse

parser = argparse.ArgumentParser(
    description="Extract transcript QTL records for a candidate manifest."
)
parser.add_argument("--manifest", required=True)
parser.add_argument("--input-dir", required=True)
parser.add_argument("--output-dir", required=True)
args = parser.parse_args()
INPUT = Path(args.input_dir)
OUT = Path(args.output_dir)
TU_OUT = OUT / "01_tuQTL"
with open(args.manifest, encoding="utf-8") as handle:
    CANDIDATES = list(csv.DictReader(handle, delimiter="\t"))
required = {
    "candidate_id",
    "tissue",
    "transcript_id",
    "parent_gene_id",
    "trait_num",
    "trait_name",
    "N_qtl",
}
if not CANDIDATES or not required.issubset(CANDIDATES[0]):
    parser.error("Candidate manifest is empty or missing required columns.")
if len({r["candidate_id"] for r in CANDIDATES}) != len(CANDIDATES):
    parser.error("Candidate IDs must be unique.")
for row in CANDIDATES:
    for key in ["candidate_id", "tissue", "transcript_id"]:
        if not row[key] or "/" in row[key] or "\\" in row[key] or row[key] in {".", ".."}:
            parser.error(f"Invalid path component in {key}.")
    try:
        n = float(row["N_qtl"])
        if not (n > 0 and n.is_integer()):
            raise ValueError()
    except ValueError:
        parser.error("N_qtl must be a positive integer.")
TU_OUT.mkdir(parents=True, exist_ok=True)
# ------------------------------------------------------------
# Manifest
# ------------------------------------------------------------
manifest = OUT / "candidate_manifest.tsv"
with open(manifest, "w", newline="") as f:
    writer = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=[
            "candidate_id",
            "tissue",
            "transcript_id",
            "parent_gene_id",
            "trait_num",
            "trait_name",
            "N_qtl",
        ],
    )
    writer.writeheader()
    writer.writerows({k: r[k] for k in writer.fieldnames} for r in CANDIDATES)
# ------------------------------------------------------------
# Group by tissue so each huge nominal file is scanned ONCE
# ------------------------------------------------------------
by_tissue = defaultdict(list)
for c in CANDIDATES:
    by_tissue[c["tissue"]].append(c)
summary = []
for tissue, candidates in by_tissue.items():
    infile = INPUT / tissue / f"{tissue}.tuQTL.sig.nominal.tsv.gz"
    if not infile.exists():
        raise RuntimeError(f"Missing nominal file:\n{infile}")
    print()
    print("=" * 90)
    print(tissue)
    print("=" * 90)
    print("Input:", infile)
    tx_to_candidates = defaultdict(list)
    for c in candidates:
        tx_to_candidates[c["transcript_id"]].append(c)
    handles = {}
    writers = {}
    counts = defaultdict(int)
    # Read header first
    with gzip.open(infile, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        fields = reader.fieldnames
        required = {
            "phenotype_id",
            "variant_id",
            "distance",
            "af",
            "pval_nominal",
            "slope",
            "slope_se",
        }
        missing = required - set(fields)
        if missing:
            raise RuntimeError(f"{tissue}: missing columns {missing}")
        # open candidate outputs
        for c in candidates:
            outfile = (
                TU_OUT / f"{c['candidate_id']}."
                f"{c['tissue']}."
                f"{c['transcript_id']}.tuQTL.tsv.gz"
            )
            h = gzip.open(outfile, "wt", newline="")
            w = csv.DictWriter(h, delimiter="\t", fieldnames=fields)
            w.writeheader()
            handles[c["candidate_id"]] = h
            writers[c["candidate_id"]] = w
        scanned = 0
        for row in reader:
            scanned += 1
            tx = row["phenotype_id"]
            if tx not in tx_to_candidates:
                continue
            for c in tx_to_candidates[tx]:
                cid = c["candidate_id"]
                writers[cid].writerow(row)
                counts[cid] += 1
            if scanned % 5000000 == 0:
                print(f"Scanned {scanned:,} rows", flush=True)
    for h in handles.values():
        h.close()
    print(f"Total scanned: {scanned:,}")
    for c in candidates:
        cid = c["candidate_id"]
        n = counts[cid]
        print(f"{cid} " f"{c['transcript_id']}: " f"{n:,} rows")
        if n == 0:
            raise RuntimeError(f"{cid}: no tuQTL rows extracted")
        summary.append(
            {
                "candidate_id": cid,
                "tissue": tissue,
                "transcript_id": c["transcript_id"],
                "n_tuQTL_rows": n,
            }
        )
summary_file = OUT / "candidate_extraction_summary.tsv"
with open(summary_file, "w", newline="") as f:
    writer = csv.DictWriter(
        f, delimiter="\t", fieldnames=["candidate_id", "tissue", "transcript_id", "n_tuQTL_rows"]
    )
    writer.writeheader()
    writer.writerows(summary)
print()
print("=" * 90)
print("Candidate extraction completed")
print("=" * 90)
print("Output:", OUT)
