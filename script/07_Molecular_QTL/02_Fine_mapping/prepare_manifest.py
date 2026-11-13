#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 07_Molecular_QTL
# Purpose: Prepare manifest.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

"""Read-only input checks; write a frozen manifest in the fine-mapping directory."""

import argparse, csv, gzip, hashlib, json
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--manifest", required=True, type=Path)
parser.add_argument("--output", required=True, type=Path)
parser.add_argument("--expected-chromosomes", type=int, default=32)
args = parser.parse_args()
if args.expected_chromosomes < 1:
    parser.error("--expected-chromosomes must be positive")
source = args.manifest.resolve()
with source.open(encoding="utf8") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    required_columns = {"tissue", "qtl_type", "fdr", "genotype", "covariates", "phenotype_dir"}
    if not required_columns <= set(reader.fieldnames or []):
        parser.error("Manifest lacks required input columns")
    rows = list(reader)
if not rows or len({(x["tissue"], x["qtl_type"]) for x in rows}) != len(rows):
    parser.error("Manifest must contain unique tissue/QTL records")


def resolve(value):
    value = Path(value)
    return value if value.is_absolute() else source.parent / value


tasks = []
for row in rows:
    t, kind = row["tissue"], row["qtl_type"]
    if kind not in {"teQTL", "tuQTL"}:
        parser.error("qtl_type must use the recorded file labels teQTL or tuQTL")
    fdr, geno, cov, beds = (
        resolve(row[key]) for key in ["fdr", "genotype", "covariates", "phenotype_dir"]
    )
    files = sorted(beds.glob("*.bed.gz"))
    assert len(files) == args.expected_chromosomes, (
        t,
        kind,
        "unexpected chromosome count",
        len(files),
    )
    required = [fdr, cov] + [Path(str(geno) + "." + e) for e in ("bed", "bim", "fam")]
    for p in required:
        assert p.is_file() and p.stat().st_size > 0, str(p)
    fam = [r.split()[1] for r in Path(str(geno) + ".fam").open()]
    c = cov.open().readline().rstrip("\r\n").split("\t")[1:]
    assert len(c) == len(set(c)) and len(fam) == len(set(fam))
    assert set(c) <= set(fam)
    for p in files:
        with gzip.open(p, "rt") as h:
            a = h.readline().rstrip("\r\n").split("\t")[4:]
        assert a == c, ("sample order", str(p))
    with fdr.open() as h:
        rd = csv.DictReader(h, delimiter="\t")
        assert {"phenotype_id", "qval"} <= set(rd.fieldnames)
        ids = [r["phenotype_id"] for r in rd if float(r["qval"]) <= 0.05]
    assert len(ids) == len(set(ids)) and ids, (t, kind, "duplicate or empty FDR ids")
    task = dict(
        tissue=t,
        qtl_type=kind,
        genotype=str(geno),
        covariates=str(cov),
        fdr=str(fdr),
        phenotype_files=[str(p) for p in files],
        samples=len(c),
        significant=len(ids),
        fdr_sha256=hashlib.sha256(fdr.read_bytes()).hexdigest(),
        cov_sha256=hashlib.sha256(cov.read_bytes()).hexdigest(),
    )
    task["files"] = {
        str(p): dict(size=p.stat().st_size, mtime_ns=p.stat().st_mtime_ns) for p in required + files
    }
    tasks.append(task)
    print(
        t,
        kind,
        "samples",
        len(c),
        "significant",
        len(ids),
        "chromosomes",
        len(files),
        flush=True,
    )
out = args.output.resolve()
payload = json.dumps(tasks, ensure_ascii=False, indent=2) + "\n"
if out.exists():
    assert (
        out.read_text(encoding="utf8") == payload
    ), "Existing manifest differs; do not replace during a run."
else:
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(payload, encoding="utf8")
print(
    f"PASS: {len(tasks)} tasks; {sum(len(t['phenotype_files']) for t in tasks)} phenotype headers checked",
    flush=True,
)
