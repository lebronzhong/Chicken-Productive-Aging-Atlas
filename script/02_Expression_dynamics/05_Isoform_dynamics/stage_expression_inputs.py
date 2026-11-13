#!/usr/bin/env python3
"""Stage declared input tables without importing historical analysis directories."""
import argparse
import csv
import shutil
from pathlib import Path


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", required=True, help="TSV columns: source, destination")
    ap.add_argument("--output-dir", required=True)
    args = ap.parse_args()
    root = Path(args.output_dir).resolve()
    manifest = Path(args.manifest).resolve()
    rows = list(csv.DictReader(manifest.open(encoding="utf-8"), delimiter="\t"))
    if not rows or not {"source", "destination"}.issubset(rows[0]):
        ap.error("Manifest must contain source and destination columns and at least one row.")
    plan = []
    for row in rows:
        src = Path(row["source"])
        if not src.is_absolute():
            src = manifest.parent / src
        dst = (root / row["destination"]).resolve()
        if not src.is_file() or root not in dst.parents:
            ap.error("Input does not exist or destination escapes the output directory.")
        if dst.exists():
            ap.error(f"Preserving existing file: {dst}")
        plan.append((src, dst))
    if len({dst for _, dst in plan}) != len(plan):
        ap.error("Duplicate destination.")
    for src, dst in plan:
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
    print(f"Staged {len(plan)} input tables.")


if __name__ == "__main__":
    main()
