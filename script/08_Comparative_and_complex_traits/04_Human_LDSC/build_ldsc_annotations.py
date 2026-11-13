#!/usr/bin/env python3
"""Audit recorded upstream assets and generate LDSC annotations from human gene sets."""
import argparse, csv, json, re, subprocess, sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--ldsc-dir", required=True, type=Path)
    parser.add_argument("--gene-coordinates", required=True, type=Path)
    parser.add_argument("--bim-pattern", required=True, help="Input BIM path with {chromosome}")
    parser.add_argument("--chromosomes", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    if "{chromosome}" not in args.bim_pattern:
        parser.error("BIM pattern requires {chromosome}")
    with args.manifest.open(encoding="utf8") as h:
        rows = list(csv.DictReader(h, delimiter="\t"))
    if not rows or len({r["annotation_id"] for r in rows}) != len(rows):
        parser.error("Annotation IDs must be nonempty and unique")
    with args.chromosomes.open(encoding="utf8") as h:
        chromosomes = [r["chromosome"] for r in csv.DictReader(h, delimiter="\t")]
    if not chromosomes or len(chromosomes) != len(set(chromosomes)):
        parser.error("Invalid chromosome list")
    for row in rows:
        label = row["annotation_id"]
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", label):
            parser.error("Invalid annotation ID")
        inputs = {}
        for key in ["t_statistic", "top_transcripts", "ortholog_table", "human_gene_set"]:
            p = Path(row[key])
            p = p if p.is_absolute() else args.manifest.resolve().parent / p
            if not args.dry_run and not p.is_file():
                raise FileNotFoundError(p)
            inputs[key] = p
        for chrom in chromosomes:
            if not re.fullmatch(r"[A-Za-z0-9_]+", chrom):
                parser.error("Invalid chromosome")
            command = [
                sys.executable,
                str(args.ldsc_dir / "make_annot.py"),
                "--gene-set-file",
                str(inputs["human_gene_set"]),
                "--gene-coord-file",
                str(args.gene_coordinates),
                "--windowsize",
                "100000",
                "--bimfile",
                args.bim_pattern.format(chromosome=chrom),
                "--annot-file",
                str(args.output_dir / f"{label}.{chrom}.annot.gz"),
            ]
            if args.dry_run:
                print(json.dumps({"annotation": label, "command": command}))
            else:
                args.output_dir.mkdir(parents=True, exist_ok=True)
                subprocess.run(command, check=True)


if __name__ == "__main__":
    main()
