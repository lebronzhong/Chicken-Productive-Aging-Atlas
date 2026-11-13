#!/usr/bin/env python3
import argparse
import gzip
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="Split TensorQTL phenotype BED.gz by chromosome.")
    parser.add_argument("--input-bed-gz", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--prefix", required=True)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    handles = {}
    header = None

    with gzip.open(args.input_bed_gz, "rt") as f:
        for line in f:
            if header is None:
                header = line
                continue
            chrom = line.split("\t", 1)[0]
            if chrom not in handles:
                out = gzip.open(outdir / f"{args.prefix}.chr{chrom}.bed.gz", "wt")
                out.write(header)
                handles[chrom] = out
            handles[chrom].write(line)

    for h in handles.values():
        h.close()

    print(f"Split into {len(handles)} chromosome files in {outdir}")


if __name__ == "__main__":
    main()
