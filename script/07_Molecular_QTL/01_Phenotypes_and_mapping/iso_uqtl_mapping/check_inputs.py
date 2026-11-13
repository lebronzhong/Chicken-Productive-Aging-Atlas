#!/usr/bin/env python3
import argparse
import gzip
import os
from pathlib import Path


def open_text(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def read_bed_samples(path: str):
    with open_text(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
    if len(header) < 5:
        raise ValueError("Phenotype BED header has too few columns")
    return header[4:]


def read_cov_samples(path: str):
    with open(path, "r") as fh:
        header = fh.readline().rstrip("\n").split("\t")
    if len(header) < 2:
        raise ValueError("Covariate file header has too few columns")
    return header[1:]


def read_fam_samples(path: str):
    samples = []
    with open(path, "r") as fh:
        for line in fh:
            if line.strip():
                fields = line.rstrip("\n").split()
                samples.append(fields[1])
    return samples


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate TensorQTL input consistency.")
    parser.add_argument("-g", "--genotype-prefix", required=True)
    parser.add_argument("-p", "--phenotype-bed", required=True)
    parser.add_argument("-c", "--covariates", required=True)
    args = parser.parse_args()

    fam = f"{args.genotype_prefix}.fam"
    for ext in ("bed", "bim", "fam"):
        path = f"{args.genotype_prefix}.{ext}"
        if not os.path.exists(path):
            raise FileNotFoundError(path)

    bed_samples = read_bed_samples(args.phenotype_bed)
    cov_samples = read_cov_samples(args.covariates)
    fam_samples = read_fam_samples(fam)

    print(f"Phenotype samples: {len(bed_samples)}")
    print(f"Covariate samples: {len(cov_samples)}")
    print(f"Genotype samples : {len(fam_samples)}")

    bed_set = set(bed_samples)
    cov_set = set(cov_samples)
    fam_set = set(fam_samples)

    print(f"Phenotype ∩ Covariates: {len(bed_set & cov_set)}")
    print(f"Phenotype ∩ Genotype  : {len(bed_set & fam_set)}")
    print(f"Covariates ∩ Genotype : {len(cov_set & fam_set)}")

    if bed_samples != cov_samples:
        print("WARNING: Phenotype and covariate sample order is not identical.")
    else:
        print("OK: Phenotype and covariate sample order is identical.")

    missing_from_cov = [x for x in bed_samples if x not in cov_set]
    missing_from_bed = [x for x in cov_samples if x not in bed_set]
    missing_from_fam = [x for x in bed_samples if x not in fam_set]

    if missing_from_cov:
        print("Missing from covariates:", ",".join(missing_from_cov[:10]))
    if missing_from_bed:
        print("Missing from phenotype :", ",".join(missing_from_bed[:10]))
    if missing_from_fam:
        print("Missing from genotype  :", ",".join(missing_from_fam[:10]))


if __name__ == "__main__":
    main()
