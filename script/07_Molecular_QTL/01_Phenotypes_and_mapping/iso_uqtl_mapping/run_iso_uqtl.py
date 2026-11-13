#!/usr/bin/env python3
"""Run the recovered iso-uQTL stages for explicit tissue input records."""
import argparse, csv, json, re, subprocess, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PREP = HERE.parent / "tuqtl_psi_preprocess_package/scripts"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--tissue", help="Optional tissue ID; omit to process all manifest rows")
    parser.add_argument(
        "--chromosomes",
        type=Path,
        help="Optional TSV with a chromosome column for nominal tests only",
    )
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    manifest = args.manifest.resolve()
    with manifest.open(encoding="utf8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {
        "tissue",
        "transcript_tpm",
        "genotype",
        "covariates",
        "salmon_dir",
        "tx2gene",
        "gtf",
        "suppa",
        "output_dir",
        "phenotype_filter",
    }
    if not rows or not all(required <= set(row) for row in rows):
        parser.error("Input manifest lacks required fields or records")
    if len({row["tissue"] for row in rows}) != len(rows):
        parser.error("Duplicate tissue rows")
    if args.tissue:
        rows = [row for row in rows if row["tissue"] == args.tissue]
        if not rows:
            parser.error("Requested tissue is absent from the manifest")
    chromosomes = None
    if args.chromosomes:
        with args.chromosomes.open(encoding="utf8") as handle:
            chromosomes = [row["chromosome"] for row in csv.DictReader(handle, delimiter="\t")]
        if not chromosomes or len(set(chromosomes)) != len(chromosomes):
            parser.error("Chromosome selection must be nonempty and unique")
    for row in rows:
        tissue = row["tissue"]
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", tissue) or tissue in {".", ".."}:
            parser.error("Tissue ID must be a safe path component")
        if row["phenotype_filter"] not in {"variance", "none"}:
            parser.error("phenotype_filter must explicitly specify variance or none")

        def path(key):
            p = Path(row[key])
            return p.resolve() if p.is_absolute() else (manifest.parent / p).resolve()

        out = path("output_dir")
        psi = out / "psi_phenotype"
        results = out / "tuqtl_mapping/results"
        split = out / "tuqtl_mapping/phenotype_by_chr"
        logs = out / "tuqtl_mapping/logs"

        def execute(command, name):
            command = [str(x) for x in command]
            if args.dry_run:
                print(json.dumps({"tissue": tissue, "stage": name, "command": command}))
                return
            logs.mkdir(parents=True, exist_ok=True)
            with (logs / (name + ".log")).open("w", encoding="utf8") as handle:
                subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, check=True)

        tpm, genotype, cov = path("transcript_tpm"), path("genotype"), path("covariates")
        execute(
            [
                "bash",
                PREP / "build_tuQTL_psi_phenotype.sh",
                "--tissues",
                tissue,
                "--transcript-tpm-path",
                tpm.parent,
                "--transcript-tpm-file",
                tpm,
                "--output-path",
                psi,
                "--gene-expr-path",
                out / "gene_expression",
                "--genotype-path",
                genotype.parent,
                "--genotype-prefix",
                genotype,
                "--salmon-out",
                path("salmon_dir"),
                "--tx2gene",
                path("tx2gene"),
                "--gtf",
                path("gtf"),
                "--suppa",
                path("suppa"),
                "--script-dir",
                PREP,
            ],
            "01_preprocess",
        )
        stem = f"salmon_pheno.ChickenGTEx.{tissue}.TPM_isoform.filtered.bimchr"
        bed = psi / (stem + ".tensorqtl.bed.gz")
        if row["phenotype_filter"] == "variance":
            filtered = psi / (stem + ".varfiltered.tensorqtl.bed.gz")
            execute(
                [
                    sys.executable,
                    HERE / "filter_tensorqtl_phenotype_bed.py",
                    "--input",
                    bed,
                    "--output",
                    filtered,
                    "--min-var",
                    "1e-4",
                    "--min-iqr",
                    "0.01",
                    "--min-nonzero-frac",
                    "0.2",
                ],
                "02_filter",
            )
            bed = filtered
        execute(
            [
                sys.executable,
                HERE / "check_inputs.py",
                "--genotype",
                genotype,
                "--phenotype",
                bed,
                "--covariates",
                cov,
            ],
            "03_check",
        )
        execute(
            [
                sys.executable,
                HERE / "split_bed_by_chr.py",
                "--input-bed-gz",
                bed,
                "--outdir",
                split,
                "--prefix",
                f"salmon_pheno.ChickenGTEx.{tissue}.PSI.bimchr",
            ],
            "04_split",
        )
        if args.dry_run:
            nominal_inputs = [
                (chrom, split / f"{{prefix}}.chr{chrom}.bed.gz")
                for chrom in (chromosomes or ["{all_input_chromosomes}"])
            ]
        else:
            nominal_inputs = []
            for p in sorted(split.glob("*.bed.gz")):
                match = re.search(r"\.chr([^.]+)\.bed\.gz$", p.name)
                if match and (chromosomes is None or match[1] in chromosomes):
                    nominal_inputs.append((match[1], p))
            if not nominal_inputs:
                raise ValueError("No nominal-test chromosomes found")
            if chromosomes is not None and set(chromosomes) != {x[0] for x in nominal_inputs}:
                raise ValueError("Requested chromosome is absent from phenotype inputs")
        for chrom, chrbed in nominal_inputs:
            dest = results / f"chr{chrom}"
            if not args.dry_run:
                dest.mkdir(parents=True, exist_ok=True)
            execute(
                [
                    sys.executable,
                    "-m",
                    "tensorqtl",
                    genotype,
                    chrbed,
                    dest / f"ChickenAging.{tissue}_PSI_chr{chrom}",
                    "--mode",
                    "cis_nominal",
                    "--covariates",
                    cov,
                    "--maf_threshold",
                    "0.05",
                ],
                f"05_nominal_{chrom}",
            )
        dest = results / "full"
        if not args.dry_run:
            dest.mkdir(parents=True, exist_ok=True)
        prefix = dest / f"ChickenAging.{tissue}_PSI"
        execute(
            [
                sys.executable,
                "-m",
                "tensorqtl",
                genotype,
                bed,
                prefix,
                "--mode",
                "cis",
                "--covariates",
                cov,
                "--maf_threshold",
                "0.05",
            ],
            "06_permutation",
        )
        execute(
            [
                sys.executable,
                HERE / "calc_tensorqtl_fdr_official.py",
                str(prefix) + ".cis_qtl.txt.gz",
                str(prefix) + ".cis_qtl.txt",
                "0.05",
            ],
            "07_fdr",
        )
        if not args.dry_run:
            import pandas as pd

            table = pd.read_csv(str(prefix) + ".cis_qtl.txt", sep="\t")
            keep = table.is_eGene.astype(str).str.lower().eq("true")
            table[keep].to_csv(str(prefix) + ".cis_qtl.fdr.txt", sep="\t", index=False)


if __name__ == "__main__":
    main()
