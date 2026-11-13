#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 01_Transcriptome
# Purpose: Gene counts.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

import csv, json, subprocess, concurrent.futures, argparse
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REF = Path("/path/to/chicken-aging/reference/galGal6/Ensemble/Gallus_gallus.GRCg6a.106.gtf")
FC = ROOT / "software/featureCounts.2.0.8"
p = argparse.ArgumentParser()
p.add_argument("platform", choices=["ont", "rna"])
a = p.parse_args()
samples = list(csv.DictReader((ROOT / "sample_manifest.tsv").open(), delimiter="\t"))


def one(s):
    out = ROOT / "gene_expression" / a.platform / s["ont_sample"]
    out.mkdir(parents=True, exist_ok=False)
    cmd = [
        str(FC),
        "-a",
        str(REF),
        "-o",
        str(out / "counts.txt"),
        "-s",
        "0",
        "-Q",
        "20",
        "-t",
        "exon",
        "-g",
        "gene_id",
    ]
    if a.platform == "ont":
        cmd += ["-L"]
        with (out / "samtools.log").open("w") as slog, (out / "featureCounts.log").open("w") as log:
            process = subprocess.Popen(
                ["samtools", "view", "-h", "-F", "2820", "-q", "20", s["ont_bam"]],
                stdout=subprocess.PIPE,
                stderr=slog,
            )
            try:
                subprocess.run(
                    cmd, stdin=process.stdout, stdout=log, stderr=subprocess.STDOUT, check=True
                )
            finally:
                process.stdout.close()
            assert process.wait() == 0
    else:
        rp = ROOT / "rnaseq" / s["ont_sample"]
        assert (rp / "COMPLETE.json").exists()
        cmd += ["-p", "--countReadPairs", "-T", "4", str(rp / "STAR_Aligned.sortedByCoord.out.bam")]
        with (out / "featureCounts.log").open("w") as log:
            subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)
    (out / "COMPLETE.json").write_text(
        json.dumps(
            {
                "sample": s["ont_sample"],
                "platform": a.platform,
                "command": cmd,
                "status": "COMPLETE",
                "scope": "Ensembl106 reference genes; unstranded exon-union assignment, ambiguous assignments excluded; ONT primary MAPQ20 reads, RNA paired fragments.",
            },
            indent=2,
        )
    )
    print("GENE_COUNTS_COMPLETE", a.platform, s["ont_sample"], flush=True)


with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
    for f in concurrent.futures.as_completed([ex.submit(one, s) for s in samples]):
        f.result()
print("GENE_COUNTS_ALL_COMPLETE", a.platform, flush=True)
