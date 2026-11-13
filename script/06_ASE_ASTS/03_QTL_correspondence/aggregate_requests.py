#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 06_ASE_ASTS
# Purpose: Aggregate requests.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

from pathlib import Path
import argparse, json, sys, tarfile
from select_inputs import add_selection_args, select_frame
import pandas as pd

R = Path(__file__).resolve().parent
L = Path("/path/to/chicken-aging/allelic")
m = pd.read_csv(L / "scripts/lorals_manifest_matched.tsv", sep="\t")
audit = []
requests = []
stats = []
parser = argparse.ArgumentParser(description="Aggregate selected allele-specific sample results")
add_selection_args(parser)
parser.add_argument("--output", required=True, type=Path)
args = parser.parse_args()
m = select_frame(m, args)
O = args.output
O.mkdir(parents=True, exist_ok=True)
for s in m.ont_sample:
    p = R / "samples" / s
    assert (p / "SUCCESS.json").exists(), s + " incomplete"
    audit.append(json.loads((p / "SUCCESS.json").read_text()))
    requests.append(pd.read_csv(p / "qtl_requests.tsv.gz", sep="\t", dtype={"contig": str}))
    stats.append(pd.read_csv(p / "site_gene_statistics.tsv.gz", sep="\t", dtype={"contig": str}))
req = pd.concat(requests, ignore_index=True).drop_duplicates()
assert not req.duplicated(["tissue", "gene", "variant_id"]).any()
req.to_csv(O / "qtl_requests.tsv.gz", sep="\t", index=False)
pd.concat(stats, ignore_index=True).to_csv(
    O / "all_site_gene_statistics.tsv.gz", sep="\t", index=False
)
pd.DataFrame(audit).to_csv(O / "sample_audit.tsv", sep="\t", index=False)
(O / "COMPLETE.json").write_text(
    json.dumps(
        dict(samples=len(m), request_pairs=len(req), sample_site_rows=sum(len(x) for x in stats)),
        indent=2,
    )
)
with tarfile.open(O / "allelic_site_bundle.tar.gz", "w:gz") as tar:
    for f in [
        "qtl_requests.tsv.gz",
        "all_site_gene_statistics.tsv.gz",
        "sample_audit.tsv",
        "COMPLETE.json",
    ]:
        tar.add(O / f, arcname=f)
    # Canonical counts retained for exact transcript effects; no raw BAM transfers.
    for s in m.ont_sample:
        tar.add(
            R / "samples" / s / "canonical_counts.tsv.gz",
            arcname="canonical_counts/" + s + ".tsv.gz",
        )
print("COMPLETE", len(req), flush=True)
