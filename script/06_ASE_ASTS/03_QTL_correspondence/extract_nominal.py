#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 06_ASE_ASTS
# Purpose: Extract nominal.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

"""Read-only streaming join against actual tested QTL rows; one tissue/type per task."""

from pathlib import Path
import argparse, sys, json, gzip, time, os
import pandas as pd, numpy as np
import pyarrow.parquet as pq

R = Path(__file__).resolve().parent
manifest = json.loads((R / "results/source_preflight.json").read_text())
parser = argparse.ArgumentParser(description="Extract nominal QTL rows for one tissue and QTL type")
parser.add_argument("--tissue", required=True)
parser.add_argument("--qtl", required=True)
args = parser.parse_args()
matches = [row for row in manifest if row["tissue"] == args.tissue and row["qtl_type"] == args.qtl]
if len(matches) != 1:
    parser.error("Expected exactly one matching tissue/QTL record")
x = matches[0]
assert x["status"] == "verified"
t = x["tissue"]
qt = x["qtl_type"]
out = R / "results/nominal_matches" / f"{t}.{qt}"
out.mkdir(parents=True, exist_ok=True)
if (out / "COMPLETE.json").exists():
    sys.exit(0)
A = Path(os.environ.get("ALLELIC_DIR", str(R / "input/allelic")))
assert (A / "COMPLETE.json").exists(), "Allelic bundle must be verified and unpacked first"
req = pd.read_csv(A / "qtl_requests.tsv.gz", sep="\t", dtype={"contig": str})
req = req[req.tissue == t]
assert len(req), "No requests for this tissue"
vids = set(req.variant_id)
mapping = pd.read_csv(R / "input/QTL_phenotype_mapping.tsv.gz", sep="\t")
mapping = mapping[mapping.gene_agrees & mapping.gene.isin(req.gene)]
full = pd.read_csv(x["actual_file"], sep="\t")
assert "pval_nominal_threshold" in full, "Need actual phenotype-specific QTL threshold"
invalid_qtl_rows = int((~np.isfinite(full.qval) | ~np.isfinite(full.pval_nominal_threshold)).sum())
full = full[np.isfinite(full.qval) & np.isfinite(full.pval_nominal_threshold)].copy()
full["phenotype_QTL_sig"] = (full.qval <= 0.05) & full.is_eGene.astype(str).str.lower().isin(
    ["true", "1"]
)
meta = full[["phenotype_id", "qval", "pval_nominal_threshold", "phenotype_QTL_sig"]].merge(
    mapping, on="phenotype_id", validate="one_to_one"
)
phs = set(meta.phenotype_id)
fm = pd.read_csv(R / "input/finemap" / f"{x['qtl_tissue']}.{qt}.finemap.tsv.gz", sep="\t")
fm = fm[fm.phenotype_id.isin(phs) & fm.variant_id.isin(vids)]
fm = fm.groupby(["phenotype_id", "variant_id"], as_index=False).agg(
    pip=("pip", "max"), cs_id=("cs_id", lambda s: ";".join(sorted(set(s.dropna().astype(str)))))
)
audit = []
parts = []
start = time.time()
for f in x["nominal_files"]:
    p = Path(f)
    dest = out / (p.name + ".matched.tsv.gz")
    if dest.exists() and Path(str(dest) + ".json").exists():
        a = json.loads(Path(str(dest) + ".json").read_text())
        audit.append(a)
        parts.append(dest)
        continue
    pf = pq.ParquetFile(p)
    cols = [
        k
        for k in ["phenotype_id", "variant_id", "slope", "slope_se", "pval_nominal", "af"]
        if k in pf.schema_arrow.names
    ]
    n = 0
    hits = 0
    header = True
    with gzip.open(str(dest) + ".partial", "wt") as fh:
        for batch in pf.iter_batches(batch_size=131072, columns=cols):
            b = batch.to_pandas()
            n += len(b)
            b = b[b.variant_id.isin(vids) & b.phenotype_id.isin(phs)]
            if not len(b):
                continue
            b = b.merge(meta, on="phenotype_id", validate="many_to_one").merge(
                req, on=["gene", "variant_id"], validate="many_to_one"
            )
            if not len(b):
                continue
            b = b.merge(fm, on=["phenotype_id", "variant_id"], how="left", validate="many_to_one")
            b["in_credible_set"] = b.cs_id.notna() & b.cs_id.fillna("").ne("") & b.phenotype_QTL_sig
            b["pip"] = b.pip.fillna(0).where(b.phenotype_QTL_sig, 0)
            b["nominal_QTL_supported"] = b.phenotype_QTL_sig & (
                b.pval_nominal <= b.pval_nominal_threshold
            )
            b["qtl_type"] = qt
            b["beta_ALT_provisional"] = b.slope * b.beta_A1_to_ALT
            b.to_csv(fh, sep="\t", index=False, header=header)
            header = False
            hits += len(b)
    Path(str(dest) + ".partial").replace(dest)
    a = dict(file=f, rows_scanned=n, matched_rows=hits)
    Path(str(dest) + ".json").write_text(json.dumps(a))
    audit.append(a)
    parts.append(dest)
    print(t, qt, a, flush=True)
frames = [
    pd.read_csv(p, sep="\t", dtype={"contig": str})
    for p, a in zip(parts, audit)
    if a["matched_rows"]
]
if frames:
    d = pd.concat(frames, ignore_index=True)
    assert not d.duplicated(["phenotype_id", "variant_id"]).any(), "Repeated QTL pair"
    d.to_csv(out / "all_matches.tsv.gz", sep="\t", index=False)
    covered = d[["gene", "variant_id"]].drop_duplicates()
else:
    covered = pd.DataFrame(columns=["gene", "variant_id"])
missing = req.merge(covered.assign(found=True), on=["gene", "variant_id"], how="left")
missing = missing[missing.found.isna()]
missing.to_csv(out / "not_confirmed_jointly_evaluable.tsv.gz", sep="\t", index=False)
(out / "COMPLETE.json").write_text(
    json.dumps(
        dict(
            tissue=t,
            qtl_type=qt,
            requested_gene_SNPs=len(req),
            confirmed_tested_gene_SNPs=len(covered),
            matched_transcript_SNPs=sum(a["matched_rows"] for a in audit),
            excluded_QTL_phenotypes_nonfinite_q_or_threshold=invalid_qtl_rows,
            seconds=time.time() - start,
            files=audit,
        ),
        indent=2,
    )
)
print("COMPLETE", t, qt, len(covered), flush=True)
