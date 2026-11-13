#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 06_ASE_ASTS
# Purpose: Build analysis tables.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

"""Run after all requested site extraction and nominal joins have passed their gates."""

from pathlib import Path
import argparse, json, sys, os
from select_inputs import listed

R = Path(__file__).resolve().parent
if (R.parent / ".analysis_deps").exists():
    sys.path.insert(0, str(R.parent / ".analysis_deps"))
import pandas as pd, numpy as np

parser = argparse.ArgumentParser(description="Build allele-effect and QTL analysis tables")
parser.add_argument("--tissues", type=Path)
parser.add_argument("--output", required=True, type=Path)
parser.add_argument("--aliases", type=Path, default=R / "tissue_aliases.json")
args = parser.parse_args()
aliases = json.loads(args.aliases.read_text(encoding="utf8"))
O = args.output
O.mkdir(parents=True, exist_ok=True)
A = Path(os.environ.get("ALLELIC_DIR", str(R / "input/allelic")))
assert (A / "COMPLETE.json").exists()
s = pd.read_csv(
    A / "all_site_gene_statistics.tsv.gz", sep="\t", dtype={"contig": str}, low_memory=False
)
assert not s.duplicated(["sample", "gene", "variant_id"]).any()
orientation = {
    (aliases.get(x["qtl_tissue"].lower(), x["qtl_tissue"].lower()), x["qtl_type"]): x
    for x in json.loads((R / "results/real_slope_validation.json").read_text())
}
manifest = json.loads((R / "results/source_preflight.json").read_text())
if args.tissues:
    tissues = set(listed(args.tissues, "tissue"))
    if tissues - {x["tissue"] for x in manifest}:
        raise ValueError("Selection contains unknown tissues")
    manifest = [x for x in manifest if x["tissue"] in tissues]
audit = []
written = {"site_enrichment_input": False, "signed_effect_input": False}
schemas = {}


def append_table(frame, name):
    # Stream each tissue/version to disk instead of holding all tissues in RAM.
    if name not in schemas:
        schemas[name] = list(frame.columns)
    assert set(frame.columns) == set(schemas[name]), "Cross-tissue schema mismatch: " + name
    frame = frame[schemas[name]]
    frame.to_csv(
        O / (name + ".partial.tsv.gz"),
        sep="\t",
        index=False,
        mode="a" if written[name] else "w",
        header=not written[name],
        compression="gzip",
    )
    written[name] = True


for spec in manifest:
    t, qt = spec["tissue"], spec["qtl_type"]
    folder = R / "results/nominal_matches" / f"{t}.{qt}"
    assert (folder / "COMPLETE.json").exists(), str(folder) + " not complete"
    complete = json.loads((folder / "COMPLETE.json").read_text())
    if not complete["matched_transcript_SNPs"]:
        audit.append(dict(tissue=t, qtl_type=qt, status="no_exact_matches"))
        continue
    q = pd.read_csv(
        folder / "all_matches.tsv.gz", sep="\t", dtype={"contig": str}, low_memory=False
    )
    assert not q.duplicated(["phenotype_id", "variant_id"]).any()
    gq = q.groupby(["gene", "variant_id"], as_index=False).agg(
        QTL_supported=("nominal_QTL_supported", "max"),
        CS_supported=("in_credible_set", "max"),
        max_PIP=("pip", "max"),
        n_QTL_tested_structures=("canonical", "nunique"),
        af=("af", "median"),
    )
    ss = s[s.tissue == t]
    multi = ss.groupby("variant_id").gene.nunique()
    oa = orientation.get((t, qt), {})
    # Direction is the required alignment check; retain strict numerical reproduction status.
    # Original archived betas are used, never substituted with audit estimates.
    orient_pass = oa.get("status") == "PASS" or (
        oa.get("anchors", 0) >= 32
        and oa.get("sign_agreement", 0) == 1
        and oa.get("median_relative_error", 1) < 0.01
        and oa.get("max_relative_error", 1) < 0.05
    )
    # Read each library only once per tissue/QTL type. Keep individual effects auditable.
    indiv = []
    qkeys = q[["gene", "variant_id", "canonical"]].drop_duplicates()
    for sample, h in ss.groupby("sample"):
        counts = pd.read_csv(
            A / "canonical_counts" / f"{sample}.tsv.gz", sep="\t", dtype={"contig": str}
        )
        c = counts.merge(
            qkeys, on=["gene", "variant_id", "canonical"], validate="one_to_many"
        ).merge(
            h[["gene", "variant_id", "gene_ref", "gene_alt", "asts_q", "excluded9", "donor"]],
            on=["gene", "variant_id"],
            validate="many_to_one",
        )
        c = c[c.refCount + c.altCount >= 10].copy()
        usage_ok = (c.gene_ref >= 10) & (c.gene_alt >= 10) & c.asts_q.notna()
        c["delta_usage"] = c.altCount / c.gene_alt - c.refCount / c.gene_ref
        c["log_OR_usage"] = np.log((c.altCount + 0.5) / (c.gene_alt - c.altCount + 0.5)) - np.log(
            (c.refCount + 0.5) / (c.gene_ref - c.refCount + 0.5)
        )
        c.loc[~usage_ok, ["delta_usage", "log_OR_usage"]] = np.nan
        c["tx_ASE_log2"] = np.log2((c.altCount + 0.5) / (c.refCount + 0.5))
        c["tx_coverage"] = c.refCount + c.altCount
        c["sample"] = sample
        indiv.append(
            c[
                [
                    "sample",
                    "donor",
                    "excluded9",
                    "gene",
                    "variant_id",
                    "canonical",
                    "refCount",
                    "altCount",
                    "gene_ref",
                    "gene_alt",
                    "delta_usage",
                    "log_OR_usage",
                    "tx_ASE_log2",
                    "tx_coverage",
                ]
            ]
        )
    individual = pd.concat(indiv, ignore_index=True)
    individual.to_csv(O / f"{t}.{qt}.individual_usage.tsv.gz", sep="\t", index=False)
    for version in ["exclude9", "all450"]:
        h = ss[~ss.excluded9] if version == "exclude9" else ss
        h = h.assign(ase_sig=h.ase_q < 0.05, asts_sig=h.asts_q < 0.05, asts_test=h.asts_q.notna())
        z = h.groupby(["gene", "variant_id"], as_index=False).agg(
            ase_donors=("donor", "nunique"),
            ase_sig_donors=("ase_sig", "sum"),
            asts_donors=("asts_test", "sum"),
            asts_sig_donors=("asts_sig", "sum"),
            median_coverage=("total_coverage", "median"),
            mean_ASE_log2=("ase_log2_ALT_REF", "mean"),
            ASE_between_donor_SD=("ase_log2_ALT_REF", "std"),
            n_observed_structures=("n_canonical", "median"),
            contig=("contig", "first"),
            position=("position", "first"),
            refAllele=("refAllele", "first"),
            altAllele=("altAllele", "first"),
        )
        z = z.merge(gq, on=["gene", "variant_id"], validate="one_to_one")
        z["single_gene_reporter"] = z.variant_id.map(multi) == 1
        z["tissue"] = t
        z["qtl_type"] = qt
        z["version"] = version
        z["block"] = z.contig.astype(str) + ":" + (z.position // 1000000).astype(str)
        append_table(z, "site_enrichment_input")
        d = individual[~individual.excluded9] if version == "exclude9" else individual
        ag = d.groupby(["gene", "variant_id", "canonical"], as_index=False).agg(
            usage_donors=("delta_usage", "count"),
            mean_delta_usage=("delta_usage", "mean"),
            SD_delta_usage=("delta_usage", "std"),
            mean_log_OR_usage=("log_OR_usage", "mean"),
            mean_tx_ASE_log2=("tx_ASE_log2", "mean"),
            tx_ASE_donors=("tx_ASE_log2", "count"),
            mean_tx_coverage=("tx_coverage", "mean"),
        )
        ef = q.merge(
            ag, on=["gene", "variant_id", "canonical"], how="left", validate="many_to_one"
        ).merge(
            z[
                [
                    "gene",
                    "variant_id",
                    "ase_donors",
                    "ase_sig_donors",
                    "asts_donors",
                    "asts_sig_donors",
                    "mean_ASE_log2",
                    "ASE_between_donor_SD",
                    "median_coverage",
                    "single_gene_reporter",
                    "block",
                ]
            ],
            on=["gene", "variant_id"],
            validate="many_to_one",
        )
        ef["version"] = version
        ef["orientation_validated"] = orient_pass
        ef["numerical_slope_reproduction_status"] = oa.get("status", "unknown")
        ef["primary_structure_eligible"] = ef.n_source_aliases == 1
        append_table(ef, "signed_effect_input")
    audit.append(
        dict(
            tissue=t,
            qtl_type=qt,
            status="ready",
            exact_QTL_pairs=len(q),
            confirmed_gene_SNPs=len(gq),
            orientation_validated=orient_pass,
            source_unique_structure_pairs=int((q.n_source_aliases == 1).sum()),
        )
    )
    print("built", t, qt, len(q), flush=True)
for name, ok in written.items():
    assert ok, "No data for " + name
    (O / (name + ".partial.tsv.gz")).replace(O / (name + ".tsv.gz"))
pd.DataFrame(audit).to_csv(O / "analysis_eligibility_audit.tsv", sep="\t", index=False)
