#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 06_ASE_ASTS
# Purpose: Extract allelic sites.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

"""Extract SNP-level allele-specific statistics from matched samples."""

from pathlib import Path
import argparse, sys, json, time, hashlib
from select_inputs import add_selection_args, select_frame
import numpy as np, pandas as pd, pysam
from scipy.stats import binomtest, chi2_contingency, fisher_exact

R = Path(__file__).resolve().parent
L = Path("/path/to/chicken-aging/allelic")
REF = "/path/to/chicken-aging/reference/galGal6/Gallus_gallus.GRCg6a.dna.toplevel.fa"
K = ["contig", "position", "refAllele", "altAllele"]


def bh(p):
    a = np.asarray(p, float)
    q = np.full(len(a), np.nan)
    ix = np.flatnonzero(np.isfinite(a))
    ix = ix[np.argsort(a[ix])]
    q[ix] = np.minimum.accumulate((a[ix] * len(ix) / np.arange(1, len(ix) + 1))[::-1])[::-1].clip(
        0, 1
    )
    return q


def process(row, cm, old, bim, fa, excluded):
    s = row.ont_sample
    out = R / "samples" / s
    out.mkdir(parents=True, exist_ok=True)
    if (out / "SUCCESS.json").exists():
        return
    start = time.time()
    raw = pd.read_csv(
        L / "asts/samples" / row.tissue / row.dataset / s / (s + "_asts_quant.tsv"),
        sep="\t",
        dtype={"contig": str},
    )
    qc = pd.read_csv(
        L / "ase/samples" / row.tissue / row.dataset / s / (s + ".ase.qc.tsv"),
        sep="\t",
        dtype={"contig": str},
    )
    assert not raw.duplicated(K + ["transcript"]).any() and not qc.duplicated(K).any()
    qc = qc[
        qc.GENOTYPE.isin(["0|1", "1|0", "0/1", "1/0"])
        & qc.refAllele.isin(list("ACGT"))
        & qc.altAllele.isin(list("ACGT"))
    ]
    for flag in [
        "GENOTYPE_WARNING",
        "BLACKLIST",
        "MULTI_MAPPING",
        "OTHER_ALLELE_WARNING",
        "HIGH_INDEL_WARNING",
    ]:
        qc = qc[qc[flag] == 0]
    d = (
        raw.merge(qc[K], on=K, validate="many_to_one")
        .merge(cm, on="transcript", how="left", validate="many_to_one")
        .merge(old, on="transcript", how="left", validate="many_to_one")
    )
    bad = d.canonical.isna() | d.gene.isna() | (d.gene != d.canonical_gene)
    bad_rows = int(bad.sum())
    d = d[~bad]
    can = d.groupby(K + ["gene", "canonical"], sort=False, as_index=False)[
        ["refCount", "altCount"]
    ].sum()
    # Reference check and one-to-one BIM identity check; never guess strand complements.
    uniq = can[K].drop_duplicates()
    vrows = []
    refbad = 0
    ambig = 0
    for x in uniq.itertuples():
        chrom = str(x.contig)
        pos = int(x.position)
        if chrom not in fa.references or fa.fetch(chrom, pos - 1, pos).upper() != x.refAllele:
            refbad += 1
            continue
        hits = bim.get((chrom, pos, frozenset([x.refAllele, x.altAllele])), [])
        if len(hits) != 1:
            if len(hits) > 1:
                ambig += 1
            continue
        vid, a1, a2 = hits[0]
        vrows.append(
            dict(
                zip(K, [chrom, pos, x.refAllele, x.altAllele]),
                variant_id=vid,
                A1=a1,
                A2=a2,
                beta_A1_to_ALT=1 if a1 == x.altAllele else -1,
            )
        )
    vm = pd.DataFrame(vrows)
    rows = []
    for key, g in can.groupby(K + ["gene"], sort=False):
        chrom, pos, ref, alt, gene = key
        r = int(g.refCount.sum())
        a = int(g.altCount.sum())
        n = r + a
        if n < 21:
            continue
        h = g[g.refCount + g.altCount >= 10]
        tab = h[["refCount", "altCount"]].to_numpy(dtype=np.int64).T
        asts_p = np.nan
        method = "not_testable"
        emin = np.nan
        efrac = np.nan
        if len(h) >= 2 and tab.sum() >= 36 and min(tab.sum(1)) >= 10:
            ex = np.outer(tab.sum(1), tab.sum(0)) / tab.sum()
            emin = float(ex.min())
            efrac = float((ex < 5).mean())
            if len(h) == 2:
                asts_p = float(fisher_exact(tab).pvalue)
                method = "Fisher_2x2"
            elif emin >= 1 and efrac <= 0.2:
                asts_p = float(chi2_contingency(tab, correction=False).pvalue)
                method = "Pearson_valid_expected"
            else:
                method = "excluded_sparse_RxC"
        rows.append(
            dict(
                contig=chrom,
                position=pos,
                refAllele=ref,
                altAllele=alt,
                gene=gene,
                gene_ref=r,
                gene_alt=a,
                total_coverage=n,
                ase_p=float(binomtest(a, n, 0.5).pvalue),
                ase_log2_ALT_REF=float(np.log2((a + 0.5) / (r + 0.5))),
                asts_p=asts_p,
                asts_method=method,
                asts_expected_min=emin,
                asts_expected_fraction_lt5=efrac,
                n_canonical=len(g),
                n_retained_ASTS=len(h),
            )
        )
    z = pd.DataFrame(rows)
    z["ase_q"] = bh(z.ase_p)
    z["asts_q"] = bh(z.asts_p)
    # FDR is calculated before restricting to QTL-panel SNPs, within each library.
    all_n = len(z)
    z = z.merge(vm, on=K, validate="many_to_one")
    z["sample"] = s
    z["tissue"] = row.tissue
    z["donor"] = row.vcf_sample
    z["age"] = str(row.age)
    z["excluded9"] = s in excluded
    z.to_csv(out / "site_gene_statistics.tsv.gz", sep="\t", index=False)
    c = can.merge(z[K + ["gene", "variant_id"]], on=K + ["gene"], validate="many_to_one")
    c.to_csv(out / "canonical_counts.tsv.gz", sep="\t", index=False)
    z[
        [
            "tissue",
            "gene",
            "variant_id",
            "contig",
            "position",
            "refAllele",
            "altAllele",
            "A1",
            "A2",
            "beta_A1_to_ALT",
        ]
    ].drop_duplicates().to_csv(out / "qtl_requests.tsv.gz", sep="\t", index=False)
    audit = dict(
        sample=s,
        tissue=row.tissue,
        all_allelic_gene_markers=all_n,
        QTL_panel_gene_markers=len(z),
        ASE_positive=int((z.ase_q < 0.05).sum()),
        ASTS_testable=int(z.asts_p.notna().sum()),
        ASTS_positive=int((z.asts_q < 0.05).sum()),
        sparse_RxC_excluded=int((z.asts_method == "excluded_sparse_RxC").sum()),
        annotation_conflict_rows=bad_rows,
        reference_mismatch_sites=refbad,
        ambiguous_BIM_sites=ambig,
        collapsed_transcript_rows=len(d) - len(can),
        seconds=time.time() - start,
    )
    (out / "SUCCESS.json").write_text(json.dumps(audit, indent=2))
    print(json.dumps(audit), flush=True)


def main():
    parser = argparse.ArgumentParser(
        description="Extract allele-specific sites for selected samples"
    )
    add_selection_args(parser)
    parser.add_argument("--shards", type=int, default=1)
    parser.add_argument("--shard-index", type=int, default=0)
    args = parser.parse_args()
    if args.shards < 1 or not 0 <= args.shard_index < args.shards:
        parser.error("Require shards >= 1 and 0 <= shard-index < shards")
    manifest = pd.read_csv(L / "scripts/lorals_manifest_matched.tsv", sep="\t")
    if manifest.ont_sample.duplicated().any():
        raise ValueError("Duplicate sample identifiers")
    manifest = (
        select_frame(manifest, args).sort_values(["tissue", "ont_sample"]).reset_index(drop=True)
    )
    cm = pd.read_csv(R / "input/canonical_transcript_mapping.tsv.gz", sep="\t", dtype=str)
    cm["source_ids"] = cm.source_ids.str.split(",")
    cm = cm.explode("source_ids")[["source_ids", "transcript_id", "gene_id"]].rename(
        columns={
            "source_ids": "transcript",
            "transcript_id": "canonical",
            "gene_id": "canonical_gene",
        }
    )
    assert not cm.transcript.duplicated().any()
    old = pd.read_csv(L / "software/lorals_ref/transcript_gene.tsv", sep="\t")
    excluded = set(pd.read_csv(R / "input/01_excluded_samples.tsv", sep="\t").ont_sample)
    b = pd.read_csv(
        R / "input/all_QTL_variants.bim.gz",
        sep=r"\s+",
        header=None,
        names=["chrom", "variant", "cm", "pos", "A1", "A2"],
        dtype={"chrom": str},
    )
    bim = {}
    for x in b.itertuples():
        if x.A1 not in list("ACGT") or x.A2 not in list("ACGT") or x.A1 == x.A2:
            continue
        bim.setdefault((str(x.chrom), int(x.pos), frozenset([x.A1, x.A2])), []).append(
            (x.variant, x.A1, x.A2)
        )
    fa = pysam.FastaFile(REF)
    ids = range(args.shard_index, len(manifest), args.shards)
    for i in ids:
        process(manifest.iloc[i], cm, old, bim, fa, excluded)


if __name__ == "__main__":
    main()
