#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 01_Transcriptome
# Purpose: Summarize expression.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

import csv, json, gzip, math, statistics, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
samples = list(csv.DictReader((ROOT / "sample_manifest.tsv").open(), delimiter="\t"))


def load(path):
    with path.open() as f:
        rows = list(csv.DictReader((x for x in f if not x.startswith("#")), delimiter="\t"))
    last = list(rows[0])[-1]
    return {r["Geneid"]: (int(r["Length"]), float(r[last])) for r in rows}


def ranks(x):
    out = [0.0] * len(x)
    ss = sorted(range(len(x)), key=lambda i: x[i])
    a = 0
    while a < len(ss):
        b = a + 1
        while b < len(ss) and x[ss[b]] == x[ss[a]]:
            b += 1
        for k in ss[a:b]:
            out[k] = (a + b - 1) / 2 + 1
        a = b
    return out


def rho(x, y):
    xx = ranks(x)
    yy = ranks(y)
    ax = statistics.mean(xx)
    ay = statistics.mean(yy)
    return sum((a - ax) * (b - ay) for a, b in zip(xx, yy)) / math.sqrt(
        sum((a - ax) ** 2 for a in xx) * sum((b - ay) ** 2 for b in yy)
    )


out = (
    ROOT
    / "reports"
    / datetime.datetime.now(datetime.timezone.utc).strftime("expression_%Y%m%dT%H%M%SZ")
)
out.mkdir()
summary = []
with gzip.open(out / "paired_gene_abundance.tsv.gz", "wt") as f:
    f.write(
        "sample\ttissue\tgene_id\tgene_union_exon_length\tONT_reads\tRNA_fragments\tONT_CPM\tRNA_length_adjusted_TPM\tRNA_CPM\tin_primary_correlation\n"
    )
    for s in samples:
        name = s["ont_sample"]
        op = ROOT / "gene_expression/ont" / name
        rp = ROOT / "gene_expression/rna" / name
        assert (op / "COMPLETE.json").exists() and (rp / "COMPLETE.json").exists()
        a = load(op / "counts.txt")
        b = load(rp / "counts.txt")
        assert set(a) == set(b)
        assert all(a[g][0] == b[g][0] for g in a)
        ta = sum(x[1] for x in a.values())
        tb = sum(x[1] for x in b.values())
        scale = sum(x[1] / x[0] for x in b.values() if x[0] > 0)
        values = []
        for g in sorted(a):
            length, oc = a[g]
            sc = b[g][1]
            ov = oc / ta * 1e6
            sv = sc / length / scale * 1e6
            rawsv = sc / tb * 1e6
            keep = oc >= 10 and sc >= 10
            f.write(
                f'{name}\t{s["tissue_en"]}\t{g}\t{length}\t{oc}\t{sc}\t{ov}\t{sv}\t{rawsv}\t{int(keep)}\n'
            )
            if keep:
                values.append((ov, sv, rawsv))
        summary.append(
            {
                "sample": name,
                "tissue": s["tissue_en"],
                "replicate": s["replicate_label"],
                "genes_at_least_10_counts_in_both": len(values),
                "spearman_ONT_CPM_vs_RNA_length_adjusted_TPM": rho(
                    [x[0] for x in values], [x[1] for x in values]
                ),
                "spearman_unadjusted_CPM_sensitivity": rho(
                    [x[0] for x in values], [x[2] for x in values]
                ),
                "ONT_gene_assigned_reads": ta,
                "RNA_gene_assigned_fragments": tb,
            }
        )
with (out / "sample_expression_concordance.tsv").open("w") as f:
    w = csv.DictWriter(f, fieldnames=list(summary[0]), delimiter="\t")
    w.writeheader()
    w.writerows(summary)
(out / "summary.json").write_text(
    json.dumps(
        {
            "status": "COMPLETE",
            "samples": summary,
            "method": "Ensembl106 exon-union gene assignment, unstranded. ONT primary MAPQ>=20 reads; RNA uniquely assigned paired fragments. Primary Spearman uses genes with >=10 counts in each platform; ONT CPM vs RNA abundance adjusted by reference gene union-exon length. CPM-vs-CPM is sensitivity analysis. This is approximate gene-level relative abundance concordance, not transcript structure validation or isoform quantification.",
            "median_spearman": statistics.median(
                x["spearman_ONT_CPM_vs_RNA_length_adjusted_TPM"] for x in summary
            ),
        },
        indent=2,
    )
)
print("REPORT_DIR=" + str(out))
print(
    "MEDIAN_SPEARMAN",
    statistics.median(x["spearman_ONT_CPM_vs_RNA_length_adjusted_TPM"] for x in summary),
)
