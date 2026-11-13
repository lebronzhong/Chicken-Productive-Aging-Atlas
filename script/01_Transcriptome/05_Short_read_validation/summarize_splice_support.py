#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 01_Transcriptome
# Purpose: Summarize.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

import csv, gzip, json, pickle, collections, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
db = pickle.load((ROOT / "prepared/database.pkl").open("rb"))
samples = list(csv.DictReader((ROOT / "sample_manifest.tsv").open(), delimiter="\t"))
chains = db["chains"]
junctions = db["junctions"]
ont = {}
rna = {}
amb = {}
qc = []
for s in samples:
    name = s["ont_sample"]
    op = ROOT / "ont" / name
    rp = ROOT / "rnaseq" / name
    oq = json.loads((op / "COMPLETE.json").read_text())
    rq = json.loads((rp / "COMPLETE.json").read_text())
    assert oq["status"] == rq["status"] == "COMPLETE"
    with gzip.open(op / "chain_counts.tsv.gz", "rt") as f:
        ont[name] = {
            r["source_chain_id"]: int(r["ONT_full_chain_reads"])
            for r in csv.DictReader(f, delimiter="\t")
        }
    with gzip.open(rp / "junction_counts.tsv.gz", "rt") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    rna[name] = {int(r["junction_id"]): int(r["unique_read_alignments"]) for r in rows}
    amb[name] = {int(r["junction_id"]): int(r["strand_unresolved_coordinate_reads"]) for r in rows}
    fp = json.loads((rp / "fastp.json").read_text())
    star = {}
    for line in (rp / "STAR_Log.final.out").read_text().splitlines():
        if "|" in line:
            k, v = line.split("|", 1)
            star[k.strip()] = v.strip()
    qc.append(
        {
            "sample": name,
            "tissue": s["tissue_en"],
            "replicate": s["replicate_label"],
            "fastp_before_reads": fp["summary"]["before_filtering"]["total_reads"],
            "fastp_after_reads": fp["summary"]["after_filtering"]["total_reads"],
            "STAR_unique_mapping_percent": star.get("Uniquely mapped reads %"),
            "ONT_primary_mapped_reads": oq.get("primary_mapped_reads", 0),
            "ONT_exact_chain_reads": oq.get("exact_query_chain_reads", 0),
            "ONT_detected_query_chains": len(ont[name]),
            "STAR": star,
        }
    )
out = (
    ROOT
    / "reports"
    / datetime.datetime.now(datetime.timezone.utc).strftime("validation_%Y%m%dT%H%M%SZ")
)
out.mkdir()


def table(name, rows):
    if not rows:
        return
    with (out / name).open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


stats = []
sample_metrics = []
tissue_stats = []
union_detected = collections.defaultdict(set)
union_supported = collections.defaultdict(set)
union_replicated = collections.defaultdict(set)
with gzip.open(out / "chain_tissue_evidence.tsv.gz", "wt") as f:
    fields = [
        "tissue",
        "source_chain_id",
        "category",
        "n_introns",
        "ONT_full_chain_reads",
        "ONT_replicates_detected",
        "minimum_pooled_SR_reads",
        "n_junctions_pooled_SR_ge3",
        "n_replicates_all_junctions_SR_ge3",
        "n_matched_replicates_ONT_and_all_SR_ge3",
        "n_junctions_with_unresolved_strand_coordinate_support_ge3",
    ]
    w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
    w.writeheader()
    for tissue in dict.fromkeys(s["tissue_en"] for s in samples):
        names = [s["ont_sample"] for s in samples if s["tissue_en"] == tissue]
        assert len(names) == 3
        ids = set().union(*(ont[name] for name in names))
        pergroup = collections.defaultdict(list)
        for cid in sorted(ids):
            c = chains[cid]
            js = c["junction_ids"]
            oc = [ont[name].get(cid, 0) for name in names]
            per = [[rna[name].get(j, 0) for j in js] for name in names]
            pooled = [sum(v) for v in zip(*per)]
            mins = [min(v) for v in per]
            ac = [sum(amb[name].get(j, 0) for name in names) for j in js]
            row = dict(
                zip(
                    fields,
                    [
                        tissue,
                        cid,
                        c["category"],
                        len(js),
                        sum(oc),
                        sum(x > 0 for x in oc),
                        min(pooled),
                        sum(x >= 3 for x in pooled),
                        sum(x >= 3 for x in mins),
                        sum(x >= 3 and y > 0 for x, y in zip(mins, oc)),
                        sum(x >= 3 for x in ac),
                    ],
                )
            )
            w.writerow(row)
            pergroup[c["category"]].append((cid, js, oc, pooled, mins, ac))
            pergroup["ALL"].append((cid, js, oc, pooled, mins, ac))
        for group, records in pergroup.items():
            for ontmin in [1, 3]:
                rr = [r for r in records if sum(r[2]) >= ontmin]
                if not rr:
                    continue
                jset = set(j for r in rr for j in r[1])
                for cutoff in [1, 3, 5]:
                    uk = (group, ontmin, cutoff)
                    union_detected[uk].update(r[0] for r in rr)
                    union_supported[uk].update(r[0] for r in rr if min(r[3]) >= cutoff)
                    union_replicated[uk].update(
                        r[0] for r in rr if sum(x >= cutoff for x in r[4]) >= 2
                    )
                    stat = {
                        "tissue": tissue,
                        "category": group,
                        "minimum_tissue_ONT_reads": ontmin,
                        "minimum_SR_reads": cutoff,
                        "detected_chain_denominator": len(rr),
                        "all_junctions_supported_pooled": sum(min(r[3]) >= cutoff for r in rr),
                        "all_junctions_supported_in_2_replicates": sum(
                            sum(x >= cutoff for x in r[4]) >= 2 for r in rr
                        ),
                        "all_junctions_supported_in_2_matched_ONT_replicates": sum(
                            sum(x >= cutoff and y > 0 for x, y in zip(r[4], r[2])) >= 2 for r in rr
                        ),
                        "unique_junction_denominator": len(jset),
                        "junctions_supported_pooled": sum(
                            sum(rna[name].get(j, 0) for name in names) >= cutoff for j in jset
                        ),
                    }
                    stat["all_junction_support_fraction"] = stat[
                        "all_junctions_supported_pooled"
                    ] / len(rr)
                    stat["junction_support_fraction"] = stat["junctions_supported_pooled"] / len(
                        jset
                    )
                    stats.append(stat)
        tissue_stats.append(
            {
                "tissue": tissue,
                "ONT_detected_chains": len(ids),
                "ONT_detected_joint_candidate_chains": sum(
                    chains[cid]["category"]
                    in ("NEW_COMBINATION_KNOWN_SITES", "AT_LEAST_ONE_NEW_SPLICE_SITE")
                    for cid in ids
                ),
            }
        )
for s in samples:
    name = s["ont_sample"]
    groups = collections.defaultdict(list)
    for cid, num in ont[name].items():
        groups["ALL"].append(cid)
        groups[chains[cid]["category"]].append(cid)
    for group, ids in groups.items():
        for cutoff in [1, 3, 5]:
            js = set(j for cid in ids for j in chains[cid]["junction_ids"])
            supported = sum(
                all(rna[name].get(j, 0) >= cutoff for j in chains[cid]["junction_ids"])
                for cid in ids
            )
            sample_metrics.append(
                {
                    "sample": name,
                    "tissue": s["tissue_en"],
                    "category": group,
                    "minimum_SR_reads": cutoff,
                    "detected_chain_denominator": len(ids),
                    "all_junctions_supported": supported,
                    "all_junction_support_fraction": supported / len(ids),
                    "unique_junction_denominator": len(js),
                    "junctions_supported": sum(rna[name].get(j, 0) >= cutoff for j in js),
                }
            )
allids = set().union(*ont.values())
candidates = {
    cid
    for cid in allids
    if chains[cid]["category"] in ("NEW_COMBINATION_KNOWN_SITES", "AT_LEAST_ONE_NEW_SPLICE_SITE")
}
summary = {
    "status": "COMPLETE",
    "samples": 21,
    "tissues": 7,
    "age_weeks": 80,
    "library": "unstranded",
    "independent_of_ONT_GTF_construction": True,
    "scope": "Chains observed as exact full intron-chain matches in primary ONT alignments with MAPQ>=20 in these 21 samples; at least one read for main scope, >=3 pooled reads sensitivity. No full-atlas accuracy claim.",
    "ONT_detected_chain_union": len(allids),
    "ONT_detected_joint_candidate_union": len(candidates),
    "short_read_policy": "NH=1 primary read alignments; >=10 aligned bases immediately adjacent to both sides of each junction. RNA strand from genomic splice motif; noncanonical strand-unresolved coordinate evidence excluded from primary support counts and retained separately. Mate reads may contribute separately.",
    "chain_support_interpretation": "Every junction observed; does not phase a full transcript or establish its ends.",
    "provenance": json.loads((ROOT / "prepared/provenance.json").read_text()),
    "tissue_summary": tissue_stats,
    "sample_qc": qc,
}
union_rows = []
for key, ids in sorted(union_detected.items()):
    group, ontmin, cutoff = key
    union_rows.append(
        {
            "category": group,
            "minimum_tissue_ONT_reads": ontmin,
            "minimum_SR_reads": cutoff,
            "unique_observed_chain_denominator": len(ids),
            "all_junctions_supported_in_at_least_one_matched_tissue": len(union_supported[key]),
            "all_junctions_supported_in_2_replicates_of_at_least_one_matched_tissue": len(
                union_replicated[key]
            ),
            "support_fraction": len(union_supported[key]) / len(ids),
        }
    )
summary["union_support_policy"] = (
    "Each chain counted once. All its junctions must meet the cutoff within at least one tissue where the chain was detected by ONT; evidence is not pooled across different tissues."
)
table("unique_chain_union_support.tsv", union_rows)
table("tissue_support_metrics.tsv", stats)
table("sample_support_metrics.tsv", sample_metrics)
table("tissue_scope.tsv", tissue_stats)
table("sample_qc.tsv", [{k: v for k, v in row.items() if k != "STAR"} for row in qc])
(out / "summary.json").write_text(json.dumps(summary, indent=2))
print("REPORT_DIR=" + str(out))
print(
    json.dumps({k: v for k, v in summary.items() if k not in ("sample_qc", "provenance")}, indent=2)
)
