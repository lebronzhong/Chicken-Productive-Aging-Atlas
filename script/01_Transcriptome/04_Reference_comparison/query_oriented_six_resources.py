#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 01_Transcriptome
# Purpose: Six-resource comparison.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

import argparse, csv, json, collections, datetime
from pathlib import Path
import gtf_utils as u
from pipeline_7b import assembly_aliases
from chain_metrics_lib import Index, covered, known_site

ROOT = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(
    description="Compare query structures against six reference resources"
)
parser.add_argument("--prepared-7b", type=Path, required=True)
parser.add_argument("--run-7b", type=Path, required=True)
parser.add_argument("--previous-comparison", type=Path, required=True)
parser.add_argument("--prepared-7w", type=Path, required=True)
parser.add_argument("--run-7w", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
parser.add_argument("--expected-models", type=int)
parser.add_argument("--expected-chains", type=int)
parser.add_argument("--expected-previous-candidates", type=int)
args = parser.parse_args()
P7B, R7B, PREVIOUS = args.prepared_7b, args.run_7b, args.previous_comparison
P7W, R7W, OUT = args.prepared_7w, args.run_7w, args.output
s7w = json.load(open(R7W / "comparison_summary.json"))
assert s7w["status"] == "COMPLETED"
OUT.mkdir(parents=True, exist_ok=False)
source = json.load(open(P7B / "source_models.json"))
assert source == json.load(open(P7W / "source_models.json"))
q6 = {m["original_id"]: m for m in source.values()}
native = {u.chain_id(u.chain_key(m)): u.chain_key(m) for m in source.values()}
assert q6 and native and None not in native.values()
if args.expected_models is not None:
    assert len(q6) == args.expected_models
if args.expected_chains is not None:
    assert len(native) == args.expected_chains
q7b, _, e = u.parse_gtf(R7B / "query7b.QC_pass.gtf")
assert not e
q7w, _, e = u.parse_gtf(R7W / "query7w.QC_pass.gtf")
assert not e
modelkeys = {
    "6a": {t: u.chain_key(x) for t, x in q6.items()},
    "7b": {t: u.chain_key(x) for t, x in q7b.items()},
    "7w": {t: u.chain_key(x) for t, x in q7w.items()},
}
eligible = {}
for assembly, run, label, models in [
    ("7b", R7B, "Ensembl114", q7b),
    ("7w", R7W, "Ensembl114_7w", q7w),
]:
    d = {}
    for row in csv.DictReader(open(run / label / "source_chain_comparison.tsv"), delimiter="\t"):
        if row["status"].startswith("UNRESOLVED"):
            continue
        ids = row["original_transcript_ids"].split(";")
        keys = {u.chain_key(models[t]) for t in ids}
        assert len(keys) == 1
        d[row["source_chain_id"]] = next(iter(keys))
    eligible[assembly] = d
downloads = ROOT.parent / "annotation_compare_20260929/work/downloads"
config = json.load(open(ROOT / "config_7b.json"))
paths = {
    "Ensembl106": (
        "6a",
        Path("/path/to/chicken-aging/reference/galGal6/Ensemble/Gallus_gallus.GRCg6a.106.gtf"),
    ),
    "Guan2022": ("6a", downloads / "Guan2022_GRCg6a.gtf"),
    "Ensembl114_7b": ("7b", P7B / "Ensembl114.normalized.gtf"),
    "RefSeq7b": ("7b", P7B / "RefSeq7b.normalized.gtf"),
    "FranceAtlas": ("7b", downloads / "France_GRCg7b_GeneEnrichedAtlas.gtf.gz"),
    "Ensembl114_7w": ("7w", P7W / "Ensembl114_7w.normalized.gtf"),
}
lengths = {f[0]: int(f[1]) for l in open(P7B / "genome7b.fa.fai") if (f := l.split("\t"))}
aliases = assembly_aliases(ROOT / "resources/GRCg7b_assembly_report.txt", lengths)
common = set(eligible["7b"]) & set(eligible["7w"])
report = {
    "denominators": {
        "all_query_models": len(q6),
        "all_query_chains": len(native),
        "eligible_7b_chains": len(eligible["7b"]),
        "eligible_7w_chains": len(eligible["7w"]),
        "common_7b_7w_chains": len(common),
    },
    "query_oriented_references": {},
    "provenance": {},
    "reference_oriented_7w_supplement": {},
    "definitions": {
        "exact": "same contig, strand and complete ordered intron chain; transcript endpoints not required to coincide",
        "partial": "query is a strict contiguous intron subchain of a reference; reference-subchain-of-query is not classified as known query",
        "candidate": "neither exact nor query-subchain match, after mapping QC and reference-error checks",
        "unresolved": "not evaluable; never included as unannotated or novel",
        "site_categories": "coordinate-level donor/acceptor novelty, not formal SQANTI gene assignments",
    },
}
cats = {}
allkeys = {"6a": set(), "7b": set(), "7w": set()}
badinterval = []
for label, (assembly, path) in paths.items():
    print("READ", label, flush=True)
    ref, summary, errs = u.parse_gtf(path)
    report["provenance"][label] = {"path": str(path), "sha256": u.digest(path), "summary": summary}
    if label == "FranceAtlas":
        assert all(reason.startswith("overlapping_or_adjacent_exons:") for _, reason in errs)
        bad = {reason.split(":", 1)[1] for _, reason in errs}
        for t, m in ref.items():
            c = m["chrom"]
            opts = {c} if c in lengths else aliases.get(c, set())
            assert len(opts) == 1
            m["chrom"] = next(iter(opts))
            assert all(1 <= a <= b <= lengths[m["chrom"]] for a, b in m["exons"])
            if t in bad:
                badinterval.append((m["chrom"], m["strand"], m["exons"][0][0], m["exons"][-1][1]))
        ref = {t: m for t, m in ref.items() if t not in bad}
    else:
        assert not errs, errs[:3]
    keys = {u.chain_key(m) for m in ref.values() if len(m["exons"]) > 1}
    allkeys[assembly].update(keys)
    idx = Index(keys)
    cache = {}

    def category(key):
        if key not in cache:
            cat, _ = idx.classify(key)
            if label == "FranceAtlas" and cat in (
                "NEW_COMBINATION_KNOWN_SITES",
                "AT_LEAST_ONE_NEW_SPLICE_SITE",
            ):
                ch, st, edges = key
                lo = min(a for a, b in edges)
                hi = max(b for a, b in edges)
                if any(ch == c and st == s and lo <= b and hi >= a for c, s, a, b in badinterval):
                    cat = "UNRESOLVED_REFERENCE_QC"
            cache[key] = cat
        return cache[key]

    chainmap = native if assembly == "6a" else eligible[assembly]
    chaincats = {
        cid: category(chainmap[cid]) if cid in chainmap else "UNRESOLVED_MAPPING" for cid in native
    }
    cats[label] = chaincats
    cc = collections.Counter(chaincats.values())
    mc = collections.Counter(
        category(modelkeys[assembly][t]) if t in modelkeys[assembly] else "UNRESOLVED_MAPPING"
        for t in q6
    )
    commoncounts = collections.Counter(chaincats[t] for t in common)
    report["query_oriented_references"][label] = {
        "assembly": assembly,
        "chain_counts": dict(cc),
        "chain_percent_all_query": {k: 100 * v / len(native) for k, v in cc.items()},
        "model_counts": dict(mc),
        "model_percent_all_query": {k: 100 * v / len(q6) for k, v in mc.items()},
        "common_chain_counts": dict(commoncounts),
        "structural_increment_chains": cc["NEW_COMBINATION_KNOWN_SITES"]
        + cc["AT_LEAST_ONE_NEW_SPLICE_SITE"],
        "structural_increment_models": mc["NEW_COMBINATION_KNOWN_SITES"]
        + mc["AT_LEAST_ONE_NEW_SPLICE_SITE"],
        "coding_scope": "query was filtered by transcript coding evidence; all reference biotypes are retained to identify previously annotated structures",
    }
    assert sum(cc.values()) == len(native) and sum(mc.values()) == len(q6)
    if assembly == "7w":
        coding = {
            u.chain_key(m)
            for m in ref.values()
            if len(m["exons"]) > 1 and m["biotype"] == "protein_coding"
        }
        report["reference_oriented_7w_supplement"] = covered(coding, set(chainmap.values()))
    print("COUNTS", label, json.dumps(report["query_oriented_references"][label]), flush=True)
    u.dump(OUT / "summary.json", report)
indices = {a: Index(k) for a, k in allkeys.items()}
joint = collections.Counter()
rows = []
step = collections.Counter()
previous = {
    r["source_chain_id"]: r["joint_category"]
    for r in csv.DictReader(open(PREVIOUS / "joint_chain_classification.tsv"), delimiter="\t")
}
for cid in sorted(native):
    if cid not in common:
        jointcat = "OUTSIDE_COMMON_MAPPABLE_SCOPE"
    else:
        states = [cats[l][cid] for l in paths]
        if "EXACT" in states:
            jointcat = "EXACT_IN_AT_LEAST_ONE_RESOURCE"
        elif "QUERY_SUBCHAIN" in states:
            jointcat = "QUERY_SUBCHAIN_IN_AT_LEAST_ONE_RESOURCE"
        elif "UNRESOLVED_REFERENCE_QC" in states:
            jointcat = "UNRESOLVED_REFERENCE_QC"
        else:
            keys = {"6a": native[cid], "7b": eligible["7b"][cid], "7w": eligible["7w"][cid]}
            n = len(keys["6a"][2])
            assert all(len(k[2]) == n for k in keys.values())
            known = all(
                any(known_site(indices[a], keys[a], keys[a][2][i], role) for a in keys)
                for i in range(n)
                for role in ("donor", "acceptor")
            )
            jointcat = "NEW_COMBINATION_KNOWN_SITES" if known else "AT_LEAST_ONE_NEW_SPLICE_SITE"
    joint[jointcat] += 1
    if previous.get(cid) in ("NEW_COMBINATION_KNOWN_SITES", "AT_LEAST_ONE_NEW_SPLICE_SITE"):
        step[jointcat] += 1
    rows.append((cid, jointcat, *[cats[l][cid] for l in paths]))
u.tab(
    OUT / "all_query_chain_classification.tsv", ["source_chain_id", "joint_category", *paths], rows
)
report["six_resource_joint_counts"] = dict(joint)
report["previous_candidates_after_7w"] = dict(step)
if args.expected_previous_candidates is not None:
    assert sum(step.values()) == args.expected_previous_candidates
report["joint_new_candidate_chains"] = (
    joint["NEW_COMBINATION_KNOWN_SITES"] + joint["AT_LEAST_ONE_NEW_SPLICE_SITE"]
)
report["six_resource_joint_percent_all_query"] = {
    k: 100 * v / len(native) for k, v in joint.items()
}
report["six_resource_joint_percent_common_scope"] = {
    k: 100 * v / len(common) for k, v in joint.items() if k != "OUTSIDE_COMMON_MAPPABLE_SCOPE"
}
report["mapping_7w_summary"] = {k: v for k, v in s7w.items() if k != "reference_results"}
report["mapping_run_manifests"] = {
    assembly: json.load(open(run / "run_manifest.json"))
    for assembly, run in [("7b", R7B), ("7w", R7W)]
}
u.dump(OUT / "summary.json", report)
u.tab(
    OUT / "query_oriented_chain_summary.tsv",
    [
        "reference",
        "query_chain_denominator",
        "exact",
        "query_subchain",
        "increment_candidates",
        "unresolved",
    ],
    [
        (
            l,
            len(native),
            r["chain_counts"].get("EXACT", 0),
            r["chain_counts"].get("QUERY_SUBCHAIN", 0),
            r["structural_increment_chains"],
            sum(v for k, v in r["chain_counts"].items() if k.startswith("UNRESOLVED")),
        )
        for l, r in report["query_oriented_references"].items()
    ],
)
u.tab(
    OUT / "query_oriented_model_summary.tsv",
    [
        "reference",
        "query_model_denominator",
        "exact",
        "query_subchain",
        "increment_candidates",
        "unresolved",
    ],
    [
        (
            l,
            len(q6),
            r["model_counts"].get("EXACT", 0),
            r["model_counts"].get("QUERY_SUBCHAIN", 0),
            r["structural_increment_models"],
            sum(v for k, v in r["model_counts"].items() if k.startswith("UNRESOLVED")),
        )
        for l, r in report["query_oriented_references"].items()
    ],
)
print("JOINT", json.dumps(dict(joint)), flush=True)
print("PREVIOUS_CANDIDATES", json.dumps(dict(step)), flush=True)
print("QUERY_ORIENTED_COMPLETED", flush=True)
