#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 07_Molecular_QTL
# Purpose: Summarize regions.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

"""Three counting levels, with separate completeness for conditional results and SuSiE CS."""

import csv, gzip, json, math, os, time, uuid
from collections import Counter
from pathlib import Path

WORK = Path(__file__).resolve().parent
OUT = Path("/path/to/chicken-aging/qtl/conditional/summary")
TISSUES = sorted(
    json.loads((WORK / "study_tissues.json").read_text(encoding="utf8"))["all_tissues"]
)
KINDS = ["teQTL", "tuQTL"]


def read(p):
    with (gzip.open if str(p).endswith(".gz") else open)(p, "rt", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write(p, rows):
    tmp = p.with_name(p.name + "." + uuid.uuid4().hex + ".tmp")
    with tmp.open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    os.replace(tmp, p)


def atomic_text(p, text):
    tmp = p.with_name(p.name + "." + uuid.uuid4().hex + ".tmp")
    tmp.write_text(text, encoding="utf-8-sig")
    os.replace(tmp, p)


def merge_count(intervals):
    chrom = None
    right = -1
    n = 0
    for ch, lo, hi in sorted(intervals):
        assert lo <= hi
        if ch != chrom or lo > right:
            n += 1
            chrom = ch
            right = hi
        else:
            right = max(right, hi)
    return n


def build():
    inventory = []
    for host in ["co", "tj"]:
        text = (WORK / (host + "_cs_aggregate_counts.log")).read_text(encoding="utf-8")
        assert "INVENTORY_DONE" in text and "Traceback" not in text
        inventory += json.loads(
            next(
                line.split("=", 1)[1]
                for line in text.splitlines()
                if line.startswith("AGGREGATE_COUNTS_JSON=")
            )
        )
    sources = {(r["tissue"], r["qtl_type"]): r for r in inventory}
    assert set(sources) == {(t, k) for t in TISSUES for k in KINDS} and len(inventory) == 36
    summary = []
    audits = []
    bykind = {k: [] for k in KINDS}
    allregions = {k: [] for k in KINDS}
    available_regions = {}
    available_allregions = {}
    for tissue in TISSUES:
        for kind in KINDS:
            source = sources[tissue, kind]
            folder = OUT / "details" / tissue / kind
            marker = folder / "COMPLETE.json"
            done = marker.exists()
            row = dict(
                tissue=tissue,
                qtl_type=kind,
                conditional_status="complete" if done else "incomplete",
                conditional_signals="NA",
                rank_ge2="NA",
                retained_transcripts="NA",
                extra_retained_signals="NA",
                credible_sets=source["credible_sets"],
                CS_status="complete",
                matched_credible_sets="NA",
                CS_match_percent="NA",
                CS_supported_conditional_signals="NA",
                CS_supported_nonredundant_regions="NA",
                all_conditional_nonredundant_regions="NA",
                host=source["host"],
            )
            if done:
                info = json.loads(marker.read_text())
                sig = read(folder / "cis_independent.tsv.gz")
                matches = read(folder / "signal_CS_matches.tsv.gz")
                edges = read(folder / "signal_CS_all_edges.tsv.gz")
                regions = [r for p in folder.glob("chr*.regions.tsv") for r in read(p)]
                checks = [r for p in folder.glob("chr*.input_check.tsv") for r in read(p)]
                tests = [r for p in folder.glob("chr*.conditional_tests.tsv") for r in read(p)]
                key = lambda r: (r["phenotype_id"], int(float(r["rank"])))
                assert len(sig) == len(matches) == len(regions) == info["signals"]
                assert len({key(r) for r in sig}) == len(sig)
                assert (
                    {key(r) for r in sig} == {key(r) for r in matches} == {key(r) for r in regions}
                )
                assert len(checks) == info["targets"] and all(r["pass"] == "True" for r in checks)
                assert tests and all(
                    math.isfinite(float(r["candidate_pval_beta"]))
                    and 0 <= float(r["candidate_pval_beta"]) <= 1
                    for r in tests
                )
                cutoff = max(float(r["threshold"]) for r in tests)
                for r in sig:
                    p = float(r.get("pval_conditional") or r["pval_beta"])
                    assert math.isfinite(p) and 0 <= p <= cutoff * (1 + 1e-6)
                    assert all(
                        math.isfinite(float(r[c])) for c in ["pval_nominal", "slope", "slope_se"]
                    )
                    assert float(r["rank"]).is_integer() and float(r["rank"]) >= 1
                matched = {
                    (r["phenotype_id"], r["cs_id"])
                    for r in edges
                    if r["direct"] == "True" or float(r["max_r2"]) >= 0.8
                }
                supported_keys = {
                    key(r) for r in edges if r["direct"] == "True" or float(r["max_r2"]) >= 0.8
                }
                assert len(matched) == info["matched_cs"] <= info["cs"] == source["credible_sets"]
                if source["matched_credible_sets"] is not None:
                    assert len(matched) == source["matched_credible_sets"]
                assert supported_keys == {key(r) for r in regions if r["CS_supported"] == "True"}
                assert supported_keys == {
                    key(r) for r in matches if r["category"] in ["direct", "LD_only", "multiple_CS"]
                }
                assert Counter(r["category"] for r in matches) == info["categories"]
                rank2 = sum(key(r)[1] >= 2 for r in sig)
                assert rank2 == info["rank_ge2"]
                perpid = Counter(r["phenotype_id"] for r in sig)
                triples = [(r["chr"], int(r["start"]), int(r["end"])) for r in regions]
                supported = [
                    (r["chr"], int(r["start"]), int(r["end"]))
                    for r in regions
                    if r["CS_supported"] == "True"
                ]
                bykind[kind].extend(supported)
                allregions[kind].extend(triples)
                available_regions[tissue, kind] = merge_count(supported)
                available_allregions[tissue, kind] = merge_count(triples)
                row.update(
                    conditional_signals=len(sig),
                    rank_ge2=rank2,
                    retained_transcripts=len(perpid),
                    extra_retained_signals=len(sig) - len(perpid),
                    matched_credible_sets=len(matched),
                    CS_match_percent=(
                        round(100 * len(matched) / info["cs"], 2) if info["cs"] else "NA"
                    ),
                    CS_supported_conditional_signals=len(supported),
                    CS_supported_nonredundant_regions=merge_count(supported),
                    all_conditional_nonredundant_regions=merge_count(triples),
                )
                repairs = folder / "numerical_repairs.jsonl"
                nr = (
                    sum(bool(line.strip()) for line in repairs.read_text().splitlines())
                    if repairs.exists()
                    else 0
                )
                audits.append(
                    dict(
                        tissue=tissue,
                        qtl_type=kind,
                        input_checks=len(checks),
                        conditional_tests=len(tests),
                        retained_signal_rows=len(sig),
                        CS_inventory_count=source["credible_sets"],
                        matched_CS=len(matched),
                        empirical_fallback_calls=nr,
                        retained_empirical_signals=sum(
                            r.get("conditional_p_method") == "empirical_fallback" for r in sig
                        ),
                        passed=True,
                    )
                )
            summary.append(row)
    totals = []
    for kind in KINDS:
        allrows = [r for r in summary if r["qtl_type"] == kind]
        done = [r for r in allrows if r["conditional_status"] == "complete"]
        n = len(done)
        complete = n == 18
        for scope, rows in [("ALL_18_TISSUES", allrows), ("COMPLETED_CONDITIONAL_SUBSET", done)]:
            usable = complete or scope == "COMPLETED_CONDITIONAL_SUBSET"
            totals.append(
                dict(
                    qtl_type=kind,
                    scope=scope,
                    tissues_in_scope=len(rows),
                    conditional_completed_tissues=n,
                    conditional_status="FINAL" if complete else "PARTIAL",
                    conditional_signals=(
                        sum(r["conditional_signals"] for r in done) if usable else "NA"
                    ),
                    rank_ge2=sum(r["rank_ge2"] for r in done) if usable else "NA",
                    extra_retained_signals=(
                        sum(r["extra_retained_signals"] for r in done) if usable else "NA"
                    ),
                    credible_sets=sum(r["credible_sets"] for r in rows),
                    CS_status=(
                        "FINAL_18_TISSUES"
                        if scope == "ALL_18_TISSUES"
                        else "COMPLETED_CONDITIONAL_SUBSET"
                    ),
                    matched_credible_sets=(
                        sum(r["matched_credible_sets"] for r in done) if usable else "NA"
                    ),
                    CS_supported_conditional_signals=(
                        sum(r["CS_supported_conditional_signals"] for r in done) if usable else "NA"
                    ),
                    CS_supported_nonredundant_regions_union=(
                        merge_count(bykind[kind]) if usable else "NA"
                    ),
                    all_conditional_nonredundant_regions_union=(
                        merge_count(allregions[kind]) if usable else "NA"
                    ),
                    missing_conditional_tissues=",".join(
                        r["tissue"] for r in allrows if r["conditional_status"] != "complete"
                    ),
                    included_tissues=",".join(r["tissue"] for r in rows),
                )
            )
    # Compare independently recomputed interval unions with the existing region membership export.
    existing = OUT / "merged_region_coordinates.tsv.gz"
    if existing.exists():
        counts = Counter((r["qtl_type"], r["scope"], r["support"]) for r in read(existing))
        old = read(OUT / "all_18tissues_summary.tsv")
        old_keys = {(r["tissue"], r["qtl_type"]) for r in old if r["status"] == "complete"}
        if old_keys == set(available_regions):
            for (t, k), n in available_regions.items():
                assert counts[k, t, "CS_supported"] == n
            for k in KINDS:
                assert counts[k, "ALL_TISSUES", "CS_supported"] == merge_count(bykind[k])
    write(OUT / "three_levels_by_tissue.tsv", summary)
    write(OUT / "three_levels_totals.tsv", totals)
    write(OUT / "three_levels_validation.tsv", audits)
    write(OUT / "three_levels_CS_source_counts.tsv", inventory)
    now = time.strftime("%Y-%m-%d %H:%M:%S %z")
    pending = [
        r["tissue"] + " " + r["qtl_type"] for r in summary if r["conditional_status"] != "complete"
    ]
    note = f"""Conditional QTL summary generated at {now}.
Completed conditional tasks: {len(audits)}/36.
Pending tasks: {'; '.join(pending) or 'none'}.
Missing results are NA, not zero; incomplete tissues do not contribute partial chromosome totals.
Conditional signals, credible sets and LD-defined regions are distinct statistical units.
Across-tissue region counts are recomputed by interval merging, not summed across tissues.
Credible-set-supported regions require matching conditional signals; they are not all SuSiE credible sets.
The completed-tissue subset is not the final 18-tissue total.
See the validation and source-count tables for the evidence supporting each completed task.
"""
    atomic_text(OUT / "summary_scope.txt", note)
    print(
        json.dumps(
            dict(updated=now, completed=len(audits), CS_groups=len(inventory), totals=totals),
            ensure_ascii=False,
        ),
        flush=True,
    )
    return len(audits) == 36


if __name__ == "__main__":
    build()
