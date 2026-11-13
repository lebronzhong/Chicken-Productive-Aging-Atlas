#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 01_Transcriptome
# Purpose: Parse transcript models and exon coordinates from GTF annotation.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

"""Stage-one, standard-library-only annotation audit and native-chain comparison.
Coordinates: GTF 1-based closed exon intervals; intron keys use flanking
exon boundaries (left_exon_end, right_exon_start), always with contig and strand.
Never interprets missing downloads or failed input validation as novelty.
"""

import argparse
import collections
import csv
import datetime
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import re
import shutil
import sys
import time
import urllib.request
import zipfile

ROOT = Path(__file__).resolve().parent
ATTR = re.compile(r'(\S+)\s+"([^"]*)"\s*(?:;|$)')


def dump(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def digest(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def optext(path):
    return (
        gzip.open(path, "rt", encoding="utf-8-sig")
        if str(path).endswith(".gz")
        else open(path, encoding="utf-8-sig")
    )


def tab(path, header, rows):
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def parse_gtf(path):
    models, errors, comments = {}, [], []
    features = collections.Counter()
    declared = set()

    def issue(line, reason):
        errors.append((line, reason))

    with optext(path) as handle:
        for n, line in enumerate(handle, 1):
            if line.startswith("#"):
                if len(comments) < 15:
                    comments.append(line.rstrip())
                continue
            if not line.strip():
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                issue(n, "not_9_columns")
                continue
            chrom, source, feature, left, right, score, strand, frame, attrs = fields
            features[feature] += 1
            if feature not in ("exon", "CDS", "transcript", "mRNA"):
                continue
            a = dict(ATTR.findall(attrs))
            tid, gid = a.get("transcript_id"), a.get("gene_id")
            if not tid or not gid:
                issue(n, "missing_GTF_transcript_id_or_gene_id")
                continue
            try:
                left, right = int(left), int(right)
                if left < 1 or right < left or strand not in ("+", "-", "."):
                    raise ValueError()
            except ValueError:
                issue(n, "invalid_coordinate_or_strand")
                continue
            if tid not in models:
                models[tid] = {
                    "gene": gid,
                    "chrom": chrom,
                    "strand": strand,
                    "exons": [],
                    "cds": [],
                    "biotype": a.get("transcript_biotype", a.get("transcript_type", "NA")),
                }
            m = models[tid]
            if (m["gene"], m["chrom"], m["strand"]) != (gid, chrom, strand):
                issue(n, "conflicting_transcript_id:" + tid)
                continue
            if feature == "exon":
                m["exons"].append((left, right))
            elif feature == "CDS":
                m["cds"].append((left, right))
            else:
                declared.add(tid)
    duplicate_exons = 0
    for tid, m in models.items():
        exons = sorted(set(m["exons"]))
        duplicate_exons += len(m["exons"]) - len(exons)
        m["exons"] = exons
        m["cds"] = sorted(set(m["cds"]))
        if not exons:
            issue(0, "no_exons:" + tid)
        if len(exons) > 1 and m["strand"] == ".":
            issue(0, "unstranded_multi_exon:" + tid)
        for previous, current in zip(exons, exons[1:]):
            if previous[1] + 1 >= current[0]:
                issue(0, "overlapping_or_adjacent_exons:" + tid)
        m["exon_length"] = sum(b - a + 1 for a, b in exons)
        m["cds_length"] = sum(b - a + 1 for a, b in m["cds"])
    summary = {
        "path": str(path),
        "transcripts": len(models),
        "genes": len({m["gene"] for m in models.values()}),
        "multi_exon": sum(len(m["exons"]) > 1 for m in models.values()),
        "single_exon": sum(len(m["exons"]) == 1 for m in models.values()),
        "unstranded_single_exon": sum(
            len(m["exons"]) == 1 and m["strand"] == "." for m in models.values()
        ),
        "unique_multi_exon_chains": len(
            {chain_key(m) for m in models.values() if len(m["exons"]) > 1}
        ),
        "duplicate_exon_rows_removed": duplicate_exons,
        "features": dict(features),
        "errors": len(errors),
        "header": comments,
        "contigs": dict(collections.Counter(m["chrom"] for m in models.values())),
    }
    return models, summary, errors


def chain_key(m):
    exons = m["exons"]
    if len(exons) < 2:
        return None
    edges = tuple((a[1], b[0]) for a, b in zip(exons, exons[1:]))
    if m["strand"] == "-":
        edges = edges[::-1]
    return m["chrom"], m["strand"], edges


def chain_id(key):
    return "IC_" + hashlib.sha256(json.dumps(key, separators=(",", ":")).encode()).hexdigest()[:20]


def fasta_lengths(path):
    lengths, alphabet, duplicates = {}, collections.Counter(), []
    current = None
    with optext(path) as handle:
        for line in handle:
            if line.startswith(">"):
                words = line[1:].split()
                if not words:
                    raise ValueError("Empty FASTA identifier")
                current = words[0]
                if current in lengths:
                    duplicates.append(current)
                lengths[current] = 0
            elif line.strip():
                if current is None:
                    raise ValueError("FASTA sequence before first header")
                sequence = "".join(line.split()).upper()
                lengths[current] += len(sequence)
                alphabet.update(sequence)
    if not lengths:
        raise ValueError("No FASTA records")
    return lengths, dict(alphabet), duplicates


def check_bounds(models, lengths):
    problems = []
    for tid, m in models.items():
        if m["chrom"] not in lengths:
            problems.append((tid, m["chrom"], "contig_not_in_genome"))
        elif any(b > lengths[m["chrom"]] for a, b in m["exons"]):
            problems.append((tid, m["chrom"], "exon_outside_genome"))
    return problems


def run_dir(work, stage):
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%S_%fZ")
    result = work / "reports" / (stage + "_" + stamp)
    result.mkdir(parents=True)
    return result


def require_file(path):
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError("Missing/empty input: " + str(path))
    return path


def provenance(paths):
    return {str(p): {"size": Path(p).stat().st_size, "sha256": digest(p)} for p in paths}


def probe(cfg):
    print("Python:", sys.version.split()[0])
    if sys.version_info < (3, 8):
        raise ValueError("Python >=3.8 required")
    for key in ("query_gtf", "provided_fasta", "genome_6a", "reference_6a"):
        path = Path(cfg[key])
        print(
            key,
            "OK" if path.is_file() and os.access(path, os.R_OK) else "MISSING/UNREADABLE",
            str(path),
            sep="\t",
        )
    print("Scope:", cfg["query_scope"])
    print("This probe does not install packages or submit jobs.")


def audit(cfg, work):
    out = run_dir(work, "audit")
    print("REPORT_DIR=" + str(out), flush=True)
    paths = [
        require_file(cfg[k]) for k in ("query_gtf", "provided_fasta", "genome_6a", "reference_6a")
    ]
    query, qs, qe = parse_gtf(paths[0])
    ref, rs, re_ = parse_gtf(paths[3])
    print("GTF parsed; scanning genome FASTA (read-only).", flush=True)
    genome, _, gd = fasta_lengths(paths[2])
    fa, alphabet, duplicates = fasta_lengths(paths[1])
    qb, rb = check_bounds(query, genome), check_bounds(ref, genome)
    common = set(query) & set(fa)
    lengths_equal = sum(fa[t] == query[t]["exon_length"] for t in common)
    cds_equal = sum(query[t]["cds_length"] > 0 and fa[t] == query[t]["cds_length"] for t in common)
    report = {
        "query": qs,
        "reference": rs,
        "query_scope": cfg["query_scope"],
        "sequence_audit": {
            "fasta_records": len(fa),
            "exact_ID_overlap": len(common),
            "query_IDs_missing_in_fasta": len(set(query) - set(fa)),
            "fasta_IDs_not_in_query": len(set(fa) - set(query)),
            "equal_exon_sum_length": lengths_equal,
            "equal_CDS_sum_length": cds_equal,
            "alphabet": alphabet,
            "duplicate_fasta_IDs": len(duplicates),
            "interpretation": "Length/ID agreement is diagnostic only; full-length sequence identity is not established.",
        },
        "query_coordinate_problems": len(qb),
        "reference_coordinate_problems": len(rb),
        "duplicate_genome_IDs": len(gd),
        "input_provenance": provenance(paths),
    }
    fatal = bool(qe or re_ or qb or rb or gd or duplicates)
    report["input_validity"] = "FAIL" if fatal else "PASS_WITH_SCOPE_REVIEW"
    dump(out / "audit_summary.json", report)
    tab(
        out / "gtf_errors.tsv",
        ["resource", "line", "reason"],
        [("query", *x) for x in qe] + [("ref6a", *x) for x in re_],
    )
    tab(
        out / "coordinate_problems.tsv",
        ["resource", "transcript_id", "contig", "reason"],
        [("query", *x) for x in qb] + [("ref6a", *x) for x in rb],
    )
    tab(out / "genome_contigs.tsv", ["contig", "length"], sorted(genome.items()))
    tab(
        out / "query_fasta_correspondence.tsv",
        [
            "transcript_id",
            "gene_id",
            "exon_sum_length",
            "CDS_sum_length",
            "fasta_length",
            "exact_ID_found",
        ],
        (
            (t, m["gene"], m["exon_length"], m["cds_length"], fa.get(t, "NA"), t in fa)
            for t, m in sorted(query.items())
        ),
    )
    tab(
        out / "fasta_only_IDs.tsv",
        ["fasta_id", "length"],
        ((t, fa[t]) for t in sorted(set(fa) - set(query))),
    )
    tab(
        out / "duplicate_fasta_IDs.tsv",
        ["resource", "id"],
        [("provided", x) for x in duplicates] + [("genome", x) for x in gd],
    )
    notes = [
        "Status: " + report["input_validity"],
        "Query models: " + str(len(query)),
        "FASTA records: " + str(len(fa)),
        "Exact ID overlap: " + str(len(common)),
        "Equal exon-sum length: " + str(lengths_equal),
        "Equal CDS-sum length (where GTF CDS exists): " + str(cds_equal),
        "Provided FASTA is NOT used for chain comparison or liftover.",
        "Confirm whether union.topo.gtf is the final manuscript universe before reporting numbers.",
        "Send audit_summary.json, gtf_errors.tsv, coordinate_problems.tsv and the log back for review.",
    ]
    (out / "READ_ME.txt").write_text("\n".join(notes) + "\n", encoding="utf-8")
    print("\n".join(notes), flush=True)
    if fatal:
        raise ValueError("Audit failed: inspect the report files before comparing.")


def validate_download(path, kind):
    if kind == "zip":
        with zipfile.ZipFile(path) as z:
            bad = z.testzip()
            if bad:
                raise ValueError("ZIP CRC failure: " + bad)
    elif kind == "gtf_gz":
        found = False
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            for line in handle:
                if not line.startswith("#") and line.strip():
                    f = line.rstrip().split("\t")
                    if len(f) != 9:
                        raise ValueError("Downloaded content is not a GTF")
                    found |= f[2] == "exon" and "transcript_id" in f[8]
        if not found:
            raise ValueError("No GTF exon records found")
    else:
        raise ValueError("Unknown download kind")


def download_one(item, folder):
    final = folder / item["filename"]
    record_path = folder / (item["id"] + ".download.json")
    if final.exists():
        if not record_path.exists():
            raise ValueError(
                "Existing file lacks provenance; move it aside before retry: " + str(final)
            )
        old = json.loads(record_path.read_text(encoding="utf-8"))
        if old.get("url") != item["url"] or old.get("sha256") != digest(final):
            raise ValueError("Existing download differs from its manifest: " + str(final))
        validate_download(final, item["kind"])
        print("Verified existing:", final.name, flush=True)
        return final
    part = final.with_name(final.name + ".part")
    for attempt in range(1, 4):
        try:
            print("Downloading:", item["id"], "attempt", attempt, flush=True)
            req = urllib.request.Request(
                item["url"], headers={"User-Agent": "ChickenAnnotationCompare/1.0"}
            )
            with urllib.request.urlopen(req, timeout=45) as response, open(part, "wb") as output:
                resolved = response.geturl()
                headers = dict(response.headers)
                shutil.copyfileobj(response, output, length=1024 * 1024)
            validate_download(part, item["kind"])
            record = dict(
                item,
                resolved_url=resolved,
                bytes=part.stat().st_size,
                sha256=digest(part),
                retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                http_headers={
                    k: v
                    for k, v in headers.items()
                    if k.lower() in ("etag", "last-modified", "content-length", "content-type")
                },
            )
            os.replace(part, final)
            dump(record_path, record)
            return final
        except Exception as exc:
            print("Download failed:", type(exc).__name__, str(exc), flush=True)
            if attempt == 3:
                raise
            time.sleep(2)


def unpack_gtf(archive, target):
    # Never extract archive paths to disk. Only copy exactly one .gtf payload
    # into the fixed target, supporting a single nested supplement ZIP.
    candidates = []
    with zipfile.ZipFile(archive) as outer:
        for info in outer.infolist():
            if info.is_dir() or info.filename.startswith("__MACOSX/"):
                continue
            if info.filename.lower().endswith(".gtf"):
                candidates.append((info.filename, None))
            elif info.filename.lower().endswith(".zip"):
                if info.file_size > 128 * 1024 * 1024:
                    raise ValueError("Nested archive exceeds safety limit")
                with zipfile.ZipFile(io.BytesIO(outer.read(info))) as inner:
                    candidates.extend(
                        (info.filename, x.filename)
                        for x in inner.infolist()
                        if not x.is_dir()
                        and x.filename.lower().endswith(".gtf")
                        and not x.filename.startswith("__MACOSX/")
                    )
        if len(candidates) != 1:
            raise ValueError("Expected one Guan GTF; found: " + repr(candidates))
        first, second = candidates[0]
        part = target.with_name(target.name + ".part")
        if second is None:
            with outer.open(first) as source, open(part, "wb") as output:
                shutil.copyfileobj(source, output)
        else:
            with zipfile.ZipFile(io.BytesIO(outer.read(first))) as inner:
                with inner.open(second) as source, open(part, "wb") as output:
                    shutil.copyfileobj(source, output)
    models, summary, errors = parse_gtf(part)
    if errors or not summary["multi_exon"]:
        raise ValueError("Guan GTF failed validation: " + repr(errors[:5]))
    os.replace(part, target)
    return {
        "archive_member": first,
        "nested_member": second,
        "sha256": digest(target),
        "gtf_summary": summary,
    }


def downloads(work):
    folder = work / "downloads"
    folder.mkdir(parents=True, exist_ok=True)
    items = json.loads((ROOT / "study_sources.json").read_text(encoding="utf-8"))
    results, errors = [], []
    for item in items:
        try:
            path = download_one(item, folder)
            if item["id"] == "guan2022":
                details = unpack_gtf(path, folder / "Guan2022_GRCg6a.gtf")
                details["source_archive_sha256"] = digest(path)
                dump(folder / "guan2022.extraction.json", details)
                print(
                    "Guan GTF transcripts:",
                    details["gtf_summary"]["transcripts"],
                    "(paper reports 55382)",
                    flush=True,
                )
            else:
                models, summary, problems = parse_gtf(path)
                dump(folder / "france2024.gtf_audit.json", dict(summary, sha256=digest(path)))
                tab(folder / "france2024.gtf_issues.tsv", ["line", "reason"], problems)
                if problems:
                    print(
                        "France structural QC issues:",
                        len(problems),
                        "; preserved unchanged. Resolve before France comparison.",
                        flush=True,
                    )
                print(
                    "France GTF transcripts:", len(models), "genes:", summary["genes"], flush=True
                )
            results.append(item["id"])
        except Exception as exc:
            errors.append({"resource": item["id"], "error": str(exc)})
    dump(
        folder / "download_summary.json",
        {
            "completed": results,
            "errors": errors,
            "france_version_note": "Website counts and paper counts differ; retain SHA256 and audit actual models before final publication comparison.",
        },
    )
    print("DOWNLOAD_DIR=" + str(folder), flush=True)
    if errors:
        raise ValueError(
            "Some downloads failed; see download_summary.json. Do not classify missing resources as absent."
        )


def native(cfg, work):
    out = run_dir(work, "native6a")
    print("REPORT_DIR=" + str(out), flush=True)
    query_path = require_file(cfg["query_gtf"])
    ref_path = require_file(cfg["reference_6a"])
    genome_path = require_file(cfg["genome_6a"])
    guan_path = require_file(work / "downloads" / "Guan2022_GRCg6a.gtf")
    query, summary, errors = parse_gtf(query_path)
    if errors:
        raise ValueError("Invalid query GTF; run audit first: " + repr(errors[:5]))
    genome, _, duplicates = fasta_lengths(genome_path)
    if duplicates or check_bounds(query, genome):
        raise ValueError("Query/genome coordinate mismatch; run audit first")
    sources, ref_summaries = {}, {}
    for label, path in [("Ensembl106_6a", ref_path), ("Guan2022_6a", guan_path)]:
        models, stats, issues = parse_gtf(path)
        bounds = check_bounds(models, genome)
        if issues or bounds:
            tab(out / (label + "_errors.tsv"), ["reason"], [[str(x)] for x in issues + bounds])
            raise ValueError(
                "Reference failed validation: "
                + label
                + ". No contig aliases applied automatically."
            )
        index = collections.defaultdict(list)
        for tid, m in models.items():
            key = chain_key(m)
            if key is not None:
                index[key].append(tid)
        sources[label] = index
        ref_summaries[label] = stats
    groups = collections.defaultdict(list)
    singles = []
    for tid, m in query.items():
        key = chain_key(m)
        if key is None:
            singles.append((tid, m["gene"], m["chrom"], m["strand"], m["exons"]))
        else:
            groups[key].append(tid)
    counts = collections.Counter()
    chain_rows, model_rows = [], []
    for key in sorted(groups):
        tids = sorted(groups[key])
        ref_hits = sorted(sources["Ensembl106_6a"].get(key, []))
        guan_hits = sorted(sources["Guan2022_6a"].get(key, []))
        category = (
            "KNOWN_ENSEMBL106"
            if ref_hits
            else "GUAN_MATCH_NOT_ENSEMBL106" if guan_hits else "UNMATCHED_IN_TWO_6A_RESOURCES_ONLY"
        )
        counts[category] += 1
        cid = chain_id(key)
        genes = sorted({query[t]["gene"] for t in tids})
        chain_rows.append(
            (
                cid,
                len(tids),
                ";".join(tids),
                ";".join(genes),
                len(genes) > 1,
                key[0],
                key[1],
                json.dumps(key[2]),
                ";".join(ref_hits) or "NA",
                ";".join(guan_hits) or "NA",
                category,
                "PENDING",
                "PENDING",
                "PENDING",
            )
        )
        model_rows.extend((t, query[t]["gene"], cid, category) for t in tids)
    tab(
        out / "chain_match_6a.tsv",
        [
            "chain_id",
            "n_query_models",
            "query_ids",
            "query_gene_ids",
            "gene_ID_conflict",
            "contig",
            "strand",
            "flanking_exon_boundaries",
            "Ensembl106_matching_ids",
            "Guan_matching_ids",
            "native6a_category",
            "GRCg7b",
            "GRCg7w",
            "France7b",
        ],
        chain_rows,
    )
    tab(
        out / "query_model_to_chain.tsv",
        ["transcript_id", "gene_id", "chain_id", "native6a_category"],
        model_rows,
    )
    tab(
        out / "single_exon_separate.tsv",
        ["transcript_id", "gene_id", "contig", "strand", "exons"],
        singles,
    )
    assert len(model_rows) + len(singles) == len(query)
    assert sum(counts.values()) == len(groups)
    dump(
        out / "native_summary.json",
        {
            "scope": "two_6a_resources_only_not_final_novelty",
            "query_scope": cfg["query_scope"],
            "query": summary,
            "references": ref_summaries,
            "chain_categories": dict(counts),
            "chain_total": len(groups),
            "single_exon_models_separate": len(singles),
            "conservation_checks": "PASS",
            "unperformed": ["7b", "7w", "France7b", "terminal_variants", "single_exon_overlap"],
            "input_provenance": provenance([query_path, ref_path, genome_path, guan_path]),
        },
    )
    print(json.dumps(dict(counts), indent=2), flush=True)
    print(
        "These are preliminary structural comparisons, NOT final manuscript novelty counts.",
        flush=True,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(ROOT / "config.json"))
    parser.add_argument("command", choices=["probe", "download", "audit", "native"])
    args = parser.parse_args()
    cfg_path = Path(args.config).resolve()
    cfg = json.loads(cfg_path.read_text(encoding="utf-8-sig"))
    work = Path(cfg["work_dir"])
    if not work.is_absolute():
        work = cfg_path.parent / work
    if args.command == "probe":
        probe(cfg)
    elif args.command == "download":
        downloads(work)
    elif args.command == "audit":
        audit(cfg, work)
    else:
        native(cfg, work)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print("ERROR: " + str(exc), file=sys.stderr, flush=True)
        sys.exit(2)
