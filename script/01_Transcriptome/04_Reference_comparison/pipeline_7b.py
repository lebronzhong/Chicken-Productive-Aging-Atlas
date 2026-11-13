#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 01_Transcriptome
# Purpose: Pipeline 7b.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

"""Liftoff pilot/full transfer with per-transcript and reciprocal-coordinate QC."""

import argparse
import collections
import datetime
import functools
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from urllib.parse import unquote
import gtf_utils as u

ROOT = Path(__file__).resolve().parent


def read_json(p):
    return json.loads(Path(p).read_text(encoding="utf-8-sig"))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def config():
    c = read_json(ROOT / "config_7b.json")
    for k in ("assembly_report", "work_dir"):
        p = Path(c[k])
        c[k] = str(p if p.is_absolute() else ROOT / p)
    return c


def stamp():
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%S_%fZ")


def checked_gtf(path):
    m, s, e = u.parse_gtf(path)
    require(m and not e, "Invalid GTF: " + str(path) + "; " + repr(e[:10]))
    return m, s


def strict_sub(a, b):
    return len(a) < len(b) and any(b[i : i + len(a)] == a for i in range(len(b) - len(a) + 1))


def assembly_aliases(report, lengths):
    aliases = collections.defaultdict(set)
    verified = set()
    for line in Path(report).read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        require(len(f) == 10, "Unexpected assembly report columns")
        name, role, molecule, kind, genbank, relation, refseq, unit, length, ucsc = f
        candidates = {x for x in (name, genbank, refseq, ucsc) if x != "na"}
        if role == "assembled-molecule":
            candidates.update((molecule, "chr" + molecule))
            if molecule in ("MT", "M"):
                candidates.update(("MT", "M", "chrM", "chrMT"))
        candidates.update("chr" + x for x in list(candidates) if not x.startswith("chr"))
        present = candidates & set(lengths)
        if not present:
            continue
        require(
            len(present) == 1,
            "Multiple representations of same sequence in FASTA: " + repr(present),
        )
        target = next(iter(present))
        require(lengths[target] == int(length), "Assembly sequence length mismatch: " + target)
        verified.add(target)
        # Only equate RefSeq/GenBank coordinates when the report says identical.
        if relation != "=":
            candidates = {target}
        if role == "assembled-molecule" and relation == "=":
            candidates.update((molecule, "chr" + molecule))
            if molecule in ("MT", "M"):
                candidates.update(("MT", "M", "chrM", "chrMT"))
        for a in list(candidates):
            candidates.add("chr" + a)
        for a in candidates:
            aliases[a].add(target)
    require(
        set(lengths) == verified,
        "FASTA contains sequences absent from the expected assembly report: "
        + repr(sorted(set(lengths) - verified)[:10]),
    )
    return aliases


def normalize_gtf(source, target, aliases, lengths):
    mapping = {}
    with u.optext(source) as f, open(target, "w", encoding="utf-8", newline="\n") as out:
        for line in f:
            if line.startswith("#") or not line.strip():
                out.write(line)
                continue
            fields = line.rstrip("\r\n").split("\t")
            require(len(fields) == 9, "Invalid reference GTF row")
            old = fields[0]
            choices = {old} if old in lengths else aliases.get(old, set())
            require(len(choices) == 1, "Unresolved/ambiguous reference contig: " + old)
            new = next(iter(choices))
            require(
                1 <= int(fields[3]) <= int(fields[4]) <= lengths[new],
                "Reference GTF coordinate outside assembly: " + line[:100],
            )
            mapping[old] = new
            fields[0] = new
            out.write("\t".join(fields) + "\n")
    return mapping


def gene_groups(models):
    groups = collections.defaultdict(list)
    for tid, m in models.items():
        groups[m["gid"]].append(tid)
    return groups


def write_gff(models, path):
    with open(path, "w", encoding="utf-8", newline="\n") as out:
        out.write("##gff-version 3\n")
        for gid, tids in sorted(gene_groups(models).items()):
            ms = [models[t] for t in tids]
            require(
                len({(m["chrom"], m["strand"]) for m in ms}) == 1,
                "Gene spans multiple contigs/strands: " + gid,
            )
            chrom, strand = ms[0]["chrom"], ms[0]["strand"]
            left = min(a for m in ms for a, b in m["exons"])
            right = max(b for m in ms for a, b in m["exons"])

            def emit(ft, a, b, attrs):
                out.write(
                    f"{chrom}\tannotation_transfer\t{ft}\t{a}\t{b}\t.\t{strand}\t.\t{attrs}\n"
                )

            emit("gene", left, right, "ID=" + gid)
            for tid in sorted(tids):
                m = models[tid]
                emit(
                    "mRNA",
                    min(a for a, b in m["exons"]),
                    max(b for a, b in m["exons"]),
                    f"ID={tid};Parent={gid}",
                )
                for eid, (a, b) in zip(m["eids"], m["exons"]):
                    emit("exon", a, b, f"ID={eid};Parent={tid}")


def parse_lifted(path):
    genes, txs, exons, seen = {}, {}, [], set()
    with open(path, encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            x = line.rstrip().split("\t")
            require(len(x) == 9, "Invalid Liftoff GFF3 output")
            attrs = dict(
                (k, unquote(v)) for k, v in (a.split("=", 1) for a in x[8].split(";") if "=" in a)
            )
            fid = attrs.get("ID")
            require(fid and fid not in seen, "Duplicate/missing output feature ID: " + str(fid))
            seen.add(fid)
            obj = dict(chrom=x[0], strand=x[6], start=int(x[3]), end=int(x[4]), attrs=attrs)
            require(obj["start"] >= 1 and obj["end"] >= obj["start"], "Invalid lifted interval")
            if x[2] == "gene":
                genes[fid] = obj
            elif x[2] in ("mRNA", "transcript"):
                obj.update(gid=attrs.get("Parent"), exons=[], eids=[])
                txs[fid] = obj
            elif x[2] == "exon":
                exons.append((fid, attrs.get("Parent"), obj))
    for eid, parent, obj in exons:
        require(parent in txs, "Orphan/multiple-parent exon: " + str(parent))
        m = txs[parent]
        require(
            (obj["chrom"], obj["strand"]) == (m["chrom"], m["strand"]),
            "Exon/transcript coordinate inconsistency",
        )
        m["exons"].append((obj["start"], obj["end"]))
        m["eids"].append(eid)
    for m in txs.values():
        pairs = sorted(zip(m["exons"], m["eids"]))
        m["exons"] = [x for x, e in pairs]
        m["eids"] = [e for x, e in pairs]
    return genes, txs


def scan_sam(folder, expected):
    alignments = collections.defaultdict(dict)
    files = sorted(Path(folder).glob("*.sam"))
    require(files, "No Liftoff SAM files: uniqueness cannot be assessed")
    for path in files:
        with open(path) as f:
            for line in f:
                if line.startswith("@"):
                    continue
                x = line.rstrip().split("\t")
                require(len(x) >= 11, "Invalid SAM row")
                gid = x[0]
                if gid not in expected:
                    # Liftoff intermediate copy suffix; copies are not requested.
                    gid = (
                        gid.removesuffix("_0")
                        if hasattr(gid, "removesuffix")
                        else gid[:-2] if gid.endswith("_0") else gid
                    )
                require(gid in expected, "Unexpected SAM query ID: " + x[0])
                flag = int(x[1])
                if flag & 4:
                    continue
                span = sum(
                    int(n) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", x[5]) if op in "MDN=X"
                )
                tags = {p.split(":", 2)[0]: p.split(":", 2)[2] for p in x[11:] if p.count(":") >= 2}
                if "AS" not in tags:
                    continue
                a = dict(
                    chrom=x[2],
                    start=int(x[3]),
                    end=int(x[3]) + span - 1,
                    reverse=bool(flag & 16),
                    primary=not bool(flag & (256 | 2048)),
                    mapq=int(x[4]),
                    score=int(tags["AS"]),
                )
                key = a["chrom"], a["start"], a["end"], a["reverse"], a["primary"]
                prev = alignments[gid].get(key)
                if prev is None or (a["score"], a["mapq"]) > (prev["score"], prev["mapq"]):
                    alignments[gid][key] = a
    return {g: list(v.values()) for g, v in alignments.items()}


def distinct_locus(a, b):
    if (a["chrom"], a["reverse"]) != (b["chrom"], b["reverse"]):
        return True
    overlap = max(0, min(a["end"], b["end"]) - max(a["start"], b["start"]) + 1)
    return overlap / min(a["end"] - a["start"] + 1, b["end"] - b["start"] + 1) < 0.8


def gene_qc(gid, gene, source, aligns, c):
    if gene is None:
        return "GENE_UNMAPPED"
    a = gene["attrs"]
    try:
        if (
            float(a["coverage"]) < c["gene_min_coverage"]
            or float(a["sequence_ID"]) < c["gene_min_identity"]
        ):
            return "GENE_LOW_COVERAGE_OR_IDENTITY"
    except (KeyError, ValueError):
        return "GENE_MISSING_ALIGNMENT_METRICS"
    if (
        a.get("extra_copy_number", "0") != "0"
        or a.get("partial_mapping") == "True"
        or a.get("low_identity") == "True"
    ):
        return "GENE_FLAGGED_MAPPING"
    candidates = aligns.get(gid, [])
    required_reverse = source["strand"] != gene["strand"]
    primary = [
        x
        for x in candidates
        if x["primary"]
        and x["mapq"] >= c["min_primary_mapq"]
        and x["mapq"] != 255
        and x["chrom"] == gene["chrom"]
        and x["reverse"] == required_reverse
        and x["start"] <= gene["start"]
        and x["end"] >= gene["end"]
    ]
    if not primary:
        return "GENE_NO_CONFIDENT_PRIMARY_SPANNING_MAPPING"
    best = max(primary, key=lambda x: x["score"])
    if any(
        x["score"] >= best["score"] * c["secondary_score_ratio"] and distinct_locus(x, best)
        for x in candidates
    ):
        return "GENE_COMPETING_ALIGNMENT"
    return "PASS"


def shape_qc(src, target):
    if target is None:
        return "TRANSCRIPT_UNMAPPED"
    if target["gid"] != src["gid"] or not target["exons"]:
        return "PARENT_OR_EXON_MISSING"
    ids = src["eids"] if src["strand"] == target["strand"] else list(reversed(src["eids"]))
    if ids != target["eids"]:
        return "EXON_ID_COUNT_OR_ORDER_CHANGED"
    if any(a[1] + 1 >= b[0] for a, b in zip(target["exons"], target["exons"][1:])):
        return "INVALID_EXON_SPACING"
    return "PASS"


def reciprocal_ok(original, back):
    return (
        back is not None
        and shape_qc(original, back) == "PASS"
        and (original["chrom"], original["strand"]) == (back["chrom"], back["strand"])
        and [tuple(x) for x in original["exons"]] == [tuple(x) for x in back["exons"]]
    )


RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq):
    return seq.translate(RC)[::-1]


def sequence_checker(f6, f7, c):
    import edlib

    @functools.lru_cache(maxsize=200000)
    def exon(ch6, a, b, ch7, x, y, reverse):
        s = str(f6[ch6][a - 1 : b]).upper()
        t = str(f7[ch7][x - 1 : y]).upper()
        if reverse:
            t = revcomp(t)
        if (
            not s
            or not t
            or sum(z not in "ACGT" for z in s) / len(s) > c["max_ambiguous_base_fraction"]
            or sum(z not in "ACGT" for z in t) / len(t) > c["max_ambiguous_base_fraction"]
        ):
            return "AMBIGUOUS_EXON_SEQUENCE", 0.0
        ratio = len(t) / len(s)
        if not c["exon_min_length_ratio"] <= ratio <= c["exon_max_length_ratio"]:
            return "EXON_LENGTH_CHANGE", 0.0
        identity = 1 - edlib.align(s, t, mode="NW", task="distance")["editDistance"] / max(
            len(s), len(t)
        )
        return ("PASS" if identity >= c["exon_min_identity"] else "EXON_LOW_IDENTITY"), identity

    def motifs(m, fa):
        seq = []
        for a, b in zip(m["exons"], m["exons"][1:]):
            donor = str(fa[m["chrom"]][a[1] : a[1] + 2]).upper()
            acceptor = str(fa[m["chrom"]][b[0] - 3 : b[0] - 1]).upper()
            seq.append(donor + acceptor)
        return seq if m["strand"] == "+" else [revcomp(x) for x in reversed(seq)]

    def check(src, tgt):
        d = dict(zip(tgt["eids"], tgt["exons"]))
        minimum = 1.0
        for eid, (a, b) in zip(src["eids"], src["exons"]):
            x, y = d[eid]
            if tgt["chrom"] not in f7 or y > len(f7[tgt["chrom"]]):
                return "EXON_OUTSIDE_TARGET", 0.0
            status, identity = exon(
                src["chrom"], a, b, tgt["chrom"], x, y, src["strand"] != tgt["strand"]
            )
            minimum = min(minimum, identity)
            if status != "PASS":
                return status, minimum
        if c["require_preserved_splice_motifs"] and motifs(src, f6) != motifs(tgt, f7):
            return "SPLICE_MOTIF_CHANGED", minimum
        return "PASS", minimum

    return check


def probe(c):
    print("Python:", sys.version.split()[0])
    failed = []
    for k in (
        "query_gtf",
        "genome_6a",
        "genome_7b",
        "reference_7b",
        "reference_7b_refseq",
        "assembly_report",
    ):
        p = Path(c[k])
        ok = p.is_file() and os.access(p, os.R_OK)
        print(k, "OK" if ok else "MISSING/UNREADABLE", str(p), sep="\t")
        if not ok:
            failed.append(k)
    require(not failed, "Fix missing inputs before submitting: " + repr(failed))
    print(
        "No alignment has been run. Liftoff/minimap2 may be installed in the isolated environment next."
    )


def prepare(c):
    from pyfaidx import Fasta

    work = Path(c["work_dir"])
    work.mkdir(parents=True, exist_ok=True)
    out = work / ("prepare_" + stamp())
    out.mkdir()
    print("PREPARE_DIR=" + str(out), flush=True)
    prov = {
        k: dict(path=c[k], sha256=u.digest(u.require_file(c[k])))
        for k in (
            "query_gtf",
            "genome_6a",
            "genome_7b",
            "reference_7b",
            "reference_7b_refseq",
            "assembly_report",
        )
    }
    for k in ("query", "genome_6a", "assembly_report"):
        pk = "query_gtf" if k == "query" else k
        require(
            prov[pk]["sha256"] == c["expected_" + k + "_sha256"],
            "Input differs from audited version: " + pk,
        )
    models, qs = checked_gtf(c["query_gtf"])
    require(
        qs["transcripts"] == c["expected_query_models"]
        and qs["unique_multi_exon_chains"] == c["expected_query_chains"],
        "Final query universe changed",
    )
    require(qs["single_exon"] == 0, "This configuration expects the audited multi-exon-only query")
    for k, name in [("genome_6a", "genome6a.fa"), ("genome_7b", "genome7b.fa")]:
        print("Staging private genome copy:", name, flush=True)
        shutil.copyfile(c[k], out / name)
    f6 = Fasta(str(out / "genome6a.fa"))
    f7 = Fasta(str(out / "genome7b.fa"))
    len6 = {k: len(f6[k]) for k in f6.keys()}
    len7 = {k: len(f7[k]) for k in f7.keys()}
    require(not u.check_bounds(models, len6), "Query/source genome coordinate mismatch")
    aliases = assembly_aliases(c["assembly_report"], len7)
    refstats = {}
    alias_counts = {}
    for label, field in [("Ensembl114", "reference_7b"), ("RefSeq7b", "reference_7b_refseq")]:
        normalized = out / (label + ".normalized.gtf")
        amap = normalize_gtf(c[field], normalized, aliases, len7)
        refs, rs = checked_gtf(normalized)
        build_headers = [x for x in rs["header"] if x.startswith("#!genome-build ")]
        require(
            not build_headers or all("GRCg7b" in x for x in build_headers),
            "Reference header reports an assembly other than GRCg7b: " + label,
        )
        require(
            not u.check_bounds(refs, len7), "Target annotation/genome coordinate mismatch: " + label
        )
        u.tab(
            out / (label + ".contig_aliases.tsv"),
            ["annotation_contig", "FASTA_contig"],
            sorted(amap.items()),
        )
        refstats[label] = rs
        alias_counts[label] = len(amap)
    genes = {g: "G%08d" % i for i, g in enumerate(sorted({m["gene"] for m in models.values()}), 1)}
    source = {}
    for i, (tid, m) in enumerate(sorted(models.items()), 1):
        sid = "T%09d" % i
        source[sid] = dict(
            original_id=tid,
            original_gene=m["gene"],
            gid=genes[m["gene"]],
            chrom=m["chrom"],
            strand=m["strand"],
            exons=m["exons"],
            eids=[sid + "E%04d" % j for j in range(1, len(m["exons"]) + 1)],
        )
    groups = gene_groups(source)
    # First represent every contig, then complex loci, then deterministic hash sampling.
    bycontig = collections.defaultdict(list)
    for gid, tids in groups.items():
        bycontig[source[tids[0]]["chrom"]].append(gid)
    rank = lambda g: hashlib.sha256(g.encode()).hexdigest()
    selected = {min(gs, key=rank) for gs in bycontig.values()}
    complex_genes = sorted(
        groups, key=lambda g: (-max(len(source[t]["exons"]) for t in groups[g]), g)
    )[:20]
    selected.update(complex_genes)
    target = min(len(groups), max(c["pilot_genes"], len(selected)))
    for g in sorted(groups, key=rank):
        if len(selected) >= target:
            break
        selected.add(g)
    pilot = {t: m for t, m in source.items() if m["gid"] in selected}
    write_gff(source, out / "all.source.gff3")
    write_gff(pilot, out / "pilot.source.gff3")
    u.dump(out / "source_models.json", source)
    u.dump(out / "pilot_ids.json", sorted(pilot))
    meta = dict(
        config=c,
        input_provenance=prov,
        query=qs,
        references_7b=refstats,
        pilot_genes=len(selected),
        pilot_transcripts=len(pilot),
        alias_counts=alias_counts,
        status="PREPARED",
        reference_scope="Ensembl114 and supplied RefSeq7b; RefSeq release must be read from header/metadata",
        generated_hashes={
            n: u.digest(out / n)
            for n in (
                "source_models.json",
                "pilot_ids.json",
                "all.source.gff3",
                "pilot.source.gff3",
                "Ensembl114.normalized.gtf",
                "RefSeq7b.normalized.gtf",
            )
        },
    )
    u.dump(out / "prepare_summary.json", meta)
    f6.close()
    f7.close()
    u.dump(work / "prepared_pointer.json", dict(directory=str(out)))
    print(
        "STATUS=PREPARED; pilot genes=%d; pilot transcripts=%d" % (len(selected), len(pilot)),
        flush=True,
    )


def lift(annotation, reference, target, out, threads, c):
    out.mkdir()
    cmd = [
        "liftoff",
        str(target),
        str(reference),
        "-g",
        str(annotation),
        "-o",
        str(out / "lifted.gff3"),
        "-u",
        str(out / "unmapped.txt"),
        "-dir",
        str(out / "intermediate"),
        "-p",
        str(threads),
        "-a",
        str(c["gene_min_coverage"]),
        "-s",
        str(c["gene_min_identity"]),
        "-exclude_partial",
    ]
    u.dump(out / "command.json", cmd)
    execution = (
        ["/usr/bin/time", "-v", "-o", str(out / "resources.txt")] + cmd
        if Path("/usr/bin/time").is_file()
        else cmd
    )
    with open(out / "liftoff.log", "w") as log:
        subprocess.run(execution, stdout=log, stderr=subprocess.STDOUT, check=True)
    u.require_file(out / "lifted.gff3")
    return parse_lifted(out / "lifted.gff3")


def compare_target(source, passed, refs, out, stage, reference_label="Ensembl114"):
    ref_index = collections.defaultdict(list)
    for tid, m in refs.items():
        k = u.chain_key(m)
        if k:
            ref_index[k].append(tid)
    source_groups = collections.defaultdict(list)
    for tid, m in source.items():
        source_groups[u.chain_key(m)].append(tid)
    counts = collections.Counter()
    rows = []
    eligible_target = set()
    for key, tids in sorted(source_groups.items()):
        ok = [t for t in tids if t in passed]
        targets = {u.chain_key(passed[t]) for t in ok}
        if len(ok) != len(tids):
            status = "UNRESOLVED_MODEL_QC"
        elif len(targets) != 1:
            status = "UNRESOLVED_DIVERGENT_PROJECTION"
        else:
            target = next(iter(targets))
            eligible_target.add(target)
            status = "MATCH" if target in ref_index else "NO_EXACT_MATCH_IN_" + reference_label
        counts[status] += 1
        rows.append(
            (
                u.chain_id(key),
                ";".join(source[t]["original_id"] for t in tids),
                len(tids),
                len(ok),
                status,
                ";".join(sorted({r for k in targets for r in ref_index.get(k, [])})) or "NA",
            )
        )
    u.tab(
        out / "source_chain_comparison.tsv",
        [
            "source_chain_id",
            "original_transcript_ids",
            "n_models",
            "n_QC_pass_models",
            "status",
            "target_reference_ids_from_passing_models",
        ],
        rows,
    )
    # Compute coverage in the target coordinate system, after deduplicating projected chains.
    contains = collections.defaultdict(list)
    starts = collections.defaultdict(list)
    for ch, st, e in eligible_target:
        starts[ch, st, e[0]].append(e)
        for junction in set(e):
            contains[ch, st, junction].append(e)
    refcounts = collections.Counter()
    reflines = []
    for key, ids in sorted(ref_index.items()):
        ch, st, e = key
        if key in eligible_target:
            cat = "EXACT"
        elif any(strict_sub(e, x) for x in contains.get((ch, st, e[0]), ())):
            cat = "REFERENCE_SUBCHAIN_OF_QUERY"
        elif any(strict_sub(x, e) for j in e for x in starts.get((ch, st, j), ())):
            cat = "QUERY_SUBCHAIN_OF_REFERENCE"
        else:
            cat = "NO_EXACT_OR_SUBCHAIN_MATCH"
        refcounts[cat] += 1
        reflines.append((u.chain_id(key), ";".join(sorted(ids)), cat))
    u.tab(
        out / "reference_chain_coverage.tsv",
        ["target_reference_chain_id", "reference_transcript_ids", "status"],
        reflines,
    )
    assert sum(counts.values()) == len(source_groups)
    assert sum(refcounts.values()) == len(ref_index)
    total = len(ref_index)
    exact = refcounts["EXACT"]
    part = refcounts["REFERENCE_SUBCHAIN_OF_QUERY"] + refcounts["QUERY_SUBCHAIN_OF_REFERENCE"]
    return dict(
        stage=stage,
        reference_label=reference_label,
        scope="PILOT_NOT_FINAL" if stage == "pilot" else reference_label + "_QC_pass_scope",
        source_chains=len(source_groups),
        source_chain_categories=dict(counts),
        eligible_distinct_target_chains=len(eligible_target),
        reference_multi_exon_chains=total,
        reference_chain_categories=dict(refcounts),
        reference_exact_percent=100 * exact / total if total else None,
        reference_exact_or_subchain_percent=100 * (exact + part) / total if total else None,
        denominator_note="Reference denominator: all deduplicated multi-exon chains in this reference. Query numerator: source chains whose ALL transcript models passed QC and projected consistently. Missing mapping is unresolved, never novelty.",
    )


def run(c, stage, threads):
    from pyfaidx import Fasta

    work = Path(c["work_dir"])
    prep = Path(read_json(work / "prepared_pointer.json")["directory"])
    meta = read_json(prep / "prepare_summary.json")
    require(
        meta["config"] == c,
        "Configuration changed after prepare; prepare again to freeze new inputs",
    )
    for n, h in meta["generated_hashes"].items():
        require(u.digest(prep / n) == h, "Prepared file changed: " + n)
    for name, k in [("genome6a.fa", "genome_6a"), ("genome7b.fa", "genome_7b")]:
        require(
            u.digest(prep / name) == meta["input_provenance"][k]["sha256"],
            "Staged genome changed: " + name,
        )
    source = read_json(prep / "source_models.json")
    if stage == "pilot":
        subset = set(read_json(prep / "pilot_ids.json"))
        source = {t: m for t, m in source.items() if t in subset}
    out = work / (stage + "_" + stamp())
    out.mkdir()
    print("RUN_DIR=" + str(out), flush=True)
    versions = {}
    for tool, args in [("liftoff", ["-V"]), ("minimap2", ["--version"])]:
        versions[tool] = subprocess.check_output(
            [tool] + args, stderr=subprocess.STDOUT, text=True
        ).strip()
    u.dump(
        out / "run_manifest.json",
        dict(
            stage=stage,
            config=c,
            prepare_directory=str(prep),
            input_provenance=meta["input_provenance"],
            tools=versions,
        ),
    )
    annotation = prep / ("pilot.source.gff3" if stage == "pilot" else "all.source.gff3")
    print("Forward Liftoff: 6a -> 7b", flush=True)
    fg, ft = lift(
        annotation, prep / "genome6a.fa", prep / "genome7b.fa", out / "forward", threads, c
    )
    require(set(ft) <= set(source), "Unexpected transcript IDs in forward mapping")
    source_groups = gene_groups(source)
    fa = scan_sam(out / "forward/intermediate", set(source_groups))
    gs = {g: gene_qc(g, fg.get(g), source[ts[0]], fa, c) for g, ts in source_groups.items()}
    f6 = Fasta(str(prep / "genome6a.fa"))
    f7 = Fasta(str(prep / "genome7b.fa"))
    seqcheck = sequence_checker(f6, f7, c)
    reasons = {}
    minimum = {}
    forward_pass = {}
    for t, s in source.items():
        why = gs[s["gid"]]
        if why == "PASS":
            why = shape_qc(s, ft.get(t))
        if why == "PASS":
            why, minimum[t] = seqcheck(s, ft[t])
        if why == "PASS":
            forward_pass[t] = ft[t]
        else:
            reasons[t] = "FORWARD_" + why
    f6.close()
    f7.close()
    u.dump(
        out / "forward_qc_summary.json",
        dict(
            selected_models=len(source),
            forward_pass_models=len(forward_pass),
            failure_counts=dict(collections.Counter(reasons.values())),
        ),
    )
    print("Forward transcript QC pass:", len(forward_pass), "/", len(source), flush=True)
    passed = {}
    if forward_pass:
        backinput = out / "reciprocal_input.gff3"
        write_gff(forward_pass, backinput)
        print("Reciprocal Liftoff: 7b -> 6a", flush=True)
        bg, bt = lift(
            backinput, prep / "genome7b.fa", prep / "genome6a.fa", out / "reciprocal", threads, c
        )
        require(set(bt) <= set(forward_pass), "Unexpected transcript IDs in reciprocal mapping")
        backgroups = gene_groups(forward_pass)
        ba = scan_sam(out / "reciprocal/intermediate", set(backgroups))
        bgs = {g: gene_qc(g, bg.get(g), forward_pass[ts[0]], ba, c) for g, ts in backgroups.items()}
        for t, m in forward_pass.items():
            why = bgs[m["gid"]]
            if why == "PASS" and not reciprocal_ok(source[t], bt.get(t)):
                why = "EXON_COORDINATES_NOT_RESTORED"
            if why == "PASS":
                passed[t] = m
                reasons[t] = "PASS"
            else:
                reasons[t] = "RECIPROCAL_" + why
    require(len(reasons) == len(source), "Missing model QC decisions")
    with open(out / "query7b.QC_pass.gtf", "w", encoding="utf-8", newline="\n") as f:
        for t, m in sorted(passed.items()):
            s = source[t]
            attrs = (
                'gene_id "' + s["original_gene"] + '"; transcript_id "' + s["original_id"] + '";'
            )
            for a, b in m["exons"]:
                f.write(f'{m["chrom"]}\tLiftoff_QC\texon\t{a}\t{b}\t.\t{m["strand"]}\t.\t{attrs}\n')
    u.tab(
        out / "transcript_qc.tsv",
        [
            "original_transcript_id",
            "original_gene_id",
            "synthetic_id",
            "status",
            "min_exon_identity",
        ],
        (
            (s["original_id"], s["original_gene"], t, reasons[t], minimum.get(t, "NA"))
            for t, s in sorted(source.items())
        ),
    )
    reference_results = {}
    for label in ("Ensembl114", "RefSeq7b"):
        refs, rs = checked_gtf(prep / (label + ".normalized.gtf"))
        refout = out / label
        refout.mkdir()
        result = compare_target(source, passed, refs, refout, stage, label)
        result["reference_annotation"] = rs
        reference_results[label] = result
    report = dict(
        stage=stage,
        scope="PILOT_NOT_FINAL" if stage == "pilot" else "Ensembl114_and_RefSeq7b_QC_pass_scope",
        status="COMPLETED",
        selected_transcript_models=len(source),
        passed_transcript_models=len(passed),
        transcript_qc_counts=dict(collections.Counter(reasons.values())),
        mapping_confidence="MAPQ>=threshold, no comparable-scoring distinct alignment among reported SAM hits, per-exon sequence/structure QC, exact reciprocal exon coordinates. This is an operational confidence screen, not proof of unique orthology.",
        forward_gene_failure_counts=dict(collections.Counter(gs.values())),
        reference_results=reference_results,
        other_resources={"France7b": "PENDING", "7w": "PENDING"},
    )
    u.dump(out / "comparison_summary.json", report)
    print(
        json.dumps(
            {
                k: v
                for k, v in report.items()
                if k not in ("reference_annotation", "reference_results")
            },
            indent=2,
        ),
        flush=True,
    )
    for label, result in reference_results.items():
        print(
            label,
            json.dumps({k: v for k, v in result.items() if k != "reference_annotation"}),
            flush=True,
        )
    print(
        "Completed. Read comparison_summary.json and transcript_qc.tsv before interpreting unmatched chains.",
        flush=True,
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("stage", choices=["probe", "prepare", "pilot", "full"])
    ap.add_argument("--threads", type=int, default=1)
    args = ap.parse_args()
    require(args.threads > 0, "Threads must be positive")
    c = config()
    if args.stage == "probe":
        probe(c)
    elif args.stage == "prepare":
        prepare(c)
    else:
        run(c, args.stage, args.threads)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print("ERROR: " + str(e), file=sys.stderr, flush=True)
        raise
