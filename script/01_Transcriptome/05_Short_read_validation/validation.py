#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 01_Transcriptome
# Purpose: Short-read validation.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

import argparse, collections, csv, gzip, json, pickle, subprocess, os, hashlib, time
from pathlib import Path
import pysam
import gtf_utils as u

ROOT = Path(__file__).resolve().parent.parent
GENOME = Path(
    "/path/to/chicken-aging/reference/galGal6/Ensemble/Gallus_gallus.GRCg6a.dna_sm.toplevel.fa"
)
REF = GENOME.parent / "Gallus_gallus.GRCg6a.106.gtf"
QUERY = Path(
    "/path/to/chicken-aging/ont/shell/10_20240903.gffcompare/9.6_final_gtf/9.9_finalgtf_annotated_rename.gtf"
)
CLASS = (
    ROOT.parent
    / "annotation_comparison/query_oriented_six_resources_20260930T100949Z/all_query_chain_classification.tsv"
)
STAR = ROOT / "software/STAR.2.7.11b"
MANIFEST = ROOT / "sample_manifest.tsv"


def dump(path, x):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(x, indent=2, ensure_ascii=False))


def read_samples():
    return list(csv.DictReader(MANIFEST.open(), delimiter="\t"))


def run(cmd, log):
    print("RUN", *[str(x) for x in cmd], flush=True)
    with log.open("w") as f:
        subprocess.run([str(x) for x in cmd], stdout=f, stderr=subprocess.STDOUT, check=True)


def load():
    return pickle.load((ROOT / "prepared/database.pkl").open("rb"))


def prepare(expected_chains=None):
    out = ROOT / "prepared"
    out.mkdir(exist_ok=False)
    q, qs, err = u.parse_gtf(QUERY)
    assert not err
    r, rs, err = u.parse_gtf(REF)
    assert not err
    groups = collections.defaultdict(list)
    for tid, m in q.items():
        groups[u.chain_key(m)].append(tid)
    assert groups and None not in groups
    if expected_chains is not None:
        assert len(groups) == expected_chains
    classes = {
        x["source_chain_id"]: x["joint_category"]
        for x in csv.DictReader(CLASS.open(), delimiter="\t")
    }
    junctions = sorted({(ch, st, a, b) for ch, st, edges in groups for a, b in edges})
    ji = {key: i for i, key in enumerate(junctions)}
    reference_junctions = {
        (m["chrom"], m["strand"], a[1], b[0])
        for m in r.values()
        for a, b in zip(m["exons"], m["exons"][1:])
    }
    genome = pysam.FastaFile(str(GENOME))
    motifs = {}
    chains = {}
    for ch, st, a, b in junctions:
        key = (ch, a, b)
        if key not in motifs:
            motif = (
                genome.fetch(ch, a, a + 2).upper() + "-" + genome.fetch(ch, b - 3, b - 1).upper()
            )
            strand = (
                "+"
                if motif in ("GT-AG", "GC-AG", "AT-AC")
                else ("-" if motif in ("CT-AC", "CT-GC", "GT-AT") else "?")
            )
            motifs[key] = (motif, strand)
    genome.close()
    for key, tids in groups.items():
        cid = u.chain_id(key)
        ch, st, edges = key
        chains[cid] = {
            "key": key,
            "junction_ids": [ji[(ch, st, a, b)] for a, b in edges],
            "transcript_ids": tids,
            "category": classes[cid],
        }
    db = {
        "chains": chains,
        "junctions": junctions,
        "motifs": motifs,
        "reference_junctions": reference_junctions,
    }
    with (out / "database.pkl").open("wb") as f:
        pickle.dump(db, f, pickle.HIGHEST_PROTOCOL)
    dump(
        out / "provenance.json",
        {
            "query_sha256": u.digest(QUERY),
            "reference_sha256": u.digest(REF),
            "classification_sha256": u.digest(CLASS),
            "query_models": len(q),
            "chains": len(chains),
            "junctions": len(junctions),
            "library": "unstranded",
            "RNAseq_used_in_ONT_construction": False,
            "sr_count_unit": "primary uniquely mapped read alignments; mates may count separately; not UMI molecules",
            "short_read_anchor_min": 10,
            "ont_mapq_min": 20,
            "ONT_scope": "exact full intron-chain matches to primary MAPQ>=20 alignments; does not establish transcript endpoints",
            "thresholds": [1, 3, 5],
            "query_not_used_to_guide_STAR": True,
        },
    )
    print("PREPARED", len(chains), len(junctions), flush=True)


def index(threads):
    out = ROOT / "star_index"
    out.mkdir(exist_ok=False)
    run(
        [
            STAR,
            "--runMode",
            "genomeGenerate",
            "--runThreadN",
            threads,
            "--genomeDir",
            out,
            "--genomeFastaFiles",
            GENOME,
            "--sjdbGTFfile",
            REF,
            "--sjdbOverhang",
            "149",
            "--genomeSAindexNbases",
            "13",
        ],
        ROOT / "logs/index_STAR.log",
    )
    dump(
        out / "COMPLETE.json",
        {
            "status": "COMPLETE",
            "STAR": subprocess.check_output([str(STAR), "--version"], text=True).strip(),
            "annotation": "Ensembl106; no ONT query annotation",
        },
    )


def edges_and_anchors(read):
    pos = read.reference_start
    out = []
    ops = read.cigartuples or []
    for i, (op, length) in enumerate(ops):
        if op == 3:
            left = ops[i - 1][1] if i > 0 and ops[i - 1][0] in (0, 7, 8) else 0
            right = ops[i + 1][1] if i + 1 < len(ops) and ops[i + 1][0] in (0, 7, 8) else 0
            out.append((pos, pos + length + 1, left, right))
        if op in (0, 2, 3, 7, 8):
            pos += length
    return out


def ont(index, threads):
    sample = read_samples()[index]
    out = ROOT / "ont" / sample["ont_sample"]
    out.mkdir(parents=True, exist_ok=False)
    db = load()
    lookup = {v["key"]: cid for cid, v in db["chains"].items()}
    counter = collections.Counter()
    qc = collections.Counter()
    run(["samtools", "quickcheck", "-v", sample["ont_bam"]], out / "quickcheck.log")
    with pysam.AlignmentFile(sample["ont_bam"], "rb", threads=min(threads, 4)) as bam:
        for read in bam.fetch(until_eof=True):
            qc["alignment_records"] += 1
            if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_qcfail:
                continue
            qc["primary_mapped_reads"] += 1
            if read.mapping_quality < 20:
                continue
            ee = edges_and_anchors(read)
            if not ee:
                continue
            qc["primary_mapq20_spliced_reads"] += 1
            strand = "-" if read.is_reverse else "+"
            edges = tuple((a, b) for a, b, l, r in ee)
            if strand == "-":
                edges = edges[::-1]
            cid = lookup.get((read.reference_name, strand, edges))
            if cid:
                counter[cid] += 1
                qc["exact_query_chain_reads"] += 1
    with gzip.open(out / "chain_counts.tsv.gz", "wt") as f:
        f.write("source_chain_id\tONT_full_chain_reads\n")
        for cid, num in sorted(counter.items()):
            f.write(f"{cid}\t{num}\n")
    dump(
        out / "COMPLETE.json",
        dict(
            qc,
            sample=sample["ont_sample"],
            distinct_detected_chains=len(counter),
            status="COMPLETE",
        ),
    )
    print("ONT_COMPLETE", sample["ont_sample"], dict(qc), flush=True)


def rna(index, threads):
    sample = read_samples()[index]
    out = ROOT / "rnaseq" / sample["ont_sample"]
    out.mkdir(parents=True, exist_ok=False)
    expected = {
        parts[1].lstrip("*"): parts[0]
        for line in (Path(sample["fastq_r1"]).parent / "data_md5.txt").read_text().splitlines()
        if len(parts := line.split()) == 2
    }
    checked = {}
    for key in ("fastq_r1", "fastq_r2"):
        path = Path(sample[key])
        h = hashlib.md5()
        with path.open("rb") as f:
            for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
                h.update(block)
        checked[path.name] = h.hexdigest()
        assert h.hexdigest() == expected[path.name], ("FASTQ_MD5_MISMATCH", path)
    dump(out / "raw_fastq_md5.json", checked)
    run(
        [
            "fastp",
            "-i",
            sample["fastq_r1"],
            "-I",
            sample["fastq_r2"],
            "-o",
            out / "clean_R1.fq.gz",
            "-O",
            out / "clean_R2.fq.gz",
            "--thread",
            min(threads, 16),
            "--detect_adapter_for_pe",
            "--length_required",
            "50",
            "--json",
            out / "fastp.json",
            "--html",
            out / "fastp.html",
        ],
        out / "fastp.log",
    )
    assert (ROOT / "star_index/COMPLETE.json").exists()
    prefix = str(out / "STAR_")
    run(
        [
            STAR,
            "--genomeDir",
            ROOT / "star_index",
            "--runThreadN",
            threads,
            "--readFilesIn",
            out / "clean_R1.fq.gz",
            out / "clean_R2.fq.gz",
            "--readFilesCommand",
            "zcat",
            "--twopassMode",
            "Basic",
            "--outFileNamePrefix",
            prefix,
            "--outSAMtype",
            "BAM",
            "SortedByCoordinate",
            "--limitBAMsortRAM",
            "24000000000",
            "--outSAMattributes",
            "NH",
            "HI",
            "AS",
            "nM",
            "NM",
            "MD",
            "jM",
            "jI",
            "--alignSJoverhangMin",
            "10",
            "--alignSJDBoverhangMin",
            "10",
            "--outFilterIntronMotifs",
            "None",
            "--quantMode",
            "GeneCounts",
        ],
        out / "STAR_run.log",
    )
    bamfile = Path(prefix + "Aligned.sortedByCoord.out.bam")
    run(["samtools", "quickcheck", "-v", bamfile], out / "quickcheck.log")
    run(["samtools", "index", "-@", min(threads, 4), bamfile], out / "index.log")
    db = load()
    wanted = set(db["motifs"])
    counts = collections.Counter()
    unknown = collections.Counter()
    qc = collections.Counter()
    with pysam.AlignmentFile(str(bamfile), "rb", threads=min(threads, 4)) as bam:
        for read in bam.fetch(until_eof=True):
            qc["alignment_records"] += 1
            if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_qcfail:
                continue
            if not read.has_tag("NH") or read.get_tag("NH") != 1:
                continue
            qc["primary_unique_reads"] += 1
            for a, b, left, right in edges_and_anchors(read):
                if min(left, right) < 10:
                    continue
                key = (read.reference_name, a, b)
                if key not in wanted:
                    continue
                motif, st = db["motifs"][key]
                if st == "?":
                    unknown[key] += 1
                    qc["query_coordinate_hits_strand_unresolved"] += 1
                else:
                    counts[(read.reference_name, st, a, b)] += 1
                    qc["query_coordinate_hits_strand_resolved"] += 1
    ji = {j: i for i, j in enumerate(db["junctions"])}
    with gzip.open(out / "junction_counts.tsv.gz", "wt") as f:
        f.write("junction_id\tunique_read_alignments\tstrand_unresolved_coordinate_reads\n")
        for key, num in sorted(ji.items(), key=lambda x: x[1]):
            ch, st, a, b = key
            v = counts[key]
            amb = unknown[(ch, a, b)]
            if v or amb:
                f.write(f"{num}\t{v}\t{amb}\n")
    dump(
        out / "COMPLETE.json",
        dict(
            qc,
            sample=sample["ont_sample"],
            status="COMPLETE",
            strand_policy="Genome splice motif; never BAM alignment flag for unstranded RNA-seq. Noncanonical unresolved-strand coordinate evidence retained separately.",
        ),
    )
    print("RNASEQ_COMPLETE", sample["ont_sample"], dict(qc), flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("stage", choices=["prepare", "index", "ont", "rna"])
    p.add_argument("--expected-chains", type=int)
    p.add_argument("--sample-index", type=int)
    p.add_argument("--threads", type=int, default=8)
    a = p.parse_args()
    if a.stage == "prepare":
        prepare(a.expected_chains)
    elif a.stage == "index":
        index(a.threads)
    elif a.stage == "ont":
        ont(a.sample_index, a.threads)
    else:
        rna(a.sample_index, a.threads)
