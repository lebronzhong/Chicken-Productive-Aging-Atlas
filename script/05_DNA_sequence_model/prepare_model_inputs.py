#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 05_DNA_sequence_model
# Purpose: Model training data preparation.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

import csv, json, pathlib, sys, shutil, hashlib, subprocess, concurrent.futures, gzip, collections
import pysam, pyBigWig

ROOT = pathlib.Path(__file__).resolve().parent
OLD = ROOT.parent
ONT = pathlib.Path("/path/to/chicken-aging/ont")
PY = sys.executable
OUT = ROOT / "training_inputs_v1"
STUDY = json.loads((ROOT / "study_samples.json").read_text(encoding="utf8"))
PREFIXES = STUDY["PREFIXES"]
NAMES = STUDY["NAMES"]
GROUPS = STUDY["GROUPS"]


def dump(name, obj):
    (ROOT / name).write_text(json.dumps(obj, indent=2, ensure_ascii=False) + "\n")


def tsv(path, rows):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def audit():
    old = list(
        csv.DictReader(open(OLD / "samples.matched.tsv", encoding="utf-8-sig"), delimiter="\t")
    )
    bysample = {s["ont_sample_id"]: s for s in old}
    byanimal = {s["animal_id"]: s for s in old}
    with gzip.open(OLD / "variants/cohort.original.vcf.gz", "rt") as f:
        for line in f:
            if line.startswith("#CHROM"):
                vcfids = set(line.rstrip().split("\t")[9:])
                break
    expected = {g + "-" + str(i) for g, (_, ids) in GROUPS.items() for i in ids}
    assert vcfids == expected, (vcfids, expected)
    ref = {
        x.split()[0]: int(x.split()[1])
        for x in open(OLD / "reference_assets/GRCg6a.training.fa.fai")
    }
    rows = []
    errors = []
    checks = []
    for prefix in PREFIXES:
        for group, (tag, ids) in GROUPS.items():
            for i in ids:
                sid = f"{prefix}{i}" if group == "80wk" else f"{prefix}{tag}-{i}"
                animal = f"{group}-{i}"
                bam = ONT / f"{group}_data" / f"{group}_ONTdata" / sid / "4_samtools/aln_sort.bam"
                row = dict(bysample.get(sid, byanimal[animal]))
                row.update(
                    ont_sample_id=sid,
                    tissue_id=NAMES[prefix],
                    bam_path=str(bam),
                    quant_sample_id=sid,
                )
                row["source_tissue_name"] = prefix
                row["count_file"] = str(
                    ONT
                    / "shell/12_20240919_TPM.count/count/5sample_over3_Tissue"
                    / f"{prefix}_merged_transcript_grouped_counts_filtered.tsv"
                )
                row["bam_index_path"] = str(bam) + ".bai"
                row["track_origin"] = "reused" if sid in bysample else "new"
                row["raw_track_dir"] = str(
                    (OLD if sid in bysample else ROOT) / "tracks_unstranded_v1" / sid
                )
                try:
                    assert bam.is_file() and bam.stat().st_size > 0
                    pysam.quickcheck(str(bam))
                    with pysam.AlignmentFile(str(bam), "rb") as b:
                        assert b.header.to_dict().get("HD", {}).get("SO") == "coordinate"
                        header = dict(zip(b.references, b.lengths))
                        assert all(
                            header.get(c) == n for c, n in ref.items()
                        ), "Reference sequence mismatch"
                        assert b.check_index()
                        next(b.fetch("1", 0, min(ref["1"], 1000000)), None)
                    if sid in bysample:
                        d = pathlib.Path(row["raw_track_dir"])
                        m = json.load(open(d / "metadata.json"))
                        assert (
                            m["sample"]["bam_path"] == str(bam)
                            and m["sample"]["animal_id"] == animal
                        )
                        assert (
                            m["independent_regions_passed"] > 0
                            and m["mode"] == "unstranded_raw_depth"
                        )
                        assert (d / "coverage.bw").stat().st_size == m["bigwig_bytes"]
                        with pyBigWig.open(str(d / "coverage.bw")) as bw:
                            assert bw.chroms() == header
                    row["include_for_training"] = "QC_pass"
                    checks.append(
                        {
                            "sample": sid,
                            "bam_bytes": bam.stat().st_size,
                            "index_older_than_bam": pathlib.Path(row["bam_index_path"])
                            .stat()
                            .st_mtime
                            < bam.stat().st_mtime,
                            "count_exists": pathlib.Path(row["count_file"]).is_file(),
                        }
                    )
                except Exception as e:
                    errors.append({"sample": sid, "error": repr(e)})
                    row["include_for_training"] = "QC_fail"
                    print("ERROR", sid, repr(e), flush=True)
                rows.append(row)
        print("AUDIT", prefix, len(rows), "errors", len(errors), flush=True)
    tsv(ROOT / "samples.matched.tsv", rows)
    tsv(
        ROOT / "tissue_sample_counts.tsv",
        [
            {
                "tissue": NAMES[p],
                "source_name": p,
                **{
                    g: sum(r["tissue_id"] == NAMES[p] and r["group"] == g for r in rows)
                    for g in GROUPS
                },
            }
            for p in PREFIXES
        ],
    )
    dump("audit_checks.json", checks)
    dump("audit_errors.json", errors)
    dump(
        "audit_summary.json",
        {
            "samples": len(rows),
            "tissues": len(PREFIXES),
            "animals": len(expected),
            "reused": sum(r["track_origin"] == "reused" for r in rows),
            "new": sum(r["track_origin"] == "new" for r in rows),
            "errors": errors,
            "identity_evidence": "Sample-ID mapping and agreement with VCF sample identifiers; no genetic-fingerprint comparison",
            "library": "ONT PCR-cDNA, unstranded coverage",
            "sequencing_batches": "20wk;80wk;50wk+100wk+M",
            "name_confirmation": "Tissue-label mapping: Vignia=vagina, Uterine=uterus/shell gland, Pectoral=pectoral muscle",
        },
    )
    assert len(rows) == 825 and not errors, "Audit failed: inspect audit_errors.json"
    OUT.mkdir(exist_ok=True)
    for name in ["chromosome_splits.tsv"]:
        shutil.copy2(OLD / "training_inputs_v1" / name, OUT / name)
    shutil.copy2(
        OLD / "downstream_20261014/inputs/data/tfrecords_v2/sequences.bed", OUT / "sequences.bed"
    )
    dump("tissue_order.json", [NAMES[p] for p in PREFIXES])
    assert (ROOT / "make_unstranded_v1.py").is_file() and (
        ROOT / "merge_training_tracks_v1.py"
    ).is_file()
    dump(
        "provenance.json",
        {
            "distributed_worker_sha256": sha(ROOT / "make_unstranded_v1.py"),
            "distributed_merge_sha256": sha(ROOT / "merge_training_tracks_v1.py"),
            "sequences_sha256": sha(OUT / "sequences.bed"),
            "split_sha256": sha(OUT / "chromosome_splits.tsv"),
            "target_transform": "sqrt(1+sum32(0.3*depth))-1; above384:383+sqrt(v-383); clip0..768; float16",
            "sequence_length": 393216,
            "pool_width": 32,
            "individual_holdout": False,
            "raw_tracks": "Primary mapped not QC-fail MAPQ>=20; duplicates retained; CIGAR N/D excluded",
        },
    )
    (ROOT / "AUDIT_PASS").write_text("825 samples checked\n")


def tracks():
    assert (ROOT / "AUDIT_PASS").exists()
    rows = list(csv.DictReader(open(ROOT / "samples.matched.tsv"), delimiter="\t"))
    todo = [i for i, r in enumerate(rows) if r["track_origin"] == "new"]
    (ROOT / "tracks_unstranded_v1").mkdir(exist_ok=True)
    # Representatives for all 15 added tissues finish before the remaining samples.
    first = []
    seen = set()
    for i in todo:
        if rows[i]["tissue_id"] not in seen:
            first.append(i)
            seen.add(rows[i]["tissue_id"])

    def run(i):
        r = rows[i]
        d = pathlib.Path(r["raw_track_dir"])
        if d.exists():
            m = json.load(open(d / "metadata.json"))
            assert (
                m["sample"]["bam_path"] == r["bam_path"]
                and m["independent_regions_passed"] > 0
                and (d / "coverage.bw").stat().st_size == m["bigwig_bytes"]
            )
            return
        with open(ROOT / "logs" / f"sample_{r['ont_sample_id']}.log", "w") as f:
            subprocess.run(
                [PY, "-u", str(ROOT / "make_unstranded_v1.py"), "--index", str(i)],
                stdout=f,
                stderr=subprocess.STDOUT,
                check=True,
            )
        print("COMPLETE", r["ont_sample_id"], flush=True)

    for batch in [first, [i for i in todo if i not in first]]:
        with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:
            for result in pool.map(run, batch):
                pass
    dump("new_tracks_summary.json", {"completed": len(todo), "representatives": len(first)})


def merge():
    assert (ROOT / "new_tracks_summary.json").exists()
    rows = list(csv.DictReader(open(ROOT / "samples.matched.tsv"), delimiter="\t"))
    chroms = {x.split()[0] for x in list(open(OUT / "chromosome_splits.tsv"))[1:]}
    inputs = {t: [] for t in json.load(open(ROOT / "tissue_order.json"))}
    for r in rows:
        d = pathlib.Path(r["raw_track_dir"])
        m = json.load(open(d / "metadata.json"))
        assert (
            m["independent_regions_passed"] > 0
            and (d / "coverage.bw").stat().st_size == m["bigwig_bytes"]
        )
        denom = sum(c["reads"] for c in m["chromosomes"] if c["chrom"] in chroms)
        assert denom > 0
        inputs[r["tissue_id"]].append(
            {
                "animal": r["animal_id"],
                "sample": r["ont_sample_id"],
                "bigwig": str(d / "coverage.bw"),
                "autosomal_primary_reads": denom,
                "scale": 1e7 / denom,
            }
        )
    for values in inputs.values():
        assert len(values) == 25 and len({v["animal"] for v in values}) == 25
    (OUT / "merge_inputs.json").write_text(json.dumps(inputs, indent=2))
    subprocess.run([PY, "-u", str(ROOT / "merge_training_tracks_v1.py")], check=True)


def bundle():
    summaries = json.load(open(OUT / "merge_summary.json"))
    assert len(summaries) == 33
    dest = ROOT / "transfer_bundle_33t25_v2"
    dest.mkdir(exist_ok=False)
    for folder in ["reference", "metadata", "tracks"]:
        (dest / folder).mkdir()

    def link(src, target):
        import os

        os.link(src, target)

    for name in [
        "GRCg6a.training.fa",
        "GRCg6a.training.fa.fai",
        "chrom.sizes",
        "assembly_gaps_N.bed",
        "transcript_exons.tsv",
    ]:
        link(OLD / "reference_assets" / name, dest / "reference" / name)
    for name in [
        "samples.matched.tsv",
        "tissue_sample_counts.tsv",
        "audit_summary.json",
        "provenance.json",
        "tissue_order.json",
    ]:
        shutil.copy2(ROOT / name, dest / "metadata" / name)
    for name in [
        "chromosome_splits.tsv",
        "sequences.bed",
        "merge_summary.json",
        "merge_inputs.json",
    ]:
        shutil.copy2(OUT / name, dest / "metadata" / name)
    targets = []
    for i, t in enumerate(json.load(open(ROOT / "tissue_order.json"))):
        d = OUT / "tissue_tracks" / t
        m = json.load(open(d / "metadata.json"))
        assert m["training_animals"] == 25 and m["checks"] > 0
        link(d / "coverage.bw", dest / "tracks" / f"{t}.bw")
        shutil.copy2(d / "metadata.json", dest / "metadata" / f"{t}.track.json")
        targets.append(
            {
                "index": i,
                "identifier": t,
                "file": f"tracks/{t}.bw",
                "strand_pair": i,
                "description": f"{t}: mean25 normalized unstranded ONT depth",
            }
        )
    tsv(dest / "targets.tsv", targets)
    (dest / "README.txt").write_text(
        "33 tissues;25 animals; reference DNA input; fixed chromosome holdouts. These are coverage inputs, NOT training-ready TFRecords. Apply the recorded target transform and validate tensor/channel order during conversion. No individual holdout. Original test chromosomes have been inspected in earlier model development. Hard-linked source files must remain immutable.\n"
    )
    files = sorted(p for p in dest.rglob("*") if p.is_file())
    with open(dest / "SHA256SUMS", "w") as f:
        for p in files:
            f.write(sha(p) + "  " + str(p.relative_to(dest)) + "\n")
    dump(
        "CPU_PREPARATION_COMPLETE.json",
        {
            "bundle": str(dest),
            "files": len(files),
            "bytes": sum(p.stat().st_size for p in files),
            "status": "INPUT_BUNDLE_READY_NOT_TFRECORDS",
        },
    )


if __name__ == "__main__":
    {"audit": audit, "tracks": tracks, "merge": merge, "bundle": bundle}[sys.argv[1]]()
