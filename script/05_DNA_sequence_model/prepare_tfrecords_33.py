#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 05_DNA_sequence_model
# Purpose: Model training data preparation.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

"""Chicken-specific orchestration of unmodified official Borzoi TFRecord tools."""

import concurrent.futures as cf
import csv, json, os, random, shutil, subprocess, sys, time
from collections import Counter
import argparse

from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--expected-split-counts", type=Path)
args = parser.parse_args()
import h5py
import numpy as np
import pyBigWig
import pysam
from baskerville import data

ROOT = Path("/path/to/chicken-aging/model")
PROJECT = ROOT / "experiment33"
CORE = PROJECT / "transfer_bundle_33t25_v2"
WORK = PROJECT / "cpu_prepare_33"
SCRIPTS = ROOT / "software/baskerville-tf/src/baskerville/scripts"
PY = str(ROOT / "software/borzoi_py310/bin/python")
FA = CORE / "reference/GRCg6a.training.fa"
LENGTH, WIDTH, SCALE, CLIP_SOFT, CLIP = 393216, 32, 0.3, 384.0, 768.0
PROCESSES = 4


def write_json(path, obj):
    path.write_text(json.dumps(obj, indent=2) + "\n")


def command(args, log):
    with log.open("w") as handle:
        subprocess.run(
            [PY] + list(map(str, args)), stdout=handle, stderr=subprocess.STDOUT, check=True
        )


def parallel(jobs):
    with cf.ThreadPoolExecutor(max_workers=PROCESSES) as pool:
        futures = [pool.submit(command, args, log) for args, log in jobs]
        for future in cf.as_completed(futures):
            future.result()


def plan():
    assert (WORK / "TRANSFER_VERIFIED.txt").is_file(), "Transfer verification must pass first"
    with (CORE / "targets.tsv").open() as handle:
        tracks = list(csv.DictReader(handle, delimiter="\t"))
    assert len(tracks) == 33
    for t in tracks:
        t["tissue"] = t["identifier"]
    with (CORE / "metadata/chromosome_splits.tsv").open() as handle:
        chroms = list(csv.DictReader(handle, delimiter="\t"))
    splitmap = {row["chrom"]: row["split"] for row in chroms}
    assert set(splitmap.values()) == {"train", "valid", "test"}
    segments = {row["chrom"]: [(0, int(row["length"]))] for row in chroms}
    with pysam.FastaFile(str(FA)) as fa:
        for c, spans in segments.items():
            assert fa.get_reference_length(c) == spans[0][1]
    for i, t in enumerate(tracks):
        assert int(t["index"]) == int(t["strand_pair"]) == i
        t["absolute_file"] = str(CORE / t["file"])
        with pyBigWig.open(t["absolute_file"]) as bw:
            for c, spans in segments.items():
                assert bw.chroms(c) == spans[0][1]
    sequences = []
    # Reuse the exact old window order, not just its split recipe.
    from collections import namedtuple

    Window = namedtuple("Window", "chr start end label")
    for line in (CORE / "metadata/sequences.bed").read_text().splitlines():
        c, a, b, label = line.split()
        sequences.append(Window(c, int(a), int(b), label))
    assert set(s.label for s in sequences) == {"train", "valid", "test"}
    if args.expected_split_counts:
        assert dict(Counter(s.label for s in sequences)) == json.loads(
            args.expected_split_counts.read_text(encoding="utf8")
        )
    assert [s.label for s in sequences] == sorted(
        [s.label for s in sequences], key=["train", "valid", "test"].index
    )
    assert all(s.end - s.start == LENGTH and splitmap[s.chr] == s.label for s in sequences)
    data.write_seqs_bed(str(WORK / "all_sequences.bed"), sequences, True)
    counts = Counter(s.label for s in sequences)
    report = {
        "counts": dict(counts),
        "windows_per_chromosome": dict(Counter(s.chr for s in sequences)),
        "chromosomes_without_windows": sorted(set(splitmap) - {s.chr for s in sequences}),
        "sequence_length": LENGTH,
        "bin_width": WIDTH,
        "crop_bp": 0,
        "seed": 20261012,
        "train_stride": 131072,
        "valid_test_stride": LENGTH,
        "split": "Fixed whole chromosomes from transferred manifest; no expression-based splitting.",
        "gaps": "No window crosses an assembly N gap. This is not a mappability blacklist.",
        "normalization": "Input is mean coverage of 25 animals normalized per 10M autosomal reads.",
        "transform": "Official latest sum_sqrt: sqrt(1+sum(0.3*depth per 32bp))-1; soft_clip=384; hard_clip=768; clip_pct=1 (no global quantile fitting).",
        "reference": "Full FASTA retained; only the fixed autosomal window list is passed to official tools.",
        "official_commits": {
            "borzoi": "f10f8c799639dd2cfc0c1eb09499bfc47ab9cf8a",
            "baskerville-tf": "b06950be1be704ab09ab23ffb6f5102d9e2a63a0",
        },
    }
    write_json(WORK / "conversion_design.json", report)
    print("WINDOW_PLAN", json.dumps(report), flush=True)
    return tracks, sequences, splitmap


def convert(out, tracks, sequences):
    assert not out.exists(), "Refusing to overwrite " + str(out)
    out.mkdir(parents=True)
    (out / "seqs_cov").mkdir()
    (out / "tfrecords").mkdir()
    (out / "logs").mkdir()
    bed = out / "sequences.bed"
    data.write_seqs_bed(str(bed), sequences, True)
    with (out / "targets.txt").open("w") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "",
                "identifier",
                "file",
                "clip",
                "clip_soft",
                "scale",
                "sum_stat",
                "strand_pair",
                "description",
            ]
        )
        for i, t in enumerate(tracks):
            writer.writerow(
                [
                    i,
                    t["tissue"].replace(" ", "_"),
                    t["absolute_file"],
                    CLIP,
                    CLIP_SOFT,
                    SCALE,
                    "sum_sqrt",
                    i,
                    "ONT_unstranded:" + t["tissue"],
                ]
            )
    jobs = []
    for i, t in enumerate(tracks):
        jobs.append(
            (
                [
                    WORK / "hound_data_read_float32.py",
                    "-w",
                    WIDTH,
                    "-u",
                    "sum_sqrt",
                    "-s",
                    SCALE,
                    "--clip_soft",
                    CLIP_SOFT,
                    "-c",
                    CLIP,
                    "--clip_pct",
                    1,
                    t["absolute_file"],
                    bed,
                    out / "seqs_cov" / ("%d.h5" % i),
                ],
                out / "logs" / ("read_%d.log" % i),
            )
        )
    parallel(jobs)
    for i in range(33):
        with h5py.File(out / "seqs_cov" / ("%d.h5" % i), "r") as h:
            assert h["targets"].shape == (len(sequences), LENGTH // WIDTH)
    jobs = []
    shardmap = []
    offset = 0
    for split in ("train", "valid", "test"):
        n = sum(s.label == split for s in sequences)
        for shard, start in enumerate(range(offset, offset + n, 64)):
            end = min(start + 64, offset + n)
            filename = "%s-%d.tfr" % (split, shard)
            jobs.append(
                (
                    [
                        SCRIPTS / "hound_data_write.py",
                        "-s",
                        start,
                        "-e",
                        end,
                        FA,
                        bed,
                        out / "seqs_cov",
                        out / "tfrecords" / filename,
                    ],
                    out / "logs" / ("write_" + filename + ".log"),
                )
            )
            shardmap.append({"file": filename, "start": start, "end": end, "split": split})
        offset += n
    parallel(jobs)
    stats = {
        "num_targets": 33,
        "seq_length": LENGTH,
        "seq_1hot": True,
        "pool_width": WIDTH,
        "crop_bp": 0,
        "target_length": LENGTH // WIDTH,
    }
    stats.update(
        {s + "_seqs": sum(x.label == s for x in sequences) for s in ("train", "valid", "test")}
    )
    write_json(out / "statistics.json", stats)
    write_json(out / "shards.json", shardmap)
    return shardmap


def expected_track(bw, seq):
    raw = bw.values(seq.chr, seq.start, seq.end, numpy=True).astype("float32")
    assert not np.isinf(raw).any() and (raw[np.isfinite(raw)] >= 0).all()
    raw = np.nan_to_num(raw, nan=0.0)
    values = -1 + np.sqrt(1 + (SCALE * raw).reshape(-1, WIDTH).sum(axis=1, dtype="float32"))
    mask = values > CLIP_SOFT
    values[mask] = CLIP_SOFT - 1 + np.sqrt(values[mask] - CLIP_SOFT + 1)
    return np.clip(values, -CLIP, CLIP).astype("float16")


def validate(out, tracks, sequences, shards, splitmap):
    import tensorflow as tf

    options = tf.data.Options()
    options.threading.private_threadpool_size = 1
    counts = Counter()
    compared = 0
    nonzero = np.zeros(33, dtype="int64")
    total_bins = 0
    lut = np.full(256, 4, dtype="uint8")
    for i, b in enumerate(b"ACGT"):
        lut[b] = i
    fa = pysam.FastaFile(str(FA))
    bws = [pyBigWig.open(t["absolute_file"]) for t in tracks]
    try:
        for shard in shards:
            count = 0
            ds = tf.data.TFRecordDataset(
                str(out / "tfrecords" / shard["file"]),
                compression_type="ZLIB",
                num_parallel_reads=1,
            ).with_options(options)
            for record in ds:
                ex = tf.train.Example.FromString(bytes(record.numpy()))
                x = np.frombuffer(
                    ex.features.feature["sequence"].bytes_list.value[0], dtype="uint8"
                )
                y = np.frombuffer(
                    ex.features.feature["target"].bytes_list.value[0], dtype="float16"
                ).reshape(LENGTH // WIDTH, 33)
                assert x.shape == (LENGTH,) and x.max() <= 4
                assert np.isfinite(y).all() and (y >= 0).all()
                seq = sequences[shard["start"] + count]
                assert seq.label == shard["split"] == splitmap[seq.chr]
                nonzero += (y > 0).sum(axis=0)
                total_bins += y.shape[0]
                if count == 0:
                    dna = fa.fetch(seq.chr, seq.start, seq.end).upper().encode("ascii")
                    assert np.array_equal(
                        x, lut[np.frombuffer(dna, dtype="uint8")]
                    ), "DNA/reference mismatch"
                    # Independent reference+bigWig reconstruction for all 33 channels, first record of every shard.
                    for i, bw in enumerate(bws):
                        expected = expected_track(bw, seq)
                        assert np.allclose(
                            y[:, i], expected, rtol=0.002, atol=0.002
                        ), "Coverage/target mismatch"
                        compared += 1
                count += 1
            assert count == shard["end"] - shard["start"], "Missing or extra records"
            counts[shard["split"]] += count
        assert all(
            counts[s] == sum(x.label == s for x in sequences) for s in ("train", "valid", "test")
        )
        assert np.all(nonzero > 0), "An entire tissue target is zero"
        report = {
            "status": "TFRECORD_CONTENT_QC_PASSED",
            "records": dict(counts),
            "shards": len(shards),
            "independent_bigwig_checks": compared,
            "nonzero_bin_fraction_by_tissue": dict(
                zip([t["tissue"] for t in tracks], (nonzero / total_bins).tolist())
            ),
            "checks": "Every record: decode, shape, finite/nonnegative targets, labels and counts. First record of each shard: reference DNA and all 33 bigWig targets.",
            "remaining": "Model input/output compatibility and real GPU training smoke test still required.",
        }
        write_json(out / "QC.json", report)
        print("QC", str(out), json.dumps(report), flush=True)
    finally:
        fa.close()
        for bw in bws:
            bw.close()


def main():
    tracks, seqs, splitmap = plan()
    pilot = [
        s for label in ("train", "valid", "test") for s in [x for x in seqs if x.label == label][:2]
    ]
    for name, selected in [("pilot_tfrecords_v2", pilot), ("tfrecords_v2", seqs)]:
        out = PROJECT / "data" / name
        start = time.monotonic()
        print("START", name, len(selected), flush=True)
        shards = convert(out, tracks, selected)
        validate(out, tracks, selected, shards, splitmap)
        print("DONE", name, "seconds", round(time.monotonic() - start), flush=True)
    write_json(
        WORK / "CPU_PREPARATION_READY_v2.json",
        {
            "status": "TFRECORDS_VALIDATED_GPU_TEST_PENDING",
            "data_dir": str(PROJECT / "data/tfrecords_v2"),
        },
    )


if __name__ == "__main__":
    main()
