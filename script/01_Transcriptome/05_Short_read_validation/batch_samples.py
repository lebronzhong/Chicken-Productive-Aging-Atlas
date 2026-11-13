#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 01_Transcriptome
# Purpose: Batch samples.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

import argparse, subprocess, sys, concurrent.futures, json
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument("stage", choices=["ont", "rna"])
p.add_argument("--workers", type=int, required=True)
p.add_argument("--threads", type=int, required=True)
p.add_argument("--indices", default="0-20")
a = p.parse_args()
indices = (
    list(range(int(a.indices.split("-")[0]), int(a.indices.split("-")[1]) + 1))
    if "-" in a.indices
    else [int(x) for x in a.indices.split(",")]
)
partition = Path(__file__).resolve().parent.parent / "rna_work_partition.json"
if a.stage == "rna" and a.indices == "3-20" and partition.exists():
    indices = json.loads(partition.read_text())["primary_remaining_indices"]
print("SCHEDULED_INDICES", a.stage, indices, flush=True)


def one(i):
    print("START_SAMPLE", a.stage, i, flush=True)
    subprocess.run(
        [
            sys.executable,
            "-u",
            str(Path(__file__).with_name("validation.py")),
            a.stage,
            "--sample-index",
            str(i),
            "--threads",
            str(a.threads),
        ],
        check=True,
    )
    print("DONE_SAMPLE", a.stage, i, flush=True)


with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as ex:
    for f in concurrent.futures.as_completed([ex.submit(one, i) for i in indices]):
        f.result()
print("SCHEDULED_SAMPLES_COMPLETE", a.stage, indices, flush=True)
