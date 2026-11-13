#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 06_ASE_ASTS
# Recovered study preparation script; inputs are supplied as arguments.
import sys
import re
from pathlib import Path

# Read inputs and matching configuration.
import argparse
import json

parser = argparse.ArgumentParser(description="Match ONT sample directories to WGS sample IDs.")
parser.add_argument("ont_root", type=Path)
parser.add_argument("wgs_sample_file", type=Path)
parser.add_argument("manifest_all", type=Path)
parser.add_argument("manifest_matched", type=Path)
parser.add_argument("--config", type=Path, default=Path(__file__).with_name("lorals_datasets.json"))
args = parser.parse_args()
ont_root, wgs_sample_file = args.ont_root, args.wgs_sample_file
manifest_all, manifest_matched = args.manifest_all, args.manifest_matched
config = json.loads(args.config.read_text(encoding="utf-8"))
for target in [manifest_all, manifest_matched]:
    target.parent.mkdir(parents=True, exist_ok=True)
wgs_samples = [x.strip() for x in wgs_sample_file.read_text().splitlines() if x.strip()]
wgs_set = set(wgs_samples)
# Configure dataset paths.
dataset_info = {k: dict(v, root=ont_root / v["root"]) for k, v in config["datasets"].items()}
# Configure tissue aliases.
tissue_alias = config["tissue_aliases"]


# Match ONT sample folder names.
def match_sample(folder_name, age, mode):
    for tissue, aliases in tissue_alias.items():
        for alias in aliases:
            if mode == "age_dash":
                pattern = "^" + re.escape(alias) + re.escape(age) + r"-(\d+)$"
            elif mode == "number_only":
                pattern = "^" + re.escape(alias) + r"(\d+)$"
            elif mode == "M_dash":
                pattern = "^" + re.escape(alias) + r"M-(\d+)$"
            else:
                continue
            m = re.match(pattern, folder_name)
            if m:
                replicate = m.group(1)
                return (tissue, alias, replicate)
    return (None, None, None)


# Match WGS donor identifiers.
def find_wgs_sample(age, replicate):
    if age == "M":
        candidates = [f"M-{replicate}", f"M_{replicate}", f"M{replicate}"]
    else:
        candidates = [
            f"{age}wk-{replicate}",
            f"{age}wk_{replicate}",
            f"{age}wk{replicate}",
            f"{age}-{replicate}",
            f"{age}_{replicate}",
        ]
    for candidate in candidates:
        if candidate in wgs_set:
            return candidate
    return "NA"


# Collect available ONT libraries.
dataset_order = {name: i for i, name in enumerate(dataset_info, 1)}
rows = []
for dataset, info in dataset_info.items():
    age = info["age"]
    data_root = info["root"]
    mode = info["mode"]
    print("================================")
    print("Dataset:", dataset)
    print("Age:", age)
    print("Mode:", mode)
    print("Root:", data_root)
    if not data_root.exists():
        print("WARNING: directory does not exist:", data_root)
        continue
    n_all_dirs = 0
    n_name_matched = 0
    n_bam_found = 0
    n_wgs_matched = 0
    for sample_dir in sorted(data_root.iterdir()):
        if not sample_dir.is_dir():
            continue
        n_all_dirs += 1
        sample = sample_dir.name
        tissue, original_tissue, replicate = match_sample(sample, age, mode)
        if tissue is None:
            continue
        n_name_matched += 1
        bam = sample_dir / "4_samtools" / "aln_sort.bam"
        if not bam.exists():
            print("WARNING BAM missing:", bam)
            continue
        n_bam_found += 1
        vcf_sample = find_wgs_sample(age, replicate)
        if vcf_sample != "NA":
            n_wgs_matched += 1
        rows.append(
            [dataset, age, tissue, original_tissue, sample, replicate, vcf_sample, str(bam)]
        )
    print("All first-level directories:", n_all_dirs)
    print("Folders matching target tissues:", n_name_matched)
    print("BAM files found:", n_bam_found)
    print("Samples matched to WGS:", n_wgs_matched)
# Sort records by tissue, dataset and replicate.
rows = sorted(rows, key=lambda x: (x[2], dataset_order[x[0]], int(x[5])))
# Write the complete manifest.
with manifest_all.open("w") as o:
    o.write(
        "task_id\t"
        "dataset\t"
        "age\t"
        "tissue\t"
        "original_tissue\t"
        "ont_sample\t"
        "replicate\t"
        "vcf_sample\t"
        "bam\n"
    )
    for i, row in enumerate(rows, 1):
        o.write(str(i) + "\t" + "\t".join(row) + "\n")
# Write the WGS-matched manifest.
matched_rows = [row for row in rows if row[6] != "NA"]
with manifest_matched.open("w") as o:
    o.write(
        "array_id\t"
        "dataset\t"
        "age\t"
        "tissue\t"
        "original_tissue\t"
        "ont_sample\t"
        "replicate\t"
        "vcf_sample\t"
        "bam\n"
    )
    for i, row in enumerate(matched_rows, 1):
        o.write(str(i) + "\t" + "\t".join(row) + "\n")
# Report sample counts.
print("================================")
print("Total target BAM samples:", len(rows))
print("Matched WGS samples:", len(matched_rows))
print("Unmatched WGS samples:", len(rows) - len(matched_rows))
# Report counts by dataset.
print("===== Samples by dataset =====")
for dataset in dataset_info:
    n_total = sum(row[0] == dataset for row in rows)
    n_matched = sum(row[0] == dataset and row[6] != "NA" for row in rows)
    print(dataset, "total=", n_total, "matched=", n_matched)
# Report counts by tissue.
print("===== Samples by tissue =====")
for tissue in tissue_alias:
    n_total = sum(row[2] == tissue for row in rows)
    n_matched = sum(row[2] == tissue and row[6] != "NA" for row in rows)
    print(tissue, "total=", n_total, "matched=", n_matched)
# 13
print("===== 80wk accepted samples =====")
for row in rows:
    if row[0] == "80wk_data":
        print(row[4], "->", row[6], "->", row[7])
# 14
print("===== Unmatched samples =====")
for row in rows:
    if row[6] == "NA":
        print("UNMATCHED:", row[0], row[1], row[2], row[4], row[7])
