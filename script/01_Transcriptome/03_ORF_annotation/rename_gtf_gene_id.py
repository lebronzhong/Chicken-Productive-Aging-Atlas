#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 01_Transcriptome
# Purpose: Rename gtf gene id.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

import re
import argparse
from collections import Counter


def parse_attributes(attr_str):
    """

    Parse the attributes in GTF column 9 into a dictionary.

    Example:

    gene_id "MSTRG.1"; transcript_id "TCONS_1"; ref_gene_id "ENSGALG...";

    """
    attrs = {}
    for item in re.finditer(r'(\S+)\s+"([^"]*)"', attr_str):
        key, value = item.group(1), item.group(2)
        attrs[key] = value
    return attrs


def rebuild_attributes(attr_str, updates):
    """

    Update selected attributes while preserving their original order and format.

    updates: dictionary mapping attribute names to replacement values.

    """
    for key, value in updates.items():
        pattern = rf'({re.escape(key)}\s+")([^"]*)(")'
        if re.search(pattern, attr_str):
            attr_str = re.sub(pattern, rf"\1{value}\3", attr_str, count=1)
        else:
            # Append attributes that were absent from the original record.
            attr_str = attr_str.strip()
            if not attr_str.endswith(";"):
                attr_str += ";"
            attr_str += f' {key} "{value}";'
    return attr_str


def first_pass_collect_transcript_map(gtf_file):
    """

    First pass:

    Collect transcript_id to ref_gene_id mappings.

    Also collect gene_id to ref_gene_id mappings from transcript records for optional gene-feature updates.

    """
    transcript_to_refgene = {}
    gene_to_refgene = {}
    stats = Counter()
    with open(gtf_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                stats["malformed_lines"] += 1
                continue
            feature = fields[2]
            attrs = parse_attributes(fields[8])
            if feature == "transcript":
                stats["transcript_lines"] += 1
                transcript_id = attrs.get("transcript_id")
                ref_gene_id = attrs.get("ref_gene_id")
                old_gene_id = attrs.get("gene_id")
                if transcript_id:
                    if ref_gene_id:
                        if (
                            transcript_id in transcript_to_refgene
                            and transcript_to_refgene[transcript_id] != ref_gene_id
                        ):
                            stats["conflicting_transcript_id"] += 1
                        transcript_to_refgene[transcript_id] = ref_gene_id
                        stats["transcripts_with_ref_gene_id"] += 1
                    else:
                        stats["transcripts_without_ref_gene_id"] += 1
                if old_gene_id and ref_gene_id:
                    if (
                        old_gene_id in gene_to_refgene
                        and gene_to_refgene[old_gene_id] != ref_gene_id
                    ):
                        stats["conflicting_gene_id_mapping"] += 1
                    gene_to_refgene[old_gene_id] = ref_gene_id
    return transcript_to_refgene, gene_to_refgene, stats


def second_pass_update_gtf(
    input_gtf,
    output_gtf,
    transcript_to_refgene,
    gene_to_refgene,
    update_gene_feature=False,
    missing_log=None,
):
    """

    Second pass:

    - Update gene_id for all features carrying a transcript_id.

    - Optionally update gene features.

    - Write logs and summary counts.

    """
    stats = Counter()
    missing_transcripts = set()
    with open(input_gtf, "r") as infile, open(output_gtf, "w") as outfile:
        for line in infile:
            if line.startswith("#"):
                outfile.write(line)
                continue
            raw_line = line.rstrip("\n")
            fields = raw_line.split("\t")
            if len(fields) != 9:
                stats["malformed_lines"] += 1
                continue
            feature = fields[2]
            attrs_str = fields[8]
            attrs = parse_attributes(attrs_str)
            updated = False
            # 1. Optionally update gene features.
            if feature == "gene" and update_gene_feature:
                old_gene_id = attrs.get("gene_id")
                if old_gene_id and old_gene_id in gene_to_refgene:
                    new_gene_id = gene_to_refgene[old_gene_id]
                    if old_gene_id != new_gene_id:
                        fields[8] = rebuild_attributes(attrs_str, {"gene_id": new_gene_id})
                        stats["updated_gene_features"] += 1
                        updated = True
            # 2. Update every feature carrying a transcript_id.
            transcript_id = attrs.get("transcript_id")
            if transcript_id:
                stats["features_with_transcript_id"] += 1
                if transcript_id in transcript_to_refgene:
                    new_gene_id = transcript_to_refgene[transcript_id]
                    old_gene_id = attrs.get("gene_id")
                    if old_gene_id:
                        if old_gene_id != new_gene_id:
                            fields[8] = rebuild_attributes(fields[8], {"gene_id": new_gene_id})
                            stats["updated_feature_gene_id"] += 1
                            updated = True
                        else:
                            stats["already_correct_feature_gene_id"] += 1
                    else:
                        fields[8] = rebuild_attributes(fields[8], {"gene_id": new_gene_id})
                        stats["added_missing_gene_id"] += 1
                        updated = True
                else:
                    stats["features_with_unmapped_transcript_id"] += 1
                    missing_transcripts.add(transcript_id)
            if not updated:
                stats["unchanged_lines"] += 1
            outfile.write("\t".join(fields) + "\n")
    # Write the list of unmapped transcripts.
    if missing_log:
        with open(missing_log, "w") as logf:
            for tid in sorted(missing_transcripts):
                logf.write(tid + "\n")
    return stats


def main():
    parser = argparse.ArgumentParser(
        description="Rename GTF gene_id attributes using ref_gene_id values from transcript features."
    )
    parser.add_argument("-i", "--input", required=True, help="Input GTF file.")
    parser.add_argument("-o", "--output", required=True, help="Output GTF file.")
    parser.add_argument(
        "--update-gene-feature",
        action="store_true",
        help="Also update gene-feature gene_id attributes (disabled by default).",
    )
    parser.add_argument("--missing-log", default=None, help="Write unmapped transcript IDs.")
    args = parser.parse_args()
    print("[1/2] Collecting transcript_id -> ref_gene_id mapping ...")
    transcript_to_refgene, gene_to_refgene, stats1 = first_pass_collect_transcript_map(args.input)
    print(f"  Total transcript lines               : {stats1['transcript_lines']}")
    print(f"  Transcripts with ref_gene_id         : {stats1['transcripts_with_ref_gene_id']}")
    print(f"  Transcripts without ref_gene_id      : {stats1['transcripts_without_ref_gene_id']}")
    print(f"  Conflicting transcript_id mappings   : {stats1['conflicting_transcript_id']}")
    print(f"  Conflicting gene_id mappings         : {stats1['conflicting_gene_id_mapping']}")
    print(f"  Malformed lines in pass 1            : {stats1['malformed_lines']}")
    print("[2/2] Updating GTF ...")
    stats2 = second_pass_update_gtf(
        args.input,
        args.output,
        transcript_to_refgene,
        gene_to_refgene,
        update_gene_feature=args.update_gene_feature,
        missing_log=args.missing_log,
    )
    print("\nDone.")
    print("Summary:")
    print(f"  Features with transcript_id          : {stats2['features_with_transcript_id']}")
    print(f"  Updated feature gene_id              : {stats2['updated_feature_gene_id']}")
    print(f"  Added missing gene_id                : {stats2['added_missing_gene_id']}")
    print(f"  Already correct feature gene_id      : {stats2['already_correct_feature_gene_id']}")
    print(
        f"  Features with unmapped transcript_id : {stats2['features_with_unmapped_transcript_id']}"
    )
    print(f"  Updated gene features                : {stats2['updated_gene_features']}")
    print(f"  Unchanged lines                      : {stats2['unchanged_lines']}")
    print(f"  Malformed lines in pass 2            : {stats2['malformed_lines']}")
    if args.missing_log:
        print(f"  Missing transcript log               : {args.missing_log}")


if __name__ == "__main__":
    main()
