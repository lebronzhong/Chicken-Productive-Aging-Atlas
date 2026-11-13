#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 01_Transcriptome
# Purpose: Rename.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

import re


def update_gtf(input_file, output_file):
    transcript_dict = {}
    # Step 1: Read the input file and update gene_id for transcript and exon lines
    with open(input_file, "r") as infile:
        lines = infile.readlines()
    # First pass: Update gene_id in transcript and exon lines
    with open(output_file, "w") as outfile:
        for line in lines:
            if line.startswith("#"):
                outfile.write(line)
                continue
            fields = line.strip().split("\t")
            if len(fields) != 9:
                continue  # Skip lines that don't have 9 columns
            feature = fields[2]
            attributes = fields[8]
            if feature == "transcript":
                ref_gene_id_match = re.search(r'ref_gene_id "([^"]+)"', attributes)
                if ref_gene_id_match:
                    ref_gene_id = ref_gene_id_match.group(1)
                    # Update gene_id
                    updated_attributes = re.sub(
                        r'gene_id "[^"]+"', f'gene_id "{ref_gene_id}"', attributes
                    )
                    # Store the updated gene_id for this transcript_id
                    transcript_id_match = re.search(r'transcript_id "([^"]+)"', attributes)
                    if transcript_id_match:
                        old_transcript_id = transcript_id_match.group(1)
                        transcript_dict[old_transcript_id] = ref_gene_id
                    fields[8] = updated_attributes
            elif feature == "exon":
                transcript_id_match = re.search(r'transcript_id "([^"]+)"', attributes)
                if transcript_id_match:
                    old_transcript_id = transcript_id_match.group(1)
                    ref_gene_id = transcript_dict.get(old_transcript_id)
                    if ref_gene_id:
                        # Update gene_id
                        updated_attributes = re.sub(
                            r'gene_id "[^"]+"', f'gene_id "{ref_gene_id}"', attributes
                        )
                        fields[8] = updated_attributes
            outfile.write("\t".join(fields) + "\n")


# Use the function to update your GTF file
input_file = "9.6_finalgtf_.annotated.gtf"
output_file = "9.9_finalgtf_annotated_rename.gtf"
update_gtf(input_file, output_file)
