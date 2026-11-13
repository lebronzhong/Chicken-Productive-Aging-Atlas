#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 06_ASE_ASTS
# Recovered study preparation script; inputs are supplied as arguments.
import sys
import os
import re
from collections import defaultdict
import pysam

# Read paths.
import argparse

parser = argparse.ArgumentParser(
    description="Prepare transcript FASTA, gene exon BED and transcript-to-gene map."
)
parser.add_argument("gtf")
parser.add_argument("reference_fasta")
parser.add_argument("output_dir")
args = parser.parse_args()
gtf, ref, outdir = args.gtf, args.reference_fasta, args.output_dir
os.makedirs(outdir, exist_ok=True)
gene_re = re.compile(r'gene_id "([^"]+)"')
tx_re = re.compile(r'transcript_id "([^"]+)"')
tx_exons = defaultdict(list)
tx_gene = {}
tx_chr = {}
tx_strand = {}
gene_exons = defaultdict(list)
# Collect transcript and gene exon coordinates.
with open(gtf) as f:
    for line in f:
        if line.startswith("#"):
            continue
        x = line.rstrip("\n").split("\t")
        if len(x) < 9:
            continue
        if x[2] != "exon":
            continue
        gm = gene_re.search(x[8])
        tm = tx_re.search(x[8])
        if gm is None or tm is None:
            continue
        gene = gm.group(1)
        tx = tm.group(1)
        chrom = x[0]
        start = int(x[3]) - 1
        end = int(x[4])
        strand = x[6]
        tx_exons[tx].append((start, end))
        tx_gene[tx] = gene
        tx_chr[tx] = chrom
        tx_strand[tx] = strand
        gene_exons[(chrom, gene)].append((start, end))
print("Transcripts:", len(tx_exons))
print("Genes:", len(set(tx_gene.values())))
# Write transcript-to-gene mapping.
txgene_out = os.path.join(outdir, "transcript_gene.tsv")
with open(txgene_out, "w") as o:
    o.write("transcript\tgene\n")
    for tx in sorted(tx_gene):
        o.write(f"{tx}\t{tx_gene[tx]}\n")
# Merge gene exon intervals and write BED.
bed_rows = []
for (chrom, gene), intervals in gene_exons.items():
    intervals = sorted(intervals)
    merged = []
    for start, end in intervals:
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    for start, end in merged:
        bed_rows.append((chrom, start, end, gene))
bed_rows.sort(key=lambda x: (x[0], x[1], x[2], x[3]))
bed_out = os.path.join(outdir, "gene_exons.bed")
with open(bed_out, "w") as o:
    for row in bed_rows:
        o.write("\t".join(map(str, row)) + "\n")
# Extract stranded transcript sequences.
fa = pysam.FastaFile(ref)
table = str.maketrans("ACGTRYMKWSBDHVNacgtrymkwsbdhvn", "TGCAYRKMWSVHDBNtgcayrkmwsvhdbn")


def revcomp(seq):
    return seq.translate(table)[::-1]


transcript_fa = os.path.join(outdir, "transcriptome.fa")
with open(transcript_fa, "w") as o:
    for tx in sorted(tx_exons):
        chrom = tx_chr[tx]
        intervals = sorted(tx_exons[tx], key=lambda x: x[0])
        seq = "".join(fa.fetch(chrom, start, end) for start, end in intervals)
        if tx_strand[tx] == "-":
            seq = revcomp(seq)
        o.write(f">{tx}\n")
        for i in range(0, len(seq), 80):
            o.write(seq[i : i + 80] + "\n")
fa.close()
print("gene_exons.bed:", bed_out)
print("transcript_gene.tsv:", txgene_out)
print("transcriptome.fa:", transcript_fa)
