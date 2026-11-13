#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 07_Molecular_QTL
# Purpose: Prepare smr input.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

from pathlib import Path
import gzip
import csv
import time
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pyarrow.compute as pc

import argparse
import glob

parser = argparse.ArgumentParser(
    description="Prepare iso-uQTL nominal and phenotype inputs for SMR."
)
parser.add_argument(
    "--manifest",
    required=True,
    help="TSV: tissue, fdr_file, nominal_parquet_glob, phenotype_bed_glob",
)
parser.add_argument("--output-dir", required=True)
args = parser.parse_args()
manifest_path = Path(args.manifest).resolve()
with manifest_path.open(encoding="utf-8") as handle:
    manifest_rows = list(csv.DictReader(handle, delimiter="\t"))
required = {"tissue", "fdr_file", "nominal_parquet_glob", "phenotype_bed_glob"}
if not manifest_rows or not required.issubset(manifest_rows[0]):
    parser.error("Manifest is empty or missing required columns.")
if len({r["tissue"] for r in manifest_rows}) != len(manifest_rows):
    parser.error("Duplicate tissue in manifest.")
if any(
    not r["tissue"] or "/" in r["tissue"] or "\\" in r["tissue"] or r["tissue"] in {".", ".."}
    for r in manifest_rows
):
    parser.error("Tissue labels must be nonempty directory names.")
settings = {r["tissue"]: r for r in manifest_rows}
TISSUES = list(settings)
OUT = Path(args.output_dir)
OUT.mkdir(parents=True, exist_ok=True)
BATCH_SIZE = 20000


def input_path(tissue, column):
    value = settings[tissue][column]
    if not value:
        raise ValueError(f"{tissue}: empty {column}")
    path = Path(value)
    return str(path if path.is_absolute() else manifest_path.parent / path)


def find_fdr(tissue):
    path = Path(input_path(tissue, "fdr_file"))
    if not path.is_file():
        raise FileNotFoundError(path)
    return path


def matching_files(tissue, column):
    files = [Path(p) for p in sorted(glob.glob(input_path(tissue, column)))]
    if not files or not all(p.is_file() for p in files):
        raise FileNotFoundError(f"{tissue}: no valid files matched {column}")
    return files


def find_parquets(tissue):
    return matching_files(tissue, "nominal_parquet_glob")


def find_psi_beds(tissue):
    return matching_files(tissue, "phenotype_bed_glob")


def load_significant_transcripts(fdr_file):
    df = pd.read_csv(fdr_file, sep="\t", low_memory=False)
    if "phenotype_id" not in df.columns:
        raise RuntimeError(f"phenotype_id absent from:\n" f"{fdr_file}")
    # --------------------------------------------
    # Prefer q-value criterion
    # --------------------------------------------
    if "qval" in df.columns:
        q = pd.to_numeric(df["qval"], errors="coerce")
        df = df[q <= 0.05]
    # --------------------------------------------
    # Fallback
    # --------------------------------------------
    elif "is_eGene" in df.columns:
        sig = df["is_eGene"].astype(str).str.strip().str.lower().isin(["true", "1", "yes"])
        df = df[sig]
    else:
        raise RuntimeError(f"{fdr_file} has neither qval nor is_eGene.")
    tx = sorted(set(df["phenotype_id"].dropna().astype(str)))
    return tx


# ============================================================
# Detect parquet schema
# ============================================================
def detect_columns(parquet_file):
    schema = pq.ParquetFile(parquet_file).schema_arrow
    names = set(schema.names)
    required = [
        "phenotype_id",
        "variant_id",
        "af",
        "pval_nominal",
        "slope",
    ]
    missing = [x for x in required if x not in names]
    if missing:
        raise RuntimeError(
            f"Missing parquet columns in:\n"
            f"{parquet_file}\n"
            f"Missing: {missing}\n"
            f"Available: {schema.names}"
        )
    # --------------------------------------------
    # tensorQTL distance field
    # --------------------------------------------
    if "start_distance" in names:
        distance = "start_distance"
    elif "tss_distance" in names:
        distance = "tss_distance"
    elif "distance" in names:
        distance = "distance"
    else:
        raise RuntimeError(
            f"No distance column in:\n" f"{parquet_file}\n" f"Available: {schema.names}"
        )
    # --------------------------------------------
    # slope SE
    # --------------------------------------------
    slope_se = "slope_se" if "slope_se" in names else None
    columns = [
        "phenotype_id",
        "variant_id",
        distance,
        "af",
        "pval_nominal",
        "slope",
    ]
    if slope_se:
        columns.append(slope_se)
    return (columns, distance, slope_se)


# ============================================================
# Extract nominal QTL associations
#
# LOW-MEMORY:
# parquet is processed batch by batch.
# ============================================================
def extract_nominal(parquet_files, sig_tx, outfile):
    sig_arrow = pa.array(sig_tx, type=pa.string())
    total_scanned = 0
    total_written = 0
    columns = None
    distance_col = None
    slope_se_col = None
    with gzip.open(outfile, "wt", newline="") as fout:
        writer = csv.writer(fout, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "phenotype_id",
                "variant_id",
                "distance",
                "af",
                "pval_nominal",
                "slope",
                "slope_se",
            ]
        )
        # ----------------------------------------
        # Chromosome by chromosome
        # ----------------------------------------
        for idx, f in enumerate(parquet_files, start=1):
            # Determine schema once
            if columns is None:
                (columns, distance_col, slope_se_col) = detect_columns(f)
                print("Parquet columns:", columns, flush=True)
            pf = pq.ParquetFile(f)
            file_written = 0
            # ------------------------------------
            # Streaming batches
            # ------------------------------------
            for batch in pf.iter_batches(batch_size=BATCH_SIZE, columns=columns, use_threads=False):
                total_scanned += batch.num_rows
                pheno_index = batch.schema.get_field_index("phenotype_id")
                pheno = batch.column(pheno_index)
                mask = pc.is_in(pheno, value_set=sig_arrow)
                filtered = batch.filter(mask)
                if filtered.num_rows == 0:
                    continue
                d = filtered.to_pydict()
                n = filtered.num_rows
                for i in range(n):
                    writer.writerow(
                        [
                            d["phenotype_id"][i],
                            d["variant_id"][i],
                            d[distance_col][i],
                            d["af"][i],
                            d["pval_nominal"][i],
                            d["slope"][i],
                            (d[slope_se_col][i] if slope_se_col else "NA"),
                        ]
                    )
                file_written += n
                total_written += n
            print(
                f"[{idx}/{len(parquet_files)}] " f"{f.name}: " f"written={file_written:,}",
                flush=True,
            )
    return (total_scanned, total_written)


# ============================================================
# Build phenotype coordinate BED
#
# Only first 4 columns are parsed.
# Sample PSI matrix is NOT loaded into memory.
# ============================================================
def build_psi_bed(psi_files, sig_set, outfile):
    seen = set()
    written = 0
    with open(outfile, "w") as out:
        out.write("chr\tstart\tend\tphenotype_id\n")
        for f in psi_files:
            with gzip.open(f, "rt") as handle:
                for line_no, line in enumerate(handle):
                    if line_no == 0:
                        continue
                    # Only split first 4 columns +
                    # remainder of sample matrix.
                    parts = line.rstrip("\n").split("\t", 4)
                    if len(parts) < 4:
                        continue
                    chrom = parts[0]
                    start = parts[1]
                    end = parts[2]
                    phenotype = parts[3]
                    if phenotype not in sig_set:
                        continue
                    if phenotype in seen:
                        continue
                    seen.add(phenotype)
                    out.write(f"{chrom}\t" f"{start}\t" f"{end}\t" f"{phenotype}\n")
                    written += 1
    return (written, seen)


# ============================================================
# Main
# ============================================================
summary = []
for tissue in TISSUES:
    print("\n" + "=" * 90, flush=True)
    print(tissue, flush=True)
    print("=" * 90, flush=True)
    t0 = time.time()
    # --------------------------------------------------------
    # Resolve paths
    # --------------------------------------------------------
    fdr_file = find_fdr(tissue)
    parquet_files = find_parquets(tissue)
    psi_files = find_psi_beds(tissue)
    print("FDR:", fdr_file, flush=True)
    print("Parquet files:", len(parquet_files), flush=True)
    print("PSI BED files:", len(psi_files), flush=True)
    # --------------------------------------------------------
    # Output directory
    # --------------------------------------------------------
    tissue_out = OUT / tissue
    tissue_out.mkdir(parents=True, exist_ok=True)
    # --------------------------------------------------------
    # Significant transcripts
    # --------------------------------------------------------
    sig_tx = load_significant_transcripts(fdr_file)
    sig_set = set(sig_tx)
    tx_file = tissue_out / f"{tissue}.significant_tu_transcripts.txt"
    with open(tx_file, "w") as out:
        for tx in sig_tx:
            out.write(tx + "\n")
    print("Significant tu transcripts:", f"{len(sig_tx):,}", flush=True)
    # --------------------------------------------------------
    # PSI phenotype BED
    # --------------------------------------------------------
    psi_out = tissue_out / f"{tissue}.tuQTL.sig.PSI.bed"
    (n_bed, bed_seen) = build_psi_bed(psi_files, sig_set, psi_out)
    missing_bed = sig_set - bed_seen
    print("PSI BED transcripts:", f"{n_bed:,}", flush=True)
    print("Significant transcripts absent from PSI BED:", f"{len(missing_bed):,}", flush=True)
    if missing_bed:
        missing_file = tissue_out / f"{tissue}.missing_from_PSI_bed.txt"
        with open(missing_file, "w") as out:
            for tx in sorted(missing_bed):
                out.write(tx + "\n")
    # --------------------------------------------------------
    # Hard QC
    # --------------------------------------------------------
    if missing_bed:
        raise RuntimeError(
            f"{tissue}: "
            f"{len(missing_bed)} significant transcripts "
            f"are missing from PSI BED."
        )
    # --------------------------------------------------------
    # Nominal associations
    # --------------------------------------------------------
    nominal_out = tissue_out / f"{tissue}.tuQTL.sig.nominal.tsv.gz"
    (scanned, written) = extract_nominal(parquet_files, sig_tx, nominal_out)
    runtime = (time.time() - t0) / 60
    # --------------------------------------------------------
    # Summary
    # --------------------------------------------------------
    summary.append(
        {
            "tissue": tissue,
            "fdr_file": str(fdr_file),
            "n_significant_transcripts": len(sig_tx),
            "n_PSI_bed_transcripts": n_bed,
            "n_missing_PSI_bed": len(missing_bed),
            "nominal_rows_scanned": scanned,
            "nominal_rows_written": written,
            "runtime_min": runtime,
        }
    )
    print("Nominal scanned:", f"{scanned:,}", flush=True)
    print("Nominal written:", f"{written:,}", flush=True)
    print("Runtime:", f"{runtime:.1f} min", flush=True)
# ============================================================
# Final summary
# ============================================================
summary_df = pd.DataFrame(summary)
summary_file = OUT / "smr_input_summary.tsv"
summary_df.to_csv(summary_file, sep="\t", index=False)
print("\n" + "=" * 100, flush=True)
print("FINAL SUMMARY", flush=True)
print("=" * 100, flush=True)
print(summary_df.to_string(index=False), flush=True)
print("\nOutput directory:", flush=True)
print(OUT, flush=True)
print("\nSMR input preparation completed.", flush=True)
