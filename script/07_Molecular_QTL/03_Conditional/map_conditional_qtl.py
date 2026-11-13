#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 07_Molecular_QTL
# Purpose: Conditional QTLs.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

import argparse, gzip, hashlib, inspect, json, os, re, time, traceback
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from scipy.stats import t as student
from tensorqtl import cis, genotypeio
from tensorqtl.core import Residualizer, calculate_corr
import importlib.metadata

ROOT = Path(__file__).resolve().parent


def atom(path, obj):
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(obj, indent=2))
    temp.replace(path)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def prepare(host):
    if host == "co":
        fm = Path("/path/to/chicken-aging/qtl/fine_mapping")
        tasks = json.loads((fm / "input_manifest.json").read_text())
        for task in tasks:
            task["susie"] = str(
                next(
                    (
                        fm
                        / "results_EstimateResidualVarianceFalse"
                        / task["tissue"]
                        / task["qtl_type"]
                    ).glob("*.SuSiE_summary.parquet")
                )
            )
    else:
        base = Path("/path/to/chicken-aging/qtl/inputs")
        tasks = []
        for tissue in json.loads((ROOT / "study_tissues.json").read_text(encoding="utf8"))[
            "partition_tissues"
        ]:
            for kind in ["teQTL", "tuQTL"]:
                if kind == "teQTL":
                    directory = base / "7.teQTL_mapping" / f"{tissue}_teQTL"
                    prefix = f"ChickenAging.{tissue}_teQTL"
                    bed = (
                        base
                        / "4.QTL_input/teQTL"
                        / f"ChickenGTEx.{tissue}.transcript.filtered.tensorqtl.bed.gz"
                    )
                else:
                    directory = base / "6.tuQTL_mapping" / f"{tissue}_PSI_qtl"
                    prefix = f"ChickenAging.{tissue}_PSI"
                    bed = (
                        base
                        / "4.QTL_input/tuQTL"
                        / f"salmon_pheno.ChickenGTEx.{tissue}.TPM_isoform.filtered.tensorqtl.bed.gz"
                    )
                tasks.append(
                    dict(
                        tissue=tissue,
                        qtl_type=kind,
                        genotype=str(base / "2.bfile" / f"ChickenGTEx.{tissue}"),
                        covariates=str(
                            base / "5.covariates_file" / f"{tissue}.combined_covariates.txt"
                        ),
                        phenotype_files=[str(bed)],
                        fdr=str(directory / "results/full" / f"{prefix}.cis_qtl.fdr.txt"),
                        susie=str(
                            next(
                                (
                                    base
                                    / "13.QTL_finemapping/results_EstimateResidualVarianceFalse"
                                    / tissue
                                    / kind
                                ).glob("*.SuSiE_summary.parquet")
                            )
                        ),
                    )
                )
    for task in tasks:
        task["full"] = task["fdr"].replace(".fdr.txt", ".txt")
        log = Path(task["fdr"]).parent
        log = next(log.glob("*.cis.log"))
        contents = log.read_text()
        assert "0.05 MAF filter" in contents and "1,000,000" in contents, str(log)
        task["mapping_log"] = str(log)
        task["samples"] = int(re.search(r"\* (\d+) samples", contents).group(1))
        task["n_covariates"] = int(re.search(r"\* (\d+) covariates", contents).group(1))
        for key, hashkey in [("fdr", "fdr_sha256"), ("covariates", "cov_sha256")]:
            if hashkey in task:
                assert digest(task[key]) == task[hashkey], task[key]
        files = [task[k] for k in ["full", "fdr", "covariates", "susie", "mapping_log"]]
        files += task["phenotype_files"] + [
            task["genotype"] + "." + ext for ext in ["bed", "bim", "fam"]
        ]
        task["current_files"] = {
            str(p): dict(size=Path(p).stat().st_size, mtime=Path(p).stat().st_mtime_ns)
            for p in files
        }
        task["hashes"] = {task[k]: digest(task[k]) for k in ["full", "fdr", "covariates"]}
        frame = pd.read_csv(task["full"], sep="\t", index_col=0)
        sig = frame[frame.qval <= 0.05]
        fdr = pd.read_csv(task["fdr"], sep="\t", index_col=0)
        assert sig.index.is_unique
        assert set(fdr.index) <= set(frame.index), task["fdr"]
        task["full_significant_count"] = len(sig)
        task["full_fdr_membership_equal"] = set(sig.index) == set(fdr.index)
        sig = fdr[fdr.qval <= 0.05]
        assert len(sig) == len(fdr)
        core_columns = [
            c for c in fdr.columns if c not in ["qval", "pval_nominal_threshold", "is_eGene"]
        ]
        pd.testing.assert_frame_equal(
            frame.loc[fdr.index, core_columns].sort_index(),
            fdr[core_columns].sort_index(),
            check_dtype=False,
        )
        task["targets"] = len(sig)
        cv = pd.read_csv(task["covariates"], sep="\t", index_col=0).T
        assert cv.shape == (task["samples"], task["n_covariates"])
        fam = pd.read_csv(task["genotype"] + ".fam", sep=r"\s+", header=None, dtype=str)
        assert set(cv.index) <= set(fam[1]) and cv.index.is_unique
        for bed in task["phenotype_files"]:
            with gzip.open(bed, "rt") as h:
                assert h.readline().rstrip().split("\t")[4:] == list(cv.index), bed
        print(task["tissue"], task["qtl_type"], len(sig), cv.shape, flush=True)
    atom(ROOT / "manifest.json", tasks)


def run(tissue, kind):
    task = next(
        t
        for t in json.loads((ROOT / "manifest.json").read_text())
        if t["tissue"] == tissue and t["qtl_type"] == kind
    )
    out = ROOT / "results" / tissue / kind
    out.mkdir(parents=True, exist_ok=True)
    if (out / "COMPLETE.json").exists():
        return
    for path, before in task["current_files"].items():
        now = Path(path).stat()
        assert (now.st_size, now.st_mtime_ns) == (before["size"], before["mtime"]), path
    for path, value in task["hashes"].items():
        assert digest(path) == value, path
    threads = int(os.environ.get("OMP_NUM_THREADS", "2"))
    torch.set_num_threads(threads)
    cv = pd.read_csv(task["covariates"], sep="\t", index_col=0).T
    full = pd.read_csv(task["fdr"], sep="\t", index_col=0)
    sig = full[full.qval <= 0.05]
    needed = set(sig.index)
    bychrom = {}
    for bed in task["phenotype_files"]:
        with gzip.open(bed, "rt") as handle:
            header = handle.readline().rstrip().split("\t")
            assert header[4:] == list(cv.index)
            for line in handle:
                row = line.rstrip().split("\t")
                if row[3] in needed:
                    bychrom.setdefault(row[0], []).append(
                        (
                            row[3],
                            int(row[1]) + 1,
                            int(row[2]),
                            np.asarray(row[4:], dtype=np.float32),
                        )
                    )
    allids = [row[0] for rows in bychrom.values() for row in rows]
    assert len(allids) == len(set(allids)) == len(sig) and set(allids) == needed
    pr = genotypeio.PlinkReader(task["genotype"], select_samples=list(cv.index), verbose=False)
    cs = pd.read_parquet(task["susie"])
    cs = cs[cs.phenotype_id.isin(needed)]
    assert cs.pip.between(0, 1.000001).all()
    cs["cs_id"] = cs.cs_id.astype(str)
    cslookup = {pid: list(frame.groupby("cs_id")) for pid, frame in cs.groupby("phenotype_id")}
    provenance = dict(
        task=task,
        version=importlib.metadata.version("tensorqtl"),
        maf=0.05,
        window=1000000,
        fdr=0.05,
        nperm=10000,
        seed=20261030,
        ld_r2=0.8,
        script_sha256=digest(__file__),
    )
    if (out / "provenance.json").exists():
        assert json.loads((out / "provenance.json").read_text()) == provenance
    else:
        atom(out / "provenance.json", provenance)
    rz = Residualizer(torch.tensor(cv.values, dtype=torch.float32))
    for chrom, rows in bychrom.items():
        stem = out / ("chr" + chrom)
        marker = Path(str(stem) + ".complete.json")
        if marker.exists():
            continue
        started = time.time()
        matrix, positions = pr.get_region(chrom)
        geno = pd.DataFrame(matrix, index=positions.index, columns=cv.index)
        variants = pd.DataFrame({"chrom": chrom, "pos": positions.values}, index=positions.index)
        y = pd.DataFrame(
            np.stack([row[3] for row in rows]), index=[row[0] for row in rows], columns=cv.index
        )
        pos = pd.DataFrame(
            {"chr": chrom, "start": [row[1] for row in rows], "end": [row[2] for row in rows]},
            index=y.index,
        )
        assert np.isfinite(y.values).all() and geno.index.is_unique

        def dose(variant):
            x = geno.loc[variant].to_numpy(dtype=np.float32).copy()
            x[x == -9] = np.nan
            x[np.isnan(x)] = np.nanmean(x)
            return x

        checks = []
        for pid in y.index:
            original = sig.loc[pid]
            variant = original.variant_id
            bp = variants.loc[variant, "pos"]
            assert bp - pos.loc[pid, "start"] == original.start_distance
            assert bp - pos.loc[pid, "end"] == original.end_distance
            r, vx, vy = calculate_corr(
                torch.tensor(dose(variant)[None, :]),
                torch.tensor(y.loc[[pid]].values),
                rz,
                return_var=True,
            )
            rho = r.item()
            p = 2 * student.sf(abs(rho) * np.sqrt(rz.dof / max(1 - rho * rho, 1e-30)), rz.dof)
            slope = (r * (vy / vx).sqrt()).item()
            p_ok = (original.pval_nominal == 0 and p < 1e-300) or np.isclose(
                p, original.pval_nominal, rtol=0.002, atol=1e-300
            )
            ok = bool(p_ok and np.isclose(slope, original.slope, rtol=0.002, atol=1e-7))
            checks.append((pid, variant, original.pval_nominal, p, original.slope, slope, ok))
        check = pd.DataFrame(
            checks,
            columns=[
                "phenotype_id",
                "variant_id",
                "original_p",
                "recomputed_p",
                "original_slope",
                "recomputed_slope",
                "pass",
            ],
        )
        check.to_csv(str(stem) + ".input_check.tsv", sep="\t", index=False)
        assert check["pass"].all(), f"Input reproduction failed {tissue} {kind} {chrom}"
        calls = []
        original_func = cis.calculate_beta_approx_pval

        def audit(*args, **kwargs):
            result = original_func(*args, **kwargs)
            loc = inspect.currentframe().f_back.f_locals
            calls.append(
                (loc.get("phenotype_id", ""), float(result[0]), float(loc["signif_threshold"]))
            )
            return result

        cis.calculate_beta_approx_pval = audit
        try:
            result = cis.map_independent(
                geno,
                variants,
                full,
                y,
                pos,
                cv,
                maf_threshold=0.05,
                fdr=0.05,
                nperm=10000,
                window=1000000,
                seed=20261030,
                verbose=False,
            )
        finally:
            cis.calculate_beta_approx_pval = original_func
        result.to_csv(str(stem) + ".signals.tsv", sep="\t", index=False)
        tests = pd.DataFrame(calls, columns=["phenotype_id", "candidate_pval_beta", "threshold"])
        tests.to_csv(str(stem) + ".conditional_tests.tsv", sep="\t", index=False)
        assert np.isfinite(tests.candidate_pval_beta).all(), "Nonfinite conditional tests"
        assert set(result.phenotype_id) == set(y.index), "Missing output phenotypes"
        assert not result.duplicated(["phenotype_id", "rank"]).any()
        assert np.isfinite(
            result[["pval_nominal", "pval_beta", "slope", "slope_se"]].to_numpy(float)
        ).all()
        matches = []
        edges = []
        regions = []
        for row in result.itertuples():
            candidate = []
            groups = cslookup.get(row.phenotype_id, [])
            x = dose(row.variant_id)
            x = x - x.mean()
            max_r2 = 0.0
            for cid, group in groups:
                direct = row.variant_id in set(group.variant_id)
                if direct:
                    score = 1.0
                    proxy = row.variant_id
                else:
                    values = np.stack([dose(v) for v in group.variant_id.unique()])
                    values -= values.mean(axis=1, keepdims=True)
                    scores = (values @ x) ** 2 / ((values * values).sum(axis=1) * (x @ x))
                    ix = int(np.nanargmax(scores))
                    score = float(scores[ix])
                    proxy = group.variant_id.unique()[ix]
                max_r2 = max(max_r2, score)
                edges.append(
                    (row.phenotype_id, row.rank, row.variant_id, cid, direct, proxy, score)
                )
                if direct or score >= 0.8:
                    candidate.append((cid, direct))
            category = (
                "no_CS"
                if not groups
                else (
                    "multiple_CS"
                    if len(candidate) > 1
                    else (
                        "direct"
                        if any(c[1] for c in candidate)
                        else "LD_only" if candidate else "no_match"
                    )
                )
            )
            chrom_positions = variants.pos.to_numpy()
            mask = (chrom_positions >= pos.loc[row.phenotype_id, "start"] - 1000000) & (
                chrom_positions <= pos.loc[row.phenotype_id, "end"] + 1000000
            )
            block = geno.to_numpy()[mask].astype(np.float32)
            block[block == -9] = np.nan
            means = np.nanmean(block, axis=1)
            block = np.where(np.isnan(block), means[:, None], block)
            af = means / 2
            block -= block.mean(axis=1, keepdims=True)
            with np.errstate(invalid="ignore", divide="ignore"):
                r2 = (block @ x) ** 2 / ((block * block).sum(axis=1) * (x @ x))
            selected = chrom_positions[mask][(r2 >= 0.8) & (af >= 0.05) & (af <= 0.95)]
            selected = np.append(selected, variants.loc[row.variant_id, "pos"])
            regions.append(
                (
                    row.phenotype_id,
                    row.rank,
                    row.variant_id,
                    chrom,
                    int(selected.min()),
                    int(selected.max()),
                    bool(candidate),
                )
            )
            matches.append(
                (
                    row.phenotype_id,
                    row.rank,
                    row.variant_id,
                    category,
                    ",".join(c[0] for c in candidate),
                    max_r2,
                )
            )
        pd.DataFrame(
            matches,
            columns=["phenotype_id", "rank", "variant_id", "category", "matched_cs", "max_r2"],
        ).to_csv(str(stem) + ".matches.tsv", sep="\t", index=False)
        pd.DataFrame(
            edges,
            columns=["phenotype_id", "rank", "variant_id", "cs_id", "direct", "proxy", "max_r2"],
        ).to_csv(str(stem) + ".edges.tsv", sep="\t", index=False)
        pd.DataFrame(
            regions,
            columns=["phenotype_id", "rank", "variant_id", "chr", "start", "end", "CS_supported"],
        ).to_csv(str(stem) + ".regions.tsv", sep="\t", index=False)
        atom(marker, dict(phenotypes=len(rows), signals=len(result), seconds=time.time() - started))
        print("DONE_CHR", tissue, kind, chrom, len(result), flush=True)
        del geno, matrix, y, result
    signals = pd.concat(
        [pd.read_csv(p, sep="\t") for p in out.glob("*.signals.tsv")], ignore_index=True
    )
    matches = pd.concat(
        [pd.read_csv(p, sep="\t") for p in out.glob("*.matches.tsv")], ignore_index=True
    )
    edges = pd.concat(
        [pd.read_csv(p, sep="\t", dtype={"cs_id": str}) for p in out.glob("*.edges.tsv")],
        ignore_index=True,
    )
    assert set(signals.phenotype_id) == needed
    for name, frame in [
        ("cis_independent", signals),
        ("signal_CS_matches", matches),
        ("signal_CS_all_edges", edges),
    ]:
        frame.to_csv(out / (name + ".tsv.gz"), sep="\t", index=False)
    paired = edges[edges.direct | (edges.max_r2 >= 0.8)]
    counts = dict(
        tissue=tissue,
        qtl_type=kind,
        targets=len(sig),
        signals=len(signals),
        rank_ge2=int((signals["rank"] >= 2).sum()),
        unique_representative_snps=signals.variant_id.nunique(),
        cs=len(cs[["phenotype_id", "cs_id"]].drop_duplicates()),
        matched_cs=len(paired[["phenotype_id", "cs_id"]].drop_duplicates()),
        categories=matches.category.value_counts().to_dict(),
    )
    atom(out / "COMPLETE.json", counts)
    print("COMPLETE", counts, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--prepare", choices=["co", "tj"])
    parser.add_argument("--tissue")
    parser.add_argument("--qtl")
    args = parser.parse_args()
    if args.prepare:
        prepare(args.prepare)
    else:
        run(args.tissue, args.qtl)
