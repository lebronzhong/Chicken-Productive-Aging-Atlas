#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 07_Molecular_QTL
# Purpose: SuSiE fine-mapping.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

"""Bounded-memory TensorQTL SuSiE; preserve original interval windows and sample order."""

import argparse, contextlib, csv, gc, gzip, hashlib, importlib.metadata, json, os, pickle, re, shutil, sys, time, traceback
from pathlib import Path
import numpy as np
import pandas as pd
import torch
import dask
from tensorqtl import genotypeio, susie
from tensorqtl.core import Residualizer, impute_mean, calculate_maf

PARAMS = dict(
    L=10,
    scaled_prior_variance=0.2,
    estimate_residual_variance=False,
    estimate_prior_variance=True,
    tol=0.001,
    coverage=0.95,
    min_abs_corr=0.5,
    max_iter=500,
)
SUMMARY = ["phenotype_id", "variant_id", "pip", "af", "cs_id"]
STATUS = [
    "phenotype_id",
    "chr",
    "start",
    "end",
    "status",
    "variants",
    "var_y",
    "var_y_resid",
    "sigma2",
    "niter",
    "seconds",
    "error",
]


def atom_json(path, obj):
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(obj, indent=2, ensure_ascii=False) + "\n")
    tmp.replace(path)


def source_identity(task):
    for name, old in task["files"].items():
        st = Path(name).stat()
        if st.st_size != old["size"] or st.st_mtime_ns != old["mtime_ns"]:
            raise ValueError("Input changed since preflight: " + name)
    for key, hashkey in [("fdr", "fdr_sha256"), ("covariates", "cov_sha256")]:
        assert hashlib.sha256(Path(task[key]).read_bytes()).hexdigest() == task[hashkey]


def tensor_numpy(x):
    return x.detach().cpu().numpy() if torch.is_tensor(x) else x


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--tissue", required=True)
    ap.add_argument("--qtl", choices=["teQTL", "tuQTL"], required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument(
        "--smoke",
        action="store_true",
        help="One locus in separate output, compare with native susie.map",
    )
    args = ap.parse_args()
    torch.set_num_threads(2)
    dask.config.set(scheduler="synchronous")
    task = next(
        x
        for x in json.loads(Path(args.manifest).read_text())
        if x["tissue"] == args.tissue and x["qtl_type"] == args.qtl
    )
    source_identity(task)
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    prefix = f"ChickenAging.{args.tissue}.{args.qtl}"
    provenance = dict(
        task=task,
        parameters=PARAMS,
        window=1000000,
        maf_threshold=0,
        fdr=0.05,
        tensorqtl=importlib.metadata.version("tensorqtl"),
        smoke=args.smoke,
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    )
    if (out / "provenance.json").exists():
        assert (
            json.loads((out / "provenance.json").read_text()) == provenance
        ), "Existing run has different inputs/code/parameters"
    else:
        atom_json(out / "provenance.json", provenance)
    if (out / "COMPLETE.json").exists():
        print("Already complete", out, flush=True)
        return
    fdr = pd.read_csv(task["fdr"], sep="\t")
    fdr = fdr.loc[fdr.qval <= 0.05].set_index("phenotype_id")
    targets = set(fdr.index)
    cov = pd.read_csv(task["covariates"], sep="\t", index_col=0).T
    assert cov.index.is_unique and np.isfinite(cov.to_numpy(dtype=float)).all()
    cx = cov.to_numpy(dtype=float)
    centered = cx - cx.mean(axis=0)
    assert (
        np.linalg.matrix_rank(centered) == cx.shape[1]
    ), "Rank-deficient covariates require review"
    residualizer = Residualizer(torch.tensor(cx, dtype=torch.float32))
    pr = genotypeio.PlinkReader(task["genotype"], select_samples=cov.index.tolist(), verbose=False)
    assert pr.bim["snp"].is_unique, "Duplicate genotype variant IDs"
    chunks = out / "chromosomes"
    chunks.mkdir(exist_ok=True)
    seen = set()
    smoke_done = False
    for bedpath in sorted(
        task["phenotype_files"], key=lambda s: int(re.search(r"\.chr(\d+)\.bed", s).group(1))
    ):
        label = re.search(r"\.chr([^.]*)\.bed", bedpath).group(1)
        stem = chunks / f"chr{label}"
        marker = Path(str(stem) + ".done.json")
        if marker.exists():
            old = json.loads(marker.read_text())
            ids = set(old["phenotypes"])
            assert not (seen & ids), "Duplicate phenotype across chromosomes"
            seen |= ids
            print("RESUME", label, len(ids), flush=True)
            continue
        if shutil.disk_usage(out).free < 3 * 1024**3:
            raise RuntimeError("Less than 3 GiB free disk; stopping safely")
        selected = []
        with gzip.open(bedpath, "rt") as h:
            header = h.readline().rstrip("\r\n").split("\t")
            assert header[4:] == cov.index.tolist()
            for line in h:
                fields = line.rstrip("\r\n").split("\t")
                if fields[3] in targets:
                    assert fields[3] not in seen, "Duplicate phenotype in BEDs"
                    seen.add(fields[3])
                    selected.append(fields)
        if args.smoke:
            selected = selected[:1]
        print("CHROMOSOME", args.tissue, args.qtl, label, "targets", len(selected), flush=True)
        output_paths = {
            k: Path(str(stem) + s)
            for k, s in [
                ("summary", ".summary.parquet"),
                ("status", ".status.tsv"),
                ("records", ".SuSiE_records.pickle.gz"),
            ]
        }
        tmp = {k: p.with_name(p.name + ".tmp") for k, p in output_paths.items()}
        summaries = []
        status_rows = []
        if selected:
            chrom = selected[0][0]
            assert all(r[0] == chrom for r in selected)
            assert chrom in pr.chrs, ("No genotype chromosome", chrom)
            geno, positions = pr.get_region(chrom)
            order = np.argsort(positions.to_numpy(), kind="stable")
            geno = geno[order]
            positions = positions.iloc[order]
            bp = positions.to_numpy()
            variant_ids = positions.index.to_numpy()
        with (
            gzip.open(tmp["records"], "wb", compresslevel=1) as records,
            tmp["status"].open("w") as st,
        ):
            writer = csv.DictWriter(st, fieldnames=STATUS, delimiter="\t")
            writer.writeheader()
            for ix, fields in enumerate(selected):
                chrom, start0, end, pid = fields[:4]
                start = int(start0) + 1
                end = int(end)
                assert 1 <= start <= end, ("Invalid BED interval", pid)
                y = np.asarray(fields[4:], dtype=np.float32)
                info = dict(
                    phenotype_id=pid,
                    chr=chrom,
                    start=start,
                    end=end,
                    status="failed",
                    variants=0,
                    var_y=None,
                    var_y_resid=None,
                    sigma2=None,
                    niter=None,
                    seconds=0,
                    error="",
                )
                began = time.time()
                try:
                    if not np.isfinite(y).all() or np.var(y) == 0:
                        raise ValueError("Non-finite or constant phenotype")
                    left = np.searchsorted(bp, start - 1000000, side="left")
                    right = np.searchsorted(bp, end + 1000000, side="right")
                    if left == right:
                        raise ValueError("No genotypes in interval window")
                    vid = variant_ids[left:right]
                    if str(fdr.loc[pid, "variant_id"]) not in vid:
                        raise ValueError("FDR lead variant absent from current window/genotypes")
                    gt = torch.tensor(geno[left:right], dtype=torch.float32)
                    impute_mean(gt)
                    keep = ~(gt == gt[:, [0]]).all(1)
                    if not keep.any():
                        raise ValueError("All variants monomorphic")
                    gt = gt[keep]
                    vid = vid[keep.numpy()]
                    if not torch.isfinite(gt).all():
                        raise ValueError("Non-finite genotype after imputation")
                    yr = residualizer.transform(torch.tensor(y).reshape(1, -1))
                    xr = residualizer.transform(gt)
                    info["variants"] = len(vid)
                    info["var_y"] = float(np.var(y, ddof=1))
                    info["var_y_resid"] = float(yr.var())
                    if not torch.isfinite(xr).all() or info["var_y_resid"] <= 0:
                        raise ValueError("Degenerate covariate-adjusted inputs")
                    with torch.no_grad():
                        result = susie.susie(xr.T, yr.T, **PARAMS)
                    pip = np.asarray(result["pip"])
                    assert np.isfinite(pip).all() and ((pip >= 0) & (pip <= 1 + 1e-6)).all()
                    af = (gt.sum(1) / (2 * gt.shape[1])).numpy()
                    pipdf = pd.DataFrame(
                        {"pip": pip, "af": af}, index=pd.Index(vid, name="variant_id")
                    )
                    info["sigma2"] = float(result["sigma2"])
                    info["niter"] = int(result["niter"])
                    sets = result["sets"]["cs"]
                    info["status"] = (
                        "nonconverged"
                        if not result["converged"]
                        else ("ok" if sets is not None else "no_credible_set")
                    )
                    if info["status"] == "ok":
                        for csid, indices in sets.items():
                            frame = pipdf.iloc[indices].reset_index()
                            frame["phenotype_id"] = pid
                            frame["cs_id"] = csid.lstrip("L")
                            summaries.append(frame[SUMMARY])
                    record = {
                        k: tensor_numpy(result[k])
                        for k in ["sets", "converged", "elbo", "niter", "lbf_variable"]
                    }
                    record.update(
                        phenotype_id=pid,
                        pip=pipdf,
                        sigma2=info["sigma2"],
                        var_y_resid=info["var_y_resid"],
                    )
                    pickle.dump(record, records, protocol=4)
                    if args.smoke:
                        native_g = pd.DataFrame(
                            geno[left:right],
                            index=pd.Index(variant_ids[left:right], name="snp"),
                            columns=cov.index,
                        )
                        native_v = pd.DataFrame(
                            {"chrom": chrom, "pos": bp[left:right]}, index=native_g.index
                        )
                        native_y = pd.DataFrame([y], index=[pid], columns=cov.index)
                        native_p = pd.DataFrame(
                            {"chr": [chrom], "start": [start], "end": [end]}, index=[pid]
                        )
                        with open(os.devnull, "w") as null, contextlib.redirect_stdout(null):
                            native_summary, native = susie.map(
                                native_g,
                                native_v,
                                native_y,
                                native_p,
                                cov,
                                summary_only=False,
                                window=1000000,
                                maf_threshold=0,
                                **PARAMS,
                            )
                        if pid in native:
                            np.testing.assert_allclose(
                                native[pid]["pip"]["pip"].to_numpy(), pip, rtol=1e-5, atol=1e-6
                            )
                            assert set(native[pid]["sets"]["cs"]) == set(sets)
                            for key in sets:
                                np.testing.assert_array_equal(
                                    native[pid]["sets"]["cs"][key], sets[key]
                                )
                        else:
                            assert sets is None
                        print(
                            "SMOKE native TensorQTL PIP/credible-set comparison PASS",
                            pid,
                            flush=True,
                        )
                        smoke_done = True
                    del gt, xr, yr, result, record, pipdf
                except (ValueError, RuntimeError, FloatingPointError) as e:
                    info["error"] = str(e).replace("\t", " ").replace("\n", " ")[:1000]
                    with (out / "errors.log").open("a") as eh:
                        eh.write(pid + "\n" + traceback.format_exc() + "\n")
                    print("FAILED", pid, info["error"], flush=True)
                info["seconds"] = round(time.time() - began, 3)
                writer.writerow(info)
                st.flush()
                status_rows.append(info)
                if ix % 25 == 0 or ix == len(selected) - 1:
                    print("PROGRESS", label, ix + 1, len(selected), pid, info["status"], flush=True)
        frame = (
            pd.concat(summaries, ignore_index=True) if summaries else pd.DataFrame(columns=SUMMARY)
        )
        frame.to_parquet(tmp["summary"], index=False)
        for key in tmp:
            tmp[key].replace(output_paths[key])
        atom_json(
            marker,
            dict(
                phenotypes=[r[3] for r in selected],
                statuses=dict(
                    pd.Series([r["status"] for r in status_rows], dtype=str)
                    .value_counts()
                    .astype(int)
                    .items()
                ),
            ),
        )
        if selected:
            del geno, positions
        del summaries, frame
        gc.collect()
        if args.smoke:
            if smoke_done:
                break
            if selected:
                raise RuntimeError("Smoke locus failed")
    if args.smoke:
        assert smoke_done, "No smoke locus found"
        atom_json(out / "SMOKE_PASS.json", {"passed": True})
        return
    missing = targets - seen
    if missing:
        atom_json(out / "missing_phenotypes.json", sorted(missing))
        raise ValueError(f"{len(missing)} FDR phenotypes absent from archived BEDs")
    frames = [pd.read_parquet(p) for p in chunks.glob("*.summary.parquet")]
    summary = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=SUMMARY)
    summary.to_parquet(out / (prefix + ".SuSiE_summary.parquet"), index=False)
    summary.to_csv(out / (prefix + ".SuSiE_summary.tsv.gz"), sep="\t", index=False)
    statuses = pd.concat(
        [pd.read_csv(p, sep="\t") for p in chunks.glob("*.status.tsv")], ignore_index=True
    )
    assert statuses.phenotype_id.is_unique and set(statuses.phenotype_id) == targets
    statuses.to_csv(out / "phenotype_status.tsv", sep="\t", index=False)
    counts = dict(
        targets=len(targets),
        statuses={str(k): int(v) for k, v in statuses.status.value_counts().items()},
        credible_sets=len(summary[["phenotype_id", "cs_id"]].drop_duplicates()),
        phenotypes_with_cs=summary.phenotype_id.nunique(),
        unique_variants=summary.variant_id.nunique(),
        parameters=PARAMS,
    )
    bad = statuses.status.isin(["failed", "nonconverged"]).sum()
    atom_json(out / ("COMPLETE_WITH_FAILURES.json" if bad else "COMPLETE.json"), counts)
    print("TASK_RESULT", json.dumps(counts), flush=True)
    if bad:
        sys.exit(2)


if __name__ == "__main__":
    main()
