#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 05_DNA_sequence_model
# Purpose: Evaluate coverage and end usage.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

import pathlib, json, hashlib, time, copy, csv
import numpy as np, pandas as pd, tensorflow as tf
from scipy.stats import pearsonr, spearmanr
from baskerville import seqnn, dataset, metrics, trainer
from test_metrics_33 import Evaluator

R = pathlib.Path("/path/to/chicken-aging/model")
P = R / "experiment33"
D = P / "data/tfrecords_v2"
O = P / "evaluations/medium33"
A = P / "gpu_scripts"
assert not O.exists(), "Preserve previous evaluation; use a new version if rerunning"
assert tf.config.list_physical_devices("GPU"), "GPU unavailable"
O.mkdir(parents=True)
tf.keras.mixed_precision.set_global_policy("float32")
tf.keras.utils.set_random_seed(20261013)
masks = json.load(open(P / "test_masks_v1/test_masks.json"))
targets = pd.read_csv(D / "targets.txt", sep="\t", index_col=0)
assert masks["tissues"] == targets.identifier.tolist() and len(targets) == 33
assert (
    hashlib.sha256((D / "sequences.bed").read_bytes()).hexdigest() == masks["sequence_bed_sha256"]
)
assert masks["n_windows"] == 223
(O / "test_masks.json").write_text(json.dumps(masks))
lossfn = trainer.parse_loss("poisson_mn", total_weight=0.2, weight_range=8, weight_exp=6)
models = [("medium33", A / "params_medium_33.json", P / "runs/medium_33_full_v1/model_best.h5", 33)]
summaries = []
provenance = []
reference_y = None


def assess(label, p, y, mm):
    ev = Evaluator(mm)
    scores, rs, arrays, counts = ev.evaluate(p, y)
    control = Evaluator(mm).evaluate(y, y)[0]
    assert all(abs(v - 1) < 1e-8 for v in control.values() if np.isfinite(v)), control
    folder = O / label
    folder.mkdir()
    np.savez_compressed(folder / "values.npz", **arrays, tissues=np.array(mm["tissues"]))
    rows = []
    for metric, values in rs.items():
        rows.extend(
            dict(metric=metric, tissue=t, value=float(v)) for t, v in zip(mm["tissues"], values)
        )
    pg = arrays["gene_prediction"]
    yg = arrays["gene_observed"]
    keep = arrays["gene_keep"]
    eps = np.array([len(g["bins"]) for g in mm["genes"]])[:, None] * 0.1
    pf = arrays["gene_prediction_fc"]
    yf = arrays["gene_observed_fc"]
    baseline = np.log2(yg.mean(1, keepdims=True) + eps)
    for ti, t in enumerate(mm["tissues"]):
        values = {
            "gene_log_expression_r": pearsonr(
                np.log2(pg[keep, ti] + eps[keep, 0]), np.log2(yg[keep, ti] + eps[keep, 0])
            )[0],
            "gene_fc_mse": np.mean((pf[keep, ti] - yf[keep, ti]) ** 2),
            "gene_zero_fc_baseline_mse": np.mean(yf[keep, ti] ** 2),
            "observed_tissue_mean_expression_baseline_r": pearsonr(
                baseline[keep, 0], np.log2(yg[keep, ti] + eps[keep, 0])
            )[0],
        }
        rows.extend(dict(metric=k, tissue=t, value=float(v)) for k, v in values.items())
    pd.DataFrame(rows).to_csv(folder / "per_tissue_metrics.tsv", sep="\t", index=False)
    gm = pd.DataFrame(
        [
            dict(
                gene=g["gene"],
                chrom=g["chrom"],
                window_index=g["window_index"],
                n_bins=len(g["bins"]),
                eligible=bool(keep[i]),
            )
            for i, g in enumerate(mm["genes"])
        ]
    )
    gm.to_csv(folder / "gene_metadata.tsv", sep="\t", index=False)
    em = pd.DataFrame(
        [
            dict(
                event=e["event"],
                gene=e["gene"],
                type=e["type"],
                window_index=e["window_index"],
                specific_bins=len(e["specific"]),
                common_bins=len(e["common"]),
                eligible=bool(arrays["event_keep"][i]),
                group_transcripts=",".join(e["group_transcripts"]),
                other_transcripts=",".join(e["other_transcripts"]),
            )
            for i, e in enumerate(mm["events"])
        ]
    )
    em.to_csv(folder / "event_metadata.tsv", sep="\t", index=False)
    for entity, meta, vals in [
        (
            "gene",
            gm,
            {
                "predicted_coverage": pg,
                "observed_coverage": yg,
                "predicted_log2_fc": pf,
                "observed_log2_fc": yf,
            },
        ),
        (
            "event",
            em,
            {
                "predicted_ratio": arrays["event_prediction_ratio"],
                "observed_ratio": arrays["event_observed_ratio"],
                "predicted_contrast": arrays["event_prediction_contrast"],
                "observed_contrast": arrays["event_observed_contrast"],
            },
        ),
    ]:
        chunks = []
        for ti, t in enumerate(mm["tissues"]):
            q = meta.copy()
            q["tissue"] = t
            for k, v in vals.items():
                q[k] = v[:, ti]
            chunks.append(q)
        pd.concat(chunks, ignore_index=True).to_csv(
            folder / f"{entity}_tissue_values.tsv.gz", sep="\t", index=False, compression="gzip"
        )
    predicted_top = pg[keep].argmax(1)
    observed_top = yg[keep].argmax(1)
    top = gm[keep].copy()
    top["predicted_top_tissue"] = np.array(mm["tissues"])[predicted_top]
    top["observed_top_tissue"] = np.array(mm["tissues"])[observed_top]
    top.to_csv(folder / "top_tissue.tsv", sep="\t", index=False)
    # Block bootstrap by held-out genomic window, not individual bins.
    rng = np.random.default_rng(20261020)
    ci = {}
    for key in ["gene_fc", "TSS_contrast", "APA_contrast"]:
        if key == "gene_fc":
            xx, yy, kk, meta = pf, yf, keep, gm
        else:
            kind = key.split("_")[0]
            xx = arrays["event_prediction_contrast"]
            yy = arrays["event_observed_contrast"]
            kk = arrays["event_keep"] & (em.type.to_numpy() == kind)
            meta = em
        blocks = np.unique(meta.window_index.to_numpy()[kk])
        indices = {b: np.flatnonzero(kk & (meta.window_index.to_numpy() == b)) for b in blocks}
        samples = []
        if kk.sum() >= 10 and len(blocks) >= 2:
            for _ in range(200):
                ix = np.concatenate(
                    [indices[b] for b in rng.choice(blocks, len(blocks), replace=True)]
                )
                samples.append(
                    np.nanmean([spearmanr(xx[ix, t], yy[ix, t])[0] for t in range(xx.shape[1])])
                )
            ci[key] = np.nanquantile(samples, [0.025, 0.975]).tolist()
        else:
            ci[key] = None
    result = dict(
        label=label,
        n_tissues=len(mm["tissues"]),
        **scores,
        counts=counts,
        window_bootstrap_ci95=ci,
        top_tissue_accuracy=float(np.mean(predicted_top == observed_top)),
        majority_top_tissue_baseline=float(
            np.bincount(observed_top, minlength=len(mm["tissues"])).max() / len(observed_top)
        ),
        identity_control=control,
        scope="Held-out chromosome coverage and annotation-defined terminal usage proxies; not isoform TPM or SNP causal validation",
    )
    (folder / "summary.json").write_text(json.dumps(result, indent=2))
    summaries.append(result)
    print("STRUCTURE_METRICS", json.dumps(result), flush=True)


for label, params, weights, nt in models:
    start = time.time()
    sha = hashlib.sha256(weights.read_bytes()).hexdigest()
    cfg = json.load(open(params))
    cfg["model"]["strand_pair"] = [np.arange(nt)]
    model = seqnn.SeqNN(cfg["model"])
    model.restore(str(weights), 0)
    ds = dataset.SeqDataset(str(D), "test", batch_size=1, mode="eval")
    predictions = np.empty((223, 12288, nt), np.float32)
    observed = np.empty_like(predictions)
    cr = metrics.PearsonR(nt)
    c2 = metrics.R2(nt)
    losses = []

    @tf.function
    def predict(x):
        return model.model(x, training=False)

    for wi, (x, y) in enumerate(ds.dataset):
        assert wi < 223
        y = y[:, :, :nt]
        q = predict(x)
        assert tuple(q.shape) == (1, 12288, nt)
        tf.debugging.assert_all_finite(q, "Nonfinite predictions")
        cr(y, q)
        c2(y, q)
        losses.append(float(lossfn(y, q).numpy()))
        predictions[wi] = q.numpy()[0]
        observed[wi] = y.numpy()[0]
        if wi % 50 == 0:
            print("PREDICT", label, wi, flush=True)
    assert wi + 1 == 223
    if reference_y is None:
        reference_y = observed.copy()
        np.save(O / "observed25_test.npy", observed.astype(np.float16))
    else:
        assert np.array_equal(observed, reference_y[:, :, :nt])
    np.save(O / f"{label}_predictions_test.npy", predictions)
    coverage_rows = []
    for ti, tissue in enumerate(masks["tissues"][:nt]):
        yy = observed[:, :, ti].reshape(-1).astype(np.float64)
        pp = predictions[:, :, ti].reshape(-1).astype(np.float64)
        coverage_rows.append(
            dict(
                tissue=tissue,
                pearson=float(pearsonr(yy, pp)[0]),
                r2=float(1 - np.sum((pp - yy) ** 2) / np.sum((yy - yy.mean()) ** 2)),
            )
        )
    pd.DataFrame(coverage_rows).to_csv(
        O / f"{label}_coverage_per_tissue.tsv", sep="\t", index=False
    )
    info = dict(
        model=label,
        weights=str(weights),
        sha256=sha,
        coverage_r=float(cr.result().numpy()),
        coverage_r2=float(c2.result().numpy()),
        loss=float(np.mean(losses)),
        seconds=time.time() - start,
        n_targets=nt,
    )
    provenance.append(info)
    print("COVERAGE", json.dumps(info), flush=True)
    mm = copy.deepcopy(masks)
    mm["tissues"] = mm["tissues"][:nt]
    assess(label, predictions, observed, mm)
    assert hashlib.sha256(weights.read_bytes()).hexdigest() == sha
    del model, ds, predictions, observed
    tf.keras.backend.clear_session()
pd.DataFrame(
    [
        {
            k: v
            for k, v in s.items()
            if k not in ["counts", "window_bootstrap_ci95", "identity_control"]
        }
        for s in summaries
    ]
).to_csv(O / "model_comparison.tsv", sep="\t", index=False)
(O / "COMPLETE.json").write_text(
    json.dumps(dict(status="COMPLETE", models=provenance, summaries=summaries), indent=2)
)
