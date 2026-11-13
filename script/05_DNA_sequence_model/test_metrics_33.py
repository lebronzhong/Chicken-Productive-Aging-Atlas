#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 05_DNA_sequence_model
# Purpose: Test metrics 33.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

"""Shared metric definitions for frozen held-out test evaluation; no fitting."""

import numpy as np
from scipy.stats import spearmanr


def inverse(x):
    x = np.asarray(x, dtype=np.float64)
    x = np.where(x > 384, (x - 383) ** 2 + 383, x)
    return ((x + 1) ** 2 - 1) / 9.6


def contrast(x):
    x = np.log2(x)
    return (x - x.mean(1, keepdims=True)) * x.shape[1] / (x.shape[1] - 1)


def rho_columns(p, y, keep):
    if np.sum(keep) < 10:
        return np.full(p.shape[1], np.nan)
    return np.array([spearmanr(p[keep, t], y[keep, t])[0] for t in range(p.shape[1])])


class Evaluator:
    def __init__(self, masks):
        self.masks = masks
        self.frozen = None

    def evaluate(self, p, y):
        assert p.shape == y.shape == (self.masks["n_windows"], 12288, len(self.masks["tissues"]))
        assert np.isfinite(p).all() and np.isfinite(y).all()
        pg = []
        yg = []
        gene_ok = []
        lengths = []
        for g in self.masks["genes"]:
            a = p[g["window_index"], g["bins"]]
            b = y[g["window_index"], g["bins"]]
            pa = inverse(a)
            ob = inverse(b)
            pg.append(pa.sum(0))
            yg.append(ob.sum(0))
            lengths.append(len(g["bins"]))
            gene_ok.append((ob.mean(0).max() >= 1) and not np.any(b >= 768))
        pg = np.array(pg)
        yg = np.array(yg)
        eps = np.array(lengths)[:, None] * 0.1

        def fc(x):
            return np.log2((x + eps) / ((x.sum(1, keepdims=True) - x) / (pg.shape[1] - 1) + eps))

        pf = fc(pg)
        yf = fc(yg)
        er = []
        ey = []
        event_ok = []
        for e in self.masks["events"]:
            wi = e["window_index"]
            s = e["specific"]
            c = e["common"]
            a = p[wi]
            b = y[wi]
            ps = inverse(a[s]).mean(0)
            pc = inverse(a[c]).mean(0)
            ys = inverse(b[s]).mean(0)
            yc = inverse(b[c]).mean(0)
            er.append((ps + 0.1) / (pc + 0.1))
            ey.append((ys + 0.1) / (yc + 0.1))
            event_ok.append(bool(np.all(yc >= 1) and not np.any(b[s + c] >= 768)))
        er = np.array(er)
        ey = np.array(ey)
        ec = contrast(er)
        yc = contrast(ey)
        frozen = (np.array(gene_ok), np.array(event_ok))
        if self.frozen is None:
            self.frozen = frozen
        else:
            assert all(
                np.array_equal(a, b) for a, b in zip(self.frozen, frozen)
            ), "Observation eligibility changed"
        gkeep, ekeep = self.frozen
        assert gkeep.sum() >= 50
        types = np.array([e["type"] for e in self.masks["events"]])
        rs = {}
        scores = {}
        rs["gene_fc"] = rho_columns(pf, yf, gkeep)
        scores["gene_fc"] = float(np.mean(rs["gene_fc"]))
        for kind in ["TSS", "APA"]:
            take = ekeep & (types == kind)
            rs[kind + "_contrast"] = rho_columns(ec, yc, take)
            scores[kind + "_contrast"] = float(np.mean(rs[kind + "_contrast"]))
            rs[kind + "_ratio"] = rho_columns(er, ey, take)
            scores[kind + "_ratio"] = float(np.mean(rs[kind + "_ratio"]))
        arrays = dict(
            gene_prediction=pg,
            gene_observed=yg,
            gene_prediction_fc=pf,
            gene_observed_fc=yf,
            gene_keep=gkeep,
            event_prediction_ratio=er,
            event_observed_ratio=ey,
            event_prediction_contrast=ec,
            event_observed_contrast=yc,
            event_keep=ekeep,
        )
        counts = {
            "genes": int(gkeep.sum()),
            "TSS": int((ekeep & (types == "TSS")).sum()),
            "APA": int((ekeep & (types == "APA")).sum()),
        }
        return scores, rs, arrays, counts
