#!/usr/bin/env python3
"""Extract matched THE-95/THE-101 score-shape diagnostics to compact JSON.

The script is intentionally read-only and writes JSON to stdout. It is designed
to run on SDCC where the completed same-source score caches are available.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.interpolate import UnivariateSpline
from sklearn.metrics import log_loss, roc_auc_score


PRODUCT = "centAsFeatBase3x3_pt15to35"
SCORE_KEY = f"score_{PRODUCT}"
CENTRALITY_BINS = ((0.0, 20.0), (20.0, 50.0), (50.0, 80.0))
PT_BINS = ((15.0, 20.0), (20.0, 25.0), (25.0, 35.0))
SCORE_BINS = np.linspace(0.0, 1.0, 51)
ETA_RANGE = (-0.7, 0.7)
N_WEIGHT_BINS = 20
ET_WEIGHT_CAP = 800.0

DEFAULT_BASE = Path(
    "/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/"
    "auauTightBDT_THE101_minbias_baseline14_extract_20260714_1630/reports"
)
DEFAULT_BEFORE = DEFAULT_BASE / (
    "model_validation_condor_the101_the95_same_source_fullstat_"
    "mlpy_recovery_20260714_2141"
)
DEFAULT_AFTER = DEFAULT_BASE / (
    "model_validation_condor_the101_candidate_same_source_fullstat_"
    "mlpy_recovery_20260714_2141"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--before-report", type=Path, default=DEFAULT_BEFORE)
    parser.add_argument("--after-report", type=Path, default=DEFAULT_AFTER)
    return parser.parse_args()


def cache_paths(report: Path) -> list[Path]:
    listed = report / "score_caches.list"
    if listed.is_file():
        paths = [Path(line.strip()) for line in listed.read_text().splitlines() if line.strip()]
    else:
        paths = sorted((report / "score_caches").glob("score_cache_*.npz"))
    if not paths:
        raise SystemExit(f"No score caches found under {report}")
    return paths


def normalize_mean_one(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype="float64")
    finite = np.isfinite(values) & (values > 0.0)
    if not finite.any():
        return np.ones(len(values), dtype="float64")
    mean = float(np.mean(values[finite]))
    out = np.ones(len(values), dtype="float64")
    out[finite] = values[finite] / mean
    return out


def inverse_pdf_weights(
    values: np.ndarray,
    *,
    fixed_range: tuple[float, float] | None,
    weight_cap: float | None,
) -> tuple[np.ndarray, dict]:
    values = np.asarray(values, dtype="float64")
    weights = np.ones(len(values), dtype="float64")
    finite = np.isfinite(values)
    report = {
        "n_entries": int(len(values)),
        "n_finite": int(finite.sum()),
        "n_bins": N_WEIGHT_BINS,
        "fixed_range": list(fixed_range) if fixed_range is not None else None,
        "weight_cap": weight_cap,
        "status": "ok",
    }
    if finite.sum() < max(10, N_WEIGHT_BINS):
        report["status"] = "insufficient_finite_values"
        return weights, report
    vals = values[finite]
    lo, hi = (float(np.min(vals)), float(np.max(vals))) if fixed_range is None else fixed_range
    report["range"] = [float(lo), float(hi)]
    edges = np.linspace(lo, hi, N_WEIGHT_BINS + 1)
    hist, bin_edges = np.histogram(vals, bins=edges, density=True)
    hist = hist.astype("float64") * N_WEIGHT_BINS
    centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    good = np.isfinite(hist)
    if good.sum() < 4:
        report["status"] = "insufficient_histogram_support"
        return weights, report
    spline = UnivariateSpline(centers[good], hist[good], s=0.0)
    pdf = np.clip(spline(vals), 1.0e-3, None)
    local = 1.0 / pdf
    if weight_cap is not None:
        local = np.clip(local, None, weight_cap)
    local = normalize_mean_one(local)
    weights[finite] = local
    report.update(
        min_weight=float(np.min(local)),
        max_weight=float(np.max(local)),
        mean_weight=float(np.mean(local)),
    )
    return weights, report


def display_weights(labels: np.ndarray, et: np.ndarray, eta: np.ndarray) -> tuple[np.ndarray, dict]:
    weights = np.ones(len(labels), dtype="float64")
    counts = {str(cls): int(np.sum(labels == cls)) for cls in (0, 1)}
    if min(counts.values()) <= 0:
        raise SystemExit(f"Both source-defined classes are required: {counts}")
    total = sum(counts.values())
    factors: dict[str, float] = {}
    eta_reports: dict[str, dict] = {}
    et_reports: dict[str, dict] = {}
    for cls in (0, 1):
        mask = labels == cls
        factor = total / (2.0 * counts[str(cls)])
        factors[str(cls)] = float(factor)
        weights[mask] *= factor
        eta_w, eta_report = inverse_pdf_weights(
            eta[mask], fixed_range=ETA_RANGE, weight_cap=None
        )
        et_w, et_report = inverse_pdf_weights(
            et[mask], fixed_range=None, weight_cap=ET_WEIGHT_CAP
        )
        weights[mask] *= eta_w * et_w
        eta_reports[str(cls)] = eta_report
        et_reports[str(cls)] = et_report
    if not np.all(np.isfinite(weights) & (weights > 0.0)):
        raise SystemExit("Display weights contain non-finite or non-positive values")
    return weights, {
        "weight_mode": "ppg12-exact-display-source-classes",
        "class_counts": counts,
        "class_weight_factors": factors,
        "eta_reweight": eta_reports,
        "et_reweight": et_reports,
        "weights_computed_before_binning": True,
        "event_weight_used": False,
        "cross_section_weight_used_for_training": False,
        "sum_weight_class0": float(np.sum(weights[labels == 0])),
        "sum_weight_class1": float(np.sum(weights[labels == 1])),
    }


def weighted_quantile(values: np.ndarray, weights: np.ndarray, q: float) -> float:
    good = np.isfinite(values) & np.isfinite(weights) & (weights > 0.0)
    if not good.any():
        return math.nan
    vals = values[good].astype("float64")
    wgt = weights[good].astype("float64")
    order = np.argsort(vals, kind="mergesort")
    vals = vals[order]
    cumulative = np.cumsum(wgt[order])
    index = min(np.searchsorted(cumulative, q * cumulative[-1], side="left"), len(vals) - 1)
    return float(vals[index])


def metrics(score: np.ndarray, labels: np.ndarray, weights: np.ndarray, mask: np.ndarray) -> dict:
    good = mask & np.isfinite(score) & np.isfinite(weights) & (weights > 0.0)
    y = labels[good]
    s = score[good].astype("float64")
    w = weights[good].astype("float64")
    if len(np.unique(y)) != 2:
        return {key: math.nan for key in ("auc", "logloss", "median_gap", "wp80_inclusive_fake", "wp80_threshold")}
    sig = y == 1
    inc = y == 0
    wp80 = weighted_quantile(s[sig], w[sig], 0.20)
    return {
        "auc": float(roc_auc_score(y, s, sample_weight=w)),
        "logloss": float(log_loss(y, np.clip(s, 1.0e-6, 1.0 - 1.0e-6), sample_weight=w, labels=[0, 1])),
        "median_signal_score": weighted_quantile(s[sig], w[sig], 0.50),
        "median_inclusive_score": weighted_quantile(s[inc], w[inc], 0.50),
        "median_gap": weighted_quantile(s[sig], w[sig], 0.50) - weighted_quantile(s[inc], w[inc], 0.50),
        "wp80_inclusive_fake": float(np.sum(w[inc & (s >= wp80)]) / np.sum(w[inc])),
        "wp80_threshold": wp80,
        "signal_entries": int(np.sum(sig)),
        "inclusive_entries": int(np.sum(inc)),
    }


def histogram(score: np.ndarray, labels: np.ndarray, mask: np.ndarray, cls: int) -> dict:
    values = score[mask & (labels == cls) & np.isfinite(score)]
    counts, _ = np.histogram(values, bins=SCORE_BINS)
    density = counts.astype("float64")
    if density.sum() > 0:
        density /= density.sum() * np.diff(SCORE_BINS)
    return {
        "counts": counts.astype(int).tolist(),
        "density": density.astype(float).tolist(),
        "entries": int(len(values)),
    }


def main() -> int:
    args = parse_args()
    before_paths = cache_paths(args.before_report)
    after_paths = cache_paths(args.after_report)
    if len(before_paths) != len(after_paths):
        raise SystemExit("Before/after score-cache counts differ")

    labels_parts: list[np.ndarray] = []
    et_parts: list[np.ndarray] = []
    eta_parts: list[np.ndarray] = []
    cent_parts: list[np.ndarray] = []
    before_parts: list[np.ndarray] = []
    after_parts: list[np.ndarray] = []

    for before_path, after_path in zip(before_paths, after_paths):
        before = np.load(before_path, allow_pickle=True)
        after = np.load(after_path, allow_pickle=True)
        for key in ("is_signal", "source_sample", "cluster_Et", "cluster_Eta", "centrality", SCORE_KEY):
            if key not in before.files or key not in after.files:
                raise SystemExit(f"Missing {key} in matched cache pair")
        for key in ("is_signal", "source_sample", "cluster_Et", "cluster_Eta", "centrality"):
            if not np.array_equal(before[key], after[key]):
                raise SystemExit(f"Matched cache rows differ for {key}: {before_path.name}")

        source = before["source_sample"].astype(str)
        truth_y = before["is_signal"].astype("int32", copy=False)
        et = before["cluster_Et"].astype("float32", copy=False)
        eta = before["cluster_Eta"].astype("float32", copy=False)
        cent = before["centrality"].astype("float32", copy=False)
        score_before = before[SCORE_KEY].astype("float32", copy=False)
        score_after = after[SCORE_KEY].astype("float32", copy=False)

        signal = (np.char.find(source, "embeddedPhoton") >= 0) & (truth_y == 1)
        inclusive = np.char.find(source, "embeddedJet") >= 0
        eligible = (
            (signal | inclusive)
            & np.isfinite(et)
            & (et >= 15.0)
            & (et < 35.0)
            & np.isfinite(eta)
            & (eta >= ETA_RANGE[0])
            & (eta < ETA_RANGE[1])
            & np.isfinite(cent)
            & (cent >= 0.0)
            & (cent < 80.0)
            & np.isfinite(score_before)
            & np.isfinite(score_after)
        )
        labels = np.where(signal, 1, 0).astype("int8")
        labels_parts.append(labels[eligible])
        et_parts.append(et[eligible])
        eta_parts.append(eta[eligible])
        cent_parts.append(cent[eligible])
        before_parts.append(score_before[eligible])
        after_parts.append(score_after[eligible])

    labels = np.concatenate(labels_parts)
    et = np.concatenate(et_parts)
    eta = np.concatenate(eta_parts)
    cent = np.concatenate(cent_parts)
    score_before = np.concatenate(before_parts)
    score_after = np.concatenate(after_parts)
    weights, weight_report = display_weights(labels, et, eta)
    all_mask = np.ones(len(labels), dtype=bool)

    models = {}
    for key, label, score in (
        ("before", "Before classifier requirement", score_before),
        ("after", "After classifier requirement", score_after),
    ):
        cells = {}
        auc_grid = {}
        for clo, chi in CENTRALITY_BINS:
            ckey = f"{clo:g}_{chi:g}"
            cmask = (cent >= clo) & (cent < chi)
            cells[ckey] = {
                "label": f"{clo:g}-{chi:g}%",
                "metrics": metrics(score, labels, weights, cmask),
                "signal": histogram(score, labels, cmask, 1),
                "inclusive": histogram(score, labels, cmask, 0),
            }
            for plo, phi in PT_BINS:
                pkey = f"{plo:g}_{phi:g}"
                pmask = (et >= plo) & (et < phi)
                auc_grid[f"pt_{pkey}_cent_{ckey}"] = metrics(
                    score, labels, weights, cmask & pmask
                )["auc"]
        models[key] = {
            "label": label,
            "inclusive_metrics": metrics(score, labels, weights, all_mask),
            "cells": cells,
            "auc_by_pt_centrality": auc_grid,
        }

    payload = {
        "schema": "THE101_MINBIAS_BEFORE_AFTER_SCORE_PAYLOAD_V1",
        "status": "READY",
        "product": PRODUCT,
        "selection": "15 <= cluster_Et < 35 GeV; 0 <= centrality < 80%; |eta| < 0.7",
        "plot_class_definition": (
            "Signal = embedded Photon12+20 rows with is_signal == 1; "
            "inclusive = all embedded Jet12+20+30+40 rows without a truth-background filter"
        ),
        "comparison_scope": "THE-95 and THE-101 evaluated on identical classifier-pass rows",
        "bin_edges": SCORE_BINS.astype(float).tolist(),
        "centrality_bins": [list(x) for x in CENTRALITY_BINS],
        "pt_bins": [list(x) for x in PT_BINS],
        "models": models,
        "weight_report": weight_report,
        "row_count": int(len(labels)),
        "signal_entries": int(np.sum(labels == 1)),
        "inclusive_entries": int(np.sum(labels == 0)),
        "source_reports": {
            "before": str(args.before_report),
            "after": str(args.after_report),
            "cache_pairs": len(before_paths),
        },
    }
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
