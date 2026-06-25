#!/usr/bin/env python3
"""Derive THE-57 WP80 slide payload from the full weighted training matrix."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

CENT_EDGES = [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0]
ET_RANGE = (15.0, 35.0)
ET_EDGES = [15.0, 17.0, 19.0, 21.0, 23.0, 25.0, 27.0, 30.0, 35.0]
TARGETS = [0.90, 0.80, 0.70]
TARGET_NAMES = {0.90: "WP90", 0.80: "WP80", 0.70: "WP70"}
FEATURES = [
    "cluster_Et",
    "cluster_weta_cogx",
    "cluster_wphi_cogx",
    "cluster_weta33_cogx",
    "cluster_wphi33_cogx",
    "vertexz",
    "cluster_Eta",
    "e11_over_e33",
    "cluster_et1",
    "cluster_et2",
    "cluster_et3",
    "cluster_et4",
    "e32_over_e35",
    "centrality",
]

PPG12_EXACT_ETA_RANGE = (-0.7, 0.7)
PPG12_EXACT_N_BINS = 20
PPG12_EXACT_ET_WEIGHT_CAP = 800.0


def normalize_mean_one(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype="float64")
    finite = np.isfinite(values)
    if not finite.any():
        return values
    mean = float(np.mean(values[finite]))
    if mean <= 0.0 or not math.isfinite(mean):
        return values
    out = values.copy()
    out[finite] = out[finite] / mean
    return out


def inverse_pdf_weights(values: np.ndarray, *, n_bins: int, fixed_range=None, weight_cap: float | None = None):
    from scipy.interpolate import UnivariateSpline

    values = np.asarray(values, dtype="float64")
    weights = np.ones(len(values), dtype="float64")
    finite = np.isfinite(values)
    report = {
        "n_entries": int(len(values)),
        "n_finite": int(finite.sum()),
        "n_bins": int(n_bins),
        "fixed_range": list(fixed_range) if fixed_range is not None else None,
        "weight_cap": float(weight_cap) if weight_cap is not None else None,
        "status": "ok",
    }
    if finite.sum() < max(10, n_bins):
        report["status"] = "insufficient_finite_values"
        return weights, report
    vals = values[finite]
    lo, hi = (float(np.min(vals)), float(np.max(vals))) if fixed_range is None else (float(fixed_range[0]), float(fixed_range[1]))
    report["range"] = [lo, hi]
    if not math.isfinite(lo) or not math.isfinite(hi) or hi <= lo:
        report["status"] = "invalid_range"
        return weights, report
    try:
        hist, bin_edges = np.histogram(vals, bins=np.linspace(lo, hi, n_bins + 1), density=True)
    except Exception as exc:  # noqa: BLE001
        report["status"] = f"histogram_failed: {exc}"
        return weights, report
    hist = hist.astype("float64") * float(n_bins)
    centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    good_hist = np.isfinite(hist)
    if good_hist.sum() < 4:
        report["status"] = "insufficient_histogram_support"
        return weights, report
    try:
        spline = UnivariateSpline(centers[good_hist], hist[good_hist], s=0.0)
        pdf = spline(vals)
    except Exception as exc:  # noqa: BLE001
        report["status"] = f"spline_failed: {exc}"
        return weights, report
    report["pdf_min_before_clip"] = float(np.nanmin(pdf)) if len(pdf) else math.nan
    report["pdf_max_before_clip"] = float(np.nanmax(pdf)) if len(pdf) else math.nan
    report["pdf_negative_fraction_before_clip"] = float(np.mean(pdf < 0.0)) if len(pdf) else math.nan
    pdf = np.clip(pdf, a_min=1.0e-3, a_max=None)
    local = 1.0 / pdf
    if weight_cap is not None:
        local = np.clip(local, a_min=None, a_max=float(weight_cap))
    local = normalize_mean_one(local)
    weights[finite] = local
    report["min_weight"] = float(np.min(local)) if len(local) else math.nan
    report["max_weight"] = float(np.max(local)) if len(local) else math.nan
    report["mean_weight"] = float(np.mean(local)) if len(local) else math.nan
    return weights, report


def compute_ppg12_exact_weights(labels: np.ndarray, et: np.ndarray, eta: np.ndarray):
    labels = np.asarray(labels, dtype="int32")
    weights = np.ones(len(labels), dtype="float64")
    report: dict[str, object] = {
        "weight_mode": "ppg12-exact",
        "event_weight_used": False,
        "vertex_reweight": False,
        "centrality_event_weight": False,
        "cross_section_weight_used_for_training": False,
        "weights_computed_before_binning": True,
        "eta_range": list(PPG12_EXACT_ETA_RANGE),
        "eta_bins": PPG12_EXACT_N_BINS,
        "et_bins": PPG12_EXACT_N_BINS,
        "et_weight_cap": PPG12_EXACT_ET_WEIGHT_CAP,
    }
    class_counts = {str(cls): int((labels == cls).sum()) for cls in (0, 1)}
    if class_counts["0"] <= 0 or class_counts["1"] <= 0:
        raise SystemExit(f"PPG12-exact weights need both classes; observed counts={class_counts}")
    n_total = class_counts["0"] + class_counts["1"]
    class_factors = {}
    eta_reports = {}
    et_reports = {}
    for cls in (0, 1):
        mask = labels == cls
        factor = float(n_total) / (2.0 * float(class_counts[str(cls)]))
        class_factors[str(cls)] = factor
        weights[mask] *= factor
        eta_w, eta_report = inverse_pdf_weights(
            eta[mask],
            n_bins=PPG12_EXACT_N_BINS,
            fixed_range=PPG12_EXACT_ETA_RANGE,
            weight_cap=None,
        )
        weights[mask] *= eta_w
        eta_reports[str(cls)] = eta_report
        et_w, et_report = inverse_pdf_weights(
            et[mask],
            n_bins=PPG12_EXACT_N_BINS,
            fixed_range=None,
            weight_cap=PPG12_EXACT_ET_WEIGHT_CAP,
        )
        weights[mask] *= et_w
        et_reports[str(cls)] = et_report
    finite_positive = np.isfinite(weights) & (weights > 0.0)
    if not finite_positive.all():
        raise SystemExit(f"PPG12-exact weights produced {int((~finite_positive).sum())} bad rows")
    report["class_counts"] = class_counts
    report["class_weight_factors"] = class_factors
    report["eta_reweight"] = eta_reports
    report["et_reweight"] = et_reports
    report["sum_weight_class0"] = float(weights[labels == 0].sum())
    report["sum_weight_class1"] = float(weights[labels == 1].sum())
    report["min_weight"] = float(np.min(weights)) if len(weights) else math.nan
    report["max_weight"] = float(np.max(weights)) if len(weights) else math.nan
    report["mean_weight"] = float(np.mean(weights)) if len(weights) else math.nan
    return weights, report


def weighted_quantile(values: np.ndarray, weights: np.ndarray, q: float) -> float:
    values = np.asarray(values, dtype="float64")
    weights = np.asarray(weights, dtype="float64")
    good = np.isfinite(values) & np.isfinite(weights) & (weights > 0.0)
    if not good.any():
        return math.nan
    values = values[good]
    weights = weights[good]
    order = np.argsort(values, kind="mergesort")
    values = values[order]
    weights = weights[order]
    cdf = np.cumsum(weights)
    total = float(cdf[-1])
    if total <= 0.0:
        return math.nan
    idx = int(np.searchsorted(cdf, max(0.0, min(1.0, q)) * total, side="left"))
    idx = max(0, min(idx, len(values) - 1))
    return float(values[idx])


def effective_entries(weights: np.ndarray) -> float:
    weights = np.asarray(weights, dtype="float64")
    good = np.isfinite(weights) & (weights > 0.0)
    if not good.any():
        return 0.0
    w = weights[good]
    denom = float(np.sum(w * w))
    return float(np.sum(w) ** 2 / denom) if denom > 0.0 else 0.0


def threshold_for_target(sig_score, sig_weight, bkg_score, bkg_weight, target: float) -> dict:
    if sig_score.size == 0 or bkg_score.size == 0:
        return {
            "threshold": math.nan,
            "threshold_stat_err": math.nan,
            "threshold_stat_err_low": math.nan,
            "threshold_stat_err_high": math.nan,
            "signal_efficiency": math.nan,
            "background_fake_rate": math.nan,
            "signal_entries": int(sig_score.size),
            "background_entries": int(bkg_score.size),
            "signal_weight_sum": 0.0,
            "background_weight_sum": 0.0,
            "signal_effective_entries": 0.0,
            "background_effective_entries": 0.0,
        }
    q = max(0.0, min(1.0, 1.0 - target))
    thr = weighted_quantile(sig_score, sig_weight, q)
    sig_sum = float(np.sum(sig_weight))
    bkg_sum = float(np.sum(bkg_weight))
    sig_eff = float(np.sum(sig_weight[sig_score > thr]) / sig_sum) if sig_sum > 0 else math.nan
    bkg_fake = float(np.sum(bkg_weight[bkg_score > thr]) / bkg_sum) if bkg_sum > 0 else math.nan
    n_eff = effective_entries(sig_weight)
    if n_eff > 1.0 and math.isfinite(thr):
        sigma_q = math.sqrt(max(0.0, q * (1.0 - q)) / n_eff)
        thr_low = weighted_quantile(sig_score, sig_weight, max(0.0, q - sigma_q))
        thr_high = weighted_quantile(sig_score, sig_weight, min(1.0, q + sigma_q))
        err_low = max(0.0, thr - thr_low)
        err_high = max(0.0, thr_high - thr)
        err = 0.5 * (err_low + err_high)
    else:
        err = err_low = err_high = math.nan
    return {
        "threshold": float(thr),
        "threshold_stat_err": float(err),
        "threshold_stat_err_low": float(err_low),
        "threshold_stat_err_high": float(err_high),
        "signal_efficiency": sig_eff,
        "background_fake_rate": bkg_fake,
        "signal_entries": int(sig_score.size),
        "background_entries": int(bkg_score.size),
        "signal_weight_sum": sig_sum,
        "background_weight_sum": bkg_sum,
        "signal_effective_entries": n_eff,
        "background_effective_entries": effective_entries(bkg_weight),
    }


def score_matrix(matrix: Path, model: Path, chunk_size: int) -> dict:
    from xgboost import XGBClassifier

    z = np.load(matrix, allow_pickle=True)
    missing = [key for key in [*FEATURES, "is_signal", "cluster_Et", "cluster_Eta", "centrality"] if key not in z.files]
    if missing:
        raise SystemExit(f"{matrix} missing keys: {missing}")
    labels = np.asarray(z["is_signal"], dtype=np.int8)
    et = np.asarray(z["cluster_Et"], dtype=np.float32)
    eta = np.asarray(z["cluster_Eta"], dtype=np.float32)
    cent = np.asarray(z["centrality"], dtype=np.float32)
    weights, weight_report = compute_ppg12_exact_weights(labels, et, eta)

    clf = XGBClassifier()
    clf.load_model(str(model))
    n = len(labels)
    score = np.empty(n, dtype=np.float32)
    for start in range(0, n, chunk_size):
        stop = min(n, start + chunk_size)
        x = np.column_stack([np.asarray(z[name][start:stop], dtype=np.float32) for name in FEATURES])
        score[start:stop] = clf.predict_proba(x)[:, 1].astype(np.float32, copy=False)
    return {
        "labels": labels,
        "et": et,
        "cent": cent,
        "score": score,
        "weights": weights,
        "weight_report": weight_report,
        "n_rows": int(n),
        "matrix_files": list(z.files),
    }


def derive_payload(args: argparse.Namespace) -> dict:
    scored = score_matrix(args.matrix, args.model, args.chunk_size)
    labels = scored["labels"]
    et = scored["et"]
    cent = scored["cent"]
    score = scored["score"]
    weights = scored["weights"]
    base = (
        np.isfinite(score)
        & np.isfinite(et)
        & np.isfinite(cent)
        & np.isin(labels, [0, 1])
        & (et >= ET_RANGE[0])
        & (et < ET_RANGE[1])
        & (cent >= CENT_EDGES[0])
        & (cent < CENT_EDGES[-1])
    )
    rows = []
    et_cells = []
    for target in TARGETS:
        for clo, chi in zip(CENT_EDGES[:-1], CENT_EDGES[1:]):
            cmask = base & (cent >= clo) & (cent < chi)
            sig = cmask & (labels == 1)
            bkg = cmask & (labels == 0)
            item = threshold_for_target(score[sig], weights[sig], score[bkg], weights[bkg], target)
            rows.append(
                {
                    "target_signal_efficiency": target,
                    "wp_label": TARGET_NAMES[target],
                    "centrality_min": clo,
                    "centrality_max": chi,
                    "centrality_center": 0.5 * (clo + chi),
                    "centrality_label": f"{int(clo)}-{int(chi)}%",
                    **item,
                }
            )
            for elo, ehi in zip(ET_EDGES[:-1], ET_EDGES[1:]):
                emask = cmask & (et >= elo) & (et < ehi)
                sig_e = emask & (labels == 1)
                bkg_e = emask & (labels == 0)
                cell = threshold_for_target(score[sig_e], weights[sig_e], score[bkg_e], weights[bkg_e], target)
                et_cells.append(
                    {
                        "target_signal_efficiency": target,
                        "wp_label": TARGET_NAMES[target],
                        "centrality_min": clo,
                        "centrality_max": chi,
                        "centrality_center": 0.5 * (clo + chi),
                        "centrality_label": f"{int(clo)}-{int(chi)}%",
                        "et_min": elo,
                        "et_max": ehi,
                        "et_center": 0.5 * (elo + ehi),
                        **cell,
                    }
                )
    flat_rows = compute_weighted_flat_rows(et_cells)
    return {
        "schema": "THE57_FULL_WEIGHTED_CENTRALITY_DEPENDENT_BDT_WP_V1",
        "matrix": str(args.matrix),
        "model": str(args.model),
        "score_key": "score_centAsFeatBase3x3_pt15to35",
        "features": FEATURES,
        "et_range": list(ET_RANGE),
        "et_edges": ET_EDGES,
        "cent_edges": CENT_EDGES,
        "targets": TARGETS,
        "files_loaded": 1,
        "rows_loaded": int(base.sum()),
        "full_matrix_rows": scored["n_rows"],
        "weighting": scored["weight_report"],
        "threshold_mode": "weighted_signal_quantile",
        "flat_fit_mode": "inverse_variance_weighted_constant_over_et_bins",
        "source_label": "THE-57 full combined weighted training matrix",
        "model_label": "default baseline Au+Au BDT with THE-58 cut",
        "training_sample": "Photon+Jet Embedded 12+20, Inclusive Jet Embedded 12+20+30+40",
        "training_inputs": "baseV3E + weta33/wphi33 + centrality; THE-58 cut; PPG12-exact ET/eta weights",
        "stat_uncertainty_method": "Weighted signal-quantile uncertainty uses signal effective entries per E_T bin; flat constants use inverse-variance weighted constant fits over E_T.",
        "rows": rows,
        "et_cells": et_cells,
        "flat_rows": flat_rows,
    }


def compute_weighted_flat_rows(cells: list[dict]) -> list[dict]:
    rows = []
    for target in TARGETS:
        for clo, chi in zip(CENT_EDGES[:-1], CENT_EDGES[1:]):
            label = f"{int(clo)}-{int(chi)}%"
            sub = [
                c
                for c in cells
                if c["centrality_label"] == label
                and abs(float(c["target_signal_efficiency"]) - target) < 1.0e-9
                and math.isfinite(float(c["threshold"]))
            ]
            if not sub:
                continue
            y = np.array([float(c["threshold"]) for c in sub], dtype="float64")
            err = np.array([float(c.get("threshold_stat_err", math.nan)) for c in sub], dtype="float64")
            good_err = np.isfinite(err) & (err > 0.0)
            if good_err.any():
                fit_w = np.where(good_err, 1.0 / np.maximum(err, 1.0e-6) ** 2, 0.0)
                const = float(np.sum(fit_w * y) / np.sum(fit_w))
            else:
                fit_w = np.ones_like(y)
                const = float(np.mean(y))
            residuals = y - const
            rows.append(
                {
                    "target_signal_efficiency": target,
                    "wp_label": TARGET_NAMES[target],
                    "centrality_min": clo,
                    "centrality_max": chi,
                    "centrality_center": 0.5 * (clo + chi),
                    "centrality_label": label,
                    "threshold": const,
                    "signal_efficiency": target,
                    "background_fake_rate": math.nan,
                    "signal_entries": int(sum(int(c["signal_entries"]) for c in sub)),
                    "background_entries": int(sum(int(c["background_entries"]) for c in sub)),
                    "signal_weight_sum": float(sum(float(c.get("signal_weight_sum", 0.0)) for c in sub)),
                    "background_weight_sum": float(sum(float(c.get("background_weight_sum", 0.0)) for c in sub)),
                    "signal_effective_entries": float(sum(float(c.get("signal_effective_entries", 0.0)) for c in sub)),
                    "background_effective_entries": float(sum(float(c.get("background_effective_entries", 0.0)) for c in sub)),
                    "n_et_points": int(len(sub)),
                    "threshold_stat_err": float(math.sqrt(1.0 / np.sum(fit_w))) if np.sum(fit_w) > 0.0 else math.nan,
                    "threshold_stat_err_low": math.nan,
                    "threshold_stat_err_high": math.nan,
                    "flat_max_abs_residual": float(np.max(np.abs(residuals))),
                    "flat_rms_residual": float(np.sqrt(np.mean(residuals**2))),
                }
            )
    return rows


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--matrix", type=Path, required=True)
    ap.add_argument("--model", type=Path, required=True)
    ap.add_argument("--json-out", type=Path, required=True)
    ap.add_argument("--chunk-size", type=int, default=500_000)
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    payload = derive_payload(args)
    if str(args.json_out) == "-":
        print("__THE57_FULL_WP_JSON_BEGIN__")
        print(json.dumps(payload, separators=(",", ":")))
        print("__THE57_FULL_WP_JSON_END__")
    else:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        print(args.json_out)


if __name__ == "__main__":
    main()
