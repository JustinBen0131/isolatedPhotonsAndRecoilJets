#!/usr/bin/env python3
"""Extract THE-37 fixed-inclusive signal train/validation matrix payload.

This script is intended to be streamed to SDCC and run against existing model
registries/matrices. It performs read-only scoring and writes JSON to stdout.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import math
from collections import Counter
from pathlib import Path

import numpy as np
import xgboost as xgb
from scipy.interpolate import UnivariateSpline
from sklearn.metrics import brier_score_loss, log_loss, roc_auc_score
from sklearn.model_selection import train_test_split


DEFAULT_BASELINE_DIR = Path(
    "/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/the57_models/a1_withcut_20260615_stagedcache_streamreduce"
)
DEFAULT_PHOTON12_DIR = Path(
    "/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/the8_models/"
    "THE8_corrected_truthdefault_20260615_sig12_bkg12_20_30_40"
)

BINS = np.linspace(0.0, 1.0, 51)
ETA_RANGE = (-0.7, 0.7)
N_BINS = 20
ET_CAP = 800.0


def safe_ratio(numer: np.ndarray, denom: np.ndarray) -> np.ndarray:
    numer = np.asarray(numer, dtype="float64")
    denom = np.asarray(denom, dtype="float64")
    out = np.full(len(numer), np.nan, dtype="float64")
    good = np.isfinite(numer) & np.isfinite(denom) & (np.abs(denom) > 1.0e-9)
    out[good] = numer[good] / denom[good]
    return out


def load_registry(path: Path) -> dict:
    data = json.loads(path.read_text())
    model = data["models"][0]
    report = model["report"]
    return {
        "model": model,
        "features": list(model["features"]),
        "pt_range": model.get("pt_range") or report.get("pt_range"),
        "cent_range": model.get("cent_range") or report.get("cent_range"),
        "split": report["split"],
        "holdout_rows": int(report["overfit_diagnostics"]["holdout_rows"]),
        "report": report,
    }


def load_frame(path: Path, features: list[str]) -> dict[str, np.ndarray]:
    matrix = np.load(path, allow_pickle=True)
    needed = set(features) | {"is_signal", "centrality", "cluster_Et", "cluster_Eta", "source_sample"}
    base_needed: set[str] = set()
    for name in needed:
        if name == "cluster_weta_over_wphi":
            base_needed.update(["cluster_weta_cogx", "cluster_wphi_cogx"])
        elif name == "cluster_weta33_over_wphi33":
            base_needed.update(["cluster_weta33_cogx", "cluster_wphi33_cogx"])
        else:
            base_needed.add(name)
    missing = sorted(col for col in base_needed if col not in matrix.files)
    if missing:
        raise SystemExit(f"{path} is missing columns: {missing}")
    frame = {name: matrix[name] for name in base_needed}
    if "cluster_weta_over_wphi" in features:
        frame["cluster_weta_over_wphi"] = safe_ratio(frame["cluster_weta_cogx"], frame["cluster_wphi_cogx"])
    if "cluster_weta33_over_wphi33" in features:
        frame["cluster_weta33_over_wphi33"] = safe_ratio(frame["cluster_weta33_cogx"], frame["cluster_wphi33_cogx"])
    return frame


def filtered_indices(frame: dict[str, np.ndarray], features: list[str], pt_range, cent_range) -> np.ndarray:
    mask = np.ones(len(frame["is_signal"]), dtype=bool)
    if pt_range is not None:
        lo, hi = float(pt_range[0]), float(pt_range[1])
        et = frame["cluster_Et"]
        mask &= np.isfinite(et) & (et >= lo) & (et < hi)
    if cent_range is not None and cent_range != "all":
        lo, hi = float(cent_range[0]), float(cent_range[1])
        cent = frame["centrality"]
        mask &= np.isfinite(cent) & (cent >= lo) & (cent < hi)
    y = frame["is_signal"].astype("int32", copy=False)
    mask &= np.isin(y, [0, 1])
    for name in features:
        mask &= np.isfinite(frame[name])
    return np.flatnonzero(mask)


def holdout_indices(frame: dict[str, np.ndarray], registry: dict, random_seed: int) -> np.ndarray:
    idx = filtered_indices(frame, registry["features"], registry["pt_range"], registry["cent_range"])
    y = frame["is_signal"][idx].astype("int32", copy=False)
    _, test_idx = train_test_split(idx, test_size=0.10, random_state=random_seed, stratify=y)
    test_idx = np.asarray(test_idx, dtype="int64")
    if len(test_idx) != registry["holdout_rows"]:
        raise SystemExit(f"holdout mismatch: got {len(test_idx)}, expected {registry['holdout_rows']}")
    return test_idx


def normalize_mean_one(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype="float64")
    finite = np.isfinite(values) & (values > 0.0)
    if not finite.any():
        return np.ones(len(values), dtype="float64")
    mean = float(np.mean(values[finite]))
    if not math.isfinite(mean) or mean <= 0.0:
        return np.ones(len(values), dtype="float64")
    out = np.ones(len(values), dtype="float64")
    out[finite] = values[finite] / mean
    return out


def inverse_pdf_weights(values: np.ndarray, *, fixed_range=None, weight_cap: float | None = None) -> tuple[np.ndarray, dict]:
    values = np.asarray(values, dtype="float64")
    weights = np.ones(len(values), dtype="float64")
    finite = np.isfinite(values)
    report = {
        "n_entries": int(len(values)),
        "n_finite": int(finite.sum()),
        "n_bins": int(N_BINS),
        "fixed_range": list(fixed_range) if fixed_range is not None else None,
        "weight_cap": float(weight_cap) if weight_cap is not None else None,
        "status": "ok",
    }
    if finite.sum() < max(10, N_BINS):
        report["status"] = "insufficient_finite_values"
        return weights, report

    vals = values[finite]
    if fixed_range is None:
        lo, hi = float(np.min(vals)), float(np.max(vals))
    else:
        lo, hi = float(fixed_range[0]), float(fixed_range[1])
    report["range"] = [lo, hi]
    if not math.isfinite(lo) or not math.isfinite(hi) or hi <= lo:
        report["status"] = "invalid_range"
        return weights, report

    edges = np.linspace(lo, hi, N_BINS + 1)
    hist, bin_edges = np.histogram(vals, bins=edges, density=True)
    hist = hist.astype("float64") * float(N_BINS)
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


def ppg12_weights_for_display(
    frame: dict[str, np.ndarray],
    source_keep_idx: np.ndarray,
    display_y: np.ndarray,
) -> tuple[np.ndarray, dict]:
    labels = display_y.astype("int32", copy=False)
    weights = np.ones(len(labels), dtype="float64")
    report: dict[str, object] = {
        "weight_mode": "ppg12-exact-display-source-classes",
        "event_weight_used": False,
        "cross_section_weight_used_for_training": False,
        "weights_computed_before_binning": True,
        "eta_range": list(ETA_RANGE),
        "eta_bins": N_BINS,
        "et_bins": N_BINS,
        "et_weight_cap": ET_CAP,
    }
    class_counts = {str(cls): int((labels == cls).sum()) for cls in (0, 1)}
    n_total = sum(class_counts.values())
    if class_counts["0"] <= 0 or class_counts["1"] <= 0:
        raise SystemExit(f"need both display classes, got {class_counts}")
    class_factors = {}
    for cls in (0, 1):
        factor = n_total / (2.0 * class_counts[str(cls)])
        class_factors[str(cls)] = float(factor)
        weights[labels == cls] *= factor

    eta_reports = {}
    et_reports = {}
    eta = frame["cluster_Eta"][source_keep_idx]
    et = frame["cluster_Et"][source_keep_idx]
    for cls in (0, 1):
        mask = labels == cls
        eta_w, eta_report = inverse_pdf_weights(eta[mask], fixed_range=ETA_RANGE, weight_cap=None)
        et_w, et_report = inverse_pdf_weights(et[mask], fixed_range=None, weight_cap=ET_CAP)
        weights[mask] *= eta_w
        weights[mask] *= et_w
        eta_reports[str(cls)] = eta_report
        et_reports[str(cls)] = et_report

    finite_positive = np.isfinite(weights) & (weights > 0.0)
    if not finite_positive.all():
        raise SystemExit(f"display weights produced {int((~finite_positive).sum())} invalid rows")
    report.update(
        {
            "class_counts": class_counts,
            "class_weight_factors": class_factors,
            "eta_reweight": eta_reports,
            "et_reweight": et_reports,
            "sum_weight_class0": float(weights[labels == 0].sum()),
            "sum_weight_class1": float(weights[labels == 1].sum()),
            "min_weight": float(np.min(weights)),
            "max_weight": float(np.max(weights)),
            "mean_weight": float(np.mean(weights)),
        }
    )
    return weights, report


def score_model(
    model_path: Path,
    frame: dict[str, np.ndarray],
    features: list[str],
    idx: np.ndarray,
    batch_size: int,
) -> np.ndarray:
    booster = xgb.Booster()
    booster.load_model(str(model_path))
    out = np.empty(len(idx), dtype="float32")
    for start in range(0, len(idx), batch_size):
        stop = min(start + batch_size, len(idx))
        local = idx[start:stop]
        x = np.column_stack([frame[name][local].astype("float32", copy=False) for name in features])
        out[start:stop] = booster.predict(xgb.DMatrix(x)).astype("float32")
    return out


def weighted_hist(values: np.ndarray, weights: np.ndarray) -> dict:
    good = np.isfinite(values) & np.isfinite(weights) & (weights > 0.0)
    counts, _ = np.histogram(values[good], bins=BINS)
    weighted_counts, _ = np.histogram(values[good], bins=BINS, weights=weights[good])
    density = weighted_counts.astype("float64")
    if density.sum() > 0:
        density = density / density.sum() / np.diff(BINS)
    return {
        "counts": counts.astype(int).tolist(),
        "weighted_counts": weighted_counts.astype(float).tolist(),
        "density": density.astype(float).tolist(),
        "entries": int(good.sum()),
        "sum_weight": float(weighted_counts.sum()),
    }


def weighted_quantile(values: np.ndarray, weights: np.ndarray, q: float) -> float:
    good = np.isfinite(values) & np.isfinite(weights) & (weights > 0.0)
    if not good.any():
        return math.nan
    v = values[good].astype("float64")
    w = weights[good].astype("float64")
    order = np.argsort(v, kind="mergesort")
    v = v[order]
    w = w[order]
    cumsum = np.cumsum(w)
    target = q * cumsum[-1]
    return float(v[min(np.searchsorted(cumsum, target, side="left"), len(v) - 1)])


def cell_record(
    train_key: str,
    train_label: str,
    validation_key: str,
    validation_label: str,
    lane: str,
    model_path: Path,
    score: np.ndarray,
    signal_mask: np.ndarray,
    inclusive_mask: np.ndarray,
    weights: np.ndarray,
) -> dict:
    sig_score = score[signal_mask]
    inc_score = score[inclusive_mask]
    sig_w = weights[signal_mask]
    inc_w = weights[inclusive_mask]
    y = np.concatenate([np.ones(len(sig_score), dtype="int32"), np.zeros(len(inc_score), dtype="int32")])
    s = np.concatenate([sig_score, inc_score]).astype("float64")
    w = np.concatenate([sig_w, inc_w]).astype("float64")
    good = np.isfinite(s) & np.isfinite(w) & (w > 0.0)
    auc = float(roc_auc_score(y[good], s[good], sample_weight=w[good])) if len(np.unique(y[good])) == 2 else math.nan
    logloss = float(log_loss(y[good], np.clip(s[good], 1e-6, 1.0 - 1e-6), sample_weight=w[good], labels=[0, 1]))
    brier = float(brier_score_loss(y[good], s[good], sample_weight=w[good]))
    med_sig = weighted_quantile(sig_score, sig_w, 0.5)
    med_inc = weighted_quantile(inc_score, inc_w, 0.5)
    wp80 = weighted_quantile(sig_score, sig_w, 0.2)
    finite_inc = np.isfinite(inc_score)
    finite_sig = np.isfinite(sig_score)
    wp80_inc_pass = float(np.sum(inc_w[finite_inc & (inc_score >= wp80)]) / np.sum(inc_w[finite_inc]))
    wp80_sig_eff = float(np.sum(sig_w[finite_sig & (sig_score >= wp80)]) / np.sum(sig_w[finite_sig]))
    return {
        "row_key": train_key,
        "row_label": train_label,
        "col_key": validation_key,
        "col_label": validation_label,
        "lane": lane,
        "model_path": str(model_path),
        "bin_edges": BINS.astype(float).tolist(),
        "signal": weighted_hist(sig_score, sig_w),
        "background": weighted_hist(inc_score, inc_w),
        "entries": {
            "signal": int(signal_mask.sum()),
            "background": int(inclusive_mask.sum()),
            "total": int(signal_mask.sum() + inclusive_mask.sum()),
        },
        "sum_weight": {
            "signal": float(np.sum(sig_w)),
            "background": float(np.sum(inc_w)),
            "total": float(np.sum(sig_w) + np.sum(inc_w)),
        },
        "metrics": {
            "auc": auc,
            "logloss": logloss,
            "brier": brier,
            "median_signal_score": med_sig,
            "median_background_score": med_inc,
            "median_gap": float(med_sig - med_inc),
            "wp80_threshold": wp80,
            "wp80_background_fake_rate": wp80_inc_pass,
            "wp80_signal_efficiency": wp80_sig_eff,
        },
    }


def build_payload(args: argparse.Namespace) -> dict:
    baseline_reg = load_registry(args.baseline_dir / "model_registry.json")
    photon12_reg = load_registry(args.photon12_dir / "model_registry.json")
    if baseline_reg["features"] != photon12_reg["features"]:
        raise SystemExit("feature mismatch between baseline and Photon12 model")

    frame = load_frame(args.baseline_dir / "training_matrix.npz", baseline_reg["features"])
    holdout = holdout_indices(frame, baseline_reg, args.random_seed)
    source = frame["source_sample"].astype(str)
    original_y = frame["is_signal"].astype("int32", copy=False)
    cent = frame["centrality"]
    cent0_20 = np.isfinite(cent) & (cent >= 0.0) & (cent < 20.0)
    idx = holdout[cent0_20[holdout]]

    signal12_global = (np.char.find(source, "embeddedPhoton12") >= 0) & (original_y == 1)
    signal1220_global = (np.char.find(source, "embeddedPhoton") >= 0) & (original_y == 1)
    inclusive_global = np.char.find(source, "embeddedJet") >= 0

    filtered = filtered_indices(frame, baseline_reg["features"], baseline_reg["pt_range"], baseline_reg["cent_range"])
    filtered_mask = np.zeros(len(source), dtype=bool)
    filtered_mask[filtered] = True

    validation_defs = [
        ("photon12", "Photon12 only", signal12_global),
        ("photon12_20", "Photon12+20", signal1220_global),
    ]
    weights_by_validation = {}
    weight_reports = {}
    for key, _label, signal_global in validation_defs:
        source_keep_idx = np.flatnonzero(filtered_mask & (signal_global | inclusive_global))
        display_y = np.zeros(len(source_keep_idx), dtype="int32")
        display_y[signal_global[source_keep_idx]] = 1
        weights, report = ppg12_weights_for_display(frame, source_keep_idx, display_y)
        source_position = {int(v): i for i, v in enumerate(source_keep_idx.tolist())}
        local_weights = np.full(len(idx), np.nan, dtype="float64")
        for local_i, global_i in enumerate(idx.tolist()):
            pos = source_position.get(int(global_i))
            if pos is not None:
                local_weights[local_i] = weights[pos]
        weights_by_validation[key] = local_weights
        weight_reports[key] = report

    models = [
        (
            "photon12",
            "Photon12 only",
            "THE8_corrected_truthdefault_20260615_sig12_bkg12_20_30_40",
            args.photon12_dir / "auau_tight_bdt_centAsFeatBase3x3_pt15to35_tmva.xgb.json",
        ),
        (
            "photon12_20",
            "Photon12+20",
            "THE57_A1_corrected_baseline_full_signal_full_inclusive",
            args.baseline_dir / "auau_tight_bdt_centAsFeatBase3x3_pt15to35_tmva.xgb.json",
        ),
    ]
    scores = {
        key: score_model(model_path, frame, baseline_reg["features"], idx, args.batch_size)
        for key, _label, _lane, model_path in models
    }

    cells = []
    for train_key, train_label, lane, model_path in models:
        score = scores[train_key]
        for validation_key, validation_label, signal_global in validation_defs:
            weights = weights_by_validation[validation_key]
            good_weight = np.isfinite(weights) & (weights > 0.0)
            cells.append(
                cell_record(
                    train_key,
                    train_label,
                    validation_key,
                    validation_label,
                    lane,
                    model_path,
                    score,
                    signal_global[idx] & good_weight,
                    inclusive_global[idx] & good_weight,
                    weights,
                )
            )

    source_counts = Counter(source[idx]).most_common(20)
    return {
        "schema": "THE37_AUAU_FIXED_INCLUSIVE_SIGNAL_TRAINVAL_SCORE_SHAPE_V1",
        "created_remote_time": dt.datetime.now(dt.timezone.utc).isoformat(),
        "auau_signal_trainval": {
            "title": "Au+Au 0-20%: fixed embedded inclusive MC, signal sample train/validation matrix",
            "row_label_prefix": "TRAIN SIGNAL",
            "col_label_prefix": "VALIDATION SIGNAL",
            "row_order": [
                {"key": "photon12", "label": "Photon12 only"},
                {"key": "photon12_20", "label": "Photon12+20"},
            ],
            "col_order": [
                {"key": "photon12", "label": "Photon12 only"},
                {"key": "photon12_20", "label": "Photon12+20"},
            ],
            "fixed_inclusive_label": "FIXED INCLUSIVE MC: Jet12+20+30+40 in training and validation",
            "validation_meta": {
                "class_definition": (
                    "Signal columns use embeddedPhoton12 truth-isolated prompt or embeddedPhoton12+20 truth-isolated prompt; "
                    "inclusive MC = embeddedJet12/20/30/40 with no truth filter."
                ),
                "fixed_inclusive_mc": (
                    "run28_embeddedJet12 + run28_embeddedJet20 + run28_embeddedJet30 + run28_embeddedJet40 "
                    "in both training and validation inclusive class."
                ),
                "validation_matrix": str(args.baseline_dir / "training_matrix.npz"),
                "validation_holdout": (
                    f"baseline 10% row holdout reconstructed with random_state={args.random_seed}, "
                    "then centrality 0-20% selected"
                ),
                "fixed_holdout_cent0_20_rows": int(len(idx)),
                "source_counts_holdout_cent0_20": [[str(k), int(v)] for k, v in source_counts],
                "background_holdout_entries": int(inclusive_global[idx].sum()),
                "signal12_holdout_entries": int(signal12_global[idx].sum()),
                "signal12_20_holdout_entries": int(signal1220_global[idx].sum()),
                "weight_reports": weight_reports,
                "train_models": {
                    "photon12": {
                        "lane": "THE8_corrected_truthdefault_20260615_sig12_bkg12_20_30_40",
                        "model_path": str(args.photon12_dir / "auau_tight_bdt_centAsFeatBase3x3_pt15to35_tmva.xgb.json"),
                        "training_signal": "run28_embeddedPhoton12",
                        "training_inclusive_mc": "run28_embeddedJet12/20/30/40",
                    },
                    "photon12_20": {
                        "lane": "THE57_A1_corrected_baseline_full_signal_full_inclusive",
                        "model_path": str(args.baseline_dir / "auau_tight_bdt_centAsFeatBase3x3_pt15to35_tmva.xgb.json"),
                        "training_signal": "run28_embeddedPhoton12 + run28_embeddedPhoton20",
                        "training_inclusive_mc": "run28_embeddedJet12/20/30/40",
                    },
                },
            },
            "cells": cells,
            "ymax": float(max(max(c["signal"]["density"] + c["background"]["density"]) for c in cells) * 1.08),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-dir", type=Path, default=DEFAULT_BASELINE_DIR)
    parser.add_argument("--photon12-dir", type=Path, default=DEFAULT_PHOTON12_DIR)
    parser.add_argument("--random-seed", type=int, default=13)
    parser.add_argument("--batch-size", type=int, default=200000)
    args = parser.parse_args()
    print(json.dumps(build_payload(args), indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
