#!/usr/bin/env python3
"""Extract same-holdout THE8 routed-BDT metric payload.

This script is designed to run read-only on SDCC against existing corrected
THE8 registries. It writes a compact JSON payload to stdout; no remote files are
created by default.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import math
from pathlib import Path

import numpy as np
import xgboost as xgb
from scipy.interpolate import UnivariateSpline
from sklearn.metrics import log_loss, roc_auc_score
from sklearn.model_selection import train_test_split


DEFAULT_BASELINE_DIR = Path(
    "/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/the57_models/"
    "a1_withcut_20260615_stagedcache_streamreduce"
)
DEFAULT_BINNED_DIR = Path(
    "/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/the8_models/"
    "THE8_corrected_truthdefault_20260615_binned14"
)
BASELINE_MODEL_ID = "centAsFeatBase3x3_pt15to35"
PRODUCTS = ["global", "cent3", "cent7", "et", "etcent3", "etcent7"]
PRODUCT_TO_REGISTRY = {
    "cent3": "base14_perCent3",
    "cent7": "base14_perCent7",
    "et": "base14_perEt",
    "etcent3": "base14_perEtCent3",
    "etcent7": "base14_perEtCent7",
}
ETA_RANGE = (-0.7, 0.7)
N_WEIGHT_BINS = 20
ET_WEIGHT_CAP = 800.0
COARSE_CENT_BINS = [(0.0, 20.0), (20.0, 50.0), (50.0, 80.0)]
FINE_CENT_BINS = [(0.0, 10.0), (10.0, 20.0), (20.0, 30.0), (30.0, 40.0), (40.0, 50.0), (50.0, 60.0), (60.0, 80.0)]
ET_BINS = [(15.0, 17.0), (17.0, 19.0), (19.0, 21.0), (21.0, 23.0), (23.0, 25.0), (25.0, 27.0), (27.0, 30.0), (30.0, 35.0)]


def range_label(lo: float, hi: float, suffix: str) -> str:
    return f"{lo:g}-{hi:g}{suffix}"


def safe_ratio(numer: np.ndarray, denom: np.ndarray) -> np.ndarray:
    numer = np.asarray(numer, dtype="float64")
    denom = np.asarray(denom, dtype="float64")
    out = np.full(len(numer), np.nan, dtype="float64")
    good = np.isfinite(numer) & np.isfinite(denom) & (np.abs(denom) > 1.0e-9)
    out[good] = numer[good] / denom[good]
    return out


def load_registry(path: Path) -> dict:
    if not path.is_file():
        raise SystemExit(f"Missing registry: {path}")
    return json.loads(path.read_text())


def model_by_id(registry: dict, model_id: str) -> dict:
    for model in registry.get("models", []):
        if model.get("model_id") == model_id:
            return model
    raise SystemExit(f"Registry has no model_id={model_id}")


def models_by_product(registry: dict, product: str) -> list[dict]:
    models = [m for m in registry.get("models", []) if m.get("product") == product]
    if not models:
        raise SystemExit(f"Registry has no product={product}")
    return sorted(models, key=lambda m: (m.get("pt_range") or [-1, -1], m.get("cent_range") or [-1, -1], m["model_id"]))


def expand_needed_columns(features: list[str]) -> set[str]:
    needed = set(features) | {"is_signal", "centrality", "cluster_Et", "cluster_Eta", "source_sample"}
    expanded: set[str] = set()
    for name in needed:
        if name == "cluster_weta_over_wphi":
            expanded.update(["cluster_weta_cogx", "cluster_wphi_cogx"])
        elif name == "cluster_weta33_over_wphi33":
            expanded.update(["cluster_weta33_cogx", "cluster_wphi33_cogx"])
        else:
            expanded.add(name)
    return expanded


def load_frame(path: Path, features: list[str]) -> dict[str, np.ndarray]:
    matrix = np.load(path, allow_pickle=True)
    needed = expand_needed_columns(features)
    missing = sorted(col for col in needed if col not in matrix.files)
    if missing:
        raise SystemExit(f"{path} is missing columns: {missing}")
    frame = {name: matrix[name] for name in needed}
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
    if cent_range is not None:
        lo, hi = float(cent_range[0]), float(cent_range[1])
        cent = frame["centrality"]
        mask &= np.isfinite(cent) & (cent >= lo) & (cent < hi)
    y = frame["is_signal"].astype("int32", copy=False)
    mask &= np.isin(y, [0, 1])
    for name in features:
        mask &= np.isfinite(frame[name])
    return np.flatnonzero(mask)


def reconstruct_holdout(frame: dict[str, np.ndarray], model: dict, random_seed: int) -> np.ndarray:
    features = list(model["features"])
    report = model.get("report", {})
    pt_range = model.get("pt_range") or report.get("pt_range")
    cent_range = model.get("cent_range") or report.get("cent_range")
    idx = filtered_indices(frame, features, pt_range, cent_range)
    y = frame["is_signal"][idx].astype("int32", copy=False)
    _, test_idx = train_test_split(idx, test_size=0.10, random_state=random_seed, stratify=y)
    expected = report.get("overfit_diagnostics", {}).get("holdout_rows")
    if expected is not None and int(expected) != len(test_idx):
        raise SystemExit(f"Holdout mismatch: got {len(test_idx)}, expected {expected}")
    return np.asarray(test_idx, dtype="int64")


def range_mask(values: np.ndarray, value_range) -> np.ndarray:
    if value_range is None:
        return np.ones(len(values), dtype=bool)
    lo, hi = float(value_range[0]), float(value_range[1])
    return np.isfinite(values) & (values >= lo) & (values < hi)


def score_model(model_path: Path, frame: dict[str, np.ndarray], features: list[str], idx: np.ndarray, batch_size: int) -> np.ndarray:
    booster = xgb.Booster()
    booster.load_model(str(model_path))
    out = np.empty(len(idx), dtype="float32")
    for start in range(0, len(idx), batch_size):
        stop = min(start + batch_size, len(idx))
        local = idx[start:stop]
        x = np.column_stack([frame[name][local].astype("float32", copy=False) for name in features])
        out[start:stop] = booster.predict(xgb.DMatrix(x)).astype("float32")
    return out


def score_routed_product(
    product: str,
    models: list[dict],
    frame: dict[str, np.ndarray],
    idx: np.ndarray,
    batch_size: int,
) -> tuple[np.ndarray, dict]:
    local_et = frame["cluster_Et"][idx]
    local_cent = frame["centrality"][idx]
    score = np.full(len(idx), np.nan, dtype="float32")
    routes = []
    for model in models:
        report = model.get("report", {})
        pt_range = model.get("pt_range") or report.get("pt_range")
        cent_range = model.get("cent_range") or report.get("cent_range")
        mask = range_mask(local_et, pt_range) & range_mask(local_cent, cent_range)
        if not np.any(mask):
            continue
        route_idx = idx[np.flatnonzero(mask)]
        route_score = score_model(Path(model["output_xgb_json"]), frame, list(model["features"]), route_idx, batch_size)
        score[mask] = route_score
        routes.append(
            {
                "model_id": model["model_id"],
                "product": product,
                "pt_range": pt_range,
                "cent_range": cent_range,
                "rows_scored": int(mask.sum()),
                "features": list(model["features"]),
                "output_xgb_json": model["output_xgb_json"],
            }
        )
    missing = int(np.sum(~np.isfinite(score)))
    if missing:
        raise SystemExit(f"{product} left {missing} rows unscored")
    return score, {"routes": routes, "route_count": len(routes)}


def class_masks(frame: dict[str, np.ndarray], idx: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    source = frame["source_sample"][idx].astype(str)
    y = frame["is_signal"][idx].astype("int32", copy=False)
    signal = (np.char.find(source, "embeddedPhoton") >= 0) & (y == 1)
    inclusive = np.char.find(source, "embeddedJet") >= 0
    overlap = signal & inclusive
    if np.any(overlap):
        raise SystemExit("Signal and inclusive masks overlap")
    return signal, inclusive


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
        "n_bins": int(N_WEIGHT_BINS),
        "fixed_range": list(fixed_range) if fixed_range is not None else None,
        "weight_cap": float(weight_cap) if weight_cap is not None else None,
        "status": "ok",
    }
    if finite.sum() < max(10, N_WEIGHT_BINS):
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
    edges = np.linspace(lo, hi, N_WEIGHT_BINS + 1)
    hist, bin_edges = np.histogram(vals, bins=edges, density=True)
    hist = hist.astype("float64") * float(N_WEIGHT_BINS)
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


def display_weights(frame: dict[str, np.ndarray], source_keep_idx: np.ndarray, display_y: np.ndarray) -> tuple[np.ndarray, dict]:
    labels = display_y.astype("int32", copy=False)
    weights = np.ones(len(labels), dtype="float64")
    report: dict[str, object] = {
        "weight_mode": "ppg12-exact-display-source-classes",
        "event_weight_used": False,
        "cross_section_weight_used_for_training": False,
        "weights_computed_before_binning": True,
        "eta_range": list(ETA_RANGE),
        "eta_bins": N_WEIGHT_BINS,
        "et_bins": N_WEIGHT_BINS,
        "et_weight_cap": ET_WEIGHT_CAP,
    }
    class_counts = {str(cls): int((labels == cls).sum()) for cls in (0, 1)}
    n_total = sum(class_counts.values())
    if class_counts["0"] <= 0 or class_counts["1"] <= 0:
        raise SystemExit(f"Need both display classes, got {class_counts}")
    class_factors = {}
    for cls in (0, 1):
        factor = n_total / (2.0 * class_counts[str(cls)])
        class_factors[str(cls)] = float(factor)
        weights[labels == cls] *= factor
    eta = frame["cluster_Eta"][source_keep_idx]
    et = frame["cluster_Et"][source_keep_idx]
    eta_reports = {}
    et_reports = {}
    for cls in (0, 1):
        mask = labels == cls
        eta_w, eta_report = inverse_pdf_weights(eta[mask], fixed_range=ETA_RANGE, weight_cap=None)
        et_w, et_report = inverse_pdf_weights(et[mask], fixed_range=None, weight_cap=ET_WEIGHT_CAP)
        weights[mask] *= eta_w
        weights[mask] *= et_w
        eta_reports[str(cls)] = eta_report
        et_reports[str(cls)] = et_report
    finite_positive = np.isfinite(weights) & (weights > 0.0)
    if not finite_positive.all():
        raise SystemExit(f"Display weights produced {int((~finite_positive).sum())} invalid rows")
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


def quantile(values: np.ndarray, q: float) -> float:
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return math.nan
    return float(np.quantile(values.astype("float64"), q))


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


def metrics_for(score: np.ndarray, signal: np.ndarray, inclusive: np.ndarray, weights: np.ndarray) -> dict:
    sig_score = score[signal & np.isfinite(score)]
    inc_score = score[inclusive & np.isfinite(score)]
    sig_w = weights[signal & np.isfinite(score)]
    inc_w = weights[inclusive & np.isfinite(score)]
    if len(sig_score) == 0 or len(inc_score) == 0:
        return {
            "auc": math.nan,
            "logloss": math.nan,
            "median_gap": math.nan,
            "wp80_inclusive_fake": math.nan,
            "wp80_threshold": math.nan,
            "signal_entries": int(len(sig_score)),
            "inclusive_entries": int(len(inc_score)),
        }
    y = np.concatenate([np.ones(len(sig_score), dtype="int32"), np.zeros(len(inc_score), dtype="int32")])
    s = np.concatenate([sig_score, inc_score]).astype("float64")
    w = np.concatenate([sig_w, inc_w]).astype("float64")
    good = np.isfinite(s) & np.isfinite(w) & (w > 0.0)
    auc = float(roc_auc_score(y[good], s[good], sample_weight=w[good])) if len(np.unique(y[good])) == 2 else math.nan
    ll = float(log_loss(y[good], np.clip(s[good], 1e-6, 1.0 - 1e-6), sample_weight=w[good], labels=[0, 1]))
    wp80 = weighted_quantile(sig_score, sig_w, 0.20)
    finite_inc = np.isfinite(inc_score) & np.isfinite(inc_w) & (inc_w > 0.0)
    fake = float(np.sum(inc_w[finite_inc & (inc_score >= wp80)]) / np.sum(inc_w[finite_inc])) if math.isfinite(wp80) else math.nan
    med_sig = weighted_quantile(sig_score, sig_w, 0.50)
    med_inc = weighted_quantile(inc_score, inc_w, 0.50)
    return {
        "auc": auc,
        "logloss": ll,
        "median_signal_score": med_sig,
        "median_inclusive_score": med_inc,
        "median_gap": med_sig - med_inc,
        "wp80_inclusive_fake": fake,
        "wp80_threshold": wp80,
        "signal_entries": int(len(sig_score)),
        "inclusive_entries": int(len(inc_score)),
        "sum_weight_signal": float(np.sum(sig_w[np.isfinite(sig_w)])),
        "sum_weight_inclusive": float(np.sum(inc_w[np.isfinite(inc_w)])),
    }


def group_mask(kind: str, item: dict, et: np.ndarray, cent: np.ndarray) -> np.ndarray:
    mask = np.ones(len(et), dtype=bool)
    if "et_range" in item:
        mask &= range_mask(et, item["et_range"])
    if "cent_range" in item:
        mask &= range_mask(cent, item["cent_range"])
    if kind == "all":
        return mask
    return mask


def build_groups() -> list[dict]:
    groups = [{"kind": "all", "key": "all", "label": "15-35 GeV, 0-80%", "sort": [0, 0]}]
    for lo, hi in COARSE_CENT_BINS:
        groups.append({"kind": "coarse_cent", "key": f"cent_{lo:g}_{hi:g}", "label": range_label(lo, hi, "%"), "cent_range": [lo, hi], "sort": [1, lo]})
    for lo, hi in FINE_CENT_BINS:
        groups.append({"kind": "fine_cent", "key": f"cent_{lo:g}_{hi:g}", "label": range_label(lo, hi, "%"), "cent_range": [lo, hi], "sort": [2, lo]})
    for lo, hi in ET_BINS:
        groups.append({"kind": "et", "key": f"et_{lo:g}_{hi:g}", "label": range_label(lo, hi, " GeV"), "et_range": [lo, hi], "sort": [3, lo]})
    for elo, ehi in ET_BINS:
        for clo, chi in COARSE_CENT_BINS:
            groups.append(
                {
                    "kind": "et_coarse_cent",
                    "key": f"et_{elo:g}_{ehi:g}_cent_{clo:g}_{chi:g}",
                    "label": f"{range_label(elo, ehi, ' GeV')} × {range_label(clo, chi, '%')}",
                    "et_range": [elo, ehi],
                    "cent_range": [clo, chi],
                    "sort": [4, elo, clo],
                }
            )
    for elo, ehi in ET_BINS:
        for clo, chi in FINE_CENT_BINS:
            groups.append(
                {
                    "kind": "et_fine_cent",
                    "key": f"et_{elo:g}_{ehi:g}_cent_{clo:g}_{chi:g}",
                    "label": f"{range_label(elo, ehi, ' GeV')} × {range_label(clo, chi, '%')}",
                    "et_range": [elo, ehi],
                    "cent_range": [clo, chi],
                    "sort": [5, elo, clo],
                }
            )
    return groups


def build_payload(args: argparse.Namespace) -> dict:
    baseline_registry = load_registry(args.baseline_dir / "model_registry.json")
    binned_registry = load_registry(args.binned_dir / "model_registry.json")
    baseline_model = model_by_id(baseline_registry, BASELINE_MODEL_ID)
    features = list(baseline_model["features"])
    for model in binned_registry.get("models", []):
        for feat in model.get("features", []):
            if feat not in features:
                features.append(feat)
    frame = load_frame(args.baseline_dir / "training_matrix.npz", features)
    holdout = reconstruct_holdout(frame, baseline_model, args.random_seed)
    et = frame["cluster_Et"][holdout]
    cent = frame["centrality"][holdout]
    analysis = np.isfinite(et) & (et >= 15.0) & (et < 35.0) & np.isfinite(cent) & (cent >= 0.0) & (cent < 80.0)
    idx = holdout[np.flatnonzero(analysis)]
    local_et = frame["cluster_Et"][idx]
    local_cent = frame["centrality"][idx]
    signal, inclusive = class_masks(frame, idx)
    filtered = filtered_indices(frame, list(baseline_model["features"]), baseline_model.get("pt_range"), [0.0, 80.0])
    filtered_mask = np.zeros(len(frame["is_signal"]), dtype=bool)
    filtered_mask[filtered] = True
    signal_all, inclusive_all = class_masks(frame, np.arange(len(frame["is_signal"]), dtype="int64"))
    source_keep_idx = np.flatnonzero(filtered_mask & (signal_all | inclusive_all))
    display_y = np.zeros(len(source_keep_idx), dtype="int32")
    display_y[signal_all[source_keep_idx]] = 1
    source_weights, weight_report = display_weights(frame, source_keep_idx, display_y)
    full_weights = np.full(len(frame["is_signal"]), np.nan, dtype="float64")
    full_weights[source_keep_idx] = source_weights
    local_weights = full_weights[idx]
    class_row = signal | inclusive
    if not np.all(np.isfinite(local_weights[class_row]) & (local_weights[class_row] > 0.0)):
        raise SystemExit("Non-finite source-class weights in analysis holdout")
    score_by_product: dict[str, np.ndarray] = {}
    route_inventory: dict[str, dict] = {}
    score_by_product["global"] = score_model(Path(baseline_model["output_xgb_json"]), frame, list(baseline_model["features"]), idx, args.batch_size)
    route_inventory["global"] = {
        "model_id": baseline_model["model_id"],
        "route_count": 1,
        "features": list(baseline_model["features"]),
        "output_xgb_json": baseline_model["output_xgb_json"],
    }
    for product_key, registry_product in PRODUCT_TO_REGISTRY.items():
        score_by_product[product_key], route_inventory[product_key] = score_routed_product(
            registry_product,
            models_by_product(binned_registry, registry_product),
            frame,
            idx,
            args.batch_size,
        )
    groups = build_groups()
    cells = []
    for group in groups:
        mask = group_mask(group["kind"], group, local_et, local_cent)
        if not np.any(mask):
            continue
        group_signal = signal & mask
        group_inclusive = inclusive & mask
        for product_key in PRODUCTS:
            record = {
                "group_kind": group["kind"],
                "group_key": group["key"],
                "group_label": group["label"],
                "product": product_key,
                "metrics": metrics_for(score_by_product[product_key], group_signal, group_inclusive, local_weights),
            }
            if "et_range" in group:
                record["et_range"] = group["et_range"]
            if "cent_range" in group:
                record["cent_range"] = group["cent_range"]
            cells.append(record)
    return {
        "schema": "THE8_ROUTED_BDT_COMPREHENSIVE_BASELINE_HOLDOUT_METRICS_V1",
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "source": "remote read-only scoring on baseline 10% row holdout reconstructed from THE57 baseline training_matrix.npz",
        "selection": "baseline 10% row holdout; 15 <= cluster_Et < 35; 0 <= centrality < 80",
        "labels": "Signal MC = embeddedPhoton source rows with is_signal==1; Inclusive MC = embeddedJet source rows, no truth-background filter",
        "metric_weighting": "PPG12-style display weights for source-class diagnostics: class balance plus eta and ET inverse-PDF weights",
        "weight_report": weight_report,
        "baseline_dir": str(args.baseline_dir),
        "binned_dir": str(args.binned_dir),
        "baseline_model": BASELINE_MODEL_ID,
        "holdout_reconstruction": {
            "mode": "row",
            "random_seed": int(args.random_seed),
            "test_size": 0.10,
            "holdout_rows_before_0_80_cent_filter": int(len(holdout)),
            "selected_source_kinematic_rows": int(len(idx)),
            "signal_rows": int(signal.sum()),
            "inclusive_rows": int(inclusive.sum()),
        },
        "products": {
            "global": BASELINE_MODEL_ID,
            "cent3": "base14_perCent3",
            "cent7": "base14_perCent7",
            "et": "base14_perEt",
            "etcent3": "base14_perEtCent3",
            "etcent7": "base14_perEtCent7",
        },
        "bins": {
            "coarse_cent": COARSE_CENT_BINS,
            "fine_cent": FINE_CENT_BINS,
            "et": ET_BINS,
        },
        "route_inventory": route_inventory,
        "cells": cells,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline-dir", type=Path, default=DEFAULT_BASELINE_DIR)
    parser.add_argument("--binned-dir", type=Path, default=DEFAULT_BINNED_DIR)
    parser.add_argument("--random-seed", type=int, default=13)
    parser.add_argument("--batch-size", type=int, default=200000)
    args = parser.parse_args()
    payload = build_payload(args)
    print(json.dumps(payload, allow_nan=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
