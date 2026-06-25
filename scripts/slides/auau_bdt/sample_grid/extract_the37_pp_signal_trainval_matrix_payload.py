#!/usr/bin/env python3
"""Extract THE-37 pp fixed-inclusive signal train/validation matrix payload.

This script is intended to be streamed to SDCC and run read-only against the
existing THE-37 pp sample-grid products. It writes compact JSON to stdout.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import math
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import xgboost as xgb
from sklearn.metrics import brier_score_loss, log_loss, roc_auc_score


REMOTE_REPO = Path("/sphenix/u/patsfan753/scratch/thesisAnalysis")
DEFAULT_GRID_BASE = (
    REMOTE_REPO
    / "dataOutput/ppPhotonMLPipeline/pp_currentian_basev3e_sample_grid_20260620_envfix2"
)
FULL_BKG_TAG = "8_12_20_30_40"
FULL_BKG_LABEL = "Jet8+12+20+30+40"
SIGNAL_SETS = [
    ("photon5", "Photon5", ("run28_photonjet5",), "sig5"),
    ("photon5_10", "Photon5+10", ("run28_photonjet5", "run28_photonjet10"), "sig5_10"),
    (
        "photon5_10_20",
        "Photon5+10+20",
        ("run28_photonjet5", "run28_photonjet10", "run28_photonjet20"),
        "sig5_10_20",
    ),
]
BKG_SETS = [
    ("jet8_12", "Jet8+12", ("run28_jet8", "run28_jet12"), "8_12"),
    ("jet8_12_20", "Jet8+12+20", ("run28_jet8", "run28_jet12", "run28_jet20"), "8_12_20"),
    (
        "jet8_12_20_30",
        "Jet8+12+20+30",
        ("run28_jet8", "run28_jet12", "run28_jet20", "run28_jet30"),
        "8_12_20_30",
    ),
    (
        "jet8_12_20_30_40",
        "Jet8+12+20+30+40",
        ("run28_jet8", "run28_jet12", "run28_jet20", "run28_jet30", "run28_jet40"),
        "8_12_20_30_40",
    ),
]
INCLUSIVE_SAMPLES = ("run28_jet8", "run28_jet12", "run28_jet20", "run28_jet30", "run28_jet40")
FULL_SIGNAL_SAMPLES = ("run28_photonjet5", "run28_photonjet10", "run28_photonjet20")
FULL_SIGNAL_TAG = "sig5_10_20"
FULL_SIGNAL_LABEL = "Photon5+10+20"
BINS = np.linspace(0.0, 1.0, 51)


def add_training_imports() -> None:
    for path in [
        REMOTE_REPO / "scripts/ml/training",
        REMOTE_REPO / "scripts/ml/validation",
        REMOTE_REPO / "scripts/ml/stacking",
    ]:
        if str(path) not in sys.path:
            sys.path.insert(0, str(path))


def json_ready(value):
    if isinstance(value, dict):
        return {str(k): json_ready(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(v) for v in value]
    if isinstance(value, np.ndarray):
        return json_ready(value.tolist())
    if isinstance(value, (np.floating, float)):
        return None if not math.isfinite(float(value)) else float(value)
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    return value


def read_registry(lane_dir: Path) -> dict:
    reg_path = lane_dir / "models/bdt_ppg12_currentIAN_basev3E/model_registry.json"
    if not reg_path.is_file():
        raise SystemExit(f"missing registry: {reg_path}")
    reg = json.loads(reg_path.read_text())
    models = reg.get("models") or []
    if len(models) != 1:
        raise SystemExit(f"expected one model in {reg_path}, got {len(models)}")
    model = models[0]
    features = [str(item) for item in model.get("features") or []]
    xgb_path = Path(model.get("output_xgb_json") or model.get("report", {}).get("output_xgb_json", ""))
    if not features:
        raise SystemExit(f"model in {reg_path} has no features")
    if not xgb_path.is_file():
        raise SystemExit(f"missing XGBoost model for {lane_dir.name}: {xgb_path}")
    return {"registry": reg_path, "model": model, "features": features, "xgb_path": xgb_path}


def expand_manifest(path: Path) -> list[Path]:
    if not path.is_file():
        raise SystemExit(f"missing manifest: {path}")
    out = []
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if line and not line.startswith("#"):
            p = Path(line)
            if not p.is_file():
                raise SystemExit(f"manifest entry missing: {p}")
            out.append(p)
    if not out:
        raise SystemExit(f"manifest is empty: {path}")
    return out


def source_subset(source: np.ndarray, samples: tuple[str, ...]) -> np.ndarray:
    mask = np.zeros(len(source), dtype=bool)
    for sample in samples:
        mask |= source == sample
    return mask


def score_model(model_path: Path, frame, features: list[str], batch_size: int) -> np.ndarray:
    booster = xgb.Booster()
    booster.load_model(str(model_path))
    out = np.empty(len(frame), dtype="float32")
    for start in range(0, len(frame), batch_size):
        stop = min(start + batch_size, len(frame))
        x = np.column_stack(
            [frame[name].to_numpy(dtype="float32", copy=False)[start:stop] for name in features]
        )
        out[start:stop] = booster.predict(xgb.DMatrix(x)).astype("float32")
    return out


def weighted_hist(values: np.ndarray, weights: np.ndarray) -> dict:
    good = np.isfinite(values) & np.isfinite(weights) & (weights > 0.0)
    counts, _ = np.histogram(values[good], bins=BINS)
    weighted, _ = np.histogram(values[good], bins=BINS, weights=weights[good])
    total = float(weighted.sum())
    widths = np.diff(BINS)
    density = np.zeros_like(weighted, dtype="float64")
    if total > 0.0:
        density = weighted / total / widths
    return {
        "entries": int(good.sum()),
        "sum_weight": total,
        "counts": counts.astype(int).tolist(),
        "weighted_counts": weighted.astype(float).tolist(),
        "density": density.astype(float).tolist(),
    }


def weighted_median(values: np.ndarray, weights: np.ndarray) -> float:
    good = np.isfinite(values) & np.isfinite(weights) & (weights > 0.0)
    if not good.any():
        return math.nan
    vals = values[good]
    w = weights[good]
    order = np.argsort(vals)
    vals = vals[order]
    w = w[order]
    cdf = np.cumsum(w)
    target = 0.5 * float(cdf[-1])
    return float(vals[int(np.searchsorted(cdf, target, side="left"))])


def wp_threshold_for_signal_efficiency(scores: np.ndarray, weights: np.ndarray, target_eff: float) -> dict:
    good = np.isfinite(scores) & np.isfinite(weights) & (weights > 0.0)
    if not good.any():
        return {"threshold": math.nan, "signal_efficiency": math.nan}
    vals = scores[good]
    w = weights[good]
    order = np.argsort(vals)[::-1]
    vals = vals[order]
    w = w[order]
    cdf = np.cumsum(w)
    total = float(cdf[-1])
    idx = int(np.searchsorted(cdf, target_eff * total, side="left"))
    idx = min(idx, len(vals) - 1)
    threshold = float(vals[idx])
    eff = float(w[vals >= threshold].sum() / total)
    return {"threshold": threshold, "signal_efficiency": eff}


def metrics_for(signal_scores: np.ndarray, inclusive_scores: np.ndarray) -> dict:
    sig_w = np.ones(len(signal_scores), dtype="float64")
    inc_w = np.ones(len(inclusive_scores), dtype="float64")
    y = np.concatenate([np.ones(len(signal_scores), dtype="int32"), np.zeros(len(inclusive_scores), dtype="int32")])
    scores = np.concatenate([signal_scores, inclusive_scores]).astype("float64", copy=False)
    weights = np.concatenate([sig_w, inc_w])
    good = np.isfinite(scores)
    y = y[good]
    scores = scores[good]
    weights = weights[good]
    clipped = np.clip(scores, 1.0e-6, 1.0 - 1.0e-6)
    wp = wp_threshold_for_signal_efficiency(signal_scores, sig_w, 0.80)
    threshold = wp["threshold"]
    inc_pass = float(np.mean(inclusive_scores[np.isfinite(inclusive_scores)] >= threshold))
    return {
        "auc": float(roc_auc_score(y, scores, sample_weight=weights)),
        "logloss": float(log_loss(y, clipped, sample_weight=weights, labels=[0, 1])),
        "brier": float(brier_score_loss(y, clipped, sample_weight=weights)),
        "median_signal_score": weighted_median(signal_scores, sig_w),
        "median_background_score": weighted_median(inclusive_scores, inc_w),
        "median_gap": weighted_median(signal_scores, sig_w) - weighted_median(inclusive_scores, inc_w),
        "wp80_threshold": float(threshold),
        "wp80_signal_efficiency": float(wp["signal_efficiency"]),
        "wp80_background_fake_rate": inc_pass,
    }


def build_cell(row, col, score: np.ndarray, source: np.ndarray, y: np.ndarray, model_path: Path, lane: str) -> dict:
    row_key, row_label, _row_samples, _row_prefix = row
    col_key, col_label, col_samples, _col_prefix = col
    signal_mask = source_subset(source, col_samples) & (y == 1)
    inclusive_mask = source_subset(source, INCLUSIVE_SAMPLES)
    signal_scores = score[signal_mask]
    inclusive_scores = score[inclusive_mask]
    sig_w = np.ones(len(signal_scores), dtype="float64")
    inc_w = np.ones(len(inclusive_scores), dtype="float64")
    return {
        "row_key": row_key,
        "row_label": row_label,
        "col_key": col_key,
        "col_label": col_label,
        "lane": lane,
        "model_path": str(model_path),
        "bin_edges": BINS.astype(float).tolist(),
        "signal": weighted_hist(signal_scores, sig_w),
        "background": weighted_hist(inclusive_scores, inc_w),
        "metrics": metrics_for(signal_scores, inclusive_scores),
        "entries": {
            "signal": int(len(signal_scores)),
            "background": int(len(inclusive_scores)),
            "total": int(len(signal_scores) + len(inclusive_scores)),
        },
        "sum_weight": {
            "signal": float(len(signal_scores)),
            "background": float(len(inclusive_scores)),
            "total": float(len(signal_scores) + len(inclusive_scores)),
        },
    }


def build_inclusive_cell(row, col, score: np.ndarray, source: np.ndarray, y: np.ndarray, model_path: Path, lane: str) -> dict:
    row_key, row_label, row_samples, _row_prefix = row
    col_key, col_label, _col_samples, _col_prefix = col
    signal_mask = source_subset(source, FULL_SIGNAL_SAMPLES) & (y == 1)
    inclusive_mask = source_subset(source, row_samples)
    signal_scores = score[signal_mask]
    inclusive_scores = score[inclusive_mask]
    sig_w = np.ones(len(signal_scores), dtype="float64")
    inc_w = np.ones(len(inclusive_scores), dtype="float64")
    return {
        "row_key": row_key,
        "row_label": row_label,
        "col_key": col_key,
        "col_label": col_label,
        "lane": lane,
        "model_path": str(model_path),
        "bin_edges": BINS.astype(float).tolist(),
        "signal": weighted_hist(signal_scores, sig_w),
        "background": weighted_hist(inclusive_scores, inc_w),
        "metrics": metrics_for(signal_scores, inclusive_scores),
        "entries": {
            "signal": int(len(signal_scores)),
            "background": int(len(inclusive_scores)),
            "total": int(len(signal_scores) + len(inclusive_scores)),
        },
        "sum_weight": {
            "signal": float(len(signal_scores)),
            "background": float(len(inclusive_scores)),
            "total": float(len(signal_scores) + len(inclusive_scores)),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grid-base", type=Path, default=DEFAULT_GRID_BASE)
    parser.add_argument("--tree", default="AuAuPhotonIDTrainingTree")
    parser.add_argument("--pt-range", default="15:35")
    parser.add_argument("--centrality-range", default="-1:0")
    parser.add_argument("--random-seed", type=int, default=42)
    parser.add_argument("--max-load-rows-per-class", type=int, default=1000000)
    parser.add_argument("--max-load-rows", type=int, default=0)
    parser.add_argument("--batch-size", type=int, default=250000)
    args = parser.parse_args()
    args.label_branch = "is_signal"
    args.max_rows = 0

    add_training_imports()
    from train_auau_photon_bdt import add_derived_features, expand_required_columns, load_frame
    from validate_pp_photon_ml_tables import parse_range, select_rows

    full_lane = f"sig5_10_20_bkg{FULL_BKG_TAG}"
    full_dir = args.grid_base / full_lane
    signal_paths = expand_manifest(full_dir / "signal_roots_currentIAN.list")
    inclusive_paths = expand_manifest(full_dir / "inclusive_roots_currentIAN.list")

    row_models = {}
    bkg_models = {}
    feature_union = {"is_signal", "cluster_Et", "cluster_Eta", "centrality", "vertexz"}
    for row in SIGNAL_SETS:
        lane = f"{row[3]}_bkg{FULL_BKG_TAG}"
        info = read_registry(args.grid_base / lane)
        row_models[row[0]] = {"lane": lane, **info}
        feature_union.update(expand_required_columns(info["features"]))
    for col in BKG_SETS:
        lane = f"{FULL_SIGNAL_TAG}_bkg{col[3]}"
        info = read_registry(args.grid_base / lane)
        bkg_models[col[0]] = {"lane": lane, **info}
        feature_union.update(expand_required_columns(info["features"]))

    paths = signal_paths + inclusive_paths
    frame, optional_seen = load_frame(
        paths,
        args.tree,
        sorted(feature_union),
        ["run", "evt"],
        None,
        None,
        skip_missing_tree=False,
        label_branch="is_signal",
        max_load_rows_per_class=args.max_load_rows_per_class,
        max_load_rows=args.max_load_rows,
        load_sample_seed=args.random_seed,
    )
    frame = add_derived_features(frame)
    frame = select_rows(frame, args)
    if len(frame) == 0:
        raise SystemExit("no validation rows selected")

    source = frame["source_sample"].astype(str).to_numpy(dtype=object)
    y = frame["is_signal"].to_numpy(dtype="int32", copy=False)
    source_counts = Counter(str(item) for item in source)
    source_signal_counts = {
        sample: int(((source == sample) & (y == 1)).sum()) for sample in sorted(source_counts)
    }
    source_all_counts = {sample: int(source_counts[sample]) for sample in sorted(source_counts)}
    inclusive_truth_signal_rows = int((source_subset(source, INCLUSIVE_SAMPLES) & (y == 1)).sum())

    cells = []
    inclusive_cells = []
    max_density = 0.0
    inclusive_max_density = 0.0
    for row in SIGNAL_SETS:
        info = row_models[row[0]]
        score = score_model(info["xgb_path"], frame, info["features"], args.batch_size)
        for col in SIGNAL_SETS:
            cell = build_cell(row, col, score, source, y, info["xgb_path"], info["lane"])
            max_density = max(
                max_density,
                max(cell["signal"]["density"] or [0.0]),
                max(cell["background"]["density"] or [0.0]),
            )
            cells.append(cell)
    for col in BKG_SETS:
        info = bkg_models[col[0]]
        score = score_model(info["xgb_path"], frame, info["features"], args.batch_size)
        for row in BKG_SETS:
            cell = build_inclusive_cell(row, col, score, source, y, info["xgb_path"], info["lane"])
            inclusive_max_density = max(
                inclusive_max_density,
                max(cell["signal"]["density"] or [0.0]),
                max(cell["background"]["density"] or [0.0]),
            )
            inclusive_cells.append(cell)

    payload = {
        "schema": "THE37_PP_SIGNAL_AND_INCLUSIVE_TRAINVAL_MATRIX_PAYLOAD_V1",
        "created_remote_time": dt.datetime.now(dt.UTC).isoformat(),
        "pp_signal_trainval": {
            "title": "pp 15-35 GeV: signal train/validation matrix with fixed inclusive MC",
            "fixed_inclusive_label": f"FIXED INCLUSIVE MC: {FULL_BKG_LABEL}",
            "row_label_prefix": "SIGNAL",
            "row_order": [{"key": key, "label": label} for key, label, _samples, _prefix in SIGNAL_SETS],
            "col_order": [{"key": key, "label": label} for key, label, _samples, _prefix in SIGNAL_SETS],
            "cells": cells,
            "ymax": min(24.0, max(18.0, float(max_density) * 1.04)),
            "validation_meta": {
                "class_definition": (
                    "Signal = selected photonjet source rows with is_signal==1; "
                    "Inclusive MC = all selected Jet8+12+20+30+40 source rows with no truth filter."
                ),
                "grid_base": str(args.grid_base),
                "fixed_signal_manifest": str(full_dir / "signal_roots_currentIAN.list"),
                "fixed_inclusive_manifest": str(full_dir / "inclusive_roots_currentIAN.list"),
                "pt_range": args.pt_range,
                "centrality_range": args.centrality_range,
                "rows_after_selection": int(len(frame)),
                "source_counts_after_selection": source_all_counts,
                "source_signal_counts_after_selection": source_signal_counts,
                "inclusive_truth_signal_rows_kept_as_inclusive": inclusive_truth_signal_rows,
                "optional_seen": sorted(optional_seen),
                "weights": "unit weights; source-aware read-only extraction from ROOT inputs",
                "max_load_rows_per_class": int(args.max_load_rows_per_class),
                "random_seed": int(args.random_seed),
            },
        },
        "pp_inclusive_trainval": {
            "title": "pp 15-35 GeV: fixed signal vs inclusive-MC train/validation sample",
            "fixed_signal_label": f"FIXED SIGNAL MC: {FULL_SIGNAL_LABEL}",
            "row_label_prefix": "INCLUSIVE",
            "row_order": [{"key": key, "label": label} for key, label, _samples, _prefix in BKG_SETS],
            "col_order": [{"key": key, "label": label} for key, label, _samples, _prefix in BKG_SETS],
            "cells": inclusive_cells,
            "ymax": min(24.0, max(18.0, float(inclusive_max_density) * 1.04)),
            "validation_meta": {
                "class_definition": (
                    "Signal = all selected Photon5+10+20 source rows with is_signal==1; "
                    "Inclusive MC = selected source-defined Jet8/12/20/30/40 rows per validation row, no truth filter."
                ),
                "grid_base": str(args.grid_base),
                "fixed_signal_manifest": str(full_dir / "signal_roots_currentIAN.list"),
                "fixed_inclusive_manifest": str(full_dir / "inclusive_roots_currentIAN.list"),
                "pt_range": args.pt_range,
                "centrality_range": args.centrality_range,
                "rows_after_selection": int(len(frame)),
                "source_counts_after_selection": source_all_counts,
                "source_signal_counts_after_selection": source_signal_counts,
                "inclusive_truth_signal_rows_kept_as_inclusive": inclusive_truth_signal_rows,
                "optional_seen": sorted(optional_seen),
                "weights": "unit weights; source-aware read-only extraction from ROOT inputs",
                "max_load_rows_per_class": int(args.max_load_rows_per_class),
                "random_seed": int(args.random_seed),
            },
        },
    }
    print("###THE37_PP_SIGNAL_TRAINVAL_PAYLOAD_BEGIN###")
    print(json.dumps(json_ready(payload), sort_keys=True))
    print("###THE37_PP_SIGNAL_TRAINVAL_PAYLOAD_END###")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
