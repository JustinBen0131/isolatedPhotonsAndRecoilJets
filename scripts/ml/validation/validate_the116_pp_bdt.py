#!/usr/bin/env python3
"""Certify the THE-116 pp BDT on its deterministic event-group holdout."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.metrics import roc_auc_score

_HERE = Path(__file__).resolve()
_TRAINING = _HERE.parents[1] / "training"
sys.path.insert(0, str(_TRAINING))
from train_auau_photon_bdt import global_event_keys, stable_seed  # noqa: E402


EXPECTED_FEATURES = [
    "cluster_Et",
    "cluster_weta_cogx",
    "cluster_wphi_cogx",
    "vertexz",
    "cluster_Eta",
    "e11_over_e33",
    "cluster_et1",
    "cluster_et2",
    "cluster_et3",
    "cluster_et4",
    "e32_over_e35",
]
SIGNAL_SOURCES = {"run28_photonjet5", "run28_photonjet10", "run28_photonjet20"}
BACKGROUND_SOURCES = {"run28_jet8", "run28_jet12", "run28_jet20", "run28_jet30"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--candidate-registry", type=Path, required=True)
    parser.add_argument("--control-registry", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--pt-edges", default="15,17,19,21,23,25,28,32,35")
    parser.add_argument("--eta-edges", default="-0.7,-0.35,0,0.35,0.7")
    parser.add_argument("--runtime-tolerance", type=float, default=2.0e-7)
    parser.add_argument("--max-auc-regression", type=float, default=0.005)
    parser.add_argument("--max-binned-auc-regression", type=float, default=0.02)
    parser.add_argument("--max-train-holdout-auc-gap", type=float, default=0.12)
    parser.add_argument("--surface-max-rms", type=float, default=0.015)
    parser.add_argument("--surface-max-residual", type=float, default=0.025)
    parser.add_argument("--surface-max-efficiency-error", type=float, default=0.02)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_registry(path: Path) -> tuple[dict, dict]:
    payload = json.loads(path.read_text())
    if payload.get("status") != "READY" or payload.get("model_count") != 1:
        raise SystemExit(f"registry is not one-model READY: {path}")
    model = payload["models"][0]
    if model.get("features") != EXPECTED_FEATURES:
        raise SystemExit(f"unexpected feature order in {path}")
    return payload, model


def parse_edges(text: str) -> np.ndarray:
    edges = np.asarray([float(item) for item in text.split(",")], dtype="float64")
    if len(edges) < 2 or not np.all(np.diff(edges) > 0):
        raise SystemExit(f"invalid edges: {text}")
    return edges


def weighted_threshold(scores: np.ndarray, weights: np.ndarray, efficiency: float) -> float:
    order = np.argsort(scores, kind="mergesort")
    sorted_scores = scores[order]
    sorted_weights = weights[order]
    cumulative = np.cumsum(sorted_weights)
    target = (1.0 - efficiency) * cumulative[-1]
    index = int(np.searchsorted(cumulative, target, side="left"))
    return float(sorted_scores[min(index, len(sorted_scores) - 1)])


def weighted_efficiency(mask: np.ndarray, weights: np.ndarray) -> float:
    denominator = float(np.sum(weights))
    return float(np.sum(weights[mask]) / denominator) if denominator > 0 else math.nan


def score_model(model: dict, matrix: np.ndarray) -> np.ndarray:
    booster = xgb.Booster()
    booster.load_model(str(model["output_xgb_json"]))
    return np.asarray(booster.predict(xgb.DMatrix(matrix)), dtype="float64")


def tmva_scores(path: Path, matrix: np.ndarray, chunk_size: int = 100_000) -> np.ndarray:
    import ROOT  # pylint: disable=import-outside-toplevel

    runtime = ROOT.TMVA.Experimental.RBDT("myBDT", str(path))
    observed = np.empty(len(matrix), dtype="float64")
    for start in range(0, len(matrix), chunk_size):
        stop = min(len(matrix), start + chunk_size)
        observed[start:stop] = np.asarray(runtime.Compute(matrix[start:stop]), dtype="float64").reshape(-1)
    return observed


def event_masks(frame: pd.DataFrame, model_id: str, seed: int = 42) -> tuple[np.ndarray, np.ndarray, dict]:
    keys, columns = global_event_keys(frame, require_file_qualified=False)
    unique = np.unique(keys)
    hashes = np.asarray([stable_seed("event50", seed, model_id, key) for key in unique], dtype=np.uint64)
    order = np.argsort(hashes, kind="mergesort")
    n_test = int(round(len(unique) * 0.5))
    test_keys = set(unique[order[:n_test]].tolist())
    test = np.asarray([key in test_keys for key in keys], dtype=bool)
    train = ~test
    overlap = len(set(keys[train].tolist()).intersection(keys[test].tolist()))
    return train, test, {
        "event_key_columns": columns,
        "unique_events": int(len(unique)),
        "train_events": int(len(set(keys[train].tolist()))),
        "test_events": int(len(set(keys[test].tolist()))),
        "train_rows": int(train.sum()),
        "test_rows": int(test.sum()),
        "event_overlap": int(overlap),
    }


def binned_auc(labels, scores, weights, values, edges) -> list[dict]:
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        mask = (values >= lo) & (values < hi)
        classes = np.unique(labels[mask])
        rows.append({
            "lo": float(lo),
            "hi": float(hi),
            "rows": int(mask.sum()),
            "signal": int(np.sum(mask & (labels == 1))),
            "background": int(np.sum(mask & (labels == 0))),
            "auc": float(roc_auc_score(labels[mask], scores[mask], sample_weight=weights[mask])) if len(classes) == 2 else None,
        })
    return rows


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    candidate_registry, candidate = load_registry(args.candidate_registry)
    control_registry, control = load_registry(args.control_registry)
    if candidate.get("pt_range") != [15.0, 35.0] or control.get("pt_range") != [5.0, 35.0]:
        raise SystemExit("candidate/control model domains are not the frozen 15-35 and 5-35 pair")

    data = np.load(args.cache, allow_pickle=True)
    required = set(EXPECTED_FEATURES + [
        "is_signal", "cluster_Et", "cluster_Eta", "source_sample", "input_file_index", "run", "evt",
        "__ppg12_exact_training_weight",
    ])
    missing = sorted(required.difference(data.files))
    if missing:
        raise SystemExit(f"weighted cache missing fields: {missing}")
    frame = pd.DataFrame({name: data[name] for name in required})
    source = frame["source_sample"].astype(str).to_numpy()
    labels = frame["is_signal"].to_numpy(dtype="int8")
    domain = (frame["cluster_Et"].to_numpy() >= 15.0) & (frame["cluster_Et"].to_numpy() < 35.0)
    finite = np.isfinite(frame[EXPECTED_FEATURES].to_numpy(dtype="float64")).all(axis=1)
    frame = frame.loc[domain & finite].reset_index(drop=True)
    source = frame["source_sample"].astype(str).to_numpy()
    labels = frame["is_signal"].to_numpy(dtype="int8")
    weights = frame["__ppg12_exact_training_weight"].to_numpy(dtype="float64")
    matrix = np.ascontiguousarray(frame[EXPECTED_FEATURES].to_numpy(dtype="float32"))

    observed_sources = set(source.tolist())
    wrong_signal = int(np.sum((labels == 1) & ~np.isin(source, sorted(SIGNAL_SOURCES))))
    wrong_background = int(np.sum((labels == 0) & ~np.isin(source, sorted(BACKGROUND_SOURCES))))
    source_gate = (
        observed_sources == (SIGNAL_SOURCES | BACKGROUND_SOURCES)
        and wrong_signal == 0
        and wrong_background == 0
    )

    train_mask, test_mask, split = event_masks(frame, str(candidate["model_id"]), seed=42)
    if split["event_overlap"] != 0:
        raise SystemExit("event-group split overlap is nonzero")
    x_test = matrix[test_mask]
    y_test = labels[test_mask]
    w_test = weights[test_mask]
    current_scores = score_model(candidate, x_test)
    control_scores = score_model(control, x_test)
    current_auc = float(roc_auc_score(y_test, current_scores, sample_weight=w_test))
    control_auc = float(roc_auc_score(y_test, control_scores, sample_weight=w_test))
    auc_delta = current_auc - control_auc

    tmva = tmva_scores(Path(candidate["output_tmva"]), x_test)
    difference = np.abs(tmva - current_scores)
    parity = {
        "rows": int(len(x_test)),
        "max_abs_difference": float(np.max(difference)),
        "mean_abs_difference": float(np.mean(difference)),
        "p99_abs_difference": float(np.quantile(difference, 0.99)),
        "nonfinite_rows": int(np.sum(~np.isfinite(tmva) | ~np.isfinite(current_scores))),
        "tolerance": args.runtime_tolerance,
    }
    parity["status"] = "PASS" if parity["nonfinite_rows"] == 0 and parity["max_abs_difference"] <= args.runtime_tolerance else "FAIL"

    report = candidate["report"]
    overfit = report["overfit_diagnostics"]
    reported_holdout_auc = float(overfit["holdout_auc"])
    holdout_reconstruction_delta = current_auc - reported_holdout_auc
    overfit_gate = float(overfit["auc_gap_train_minus_holdout"]) <= args.max_train_holdout_auc_gap

    pt_edges = parse_edges(args.pt_edges)
    eta_edges = parse_edges(args.eta_edges)
    pt = frame.loc[test_mask, "cluster_Et"].to_numpy(dtype="float64")
    eta = frame.loc[test_mask, "cluster_Eta"].to_numpy(dtype="float64")
    et_current = binned_auc(y_test, current_scores, w_test, pt, pt_edges)
    et_control = binned_auc(y_test, control_scores, w_test, pt, pt_edges)
    eta_current = binned_auc(y_test, current_scores, w_test, eta, eta_edges)
    eta_control = binned_auc(y_test, control_scores, w_test, eta, eta_edges)
    pt_auc_deltas = [
        current["auc"] - control_row["auc"]
        for current, control_row in zip(et_current, et_control)
        if current["auc"] is not None and control_row["auc"] is not None
    ]
    eta_auc_deltas = [
        current["auc"] - control_row["auc"]
        for current, control_row in zip(eta_current, eta_control)
        if current["auc"] is not None and control_row["auc"] is not None
    ]

    wp_rows = []
    for efficiency, label in ((0.90, "WP90"), (0.80, "WP80"), (0.70, "WP70")):
        for lo, hi in zip(pt_edges[:-1], pt_edges[1:]):
            local = (pt >= lo) & (pt < hi)
            signal = local & (y_test == 1)
            background = local & (y_test == 0)
            threshold = weighted_threshold(current_scores[signal], w_test[signal], efficiency)
            wp_rows.append({
                "wp": label,
                "target_signal_efficiency": efficiency,
                "pt_lo": float(lo),
                "pt_hi": float(hi),
                "pt_center": float(0.5 * (lo + hi)),
                "threshold": threshold,
                "achieved_signal_efficiency": weighted_efficiency(current_scores[signal] > threshold, w_test[signal]),
                "background_acceptance": weighted_efficiency(current_scores[background] > threshold, w_test[background]),
                "signal_rows": int(signal.sum()),
                "background_rows": int(background.sum()),
            })

    surfaces = {}
    for label in ("WP90", "WP80", "WP70"):
        rows = [row for row in wp_rows if row["wp"] == label]
        centers = np.asarray([row["pt_center"] for row in rows])
        thresholds = np.asarray([row["threshold"] for row in rows])
        slope, intercept = np.polyfit(centers, thresholds, 1)
        fitted = intercept + slope * centers
        residual = thresholds - fitted
        target = rows[0]["target_signal_efficiency"]
        applied = intercept + slope * pt
        signal_all = y_test == 1
        achieved = []
        for lo, hi in zip(pt_edges[:-1], pt_edges[1:]):
            local = signal_all & (pt >= lo) & (pt < hi)
            achieved.append(weighted_efficiency(current_scores[local] > applied[local], w_test[local]))
        max_eff_error = float(np.max(np.abs(np.asarray(achieved) - target)))
        rms = float(np.sqrt(np.mean(np.square(residual))))
        max_residual = float(np.max(np.abs(residual)))
        accepted = rms <= args.surface_max_rms and max_residual <= args.surface_max_residual and max_eff_error <= args.surface_max_efficiency_error
        surfaces[label] = {
            "status": "ACCEPTED" if accepted else "BINNED_THRESHOLDS_ONLY",
            "intercept": float(intercept),
            "slope_per_gev": float(slope),
            "rms_residual": rms,
            "max_abs_residual": max_residual,
            "max_abs_efficiency_error": max_eff_error,
            "achieved_efficiency_by_bin": achieved,
        }

    critical_gates = {
        "candidate_registry_ready": candidate_registry.get("status") == "READY",
        "control_registry_ready": control_registry.get("status") == "READY",
        "feature_order": candidate.get("features") == EXPECTED_FEATURES,
        "source_role_label_closure": bool(source_gate),
        "finite_model_inputs": bool(np.isfinite(matrix).all()),
        "finite_positive_weights": bool(np.isfinite(weights).all() and np.all(weights > 0)),
        "event_group_overlap_zero": split["event_overlap"] == 0,
        "reported_holdout_reconstructed": abs(holdout_reconstruction_delta) <= 1.0e-10,
        "runtime_parity": parity["status"] == "PASS",
        "no_material_auc_regression": auc_delta >= -args.max_auc_regression,
        "overtraining_gate": bool(overfit_gate),
        "all_pt_bins_have_both_classes": all(row["signal"] > 0 and row["background"] > 0 for row in et_current),
        "all_eta_bins_have_both_classes": all(row["signal"] > 0 and row["background"] > 0 for row in eta_current),
        "pt_resolved_stability": bool(
            pt_auc_deltas and min(pt_auc_deltas) >= -args.max_binned_auc_regression
        ),
        "eta_resolved_stability": bool(
            eta_auc_deltas and min(eta_auc_deltas) >= -args.max_binned_auc_regression
        ),
    }
    status = "PASS" if all(critical_gates.values()) else "FAIL"
    payload = {
        "schema": "THE116_PP_BDT_VALIDATION_V1",
        "status": status,
        "promotion_status": "NOT_PROMOTED",
        "model_domain_gev": [15.0, 35.0],
        "feature_order": EXPECTED_FEATURES,
        "critical_gates": critical_gates,
        "source_label_closure": {
            "observed_sources": sorted(observed_sources),
            "wrong_signal_rows": wrong_signal,
            "wrong_background_rows": wrong_background,
        },
        "split": split,
        "performance": {
            "candidate_weighted_holdout_auc": current_auc,
            "control_weighted_same_holdout_auc": control_auc,
            "candidate_minus_control_auc": auc_delta,
            "reported_holdout_auc": reported_holdout_auc,
            "reported_holdout_reconstruction_delta": holdout_reconstruction_delta,
            "train_auc": float(overfit["train_auc"]),
            "train_minus_holdout_auc": float(overfit["auc_gap_train_minus_holdout"]),
            "et_candidate": et_current,
            "et_control": et_control,
            "eta_candidate": eta_current,
            "eta_control": eta_control,
            "pt_candidate_minus_control_auc": pt_auc_deltas,
            "eta_candidate_minus_control_auc": eta_auc_deltas,
        },
        "runtime_parity": parity,
        "working_points": wp_rows,
        "runtime_surfaces": surfaces,
        "threshold_authority": "exact binned weighted holdout thresholds",
        "provenance": {
            "cache": str(args.cache), "cache_sha256": sha256(args.cache),
            "candidate_registry": str(args.candidate_registry), "candidate_registry_sha256": sha256(args.candidate_registry),
            "control_registry": str(args.control_registry), "control_registry_sha256": sha256(args.control_registry),
            "candidate_xgb_sha256": sha256(Path(candidate["output_xgb_json"])),
            "candidate_tmva_sha256": sha256(Path(candidate["output_tmva"])),
            "control_xgb_sha256": sha256(Path(control["output_xgb_json"])),
        },
        "gates": {
            "runtime_tolerance": args.runtime_tolerance,
            "max_auc_regression": args.max_auc_regression,
            "max_binned_auc_regression": args.max_binned_auc_regression,
            "max_train_holdout_auc_gap": args.max_train_holdout_auc_gap,
            "surface_max_rms": args.surface_max_rms,
            "surface_max_residual": args.surface_max_residual,
            "surface_max_efficiency_error": args.surface_max_efficiency_error,
        },
    }
    json_path = args.outdir / "the116_pp_bdt_validation.json"
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    csv_path = args.outdir / "the116_pp_bdt_working_points.csv"
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(wp_rows[0]))
        writer.writeheader()
        writer.writerows(wp_rows)
    print(json.dumps({
        "status": status,
        "candidate_auc": current_auc,
        "control_auc": control_auc,
        "auc_delta": auc_delta,
        "runtime_max_abs": parity["max_abs_difference"],
        "train_holdout_auc_gap": float(overfit["auc_gap_train_minus_holdout"]),
        "surfaces": {key: value["status"] for key, value in surfaces.items()},
        "json": str(json_path),
        "working_points": str(csv_path),
    }, indent=2, sort_keys=True))
    if status != "PASS":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
