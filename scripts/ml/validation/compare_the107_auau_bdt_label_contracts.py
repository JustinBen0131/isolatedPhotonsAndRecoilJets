#!/usr/bin/env python3
"""Validate and render the paired THE-107 Au+Au BDT label contracts.

The inputs are the exact row-level holdout caches written by the trainer for
the nominal isolated-prompt Au+Au model and the PPG12 source-role variant.
This script verifies model/config provenance, Python/TMVA score parity, native
held-out performance, cross-label behavior on the intersection of both native
holdouts, and prompt-photon score dependence on truth isolation.  Its primary
reader-facing artifact is the requested 2x3 centrality comparison.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


EXPECTED_FEATURES = [
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
EXPECTED_XGBOOST = {
    "n_estimators": 450,
    "max_depth": 4,
    "learning_rate": 0.035,
    "subsample": 0.85,
    "colsample_bytree": 0.85,
    "reg_alpha": 5.0,
    "reg_lambda": 0.3,
    "grow_policy": "lossguide",
    "max_bin": 256,
    "tree_method": "hist",
    "objective": "binary:logistic",
    "eval_metric": ["auc", "logloss"],
    "random_state": 13,
}
CENTRALITY_BINS = [(0.0, 20.0), (20.0, 50.0), (50.0, 80.0)]
SCORE_EDGES = np.linspace(0.0, 1.0, 51)
GREEN = "#159947"
ORANGE = "#D55E00"
BLUE = "#1479C9"
INK = "#18212B"
MUTED = "#5C6670"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def json_safe(value):
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def load_npz(path: Path) -> dict[str, np.ndarray]:
    if not path.is_file() or path.stat().st_size <= 0:
        raise SystemExit(f"missing holdout cache: {path}")
    with np.load(path, allow_pickle=True) as payload:
        return {name: np.asarray(payload[name]) for name in payload.files}


def required(payload: dict[str, np.ndarray], names: list[str], source: Path) -> None:
    missing = [name for name in names if name not in payload]
    if missing:
        raise SystemExit(f"{source} is missing required holdout fields: {missing}")


def weighted_quantile(values: np.ndarray, weights: np.ndarray, q: float) -> float:
    good = np.isfinite(values) & np.isfinite(weights) & (weights > 0.0)
    if not good.any():
        return math.nan
    values = np.asarray(values[good], dtype="float64")
    weights = np.asarray(weights[good], dtype="float64")
    order = np.argsort(values, kind="mergesort")
    values = values[order]
    cdf = np.cumsum(weights[order])
    target = float(np.clip(q, 0.0, 1.0)) * float(cdf[-1])
    index = int(np.searchsorted(cdf, target, side="left"))
    return float(values[min(max(index, 0), len(values) - 1)])


def weighted_auc(labels: np.ndarray, scores: np.ndarray, weights: np.ndarray | None) -> float:
    from sklearn.metrics import roc_auc_score

    good = np.isfinite(scores) & np.isin(labels, [0, 1])
    if weights is not None:
        good &= np.isfinite(weights) & (weights > 0.0)
    labels = labels[good]
    scores = scores[good]
    if len(np.unique(labels)) != 2:
        return math.nan
    return float(
        roc_auc_score(
            labels,
            scores,
            sample_weight=None if weights is None else weights[good],
        )
    )


def finite_mean(values: np.ndarray) -> float:
    values = np.asarray(values, dtype="float64")
    values = values[np.isfinite(values)]
    return float(np.mean(values)) if len(values) else math.nan


def wp80_metrics(labels: np.ndarray, scores: np.ndarray, weights: np.ndarray) -> dict:
    sig = labels == 1
    bkg = labels == 0
    if not sig.any() or not bkg.any():
        return {
            "threshold": math.nan,
            "signal_efficiency": math.nan,
            "background_acceptance": math.nan,
        }
    threshold = weighted_quantile(scores[sig], weights[sig], 0.20)
    sig_sum = float(np.sum(weights[sig]))
    bkg_sum = float(np.sum(weights[bkg]))
    return {
        "threshold": threshold,
        "signal_efficiency": float(np.sum(weights[sig & (scores > threshold)]) / sig_sum),
        "background_acceptance": float(np.sum(weights[bkg & (scores > threshold)]) / bkg_sum),
    }


def candidate_keys(payload: dict[str, np.ndarray], source: Path) -> np.ndarray:
    columns = [
        "source_role",
        "source_sample_code",
        "input_file_index",
        "run",
        "evt",
        "input_tree_entry",
    ]
    required(payload, columns, source)
    arrays = [np.asarray(payload[name], dtype="int64") for name in columns]
    key = np.rec.fromarrays(arrays, names=",".join(columns))
    if len(np.unique(key)) != len(key):
        raise SystemExit(f"candidate identity is not unique in {source}")
    return key


def score_xgboost(model_path: Path, x: np.ndarray) -> np.ndarray:
    import xgboost as xgb

    model = xgb.Booster()
    model.load_model(str(model_path))
    feature_names = model.feature_names
    matrix = xgb.DMatrix(
        np.asarray(x, dtype="float32"),
        feature_names=feature_names if feature_names else None,
    )
    return np.asarray(model.predict(matrix), dtype="float32")


def runtime_parity(tmva_path: Path, x: np.ndarray, xgb_score: np.ndarray, sample_size: int) -> dict:
    import ROOT

    n = min(int(sample_size), len(x))
    if n <= 0:
        raise SystemExit("runtime parity received an empty holdout")
    indices = np.unique(np.linspace(0, len(x) - 1, n, dtype="int64"))
    model = ROOT.TMVA.Experimental.RBDT("myBDT", str(tmva_path))
    tmva = np.full(len(indices), np.nan, dtype="float64")
    for out_index, row_index in enumerate(indices):
        vector = ROOT.std.vector("float")()
        for value in x[row_index]:
            vector.push_back(float(value))
        result = model.Compute(vector)
        if len(result):
            tmva[out_index] = float(result[0])
    reference = np.asarray(xgb_score[indices], dtype="float64")
    finite = np.isfinite(reference) & np.isfinite(tmva)
    if not finite.all():
        raise SystemExit(
            f"runtime parity produced non-finite scores for {int((~finite).sum())}/{len(finite)} rows"
        )
    delta = tmva - reference
    report = {
        "rows": int(len(indices)),
        "max_abs_difference": float(np.max(np.abs(delta))),
        "rms_difference": float(np.sqrt(np.mean(delta * delta))),
        "mean_difference": float(np.mean(delta)),
        "tolerance": 2.0e-6,
    }
    report["passed"] = report["max_abs_difference"] <= report["tolerance"]
    if not report["passed"]:
        raise SystemExit(f"Python/TMVA runtime parity failed for {tmva_path}: {report}")
    return report


def load_lane(
    name: str,
    holdout_path: Path,
    model_path: Path,
    tmva_path: Path,
    metadata_path: Path,
    runtime_rows: int,
) -> dict:
    payload = load_npz(holdout_path)
    required(
        payload,
        [
            "features",
            "x",
            "is_signal",
            "training_weight",
            "score_xgboost",
            "centrality",
            "cluster_Et",
            "truth_is_prompt",
            "truth_iso_et",
            "truth_iso_pass",
            "ppg12_source_role_label",
            "nominal_is_signal",
        ],
        holdout_path,
    )
    features = [str(item) for item in payload["features"].tolist()]
    if features != EXPECTED_FEATURES:
        raise SystemExit(f"{name} feature order differs from the frozen 14-feature contract: {features}")
    x = np.asarray(payload["x"], dtype="float32")
    if x.ndim != 2 or x.shape[1] != len(EXPECTED_FEATURES):
        raise SystemExit(f"{name} holdout matrix has unexpected shape {x.shape}")
    recomputed = score_xgboost(model_path, x)
    stored = np.asarray(payload["score_xgboost"], dtype="float32")
    score_delta = float(np.max(np.abs(recomputed - stored)))
    if score_delta > 2.0e-6:
        raise SystemExit(f"{name} stored/recomputed XGBoost score mismatch: {score_delta}")
    metadata = json.loads(metadata_path.read_text())
    runtime = runtime_parity(tmva_path, x, recomputed, runtime_rows)
    return {
        "name": name,
        "payload": payload,
        "features": features,
        "x": x,
        "labels": np.asarray(payload["is_signal"], dtype="int8"),
        "weights": np.asarray(payload["training_weight"], dtype="float64"),
        "scores": recomputed,
        "keys": candidate_keys(payload, holdout_path),
        "metadata": metadata,
        "stored_score_max_abs_difference": score_delta,
        "runtime_parity": runtime,
        "paths": {
            "holdout": str(holdout_path),
            "xgboost": str(model_path),
            "tmva": str(tmva_path),
            "metadata": str(metadata_path),
        },
        "sha256": {
            "holdout": sha256(holdout_path),
            "xgboost": sha256(model_path),
            "tmva": sha256(tmva_path),
            "metadata": sha256(metadata_path),
        },
    }


def native_metrics(lane: dict) -> dict:
    labels = lane["labels"]
    scores = lane["scores"]
    weights = lane["weights"]
    centrality = np.asarray(lane["payload"]["centrality"], dtype="float64")
    bins = []
    for low, high in CENTRALITY_BINS:
        mask = (centrality >= low) & (centrality < high)
        local_y = labels[mask]
        local_s = scores[mask]
        local_w = weights[mask]
        bins.append(
            {
                "centrality": [low, high],
                "rows": int(mask.sum()),
                "signal_rows": int(np.sum(local_y == 1)),
                "background_rows": int(np.sum(local_y == 0)),
                "weighted_auc": weighted_auc(local_y, local_s, local_w),
                "wp80": wp80_metrics(local_y, local_s, local_w),
            }
        )
    return {
        "rows": int(len(labels)),
        "signal_rows": int(np.sum(labels == 1)),
        "background_rows": int(np.sum(labels == 0)),
        "weighted_auc": weighted_auc(labels, scores, weights),
        "centrality_bins": bins,
    }


def cross_label_metrics(nominal: dict, variant: dict, min_common_rows: int = 1000) -> dict:
    common, nominal_index, variant_index = np.intersect1d(
        nominal["keys"], variant["keys"], return_indices=True
    )
    if len(common) < int(min_common_rows):
        raise SystemExit(
            "native holdout intersection is unexpectedly small: "
            f"{len(common)} < required {int(min_common_rows)}"
        )
    x_nom = nominal["x"][nominal_index]
    x_var = variant["x"][variant_index]
    max_feature_difference = float(np.max(np.abs(x_nom - x_var)))
    if max_feature_difference > 1.0e-6:
        raise SystemExit(f"common holdout feature rows disagree: max delta={max_feature_difference}")

    nominal_label = np.asarray(
        nominal["payload"]["nominal_is_signal"][nominal_index], dtype="int8"
    )
    ppg12_label_nom = np.asarray(
        nominal["payload"]["ppg12_source_role_label"][nominal_index], dtype="int8"
    )
    ppg12_label_var = np.asarray(variant["labels"][variant_index], dtype="int8")
    if not np.array_equal(ppg12_label_nom, ppg12_label_var):
        raise SystemExit("common holdout PPG12 labels disagree between lane caches")
    nominal_label_var = np.asarray(
        variant["payload"]["nominal_is_signal"][variant_index], dtype="int8"
    )
    if not np.array_equal(nominal_label, nominal_label_var):
        raise SystemExit("common holdout nominal labels disagree between lane caches")

    score_nominal = nominal["scores"][nominal_index]
    score_variant = variant["scores"][variant_index]
    matrix = {}
    for model_name, score in (("nominal_model", score_nominal), ("ppg12_variant_model", score_variant)):
        matrix[model_name] = {
            "nominal_label_auc_unweighted": weighted_auc(nominal_label, score, None),
            "ppg12_label_auc_unweighted": weighted_auc(ppg12_label_nom, score, None),
        }

    prompt = np.asarray(nominal["payload"]["truth_is_prompt"][nominal_index], dtype="int8") == 1
    truth_iso_pass = np.asarray(
        nominal["payload"]["truth_iso_pass"][nominal_index], dtype="int8"
    )
    truth_iso_et = np.asarray(
        nominal["payload"]["truth_iso_et"][nominal_index], dtype="float64"
    )
    conditional = {}
    from scipy.stats import spearmanr

    for model_name, score in (("nominal_model", score_nominal), ("ppg12_variant_model", score_variant)):
        valid = prompt & np.isfinite(truth_iso_et) & np.isfinite(score)
        isolated = valid & (truth_iso_pass == 1)
        nonisolated = valid & (truth_iso_pass == 0)
        correlation = math.nan
        if (
            int(valid.sum()) >= 3
            and float(np.nanstd(score[valid])) > 0.0
            and float(np.nanstd(truth_iso_et[valid])) > 0.0
        ):
            correlation = float(spearmanr(score[valid], truth_iso_et[valid]).statistic)
        conditional[model_name] = {
            "prompt_rows": int(valid.sum()),
            "isolated_prompt_rows": int(isolated.sum()),
            "nonisolated_prompt_rows": int(nonisolated.sum()),
            "mean_score_isolated_prompt": finite_mean(score[isolated]),
            "mean_score_nonisolated_prompt": finite_mean(score[nonisolated]),
            "score_auc_for_truth_isolation_within_prompt": weighted_auc(
                truth_iso_pass[valid], score[valid], None
            ),
            "spearman_score_vs_truth_iso_et_within_prompt": float(correlation),
        }
    return {
        "rows": int(len(common)),
        "fraction_of_nominal_native_holdout": float(len(common) / len(nominal["keys"])),
        "fraction_of_variant_native_holdout": float(len(common) / len(variant["keys"])),
        "max_abs_feature_difference": max_feature_difference,
        "auc_matrix": matrix,
        "conditional_truth_isolation": conditional,
    }


def draw_weighted_density(ax, scores, weights, color, label, fill_alpha) -> None:
    counts, _ = np.histogram(scores, bins=SCORE_EDGES, weights=weights)
    total = float(np.sum(counts))
    density = counts / total if total > 0.0 else counts
    y = np.r_[density, density[-1]]
    ax.fill_between(SCORE_EDGES, y, step="post", color=color, alpha=fill_alpha, linewidth=0)
    ax.step(SCORE_EDGES, y, where="post", color=color, lw=2.0, label=label)


def render_target(nominal: dict, variant: dict, metrics: dict, png: Path, pdf: Path) -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.size": 10.5,
            "axes.linewidth": 1.0,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )
    fig, axes = plt.subplots(2, 3, figsize=(15.6, 9.2), sharex=True, sharey=False)
    lanes = [
        (
            nominal,
            "Current Au+Au label",
            r"signal = prompt $\cap$ truth-isolated; all other eligible rows = background",
        ),
        (
            variant,
            "PPG12-equivalent source-role label",
            "photon sources: prompt signal; jet sources: non-prompt background; cross-role rows discarded",
        ),
    ]
    for row, (lane, row_title, definition) in enumerate(lanes):
        labels = lane["labels"]
        scores = lane["scores"]
        weights = lane["weights"]
        centrality = np.asarray(lane["payload"]["centrality"], dtype="float64")
        row_metrics = metrics["native"][lane["name"]]["centrality_bins"]
        for col, ((low, high), local_metrics) in enumerate(zip(CENTRALITY_BINS, row_metrics)):
            ax = axes[row, col]
            mask = (centrality >= low) & (centrality < high)
            signal = mask & (labels == 1)
            background = mask & (labels == 0)
            draw_weighted_density(ax, scores[background], weights[background], ORANGE, "Background", 0.11)
            draw_weighted_density(ax, scores[signal], weights[signal], GREEN, "Signal", 0.13)
            threshold = float(local_metrics["wp80"]["threshold"])
            if math.isfinite(threshold):
                ax.axvline(threshold, color=BLUE, lw=1.7, ls="--", alpha=0.9)
            ax.set_xlim(0.0, 1.0)
            ax.set_ylim(bottom=0.0)
            ax.grid(axis="y", color="0.90", lw=0.7)
            ax.set_title(f"{int(low)}–{int(high)}% centrality", fontsize=12.5, weight="bold", pad=7)
            ax.text(
                0.035,
                0.955,
                f"weighted AUC = {local_metrics['weighted_auc']:.4f}\n"
                f"WP80 = {threshold:.3f}   bkg acc. = {100.0 * local_metrics['wp80']['background_acceptance']:.1f}%\n"
                f"holdout rows: S = {local_metrics['signal_rows']:,}, B = {local_metrics['background_rows']:,}",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=9.2,
                color=INK,
                bbox={"boxstyle": "round,pad=0.28", "facecolor": "white", "edgecolor": "0.80", "alpha": 0.92},
            )
            if col == 0:
                ax.set_ylabel("Weighted fraction / 0.02", fontsize=11.2)
            if row == 1:
                ax.set_xlabel("BDT score", fontsize=11.5)
            if row == 0 and col == 2:
                ax.legend(loc="center right", frameon=False, fontsize=9.6)

    fig.text(0.055, 0.982, r"$\bf{sPHENIX}$ Simulation", fontsize=14.5, color=INK, va="top")
    fig.text(
        0.055,
        0.947,
        "Au+Au photon-ID BDT: isolated-prompt versus PPG12 source-role class construction",
        fontsize=18.2,
        weight="bold",
        color=INK,
        va="top",
    )
    fig.text(
        0.055,
        0.912,
        r"Exact candidate-row holdouts; $15\leq E_T^\gamma<35$ GeV, $|\eta|<0.7$, $|z_{\rm vtx}|<10$ cm; same 14 inputs, weights, split seed, and XGBoost configuration",
        fontsize=11.2,
        color=MUTED,
        va="top",
    )
    row_header_y = (0.858, 0.452)
    row_colors = (GREEN, BLUE)
    for y, color, (_lane, row_title, definition) in zip(row_header_y, row_colors, lanes):
        fig.text(0.055, y, row_title, fontsize=11.7, weight="bold", color=color, va="center")
        fig.text(0.335, y, definition, fontsize=9.8, color=MUTED, va="center")
    cross = metrics["cross_label"]["auc_matrix"]
    fig.text(
        0.055,
        0.025,
        "Cross-label AUC on the common native-holdout intersection:  "
        f"nominal model = {cross['nominal_model']['nominal_label_auc_unweighted']:.4f} (nominal), "
        f"{cross['nominal_model']['ppg12_label_auc_unweighted']:.4f} (PPG12);  "
        f"variant model = {cross['ppg12_variant_model']['nominal_label_auc_unweighted']:.4f} (nominal), "
        f"{cross['ppg12_variant_model']['ppg12_label_auc_unweighted']:.4f} (PPG12).  "
        "No data application or model promotion.",
        fontsize=9.8,
        color=MUTED,
    )
    fig.subplots_adjust(left=0.075, right=0.985, top=0.815, bottom=0.105, hspace=0.50, wspace=0.18)
    png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png, dpi=240, facecolor="white")
    fig.savefig(pdf, facecolor="white")
    plt.close(fig)


def compare_frozen_configuration(nominal: dict, variant: dict, *, strict_production: bool) -> dict:
    keys = [
        "features",
        "pt_range",
        "cent_range",
        "weight_mode",
        "majority_cap_ratio",
        "majority_class_optimization",
        "background_subsampling",
        "xgboost",
    ]
    differences = {}
    for key in keys:
        left = nominal["metadata"].get(key)
        right = variant["metadata"].get(key)
        if left != right:
            differences[key] = {"nominal": left, "variant": right}
    nominal_split = nominal["metadata"].get("split") or {}
    variant_split = variant["metadata"].get("split") or {}
    split_contract = {
        "mode": (nominal_split.get("mode"), variant_split.get("mode")),
        "random_seed": (nominal_split.get("random_seed"), variant_split.get("random_seed")),
        "stratified_by_class": (
            nominal_split.get("stratified_by_class"),
            variant_split.get("stratified_by_class"),
        ),
        "test_fraction_requested": (
            nominal_split.get("test_fraction_requested"),
            variant_split.get("test_fraction_requested"),
        ),
    }
    for field, pair in split_contract.items():
        if pair[0] != pair[1]:
            differences[f"split.{field}"] = {"nominal": pair[0], "variant": pair[1]}
    # Row counts, labels, label audit, and input cache paths are expected to
    # differ.  Every model/phase-space/hyperparameter field above must match.
    if differences:
        raise SystemExit(f"frozen training configuration differs between lanes: {differences}")

    expected_failures = {}
    if strict_production:
        nominal_input_cache = nominal["metadata"].get("cache_input_provenance") or {}
        nominal_final_cache = nominal["metadata"].get("cache_final_provenance") or {}
        variant_input_cache = variant["metadata"].get("cache_input_provenance") or {}
        variant_final_cache = variant["metadata"].get("cache_final_provenance") or {}
        cache_copy_checks = {
            "nominal input/final cache SHA": (
                nominal_input_cache.get("sha256"),
                nominal_final_cache.get("sha256"),
            ),
            "nominal final/variant input cache SHA": (
                nominal_final_cache.get("sha256"),
                variant_input_cache.get("sha256"),
            ),
            "nominal final/variant input cache size": (
                nominal_final_cache.get("size_bytes"),
                variant_input_cache.get("size_bytes"),
            ),
        }
        for field, (left, right) in cache_copy_checks.items():
            if left is None or right is None or left != right:
                expected_failures[field] = {"left": left, "right": right}
        if nominal_final_cache.get("path") == variant_input_cache.get("path"):
            expected_failures["separate copied cache paths"] = {
                "nominal": nominal_final_cache.get("path"),
                "variant": variant_input_cache.get("path"),
            }
        if variant_input_cache.get("sha256") == variant_final_cache.get("sha256"):
            expected_failures["variant input/final cache transformation"] = {
                "input": variant_input_cache.get("sha256"),
                "final": variant_final_cache.get("sha256"),
                "expected": "different after cross-role row removal and weight recomputation",
            }
        lane_expectations = {
            "nominal": "nominal-isolated-prompt",
            "ppg12_variant": "ppg12-source-role",
        }
        for lane in (nominal, variant):
            name = lane["name"]
            metadata = lane["metadata"]
            checks = {
                "status": (metadata.get("status"), "trained"),
                "features": (metadata.get("features"), EXPECTED_FEATURES),
                "pt_range": (metadata.get("pt_range"), [15.0, 35.0]),
                "cent_range": (metadata.get("cent_range"), None),
                "weight_mode": (metadata.get("weight_mode"), "ppg12-exact"),
                "majority_cap_ratio": (metadata.get("majority_cap_ratio"), 0.0),
                "majority_class_optimization.enabled": (
                    metadata.get("majority_class_optimization", {}).get("enabled"),
                    False,
                ),
                "background_subsampling.enabled": (
                    metadata.get("background_subsampling", {}).get("enabled"),
                    False,
                ),
                "weighting.weight_mode": (
                    metadata.get("weighting", {}).get("weight_mode"),
                    "ppg12-exact",
                ),
                "weighting.source": (
                    metadata.get("weighting", {}).get("source"),
                    "__ppg12_exact_training_weight",
                ),
                "weighting.weights_computed_before_binning": (
                    metadata.get("weighting", {}).get("weights_computed_before_binning"),
                    True,
                ),
                "weighting.event_weight_used": (
                    metadata.get("weighting", {}).get("event_weight_used"),
                    False,
                ),
                "weighting.cross_section_weight_used_for_training": (
                    metadata.get("weighting", {}).get("cross_section_weight_used_for_training"),
                    False,
                ),
                "weighting.centrality_event_weight": (
                    metadata.get("weighting", {}).get("centrality_event_weight"),
                    False,
                ),
                "weighting.vertex_reweight": (
                    metadata.get("weighting", {}).get("vertex_reweight"),
                    False,
                ),
                "split.mode": (metadata.get("split", {}).get("mode"), "row"),
                "split.random_seed": (metadata.get("split", {}).get("random_seed"), 13),
                "split.stratified_by_class": (
                    metadata.get("split", {}).get("stratified_by_class"),
                    True,
                ),
                "split.test_fraction_requested": (
                    metadata.get("split", {}).get("test_fraction_requested"),
                    0.10,
                ),
                "label_contract.contract": (
                    metadata.get("label_contract", {}).get("contract"),
                    lane_expectations[name],
                ),
                "label_contract.nominal_label_closure_mismatches": (
                    metadata.get("label_contract", {}).get("nominal_label_closure_mismatches"),
                    0,
                ),
                "label_contract.ppg12_label_closure_mismatches": (
                    metadata.get("label_contract", {}).get("ppg12_label_closure_mismatches"),
                    0,
                ),
                "label_contract.source_role_closure_mismatches": (
                    metadata.get("label_contract", {}).get("source_role_closure_mismatches"),
                    0,
                ),
                "label_contract.truth_isolation_used_by_label": (
                    metadata.get("label_contract", {}).get("truth_isolation_used_by_label"),
                    name == "nominal",
                ),
                "label_contract.source_role_used_by_label": (
                    metadata.get("label_contract", {}).get("source_role_used_by_label"),
                    name == "ppg12_variant",
                ),
            }
            if name == "nominal":
                checks["label_contract.rows_discarded"] = (
                    metadata.get("label_contract", {}).get("rows_discarded"),
                    0,
                )
                checks["label_contract.discarded_stale_precomputed_weight"] = (
                    metadata.get("label_contract", {}).get("discarded_stale_precomputed_weight"),
                    False,
                )
            else:
                discarded = metadata.get("label_contract", {}).get("rows_discarded")
                if not isinstance(discarded, int) or discarded <= 0:
                    expected_failures[
                        f"{name}.label_contract.rows_discarded"
                    ] = {"observed": discarded, "expected": "> 0"}
                checks["label_contract.discarded_stale_precomputed_weight"] = (
                    metadata.get("label_contract", {}).get("discarded_stale_precomputed_weight"),
                    True,
                )
            for field, (observed, expected) in checks.items():
                if observed != expected:
                    expected_failures[f"{name}.{field}"] = {
                        "observed": observed,
                        "expected": expected,
                    }
            xgboost = metadata.get("xgboost") or {}
            for field, expected in EXPECTED_XGBOOST.items():
                observed = xgboost.get(field)
                if observed != expected:
                    expected_failures[f"{name}.xgboost.{field}"] = {
                        "observed": observed,
                        "expected": expected,
                    }
        if expected_failures:
            raise SystemExit(
                "paired models do not match the frozen THE-107 production contract: "
                f"{expected_failures}"
            )
    return {
        "passed": True,
        "strict_production_contract": bool(strict_production),
        "checked_fields": keys
        + [
            "split.mode",
            "split.random_seed",
            "split.stratified_by_class",
            "split.test_fraction_requested",
            "label_contract.contract",
            "label_contract closure and lane-semantics fields",
            "majority/downsampling disabled",
            "global class/eta/ET weighting without event/cross-section/centrality/vertex weights",
            "byte-identical nominal-final to variant-input cache copy on separate paths",
            "variant cache transformed exactly once before WP derivation",
        ]
        + [f"xgboost.{field}" for field in EXPECTED_XGBOOST],
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    for lane in ("nominal", "variant"):
        parser.add_argument(f"--{lane}-holdout", type=Path, required=True)
        parser.add_argument(f"--{lane}-xgboost", type=Path, required=True)
        parser.add_argument(f"--{lane}-tmva", type=Path, required=True)
        parser.add_argument(f"--{lane}-metadata", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument(
        "--campaign-label",
        default="THE-107",
        help="Campaign identifier recorded in the validation payload.",
    )
    parser.add_argument(
        "--artifact-prefix",
        default="the107_auau_bdt_label_contract",
        help="Filename prefix for the PNG, PDF, and validation JSON outputs.",
    )
    parser.add_argument("--runtime-parity-rows", type=int, default=5000)
    parser.add_argument(
        "--allow-nonproduction-config",
        action="store_true",
        help="Retain pair-equality checks but skip the exact 15--35 GeV/450-tree production contract (canaries only).",
    )
    parser.add_argument(
        "--min-common-holdout-rows",
        type=int,
        default=1000,
        help=(
            "Minimum candidate-identity intersection required for cross-label metrics. "
            "The production default is 1000; use a lower value only for extraction/training canaries."
        ),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    nominal = load_lane(
        "nominal",
        args.nominal_holdout,
        args.nominal_xgboost,
        args.nominal_tmva,
        args.nominal_metadata,
        args.runtime_parity_rows,
    )
    variant = load_lane(
        "ppg12_variant",
        args.variant_holdout,
        args.variant_xgboost,
        args.variant_tmva,
        args.variant_metadata,
        args.runtime_parity_rows,
    )
    metrics = {
        "schema": "THE107_AUAU_BDT_LABEL_CONTRACT_VALIDATION_V1",
        "campaign": args.campaign_label,
        "status": "VALIDATED_SIMULATION_ONLY_NO_PROMOTION",
        "frozen_configuration": compare_frozen_configuration(
            nominal,
            variant,
            strict_production=not args.allow_nonproduction_config,
        ),
        "native": {
            nominal["name"]: native_metrics(nominal),
            variant["name"]: native_metrics(variant),
        },
        "cross_label": cross_label_metrics(
            nominal, variant, min_common_rows=args.min_common_holdout_rows
        ),
        "runtime_parity": {
            nominal["name"]: nominal["runtime_parity"],
            variant["name"]: variant["runtime_parity"],
        },
        "stored_score_reproduction": {
            nominal["name"]: nominal["stored_score_max_abs_difference"],
            variant["name"]: variant["stored_score_max_abs_difference"],
        },
        "inputs": {
            nominal["name"]: {"paths": nominal["paths"], "sha256": nominal["sha256"]},
            variant["name"]: {"paths": variant["paths"], "sha256": variant["sha256"]},
        },
        "boundaries": [
            "candidate-row holdout, not event-grouped",
            "simulation-only model comparison",
            "no data scoring",
            "no canonical model promotion",
        ],
    }
    args.outdir.mkdir(parents=True, exist_ok=True)
    png = args.outdir / f"{args.artifact_prefix}_score_separation_2x3.png"
    pdf = args.outdir / f"{args.artifact_prefix}_score_separation_2x3.pdf"
    render_target(nominal, variant, metrics, png, pdf)
    metrics["outputs"] = {
        "png": str(png),
        "png_sha256": sha256(png),
        "pdf": str(pdf),
        "pdf_sha256": sha256(pdf),
    }
    output = args.outdir / f"{args.artifact_prefix}_validation.json"
    output.write_text(json.dumps(json_safe(metrics), indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": metrics["status"], "png": str(png), "metrics": str(output)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
