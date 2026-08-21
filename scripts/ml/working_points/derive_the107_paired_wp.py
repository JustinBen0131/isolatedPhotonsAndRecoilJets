#!/usr/bin/env python3
"""Derive source-verifiable WP surfaces for the two THE-107 label lanes.

The numerical threshold construction intentionally reuses the accepted
full-matrix weighted-quantile implementation.  This wrapper replaces its
historical campaign labels, records the exact label contract and input hashes,
and fits the seven per-centrality constants with a continuous linear runtime
surface for WP90, WP80, and WP70.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from derive_the57_full_weighted_wp80 import FEATURES, derive_payload


LANE_LABELS = {
    "nominal-isolated-prompt": {
        "source_label": "Au+Au candidate-level isolated-prompt label",
        "model_label": "current Au+Au label construction control",
        "truth_isolation_used_by_label": True,
        "source_role_used_by_label": False,
    },
    "ppg12-source-role": {
        "source_label": "PPG12-equivalent source-role prompt/background label",
        "model_label": "PPG12 source-role label-construction variant",
        "truth_isolation_used_by_label": False,
        "source_role_used_by_label": True,
    },
}

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


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def linear_fit(flat_rows: list[dict], target: float) -> dict:
    selected = [
        row
        for row in flat_rows
        if abs(float(row["target_signal_efficiency"]) - float(target)) < 1.0e-9
        and math.isfinite(float(row["threshold"]))
    ]
    if len(selected) < 2:
        raise SystemExit(f"insufficient centrality constants for WP{int(100 * target)}")
    x = np.asarray([row["centrality_center"] for row in selected], dtype="float64")
    y = np.asarray([row["threshold"] for row in selected], dtype="float64")
    err = np.asarray(
        [row.get("threshold_stat_err", math.nan) for row in selected], dtype="float64"
    )
    good = np.isfinite(err) & (err > 0.0)
    weights = np.ones_like(err, dtype="float64")
    weights[good] = 1.0 / np.maximum(err[good], 1.0e-8)
    slope, intercept = np.polyfit(x, y, 1, w=weights)
    predicted = intercept + slope * x
    residual = y - predicted
    return {
        "target_signal_efficiency": float(target),
        "wp_label": f"WP{int(round(100 * target))}",
        "intercept": float(intercept),
        "slope_per_centrality_percentile": float(slope),
        "formula": f"T(c) = {intercept:.10f} + {slope:.10f} c",
        "fit_mode": "inverse-threshold-error weighted linear fit",
        "n_centrality_points": int(len(x)),
        "rms_residual": float(np.sqrt(np.mean(residual * residual))),
        "max_abs_residual": float(np.max(np.abs(residual))),
        "centrality_centers": x.tolist(),
        "thresholds": y.tolist(),
        "residuals": residual.tolist(),
    }


def same_path(left: str | Path, right: str | Path) -> bool:
    return Path(left).expanduser().resolve() == Path(right).expanduser().resolve()


def verify_metadata(args: argparse.Namespace) -> dict:
    metadata = json.loads(args.metadata.read_text())
    failures = {}
    matrix_sha256 = sha256(args.matrix)
    checks = {
        "status": (metadata.get("status"), "trained"),
        "features": (metadata.get("features"), FEATURES),
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
        "label_contract.contract": (
            metadata.get("label_contract", {}).get("contract"),
            args.lane,
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
            args.lane == "nominal-isolated-prompt",
        ),
        "label_contract.source_role_used_by_label": (
            metadata.get("label_contract", {}).get("source_role_used_by_label"),
            args.lane == "ppg12-source-role",
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
    }
    if args.lane == "nominal-isolated-prompt":
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
            failures["label_contract.rows_discarded"] = {
                "observed": discarded,
                "expected": "> 0",
            }
        checks["label_contract.discarded_stale_precomputed_weight"] = (
            metadata.get("label_contract", {}).get("discarded_stale_precomputed_weight"),
            True,
        )
    for field, (observed, expected) in checks.items():
        if observed != expected:
            failures[field] = {"observed": observed, "expected": expected}
    for field, expected in EXPECTED_XGBOOST.items():
        observed = (metadata.get("xgboost") or {}).get(field)
        if observed != expected:
            failures[f"xgboost.{field}"] = {"observed": observed, "expected": expected}
    if not same_path(metadata.get("cache_file", ""), args.matrix):
        failures["cache_file"] = {
            "observed": metadata.get("cache_file"),
            "expected": str(args.matrix),
        }
    if not same_path(metadata.get("output_xgb_json", ""), args.model):
        failures["output_xgb_json"] = {
            "observed": metadata.get("output_xgb_json"),
            "expected": str(args.model),
        }
    final_cache = metadata.get("cache_final_provenance") or {}
    if final_cache.get("sha256") != matrix_sha256:
        failures["cache_final_provenance.sha256"] = {
            "observed": final_cache.get("sha256"),
            "expected": matrix_sha256,
        }
    if final_cache.get("size_bytes") != args.matrix.stat().st_size:
        failures["cache_final_provenance.size_bytes"] = {
            "observed": final_cache.get("size_bytes"),
            "expected": args.matrix.stat().st_size,
        }
    if failures:
        raise SystemExit(f"WP input metadata does not match the frozen THE-107 lane: {failures}")
    return {
        "status": metadata.get("status"),
        "model_id": metadata.get("model_id"),
        "product": metadata.get("product"),
        "label_contract": metadata.get("label_contract"),
        "features": metadata.get("features"),
        "pt_range": metadata.get("pt_range"),
        "cent_range": metadata.get("cent_range"),
        "weight_mode": metadata.get("weight_mode"),
        "majority_cap_ratio": metadata.get("majority_cap_ratio"),
        "majority_class_optimization": metadata.get("majority_class_optimization"),
        "background_subsampling": metadata.get("background_subsampling"),
        "weighting": metadata.get("weighting"),
        "cache_input_provenance": metadata.get("cache_input_provenance"),
        "cache_final_provenance": final_cache,
        "split": metadata.get("split"),
        "xgboost": metadata.get("xgboost"),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--lane", choices=sorted(LANE_LABELS), required=True)
    parser.add_argument("--json-out", type=Path, required=True)
    parser.add_argument(
        "--campaign-label",
        default="THE-107",
        help="Campaign identifier recorded in the working-point payload.",
    )
    parser.add_argument("--chunk-size", type=int, default=500_000)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not args.matrix.is_file() or not args.model.is_file() or not args.metadata.is_file():
        raise SystemExit(
            f"missing matrix/model/metadata: {args.matrix} ; {args.model} ; {args.metadata}"
        )
    verified_metadata = verify_metadata(args)
    payload = derive_payload(args)
    lane = LANE_LABELS[args.lane]
    payload.update(
        {
            "schema": "THE107_AUAU_BDT_PAIRED_LABEL_WP_V1",
            "campaign": args.campaign_label,
            "status": "SIMULATION_VALIDATION_ONLY_NO_PROMOTION",
            "label_contract": args.lane,
            "source_label": lane["source_label"],
            "model_label": lane["model_label"],
            "truth_isolation_used_by_label": lane["truth_isolation_used_by_label"],
            "source_role_used_by_label": lane["source_role_used_by_label"],
            "training_inputs": (
                "frozen 14-feature Au+Au baseline; 15 <= cluster E_T < 35 GeV; "
                "PPG12-exact global class/eta/E_T weights"
            ),
            "input_sha256": {
                "matrix": verified_metadata["cache_final_provenance"]["sha256"],
                "xgboost_model": sha256(args.model),
                "metadata": sha256(args.metadata),
            },
            "verified_model_metadata": verified_metadata,
            "continuous_fits": [linear_fit(payload["flat_rows"], target) for target in payload["targets"]],
            "boundaries": [
                "derived from simulation training population",
                "not a detector-level closure result",
                "no data application",
                "no canonical-model promotion",
            ],
        }
    )
    args.json_out.parent.mkdir(parents=True, exist_ok=True)
    args.json_out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(args.json_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
