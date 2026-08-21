#!/usr/bin/env python3
"""Derive exact weighted WP70/WP80/WP90 for a THE-134 factorial view.

The exact binned thresholds are always authoritative.  A compact linear
surface is emitted only as an optional runtime candidate and is marked
accepted only when both residual and achieved-efficiency gates pass.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve()
_CONTRACTS = _HERE.parents[1] / "contracts"
sys.path.insert(0, str(_CONTRACTS))
from the134_h70_contract import (  # noqa: E402
    ALL_SHOWER_VIEWS,
    AUAU_CENTRALITY_EDGES,
    FEATURES_BY_SYSTEM,
    MAX_SURFACE_ABS_RESIDUAL,
    MAX_SURFACE_EFFICIENCY_ERROR,
    MAX_SURFACE_RMS,
    MAX_WP_EXACT_EFFICIENCY_ERROR,
    MODEL_DOMAIN_GEV,
    READER_ET_EDGES,
    TARGET_EFFICIENCIES,
    VIEW_NAME,
    expected_model_origin,
    expected_model_split,
    is_sha256,
    model_extraction_authority_checks,
    sha256_file,
    shower_semantic_sha256,
    valid_reuse_pinned_hashes,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--system", choices=sorted(FEATURES_BY_SYSTEM), required=True)
    parser.add_argument("--view", choices=ALL_SHOWER_VIEWS, default=VIEW_NAME)
    parser.add_argument("--score-sample", type=Path, required=True)
    parser.add_argument("--population-certificate", type=Path, required=True)
    parser.add_argument("--model-metadata", type=Path, required=True)
    parser.add_argument("--json-out", type=Path, required=True)
    parser.add_argument("--surface-max-rms", type=float, default=MAX_SURFACE_RMS)
    parser.add_argument(
        "--surface-max-residual", type=float, default=MAX_SURFACE_ABS_RESIDUAL
    )
    parser.add_argument(
        "--surface-max-efficiency-error",
        type=float,
        default=MAX_SURFACE_EFFICIENCY_ERROR,
    )
    parser.add_argument(
        "--exact-max-efficiency-error",
        type=float,
        default=MAX_WP_EXACT_EFFICIENCY_ERROR,
    )
    return parser.parse_args()


def enforce_frozen_gate_maxima(args: argparse.Namespace) -> None:
    gates = {
        "surface_max_rms": (args.surface_max_rms, MAX_SURFACE_RMS),
        "surface_max_residual": (
            args.surface_max_residual,
            MAX_SURFACE_ABS_RESIDUAL,
        ),
        "surface_max_efficiency_error": (
            args.surface_max_efficiency_error,
            MAX_SURFACE_EFFICIENCY_ERROR,
        ),
        "exact_max_efficiency_error": (
            args.exact_max_efficiency_error,
            MAX_WP_EXACT_EFFICIENCY_ERROR,
        ),
    }
    invalid = {
        name: {"requested": value, "frozen_maximum": maximum}
        for name, (value, maximum) in gates.items()
        if not math.isfinite(value) or value < 0.0 or value > maximum
    }
    if invalid:
        raise SystemExit("scientific gate widening is forbidden: " + json.dumps(invalid, sort_keys=True))


def weighted_threshold(scores: np.ndarray, weights: np.ndarray, efficiency: float) -> float:
    scores = np.asarray(scores, dtype=np.float64)
    weights = np.asarray(weights, dtype=np.float64)
    good = np.isfinite(scores) & np.isfinite(weights) & (weights > 0.0)
    if not np.any(good):
        return math.nan
    scores = scores[good]
    weights = weights[good]
    order = np.argsort(scores, kind="mergesort")
    scores = scores[order]
    weights = weights[order]
    cumulative = np.cumsum(weights)
    target = (1.0 - efficiency) * float(cumulative[-1])
    index = int(np.searchsorted(cumulative, target, side="left"))
    return float(scores[min(index, len(scores) - 1)])


def weighted_efficiency(mask: np.ndarray, weights: np.ndarray) -> float:
    weights = np.asarray(weights, dtype=np.float64)
    denominator = float(np.sum(weights))
    return float(np.sum(weights[np.asarray(mask, dtype=bool)]) / denominator) if denominator > 0 else math.nan


def require_model_metadata(metadata: dict, system: str, view_name: str = VIEW_NAME) -> None:
    expected_features = list(FEATURES_BY_SYSTEM[system])
    nested_contract = metadata.get("the134_view_contract", {})
    nested_view = nested_contract.get("shower_definition")
    nested_semantic = nested_contract.get("shower_semantic_sha256")
    top_view = metadata.get("shower_definition")
    top_semantic = metadata.get("shower_semantic_sha256")
    observed_view = top_view or nested_view
    observed_semantic = top_semantic or nested_semantic
    weighting = metadata.get("weighting", {})
    expected_origin = expected_model_origin(system, view_name)
    expected_split = expected_model_split(system, view_name)
    reuse_pins = metadata.get("reuse_pinned_hashes", {})
    split = metadata.get("split", {})
    checks = {
        "features": (metadata.get("features"), expected_features),
        "pt_range": (metadata.get("pt_range"), list(MODEL_DOMAIN_GEV)),
        "weight_mode": (metadata.get("weight_mode"), "ppg12-exact"),
        "shower_definition": (observed_view, view_name),
        "shower_semantic_sha256": (
            observed_semantic,
            shower_semantic_sha256(view_name),
        ),
        "split.mode": (
            metadata.get("split", {}).get("mode"),
            expected_split["mode"],
        ),
        "split.test_fraction": (
            metadata.get("split", {}).get("test_fraction_requested"),
            expected_split["test_fraction_requested"],
        ),
        "split.seed": (
            metadata.get("split", {}).get("random_seed"),
            expected_split["random_seed"],
        ),
        "split.positive_rows": (
            int(split.get("train_rows", 0)) > 0 and int(split.get("test_rows", 0)) > 0,
            True,
        ),
        "split.positive_events": (
            (
                int(split.get("train_events", 0)) > 0
                and int(split.get("test_events", 0)) > 0
            )
            if expected_split["mode"] == "event50"
            else True,
            True,
        ),
        "split.accepted_legacy_boundary": (
            metadata.get("accepted_split_contract") == expected_split
            if expected_split["mode"] == "row"
            else True,
            True,
        ),
        "event_weight_unused": (weighting.get("event_weight_used"), False),
        "vertex_weight_unused": (weighting.get("vertex_reweight"), False),
        "centrality_weight_unused": (weighting.get("centrality_event_weight"), False),
        "cross_section_weight_unused": (
            weighting.get("cross_section_weight_used_for_training"),
            False,
        ),
        "weights_computed_before_binning": (
            weighting.get("weights_computed_before_binning"),
            True,
        ),
        "old_low_calo_veto_disabled": (
            bool(metadata.get("event_quality_filter", {}).get("enabled", False)),
            False,
        ),
        "model_origin": (metadata.get("model_origin"), expected_origin),
        "reuse_pins": (
            valid_reuse_pinned_hashes(reuse_pins, system, view_name)
            if expected_origin.startswith("REUSED_")
            else not reuse_pins,
            True,
        ),
        "full_extraction_authority": (
            all(model_extraction_authority_checks(metadata).values()),
            True,
        ),
        "view_contract_schema": (
            nested_contract.get("schema"),
            "THE134_FACTORIAL_VIEW_MODEL_VIEW_CONTRACT_V1",
        ),
        "view_contract_system": (nested_contract.get("system"), system),
        "top_nested_view_consistency": (top_view in {None, nested_view}, True),
        "top_nested_semantic_consistency": (
            top_semantic in {None, nested_semantic},
            True,
        ),
    }
    failures = {
        name: {"observed": observed, "expected": expected}
        for name, (observed, expected) in checks.items()
        if observed != expected
    }
    if failures:
        raise SystemExit("model metadata violates THE-134 view contract: " + json.dumps(failures, sort_keys=True))


def derive_rows(
    scores: np.ndarray,
    labels: np.ndarray,
    weights: np.ndarray,
    axis: np.ndarray,
    edges: np.ndarray,
) -> list[dict]:
    rows: list[dict] = []
    for efficiency in TARGET_EFFICIENCIES:
        for lo, hi in zip(edges[:-1], edges[1:]):
            local = np.isfinite(axis) & (axis >= lo) & (axis < hi)
            signal = local & (labels == 1)
            background = local & (labels == 0)
            threshold = weighted_threshold(scores[signal], weights[signal], efficiency)
            achieved_signal_efficiency = weighted_efficiency(
                scores[signal] > threshold, weights[signal]
            )
            tied_signal_weight_fraction = weighted_efficiency(
                scores[signal] == threshold, weights[signal]
            )
            rows.append(
                {
                    "wp": f"WP{int(round(100 * efficiency))}",
                    "target_signal_efficiency": float(efficiency),
                    "bin_lo": float(lo),
                    "bin_hi": float(hi),
                    "bin_center": float(0.5 * (lo + hi)),
                    "threshold": threshold,
                    "achieved_signal_efficiency": achieved_signal_efficiency,
                    "achieved_minus_target": achieved_signal_efficiency - efficiency,
                    "abs_efficiency_error": abs(achieved_signal_efficiency - efficiency),
                    "signal_weight_fraction_tied_at_threshold": tied_signal_weight_fraction,
                    "background_acceptance": weighted_efficiency(
                        scores[background] > threshold, weights[background]
                    ),
                    "signal_rows": int(np.sum(signal)),
                    "background_rows": int(np.sum(background)),
                }
            )
    return rows


def fit_surfaces(
    rows: list[dict],
    scores: np.ndarray,
    labels: np.ndarray,
    weights: np.ndarray,
    axis: np.ndarray,
    edges: np.ndarray,
    args: argparse.Namespace,
) -> dict:
    surfaces = {}
    for efficiency in TARGET_EFFICIENCIES:
        label = f"WP{int(round(100 * efficiency))}"
        selected = [row for row in rows if row["wp"] == label]
        centers = np.asarray([row["bin_center"] for row in selected], dtype=np.float64)
        thresholds = np.asarray([row["threshold"] for row in selected], dtype=np.float64)
        finite = np.isfinite(thresholds)
        if int(np.sum(finite)) < 2:
            surfaces[label] = {"status": "BINNED_THRESHOLDS_ONLY", "reason": "insufficient finite bins"}
            continue
        slope, intercept = np.polyfit(centers[finite], thresholds[finite], 1)
        residual = thresholds[finite] - (intercept + slope * centers[finite])
        applied = intercept + slope * axis
        achieved: list[float] = []
        for lo, hi in zip(edges[:-1], edges[1:]):
            local = (labels == 1) & np.isfinite(axis) & (axis >= lo) & (axis < hi)
            achieved.append(weighted_efficiency(scores[local] > applied[local], weights[local]))
        finite_achieved = np.asarray([value for value in achieved if math.isfinite(value)], dtype=np.float64)
        rms = float(np.sqrt(np.mean(np.square(residual))))
        max_residual = float(np.max(np.abs(residual)))
        max_eff_error = (
            float(np.max(np.abs(finite_achieved - efficiency))) if len(finite_achieved) else math.inf
        )
        accepted = (
            rms <= args.surface_max_rms
            and max_residual <= args.surface_max_residual
            and max_eff_error <= args.surface_max_efficiency_error
        )
        surfaces[label] = {
            "status": "ACCEPTED" if accepted else "BINNED_THRESHOLDS_ONLY",
            "intercept": float(intercept),
            "slope_per_axis_unit": float(slope),
            "rms_residual": rms,
            "max_abs_residual": max_residual,
            "max_abs_efficiency_error": max_eff_error,
            "achieved_efficiency_by_bin": achieved,
            "gates": {
                "max_rms": args.surface_max_rms,
                "max_abs_residual": args.surface_max_residual,
                "max_abs_efficiency_error": args.surface_max_efficiency_error,
            },
        }
    return surfaces


def exact_efficiency_gate(rows: list[dict], tolerance: float) -> bool:
    return all(
        math.isfinite(row["abs_efficiency_error"])
        and row["abs_efficiency_error"] <= tolerance
        for row in rows
    )


def main() -> int:
    args = parse_args()
    enforce_frozen_gate_maxima(args)
    if not all(
        path.is_file()
        for path in (args.score_sample, args.population_certificate, args.model_metadata)
    ):
        raise SystemExit("missing score sample, population certificate, or model metadata")
    metadata = json.loads(args.model_metadata.read_text())
    require_model_metadata(metadata, args.system, args.view)
    population_certificate = json.loads(args.population_certificate.read_text())
    expected_certificate_schema = (
        "THE134_FACTORIAL_VIEW_EXACT_HOLDOUT_CERTIFICATE_V1"
        if args.system == "pp"
        else "THE134_FACTORIAL_VIEW_FULL_WEIGHTED_WP_SAMPLE_CERTIFICATE_V1"
    )
    population_provenance = population_certificate.get("provenance", {})
    observed_sample_sha = population_provenance.get(
        "holdout_sha256" if args.system == "pp" else "score_sample_sha256"
    )
    population_certificate_checks = {
        "schema": population_certificate.get("schema") == expected_certificate_schema,
        "status": population_certificate.get("status") == "PASS",
        "system": population_certificate.get("system") == args.system,
        "shower_definition": population_certificate.get("shower_definition") == args.view,
        "shower_semantic_sha256": population_certificate.get("shower_semantic_sha256")
        == shower_semantic_sha256(args.view),
        "feature_order": population_certificate.get("feature_order")
        == list(FEATURES_BY_SYSTEM[args.system]),
        "score_sample_sha256": observed_sample_sha == sha256_file(args.score_sample),
        "model_metadata_sha256": population_provenance.get("model_metadata_sha256")
        == sha256_file(args.model_metadata),
        "model_xgb_sha256": is_sha256(
            population_provenance.get("model_xgb_sha256")
        ),
        "model_receipt_sha256": is_sha256(
            population_provenance.get("model_receipt_sha256")
        ),
        "model_origin": population_certificate.get("model_origin")
        == metadata.get("model_origin"),
        "reuse_pinned_hashes": population_certificate.get("reuse_pinned_hashes")
        == metadata.get("reuse_pinned_hashes"),
    }
    if args.system == "pp":
        population_certificate_checks["zero_event_group_overlap"] = (
            population_certificate.get("split", {}).get("event_overlap") == 0
        )
    else:
        population_certificate_checks["full_combined_signal_population"] = (
            population_certificate.get("population")
            == "full combined eligible Photon12+20 and Jet12+20+30+40"
            and population_certificate.get("population_checks", {}).get(
                "signal_source_complete"
            )
            is True
        )
    if not all(population_certificate_checks.values()):
        raise SystemExit(
            "WP population certificate violates THE-134 contract: "
            + json.dumps(population_certificate_checks, sort_keys=True)
        )
    data = np.load(args.score_sample, allow_pickle=True)
    required = {"features", "x", "is_signal", "training_weight", "score_xgboost", "cluster_Et"}
    if args.system == "auau":
        required.add("centrality")
    missing = sorted(required - set(data.files))
    if missing:
        raise SystemExit(f"score sample missing fields: {missing}")
    observed_features = [str(item) for item in data["features"].tolist()]
    if observed_features != list(FEATURES_BY_SYSTEM[args.system]):
        raise SystemExit(f"score-sample feature order mismatch: {observed_features}")
    labels = np.asarray(data["is_signal"], dtype=np.int8)
    weights = np.asarray(data["training_weight"], dtype=np.float64)
    scores = np.asarray(data["score_xgboost"], dtype=np.float64)
    et = np.asarray(data["cluster_Et"], dtype=np.float64)
    base = (
        np.isfinite(scores)
        & np.isfinite(weights)
        & (weights > 0.0)
        & np.isin(labels, [0, 1])
        & np.isfinite(et)
        & (et >= MODEL_DOMAIN_GEV[0])
        & (et < MODEL_DOMAIN_GEV[1])
    )
    if not np.all(base):
        raise SystemExit(f"score sample contains {int(np.sum(~base))} invalid/out-of-domain rows")
    if args.system == "pp":
        axis_name = "cluster_Et"
        axis = et
        edges = np.asarray(READER_ET_EDGES, dtype=np.float64)
    else:
        axis_name = "centrality"
        axis = np.asarray(data["centrality"], dtype=np.float64)
        edges = np.asarray(AUAU_CENTRALITY_EDGES, dtype=np.float64)
    rows = derive_rows(scores, labels, weights, axis, edges)
    bin_population_gate = all(row["signal_rows"] > 0 and row["background_rows"] > 0 for row in rows)
    finite_threshold_gate = all(math.isfinite(row["threshold"]) for row in rows)
    exact_efficiency_pass = exact_efficiency_gate(rows, args.exact_max_efficiency_error)
    surfaces = fit_surfaces(rows, scores, labels, weights, axis, edges, args)
    gates = {
        "model_metadata_identity": True,
        "certified_wp_population": all(population_certificate_checks.values()),
        "holdout_feature_order": True,
        "finite_positive_weights": bool(np.all(np.isfinite(weights) & (weights > 0.0))),
        "both_classes_in_every_bin": bin_population_gate,
        "all_exact_thresholds_finite": finite_threshold_gate,
        "exact_weighted_efficiency_within_tolerance": exact_efficiency_pass,
    }
    payload = {
        "schema": "THE134_FACTORIAL_VIEW_WEIGHTED_WORKING_POINTS_V1",
        "status": "PASS" if all(gates.values()) else "FAIL",
        "promotion_status": "CANDIDATE_CURRENT_NOT_CANONICAL",
        "system": args.system,
        "shower_definition": args.view,
        "shower_semantic_sha256": shower_semantic_sha256(args.view),
        "model_domain_gev": list(MODEL_DOMAIN_GEV),
        "feature_order": list(FEATURES_BY_SYSTEM[args.system]),
        "model_origin": metadata.get("model_origin"),
        "reuse_pinned_hashes": metadata.get("reuse_pinned_hashes"),
        "extraction_authority": metadata.get("the134_view_contract", {}).get(
            "extraction_authority"
        ),
        "axis": axis_name,
        "bin_edges": edges.tolist(),
        "targets": list(TARGET_EFFICIENCIES),
        "threshold_authority": "exact binned weighted holdout signal quantiles",
        "comparison_operator": "score > threshold",
        "tie_policy": (
            "stable mergesort weighted quantile; score ties are never split; exact binned "
            "authority fails if strict score>threshold efficiency misses target tolerance"
        ),
        "exact_max_efficiency_error": args.exact_max_efficiency_error,
        "gates": gates,
        "population_certificate_gates": population_certificate_checks,
        "exact_thresholds": rows,
        "runtime_surfaces": surfaces,
        "provenance": {
            "score_sample": str(args.score_sample),
            "score_sample_sha256": sha256_file(args.score_sample),
            "population_certificate": str(args.population_certificate),
            "population_certificate_sha256": sha256_file(args.population_certificate),
            "model_metadata": str(args.model_metadata),
            "model_metadata_sha256": sha256_file(args.model_metadata),
            "model_xgb": population_provenance.get("model_xgb"),
            "model_xgb_sha256": population_provenance.get(
                "model_xgb_sha256"
            ),
            "model_receipt": population_provenance.get("model_receipt"),
            "model_receipt_sha256": population_provenance.get(
                "model_receipt_sha256"
            ),
        },
    }
    args.json_out.parent.mkdir(parents=True, exist_ok=True)
    args.json_out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(args.json_out)
    return 0 if payload["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
