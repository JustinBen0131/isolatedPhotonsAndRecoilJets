#!/usr/bin/env python3
"""Certify one THE-134 factorial-view model on its exact event-group holdout."""

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
    BACKGROUND_SOURCES,
    CONTROL_MODEL_ORIGIN_BY_SYSTEM,
    CONTROL_TMVA_SHA256_BY_SYSTEM,
    CONTROL_VIEW_BY_SYSTEM,
    FEATURES_BY_SYSTEM,
    MAX_BINNED_AUC_REGRESSION,
    MAX_OVERALL_AUC_REGRESSION,
    MAX_RUNTIME_SCORE_ABS_DIFFERENCE,
    MAX_TRAIN_HOLDOUT_AUC_GAP,
    MODEL_DOMAIN_GEV,
    READER_ET_EDGES,
    SIGNAL_SOURCES,
    VIEW_NAME,
    extraction_authority_binding,
    expected_model_origin,
    expected_model_split,
    full_extraction_authority_checks,
    is_sha256,
    model_extraction_authority_checks,
    model_artifact_receipt_checks,
    sha256_file,
    shower_semantic_sha256,
    valid_reuse_pinned_hashes,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--system", choices=sorted(FEATURES_BY_SYSTEM), required=True)
    parser.add_argument("--view", choices=ALL_SHOWER_VIEWS, default=VIEW_NAME)
    parser.add_argument("--control-view", choices=ALL_SHOWER_VIEWS, required=True)
    parser.add_argument("--holdout", type=Path, required=True)
    parser.add_argument("--holdout-certificate", type=Path, required=True)
    parser.add_argument("--model-xgb", type=Path, required=True)
    parser.add_argument("--model-tmva", type=Path, required=True)
    parser.add_argument("--model-metadata", type=Path, required=True)
    parser.add_argument("--model-receipt", type=Path, required=True)
    parser.add_argument("--extraction-audit", type=Path, required=True)
    parser.add_argument("--control-tmva", type=Path, required=True)
    parser.add_argument("--control-metadata", type=Path, required=True)
    parser.add_argument("--control-view-certificate", type=Path, required=True)
    parser.add_argument("--json-out", type=Path, required=True)
    parser.add_argument("--eta-edges", default="-0.7,-0.35,0,0.35,0.7")
    parser.add_argument(
        "--runtime-tolerance", type=float, default=MAX_RUNTIME_SCORE_ABS_DIFFERENCE
    )
    parser.add_argument(
        "--max-auc-regression", type=float, default=MAX_OVERALL_AUC_REGRESSION
    )
    parser.add_argument(
        "--max-binned-auc-regression", type=float, default=MAX_BINNED_AUC_REGRESSION
    )
    parser.add_argument(
        "--max-train-holdout-auc-gap", type=float, default=MAX_TRAIN_HOLDOUT_AUC_GAP
    )
    return parser.parse_args()


def enforce_frozen_gate_maxima(args: argparse.Namespace) -> None:
    gates = {
        "runtime_tolerance": (args.runtime_tolerance, MAX_RUNTIME_SCORE_ABS_DIFFERENCE),
        "max_auc_regression": (args.max_auc_regression, MAX_OVERALL_AUC_REGRESSION),
        "max_binned_auc_regression": (
            args.max_binned_auc_regression,
            MAX_BINNED_AUC_REGRESSION,
        ),
        "max_train_holdout_auc_gap": (
            args.max_train_holdout_auc_gap,
            MAX_TRAIN_HOLDOUT_AUC_GAP,
        ),
    }
    invalid = {
        name: {"requested": value, "frozen_maximum": maximum}
        for name, (value, maximum) in gates.items()
        if not math.isfinite(value) or value < 0.0 or value > maximum
    }
    if invalid:
        raise SystemExit("scientific gate widening is forbidden: " + json.dumps(invalid, sort_keys=True))


def parse_edges(text: str) -> np.ndarray:
    edges = np.asarray([float(item) for item in text.split(",")], dtype=np.float64)
    if len(edges) < 2 or not np.all(np.diff(edges) > 0):
        raise SystemExit(f"invalid edges: {text}")
    return edges


def metadata_view(metadata: dict) -> tuple[object, object]:
    contract = metadata.get("the134_view_contract", {})
    return (
        metadata.get("shower_definition") or contract.get("shower_definition"),
        metadata.get("shower_semantic_sha256") or contract.get("shower_semantic_sha256"),
    )


def validate_metadata(metadata: dict, system: str, *, expected_view: str) -> dict:
    view, semantic = metadata_view(metadata)
    nested_contract = metadata.get("the134_view_contract", {})
    top_view = metadata.get("shower_definition")
    top_semantic = metadata.get("shower_semantic_sha256")
    nested_view = nested_contract.get("shower_definition")
    nested_semantic = nested_contract.get("shower_semantic_sha256")
    weighting = metadata.get("weighting", {})
    expected_split = expected_model_split(system, expected_view)
    expected_origin = expected_model_origin(system, expected_view)
    reuse_pins = metadata.get("reuse_pinned_hashes", {})
    split = metadata.get("split", {})
    positive_rows = int(split.get("train_rows", 0)) > 0 and int(
        split.get("test_rows", 0)
    ) > 0
    if expected_split["mode"] == "event50":
        split_population = (
            positive_rows
            and int(split.get("train_events", 0)) > 0
            and int(split.get("test_events", 0)) > 0
        )
        accepted_boundary = True
    else:
        split_population = positive_rows
        accepted_boundary = metadata.get("accepted_split_contract") == expected_split
    checks = {
        "features": metadata.get("features") == list(FEATURES_BY_SYSTEM[system]),
        "pt_range": metadata.get("pt_range") == list(MODEL_DOMAIN_GEV),
        "weight_mode": metadata.get("weight_mode") == "ppg12-exact",
        "split_mode": metadata.get("split", {}).get("mode")
        == expected_split["mode"],
        "split_fraction": metadata.get("split", {}).get("test_fraction_requested")
        == expected_split["test_fraction_requested"],
        "split_seed": metadata.get("split", {}).get("random_seed")
        == expected_split["random_seed"],
        "split_population": split_population,
        "accepted_legacy_split_boundary": accepted_boundary,
        "event_weight_unused": weighting.get("event_weight_used") is False,
        "vertex_weight_unused": weighting.get("vertex_reweight") is False,
        "centrality_weight_unused": weighting.get("centrality_event_weight") is False,
        "cross_section_weight_unused": weighting.get("cross_section_weight_used_for_training")
        is False,
        "weights_computed_before_binning": weighting.get("weights_computed_before_binning")
        is True,
        "old_low_calo_veto_disabled": not bool(
            metadata.get("event_quality_filter", {}).get("enabled", False)
        ),
        "model_origin": metadata.get("model_origin") == expected_origin,
        "reuse_pins": (
            valid_reuse_pinned_hashes(reuse_pins, system, expected_view)
            if expected_origin.startswith("REUSED_")
            else not reuse_pins
        ),
        "view_contract_schema": nested_contract.get("schema")
        == "THE134_FACTORIAL_VIEW_MODEL_VIEW_CONTRACT_V1",
        "view_contract_system": nested_contract.get("system") == system,
        "top_nested_view_consistency": top_view in {None, nested_view},
        "top_nested_semantic_consistency": top_semantic in {
            None,
            nested_semantic,
        },
    }
    checks.update(
        {
            "shower_definition": view == expected_view,
            "shower_semantic_sha256": semantic
            == shower_semantic_sha256(expected_view),
        }
    )
    return checks


def control_view_certificate_checks(
    certificate: dict,
    system: str,
    *,
    expected_control_view: str,
    model_tmva_sha256: str,
    model_metadata_sha256: str,
) -> dict:
    artifacts = certificate.get("artifacts", {})
    observed_tmva_sha256 = certificate.get("model_tmva_sha256") or artifacts.get(
        "tmva_sha256"
    )
    observed_metadata_sha256 = certificate.get(
        "model_metadata_sha256"
    ) or artifacts.get("metadata_sha256")
    return {
        "schema": certificate.get("schema")
        in {
            "THE134_CONTROL_VIEW_AUTHORITY_V1",
            "THE134_FACTORIAL_VIEW_MODEL_REUSE_AUDIT_V1",
        },
        "status": certificate.get("status") == "PASS",
        "system": certificate.get("system") == system,
        "shower_definition": certificate.get("shower_definition")
        == expected_control_view,
        "shower_semantic_sha256": certificate.get("shower_semantic_sha256")
        == shower_semantic_sha256(expected_control_view),
        "feature_order": certificate.get("feature_order")
        == list(FEATURES_BY_SYSTEM[system]),
        "model_origin": certificate.get("model_origin")
        == CONTROL_MODEL_ORIGIN_BY_SYSTEM[system],
        "model_tmva_sha256": observed_tmva_sha256
        == model_tmva_sha256
        == CONTROL_TMVA_SHA256_BY_SYSTEM[system],
        "model_metadata_sha256": observed_metadata_sha256 == model_metadata_sha256,
    }


def score_xgb(path: Path, matrix: np.ndarray) -> np.ndarray:
    try:
        import xgboost as xgb
    except ImportError as exc:
        raise SystemExit("validate_the134_h70_model.py requires xgboost") from exc
    booster = xgb.Booster()
    booster.load_model(str(path))
    return np.asarray(booster.predict(xgb.DMatrix(matrix)), dtype=np.float64)


def score_tmva(path: Path, matrix: np.ndarray, chunk_size: int = 100_000) -> np.ndarray:
    try:
        import ROOT
    except ImportError as exc:
        raise SystemExit("validate_the134_h70_model.py requires PyROOT for runtime parity") from exc
    runtime = ROOT.TMVA.Experimental.RBDT("myBDT", str(path))
    observed = np.empty(len(matrix), dtype=np.float64)
    for start in range(0, len(matrix), chunk_size):
        stop = min(len(matrix), start + chunk_size)
        observed[start:stop] = np.asarray(
            runtime.Compute(matrix[start:stop]), dtype=np.float64
        ).reshape(-1)
    return observed


def binned_auc(labels, scores, weights, values, edges) -> list[dict]:
    try:
        from sklearn.metrics import roc_auc_score
    except ImportError as exc:
        raise SystemExit("validate_the134_h70_model.py requires scikit-learn") from exc
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        mask = np.isfinite(values) & (values >= lo) & (values < hi)
        signal = int(np.sum(mask & (labels == 1)))
        background = int(np.sum(mask & (labels == 0)))
        auc = (
            float(roc_auc_score(labels[mask], scores[mask], sample_weight=weights[mask]))
            if signal > 0 and background > 0
            else None
        )
        rows.append(
            {
                "lo": float(lo),
                "hi": float(hi),
                "rows": int(np.sum(mask)),
                "signal": signal,
                "background": background,
                "auc": auc,
            }
        )
    return rows


def main() -> int:
    args = parse_args()
    enforce_frozen_gate_maxima(args)
    if args.control_view != CONTROL_VIEW_BY_SYSTEM[args.system]:
        raise SystemExit(
            f"control view for {args.system} is frozen to "
            f"{CONTROL_VIEW_BY_SYSTEM[args.system]}, not {args.control_view}"
        )
    paths = (
        args.holdout,
        args.holdout_certificate,
        args.model_xgb,
        args.model_tmva,
        args.model_metadata,
        args.model_receipt,
        args.extraction_audit,
        args.control_tmva,
        args.control_metadata,
        args.control_view_certificate,
    )
    missing_paths = [str(path) for path in paths if not path.is_file()]
    if missing_paths:
        raise SystemExit(f"missing inputs: {missing_paths}")
    metadata = json.loads(args.model_metadata.read_text())
    model_receipt = json.loads(args.model_receipt.read_text())
    control_metadata = json.loads(args.control_metadata.read_text())
    control_view_certificate = json.loads(args.control_view_certificate.read_text())
    extraction = json.loads(args.extraction_audit.read_text())
    holdout_certificate = json.loads(args.holdout_certificate.read_text())
    metadata_gates = validate_metadata(metadata, args.system, expected_view=args.view)
    control_metadata_gates = validate_metadata(
        control_metadata, args.system, expected_view=args.control_view
    )
    expected_control_view = args.control_view
    control_view_gates = control_view_certificate_checks(
        control_view_certificate,
        args.system,
        expected_control_view=expected_control_view,
        model_tmva_sha256=sha256_file(args.control_tmva),
        model_metadata_sha256=sha256_file(args.control_metadata),
    )
    extraction_gates = {
        "status": extraction.get("status") == "PASS",
        "system": extraction.get("system") == args.system,
        "shower_definition": extraction.get("shower_definition") == args.view,
        "shower_semantic_sha256": extraction.get("shower_semantic_sha256")
        == shower_semantic_sha256(args.view),
        "full_authority": all(
            full_extraction_authority_checks(
                extraction, args.system, args.view
            ).values()
        ),
        "model_authority_binding": all(
            model_extraction_authority_checks(metadata).values()
        )
        and metadata.get("the134_view_contract", {}).get(
            "extraction_authority"
        )
        == extraction_authority_binding(extraction),
    }
    holdout_certificate_gates = {
        "schema": holdout_certificate.get("schema")
        == "THE134_FACTORIAL_VIEW_EXACT_HOLDOUT_CERTIFICATE_V1",
        "status": holdout_certificate.get("status") == "PASS",
        "system": holdout_certificate.get("system") == args.system,
        "shower_definition": holdout_certificate.get("shower_definition") == args.view,
        "shower_semantic_sha256": holdout_certificate.get("shower_semantic_sha256")
        == shower_semantic_sha256(args.view),
        "feature_order": holdout_certificate.get("feature_order")
        == list(FEATURES_BY_SYSTEM[args.system]),
        "model_origin": holdout_certificate.get("model_origin")
        == metadata.get("model_origin"),
        "reuse_pinned_hashes": holdout_certificate.get("reuse_pinned_hashes")
        == metadata.get("reuse_pinned_hashes"),
        "control_shower_definition": holdout_certificate.get(
            "control_shower_definition"
        )
        == expected_control_view,
        "control_shower_semantic_sha256": holdout_certificate.get(
            "control_shower_semantic_sha256"
        )
        == shower_semantic_sha256(expected_control_view),
        "zero_event_group_overlap": holdout_certificate.get("split", {}).get("event_overlap")
        == 0,
        "holdout_sha256": holdout_certificate.get("provenance", {}).get("holdout_sha256")
        == sha256_file(args.holdout),
        "model_xgb_sha256": holdout_certificate.get("provenance", {}).get("model_xgb_sha256")
        == sha256_file(args.model_xgb),
        "model_metadata_sha256": holdout_certificate.get("provenance", {}).get(
            "model_metadata_sha256"
        )
        == sha256_file(args.model_metadata),
        "extraction_audit_sha256": holdout_certificate.get("provenance", {}).get(
            "extraction_audit_sha256"
        )
        == sha256_file(args.extraction_audit),
        "extraction_authority_binding": holdout_certificate.get(
            "extraction_authority"
        )
        == extraction_authority_binding(extraction),
        "model_receipt_sha256": holdout_certificate.get("provenance", {}).get(
            "model_receipt_sha256"
        )
        == sha256_file(args.model_receipt),
        "model_artifact_receipt": all(
            model_artifact_receipt_checks(
                model_receipt,
                args.system,
                args.view,
                model_xgb_sha256=sha256_file(args.model_xgb),
                model_metadata_sha256=sha256_file(args.model_metadata),
            ).values()
        ),
    }

    data = np.load(args.holdout, allow_pickle=True)
    required = {
        "features",
        "x",
        "is_signal",
        "training_weight",
        "score_xgboost",
        "cluster_Et",
        "cluster_Eta",
        "source_sample",
        "x_control",
        "control_features",
        "control_shower_definition",
        "control_shower_semantic_sha256",
    }
    missing = sorted(required - set(data.files))
    if missing:
        raise SystemExit(f"holdout missing fields: {missing}")
    features = [str(item) for item in data["features"].tolist()]
    matrix = np.ascontiguousarray(data["x"], dtype=np.float32)
    control_matrix = np.ascontiguousarray(data["x_control"], dtype=np.float32)
    labels = np.asarray(data["is_signal"], dtype=np.int8)
    weights = np.asarray(data["training_weight"], dtype=np.float64)
    cached_scores = np.asarray(data["score_xgboost"], dtype=np.float64)
    et = np.asarray(data["cluster_Et"], dtype=np.float64)
    eta = np.asarray(data["cluster_Eta"], dtype=np.float64)
    sources = np.asarray([str(item) for item in data["source_sample"].tolist()])
    if len(matrix) == 0:
        raise SystemExit("empty holdout")
    if matrix.shape != (len(labels), len(FEATURES_BY_SYSTEM[args.system])):
        raise SystemExit(f"holdout matrix shape mismatch: {matrix.shape}")
    if control_matrix.shape != matrix.shape:
        raise SystemExit(f"control holdout matrix shape mismatch: {control_matrix.shape}")
    observed_control_features = [str(item) for item in data["control_features"].tolist()]
    observed_control_view = str(np.asarray(data["control_shower_definition"]).item())
    observed_control_semantic = str(
        np.asarray(data["control_shower_semantic_sha256"]).item()
    )

    candidate_scores = score_xgb(args.model_xgb, matrix)
    control_scores = score_tmva(args.control_tmva, control_matrix)
    tmva_scores = score_tmva(args.model_tmva, matrix)
    candidate_tmva_delta = np.abs(candidate_scores - tmva_scores)
    cache_delta = np.abs(candidate_scores - cached_scores)
    finite = (
        np.isfinite(matrix).all(axis=1)
        & np.isfinite(control_matrix).all(axis=1)
        & np.isfinite(candidate_scores)
        & np.isfinite(control_scores)
        & np.isfinite(tmva_scores)
        & np.isfinite(weights)
        & (weights > 0.0)
        & np.isin(labels, [0, 1])
    )
    expected_source_labels = {
        **{source: 1 for source in SIGNAL_SOURCES[args.system]},
        **{source: 0 for source in BACKGROUND_SOURCES[args.system]},
    }
    wrong_source = ~np.isin(sources, sorted(expected_source_labels))
    wrong_label = np.asarray(
        [expected_source_labels.get(source, -999) != int(label) for source, label in zip(sources, labels)]
    )

    try:
        from sklearn.metrics import roc_auc_score
    except ImportError as exc:
        raise SystemExit("validate_the134_h70_model.py requires scikit-learn") from exc
    candidate_auc = float(roc_auc_score(labels, candidate_scores, sample_weight=weights))
    control_auc = float(roc_auc_score(labels, control_scores, sample_weight=weights))
    et_edges = np.asarray(READER_ET_EDGES, dtype=np.float64)
    eta_edges = parse_edges(args.eta_edges)
    candidate_et = binned_auc(labels, candidate_scores, weights, et, et_edges)
    control_et = binned_auc(labels, control_scores, weights, et, et_edges)
    candidate_eta = binned_auc(labels, candidate_scores, weights, eta, eta_edges)
    control_eta = binned_auc(labels, control_scores, weights, eta, eta_edges)
    et_deltas = [
        left["auc"] - right["auc"]
        for left, right in zip(candidate_et, control_et)
        if left["auc"] is not None and right["auc"] is not None
    ]
    eta_deltas = [
        left["auc"] - right["auc"]
        for left, right in zip(candidate_eta, control_eta)
        if left["auc"] is not None and right["auc"] is not None
    ]
    overfit = metadata.get("overfit_diagnostics", {})
    auc_gap = overfit.get("auc_gap_train_minus_holdout")

    critical_gates = {
        "candidate_metadata": all(metadata_gates.values()),
        "control_metadata": all(control_metadata_gates.values()),
        "control_view_authority": all(control_view_gates.values()),
        "extraction_certificate": all(extraction_gates.values()),
        "exact_event_group_holdout_certificate": all(holdout_certificate_gates.values()),
        "feature_order": features == list(FEATURES_BY_SYSTEM[args.system]),
        "control_feature_order": observed_control_features
        == list(FEATURES_BY_SYSTEM[args.system]),
        "control_view_identity": observed_control_view == expected_control_view
        and observed_control_semantic == shower_semantic_sha256(expected_control_view),
        "source_role_label_closure": int(np.sum(wrong_source | wrong_label)) == 0,
        "finite_inputs_scores_positive_weights": bool(np.all(finite)),
        "cached_score_reconstruction": float(np.max(cache_delta)) <= 1.0e-10,
        "python_tmva_runtime_parity": float(np.max(candidate_tmva_delta))
        <= args.runtime_tolerance,
        "no_material_overall_auc_regression": candidate_auc - control_auc
        >= -args.max_auc_regression,
        "et_bins_populated": all(row["signal"] > 0 and row["background"] > 0 for row in candidate_et),
        "eta_bins_populated": all(row["signal"] > 0 and row["background"] > 0 for row in candidate_eta),
        "et_resolved_stability": bool(et_deltas)
        and min(et_deltas) >= -args.max_binned_auc_regression,
        "eta_resolved_stability": bool(eta_deltas)
        and min(eta_deltas) >= -args.max_binned_auc_regression,
        "overtraining_gate": isinstance(auc_gap, (float, int))
        and float(auc_gap) <= args.max_train_holdout_auc_gap,
    }
    payload = {
        "schema": "THE134_FACTORIAL_VIEW_MODEL_VALIDATION_V1",
        "status": "PASS" if all(critical_gates.values()) else "FAIL",
        "promotion_status": "CANDIDATE_CURRENT_NOT_CANONICAL",
        "system": args.system,
        "shower_definition": args.view,
        "shower_semantic_sha256": shower_semantic_sha256(args.view),
        "feature_order": list(FEATURES_BY_SYSTEM[args.system]),
        "model_domain_gev": list(MODEL_DOMAIN_GEV),
        "model_origin": metadata.get("model_origin"),
        "reuse_pinned_hashes": metadata.get("reuse_pinned_hashes"),
        "extraction_authority": extraction_authority_binding(extraction),
        "critical_gates": critical_gates,
        "metadata_gates": metadata_gates,
        "control_metadata_gates": control_metadata_gates,
        "control_view_gates": control_view_gates,
        "extraction_gates": extraction_gates,
        "holdout_certificate_gates": holdout_certificate_gates,
        "runtime_parity": {
            "rows": int(len(matrix)),
            "max_python_tmva_abs_difference": float(np.max(candidate_tmva_delta)),
            "max_cached_python_abs_difference": float(np.max(cache_delta)),
            "tolerance": args.runtime_tolerance,
        },
        "source_label_closure": {
            "observed_sources": sorted(set(sources.tolist())),
            "wrong_source_rows": int(np.sum(wrong_source)),
            "wrong_label_rows": int(np.sum(wrong_label)),
        },
        "performance": {
            "candidate_weighted_holdout_auc": candidate_auc,
            "control_weighted_same_holdout_auc": control_auc,
            "candidate_minus_control_auc": candidate_auc - control_auc,
            "train_minus_holdout_auc": auc_gap,
            "et_candidate": candidate_et,
            "et_control": control_et,
            "eta_candidate": candidate_eta,
            "eta_control": control_eta,
            "et_candidate_minus_control_auc": et_deltas,
            "eta_candidate_minus_control_auc": eta_deltas,
        },
        "provenance": {
            "holdout": str(args.holdout),
            "holdout_sha256": sha256_file(args.holdout),
            "holdout_certificate": str(args.holdout_certificate),
            "holdout_certificate_sha256": sha256_file(args.holdout_certificate),
            "model_xgb": str(args.model_xgb),
            "model_xgb_sha256": sha256_file(args.model_xgb),
            "model_tmva": str(args.model_tmva),
            "model_tmva_sha256": sha256_file(args.model_tmva),
            "model_metadata": str(args.model_metadata),
            "model_metadata_sha256": sha256_file(args.model_metadata),
            "model_receipt": str(args.model_receipt),
            "model_receipt_sha256": sha256_file(args.model_receipt),
            "extraction_audit": str(args.extraction_audit),
            "extraction_audit_sha256": sha256_file(args.extraction_audit),
            "control_tmva": str(args.control_tmva),
            "control_tmva_sha256": sha256_file(args.control_tmva),
            "control_metadata": str(args.control_metadata),
            "control_metadata_sha256": sha256_file(args.control_metadata),
            "control_view_certificate": str(args.control_view_certificate),
            "control_view_certificate_sha256": sha256_file(args.control_view_certificate),
        },
    }
    args.json_out.parent.mkdir(parents=True, exist_ok=True)
    args.json_out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(args.json_out)
    return 0 if payload["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
