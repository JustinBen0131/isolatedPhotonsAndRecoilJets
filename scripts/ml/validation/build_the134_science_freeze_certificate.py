#!/usr/bin/env python3
"""Build the training-source-complete seven-view THE-134 freeze certificate."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
from the134_aggregate_contract import (  # noqa: E402
    MODEL_ARTIFACT_KEYS,
    MODEL_CERTIFICATE_ARTIFACT_PAIRS,
    MODEL_CERTIFICATE_KEYS,
    MODEL_REGISTRY_KEYS,
    MODEL_VALIDATION_DETAIL_GATE_FIELDS,
    MODEL_VALIDATION_GATES,
    MODEL_VALIDATION_KEYS,
    MODEL_VALIDATION_SCHEMA,
    PAIRED_REGISTRY_SCHEMA,
    REPLAY_CERTIFICATE_GATES,
    REPLAY_CERTIFICATE_SCHEMA,
    REQUIRED_SOURCE_KEYS,
    REQUIRED_SOURCES_BY_SYSTEM,
    REQUIRED_SYSTEM_VIEW_KEYS,
    REQUIRED_VIEW_KEYS,
    SCIENCE_CERTIFICATE_SCHEMA,
    SCIENCE_MANIFEST_SCHEMA,
    SOURCE_WITNESS_GATES,
    SYSTEMS,
    WORKING_POINT_AXIS_BY_SYSTEM,
    WORKING_POINT_COMPARISON_OPERATOR,
    WORKING_POINT_EDGES_BY_SYSTEM,
    WORKING_POINT_GATES,
    WORKING_POINT_LABEL_BY_TARGET,
    WORKING_POINT_THRESHOLD_AUTHORITY,
    WORKING_POINT_THRESHOLD_ROW_KEYS,
    WORKING_POINT_TIE_POLICY,
    WORKING_POINT_RECEIPT_KEYS,
    WORKING_POINTS_SCHEMA,
    ContractError,
    expected_view_semantic,
    require_exact_key_inventory,
    require_exact_true_gates,
    require_git_commit,
    require_sha256,
    read_json,
    semantic_receipt,
    verify_registered_file,
    verify_registered_json,
)

_CONTRACTS = _HERE.parents[1] / "contracts"
sys.path.insert(0, str(_CONTRACTS))
from the134_h70_contract import (  # noqa: E402
    FEATURES_BY_SYSTEM,
    MAX_RUNTIME_SCORE_ABS_DIFFERENCE,
    MAX_SURFACE_ABS_RESIDUAL,
    MAX_SURFACE_EFFICIENCY_ERROR,
    MAX_SURFACE_RMS,
    MAX_WP_EXACT_EFFICIENCY_ERROR,
    MODEL_DOMAIN_GEV,
    TARGET_EFFICIENCIES,
    canonical_json_sha256,
    expected_model_origin,
    extraction_authority_binding,
    full_extraction_authority_checks,
    model_artifact_receipt_checks,
    valid_reuse_pinned_hashes,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--json-out", type=Path, required=True)
    return parser.parse_args()


def require_exact_keys(payload: object, expected: tuple[str, ...], label: str) -> dict:
    if not isinstance(payload, dict):
        raise ContractError(f"{label} must be an object")
    missing = set(expected) - set(payload)
    extra = set(payload) - set(expected)
    if missing or extra:
        raise ContractError(
            f"{label} key inventory mismatch: missing={sorted(missing)} "
            f"extra={sorted(extra)}"
        )
    return payload


def require_finite_number(value: object, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ContractError(f"{label} must be a finite number")
    result = float(value)
    if not math.isfinite(result):
        raise ContractError(f"{label} must be a finite number")
    return result


def validate_exact_thresholds(
    system: str,
    rows: object,
    label: str,
    *,
    max_efficiency_error: float = MAX_WP_EXACT_EFFICIENCY_ERROR,
) -> list[dict]:
    if not isinstance(rows, list):
        raise ContractError(f"{label} must be a list")
    edges = WORKING_POINT_EDGES_BY_SYSTEM[system]
    targets = tuple(float(value) for value in TARGET_EFFICIENCIES)
    expected_count = len(targets) * (len(edges) - 1)
    if len(rows) != expected_count:
        raise ContractError(
            f"{label} row count={len(rows)}, expected={expected_count}"
        )
    normalized = []
    for index, row in enumerate(rows):
        row = require_exact_keys(
            row, WORKING_POINT_THRESHOLD_ROW_KEYS, f"{label}[{index}]"
        )
        target_index = index // (len(edges) - 1)
        bin_index = index % (len(edges) - 1)
        target = targets[target_index]
        bin_lo, bin_hi = edges[bin_index], edges[bin_index + 1]
        expected_wp = WORKING_POINT_LABEL_BY_TARGET[target]
        expected = {
            "wp": expected_wp,
            "target_signal_efficiency": target,
            "bin_lo": bin_lo,
            "bin_hi": bin_hi,
            "bin_center": 0.5 * (bin_lo + bin_hi),
        }
        for field, expected_value in expected.items():
            observed = row.get(field)
            if field == "wp":
                if observed != expected_value:
                    raise ContractError(
                        f"{label}[{index}].wp={observed!r}, "
                        f"expected={expected_value!r}"
                    )
            elif require_finite_number(
                observed, f"{label}[{index}].{field}"
            ) != expected_value:
                raise ContractError(
                    f"{label}[{index}].{field}={observed!r}, "
                    f"expected={expected_value!r}"
                )
        for field in (
            "threshold",
            "achieved_signal_efficiency",
            "achieved_minus_target",
            "abs_efficiency_error",
            "signal_weight_fraction_tied_at_threshold",
            "background_acceptance",
        ):
            require_finite_number(row.get(field), f"{label}[{index}].{field}")
        achieved = float(row["achieved_signal_efficiency"])
        achieved_minus_target = float(row["achieved_minus_target"])
        abs_error = float(row["abs_efficiency_error"])
        if not math.isclose(
            achieved_minus_target,
            achieved - target,
            rel_tol=0.0,
            abs_tol=1.0e-15,
        ):
            raise ContractError(
                f"{label}[{index}].achieved_minus_target is internally inconsistent"
            )
        if not math.isclose(
            abs_error,
            abs(achieved_minus_target),
            rel_tol=0.0,
            abs_tol=1.0e-15,
        ):
            raise ContractError(
                f"{label}[{index}].abs_efficiency_error is internally inconsistent"
            )
        for field in (
            "achieved_signal_efficiency",
            "signal_weight_fraction_tied_at_threshold",
            "background_acceptance",
        ):
            if not 0.0 <= float(row[field]) <= 1.0:
                raise ContractError(
                    f"{label}[{index}].{field} must be within [0, 1]"
                )
        if require_finite_number(
            row.get("abs_efficiency_error"),
            f"{label}[{index}].abs_efficiency_error",
        ) > max_efficiency_error:
            raise ContractError(
                f"{label}[{index}] exceeds frozen exact-efficiency tolerance"
            )
        for field in ("signal_rows", "background_rows"):
            value = row.get(field)
            if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
                raise ContractError(f"{label}[{index}].{field} must be positive")
        normalized.append(dict(row))
    return normalized


def validate_runtime_surfaces(
    system: str,
    value: object,
    label: str,
) -> dict:
    if not isinstance(value, dict):
        raise ContractError(f"{label} must be an object")
    expected_wp = set(WORKING_POINT_LABEL_BY_TARGET.values())
    if set(value) != expected_wp:
        raise ContractError(
            f"{label} inventory mismatch: missing={sorted(expected_wp-set(value))} "
            f"extra={sorted(set(value)-expected_wp)}"
        )
    for wp, surface in value.items():
        surface = require_exact_keys(
            surface,
            (
                "status",
                "intercept",
                "slope_per_axis_unit",
                "rms_residual",
                "max_abs_residual",
                "max_abs_efficiency_error",
                "achieved_efficiency_by_bin",
                "gates",
            ),
            f"{label}.{wp}",
        )
        if surface.get("status") not in {
            "ACCEPTED",
            "BINNED_THRESHOLDS_ONLY",
        }:
            raise ContractError(
                f"{label}.{wp}.status={surface.get('status')!r} is invalid"
            )
        for field in (
            "intercept",
            "slope_per_axis_unit",
            "rms_residual",
            "max_abs_residual",
            "max_abs_efficiency_error",
        ):
            require_finite_number(surface.get(field), f"{label}.{wp}.{field}")
        achieved = surface.get("achieved_efficiency_by_bin")
        expected_bins = len(WORKING_POINT_EDGES_BY_SYSTEM[system]) - 1
        if not isinstance(achieved, list) or len(achieved) != expected_bins:
            raise ContractError(
                f"{label}.{wp}.achieved_efficiency_by_bin must contain "
                f"{expected_bins} bins"
            )
        for index, efficiency in enumerate(achieved):
            finite = require_finite_number(
                efficiency,
                f"{label}.{wp}.achieved_efficiency_by_bin[{index}]",
            )
            if not 0.0 <= finite <= 1.0:
                raise ContractError(
                    f"{label}.{wp}.achieved_efficiency_by_bin[{index}] "
                    "must be within [0, 1]"
                )
        gates = require_exact_keys(
            surface.get("gates"),
            ("max_rms", "max_abs_residual", "max_abs_efficiency_error"),
            f"{label}.{wp}.gates",
        )
        maximum_gates = {
            "max_rms": MAX_SURFACE_RMS,
            "max_abs_residual": MAX_SURFACE_ABS_RESIDUAL,
            "max_abs_efficiency_error": MAX_SURFACE_EFFICIENCY_ERROR,
        }
        normalized_gates = {}
        for gate, frozen_maximum in maximum_gates.items():
            observed = require_finite_number(
                gates.get(gate), f"{label}.{wp}.gates.{gate}"
            )
            if observed < 0.0 or observed > frozen_maximum:
                raise ContractError(
                    f"{label}.{wp}.gates.{gate}={observed} exceeds its "
                    f"frozen maximum={frozen_maximum}"
                )
            normalized_gates[gate] = observed
        accepted = (
            float(surface["rms_residual"]) <= normalized_gates["max_rms"]
            and float(surface["max_abs_residual"])
            <= normalized_gates["max_abs_residual"]
            and float(surface["max_abs_efficiency_error"])
            <= normalized_gates["max_abs_efficiency_error"]
        )
        expected_status = "ACCEPTED" if accepted else "BINNED_THRESHOLDS_ONLY"
        if surface["status"] != expected_status:
            raise ContractError(
                f"{label}.{wp}.status={surface['status']!r}, "
                f"expected={expected_status!r} from frozen gates"
            )
    return value


def validate_validation_receipt(
    system: str,
    view: str,
    receipt: dict,
    model: dict,
    model_origin: str,
    reuse_pinned_hashes: object,
    label: str,
) -> None:
    receipt = require_exact_keys(receipt, MODEL_VALIDATION_KEYS, label)
    failures = []
    expected = {
        "schema": MODEL_VALIDATION_SCHEMA,
        "status": "PASS",
        "system": system,
        "shower_definition": view,
        "shower_semantic_sha256": expected_view_semantic(view),
        "feature_order": list(FEATURES_BY_SYSTEM[system]),
        "model_domain_gev": list(MODEL_DOMAIN_GEV),
        "model_origin": model_origin,
        "reuse_pinned_hashes": reuse_pinned_hashes,
        "promotion_status": "CANDIDATE_CURRENT_NOT_CANONICAL",
    }
    for field, wanted in expected.items():
        if receipt.get(field) != wanted:
            failures.append(f"{field} mismatch")
    if failures:
        raise ContractError(f"{label} identity failed: {failures}")
    require_exact_true_gates(
        receipt.get("critical_gates"), MODEL_VALIDATION_GATES, label
    )
    for field in MODEL_VALIDATION_DETAIL_GATE_FIELDS:
        gates = receipt.get(field)
        if not isinstance(gates, dict) or not gates:
            raise ContractError(f"{label}.{field} must be a non-empty object")
        if any(value is not True for value in gates.values()):
            raise ContractError(f"{label}.{field} is not all true")
    source_closure = receipt.get("source_label_closure")
    if not isinstance(source_closure, dict):
        raise ContractError(f"{label}.source_label_closure must be an object")
    expected_sources = sorted(REQUIRED_SOURCES_BY_SYSTEM[system])
    if (
        source_closure.get("observed_sources") != expected_sources
        or source_closure.get("wrong_source_rows") != 0
        or source_closure.get("wrong_label_rows") != 0
    ):
        raise ContractError(f"{label}.source_label_closure failed")
    if not isinstance(receipt.get("performance"), dict) or not receipt["performance"]:
        raise ContractError(f"{label}.performance must be a non-empty object")
    runtime = receipt.get("runtime_parity")
    if not isinstance(runtime, dict):
        raise ContractError(f"{label}.runtime_parity must be an object")
    rows = runtime.get("rows")
    if isinstance(rows, bool) or not isinstance(rows, int) or rows <= 0:
        raise ContractError(f"{label}.runtime_parity.rows must be positive")
    tolerance = require_finite_number(
        runtime.get("tolerance"), f"{label}.runtime_parity.tolerance"
    )
    if tolerance < 0.0 or tolerance > MAX_RUNTIME_SCORE_ABS_DIFFERENCE:
        raise ContractError(
            f"{label}.runtime_parity.tolerance={tolerance}, "
            f"frozen maximum={MAX_RUNTIME_SCORE_ABS_DIFFERENCE}"
        )
    runtime_delta = require_finite_number(
        runtime.get("max_python_tmva_abs_difference"),
        f"{label}.runtime_parity.max_python_tmva_abs_difference",
    )
    if runtime_delta > tolerance:
        raise ContractError(f"{label} Python/TMVA runtime parity failed")
    cache_delta = require_finite_number(
        runtime.get("max_cached_python_abs_difference"),
        f"{label}.runtime_parity.max_cached_python_abs_difference",
    )
    if cache_delta > 1.0e-10:
        raise ContractError(f"{label} cached/Python score identity failed")
    provenance = receipt.get("provenance")
    if not isinstance(provenance, dict):
        raise ContractError(f"{label}.provenance must be an object")
    for registry_key, receipt_key in (
        ("xgboost_sha256", "model_xgb_sha256"),
        ("tmva_sha256", "model_tmva_sha256"),
        ("metadata_sha256", "model_metadata_sha256"),
    ):
        if provenance.get(receipt_key) != model.get(registry_key):
            raise ContractError(
                f"{label}.{receipt_key} does not match registry model"
            )


def validate_working_point_receipt(
    system: str,
    view: str,
    receipt: dict,
    model: dict,
    model_origin: str,
    reuse_pinned_hashes: object,
    label: str,
) -> list[dict]:
    receipt = require_exact_keys(receipt, WORKING_POINT_RECEIPT_KEYS, label)
    expected = {
        "schema": WORKING_POINTS_SCHEMA,
        "status": "PASS",
        "system": system,
        "shower_definition": view,
        "shower_semantic_sha256": expected_view_semantic(view),
        "feature_order": list(FEATURES_BY_SYSTEM[system]),
        "model_domain_gev": list(MODEL_DOMAIN_GEV),
        "model_origin": model_origin,
        "axis": WORKING_POINT_AXIS_BY_SYSTEM[system],
        "bin_edges": list(WORKING_POINT_EDGES_BY_SYSTEM[system]),
        "targets": list(TARGET_EFFICIENCIES),
        "threshold_authority": WORKING_POINT_THRESHOLD_AUTHORITY,
        "comparison_operator": WORKING_POINT_COMPARISON_OPERATOR,
        "tie_policy": WORKING_POINT_TIE_POLICY,
        "reuse_pinned_hashes": reuse_pinned_hashes,
        "promotion_status": "CANDIDATE_CURRENT_NOT_CANONICAL",
    }
    failures = [
        field for field, wanted in expected.items() if receipt.get(field) != wanted
    ]
    if failures:
        raise ContractError(f"{label} identity failed: {failures}")
    require_exact_true_gates(receipt.get("gates"), WORKING_POINT_GATES, label)
    population_gates = receipt.get("population_certificate_gates")
    if not isinstance(population_gates, dict) or not population_gates:
        raise ContractError(
            f"{label}.population_certificate_gates must be a non-empty object"
        )
    if any(value is not True for value in population_gates.values()):
        raise ContractError(f"{label}.population_certificate_gates is not all true")
    exact_tolerance = require_finite_number(
        receipt.get("exact_max_efficiency_error"),
        f"{label}.exact_max_efficiency_error",
    )
    if exact_tolerance < 0.0 or exact_tolerance > MAX_WP_EXACT_EFFICIENCY_ERROR:
        raise ContractError(
            f"{label}.exact_max_efficiency_error={exact_tolerance}, "
            f"frozen maximum={MAX_WP_EXACT_EFFICIENCY_ERROR}"
        )
    rows = validate_exact_thresholds(
        system,
        receipt.get("exact_thresholds"),
        f"{label}.exact_thresholds",
        max_efficiency_error=exact_tolerance,
    )
    validate_runtime_surfaces(
        system, receipt.get("runtime_surfaces"), f"{label}.runtime_surfaces"
    )
    provenance = receipt.get("provenance")
    if not isinstance(provenance, dict):
        raise ContractError(f"{label}.provenance must be an object")
    for registry_key, receipt_key in (
        ("xgboost_sha256", "model_xgb_sha256"),
        ("metadata_sha256", "model_metadata_sha256"),
    ):
        if provenance.get(receipt_key) != model.get(registry_key):
            raise ContractError(
                f"{label}.{receipt_key} does not match registry model"
            )
    return rows


def validate_registry(
    view: str,
    record: dict,
    public_commit: str,
    code_sha256: str,
) -> dict:
    path, registry, observed_sha = verify_registered_json(
        record, f"{view} paired registry"
    )
    failures = []
    if registry.get("schema") != PAIRED_REGISTRY_SCHEMA:
        failures.append(f"schema={registry.get('schema')!r}")
    if registry.get("status") != "READY":
        failures.append(f"status={registry.get('status')!r}")
    if registry.get("promotion_status") != "CANDIDATE_CURRENT_NOT_CANONICAL":
        failures.append(f"promotion_status={registry.get('promotion_status')!r}")
    if registry.get("shower_definition") != view:
        failures.append("shower_definition mismatch")
    if registry.get("shower_semantic_sha256") != expected_view_semantic(view):
        failures.append("shower_semantic_sha256 mismatch")
    if registry.get("public_commit") != public_commit:
        failures.append("public_commit mismatch")
    if registry.get("code_sha256") != code_sha256:
        failures.append("code_sha256 mismatch")
    if registry.get("model_count") != 2:
        failures.append("model_count must equal 2")
    if registry.get("systems") != list(SYSTEMS):
        failures.append("systems must be ordered exactly as pp, auau")
    models = registry.get("models")
    if not isinstance(models, list):
        failures.append("models must be a list")
        models = []
    model_by_system = {}
    for model in models:
        if not isinstance(model, dict):
            failures.append("model entries must be objects")
            continue
        system = model.get("system")
        if system in model_by_system:
            failures.append(f"duplicate model system={system!r}")
        model_by_system[system] = model
    if set(model_by_system) != set(SYSTEMS):
        failures.append("models must contain exactly one pp and one auau model")
    if [
        model.get("system") for model in models if isinstance(model, dict)
    ] != list(SYSTEMS):
        failures.append("models must be ordered exactly as pp, auau")
    recorded_semantic = registry.get("registry_semantic_sha256")
    registry_without_semantic = dict(registry)
    registry_without_semantic.pop("registry_semantic_sha256", None)
    if recorded_semantic != canonical_json_sha256(registry_without_semantic):
        failures.append("registry_semantic_sha256 mismatch")
    if failures:
        raise ContractError(f"{view} paired registry failed: {failures}")

    normalized_models = []
    for system in SYSTEMS:
        model = require_exact_keys(
            model_by_system[system],
            MODEL_REGISTRY_KEYS,
            f"{view}/{system} registry model",
        )
        model_origin = expected_model_origin(system, view)
        expected_identity = {
            "system": system,
            "status": "READY",
            "model_origin": model_origin,
            "shower_definition": view,
            "shower_semantic_sha256": expected_view_semantic(view),
            "feature_order": list(FEATURES_BY_SYSTEM[system]),
            "model_domain_gev": list(MODEL_DOMAIN_GEV),
            "working_point_axis": WORKING_POINT_AXIS_BY_SYSTEM[system],
            "working_point_bin_edges": list(
                WORKING_POINT_EDGES_BY_SYSTEM[system]
            ),
            "threshold_authority": WORKING_POINT_THRESHOLD_AUTHORITY,
        }
        identity_failures = [
            field
            for field, expected in expected_identity.items()
            if model.get(field) != expected
        ]
        if identity_failures:
            raise ContractError(
                f"{view}/{system} registry model identity failed: "
                f"{identity_failures}"
            )
        reuse_pins = model.get("reuse_pinned_hashes")
        if model_origin.startswith("REUSED_"):
            if not isinstance(reuse_pins, dict) or not valid_reuse_pinned_hashes(
                reuse_pins, system, view
            ):
                raise ContractError(
                    f"{view}/{system} reused-model pins do not match frozen authority"
                )
        elif reuse_pins is not None:
            raise ContractError(
                f"{view}/{system} trained model unexpectedly carries reuse pins"
            )
        validate_exact_thresholds(
            system,
            model.get("exact_working_points"),
            f"{view}/{system} registry exact_working_points",
        )
        validate_runtime_surfaces(
            system,
            model.get("runtime_surfaces"),
            f"{view}/{system} registry runtime_surfaces",
        )

        model_artifacts = require_exact_keys(
            model.get("model"),
            MODEL_ARTIFACT_KEYS,
            f"{view}/{system} model artifacts",
        )
        for path_key, sha_key in (
            ("xgboost", "xgboost_sha256"),
            ("tmva", "tmva_sha256"),
            ("metadata", "metadata_sha256"),
        ):
            verify_registered_file(
                model_artifacts[path_key],
                model_artifacts[sha_key],
                f"{view}/{system} model {path_key}",
            )

        certificates = require_exact_keys(
            model.get("certificates"),
            MODEL_CERTIFICATE_KEYS,
            f"{view}/{system} model certificates",
        )
        for path_key, sha_key in MODEL_CERTIFICATE_ARTIFACT_PAIRS:
            verify_registered_file(
                certificates[path_key],
                certificates[sha_key],
                f"{view}/{system} certificate {path_key}",
            )

        _, extraction_receipt, _ = verify_registered_json(
            {
                "path": certificates["extraction"],
                "sha256": certificates["extraction_sha256"],
            },
            f"{view}/{system} extraction receipt",
        )
        extraction_checks = full_extraction_authority_checks(
            extraction_receipt, system, view
        )
        if not all(extraction_checks.values()):
            failed = sorted(
                key for key, value in extraction_checks.items() if not value
            )
            raise ContractError(
                f"{view}/{system} extraction receipt authority failed: {failed}"
            )
        extraction_binding = extraction_authority_binding(extraction_receipt)

        _, model_receipt, _ = verify_registered_json(
            {
                "path": certificates["model_receipt"],
                "sha256": certificates["model_receipt_sha256"],
            },
            f"{view}/{system} model completion receipt",
        )
        model_receipt_checks = model_artifact_receipt_checks(
            model_receipt,
            system,
            view,
            model_xgb_sha256=model_artifacts["xgboost_sha256"],
            model_metadata_sha256=model_artifacts["metadata_sha256"],
        )
        if not all(model_receipt_checks.values()):
            failed = sorted(
                key for key, value in model_receipt_checks.items() if not value
            )
            raise ContractError(
                f"{view}/{system} model completion receipt failed: {failed}"
            )
        if model_receipt.get("extraction_authority") != extraction_binding:
            raise ContractError(
                f"{view}/{system} model completion receipt extraction authority "
                "does not match the registered extraction"
            )

        _, validation, validation_sha256 = verify_registered_json(
            {
                "path": certificates["validation"],
                "sha256": certificates["validation_sha256"],
            },
            f"{view}/{system} validation receipt",
        )
        _, working_points, working_points_sha256 = verify_registered_json(
            {
                "path": certificates["working_points"],
                "sha256": certificates["working_points_sha256"],
            },
            f"{view}/{system} working-point receipt",
        )
        validate_validation_receipt(
            system,
            view,
            validation,
            model_artifacts,
            model_origin,
            reuse_pins,
            f"{view}/{system} validation receipt",
        )
        exact_thresholds = validate_working_point_receipt(
            system,
            view,
            working_points,
            model_artifacts,
            model_origin,
            reuse_pins,
            f"{view}/{system} working-point receipt",
        )
        if validation.get("extraction_authority") != extraction_binding:
            raise ContractError(
                f"{view}/{system} validation extraction authority does not "
                "match the registered extraction"
            )
        if working_points.get("extraction_authority") != extraction_binding:
            raise ContractError(
                f"{view}/{system} working-point extraction authority does not "
                "match the registered extraction"
            )

        expected_registry_bindings = {
            "exact_working_points": exact_thresholds,
            "working_point_axis": working_points["axis"],
            "working_point_bin_edges": working_points["bin_edges"],
            "threshold_authority": working_points["threshold_authority"],
            "runtime_surfaces": working_points["runtime_surfaces"],
        }
        binding_failures = [
            field
            for field, expected in expected_registry_bindings.items()
            if model.get(field) != expected
        ]
        if binding_failures:
            raise ContractError(
                f"{view}/{system} registry/working-point receipt binding failed: "
                f"{binding_failures}"
            )
        validation_provenance = validation["provenance"]
        working_point_provenance = working_points["provenance"]
        path_and_hash_bindings = (
            (
                validation_provenance,
                "model_xgb",
                "model_xgb_sha256",
                model_artifacts,
                "xgboost",
                "xgboost_sha256",
            ),
            (
                validation_provenance,
                "model_tmva",
                "model_tmva_sha256",
                model_artifacts,
                "tmva",
                "tmva_sha256",
            ),
            (
                validation_provenance,
                "model_metadata",
                "model_metadata_sha256",
                model_artifacts,
                "metadata",
                "metadata_sha256",
            ),
            (
                validation_provenance,
                "model_receipt",
                "model_receipt_sha256",
                certificates,
                "model_receipt",
                "model_receipt_sha256",
            ),
            (
                validation_provenance,
                "holdout",
                "holdout_sha256",
                certificates,
                "holdout",
                "holdout_sha256",
            ),
            (
                validation_provenance,
                "holdout_certificate",
                "holdout_certificate_sha256",
                certificates,
                "holdout_certificate",
                "holdout_certificate_sha256",
            ),
            (
                validation_provenance,
                "extraction_audit",
                "extraction_audit_sha256",
                certificates,
                "extraction",
                "extraction_sha256",
            ),
            (
                working_point_provenance,
                "score_sample",
                "score_sample_sha256",
                certificates,
                "working_point_score_sample",
                "working_point_score_sample_sha256",
            ),
            (
                working_point_provenance,
                "population_certificate",
                "population_certificate_sha256",
                certificates,
                "working_point_population_certificate",
                "working_point_population_certificate_sha256",
            ),
            (
                working_point_provenance,
                "model_xgb",
                "model_xgb_sha256",
                model_artifacts,
                "xgboost",
                "xgboost_sha256",
            ),
            (
                working_point_provenance,
                "model_metadata",
                "model_metadata_sha256",
                model_artifacts,
                "metadata",
                "metadata_sha256",
            ),
            (
                working_point_provenance,
                "model_receipt",
                "model_receipt_sha256",
                certificates,
                "model_receipt",
                "model_receipt_sha256",
            ),
        )
        for (
            receipt_payload,
            receipt_path_key,
            receipt_sha_key,
            registry_payload,
            registry_path_key,
            registry_sha_key,
        ) in path_and_hash_bindings:
            if (
                receipt_payload.get(receipt_path_key)
                != registry_payload.get(registry_path_key)
                or receipt_payload.get(receipt_sha_key)
                != registry_payload.get(registry_sha_key)
            ):
                raise ContractError(
                    f"{view}/{system} {receipt_path_key} path/hash binding "
                    "does not match the paired registry"
                )

        normalized_models.append(
            {
                "system": system,
                "model_origin": model_origin,
                "model_xgboost_sha256": model_artifacts["xgboost_sha256"],
                "model_tmva_sha256": model_artifacts["tmva_sha256"],
                "model_metadata_sha256": model_artifacts["metadata_sha256"],
                "validation_certificate_sha256": validation_sha256,
                "working_points_certificate_sha256": working_points_sha256,
                "exact_working_point_count": len(exact_thresholds),
            }
        )
    return {
        "view": view,
        "path": str(path),
        "sha256": observed_sha,
        "registry_semantic_sha256": recorded_semantic,
        "shower_semantic_sha256": expected_view_semantic(view),
        "systems": list(SYSTEMS),
        "models": normalized_models,
    }


def validate_source_witnesses(system: str, view: str, certificate: dict) -> list[dict]:
    witnesses = require_exact_key_inventory(
        certificate.get("source_witnesses"),
        REQUIRED_SOURCES_BY_SYSTEM[system],
        lambda record: record.get("source"),
        f"{system}/{view} source witnesses",
    )
    normalized = []
    for source in REQUIRED_SOURCES_BY_SYSTEM[system]:
        witness = witnesses[source]
        if witness.get("system") != system:
            raise ContractError(f"{system}/{view}/{source} system mismatch")
        if witness.get("shower_definition") != view:
            raise ContractError(
                f"{system}/{view}/{source} shower_definition mismatch"
            )
        hashes = {
            field: require_sha256(
                witness.get(field), f"{system}/{view}/{source}.{field}"
            )
            for field in (
                "source_manifest_sha256",
                "direct_artifact_sha256",
                "writer_artifact_sha256",
                "cache_receipt_sha256",
            )
        }
        require_exact_true_gates(
            witness.get("gates"),
            SOURCE_WITNESS_GATES,
            f"{system}/{view}/{source}",
        )
        normalized.append(
            {
                "system": system,
                "view": view,
                "source": source,
                **hashes,
            }
        )
    return normalized


def validate_replay_certificate(
    system: str,
    view: str,
    record: dict,
    registry_sha256: str,
    public_commit: str,
    code_sha256: str,
) -> dict:
    path, certificate, observed_sha = verify_registered_json(
        record, f"{system}/{view} replay certificate"
    )
    failures = []
    if certificate.get("schema") != REPLAY_CERTIFICATE_SCHEMA:
        failures.append(f"schema={certificate.get('schema')!r}")
    if certificate.get("status") != "PASS":
        failures.append(f"status={certificate.get('status')!r}")
    if certificate.get("promotion_status") != "CANDIDATE_CURRENT_NOT_CANONICAL":
        failures.append(f"promotion_status={certificate.get('promotion_status')!r}")
    if certificate.get("system") != system:
        failures.append("system mismatch")
    if certificate.get("shower_definition") != view:
        failures.append("shower_definition mismatch")
    if certificate.get("shower_semantic_sha256") != expected_view_semantic(view):
        failures.append("shower_semantic_sha256 mismatch")
    if certificate.get("paired_registry_sha256") != registry_sha256:
        failures.append("paired_registry_sha256 mismatch")
    if certificate.get("public_commit") != public_commit:
        failures.append("public_commit mismatch")
    if certificate.get("code_sha256") != code_sha256:
        failures.append("code_sha256 mismatch")
    if failures:
        raise ContractError(
            f"{system}/{view} replay certificate failed: {failures}"
        )
    require_exact_true_gates(
        certificate.get("gates"),
        REPLAY_CERTIFICATE_GATES,
        f"{system}/{view} replay certificate",
    )
    witnesses = validate_source_witnesses(system, view, certificate)
    return {
        "system": system,
        "view": view,
        "path": str(path),
        "sha256": observed_sha,
        "paired_registry_sha256": registry_sha256,
        "source_witness_count": len(witnesses),
        "source_witnesses": witnesses,
    }


def build_certificate(manifest_path: Path) -> dict:
    manifest = read_json(manifest_path, "science-freeze aggregate manifest")
    if manifest.get("schema") != SCIENCE_MANIFEST_SCHEMA:
        raise ContractError(
            f"manifest.schema={manifest.get('schema')!r}, "
            f"expected={SCIENCE_MANIFEST_SCHEMA!r}"
        )
    public_commit = require_git_commit(manifest.get("public_commit"))
    code_sha256 = require_sha256(manifest.get("code_sha256"), "code_sha256")
    registry_records = require_exact_key_inventory(
        manifest.get("paired_registries"),
        REQUIRED_VIEW_KEYS,
        lambda record: record.get("view"),
        "paired registries",
    )
    certificate_records = require_exact_key_inventory(
        manifest.get("replay_certificates"),
        REQUIRED_SYSTEM_VIEW_KEYS,
        lambda record: (record.get("system"), record.get("view")),
        "replay certificates",
    )

    registries = {}
    normalized_registries = []
    for view in REQUIRED_VIEW_KEYS:
        normalized = validate_registry(
            view,
            registry_records[view],
            public_commit,
            code_sha256,
        )
        registries[view] = normalized
        normalized_registries.append(normalized)

    normalized_certificates = []
    training_source_bindings: dict[tuple[str, str], dict[str, str]] = {}
    for system, view in REQUIRED_SYSTEM_VIEW_KEYS:
        normalized = validate_replay_certificate(
            system,
            view,
            certificate_records[(system, view)],
            registries[view]["sha256"],
            public_commit,
            code_sha256,
        )
        normalized_certificates.append(normalized)
        for witness in normalized["source_witnesses"]:
            key = (system, witness["source"])
            stable = {
                "source_manifest_sha256": witness["source_manifest_sha256"],
                "direct_artifact_sha256": witness["direct_artifact_sha256"],
                "writer_artifact_sha256": witness["writer_artifact_sha256"],
            }
            previous = training_source_bindings.setdefault(key, stable)
            if previous != stable:
                raise ContractError(
                    f"{system}/{witness['source']} artifact/source binding "
                    "changes across shower views"
                )

    if set(training_source_bindings) != set(REQUIRED_SOURCE_KEYS):
        raise ContractError(
            "aggregate training-source binding inventory is incomplete"
        )

    payload = {
        "schema": SCIENCE_CERTIFICATE_SCHEMA,
        "status": "PASS",
        "gate": "P5A-S",
        "promotion_status": "CANDIDATE_CURRENT_NOT_CANONICAL",
        "public_commit": public_commit,
        "code_sha256": code_sha256,
        "manifest": str(manifest_path),
        "manifest_sha256": hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
        "view_count": len(REQUIRED_VIEW_KEYS),
        "system_count": len(SYSTEMS),
        "paired_registry_count": len(normalized_registries),
        "replay_certificate_count": len(normalized_certificates),
        "training_source_view_witness_count": sum(
            item["source_witness_count"] for item in normalized_certificates
        ),
        "training_source_count": len(REQUIRED_SOURCE_KEYS),
        "paired_registries": normalized_registries,
        "replay_certificates": normalized_certificates,
        "training_source_bindings": [
            {
                "system": system,
                "source": source,
                **training_source_bindings[(system, source)],
            }
            for system, source in REQUIRED_SOURCE_KEYS
        ],
        "gates": {
            "all_seven_paired_registries_ready": True,
            "all_fourteen_view_system_certificates_pass": True,
            "training_source_complete_direct_writer_cache_witnesses": True,
            "cross_view_training_source_and_artifact_bindings_stable": True,
            "exact_wp70_wp80_wp90_present": True,
            "no_automatic_canonical_promotion": True,
        },
        "boundaries": [
            "This certificate freezes THE-134 P5A-S scientific authority only.",
            "THE-121/THE-122 remain closed pending THE-131 P5B and joint P5C.",
            "No model or physics output is promoted to CANONICAL.",
        ],
    }
    return semantic_receipt(payload, "aggregate_semantic_sha256")


def main() -> int:
    args = parse_args()
    try:
        payload = build_certificate(args.manifest)
    except ContractError as exc:
        print(f"THE-134 aggregate certification failed: {exc}", file=sys.stderr)
        return 2
    args.json_out.parent.mkdir(parents=True, exist_ok=True)
    args.json_out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(args.json_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
