#!/usr/bin/env python3
"""Detector-neutral identities for THE-134 aggregate release certificates.

This module deliberately contains no ROOT, scheduler, or filesystem mutation
logic. It defines the complete view/system/training-source Cartesian products
required by the THE-134 science-freeze and storage gates so both receipt
builders fail closed on omissions and duplicates.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Iterable

_HERE = Path(__file__).resolve()
_CONTRACTS = _HERE.parents[1] / "contracts"
sys.path.insert(0, str(_CONTRACTS))
from the134_h70_contract import (  # noqa: E402
    AUAU_CENTRALITY_EDGES,
    ALL_SHOWER_VIEWS,
    BACKGROUND_SOURCES,
    FEATURES_BY_SYSTEM,
    MAX_RUNTIME_SCORE_ABS_DIFFERENCE,
    MAX_WP_EXACT_EFFICIENCY_ERROR,
    MODEL_DOMAIN_GEV,
    READER_ET_EDGES,
    SIGNAL_SOURCES,
    TARGET_EFFICIENCIES,
    canonical_json_sha256,
    expected_model_origin,
    is_sha256,
    sha256_file,
    shower_semantic_sha256,
)


SYSTEMS = ("pp", "auau")
REQUIRED_SOURCES_BY_SYSTEM = {
    system: tuple(SIGNAL_SOURCES[system]) + tuple(BACKGROUND_SOURCES[system])
    for system in SYSTEMS
}
REQUIRED_SOURCE_KEYS = tuple(
    (system, source)
    for system in SYSTEMS
    for source in REQUIRED_SOURCES_BY_SYSTEM[system]
)
REQUIRED_VIEW_KEYS = tuple(ALL_SHOWER_VIEWS)
REQUIRED_SYSTEM_VIEW_KEYS = tuple(
    (system, view) for system in SYSTEMS for view in ALL_SHOWER_VIEWS
)

SCIENCE_MANIFEST_SCHEMA = "THE134_SCIENCE_FREEZE_AGGREGATE_MANIFEST_V1"
PAIRED_REGISTRY_SCHEMA = "THE134_FACTORIAL_VIEW_PAIRED_MODEL_REGISTRY_V1"
REPLAY_CERTIFICATE_SCHEMA = "THE134_FACTORIAL_VIEW_REPLAY_CERTIFICATE_V1"
SCIENCE_CERTIFICATE_SCHEMA = (
    "THE134_FACTORIAL_SCIENCE_FREEZE_AGGREGATE_CERTIFICATE_V1"
)
STORAGE_MANIFEST_SCHEMA = "THE134_STORAGE_QUOTA_PROJECTION_MANIFEST_V1"
STORAGE_CERTIFICATE_SCHEMA = "THE134_STORAGE_QUOTA_PROJECTION_CERTIFICATE_V1"
MODEL_VALIDATION_SCHEMA = "THE134_FACTORIAL_VIEW_MODEL_VALIDATION_V1"
WORKING_POINTS_SCHEMA = "THE134_FACTORIAL_VIEW_WEIGHTED_WORKING_POINTS_V1"

WORKING_POINT_AXIS_BY_SYSTEM = {
    "pp": "cluster_Et",
    "auau": "centrality",
}
WORKING_POINT_EDGES_BY_SYSTEM = {
    "pp": tuple(float(value) for value in READER_ET_EDGES),
    "auau": tuple(float(value) for value in AUAU_CENTRALITY_EDGES),
}
WORKING_POINT_LABEL_BY_TARGET = {
    0.90: "WP90",
    0.80: "WP80",
    0.70: "WP70",
}
WORKING_POINT_THRESHOLD_AUTHORITY = (
    "exact binned weighted holdout signal quantiles"
)
WORKING_POINT_COMPARISON_OPERATOR = "score > threshold"
WORKING_POINT_TIE_POLICY = (
    "stable mergesort weighted quantile; score ties are never split; exact "
    "binned authority fails if strict score>threshold efficiency misses "
    "target tolerance"
)
WORKING_POINT_THRESHOLD_ROW_KEYS = (
    "wp",
    "target_signal_efficiency",
    "bin_lo",
    "bin_hi",
    "bin_center",
    "threshold",
    "achieved_signal_efficiency",
    "achieved_minus_target",
    "abs_efficiency_error",
    "signal_weight_fraction_tied_at_threshold",
    "background_acceptance",
    "signal_rows",
    "background_rows",
)

MODEL_VALIDATION_GATES = (
    "candidate_metadata",
    "control_metadata",
    "control_view_authority",
    "extraction_certificate",
    "exact_event_group_holdout_certificate",
    "feature_order",
    "control_feature_order",
    "control_view_identity",
    "source_role_label_closure",
    "finite_inputs_scores_positive_weights",
    "cached_score_reconstruction",
    "python_tmva_runtime_parity",
    "no_material_overall_auc_regression",
    "et_bins_populated",
    "eta_bins_populated",
    "et_resolved_stability",
    "eta_resolved_stability",
    "overtraining_gate",
)
MODEL_VALIDATION_KEYS = (
    "schema",
    "status",
    "promotion_status",
    "system",
    "shower_definition",
    "shower_semantic_sha256",
    "feature_order",
    "model_domain_gev",
    "model_origin",
    "reuse_pinned_hashes",
    "extraction_authority",
    "critical_gates",
    "metadata_gates",
    "control_metadata_gates",
    "control_view_gates",
    "extraction_gates",
    "holdout_certificate_gates",
    "runtime_parity",
    "source_label_closure",
    "performance",
    "provenance",
)
MODEL_VALIDATION_DETAIL_GATE_FIELDS = (
    "metadata_gates",
    "control_metadata_gates",
    "control_view_gates",
    "extraction_gates",
    "holdout_certificate_gates",
)
WORKING_POINT_GATES = (
    "model_metadata_identity",
    "certified_wp_population",
    "holdout_feature_order",
    "finite_positive_weights",
    "both_classes_in_every_bin",
    "all_exact_thresholds_finite",
    "exact_weighted_efficiency_within_tolerance",
)
WORKING_POINT_RECEIPT_KEYS = (
    "schema",
    "status",
    "promotion_status",
    "system",
    "shower_definition",
    "shower_semantic_sha256",
    "model_domain_gev",
    "feature_order",
    "model_origin",
    "reuse_pinned_hashes",
    "extraction_authority",
    "axis",
    "bin_edges",
    "targets",
    "threshold_authority",
    "comparison_operator",
    "tie_policy",
    "exact_max_efficiency_error",
    "gates",
    "population_certificate_gates",
    "exact_thresholds",
    "runtime_surfaces",
    "provenance",
)
MODEL_REGISTRY_KEYS = (
    "system",
    "status",
    "model_origin",
    "reuse_pinned_hashes",
    "shower_definition",
    "shower_semantic_sha256",
    "feature_order",
    "model_domain_gev",
    "model",
    "exact_working_points",
    "working_point_axis",
    "working_point_bin_edges",
    "threshold_authority",
    "runtime_surfaces",
    "certificates",
)
MODEL_ARTIFACT_KEYS = (
    "xgboost",
    "xgboost_sha256",
    "tmva",
    "tmva_sha256",
    "metadata",
    "metadata_sha256",
)
MODEL_CERTIFICATE_KEYS = (
    "extraction",
    "extraction_sha256",
    "validation",
    "validation_sha256",
    "working_points",
    "working_points_sha256",
    "model_receipt",
    "model_receipt_sha256",
    "holdout",
    "holdout_sha256",
    "holdout_certificate",
    "holdout_certificate_sha256",
    "working_point_score_sample",
    "working_point_score_sample_sha256",
    "working_point_population_certificate",
    "working_point_population_certificate_sha256",
)
MODEL_CERTIFICATE_ARTIFACT_PAIRS = (
    ("extraction", "extraction_sha256"),
    ("validation", "validation_sha256"),
    ("working_points", "working_points_sha256"),
    ("model_receipt", "model_receipt_sha256"),
    ("holdout", "holdout_sha256"),
    ("holdout_certificate", "holdout_certificate_sha256"),
    ("working_point_score_sample", "working_point_score_sample_sha256"),
    (
        "working_point_population_certificate",
        "working_point_population_certificate_sha256",
    ),
)

REPLAY_CERTIFICATE_GATES = (
    "artifact_health",
    "direct_writer_histogram_sumw2_neutrality",
    "writer_cache_replay_closure",
    "normalized_tree_identity_population_closure",
    "model_input_runtime_score_parity",
    "below15_diagnostic_null_wp_safety",
    "source_weight_once_closure",
    "response_accounting_closure",
)

SOURCE_WITNESS_GATES = (
    "direct_health",
    "writer_health",
    "cache_health",
    "direct_writer_histogram_sumw2_neutrality",
    "writer_cache_replay_closure",
    "normalized_tree_identity_population_closure",
    "model_input_runtime_score_parity",
    "below15_diagnostic_null_wp_safety",
    "source_weight_once_closure",
    "response_accounting_closure",
)

SCIENCE_CERTIFICATE_GATES = (
    "all_seven_paired_registries_ready",
    "all_fourteen_view_system_certificates_pass",
    "training_source_complete_direct_writer_cache_witnesses",
    "cross_view_training_source_and_artifact_bindings_stable",
    "exact_wp70_wp80_wp90_present",
    "no_automatic_canonical_promotion",
)

SCIENCE_CERTIFICATE_KEYS = (
    "schema",
    "status",
    "gate",
    "promotion_status",
    "public_commit",
    "code_sha256",
    "manifest",
    "manifest_sha256",
    "view_count",
    "system_count",
    "paired_registry_count",
    "replay_certificate_count",
    "training_source_view_witness_count",
    "training_source_count",
    "paired_registries",
    "replay_certificates",
    "training_source_bindings",
    "gates",
    "boundaries",
    "aggregate_semantic_sha256",
)

MINIMUM_FREE_FRACTION = 0.20
SIMULTANEOUS_MERGE_WORKSPACE_COUNT = 2


class ContractError(ValueError):
    """Raised when a THE-134 aggregate input violates its frozen contract."""


def read_json(path: Path, label: str) -> dict:
    if not path.is_file():
        raise ContractError(f"missing {label} JSON: {path}")
    try:
        payload = json.loads(path.read_text())
    except json.JSONDecodeError as exc:
        raise ContractError(f"malformed {label} JSON {path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise ContractError(f"{label} JSON must contain one object: {path}")
    return payload


def require_sha256(value: object, label: str) -> str:
    if not is_sha256(value):
        raise ContractError(f"{label} is not an exact 64-character SHA-256")
    return str(value).lower()


def require_git_commit(value: object, label: str = "public_commit") -> str:
    text = str(value).lower()
    if len(text) != 40 or any(
        character not in "0123456789abcdef" for character in text
    ):
        raise ContractError(f"{label} is not an exact 40-character Git commit")
    return text


def require_exact_true_gates(
    gates: object,
    expected: Iterable[str],
    label: str,
) -> None:
    expected_set = set(expected)
    if not isinstance(gates, dict):
        raise ContractError(f"{label} must be an object")
    observed_set = set(gates)
    if observed_set != expected_set:
        missing = sorted(expected_set - observed_set)
        extra = sorted(observed_set - expected_set)
        raise ContractError(
            f"{label} gate inventory mismatch: missing={missing} extra={extra}"
        )
    failed = sorted(key for key, value in gates.items() if value is not True)
    if failed:
        raise ContractError(f"{label} gates are not all true: {failed}")


def require_exact_key_inventory(
    records: object,
    expected_keys: Iterable[object],
    key_function,
    label: str,
) -> dict:
    if not isinstance(records, list):
        raise ContractError(f"{label} must be a list")
    expected = tuple(expected_keys)
    keyed: dict[object, dict] = {}
    duplicates = []
    for record in records:
        if not isinstance(record, dict):
            raise ContractError(f"{label} records must be objects")
        key = key_function(record)
        if key in keyed:
            duplicates.append(key)
        keyed[key] = record
    if duplicates:
        raise ContractError(f"{label} contains duplicate keys: {sorted(duplicates)}")
    missing = set(expected) - set(keyed)
    extra = set(keyed) - set(expected)
    if missing or extra:
        raise ContractError(
            f"{label} inventory mismatch: missing={sorted(missing)} "
            f"extra={sorted(extra)}"
        )
    return keyed


def verify_registered_json(
    record: dict,
    label: str,
) -> tuple[Path, dict, str]:
    path = Path(str(record.get("path", "")))
    expected_sha = require_sha256(record.get("sha256"), f"{label}.sha256")
    payload = read_json(path, label)
    observed_sha = sha256_file(path)
    if observed_sha != expected_sha:
        raise ContractError(
            f"{label} SHA-256 mismatch: expected={expected_sha} "
            f"observed={observed_sha}"
        )
    return path, payload, observed_sha


def verify_registered_file(
    path_value: object,
    sha_value: object,
    label: str,
) -> tuple[Path, str]:
    """Rehash one producer-registered file at the aggregate boundary."""

    path = Path(str(path_value or ""))
    expected_sha = require_sha256(sha_value, f"{label}.sha256")
    if not path.is_file() or path.stat().st_size <= 0:
        raise ContractError(f"missing or empty {label}: {path}")
    observed_sha = sha256_file(path)
    if observed_sha != expected_sha:
        raise ContractError(
            f"{label} SHA-256 mismatch: expected={expected_sha} "
            f"observed={observed_sha}"
        )
    return path, observed_sha


def semantic_receipt(payload: dict, field: str) -> dict:
    result = dict(payload)
    result[field] = canonical_json_sha256(payload)
    return result


def verify_semantic_receipt(payload: dict, field: str, label: str) -> str:
    if field not in payload:
        raise ContractError(f"{label}.{field} is missing")
    semantic_payload = dict(payload)
    recorded = require_sha256(
        semantic_payload.pop(field), f"{label}.{field}"
    )
    expected = canonical_json_sha256(semantic_payload)
    if recorded != expected:
        raise ContractError(
            f"{label}.{field} mismatch: expected={expected} "
            f"observed={recorded}"
        )
    return recorded


def expected_view_semantic(view: str) -> str:
    if view not in ALL_SHOWER_VIEWS:
        raise ContractError(f"unknown shower view: {view!r}")
    return shower_semantic_sha256(view)
