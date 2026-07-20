#!/usr/bin/env python3
"""Fail-closed gate for the complete PPG12 stitched-purity estimator.

The gate is deliberately ROOT-independent.  Extraction code writes compact
JSON manifests containing the audited candidate-level summaries and weighted
histogram cells; this program validates and compares those manifests.  It
never submits jobs, merges ROOT files, changes current pointers, or promotes
artifacts.

Commands
--------
diagnose
    Compare an executable-PPG12 reference manifest with a RecoilJets manifest
    and always write structured diagnostic evidence.
admit
    Run the same comparison plus the additive-merge audit.  Only a complete
    pass emits ``admission_manifest.json``.
verify-production
    Require an unchanged admission, rerun the full comparison on the
    production reference/candidate pair, validate merge arithmetic and the
    historical-archive comparison, and emit ``production_gate_report.json``
    only on pass.

The manifest schemas are documented in
``agent_context/analysis_contracts/ppg12_stitched_purity_closure.yaml``.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import sys
from pathlib import Path
from typing import Any, Iterable


REPO = Path(__file__).resolve().parents[3]
DEFAULT_CONTRACT = (
    REPO
    / "agent_context/analysis_contracts/ppg12_stitched_purity_closure.yaml"
)
REQUIRED_HISTORICAL_SOURCE_ROLES = {
    "historical_purity",
    "historical_abcd",
    "historical_leakage",
    "candidate_purity",
    "candidate_inclusive",
    "candidate_photon",
}


class GateInputError(RuntimeError):
    """Input could not be parsed well enough to run the scientific gate."""

    def __init__(self, code: str, message: str):
        super().__init__(message)
        self.code = code


def _read_json(path: Path, label: str) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text())
    except FileNotFoundError as exc:
        raise GateInputError("manifest_invalid", f"missing {label}: {path}") from exc
    except json.JSONDecodeError as exc:
        raise GateInputError(
            "manifest_invalid", f"invalid JSON-compatible YAML in {label} {path}: {exc}"
        ) from exc
    if not isinstance(payload, dict):
        raise GateInputError("manifest_invalid", f"{label} must contain an object: {path}")
    return payload


def _canonical_bytes(payload: Any) -> bytes:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def _payload_sha256(payload: Any) -> str:
    return hashlib.sha256(_canonical_bytes(payload)).hexdigest()


def _canonical_groups(
    lane_id: str, groups: Any, group_size: int
) -> tuple[list[dict[str, Any]], int]:
    if not isinstance(groups, list):
        raise GateInputError("canary_underpopulated", f"{lane_id} groups must be a list")
    normalized: list[dict[str, Any]] = []
    seen_rows: set[str] = set()
    for group_index, row in enumerate(groups):
        if not isinstance(row, dict):
            raise GateInputError(
                "canary_underpopulated", f"{lane_id} group {group_index} is malformed"
            )
        event_rows = row.get("event_rows")
        if (
            not isinstance(event_rows, list)
            or len(event_rows) != group_size
            or not all(isinstance(item, str) and item for item in event_rows)
        ):
            raise GateInputError(
                "canary_underpopulated",
                f"{lane_id} group {group_index} is not an exact five-file group",
            )
        if any(item in seen_rows for item in event_rows):
            raise GateInputError(
                "canary_underpopulated", f"{lane_id} group inputs are duplicated"
            )
        seen_rows.update(event_rows)
        group_sha = _payload_sha256(
            {
                "lane_id": lane_id,
                "group_index": group_index,
                "event_rows": event_rows,
            }
        )
        canonical = {
            "group_index": group_index,
            "group_id": f"group-{group_index:05d}-{group_sha[:16]}",
            "event_rows": event_rows,
            "group_sha256": group_sha,
        }
        if row != canonical:
            raise GateInputError(
                "canary_underpopulated",
                f"{lane_id} group {group_index} identity/hash is stale",
            )
        normalized.append(canonical)
    return normalized, len(seen_rows)


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n")


def _write_csv(path: Path, rows: list[dict[str, Any]], preferred: Iterable[str]) -> None:
    fields = list(preferred)
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _failure(code: str, message: str, **context: Any) -> dict[str, Any]:
    row: dict[str, Any] = {"code": code, "message": message}
    if context:
        row["context"] = context
    return row


def _finite_number(value: Any) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(float(value))


def _count_or_zero(value: Any) -> int:
    return value if isinstance(value, int) and not isinstance(value, bool) and value >= 0 else 0


def _as_finite_vector(value: Any, label: str) -> list[float]:
    if not isinstance(value, list) or not value:
        raise GateInputError("manifest_invalid", f"{label} must be a non-empty list")
    if not all(_finite_number(item) for item in value):
        raise GateInputError("manifest_invalid", f"{label} contains a non-finite value")
    return [float(item) for item in value]


def _lane_id(family: str, sample: str, period: str, interaction: str) -> str:
    return f"{family}:{sample}:{period}:{interaction}"


def _expected_lanes(contract: dict[str, Any]) -> dict[str, dict[str, str]]:
    expected: dict[str, dict[str, str]] = {}
    for family, spec in contract["families"].items():
        for sample in spec["samples"]:
            for period in contract["periods"]:
                for interaction in contract["interactions"]:
                    lane_id = _lane_id(family, sample, period, interaction)
                    expected[lane_id] = {
                        "family": family,
                        "sample": sample,
                        "period": period,
                        "interaction": interaction,
                    }
    return expected


def _expected_candidate_parity_lanes(contract: dict[str, Any]) -> set[str]:
    """Return the exact photon-lane set covered by same-event comparisons."""
    return {
        _lane_id("photon", sample, period, interaction)
        for sample in contract["families"]["photon"]["samples"]
        for period in contract["periods"]
        for interaction in contract["interactions"]
    }


def _is_sha256(value: Any) -> bool:
    return (
        isinstance(value, str)
        and len(value) == 64
        and all(character in "0123456789abcdef" for character in value)
    )


def _canonical_assembler(contract: dict[str, Any]) -> Any:
    path = (REPO / contract["canonical_tools"]["assembler"]).resolve()
    spec = importlib.util.spec_from_file_location("ppg12_closure_assembler", path)
    if spec is None or spec.loader is None:
        raise GateInputError("manifest_invalid", f"cannot load canonical assembler: {path}")
    module = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(module)
    except Exception as exc:
        raise GateInputError(
            "manifest_invalid", f"cannot import canonical assembler {path}: {exc}"
        ) from exc
    return module


def _canonical_manifest_gate(
    path: Path,
    contract: dict[str, Any],
    roles: set[str],
    label: str,
    failures: list[dict[str, Any]],
) -> None:
    try:
        _canonical_assembler(contract).verify_manifest_evidence(path, contract, roles)
    except Exception as exc:
        failures.append(
            _failure(
                "manifest_invalid",
                f"{label} is not a canonical hash-bound assembly: {exc}",
            )
        )


def _index_lanes(
    manifest: dict[str, Any],
    label: str,
    contract: dict[str, Any],
    failures: list[dict[str, Any]],
) -> dict[str, dict[str, Any]]:
    lanes = manifest.get("lanes")
    if not isinstance(lanes, list):
        failures.append(_failure("manifest_invalid", f"{label}.lanes must be a list"))
        return {}
    indexed: dict[str, dict[str, Any]] = {}
    for index, lane in enumerate(lanes):
        if not isinstance(lane, dict):
            failures.append(_failure("manifest_invalid", f"{label}.lanes[{index}] is not an object"))
            continue
        keys = ("family", "sample", "period", "interaction")
        if any(not isinstance(lane.get(key), str) for key in keys):
            failures.append(_failure("manifest_invalid", f"{label}.lanes[{index}] has incomplete identity"))
            continue
        canonical = _lane_id(*(lane[key] for key in keys))
        declared = lane.get("lane_id")
        if declared != canonical:
            failures.append(
                _failure(
                    "source_contract_failed",
                    f"{label} lane_id does not match its identity",
                    declared=declared,
                    expected=canonical,
                )
            )
        if canonical in indexed:
            failures.append(_failure("source_contract_failed", f"duplicate lane in {label}: {canonical}"))
            continue
        indexed[canonical] = lane

    expected = set(_expected_lanes(contract))
    observed = set(indexed)
    for lane_id in sorted(expected - observed):
        failures.append(_failure("source_contract_failed", f"missing lane in {label}: {lane_id}"))
    for lane_id in sorted(observed - expected):
        failures.append(_failure("source_contract_failed", f"unexpected lane in {label}: {lane_id}"))
    return indexed


def _validate_manifest_header(
    manifest: dict[str, Any],
    label: str,
    allowed_roles: set[str],
    contract: dict[str, Any],
    failures: list[dict[str, Any]],
) -> None:
    if manifest.get("schema") != "ppg12-stitched-purity-manifest/v1":
        failures.append(_failure("manifest_invalid", f"unsupported {label} schema"))
    if manifest.get("role") not in allowed_roles:
        failures.append(
            _failure("manifest_invalid", f"invalid {label} role", role=manifest.get("role"))
        )
    if manifest.get("abcd_population") != contract["abcd_population"]:
        failures.append(
            _failure(
                "candidate_contract_failed" if label != "reference" else "reference_invalid",
                f"{label} must use unsuffixed inclusive A/B/C/D",
                observed=manifest.get("abcd_population"),
            )
        )
    scale = manifest.get("external_scale")
    if not _finite_number(scale) or abs(float(scale) - float(contract["external_scale"])) > 1e-12:
        failures.append(
            _failure(
                "candidate_contract_failed" if label != "reference" else "reference_invalid",
                f"{label} uses a forbidden external normalization",
                observed=scale,
            )
        )
    if manifest.get("random_seed") != contract["random_seed"]:
        failures.append(_failure("purity_solver_failed", f"{label} random seed is not frozen seed 42"))
    if manifest.get("toy_count") != contract["toy_count"]:
        failures.append(_failure("purity_solver_failed", f"{label} toy count is not frozen at 20000"))
    provenance = manifest.get("provenance")
    if not isinstance(provenance, dict):
        failures.append(_failure("manifest_invalid", f"{label}.provenance must be an object"))
        return
    for field in contract["frozen_provenance_fields"]:
        value = provenance.get(field)
        if not isinstance(value, str) or not value:
            failures.append(_failure("manifest_invalid", f"{label}.provenance.{field} is required"))


def _compare_headers(
    reference: dict[str, Any],
    candidate: dict[str, Any],
    contract: dict[str, Any],
    failures: list[dict[str, Any]],
) -> None:
    comparable = [field for field in contract["frozen_provenance_fields"] if field != "implementation_sha256"]
    ref = reference.get("provenance", {})
    cand = candidate.get("provenance", {})
    for field in comparable:
        if ref.get(field) != cand.get(field):
            failures.append(
                _failure(
                    "source_contract_failed",
                    f"reference/candidate provenance mismatch: {field}",
                    reference=ref.get(field),
                    candidate=cand.get(field),
                )
            )


def _feature_and_tag_gate(
    candidate: dict[str, Any], contract: dict[str, Any], failures: list[dict[str, Any]]
) -> dict[str, Any]:
    parity = candidate.get("candidate_parity")
    summary: dict[str, Any] = {}
    if not isinstance(parity, dict):
        failures.append(_failure("candidate_contract_failed", "candidate_parity summary is required"))
        return summary

    declared_tolerances = parity.get("tolerances")
    expected_tolerances = contract["candidate_parity_tolerances"]
    if not isinstance(declared_tolerances, dict) or set(declared_tolerances) != set(expected_tolerances):
        failures.append(
            _failure(
                "candidate_contract_failed",
                "candidate-parity tolerance contract is missing or incomplete",
                observed=declared_tolerances,
                expected=expected_tolerances,
            )
        )
    else:
        for name, expected in expected_tolerances.items():
            value = declared_tolerances.get(name)
            if not _finite_number(value) or float(value) != float(expected):
                failures.append(
                    _failure(
                        "candidate_contract_failed",
                        f"candidate-parity tolerance drift: {name}",
                        observed=value,
                        expected=expected,
                    )
                )

    expected_lane_ids = _expected_candidate_parity_lanes(contract)
    lane_rows = parity.get("lane_coverage")
    indexed_lane_rows: dict[str, dict[str, Any]] = {}
    if not isinstance(lane_rows, list):
        failures.append(
            _failure(
                "candidate_contract_failed",
                "candidate-parity lane_coverage must contain the exact 12 photon lanes",
            )
        )
    else:
        for index, row in enumerate(lane_rows):
            lane_id = row.get("lane_id") if isinstance(row, dict) else None
            if not isinstance(lane_id, str) or lane_id in indexed_lane_rows:
                failures.append(
                    _failure(
                        "candidate_contract_failed",
                        f"invalid or duplicate candidate-parity lane_coverage row {index}",
                        lane_id=lane_id,
                    )
                )
                continue
            indexed_lane_rows[lane_id] = row
        if set(indexed_lane_rows) != expected_lane_ids:
            failures.append(
                _failure(
                    "candidate_contract_failed",
                    "candidate-parity lane coverage is not the exact 12-lane photon contract",
                    missing=sorted(expected_lane_ids - set(indexed_lane_rows)),
                    unexpected=sorted(set(indexed_lane_rows) - expected_lane_ids),
                )
            )

    score_tol = float(contract["tolerances"]["score_abs"])
    for lane_id, row in sorted(indexed_lane_rows.items()):
        population_row = row.get("population", {})
        ref_count = population_row.get("reference_count") if isinstance(population_row, dict) else None
        cand_count = population_row.get("candidate_count") if isinstance(population_row, dict) else None
        unmatched_ref = population_row.get("unmatched_reference") if isinstance(population_row, dict) else None
        unmatched_cand = population_row.get("unmatched_candidate") if isinstance(population_row, dict) else None
        if (
            not isinstance(ref_count, int)
            or isinstance(ref_count, bool)
            or ref_count <= 0
            or cand_count != ref_count
            or unmatched_ref != 0
            or unmatched_cand != 0
        ):
            failures.append(
                _failure(
                    "candidate_contract_failed",
                    f"candidate populations do not match in photon lane {lane_id}",
                    population=population_row,
                )
            )
        feature_row = row.get("features", {})
        if (
            not isinstance(feature_row, dict)
            or not isinstance(feature_row.get("comparisons"), int)
            or feature_row.get("comparisons", 0) <= 0
            or feature_row.get("violations") != 0
            or not _finite_number(feature_row.get("max_normalized_delta"))
            or float(feature_row["max_normalized_delta"]) > 1.0
        ):
            failures.append(
                _failure(
                    "candidate_contract_failed",
                    f"feature tolerance failed in photon lane {lane_id}",
                    observed=feature_row,
                )
            )
        score_row = row.get("scores", {})
        if (
            not isinstance(score_row, dict)
            or not isinstance(score_row.get("comparisons"), int)
            or score_row.get("comparisons", 0) <= 0
            or score_row.get("mismatches") != 0
            or not _finite_number(score_row.get("max_abs_delta"))
            or float(score_row["max_abs_delta"]) > score_tol
        ):
            failures.append(
                _failure(
                    "candidate_contract_failed",
                    f"score tolerance failed in photon lane {lane_id}",
                    observed=score_row,
                )
            )
        route_row = row.get("model_route", {})
        if (
            not isinstance(route_row, dict)
            or not isinstance(route_row.get("comparisons"), int)
            or route_row.get("comparisons", 0) <= 0
            or route_row.get("mismatches") != 0
        ):
            failures.append(
                _failure(
                    "candidate_contract_failed",
                    f"model-route comparison failed in photon lane {lane_id}",
                    observed=route_row,
                )
            )
        tag_rows = row.get("tags", {})
        if not isinstance(tag_rows, dict) or set(tag_rows) != set(contract["required_tag_checks"]):
            failures.append(
                _failure(
                    "candidate_contract_failed",
                    f"tag coverage is incomplete in photon lane {lane_id}",
                )
            )
        else:
            for name in contract["required_tag_checks"]:
                tag_row = tag_rows[name]
                if (
                    not isinstance(tag_row, dict)
                    or not isinstance(tag_row.get("comparisons"), int)
                    or tag_row.get("comparisons", 0) <= 0
                    or tag_row.get("mismatches") != 0
                ):
                    failures.append(
                        _failure(
                            "candidate_contract_failed",
                            f"tag comparison failed in photon lane {lane_id}: {name}",
                            observed=tag_row,
                        )
                    )

    population = parity.get("population", {})
    ref_count = population.get("reference_count")
    cand_count = population.get("candidate_count")
    unmatched_ref = population.get("unmatched_reference")
    unmatched_cand = population.get("unmatched_candidate")
    if not all(isinstance(value, int) and value >= 0 for value in (ref_count, cand_count, unmatched_ref, unmatched_cand)):
        failures.append(_failure("candidate_contract_failed", "invalid candidate population summary"))
    elif ref_count == 0 or cand_count != ref_count or unmatched_ref or unmatched_cand:
        failures.append(
            _failure(
                "candidate_contract_failed",
                "candidate populations do not match one-to-one",
                reference_count=ref_count,
                candidate_count=cand_count,
                unmatched_reference=unmatched_ref,
                unmatched_candidate=unmatched_cand,
            )
        )
    summary["population"] = population

    feature_rows = parity.get("features")
    feature_map = {
        row.get("name"): row
        for row in feature_rows
        if isinstance(row, dict) and isinstance(row.get("name"), str)
    } if isinstance(feature_rows, list) else {}
    for feature in contract["required_features"]:
        row = feature_map.get(feature)
        if row is None:
            failures.append(_failure("candidate_contract_failed", f"missing feature comparison: {feature}"))
            continue
        comparisons = row.get("comparisons")
        violations = row.get("violations")
        max_normalized = row.get("max_normalized_delta")
        if not isinstance(comparisons, int) or comparisons <= 0:
            failures.append(_failure("candidate_contract_failed", f"feature {feature} has no comparisons"))
        if violations != 0 or not _finite_number(max_normalized) or float(max_normalized) > 1.0:
            failures.append(
                _failure(
                    "candidate_contract_failed",
                    f"feature tolerance failed: {feature}",
                    violations=violations,
                    max_normalized_delta=max_normalized,
                )
            )

    score_rows = parity.get("scores")
    score_map = {
        row.get("name"): row
        for row in score_rows
        if isinstance(row, dict) and isinstance(row.get("name"), str)
    } if isinstance(score_rows, list) else {}
    score_tol = float(contract["tolerances"]["score_abs"])
    for model in contract["required_score_models"]:
        row = score_map.get(model)
        if row is None:
            failures.append(_failure("candidate_contract_failed", f"missing score comparison: {model}"))
            continue
        if (
            not isinstance(row.get("comparisons"), int)
            or row.get("comparisons", 0) <= 0
            or row.get("mismatches") != 0
            or not _finite_number(row.get("max_abs_delta"))
            or float(row["max_abs_delta"]) > score_tol
        ):
            failures.append(_failure("candidate_contract_failed", f"score tolerance failed: {model}", **row))

    route = parity.get("model_route", {})
    if not isinstance(route.get("comparisons"), int) or route.get("comparisons", 0) <= 0 or route.get("mismatches") != 0:
        failures.append(_failure("candidate_contract_failed", "model-route comparison failed", **route))

    tags = parity.get("tags", {})
    for name in contract["required_tag_checks"]:
        row = tags.get(name, {}) if isinstance(tags, dict) else {}
        if not isinstance(row.get("comparisons"), int) or row.get("comparisons", 0) <= 0 or row.get("mismatches") != 0:
            failures.append(_failure("candidate_contract_failed", f"tag comparison failed: {name}", **row))

    if set(indexed_lane_rows) == expected_lane_ids:
        lane_population_rows = [
            row.get("population") if isinstance(row.get("population"), dict) else {}
            for row in indexed_lane_rows.values()
        ]
        lane_reference_total = sum(_count_or_zero(row.get("reference_count")) for row in lane_population_rows)
        lane_candidate_total = sum(_count_or_zero(row.get("candidate_count")) for row in lane_population_rows)
        lane_unmatched_reference = sum(_count_or_zero(row.get("unmatched_reference")) for row in lane_population_rows)
        lane_unmatched_candidate = sum(_count_or_zero(row.get("unmatched_candidate")) for row in lane_population_rows)
        if (
            ref_count != lane_reference_total
            or cand_count != lane_candidate_total
            or unmatched_ref != lane_unmatched_reference
            or unmatched_cand != lane_unmatched_candidate
        ):
            failures.append(
                _failure(
                    "candidate_contract_failed",
                    "global candidate population does not equal the exact 12-lane sum",
                    global_population=population,
                    lane_population={
                        "reference_count": lane_reference_total,
                        "candidate_count": lane_candidate_total,
                        "unmatched_reference": lane_unmatched_reference,
                        "unmatched_candidate": lane_unmatched_candidate,
                    },
                )
            )

        lane_feature_rows = [
            row.get("features") if isinstance(row.get("features"), dict) else {}
            for row in indexed_lane_rows.values()
        ]
        feature_comparison_total = sum(_count_or_zero(row.get("comparisons")) for row in lane_feature_rows)
        feature_violation_total = sum(_count_or_zero(row.get("violations")) for row in lane_feature_rows)
        global_feature_comparisons = sum(_count_or_zero(row.get("comparisons")) for row in feature_map.values())
        global_feature_violations = sum(_count_or_zero(row.get("violations")) for row in feature_map.values())
        if feature_comparison_total != global_feature_comparisons or feature_violation_total != global_feature_violations:
            failures.append(
                _failure(
                    "candidate_contract_failed",
                    "global feature totals do not equal the exact 12-lane sum",
                    lane_comparisons=feature_comparison_total,
                    global_comparisons=global_feature_comparisons,
                    lane_violations=feature_violation_total,
                    global_violations=global_feature_violations,
                )
            )

        lane_score_rows = [
            row.get("scores") if isinstance(row.get("scores"), dict) else {}
            for row in indexed_lane_rows.values()
        ]
        score_comparison_total = sum(_count_or_zero(row.get("comparisons")) for row in lane_score_rows)
        score_mismatch_total = sum(_count_or_zero(row.get("mismatches")) for row in lane_score_rows)
        global_score_comparisons = sum(_count_or_zero(row.get("comparisons")) for row in score_map.values())
        global_score_mismatches = sum(_count_or_zero(row.get("mismatches")) for row in score_map.values())
        if score_comparison_total != global_score_comparisons or score_mismatch_total != global_score_mismatches:
            failures.append(
                _failure(
                    "candidate_contract_failed",
                    "global score totals do not equal the exact 12-lane sum",
                    lane_comparisons=score_comparison_total,
                    global_comparisons=global_score_comparisons,
                    lane_mismatches=score_mismatch_total,
                    global_mismatches=global_score_mismatches,
                )
            )

        lane_route_rows = [
            row.get("model_route") if isinstance(row.get("model_route"), dict) else {}
            for row in indexed_lane_rows.values()
        ]
        lane_route_comparisons = sum(_count_or_zero(row.get("comparisons")) for row in lane_route_rows)
        lane_route_mismatches = sum(_count_or_zero(row.get("mismatches")) for row in lane_route_rows)
        if route.get("comparisons") != lane_route_comparisons or route.get("mismatches") != lane_route_mismatches:
            failures.append(
                _failure(
                    "candidate_contract_failed",
                    "global model-route totals do not equal the exact 12-lane sum",
                )
            )
        for name in contract["required_tag_checks"]:
            lane_tag_rows = [
                row.get("tags", {}).get(name, {})
                if isinstance(row.get("tags"), dict)
                and isinstance(row.get("tags", {}).get(name), dict)
                else {}
                for row in indexed_lane_rows.values()
            ]
            lane_comparisons = sum(_count_or_zero(row.get("comparisons")) for row in lane_tag_rows)
            lane_mismatches = sum(_count_or_zero(row.get("mismatches")) for row in lane_tag_rows)
            global_row = tags.get(name, {}) if isinstance(tags, dict) else {}
            if global_row.get("comparisons") != lane_comparisons or global_row.get("mismatches") != lane_mismatches:
                failures.append(
                    _failure(
                        "candidate_contract_failed",
                        f"global tag totals do not equal the exact 12-lane sum: {name}",
                    )
                )

    summary["feature_count"] = len(feature_map)
    summary["score_model_count"] = len(score_map)
    summary["candidate_parity_lane_count"] = len(indexed_lane_rows)
    summary["candidate_parity_tolerances"] = declared_tolerances
    summary["model_route"] = route
    summary["tags"] = tags
    return summary


def _histogram_vectors(hist: Any, label: str) -> tuple[list[float], list[float], list[float], list[float]]:
    if not isinstance(hist, dict):
        raise GateInputError("root_family_missing", f"missing histogram payload: {label}")
    edges = _as_finite_vector(hist.get("bin_edges"), f"{label}.bin_edges")
    sumw = _as_finite_vector(hist.get("sumw"), f"{label}.sumw")
    sumw2 = _as_finite_vector(hist.get("sumw2"), f"{label}.sumw2")
    fills = _as_finite_vector(hist.get("fills"), f"{label}.fills")
    if len(edges) != len(sumw) + 1 or len(sumw2) != len(sumw) or len(fills) != len(sumw):
        raise GateInputError("binning_mismatch", f"inconsistent vector lengths: {label}")
    if any(right <= left for left, right in zip(edges, edges[1:])):
        raise GateInputError("binning_mismatch", f"non-increasing bin edges: {label}")
    if any(value < 0 for value in sumw2) or any(value < 0 for value in fills):
        raise GateInputError("manifest_invalid", f"negative sumw2/fill count: {label}")
    return edges, sumw, sumw2, fills


def _cell_delta(candidate: float, reference: float, abs_floor: float) -> tuple[float, float, float]:
    absolute = candidate - reference
    denominator = max(abs(reference), abs_floor)
    relative = absolute / denominator
    ratio = candidate / reference if reference != 0.0 else (1.0 if candidate == 0.0 else math.inf)
    return absolute, relative, ratio


def _compare_lanes(
    reference_lanes: dict[str, dict[str, Any]],
    candidate_lanes: dict[str, dict[str, Any]],
    contract: dict[str, Any],
    failures: list[dict[str, Any]],
) -> tuple[list[dict[str, Any]], dict[str, dict[str, tuple[list[float], list[float]]]]]:
    rows: list[dict[str, Any]] = []
    aggregate: dict[str, dict[str, tuple[list[float], list[float]]]] = {"reference": {}, "candidate": {}}
    expected = _expected_lanes(contract)
    cell_tol = float(contract["tolerances"]["weighted_cell_relative"])
    abs_floor = float(contract["tolerances"]["weighted_cell_abs_floor"])

    aggregate_sums: dict[str, dict[str, list[float]]] = {"reference": {}, "candidate": {}}
    aggregate_edges: dict[str, dict[str, list[float]]] = {"reference": {}, "candidate": {}}

    for lane_id, identity in sorted(expected.items()):
        ref_lane = reference_lanes.get(lane_id)
        cand_lane = candidate_lanes.get(lane_id)
        if ref_lane is None or cand_lane is None:
            continue
        for field in contract["paired_lane_fields"]:
            if not isinstance(ref_lane.get(field), str) or not ref_lane.get(field):
                failures.append(_failure("manifest_invalid", f"reference {lane_id} missing {field}"))
            if not isinstance(cand_lane.get(field), str) or not cand_lane.get(field):
                failures.append(_failure("manifest_invalid", f"candidate {lane_id} missing {field}"))
            if ref_lane.get(field) != cand_lane.get(field):
                failures.append(
                    _failure(
                        "source_contract_failed",
                        f"lane contract mismatch for {lane_id}: {field}",
                        reference=ref_lane.get(field),
                        candidate=cand_lane.get(field),
                    )
                )
        for lane, label in ((ref_lane, "reference"), (cand_lane, "candidate")):
            group_contract = contract["admission_group_contract"]
            group_count = lane.get("group_count")
            group_size = lane.get("group_size")
            if (
                not isinstance(group_count, int)
                or isinstance(group_count, bool)
                or group_count < int(group_contract["minimum_group_count"])
                or group_size != int(group_contract["group_size"])
            ):
                failures.append(
                    _failure(
                        "canary_underpopulated",
                        f"{label} {lane_id} lacks complete five-file groups",
                    )
                )
            else:
                try:
                    groups, event_row_count = _canonical_groups(
                        lane_id, lane.get("groups"), group_size
                    )
                    if (
                        len(groups) != group_count
                        or lane.get("event_set_row_count") != event_row_count
                        or event_row_count != group_count * group_size
                        or lane.get("group_set_sha256") != _payload_sha256(groups)
                    ):
                        raise GateInputError(
                            "canary_underpopulated",
                            f"{label} {lane_id} group coverage/hash is stale",
                        )
                except GateInputError as exc:
                    failures.append(_failure(exc.code, str(exc)))
            lane_scale = lane.get("external_scale")
            if not _finite_number(lane_scale) or abs(float(lane_scale) - 1.0) > 1e-12:
                failures.append(
                    _failure(
                        "candidate_contract_failed" if label == "candidate" else "reference_invalid",
                        f"{label} {lane_id} uses a forbidden lane-level normalization",
                        observed=lane_scale,
                    )
                )

        required = contract["families"][identity["family"]]["required_observables"]
        ref_obs = ref_lane.get("observables", {})
        cand_obs = cand_lane.get("observables", {})
        if not isinstance(ref_obs, dict) or not isinstance(cand_obs, dict):
            failures.append(_failure("root_family_missing", f"missing observables for lane {lane_id}"))
            continue
        for observable in required:
            if observable not in ref_obs or observable not in cand_obs:
                failures.append(
                    _failure(
                        "root_family_missing",
                        f"missing required observable {lane_id}:{observable}",
                        reference_present=observable in ref_obs,
                        candidate_present=observable in cand_obs,
                    )
                )
                continue
            try:
                ref_edges, ref_sumw, ref_sumw2, ref_fills = _histogram_vectors(
                    ref_obs[observable], f"reference:{lane_id}:{observable}"
                )
                cand_edges, cand_sumw, cand_sumw2, cand_fills = _histogram_vectors(
                    cand_obs[observable], f"candidate:{lane_id}:{observable}"
                )
            except GateInputError as exc:
                failures.append(_failure(exc.code, str(exc)))
                continue
            if ref_edges != cand_edges:
                failures.append(_failure("binning_mismatch", f"bin edges differ for {lane_id}:{observable}"))
                continue
            category = "photon_leakage_failed" if identity["family"] == "photon" else "inclusive_abcd_failed"
            for index, (r_sumw, c_sumw, r_sumw2, c_sumw2, r_fills, c_fills) in enumerate(
                zip(ref_sumw, cand_sumw, ref_sumw2, cand_sumw2, ref_fills, cand_fills), start=1
            ):
                abs_delta, rel_delta, ratio = _cell_delta(c_sumw, r_sumw, abs_floor)
                _, sumw2_rel, _ = _cell_delta(c_sumw2, r_sumw2, abs_floor)
                fill_delta = c_fills - r_fills
                rows.append(
                    {
                        "lane_id": lane_id,
                        "family": identity["family"],
                        "sample": identity["sample"],
                        "period": identity["period"],
                        "interaction": identity["interaction"],
                        "observable": observable,
                        "bin": index,
                        "x_low": ref_edges[index - 1],
                        "x_high": ref_edges[index],
                        "reference_sumw": r_sumw,
                        "candidate_sumw": c_sumw,
                        "candidate_over_reference": ratio,
                        "absolute_delta": abs_delta,
                        "relative_delta": rel_delta,
                        "reference_sumw2": r_sumw2,
                        "candidate_sumw2": c_sumw2,
                        "sumw2_relative_delta": sumw2_rel,
                        "reference_fills": r_fills,
                        "candidate_fills": c_fills,
                        "fill_delta": fill_delta,
                    }
                )
                if abs(rel_delta) > cell_tol or abs(sumw2_rel) > cell_tol or fill_delta != 0.0:
                    failures.append(
                        _failure(
                            category,
                            f"weighted cell mismatch {lane_id}:{observable}:bin{index}",
                            relative_delta=rel_delta,
                            sumw2_relative_delta=sumw2_rel,
                            fill_delta=fill_delta,
                        )
                    )

            aggregate_key = f"{identity['family']}:{observable}"
            for label, edges, values in (
                ("reference", ref_edges, ref_sumw), ("candidate", cand_edges, cand_sumw)
            ):
                prior_edges = aggregate_edges[label].get(aggregate_key)
                if prior_edges is not None and prior_edges != edges:
                    failures.append(_failure("binning_mismatch", f"aggregate binning mismatch: {aggregate_key}"))
                else:
                    aggregate_edges[label][aggregate_key] = edges
                    prior = aggregate_sums[label].get(aggregate_key)
                    aggregate_sums[label][aggregate_key] = (
                        list(values) if prior is None else [left + right for left, right in zip(prior, values)]
                    )

    for label in ("reference", "candidate"):
        for key, values in aggregate_sums[label].items():
            aggregate[label][key] = (aggregate_edges[label][key], values)
    return rows, aggregate


def _estimator_inputs(
    aggregate: dict[str, dict[str, tuple[list[float], list[float]]]],
    contract: dict[str, Any],
    failures: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    required = [f"inclusive:{name}" for name in ("A", "B", "C", "D")] + [
        f"photon:{name}_signal" for name in ("A", "B", "C", "D")
    ]
    if any(key not in aggregate[side] for side in ("reference", "candidate") for key in required):
        return []
    edges = aggregate["reference"]["inclusive:A"][0]
    if any(aggregate[side][key][0] != edges for side in ("reference", "candidate") for key in required):
        failures.append(_failure("binning_mismatch", "stitched estimator inputs do not share one binning"))
        return []
    rows: list[dict[str, Any]] = []
    tol = float(contract["tolerances"]["weighted_cell_relative"])
    abs_floor = float(contract["tolerances"]["weighted_cell_abs_floor"])
    for index in range(len(edges) - 1):
        row: dict[str, Any] = {"bin": index + 1, "x_low": edges[index], "x_high": edges[index + 1]}
        for side in ("reference", "candidate"):
            for region in ("A", "B", "C", "D"):
                row[f"{side}_{region}"] = aggregate[side][f"inclusive:{region}"][1][index]
            signal_a = aggregate[side]["photon:A_signal"][1][index]
            for region in ("B", "C", "D"):
                numerator = aggregate[side][f"photon:{region}_signal"][1][index]
                row[f"{side}_c{region}"] = numerator / signal_a if signal_a != 0.0 else math.nan
        if any(not _finite_number(row[f"{side}_{region}"]) or row[f"{side}_{region}"] <= 0.0
               for side in ("reference", "candidate") for region in ("A", "B", "C", "D")):
            failures.append(_failure("canary_underpopulated", f"ABCD estimator bin {index + 1} is not populated"))
        if any(
            not _finite_number(row[f"{side}_c{region}"])
            or row[f"{side}_c{region}"] <= 0.0
            for side in ("reference", "candidate")
            for region in ("B", "C", "D")
        ):
            failures.append(_failure("canary_underpopulated", f"leakage estimator bin {index + 1} is not populated"))
        for key in ("A", "B", "C", "D", "cB", "cC", "cD"):
            absolute, relative, ratio = _cell_delta(
                float(row[f"candidate_{key}"]), float(row[f"reference_{key}"]), abs_floor
            )
            row[f"{key}_candidate_over_reference"] = ratio
            row[f"{key}_relative_delta"] = relative
            category = "photon_leakage_failed" if key.startswith("c") else "inclusive_abcd_failed"
            if abs(relative) > tol:
                failures.append(
                    _failure(category, f"stitched estimator input mismatch {key}:bin{index + 1}", relative_delta=relative)
                )
        rows.append(row)
    return rows


def _purity_gate(
    reference: dict[str, Any],
    candidate: dict[str, Any],
    contract: dict[str, Any],
    failures: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    ref_purity = reference.get("purity")
    cand_purity = candidate.get("purity")
    if not isinstance(ref_purity, dict) or not isinstance(cand_purity, dict):
        failures.append(_failure("purity_solver_failed", "truth/raw/corrected purity payloads are required"))
        return []

    for label, purity in (("reference", ref_purity), ("candidate", cand_purity)):
        repetition = purity.get("fixed_seed_repetition")
        series_payload = {
            "bin_edges": purity.get("bin_edges"),
            **{name: purity.get(name) for name in ("truth", "raw", "corrected")},
        }
        try:
            observed_payload_sha = _payload_sha256(series_payload)
        except (TypeError, ValueError):
            observed_payload_sha = None
        first_sha = repetition.get("first_output_sha256") if isinstance(repetition, dict) else None
        repeated_sha = repetition.get("repeated_output_sha256") if isinstance(repetition, dict) else None
        if (
            not _is_sha256(first_sha)
            or not _is_sha256(repeated_sha)
            or first_sha != repeated_sha
            or observed_payload_sha is None
            or first_sha != observed_payload_sha
        ):
            failures.append(
                _failure(
                    "purity_solver_failed",
                    f"{label} fixed-seed estimator repetition is missing or non-deterministic",
                    first_output_sha256=first_sha,
                    repeated_output_sha256=repeated_sha,
                    manifest_payload_sha256=observed_payload_sha,
                )
            )
    try:
        ref_edges = _as_finite_vector(ref_purity.get("bin_edges"), "reference.purity.bin_edges")
        cand_edges = _as_finite_vector(cand_purity.get("bin_edges"), "candidate.purity.bin_edges")
    except GateInputError as exc:
        failures.append(_failure(exc.code, str(exc)))
        return []
    if ref_edges != cand_edges:
        failures.append(_failure("binning_mismatch", "purity bin edges differ"))
        return []
    rows: list[dict[str, Any]] = []
    max_tol = float(contract["tolerances"]["purity_ratio_max_abs"])
    rms_tol = float(contract["tolerances"]["purity_ratio_rms"])
    abs_floor = float(contract["tolerances"]["weighted_cell_abs_floor"])
    for series in ("truth", "raw", "corrected"):
        ref_series = ref_purity.get(series, {})
        cand_series = cand_purity.get(series, {})
        try:
            ref_values = _as_finite_vector(ref_series.get("value"), f"reference.purity.{series}.value")
            cand_values = _as_finite_vector(cand_series.get("value"), f"candidate.purity.{series}.value")
            ref_errors = _as_finite_vector(ref_series.get("error"), f"reference.purity.{series}.error")
            cand_errors = _as_finite_vector(cand_series.get("error"), f"candidate.purity.{series}.error")
        except GateInputError as exc:
            failures.append(_failure(exc.code, str(exc)))
            continue
        if any(len(values) != len(ref_edges) - 1 for values in (ref_values, cand_values, ref_errors, cand_errors)):
            failures.append(_failure("binning_mismatch", f"purity vector length differs: {series}"))
            continue
        relative_values: list[float] = []
        relative_errors: list[float] = []
        for index, (ref_value, cand_value, ref_error, cand_error) in enumerate(
            zip(ref_values, cand_values, ref_errors, cand_errors), start=1
        ):
            absolute, relative, ratio = _cell_delta(cand_value, ref_value, abs_floor)
            error_absolute, error_relative, error_ratio = _cell_delta(
                cand_error, ref_error, abs_floor
            )
            relative_values.append(relative)
            relative_errors.append(error_relative)
            rows.append(
                {
                    "series": series,
                    "bin": index,
                    "x_low": ref_edges[index - 1],
                    "x_high": ref_edges[index],
                    "reference_value": ref_value,
                    "candidate_value": cand_value,
                    "candidate_over_reference": ratio,
                    "absolute_delta": absolute,
                    "relative_delta": relative,
                    "reference_error": ref_error,
                    "candidate_error": cand_error,
                    "error_candidate_over_reference": error_ratio,
                    "error_absolute_delta": error_absolute,
                    "error_relative_delta": error_relative,
                }
            )
        maximum = max((abs(value) for value in relative_values), default=math.inf)
        rms = math.sqrt(sum(value * value for value in relative_values) / len(relative_values)) if relative_values else math.inf
        if maximum > max_tol or rms > rms_tol:
            failures.append(
                _failure(
                    "stitched_purity_failed",
                    f"{series} purity parity failed",
                    max_abs_ratio_minus_one=maximum,
                    rms_ratio_minus_one=rms,
                    max_allowed=max_tol,
                    rms_allowed=rms_tol,
                )
            )
        error_maximum = max(
            (abs(value) for value in relative_errors), default=math.inf
        )
        error_rms = (
            math.sqrt(
                sum(value * value for value in relative_errors)
                / len(relative_errors)
            )
            if relative_errors
            else math.inf
        )
        if error_maximum > max_tol or error_rms > rms_tol:
            failures.append(
                _failure(
                    "stitched_purity_failed",
                    f"{series} purity uncertainty parity failed",
                    max_abs_error_ratio_minus_one=error_maximum,
                    rms_error_ratio_minus_one=error_rms,
                    max_allowed=max_tol,
                    rms_allowed=rms_tol,
                )
            )
    return rows


def _merge_audit_gate(
    payload: dict[str, Any],
    merge_audit_path: Path,
    candidate_manifest_path: Path,
    contract: dict[str, Any],
    failures: list[dict[str, Any]],
) -> dict[str, Any]:
    canonical_outputs: list[dict[str, str]] = []
    try:
        _, canonical_outputs = _canonical_assembler(contract).verify_merge_audit_evidence(
            merge_audit_path, candidate_manifest_path, contract
        )
    except Exception as exc:
        failures.append(
            _failure(
                "merge_arithmetic_failed",
                f"merge audit is not canonical ROOT-native evidence: {exc}",
            )
        )
    if payload.get("schema") != "ppg12-stitched-purity-merge-audit/v1":
        failures.append(_failure("merge_arithmetic_failed", "unsupported merge-audit schema"))
    candidate_link = payload.get("candidate_manifest")
    expected_candidate_sha = _file_sha256(candidate_manifest_path)
    if not isinstance(candidate_link, dict):
        failures.append(
            _failure(
                "merge_arithmetic_failed",
                "merge audit is not bound to the exact candidate manifest",
            )
        )
    else:
        raw_path = candidate_link.get("path")
        linked_sha = candidate_link.get("sha256")
        linked_path = Path(raw_path) if isinstance(raw_path, str) else None
        if linked_path is not None and not linked_path.is_absolute():
            linked_path = (merge_audit_path.parent / linked_path).resolve()
        if (
            linked_path is None
            or linked_path != candidate_manifest_path.resolve()
            or linked_sha != expected_candidate_sha
        ):
            failures.append(
                _failure(
                    "merge_arithmetic_failed",
                    "merge audit candidate-manifest binding is stale or mismatched",
                    linked_path=str(linked_path) if linked_path is not None else raw_path,
                    linked_sha256=linked_sha,
                    expected_path=str(candidate_manifest_path.resolve()),
                    expected_sha256=expected_candidate_sha,
                )
            )
    audits = payload.get("audits") if isinstance(payload.get("audits"), list) else [payload]
    expected_counts = {"inclusive": 20, "photon": 12}
    observed_counts: dict[str, int] = {}
    all_inputs: list[str] = []
    merge_tol = float(contract["tolerances"]["merge_abs"])
    for index, audit in enumerate(audits):
        if not isinstance(audit, dict):
            failures.append(_failure("merge_arithmetic_failed", f"merge audit {index} is not an object"))
            continue
        family = audit.get("family")
        if len(audits) == 1 and family is None:
            family = "combined"
        status = str(audit.get("status", "")).upper()
        if status != "PASS" or audit.get("failures") not in ([], None):
            failures.append(_failure("merge_arithmetic_failed", f"merge audit failed: {family}", status=status))
        for field in ("max_content_delta", "max_sumw2_delta"):
            value = audit.get(field)
            if not _finite_number(value) or abs(float(value)) > merge_tol:
                failures.append(_failure("merge_arithmetic_failed", f"{family} {field} exceeds tolerance", observed=value))
        inputs = audit.get("inputs_fixed_order") or audit.get("input_ids")
        if not isinstance(inputs, list) or not inputs or not all(isinstance(item, str) and item for item in inputs):
            failures.append(_failure("merge_arithmetic_failed", f"{family} merge audit lacks exact ordered inputs"))
            inputs = []
        declared_count = audit.get("inputs_fixed_order_count", len(inputs))
        if declared_count != len(inputs):
            failures.append(_failure("merge_arithmetic_failed", f"{family} input count does not match exact list"))
        if len(set(inputs)) != len(inputs):
            failures.append(_failure("merge_arithmetic_failed", f"{family} merge input is duplicated"))
        all_inputs.extend(inputs)
        if family in expected_counts:
            observed_counts[family] = len(inputs)
            if len(inputs) != expected_counts[family]:
                failures.append(
                    _failure(
                        "merge_arithmetic_failed",
                        f"{family} merge must consume exactly {expected_counts[family]} lane inputs",
                        observed=len(inputs),
                    )
                )
    if len(audits) == 1 and audits[0].get("family") is None and len(all_inputs) != 32:
        failures.append(_failure("merge_arithmetic_failed", "combined merge must consume exactly 32 lane inputs"))
    if len(set(all_inputs)) != len(all_inputs):
        failures.append(_failure("merge_arithmetic_failed", "an input appears in more than one merge audit"))
    if len(audits) > 1 and observed_counts != expected_counts:
        failures.append(_failure("merge_arithmetic_failed", "merge audit must contain inclusive=20 and photon=12 inputs"))
    return {
        "audit_count": len(audits),
        "input_count": len(all_inputs),
        "input_sha256": _payload_sha256(all_inputs),
        "output_artifacts": canonical_outputs,
    }


def _finalize_report(
    mode: str,
    contract_path: Path,
    contract: dict[str, Any],
    reference_path: Path | None,
    candidate_path: Path | None,
    failures: list[dict[str, Any]],
    summary: dict[str, Any],
) -> dict[str, Any]:
    unique: list[dict[str, Any]] = []
    seen: set[bytes] = set()
    for failure in failures:
        key = _canonical_bytes(failure)
        if key not in seen:
            seen.add(key)
            unique.append(failure)
    return {
        "schema": "ppg12-stitched-purity-gate-report/v1",
        "mode": mode,
        "status": "PASS" if not unique else "FAIL",
        "contract": {"path": str(contract_path), "sha256": _payload_sha256(contract)},
        "reference_manifest": str(reference_path) if reference_path else None,
        "candidate_manifest": str(candidate_path) if candidate_path else None,
        "failure_count": len(unique),
        "failures": unique,
        "summary": summary,
    }


def _run_pair_gate(
    mode: str,
    reference_path: Path,
    candidate_path: Path,
    contract_path: Path,
    merge_audit_path: Path | None = None,
) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], dict[str, Any], dict[str, Any]]:
    contract = _read_json(contract_path, "contract")
    reference = _read_json(reference_path, "reference manifest")
    candidate = _read_json(candidate_path, "candidate manifest")
    failures: list[dict[str, Any]] = []
    _canonical_manifest_gate(
        reference_path, contract, {"reference"}, "reference manifest", failures
    )
    _canonical_manifest_gate(
        candidate_path,
        contract,
        {"candidate", "production"},
        "candidate manifest",
        failures,
    )
    _validate_manifest_header(reference, "reference", {"reference"}, contract, failures)
    _validate_manifest_header(candidate, "candidate", {"candidate", "production"}, contract, failures)
    _compare_headers(reference, candidate, contract, failures)
    reference_lanes = _index_lanes(reference, "reference", contract, failures)
    candidate_lanes = _index_lanes(candidate, "candidate", contract, failures)
    parity_summary = _feature_and_tag_gate(candidate, contract, failures)
    delta_rows, aggregate = _compare_lanes(reference_lanes, candidate_lanes, contract, failures)
    estimator_rows = _estimator_inputs(aggregate, contract, failures)
    purity_rows = _purity_gate(reference, candidate, contract, failures)
    merge_summary: dict[str, Any] = {}
    if merge_audit_path is not None:
        merge_summary = _merge_audit_gate(
            _read_json(merge_audit_path, "merge audit"),
            merge_audit_path,
            candidate_path,
            contract,
            failures,
        )
    summary = {
        "expected_lane_count": len(_expected_lanes(contract)),
        "reference_lane_count": len(reference_lanes),
        "candidate_lane_count": len(candidate_lanes),
        "lane_observable_cell_count": len(delta_rows),
        "estimator_bin_count": len(estimator_rows),
        "purity_point_count": len(purity_rows),
        "candidate_parity": parity_summary,
        "merge": merge_summary,
    }
    report = _finalize_report(mode, contract_path, contract, reference_path, candidate_path, failures, summary)
    return report, delta_rows, estimator_rows, purity_rows, reference, candidate


def _write_pair_outputs(
    outdir: Path,
    report: dict[str, Any],
    delta_rows: list[dict[str, Any]],
    estimator_rows: list[dict[str, Any]],
    purity_rows: list[dict[str, Any]],
) -> None:
    _write_json(outdir / "gate_report.json", report)
    _write_csv(
        outdir / "lane_observable_deltas.csv",
        delta_rows,
        ("lane_id", "family", "sample", "period", "interaction", "observable", "bin"),
    )
    _write_csv(outdir / "estimator_inputs.csv", estimator_rows, ("bin", "x_low", "x_high"))
    _write_csv(outdir / "purity_points.csv", purity_rows, ("series", "bin", "x_low", "x_high"))


def _frozen_snapshot(manifest: dict[str, Any], contract: dict[str, Any]) -> dict[str, Any]:
    # The canary and broad production intentionally process different event
    # subsets, so event_set_sha256 must agree only within each PPG12/RecoilJets
    # pair.  The full source-list and physics contracts remain frozen across
    # admission and production.
    lane_fields = contract["production_frozen_lane_fields"]
    provenance = manifest.get("provenance", {})
    lanes = manifest.get("lanes", [])
    return {
        "provenance": {field: provenance.get(field) for field in contract["frozen_provenance_fields"]},
        "lanes": {
            lane.get("lane_id", f"invalid:{index}"): {field: lane.get(field) for field in lane_fields}
            for index, lane in enumerate(sorted(
                (row for row in lanes if isinstance(row, dict)), key=lambda row: str(row.get("lane_id"))
            ))
        },
        "abcd_population": manifest.get("abcd_population"),
        "external_scale": manifest.get("external_scale"),
        "random_seed": manifest.get("random_seed"),
        "toy_count": manifest.get("toy_count"),
    }


def _admission_payload(
    report: dict[str, Any],
    reference_path: Path,
    candidate_path: Path,
    merge_audit_path: Path,
    reference: dict[str, Any],
    candidate: dict[str, Any],
    contract: dict[str, Any],
) -> dict[str, Any]:
    return {
        "schema": "ppg12-stitched-purity-admission/v1",
        "status": "PASS",
        "contract_sha256": _payload_sha256(contract),
        "reference_manifest": {"path": str(reference_path), "sha256": _file_sha256(reference_path)},
        "candidate_manifest": {"path": str(candidate_path), "sha256": _file_sha256(candidate_path)},
        "merge_audit": {"path": str(merge_audit_path), "sha256": _file_sha256(merge_audit_path)},
        "gate_report_payload_sha256": _payload_sha256(report),
        "reference_frozen": _frozen_snapshot(reference, contract),
        "candidate_frozen": _frozen_snapshot(candidate, contract),
    }


def _resolve_link(wrapper_path: Path, link: Any, label: str) -> Path:
    if not isinstance(link, dict) or not isinstance(link.get("path"), str) or not isinstance(link.get("sha256"), str):
        raise GateInputError("manifest_invalid", f"production {label} link must contain path and sha256")
    path = Path(link["path"])
    if not path.is_absolute():
        path = (wrapper_path.parent / path).resolve()
    if not path.exists() or _file_sha256(path) != link["sha256"]:
        raise GateInputError("hash_drift", f"production {label} hash does not match: {path}")
    return path


def _check_frozen(
    label: str,
    frozen: dict[str, Any],
    manifest: dict[str, Any],
    contract: dict[str, Any],
    failures: list[dict[str, Any]],
) -> None:
    current = _frozen_snapshot(manifest, contract)
    if frozen != current:
        failures.append(_failure("hash_drift", f"{label} contract drifted after admission"))


def _candidate_artifacts_gate(
    wrapper_path: Path,
    payload: Any,
    contract: dict[str, Any],
    failures: list[dict[str, Any]],
) -> list[dict[str, str]]:
    if not isinstance(payload, list):
        failures.append(_failure("promotion_blocked", "production candidate_artifacts must be a list"))
        return []
    expected = set(contract["required_candidate_artifact_families"])
    observed: dict[str, dict[str, str]] = {}
    for index, row in enumerate(payload):
        if not isinstance(row, dict):
            failures.append(_failure("manifest_invalid", f"candidate_artifacts[{index}] is not an object"))
            continue
        family = row.get("family")
        raw_path = row.get("path")
        expected_hash = row.get("sha256")
        if family not in expected or not isinstance(raw_path, str) or not isinstance(expected_hash, str):
            failures.append(_failure("manifest_invalid", f"invalid candidate artifact entry {index}"))
            continue
        if family in observed:
            failures.append(_failure("promotion_blocked", f"duplicate candidate artifact family: {family}"))
            continue
        path = Path(raw_path)
        if not path.is_absolute():
            path = (wrapper_path.parent / path).resolve()
        if not path.exists() or _file_sha256(path) != expected_hash:
            failures.append(
                _failure(
                    "hash_drift",
                    f"candidate artifact hash does not match: {family}",
                    path=str(path),
                )
            )
            continue
        observed[family] = {"family": family, "path": str(path), "sha256": expected_hash}
    for family in sorted(expected - set(observed)):
        failures.append(_failure("promotion_blocked", f"missing candidate artifact family: {family}"))
    return [observed[family] for family in sorted(observed)]


def _historical_gate(payload: Any, contract: dict[str, Any], failures: list[dict[str, Any]]) -> dict[str, Any]:
    if not isinstance(payload, dict):
        failures.append(_failure("promotion_blocked", "historical archive comparison is required"))
        return {}
    tolerances = contract["tolerances"]
    checks = (
        ("chi2_ndf", "historical_chi2_ndf_max", lambda value, limit: value <= limit),
        ("max_abs_pull", "historical_pull_max_abs", lambda value, limit: value <= limit),
        (
            "abcd_coherent_trend_max_sigma",
            "historical_coherent_trend_max_sigma",
            lambda value, limit: value <= limit,
        ),
        (
            "leakage_coherent_trend_max_sigma",
            "historical_coherent_trend_max_sigma",
            lambda value, limit: value <= limit,
        ),
    )
    for field, limit_field, predicate in checks:
        value = payload.get(field)
        limit = float(tolerances[limit_field])
        if not _finite_number(value) or float(value) < 0.0 or not predicate(float(value), limit):
            failures.append(_failure("promotion_blocked", f"historical comparison failed: {field}", observed=value, limit=limit))
    mean = payload.get("weighted_mean_ratio")
    error = payload.get("weighted_mean_ratio_error")
    sigma = float(contract["production_historical_requirements"]["weighted_mean_ratio_compatible_with_unity_sigma"])
    if not _finite_number(mean) or not _finite_number(error) or float(error) <= 0.0 or abs(float(mean) - 1.0) > sigma * float(error):
        failures.append(
            _failure(
                "promotion_blocked",
                "historical weighted mean ratio is not compatible with unity",
                weighted_mean_ratio=mean,
                weighted_mean_ratio_error=error,
            )
        )
    return payload


def _historical_evidence_gate(
    wrapper_path: Path,
    wrapper: dict[str, Any],
    embedded: dict[str, Any],
    candidate_artifacts: list[dict[str, str]],
    candidate_manifest: dict[str, Any],
    contract: dict[str, Any],
    failures: list[dict[str, Any]],
) -> None:
    """Bind historical metrics to their exact evidence and candidate artifacts."""
    assembly = wrapper.get("assembly")
    if not isinstance(assembly, dict):
        failures.append(_failure("promotion_blocked", "production wrapper assembly evidence is required"))
        return
    if assembly.get("contract_sha256") != _payload_sha256(contract):
        failures.append(_failure("hash_drift", "production wrapper contract hash is stale"))
    expected_artifact_sha = _payload_sha256(candidate_artifacts)
    if assembly.get("candidate_artifact_set_sha256") != expected_artifact_sha:
        failures.append(_failure("hash_drift", "production wrapper candidate artifact set is stale"))
    try:
        evidence_path = _resolve_link(
            wrapper_path,
            assembly.get("historical_comparison"),
            "historical_comparison",
        )
        evidence = _read_json(evidence_path, "historical comparison evidence")
        if evidence.get("schema") != "ppg12-stitched-purity-historical-comparison/v1":
            failures.append(_failure("promotion_blocked", "unsupported historical-comparison evidence schema"))
            return
        if evidence.get("final_purity_series") != "corrected":
            failures.append(_failure("promotion_blocked", "historical comparison is not the corrected-purity series"))
        input_evidence = evidence.get("input")
        if not isinstance(input_evidence, dict) or input_evidence.get("mode") != "direct_root_objects":
            failures.append(_failure("promotion_blocked", "historical metrics were not derived directly from ROOT objects"))
        evidence_metrics = {
            key: evidence.get(key)
            for key in (
                "chi2_ndf",
                "max_abs_pull",
                "weighted_mean_ratio",
                "weighted_mean_ratio_error",
                "abcd_coherent_trend_max_sigma",
                "leakage_coherent_trend_max_sigma",
            )
        }
        if evidence_metrics != embedded:
            failures.append(
                _failure(
                    "hash_drift",
                    "embedded historical metrics differ from their hash-bound evidence",
                )
            )
        source_links = evidence.get("source_links")
        if not isinstance(source_links, list):
            failures.append(_failure("promotion_blocked", "historical source_links are required"))
            return
        normalized_links: list[dict[str, str]] = []
        observed_roles: set[str] = set()
        for index, link in enumerate(source_links):
            if not isinstance(link, dict) or link.get("role") not in REQUIRED_HISTORICAL_SOURCE_ROLES:
                failures.append(_failure("promotion_blocked", f"invalid historical source link {index}"))
                continue
            role = str(link["role"])
            if role in observed_roles:
                failures.append(_failure("promotion_blocked", f"duplicate historical source role: {role}"))
                continue
            observed_roles.add(role)
            try:
                source_path = _resolve_link(evidence_path, link, f"historical source {role}")
            except GateInputError as exc:
                failures.append(_failure(exc.code, str(exc)))
                continue
            normalized_links.append(
                {"role": role, "path": str(source_path), "sha256": str(link["sha256"])}
            )
        if observed_roles != REQUIRED_HISTORICAL_SOURCE_ROLES:
            failures.append(
                _failure(
                    "promotion_blocked",
                    "historical evidence does not bind the complete source set",
                    missing=sorted(REQUIRED_HISTORICAL_SOURCE_ROLES - observed_roles),
                )
            )
        normalized_links.sort(key=lambda row: row["role"])
        if evidence.get("source_set_sha256") != _payload_sha256(normalized_links):
            failures.append(_failure("hash_drift", "historical source-set hash is stale"))

        indexed_sources = {row["role"]: row for row in normalized_links}
        indexed_artifacts = {row["family"]: row for row in candidate_artifacts}
        for role, family in (
            ("candidate_inclusive", "inclusive"),
            ("candidate_photon", "photon"),
        ):
            source = indexed_sources.get(role)
            artifact = indexed_artifacts.get(family)
            if source is not None and artifact is not None and (
                source["sha256"] != artifact["sha256"]
            ):
                failures.append(
                    _failure(
                        "hash_drift",
                        f"historical evidence {role} differs from promoted candidate artifact",
                    )
                )
        candidate_purity_link = (
            candidate_manifest.get("assembly", {}).get("purity")
            if isinstance(candidate_manifest.get("assembly"), dict)
            else None
        )
        candidate_purity_source = indexed_sources.get("candidate_purity")
        if not isinstance(candidate_purity_link, dict):
            failures.append(
                _failure(
                    "manifest_invalid",
                    "production candidate manifest lacks hash-bound purity evidence",
                )
            )
        elif (
            candidate_purity_source is not None
            and candidate_purity_source["sha256"]
            != candidate_purity_link.get("sha256")
        ):
            failures.append(
                _failure(
                    "hash_drift",
                    "historical evidence candidate purity differs from production manifest purity",
                )
            )
    except GateInputError as exc:
        failures.append(_failure(exc.code, str(exc)))


def _exit_code(report: dict[str, Any], contract: dict[str, Any]) -> int:
    if report["status"] == "PASS":
        return 0
    codes = contract["terminal_codes"]
    observed = [codes.get(failure["code"], codes["manifest_invalid"]) for failure in report["failures"]]
    return min(observed) if observed else codes["manifest_invalid"]


def _run_diagnose_or_admit(args: argparse.Namespace) -> int:
    contract_path = Path(args.contract).resolve()
    reference_path = Path(args.reference_manifest).resolve()
    candidate_path = Path(args.candidate_manifest).resolve()
    outdir = Path(args.outdir).resolve()
    contract = _read_json(contract_path, "contract")
    merge = Path(args.merge_audit).resolve() if getattr(args, "merge_audit", None) else None
    try:
        report, deltas, estimator, purity, reference, candidate = _run_pair_gate(
            args.command, reference_path, candidate_path, contract_path, merge
        )
    except GateInputError as exc:
        report = _finalize_report(
            args.command,
            contract_path,
            contract,
            reference_path,
            candidate_path,
            [_failure(exc.code, str(exc))],
            {},
        )
        deltas, estimator, purity = [], [], []
        reference, candidate = {}, {}
    _write_pair_outputs(outdir, report, deltas, estimator, purity)
    admission = outdir / "admission_manifest.json"
    if admission.exists():
        admission.unlink()
    if args.command == "admit" and report["status"] == "PASS":
        assert merge is not None
        payload = _admission_payload(report, reference_path, candidate_path, merge, reference, candidate, contract)
        _write_json(admission, payload)
    print(json.dumps({"status": report["status"], "report": str(outdir / "gate_report.json"), "failures": report["failure_count"]}))
    return _exit_code(report, contract)


def _run_verify_production(args: argparse.Namespace) -> int:
    contract_path = Path(args.contract).resolve()
    contract = _read_json(contract_path, "contract")
    admission_path = Path(args.admission_manifest).resolve()
    production_path = Path(args.production_manifest).resolve()
    merge_path = Path(args.merge_audit).resolve()
    outdir = Path(args.outdir).resolve()
    failures: list[dict[str, Any]] = []
    report: dict[str, Any]
    deltas: list[dict[str, Any]] = []
    estimator: list[dict[str, Any]] = []
    purity: list[dict[str, Any]] = []
    production_report = outdir / "production_gate_report.json"
    if production_report.exists():
        production_report.unlink()
    try:
        admission = _read_json(admission_path, "admission manifest")
        wrapper = _read_json(production_path, "production manifest")
        if admission.get("schema") != "ppg12-stitched-purity-admission/v1" or admission.get("status") != "PASS":
            failures.append(_failure("promotion_blocked", "admission manifest is not a passing v1 admission"))
        if admission.get("contract_sha256") != _payload_sha256(contract):
            failures.append(_failure("hash_drift", "contract changed after admission"))
        admission_reference_path = _resolve_link(
            admission_path,
            admission.get("reference_manifest"),
            "admission reference_manifest",
        )
        admission_candidate_path = _resolve_link(
            admission_path,
            admission.get("candidate_manifest"),
            "admission candidate_manifest",
        )
        admission_merge_path = _resolve_link(
            admission_path,
            admission.get("merge_audit"),
            "admission merge_audit",
        )

        # An admission JSON is a receipt, not an authority.  Re-run the exact
        # admission pair gate from its linked immutable evidence and rebuild
        # the canonical receipt.  This rejects a hand-authored PASS-shaped
        # document even when all of its individual path/SHA links exist.
        (
            admission_report,
            _,
            _,
            _,
            admission_reference,
            admission_candidate,
        ) = _run_pair_gate(
            "admit",
            admission_reference_path,
            admission_candidate_path,
            contract_path,
            admission_merge_path,
        )
        failures.extend(admission_report["failures"])
        if admission_report["status"] == "PASS":
            replayed_admission = _admission_payload(
                admission_report,
                admission_reference_path,
                admission_candidate_path,
                admission_merge_path,
                admission_reference,
                admission_candidate,
                contract,
            )
            if admission != replayed_admission:
                failures.append(
                    _failure(
                        "hash_drift",
                        "admission manifest is not the canonical replay-derived receipt",
                    )
                )
        else:
            failures.append(
                _failure(
                    "promotion_blocked",
                    "linked admission evidence no longer passes the admission gate",
                )
            )
        if wrapper.get("schema") != "ppg12-stitched-purity-production/v1":
            failures.append(_failure("manifest_invalid", "unsupported production wrapper schema"))
        if wrapper.get("admission_sha256") != _file_sha256(admission_path):
            failures.append(_failure("hash_drift", "production wrapper does not reference this exact admission"))
        reference_path = _resolve_link(production_path, wrapper.get("reference_manifest"), "reference_manifest")
        candidate_path = _resolve_link(production_path, wrapper.get("candidate_manifest"), "candidate_manifest")
        wrapper_merge_path = _resolve_link(
            production_path, wrapper.get("merge_audit"), "merge_audit"
        )
        if wrapper_merge_path != merge_path:
            failures.append(
                _failure("hash_drift", "production wrapper references a different merge audit")
            )
        pair_report, deltas, estimator, purity, reference, candidate = _run_pair_gate(
            "verify-production", reference_path, candidate_path, contract_path, merge_path
        )
        failures.extend(pair_report["failures"])
        _check_frozen("reference", admission.get("reference_frozen", {}), reference, contract, failures)
        _check_frozen("candidate", admission.get("candidate_frozen", {}), candidate, contract, failures)
        candidate_artifacts = _candidate_artifacts_gate(
            production_path, wrapper.get("candidate_artifacts"), contract, failures
        )
        merge_artifacts = pair_report.get("summary", {}).get("merge", {}).get(
            "output_artifacts", []
        )
        if candidate_artifacts != merge_artifacts:
            failures.append(
                _failure(
                    "merge_arithmetic_failed",
                    "promoted candidate artifacts are not the exact ROOT-native merge outputs",
                )
            )
        historical = _historical_gate(wrapper.get("historical_archive_comparison"), contract, failures)
        _historical_evidence_gate(
            production_path,
            wrapper,
            historical,
            candidate_artifacts,
            candidate,
            contract,
            failures,
        )
        report = _finalize_report(
            "verify-production",
            contract_path,
            contract,
            reference_path,
            candidate_path,
            failures,
            {
                **pair_report["summary"],
                "admission_manifest": str(admission_path),
                "admission_sha256": _file_sha256(admission_path),
                "production_wrapper": str(production_path),
                "production_wrapper_sha256": _file_sha256(production_path),
                "historical_archive_comparison": historical,
                "candidate_artifacts": candidate_artifacts,
                "candidate_artifact_set_sha256": _payload_sha256(candidate_artifacts),
            },
        )
    except GateInputError as exc:
        report = _finalize_report(
            "verify-production",
            contract_path,
            contract,
            None,
            None,
            failures + [_failure(exc.code, str(exc))],
            {"admission_manifest": str(admission_path), "production_wrapper": str(production_path)},
        )
    _write_pair_outputs(outdir, report, deltas, estimator, purity)
    if report["status"] == "PASS":
        _write_json(
            production_report,
            {
                "schema": "ppg12-stitched-purity-production-gate/v1",
                "status": "PASS",
                "contract_sha256": _payload_sha256(contract),
                "gate_report_payload_sha256": _payload_sha256(report),
                "admission_sha256": _file_sha256(admission_path),
                "production_manifest_sha256": _file_sha256(production_path),
                "reference_manifest_sha256": _file_sha256(reference_path),
                "candidate_manifest_sha256": _file_sha256(candidate_path),
                "merge_audit_sha256": _file_sha256(merge_path),
                "candidate_artifacts": report["summary"]["candidate_artifacts"],
                "candidate_artifact_set_sha256": report["summary"]["candidate_artifact_set_sha256"],
            },
        )
    print(json.dumps({"status": report["status"], "report": str(outdir / "gate_report.json"), "failures": report["failure_count"]}))
    return _exit_code(report, contract)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", default=str(DEFAULT_CONTRACT), help="JSON-compatible YAML closure contract")
    subparsers = parser.add_subparsers(dest="command", required=True)
    for name in ("diagnose", "admit"):
        command = subparsers.add_parser(name)
        command.add_argument("--reference-manifest", required=True)
        command.add_argument("--candidate-manifest", required=True)
        command.add_argument("--outdir", required=True)
        if name == "admit":
            command.add_argument("--merge-audit", required=True)
    verify = subparsers.add_parser("verify-production")
    verify.add_argument("--admission-manifest", required=True)
    verify.add_argument("--production-manifest", required=True)
    verify.add_argument("--merge-audit", required=True)
    verify.add_argument("--outdir", required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _parser()
    args = parser.parse_args(argv)
    if args.command in {"diagnose", "admit"}:
        return _run_diagnose_or_admit(args)
    return _run_verify_production(args)


if __name__ == "__main__":
    sys.exit(main())
