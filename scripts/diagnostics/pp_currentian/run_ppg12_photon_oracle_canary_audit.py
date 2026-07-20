#!/usr/bin/env python3
"""Audit the exact 12-lane source-locked paired executable photon canary.

This is a terminal-only consumer.  It never queries or controls Condor and it
never submits, merges, promotes, or patches physics source.  It fails closed
unless the plan and PASS receipts are exact and each physical lane contains a
unique, contract-bound preserved-PPG12 trace, paired candidate CSV, and
executable aggregate. Archived PPG12 ROOTs and Python selection shadows are
not accepted as oracle evidence.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

try:
    import ROOT  # type: ignore
except ImportError:
    ROOT = None  # type: ignore


PLAN_SCHEMA = "ppg12-stitched-purity-photon-canary-plan/v5"
CANDIDATE_PARITY_SCHEMA = "ppg12-stitched-purity-candidate-parity/v1"
MIN_ROOT_BYTES = 50_000
BASE_E_MODEL = Path(
    "/sphenix/user/shuhangli/ppg12/FunWithxgboost/"
    "binned_models/model_base_E_split_single_tmva.root"
)
BASE_V3E_MODEL = Path(
    "/sphenix/user/shuhangli/ppg12/FunWithxgboost/"
    "binned_models/model_base_v3E_split_single_tmva.root"
)
MASK_ROOT = Path(
    "/sphenix/user/shuhangli/ppg12/efficiencytool/tower_masks_bdt_nom.root"
)

FEATURE_NAMES = (
    "cluster_Et_score_input",
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
)
FEATURE_NAME_ALIASES = {"cluster_Et_score_input": "cluster_Et"}
SCORE_NAME_ALIASES = {"base_v3E": "baseV3E"}
CANDIDATE_PARITY_TOLERANCES = {
    "feature_abs_floor": 1.0e-6,
    "feature_relative": 1.0e-5,
    "score_abs": 1.0e-6,
}
SCORE_PAIRS = (
    ("base_E", "rj_bdt_base_E", "ppg12_bdt_base_E"),
    ("base_v3E", "rj_bdt_base_v3E", "ppg12_bdt_base_v3E"),
    ("selected_recomputed", "rj_selected_bdt_score", "ppg12_selected_bdt_score"),
    ("selected_stored", "rj_stored_bdt_score", "ppg12_selected_bdt_score"),
)
ISOLATION_PAIRS = (
    ("raw_eiso", "rj_raw_eiso", "ppg12_raw_eiso"),
    ("corrected_eiso", "rj_corrected_eiso", "ppg12_corrected_eiso"),
    ("iso_threshold", "rj_iso_threshold", "ppg12_iso_threshold"),
    ("noniso_threshold", "rj_noniso_threshold", "ppg12_noniso_threshold"),
)
FIRST_DIVERGENCE_ORDER = (
    "candidate_population",
    "features",
    "scores",
    "model_route",
    "tags",
    "isolation_abcd",
    "weight_metadata",
)

REQUIRED_RJ_BRANCHES = {
    "evt",
    "eventnumber",
    "is_signal",
    "truth_track_id",
    "cluster_Et",
    "cluster_Et_score_input",
    "cluster_Eta",
    "cluster_Phi",
    "ppg12_kin_vertexz",
    "cluster_weta_cogx",
    "cluster_wphi_cogx",
    "e11_over_e33",
    "e32_over_e35",
    "cluster_et1",
    "cluster_et2",
    "cluster_et3",
    "cluster_et4",
    "cluster_prob",
    "npb_score",
    "tight_bdt_score",
    "ppg12_common_pass",
    "ppg12_tight_tag",
    "event_weight",
    "ppg12_raw_eiso",
    "ppg12_reco_eiso",
    "ppg12_iso_threshold",
    "ppg12_noniso_threshold",
    "ppg12_is_iso",
    "ppg12_is_noniso",
    "cluster_index",
    "ppg12_sample_bin",
    "ppg12_xsec_pb",
    "ppg12_xsec_weight",
    "ppg12_window_low",
    "ppg12_window_high",
    "ppg12_truth_window_pass_r04",
}


def expected_lane_ids() -> list[str]:
    return [
        f"photon:photon{photon}:{period}:{interaction}"
        for photon in (5, 10, 20)
        for period in ("0mrad", "1p5mrad")
        for interaction in ("si", "di")
    ]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while chunk := stream.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_payload_sha256(payload: Any) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode()
    return hashlib.sha256(encoded).hexdigest()


def _bound_file(item: Any, path: Path, label: str) -> None:
    if not isinstance(item, dict):
        raise RuntimeError(f"{label} provenance is missing")
    if Path(str(item.get("path", ""))).resolve() != path.resolve():
        raise RuntimeError(f"{label} provenance path differs")
    if item.get("sha256") != sha256_file(path):
        raise RuntimeError(f"{label} provenance hash differs")


def validate_lane_evidence_bundle(lane_id: str, bundle: dict[str, Path]) -> dict[str, Any]:
    required = {
        "runtime_contract", "candidate_csv", "executable_aggregate",
        "trace_csv", "response_trace_csv",
    }
    if set(bundle) != required:
        raise RuntimeError(f"{lane_id}: paired evidence bundle role set differs")
    paths = {name: Path(path).resolve() for name, path in bundle.items()}
    for name, path in paths.items():
        if not path.is_file() or path.stat().st_size <= 0:
            raise RuntimeError(f"{lane_id}: missing paired evidence {name}: {path}")
    contract = json.loads(paths["runtime_contract"].read_text())
    lane = contract.get("lane", {})
    if lane.get("lane_id") != lane_id:
        raise RuntimeError(f"{lane_id}: contract embeds another physical lane")
    contract_hash = sha256_file(paths["runtime_contract"])
    contract_paths = contract.get("paths", {})
    for bundle_name, contract_name in (
        ("candidate_csv", "candidate_csv"),
        ("executable_aggregate", "ppg12_executable_aggregate"),
        ("trace_csv", "ppg12_executable_trace"),
        ("response_trace_csv", "ppg12_executable_response_trace"),
    ):
        if Path(str(contract_paths.get(contract_name, ""))).resolve() != paths[bundle_name]:
            raise RuntimeError(f"{lane_id}: {bundle_name} differs from runtime contract")
    rows = read_csv_rows(paths["candidate_csv"])
    if not rows:
        raise RuntimeError(f"{lane_id}: candidate CSV contains no rows")
    for row_number, row in enumerate(rows, start=2):
        if row.get("lane_id") != lane_id:
            raise RuntimeError(f"{lane_id}: candidate row {row_number} belongs to another lane")
        if row.get("runtime_contract_sha256") != contract_hash:
            raise RuntimeError(f"{lane_id}: candidate row {row_number} has stale contract binding")
    aggregate = json.loads(paths["executable_aggregate"].read_text())
    if aggregate.get("evidence_source") != "preserved_ppg12_executable_aggregate":
        raise RuntimeError(f"{lane_id}: aggregate is not preserved-executable evidence")
    if aggregate.get("status") != "PASS" or aggregate.get("mode") != "full":
        raise RuntimeError(f"{lane_id}: executable aggregate did not pass")
    identity = aggregate.get("lane_identity", {})
    if identity != {"lane_id": lane_id, "runtime_contract_sha256": contract_hash}:
        raise RuntimeError(f"{lane_id}: aggregate lane/contract identity differs")
    provenance = aggregate.get("provenance", {})
    for aggregate_role, bundle_name in (
        ("runtime_contract", "runtime_contract"),
        ("candidate_csv", "candidate_csv"),
        ("trace_csv", "trace_csv"),
        ("response_trace_csv", "response_trace_csv"),
    ):
        _bound_file(provenance.get(aggregate_role), paths[bundle_name], f"{lane_id} {aggregate_role}")
    return {
        "paths": paths,
        "contract_sha256": contract_hash,
        "candidate_rows": rows,
        "aggregate": aggregate,
    }


def build_trace_evidence(
    trace_bundle_by_lane: dict[str, dict[str, Path]],
) -> dict[str, Any]:
    expected = sorted(expected_lane_ids())
    if set(trace_bundle_by_lane) != set(expected):
        raise RuntimeError("paired trace bundles are not the exact 12 photon lanes")
    auditor_path = Path(__file__).with_name("audit_ppg12_recoiljets_paired_oracle.py").resolve()
    spec = importlib.util.spec_from_file_location("ppg12_trace_auditor", auditor_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load candidate trace auditor: {auditor_path}")
    auditor = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(auditor)
    lanes: list[dict[str, Any]] = []
    uniqueness: dict[str, set[str]] = {
        role: set() for role in (
            "runtime_contract", "candidate_csv", "executable_aggregate",
            "trace_csv", "response_trace_csv",
        )
    }
    digest_uniqueness = {role: set() for role in uniqueness}
    for lane_id in expected:
        validated = validate_lane_evidence_bundle(lane_id, trace_bundle_by_lane[lane_id])
        paths = validated["paths"]
        path = paths["candidate_csv"]
        rows = validated["candidate_rows"]
        report = auditor.analyze_candidate_report(
            path,
            expected_lane_id=lane_id,
            expected_runtime_contract_sha256=validated["contract_sha256"],
        )
        if report.get("status") != "PASS" or report.get("first_divergence") is not None:
            raise RuntimeError(f"candidate trace is not an executable parity PASS: {lane_id}")
        evidence = {}
        for role, item_path in paths.items():
            resolved = str(item_path)
            item_hash = sha256_file(item_path)
            if resolved in uniqueness[role] or item_hash in digest_uniqueness[role]:
                raise RuntimeError(f"{lane_id}: {role} path/hash is reused by another physical lane")
            uniqueness[role].add(resolved)
            digest_uniqueness[role].add(item_hash)
            evidence[role] = {"path": resolved, "sha256": item_hash}
        lanes.append({
            "lane_id": lane_id,
            "paired_evidence": evidence,
            "row_count": len(rows),
            "audit_payload_sha256": canonical_payload_sha256(report),
        })
    lanes.sort(key=lambda row: row["lane_id"])
    return {
        "auditor": {"path": str(auditor_path), "sha256": sha256_file(auditor_path)},
        "summary_builder": {
            "path": str(Path(__file__).resolve()),
            "sha256": sha256_file(Path(__file__).resolve()),
        },
        "lanes": lanes,
        "lane_set_sha256": canonical_payload_sha256(lanes),
    }


def finite_float(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def close_feature(left: Any, right: Any) -> bool:
    lval = finite_float(left)
    rval = finite_float(right)
    if lval is None or rval is None:
        return False
    return abs(lval - rval) <= max(1.0e-6, 1.0e-5 * abs(rval))


def close_score(left: Any, right: Any) -> bool:
    lval = finite_float(left)
    rval = finite_float(right)
    return lval is not None and rval is not None and abs(lval - rval) <= 1.0e-6


def _comparison_delta(left: Any, right: Any) -> tuple[float, float, bool]:
    """Return absolute/normalized feature deltas and the tolerance decision."""
    lval = finite_float(left)
    rval = finite_float(right)
    if lval is None or rval is None:
        # Keep the artifact valid JSON while forcing the closure gate to fail.
        return 0.0, 2.0, False
    absolute = abs(lval - rval)
    denominator = max(
        CANDIDATE_PARITY_TOLERANCES["feature_abs_floor"],
        CANDIDATE_PARITY_TOLERANCES["feature_relative"] * abs(rval),
    )
    normalized = absolute / denominator
    return absolute, normalized, normalized <= 1.0


def _score_delta(left: Any, right: Any) -> tuple[float, bool]:
    lval = finite_float(left)
    rval = finite_float(right)
    if lval is None or rval is None:
        return 2.0 * CANDIDATE_PARITY_TOLERANCES["score_abs"], False
    absolute = abs(lval - rval)
    return absolute, absolute <= CANDIDATE_PARITY_TOLERANCES["score_abs"]


def _tag_value(row: dict[str, str], name: str) -> int | None:
    try:
        return int(row.get(name, ""))
    except (TypeError, ValueError):
        return None


def _candidate_parity_lane(lane_id: str, rows: list[dict[str, str]]) -> dict[str, Any]:
    statuses = Counter(row.get("match_status", "") for row in rows)
    matched = [row for row in rows if row.get("match_status") == "matched"]

    feature_metrics = {
        FEATURE_NAME_ALIASES.get(source_name, source_name): {
            "comparisons": 0,
            "violations": 0,
            "max_abs_delta": 0.0,
            "max_normalized_delta": 0.0,
        }
        for source_name in FEATURE_NAMES
    }
    score_pairs = {
        SCORE_NAME_ALIASES.get(label, label): (left_name, right_name)
        for label, left_name, right_name in SCORE_PAIRS
        if label in {"base_E", "base_v3E"}
    }
    score_metrics = {
        name: {"comparisons": 0, "mismatches": 0, "max_abs_delta": 0.0}
        for name in score_pairs
    }
    route = {"comparisons": 0, "mismatches": 0}
    tags = {
        name: {"comparisons": 0, "mismatches": 0}
        for name in ("common", "tight", "non_tight", "abcd")
    }

    for row in matched:
        for source_name in FEATURE_NAMES:
            output_name = FEATURE_NAME_ALIASES.get(source_name, source_name)
            metric = feature_metrics[output_name]
            absolute, normalized, passes = _comparison_delta(
                row.get(f"{source_name}_rj"), row.get(f"{source_name}_ppg12")
            )
            metric["comparisons"] += 1
            metric["violations"] += int(not passes)
            metric["max_abs_delta"] = max(metric["max_abs_delta"], absolute)
            metric["max_normalized_delta"] = max(
                metric["max_normalized_delta"], normalized
            )

        for name, (left_name, right_name) in score_pairs.items():
            metric = score_metrics[name]
            absolute, passes = _score_delta(row.get(left_name), row.get(right_name))
            metric["comparisons"] += 1
            metric["mismatches"] += int(not passes)
            metric["max_abs_delta"] = max(metric["max_abs_delta"], absolute)

        route["comparisons"] += 1
        route["mismatches"] += int(
            row.get("model_route_agree") != "1"
            or row.get("rj_inferred_stored_model") != row.get("ppg12_selected_model")
        )

        reference_tag = _tag_value(row, "ppg12_recomputed_tag")
        recomputed_tag = _tag_value(row, "rj_recomputed_tag")
        stored_tag = _tag_value(row, "rj_stored_tag")
        tag_mismatches = {
            "common": (
                row.get("common_agree") != "1"
                or row.get("stored_common_agree") != "1"
            ),
            "tight": (
                reference_tag is None
                or recomputed_tag is None
                or stored_tag is None
                or (recomputed_tag == 1) != (reference_tag == 1)
                or (stored_tag == 1) != (reference_tag == 1)
            ),
            "non_tight": (
                reference_tag is None
                or recomputed_tag is None
                or stored_tag is None
                or (recomputed_tag == 2) != (reference_tag == 2)
                or (stored_tag == 2) != (reference_tag == 2)
            ),
            "abcd": (
                row.get("abcd_agree") != "1"
                or row.get("stored_abcd_agree") != "1"
            ),
        }
        for name, mismatch in tag_mismatches.items():
            tags[name]["comparisons"] += 1
            tags[name]["mismatches"] += int(mismatch)

    return {
        "lane_id": lane_id,
        "population": {
            "reference_count": statuses.get("matched", 0) + statuses.get("ppg12_only", 0),
            "candidate_count": statuses.get("matched", 0) + statuses.get("rj_only", 0),
            "unmatched_reference": statuses.get("ppg12_only", 0),
            "unmatched_candidate": statuses.get("rj_only", 0),
        },
        "features": {
            "comparisons": sum(row["comparisons"] for row in feature_metrics.values()),
            "violations": sum(row["violations"] for row in feature_metrics.values()),
            "max_normalized_delta": max(
                (row["max_normalized_delta"] for row in feature_metrics.values()),
                default=0.0,
            ),
        },
        "scores": {
            "comparisons": sum(row["comparisons"] for row in score_metrics.values()),
            "mismatches": sum(row["mismatches"] for row in score_metrics.values()),
            "max_abs_delta": max(
                (row["max_abs_delta"] for row in score_metrics.values()), default=0.0
            ),
        },
        "model_route": route,
        "tags": tags,
        "_feature_metrics": feature_metrics,
        "_score_metrics": score_metrics,
    }


def build_candidate_parity(
    lane_rows: dict[str, list[dict[str, str]]],
    *,
    created_utc: str | None = None,
    candidate_row_paths: dict[str, Path] | None = None,
    trace_bundle_by_lane: dict[str, dict[str, Path]] | None = None,
) -> dict[str, Any]:
    """Build the exact gate-compatible candidate-parity artifact."""
    expected = sorted(expected_lane_ids())
    if set(lane_rows) != set(expected) or len(lane_rows) != len(expected):
        raise RuntimeError(
            "candidate-parity input is not the exact 12 photon lanes: "
            f"missing={sorted(set(expected) - set(lane_rows))}, "
            f"unexpected={sorted(set(lane_rows) - set(expected))}"
        )

    detailed_lanes = [_candidate_parity_lane(lane_id, lane_rows[lane_id]) for lane_id in expected]
    lane_coverage = [
        {key: value for key, value in lane.items() if not key.startswith("_")}
        for lane in detailed_lanes
    ]

    population_names = (
        "reference_count",
        "candidate_count",
        "unmatched_reference",
        "unmatched_candidate",
    )
    population = {
        name: sum(lane["population"][name] for lane in detailed_lanes)
        for name in population_names
    }
    features = []
    for source_name in FEATURE_NAMES:
        name = FEATURE_NAME_ALIASES.get(source_name, source_name)
        rows = [lane["_feature_metrics"][name] for lane in detailed_lanes]
        features.append(
            {
                "name": name,
                "comparisons": sum(row["comparisons"] for row in rows),
                "violations": sum(row["violations"] for row in rows),
                "max_abs_delta": max(
                    (row["max_abs_delta"] for row in rows), default=0.0
                ),
                "max_normalized_delta": max(
                    (row["max_normalized_delta"] for row in rows), default=0.0
                ),
            }
        )
    scores = []
    for source_name in ("base_v3E", "base_E"):
        name = SCORE_NAME_ALIASES.get(source_name, source_name)
        rows = [lane["_score_metrics"][name] for lane in detailed_lanes]
        scores.append(
            {
                "name": name,
                "comparisons": sum(row["comparisons"] for row in rows),
                "mismatches": sum(row["mismatches"] for row in rows),
                "max_abs_delta": max(
                    (row["max_abs_delta"] for row in rows), default=0.0
                ),
            }
        )
    route = {
        name: sum(lane["model_route"][name] for lane in detailed_lanes)
        for name in ("comparisons", "mismatches")
    }
    tags = {
        tag_name: {
            name: sum(lane["tags"][tag_name][name] for lane in detailed_lanes)
            for name in ("comparisons", "mismatches")
        }
        for tag_name in ("common", "tight", "non_tight", "abcd")
    }
    result = {
        "schema": CANDIDATE_PARITY_SCHEMA,
        "created_utc": created_utc or datetime.now(timezone.utc).isoformat(),
        "name_aliases": {
            "features": dict(FEATURE_NAME_ALIASES),
            "scores": dict(SCORE_NAME_ALIASES),
        },
        "candidate_parity": {
            "tolerances": dict(CANDIDATE_PARITY_TOLERANCES),
            "lane_coverage": lane_coverage,
            "population": population,
            "features": features,
            "scores": scores,
            "model_route": route,
            "tags": tags,
        },
    }
    if candidate_row_paths is not None and trace_bundle_by_lane is None:
        raise RuntimeError(
            "candidate-row paths alone are inadmissible; provide the bound paired "
            "contract/trace/aggregate bundle"
        )
    if trace_bundle_by_lane is not None:
        result["trace_evidence"] = build_trace_evidence(trace_bundle_by_lane)
    return result


def candidate_identity(row: dict[str, str]) -> dict[str, Any]:
    return {
        "segment": row.get("segment", ""),
        "eventnumber": row.get("eventnumber", ""),
        "truth_track_id": row.get("truth_track_id", ""),
        "rj_cluster_index": row.get("rj_cluster_index", ""),
        "ppg12_cluster_index": row.get("ppg12_cluster_index", ""),
    }


def summarize_candidate_rows(rows: Iterable[dict[str, str]]) -> dict[str, Any]:
    rows = list(rows)
    statuses = Counter(row.get("match_status", "") for row in rows)
    matched = [row for row in rows if row.get("match_status") == "matched"]
    counts: dict[str, int] = {stage: 0 for stage in FIRST_DIVERGENCE_ORDER}
    examples: dict[str, dict[str, Any]] = {}

    population_failures = statuses.get("rj_only", 0) + statuses.get("ppg12_only", 0)
    if not matched:
        population_failures += 1
    counts["candidate_population"] = population_failures
    if population_failures:
        first = next(
            (row for row in rows if row.get("match_status") != "matched"),
            rows[0] if rows else {},
        )
        examples["candidate_population"] = {
            **candidate_identity(first),
            "match_status": first.get("match_status", "no_candidate_rows"),
        }

    for row in matched:
        for name in FEATURE_NAMES:
            left = row.get(f"{name}_rj")
            right = row.get(f"{name}_ppg12")
            if not close_feature(left, right):
                counts["features"] += 1
                examples.setdefault(
                    "features",
                    {
                        **candidate_identity(row),
                        "quantity": name,
                        "recoiljets": left,
                        "ppg12": right,
                    },
                )

        for label, left_name, right_name in SCORE_PAIRS:
            left = row.get(left_name)
            right = row.get(right_name)
            if not close_score(left, right):
                counts["scores"] += 1
                examples.setdefault(
                    "scores",
                    {
                        **candidate_identity(row),
                        "quantity": label,
                        "recoiljets": left,
                        "ppg12": right,
                    },
                )

        inferred = row.get("rj_inferred_stored_model", "")
        expected = row.get("ppg12_selected_model", "")
        exact_route_agree = row.get("model_route_agree") == "1"
        stored_route_agree = inferred == expected
        if not exact_route_agree or not stored_route_agree:
            counts["model_route"] += 1
            examples.setdefault(
                "model_route",
                {
                    **candidate_identity(row),
                    "recomputed_route": row.get("rj_selected_model", ""),
                    "stored_route_inferred": inferred,
                    "ppg12_route": expected,
                },
            )

        tag_flags = (
            row.get("common_agree") == "1",
            row.get("stored_common_agree") == "1",
            row.get("tag_agree") == "1",
            row.get("stored_tag_agree") == "1",
        )
        if not all(tag_flags):
            counts["tags"] += 1
            examples.setdefault(
                "tags",
                {
                    **candidate_identity(row),
                    "rj_recomputed_common": row.get("rj_recomputed_common_pass", ""),
                    "rj_stored_common": row.get("rj_stored_common_pass", ""),
                    "ppg12_common": row.get("ppg12_common_pass", ""),
                    "rj_recomputed_tag": row.get("rj_recomputed_tag", ""),
                    "rj_stored_tag": row.get("rj_stored_tag", ""),
                    "ppg12_tag": row.get("ppg12_recomputed_tag", ""),
                },
            )

        isolation_ok = True
        for label, left_name, right_name in ISOLATION_PAIRS:
            left = row.get(left_name)
            right = row.get(right_name)
            if not close_feature(left, right):
                isolation_ok = False
                examples.setdefault(
                    "isolation_abcd",
                    {
                        **candidate_identity(row),
                        "quantity": label,
                        "recoiljets": left,
                        "ppg12": right,
                    },
                )
        for left_name, right_name in (
            ("rj_is_iso", "ppg12_is_iso"),
            ("rj_is_noniso", "ppg12_is_noniso"),
        ):
            if row.get(left_name) != row.get(right_name):
                isolation_ok = False
                examples.setdefault(
                    "isolation_abcd",
                    {
                        **candidate_identity(row),
                        "quantity": f"{left_name}_vs_{right_name}",
                        "recoiljets": row.get(left_name, ""),
                        "ppg12": row.get(right_name, ""),
                    },
                )
        if row.get("abcd_agree") != "1" or row.get("stored_abcd_agree") != "1":
            isolation_ok = False
            examples.setdefault(
                "isolation_abcd",
                {
                    **candidate_identity(row),
                    "quantity": "ABCD region",
                    "rj_recomputed": row.get("rj_abcd_region", ""),
                    "rj_stored": row.get("rj_stored_abcd_region", ""),
                    "ppg12": row.get("ppg12_abcd_region", ""),
                },
            )
        if not isolation_ok:
            counts["isolation_abcd"] += 1

        finite_weight_values = (
            row.get("rj_event_weight"),
            row.get("rj_xsec_pb"),
            row.get("rj_xsec_weight"),
            row.get("rj_window_low"),
        )
        try:
            window_high = float(row.get("rj_window_high", "nan"))
        except (TypeError, ValueError):
            window_high = float("nan")
        # Photon20 intentionally owns an open-ended truth-response window, so
        # +infinity is valid for its upper edge while NaN is never valid.
        if any(finite_float(value) is None for value in finite_weight_values) or math.isnan(window_high):
            counts["weight_metadata"] += 1
            examples.setdefault(
                "weight_metadata",
                {
                    **candidate_identity(row),
                    "event_weight": row.get("rj_event_weight", ""),
                    "xsec_pb": row.get("rj_xsec_pb", ""),
                    "xsec_weight": row.get("rj_xsec_weight", ""),
                    "window_low": row.get("rj_window_low", ""),
                    "window_high": row.get("rj_window_high", ""),
                },
            )

    first_divergence = next(
        (stage for stage in FIRST_DIVERGENCE_ORDER if counts[stage] > 0),
        None,
    )
    return {
        "candidate_rows": len(rows),
        "match_status_counts": dict(sorted(statuses.items())),
        "stage_failure_counts": counts,
        "stage_examples": examples,
        "first_divergence": first_divergence,
        "weight_parity_scope": (
            "metadata completeness only; PPG12 stitched event weights are evaluated "
            "by the 32-lane closure gate"
        ),
        "fill_multiplicity_scope": (
            "PPG12 oracle multiplicity is exported; RecoilJets fill multiplicity is "
            "not encoded as a candidate-tree branch and remains a histogram-level gate"
        ),
    }


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def validate_plan_and_receipts(plan_path: Path, receipts_path: Path) -> list[dict[str, Any]]:
    plan = json.loads(plan_path.read_text())
    if plan.get("schema") != PLAN_SCHEMA:
        raise RuntimeError(f"unexpected plan schema: {plan.get('schema')}")
    lanes = plan.get("lanes")
    if not isinstance(lanes, list):
        raise RuntimeError("plan lanes must be a list")
    lane_ids = [lane.get("lane_id") for lane in lanes]
    if lane_ids != expected_lane_ids():
        raise RuntimeError(f"lane order/coverage mismatch: {lane_ids}")

    auth_payload = {
        key: plan[key]
        for key in (
            "schema", "campaign_tag", "output_root", "source_manifest",
            "paired_driver", "common", "lanes",
        )
    }
    expected_token = "RUN_PPG12_PAIRED_" + canonical_payload_sha256(auth_payload)
    if plan.get("submission_token") != expected_token:
        raise RuntimeError("plan authorization token does not bind its exact lane/assets")
    for label, item in (
        ("source_manifest", plan.get("source_manifest")),
        ("paired_driver", plan.get("paired_driver")),
    ):
        path = Path(str((item or {}).get("path", "")))
        _bound_file(item, path, f"plan {label}")
    common = plan.get("common", {})
    expected_common = {
        "setup_script", "apply_bdt", "apply_config", "base_e_model",
        "base_v3e_model", "npb_model", "tower_mask",
        "recoil_runtime_manifest", "recoil_config",
    }
    if set(common) != expected_common:
        raise RuntimeError("plan common asset role set differs")
    for role, item in common.items():
        path = Path(str(item.get("path", "")))
        _bound_file(item, path, f"plan common {role}")
    with receipts_path.open(newline="") as stream:
        receipts = list(csv.DictReader(stream, delimiter="\t"))
    if [row.get("lane_id") for row in receipts] != expected_lane_ids():
        raise RuntimeError("submission receipt lane order/coverage mismatch")
    receipt_by_lane = {row["lane_id"]: row for row in receipts}

    unique_paths: dict[str, set[str]] = {
        role: set() for role in (
            "runtime_contract", "candidate_csv", "executable_aggregate",
            "trace_csv", "response_trace_csv",
        )
    }
    unique_hashes = {role: set() for role in unique_paths}
    parsed: list[dict[str, Any]] = []
    for lane in lanes:
        lane_id = lane["lane_id"]
        receipt = receipt_by_lane[lane_id]
        if lane.get("execution") != "source_locked_paired_executable" or lane.get("rows") != 5:
            raise RuntimeError(f"{lane_id}: plan is not a five-row paired executable lane")
        _, sample_key, period, interaction = lane_id.split(":")
        if lane.get("sample") != "Photon" + sample_key.removeprefix("photon"):
            raise RuntimeError(f"{lane_id}: sample identity differs")
        if lane.get("period") != period or lane.get("interaction") != interaction.upper():
            raise RuntimeError(f"{lane_id}: period/interaction identity differs")
        sources = lane.get("sources", {})
        if set(sources) != {"ppg_macro", "g4_full_list", "truthjet_full_list"}:
            raise RuntimeError(f"{lane_id}: source role set differs")
        for role, item in sources.items():
            path = Path(str(item.get("path", "")))
            _bound_file(item, path, f"{lane_id} source {role}")
        if receipt.get("status") != "PASS":
            raise RuntimeError(f"{lane_id}: receipt is not PASS")
        if receipt.get("output_base") != lane.get("output_base"):
            raise RuntimeError(f"{lane_id}: receipt output_base differs from plan")
        expected_evidence = lane.get("expected_evidence", {})
        for receipt_key, expected_key in (
            ("runtime_contract", "runtime_contract"),
            ("candidate_csv", "candidate_csv"),
            ("executable_aggregate", "executable_aggregate"),
        ):
            if receipt.get(receipt_key) != expected_evidence.get(expected_key):
                raise RuntimeError(f"{lane_id}: {receipt_key} path differs from plan")
            path = Path(str(receipt[receipt_key])).resolve()
            if not path.is_file() or path.stat().st_size <= 0:
                raise RuntimeError(f"{lane_id}: missing receipt evidence {receipt_key}")
            if receipt.get(receipt_key + "_sha256") != sha256_file(path):
                raise RuntimeError(f"{lane_id}: {receipt_key} receipt hash differs")
        run_state = Path(str(expected_evidence.get("run_state", "")))
        if not run_state.is_file() or run_state.read_text().strip() != "PASS":
            raise RuntimeError(f"{lane_id}: paired run state is not PASS")
        contract_path = Path(receipt["runtime_contract"]).resolve()
        contract = json.loads(contract_path.read_text())
        contract_paths = contract.get("paths", {})
        bundle = {
            "runtime_contract": contract_path,
            "candidate_csv": Path(receipt["candidate_csv"]).resolve(),
            "executable_aggregate": Path(receipt["executable_aggregate"]).resolve(),
            "trace_csv": Path(str(contract_paths.get("ppg12_executable_trace", ""))).resolve(),
            "response_trace_csv": Path(str(contract_paths.get("ppg12_executable_response_trace", ""))).resolve(),
        }
        validate_lane_evidence_bundle(lane_id, bundle)
        for role, path in bundle.items():
            resolved, item_hash = str(path), sha256_file(path)
            if resolved in unique_paths[role] or item_hash in unique_hashes[role]:
                raise RuntimeError(f"{lane_id}: {role} is reused by another physical lane")
            unique_paths[role].add(resolved)
            unique_hashes[role].add(item_hash)
        enriched = dict(lane)
        enriched["trace_bundle"] = bundle
        parsed.append(enriched)
    return parsed


def audit_root(path: Path, tree_name: str) -> dict[str, Any]:
    if ROOT is None:
        raise RuntimeError("PyROOT is required; run under the sPHENIX analysis environment")
    size = path.stat().st_size
    if size <= MIN_ROOT_BYTES:
        raise RuntimeError(f"ROOT is not terminal/non-tiny: {size} bytes")
    root_file = ROOT.TFile.Open(str(path))
    if not root_file or root_file.IsZombie():
        raise RuntimeError("ROOT is unreadable or zombie")
    recovered = bool(root_file.TestBit(ROOT.TFile.kRecovered))
    if recovered:
        root_file.Close()
        raise RuntimeError("ROOT carries TFile::kRecovered")
    tree = root_file.Get(tree_name)
    if not tree:
        root_file.Close()
        raise RuntimeError(f"missing {tree_name}")
    branches = {branch.GetName() for branch in tree.GetListOfBranches()}
    missing = sorted(REQUIRED_RJ_BRANCHES - branches)
    entries = int(tree.GetEntries())
    root_file.Close()
    if missing:
        raise RuntimeError(f"missing candidate branches: {missing}")
    if entries <= 0:
        raise RuntimeError(f"{tree_name} contains zero rows")
    return {
        "path": str(path),
        "bytes": size,
        "sha256": sha256_file(path),
        "tree": tree_name,
        "tree_entries": entries,
        "required_branch_count": len(REQUIRED_RJ_BRANCHES),
        "missing_branches": [],
        "zombie": False,
        "recovered": False,
    }


def audit_large_oracle_root(path: Path) -> dict[str, Any]:
    """Identify the multi-GB preserved oracle without a 100+ GB hash pass."""
    if ROOT is None:
        raise RuntimeError("PyROOT is required; run under the sPHENIX analysis environment")
    stat = path.stat()
    root_file = ROOT.TFile.Open(str(path))
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"PPG12 oracle ROOT is unreadable or zombie: {path}")
    if root_file.TestBit(ROOT.TFile.kRecovered):
        root_file.Close()
        raise RuntimeError(f"PPG12 oracle ROOT carries TFile::kRecovered: {path}")
    tree = root_file.Get("slimtree")
    if not tree:
        root_file.Close()
        raise RuntimeError(f"PPG12 oracle ROOT is missing slimtree: {path}")
    entries = int(tree.GetEntries())
    uuid = str(root_file.GetUUID().AsString())
    root_file.Close()
    return {
        "path": str(path),
        "bytes": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "root_uuid": uuid,
        "tree": "slimtree",
        "tree_entries": entries,
        "sha256": None,
        "identity_note": (
            "Frozen by absolute path, byte size, mtime_ns, ROOT UUID, and tree entries; "
            "a full SHA-256 pass over the six 11-26 GB preserved files is deliberately "
            "outside this time-critical candidate audit."
        ),
    }


def find_single_root(output_base: Path) -> Path:
    roots = sorted(output_base.rglob("*.root"))
    if len(roots) != 1:
        raise RuntimeError(f"expected exactly one ROOT, found {len(roots)}")
    return roots[0]


def load_comparator(path: Path) -> Any:
    spec = importlib.util.spec_from_file_location("ppg12_same_cluster_oracle", path)
    if not spec or not spec.loader:
        raise RuntimeError(f"cannot import comparator: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def write_markdown(summary: dict[str, Any], path: Path) -> None:
    lines = [
        "# THE-97 12-lane photon oracle canary audit",
        "",
        f"- Status: **{summary['status']}**",
        f"- First divergence: `{summary.get('first_divergence') or 'none in same-event scope'}`",
        f"- Plan: `{summary['plan_json']}`",
        f"- Receipts: `{summary['receipt_tsv']}`",
        "- Weight parity and the final purity estimator remain the downstream 32-lane gate.",
        "",
        "| lane | ROOT rows | matched | RJ-only | PPG12-only | population | features | scores | route | tags | isolation/ABCD | first |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for lane in summary.get("lanes", []):
        candidate = lane.get("candidate_summary", {})
        statuses = candidate.get("match_status_counts", {})
        counts = candidate.get("stage_failure_counts", {})
        lines.append(
            f"| {lane['lane_id']} | {lane.get('root_audit', {}).get('tree_entries', 0)} | "
            f"{statuses.get('matched', 0)} | {statuses.get('rj_only', 0)} | "
            f"{statuses.get('ppg12_only', 0)} | {counts.get('candidate_population', 0)} | "
            f"{counts.get('features', 0)} | {counts.get('scores', 0)} | "
            f"{counts.get('model_route', 0)} | {counts.get('tags', 0)} | "
            f"{counts.get('isolation_abcd', 0)} | {candidate.get('first_divergence') or 'none'} |"
        )
    if summary.get("first_divergence_example"):
        lines.extend(
            [
                "",
                "## First divergence example",
                "",
                "```json",
                json.dumps(summary["first_divergence_example"], indent=2, sort_keys=True),
                "```",
            ]
        )
    if summary.get("errors"):
        lines.extend(["", "## Errors", ""])
        for error in summary["errors"]:
            lines.append(f"- `{error['lane_id']}`: {error['error']}")
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-json", required=True, type=Path)
    parser.add_argument("--receipt-tsv", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    if args.output_dir.exists():
        raise RuntimeError(f"output directory already exists: {args.output_dir}")
    lanes = validate_plan_and_receipts(args.plan_json, args.receipt_tsv)
    args.output_dir.mkdir(parents=True)

    summary: dict[str, Any] = {
        "schema": "ppg12-photon-oracle-canary-audit/v2",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "AUDITING",
        "plan_json": str(args.plan_json),
        "receipt_tsv": str(args.receipt_tsv),
        "lanes": [],
        "errors": [],
    }
    auditor_path = Path(__file__).with_name("audit_ppg12_recoiljets_paired_oracle.py")
    spec = importlib.util.spec_from_file_location("paired_postrun_auditor", auditor_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load paired auditor: {auditor_path}")
    paired_auditor = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = paired_auditor
    spec.loader.exec_module(paired_auditor)
    lane_rows: dict[str, list[dict[str, str]]] = {}
    bundles: dict[str, dict[str, Path]] = {}
    for lane in lanes:
        result: dict[str, Any] = {
            "lane_id": lane["lane_id"],
            "period": lane["period"],
            "interaction": lane["interaction"],
        }
        try:
            bundle = lane["trace_bundle"]
            paired_auditor.validate_postrun(bundle["runtime_contract"])
            rows = read_csv_rows(bundle["candidate_csv"])
            lane_rows[lane["lane_id"]] = rows
            bundles[lane["lane_id"]] = bundle
            result["candidate_summary"] = summarize_candidate_rows(rows)
            result["paired_evidence"] = {
                role: {"path": str(path), "sha256": sha256_file(path)}
                for role, path in bundle.items()
            }
        except Exception as exc:
            result["audit_error"] = str(exc)
            summary["errors"].append({"lane_id": lane["lane_id"], "error": str(exc)})
        summary["lanes"].append(result)

    if summary["errors"]:
        summary["status"] = "PAIRED_EVIDENCE_AUDIT_FAILED"
        summary["first_divergence"] = "paired_evidence"
        (args.output_dir / "audit_summary.json").write_text(
            json.dumps(summary, indent=2, sort_keys=True) + "\n"
        )
        write_markdown(summary, args.output_dir / "audit_summary.md")
        print(f"status=PAIRED_EVIDENCE_AUDIT_FAILED output={args.output_dir}")
        return 2
    candidate_parity_path = args.output_dir / "candidate_parity.json"
    candidate_parity = build_candidate_parity(
        lane_rows,
        created_utc=summary["created_utc"],
        trace_bundle_by_lane=bundles,
    )
    candidate_parity_path.write_text(
        json.dumps(candidate_parity, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    summary["candidate_parity_json"] = str(candidate_parity_path)
    summary["candidate_parity_sha256"] = sha256_file(candidate_parity_path)
    global_counts = Counter()
    for result in summary["lanes"]:
        global_counts.update(result["candidate_summary"]["stage_failure_counts"])
    summary["global_stage_failure_counts"] = dict(global_counts)
    first = next(
        (stage for stage in FIRST_DIVERGENCE_ORDER if global_counts[stage] > 0), None
    )
    summary["first_divergence"] = first
    summary["status"] = "SAME_EVENT_SCOPE_PASSED" if first is None else "DIVERGENCE_FOUND"
    exit_code = 0 if first is None else 1

    (args.output_dir / "audit_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    write_markdown(summary, args.output_dir / "audit_summary.md")
    print(
        f"status={summary['status']} first_divergence={summary.get('first_divergence')} "
        f"output={args.output_dir}"
    )
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
