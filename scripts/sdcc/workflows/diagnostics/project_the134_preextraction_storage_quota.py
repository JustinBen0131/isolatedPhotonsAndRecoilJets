#!/usr/bin/env python3
"""Fail-closed storage admission for the full THE-134 source extraction.

This controller is intentionally non-submitting.  It has three actions:

``measure``
    Rehash the corrected full-extraction evidence chain and turn thirteen
    explicitly declared per-row artifact ceilings into one deterministic
    storage manifest.

``snapshot``
    Read live Lustre quota state for the four logical storage domains bound by
    the corrected extraction plan.  Intended namespaces are checked for
    freshness, while quota queries use the nearest existing parent.

``project``
    Combine a passed measurement manifest with a fresh authoritative quota
    snapshot.  Logical domains sharing one physical quota domain are grouped
    before applying the >=20% byte and inode headroom gates.

None of these actions invokes Condor, creates a production namespace, grants
training/science/production authority, or promotes an artifact to CANONICAL.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping, Sequence


HERE = Path(__file__).resolve().parent
RESOLVER_PATH = HERE / "resolve_the134_full_multiview_extraction.py"
AMENDMENT_PATH = HERE / "build_the134_capacity_count_amendment.py"
CAPACITY_BINDING_PATH = (
    HERE / "build_the134_capacity_partition_binding.py"
)


def _load_local_module(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load local module: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


resolver = _load_local_module("the134_preextraction_resolver", RESOLVER_PATH)
amendment_tool = _load_local_module(
    "the134_preextraction_count_amendment", AMENDMENT_PATH
)
capacity_binding_tool = _load_local_module(
    "the134_preextraction_capacity_binding", CAPACITY_BINDING_PATH
)


MEASUREMENT_SPEC_SCHEMA = (
    "THE134_PREEXTRACTION_ARTIFACT_MEASUREMENT_SPEC_V1"
)
MEASUREMENT_SCHEMA = "THE134_PREEXTRACTION_ARTIFACT_MEASUREMENT_V1"
QUOTA_SNAPSHOT_SCHEMA = "THE134_PREEXTRACTION_LIVE_QUOTA_SNAPSHOT_V1"
STORAGE_MANIFEST_SCHEMA = "THE134_PREEXTRACTION_STORAGE_MANIFEST_V1"
STORAGE_CERTIFICATE_SCHEMA = (
    "THE134_PREEXTRACTION_STORAGE_CERTIFICATE_V1"
)
CONTROLLER_BUDGET_SCHEMA = (
    "THE134_FULL_EXTRACTION_CONTROLLER_DRY_MATERIALIZATION_BUDGET_V1"
)
CONTROLLER_DERIVATION_SCHEMA = (
    "THE134_FULL_EXTRACTION_CONTROLLER_DRY_MATERIALIZATION_BUDGET_DERIVATION_V1"
)
CONTROLLER_DERIVATION_STATUS = "PASS_EXACT_SERIALIZER_MEASUREMENT"
CONTROLLER_SERIALIZER_PATH = (
    HERE / "materialize_the134_full_multiview_extraction.py"
)
CONTROLLER_BUDGET_BUILDER_PATH = (
    HERE / "build_the134_controller_dry_materialization_budget.py"
)
CONTROLLER_ARTIFACT_SIZE_WIDTH_CEILING = 9_999_999_999_999_999_999
GATE = "PRE_P5A_SOURCE_EXTRACTION_STORAGE_ADMISSION"
AUTHORITY_STATE = "NON_SUBMITTING_EXTRACTION_STORAGE_INPUT_ONLY"

PASS_STATUS = "PASS_NON_SUBMITTING_PREEXTRACTION_ONLY"
FAIL_STATUS = "FAIL_CAPACITY"
BLOCKED_STATUS = "BLOCKED_INCOMPLETE_EVIDENCE"

EXIT_PASS = 0
EXIT_MALFORMED = 2
EXIT_CAPACITY = 3
EXIT_BLOCKED = 4

EXPECTED_ROW_COUNT = 13
EXPECTED_SOURCE_TUPLES = 129_998
EXPECTED_GROUP_SIZE = 7
EXPECTED_JOBS = 18_577
EXPECTED_OUTPUT_PAIRS = 18_577
EXPECTED_ANALYSIS_OUTPUTS = 18_577
EXPECTED_SIDECAR_OUTPUTS = 18_577
EXPECTED_ROOT_ARTIFACTS = 37_154
EXPECTED_RETAINED_ANALYSIS_OUTPUTS = 0
EXPECTED_DURABLE_ROOT_ARTIFACTS = EXPECTED_SIDECAR_OUTPUTS
EXPECTED_SOURCE_OCCURRENCES = 18_577
EXPECTED_REQUEST_MEMORY_MB = 8_000

MIN_HEADROOM_NUMERATOR = 1
MIN_HEADROOM_DENOMINATOR = 5
MIN_RETRY_RESERVE_NUMERATOR = 1
MIN_RETRY_RESERVE_DENOMINATOR = 10
MIN_WITNESS_CEILING_NUMERATOR = 5
MIN_WITNESS_CEILING_DENOMINATOR = 4
DEFAULT_SNAPSHOT_MAX_AGE_SECONDS = 900

STORAGE_DOMAINS = (
    "bulk_science",
    "scratch_control",
    "scratch_evidence",
    "scheduler_streams",
)
ARTIFACT_CLASS_CONTRACT = {
    "analysis_root": {
        "domain": "bulk_science",
        "count_per_job": 0,
        "basis": "EPHEMERAL_CONDOR_SCRATCH_VALIDATED_NOT_RETAINED",
        "requires_witness": False,
    },
    "training_sidecar_root": {
        "domain": "bulk_science",
        "count_per_job": 1,
        "basis": "CONSERVATIVE_CAPACITY_WITNESS_ENVELOPE",
        "requires_witness": True,
    },
    "evidence_record": {
        "domain": "scratch_evidence",
        "count_per_job": 1,
        "basis": "DECLARED_CONSERVATIVE_CONTROL_BOUND",
        "requires_witness": False,
    },
    "scheduler_stdout": {
        "domain": "scheduler_streams",
        "count_per_job": 1,
        "basis": "DECLARED_CONSERVATIVE_CONTROL_BOUND",
        "requires_witness": False,
    },
    "scheduler_stderr": {
        "domain": "scheduler_streams",
        "count_per_job": 1,
        "basis": "DECLARED_CONSERVATIVE_CONTROL_BOUND",
        "requires_witness": False,
    },
    "scheduler_event_log": {
        "domain": "scheduler_streams",
        "count_per_job": 1,
        "basis": "DECLARED_CONSERVATIVE_CONTROL_BOUND",
        "requires_witness": False,
    },
}
ANALYSIS_HEALTH_PROFILE = "THE134_ANALYSIS_ROOT_ARTIFACT_HEALTH_V1"
SIDECAR_HEALTH_PROFILE = "THE134_TRAINING_SIDECAR_ARTIFACT_HEALTH_V1"

EXPECTED_AUTHORITY = {
    "state": AUTHORITY_STATE,
    "submission_performed": False,
    "full_training_authority": 0,
    "full_extraction_authority": False,
    "science_freeze_authority": False,
    "broad_production_authority": False,
    "canonical_promotion": False,
}
AUTHORITY_FIELDS = {
    key: value for key, value in EXPECTED_AUTHORITY.items() if key != "state"
}

SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
PRINCIPAL_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_.-]{0,63}$")
SIZE_RE = re.compile(r"^([0-9]+(?:\.[0-9]+)?)([KMGTPEkmgtpe]?)$")
P5A_S_RE = re.compile(r"(?<![A-Z0-9])P5A[\s_-]*S(?![A-Z0-9])", re.IGNORECASE)

BLOCKER_ORDER = (
    "P5A_S_INPUT_FORBIDDEN",
    "PLAN_BINDING_DRIFT",
    "COUNT_CONTRACT_DRIFT",
    "ROW_ENVELOPE_MISSING",
    "CAPACITY_WITNESS_NOT_STORAGE_CEILING",
    "CONTROLLER_BUDGET_MISSING",
    "CONTROLLER_DERIVATION_MISSING",
    "ARTIFACT_CLASS_MISSING",
    "QUOTA_DOMAIN_UNBOUND",
    "QUOTA_SNAPSHOT_STALE",
    "QUOTA_USAGE_UNAUTHORITATIVE",
    "NAMESPACE_NOT_FRESH",
    "BYTE_HEADROOM_LT_20PCT",
    "INODE_HEADROOM_LT_20PCT",
)


class ProjectionError(RuntimeError):
    """Malformed, cryptographically inconsistent, or tampered input."""


def canonical_json_bytes(payload: Any) -> bytes:
    return (
        json.dumps(
            payload,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
        )
        + "\n"
    ).encode("utf-8")


def semantic_sha256(payload: Any) -> str:
    return hashlib.sha256(canonical_json_bytes(payload)).hexdigest()


def file_sha256(path: Path) -> str:
    try:
        with path.open("rb") as stream:
            digest = hashlib.sha256()
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
    except OSError as exc:
        raise ProjectionError(f"cannot read artifact: {path}") from exc
    return digest.hexdigest()


def require_sha256(value: object, label: str) -> str:
    if not isinstance(value, str) or SHA256_RE.fullmatch(value) is None:
        raise ProjectionError(f"{label} must be a lowercase SHA-256")
    return value


def require_mapping(value: object, label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ProjectionError(f"{label} must be an object")
    return value


def require_sequence(value: object, label: str) -> list[Any]:
    if not isinstance(value, list):
        raise ProjectionError(f"{label} must be an array")
    return value


def require_exact_keys(
    payload: Mapping[str, Any], expected: Iterable[str], label: str
) -> None:
    observed = set(payload)
    expected_set = set(expected)
    if observed != expected_set:
        raise ProjectionError(
            f"{label} key inventory differs: "
            f"missing={sorted(expected_set - observed)} "
            f"extra={sorted(observed - expected_set)}"
        )


def require_nonnegative_int(value: object, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ProjectionError(f"{label} must be a nonnegative integer")
    return value


def require_positive_int(value: object, label: str) -> int:
    result = require_nonnegative_int(value, label)
    if result == 0:
        raise ProjectionError(f"{label} must be greater than zero")
    return result


def require_exact_int(value: object, expected: int, label: str) -> None:
    observed = require_nonnegative_int(value, label)
    if observed != expected:
        raise ProjectionError(f"{label}={observed}, expected={expected}")


def require_absolute_remote_root(value: object, label: str) -> str:
    if (
        not isinstance(value, str)
        or not value.startswith("/")
        or "\x00" in value
        or ".." in Path(value).parts
    ):
        raise ProjectionError(f"{label} must be a safe absolute path")
    return value.rstrip("/") or "/"


def ordered_blockers(blockers: Iterable[str]) -> list[str]:
    unique = set(blockers)
    return [code for code in BLOCKER_ORDER if code in unique] + sorted(
        unique - set(BLOCKER_ORDER)
    )


def authority_surface() -> dict[str, Any]:
    return dict(EXPECTED_AUTHORITY)


def apply_authority(payload: dict[str, Any]) -> dict[str, Any]:
    payload["authority"] = authority_surface()
    payload.update(AUTHORITY_FIELDS)
    return payload


def validate_authority(payload: Mapping[str, Any], label: str) -> None:
    authority = require_mapping(payload.get("authority"), f"{label}.authority")
    if authority != EXPECTED_AUTHORITY:
        raise ProjectionError(f"{label} authority boundary differs")
    for field, expected in AUTHORITY_FIELDS.items():
        if payload.get(field) != expected:
            raise ProjectionError(f"{label}.{field} must remain {expected!r}")


def validate_all_authority_occurrences(payload: Any, label: str) -> None:
    if isinstance(payload, dict):
        for key, value in payload.items():
            if key in AUTHORITY_FIELDS and value != AUTHORITY_FIELDS[key]:
                raise ProjectionError(
                    f"{label}.{key} escalates non-submitting authority"
                )
            validate_all_authority_occurrences(value, f"{label}.{key}")
    elif isinstance(payload, list):
        for index, value in enumerate(payload):
            validate_all_authority_occurrences(value, f"{label}[{index}]")


def contains_p5a_s_input(payload: Any) -> bool:
    if isinstance(payload, dict):
        for key, value in payload.items():
            key_text = str(key)
            lowered = key_text.lower()
            if (
                "science_freeze_certificate" in lowered
                or P5A_S_RE.search(key_text) is not None
            ):
                return True
            if contains_p5a_s_input(value):
                return True
        return False
    if isinstance(payload, list):
        return any(contains_p5a_s_input(value) for value in payload)
    if isinstance(payload, str):
        return P5A_S_RE.search(payload) is not None
    return False


def artifact_ref(value: object, label: str) -> dict[str, str]:
    record = require_mapping(value, label)
    require_exact_keys(record, {"path", "sha256"}, label)
    path_text = record.get("path")
    if not isinstance(path_text, str) or not path_text:
        raise ProjectionError(f"{label}.path must be nonempty")
    return {
        "path": path_text,
        "sha256": require_sha256(record.get("sha256"), f"{label}.sha256"),
    }


def load_artifact(
    value: object, label: str, *, expected_schema: str | None = None
) -> tuple[dict[str, Any], dict[str, Any]]:
    record = artifact_ref(value, label)
    path = Path(record["path"])
    observed_sha256 = file_sha256(path)
    if observed_sha256 != record["sha256"]:
        raise ProjectionError(f"{label} SHA-256 differs")
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ProjectionError(f"{label} is not readable JSON") from exc
    payload = require_mapping(payload, label)
    if expected_schema is not None and payload.get("schema") != expected_schema:
        raise ProjectionError(
            f"{label}.schema={payload.get('schema')!r}, "
            f"expected={expected_schema!r}"
        )
    return payload, {
        "path": str(path),
        "sha256": observed_sha256,
        "size_bytes": path.stat().st_size,
    }


def verify_semantic_receipt(
    payload: Mapping[str, Any], field: str, label: str
) -> None:
    expected = require_sha256(payload.get(field), f"{label}.{field}")
    semantic_payload = dict(payload)
    del semantic_payload[field]
    if semantic_sha256(semantic_payload) != expected:
        raise ProjectionError(f"{label}.{field} differs")


def ceil_ratio(value: int, numerator: int, denominator: int) -> int:
    return (value * numerator + denominator - 1) // denominator


def expected_row_ids() -> list[str]:
    return [record["row_id"] for record in resolver.inventory_rows()]


def validate_plan_counts(plan: Mapping[str, Any]) -> tuple[list[dict[str, Any]], list[str]]:
    blockers: list[str] = []
    if (
        plan.get("schema") != resolver.PLAN_SCHEMA
        or plan.get("status") != "PREFLIGHT_PASS"
        or plan.get("submission_performed") is not False
    ):
        raise ProjectionError("corrected extraction plan schema/status differs")
    resolver.validate_preflight_authority_payload(plan, label="plan")
    rows = require_sequence(plan.get("rows"), "plan.rows")
    row_ids = [
        row.get("row_id") for row in rows if isinstance(row, dict)
    ]
    if row_ids != expected_row_ids() or len(rows) != EXPECTED_ROW_COUNT:
        blockers.append("COUNT_CONTRACT_DRIFT")

    normalized_rows: list[dict[str, Any]] = []
    for raw in rows:
        row = require_mapping(raw, "plan row")
        row_id = str(row.get("row_id", ""))
        input_contract = require_mapping(
            row.get("input_contract"), f"{row_id}.input_contract"
        )
        execution_contract = require_mapping(
            row.get("execution_contract"), f"{row_id}.execution_contract"
        )
        source_tuples = require_positive_int(
            input_contract.get("source_tuple_count"),
            f"{row_id}.source_tuple_count",
        )
        jobs = require_positive_int(
            input_contract.get("expected_job_count"),
            f"{row_id}.expected_job_count",
        )
        expected_jobs = ceil_ratio(
            source_tuples, 1, EXPECTED_GROUP_SIZE
        )
        count_fields = (
            "expected_chunk_count",
            "expected_output_pair_count",
            "expected_analysis_output_count",
            "expected_sidecar_output_count",
            "expected_source_occurrence_count",
        )
        if (
            input_contract.get("group_size") != EXPECTED_GROUP_SIZE
            or jobs != expected_jobs
            or any(input_contract.get(field) != jobs for field in count_fields)
            or input_contract.get("source_occurrences_per_output_pair") != 1
        ):
            blockers.append("COUNT_CONTRACT_DRIFT")
        materialization_environment = require_mapping(
            execution_contract.get("materialization_environment"),
            f"{row_id}.materialization_environment",
        )
        worker_environment = require_mapping(
            execution_contract.get("worker_environment"),
            f"{row_id}.worker_environment",
        )
        if (
            materialization_environment.get("RJ_REQUEST_MEMORY")
            != f"{EXPECTED_REQUEST_MEMORY_MB}MB"
            or worker_environment.get("RJ_THE134_EPHEMERAL_ANALYSIS_OUTPUT")
            != "1"
        ):
            blockers.append("COUNT_CONTRACT_DRIFT")
        system = row.get("system")
        if system not in {"pp", "auau"}:
            raise ProjectionError(f"{row_id}.system differs")
        normalized_rows.append(
            {
                "row_id": row_id,
                "system": system,
                "source_tuple_count": source_tuples,
                "expected_job_count": jobs,
            }
        )

    partition = require_mapping(
        plan.get("execution_partition"), "plan.execution_partition"
    )
    exact_partition = {
        "group_size": EXPECTED_GROUP_SIZE,
        "source_tuple_count": EXPECTED_SOURCE_TUPLES,
        "expected_chunk_count": EXPECTED_JOBS,
        "expected_job_count": EXPECTED_JOBS,
        "expected_output_pair_count": EXPECTED_OUTPUT_PAIRS,
        "expected_analysis_output_count": EXPECTED_ANALYSIS_OUTPUTS,
        "expected_sidecar_output_count": EXPECTED_SIDECAR_OUTPUTS,
        "expected_physical_root_artifact_count": EXPECTED_ROOT_ARTIFACTS,
        "expected_retained_analysis_output_count": (
            EXPECTED_RETAINED_ANALYSIS_OUTPUTS
        ),
        "expected_durable_root_artifact_count": (
            EXPECTED_DURABLE_ROOT_ARTIFACTS
        ),
        "expected_source_occurrence_count": EXPECTED_SOURCE_OCCURRENCES,
        "source_occurrences_per_output_pair": 1,
    }
    if (
        partition.get("schema") != resolver.PARTITION_SCHEMA
        or any(partition.get(key) != value for key, value in exact_partition.items())
        or sum(row["source_tuple_count"] for row in normalized_rows)
        != EXPECTED_SOURCE_TUPLES
        or sum(row["expected_job_count"] for row in normalized_rows)
        != EXPECTED_JOBS
    ):
        blockers.append("COUNT_CONTRACT_DRIFT")
    return normalized_rows, ordered_blockers(blockers)


def validate_preflight_receipt(
    payload: Mapping[str, Any],
    *,
    plan_sha256: str,
) -> None:
    if (
        payload.get("schema") != resolver.RECEIPT_SCHEMA
        or payload.get("status") != "PASS"
        or payload.get("submission_performed") is not False
        or payload.get("row_count") != EXPECTED_ROW_COUNT
        or payload.get("full_training_authority") != 0
    ):
        raise ProjectionError("corrected preflight receipt differs")
    resolver.validate_preflight_authority_payload(payload, label="preflight receipt")
    artifacts = require_mapping(
        payload.get("artifacts"), "preflight receipt.artifacts"
    )
    plan_record = require_mapping(artifacts.get("plan"), "preflight receipt plan")
    if plan_record.get("sha256") != plan_sha256:
        raise ProjectionError("preflight receipt plan SHA-256 differs")


def validate_amendment_readback(
    payload: Mapping[str, Any],
    *,
    amendment_artifact: Mapping[str, Any],
    amendment: Mapping[str, Any],
) -> None:
    require_exact_keys(
        payload,
        {
            "schema",
            "status",
            "authority_state",
            "amendment",
            "amendment_sha256",
            "amendment_semantic_sha256",
            "all_evidence_rehashed",
            "byte_exact_rebuild",
            "submission_performed",
            "full_training_authority",
            "full_extraction_authority",
        },
        "capacity-count amendment readback",
    )
    amendment_path = Path(str(amendment_artifact.get("path", ""))).absolute()
    readback_path = Path(str(payload.get("amendment", ""))).absolute()
    if (
        payload.get("schema") != amendment_tool.READBACK_SCHEMA
        or payload.get("status") != "PASS"
        or payload.get("authority_state") != amendment_tool.AUTHORITY_STATE
        or readback_path != amendment_path
        or payload.get("amendment_sha256") != amendment_artifact.get("sha256")
        or payload.get("amendment_semantic_sha256")
        != amendment.get("amendment_semantic_sha256")
        or payload.get("all_evidence_rehashed") is not True
        or payload.get("byte_exact_rebuild") is not True
        or payload.get("submission_performed") is not False
        or payload.get("full_training_authority") != 0
        or payload.get("full_extraction_authority") is not False
    ):
        raise ProjectionError("capacity-count amendment readback differs")


def validate_capacity_binding_readback(
    payload: Mapping[str, Any],
    *,
    binding_artifact: Mapping[str, Any],
    binding: Mapping[str, Any],
) -> None:
    require_exact_keys(
        payload,
        {
            "schema",
            "status",
            "authority_state",
            "binding",
            "binding_sha256",
            "binding_semantic_sha256",
            "all_evidence_rehashed",
            "byte_exact_rebuild",
            "legacy_evidence_consumed",
            "submission_performed",
            "full_training_authority",
            "full_extraction_authority",
        },
        "capacity partition binding readback",
    )
    binding_path = Path(str(binding_artifact.get("path", ""))).absolute()
    readback_path = Path(str(payload.get("binding", ""))).absolute()
    if (
        payload.get("schema") != capacity_binding_tool.READBACK_SCHEMA
        or payload.get("status") != "PASS"
        or payload.get("authority_state")
        != capacity_binding_tool.AUTHORITY_STATE
        or readback_path != binding_path
        or payload.get("binding_sha256") != binding_artifact.get("sha256")
        or payload.get("binding_semantic_sha256")
        != binding.get("binding_semantic_sha256")
        or payload.get("all_evidence_rehashed") is not True
        or payload.get("byte_exact_rebuild") is not True
        or payload.get("legacy_evidence_consumed") is not False
        or payload.get("submission_performed") is not False
        or payload.get("full_training_authority") != 0
        or payload.get("full_extraction_authority") is not False
    ):
        raise ProjectionError("capacity partition binding readback differs")


def derive_capacity_witnesses(
    amendment: Mapping[str, Any],
) -> dict[str, dict[str, Any]]:
    if amendment.get("schema") == capacity_binding_tool.BINDING_SCHEMA:
        binding = require_mapping(
            amendment.get("capacity_binding"),
            "capacity binding.capacity_binding",
        )
    else:
        binding = require_mapping(
            amendment.get("capacity_to_partition_binding"),
            "amendment.capacity_to_partition_binding",
        )
    current_ephemeral_binding = (
        amendment.get("schema") == capacity_binding_tool.BINDING_SCHEMA
    )
    capacity_rows = require_sequence(
        binding.get("selected_rows"), "capacity binding selected rows"
    )
    witnesses: dict[str, dict[str, Any]] = {}
    for raw in capacity_rows:
        row = require_mapping(raw, "capacity witness")
        row_id = str(row.get("row_id", ""))
        system = row.get("system")
        if system not in {"pp", "auau"} or row_id in witnesses:
            raise ProjectionError("capacity witness identity differs")
        analysis_field = (
            "analysis_health"
            if current_ephemeral_binding
            else "analysis_output"
        )
        analysis = require_mapping(
            row.get(analysis_field), f"{row_id}.{analysis_field}"
        )
        sidecar = require_mapping(
            row.get("sidecar_output"), f"{row_id}.sidecar_output"
        )
        witnesses[row_id] = {
            "witness_id": row_id,
            "row_id": row_id,
            "system": system,
            "tuple_count": EXPECTED_GROUP_SIZE,
            "request_memory_mb": EXPECTED_REQUEST_MEMORY_MB,
            "wall_time_seconds": require_positive_int(
                row.get("remote_wall_clock_seconds"),
                f"{row_id}.remote_wall_clock_seconds",
            ),
            "max_rss_mb": require_positive_int(
                row.get("memory_usage_mb"), f"{row_id}.memory_usage_mb"
            ),
            "storage_ceiling_authority": False,
            "analysis_root": (
                {
                    "path": str(analysis.get("path", "")),
                    "artifact_state": str(
                        analysis.get("artifact_state", "")
                    ),
                    "health_receipt_path": str(
                        analysis.get("health_receipt_path", "")
                    ),
                    "health_receipt_sha256": require_sha256(
                        analysis.get("health_receipt_sha256"),
                        f"{row_id}.analysis_health.health_receipt_sha256",
                    ),
                    "key_inventory_sha256": require_sha256(
                        analysis.get("key_inventory_sha256"),
                        f"{row_id}.analysis_health.key_inventory_sha256",
                    ),
                    "size_bytes": require_positive_int(
                        analysis.get("size_bytes"),
                        f"{row_id}.analysis_health.size_bytes",
                    ),
                    "health_profile": ANALYSIS_HEALTH_PROFILE,
                    "health_status": "PASS",
                }
                if current_ephemeral_binding
                else {
                    "path": str(analysis.get("path", "")),
                    "sha256": require_sha256(
                        analysis.get("sha256"),
                        f"{row_id}.analysis_output.sha256",
                    ),
                    "size_bytes": require_positive_int(
                        analysis.get("size_bytes"),
                        f"{row_id}.analysis_output.size_bytes",
                    ),
                    "health_profile": ANALYSIS_HEALTH_PROFILE,
                    "health_status": "PASS",
                }
            ),
            "training_sidecar_root": {
                "path": str(sidecar.get("path", "")),
                "sha256": require_sha256(
                    sidecar.get("sha256"), f"{row_id}.sidecar_output.sha256"
                ),
                "size_bytes": require_positive_int(
                    sidecar.get("size_bytes"),
                    f"{row_id}.sidecar_output.size_bytes",
                ),
                "health_profile": SIDECAR_HEALTH_PROFILE,
                "health_status": "PASS",
            },
        }
    if (
        set(witnesses) != set(amendment_tool.SELECTED_CAPACITY_ROWS)
        or {record["system"] for record in witnesses.values()} != {"pp", "auau"}
    ):
        raise ProjectionError("capacity witness inventory differs")
    return witnesses


def load_revalidated_amendment_chain(
    amendment_ref: object,
    readback_ref: object,
) -> tuple[
    dict[str, Any],
    dict[str, Any],
    dict[str, Any],
    dict[str, dict[str, Any]],
]:
    amendment, amendment_artifact = load_artifact(
        amendment_ref, "capacity evidence binding"
    )
    schema = amendment.get("schema")
    if schema == amendment_tool.AMENDMENT_SCHEMA:
        try:
            amendment = amendment_tool.validate_amendment_payload(amendment)
        except Exception as exc:
            raise ProjectionError(
                f"capacity-count amendment validation failed: {exc}"
            ) from exc
        readback, readback_artifact = load_artifact(
            readback_ref,
            "capacity-count amendment readback",
            expected_schema=amendment_tool.READBACK_SCHEMA,
        )
        validate_amendment_readback(
            readback,
            amendment_artifact=amendment_artifact,
            amendment=amendment,
        )
    elif schema == capacity_binding_tool.BINDING_SCHEMA:
        try:
            amendment = capacity_binding_tool.validate_binding_payload(
                amendment
            )
        except Exception as exc:
            raise ProjectionError(
                f"capacity partition binding validation failed: {exc}"
            ) from exc
        readback, readback_artifact = load_artifact(
            readback_ref,
            "capacity partition binding readback",
            expected_schema=capacity_binding_tool.READBACK_SCHEMA,
        )
        validate_capacity_binding_readback(
            readback,
            binding_artifact=amendment_artifact,
            binding=amendment,
        )
    else:
        raise ProjectionError(
            "capacity evidence binding schema is neither the preserved "
            "legacy amendment nor the current partition binding"
        )
    return (
        amendment,
        amendment_artifact,
        readback_artifact,
        derive_capacity_witnesses(amendment),
    )


def validate_chain(
    spec: Mapping[str, Any],
) -> tuple[
    dict[str, Any],
    list[dict[str, Any]],
    dict[str, dict[str, Any]],
    dict[str, Any],
    list[str],
]:
    bindings = require_mapping(spec.get("bindings"), "measurement bindings")
    require_exact_keys(
        bindings,
        {
            "plan",
            "preflight_receipt",
            "capacity_count_amendment",
            "capacity_count_amendment_readback",
            "bundle_manifest",
            "materialization_receipt",
        },
        "measurement bindings",
    )
    plan, plan_artifact = load_artifact(
        bindings["plan"], "corrected extraction plan", expected_schema=resolver.PLAN_SCHEMA
    )
    rows, blockers = validate_plan_counts(plan)
    preflight, preflight_artifact = load_artifact(
        bindings["preflight_receipt"],
        "corrected preflight receipt",
        expected_schema=resolver.RECEIPT_SCHEMA,
    )
    validate_preflight_receipt(
        preflight, plan_sha256=plan_artifact["sha256"]
    )
    (
        amendment,
        amendment_artifact,
        readback_artifact,
        witnesses,
    ) = load_revalidated_amendment_chain(
        bindings["capacity_count_amendment"],
        bindings["capacity_count_amendment_readback"],
    )
    _bundle, bundle_artifact = load_artifact(
        bindings["bundle_manifest"], "immutable bundle manifest"
    )
    _materialization, materialization_artifact = load_artifact(
        bindings["materialization_receipt"], "materialization receipt"
    )

    is_current_binding = (
        amendment.get("schema") == capacity_binding_tool.BINDING_SCHEMA
    )
    corrected = require_mapping(
        amendment.get("preflight")
        if is_current_binding
        else amendment.get("corrected_preflight"),
        "capacity binding preflight",
    )
    immutable = require_mapping(
        amendment.get("immutable_authority"),
        "capacity binding immutable_authority",
    )
    expected_bindings = {
        "plan": plan_artifact["sha256"],
        "preflight_receipt": preflight_artifact["sha256"],
    }
    for field, expected_sha in expected_bindings.items():
        record = require_mapping(
            corrected.get(field), f"amendment.corrected_preflight.{field}"
        )
        if record.get("sha256") != expected_sha:
            blockers.append("PLAN_BINDING_DRIFT")
    immutable_bindings = {
        "bundle_manifest": bundle_artifact["sha256"],
        "materialization_receipt": materialization_artifact["sha256"],
    }
    for field, expected_sha in immutable_bindings.items():
        record = require_mapping(
            immutable.get(field), f"amendment.immutable_authority.{field}"
        )
        if record.get("sha256") != expected_sha:
            blockers.append("PLAN_BINDING_DRIFT")

    if is_current_binding:
        count = require_mapping(
            amendment.get("count_contract"),
            "capacity binding count_contract",
        )
        exact_count_values = {
            "row_count": EXPECTED_ROW_COUNT,
            "group_size": EXPECTED_GROUP_SIZE,
            "source_tuple_count": EXPECTED_SOURCE_TUPLES,
            "expected_chunk_count": EXPECTED_JOBS,
            "expected_job_count": EXPECTED_JOBS,
            "expected_output_pair_count": EXPECTED_OUTPUT_PAIRS,
            "expected_analysis_output_count": EXPECTED_ANALYSIS_OUTPUTS,
            "expected_sidecar_output_count": EXPECTED_SIDECAR_OUTPUTS,
            "expected_physical_root_artifact_count": EXPECTED_ROOT_ARTIFACTS,
            "expected_source_occurrence_count": EXPECTED_SOURCE_OCCURRENCES,
            "source_occurrences_per_output_pair": 1,
            "request_memory_mb": EXPECTED_REQUEST_MEMORY_MB,
        }
    else:
        count = require_mapping(
            amendment.get("count_correction"),
            "amendment.count_correction",
        )
        exact_count_values = {
            "group_size": EXPECTED_GROUP_SIZE,
            "source_tuple_count": EXPECTED_SOURCE_TUPLES,
            "corrected_chunk_count": EXPECTED_JOBS,
            "corrected_job_count": EXPECTED_JOBS,
            "corrected_output_pair_count": EXPECTED_OUTPUT_PAIRS,
            "corrected_analysis_output_count": EXPECTED_ANALYSIS_OUTPUTS,
            "corrected_sidecar_output_count": EXPECTED_SIDECAR_OUTPUTS,
            "physical_root_artifact_count": EXPECTED_ROOT_ARTIFACTS,
            "corrected_source_occurrence_count": EXPECTED_SOURCE_OCCURRENCES,
            "source_occurrences_per_output_pair": 1,
        }
    if any(count.get(key) != value for key, value in exact_count_values.items()):
        blockers.append("COUNT_CONTRACT_DRIFT")

    capacity_partition = (
        require_mapping(
            amendment.get("capacity_binding"),
            "capacity binding.capacity_binding",
        )
        if is_current_binding
        else None
    )
    if capacity_partition is not None and (
        capacity_partition.get("execution_partition_sha256")
        != corrected.get("execution_partition_sha256")
        or capacity_partition.get("partition_artifact_sha256")
        != require_mapping(
            corrected.get("partition"),
            "capacity binding preflight partition",
        ).get("sha256")
    ):
        blockers.append("PLAN_BINDING_DRIFT")
    bindings_out = {
        "plan": plan_artifact,
        "preflight_receipt": preflight_artifact,
        "capacity_count_amendment": amendment_artifact,
        "capacity_count_amendment_readback": readback_artifact,
        "bundle_manifest": bundle_artifact,
        "materialization_receipt": materialization_artifact,
        "execution_partition_sha256": require_sha256(
            (
                capacity_partition.get("execution_partition_sha256")
                if capacity_partition is not None
                else corrected.get("execution_partition_sha256")
            ),
            "capacity binding execution_partition_sha256",
        ),
        "partition_artifact_sha256": require_sha256(
            (
                capacity_partition.get("partition_artifact_sha256")
                if capacity_partition is not None
                else require_mapping(
                    corrected.get("partition"),
                    "amendment corrected partition",
                ).get("sha256")
            ),
            "capacity binding partition SHA-256",
        ),
    }
    return plan, rows, witnesses, bindings_out, ordered_blockers(blockers)


def validate_controller_budget(
    value: object,
) -> tuple[dict[str, Any] | None, list[str]]:
    if value is None:
        return None, ["CONTROLLER_BUDGET_MISSING"]
    try:
        payload, artifact = load_artifact(
            value,
            "controller dry-materialization budget",
            expected_schema=CONTROLLER_BUDGET_SCHEMA,
        )
    except ProjectionError as exc:
        if "cannot read artifact" in str(exc):
            return None, ["CONTROLLER_BUDGET_MISSING"]
        raise
    expected_keys = {
        "schema",
        "status",
        "submission_performed",
        "expected_job_count",
        "storage_budget",
        "full_training_authority",
        "full_extraction_authority",
    }
    require_exact_keys(payload, expected_keys, "controller budget receipt")
    if (
        payload.get("status") != "PASS"
        or payload.get("submission_performed") is not False
        or payload.get("full_training_authority") != 0
        or payload.get("full_extraction_authority") is not False
        or payload.get("expected_job_count") != EXPECTED_JOBS
    ):
        return None, ["CONTROLLER_BUDGET_MISSING"]
    budget = require_mapping(payload.get("storage_budget"), "controller storage budget")
    require_exact_keys(
        budget,
        {"fixed_bytes", "fixed_inodes", "bytes_per_job", "inodes_per_job"},
        "controller storage budget",
    )
    normalized = {
        field: require_nonnegative_int(
            budget.get(field), f"controller storage budget.{field}"
        )
        for field in (
            "fixed_bytes",
            "fixed_inodes",
            "bytes_per_job",
            "inodes_per_job",
        )
    }
    if normalized["bytes_per_job"] == 0 or normalized["inodes_per_job"] == 0:
        return None, ["CONTROLLER_BUDGET_MISSING"]
    return {
        "receipt": artifact,
        "expected_job_count": EXPECTED_JOBS,
        "storage_domain": "scratch_control",
        **normalized,
        "projected_bytes": (
            normalized["fixed_bytes"]
            + EXPECTED_JOBS * normalized["bytes_per_job"]
        ),
        "projected_inodes": (
            normalized["fixed_inodes"]
            + EXPECTED_JOBS * normalized["inodes_per_job"]
        ),
    }, []


def _require_matching_artifact(
    value: object,
    expected: Mapping[str, Any],
    label: str,
) -> dict[str, Any]:
    observed = validate_materialized_artifact_record(value, label)
    expected_path = Path(str(expected.get("path", "")))
    if (
        Path(observed["path"]) != expected_path
        or observed["sha256"] != expected.get("sha256")
        or (
            "size_bytes" in expected
            and observed["size_bytes"] != expected.get("size_bytes")
        )
    ):
        raise ProjectionError(f"{label} differs from revalidated chain")
    return observed


def validate_controller_derivation(
    value: object,
    *,
    controller: Mapping[str, Any] | None,
    plan: Mapping[str, Any],
    bindings: Mapping[str, Any],
) -> tuple[dict[str, Any] | None, list[str]]:
    if value is None:
        return None, ["CONTROLLER_DERIVATION_MISSING"]
    if controller is None:
        return None, ["CONTROLLER_BUDGET_MISSING"]
    try:
        payload, artifact = load_artifact(
            value,
            "controller dry-materialization derivation",
            expected_schema=CONTROLLER_DERIVATION_SCHEMA,
        )
    except ProjectionError as exc:
        if "cannot read artifact" in str(exc):
            return None, ["CONTROLLER_DERIVATION_MISSING"]
        raise
    require_exact_keys(
        payload,
        {
            "schema",
            "status",
            "submission_performed",
            "authority",
            "serializer",
            "inputs",
            "measurements",
            "budget_derivation",
            "future_bindings",
            "budget_receipt",
            "boundaries",
            "derivation_semantic_sha256",
        },
        "controller derivation",
    )
    verify_semantic_receipt(
        payload,
        "derivation_semantic_sha256",
        "controller derivation",
    )
    authority = require_mapping(
        payload.get("authority"), "controller derivation authority"
    )
    if (
        payload.get("status") != CONTROLLER_DERIVATION_STATUS
        or payload.get("submission_performed") is not False
        or authority
        != {
            "state": "NON_SUBMITTING_CONTROLLER_STORAGE_INPUT_ONLY",
            "full_training_authority": 0,
            "full_extraction_authority": False,
            "science_freeze_authority": False,
            "broad_production_authority": False,
            "canonical_promotion": False,
        }
    ):
        raise ProjectionError("controller derivation authority/status differs")

    receipt = require_mapping(
        controller.get("receipt"), "normalized controller budget receipt"
    )
    _require_matching_artifact(
        payload.get("budget_receipt"),
        receipt,
        "controller derivation budget receipt",
    )

    serializer = require_mapping(
        payload.get("serializer"), "controller derivation serializer"
    )
    require_exact_keys(
        serializer,
        {"path", "sha256", "row_function", "job_function", "manifest_function"},
        "controller derivation serializer",
    )
    serializer_path = CONTROLLER_SERIALIZER_PATH.absolute()
    if (
        Path(str(serializer.get("path", ""))) != serializer_path
        or serializer.get("sha256") != file_sha256(serializer_path)
        or serializer.get("row_function") != "staged_row"
        or serializer.get("job_function") != "staged_job"
        or serializer.get("manifest_function") != "build_staged_artifacts"
    ):
        raise ProjectionError("controller derivation serializer differs")

    inputs = require_mapping(
        payload.get("inputs"), "controller derivation inputs"
    )
    require_exact_keys(
        inputs,
        {
            "plan",
            "preflight_receipt",
            "immutable_bundle",
            "immutable_materialization",
            "source_authority",
            "source_partition",
            "duplicate_fingerprint_sha256",
            "execution_fingerprint_sha256",
        },
        "controller derivation inputs",
    )
    for derivation_key, binding_key in (
        ("plan", "plan"),
        ("preflight_receipt", "preflight_receipt"),
        ("immutable_bundle", "bundle_manifest"),
        ("immutable_materialization", "materialization_receipt"),
    ):
        expected = require_mapping(
            bindings.get(binding_key), f"revalidated {binding_key} binding"
        )
        _require_matching_artifact(
            inputs.get(derivation_key),
            expected,
            f"controller derivation {derivation_key}",
        )

    plan_inputs = require_mapping(
        plan.get("input_manifests"), "revalidated plan input manifests"
    )
    source_ref = require_mapping(
        plan_inputs.get("sources"), "revalidated plan source authority"
    )
    source_payload, source_artifact = load_artifact(
        {
            "path": source_ref.get("path"),
            "sha256": source_ref.get("sha256"),
        },
        "revalidated plan source authority",
    )
    if (
        source_payload.get("schema") != resolver.SOURCE_SCHEMA
        or source_payload.get("status") != "PASS"
    ):
        raise ProjectionError("revalidated plan source authority differs")
    _require_matching_artifact(
        inputs.get("source_authority"),
        source_artifact,
        "controller derivation source authority",
    )

    partition_record = _require_matching_artifact(
        inputs.get("source_partition"),
        require_mapping(
            inputs.get("source_partition"),
            "controller derivation partition",
        ),
        "controller derivation source partition",
    )
    plan_artifact = require_mapping(
        bindings.get("plan"), "revalidated plan binding"
    )
    partition_contract = require_mapping(
        plan.get("execution_partition"), "revalidated execution partition"
    )
    plan_partition = require_mapping(
        partition_contract.get("partition_artifact"),
        "revalidated plan partition artifact",
    )
    expected_partition_path = (
        Path(str(plan_artifact.get("path", ""))).absolute().parent
        / str(plan_partition.get("name", ""))
    )
    if (
        Path(partition_record["path"]) != expected_partition_path
        or partition_record["sha256"]
        != bindings.get("partition_artifact_sha256")
        or partition_record["sha256"] != plan_partition.get("sha256")
        or partition_contract.get("expected_job_count") != EXPECTED_JOBS
        or resolver.canonical_sha256(partition_contract)
        != bindings.get("execution_partition_sha256")
    ):
        raise ProjectionError("controller derivation source partition differs")

    if (
        inputs.get("duplicate_fingerprint_sha256")
        != plan.get("duplicate_fingerprint_sha256")
        or inputs.get("execution_fingerprint_sha256")
        != plan.get("execution_fingerprint_sha256")
    ):
        raise ProjectionError("controller derivation plan fingerprints differ")
    require_sha256(
        inputs.get("duplicate_fingerprint_sha256"),
        "controller derivation duplicate fingerprint",
    )
    require_sha256(
        inputs.get("execution_fingerprint_sha256"),
        "controller derivation execution fingerprint",
    )

    measurements = require_mapping(
        payload.get("measurements"), "controller derivation measurements"
    )
    require_exact_keys(
        measurements,
        {
            "row_record_count",
            "row_manifest_bytes",
            "job_record_count",
            "job_manifest_exact_bytes",
            "job_record_min_bytes",
            "job_record_max_bytes",
            "job_record_ceiling_bytes",
            "submit_description_bytes",
            "manifest_envelope_bytes",
            "exact_envelope_total_bytes",
            "exact_envelope_total_inodes",
            "fixed_point_iterations",
        },
        "controller derivation measurements",
    )
    numeric = {
        key: require_nonnegative_int(
            measurements.get(key), f"controller derivation measurements.{key}"
        )
        for key in measurements
    }
    if (
        numeric["row_record_count"] != EXPECTED_ROW_COUNT
        or numeric["job_record_count"] != EXPECTED_JOBS
        or numeric["job_record_max_bytes"] != controller["bytes_per_job"]
        or numeric["job_record_min_bytes"] <= 0
        or numeric["job_record_min_bytes"] > numeric["job_record_max_bytes"]
        or numeric["job_record_ceiling_bytes"]
        != EXPECTED_JOBS * controller["bytes_per_job"]
        or numeric["job_manifest_exact_bytes"]
        > numeric["job_record_ceiling_bytes"]
        or numeric["fixed_point_iterations"] <= 0
        or numeric["fixed_point_iterations"] > 32
        or numeric["exact_envelope_total_inodes"] != 4
    ):
        raise ProjectionError("controller derivation measurements differ")

    derived = require_mapping(
        payload.get("budget_derivation"), "controller budget derivation"
    )
    require_exact_keys(
        derived,
        {
            "fixed_bytes",
            "fixed_inodes",
            "bytes_per_job",
            "inodes_per_job",
            "expected_job_count",
            "projected_bytes",
            "projected_inodes",
        },
        "controller budget derivation",
    )
    bases = {
        "fixed_bytes": (
            "exact row manifest + exact submit description + "
            "maximum-width future-binding materialization manifest"
        ),
        "fixed_inodes": "exact dry-materialization fixed artifact inventory",
        "bytes_per_job": "maximum exact staged_job canonical JSON record bytes",
        "inodes_per_job": "conservative one logical controller record per job",
        "expected_job_count": "validated resolver partition record count",
    }
    expected_values = {
        "fixed_bytes": controller["fixed_bytes"],
        "fixed_inodes": controller["fixed_inodes"],
        "bytes_per_job": controller["bytes_per_job"],
        "inodes_per_job": controller["inodes_per_job"],
        "expected_job_count": EXPECTED_JOBS,
    }
    for key, expected_value in expected_values.items():
        record = require_mapping(derived.get(key), f"controller derivation {key}")
        require_exact_keys(record, {"value", "basis"}, f"controller derivation {key}")
        if record.get("value") != expected_value or record.get("basis") != bases[key]:
            raise ProjectionError(f"controller derivation {key} differs")
    if (
        derived.get("projected_bytes") != controller["projected_bytes"]
        or derived.get("projected_inodes") != controller["projected_inodes"]
        or controller["fixed_bytes"]
        != (
            numeric["row_manifest_bytes"]
            + numeric["submit_description_bytes"]
            + numeric["manifest_envelope_bytes"]
        )
        or numeric["exact_envelope_total_bytes"]
        != numeric["job_manifest_exact_bytes"] + controller["fixed_bytes"]
        or numeric["exact_envelope_total_bytes"] > controller["projected_bytes"]
        or numeric["exact_envelope_total_inodes"] > controller["projected_inodes"]
    ):
        raise ProjectionError("controller derivation arithmetic differs")

    future = require_mapping(
        payload.get("future_bindings"), "controller derivation future bindings"
    )
    require_exact_keys(
        future,
        {
            "budget_output_path",
            "storage_certificate_path",
            "artifact_size_width_ceiling",
        },
        "controller derivation future bindings",
    )
    storage_path = Path(str(future.get("storage_certificate_path", "")))
    if (
        Path(str(future.get("budget_output_path", "")))
        != Path(str(receipt.get("path", "")))
        or not storage_path.is_absolute()
        or future.get("artifact_size_width_ceiling")
        != CONTROLLER_ARTIFACT_SIZE_WIDTH_CEILING
    ):
        raise ProjectionError("controller derivation future bindings differ")
    if payload.get("boundaries") != [
        "No Condor submission or job control.",
        "No output, evidence, or submit namespace creation.",
        "No scientific, model, working-point, production, or CANONICAL authority.",
    ]:
        raise ProjectionError("controller derivation boundaries differ")
    verify_exact_controller_derivation(
        payload,
        controller=controller,
        plan=plan,
        bindings=bindings,
    )
    return artifact, []


def verify_exact_controller_derivation(
    payload: Mapping[str, Any],
    *,
    controller: Mapping[str, Any],
    plan: Mapping[str, Any],
    bindings: Mapping[str, Any],
) -> None:
    """Replay the exact production serializers and reject self-authored claims.

    The strict budget and its derivation are produced together, so arithmetic
    consistency and a semantic seal alone cannot establish that the claimed
    measurements came from the current serializers.  Reopen the pinned input
    chain with the dedicated non-submitting builder, execute all 18,577 exact
    staged-job serializations plus the fixed-point manifest envelope, and
    require byte-for-byte equality with both receipts.
    """

    builder = _load_local_module(
        "the134_controller_budget_revalidator",
        CONTROLLER_BUDGET_BUILDER_PATH,
    )
    plan_binding = require_mapping(
        bindings.get("plan"), "revalidated plan binding"
    )
    preflight_binding = require_mapping(
        bindings.get("preflight_receipt"),
        "revalidated preflight receipt binding",
    )
    try:
        replay_plan, replay_plan_artifact = builder.load_pinned_json(
            Path(str(plan_binding.get("path", ""))),
            str(plan_binding.get("sha256", "")),
            "revalidated extraction plan",
        )
        replay_receipt, replay_receipt_artifact = builder.load_pinned_json(
            Path(str(preflight_binding.get("path", ""))),
            str(preflight_binding.get("sha256", "")),
            "revalidated resolver preflight receipt",
        )
        if replay_plan != dict(plan):
            raise ProjectionError(
                "controller derivation replay plan differs from validated plan"
            )
        context = builder.validate_input_chain(
            plan=replay_plan,
            plan_path=Path(str(plan_binding["path"])),
            plan_artifact=replay_plan_artifact,
            receipt=replay_receipt,
            receipt_artifact=replay_receipt_artifact,
        )
        future = require_mapping(
            payload.get("future_bindings"),
            "controller derivation future bindings",
        )
        expected_budget, expected_derivation = builder.derive_budget(
            context,
            budget_output=Path(
                str(
                    require_mapping(
                        controller.get("receipt"),
                        "normalized controller budget receipt",
                    )["path"]
                )
            ),
            future_storage_certificate=Path(
                str(future["storage_certificate_path"])
            ),
        )
    except ProjectionError:
        raise
    except (builder.BudgetError, FileNotFoundError) as exc:
        raise ProjectionError(
            f"controller derivation exact serializer replay failed: {exc}"
        ) from exc

    controller_receipt = require_mapping(
        controller.get("receipt"), "normalized controller budget receipt"
    )
    observed_budget, _ = load_artifact(
        {
            "path": controller_receipt.get("path"),
            "sha256": controller_receipt.get("sha256"),
        },
        "controller budget exact serializer replay",
        expected_schema=CONTROLLER_BUDGET_SCHEMA,
    )
    if observed_budget != expected_budget:
        raise ProjectionError(
            "controller budget differs from exact serializer replay"
        )
    if dict(payload) != expected_derivation:
        raise ProjectionError(
            "controller derivation differs from exact serializer replay"
        )


def normalize_row_envelopes(
    raw_envelopes: object,
    *,
    rows: Sequence[Mapping[str, Any]],
    witnesses: Mapping[str, Mapping[str, Any]],
) -> tuple[list[dict[str, Any]], dict[str, dict[str, int]], list[str]]:
    blockers: list[str] = []
    raw_records = require_sequence(raw_envelopes, "row_envelopes")
    raw_by_id: dict[str, dict[str, Any]] = {}
    for raw in raw_records:
        record = require_mapping(raw, "row envelope")
        require_exact_keys(record, {"row_id", "artifact_classes"}, "row envelope")
        row_id = str(record.get("row_id", ""))
        if not row_id or row_id in raw_by_id:
            raise ProjectionError("row envelope identity is missing or duplicated")
        raw_by_id[row_id] = record

    expected_ids = [str(row["row_id"]) for row in rows]
    if set(raw_by_id) != set(expected_ids):
        blockers.append("ROW_ENVELOPE_MISSING")

    domain_totals = {
        domain: {"base_bytes": 0, "base_inodes": 0}
        for domain in STORAGE_DOMAINS
    }
    normalized: list[dict[str, Any]] = []
    for row in rows:
        row_id = str(row["row_id"])
        raw = raw_by_id.get(row_id)
        if raw is None:
            continue
        classes_raw = require_sequence(
            raw.get("artifact_classes"), f"{row_id}.artifact_classes"
        )
        classes_by_name: dict[str, dict[str, Any]] = {}
        for class_raw in classes_raw:
            class_record = require_mapping(
                class_raw, f"{row_id} artifact class"
            )
            require_exact_keys(
                class_record,
                {
                    "artifact_class",
                    "storage_domain",
                    "count_per_job",
                    "bytes_per_artifact_ceiling",
                    "basis",
                    "witness_row_ids",
                },
                f"{row_id} artifact class",
            )
            name = str(class_record.get("artifact_class", ""))
            if name in classes_by_name:
                raise ProjectionError(f"{row_id} duplicates artifact class {name}")
            classes_by_name[name] = class_record
        if set(classes_by_name) != set(ARTIFACT_CLASS_CONTRACT):
            blockers.append("ARTIFACT_CLASS_MISSING")

        normalized_classes: list[dict[str, Any]] = []
        for name, contract in ARTIFACT_CLASS_CONTRACT.items():
            raw_class = classes_by_name.get(name)
            if raw_class is None:
                continue
            domain = raw_class.get("storage_domain")
            count_per_job = raw_class.get("count_per_job")
            basis = raw_class.get("basis")
            ceiling = require_positive_int(
                raw_class.get("bytes_per_artifact_ceiling"),
                f"{row_id}.{name}.bytes_per_artifact_ceiling",
            )
            witness_ids = require_sequence(
                raw_class.get("witness_row_ids"),
                f"{row_id}.{name}.witness_row_ids",
            )
            if (
                domain != contract["domain"]
                or count_per_job != contract["count_per_job"]
                or basis != contract["basis"]
                or not all(isinstance(value, str) for value in witness_ids)
            ):
                blockers.append("ARTIFACT_CLASS_MISSING")
            if contract["requires_witness"]:
                compatible = [
                    witnesses[witness_id]
                    for witness_id in witness_ids
                    if witness_id in witnesses
                    and witnesses[witness_id]["system"] == row["system"]
                ]
                if not compatible:
                    blockers.append("CAPACITY_WITNESS_NOT_STORAGE_CEILING")
                else:
                    witnessed_max = max(
                        record[name]["size_bytes"] for record in compatible
                    )
                    required_ceiling = ceil_ratio(
                        witnessed_max,
                        MIN_WITNESS_CEILING_NUMERATOR,
                        MIN_WITNESS_CEILING_DENOMINATOR,
                    )
                    if ceiling < required_ceiling:
                        blockers.append(
                            "CAPACITY_WITNESS_NOT_STORAGE_CEILING"
                        )
            elif witness_ids:
                blockers.append("ARTIFACT_CLASS_MISSING")

            job_count = int(row["expected_job_count"])
            artifact_count = job_count * int(contract["count_per_job"])
            projected_bytes = artifact_count * ceiling
            if domain in domain_totals:
                domain_totals[domain]["base_bytes"] += projected_bytes
                domain_totals[domain]["base_inodes"] += artifact_count
            normalized_classes.append(
                {
                    "artifact_class": name,
                    "storage_domain": domain,
                    "count_per_job": count_per_job,
                    "bytes_per_artifact_ceiling": ceiling,
                    "basis": basis,
                    "witness_row_ids": list(witness_ids),
                    "projected_artifact_count": artifact_count,
                    "projected_bytes": projected_bytes,
                    "projected_inodes": artifact_count,
                }
            )
        normalized.append(
            {
                **dict(row),
                "artifact_classes": normalized_classes,
            }
        )
    return normalized, domain_totals, ordered_blockers(blockers)


def build_storage_manifest_from_validated(
    spec: Mapping[str, Any],
    *,
    plan: Mapping[str, Any],
    rows: Sequence[Mapping[str, Any]],
    witnesses: Mapping[str, Mapping[str, Any]],
    bindings: Mapping[str, Any],
    chain_blockers: Sequence[str] = (),
) -> dict[str, Any]:
    blockers = list(chain_blockers)
    controller, controller_blockers = validate_controller_budget(
        spec.get("controller_dry_materialization_budget")
    )
    blockers.extend(controller_blockers)
    controller_derivation, derivation_blockers = validate_controller_derivation(
        spec.get("controller_dry_materialization_derivation"),
        controller=controller,
        plan=plan,
        bindings=bindings,
    )
    blockers.extend(derivation_blockers)
    retry = require_mapping(spec.get("retry_reserve"), "retry_reserve")
    require_exact_keys(retry, {"numerator", "denominator"}, "retry_reserve")
    retry_numerator = require_positive_int(
        retry.get("numerator"), "retry_reserve.numerator"
    )
    retry_denominator = require_positive_int(
        retry.get("denominator"), "retry_reserve.denominator"
    )
    if (
        retry_numerator >= retry_denominator
        or retry_numerator * MIN_RETRY_RESERVE_DENOMINATOR
        < retry_denominator * MIN_RETRY_RESERVE_NUMERATOR
    ):
        raise ProjectionError(
            "retry reserve must be at least 10 percent and below 100 percent"
        )

    normalized_rows, domain_totals, envelope_blockers = normalize_row_envelopes(
        spec.get("row_envelopes"),
        rows=rows,
        witnesses=witnesses,
    )
    blockers.extend(envelope_blockers)
    for domain, totals in domain_totals.items():
        totals["retry_reserve_bytes"] = ceil_ratio(
            totals["base_bytes"], retry_numerator, retry_denominator
        )
        totals["retry_reserve_inodes"] = ceil_ratio(
            totals["base_inodes"], retry_numerator, retry_denominator
        )
        totals["controller_bytes"] = 0
        totals["controller_inodes"] = 0
        totals["projected_increment_bytes"] = (
            totals["base_bytes"] + totals["retry_reserve_bytes"]
        )
        totals["projected_increment_inodes"] = (
            totals["base_inodes"] + totals["retry_reserve_inodes"]
        )
    if controller is not None:
        control = domain_totals["scratch_control"]
        control["controller_bytes"] = controller["projected_bytes"]
        control["controller_inodes"] = controller["projected_inodes"]
        control["projected_increment_bytes"] += controller["projected_bytes"]
        control["projected_increment_inodes"] += controller["projected_inodes"]

    blocker_codes = ordered_blockers(blockers)
    status = PASS_STATUS if not blocker_codes else BLOCKED_STATUS
    campaign = require_mapping(plan.get("campaign"), "plan.campaign")
    plan_roots = {
        "bulk_science": require_absolute_remote_root(
            campaign.get("output_root"), "campaign.output_root"
        ),
        "scratch_control": require_absolute_remote_root(
            campaign.get("submit_root"), "campaign.submit_root"
        ),
        "scratch_evidence": require_absolute_remote_root(
            campaign.get("evidence_root"), "campaign.evidence_root"
        ),
        "scheduler_streams": require_absolute_remote_root(
            campaign.get("submit_root"), "campaign.submit_root"
        ),
    }
    measurement = apply_authority(
        {
            "schema": MEASUREMENT_SCHEMA,
            "status": status,
            "gate": GATE,
            "bindings": dict(bindings),
            "count_contract": {
                "row_count": EXPECTED_ROW_COUNT,
                "source_tuple_count": EXPECTED_SOURCE_TUPLES,
                "group_size": EXPECTED_GROUP_SIZE,
                "expected_job_count": EXPECTED_JOBS,
                "expected_output_pair_count": EXPECTED_OUTPUT_PAIRS,
                "expected_analysis_output_count": EXPECTED_ANALYSIS_OUTPUTS,
                "expected_sidecar_output_count": EXPECTED_SIDECAR_OUTPUTS,
                "expected_physical_root_artifact_count": EXPECTED_ROOT_ARTIFACTS,
                "expected_retained_analysis_output_count": (
                    EXPECTED_RETAINED_ANALYSIS_OUTPUTS
                ),
                "expected_durable_root_artifact_count": (
                    EXPECTED_DURABLE_ROOT_ARTIFACTS
                ),
                "expected_source_occurrence_count": EXPECTED_SOURCE_OCCURRENCES,
                "request_memory_mb": EXPECTED_REQUEST_MEMORY_MB,
            },
            "capacity_witnesses": [
                dict(witnesses[row_id]) for row_id in sorted(witnesses)
            ],
            "row_envelopes": normalized_rows,
            "controller_dry_materialization_budget": controller,
            "controller_dry_materialization_derivation": controller_derivation,
            "retry_reserve": {
                "numerator": retry_numerator,
                "denominator": retry_denominator,
            },
            "storage_domain_roots": plan_roots,
            "storage_domain_totals": domain_totals,
            "blocker_codes": blocker_codes,
        }
    )
    measurement["measurement_semantic_sha256"] = semantic_sha256(measurement)
    manifest = apply_authority(
        {
            "schema": STORAGE_MANIFEST_SCHEMA,
            "status": status,
            "gate": GATE,
            "artifact_measurement": measurement,
            "blocker_codes": blocker_codes,
        }
    )
    manifest["manifest_semantic_sha256"] = semantic_sha256(manifest)
    return manifest


def build_storage_manifest(spec_path: Path) -> dict[str, Any]:
    try:
        spec = json.loads(spec_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ProjectionError(f"cannot read measurement spec: {spec_path}") from exc
    spec = require_mapping(spec, "measurement spec")
    require_exact_keys(
        spec,
        {
            "schema",
            "authority",
            *AUTHORITY_FIELDS,
            "bindings",
            "controller_dry_materialization_budget",
            "controller_dry_materialization_derivation",
            "retry_reserve",
            "row_envelopes",
        },
        "measurement spec",
    )
    if spec.get("schema") != MEASUREMENT_SPEC_SCHEMA:
        raise ProjectionError("measurement spec schema differs")
    validate_authority(spec, "measurement spec")
    validate_all_authority_occurrences(spec, "measurement spec")
    plan, rows, witnesses, bindings, blockers = validate_chain(spec)
    if contains_p5a_s_input(spec):
        blockers.append("P5A_S_INPUT_FORBIDDEN")
    return build_storage_manifest_from_validated(
        spec,
        plan=plan,
        rows=rows,
        witnesses=witnesses,
        bindings=bindings,
        chain_blockers=blockers,
    )


def intended_storage_roots(plan: Mapping[str, Any]) -> dict[str, str]:
    campaign = require_mapping(plan.get("campaign"), "plan.campaign")
    return {
        "bulk_science": require_absolute_remote_root(
            campaign.get("output_root"), "campaign.output_root"
        ),
        "scratch_control": require_absolute_remote_root(
            campaign.get("submit_root"), "campaign.submit_root"
        ),
        "scratch_evidence": require_absolute_remote_root(
            campaign.get("evidence_root"), "campaign.evidence_root"
        ),
        "scheduler_streams": require_absolute_remote_root(
            campaign.get("submit_root"), "campaign.submit_root"
        ),
    }


def nearest_existing_parent(path: Path) -> Path:
    candidate = path
    while not candidate.exists():
        if candidate.parent == candidate:
            raise ProjectionError(f"no existing parent for quota path: {path}")
        candidate = candidate.parent
    return candidate


def parse_human_size(value: str, label: str) -> int:
    token = value.strip().rstrip("*")
    token = token.strip("[]")
    if token in {"-", "none", "None"}:
        return 0
    match = SIZE_RE.fullmatch(token)
    if match is None:
        raise ProjectionError(f"{label} has unsupported size token {value!r}")
    number_text, suffix = match.groups()
    exponent = {"": 0, "K": 1, "M": 2, "G": 3, "T": 4, "P": 5, "E": 6}[
        suffix.upper()
    ]
    # Lustre's default unsuffixed block unit is KiB.  Human-readable output
    # normally carries a suffix; preserving the default here is fail-safe.
    multiplier = 1024 ** (exponent if suffix else 1)
    if "." in number_text:
        whole, fraction = number_text.split(".", 1)
        scale = 10 ** len(fraction)
        numerator = int(whole) * scale + int(fraction)
        return (numerator * multiplier + scale - 1) // scale
    return int(number_text) * multiplier


def parse_lfs_quota_output(
    stdout: str,
    stderr: str,
    returncode: int,
) -> dict[str, Any]:
    inaccurate = (
        "[" in stdout
        or "]" in stdout
        or "inaccurate" in stdout.lower()
        or "inaccurate" in stderr.lower()
    )
    lines = [line.strip() for line in stdout.splitlines() if line.strip()]
    tokens: list[str] | None = None
    for index, line in enumerate(lines):
        if not line.startswith("/"):
            continue
        candidate = line.split()
        if len(candidate) < 8 and index + 1 < len(lines):
            candidate.extend(lines[index + 1].split())
        if len(candidate) >= 8:
            tokens = candidate
            break
    if tokens is None:
        return {
            "usage_authoritative": False,
            "blocker_codes": ["QUOTA_USAGE_UNAUTHORITATIVE"],
        }
    filesystem = tokens[0]
    try:
        used_bytes = parse_human_size(tokens[1], "quota used bytes")
        soft_bytes = parse_human_size(tokens[2], "quota soft bytes")
        hard_bytes = parse_human_size(tokens[3], "quota hard bytes")
        used_inodes = int(tokens[5].rstrip("*[]"))
        soft_inodes = int(tokens[6].rstrip("*[]"))
        hard_inodes = int(tokens[7].rstrip("*[]"))
    except (ProjectionError, ValueError):
        return {
            "usage_authoritative": False,
            "blocker_codes": ["QUOTA_USAGE_UNAUTHORITATIVE"],
        }
    byte_limits = [value for value in (soft_bytes, hard_bytes) if value > 0]
    inode_limits = [value for value in (soft_inodes, hard_inodes) if value > 0]
    authoritative = (
        returncode == 0
        and not inaccurate
        and bool(byte_limits)
        and bool(inode_limits)
        and used_bytes >= 0
        and used_inodes >= 0
    )
    return {
        "filesystem": filesystem,
        "used_bytes": used_bytes,
        "quota_bytes": min(byte_limits) if byte_limits else None,
        "used_inodes": used_inodes,
        "quota_inodes": min(inode_limits) if inode_limits else None,
        "usage_authoritative": authoritative,
        "blocker_codes": [] if authoritative else ["QUOTA_USAGE_UNAUTHORITATIVE"],
    }


def build_quota_snapshot(
    plan_path: Path,
    principal: str,
    *,
    max_age_seconds: int = DEFAULT_SNAPSHOT_MAX_AGE_SECONDS,
    runner: Callable[..., subprocess.CompletedProcess[str]] = subprocess.run,
    now_seconds: int | None = None,
    lfs_binary: str | None = None,
) -> dict[str, Any]:
    if PRINCIPAL_RE.fullmatch(principal) is None:
        raise ProjectionError("principal contains unsupported characters")
    max_age_seconds = require_positive_int(max_age_seconds, "max_age_seconds")
    try:
        plan = json.loads(plan_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ProjectionError(f"cannot read extraction plan: {plan_path}") from exc
    plan = require_mapping(plan, "extraction plan")
    rows, count_blockers = validate_plan_counts(plan)
    del rows
    plan_sha256 = file_sha256(plan_path)
    roots = intended_storage_roots(plan)
    binary = lfs_binary or shutil.which("lfs")
    observed_seconds = int(time.time()) if now_seconds is None else now_seconds
    require_nonnegative_int(observed_seconds, "observed time")
    domains: list[dict[str, Any]] = []
    blockers = list(count_blockers)
    for domain in STORAGE_DOMAINS:
        root = roots[domain]
        namespace_fresh = not os.path.lexists(root)
        if not namespace_fresh:
            blockers.append("NAMESPACE_NOT_FRESH")
        probe_path = nearest_existing_parent(Path(root))
        argv = [
            binary or "lfs",
            "quota",
            "-h",
            "-u",
            principal,
            str(probe_path),
        ]
        if binary is None:
            result = subprocess.CompletedProcess(argv, 127, "", "lfs not found")
        else:
            try:
                result = runner(
                    argv,
                    text=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    check=False,
                    timeout=30,
                )
            except (OSError, subprocess.SubprocessError) as exc:
                result = subprocess.CompletedProcess(argv, 126, "", str(exc))
        parsed = parse_lfs_quota_output(
            result.stdout or "", result.stderr or "", result.returncode
        )
        blockers.extend(parsed.get("blocker_codes", []))
        filesystem = parsed.get("filesystem")
        quota_domain_id = (
            semantic_sha256(
                {
                    "principal": principal,
                    "filesystem": filesystem,
                }
            )
            if isinstance(filesystem, str) and filesystem
            else None
        )
        domains.append(
            {
                "storage_domain": domain,
                "intended_root": root,
                "namespace_fresh": namespace_fresh,
                "quota_probe_path": str(probe_path),
                "quota_domain_id": quota_domain_id,
                "filesystem": filesystem,
                "query_argv": argv,
                "query_exit_code": result.returncode,
                "stdout": result.stdout or "",
                "stderr": result.stderr or "",
                "stdout_sha256": hashlib.sha256(
                    (result.stdout or "").encode("utf-8")
                ).hexdigest(),
                "stderr_sha256": hashlib.sha256(
                    (result.stderr or "").encode("utf-8")
                ).hexdigest(),
                "used_bytes": parsed.get("used_bytes"),
                "quota_bytes": parsed.get("quota_bytes"),
                "used_inodes": parsed.get("used_inodes"),
                "quota_inodes": parsed.get("quota_inodes"),
                "usage_authoritative": parsed.get(
                    "usage_authoritative", False
                ),
                "blocker_codes": ordered_blockers(
                    parsed.get("blocker_codes", [])
                    + ([] if namespace_fresh else ["NAMESPACE_NOT_FRESH"])
                ),
            }
        )
    blocker_codes = ordered_blockers(blockers)
    status = PASS_STATUS if not blocker_codes else BLOCKED_STATUS
    payload = apply_authority(
        {
            "schema": QUOTA_SNAPSHOT_SCHEMA,
            "status": status,
            "gate": GATE,
            "plan": {
                "path": str(plan_path),
                "sha256": plan_sha256,
            },
            "principal": principal,
            "observed_at": datetime.fromtimestamp(
                observed_seconds, tz=timezone.utc
            ).isoformat(),
            "observed_unix_seconds": observed_seconds,
            "max_age_seconds": max_age_seconds,
            "storage_domains": domains,
            "blocker_codes": blocker_codes,
        }
    )
    payload["snapshot_semantic_sha256"] = semantic_sha256(payload)
    return payload


def validate_materialized_artifact_record(
    value: object,
    label: str,
) -> dict[str, Any]:
    record = require_mapping(value, label)
    require_exact_keys(record, {"path", "sha256", "size_bytes"}, label)
    path_text = record.get("path")
    if not isinstance(path_text, str) or not path_text:
        raise ProjectionError(f"{label}.path must be nonempty")
    path = Path(path_text)
    if not path.is_absolute() or not path.is_file() or path.is_symlink():
        raise ProjectionError(
            f"{label}.path must be an existing absolute regular file"
        )
    expected_sha256 = require_sha256(record.get("sha256"), f"{label}.sha256")
    observed_sha256 = file_sha256(path)
    if observed_sha256 != expected_sha256:
        raise ProjectionError(f"{label} SHA-256 differs")
    expected_size = require_nonnegative_int(
        record.get("size_bytes"), f"{label}.size_bytes"
    )
    if path.stat().st_size != expected_size:
        raise ProjectionError(f"{label}.size_bytes differs")
    return {
        "path": str(path),
        "sha256": observed_sha256,
        "size_bytes": expected_size,
    }


def validate_capacity_witnesses(
    value: object,
) -> dict[str, dict[str, Any]]:
    records = require_sequence(value, "measurement capacity witnesses")
    by_id: dict[str, dict[str, Any]] = {}
    expected_witness_ids = set(amendment_tool.SELECTED_CAPACITY_ROWS)
    expected_witness_keys = {
        "witness_id",
        "row_id",
        "system",
        "tuple_count",
        "request_memory_mb",
        "wall_time_seconds",
        "max_rss_mb",
        "storage_ceiling_authority",
        "analysis_root",
        "training_sidecar_root",
    }
    expected_sidecar_artifact_keys = {
        "path",
        "sha256",
        "size_bytes",
        "health_profile",
        "health_status",
    }
    expected_analysis_artifact_keys = {
        "path",
        "artifact_state",
        "health_receipt_path",
        "health_receipt_sha256",
        "key_inventory_sha256",
        "size_bytes",
        "health_profile",
        "health_status",
    }
    for raw in records:
        record = require_mapping(raw, "measurement capacity witness")
        require_exact_keys(
            record, expected_witness_keys, "measurement capacity witness"
        )
        witness_id = str(record.get("witness_id", ""))
        if (
            not witness_id
            or witness_id in by_id
            or record.get("row_id") != witness_id
            or record.get("system")
            != amendment_tool.EXPECTED_SYSTEM.get(witness_id)
        ):
            raise ProjectionError("measurement capacity witness identity differs")
        require_exact_int(
            record.get("tuple_count"),
            EXPECTED_GROUP_SIZE,
            f"{witness_id}.tuple_count",
        )
        require_exact_int(
            record.get("request_memory_mb"),
            EXPECTED_REQUEST_MEMORY_MB,
            f"{witness_id}.request_memory_mb",
        )
        require_positive_int(
            record.get("wall_time_seconds"),
            f"{witness_id}.wall_time_seconds",
        )
        require_positive_int(
            record.get("max_rss_mb"), f"{witness_id}.max_rss_mb"
        )
        if record.get("storage_ceiling_authority") is not False:
            raise ProjectionError(
                f"{witness_id}.storage_ceiling_authority must remain false"
            )
        normalized = dict(record)
        for artifact_class, health_profile in (
            ("analysis_root", ANALYSIS_HEALTH_PROFILE),
            ("training_sidecar_root", SIDECAR_HEALTH_PROFILE),
        ):
            artifact = require_mapping(
                record.get(artifact_class),
                f"{witness_id}.{artifact_class}",
            )
            current_ephemeral_analysis = (
                artifact_class == "analysis_root"
                and "artifact_state" in artifact
            )
            require_exact_keys(
                artifact,
                (
                    expected_analysis_artifact_keys
                    if current_ephemeral_analysis
                    else expected_sidecar_artifact_keys
                ),
                f"{witness_id}.{artifact_class}",
            )
            path_text = artifact.get("path")
            if not isinstance(path_text, str) or not path_text:
                raise ProjectionError(
                    f"{witness_id}.{artifact_class}.path must be nonempty"
                )
            if current_ephemeral_analysis:
                if (
                    artifact.get("artifact_state")
                    != "EPHEMERAL_VALIDATED_NOT_RETAINED"
                    or not isinstance(
                        artifact.get("health_receipt_path"), str
                    )
                    or not artifact.get("health_receipt_path")
                ):
                    raise ProjectionError(
                        f"{witness_id}.analysis_root receipt binding differs"
                    )
                require_sha256(
                    artifact.get("health_receipt_sha256"),
                    f"{witness_id}.analysis_root.health_receipt_sha256",
                )
                require_sha256(
                    artifact.get("key_inventory_sha256"),
                    f"{witness_id}.analysis_root.key_inventory_sha256",
                )
            else:
                require_sha256(
                    artifact.get("sha256"),
                    f"{witness_id}.{artifact_class}.sha256",
                )
            require_positive_int(
                artifact.get("size_bytes"),
                f"{witness_id}.{artifact_class}.size_bytes",
            )
            if (
                artifact.get("health_profile") != health_profile
                or artifact.get("health_status") != "PASS"
            ):
                raise ProjectionError(
                    f"{witness_id}.{artifact_class} health binding differs"
                )
        by_id[witness_id] = normalized
    if set(by_id) != expected_witness_ids:
        raise ProjectionError("measurement capacity witness inventory differs")
    return by_id


def validate_measurement_derivation(
    measurement: Mapping[str, Any],
) -> tuple[dict[str, str], dict[str, dict[str, int]]]:
    bindings = require_mapping(
        measurement.get("bindings"), "measurement bindings"
    )
    require_exact_keys(
        bindings,
        {
            "plan",
            "preflight_receipt",
            "capacity_count_amendment",
            "capacity_count_amendment_readback",
            "bundle_manifest",
            "materialization_receipt",
            "execution_partition_sha256",
            "partition_artifact_sha256",
        },
        "measurement bindings",
    )
    validated_binding_records: dict[str, dict[str, Any]] = {}
    for field in (
        "plan",
        "preflight_receipt",
        "capacity_count_amendment",
        "capacity_count_amendment_readback",
        "bundle_manifest",
        "materialization_receipt",
    ):
        validated_binding_records[field] = validate_materialized_artifact_record(
            bindings.get(field), f"measurement bindings.{field}"
        )
    for field in (
        "execution_partition_sha256",
        "partition_artifact_sha256",
    ):
        require_sha256(
            bindings.get(field), f"measurement bindings.{field}"
        )

    (
        _amendment,
        revalidated_amendment_artifact,
        revalidated_readback_artifact,
        expected_witnesses,
    ) = load_revalidated_amendment_chain(
        {
            "path": validated_binding_records["capacity_count_amendment"]["path"],
            "sha256": validated_binding_records[
                "capacity_count_amendment"
            ]["sha256"],
        },
        {
            "path": validated_binding_records[
                "capacity_count_amendment_readback"
            ]["path"],
            "sha256": validated_binding_records[
                "capacity_count_amendment_readback"
            ]["sha256"],
        },
    )
    if (
        revalidated_amendment_artifact
        != validated_binding_records["capacity_count_amendment"]
        or revalidated_readback_artifact
        != validated_binding_records["capacity_count_amendment_readback"]
    ):
        raise ProjectionError(
            "measurement capacity amendment artifact binding differs"
        )

    raw_witnesses = require_sequence(
        measurement.get("capacity_witnesses"),
        "measurement capacity witnesses",
    )
    expected_witness_order = sorted(expected_witnesses)
    if [
        record.get("witness_id") if isinstance(record, dict) else None
        for record in raw_witnesses
    ] != expected_witness_order:
        raise ProjectionError("measurement capacity witness ordering differs")
    witnesses = validate_capacity_witnesses(
        measurement.get("capacity_witnesses")
    )
    if witnesses != expected_witnesses:
        raise ProjectionError(
            "measurement capacity witnesses differ from the bound amendment"
        )
    rows = require_sequence(
        measurement.get("row_envelopes"), "measurement row envelopes"
    )
    expected_ids = expected_row_ids()
    observed_ids = [
        row.get("row_id") for row in rows if isinstance(row, dict)
    ]
    if observed_ids != expected_ids or len(rows) != EXPECTED_ROW_COUNT:
        raise ProjectionError("measurement row envelope inventory differs")

    derived_totals = {
        domain: {"base_bytes": 0, "base_inodes": 0}
        for domain in STORAGE_DOMAINS
    }
    derived_artifact_counts = {
        name: 0 for name in ARTIFACT_CLASS_CONTRACT
    }
    source_tuple_total = 0
    job_total = 0
    expected_row_keys = {
        "row_id",
        "system",
        "source_tuple_count",
        "expected_job_count",
        "artifact_classes",
    }
    expected_class_keys = {
        "artifact_class",
        "storage_domain",
        "count_per_job",
        "bytes_per_artifact_ceiling",
        "basis",
        "witness_row_ids",
        "projected_artifact_count",
        "projected_bytes",
        "projected_inodes",
    }
    for raw_row in rows:
        row = require_mapping(raw_row, "measurement row envelope")
        require_exact_keys(row, expected_row_keys, "measurement row envelope")
        row_id = str(row.get("row_id", ""))
        system = row.get("system")
        if system not in {"pp", "auau"}:
            raise ProjectionError(f"{row_id}.system differs")
        source_tuples = require_positive_int(
            row.get("source_tuple_count"), f"{row_id}.source_tuple_count"
        )
        jobs = require_positive_int(
            row.get("expected_job_count"), f"{row_id}.expected_job_count"
        )
        if jobs != ceil_ratio(source_tuples, 1, EXPECTED_GROUP_SIZE):
            raise ProjectionError(f"{row_id}.expected_job_count differs")
        source_tuple_total += source_tuples
        job_total += jobs

        class_records = require_sequence(
            row.get("artifact_classes"), f"{row_id}.artifact_classes"
        )
        if [
            record.get("artifact_class")
            for record in class_records
            if isinstance(record, dict)
        ] != list(ARTIFACT_CLASS_CONTRACT):
            raise ProjectionError(f"{row_id} artifact class inventory differs")
        for raw_class in class_records:
            record = require_mapping(
                raw_class, f"{row_id} artifact class"
            )
            require_exact_keys(
                record, expected_class_keys, f"{row_id} artifact class"
            )
            name = str(record.get("artifact_class", ""))
            contract = ARTIFACT_CLASS_CONTRACT.get(name)
            if contract is None:
                raise ProjectionError(
                    f"{row_id} has unknown artifact class {name!r}"
                )
            if (
                record.get("storage_domain") != contract["domain"]
                or record.get("count_per_job") != contract["count_per_job"]
                or record.get("basis") != contract["basis"]
            ):
                raise ProjectionError(
                    f"{row_id}.{name} artifact contract differs"
                )
            ceiling = require_positive_int(
                record.get("bytes_per_artifact_ceiling"),
                f"{row_id}.{name}.bytes_per_artifact_ceiling",
            )
            witness_ids = require_sequence(
                record.get("witness_row_ids"),
                f"{row_id}.{name}.witness_row_ids",
            )
            if not all(isinstance(value, str) for value in witness_ids):
                raise ProjectionError(f"{row_id}.{name} witness IDs differ")
            if contract["requires_witness"]:
                compatible = [
                    witnesses[witness_id]
                    for witness_id in witness_ids
                    if witness_id in witnesses
                    and witnesses[witness_id]["system"] == system
                ]
                if not compatible:
                    raise ProjectionError(
                        f"{row_id}.{name} lacks a compatible capacity witness"
                    )
                witnessed_max = max(
                    witness[name]["size_bytes"] for witness in compatible
                )
                minimum_ceiling = ceil_ratio(
                    witnessed_max,
                    MIN_WITNESS_CEILING_NUMERATOR,
                    MIN_WITNESS_CEILING_DENOMINATOR,
                )
                if ceiling < minimum_ceiling:
                    raise ProjectionError(
                        f"{row_id}.{name} storage ceiling is below its witness bound"
                    )
            elif witness_ids:
                raise ProjectionError(
                    f"{row_id}.{name} must not bind capacity witnesses"
                )

            artifact_count = jobs * int(contract["count_per_job"])
            projected_bytes = artifact_count * ceiling
            if (
                record.get("projected_artifact_count") != artifact_count
                or record.get("projected_inodes") != artifact_count
                or record.get("projected_bytes") != projected_bytes
            ):
                raise ProjectionError(
                    f"{row_id}.{name} projected artifact arithmetic differs"
                )
            domain = str(contract["domain"])
            derived_artifact_counts[name] += artifact_count
            derived_totals[domain]["base_bytes"] += projected_bytes
            derived_totals[domain]["base_inodes"] += artifact_count

    if (
        source_tuple_total != EXPECTED_SOURCE_TUPLES
        or job_total != EXPECTED_JOBS
    ):
        raise ProjectionError("measurement row aggregate counts differ")
    if (
        derived_artifact_counts["analysis_root"]
        != EXPECTED_RETAINED_ANALYSIS_OUTPUTS
        or derived_artifact_counts["training_sidecar_root"]
        != EXPECTED_SIDECAR_OUTPUTS
        or (
            derived_artifact_counts["analysis_root"]
            + derived_artifact_counts["training_sidecar_root"]
        )
        != EXPECTED_DURABLE_ROOT_ARTIFACTS
        or any(
            derived_artifact_counts[name] != EXPECTED_JOBS
            for name in (
                "evidence_record",
                "scheduler_stdout",
                "scheduler_stderr",
                "scheduler_event_log",
            )
        )
    ):
        raise ProjectionError("measurement artifact aggregate counts differ")

    retry = require_mapping(
        measurement.get("retry_reserve"), "measurement retry reserve"
    )
    require_exact_keys(
        retry, {"numerator", "denominator"}, "measurement retry reserve"
    )
    retry_numerator = require_positive_int(
        retry.get("numerator"), "measurement retry reserve.numerator"
    )
    retry_denominator = require_positive_int(
        retry.get("denominator"), "measurement retry reserve.denominator"
    )
    if (
        retry_numerator >= retry_denominator
        or retry_numerator * MIN_RETRY_RESERVE_DENOMINATOR
        < retry_denominator * MIN_RETRY_RESERVE_NUMERATOR
    ):
        raise ProjectionError("measurement retry reserve differs")

    controller = require_mapping(
        measurement.get("controller_dry_materialization_budget"),
        "measurement controller budget",
    )
    require_exact_keys(
        controller,
        {
            "receipt",
            "expected_job_count",
            "storage_domain",
            "fixed_bytes",
            "fixed_inodes",
            "bytes_per_job",
            "inodes_per_job",
            "projected_bytes",
            "projected_inodes",
        },
        "measurement controller budget",
    )
    validate_materialized_artifact_record(
        controller.get("receipt"), "measurement controller budget receipt"
    )
    if (
        controller.get("expected_job_count") != EXPECTED_JOBS
        or controller.get("storage_domain") != "scratch_control"
    ):
        raise ProjectionError("measurement controller count/domain differs")
    controller_values = {
        field: require_nonnegative_int(
            controller.get(field), f"measurement controller budget.{field}"
        )
        for field in (
            "fixed_bytes",
            "fixed_inodes",
            "bytes_per_job",
            "inodes_per_job",
        )
    }
    if (
        controller_values["bytes_per_job"] == 0
        or controller_values["inodes_per_job"] == 0
    ):
        raise ProjectionError("measurement controller per-job budget differs")
    controller_bytes = (
        controller_values["fixed_bytes"]
        + EXPECTED_JOBS * controller_values["bytes_per_job"]
    )
    controller_inodes = (
        controller_values["fixed_inodes"]
        + EXPECTED_JOBS * controller_values["inodes_per_job"]
    )
    if (
        controller.get("projected_bytes") != controller_bytes
        or controller.get("projected_inodes") != controller_inodes
    ):
        raise ProjectionError("measurement controller arithmetic differs")
    plan_payload, _ = load_artifact(
        {
            "path": validated_binding_records["plan"]["path"],
            "sha256": validated_binding_records["plan"]["sha256"],
        },
        "measurement derivation plan",
        expected_schema=resolver.PLAN_SCHEMA,
    )
    derivation_record = validate_materialized_artifact_record(
        measurement.get("controller_dry_materialization_derivation"),
        "measurement controller derivation receipt",
    )
    revalidated_derivation, derivation_blockers = validate_controller_derivation(
        {
            "path": derivation_record["path"],
            "sha256": derivation_record["sha256"],
        },
        controller=controller,
        plan=plan_payload,
        bindings=bindings,
    )
    if derivation_blockers or revalidated_derivation != derivation_record:
        raise ProjectionError("measurement controller derivation differs")

    for domain, totals in derived_totals.items():
        totals["retry_reserve_bytes"] = ceil_ratio(
            totals["base_bytes"], retry_numerator, retry_denominator
        )
        totals["retry_reserve_inodes"] = ceil_ratio(
            totals["base_inodes"], retry_numerator, retry_denominator
        )
        totals["controller_bytes"] = (
            controller_bytes if domain == "scratch_control" else 0
        )
        totals["controller_inodes"] = (
            controller_inodes if domain == "scratch_control" else 0
        )
        totals["projected_increment_bytes"] = (
            totals["base_bytes"]
            + totals["retry_reserve_bytes"]
            + totals["controller_bytes"]
        )
        totals["projected_increment_inodes"] = (
            totals["base_inodes"]
            + totals["retry_reserve_inodes"]
            + totals["controller_inodes"]
        )

    observed_totals = require_mapping(
        measurement.get("storage_domain_totals"),
        "measurement storage-domain totals",
    )
    if observed_totals != derived_totals:
        raise ProjectionError("measurement storage-domain totals differ")

    roots_raw = require_mapping(
        measurement.get("storage_domain_roots"),
        "measurement storage-domain roots",
    )
    if set(roots_raw) != set(STORAGE_DOMAINS):
        raise ProjectionError("measurement storage-domain root inventory differs")
    roots = {
        domain: require_absolute_remote_root(
            roots_raw.get(domain), f"measurement root {domain}"
        )
        for domain in STORAGE_DOMAINS
    }
    return roots, derived_totals


def validate_storage_manifest(payload: Mapping[str, Any]) -> dict[str, Any]:
    if payload.get("schema") != STORAGE_MANIFEST_SCHEMA:
        raise ProjectionError("storage manifest schema differs")
    require_exact_keys(
        payload,
        {
            "schema",
            "status",
            "gate",
            "artifact_measurement",
            "blocker_codes",
            "authority",
            *AUTHORITY_FIELDS,
            "manifest_semantic_sha256",
        },
        "storage manifest",
    )
    validate_authority(payload, "storage manifest")
    validate_all_authority_occurrences(payload, "storage manifest")
    verify_semantic_receipt(
        payload, "manifest_semantic_sha256", "storage manifest"
    )
    if payload.get("gate") != GATE:
        raise ProjectionError("storage manifest gate differs")
    blocker_codes = require_sequence(
        payload.get("blocker_codes"), "storage manifest blocker codes"
    )
    if blocker_codes != ordered_blockers(blocker_codes):
        raise ProjectionError("storage manifest blocker ordering differs")
    expected_status = PASS_STATUS if not blocker_codes else BLOCKED_STATUS
    if payload.get("status") != expected_status:
        raise ProjectionError("storage manifest status/blockers differ")
    measurement = require_mapping(
        payload.get("artifact_measurement"), "artifact measurement"
    )
    if measurement.get("schema") != MEASUREMENT_SCHEMA:
        raise ProjectionError("artifact measurement schema differs")
    require_exact_keys(
        measurement,
        {
            "schema",
            "status",
            "gate",
            "bindings",
            "count_contract",
            "capacity_witnesses",
            "row_envelopes",
            "controller_dry_materialization_budget",
            "controller_dry_materialization_derivation",
            "retry_reserve",
            "storage_domain_roots",
            "storage_domain_totals",
            "blocker_codes",
            "authority",
            *AUTHORITY_FIELDS,
            "measurement_semantic_sha256",
        },
        "artifact measurement",
    )
    validate_authority(measurement, "artifact measurement")
    verify_semantic_receipt(
        measurement,
        "measurement_semantic_sha256",
        "artifact measurement",
    )
    if (
        measurement.get("gate") != GATE
        or measurement.get("status") != expected_status
        or measurement.get("blocker_codes") != blocker_codes
    ):
        raise ProjectionError("artifact measurement status/gate differs")
    count_contract = require_mapping(
        measurement.get("count_contract"), "measurement count contract"
    )
    exact_counts = {
        "row_count": EXPECTED_ROW_COUNT,
        "source_tuple_count": EXPECTED_SOURCE_TUPLES,
        "group_size": EXPECTED_GROUP_SIZE,
        "expected_job_count": EXPECTED_JOBS,
        "expected_output_pair_count": EXPECTED_OUTPUT_PAIRS,
        "expected_analysis_output_count": EXPECTED_ANALYSIS_OUTPUTS,
        "expected_sidecar_output_count": EXPECTED_SIDECAR_OUTPUTS,
        "expected_physical_root_artifact_count": EXPECTED_ROOT_ARTIFACTS,
        "expected_retained_analysis_output_count": (
            EXPECTED_RETAINED_ANALYSIS_OUTPUTS
        ),
        "expected_durable_root_artifact_count": (
            EXPECTED_DURABLE_ROOT_ARTIFACTS
        ),
        "expected_source_occurrence_count": EXPECTED_SOURCE_OCCURRENCES,
        "request_memory_mb": EXPECTED_REQUEST_MEMORY_MB,
    }
    if count_contract != exact_counts:
        raise ProjectionError("measurement count contract differs")
    validate_measurement_derivation(measurement)
    if contains_p5a_s_input(payload):
        raise ProjectionError("P5A-S input is forbidden")
    return dict(payload)


def validate_quota_snapshot(payload: Mapping[str, Any]) -> dict[str, Any]:
    if payload.get("schema") != QUOTA_SNAPSHOT_SCHEMA:
        raise ProjectionError("quota snapshot schema differs")
    require_exact_keys(
        payload,
        {
            "schema",
            "status",
            "gate",
            "plan",
            "principal",
            "observed_at",
            "observed_unix_seconds",
            "max_age_seconds",
            "storage_domains",
            "blocker_codes",
            "authority",
            *AUTHORITY_FIELDS,
            "snapshot_semantic_sha256",
        },
        "quota snapshot",
    )
    validate_authority(payload, "quota snapshot")
    validate_all_authority_occurrences(payload, "quota snapshot")
    verify_semantic_receipt(
        payload, "snapshot_semantic_sha256", "quota snapshot"
    )
    if payload.get("gate") != GATE:
        raise ProjectionError("quota snapshot gate differs")
    principal = payload.get("principal")
    if not isinstance(principal, str) or PRINCIPAL_RE.fullmatch(principal) is None:
        raise ProjectionError("quota snapshot principal differs")
    blocker_codes = require_sequence(
        payload.get("blocker_codes"), "quota snapshot blocker codes"
    )
    if blocker_codes != ordered_blockers(blocker_codes):
        raise ProjectionError("quota snapshot blocker ordering differs")
    expected_status = PASS_STATUS if not blocker_codes else BLOCKED_STATUS
    if payload.get("status") != expected_status:
        raise ProjectionError("quota snapshot status/blockers differ")
    plan = artifact_ref(payload.get("plan"), "quota snapshot plan")
    del plan
    domains = require_sequence(
        payload.get("storage_domains"), "quota snapshot storage domains"
    )
    by_name: dict[str, dict[str, Any]] = {}
    expected_domain_keys = {
        "storage_domain",
        "intended_root",
        "namespace_fresh",
        "quota_probe_path",
        "quota_domain_id",
        "filesystem",
        "query_argv",
        "query_exit_code",
        "stdout",
        "stderr",
        "stdout_sha256",
        "stderr_sha256",
        "used_bytes",
        "quota_bytes",
        "used_inodes",
        "quota_inodes",
        "usage_authoritative",
        "blocker_codes",
    }
    for raw in domains:
        record = require_mapping(raw, "quota snapshot domain")
        require_exact_keys(record, expected_domain_keys, "quota snapshot domain")
        name = str(record.get("storage_domain", ""))
        if name in by_name:
            raise ProjectionError(f"quota snapshot duplicates domain {name}")
        by_name[name] = record
        require_absolute_remote_root(
            record.get("intended_root"), f"{name}.intended_root"
        )
        require_absolute_remote_root(
            record.get("quota_probe_path"), f"{name}.quota_probe_path"
        )
        require_sha256(record.get("stdout_sha256"), f"{name}.stdout_sha256")
        require_sha256(record.get("stderr_sha256"), f"{name}.stderr_sha256")
        query_argv = require_sequence(record.get("query_argv"), f"{name}.query_argv")
        if not query_argv or not all(
            isinstance(value, str) and value for value in query_argv
        ):
            raise ProjectionError(f"{name}.query_argv differs")
        intended_root = require_absolute_remote_root(
            record.get("intended_root"), f"{name}.intended_root"
        )
        quota_probe_path = require_absolute_remote_root(
            record.get("quota_probe_path"), f"{name}.quota_probe_path"
        )
        if (
            len(query_argv) != 6
            or Path(query_argv[0]).name != "lfs"
            or query_argv[1:]
            != [
                "quota",
                "-h",
                "-u",
                principal,
                quota_probe_path,
            ]
        ):
            raise ProjectionError(f"{name}.query_argv is not canonical")
        expected_probe = nearest_existing_parent(Path(intended_root))
        if expected_probe != Path(quota_probe_path):
            raise ProjectionError(
                f"{name}.quota_probe_path is not the nearest existing parent"
            )
        authoritative = record.get("usage_authoritative")
        namespace_fresh = record.get("namespace_fresh")
        if not isinstance(authoritative, bool) or not isinstance(
            namespace_fresh, bool
        ):
            raise ProjectionError(f"{name} quota booleans differ")
        if namespace_fresh != (not os.path.lexists(intended_root)):
            raise ProjectionError(f"{name}.namespace_fresh differs from readback")
        stdout = record.get("stdout")
        stderr = record.get("stderr")
        if not isinstance(stdout, str) or not isinstance(stderr, str):
            raise ProjectionError(f"{name} quota output must be text")
        if hashlib.sha256(stdout.encode("utf-8")).hexdigest() != record.get(
            "stdout_sha256"
        ):
            raise ProjectionError(f"{name}.stdout_sha256 differs")
        if hashlib.sha256(stderr.encode("utf-8")).hexdigest() != record.get(
            "stderr_sha256"
        ):
            raise ProjectionError(f"{name}.stderr_sha256 differs")
        query_exit_code = record.get("query_exit_code")
        if isinstance(query_exit_code, bool) or not isinstance(
            query_exit_code, int
        ):
            raise ProjectionError(f"{name}.query_exit_code differs")
        reparsed = parse_lfs_quota_output(
            stdout, stderr, query_exit_code
        )
        for field in (
            "filesystem",
            "used_bytes",
            "quota_bytes",
            "used_inodes",
            "quota_inodes",
            "usage_authoritative",
        ):
            if record.get(field) != reparsed.get(field):
                raise ProjectionError(
                    f"{name}.{field} differs from raw quota output"
                )
        domain_blockers = require_sequence(
            record.get("blocker_codes"), f"{name}.blocker_codes"
        )
        if domain_blockers != ordered_blockers(domain_blockers):
            raise ProjectionError(f"{name} blocker ordering differs")
        expected_domain_blockers = ordered_blockers(
            list(reparsed.get("blocker_codes", []))
            + ([] if namespace_fresh else ["NAMESPACE_NOT_FRESH"])
        )
        if domain_blockers != expected_domain_blockers:
            raise ProjectionError(
                f"{name} blockers differ from raw quota output"
            )
        if authoritative:
            if record.get("query_exit_code") != 0:
                raise ProjectionError(
                    f"{name} claims authoritative nonzero quota query"
                )
            filesystem = record.get("filesystem")
            if not isinstance(filesystem, str) or not filesystem:
                raise ProjectionError(f"{name}.filesystem differs")
            expected_id = semantic_sha256(
                {"principal": principal, "filesystem": filesystem}
            )
            if record.get("quota_domain_id") != expected_id:
                raise ProjectionError(f"{name}.quota_domain_id differs")
            for field in ("used_bytes", "used_inodes"):
                require_nonnegative_int(record.get(field), f"{name}.{field}")
            for field in ("quota_bytes", "quota_inodes"):
                require_positive_int(record.get(field), f"{name}.{field}")
            if "QUOTA_USAGE_UNAUTHORITATIVE" in domain_blockers:
                raise ProjectionError(f"{name} authority/blockers differ")
        elif "QUOTA_USAGE_UNAUTHORITATIVE" not in domain_blockers:
            raise ProjectionError(f"{name} missing quota authority blocker")
        if namespace_fresh == (
            "NAMESPACE_NOT_FRESH" in domain_blockers
        ):
            raise ProjectionError(f"{name} namespace/blockers differ")
    if set(by_name) != set(STORAGE_DOMAINS):
        raise ProjectionError("quota snapshot storage-domain inventory differs")
    return dict(payload)


def build_projection(
    manifest_path: Path,
    snapshot_path: Path,
    *,
    now_seconds: int | None = None,
) -> dict[str, Any]:
    try:
        manifest_raw = json.loads(manifest_path.read_text(encoding="utf-8"))
        snapshot_raw = json.loads(snapshot_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ProjectionError("cannot read storage manifest or quota snapshot") from exc
    manifest = validate_storage_manifest(
        require_mapping(manifest_raw, "storage manifest")
    )
    snapshot = validate_quota_snapshot(
        require_mapping(snapshot_raw, "quota snapshot")
    )
    measurement = require_mapping(
        manifest.get("artifact_measurement"), "artifact measurement"
    )
    blockers = list(manifest.get("blocker_codes", []))
    blockers.extend(snapshot.get("blocker_codes", []))
    plan_binding = require_mapping(
        require_mapping(
            measurement.get("bindings"), "measurement.bindings"
        ).get("plan"),
        "measurement plan binding",
    )
    snapshot_plan = require_mapping(snapshot.get("plan"), "snapshot plan binding")
    if (
        snapshot_plan.get("sha256") != plan_binding.get("sha256")
        or snapshot_plan.get("path") != plan_binding.get("path")
    ):
        blockers.append("PLAN_BINDING_DRIFT")

    observed = require_nonnegative_int(
        snapshot.get("observed_unix_seconds"), "snapshot observed time"
    )
    max_age = require_positive_int(
        snapshot.get("max_age_seconds"), "snapshot max age"
    )
    now = int(time.time()) if now_seconds is None else now_seconds
    require_nonnegative_int(now, "projection time")
    if now < observed or now - observed > max_age:
        blockers.append("QUOTA_SNAPSHOT_STALE")

    domain_totals = require_mapping(
        measurement.get("storage_domain_totals"),
        "measurement storage-domain totals",
    )
    manifest_roots = require_mapping(
        measurement.get("storage_domain_roots"),
        "measurement storage-domain roots",
    )
    if set(domain_totals) != set(STORAGE_DOMAINS):
        blockers.append("QUOTA_DOMAIN_UNBOUND")
    snapshot_domains_raw = require_sequence(
        snapshot.get("storage_domains"), "quota snapshot storage domains"
    )
    snapshot_by_domain: dict[str, dict[str, Any]] = {}
    for raw in snapshot_domains_raw:
        record = require_mapping(raw, "quota snapshot domain")
        domain = str(record.get("storage_domain", ""))
        if domain in snapshot_by_domain:
            raise ProjectionError(f"quota snapshot duplicates domain {domain}")
        snapshot_by_domain[domain] = record
    if set(snapshot_by_domain) != set(STORAGE_DOMAINS):
        blockers.append("QUOTA_DOMAIN_UNBOUND")

    quota_groups: dict[str, dict[str, Any]] = {}
    for domain in STORAGE_DOMAINS:
        quota_record = snapshot_by_domain.get(domain)
        totals = domain_totals.get(domain)
        if quota_record is None or not isinstance(totals, dict):
            continue
        if quota_record.get("intended_root") != manifest_roots.get(domain):
            blockers.append("PLAN_BINDING_DRIFT")
        if quota_record.get("usage_authoritative") is not True:
            blockers.append("QUOTA_USAGE_UNAUTHORITATIVE")
        if quota_record.get("namespace_fresh") is not True:
            blockers.append("NAMESPACE_NOT_FRESH")
        quota_domain_id = quota_record.get("quota_domain_id")
        if not isinstance(quota_domain_id, str) or not quota_domain_id:
            blockers.append("QUOTA_DOMAIN_UNBOUND")
            continue
        used_bytes = quota_record.get("used_bytes")
        quota_bytes = quota_record.get("quota_bytes")
        used_inodes = quota_record.get("used_inodes")
        quota_inodes = quota_record.get("quota_inodes")
        if not all(
            isinstance(value, int) and not isinstance(value, bool) and value >= 0
            for value in (used_bytes, quota_bytes, used_inodes, quota_inodes)
        ) or quota_bytes == 0 or quota_inodes == 0:
            blockers.append("QUOTA_USAGE_UNAUTHORITATIVE")
            continue
        group = quota_groups.get(quota_domain_id)
        if group is None:
            group = {
                "quota_domain_id": quota_domain_id,
                "filesystem": quota_record.get("filesystem"),
                "storage_domains": [],
                "used_bytes": used_bytes,
                "quota_bytes": quota_bytes,
                "used_inodes": used_inodes,
                "quota_inodes": quota_inodes,
                "quota_observation_count": 0,
                "projected_increment_bytes": 0,
                "projected_increment_inodes": 0,
            }
            quota_groups[quota_domain_id] = group
        elif group["filesystem"] != quota_record.get("filesystem"):
            raise ProjectionError(
                f"shared quota domain {quota_domain_id} filesystem differs"
            )
        group["used_bytes"] = max(group["used_bytes"], used_bytes)
        group["quota_bytes"] = min(group["quota_bytes"], quota_bytes)
        group["used_inodes"] = max(group["used_inodes"], used_inodes)
        group["quota_inodes"] = min(group["quota_inodes"], quota_inodes)
        group["quota_observation_count"] += 1
        group["storage_domains"].append(domain)
        group["projected_increment_bytes"] += require_nonnegative_int(
            totals.get("projected_increment_bytes"),
            f"{domain}.projected_increment_bytes",
        )
        group["projected_increment_inodes"] += require_nonnegative_int(
            totals.get("projected_increment_inodes"),
            f"{domain}.projected_increment_inodes",
        )

    projections: list[dict[str, Any]] = []
    for quota_domain_id in sorted(quota_groups):
        group = quota_groups[quota_domain_id]
        projected_used_bytes = (
            group["used_bytes"] + group["projected_increment_bytes"]
        )
        projected_free_bytes = group["quota_bytes"] - projected_used_bytes
        projected_used_inodes = (
            group["used_inodes"] + group["projected_increment_inodes"]
        )
        projected_free_inodes = group["quota_inodes"] - projected_used_inodes
        byte_headroom_pass = (
            projected_free_bytes >= 0
            and projected_free_bytes * MIN_HEADROOM_DENOMINATOR
            >= group["quota_bytes"] * MIN_HEADROOM_NUMERATOR
        )
        inode_headroom_pass = (
            projected_free_inodes >= 0
            and projected_free_inodes * MIN_HEADROOM_DENOMINATOR
            >= group["quota_inodes"] * MIN_HEADROOM_NUMERATOR
        )
        if not byte_headroom_pass:
            blockers.append("BYTE_HEADROOM_LT_20PCT")
        if not inode_headroom_pass:
            blockers.append("INODE_HEADROOM_LT_20PCT")
        projections.append(
            {
                **group,
                "storage_domains": sorted(group["storage_domains"]),
                "projected_used_bytes": projected_used_bytes,
                "projected_free_bytes": projected_free_bytes,
                "projected_used_inodes": projected_used_inodes,
                "projected_free_inodes": projected_free_inodes,
                "byte_headroom_at_least_20_percent": byte_headroom_pass,
                "inode_headroom_at_least_20_percent": inode_headroom_pass,
            }
        )

    blocker_codes = ordered_blockers(blockers)
    evidence_blockers = [
        code
        for code in blocker_codes
        if code
        not in {"BYTE_HEADROOM_LT_20PCT", "INODE_HEADROOM_LT_20PCT"}
    ]
    capacity_blockers = [
        code
        for code in blocker_codes
        if code in {"BYTE_HEADROOM_LT_20PCT", "INODE_HEADROOM_LT_20PCT"}
    ]
    if evidence_blockers:
        status = BLOCKED_STATUS
    elif capacity_blockers:
        status = FAIL_STATUS
    else:
        status = PASS_STATUS
    payload = apply_authority(
        {
            "schema": STORAGE_CERTIFICATE_SCHEMA,
            "status": status,
            "gate": GATE,
            "manifest": {
                "path": str(manifest_path),
                "sha256": file_sha256(manifest_path),
                "manifest_semantic_sha256": manifest[
                    "manifest_semantic_sha256"
                ],
            },
            "quota_snapshot": {
                "path": str(snapshot_path),
                "sha256": file_sha256(snapshot_path),
                "snapshot_semantic_sha256": snapshot[
                    "snapshot_semantic_sha256"
                ],
            },
            "minimum_headroom": {
                "numerator": MIN_HEADROOM_NUMERATOR,
                "denominator": MIN_HEADROOM_DENOMINATOR,
            },
            "quota_domain_projections": projections,
            "blocker_codes": blocker_codes,
            "boundaries": [
                "This certificate grants storage admission only.",
                "It does not submit or materialize Condor jobs.",
                "It does not grant full extraction, science-freeze, broad-production, or CANONICAL authority.",
                "P5A-S science-freeze evidence is forbidden as an input to this pre-extraction gate.",
            ],
        }
    )
    payload["certificate_semantic_sha256"] = semantic_sha256(payload)
    return payload


def canonical_path_with_missing_tail(path: Path) -> Path:
    candidate = Path(os.path.abspath(path))
    missing_tail: list[str] = []
    while not os.path.lexists(candidate):
        if candidate.parent == candidate:
            raise ProjectionError(f"cannot resolve output boundary: {path}")
        missing_tail.append(candidate.name)
        candidate = candidate.parent
    try:
        resolved = candidate.resolve(strict=True)
    except OSError as exc:
        raise ProjectionError(f"cannot resolve output boundary: {path}") from exc
    for part in reversed(missing_tail):
        resolved /= part
    return resolved


def path_is_equal_or_within(path: Path, root: Path) -> bool:
    candidate = canonical_path_with_missing_tail(path)
    boundary = canonical_path_with_missing_tail(root)
    try:
        candidate.relative_to(boundary)
    except ValueError:
        return False
    return True


def write_json(
    path: Path,
    payload: Mapping[str, Any],
    *,
    forbidden_roots: Iterable[str] = (),
) -> None:
    if not path.is_absolute() or ".." in path.parts:
        raise ProjectionError("JSON output must be one absolute non-traversing path")
    for root_text in forbidden_roots:
        root = Path(require_absolute_remote_root(root_text, "forbidden output root"))
        if path_is_equal_or_within(path, root):
            raise ProjectionError(
                f"JSON output must remain outside intended campaign root: {root}"
            )
    parent = path.parent
    if not parent.exists() or not parent.is_dir() or parent.is_symlink():
        raise ProjectionError(
            "JSON output parent must be one existing non-symlink directory"
        )
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    try:
        descriptor = os.open(path, flags, 0o444)
    except OSError as exc:
        raise ProjectionError(
            f"JSON output must be new and non-symlinked: {path}"
        ) from exc
    try:
        data = canonical_json_bytes(payload)
        offset = 0
        while offset < len(data):
            written = os.write(descriptor, data[offset:])
            if written <= 0:
                raise ProjectionError(f"short write for JSON output: {path}")
            offset += written
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def forbidden_roots_for_output(
    action: str,
    payload: Mapping[str, Any],
    args: argparse.Namespace,
) -> list[str]:
    if action == "measure":
        measurement = require_mapping(
            payload.get("artifact_measurement"), "artifact measurement"
        )
        roots = require_mapping(
            measurement.get("storage_domain_roots"),
            "measurement storage-domain roots",
        )
        return [str(roots[domain]) for domain in STORAGE_DOMAINS]
    if action == "snapshot":
        domains = require_sequence(
            payload.get("storage_domains"), "quota snapshot storage domains"
        )
        return [
            str(
                require_mapping(record, "quota snapshot domain").get(
                    "intended_root"
                )
            )
            for record in domains
        ]
    try:
        manifest_raw = json.loads(
            args.manifest.read_text(encoding="utf-8")
        )
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ProjectionError(
            "cannot reread storage manifest for output-path safety"
        ) from exc
    manifest = validate_storage_manifest(
        require_mapping(manifest_raw, "storage manifest")
    )
    measurement = require_mapping(
        manifest.get("artifact_measurement"), "artifact measurement"
    )
    roots = require_mapping(
        measurement.get("storage_domain_roots"),
        "measurement storage-domain roots",
    )
    return [str(roots[domain]) for domain in STORAGE_DOMAINS]


def exit_for_status(status: str) -> int:
    if status == PASS_STATUS:
        return EXIT_PASS
    if status == FAIL_STATUS:
        return EXIT_CAPACITY
    if status == BLOCKED_STATUS:
        return EXIT_BLOCKED
    raise ProjectionError(f"unknown result status {status!r}")


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="action", required=True)
    measure = subparsers.add_parser(
        "measure", help="rehash evidence and build the storage manifest"
    )
    measure.add_argument("--manifest", type=Path, required=True)
    measure.add_argument("--json-out", type=Path, required=True)
    snapshot = subparsers.add_parser(
        "snapshot", help="capture live authoritative Lustre quota state"
    )
    snapshot.add_argument("--plan", type=Path, required=True)
    snapshot.add_argument("--principal", required=True)
    snapshot.add_argument("--json-out", type=Path, required=True)
    snapshot.add_argument(
        "--max-age-seconds",
        type=int,
        default=DEFAULT_SNAPSHOT_MAX_AGE_SECONDS,
    )
    project = subparsers.add_parser(
        "project", help="apply byte and inode admission gates"
    )
    project.add_argument("--manifest", type=Path, required=True)
    project.add_argument("--quota-snapshot", type=Path, required=True)
    project.add_argument("--json-out", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        if args.action == "measure":
            payload = build_storage_manifest(args.manifest)
        elif args.action == "snapshot":
            payload = build_quota_snapshot(
                args.plan,
                args.principal,
                max_age_seconds=args.max_age_seconds,
            )
        else:
            payload = build_projection(args.manifest, args.quota_snapshot)
        write_json(
            args.json_out,
            payload,
            forbidden_roots=forbidden_roots_for_output(
                args.action, payload, args
            ),
        )
        print(args.json_out)
        return exit_for_status(str(payload["status"]))
    except ProjectionError as exc:
        print(f"THE-134 pre-extraction storage gate failed: {exc}", file=sys.stderr)
        return EXIT_MALFORMED


if __name__ == "__main__":
    raise SystemExit(main())
