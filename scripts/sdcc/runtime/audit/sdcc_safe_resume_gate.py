#!/usr/bin/env python3
"""Fail-closed SDCC shared-resource admission certificate.

The gate is deliberately detector- and campaign-semantics neutral.  It binds
already-produced scientific, permission, duplicate, resource, and execution
evidence to one small SDCC operation.  It performs no network access, queue
query, filesystem walk, submission, or job control.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import socket
import sys
import time
from pathlib import Path
from typing import Any, Iterable, Mapping


CONTRACT_SCHEMA = "SDCC_SAFE_RESUME_CONTRACT_V1"
CERTIFICATE_SCHEMA = "SDCC_SAFE_RESUME_CERTIFICATE_V1"
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
SAFE_TEXT_RE = re.compile(r"^[A-Za-z0-9_.:/+-]+$")
ALLOWED_HOST = "sphnxuser05"
ALLOWED_SERVICES = frozenset(
    {"sdcc_ssh", "sdcc_sftp", "condor_schedd", "sphenix_cvmfs"}
)
MAX_EXPLICIT_PATHS = 256
MAX_MANIFEST_VALIDATION_ROWS = 256
NORMAL_MANIFEST_VALIDATION_ROWS = 32
MAX_GROUP_FILES_PER_ROW = 2048
MAX_SUBMITTER_SECONDS = 600
MAX_CAPTURE_BYTES = 256 * 1024
MIN_MEMORY_HEADROOM = 1.5
MAX_MEMORY_HEADROOM = 3.0
MAX_CERTIFICATE_LIFETIME_SECONDS = 3600


class SafeResumeError(RuntimeError):
    """A contract cannot safely authorize the declared SDCC operation."""


def canonical_json_bytes(payload: Any) -> bytes:
    return (
        json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
        + "\n"
    ).encode("utf-8")


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
    except OSError as exc:
        raise SafeResumeError(f"cannot read artifact: {path}") from exc
    return digest.hexdigest()


def verify_artifact_record(
    record: object,
    label: str,
    *,
    require_assertion_count: bool,
) -> dict[str, Any]:
    """Reopen one certificate-bound artifact without traversing its tree."""

    item = require_mapping(record, label)
    expected_keys = {"path", "sha256", "size_bytes"}
    if require_assertion_count:
        expected_keys.add("assertion_count")
    require_exact_keys(item, expected_keys, label)
    raw_path = item.get("path")
    if not isinstance(raw_path, str) or not Path(raw_path).is_absolute():
        raise SafeResumeError(f"{label}.path must be absolute")
    path = Path(raw_path)
    if not path.is_file() or path.is_symlink():
        raise SafeResumeError(f"{label} must be a regular non-symlink file")
    expected_sha256 = require_sha256(item.get("sha256"), f"{label}.sha256")
    size_bytes = require_int(item.get("size_bytes"), f"{label}.size_bytes")
    if path.stat().st_size != size_bytes or file_sha256(path) != expected_sha256:
        raise SafeResumeError(f"{label} artifact differs")
    if require_assertion_count:
        require_int(
            item.get("assertion_count"),
            f"{label}.assertion_count",
            minimum=1,
        )
    return dict(item)


def require_sha256(value: object, label: str) -> str:
    if not isinstance(value, str) or SHA256_RE.fullmatch(value) is None:
        raise SafeResumeError(f"{label} must be a lowercase SHA-256")
    return value


def require_mapping(value: object, label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise SafeResumeError(f"{label} must be an object")
    return value


def require_exact_keys(
    payload: Mapping[str, Any], expected: Iterable[str], label: str
) -> None:
    wanted = set(expected)
    observed = set(payload)
    if observed != wanted:
        raise SafeResumeError(
            f"{label} key inventory differs: "
            f"missing={sorted(wanted - observed)} extra={sorted(observed - wanted)}"
        )


def require_int(value: object, label: str, *, minimum: int = 0) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise SafeResumeError(f"{label} must be an integer >= {minimum}")
    return value


def safe_text(value: object, label: str) -> str:
    if (
        not isinstance(value, str)
        or not value
        or value != value.strip()
        or SAFE_TEXT_RE.fullmatch(value) is None
    ):
        raise SafeResumeError(f"{label} is unsafe")
    return value


def load_json_artifact(record: object, label: str) -> tuple[dict[str, Any], dict[str, Any]]:
    item = require_mapping(record, label)
    require_exact_keys(item, {"path", "sha256", "assertions"}, label)
    raw_path = item.get("path")
    if not isinstance(raw_path, str) or not Path(raw_path).is_absolute():
        raise SafeResumeError(f"{label}.path must be absolute")
    path = Path(raw_path)
    if not path.is_file() or path.is_symlink():
        raise SafeResumeError(f"{label} must be a regular non-symlink file")
    expected = require_sha256(item.get("sha256"), f"{label}.sha256")
    observed = file_sha256(path)
    if observed != expected:
        raise SafeResumeError(f"{label} SHA-256 differs")
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise SafeResumeError(f"{label} is not readable JSON") from exc
    payload = require_mapping(payload, f"{label} payload")
    assertions = item.get("assertions")
    if not isinstance(assertions, list) or not assertions:
        raise SafeResumeError(f"{label}.assertions must be a non-empty array")
    for index, raw in enumerate(assertions):
        assertion = require_mapping(raw, f"{label}.assertions[{index}]")
        require_exact_keys(assertion, {"path", "equals"}, "artifact assertion")
        dotted = assertion.get("path")
        if not isinstance(dotted, str) or not dotted or dotted.startswith("."):
            raise SafeResumeError(f"{label} assertion path is invalid")
        cursor: object = payload
        for component in dotted.split("."):
            if not isinstance(cursor, dict) or component not in cursor:
                raise SafeResumeError(f"{label} assertion path is missing: {dotted}")
            cursor = cursor[component]
        if cursor != assertion.get("equals"):
            raise SafeResumeError(f"{label} assertion differs: {dotted}")
    return payload, {
        "path": str(path),
        "sha256": observed,
        "size_bytes": path.stat().st_size,
        "assertion_count": len(assertions),
    }


def validate_operation(operation: Mapping[str, Any]) -> dict[str, Any]:
    require_exact_keys(
        operation,
        {
            "class",
            "campaign_tag",
            "submit_host",
            "row_count",
            "logical_job_count",
            "cluster_count",
            "group_size",
            "request_memory_mb",
            "max_materialize_per_cluster",
            "max_idle_per_cluster",
            "attended_submission",
            "automatic_job_control",
            "auto_memory_retry",
            "hold_failed_workers",
        },
        "operation",
    )
    operation_class = operation.get("class")
    if operation_class not in {"watched_capacity_canary", "full_extraction"}:
        raise SafeResumeError("operation.class is unsupported")
    campaign_tag = safe_text(operation.get("campaign_tag"), "campaign tag")
    if operation.get("submit_host") != ALLOWED_HOST:
        raise SafeResumeError("submit host differs from the administrator-approved host")
    if operation.get("attended_submission") is not True:
        raise SafeResumeError("submission must be attended")
    if operation.get("automatic_job_control") is not False:
        raise SafeResumeError("automatic job control is forbidden")
    if operation.get("auto_memory_retry") is not False:
        raise SafeResumeError("automatic memory retry is forbidden")
    if operation.get("hold_failed_workers") is not True:
        raise SafeResumeError("failed workers must remain held and visible")
    limits = {
        key: require_int(operation.get(key), f"operation.{key}", minimum=1)
        for key in (
            "row_count",
            "logical_job_count",
            "cluster_count",
            "group_size",
            "request_memory_mb",
            "max_materialize_per_cluster",
            "max_idle_per_cluster",
        )
    }
    if limits["max_idle_per_cluster"] > limits["max_materialize_per_cluster"]:
        raise SafeResumeError("max idle exceeds max materialize")
    if operation_class == "watched_capacity_canary":
        expected = {
            "row_count": 2,
            "logical_job_count": 2,
            "cluster_count": 2,
            "group_size": 7,
            "request_memory_mb": 3000,
            "max_materialize_per_cluster": 1,
            "max_idle_per_cluster": 1,
        }
    else:
        expected = {
            "row_count": 13,
            "logical_job_count": 18577,
            "cluster_count": 13,
            "group_size": 7,
            "request_memory_mb": 3000,
            "max_materialize_per_cluster": 20,
            "max_idle_per_cluster": 5,
        }
    if limits != expected:
        raise SafeResumeError(f"{operation_class} limits differ")
    return {"class": operation_class, "campaign_tag": campaign_tag, **limits}


def validate_login_node(payload: Mapping[str, Any]) -> dict[str, Any]:
    require_exact_keys(
        payload,
        {
            "explicit_path_count",
            "max_explicit_paths",
            "manifest_validation_rows",
            "hard_manifest_validation_limit",
            "max_group_files_per_row",
            "tree_walk",
            "unbounded_glob",
            "per_manifest_row_subprocess",
            "max_concurrent_submitters",
            "submitter_timeout_seconds",
            "capture_limit_bytes",
            "high_cardinality_location",
        },
        "login_node",
    )
    if payload.get("max_explicit_paths") != MAX_EXPLICIT_PATHS:
        raise SafeResumeError("explicit path hard limit differs")
    if require_int(payload.get("explicit_path_count"), "explicit_path_count") > MAX_EXPLICIT_PATHS:
        raise SafeResumeError("explicit path count exceeds the hard limit")
    if payload.get("hard_manifest_validation_limit") != MAX_MANIFEST_VALIDATION_ROWS:
        raise SafeResumeError("manifest validation hard limit differs")
    rows = require_int(payload.get("manifest_validation_rows"), "manifest_validation_rows")
    if rows > NORMAL_MANIFEST_VALIDATION_ROWS:
        raise SafeResumeError("routine manifest validation exceeds 32 rows")
    if payload.get("max_group_files_per_row") != MAX_GROUP_FILES_PER_ROW:
        raise SafeResumeError("group-file hard limit differs")
    for field in ("tree_walk", "unbounded_glob", "per_manifest_row_subprocess"):
        if payload.get(field) is not False:
            raise SafeResumeError(f"login-node {field} must be false")
    if payload.get("max_concurrent_submitters") != 1:
        raise SafeResumeError("submitters must execute serially")
    timeout = require_int(payload.get("submitter_timeout_seconds"), "submitter timeout", minimum=1)
    if timeout > MAX_SUBMITTER_SECONDS:
        raise SafeResumeError("submitter timeout exceeds the hard limit")
    capture = require_int(payload.get("capture_limit_bytes"), "capture limit", minimum=1)
    if capture > MAX_CAPTURE_BYTES:
        raise SafeResumeError("capture limit exceeds the hard limit")
    if payload.get("high_cardinality_location") != "condor":
        raise SafeResumeError("high-cardinality work must run through Condor")
    return {
        "explicit_path_count": payload["explicit_path_count"],
        "manifest_validation_rows": rows,
        "max_concurrent_submitters": 1,
        "submitter_timeout_seconds": timeout,
        "capture_limit_bytes": capture,
        "tree_walk_performed": False,
    }


def validate_network(payload: Mapping[str, Any]) -> dict[str, Any]:
    require_exact_keys(
        payload,
        {
            "allowed_services",
            "scanning",
            "arbitrary_upload",
            "host_switch",
            "model_directed_exploration",
        },
        "network",
    )
    services = payload.get("allowed_services")
    if (
        not isinstance(services, list)
        or not services
        or len(services) != len(set(services))
        or not set(services).issubset(ALLOWED_SERVICES)
    ):
        raise SafeResumeError("network allowlist differs")
    for field in ("scanning", "arbitrary_upload", "host_switch", "model_directed_exploration"):
        if payload.get(field) is not False:
            raise SafeResumeError(f"network {field} must be false")
    return {"allowed_services": sorted(services), "all_exploration_disabled": True}


def validate_resources(payload: Mapping[str, Any], operation: Mapping[str, Any]) -> dict[str, Any]:
    require_exact_keys(
        payload,
        {
            "measured_peak_memory_mb",
            "measured_witness_count",
            "request_memory_mb",
            "minimum_headroom_ratio",
            "maximum_headroom_ratio",
            "automatic_widening",
            "measurement_receipt",
        },
        "resources",
    )
    measured = require_int(payload.get("measured_peak_memory_mb"), "measured peak memory", minimum=1)
    witnesses = require_int(payload.get("measured_witness_count"), "measured witness count", minimum=2)
    requested = require_int(payload.get("request_memory_mb"), "requested memory", minimum=1)
    if requested != operation["request_memory_mb"]:
        raise SafeResumeError("resource and operation memory requests differ")
    if payload.get("minimum_headroom_ratio") != MIN_MEMORY_HEADROOM:
        raise SafeResumeError("minimum memory headroom contract differs")
    if payload.get("maximum_headroom_ratio") != MAX_MEMORY_HEADROOM:
        raise SafeResumeError("maximum memory headroom contract differs")
    ratio = requested / measured
    if ratio < MIN_MEMORY_HEADROOM or ratio > MAX_MEMORY_HEADROOM:
        raise SafeResumeError(f"memory request ratio is outside [{MIN_MEMORY_HEADROOM},{MAX_MEMORY_HEADROOM}]")
    if payload.get("automatic_widening") is not False:
        raise SafeResumeError("automatic resource widening is forbidden")
    _, artifact = load_json_artifact(payload.get("measurement_receipt"), "measurement receipt")
    return {
        "measured_peak_memory_mb": measured,
        "measured_witness_count": witnesses,
        "request_memory_mb": requested,
        "headroom_ratio": ratio,
        "measurement_receipt": artifact,
    }


def validate_authority(payload: Mapping[str, Any], operation_class: str) -> dict[str, Any]:
    require_exact_keys(
        payload,
        {
            "site_admin_acknowledged",
            "user_approved",
            "exact_scope",
            "submission_authority",
            "full_extraction_authority",
            "the121_authority",
            "the122_authority",
            "canonical_promotion",
        },
        "authority",
    )
    if payload.get("site_admin_acknowledged") is not True or payload.get("user_approved") is not True:
        raise SafeResumeError("site and user resume authority are required")
    if safe_text(payload.get("exact_scope"), "exact authority scope") != operation_class:
        raise SafeResumeError("authority scope differs from operation")
    if payload.get("submission_authority") is not True:
        raise SafeResumeError("exact submission authority is missing")
    if payload.get("full_extraction_authority") is not (operation_class == "full_extraction"):
        raise SafeResumeError("full-extraction authority differs")
    for field in ("the121_authority", "the122_authority", "canonical_promotion"):
        if payload.get(field) is not False:
            raise SafeResumeError(f"{field} must remain false")
    return dict(payload)


def validate_lifecycle(payload: Mapping[str, Any]) -> dict[str, Any]:
    expected = {
        "LOCAL_CHECK": "PASS",
        "SCIENTIFIC_CERTIFICATE": "PASS",
        "SITE_ADMISSION": "PENDING",
        "SUBMISSION": "NOT_STARTED",
        "RUNNING": "NOT_STARTED",
        "TERMINAL": "NOT_STARTED",
        "PRODUCTION_AUTHORIZED": False,
    }
    require_exact_keys(payload, expected, "lifecycle")
    if dict(payload) != expected:
        raise SafeResumeError("lifecycle states are conflated or premature")
    return expected


def certify_contract(contract_path: Path, *, now_seconds: int | None = None) -> dict[str, Any]:
    if not contract_path.is_file() or contract_path.is_symlink():
        raise SafeResumeError("contract must be a regular non-symlink file")
    try:
        contract = json.loads(contract_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise SafeResumeError("contract is not readable JSON") from exc
    contract = require_mapping(contract, "contract")
    require_exact_keys(
        contract,
        {
            "schema",
            "status",
            "operation",
            "login_node",
            "network",
            "resources",
            "evidence",
            "authority",
            "lifecycle",
            "expires_at_unix",
        },
        "contract",
    )
    if contract.get("schema") != CONTRACT_SCHEMA or contract.get("status") != "PREPARED_NOT_SUBMITTED":
        raise SafeResumeError("contract schema/status differs")
    now = int(time.time()) if now_seconds is None else now_seconds
    expires = require_int(contract.get("expires_at_unix"), "contract expiry", minimum=1)
    if expires <= now or expires > now + MAX_CERTIFICATE_LIFETIME_SECONDS:
        raise SafeResumeError("contract expiry is stale or too far in the future")
    operation = validate_operation(require_mapping(contract.get("operation"), "operation"))
    login_node = validate_login_node(require_mapping(contract.get("login_node"), "login_node"))
    network = validate_network(require_mapping(contract.get("network"), "network"))
    resources = validate_resources(require_mapping(contract.get("resources"), "resources"), operation)
    evidence = require_mapping(contract.get("evidence"), "evidence")
    expected_evidence = {
        "scientific_certificate",
        "permission_receipt",
        "duplicate_receipt",
        "execution_binding",
        "capacity_certificate",
        "quota_certificate",
    }
    require_exact_keys(evidence, expected_evidence, "evidence")
    artifacts: dict[str, Any] = {}
    for name in sorted(expected_evidence - {"quota_certificate"}):
        _, artifacts[name] = load_json_artifact(evidence[name], name.replace("_", " "))
    if operation["class"] == "full_extraction":
        _, artifacts["quota_certificate"] = load_json_artifact(
            evidence["quota_certificate"], "quota certificate"
        )
    elif evidence["quota_certificate"] is not None:
        raise SafeResumeError("capacity canary must not claim full quota authority")
    authority = validate_authority(require_mapping(contract.get("authority"), "authority"), operation["class"])
    lifecycle = validate_lifecycle(require_mapping(contract.get("lifecycle"), "lifecycle"))
    gate_path = Path(__file__).resolve(strict=True)
    return {
        "schema": CERTIFICATE_SCHEMA,
        "status": "SITE_ADMISSION_PASS",
        "operation": operation,
        "submit_host": ALLOWED_HOST,
        "login_node": login_node,
        "network": network,
        "resources": resources,
        "evidence": artifacts,
        "authority": authority,
        "lifecycle": {**lifecycle, "SITE_ADMISSION": "PASS"},
        "contract": {
            "path": str(contract_path.resolve(strict=True)),
            "sha256": file_sha256(contract_path),
            "size_bytes": contract_path.stat().st_size,
        },
        "gate": {
            "path": str(gate_path),
            "sha256": file_sha256(gate_path),
            "size_bytes": gate_path.stat().st_size,
        },
        "issued_at_unix": now,
        "expires_at_unix": expires,
        "submission_performed": False,
    }


def verify_certificate(
    path: Path,
    expected_sha256: str,
    *,
    operation_class: str,
    campaign_tag: str,
    submit_host: str,
    now_seconds: int | None = None,
) -> dict[str, Any]:
    expected = require_sha256(expected_sha256, "certificate SHA-256")
    if not path.is_file() or path.is_symlink() or file_sha256(path) != expected:
        raise SafeResumeError("certificate artifact differs")
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise SafeResumeError("certificate is not readable JSON") from exc
    payload = require_mapping(payload, "certificate")
    require_exact_keys(
        payload,
        {
            "schema",
            "status",
            "operation",
            "submit_host",
            "login_node",
            "network",
            "resources",
            "evidence",
            "authority",
            "lifecycle",
            "contract",
            "gate",
            "issued_at_unix",
            "expires_at_unix",
            "submission_performed",
        },
        "certificate",
    )
    if payload.get("schema") != CERTIFICATE_SCHEMA or payload.get("status") != "SITE_ADMISSION_PASS":
        raise SafeResumeError("certificate schema/status differs")
    operation = require_mapping(payload.get("operation"), "certificate operation")
    require_exact_keys(
        operation,
        {
            "class",
            "campaign_tag",
            "row_count",
            "logical_job_count",
            "cluster_count",
            "group_size",
            "request_memory_mb",
            "max_materialize_per_cluster",
            "max_idle_per_cluster",
        },
        "certificate operation",
    )
    if operation.get("class") != operation_class or operation.get("campaign_tag") != campaign_tag:
        raise SafeResumeError("certificate operation binding differs")
    expected_limits = (
        {
            "row_count": 2,
            "logical_job_count": 2,
            "cluster_count": 2,
            "group_size": 7,
            "request_memory_mb": 3000,
            "max_materialize_per_cluster": 1,
            "max_idle_per_cluster": 1,
        }
        if operation_class == "watched_capacity_canary"
        else {
            "row_count": 13,
            "logical_job_count": 18577,
            "cluster_count": 13,
            "group_size": 7,
            "request_memory_mb": 3000,
            "max_materialize_per_cluster": 20,
            "max_idle_per_cluster": 5,
        }
    )
    if operation_class not in {"watched_capacity_canary", "full_extraction"}:
        raise SafeResumeError("certificate operation class is unsupported")
    if {key: operation.get(key) for key in expected_limits} != expected_limits:
        raise SafeResumeError("certificate operation limits differ")
    if payload.get("submit_host") != submit_host or submit_host != ALLOWED_HOST:
        raise SafeResumeError("certificate submit-host binding differs")
    login_node = require_mapping(payload.get("login_node"), "certificate login node")
    require_exact_keys(
        login_node,
        {
            "explicit_path_count",
            "manifest_validation_rows",
            "max_concurrent_submitters",
            "submitter_timeout_seconds",
            "capture_limit_bytes",
            "tree_walk_performed",
        },
        "certificate login node",
    )
    if (
        require_int(login_node.get("explicit_path_count"), "certificate explicit paths")
        > MAX_EXPLICIT_PATHS
        or require_int(
            login_node.get("manifest_validation_rows"),
            "certificate manifest rows",
        )
        > NORMAL_MANIFEST_VALIDATION_ROWS
        or login_node.get("max_concurrent_submitters") != 1
        or require_int(
            login_node.get("submitter_timeout_seconds"),
            "certificate submitter timeout",
            minimum=1,
        )
        > MAX_SUBMITTER_SECONDS
        or require_int(
            login_node.get("capture_limit_bytes"),
            "certificate capture limit",
            minimum=1,
        )
        > MAX_CAPTURE_BYTES
        or login_node.get("tree_walk_performed") is not False
    ):
        raise SafeResumeError("certificate login-node budget differs")
    network = require_mapping(payload.get("network"), "certificate network")
    if network != {
        "allowed_services": sorted(network.get("allowed_services", [])),
        "all_exploration_disabled": True,
    }:
        raise SafeResumeError("certificate network boundary differs")
    if not set(network["allowed_services"]).issubset(ALLOWED_SERVICES):
        raise SafeResumeError("certificate network allowlist differs")
    resources = require_mapping(payload.get("resources"), "certificate resources")
    require_exact_keys(
        resources,
        {
            "measured_peak_memory_mb",
            "measured_witness_count",
            "request_memory_mb",
            "headroom_ratio",
            "measurement_receipt",
        },
        "certificate resources",
    )
    measured = require_int(
        resources.get("measured_peak_memory_mb"),
        "certificate measured memory",
        minimum=1,
    )
    requested = require_int(
        resources.get("request_memory_mb"),
        "certificate requested memory",
        minimum=1,
    )
    ratio = requested / measured
    if (
        require_int(
            resources.get("measured_witness_count"),
            "certificate witness count",
            minimum=2,
        )
        < 2
        or requested != operation["request_memory_mb"]
        or not isinstance(resources.get("headroom_ratio"), (int, float))
        or abs(float(resources["headroom_ratio"]) - ratio) > 1e-12
        or ratio < MIN_MEMORY_HEADROOM
        or ratio > MAX_MEMORY_HEADROOM
    ):
        raise SafeResumeError("certificate resource authority differs")
    verify_artifact_record(
        resources.get("measurement_receipt"),
        "certificate measurement receipt",
        require_assertion_count=True,
    )
    evidence = require_mapping(payload.get("evidence"), "certificate evidence")
    expected_evidence = {
        "scientific_certificate",
        "permission_receipt",
        "duplicate_receipt",
        "execution_binding",
        "capacity_certificate",
    }
    if operation_class == "full_extraction":
        expected_evidence.add("quota_certificate")
    require_exact_keys(evidence, expected_evidence, "certificate evidence")
    for name in sorted(expected_evidence):
        verify_artifact_record(
            evidence[name],
            f"certificate {name.replace('_', ' ')}",
            require_assertion_count=True,
        )
    authority = validate_authority(
        require_mapping(payload.get("authority"), "certificate authority"),
        operation_class,
    )
    if authority != payload["authority"]:
        raise SafeResumeError("certificate authority differs")
    lifecycle = require_mapping(payload.get("lifecycle"), "certificate lifecycle")
    expected_lifecycle = {
        "LOCAL_CHECK": "PASS",
        "SCIENTIFIC_CERTIFICATE": "PASS",
        "SITE_ADMISSION": "PASS",
        "SUBMISSION": "NOT_STARTED",
        "RUNNING": "NOT_STARTED",
        "TERMINAL": "NOT_STARTED",
        "PRODUCTION_AUTHORIZED": False,
    }
    if dict(lifecycle) != expected_lifecycle or payload.get("submission_performed") is not False:
        raise SafeResumeError("certificate lifecycle authority differs")
    now = int(time.time()) if now_seconds is None else now_seconds
    issued = require_int(payload.get("issued_at_unix"), "certificate issue time", minimum=1)
    expires = require_int(payload.get("expires_at_unix"), "certificate expiry", minimum=1)
    if issued > now or expires <= now or expires > issued + MAX_CERTIFICATE_LIFETIME_SECONDS:
        raise SafeResumeError("certificate is expired")
    verify_artifact_record(
        payload.get("contract"),
        "certificate contract",
        require_assertion_count=False,
    )
    gate_record = verify_artifact_record(
        payload.get("gate"),
        "certificate gate",
        require_assertion_count=False,
    )
    current_gate = Path(__file__).resolve(strict=True)
    if (
        gate_record["path"] != str(current_gate)
        or gate_record["sha256"] != file_sha256(current_gate)
        or gate_record["size_bytes"] != current_gate.stat().st_size
    ):
        raise SafeResumeError("certificate gate identity differs")
    return payload


def write_new_json(path: Path, payload: Mapping[str, Any]) -> None:
    if not path.is_absolute() or os.path.lexists(path):
        raise SafeResumeError("output must be a fresh absolute path")
    parent = path.parent
    if not parent.is_dir() or parent.is_symlink():
        raise SafeResumeError("output parent must be an existing real directory")
    data = canonical_json_bytes(payload)
    try:
        descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
        os.fchmod(descriptor, 0o644)
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
    except OSError as exc:
        raise SafeResumeError(f"cannot write certificate: {path}") from exc
    if path.read_bytes() != data:
        raise SafeResumeError("certificate readback differs")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    certify = subparsers.add_parser("certify")
    certify.add_argument("--contract", type=Path, required=True)
    certify.add_argument("--output", type=Path, required=True)
    certify.add_argument("--now-unix", type=int)
    verify = subparsers.add_parser("verify-certificate")
    verify.add_argument("--certificate", type=Path, required=True)
    verify.add_argument("--certificate-sha256", required=True)
    verify.add_argument("--operation-class", required=True)
    verify.add_argument("--campaign-tag", required=True)
    verify.add_argument("--submit-host", default=socket.gethostname().split(".")[0])
    verify.add_argument("--now-unix", type=int)
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.command == "certify":
            certificate = certify_contract(args.contract, now_seconds=args.now_unix)
            write_new_json(args.output, certificate)
            summary = {
                "status": certificate["status"],
                "operation_class": certificate["operation"]["class"],
                "campaign_tag": certificate["operation"]["campaign_tag"],
                "certificate": str(args.output),
                "certificate_sha256": file_sha256(args.output),
                "submission_performed": False,
            }
        else:
            certificate = verify_certificate(
                args.certificate,
                args.certificate_sha256,
                operation_class=args.operation_class,
                campaign_tag=args.campaign_tag,
                submit_host=args.submit_host,
                now_seconds=args.now_unix,
            )
            summary = {
                "status": "SITE_ADMISSION_PASS",
                "operation_class": certificate["operation"]["class"],
                "campaign_tag": certificate["operation"]["campaign_tag"],
                "submission_performed": False,
            }
    except (SafeResumeError, OSError) as exc:
        print(f"[SDCC-SAFE-RESUME][ERROR] {exc}", file=sys.stderr)
        return 2
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
