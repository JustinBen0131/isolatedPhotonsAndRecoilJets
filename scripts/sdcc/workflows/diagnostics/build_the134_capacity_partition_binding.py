#!/usr/bin/env python3
"""Bind a fresh THE-134 group-of-seven capacity canary to its exact plan.

This tool is the current-head successor to the historical count-amendment
receipt.  It deliberately consumes no legacy V1 plan or legacy count alias.
Instead it rehashes and validates one internally consistent evidence chain:

* immutable resolver bundle and materialization receipt;
* current source manifest, plan, rows, duplicate fingerprint, and partition;
* current capacity resource and ROOT/identity certificates; and
* current p+p and Au+Au multiview audits.

The resulting receipt proves that the passed two-row capacity canary is bound
to the exact 13-row, 129,998-tuple, 18,577-job group-of-seven partition.  It
does not submit jobs and grants no extraction, training, broad-production, or
CANONICAL authority.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
from pathlib import Path
from typing import Any, Iterable, Mapping


HERE = Path(__file__).resolve().parent
LEGACY_VALIDATOR_PATH = HERE / "build_the134_capacity_count_amendment.py"


def _load_local_module(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load local module: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# The historical tool contains the already-adversarially-tested validators for
# the current resolver, immutable bundle, runtime snapshot, ROOT health, source
# provenance, and multiview evidence.  Reuse those validators, not its legacy
# V1 amendment schema or its legacy evidence requirements.
evidence = _load_local_module(
    "the134_capacity_evidence_validators", LEGACY_VALIDATOR_PATH
)
resolver = evidence.resolver
source_builder = evidence.source_builder


SPEC_SCHEMA = "THE134_GROUP7_CAPACITY_PARTITION_BINDING_SPEC_V1"
BINDING_SCHEMA = "THE134_GROUP7_CAPACITY_PARTITION_BINDING_V1"
READBACK_SCHEMA = "THE134_GROUP7_CAPACITY_PARTITION_BINDING_READBACK_V1"
AUTHORITY_STATE = "CURRENT_HEAD_CAPACITY_BOUND_NOT_EXTRACTION_AUTHORITY"
GENERATED_BY = "build_the134_capacity_partition_binding.py"

EXPECTED_ROW_COUNT = 13
EXPECTED_SOURCE_TUPLES = 129_998
EXPECTED_GROUP_SIZE = 7
EXPECTED_JOBS = 18_577
EXPECTED_OUTPUT_PAIRS = 18_577
EXPECTED_ANALYSIS_OUTPUTS = 18_577
EXPECTED_SIDECAR_OUTPUTS = 18_577
EXPECTED_ROOT_ARTIFACTS = 37_154
EXPECTED_SOURCE_OCCURRENCES = 18_577
EXPECTED_REQUEST_MEMORY_MB = 8_000

SPEC_KEYS = frozenset(
    {"schema", "immutable_authority", "preflight", "capacity"}
)
PREFLIGHT_SPEC_KEYS = evidence.CORRECTED_SPEC_KEYS
CAPACITY_SPEC_KEYS = evidence.CAPACITY_SPEC_KEYS
IMMUTABLE_SPEC_KEYS = evidence.IMMUTABLE_SPEC_KEYS

BINDING_KEYS = frozenset(
    {
        "schema",
        "status",
        "authority_state",
        "generated_by",
        "immutable_authority",
        "preflight",
        "capacity_evidence",
        "count_contract",
        "capacity_binding",
        "checks",
        "submission_performed",
        "jobs_rerun",
        "full_training_authority",
        "full_extraction_authority",
        "broad_production_authority",
        "canonical_promotion",
        "binding_semantic_sha256",
    }
)
COUNT_KEYS = frozenset(
    {
        "basis",
        "row_count",
        "source_tuple_count",
        "group_size",
        "expected_chunk_count",
        "expected_job_count",
        "expected_output_pair_count",
        "expected_analysis_output_count",
        "expected_sidecar_output_count",
        "expected_physical_root_artifact_count",
        "expected_source_occurrence_count",
        "source_occurrences_per_output_pair",
        "request_memory_mb",
    }
)
CAPACITY_BINDING_KEYS = frozenset(
    {
        "status",
        "reuse_scope",
        "execution_partition_sha256",
        "partition_artifact_sha256",
        "selected_rows",
        "capacity_authority_earned",
        "full_training_authority",
        "full_extraction_authority",
        "broad_production_authority",
    }
)
CHECK_KEYS = frozenset(
    {
        "current_preflight_evidence_rehashed",
        "rowwise_partition_reconstructed",
        "partition_jsonl_byte_exact",
        "aggregate_counts_exact",
        "capacity_certificate_chain_rehashed",
        "capacity_rows_terminal_single_start_no_hold",
        "capacity_resource_request_8000_mb",
        "valid_empty_and_populated_witnesses_present",
        "capacity_chunks_equal_current_partition_chunk0",
        "resolved_and_materialized_config_chain_revalidated",
        "frozen_snapshot_inventory_and_loader_revalidated",
        "manifest_runtime_identities_bound_to_immutable_authority",
        "audit_validation_authority_cross_bound",
        "output_hashes_revalidated",
        "legacy_evidence_consumed",
        "no_science_change",
        "no_tolerance_change",
        "no_submission_or_job_control",
    }
)


class BindingError(RuntimeError):
    """Malformed, inconsistent, stale, or tampered capacity evidence."""


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
        raise BindingError(f"cannot read artifact: {path}") from exc
    return digest.hexdigest()


def require_mapping(value: object, label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise BindingError(f"{label} must be an object")
    return value


def require_exact_keys(
    payload: Mapping[str, Any], expected: Iterable[str], label: str
) -> None:
    observed = set(payload)
    expected_set = set(expected)
    if observed != expected_set:
        raise BindingError(
            f"{label} key inventory differs: "
            f"missing={sorted(expected_set - observed)} "
            f"extra={sorted(observed - expected_set)}"
        )


def artifact_input(
    record: object,
    label: str,
) -> dict[str, str]:
    payload = require_mapping(record, label)
    path = payload.get("path")
    sha256 = payload.get("sha256")
    if not isinstance(path, str) or not path:
        raise BindingError(f"{label}.path must be nonempty")
    try:
        sha256 = evidence.require_sha256(f"{label}.sha256", sha256)
    except evidence.AmendmentError as exc:
        raise BindingError(str(exc)) from exc
    return {"path": path, "sha256": sha256}


def load_spec(path: Path, expected_sha256: str) -> dict[str, Any]:
    try:
        expected = evidence.require_sha256(
            "expected spec SHA-256", expected_sha256
        )
    except evidence.AmendmentError as exc:
        raise BindingError(str(exc)) from exc
    if file_sha256(path) != expected:
        raise BindingError("capacity-binding spec SHA-256 differs")
    try:
        payload = evidence.strict_json_loads(
            path.read_text(encoding="utf-8")
        )
    except evidence.AmendmentError as exc:
        raise BindingError(
            f"capacity-binding spec has ambiguous JSON: {exc}"
        ) from exc
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise BindingError("capacity-binding spec is not valid JSON") from exc
    payload = require_mapping(payload, "capacity-binding spec")
    require_exact_keys(payload, SPEC_KEYS, "capacity-binding spec")
    if payload.get("schema") != SPEC_SCHEMA:
        raise BindingError(f"capacity-binding spec schema must be {SPEC_SCHEMA}")
    for field, keys in (
        ("immutable_authority", IMMUTABLE_SPEC_KEYS),
        ("preflight", PREFLIGHT_SPEC_KEYS),
        ("capacity", CAPACITY_SPEC_KEYS),
    ):
        require_exact_keys(
            require_mapping(payload.get(field), f"spec.{field}"),
            keys,
            f"spec.{field}",
        )
    return payload


def _load_current_preflight(
    spec: Mapping[str, Any],
    immutable: Mapping[str, Any],
) -> tuple[dict[str, Any], list[dict[str, Any]], dict[str, Any]]:
    plan, plan_artifact = evidence.load_json_artifact(
        "current plan",
        spec["plan"],
        expected_schema=evidence.CORRECTED_PLAN_SCHEMA,
    )
    receipt, receipt_artifact = evidence.load_json_artifact(
        "current preflight receipt",
        spec["preflight_receipt"],
        expected_schema=evidence.CORRECTED_RECEIPT_SCHEMA,
    )
    source, source_artifact = evidence.load_source_manifest_artifact(
        "current source manifest",
        spec["source_manifest"],
        expected_schema=evidence.CORRECTED_SOURCE_SCHEMA,
    )
    try:
        source_records = resolver.validate_source_manifest(source)
    except resolver.ControllerError as exc:
        raise BindingError(
            f"current source partition reconstruction failed: {exc}"
        ) from exc
    rows, rows_artifact, _rows_bytes = evidence.load_jsonl_artifact(
        "current rows", spec["rows"]
    )
    duplicate_text, duplicate_artifact = evidence.load_text_artifact(
        "current duplicate fingerprint", spec["duplicate_fingerprint"]
    )
    partition, partition_artifact, partition_bytes = (
        evidence.load_jsonl_artifact("current partition", spec["partition"])
    )
    artifacts = {
        "plan": plan_artifact,
        "preflight_receipt": receipt_artifact,
        "source_manifest": source_artifact,
        "rows": rows_artifact,
        "duplicate_fingerprint": duplicate_artifact,
        "partition": partition_artifact,
    }
    try:
        current = evidence.validate_corrected_preflight(
            plan,
            receipt,
            source,
            source_records,
            rows,
            duplicate_text,
            partition,
            partition_bytes,
            artifacts,
        )
        evidence.validate_preflight_immutable_binding(
            "current plan", plan, receipt, immutable
        )
        evidence.validate_corrected_bundle_authority(plan, immutable)
    except evidence.AmendmentError as exc:
        raise BindingError(str(exc)) from exc
    artifact_profile = require_mapping(
        plan.get("artifact_profile"), "current plan.artifact_profile"
    )
    if artifact_profile != resolver.SIDECAR_ONLY_ARTIFACT_PROFILE:
        raise BindingError("current plan artifact profile differs")
    return current, source_records, artifact_profile


def _load_capacity(
    spec: Mapping[str, Any],
    *,
    source_records: list[dict[str, Any]],
    current: Mapping[str, Any],
    immutable: Mapping[str, Any],
    artifact_profile: Mapping[str, Any],
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    resource, resource_artifact = evidence.load_json_artifact(
        "current capacity resource certificate",
        spec["resource_certificate"],
        expected_schema=evidence.CAPACITY_SCHEMA,
    )
    root_join, root_join_artifact = evidence.load_json_artifact(
        "current capacity root/join certificate",
        spec["root_join_certificate"],
        expected_schema=evidence.ROOT_JOIN_SCHEMA,
    )
    pp_audit, pp_audit_artifact = evidence.load_json_artifact(
        "current p+p capacity audit",
        spec["pp_audit"],
        expected_schema=evidence.CAPACITY_AUDIT_SCHEMA,
    )
    auau_audit, auau_audit_artifact = evidence.load_json_artifact(
        "current Au+Au capacity audit",
        spec["auau_audit"],
        expected_schema=evidence.CAPACITY_AUDIT_SCHEMA,
    )
    artifacts = {
        "resource_certificate": resource_artifact,
        "root_join_certificate": root_join_artifact,
        "pp_audit": pp_audit_artifact,
        "auau_audit": auau_audit_artifact,
    }
    # The current sidecar-only root-health producer records the complete
    # artifact profile.  The historical evidence validator predates that
    # additive field, so validate it here against the already-rehashed plan
    # and pass only a compatibility projection to the historical validator.
    # The pinned artifact record below still retains the hash of the complete
    # unmodified certificate.
    require_exact_keys(
        root_join,
        evidence.ROOT_JOIN_KEYS | frozenset({"artifact_profile"}),
        "current capacity root/join certificate",
    )
    observed_artifact_profile = require_mapping(
        root_join.get("artifact_profile"),
        "current capacity root/join certificate.artifact_profile",
    )
    if observed_artifact_profile != artifact_profile:
        raise BindingError(
            "capacity root/join artifact profile differs from current plan"
        )
    legacy_root_join = dict(root_join)
    legacy_root_join.pop("artifact_profile")
    try:
        # Passing the current preflight in both positions deliberately binds
        # the resource certificate and the staged chunks to one current plan.
        # No historical V1 receipt is fabricated or consumed.
        return evidence.validate_capacity_evidence(
            resource,
            legacy_root_join,
            {
                "pp_background_jet8": pp_audit,
                "auau_background_jet12": auau_audit,
            },
            artifacts,
            source_records,
            current,
            current,
            immutable,
        )
    except evidence.AmendmentError as exc:
        raise BindingError(str(exc)) from exc


def _assemble_binding(spec: Mapping[str, Any]) -> dict[str, Any]:
    require_exact_keys(spec, SPEC_KEYS, "capacity-binding spec")
    if spec.get("schema") != SPEC_SCHEMA:
        raise BindingError(f"capacity-binding spec schema must be {SPEC_SCHEMA}")
    immutable_spec = require_mapping(
        spec.get("immutable_authority"), "spec.immutable_authority"
    )
    preflight_spec = require_mapping(spec.get("preflight"), "spec.preflight")
    capacity_spec = require_mapping(spec.get("capacity"), "spec.capacity")
    require_exact_keys(
        immutable_spec, IMMUTABLE_SPEC_KEYS, "spec.immutable_authority"
    )
    require_exact_keys(
        preflight_spec, PREFLIGHT_SPEC_KEYS, "spec.preflight"
    )
    require_exact_keys(capacity_spec, CAPACITY_SPEC_KEYS, "spec.capacity")

    try:
        immutable = evidence.validate_immutable_authority(immutable_spec)
    except evidence.AmendmentError as exc:
        raise BindingError(str(exc)) from exc
    current, source_records, artifact_profile = _load_current_preflight(
        preflight_spec, immutable
    )
    capacity_evidence, selected_rows = _load_capacity(
        capacity_spec,
        source_records=source_records,
        current=current,
        immutable=immutable,
        artifact_profile=artifact_profile,
    )

    count_contract = {
        "basis": (
            "ordered_disjoint_group_of_seven_partition_one_execution_and_"
            "source_occurrence_per_output_pair"
        ),
        "row_count": EXPECTED_ROW_COUNT,
        "source_tuple_count": EXPECTED_SOURCE_TUPLES,
        "group_size": EXPECTED_GROUP_SIZE,
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
    capacity_binding = {
        "status": "PASS",
        "reuse_scope": "exact_current_partition_group_size_7_capacity_only",
        "execution_partition_sha256": current[
            "execution_partition_sha256"
        ],
        "partition_artifact_sha256": current["partition"]["sha256"],
        "selected_rows": selected_rows,
        "capacity_authority_earned": True,
        "full_training_authority": 0,
        "full_extraction_authority": False,
        "broad_production_authority": False,
    }
    payload: dict[str, Any] = {
        "schema": BINDING_SCHEMA,
        "status": "PASS",
        "authority_state": AUTHORITY_STATE,
        "generated_by": GENERATED_BY,
        "immutable_authority": immutable,
        "preflight": current,
        "capacity_evidence": capacity_evidence,
        "count_contract": count_contract,
        "capacity_binding": capacity_binding,
        "checks": {
            "current_preflight_evidence_rehashed": True,
            "rowwise_partition_reconstructed": True,
            "partition_jsonl_byte_exact": True,
            "aggregate_counts_exact": True,
            "capacity_certificate_chain_rehashed": True,
            "capacity_rows_terminal_single_start_no_hold": True,
            "capacity_resource_request_8000_mb": True,
            "valid_empty_and_populated_witnesses_present": True,
            "capacity_chunks_equal_current_partition_chunk0": True,
            "resolved_and_materialized_config_chain_revalidated": True,
            "frozen_snapshot_inventory_and_loader_revalidated": True,
            "manifest_runtime_identities_bound_to_immutable_authority": True,
            "audit_validation_authority_cross_bound": True,
            "output_hashes_revalidated": True,
            "legacy_evidence_consumed": False,
            "no_science_change": True,
            "no_tolerance_change": True,
            "no_submission_or_job_control": True,
        },
        "submission_performed": False,
        "jobs_rerun": 0,
        "full_training_authority": 0,
        "full_extraction_authority": False,
        "broad_production_authority": False,
        "canonical_promotion": False,
    }
    payload["binding_semantic_sha256"] = semantic_sha256(payload)
    return payload


def reconstruct_spec(payload: Mapping[str, Any]) -> dict[str, Any]:
    immutable = require_mapping(
        payload.get("immutable_authority"), "binding.immutable_authority"
    )
    preflight = require_mapping(payload.get("preflight"), "binding.preflight")
    capacity = require_mapping(
        payload.get("capacity_evidence"), "binding.capacity_evidence"
    )
    return {
        "schema": SPEC_SCHEMA,
        "immutable_authority": {
            "bundle_manifest": artifact_input(
                immutable.get("bundle_manifest"),
                "binding immutable bundle manifest",
            ),
            "materialization_receipt": artifact_input(
                immutable.get("materialization_receipt"),
                "binding immutable materialization receipt",
            ),
        },
        "preflight": {
            "plan": artifact_input(preflight.get("plan"), "binding plan"),
            "preflight_receipt": artifact_input(
                preflight.get("preflight_receipt"),
                "binding preflight receipt",
            ),
            "source_manifest": artifact_input(
                preflight.get("source_manifest"),
                "binding source manifest",
            ),
            "rows": artifact_input(preflight.get("rows"), "binding rows"),
            "duplicate_fingerprint": artifact_input(
                preflight.get("duplicate_fingerprint"),
                "binding duplicate fingerprint",
            ),
            "partition": artifact_input(
                preflight.get("partition"), "binding partition"
            ),
        },
        "capacity": {
            "resource_certificate": artifact_input(
                capacity.get("resource_certificate"),
                "binding capacity resource certificate",
            ),
            "root_join_certificate": artifact_input(
                capacity.get("root_join_certificate"),
                "binding capacity root/join certificate",
            ),
            "pp_audit": artifact_input(
                capacity.get("pp_audit"), "binding p+p capacity audit"
            ),
            "auau_audit": artifact_input(
                capacity.get("auau_audit"), "binding Au+Au capacity audit"
            ),
        },
    }


def _validate_surface(payload: Mapping[str, Any]) -> None:
    require_exact_keys(payload, BINDING_KEYS, "capacity binding")
    if (
        payload.get("schema") != BINDING_SCHEMA
        or payload.get("status") != "PASS"
        or payload.get("authority_state") != AUTHORITY_STATE
        or payload.get("generated_by") != GENERATED_BY
        or payload.get("submission_performed") is not False
        or payload.get("jobs_rerun") != 0
        or payload.get("full_training_authority") != 0
        or payload.get("full_extraction_authority") is not False
        or payload.get("broad_production_authority") is not False
        or payload.get("canonical_promotion") is not False
    ):
        raise BindingError("capacity-binding authority surface differs")
    counts = require_mapping(
        payload.get("count_contract"), "binding.count_contract"
    )
    require_exact_keys(counts, COUNT_KEYS, "binding.count_contract")
    exact_counts = {
        "row_count": EXPECTED_ROW_COUNT,
        "source_tuple_count": EXPECTED_SOURCE_TUPLES,
        "group_size": EXPECTED_GROUP_SIZE,
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
    if any(counts.get(key) != value for key, value in exact_counts.items()):
        raise BindingError("capacity-binding exact count contract differs")
    binding = require_mapping(
        payload.get("capacity_binding"), "binding.capacity_binding"
    )
    require_exact_keys(
        binding, CAPACITY_BINDING_KEYS, "binding.capacity_binding"
    )
    if (
        binding.get("status") != "PASS"
        or binding.get("reuse_scope")
        != "exact_current_partition_group_size_7_capacity_only"
        or binding.get("capacity_authority_earned") is not True
        or binding.get("full_training_authority") != 0
        or binding.get("full_extraction_authority") is not False
        or binding.get("broad_production_authority") is not False
    ):
        raise BindingError("capacity-binding scope or authority differs")
    checks = require_mapping(payload.get("checks"), "binding.checks")
    require_exact_keys(checks, CHECK_KEYS, "binding.checks")
    expected_checks = {key: True for key in CHECK_KEYS}
    expected_checks["legacy_evidence_consumed"] = False
    if checks != expected_checks:
        raise BindingError("capacity-binding check surface differs")
    observed_semantic = payload.get("binding_semantic_sha256")
    try:
        observed_semantic = evidence.require_sha256(
            "binding semantic SHA-256", observed_semantic
        )
    except evidence.AmendmentError as exc:
        raise BindingError(str(exc)) from exc
    semantic_payload = dict(payload)
    semantic_payload.pop("binding_semantic_sha256")
    if observed_semantic != semantic_sha256(semantic_payload):
        raise BindingError("capacity-binding semantic SHA-256 differs")


def validate_binding_payload(payload: Mapping[str, Any]) -> dict[str, Any]:
    """Rehash every pinned input and prove a byte-exact deterministic rebuild."""

    _validate_surface(payload)
    rebuilt = _assemble_binding(reconstruct_spec(payload))
    if canonical_json_bytes(rebuilt) != canonical_json_bytes(payload):
        raise BindingError(
            "capacity binding differs from current pinned evidence rebuild"
        )
    return dict(payload)


def build_binding(spec: Mapping[str, Any]) -> dict[str, Any]:
    payload = _assemble_binding(spec)
    _validate_surface(payload)
    return payload


def write_once_readonly_json(
    path: Path,
    payload: Mapping[str, Any],
    *,
    allowed_output_root: Path,
) -> str:
    try:
        return evidence.write_once_readonly_json(
            path,
            payload,
            allowed_output_root=allowed_output_root,
        )
    except evidence.AmendmentError as exc:
        raise BindingError(str(exc)) from exc


def load_binding(path: Path, expected_sha256: str) -> dict[str, Any]:
    try:
        payload, _artifact = evidence.load_json_artifact(
            "capacity partition binding",
            {"path": str(path), "sha256": expected_sha256},
            expected_schema=BINDING_SCHEMA,
        )
    except evidence.AmendmentError as exc:
        raise BindingError(str(exc)) from exc
    return validate_binding_payload(payload)


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="action", required=True)
    build = subparsers.add_parser("build", help="build the immutable binding")
    build.add_argument("--spec", type=Path, required=True)
    build.add_argument("--expected-spec-sha256", required=True)
    build.add_argument("--output", type=Path, required=True)
    build.add_argument("--allowed-output-root", type=Path, required=True)
    verify = subparsers.add_parser(
        "verify", help="rehash all evidence and byte-compare the binding"
    )
    verify.add_argument("--spec", type=Path, required=True)
    verify.add_argument("--expected-spec-sha256", required=True)
    verify.add_argument("--binding", type=Path, required=True)
    verify.add_argument("--expected-binding-sha256", required=True)
    verify.add_argument("--output-json", type=Path)
    verify.add_argument("--allowed-output-root", type=Path)
    return parser.parse_args(argv)


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        spec = load_spec(args.spec, args.expected_spec_sha256)
        rebuilt = build_binding(spec)
        if args.action == "build":
            disposition = write_once_readonly_json(
                args.output,
                rebuilt,
                allowed_output_root=args.allowed_output_root,
            )
            summary = {
                "schema": BINDING_SCHEMA,
                "status": "PASS",
                "authority_state": AUTHORITY_STATE,
                "path": str(args.output),
                "sha256": file_sha256(args.output),
                "binding_semantic_sha256": rebuilt[
                    "binding_semantic_sha256"
                ],
                "write_disposition": disposition,
                "submission_performed": False,
                "full_training_authority": 0,
                "full_extraction_authority": False,
            }
        else:
            observed = load_binding(
                args.binding, args.expected_binding_sha256
            )
            if canonical_json_bytes(observed) != canonical_json_bytes(rebuilt):
                raise BindingError(
                    "binding readback differs from current pinned evidence"
                )
            summary = {
                "schema": READBACK_SCHEMA,
                "status": "PASS",
                "authority_state": AUTHORITY_STATE,
                "binding": str(args.binding),
                "binding_sha256": args.expected_binding_sha256,
                "binding_semantic_sha256": observed[
                    "binding_semantic_sha256"
                ],
                "all_evidence_rehashed": True,
                "byte_exact_rebuild": True,
                "legacy_evidence_consumed": False,
                "submission_performed": False,
                "full_training_authority": 0,
                "full_extraction_authority": False,
            }
            if args.output_json is not None:
                if args.allowed_output_root is None:
                    raise BindingError(
                        "--allowed-output-root is required with --output-json"
                    )
                write_once_readonly_json(
                    args.output_json,
                    summary,
                    allowed_output_root=args.allowed_output_root,
                )
        print(json.dumps(summary, sort_keys=True))
        return 0
    except (
        BindingError,
        evidence.AmendmentError,
        resolver.ControllerError,
        source_builder.ManifestError,
    ) as exc:
        print(
            json.dumps(
                {
                    "schema": BINDING_SCHEMA,
                    "status": "FAIL",
                    "authority_state": AUTHORITY_STATE,
                    "error": str(exc),
                    "submission_performed": False,
                    "full_training_authority": 0,
                    "full_extraction_authority": False,
                },
                sort_keys=True,
            ),
            file=sys.stderr,
        )
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
