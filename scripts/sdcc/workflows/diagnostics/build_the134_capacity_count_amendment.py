#!/usr/bin/env python3
"""Bind the passed THE-134 group-of-seven canary to corrected V2 counts.

The original non-submitting V1 extraction plan correctly projected 18,577
group-of-seven jobs, but its source-occurrence field retained the unpartitioned
129,998-tuple count.  The V2 resolver replaces that ambiguous alias with an
explicit, ordered, disjoint partition:

    129,998 source tuples
      -> 18,577 chunks/jobs/output pairs/source occurrences
      -> 37,154 physical ROOT artifacts

This tool creates a deterministic mechanical amendment.  It rehashes and
cross-validates the preserved V1 plan/source/receipt, corrected V2
plan/source/receipt/partition, and terminal capacity certificate chain.  For
the two capacity rows it also proves that the actually staged seven-tuple file
is byte-for-byte the first V2 execution chunk for that source.

The amendment does not submit, copy, remove, release, requeue, merge, or
promote anything.  It grants no training, extraction, broad-production, or
CANONICAL authority.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import os
import re
import stat
import sys
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


HERE = Path(__file__).resolve().parent
SOURCE_BUILDER_PATH = HERE / "build_the134_source_authority_manifest.py"
RESOLVER_PATH = HERE / "resolve_the134_full_multiview_extraction.py"


def _load_local_module(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load local module: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


source_builder = _load_local_module("the134_source_builder", SOURCE_BUILDER_PATH)
resolver = _load_local_module("the134_full_resolver", RESOLVER_PATH)


SPEC_SCHEMA = "THE134_GROUP7_CAPACITY_COUNT_AMENDMENT_SPEC_V1"
AMENDMENT_SCHEMA = "THE134_GROUP7_CAPACITY_COUNT_AMENDMENT_V1"
READBACK_SCHEMA = "THE134_GROUP7_CAPACITY_COUNT_AMENDMENT_READBACK_V1"
AUTHORITY_STATE = "MECHANICAL_COUNT_SEMANTICS_AMENDED_NOT_EXTRACTION_AUTHORITY"
GENERATED_BY = "build_the134_capacity_count_amendment.py"

LEGACY_PLAN_SCHEMA = "THE134_FULL_MULTIVIEW_EXTRACTION_PLAN_V1"
LEGACY_RECEIPT_SCHEMA = "THE134_FULL_MULTIVIEW_EXTRACTION_PREFLIGHT_RECEIPT_V1"
LEGACY_PARTITION_SCHEMA = "THE134_FULL_EXTRACTION_PARTITION_CONTRACT_V1"
CORRECTED_PLAN_SCHEMA = resolver.PLAN_SCHEMA
CORRECTED_RECEIPT_SCHEMA = resolver.RECEIPT_SCHEMA
CORRECTED_SOURCE_SCHEMA = resolver.SOURCE_SCHEMA
CORRECTED_PARTITION_SCHEMA = resolver.PARTITION_SCHEMA
CAPACITY_SCHEMA = "THE134_GROUP7_PARTITION_CAPACITY_CERTIFICATE_V1"
ROOT_JOIN_SCHEMA = "THE134_SMOKE_ROOT_HEALTH_IDENTITY_JOIN_V1"
CAPACITY_AUDIT_SCHEMA = "THE134_CAPACITY_MULTIVIEW_MATRIX_AUDIT_V1"

GROUP_SIZE = 7
REQUEST_MEMORY_MB = 8000
EXPECTED_SOURCE_TUPLES = 129_998
EXPECTED_CHUNKS = 18_577
EXPECTED_ROOT_ARTIFACTS = 37_154
LEGACY_OVERSTATEMENT = EXPECTED_SOURCE_TUPLES - EXPECTED_CHUNKS
SELECTED_CAPACITY_ROWS = (
    "pp_background_jet8",
    "auau_background_jet12",
)
EXPECTED_POPULATION_STATE = {
    "pp_background_jet8": "VALID_EMPTY",
    "auau_background_jet12": "POPULATED",
}
EXPECTED_SYSTEM = {
    "pp_background_jet8": "pp",
    "auau_background_jet12": "auau",
}

SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
GIT_COMMIT_RE = re.compile(r"^[0-9a-f]{40}$")
CLUSTER_PROC_RE = re.compile(r"^[1-9][0-9]*\.[0-9]+$")

ARTIFACT_REF_KEYS = frozenset({"path", "sha256"})
SPEC_KEYS = frozenset(
    {"schema", "immutable_authority", "legacy", "corrected", "capacity"}
)
IMMUTABLE_SPEC_KEYS = frozenset(
    {"bundle_manifest", "materialization_receipt"}
)
LEGACY_SPEC_KEYS = frozenset(
    {"plan", "preflight_receipt", "source_manifest"}
)
CORRECTED_SPEC_KEYS = frozenset(
    {
        "plan",
        "preflight_receipt",
        "source_manifest",
        "rows",
        "duplicate_fingerprint",
        "partition",
    }
)
CAPACITY_SPEC_KEYS = frozenset(
    {
        "resource_certificate",
        "root_join_certificate",
        "pp_audit",
        "auau_audit",
    }
)
COUNT_CORRECTION_KEYS = frozenset(
    {
        "basis",
        "group_size",
        "source_tuple_count",
        "legacy_expected_occurrence_count",
        "corrected_chunk_count",
        "corrected_job_count",
        "corrected_output_pair_count",
        "corrected_analysis_output_count",
        "corrected_sidecar_output_count",
        "physical_root_artifact_count",
        "corrected_source_occurrence_count",
        "source_occurrences_per_output_pair",
        "legacy_occurrence_overstatement",
        "science_semantics_changed",
        "tolerances_changed",
    }
)
CAPACITY_BINDING_KEYS = frozenset(
    {
        "status",
        "reuse_scope",
        "corrected_execution_partition_sha256",
        "corrected_partition_artifact_sha256",
        "selected_rows",
        "capacity_authority_earned",
        "full_training_authority",
        "full_extraction_authority",
        "broad_production_authority",
    }
)
CHECK_KEYS = frozenset(
    {
        "legacy_v1_evidence_rehashed",
        "legacy_tuple_aliases_preserved_as_history_only",
        "corrected_v2_evidence_rehashed",
        "source_semantics_identical_v1_to_v2",
        "rowwise_partition_reconstructed",
        "partition_jsonl_byte_exact",
        "aggregate_counts_exact",
        "capacity_certificate_chain_rehashed",
        "capacity_rows_terminal_single_start_no_hold",
        "capacity_resource_request_8000_mb",
        "valid_empty_and_populated_witnesses_present",
        "capacity_chunks_equal_corrected_chunk0",
        "resolved_and_materialized_config_chain_revalidated",
        "frozen_snapshot_inventory_and_loader_revalidated",
        "manifest_runtime_identities_bound_to_immutable_authority",
        "audit_validation_authority_cross_bound",
        "output_hashes_revalidated",
        "no_science_change",
        "no_tolerance_change",
        "no_submission_or_job_control",
    }
)

LEGACY_PLAN_KEYS = frozenset(
    {
        "schema",
        "status",
        "execution_state",
        "submission_performed",
        "campaign",
        "training_period_si_contract",
        "training_period_si_contract_sha256",
        "authority",
        "closure_witness_boundary",
        "source_family_closure",
        "execution_partition",
        "input_manifests",
        "duplicate_contract",
        "duplicate_fingerprint_sha256",
        "execution_fingerprint_sha256",
        "rows",
    }
)
LEGACY_RECEIPT_KEYS = frozenset(
    {
        "schema",
        "status",
        "submission_performed",
        "row_count",
        "requested_scope",
        "full_training_authority",
        "authority_state",
        "training_period_si_contract_sha256",
        "closure_witness_boundary_sha256",
        "execution_partition_sha256",
        "bundle_manifest_sha256",
        "materialization_receipt_sha256",
        "source_manifest_sha256",
        "duplicate_fingerprint_sha256",
        "execution_fingerprint_sha256",
        "artifacts",
    }
)
CORRECTED_PLAN_KEYS = LEGACY_PLAN_KEYS
CORRECTED_RECEIPT_KEYS = LEGACY_RECEIPT_KEYS
LEGACY_RECEIPT_ARTIFACT_KEYS = frozenset(
    {"plan", "rows", "duplicate_fingerprint"}
)
CORRECTED_RECEIPT_ARTIFACT_KEYS = frozenset(
    {"plan", "rows", "duplicate_fingerprint", "partition"}
)
RECEIPT_ARTIFACT_KEYS = frozenset({"name", "sha256"})
PARTITION_RECEIPT_ARTIFACT_KEYS = frozenset(
    {"name", "sha256", "record_count"}
)
SOURCE_PROVENANCE_KEYS = frozenset({"schema", "inputs"})
SOURCE_PROVENANCE_RECORD_KEYS = frozenset(
    {
        "code_sha256",
        "config_sha256",
        "input_file_sha256",
        "input_uri_sha256",
        "path",
        "row_id",
        "source_manifest_sha256",
        "source_sample",
        "system",
    }
)
AUDIT_SOURCE_PROVENANCE_KEYS = frozenset({"path", "record", "sha256"})
SOURCE_CONTRACT_KEYS = frozenset(
    {
        "dataset",
        "input_file_sha256",
        "input_uri_hash",
        "lane",
        "ownership_state",
        "period",
        "run",
        "sample",
        "segment",
        "si_di_role",
        "source_manifest_sha256",
    }
)
SOURCE_EXECUTION_CONTRACT_KEYS = frozenset(
    {
        "args_file",
        "chunk_index",
        "run",
        "staged_chunk_list",
        "submitted_args_sha256",
    }
)
RESOURCE_ROW_KEYS = frozenset(
    {
        "row_id",
        "cluster_proc",
        "job_status",
        "exit_code",
        "num_holds",
        "num_job_starts",
        "request_memory_mb",
        "memory_usage_mb",
        "memory_usage_resolution",
        "resident_set_size_kb",
        "remote_wall_clock_seconds",
    }
)
IMMUTABLE_AUTHORITY_OUTPUT_KEYS = frozenset(
    {
        "status",
        "bundle_manifest",
        "materialization_receipt",
        "bundle_identity_sha256",
        "bundle_semantic_fingerprint_sha256",
        "public_commit",
        "code_sha256",
        "replay_schema_sha256",
        "training_schema_sha256",
        "semantic_sha256",
        "runtime",
        "bundle_artifacts",
        "materialized_bundle_path",
        "artifact_count",
        "total_artifact_bytes",
        "readback",
    }
)
IMMUTABLE_RUNTIME_OUTPUT_KEYS = frozenset(
    {
        "release",
        "offline_main",
        "calo_reco_soname",
        "request_memory_mb",
        "release_core_lib_dir",
        "release_core_lib64_dir",
    }
)
BUNDLE_ARTIFACT_OUTPUT_KEYS = frozenset(
    {"role", "path", "resolved_path", "sha256", "size_bytes"}
)
INPUT_MANIFEST_KEYS = frozenset({"materialization", "bundle", "sources"})
BUNDLE_INPUT_KEYS = frozenset(
    {"path", "sha256", "semantic_fingerprint_sha256"}
)
MATERIALIZATION_INPUT_KEYS = frozenset(
    {
        "path",
        "sha256",
        "bundle_identity_sha256",
        "digest_named_bundle_path",
        "readback",
    }
)
SOURCE_INPUT_KEYS = frozenset({"path", "sha256"})
PINNED_ARTIFACT_OUTPUT_KEYS = frozenset(
    {"label", "path", "resolved_path", "sha256", "size_bytes"}
)
PINNED_SCHEMA_ARTIFACT_OUTPUT_KEYS = frozenset(
    {*PINNED_ARTIFACT_OUTPUT_KEYS, "schema"}
)
LEGACY_PREFLIGHT_OUTPUT_KEYS = frozenset(
    {
        "status",
        "plan",
        "preflight_receipt",
        "source_manifest",
        "submission_performed",
        "full_training_authority",
        "execution_partition_sha256",
        "bundle_manifest_sha256",
        "materialization_receipt_sha256",
        "source_manifest_sha256",
        "source_tuple_count",
        "expected_job_count",
        "expected_output_count",
        "legacy_expected_occurrence_count",
        "rows",
    }
)
CORRECTED_PREFLIGHT_OUTPUT_KEYS = frozenset(
    {
        "status",
        "plan",
        "preflight_receipt",
        "source_manifest",
        "rows",
        "duplicate_fingerprint",
        "partition",
        "submission_performed",
        "full_training_authority",
        "bundle_manifest_sha256",
        "materialization_receipt_sha256",
        "source_manifest_sha256",
        "duplicate_fingerprint_sha256",
        "execution_partition_sha256",
        "execution_fingerprint_sha256",
        "source_tuple_count",
        "expected_chunk_count",
        "expected_job_count",
        "expected_output_pair_count",
        "expected_analysis_output_count",
        "expected_sidecar_output_count",
        "expected_physical_root_artifact_count",
        "expected_source_occurrence_count",
        "source_occurrences_per_output_pair",
        "rows_proof",
        "capacity_source_semantics",
        "bundle_contract_semantic_sha256_by_system",
        "science_contract_semantic_sha256",
    }
)
SOURCE_CROSSWALK_OUTPUT_KEYS = frozenset(
    {
        "status",
        "row_count",
        "source_tuple_count",
        "global_unique_tuple_count",
        "source_rows_semantic_sha256",
        "manifest_semantic_sha256",
        "global_tuple_identity_records_sha256",
        "rows",
    }
)
CAPACITY_EVIDENCE_OUTPUT_KEYS = frozenset(
    {
        "resource_certificate",
        "root_join_certificate",
        "pp_audit",
        "auau_audit",
        "validation_authority",
    }
)
VALIDATION_AUTHORITY_OUTPUT_KEYS = frozenset(
    {
        "controller",
        "runtime_authority",
        "submission_journal",
        "submission_manifest",
        "submission_receipt",
    }
)
ROW_PROOF_OUTPUT_KEYS = frozenset(
    {
        "row_id",
        "source_tuple_count",
        "expected_chunk_count",
        "expected_job_count",
        "expected_output_pair_count",
        "expected_source_occurrence_count",
        "tail_tuple_count",
        "row_partition_sha256",
        "chunk_records_sha256",
        "first_execution_chunk_sha256",
        "row_fingerprint_sha256",
    }
)
CAPACITY_SOURCE_SEMANTIC_OUTPUT_KEYS = frozenset(
    {
        "row_id",
        "system",
        "lane",
        "dataset",
        "sample",
        "source_role",
        "minimum_bias_gate",
        "period",
        "run",
        "si_di_role",
        "ownership_state",
        "code_sha256",
        "replay_schema_sha256",
        "training_schema_sha256",
        "semantic_sha256",
        "library_sha256",
        "model_sha256",
        "base_config_sha256",
    }
)
CAPACITY_ROW_OUTPUT_KEYS = frozenset(
    {
        "row_id",
        "system",
        "population_state",
        "cluster_proc",
        "exit_code",
        "num_job_starts",
        "num_holds",
        "request_memory_mb",
        "memory_usage_mb",
        "resident_set_size_kb",
        "remote_wall_clock_seconds",
        "capacity_chunk_index",
        "corrected_partition_chunk_index",
        "segment",
        "tuple_count",
        "staged_chunk_list",
        "staged_chunk_sha256",
        "corrected_chunk_fingerprint_sha256",
        "corrected_execution_chunk_sha256",
        "full_source_manifest_sha256",
        "source_contract",
        "source_identity_canonical_sha256",
        "source_occurrence_id_hex",
        "base_config_sha256",
        "capacity_runtime_config",
        "capacity_materialized_config",
        "capacity_fanout_contract",
        "runtime_config_chain_semantic_sha256",
        "snapshot_loader_receipt",
        "snapshot_manifest",
        "capacity_submission_manifest_row_sha256",
        "analysis_output",
        "sidecar_output",
        "source_provenance",
        "staged_chunk_equals_corrected_partition_chunk0",
    }
)
OUTPUT_ARTIFACT_KEYS = frozenset({"path", "sha256", "size_bytes"})
LEGACY_ROW_OUTPUT_KEYS = frozenset(
    {
        "row_id",
        "source_tuple_count",
        "expected_chunk_count",
        "legacy_expected_occurrence_count",
    }
)
SOURCE_CROSSWALK_ROW_KEYS = frozenset(
    {
        "row_id",
        "tuple_count",
        "tuple_records_sha256",
        "full_source_manifest_sha256",
        "first_tuple_sha256",
        "last_tuple_sha256",
        "input_lists",
        "legacy_tuple_aliases_preserved_as_history_only",
        "corrected_execution_aliases_absent",
    }
)
CROSSWALK_INPUT_LIST_KEYS = frozenset(
    {
        "role",
        "path",
        "resolved_path",
        "sha256",
        "size_bytes",
        "line_count",
        "executable_count",
    }
)
CAPACITY_KEYS = frozenset(
    {
        "schema",
        "status",
        "capacity_canary_id",
        "full_plan",
        "full_plan_sha256",
        "preflight_receipt",
        "preflight_receipt_sha256",
        "bundle_manifest_sha256",
        "materialization_receipt_sha256",
        "execution_partition_sha256",
        "selected_rows",
        "rows",
        "root_health_identity_join_certificate",
        "root_health_identity_join_certificate_sha256",
        "multiview_audits",
        "execution_group_size",
        "source_occurrences_per_output",
        "capacity_authority_earned",
        "full_training_authority",
        "submission_performed",
        "validation_authority",
    }
)
ROOT_JOIN_KEYS = frozenset(
    {
        "schema",
        "status",
        "scope",
        "analysis_health_profile",
        "sidecar_health_profile",
        "row_count",
        "populated_row_count",
        "valid_empty_row_count",
        "populated_rows",
        "valid_empty_rows",
        "rows",
        "failures",
        "execution_group_size",
        "source_occurrences_per_output",
        "capacity_authority_earned",
        "full_training_authority",
        "validation_authority",
    }
)
AUDIT_KEYS = frozenset(
    {
        "schema",
        "status",
        "scope",
        "row_id",
        "system",
        "selected_rows",
        "authority_state",
        "execution_group_size",
        "source_occurrences_per_output",
        "source_occurrence_count",
        "source_occurrence_id_hex",
        "source_contract",
        "source_execution_contract",
        "source_identity_canonical_sha256",
        "source_provenance",
        "submission_manifest",
        "submission_receipt",
        "root_health_identity_join_certificate",
        "sidecar",
        "file_report",
        "file_validator",
        "checks",
        "population_state",
        "class_balance_state",
        "matrix_materialized",
        "matrix_out",
        "full_training_authority",
        "training_matrix_authority_earned",
        "failures",
    }
)
PLAN_AUTHORITY_KEYS = frozenset(
    {
        "requested_scope",
        "full_training_authority",
        "authority_state",
    }
)
CAPACITY_AUDIT_CHECK_KEYS = frozenset(
    {
        "exact_branch_inventory",
        "exact_manifest_receipt_source_binding",
        "file_level_semantics",
        "matrix_absent",
        "root_identity_join_binding",
        "sidecar_content_hash_binding",
    }
)


class AmendmentError(ValueError):
    """Fail-closed count-amendment contract violation."""


def canonical_json_bytes(payload: Any) -> bytes:
    return (
        json.dumps(
            payload,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        )
        + "\n"
    ).encode("utf-8")


def semantic_sha256(payload: Any) -> str:
    return hashlib.sha256(
        json.dumps(
            payload,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        ).encode("utf-8")
    ).hexdigest()


def strict_json_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    """Reject duplicate JSON keys instead of silently accepting the last."""

    payload: dict[str, Any] = {}
    for key, value in pairs:
        if key in payload:
            raise AmendmentError(f"JSON object duplicates key: {key}")
        payload[key] = value
    return payload


def strict_json_loads(text: str) -> Any:
    def reject_constant(value: str) -> None:
        raise AmendmentError(f"JSON contains non-finite constant: {value}")

    return json.loads(
        text,
        object_pairs_hook=strict_json_object,
        parse_constant=reject_constant,
    )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require_exact_keys(
    label: str,
    payload: Mapping[str, Any],
    keys: set[str] | frozenset[str],
) -> None:
    observed = set(payload)
    expected = set(keys)
    if observed != expected:
        raise AmendmentError(
            f"{label} field inventory differs: "
            f"missing={sorted(expected - observed)} "
            f"extra={sorted(observed - expected)}"
        )


def require_mapping(label: str, value: object) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise AmendmentError(f"{label} must be one object")
    return dict(value)


def require_sequence(label: str, value: object) -> list[Any]:
    if not isinstance(value, list):
        raise AmendmentError(f"{label} must be one list")
    return list(value)


def require_sha256(label: str, value: object) -> str:
    text = str(value)
    if not SHA256_RE.fullmatch(text):
        raise AmendmentError(f"{label} must be a lowercase 64-character SHA-256")
    return text


def require_bool(label: str, value: object, expected: bool) -> bool:
    if value is not expected:
        raise AmendmentError(f"{label} must be {expected}")
    return expected


def require_exact_int(label: str, value: object, expected: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value != expected:
        raise AmendmentError(f"{label} must be exactly {expected}")
    return expected


def require_nonnegative_int(label: str, value: object) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise AmendmentError(f"{label} must be a nonnegative integer")
    return value


def require_nonempty_text(label: str, value: object) -> str:
    text = str(value)
    if not text or any(character in text for character in "\x00\n\r\t"):
        raise AmendmentError(f"{label} is empty or contains a control character")
    return text


def require_decimal_int_text(label: str, value: object) -> int:
    text = require_nonempty_text(label, value)
    if not re.fullmatch(r"(?:0|[1-9][0-9]*)", text):
        raise AmendmentError(f"{label} must be a canonical nonnegative integer")
    return int(text)


def require_absolute_file(label: str, value: object) -> Path:
    path = Path(require_nonempty_text(label, value)).expanduser()
    if not path.is_absolute() or not path.is_file():
        raise AmendmentError(f"{label} must be an existing absolute file: {path}")
    if path.is_symlink():
        raise AmendmentError(f"{label} must not be a symlink: {path}")
    return path


def require_runtime_file(label: str, value: object) -> Path:
    """Require an absolute runtime file while permitting declared symlinks."""

    path = Path(require_nonempty_text(label, value)).expanduser()
    if not path.is_absolute() or not path.is_file():
        raise AmendmentError(f"{label} must be an existing absolute file: {path}")
    return path


def require_runtime_directory(label: str, value: object) -> Path:
    """Require an absolute runtime directory while permitting site aliases."""

    path = Path(require_nonempty_text(label, value)).expanduser()
    if not path.is_absolute() or not path.is_dir():
        raise AmendmentError(
            f"{label} must be an existing absolute directory: {path}"
        )
    return path


def normalize_runtime_artifact_ref(
    label: str, path_value: object, sha256_value: object
) -> dict[str, Any]:
    """Rehash a runtime artifact, including a declared CVMFS symlink."""

    path = require_runtime_file(f"{label}.path", path_value)
    expected = require_sha256(f"{label}.sha256", sha256_value)
    observed = sha256_file(path)
    if observed != expected:
        raise AmendmentError(
            f"{label} hash drift: expected={expected} observed={observed} "
            f"path={path}"
        )
    return {
        "path": str(path),
        "resolved_path": str(path.resolve(strict=True)),
        "sha256": observed,
        "size_bytes": path.stat().st_size,
    }


def require_same_resolved_file(
    label: str, left: object, right: object
) -> tuple[Path, Path]:
    """Require two absolute aliases to identify the same existing file."""

    left_path = require_absolute_file(f"{label}.left", left)
    right_path = require_absolute_file(f"{label}.right", right)
    if not os.path.samefile(left_path, right_path):
        raise AmendmentError(
            f"{label} aliases identify different files: "
            f"left={left_path} right={right_path}"
        )
    return left_path, right_path


def require_same_runtime_file(
    label: str, left: object, right: object
) -> tuple[Path, Path]:
    """Alias-safe file identity for runtime paths that may be symlinks."""

    left_path = require_runtime_file(f"{label}.left", left)
    right_path = require_runtime_file(f"{label}.right", right)
    if not os.path.samefile(left_path, right_path):
        raise AmendmentError(
            f"{label} aliases identify different files: "
            f"left={left_path} right={right_path}"
        )
    return left_path, right_path


def require_same_runtime_directory(
    label: str, left: object, right: object
) -> tuple[Path, Path]:
    """Alias-safe directory identity for CVMFS and site-path aliases."""

    left_path = require_runtime_directory(f"{label}.left", left)
    right_path = require_runtime_directory(f"{label}.right", right)
    if not os.path.samefile(left_path, right_path):
        raise AmendmentError(
            f"{label} aliases identify different directories: "
            f"left={left_path} right={right_path}"
        )
    return left_path, right_path


def reject_symlinked_existing_ancestors(label: str, path: Path) -> None:
    """Reject write redirection while allowing normal unresolved new suffixes."""

    current = Path(path.anchor)
    for part in path.parts[1:]:
        current = current / part
        if current.is_symlink():
            raise AmendmentError(
                f"{label} traverses a symlinked path component: {current}"
            )
        if not current.exists():
            break


def validate_output_parent(
    label: str,
    path: Path,
    *,
    allowed_output_root: Path | None,
) -> None:
    """Protect writes while permitting a declared site alias root."""

    if ".." in path.parts:
        raise AmendmentError(f"{label} must not contain parent traversal")
    if allowed_output_root is None:
        reject_symlinked_existing_ancestors(label, path)
        return
    root = allowed_output_root.expanduser()
    if (
        not root.is_absolute()
        or not root.exists()
        or not root.is_dir()
    ):
        raise AmendmentError(
            f"{label} allowed output root must be an existing absolute "
            f"directory: {root}"
        )
    lexical_root = Path(os.path.abspath(root))
    lexical_path = Path(os.path.abspath(path))
    try:
        lexical_path.relative_to(lexical_root)
    except ValueError as exc:
        raise AmendmentError(
            f"{label} escapes its declared output root: "
            f"path={path} root={root}"
        ) from exc
    resolved_root = root.resolve(strict=True)
    existing = lexical_path
    while not existing.exists():
        if existing == existing.parent:
            raise AmendmentError(
                f"{label} has no existing ancestor: {path}"
            )
        existing = existing.parent
    resolved_existing = existing.resolve(strict=True)
    try:
        resolved_existing.relative_to(resolved_root)
    except ValueError as exc:
        raise AmendmentError(
            f"{label} resolves outside its declared output root: "
            f"path={path} root={root} resolved={resolved_existing}"
        ) from exc


def normalize_artifact_ref(label: str, raw: object) -> dict[str, Any]:
    payload = require_mapping(label, raw)
    require_exact_keys(label, payload, ARTIFACT_REF_KEYS)
    path = require_absolute_file(f"{label}.path", payload["path"])
    expected = require_sha256(f"{label}.sha256", payload["sha256"])
    observed = sha256_file(path)
    if observed != expected:
        raise AmendmentError(
            f"{label} hash drift: expected={expected} observed={observed} "
            f"path={path}"
        )
    return {
        "path": str(path),
        "resolved_path": str(path.resolve(strict=True)),
        "sha256": observed,
        "size_bytes": path.stat().st_size,
    }


def load_json_artifact(
    label: str,
    raw: object,
    *,
    expected_schema: str | None = None,
) -> tuple[dict[str, Any], dict[str, Any]]:
    artifact = normalize_artifact_ref(label, raw)
    path = Path(artifact["path"])
    try:
        payload = strict_json_loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise AmendmentError(f"{label} is not valid JSON") from exc
    if not isinstance(payload, dict):
        raise AmendmentError(f"{label} must contain one JSON object")
    if expected_schema is not None and payload.get("schema") != expected_schema:
        raise AmendmentError(
            f"{label} schema must be {expected_schema}, "
            f"observed={payload.get('schema')!r}"
        )
    return payload, artifact


def load_jsonl_artifact(
    label: str,
    raw: object,
) -> tuple[list[dict[str, Any]], dict[str, Any], bytes]:
    artifact = normalize_artifact_ref(label, raw)
    path = Path(artifact["path"])
    records: list[dict[str, Any]] = []
    try:
        for line_number, line in enumerate(
            path.read_text(encoding="utf-8").splitlines(), start=1
        ):
            if not line:
                raise AmendmentError(f"{label} contains an empty line")
            record = strict_json_loads(line)
            if not isinstance(record, dict):
                raise AmendmentError(
                    f"{label} line {line_number} is not one JSON object"
                )
            records.append(record)
    except (OSError, json.JSONDecodeError) as exc:
        raise AmendmentError(f"{label} is not valid JSONL") from exc
    if not records:
        raise AmendmentError(f"{label} must contain at least one record")
    return records, artifact, path.read_bytes()


def load_text_artifact(label: str, raw: object) -> tuple[str, dict[str, Any]]:
    artifact = normalize_artifact_ref(label, raw)
    try:
        text = Path(artifact["path"]).read_text(encoding="utf-8")
    except OSError as exc:
        raise AmendmentError(f"{label} cannot be read") from exc
    return text, artifact


def load_spec(path: Path, expected_sha256: str) -> dict[str, Any]:
    payload, _artifact = load_json_artifact(
        "amendment spec",
        {"path": str(path), "sha256": expected_sha256},
        expected_schema=SPEC_SCHEMA,
    )
    require_exact_keys("amendment spec", payload, SPEC_KEYS)
    immutable = require_mapping(
        "spec.immutable_authority", payload["immutable_authority"]
    )
    legacy = require_mapping("spec.legacy", payload["legacy"])
    corrected = require_mapping("spec.corrected", payload["corrected"])
    capacity = require_mapping("spec.capacity", payload["capacity"])
    require_exact_keys(
        "spec.immutable_authority", immutable, IMMUTABLE_SPEC_KEYS
    )
    require_exact_keys("spec.legacy", legacy, LEGACY_SPEC_KEYS)
    require_exact_keys("spec.corrected", corrected, CORRECTED_SPEC_KEYS)
    require_exact_keys("spec.capacity", capacity, CAPACITY_SPEC_KEYS)
    return payload


def pinned_artifact_record(
    label: str,
    artifact: Mapping[str, Any],
    *,
    schema: str | None = None,
) -> dict[str, Any]:
    record = {
        "label": label,
        "path": artifact["path"],
        "resolved_path": artifact["resolved_path"],
        "sha256": artifact["sha256"],
        "size_bytes": artifact["size_bytes"],
    }
    if schema is not None:
        record["schema"] = schema
    return record


def validate_immutable_authority(
    spec: Mapping[str, Any],
) -> dict[str, Any]:
    bundle_payload, bundle_artifact = load_json_artifact(
        "immutable bundle manifest",
        spec["bundle_manifest"],
        expected_schema=resolver.BUNDLE_SCHEMA,
    )
    try:
        bundle = resolver.validate_bundle(bundle_payload)
    except resolver.ControllerError as exc:
        raise AmendmentError(
            f"immutable bundle validation failed: {exc}"
        ) from exc

    materialization_payload, materialization_artifact = load_json_artifact(
        "immutable materialization receipt",
        spec["materialization_receipt"],
        expected_schema=resolver.MATERIALIZATION_SCHEMA,
    )
    try:
        materialization = resolver.validate_materialization_binding(
            materialization_payload,
            materialization_path=Path(materialization_artifact["path"]),
            bundle_path=Path(bundle_artifact["path"]).resolve(strict=True),
            bundle_file_sha256=bundle_artifact["sha256"],
            bundle=bundle,
        )
    except resolver.ControllerError as exc:
        raise AmendmentError(
            f"immutable materialization validation failed: {exc}"
        ) from exc

    bundle_artifacts = {
        role: dict(bundle["artifact_by_role"][role])
        for role in sorted(bundle["artifact_by_role"])
    }
    return {
        "status": "PASS",
        "bundle_manifest": pinned_artifact_record(
            "immutable_bundle_manifest",
            bundle_artifact,
            schema=resolver.BUNDLE_SCHEMA,
        ),
        "materialization_receipt": pinned_artifact_record(
            "immutable_materialization_receipt",
            materialization_artifact,
            schema=resolver.MATERIALIZATION_SCHEMA,
        ),
        "bundle_identity_sha256": require_sha256(
            "immutable bundle identity",
            bundle["bundle_identity_sha256"],
        ),
        "bundle_semantic_fingerprint_sha256": require_sha256(
            "immutable bundle semantic fingerprint",
            bundle["semantic_fingerprint_sha256"],
        ),
        "public_commit": require_nonempty_text(
            "immutable public commit", bundle["public_commit"]
        ),
        "code_sha256": require_sha256(
            "immutable code SHA-256", bundle["code_sha256"]
        ),
        "replay_schema_sha256": require_sha256(
            "immutable replay schema SHA-256",
            bundle["replay_schema_sha256"],
        ),
        "training_schema_sha256": require_sha256(
            "immutable training schema SHA-256",
            bundle["training_schema_sha256"],
        ),
        "semantic_sha256": require_sha256(
            "immutable semantic SHA-256", bundle["semantic_sha256"]
        ),
        "runtime": dict(bundle["runtime"]),
        "bundle_artifacts": bundle_artifacts,
        "materialized_bundle_path": require_nonempty_text(
            "immutable materialized bundle path",
            materialization["digest_named_bundle_path"],
        ),
        "artifact_count": require_nonnegative_int(
            "immutable artifact count", materialization["artifact_count"]
        ),
        "total_artifact_bytes": require_nonnegative_int(
            "immutable total artifact bytes",
            materialization["total_artifact_bytes"],
        ),
        "readback": require_nonempty_text(
            "immutable materialization readback", materialization["readback"]
        ),
    }


def validate_preflight_immutable_binding(
    label: str,
    plan: Mapping[str, Any],
    receipt: Mapping[str, Any],
    immutable: Mapping[str, Any],
) -> None:
    input_manifests = require_mapping(
        f"{label}.input_manifests", plan.get("input_manifests")
    )
    require_exact_keys(
        f"{label}.input_manifests", input_manifests, INPUT_MANIFEST_KEYS
    )
    bundle = require_mapping(
        f"{label}.input_manifests.bundle", input_manifests["bundle"]
    )
    materialization = require_mapping(
        f"{label}.input_manifests.materialization",
        input_manifests["materialization"],
    )
    sources = require_mapping(
        f"{label}.input_manifests.sources", input_manifests["sources"]
    )
    require_exact_keys(
        f"{label}.input_manifests.bundle", bundle, BUNDLE_INPUT_KEYS
    )
    require_exact_keys(
        f"{label}.input_manifests.materialization",
        materialization,
        MATERIALIZATION_INPUT_KEYS,
    )
    require_exact_keys(
        f"{label}.input_manifests.sources", sources, SOURCE_INPUT_KEYS
    )

    authority_bundle = require_mapping(
        "immutable bundle manifest", immutable["bundle_manifest"]
    )
    authority_materialization = require_mapping(
        "immutable materialization receipt",
        immutable["materialization_receipt"],
    )
    bundle_path = require_absolute_file(
        f"{label}.input_manifests.bundle.path", bundle["path"]
    )
    materialization_path = require_absolute_file(
        f"{label}.input_manifests.materialization.path",
        materialization["path"],
    )
    if (
        bundle.get("sha256") != authority_bundle["sha256"]
        or receipt.get("bundle_manifest_sha256")
        != authority_bundle["sha256"]
    ):
        raise AmendmentError(
            f"{label} immutable bundle authority binding differs"
        )
    require_same_resolved_file(
        f"{label} immutable bundle authority",
        bundle_path,
        authority_bundle["path"],
    )
    if (
        materialization.get("sha256")
        != authority_materialization["sha256"]
        or receipt.get("materialization_receipt_sha256")
        != authority_materialization["sha256"]
        or materialization.get("bundle_identity_sha256")
        != immutable["bundle_identity_sha256"]
        or materialization.get("digest_named_bundle_path")
        != immutable["materialized_bundle_path"]
        or materialization.get("readback") != immutable["readback"]
    ):
        raise AmendmentError(
            f"{label} immutable materialization authority binding differs"
        )
    require_same_resolved_file(
        f"{label} immutable materialization authority",
        materialization_path,
        authority_materialization["path"],
    )
    if (
        bundle.get("semantic_fingerprint_sha256")
        != immutable["bundle_semantic_fingerprint_sha256"]
    ):
        raise AmendmentError(
            f"{label} immutable bundle semantic fingerprint differs"
        )


def validate_corrected_bundle_authority(
    plan: Mapping[str, Any],
    immutable: Mapping[str, Any],
) -> None:
    """Bind every corrected row bundle contract to the rehashed bundle."""

    authority_artifacts = require_mapping(
        "immutable bundle artifacts", immutable.get("bundle_artifacts")
    )
    require_exact_keys(
        "immutable bundle artifacts",
        authority_artifacts,
        resolver.REQUIRED_ARTIFACT_ROLES,
    )
    for row in require_sequence("corrected plan.rows", plan.get("rows")):
        descriptor = require_mapping("corrected plan row", row)
        row_id = require_nonempty_text(
            "corrected plan row_id", descriptor.get("row_id")
        )
        system = require_nonempty_text(
            f"{row_id}.system", descriptor.get("system")
        )
        if system not in {"pp", "auau"}:
            raise AmendmentError(f"{row_id} corrected system differs")
        contract = require_mapping(
            f"{row_id}.bundle_contract", descriptor.get("bundle_contract")
        )
        for field in (
            "public_commit",
            "code_sha256",
            "replay_schema_sha256",
            "training_schema_sha256",
            "semantic_sha256",
        ):
            if contract.get(field) != immutable[field]:
                raise AmendmentError(
                    f"{row_id} corrected immutable {field} binding differs"
                )
        expected_roles = {
            "library": f"{system}_library",
            "model": f"{system}_model",
            "config": f"{system}_config",
            "submitter": "submitter",
            "executor": f"{system}_executor",
        }
        for field, role in expected_roles.items():
            observed = require_mapping(
                f"{row_id}.bundle_contract.{field}",
                contract.get(field),
            )
            require_exact_keys(
                f"{row_id}.bundle_contract.{field}",
                observed,
                BUNDLE_ARTIFACT_OUTPUT_KEYS,
            )
            expected = require_mapping(
                f"immutable bundle artifact {role}",
                authority_artifacts.get(role),
            )
            require_exact_keys(
                f"immutable bundle artifact {role}",
                expected,
                BUNDLE_ARTIFACT_OUTPUT_KEYS,
            )
            for key in ("role", "sha256", "size_bytes"):
                if observed.get(key) != expected.get(key):
                    raise AmendmentError(
                        f"{row_id} corrected bundle artifact differs: "
                        f"{field}.{key}"
                    )
            require_same_runtime_file(
                f"{row_id} corrected bundle artifact {field}",
                observed.get("path"),
                expected.get("path"),
            )
            require_same_runtime_file(
                f"{row_id} corrected resolved bundle artifact {field}",
                observed.get("resolved_path"),
                expected.get("resolved_path"),
            )


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    try:
        with path.open("r", encoding="utf-8", newline="") as stream:
            reader = csv.DictReader(stream, delimiter="\t")
            if reader.fieldnames is None:
                raise AmendmentError(f"TSV has no header: {path}")
            if (
                any(not isinstance(field, str) or not field for field in reader.fieldnames)
                or len(reader.fieldnames) != len(set(reader.fieldnames))
            ):
                raise AmendmentError(
                    f"TSV header has empty or duplicate columns: {path}"
                )
            rows = list(reader)
    except OSError as exc:
        raise AmendmentError(f"cannot read TSV: {path}") from exc
    for index, row in enumerate(rows, start=2):
        if None in row:
            raise AmendmentError(
                f"TSV row has more fields than its header: {path}:{index}"
            )
        for field, value in row.items():
            if not isinstance(field, str) or not isinstance(value, str) or not value:
                raise AmendmentError(
                    f"TSV row has an empty or malformed cell: "
                    f"{path}:{index}:{field}"
                )
    return list(reader.fieldnames), rows


def expected_row_ids() -> list[str]:
    return [row["row_id"] for row in resolver.inventory_rows()]


def source_rows_by_id(
    label: str, manifest: Mapping[str, Any]
) -> dict[str, dict[str, Any]]:
    rows = require_sequence(f"{label}.rows", manifest.get("rows"))
    by_id: dict[str, dict[str, Any]] = {}
    for raw in rows:
        row = require_mapping(f"{label}.row", raw)
        row_id = require_nonempty_text(f"{label}.row_id", row.get("row_id", ""))
        if row_id in by_id:
            raise AmendmentError(f"{label} duplicates row_id={row_id}")
        by_id[row_id] = row
    if list(by_id) != expected_row_ids():
        raise AmendmentError(f"{label} row order/identity differs")
    return by_id


def compare_source_manifests(
    legacy: Mapping[str, Any],
    corrected: Mapping[str, Any],
) -> dict[str, Any]:
    if legacy.get("schema") != source_builder.SOURCE_SCHEMA_V1:
        raise AmendmentError("legacy source manifest is not V1")
    if corrected.get("schema") != source_builder.SOURCE_SCHEMA:
        raise AmendmentError("corrected source manifest is not V2")
    for field in (
        "status",
        "authority_state",
        "authority",
        "sim_list_root",
        "resolved_sim_list_root",
        "row_count",
        "global_five_tuple_check",
        "source_rows_semantic_sha256",
        "manifest_semantic_sha256",
    ):
        if legacy.get(field) != corrected.get(field):
            raise AmendmentError(
                f"V1/V2 source crosswalk field differs: {field}"
            )
    global_check = require_mapping(
        "corrected.global_five_tuple_check",
        corrected["global_five_tuple_check"],
    )
    require_exact_int(
        "source global total",
        global_check.get("total"),
        EXPECTED_SOURCE_TUPLES,
    )
    require_exact_int(
        "source global unique",
        global_check.get("unique"),
        EXPECTED_SOURCE_TUPLES,
    )
    require_exact_int(
        "source global duplicate_count",
        global_check.get("duplicate_count"),
        0,
    )
    legacy_by_id = source_rows_by_id("legacy source", legacy)
    corrected_by_id = source_rows_by_id("corrected source", corrected)
    crosswalk_rows: list[dict[str, Any]] = []
    for row_id in expected_row_ids():
        old = legacy_by_id[row_id]
        new = corrected_by_id[row_id]
        tuple_count = require_nonnegative_int(
            f"{row_id}.tuple_count", new.get("tuple_count")
        )
        require_exact_int(
            f"{row_id}.legacy.expected_input_count",
            old.get("expected_input_count"),
            tuple_count,
        )
        require_exact_int(
            f"{row_id}.legacy.expected_occurrence_count",
            old.get("expected_occurrence_count"),
            tuple_count,
        )
        if any(
            alias in new
            for alias in ("expected_input_count", "expected_occurrence_count")
        ):
            raise AmendmentError(
                f"{row_id} corrected V2 source leaks V1 count aliases"
            )
        comparable_old = {
            key: value
            for key, value in old.items()
            if key not in ("expected_input_count", "expected_occurrence_count")
        }
        if comparable_old != new:
            raise AmendmentError(f"{row_id} V1/V2 source payload differs")
        input_lists = require_sequence(
            f"{row_id}.input_lists", new.get("input_lists")
        )
        list_bindings = []
        for record in input_lists:
            item = require_mapping(f"{row_id}.input_list", record)
            list_bindings.append(
                {
                    key: item[key]
                    for key in (
                        "role",
                        "path",
                        "resolved_path",
                        "sha256",
                        "size_bytes",
                        "line_count",
                        "executable_count",
                    )
                }
            )
        crosswalk_rows.append(
            {
                "row_id": row_id,
                "tuple_count": tuple_count,
                "tuple_records_sha256": require_sha256(
                    f"{row_id}.tuple_records_sha256",
                    new.get("tuple_records_sha256"),
                ),
                "full_source_manifest_sha256": require_sha256(
                    f"{row_id}.full_source_manifest_sha256",
                    new.get("full_source_manifest_sha256"),
                ),
                "first_tuple_sha256": require_sha256(
                    f"{row_id}.first_tuple_sha256",
                    new.get("first_tuple_sha256"),
                ),
                "last_tuple_sha256": require_sha256(
                    f"{row_id}.last_tuple_sha256",
                    new.get("last_tuple_sha256"),
                ),
                "input_lists": list_bindings,
                "legacy_tuple_aliases_preserved_as_history_only": True,
                "corrected_execution_aliases_absent": True,
            }
        )
    if sum(row["tuple_count"] for row in crosswalk_rows) != EXPECTED_SOURCE_TUPLES:
        raise AmendmentError("source row tuple-count sum differs")
    return {
        "status": "PASS",
        "row_count": len(crosswalk_rows),
        "source_tuple_count": EXPECTED_SOURCE_TUPLES,
        "global_unique_tuple_count": EXPECTED_SOURCE_TUPLES,
        "source_rows_semantic_sha256": require_sha256(
            "source_rows_semantic_sha256",
            corrected["source_rows_semantic_sha256"],
        ),
        "manifest_semantic_sha256": require_sha256(
            "manifest_semantic_sha256",
            corrected["manifest_semantic_sha256"],
        ),
        "global_tuple_identity_records_sha256": require_sha256(
            "global tuple identity records SHA-256",
            global_check["tuple_identity_records_sha256"],
        ),
        "rows": crosswalk_rows,
    }


def require_preflight_authority(
    label: str,
    payload: Mapping[str, Any],
) -> None:
    require_bool(f"{label}.submission_performed", payload.get("submission_performed"), False)
    authority = payload.get("authority")
    if isinstance(authority, dict):
        require_exact_keys(
            f"{label}.authority",
            authority,
            PLAN_AUTHORITY_KEYS,
        )
        require_exact_int(
            f"{label}.authority.full_training_authority",
            authority.get("full_training_authority"),
            0,
        )
        if authority.get("requested_scope") != "full":
            raise AmendmentError(f"{label} requested_scope must remain full")
        if authority.get("authority_state") != "PREFLIGHT_RESOLVED_NOT_EARNED":
            raise AmendmentError(f"{label} authority_state differs")
    else:
        require_exact_int(
            f"{label}.full_training_authority",
            payload.get("full_training_authority"),
            0,
        )
        if payload.get("requested_scope") != "full":
            raise AmendmentError(f"{label} requested_scope must remain full")
        if payload.get("authority_state") != "PREFLIGHT_RESOLVED_NOT_EARNED":
            raise AmendmentError(f"{label} authority_state differs")


def validate_legacy_preflight(
    plan: Mapping[str, Any],
    receipt: Mapping[str, Any],
    source_manifest: Mapping[str, Any],
    artifacts: Mapping[str, Mapping[str, Any]],
) -> dict[str, Any]:
    require_exact_keys("legacy plan", plan, LEGACY_PLAN_KEYS)
    require_exact_keys("legacy receipt", receipt, LEGACY_RECEIPT_KEYS)
    if (
        plan.get("schema") != LEGACY_PLAN_SCHEMA
        or plan.get("status") != "PREFLIGHT_PASS"
        or plan.get("execution_state") != "PREFLIGHT_ONLY_NO_CONDOR_MUTATION"
    ):
        raise AmendmentError("legacy plan schema/status/state differs")
    if receipt.get("schema") != LEGACY_RECEIPT_SCHEMA or receipt.get("status") != "PASS":
        raise AmendmentError("legacy receipt schema/status differs")
    require_preflight_authority("legacy plan", plan)
    require_preflight_authority("legacy receipt", receipt)
    require_exact_int("legacy receipt row_count", receipt.get("row_count"), 13)

    receipt_artifacts = require_mapping(
        "legacy receipt.artifacts", receipt.get("artifacts")
    )
    require_exact_keys(
        "legacy receipt.artifacts",
        receipt_artifacts,
        LEGACY_RECEIPT_ARTIFACT_KEYS,
    )
    for artifact_name in LEGACY_RECEIPT_ARTIFACT_KEYS:
        require_exact_keys(
            f"legacy receipt.artifacts.{artifact_name}",
            require_mapping(
                f"legacy receipt.artifacts.{artifact_name}",
                receipt_artifacts[artifact_name],
            ),
            RECEIPT_ARTIFACT_KEYS,
        )
    plan_artifact = require_mapping(
        "legacy receipt.artifacts.plan", receipt_artifacts.get("plan")
    )
    if (
        plan_artifact.get("name") != Path(artifacts["plan"]["path"]).name
        or plan_artifact.get("sha256") != artifacts["plan"]["sha256"]
    ):
        raise AmendmentError("legacy receipt does not pin the supplied plan")
    if receipt.get("source_manifest_sha256") != artifacts["source_manifest"]["sha256"]:
        raise AmendmentError("legacy receipt source-manifest hash differs")
    plan_sources = require_mapping(
        "legacy plan.input_manifests.sources",
        require_mapping("legacy plan.input_manifests", plan["input_manifests"]).get(
            "sources"
        ),
    )
    if (
        plan_sources.get("sha256") != artifacts["source_manifest"]["sha256"]
        or receipt.get("source_manifest_sha256") != plan_sources.get("sha256")
    ):
        raise AmendmentError("legacy source-manifest binding differs")

    execution_partition = require_mapping(
        "legacy execution_partition", plan["execution_partition"]
    )
    if set(execution_partition) != {
        "schema",
        "group_size",
        "tuple_count",
        "expected_job_count",
        "expected_output_count",
        "basis",
        "capacity_canary_required_before_submission",
        "capacity_authority_earned",
    }:
        raise AmendmentError("legacy execution-partition inventory differs")
    if execution_partition.get("schema") != LEGACY_PARTITION_SCHEMA:
        raise AmendmentError("legacy execution-partition schema differs")
    require_exact_int(
        "legacy execution group size",
        execution_partition.get("group_size"),
        GROUP_SIZE,
    )
    require_exact_int(
        "legacy execution tuple_count",
        execution_partition.get("tuple_count"),
        EXPECTED_SOURCE_TUPLES,
    )
    require_exact_int(
        "legacy execution expected_job_count",
        execution_partition.get("expected_job_count"),
        EXPECTED_CHUNKS,
    )
    require_exact_int(
        "legacy execution expected_output_count",
        execution_partition.get("expected_output_count"),
        EXPECTED_CHUNKS,
    )
    require_bool(
        "legacy execution capacity_authority_earned",
        execution_partition.get("capacity_authority_earned"),
        False,
    )
    require_bool(
        "legacy execution capacity_canary_required_before_submission",
        execution_partition.get("capacity_canary_required_before_submission"),
        True,
    )
    execution_partition_sha = semantic_sha256(execution_partition)
    if receipt.get("execution_partition_sha256") != execution_partition_sha:
        raise AmendmentError("legacy execution-partition SHA-256 differs")

    source_by_id = source_rows_by_id("legacy source", source_manifest)
    plan_rows = require_sequence("legacy plan.rows", plan["rows"])
    if [row.get("row_id") for row in plan_rows if isinstance(row, dict)] != expected_row_ids():
        raise AmendmentError("legacy plan row order/identity differs")
    row_counts = []
    bundle_semantics: dict[str, set[str]] = {"pp": set(), "auau": set()}
    for raw in plan_rows:
        row = require_mapping("legacy plan row", raw)
        row_id = row["row_id"]
        source = source_by_id[row_id]
        contract = require_mapping(f"{row_id}.input_contract", row.get("input_contract"))
        tuple_count = require_nonnegative_int(
            f"{row_id}.tuple_count", contract.get("tuple_count")
        )
        if tuple_count != source.get("tuple_count"):
            raise AmendmentError(f"{row_id} legacy plan/source tuple count differs")
        expected_chunks = (tuple_count + GROUP_SIZE - 1) // GROUP_SIZE
        require_exact_int(
            f"{row_id}.expected_job_count",
            contract.get("expected_job_count"),
            expected_chunks,
        )
        require_exact_int(
            f"{row_id}.expected_output_count",
            contract.get("expected_output_count"),
            expected_chunks,
        )
        require_exact_int(
            f"{row_id}.expected_occurrence_count",
            contract.get("expected_occurrence_count"),
            tuple_count,
        )
        for field in (
            "tuple_records_sha256",
            "full_source_manifest_sha256",
            "first_tuple_sha256",
            "last_tuple_sha256",
        ):
            if contract.get(field) != source.get(field):
                raise AmendmentError(f"{row_id} legacy {field} differs")
        system = str(row.get("system", ""))
        if system not in bundle_semantics:
            raise AmendmentError(f"{row_id} legacy system differs")
        bundle_semantics[system].add(
            semantic_sha256(require_mapping(f"{row_id}.bundle_contract", row["bundle_contract"]))
        )
        row_counts.append(
            {
                "row_id": row_id,
                "source_tuple_count": tuple_count,
                "expected_chunk_count": expected_chunks,
                "legacy_expected_occurrence_count": tuple_count,
            }
        )
    if any(len(values) != 1 for values in bundle_semantics.values()):
        raise AmendmentError(
            "legacy plan rows do not share one bundle contract per system"
        )
    return {
        "status": "PASS",
        "plan": pinned_artifact_record(
            "legacy_plan", artifacts["plan"], schema=LEGACY_PLAN_SCHEMA
        ),
        "preflight_receipt": pinned_artifact_record(
            "legacy_preflight_receipt",
            artifacts["preflight_receipt"],
            schema=LEGACY_RECEIPT_SCHEMA,
        ),
        "source_manifest": pinned_artifact_record(
            "legacy_source_manifest",
            artifacts["source_manifest"],
            schema=source_builder.SOURCE_SCHEMA_V1,
        ),
        "submission_performed": False,
        "full_training_authority": 0,
        "execution_partition_sha256": execution_partition_sha,
        "bundle_manifest_sha256": require_sha256(
            "legacy bundle manifest SHA-256",
            receipt["bundle_manifest_sha256"],
        ),
        "materialization_receipt_sha256": require_sha256(
            "legacy materialization receipt SHA-256",
            receipt["materialization_receipt_sha256"],
        ),
        "source_manifest_sha256": artifacts["source_manifest"]["sha256"],
        "source_tuple_count": EXPECTED_SOURCE_TUPLES,
        "expected_job_count": EXPECTED_CHUNKS,
        "expected_output_count": EXPECTED_CHUNKS,
        "legacy_expected_occurrence_count": EXPECTED_SOURCE_TUPLES,
        "rows": row_counts,
    }


def canonical_jsonl_bytes(records: Sequence[Mapping[str, Any]]) -> bytes:
    return b"".join(canonical_json_bytes(record) for record in records)


def validate_corrected_preflight(
    plan: Mapping[str, Any],
    receipt: Mapping[str, Any],
    source_manifest: Mapping[str, Any],
    source_records: list[dict[str, Any]],
    rows_jsonl: Sequence[Mapping[str, Any]],
    duplicate_text: str,
    partition_records: Sequence[Mapping[str, Any]],
    partition_bytes: bytes,
    artifacts: Mapping[str, Mapping[str, Any]],
) -> dict[str, Any]:
    require_exact_keys("corrected plan", plan, CORRECTED_PLAN_KEYS)
    require_exact_keys(
        "corrected receipt", receipt, CORRECTED_RECEIPT_KEYS
    )
    if (
        plan.get("schema") != CORRECTED_PLAN_SCHEMA
        or plan.get("status") != "PREFLIGHT_PASS"
        or plan.get("execution_state") != "PREFLIGHT_ONLY_NO_CONDOR_MUTATION"
    ):
        raise AmendmentError("corrected plan schema/status/state differs")
    if (
        receipt.get("schema") != CORRECTED_RECEIPT_SCHEMA
        or receipt.get("status") != "PASS"
    ):
        raise AmendmentError("corrected receipt schema/status differs")
    require_preflight_authority("corrected plan", plan)
    require_preflight_authority("corrected receipt", receipt)
    require_exact_int("corrected receipt row_count", receipt.get("row_count"), 13)
    if source_manifest.get("schema") != CORRECTED_SOURCE_SCHEMA:
        raise AmendmentError("corrected source schema differs")

    plan_rows = require_sequence("corrected plan.rows", plan.get("rows"))
    if [row.get("row_id") for row in plan_rows if isinstance(row, dict)] != expected_row_ids():
        raise AmendmentError("corrected plan row order/identity differs")
    if list(rows_jsonl) != plan_rows:
        raise AmendmentError("corrected rows JSONL differs from plan rows")
    if partition_bytes != canonical_jsonl_bytes(partition_records):
        raise AmendmentError("corrected partition JSONL is not canonical")

    receipt_artifacts = require_mapping(
        "corrected receipt.artifacts", receipt.get("artifacts")
    )
    require_exact_keys(
        "corrected receipt.artifacts",
        receipt_artifacts,
        CORRECTED_RECEIPT_ARTIFACT_KEYS,
    )
    for artifact_name in (
        "plan",
        "rows",
        "duplicate_fingerprint",
    ):
        require_exact_keys(
            f"corrected receipt.artifacts.{artifact_name}",
            require_mapping(
                f"corrected receipt.artifacts.{artifact_name}",
                receipt_artifacts[artifact_name],
            ),
            RECEIPT_ARTIFACT_KEYS,
        )
    require_exact_keys(
        "corrected receipt.artifacts.partition",
        require_mapping(
            "corrected receipt.artifacts.partition",
            receipt_artifacts["partition"],
        ),
        PARTITION_RECEIPT_ARTIFACT_KEYS,
    )
    expected_artifact_names = {
        "plan": "plan",
        "rows": "rows",
        "duplicate_fingerprint": "duplicate_fingerprint",
        "partition": "partition",
    }
    for receipt_key, spec_key in expected_artifact_names.items():
        record = require_mapping(
            f"corrected receipt.artifacts.{receipt_key}",
            receipt_artifacts.get(receipt_key),
        )
        artifact = artifacts[spec_key]
        if (
            record.get("name") != Path(artifact["path"]).name
            or record.get("sha256") != artifact["sha256"]
        ):
            raise AmendmentError(
                f"corrected receipt does not pin supplied {spec_key}"
            )
    require_exact_int(
        "corrected partition record_count",
        require_mapping(
            "corrected receipt.artifacts.partition",
            receipt_artifacts["partition"],
        ).get("record_count"),
        EXPECTED_CHUNKS,
    )
    if receipt.get("source_manifest_sha256") != artifacts["source_manifest"]["sha256"]:
        raise AmendmentError("corrected receipt source-manifest hash differs")

    input_manifests = require_mapping(
        "corrected plan.input_manifests", plan.get("input_manifests")
    )
    require_exact_keys(
        "corrected plan.input_manifests",
        input_manifests,
        INPUT_MANIFEST_KEYS,
    )
    require_exact_keys(
        "corrected input_manifests.sources",
        require_mapping(
            "corrected input_manifests.sources", input_manifests["sources"]
        ),
        SOURCE_INPUT_KEYS,
    )
    require_exact_keys(
        "corrected input_manifests.bundle",
        require_mapping(
            "corrected input_manifests.bundle", input_manifests["bundle"]
        ),
        BUNDLE_INPUT_KEYS,
    )
    require_exact_keys(
        "corrected input_manifests.materialization",
        require_mapping(
            "corrected input_manifests.materialization",
            input_manifests["materialization"],
        ),
        MATERIALIZATION_INPUT_KEYS,
    )
    if input_manifests["sources"].get("sha256") != artifacts["source_manifest"]["sha256"]:
        raise AmendmentError("corrected plan source-manifest hash differs")
    if receipt.get("bundle_manifest_sha256") != input_manifests["bundle"].get("sha256"):
        raise AmendmentError("corrected bundle-manifest binding differs")
    if (
        receipt.get("materialization_receipt_sha256")
        != input_manifests["materialization"].get("sha256")
    ):
        raise AmendmentError("corrected materialization binding differs")

    duplicate_fingerprint = require_sha256(
        "corrected duplicate fingerprint",
        plan.get("duplicate_fingerprint_sha256"),
    )
    if duplicate_text != duplicate_fingerprint + "\n":
        raise AmendmentError("corrected duplicate-fingerprint file differs")
    if receipt.get("duplicate_fingerprint_sha256") != duplicate_fingerprint:
        raise AmendmentError("corrected duplicate fingerprint receipt differs")

    expected_partition_records = resolver.ordered_partition_chunks(source_records)
    if list(partition_records) != expected_partition_records:
        raise AmendmentError(
            "corrected partition records differ from source reconstruction"
        )
    partition_artifact_sha = artifacts["partition"]["sha256"]
    execution_partition = resolver.aggregate_partition_contract(
        source_records,
        partition_artifact_name=Path(artifacts["partition"]["path"]).name,
        partition_artifact_sha256=partition_artifact_sha,
    )
    if plan.get("execution_partition") != execution_partition:
        raise AmendmentError("corrected aggregate execution partition differs")
    execution_partition_sha = semantic_sha256(execution_partition)
    if receipt.get("execution_partition_sha256") != execution_partition_sha:
        raise AmendmentError("corrected execution-partition SHA-256 differs")

    exact_counts = {
        "source_tuple_count": EXPECTED_SOURCE_TUPLES,
        "expected_chunk_count": EXPECTED_CHUNKS,
        "expected_job_count": EXPECTED_CHUNKS,
        "expected_output_pair_count": EXPECTED_CHUNKS,
        "expected_analysis_output_count": EXPECTED_CHUNKS,
        "expected_sidecar_output_count": EXPECTED_CHUNKS,
        "expected_physical_root_artifact_count": EXPECTED_ROOT_ARTIFACTS,
        "expected_source_occurrence_count": EXPECTED_CHUNKS,
        "source_occurrences_per_output_pair": 1,
    }
    for field, expected in exact_counts.items():
        require_exact_int(
            f"corrected execution_partition.{field}",
            execution_partition.get(field),
            expected,
        )
    require_exact_int(
        "corrected execution_partition.group_size",
        execution_partition.get("group_size"),
        GROUP_SIZE,
    )
    require_bool(
        "corrected execution_partition.capacity_authority_earned",
        execution_partition.get("capacity_authority_earned"),
        False,
    )
    require_bool(
        "corrected execution_partition.capacity_canary_required_before_submission",
        execution_partition.get("capacity_canary_required_before_submission"),
        True,
    )

    source_by_id = {record["row_id"]: record for record in source_records}
    inventory_by_id = {
        record["row_id"]: record for record in resolver.inventory_rows()
    }
    row_proofs: list[dict[str, Any]] = []
    capacity_source_semantics: dict[str, dict[str, Any]] = {}
    bundle_semantics: dict[str, set[str]] = {"pp": set(), "auau": set()}
    science_semantics: set[str] = set()
    for descriptor in plan_rows:
        row_id = descriptor["row_id"]
        source = source_by_id[row_id]
        inventory = inventory_by_id[row_id]
        for field in (
            "system",
            "lane",
            "dataset",
            "sample",
            "source_role",
            "minimum_bias_gate",
            "photon_id_row_match",
        ):
            if descriptor.get(field) != inventory[field]:
                raise AmendmentError(
                    f"{row_id} corrected inventory field differs: {field}"
                )
        fingerprint_payload = {
            key: value
            for key, value in descriptor.items()
            if key != "row_fingerprint_sha256"
        }
        if descriptor.get("row_fingerprint_sha256") != semantic_sha256(
            fingerprint_payload
        ):
            raise AmendmentError(f"{row_id} row fingerprint differs")
        contract = require_mapping(
            f"{row_id}.input_contract", descriptor.get("input_contract")
        )
        partition = source["partition_contract"]
        expected_bindings = {
            "group_size": GROUP_SIZE,
            "source_tuple_count": source["tuple_count"],
            "expected_chunk_count": partition["expected_chunk_count"],
            "expected_job_count": partition["expected_job_count"],
            "expected_output_pair_count": partition["expected_output_pair_count"],
            "expected_analysis_output_count": partition[
                "expected_analysis_output_count"
            ],
            "expected_sidecar_output_count": partition[
                "expected_sidecar_output_count"
            ],
            "expected_source_occurrence_count": partition[
                "expected_source_occurrence_count"
            ],
            "source_occurrences_per_output_pair": 1,
            "row_partition_sha256": partition["partition_sha256"],
            "chunk_records_sha256": partition["chunk_records_sha256"],
            "partition_artifact_name": Path(artifacts["partition"]["path"]).name,
            "partition_artifact_sha256": partition_artifact_sha,
            "full_source_manifest_sha256": source[
                "full_source_manifest_sha256"
            ],
            "tuple_records_sha256": source["tuple_records_sha256"],
            "first_tuple_sha256": source["first_tuple_sha256"],
            "last_tuple_sha256": source["last_tuple_sha256"],
        }
        for field, expected in expected_bindings.items():
            if contract.get(field) != expected:
                raise AmendmentError(f"{row_id} corrected {field} differs")
        execution_contract = require_mapping(
            f"{row_id}.execution_contract",
            descriptor.get("execution_contract"),
        )
        if (
            execution_contract.get("schema") != resolver.EXECUTION_SCHEMA
            or execution_contract.get("state")
            != "PREFLIGHT_ONLY_NO_CONDOR_MUTATION"
            or execution_contract.get("partition_artifact_sha256")
            != partition_artifact_sha
            or execution_contract.get("row_partition_sha256")
            != partition["partition_sha256"]
        ):
            raise AmendmentError(
                f"{row_id} corrected execution-contract binding differs"
            )
        system = str(descriptor.get("system", ""))
        if system not in bundle_semantics:
            raise AmendmentError(f"{row_id} corrected system differs")
        bundle_contract = require_mapping(
            f"{row_id}.bundle_contract",
            descriptor.get("bundle_contract"),
        )
        bundle_semantics[system].add(semantic_sha256(bundle_contract))
        science_semantics.add(
            semantic_sha256(
                require_mapping(
                    f"{row_id}.science_contract",
                    descriptor.get("science_contract"),
                )
            )
        )
        sample = require_nonempty_text(
            f"{row_id}.sample", descriptor.get("sample")
        )
        run_match = re.match(r"^run([0-9]+)_", sample)
        if run_match is None:
            raise AmendmentError(
                f"{row_id} sample does not encode one frozen run"
            )
        config = require_mapping(
            f"{row_id}.bundle_contract.config",
            bundle_contract.get("config"),
        )
        library = require_mapping(
            f"{row_id}.bundle_contract.library",
            bundle_contract.get("library"),
        )
        model = require_mapping(
            f"{row_id}.bundle_contract.model",
            bundle_contract.get("model"),
        )
        capacity_source_semantics[row_id] = {
            "row_id": row_id,
            "system": system,
            "lane": descriptor["lane"],
            "dataset": descriptor["dataset"],
            "sample": sample,
            "source_role": descriptor["source_role"],
            "minimum_bias_gate": descriptor["minimum_bias_gate"],
            "period": descriptor["source_period"],
            "run": int(run_match.group(1)),
            "si_di_role": descriptor["source_si_di_role"],
            "ownership_state": "source_role_frozen",
            "code_sha256": require_sha256(
                f"{row_id}.bundle_contract.code_sha256",
                bundle_contract.get("code_sha256"),
            ),
            "replay_schema_sha256": require_sha256(
                f"{row_id}.bundle_contract.replay_schema_sha256",
                bundle_contract.get("replay_schema_sha256"),
            ),
            "training_schema_sha256": require_sha256(
                f"{row_id}.bundle_contract.training_schema_sha256",
                bundle_contract.get("training_schema_sha256"),
            ),
            "semantic_sha256": require_sha256(
                f"{row_id}.bundle_contract.semantic_sha256",
                bundle_contract.get("semantic_sha256"),
            ),
            "library_sha256": require_sha256(
                f"{row_id}.bundle_contract.library.sha256",
                library.get("sha256"),
            ),
            "model_sha256": require_sha256(
                f"{row_id}.bundle_contract.model.sha256",
                model.get("sha256"),
            ),
            "base_config_sha256": require_sha256(
                f"{row_id}.bundle_contract.config.sha256",
                config.get("sha256"),
            ),
        }
        row_proofs.append(
            {
                "row_id": row_id,
                "source_tuple_count": source["tuple_count"],
                "expected_chunk_count": partition["expected_chunk_count"],
                "expected_job_count": partition["expected_job_count"],
                "expected_output_pair_count": partition[
                    "expected_output_pair_count"
                ],
                "expected_source_occurrence_count": partition[
                    "expected_source_occurrence_count"
                ],
                "tail_tuple_count": partition["tail_tuple_count"],
                "row_partition_sha256": partition["partition_sha256"],
                "chunk_records_sha256": partition["chunk_records_sha256"],
                "first_execution_chunk_sha256": next(
                    record["execution_chunk_sha256"]
                    for record in partition_records
                    if record["row_id"] == row_id
                    and record["chunk_index"] == 0
                ),
                "row_fingerprint_sha256": descriptor[
                    "row_fingerprint_sha256"
                ],
            }
        )
    if (
        any(len(values) != 1 for values in bundle_semantics.values())
        or len(science_semantics) != 1
    ):
        raise AmendmentError(
            "corrected rows do not share one bundle/science contract"
        )
    if sum(row["expected_chunk_count"] for row in row_proofs) != EXPECTED_CHUNKS:
        raise AmendmentError("corrected rowwise chunk sum differs")

    execution_fingerprint = semantic_sha256(
        {
            "schema": resolver.EXECUTION_SCHEMA,
            "tag": plan["campaign"]["tag"],
            "output_root": plan["campaign"]["output_root"],
            "evidence_root": plan["campaign"]["evidence_root"],
            "submit_root": plan["campaign"]["submit_root"],
            "materialization_receipt_sha256": receipt[
                "materialization_receipt_sha256"
            ],
            "bundle_manifest_sha256": receipt["bundle_manifest_sha256"],
            "source_manifest_sha256": receipt["source_manifest_sha256"],
            "duplicate_fingerprint_sha256": duplicate_fingerprint,
            "partition_artifact_sha256": partition_artifact_sha,
            "execution_partition_sha256": execution_partition_sha,
            "row_fingerprints": [
                descriptor["row_fingerprint_sha256"]
                for descriptor in plan_rows
            ],
        }
    )
    if (
        plan.get("execution_fingerprint_sha256") != execution_fingerprint
        or receipt.get("execution_fingerprint_sha256") != execution_fingerprint
    ):
        raise AmendmentError("corrected execution fingerprint differs")

    return {
        "status": "PASS",
        "plan": pinned_artifact_record(
            "corrected_plan", artifacts["plan"], schema=CORRECTED_PLAN_SCHEMA
        ),
        "preflight_receipt": pinned_artifact_record(
            "corrected_preflight_receipt",
            artifacts["preflight_receipt"],
            schema=CORRECTED_RECEIPT_SCHEMA,
        ),
        "source_manifest": pinned_artifact_record(
            "corrected_source_manifest",
            artifacts["source_manifest"],
            schema=CORRECTED_SOURCE_SCHEMA,
        ),
        "rows": pinned_artifact_record("corrected_rows", artifacts["rows"]),
        "duplicate_fingerprint": pinned_artifact_record(
            "corrected_duplicate_fingerprint",
            artifacts["duplicate_fingerprint"],
        ),
        "partition": {
            **pinned_artifact_record(
                "corrected_partition", artifacts["partition"]
            ),
            "record_count": len(partition_records),
        },
        "submission_performed": False,
        "full_training_authority": 0,
        "bundle_manifest_sha256": require_sha256(
            "corrected bundle manifest SHA-256",
            receipt["bundle_manifest_sha256"],
        ),
        "materialization_receipt_sha256": require_sha256(
            "corrected materialization receipt SHA-256",
            receipt["materialization_receipt_sha256"],
        ),
        "source_manifest_sha256": artifacts["source_manifest"]["sha256"],
        "duplicate_fingerprint_sha256": duplicate_fingerprint,
        "execution_partition_sha256": execution_partition_sha,
        "execution_fingerprint_sha256": execution_fingerprint,
        **exact_counts,
        "rows_proof": row_proofs,
        "capacity_source_semantics": [
            capacity_source_semantics[row_id]
            for row_id in expected_row_ids()
        ],
        "bundle_contract_semantic_sha256_by_system": {
            system: next(iter(values))
            for system, values in bundle_semantics.items()
        },
        "science_contract_semantic_sha256": next(iter(science_semantics)),
    }


SUBMISSION_RECEIPT_FIELDS = (
    "row_id",
    "cluster_proc",
    "submit_file",
    "args_file",
    "staged_chunk_list",
    "staged_chunk_sha256",
    "fanout_contract_file",
    "fanout_contract_sha256",
    "condor_log",
    "condor_stdout",
    "condor_stderr",
    "analysis_output_root",
    "multiview_sidecar",
    "sidecar_owner_count",
    "submitted_args_sha256",
    "snapshot_dir",
    "snapshot_builder_header",
    "snapshot_builder_header_sha256",
    "snapshot_calo_reco_library",
    "snapshot_calo_reco_library_sha256",
    "snapshot_analysis_library",
    "snapshot_analysis_library_sha256",
    "materialized_config",
    "materialized_config_sha256",
    "snapshot_loader_receipt",
    "snapshot_loader_receipt_sha256",
    "snapshot_manifest",
    "snapshot_manifest_sha256",
)

SUBMISSION_MANIFEST_FIELDS = (
    "row_id",
    "system",
    "lane",
    "dataset",
    "sample",
    "source_role",
    "minimum_bias_gate",
    "input_files",
    "input_jobs",
    "source_manifest_sha256",
    "first_input_tuple_sha256",
    "resolved_config",
    "resolved_config_sha256",
    "library",
    "library_sha256",
    "model",
    "model_sha256",
    "code_sha256",
    "replay_schema_sha256",
    "training_schema_sha256",
    "semantic_sha256",
    "nominal_et_min_gev",
    "nominal_et_max_gev_exclusive",
    "loose_capture_et_min_gev",
    "legacy_training_tree_max_entries",
    "event_limit_per_job",
    "analysis_output_namespace",
    "multiview_sidecar",
    "submit_namespace",
    "scheduler_log_dir",
    "scheduler_stdout_dir",
    "scheduler_stderr_dir",
    "executable_input_tuple_count",
    "full_training_authority",
    "photon_cluster_builder_header",
    "photon_cluster_builder_header_sha256",
    "calo_reco_library",
    "calo_reco_library_sha256",
    "release_core_lib_dir",
    "release_core_lib64_dir",
    "release_calo_io",
    "release_calo_io_sha256",
    "release_clusteriso",
    "release_clusteriso_sha256",
    "release_jetbase",
    "release_jetbase_sha256",
)

BASE_TO_RESOLVED_KEYS = {
    "pp": frozenset(
        {
            "coneR",
            "pp_photonid_extract_only",
            "pp_photonid_training_tree",
            "pp_photonid_training_tree_max_entries",
            "pp_photonid_source_role",
            "pp_photonid_ppg12_filter",
            "pp_photonid_require_preselection",
        }
    ),
    "auau": frozenset({"coneR", "auau_bdt_training_tree"}),
}
EXPECTED_BASE_TO_RESOLVED_VALUES = {
    "pp": {
        "base": {"coneR": "[0.40, 0.30]"},
        "resolved": {
            "coneR": "[0.40]",
            "pp_photonid_extract_only": "true",
            "pp_photonid_training_tree": "true",
            "pp_photonid_training_tree_max_entries": "0",
            "pp_photonid_source_role": "background",
            "pp_photonid_ppg12_filter": "true",
            "pp_photonid_require_preselection": "false",
        },
    },
    "auau": {
        "base": {
            "coneR": "[0.40, 0.30]",
            "auau_bdt_training_tree": "false",
        },
        "resolved": {
            "coneR": "[0.40]",
            "auau_bdt_training_tree": "true",
        },
    },
}
MATERIALIZED_CONFIG_KEYS = frozenset(
    {
        "jet_pt_min",
        "back_to_back_dphi_min_pi_fraction",
        "vz_cut_cm",
        "clusterUEpipeline",
        "fixedGeV",
        "coneR",
        "preselection",
        "tight",
        "nonTight",
    }
)
EXPECTED_RESOLVED_MATERIALIZATION_VALUES = {
    "pp": {
        "resolved": {
            "jet_pt_min": "[5.0, 7.0, 10.0, 12.0]",
            "back_to_back_dphi_min_pi_fraction": "[0.5, 0.875]",
            "vz_cut_cm": "[10, 30, 60]",
            "clusterUEpipeline": "[noSub]",
            "fixedGeV": "[0.0]",
            "coneR": "[0.40]",
        },
        "materialized": {
            "jet_pt_min": "5.0",
            "back_to_back_dphi_min_pi_fraction": "0.875",
            "vz_cut_cm": "60",
            "clusterUEpipeline": "noSub",
            "fixedGeV": "2.0",
            "coneR": "0.40",
            "preselection": "newPPG12",
            "tight": "newPPG12",
            "nonTight": "newPPG12",
        },
    },
    "auau": {
        "resolved": {
            "jet_pt_min": "[5.0, 7.0, 10.0, 12.0]",
            "back_to_back_dphi_min_pi_fraction": "[0.5, 0.875]",
            "vz_cut_cm": "[10, 30, 60]",
            "clusterUEpipeline": "[baseVariant]",
            "fixedGeV": "0.0",
            "coneR": "[0.40]",
        },
        "materialized": {
            "jet_pt_min": "5.0",
            "back_to_back_dphi_min_pi_fraction": "0.875",
            "vz_cut_cm": "10",
            "clusterUEpipeline": "baseVariant",
            "fixedGeV": "0.0",
            "coneR": "0.40",
            "preselection": "newPPG12",
            "tight": "auauCentInputBase3x3BDT",
            "nonTight": "auauBDTSideband",
        },
    },
}
SNAPSHOT_LOADER_SCHEMA = "RJ_SNAPSHOT_LOADER_RECEIPT_V1"
SNAPSHOT_MANIFEST_SCHEMA = "RJ_FROZEN_SNAPSHOT_MANIFEST_V1"
SNAPSHOT_LOADER_KEYS = frozenset(
    {
        "schema",
        "status",
        "mode",
        "snapshot_lib",
        "release_roots",
        "providers",
        "targets_inspected",
        "pinned_calo_reco_release_companions",
        "ldd_report",
    }
)
SNAPSHOT_MANIFEST_KEYS = frozenset({"schema", "status", "root", "entries"})
SNAPSHOT_PROVIDER_FAMILIES = (
    "libcalo_reco.so",
    "libcalo_io.so",
    "libclusteriso.so",
    "libjetbase.so",
)
RUNTIME_AUTHORITY_SCHEMA = "THE134_SINGLE_PROVIDER_RUNTIME_AUTHORITY_V2"
RUNTIME_AUTHORITY_KEYS = frozenset(
    {
        "calo_reco",
        "pp_sim_weight_contract",
        "release",
        "schema",
        "status",
    }
)
RUNTIME_CALO_RECO_KEYS = frozenset(
    {
        "build_receipt",
        "build_receipt_sha256",
        "library",
        "library_sha256",
        "soname",
        "source_manifest",
        "source_manifest_sha256",
    }
)
RUNTIME_RELEASE_KEYS = frozenset(
    {"lib", "lib64", "name", "offline_main", "providers"}
)
RUNTIME_PROVIDER_KEYS = frozenset({"path", "sha256"})
RUNTIME_WEIGHT_KEYS = frozenset(
    {
        "interaction",
        "mix_weight",
        "period",
        "period_lumi_weight",
        "vertex_reweight",
    }
)
EXPECTED_PP_SIM_WEIGHT_CONTRACT = {
    "interaction": "SI",
    "mix_weight": "period_auto",
    "period": "0mrad",
    "period_lumi_weight": True,
    "vertex_reweight": "period_auto",
}
RUNTIME_PROVIDER_ROLE_BY_FAMILY = {
    "libcalo_io.so": "release_calo_io",
    "libclusteriso.so": "release_clusteriso",
    "libjetbase.so": "release_jetbase",
}
MUTABLE_PRIVATE_PROVIDER_RE = re.compile(
    r"^/sphenix/(?:u|user)/[^/]+/"
    r"(?:(?:thesisAnalysis|thesisAnalysis_auau)/)?install/lib(?:64)?/"
)

ROOT_ROW_KEYS = frozenset(
    {
        "row_id",
        "status",
        "analysis_output_root",
        "analysis_bytes",
        "analysis_sha256",
        "sidecar",
        "sidecar_bytes",
        "sidecar_sha256",
        "sidecar_size_diagnostic",
        "replay_tree_count",
        "replay_sources",
        "source_contract",
        "source_execution_contract",
        "source_identity_canonical_sha256",
        "source_occurrence_id_hex",
        "replay_events",
        "replay_candidates",
        "sidecar_entries",
        "sidecar_candidates",
        "joined_sidecar_rows",
        "population_state",
        "snapshot_loader_receipt_sha256",
        "snapshot_manifest_sha256",
        "full_training_authority",
    }
)


def verify_nested_file_reference(
    label: str,
    raw: object,
    *,
    allowed_extra: Iterable[str] = (),
) -> dict[str, Any]:
    payload = require_mapping(label, raw)
    expected_keys = {"path", "sha256", *allowed_extra}
    require_exact_keys(label, payload, expected_keys)
    artifact = normalize_artifact_ref(
        label,
        {"path": payload["path"], "sha256": payload["sha256"]},
    )
    result = pinned_artifact_record(label, artifact)
    for field in allowed_extra:
        result[field] = payload[field]
    return result


def validate_runtime_authority_artifact(
    artifact: Mapping[str, Any],
    immutable: Mapping[str, Any],
) -> dict[str, Any]:
    """Reparse the capacity runtime authority and bind every provider."""

    bundle_artifacts = require_mapping(
        "immutable bundle artifacts", immutable.get("bundle_artifacts")
    )
    bundle_authority = require_mapping(
        "immutable runtime authority manifest",
        bundle_artifacts.get("runtime_authority_manifest"),
    )
    if artifact.get("sha256") != bundle_authority.get("sha256"):
        raise AmendmentError(
            "capacity runtime authority differs from the immutable bundle"
        )
    authority_path = require_absolute_file(
        "capacity runtime authority path", artifact.get("path")
    )
    try:
        payload = strict_json_loads(
            authority_path.read_text(encoding="utf-8")
        )
    except (OSError, json.JSONDecodeError) as exc:
        raise AmendmentError(
            "capacity runtime authority is not valid JSON"
        ) from exc
    payload = require_mapping("capacity runtime authority payload", payload)
    require_exact_keys(
        "capacity runtime authority payload",
        payload,
        RUNTIME_AUTHORITY_KEYS,
    )
    if (
        payload.get("schema") != RUNTIME_AUTHORITY_SCHEMA
        or payload.get("status") != "PASS"
    ):
        raise AmendmentError(
            "capacity runtime authority schema/status differs"
        )

    runtime = require_mapping("immutable runtime", immutable.get("runtime"))
    calo_reco = require_mapping(
        "capacity runtime authority calo_reco", payload.get("calo_reco")
    )
    require_exact_keys(
        "capacity runtime authority calo_reco",
        calo_reco,
        RUNTIME_CALO_RECO_KEYS,
    )
    if calo_reco.get("soname") != runtime["calo_reco_soname"]:
        raise AmendmentError("capacity runtime CaloReco SONAME differs")
    calo_roles = {
        "build_receipt": "calo_reco_build_receipt",
        "library": "calo_reco_library",
        "source_manifest": "calo_reco_source_manifest",
    }
    normalized_calo: dict[str, dict[str, Any]] = {}
    for field, role in calo_roles.items():
        observed = normalize_runtime_artifact_ref(
            f"capacity runtime authority calo_reco.{field}",
            calo_reco.get(field),
            calo_reco.get(f"{field}_sha256"),
        )
        expected = require_mapping(
            f"immutable bundle artifact {role}",
            bundle_artifacts.get(role),
        )
        if (
            observed["sha256"] != expected.get("sha256")
            or observed["size_bytes"] != expected.get("size_bytes")
        ):
            raise AmendmentError(
                f"capacity runtime authority CaloReco copy differs: {field}"
            )
        normalized_calo[field] = observed

    release = require_mapping(
        "capacity runtime authority release", payload.get("release")
    )
    require_exact_keys(
        "capacity runtime authority release", release, RUNTIME_RELEASE_KEYS
    )
    if release.get("name") != runtime["release"]:
        raise AmendmentError("capacity runtime release identity differs")
    require_same_runtime_directory(
        "capacity runtime authority offline_main",
        release.get("offline_main"),
        runtime["offline_main"],
    )
    release_dirs = {
        "lib": require_runtime_directory(
            "capacity runtime authority release.lib", release.get("lib")
        ),
        "lib64": require_runtime_directory(
            "capacity runtime authority release.lib64", release.get("lib64")
        ),
    }
    require_same_runtime_directory(
        "capacity runtime authority release lib",
        release_dirs["lib"],
        runtime["release_core_lib_dir"],
    )
    require_same_runtime_directory(
        "capacity runtime authority release lib64",
        release_dirs["lib64"],
        runtime["release_core_lib64_dir"],
    )

    providers = require_mapping(
        "capacity runtime authority release.providers",
        release.get("providers"),
    )
    require_exact_keys(
        "capacity runtime authority release.providers",
        providers,
        RUNTIME_PROVIDER_ROLE_BY_FAMILY,
    )
    normalized_providers: dict[str, dict[str, Any]] = {}
    for family, role in RUNTIME_PROVIDER_ROLE_BY_FAMILY.items():
        provider = require_mapping(
            f"capacity runtime authority provider {family}",
            providers.get(family),
        )
        require_exact_keys(
            f"capacity runtime authority provider {family}",
            provider,
            RUNTIME_PROVIDER_KEYS,
        )
        observed = normalize_runtime_artifact_ref(
            f"capacity runtime authority provider {family}",
            provider.get("path"),
            provider.get("sha256"),
        )
        expected = require_mapping(
            f"immutable bundle artifact {role}",
            bundle_artifacts.get(role),
        )
        if (
            observed["sha256"] != expected.get("sha256")
            or observed["size_bytes"] != expected.get("size_bytes")
        ):
            raise AmendmentError(
                f"capacity runtime authority provider copy differs: {family}"
            )
        if not runtime_path_beneath(
            Path(observed["resolved_path"]),
            release_dirs.values(),
        ):
            raise AmendmentError(
                f"capacity runtime provider escaped the pinned release: "
                f"{family}"
            )
        normalized_providers[family] = observed

    weight_contract = require_mapping(
        "capacity runtime authority pp_sim_weight_contract",
        payload.get("pp_sim_weight_contract"),
    )
    require_exact_keys(
        "capacity runtime authority pp_sim_weight_contract",
        weight_contract,
        RUNTIME_WEIGHT_KEYS,
    )
    if weight_contract != EXPECTED_PP_SIM_WEIGHT_CONTRACT:
        raise AmendmentError(
            "capacity runtime p+p simulation weight contract differs"
        )
    return {
        "payload": payload,
        "calo_reco": normalized_calo,
        "release_dirs": {
            name: str(path) for name, path in release_dirs.items()
        },
        "providers": normalized_providers,
    }


def config_values_and_residual(
    label: str,
    text: str,
    mutable_keys: frozenset[str],
) -> tuple[dict[str, str], list[str]]:
    """Split exact top-level YAML scalar overrides from untouched lines."""

    values: dict[str, str] = {}
    residual: list[str] = []
    for raw_line in text.splitlines():
        if not raw_line.strip():
            continue
        match = re.match(r"^([A-Za-z_][A-Za-z0-9_]*):[ \t]*(.*)$", raw_line)
        if match is None or match.group(1) not in mutable_keys:
            residual.append(raw_line)
            continue
        key = match.group(1)
        if key in values:
            raise AmendmentError(f"{label} duplicates mutable key {key}")
        value = re.sub(r"[ \t]+#.*$", "", match.group(2)).strip()
        if not value:
            raise AmendmentError(f"{label} has an empty mutable key {key}")
        values[key] = value
    return values, residual


def validate_runtime_config_chain(
    row_id: str,
    system: str,
    *,
    base_config: Mapping[str, Any],
    resolved_config: Mapping[str, Any],
    materialized_config: Mapping[str, Any],
    fanout_contract: Mapping[str, Any],
    analysis_output_namespace: str,
) -> str:
    """Prove the exact immutable-base -> resolved -> materialized chain."""

    if system not in {"pp", "auau"}:
        raise AmendmentError(f"{row_id} runtime-config system differs")
    base_text = Path(str(base_config["path"])).read_text(encoding="utf-8")
    resolved_text = Path(str(resolved_config["path"])).read_text(
        encoding="utf-8"
    )
    materialized_text = Path(str(materialized_config["path"])).read_text(
        encoding="utf-8"
    )

    transition_keys = BASE_TO_RESOLVED_KEYS[system]
    base_values, base_residual = config_values_and_residual(
        f"{row_id}.base_config", base_text, transition_keys
    )
    resolved_transition_values, resolved_transition_residual = (
        config_values_and_residual(
            f"{row_id}.resolved_config", resolved_text, transition_keys
        )
    )
    expected_transition = EXPECTED_BASE_TO_RESOLVED_VALUES[system]
    if (
        base_values != expected_transition["base"]
        or resolved_transition_values != expected_transition["resolved"]
        or base_residual != resolved_transition_residual
    ):
        raise AmendmentError(
            f"{row_id} base-to-resolved config scalarization differs"
        )

    resolved_values, resolved_residual = config_values_and_residual(
        f"{row_id}.resolved_materialization",
        resolved_text,
        MATERIALIZED_CONFIG_KEYS,
    )
    materialized_values, materialized_residual = config_values_and_residual(
        f"{row_id}.materialized_config",
        materialized_text,
        MATERIALIZED_CONFIG_KEYS,
    )
    expected_materialization = EXPECTED_RESOLVED_MATERIALIZATION_VALUES[
        system
    ]
    if (
        resolved_values != expected_materialization["resolved"]
        or materialized_values != expected_materialization["materialized"]
        or resolved_residual != materialized_residual
    ):
        raise AmendmentError(
            f"{row_id} resolved-to-materialized config scalarization differs"
        )

    fanout_lines = [
        line
        for line in Path(str(fanout_contract["path"]))
        .read_text(encoding="utf-8")
        .splitlines()
        if line.strip()
    ]
    if len(fanout_lines) != 1:
        raise AmendmentError(
            f"{row_id} fanout contract must contain exactly one row"
        )
    fanout_fields = fanout_lines[0].split("|")
    if len(fanout_fields) != 5 or any(not field for field in fanout_fields):
        raise AmendmentError(
            f"{row_id} fanout contract must contain five nonempty fields"
        )
    output_path = Path(fanout_fields[0])
    if (
        not output_path.is_absolute()
        or ".." in output_path.parts
        or output_path.name != fanout_fields[1]
        or str(output_path.parent)
        != require_nonempty_text(
            f"{row_id}.analysis_output_namespace",
            analysis_output_namespace,
        ).rstrip("/")
    ):
        raise AmendmentError(
            f"{row_id} fanout output namespace/identity differs"
        )
    expected_tags = expected_materialization["materialized"]
    if fanout_fields[2:] != [
        expected_tags["preselection"],
        expected_tags["tight"],
        expected_tags["nonTight"],
    ]:
        raise AmendmentError(f"{row_id} fanout selection identity differs")

    return semantic_sha256(
        {
            "schema": "THE134_RUNTIME_CONFIG_CHAIN_V1",
            "row_id": row_id,
            "system": system,
            "base_config_sha256": base_config["sha256"],
            "resolved_config_sha256": resolved_config["sha256"],
            "materialized_config_sha256": materialized_config["sha256"],
            "fanout_contract_sha256": fanout_contract["sha256"],
            "base_to_resolved_values": {
                "base": base_values,
                "resolved": resolved_transition_values,
            },
            "resolved_to_materialized_values": {
                "resolved": resolved_values,
                "materialized": materialized_values,
            },
            "fanout_fields": fanout_fields,
        }
    )


def require_snapshot_member(
    label: str, path: Path, snapshot_root: Path
) -> Path:
    try:
        resolved = path.resolve(strict=True)
    except (FileNotFoundError, OSError, RuntimeError) as exc:
        raise AmendmentError(
            f"{label} does not resolve to an existing snapshot member: {path}"
        ) from exc
    if snapshot_root != resolved and snapshot_root not in resolved.parents:
        raise AmendmentError(
            f"{label} is outside the frozen snapshot: {path}"
        )
    return resolved


def runtime_path_beneath(path: Path, roots: Iterable[Path]) -> bool:
    resolved = path.resolve(strict=True)
    for root in roots:
        resolved_root = root.resolve(strict=True)
        try:
            resolved.relative_to(resolved_root)
            return True
        except ValueError:
            continue
    return False


def snapshot_provider_family(name: str) -> str | None:
    for family in SNAPSHOT_PROVIDER_FAMILIES:
        if name == family or name.startswith(f"{family}."):
            return family
    return None


def validate_snapshot_ldd_report(
    row_id: str,
    system: str,
    *,
    ldd_path: Path,
    snapshot_root: Path,
    expected_providers: Mapping[str, Path],
    loader_providers: Mapping[str, Any],
    release_roots: Sequence[Path],
) -> None:
    """Replay the sealed loader report without executing mutable tooling."""

    snapshot_lib = (snapshot_root / "lib").resolve(strict=True)
    analysis_library = snapshot_lib / (
        "libRecoilJets.so" if system == "pp" else "libRecoilJetsAuAu.so"
    )
    expected_targets = [
        analysis_library,
        *[
            require_runtime_file(
                f"{row_id}.expected ldd target {family}",
                expected_providers[family],
            )
            for family in SNAPSHOT_PROVIDER_FAMILIES
        ],
    ]
    try:
        report_lines = ldd_path.read_text(
            encoding="utf-8", errors="replace"
        ).splitlines()
    except OSError as exc:
        raise AmendmentError(
            f"{row_id} snapshot ldd report cannot be read"
        ) from exc

    staged_names = {entry.name for entry in snapshot_lib.iterdir()}
    staged_families = {
        name.split(".so", 1)[0] + ".so"
        for name in staged_names
        if ".so" in name
    }
    observed_dependencies: dict[str, list[Path]] = {
        family: [] for family in SNAPSHOT_PROVIDER_FAMILIES
    }
    section_observations = [0 for _target in expected_targets]
    target_index = 0
    current_target_index: int | None = None
    for raw in report_lines:
        line = raw.strip()
        if line.startswith("@@TARGET"):
            match = re.fullmatch(r"@@TARGET (/.+)", line)
            if match is None or target_index >= len(expected_targets):
                raise AmendmentError(
                    f"{row_id} malformed or extra ldd target marker"
                )
            target = require_runtime_file(
                f"{row_id}.ldd target {target_index}",
                match.group(1),
            )
            if not os.path.samefile(target, expected_targets[target_index]):
                raise AmendmentError(
                    f"{row_id} ldd target order/identity differs at "
                    f"index {target_index}"
                )
            current_target_index = target_index
            target_index += 1
            continue
        if "=> not found" in line:
            raise AmendmentError(
                f"{row_id} ldd report contains an unresolved dependency"
            )
        match = re.match(r"(\S+)\s+=>\s+(\S+)\s+", line)
        if match is None:
            if (
                current_target_index is not None
                and re.fullmatch(r"\S.*\s+\(0x[0-9A-Fa-f]+\)", line)
            ):
                section_observations[current_target_index] += 1
            continue
        if current_target_index is None:
            raise AmendmentError(
                f"{row_id} ldd dependency precedes its target marker"
            )
        section_observations[current_target_index] += 1
        name, resolved_text = match.groups()
        if MUTABLE_PRIVATE_PROVIDER_RE.match(resolved_text):
            raise AmendmentError(
                f"{row_id} ldd report uses a mutable private provider"
            )
        resolved = require_runtime_file(
            f"{row_id}.ldd dependency {name}", resolved_text
        )
        staged_family = (
            name.split(".so", 1)[0] + ".so" if ".so" in name else name
        )
        if (
            name in staged_names or staged_family in staged_families
        ) and not runtime_path_beneath(resolved, (snapshot_lib,)):
            raise AmendmentError(
                f"{row_id} staged ldd dependency escaped the snapshot: "
                f"{name}"
            )
        family = snapshot_provider_family(name)
        if family is None:
            continue
        expected = expected_providers[family]
        if not os.path.samefile(resolved, expected):
            raise AmendmentError(
                f"{row_id} ldd provider identity differs: {family}"
            )
        if family != "libcalo_reco.so" and not runtime_path_beneath(
            resolved, release_roots
        ):
            raise AmendmentError(
                f"{row_id} ldd provider escaped the pinned release: "
                f"{family}"
            )
        observed_dependencies[family].append(resolved)

    if target_index != len(expected_targets):
        raise AmendmentError(
            f"{row_id} ldd report target count differs: "
            f"expected={len(expected_targets)} observed={target_index}"
        )
    if any(count == 0 for count in section_observations):
        raise AmendmentError(
            f"{row_id} ldd report contains an empty target section"
        )
    for family, expected in expected_providers.items():
        if any(
            not os.path.samefile(path, expected)
            for path in observed_dependencies[family]
        ):
            raise AmendmentError(
                f"{row_id} ldd provider closure differs: {family}"
            )
        provider = require_mapping(
            f"{row_id}.snapshot provider {family}",
            loader_providers.get(family),
        )
        observed_resolutions = require_sequence(
            f"{row_id}.snapshot provider {family}.observed_resolutions",
            provider.get("observed_resolutions"),
        )
        if len(observed_resolutions) != 1:
            raise AmendmentError(
                f"{row_id} loader provider is not singleton: {family}"
            )
        require_same_runtime_file(
            f"{row_id} loader observed provider {family}",
            observed_resolutions[0],
            expected,
        )
        require_same_runtime_file(
            f"{row_id} loader real provider {family}",
            provider.get("realpath"),
            expected,
        )
        if family != "libcalo_reco.so":
            require_same_runtime_file(
                f"{row_id} loader declared provider {family}",
                provider.get("declared_path"),
                expected,
            )


def validate_frozen_snapshot_receipts(
    row_id: str,
    system: str,
    receipt_row: Mapping[str, str],
    manifest_row: Mapping[str, str],
    immutable: Mapping[str, Any],
    runtime_authority: Mapping[str, Any],
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Replay the sealed snapshot and loader-receipt authority checks."""

    snapshot_dir = Path(
        require_nonempty_text(
            f"{row_id}.receipt.snapshot_dir", receipt_row.get("snapshot_dir")
        )
    )
    if (
        not snapshot_dir.is_absolute()
        or ".." in snapshot_dir.parts
        or not snapshot_dir.is_dir()
        or snapshot_dir.is_symlink()
    ):
        raise AmendmentError(
            f"{row_id} snapshot_dir must be an absolute real directory"
        )
    snapshot_root = snapshot_dir.resolve(strict=True)
    if snapshot_root.stat().st_mode & 0o222:
        raise AmendmentError(f"{row_id} frozen snapshot root is writable")

    loader_payload, loader_artifact = load_json_artifact(
        f"{row_id}.snapshot_loader_receipt",
        {
            "path": receipt_row.get("snapshot_loader_receipt"),
            "sha256": receipt_row.get("snapshot_loader_receipt_sha256"),
        },
        expected_schema=SNAPSHOT_LOADER_SCHEMA,
    )
    manifest_payload, manifest_artifact = load_json_artifact(
        f"{row_id}.snapshot_manifest",
        {
            "path": receipt_row.get("snapshot_manifest"),
            "sha256": receipt_row.get("snapshot_manifest_sha256"),
        },
        expected_schema=SNAPSHOT_MANIFEST_SCHEMA,
    )
    declared_loader_path = Path(loader_artifact["path"])
    declared_manifest_path = Path(manifest_artifact["path"])
    if (
        declared_loader_path
        != snapshot_dir / "snapshot_loader_receipt.json"
        or declared_manifest_path != snapshot_dir / "snapshot_manifest.json"
    ):
        raise AmendmentError(
            f"{row_id} snapshot receipt/manifest paths are not lexically "
            "canonical"
        )
    loader_path = require_snapshot_member(
        f"{row_id}.snapshot_loader_receipt",
        declared_loader_path,
        snapshot_root,
    )
    snapshot_manifest_path = require_snapshot_member(
        f"{row_id}.snapshot_manifest",
        declared_manifest_path,
        snapshot_root,
    )
    lexical_loader_path = snapshot_root / "snapshot_loader_receipt.json"
    lexical_manifest_path = snapshot_root / "snapshot_manifest.json"
    for label, lexical_path in (
        ("snapshot loader receipt", lexical_loader_path),
        ("snapshot manifest", lexical_manifest_path),
    ):
        if (
            lexical_path.is_symlink()
            or not lexical_path.is_file()
            or lexical_path.lstat().st_nlink != 1
        ):
            raise AmendmentError(
                f"{row_id} canonical {label} must be an unlinked regular "
                "non-symlink file"
            )
    expected_loader_path = lexical_loader_path.resolve(strict=True)
    expected_manifest_path = lexical_manifest_path.resolve(strict=True)
    if (
        loader_path != expected_loader_path
        or snapshot_manifest_path != expected_manifest_path
    ):
        raise AmendmentError(
            f"{row_id} snapshot receipt/manifest names are not canonical"
        )
    if (
        loader_path.stat().st_mode & 0o222
        or snapshot_manifest_path.stat().st_mode & 0o222
    ):
        raise AmendmentError(f"{row_id} frozen snapshot receipt is writable")

    require_exact_keys(
        f"{row_id}.snapshot_loader_receipt",
        loader_payload,
        SNAPSHOT_LOADER_KEYS,
    )
    if (
        loader_payload.get("status") != "PASS"
        or loader_payload.get("mode") != system
        or loader_payload.get("pinned_calo_reco_release_companions") is not True
    ):
        raise AmendmentError(f"{row_id} snapshot loader contract differs")
    require_exact_int(
        f"{row_id}.snapshot_loader_receipt.targets_inspected",
        loader_payload.get("targets_inspected"),
        5,
    )
    runtime = require_mapping(
        "immutable runtime", immutable.get("runtime")
    )
    release_roots = require_sequence(
        f"{row_id}.snapshot release_roots",
        loader_payload.get("release_roots"),
    )
    if len(release_roots) != 2:
        raise AmendmentError(f"{row_id} snapshot release roots differ")
    require_same_runtime_directory(
        f"{row_id} snapshot release root lib64",
        release_roots[0],
        runtime["release_core_lib64_dir"],
    )
    require_same_runtime_directory(
        f"{row_id} snapshot release root lib",
        release_roots[1],
        runtime["release_core_lib_dir"],
    )
    snapshot_lib = Path(
        require_nonempty_text(
            f"{row_id}.snapshot_loader_receipt.snapshot_lib",
            loader_payload.get("snapshot_lib"),
        )
    )
    if (
        snapshot_lib != snapshot_dir / "lib"
        or not snapshot_lib.is_dir()
        or not os.path.samefile(snapshot_lib, snapshot_root / "lib")
    ):
        raise AmendmentError(f"{row_id} snapshot library root differs")

    ldd_report = require_mapping(
        f"{row_id}.snapshot_loader_receipt.ldd_report",
        loader_payload.get("ldd_report"),
    )
    require_exact_keys(
        f"{row_id}.snapshot_loader_receipt.ldd_report",
        ldd_report,
        {"path", "sha256"},
    )
    ldd_relative = Path(
        require_nonempty_text(
            f"{row_id}.snapshot_loader_receipt.ldd_report.path",
            ldd_report.get("path"),
        )
    )
    if ldd_relative.parts != ("snapshot_loader_receipt.ldd.txt",):
        raise AmendmentError(f"{row_id} snapshot ldd report path is unsafe")
    lexical_ldd_path = snapshot_root / "snapshot_loader_receipt.ldd.txt"
    if (
        lexical_ldd_path.is_symlink()
        or not lexical_ldd_path.is_file()
        or lexical_ldd_path.lstat().st_nlink != 1
    ):
        raise AmendmentError(
            f"{row_id} canonical snapshot ldd report must be an unlinked "
            "regular non-symlink file"
        )
    ldd_path = require_snapshot_member(
        f"{row_id}.snapshot_loader_receipt.ldd_report",
        lexical_ldd_path,
        snapshot_root,
    )
    if ldd_path != lexical_ldd_path.resolve(strict=True):
        raise AmendmentError(f"{row_id} snapshot ldd report name differs")
    if sha256_file(ldd_path) != require_sha256(
        f"{row_id}.snapshot_loader_receipt.ldd_report.sha256",
        ldd_report.get("sha256"),
    ):
        raise AmendmentError(f"{row_id} snapshot ldd report hash drift")

    bundle_artifacts = require_mapping(
        "immutable bundle artifacts", immutable.get("bundle_artifacts")
    )
    role_by_manifest_field = {
        "photon_cluster_builder_header": "photon_cluster_builder_header",
        "calo_reco_library": "calo_reco_library",
        "release_calo_io": "release_calo_io",
        "release_clusteriso": "release_clusteriso",
        "release_jetbase": "release_jetbase",
    }
    for manifest_field, role in role_by_manifest_field.items():
        expected = require_mapping(
            f"immutable bundle artifact {role}",
            bundle_artifacts.get(role),
        )
        if manifest_row.get(f"{manifest_field}_sha256") != expected["sha256"]:
            raise AmendmentError(
                f"{row_id} manifest immutable {manifest_field} hash differs"
            )

    expected_snapshot_authorities = (
        (
            "snapshot_builder_header",
            snapshot_dir / "PhotonClusterBuilder.h",
            snapshot_root / "PhotonClusterBuilder.h",
            manifest_row["photon_cluster_builder_header_sha256"],
        ),
        (
            "snapshot_calo_reco_library",
            snapshot_dir / "lib/libcalo_reco.so",
            snapshot_root / "lib/libcalo_reco.so",
            manifest_row["calo_reco_library_sha256"],
        ),
        (
            "snapshot_analysis_library",
            snapshot_dir
            / "lib"
            / (
                "libRecoilJets.so"
                if system == "pp"
                else "libRecoilJetsAuAu.so"
            ),
            snapshot_root
            / "lib"
            / (
                "libRecoilJets.so"
                if system == "pp"
                else "libRecoilJetsAuAu.so"
            ),
            manifest_row["library_sha256"],
        ),
    )
    for receipt_field, expected_lexical_path, expected_path, expected_sha in (
        expected_snapshot_authorities
    ):
        observed = normalize_runtime_artifact_ref(
            f"{row_id}.{receipt_field}",
            receipt_row.get(receipt_field),
            receipt_row.get(f"{receipt_field}_sha256"),
        )
        if Path(observed["path"]) != expected_lexical_path:
            raise AmendmentError(
                f"{row_id}.{receipt_field} path is not lexically canonical"
            )
        require_same_runtime_file(
            f"{row_id}.{receipt_field} ownership",
            observed["path"],
            expected_path,
        )
        if observed["sha256"] != expected_sha:
            raise AmendmentError(
                f"{row_id} snapshotted authority differs: {receipt_field}"
            )
    forbidden_snapshot_patterns = (
        "libphoton_cluster_builder_override.so*",
        "libcalo_io.so*",
        "libclusteriso.so*",
        "libjetbase.so*",
    )
    for pattern in forbidden_snapshot_patterns:
        if any((snapshot_root / "lib").glob(pattern)):
            raise AmendmentError(
                f"{row_id} forbidden snapshot-local provider is present: "
                f"{pattern}"
            )

    calo_api = snapshot_root / "lib/libcalo_reco.so"
    calo_soname = snapshot_root / "lib" / runtime["calo_reco_soname"]
    if (
        not calo_api.is_file()
        or not calo_soname.is_symlink()
        or os.readlink(calo_soname) != "libcalo_reco.so"
        or not os.path.samefile(calo_soname, calo_api)
    ):
        raise AmendmentError(f"{row_id} snapshot CaloReco SONAME differs")

    runtime_providers = require_mapping(
        "capacity runtime authority providers",
        runtime_authority.get("providers"),
    )
    expected_provider_paths = {
        "libcalo_reco.so": calo_api,
        "libcalo_io.so": require_runtime_file(
            f"{row_id}.runtime libcalo_io.so",
            require_mapping(
                f"{row_id}.runtime libcalo_io.so",
                runtime_providers.get("libcalo_io.so"),
            ).get("path"),
        ),
        "libclusteriso.so": require_runtime_file(
            f"{row_id}.runtime libclusteriso.so",
            require_mapping(
                f"{row_id}.runtime libclusteriso.so",
                runtime_providers.get("libclusteriso.so"),
            ).get("path"),
        ),
        "libjetbase.so": require_runtime_file(
            f"{row_id}.runtime libjetbase.so",
            require_mapping(
                f"{row_id}.runtime libjetbase.so",
                runtime_providers.get("libjetbase.so"),
            ).get("path"),
        ),
    }
    expected_provider_hashes = {
        "libcalo_reco.so": manifest_row["calo_reco_library_sha256"],
        "libcalo_io.so": manifest_row["release_calo_io_sha256"],
        "libclusteriso.so": manifest_row["release_clusteriso_sha256"],
        "libjetbase.so": manifest_row["release_jetbase_sha256"],
    }
    providers = require_mapping(
        f"{row_id}.snapshot_loader_receipt.providers",
        loader_payload.get("providers"),
    )
    require_exact_keys(
        f"{row_id}.snapshot_loader_receipt.providers",
        providers,
        SNAPSHOT_PROVIDER_FAMILIES,
    )
    for family in SNAPSHOT_PROVIDER_FAMILIES:
        expected_path = expected_provider_paths[family]
        expected_sha = expected_provider_hashes[family]
        provider = require_mapping(
            f"{row_id}.snapshot provider {family}", providers.get(family)
        )
        provider_keys = {"observed_resolutions", "realpath", "sha256"}
        if family != "libcalo_reco.so":
            provider_keys.add("declared_path")
        require_exact_keys(
            f"{row_id}.snapshot provider {family}",
            provider,
            provider_keys,
        )
        observed_resolutions = require_sequence(
            f"{row_id}.snapshot provider {family}.observed_resolutions",
            provider.get("observed_resolutions"),
        )
        if (
            len(observed_resolutions) != 1
            or provider.get("sha256") != expected_sha
            or sha256_file(expected_path) != expected_sha
        ):
            raise AmendmentError(
                f"{row_id} snapshot loader provider differs: {family}"
            )
        require_same_runtime_file(
            f"{row_id} snapshot loader real provider {family}",
            provider.get("realpath"),
            expected_path,
        )
        require_same_runtime_file(
            f"{row_id} snapshot loader observed provider {family}",
            observed_resolutions[0],
            expected_path,
        )
        if family != "libcalo_reco.so":
            require_same_runtime_file(
                f"{row_id} snapshot declared provider {family}",
                provider.get("declared_path"),
                expected_path,
            )
    validate_snapshot_ldd_report(
        row_id,
        system,
        ldd_path=ldd_path,
        snapshot_root=snapshot_root,
        expected_providers=expected_provider_paths,
        loader_providers=providers,
        release_roots=[
            require_runtime_directory(
                f"{row_id}.snapshot release root {index}", root
            )
            for index, root in enumerate(release_roots)
        ],
    )

    require_exact_keys(
        f"{row_id}.snapshot_manifest",
        manifest_payload,
        SNAPSHOT_MANIFEST_KEYS,
    )
    if manifest_payload.get("status") != "PASS":
        raise AmendmentError(f"{row_id} snapshot manifest is not PASS")
    manifest_root = Path(
        require_nonempty_text(
            f"{row_id}.snapshot_manifest.root",
            manifest_payload.get("root"),
        )
    )
    if (
        manifest_root != snapshot_dir
        or not manifest_root.is_dir()
        or not os.path.samefile(manifest_root, snapshot_root)
    ):
        raise AmendmentError(f"{row_id} snapshot manifest root differs")

    declared_entries = require_sequence(
        f"{row_id}.snapshot_manifest.entries",
        manifest_payload.get("entries"),
    )
    declared_paths: set[str] = set()
    for index, raw_entry in enumerate(declared_entries):
        entry = require_mapping(
            f"{row_id}.snapshot_manifest.entries[{index}]", raw_entry
        )
        entry_type = entry.get("type")
        expected_entry_keys = {
            "directory": {"path", "type"},
            "file": {"path", "sha256", "size", "type"},
            "symlink": {
                "path",
                "sha256",
                "symlink_target",
                "type",
            },
        }.get(str(entry_type))
        if expected_entry_keys is None:
            raise AmendmentError(
                f"{row_id} snapshot manifest entry type differs"
            )
        require_exact_keys(
            f"{row_id}.snapshot_manifest.entries[{index}]",
            entry,
            expected_entry_keys,
        )
        relative = require_nonempty_text(
            f"{row_id}.snapshot_manifest.entries[{index}].path",
            entry.get("path"),
        )
        relative_path = Path(relative)
        if (
            relative_path.is_absolute()
            or ".." in relative_path.parts
            or relative in declared_paths
            or relative == "snapshot_manifest.json"
        ):
            raise AmendmentError(
                f"{row_id} snapshot manifest entry path differs: {relative}"
            )
        declared_paths.add(relative)
        if entry_type in {"file", "symlink"}:
            require_sha256(
                f"{row_id}.snapshot_manifest.entries[{index}].sha256",
                entry.get("sha256"),
            )
        if entry_type == "file":
            size = entry.get("size")
            if (
                isinstance(size, bool)
                or not isinstance(size, int)
                or size < 0
            ):
                raise AmendmentError(
                    f"{row_id} snapshot manifest file size differs"
                )

    observed_entries: list[dict[str, Any]] = []
    for path in sorted(
        snapshot_root.rglob("*"),
        key=lambda item: item.relative_to(snapshot_root).as_posix(),
    ):
        if path == snapshot_root / "snapshot_manifest.json":
            continue
        relative = path.relative_to(snapshot_root).as_posix()
        if path.is_symlink():
            target = os.readlink(path)
            if os.path.isabs(target):
                raise AmendmentError(
                    f"{row_id} snapshot symlink must be relative: {relative}"
                )
            resolved = require_snapshot_member(
                f"{row_id}.snapshot_symlink.{relative}",
                path,
                snapshot_root,
            )
            if (
                not resolved.is_file()
                or
                path.parent.stat().st_mode & 0o222
                or resolved.stat().st_mode & 0o222
            ):
                raise AmendmentError(
                    f"{row_id} writable snapshot symlink survived: {relative}"
                )
            observed_entries.append(
                {
                    "path": relative,
                    "sha256": sha256_file(resolved),
                    "symlink_target": target,
                    "type": "symlink",
                }
            )
        elif path.is_file():
            if (
                path.lstat().st_mode & 0o222
                or path.lstat().st_nlink != 1
            ):
                raise AmendmentError(
                    f"{row_id} writable or externally linked snapshot file "
                    f"survived: {relative}"
                )
            observed_entries.append(
                {
                    "path": relative,
                    "sha256": sha256_file(path),
                    "size": path.stat().st_size,
                    "type": "file",
                }
            )
        elif path.is_dir():
            if path.lstat().st_mode & 0o222:
                raise AmendmentError(
                    f"{row_id} writable snapshot directory survived: "
                    f"{relative}"
                )
            observed_entries.append({"path": relative, "type": "directory"})
        else:
            raise AmendmentError(
                f"{row_id} unsupported snapshot artifact: {relative}"
            )
    if declared_entries != observed_entries:
        raise AmendmentError(
            f"{row_id} complete snapshot inventory/hash closure differs"
        )
    return (
        pinned_artifact_record(
            f"{row_id}_snapshot_loader_receipt",
            loader_artifact,
            schema=SNAPSHOT_LOADER_SCHEMA,
        ),
        pinned_artifact_record(
            f"{row_id}_snapshot_manifest",
            manifest_artifact,
            schema=SNAPSHOT_MANIFEST_SCHEMA,
        ),
    )


def materialize_first_chunk_bytes(
    source: Mapping[str, Any],
) -> tuple[bytes, dict[str, Any]]:
    chunks = source.get("_chunks")
    if not isinstance(chunks, list) or not chunks:
        raise AmendmentError(
            f"{source.get('row_id')} corrected chunks are unavailable"
        )
    chunk = require_mapping(f"{source['row_id']}.chunk0", chunks[0])
    require_exact_int(
        f"{source['row_id']}.chunk0.index", chunk.get("chunk_index"), 0
    )
    require_exact_int(
        f"{source['row_id']}.chunk0.segment", chunk.get("segment"), 1
    )
    require_exact_int(
        f"{source['row_id']}.chunk0.tuple_count",
        chunk.get("tuple_count"),
        GROUP_SIZE,
    )
    input_lists = require_sequence(
        f"{source['row_id']}.input_lists", source.get("input_lists")
    )
    by_role = {
        record["role"]: require_mapping(
            f"{source['row_id']}.{record.get('role')}", record
        )
        for record in input_lists
        if isinstance(record, dict)
    }
    if set(by_role) != set(resolver.LIST_ROLES):
        raise AmendmentError(
            f"{source['row_id']} corrected source-list roles differ"
        )
    lines_by_role: dict[str, list[str]] = {}
    for role in resolver.LIST_ROLES:
        path = require_absolute_file(
            f"{source['row_id']}.{role}.path", by_role[role]["path"]
        )
        lines_by_role[role] = path.read_text(encoding="utf-8").splitlines()
    output_lines: list[str] = []
    reconstructed_inputs: list[dict[str, str]] = []
    for physical_line in chunk["physical_lines"]:
        if (
            isinstance(physical_line, bool)
            or not isinstance(physical_line, int)
            or physical_line <= 0
        ):
            raise AmendmentError(
                f"{source['row_id']} chunk physical line is invalid"
            )
        values: dict[str, str] = {}
        for role in resolver.LIST_ROLES:
            try:
                value = resolver.executable_line(
                    lines_by_role[role][physical_line - 1]
                )
            except IndexError as exc:
                raise AmendmentError(
                    f"{source['row_id']} chunk physical line is out of range"
                ) from exc
            if value is None:
                raise AmendmentError(
                    f"{source['row_id']} chunk references a non-executable line"
                )
            values[role] = value
        reconstructed_inputs.append(values)
        output_lines.append("\t".join(values[role] for role in resolver.LIST_ROLES))
    expected_tuple_hashes = [
        resolver.canonical_sha256(inputs) for inputs in reconstructed_inputs
    ]
    if chunk.get("tuple_input_sha256s") != expected_tuple_hashes:
        raise AmendmentError(
            f"{source['row_id']} reconstructed first-chunk tuple hashes differ"
        )
    return ("\n".join(output_lines) + "\n").encode("utf-8"), chunk


def canonical_source_identity_sha256(
    source_contract: Mapping[str, Any],
) -> str:
    fields = (
        "lane",
        "dataset",
        "sample",
        "period",
        "run",
        "segment",
        "input_uri_hash",
        "input_file_sha256",
        "source_manifest_sha256",
    )
    canonical = "|".join(str(source_contract[field]) for field in fields)
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def validate_capacity_resource_row(
    row_id: str,
    raw: Mapping[str, Any],
) -> dict[str, Any]:
    row = dict(raw)
    require_exact_keys(f"{row_id} capacity resource row", row, RESOURCE_ROW_KEYS)
    if row.get("row_id") != row_id:
        raise AmendmentError(f"{row_id} capacity resource row identity differs")
    if not CLUSTER_PROC_RE.fullmatch(str(row.get("cluster_proc", ""))):
        raise AmendmentError(f"{row_id} capacity cluster.proc is malformed")
    for field, expected in (
        ("exit_code", 0),
        ("job_status", 4),
        ("num_holds", 0),
        ("num_job_starts", 1),
        ("request_memory_mb", REQUEST_MEMORY_MB),
    ):
        require_exact_int(f"{row_id}.{field}", row.get(field), expected)
    memory_usage_mb = row.get("memory_usage_mb")
    if (
        isinstance(memory_usage_mb, bool)
        or not isinstance(memory_usage_mb, int)
        or memory_usage_mb <= 0
        or memory_usage_mb > REQUEST_MEMORY_MB
    ):
        raise AmendmentError(
            f"{row_id}.memory_usage_mb must be within 1..{REQUEST_MEMORY_MB}"
        )
    resident_set_size_kb = row.get("resident_set_size_kb")
    maximum_rss_kb = REQUEST_MEMORY_MB * 1024
    if (
        isinstance(resident_set_size_kb, bool)
        or not isinstance(resident_set_size_kb, int)
        or resident_set_size_kb <= 0
        or resident_set_size_kb > maximum_rss_kb
    ):
        raise AmendmentError(
            f"{row_id}.resident_set_size_kb must be within "
            f"1..{maximum_rss_kb}"
        )
    if (
        row.get("memory_usage_resolution")
        != "EVALUATED_FROM_RESIDENT_SET_SIZE_KB"
    ):
        raise AmendmentError(f"{row_id}.memory_usage_resolution differs")
    wall = row.get("remote_wall_clock_seconds")
    if (
        isinstance(wall, bool)
        or not isinstance(wall, (int, float))
        or wall <= 0
    ):
        raise AmendmentError(
            f"{row_id}.remote_wall_clock_seconds must be positive"
        )
    return row


def validate_capacity_source_binding(
    row_id: str,
    source_contract: Mapping[str, Any],
    source_execution: Mapping[str, Any],
    corrected_semantics: Mapping[str, Any],
) -> None:
    require_exact_keys(
        f"{row_id}.source_contract",
        source_contract,
        SOURCE_CONTRACT_KEYS,
    )
    require_exact_keys(
        f"{row_id}.source_execution_contract",
        source_execution,
        SOURCE_EXECUTION_CONTRACT_KEYS,
    )
    for field in (
        "lane",
        "dataset",
        "sample",
        "period",
        "run",
        "si_di_role",
        "ownership_state",
    ):
        if source_contract.get(field) != corrected_semantics[field]:
            raise AmendmentError(
                f"{row_id} source contract differs from corrected plan: "
                f"{field}"
            )
    if source_execution.get("run") != corrected_semantics["run"]:
        raise AmendmentError(
            f"{row_id} source execution run differs from corrected plan"
        )
    for hash_field in (
        "input_uri_hash",
        "input_file_sha256",
        "source_manifest_sha256",
    ):
        require_sha256(
            f"{row_id}.{hash_field}", source_contract.get(hash_field)
        )


def validate_submission_manifest_record(
    row_id: str,
    manifest_row: Mapping[str, str],
    *,
    source_contract: Mapping[str, Any],
    corrected_semantics: Mapping[str, Any],
    corrected_source: Mapping[str, Any],
    immutable: Mapping[str, Any],
    runtime_authority: Mapping[str, Any],
) -> dict[str, Any]:
    require_exact_keys(
        f"{row_id}.submission_manifest_row",
        manifest_row,
        SUBMISSION_MANIFEST_FIELDS,
    )
    for field in (
        "system",
        "lane",
        "dataset",
        "sample",
        "source_role",
        "minimum_bias_gate",
    ):
        if manifest_row.get(field) != str(corrected_semantics[field]):
            raise AmendmentError(
                f"{row_id} submission manifest differs from corrected "
                f"source semantics: {field}"
            )
    if (
        manifest_row.get("source_manifest_sha256")
        != source_contract.get("source_manifest_sha256")
    ):
        raise AmendmentError(
            f"{row_id} submission manifest source authority differs"
        )
    if (
        require_sha256(
            f"{row_id}.submission_manifest.code_sha256",
            manifest_row.get("code_sha256"),
        )
        != corrected_semantics["code_sha256"]
    ):
        raise AmendmentError(
            f"{row_id} submission manifest code authority differs"
        )
    runtime_config = normalize_artifact_ref(
        f"{row_id}.submission_manifest.resolved_config",
        {
            "path": manifest_row.get("resolved_config"),
            "sha256": manifest_row.get("resolved_config_sha256"),
        },
    )
    if (
        manifest_row.get("first_input_tuple_sha256")
        != corrected_source.get("first_tuple_sha256")
    ):
        raise AmendmentError(
            f"{row_id} submission manifest first-tuple authority differs"
        )
    require_exact_int(
        f"{row_id}.submission_manifest.input_files",
        require_decimal_int_text(
            f"{row_id}.submission_manifest.input_files",
            manifest_row.get("input_files"),
        ),
        GROUP_SIZE,
    )
    require_exact_int(
        f"{row_id}.submission_manifest.input_jobs",
        require_decimal_int_text(
            f"{row_id}.submission_manifest.input_jobs",
            manifest_row.get("input_jobs"),
        ),
        1,
    )
    require_exact_int(
        f"{row_id}.submission_manifest.executable_input_tuple_count",
        require_decimal_int_text(
            f"{row_id}.submission_manifest.executable_input_tuple_count",
            manifest_row.get("executable_input_tuple_count"),
        ),
        require_nonnegative_int(
            f"{row_id}.corrected_source.tuple_count",
            corrected_source.get("tuple_count"),
        ),
    )
    require_exact_int(
        f"{row_id}.submission_manifest.full_training_authority",
        require_decimal_int_text(
            f"{row_id}.submission_manifest.full_training_authority",
            manifest_row.get("full_training_authority"),
        ),
        0,
    )
    exact_science_values = {
        "nominal_et_min_gev": "15.0",
        "nominal_et_max_gev_exclusive": "35.0",
        "loose_capture_et_min_gev": "5.0",
        "legacy_training_tree_max_entries": "0",
        "event_limit_per_job": "0",
    }
    for field, expected in exact_science_values.items():
        if manifest_row.get(field) != expected:
            raise AmendmentError(
                f"{row_id} submission manifest scientific field differs: "
                f"{field}"
            )
    for field in (
        "source_manifest_sha256",
        "first_input_tuple_sha256",
        "resolved_config_sha256",
        "library_sha256",
        "model_sha256",
        "code_sha256",
        "replay_schema_sha256",
        "training_schema_sha256",
        "semantic_sha256",
        "photon_cluster_builder_header_sha256",
        "calo_reco_library_sha256",
        "release_calo_io_sha256",
        "release_clusteriso_sha256",
        "release_jetbase_sha256",
    ):
        require_sha256(
            f"{row_id}.submission_manifest.{field}",
            manifest_row.get(field),
        )

    for field in (
        "replay_schema_sha256",
        "training_schema_sha256",
        "semantic_sha256",
    ):
        if manifest_row.get(field) != immutable[field]:
            raise AmendmentError(
                f"{row_id} submission manifest immutable {field} differs"
            )
    bundle_artifacts = require_mapping(
        "immutable bundle artifacts", immutable.get("bundle_artifacts")
    )
    system = str(corrected_semantics["system"])
    expected_roles = {
        "library": f"{system}_library",
        "model": f"{system}_model",
        "photon_cluster_builder_header": "photon_cluster_builder_header",
        "calo_reco_library": "calo_reco_library",
        "release_calo_io": "release_calo_io",
        "release_clusteriso": "release_clusteriso",
        "release_jetbase": "release_jetbase",
    }
    runtime_artifacts: dict[str, dict[str, Any]] = {}
    for field, role in expected_roles.items():
        artifact = normalize_runtime_artifact_ref(
            f"{row_id}.submission_manifest.{field}",
            manifest_row.get(field),
            manifest_row.get(f"{field}_sha256"),
        )
        authority = require_mapping(
            f"immutable bundle artifact {role}",
            bundle_artifacts.get(role),
        )
        if (
            artifact["sha256"] != authority.get("sha256")
            or artifact["size_bytes"] != authority.get("size_bytes")
        ):
            raise AmendmentError(
                f"{row_id} submission manifest {field} authority differs"
            )
        if field in {"library", "photon_cluster_builder_header"}:
            require_same_runtime_file(
                f"{row_id} submission manifest {field}",
                artifact["path"],
                authority.get("path"),
            )
        runtime_artifacts[field] = artifact
    runtime_copy_bindings = {
        "calo_reco_library": require_mapping(
            "capacity runtime authority CaloReco library",
            require_mapping(
                "capacity runtime authority CaloReco",
                runtime_authority.get("calo_reco"),
            ).get("library"),
        ),
        "release_calo_io": require_mapping(
            "capacity runtime authority libcalo_io.so",
            require_mapping(
                "capacity runtime authority providers",
                runtime_authority.get("providers"),
            ).get("libcalo_io.so"),
        ),
        "release_clusteriso": require_mapping(
            "capacity runtime authority libclusteriso.so",
            require_mapping(
                "capacity runtime authority providers",
                runtime_authority.get("providers"),
            ).get("libclusteriso.so"),
        ),
        "release_jetbase": require_mapping(
            "capacity runtime authority libjetbase.so",
            require_mapping(
                "capacity runtime authority providers",
                runtime_authority.get("providers"),
            ).get("libjetbase.so"),
        ),
    }
    for field, authority in runtime_copy_bindings.items():
        require_same_runtime_file(
            f"{row_id} submission manifest runtime authority {field}",
            runtime_artifacts[field]["path"],
            authority["path"],
        )
    runtime = require_mapping("immutable runtime", immutable.get("runtime"))
    runtime_release_dirs = require_mapping(
        "capacity runtime authority release directories",
        runtime_authority.get("release_dirs"),
    )
    release_dir_bindings = {
        "release_core_lib_dir": "lib",
        "release_core_lib64_dir": "lib64",
    }
    for field, runtime_field in release_dir_bindings.items():
        require_same_runtime_directory(
            f"{row_id} submission manifest immutable {field}",
            manifest_row.get(field),
            runtime[field],
        )
        require_same_runtime_directory(
            f"{row_id} submission manifest runtime authority {field}",
            manifest_row.get(field),
            runtime_release_dirs[runtime_field],
        )
    return {
        "runtime_config": runtime_config,
        "runtime_artifacts": runtime_artifacts,
        "row_semantic_sha256": semantic_sha256(dict(manifest_row)),
    }


def validate_source_provenance_artifact(
    row_id: str,
    source_provenance: Mapping[str, Any],
    *,
    sidecar_path: Path,
    staged_sha256: str,
    source_contract: Mapping[str, Any],
    corrected_semantics: Mapping[str, Any],
    submission_manifest_row: Mapping[str, str],
) -> dict[str, Any]:
    require_exact_keys(
        f"{row_id}.source_provenance",
        source_provenance,
        AUDIT_SOURCE_PROVENANCE_KEYS,
    )
    provenance_artifact = normalize_artifact_ref(
        f"{row_id}.source_provenance",
        {
            "path": source_provenance.get("path"),
            "sha256": source_provenance.get("sha256"),
        },
    )
    try:
        provenance_payload = strict_json_loads(
            Path(provenance_artifact["path"]).read_text(encoding="utf-8")
        )
    except (OSError, json.JSONDecodeError) as exc:
        raise AmendmentError(
            f"{row_id} source provenance is not valid JSON"
        ) from exc
    provenance_payload = require_mapping(
        f"{row_id}.source_provenance.payload", provenance_payload
    )
    require_exact_keys(
        f"{row_id}.source_provenance.payload",
        provenance_payload,
        SOURCE_PROVENANCE_KEYS,
    )
    if provenance_payload.get("schema") != "THE134_SOURCE_PROVENANCE_V1":
        raise AmendmentError(f"{row_id} source provenance schema differs")
    provenance_inputs = require_sequence(
        f"{row_id}.source_provenance.inputs",
        provenance_payload.get("inputs"),
    )
    if len(provenance_inputs) != 1:
        raise AmendmentError(
            f"{row_id} source provenance must contain exactly one input"
        )
    provenance_record = require_mapping(
        f"{row_id}.source_provenance.input", provenance_inputs[0]
    )
    require_exact_keys(
        f"{row_id}.source_provenance.input",
        provenance_record,
        SOURCE_PROVENANCE_RECORD_KEYS,
    )
    embedded_record = require_mapping(
        f"{row_id}.source_provenance.record",
        source_provenance.get("record"),
    )
    require_exact_keys(
        f"{row_id}.source_provenance.record",
        embedded_record,
        SOURCE_PROVENANCE_RECORD_KEYS,
    )
    expected_provenance_record = {
        "path": str(sidecar_path),
        "row_id": row_id,
        "system": corrected_semantics["system"],
        "source_sample": corrected_semantics["sample"],
        "input_uri_sha256": staged_sha256,
        "input_file_sha256": staged_sha256,
        "source_manifest_sha256": source_contract["source_manifest_sha256"],
        "config_sha256": submission_manifest_row["resolved_config_sha256"],
        "code_sha256": submission_manifest_row["code_sha256"],
    }
    for observed in (provenance_record, embedded_record):
        require_same_resolved_file(
            f"{row_id} source provenance sidecar",
            observed.get("path"),
            sidecar_path,
        )
    without_path = lambda record: {
        key: value for key, value in record.items() if key != "path"
    }
    if (
        without_path(provenance_record) != without_path(embedded_record)
        or without_path(provenance_record)
        != without_path(expected_provenance_record)
    ):
        raise AmendmentError(
            f"{row_id} source provenance content differs from corrected "
            "plan and capacity source authority"
        )
    return provenance_artifact


def validate_capacity_non_training_audit(
    row_id: str,
    payload: Mapping[str, Any],
    artifact: Mapping[str, Any],
) -> None:
    """Replay the exact capacity-only, non-training audit authority."""

    expected_system = EXPECTED_SYSTEM[row_id]
    if (
        payload.get("scope") != "capacity"
        or payload.get("authority_state") != "CAPACITY_NON_TRAINING_ONLY"
        or payload.get("row_id") != row_id
        or payload.get("system") != expected_system
        or payload.get("class_balance_state")
        != "NOT_APPLICABLE_CAPACITY_NON_TRAINING"
        or payload.get("failures") != []
    ):
        raise AmendmentError(
            f"{row_id} capacity audit non-training authority differs"
        )
    require_exact_int(
        f"{row_id}.source_occurrence_count",
        payload.get("source_occurrence_count"),
        1,
    )
    require_exact_int(
        f"{row_id}.full_training_authority",
        payload.get("full_training_authority"),
        0,
    )
    require_bool(
        f"{row_id}.matrix_materialized",
        payload.get("matrix_materialized"),
        False,
    )
    require_bool(
        f"{row_id}.training_matrix_authority_earned",
        payload.get("training_matrix_authority_earned"),
        False,
    )
    checks = require_mapping(f"{row_id}.checks", payload.get("checks"))
    require_exact_keys(
        f"{row_id}.checks", checks, CAPACITY_AUDIT_CHECK_KEYS
    )
    for check_name in CAPACITY_AUDIT_CHECK_KEYS:
        require_bool(
            f"{row_id}.checks.{check_name}",
            checks.get(check_name),
            True,
        )
    audit_path = require_absolute_file(
        f"{row_id} capacity audit artifact",
        artifact["path"],
    )
    expected_matrix_path = audit_path.parent / (
        f"{expected_system}_h70_source_complete_smoke_matrix.npz"
    )
    matrix_out = Path(
        require_nonempty_text(
            f"{row_id}.matrix_out", payload.get("matrix_out")
        )
    )
    if (
        matrix_out != expected_matrix_path
        or expected_matrix_path.exists()
        or expected_matrix_path.is_symlink()
    ):
        raise AmendmentError(
            f"{row_id} capacity audit matrix-absence contract differs"
        )


def validate_capacity_evidence(
    resource: Mapping[str, Any],
    root_join: Mapping[str, Any],
    audits: Mapping[str, Mapping[str, Any]],
    artifacts: Mapping[str, Mapping[str, Any]],
    corrected_sources: Sequence[Mapping[str, Any]],
    legacy_preflight: Mapping[str, Any],
    corrected_preflight: Mapping[str, Any],
    immutable: Mapping[str, Any],
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    require_exact_keys("capacity resource certificate", resource, CAPACITY_KEYS)
    require_exact_keys("capacity root/join certificate", root_join, ROOT_JOIN_KEYS)
    if resource.get("schema") != CAPACITY_SCHEMA or resource.get("status") != "PASS":
        raise AmendmentError("capacity resource certificate schema/status differs")
    if (
        root_join.get("schema") != ROOT_JOIN_SCHEMA
        or root_join.get("status") != "PASS"
        or root_join.get("scope") != "capacity"
    ):
        raise AmendmentError("capacity root/join schema/status/scope differs")
    for label, payload in audits.items():
        require_exact_keys(f"{label} capacity audit", payload, AUDIT_KEYS)
        if (
            payload.get("schema") != CAPACITY_AUDIT_SCHEMA
            or payload.get("status") != "PASS"
        ):
            raise AmendmentError(f"{label} capacity audit schema/status differs")
        audit_artifact = artifacts[
            "pp_audit" if label.startswith("pp_") else "auau_audit"
        ]
        validate_capacity_non_training_audit(
            label, payload, audit_artifact
        )

    if tuple(resource.get("selected_rows", [])) != SELECTED_CAPACITY_ROWS:
        raise AmendmentError("capacity selected-row order/identity differs")
    for label, payload in (
        ("capacity resource", resource),
        ("capacity root/join", root_join),
        *[(f"{row_id} audit", audits[row_id]) for row_id in SELECTED_CAPACITY_ROWS],
    ):
        require_exact_int(
            f"{label}.execution_group_size",
            payload.get("execution_group_size"),
            GROUP_SIZE,
        )
        require_exact_int(
            f"{label}.source_occurrences_per_output",
            payload.get("source_occurrences_per_output"),
            1,
        )
        require_exact_int(
            f"{label}.full_training_authority",
            payload.get("full_training_authority"),
            0,
        )
    require_bool(
        "capacity resource.capacity_authority_earned",
        resource.get("capacity_authority_earned"),
        True,
    )
    require_bool(
        "capacity resource.submission_performed",
        resource.get("submission_performed"),
        True,
    )
    require_bool(
        "capacity root/join.capacity_authority_earned",
        root_join.get("capacity_authority_earned"),
        True,
    )
    legacy_bindings = {
        "full_plan_sha256": legacy_preflight["plan"]["sha256"],
        "preflight_receipt_sha256": legacy_preflight[
            "preflight_receipt"
        ]["sha256"],
        "bundle_manifest_sha256": legacy_preflight[
            "bundle_manifest_sha256"
        ],
        "materialization_receipt_sha256": legacy_preflight[
            "materialization_receipt_sha256"
        ],
        "execution_partition_sha256": legacy_preflight[
            "execution_partition_sha256"
        ],
    }
    for field, expected in legacy_bindings.items():
        if resource.get(field) != expected:
            raise AmendmentError(f"capacity resource {field} binding differs")

    if (
        resource.get("root_health_identity_join_certificate_sha256")
        != artifacts["root_join_certificate"]["sha256"]
    ):
        raise AmendmentError("capacity resource root/join binding differs")
    require_same_resolved_file(
        "capacity resource root/join binding",
        resource.get("root_health_identity_join_certificate"),
        artifacts["root_join_certificate"]["path"],
    )
    resource_audits = require_sequence(
        "capacity resource.multiview_audits",
        resource.get("multiview_audits"),
    )
    resource_audit_by_row: dict[str, dict[str, Any]] = {}
    for raw in resource_audits:
        item = require_mapping("capacity resource audit", raw)
        row_id = require_nonempty_text(
            "capacity resource audit row_id", item.get("row_id")
        )
        if row_id in resource_audit_by_row:
            raise AmendmentError(
                f"capacity resource duplicates audit row {row_id}"
            )
        resource_audit_by_row[row_id] = item
    if list(resource_audit_by_row) != list(SELECTED_CAPACITY_ROWS):
        raise AmendmentError("capacity resource audit order/identity differs")
    for row_id in SELECTED_CAPACITY_ROWS:
        expected_artifact = artifacts[
            "pp_audit" if row_id.startswith("pp_") else "auau_audit"
        ]
        resource_audit = resource_audit_by_row[row_id]
        require_same_resolved_file(
            f"capacity resource audit binding {row_id}",
            resource_audit.get("audit"),
            expected_artifact["path"],
        )
        if (
            resource_audit.get("audit_sha256") != expected_artifact["sha256"]
            or resource_audit.get("population_state")
            != EXPECTED_POPULATION_STATE[row_id]
            or resource_audit.get("source_occurrence_id_hex")
            != audits[row_id].get("source_occurrence_id_hex")
        ):
            raise AmendmentError(
                f"capacity resource audit binding differs for {row_id}"
            )

    root_rows = require_sequence("capacity root/join.rows", root_join.get("rows"))
    root_by_id: dict[str, dict[str, Any]] = {}
    for raw in root_rows:
        row = require_mapping("capacity root/join row", raw)
        require_exact_keys(
            f"capacity root/join row {row.get('row_id')}",
            row,
            ROOT_ROW_KEYS,
        )
        row_id = str(row.get("row_id", ""))
        if row_id in root_by_id:
            raise AmendmentError(f"capacity root/join duplicates {row_id}")
        root_by_id[row_id] = row
    if list(root_by_id) != list(SELECTED_CAPACITY_ROWS):
        raise AmendmentError("capacity root/join row order/identity differs")
    if (
        root_join.get("populated_rows") != ["auau_background_jet12"]
        or root_join.get("valid_empty_rows") != ["pp_background_jet8"]
    ):
        raise AmendmentError("capacity populated/valid-empty matrix differs")
    require_exact_int("capacity root row_count", root_join.get("row_count"), 2)
    require_exact_int(
        "capacity root populated_row_count",
        root_join.get("populated_row_count"),
        1,
    )
    require_exact_int(
        "capacity root valid_empty_row_count",
        root_join.get("valid_empty_row_count"),
        1,
    )
    if root_join.get("failures") != []:
        raise AmendmentError("capacity root/join failures are not empty")

    resource_rows = require_sequence("capacity resource.rows", resource.get("rows"))
    resource_by_id: dict[str, dict[str, Any]] = {}
    for raw in resource_rows:
        row = require_mapping("capacity resource row", raw)
        row_id = str(row.get("row_id", ""))
        if row_id in resource_by_id:
            raise AmendmentError(f"capacity resource duplicates {row_id}")
        resource_by_id[row_id] = validate_capacity_resource_row(row_id, row)
    if list(resource_by_id) != list(SELECTED_CAPACITY_ROWS):
        raise AmendmentError("capacity resource row order/identity differs")

    validation_authority = require_mapping(
        "capacity resource.validation_authority",
        resource.get("validation_authority"),
    )
    if set(validation_authority) != {
        "controller",
        "runtime_authority",
        "submission_journal",
        "submission_manifest",
        "submission_receipt",
    }:
        raise AmendmentError("capacity validation-authority inventory differs")
    validation_artifacts = {
        "controller": verify_nested_file_reference(
            "capacity controller",
            validation_authority["controller"],
            allowed_extra=("validation_commit",),
        ),
        "runtime_authority": verify_nested_file_reference(
            "capacity runtime authority",
            validation_authority["runtime_authority"],
        ),
        "submission_journal": verify_nested_file_reference(
            "capacity submission journal",
            validation_authority["submission_journal"],
        ),
        "submission_manifest": verify_nested_file_reference(
            "capacity submission manifest",
            validation_authority["submission_manifest"],
        ),
        "submission_receipt": verify_nested_file_reference(
            "capacity submission receipt",
            validation_authority["submission_receipt"],
        ),
    }
    root_validation = require_mapping(
        "capacity root/join.validation_authority",
        root_join.get("validation_authority"),
    )
    require_exact_keys(
        "capacity root/join.validation_authority",
        root_validation,
        VALIDATION_AUTHORITY_OUTPUT_KEYS,
    )
    for field, authority_reference in validation_artifacts.items():
        root_reference = verify_nested_file_reference(
            f"capacity root/join {field}",
            root_validation.get(field),
            allowed_extra=("validation_commit",)
            if field == "controller"
            else (),
        )
        if root_reference["sha256"] != authority_reference["sha256"]:
            raise AmendmentError(
                f"capacity root/join validation authority differs: {field}"
            )
        require_same_resolved_file(
            f"capacity root/join validation authority {field}",
            root_reference["path"],
            authority_reference["path"],
        )
        if (
            field == "controller"
            and root_reference["validation_commit"]
            != authority_reference["validation_commit"]
        ):
            raise AmendmentError(
                "capacity root/join validation commit differs"
            )
    runtime_authority_binding = validate_runtime_authority_artifact(
        validation_artifacts["runtime_authority"],
        immutable,
    )
    for row_id in SELECTED_CAPACITY_ROWS:
        audit = audits[row_id]
        for field in ("submission_manifest", "submission_receipt"):
            audit_reference = verify_nested_file_reference(
                f"{row_id} audit {field}", audit.get(field)
            )
            authority_reference = validation_artifacts[field]
            if audit_reference["sha256"] != authority_reference["sha256"]:
                raise AmendmentError(
                    f"{row_id} audit {field} hash authority differs"
                )
            require_same_resolved_file(
                f"{row_id} audit {field} authority",
                audit_reference["path"],
                authority_reference["path"],
            )

    submission_manifest_path = Path(
        validation_artifacts["submission_manifest"]["path"]
    )
    manifest_header, manifest_rows = read_tsv(submission_manifest_path)
    if tuple(manifest_header) != SUBMISSION_MANIFEST_FIELDS:
        raise AmendmentError("capacity submission-manifest columns differ")
    manifest_by_id: dict[str, dict[str, str]] = {}
    for row in manifest_rows:
        row_id = require_nonempty_text(
            "capacity submission manifest row_id", row.get("row_id")
        )
        if row_id in manifest_by_id:
            raise AmendmentError(
                f"capacity submission manifest duplicates {row_id}"
            )
        manifest_by_id[row_id] = row
    if list(manifest_by_id) != expected_row_ids():
        raise AmendmentError("capacity submission-manifest row closure differs")

    receipt_path = Path(validation_artifacts["submission_receipt"]["path"])
    header, receipt_rows = read_tsv(receipt_path)
    if tuple(header) != SUBMISSION_RECEIPT_FIELDS:
        raise AmendmentError("capacity submission-receipt columns differ")
    receipt_by_id: dict[str, dict[str, str]] = {}
    for row in receipt_rows:
        row_id = require_nonempty_text(
            "capacity submission receipt row_id", row.get("row_id")
        )
        if row_id in receipt_by_id:
            raise AmendmentError(
                f"capacity submission receipt duplicates {row_id}"
            )
        receipt_by_id[row_id] = row
    if list(receipt_by_id) != list(SELECTED_CAPACITY_ROWS):
        raise AmendmentError("capacity submission-receipt row closure differs")

    corrected_by_id = {source["row_id"]: source for source in corrected_sources}
    corrected_semantics_by_id: dict[str, dict[str, Any]] = {}
    for raw in require_sequence(
        "corrected capacity source semantics",
        corrected_preflight.get("capacity_source_semantics"),
    ):
        record = require_mapping("corrected capacity source semantics", raw)
        row_id = require_nonempty_text(
            "corrected capacity source semantics row_id",
            record.get("row_id"),
        )
        if row_id in corrected_semantics_by_id:
            raise AmendmentError(
                f"corrected capacity source semantics duplicates {row_id}"
            )
        corrected_semantics_by_id[row_id] = record
    if list(corrected_semantics_by_id) != expected_row_ids():
        raise AmendmentError(
            "corrected capacity source-semantics order/identity differs"
        )
    capacity_rows: list[dict[str, Any]] = []
    for row_id in SELECTED_CAPACITY_ROWS:
        audit = audits[row_id]
        root_row = root_by_id[row_id]
        receipt_row = receipt_by_id[row_id]
        if audit.get("row_id") != row_id or audit.get("system") != EXPECTED_SYSTEM[row_id]:
            raise AmendmentError(f"{row_id} capacity audit identity differs")
        if audit.get("population_state") != EXPECTED_POPULATION_STATE[row_id]:
            raise AmendmentError(f"{row_id} capacity population state differs")
        if audit.get("failures") != []:
            raise AmendmentError(f"{row_id} capacity audit failures are not empty")
        if root_row.get("status") != "PASS":
            raise AmendmentError(f"{row_id} root/join row is not PASS")
        if root_row.get("population_state") != EXPECTED_POPULATION_STATE[row_id]:
            raise AmendmentError(f"{row_id} root population state differs")
        for field in (
            "source_contract",
            "source_execution_contract",
            "source_identity_canonical_sha256",
            "source_occurrence_id_hex",
        ):
            if audit.get(field) != root_row.get(field):
                raise AmendmentError(
                    f"{row_id} audit/root source binding differs: {field}"
                )
        source_contract = require_mapping(
            f"{row_id}.source_contract", audit.get("source_contract")
        )
        source_execution = require_mapping(
            f"{row_id}.source_execution_contract",
            audit.get("source_execution_contract"),
        )
        corrected_semantics = corrected_semantics_by_id[row_id]
        validate_capacity_source_binding(
            row_id,
            source_contract,
            source_execution,
            corrected_semantics,
        )
        expected_identity = canonical_source_identity_sha256(source_contract)
        if audit.get("source_identity_canonical_sha256") != expected_identity:
            raise AmendmentError(f"{row_id} canonical source identity differs")
        occurrence_id = expected_identity[:32]
        if occurrence_id == "0" * 32:
            occurrence_id = "0" * 31 + "1"
        if audit.get("source_occurrence_id_hex") != occurrence_id:
            raise AmendmentError(f"{row_id} source occurrence ID differs")

        corrected_source = corrected_by_id[row_id]
        submission_manifest_row = manifest_by_id[row_id]
        submission_manifest_binding = validate_submission_manifest_record(
            row_id,
            submission_manifest_row,
            source_contract=source_contract,
            corrected_semantics=corrected_semantics,
            corrected_source=corrected_source,
            immutable=immutable,
            runtime_authority=runtime_authority_binding,
        )
        materialized_config = normalize_artifact_ref(
            f"{row_id}.receipt.materialized_config",
            {
                "path": receipt_row.get("materialized_config"),
                "sha256": receipt_row.get("materialized_config_sha256"),
            },
        )
        fanout_contract = normalize_artifact_ref(
            f"{row_id}.receipt.fanout_contract",
            {
                "path": receipt_row.get("fanout_contract_file"),
                "sha256": receipt_row.get("fanout_contract_sha256"),
            },
        )
        base_config = require_mapping(
            f"immutable {corrected_semantics['system']} config",
            require_mapping(
                "immutable bundle artifacts",
                immutable.get("bundle_artifacts"),
            ).get(f"{corrected_semantics['system']}_config"),
        )
        config_chain_sha256 = validate_runtime_config_chain(
            row_id,
            str(corrected_semantics["system"]),
            base_config=base_config,
            resolved_config=submission_manifest_binding["runtime_config"],
            materialized_config=materialized_config,
            fanout_contract=fanout_contract,
            analysis_output_namespace=submission_manifest_row[
                "analysis_output_namespace"
            ],
        )
        snapshot_loader, snapshot_manifest = (
            validate_frozen_snapshot_receipts(
                row_id,
                str(corrected_semantics["system"]),
                receipt_row,
                submission_manifest_row,
                immutable,
                runtime_authority_binding,
            )
        )
        if (
            root_row.get("snapshot_loader_receipt_sha256")
            != snapshot_loader["sha256"]
            or root_row.get("snapshot_manifest_sha256")
            != snapshot_manifest["sha256"]
        ):
            raise AmendmentError(
                f"{row_id} root/snapshot receipt authority differs"
            )
        expected_chunk_bytes, corrected_chunk = materialize_first_chunk_bytes(
            corrected_source
        )
        staged_path = require_absolute_file(
            f"{row_id}.staged_chunk_list",
            source_execution.get("staged_chunk_list"),
        )
        staged_sha = sha256_file(staged_path)
        expected_staged_sha = hashlib.sha256(expected_chunk_bytes).hexdigest()
        if staged_path.read_bytes() != expected_chunk_bytes:
            raise AmendmentError(
                f"{row_id} staged capacity chunk differs from corrected chunk 0"
            )
        if staged_sha != expected_staged_sha:
            raise AmendmentError(
                f"{row_id} staged capacity chunk hash reconstruction differs"
            )
        if (
            source_contract.get("input_uri_hash") != staged_sha
            or source_contract.get("input_file_sha256") != staged_sha
            or receipt_row.get("staged_chunk_sha256") != staged_sha
        ):
            raise AmendmentError(
                f"{row_id} staged capacity chunk receipt/source binding differs"
            )
        require_same_resolved_file(
            f"{row_id} staged capacity chunk receipt",
            staged_path,
            receipt_row.get("staged_chunk_list"),
        )
        require_exact_int(
            f"{row_id}.source_execution.chunk_index",
            source_execution.get("chunk_index"),
            1,
        )
        require_exact_int(
            f"{row_id}.source_contract.segment",
            source_contract.get("segment"),
            1,
        )
        args_path = require_absolute_file(
            f"{row_id}.source_execution.args_file",
            source_execution.get("args_file"),
        )
        submitted_args_sha256 = sha256_file(args_path)
        if (
            require_sha256(
                f"{row_id}.source_execution.submitted_args_sha256",
                source_execution.get("submitted_args_sha256"),
            )
            != submitted_args_sha256
            or require_sha256(
                f"{row_id}.receipt.submitted_args_sha256",
                receipt_row.get("submitted_args_sha256"),
            )
            != submitted_args_sha256
            or receipt_row.get("cluster_proc")
            != resource_by_id[row_id].get("cluster_proc")
        ):
            raise AmendmentError(f"{row_id} args/cluster binding differs")
        require_same_resolved_file(
            f"{row_id} args receipt",
            args_path,
            receipt_row.get("args_file"),
        )

        analysis_path = require_absolute_file(
            f"{row_id}.analysis_output_root",
            root_row.get("analysis_output_root"),
        )
        sidecar_path = require_absolute_file(
            f"{row_id}.sidecar", root_row.get("sidecar")
        )
        require_same_resolved_file(
            f"{row_id} analysis output ownership",
            receipt_row.get("analysis_output_root"),
            analysis_path,
        )
        require_same_resolved_file(
            f"{row_id} sidecar output ownership",
            receipt_row.get("multiview_sidecar"),
            sidecar_path,
        )
        require_same_resolved_file(
            f"{row_id} submission-manifest sidecar binding",
            submission_manifest_row.get("multiview_sidecar"),
            sidecar_path,
        )
        if (
            analysis_path.stat().st_size != root_row.get("analysis_bytes")
            or sha256_file(analysis_path) != root_row.get("analysis_sha256")
            or sidecar_path.stat().st_size != root_row.get("sidecar_bytes")
            or sha256_file(sidecar_path) != root_row.get("sidecar_sha256")
        ):
            raise AmendmentError(f"{row_id} capacity output hash/size drift")
        source_provenance = require_mapping(
            f"{row_id}.source_provenance", audit.get("source_provenance")
        )
        provenance_artifact = validate_source_provenance_artifact(
            row_id,
            source_provenance,
            sidecar_path=sidecar_path,
            staged_sha256=staged_sha,
            source_contract=source_contract,
            corrected_semantics=corrected_semantics,
            submission_manifest_row=submission_manifest_row,
        )

        capacity_rows.append(
            {
                "row_id": row_id,
                "system": EXPECTED_SYSTEM[row_id],
                "population_state": EXPECTED_POPULATION_STATE[row_id],
                "cluster_proc": resource_by_id[row_id]["cluster_proc"],
                "exit_code": 0,
                "num_job_starts": 1,
                "num_holds": 0,
                "request_memory_mb": REQUEST_MEMORY_MB,
                "memory_usage_mb": resource_by_id[row_id]["memory_usage_mb"],
                "resident_set_size_kb": resource_by_id[row_id][
                    "resident_set_size_kb"
                ],
                "remote_wall_clock_seconds": resource_by_id[row_id][
                    "remote_wall_clock_seconds"
                ],
                "capacity_chunk_index": 1,
                "corrected_partition_chunk_index": 0,
                "segment": 1,
                "tuple_count": GROUP_SIZE,
                "staged_chunk_list": str(staged_path),
                "staged_chunk_sha256": staged_sha,
                "corrected_chunk_fingerprint_sha256": corrected_chunk[
                    "chunk_fingerprint_sha256"
                ],
                "corrected_execution_chunk_sha256": corrected_preflight[
                    "rows_proof"
                ][expected_row_ids().index(row_id)][
                    "first_execution_chunk_sha256"
                ],
                "full_source_manifest_sha256": corrected_source[
                    "full_source_manifest_sha256"
                ],
                "source_contract": source_contract,
                "source_identity_canonical_sha256": expected_identity,
                "source_occurrence_id_hex": occurrence_id,
                "base_config_sha256": corrected_semantics[
                    "base_config_sha256"
                ],
                "capacity_runtime_config": pinned_artifact_record(
                    f"{row_id}_capacity_runtime_config",
                    submission_manifest_binding["runtime_config"],
                ),
                "capacity_materialized_config": pinned_artifact_record(
                    f"{row_id}_capacity_materialized_config",
                    materialized_config,
                ),
                "capacity_fanout_contract": pinned_artifact_record(
                    f"{row_id}_capacity_fanout_contract",
                    fanout_contract,
                ),
                "runtime_config_chain_semantic_sha256":
                    config_chain_sha256,
                "snapshot_loader_receipt": snapshot_loader,
                "snapshot_manifest": snapshot_manifest,
                "capacity_submission_manifest_row_sha256":
                    submission_manifest_binding["row_semantic_sha256"],
                "analysis_output": {
                    "path": str(analysis_path),
                    "sha256": root_row["analysis_sha256"],
                    "size_bytes": root_row["analysis_bytes"],
                },
                "sidecar_output": {
                    "path": str(sidecar_path),
                    "sha256": root_row["sidecar_sha256"],
                    "size_bytes": root_row["sidecar_bytes"],
                },
                "source_provenance": pinned_artifact_record(
                    f"{row_id}_source_provenance", provenance_artifact
                ),
                "staged_chunk_equals_corrected_partition_chunk0": True,
            }
        )

    capacity_artifacts = {
        "resource_certificate": pinned_artifact_record(
            "capacity_resource_certificate",
            artifacts["resource_certificate"],
            schema=CAPACITY_SCHEMA,
        ),
        "root_join_certificate": pinned_artifact_record(
            "capacity_root_join_certificate",
            artifacts["root_join_certificate"],
            schema=ROOT_JOIN_SCHEMA,
        ),
        "pp_audit": pinned_artifact_record(
            "pp_capacity_audit",
            artifacts["pp_audit"],
            schema=CAPACITY_AUDIT_SCHEMA,
        ),
        "auau_audit": pinned_artifact_record(
            "auau_capacity_audit",
            artifacts["auau_audit"],
            schema=CAPACITY_AUDIT_SCHEMA,
        ),
        "validation_authority": validation_artifacts,
    }
    return capacity_artifacts, capacity_rows


def load_source_manifest_artifact(
    label: str,
    raw: object,
    *,
    expected_schema: str,
) -> tuple[dict[str, Any], dict[str, Any]]:
    artifact = normalize_artifact_ref(label, raw)
    try:
        payload = source_builder.load_manifest(
            Path(artifact["path"]), artifact["sha256"]
        )
    except (source_builder.ManifestError, OSError) as exc:
        raise AmendmentError(f"{label} validation failed: {exc}") from exc
    if payload.get("schema") != expected_schema:
        raise AmendmentError(
            f"{label} schema must be {expected_schema}, "
            f"observed={payload.get('schema')!r}"
        )
    return payload, artifact


def validate_pinned_output_artifact(
    label: str,
    raw: object,
    *,
    expected_schema: str | None = None,
    extra_keys: Iterable[str] = (),
) -> dict[str, Any]:
    payload = require_mapping(label, raw)
    expected_keys = set(PINNED_ARTIFACT_OUTPUT_KEYS) | set(extra_keys)
    if expected_schema is not None:
        expected_keys.add("schema")
    require_exact_keys(label, payload, expected_keys)
    require_nonempty_text(f"{label}.label", payload.get("label"))
    for field in ("path", "resolved_path"):
        path = Path(require_nonempty_text(f"{label}.{field}", payload.get(field)))
        if not path.is_absolute():
            raise AmendmentError(f"{label}.{field} must be absolute")
    require_sha256(f"{label}.sha256", payload.get("sha256"))
    size_bytes = require_nonnegative_int(
        f"{label}.size_bytes", payload.get("size_bytes")
    )
    if size_bytes == 0:
        raise AmendmentError(f"{label}.size_bytes must be positive")
    if expected_schema is not None and payload.get("schema") != expected_schema:
        raise AmendmentError(f"{label}.schema differs")
    for field in extra_keys:
        require_nonempty_text(f"{label}.{field}", payload.get(field))
    return payload


def validate_nested_amendment_evidence(payload: Mapping[str, Any]) -> None:
    immutable = require_mapping(
        "amendment.immutable_authority",
        payload.get("immutable_authority"),
    )
    require_exact_keys(
        "amendment.immutable_authority",
        immutable,
        IMMUTABLE_AUTHORITY_OUTPUT_KEYS,
    )
    if immutable.get("status") != "PASS":
        raise AmendmentError("amendment immutable authority is not PASS")
    validate_pinned_output_artifact(
        "amendment immutable bundle",
        immutable.get("bundle_manifest"),
        expected_schema=resolver.BUNDLE_SCHEMA,
    )
    validate_pinned_output_artifact(
        "amendment immutable materialization",
        immutable.get("materialization_receipt"),
        expected_schema=resolver.MATERIALIZATION_SCHEMA,
    )
    for field in (
        "bundle_identity_sha256",
        "bundle_semantic_fingerprint_sha256",
        "code_sha256",
        "replay_schema_sha256",
        "training_schema_sha256",
        "semantic_sha256",
    ):
        require_sha256(f"amendment immutable {field}", immutable.get(field))
    immutable_runtime = require_mapping(
        "amendment immutable runtime", immutable.get("runtime")
    )
    require_exact_keys(
        "amendment immutable runtime",
        immutable_runtime,
        IMMUTABLE_RUNTIME_OUTPUT_KEYS,
    )
    if (
        immutable_runtime.get("release") != "ana.560"
        or immutable_runtime.get("calo_reco_soname")
        != "libcalo_reco.so.0"
    ):
        raise AmendmentError("amendment immutable runtime identity differs")
    require_exact_int(
        "amendment immutable runtime request_memory_mb",
        immutable_runtime.get("request_memory_mb"),
        REQUEST_MEMORY_MB,
    )
    for field in (
        "offline_main",
        "release_core_lib_dir",
        "release_core_lib64_dir",
    ):
        path = Path(
            require_nonempty_text(
                f"amendment immutable runtime {field}",
                immutable_runtime.get(field),
            )
        )
        if not path.is_absolute():
            raise AmendmentError(
                f"amendment immutable runtime {field} must be absolute"
            )
    immutable_artifacts = require_mapping(
        "amendment immutable bundle artifacts",
        immutable.get("bundle_artifacts"),
    )
    require_exact_keys(
        "amendment immutable bundle artifacts",
        immutable_artifacts,
        resolver.REQUIRED_ARTIFACT_ROLES,
    )
    for role, raw in immutable_artifacts.items():
        artifact = require_mapping(
            f"amendment immutable bundle artifact {role}", raw
        )
        require_exact_keys(
            f"amendment immutable bundle artifact {role}",
            artifact,
            BUNDLE_ARTIFACT_OUTPUT_KEYS,
        )
        if artifact.get("role") != role:
            raise AmendmentError(
                f"amendment immutable bundle artifact role differs: {role}"
            )
        require_sha256(
            f"amendment immutable bundle artifact {role} sha256",
            artifact.get("sha256"),
        )
        if require_nonnegative_int(
            f"amendment immutable bundle artifact {role} size_bytes",
            artifact.get("size_bytes"),
        ) <= 0:
            raise AmendmentError(
                f"amendment immutable bundle artifact {role} is empty"
            )
        for field in ("path", "resolved_path"):
            path = Path(
                require_nonempty_text(
                    f"amendment immutable bundle artifact {role} {field}",
                    artifact.get(field),
                )
            )
            if not path.is_absolute():
                raise AmendmentError(
                    f"amendment immutable bundle artifact {role} "
                    f"{field} must be absolute"
                )
    if not GIT_COMMIT_RE.fullmatch(str(immutable.get("public_commit", ""))):
        raise AmendmentError("amendment immutable public_commit differs")
    require_nonempty_text(
        "amendment immutable materialized_bundle_path",
        immutable.get("materialized_bundle_path"),
    )
    if (
        require_nonnegative_int(
            "amendment immutable artifact_count",
            immutable.get("artifact_count"),
        )
        <= 0
        or require_nonnegative_int(
            "amendment immutable total_artifact_bytes",
            immutable.get("total_artifact_bytes"),
        )
        <= 0
        or immutable.get("readback")
        != "PASS_SYMLINK_FREE_READONLY_CONTENT_EXACT"
    ):
        raise AmendmentError("amendment immutable readback evidence differs")

    legacy = require_mapping(
        "amendment.legacy_preflight", payload.get("legacy_preflight")
    )
    require_exact_keys(
        "amendment.legacy_preflight",
        legacy,
        LEGACY_PREFLIGHT_OUTPUT_KEYS,
    )
    if legacy.get("status") != "PASS":
        raise AmendmentError("amendment legacy preflight is not PASS")
    require_bool(
        "amendment legacy submission_performed",
        legacy.get("submission_performed"),
        False,
    )
    require_exact_int(
        "amendment legacy full_training_authority",
        legacy.get("full_training_authority"),
        0,
    )
    for field, schema in (
        ("plan", LEGACY_PLAN_SCHEMA),
        ("preflight_receipt", LEGACY_RECEIPT_SCHEMA),
        ("source_manifest", source_builder.SOURCE_SCHEMA_V1),
    ):
        validate_pinned_output_artifact(
            f"amendment legacy {field}",
            legacy.get(field),
            expected_schema=schema,
        )
    for field in (
        "execution_partition_sha256",
        "bundle_manifest_sha256",
        "materialization_receipt_sha256",
        "source_manifest_sha256",
    ):
        require_sha256(f"amendment legacy {field}", legacy.get(field))
    for field, expected in (
        ("source_tuple_count", EXPECTED_SOURCE_TUPLES),
        ("expected_job_count", EXPECTED_CHUNKS),
        ("expected_output_count", EXPECTED_CHUNKS),
        ("legacy_expected_occurrence_count", EXPECTED_SOURCE_TUPLES),
    ):
        require_exact_int(f"amendment legacy {field}", legacy.get(field), expected)
    legacy_rows = require_sequence(
        "amendment legacy rows", legacy.get("rows")
    )
    if [row.get("row_id") for row in legacy_rows if isinstance(row, dict)] != expected_row_ids():
        raise AmendmentError("amendment legacy row order/identity differs")
    for raw in legacy_rows:
        row = require_mapping("amendment legacy row", raw)
        require_exact_keys(
            f"amendment legacy row {row.get('row_id')}",
            row,
            LEGACY_ROW_OUTPUT_KEYS,
        )
    if (
        sum(row["source_tuple_count"] for row in legacy_rows)
        != EXPECTED_SOURCE_TUPLES
        or sum(row["expected_chunk_count"] for row in legacy_rows)
        != EXPECTED_CHUNKS
        or sum(
            row["legacy_expected_occurrence_count"]
            for row in legacy_rows
        )
        != EXPECTED_SOURCE_TUPLES
    ):
        raise AmendmentError("amendment legacy rowwise totals differ")

    corrected = require_mapping(
        "amendment.corrected_preflight",
        payload.get("corrected_preflight"),
    )
    require_exact_keys(
        "amendment.corrected_preflight",
        corrected,
        CORRECTED_PREFLIGHT_OUTPUT_KEYS,
    )
    if corrected.get("status") != "PASS":
        raise AmendmentError("amendment corrected preflight is not PASS")
    require_bool(
        "amendment corrected submission_performed",
        corrected.get("submission_performed"),
        False,
    )
    require_exact_int(
        "amendment corrected full_training_authority",
        corrected.get("full_training_authority"),
        0,
    )
    for field, schema in (
        ("plan", CORRECTED_PLAN_SCHEMA),
        ("preflight_receipt", CORRECTED_RECEIPT_SCHEMA),
        ("source_manifest", CORRECTED_SOURCE_SCHEMA),
    ):
        validate_pinned_output_artifact(
            f"amendment corrected {field}",
            corrected.get(field),
            expected_schema=schema,
        )
    for field in ("rows", "duplicate_fingerprint"):
        validate_pinned_output_artifact(
            f"amendment corrected {field}", corrected.get(field)
        )
    partition = validate_pinned_output_artifact(
        "amendment corrected partition",
        corrected.get("partition"),
        extra_keys=("record_count",),
    )
    require_exact_int(
        "amendment corrected partition record_count",
        partition.get("record_count"),
        EXPECTED_CHUNKS,
    )
    for field in (
        "bundle_manifest_sha256",
        "materialization_receipt_sha256",
        "source_manifest_sha256",
        "duplicate_fingerprint_sha256",
        "execution_partition_sha256",
        "execution_fingerprint_sha256",
        "science_contract_semantic_sha256",
    ):
        require_sha256(f"amendment corrected {field}", corrected.get(field))
    corrected_counts = {
        "source_tuple_count": EXPECTED_SOURCE_TUPLES,
        "expected_chunk_count": EXPECTED_CHUNKS,
        "expected_job_count": EXPECTED_CHUNKS,
        "expected_output_pair_count": EXPECTED_CHUNKS,
        "expected_analysis_output_count": EXPECTED_CHUNKS,
        "expected_sidecar_output_count": EXPECTED_CHUNKS,
        "expected_physical_root_artifact_count": EXPECTED_ROOT_ARTIFACTS,
        "expected_source_occurrence_count": EXPECTED_CHUNKS,
        "source_occurrences_per_output_pair": 1,
    }
    for field, expected in corrected_counts.items():
        require_exact_int(
            f"amendment corrected {field}", corrected.get(field), expected
        )
    row_proofs = require_sequence(
        "amendment corrected rows_proof", corrected.get("rows_proof")
    )
    if [row.get("row_id") for row in row_proofs if isinstance(row, dict)] != expected_row_ids():
        raise AmendmentError("amendment corrected row proof order differs")
    for raw in row_proofs:
        row = require_mapping("amendment corrected row proof", raw)
        require_exact_keys(
            f"amendment corrected row proof {row.get('row_id')}",
            row,
            ROW_PROOF_OUTPUT_KEYS,
        )
    if (
        sum(row["source_tuple_count"] for row in row_proofs)
        != EXPECTED_SOURCE_TUPLES
        or sum(row["expected_chunk_count"] for row in row_proofs)
        != EXPECTED_CHUNKS
        or sum(row["expected_job_count"] for row in row_proofs)
        != EXPECTED_CHUNKS
        or sum(row["expected_output_pair_count"] for row in row_proofs)
        != EXPECTED_CHUNKS
        or sum(
            row["expected_source_occurrence_count"]
            for row in row_proofs
        )
        != EXPECTED_CHUNKS
    ):
        raise AmendmentError("amendment corrected rowwise totals differ")
    source_semantics = require_sequence(
        "amendment corrected capacity_source_semantics",
        corrected.get("capacity_source_semantics"),
    )
    if [
        row.get("row_id")
        for row in source_semantics
        if isinstance(row, dict)
    ] != expected_row_ids():
        raise AmendmentError(
            "amendment corrected source-semantics order differs"
        )
    for raw in source_semantics:
        row = require_mapping("amendment corrected source semantics", raw)
        require_exact_keys(
            f"amendment corrected source semantics {row.get('row_id')}",
            row,
            CAPACITY_SOURCE_SEMANTIC_OUTPUT_KEYS,
        )
        require_sha256(
            f"{row.get('row_id')} corrected code_sha256",
            row.get("code_sha256"),
        )
        require_sha256(
            f"{row.get('row_id')} corrected base_config_sha256",
            row.get("base_config_sha256"),
        )
        for field in (
            "replay_schema_sha256",
            "training_schema_sha256",
            "semantic_sha256",
            "library_sha256",
            "model_sha256",
        ):
            require_sha256(
                f"{row.get('row_id')} corrected {field}",
                row.get(field),
            )
    bundle_semantics = require_mapping(
        "amendment corrected bundle semantics",
        corrected.get("bundle_contract_semantic_sha256_by_system"),
    )
    require_exact_keys(
        "amendment corrected bundle semantics",
        bundle_semantics,
        {"pp", "auau"},
    )
    for system in ("pp", "auau"):
        require_sha256(
            f"amendment corrected bundle semantics {system}",
            bundle_semantics.get(system),
        )

    crosswalk = require_mapping(
        "amendment.source_semantic_crosswalk",
        payload.get("source_semantic_crosswalk"),
    )
    require_exact_keys(
        "amendment.source_semantic_crosswalk",
        crosswalk,
        SOURCE_CROSSWALK_OUTPUT_KEYS,
    )
    if crosswalk.get("status") != "PASS":
        raise AmendmentError("amendment source crosswalk is not PASS")
    require_exact_int("amendment crosswalk row_count", crosswalk.get("row_count"), 13)
    require_exact_int(
        "amendment crosswalk source_tuple_count",
        crosswalk.get("source_tuple_count"),
        EXPECTED_SOURCE_TUPLES,
    )
    require_exact_int(
        "amendment crosswalk global_unique_tuple_count",
        crosswalk.get("global_unique_tuple_count"),
        EXPECTED_SOURCE_TUPLES,
    )
    for field in (
        "source_rows_semantic_sha256",
        "manifest_semantic_sha256",
        "global_tuple_identity_records_sha256",
    ):
        require_sha256(f"amendment crosswalk {field}", crosswalk.get(field))
    crosswalk_rows = require_sequence(
        "amendment crosswalk rows", crosswalk.get("rows")
    )
    if [row.get("row_id") for row in crosswalk_rows if isinstance(row, dict)] != expected_row_ids():
        raise AmendmentError("amendment source crosswalk row order differs")
    for raw in crosswalk_rows:
        row = require_mapping("amendment source crosswalk row", raw)
        require_exact_keys(
            f"amendment source crosswalk row {row.get('row_id')}",
            row,
            SOURCE_CROSSWALK_ROW_KEYS,
        )
        require_bool(
            f"{row.get('row_id')} legacy alias history flag",
            row.get("legacy_tuple_aliases_preserved_as_history_only"),
            True,
        )
        require_bool(
            f"{row.get('row_id')} corrected alias absence flag",
            row.get("corrected_execution_aliases_absent"),
            True,
        )
        input_lists = require_sequence(
            f"{row.get('row_id')} crosswalk input_lists",
            row.get("input_lists"),
        )
        if [
            record.get("role")
            for record in input_lists
            if isinstance(record, dict)
        ] != list(resolver.LIST_ROLES):
            raise AmendmentError(
                f"{row.get('row_id')} crosswalk input-list order differs"
            )
        for raw_list in input_lists:
            record = require_mapping(
                f"{row.get('row_id')} crosswalk input list", raw_list
            )
            require_exact_keys(
                f"{row.get('row_id')} crosswalk input list "
                f"{record.get('role')}",
                record,
                CROSSWALK_INPUT_LIST_KEYS,
            )
            require_sha256(
                f"{row.get('row_id')} {record.get('role')} SHA-256",
                record.get("sha256"),
            )
            for field in (
                "size_bytes",
                "line_count",
                "executable_count",
            ):
                value = require_nonnegative_int(
                    f"{row.get('row_id')} {record.get('role')} {field}",
                    record.get(field),
                )
                if value <= 0:
                    raise AmendmentError(
                        f"{row.get('row_id')} {record.get('role')} {field} "
                        "must be positive"
                    )
                if field == "executable_count" and value != row["tuple_count"]:
                    raise AmendmentError(
                        f"{row.get('row_id')} {record.get('role')} {field} "
                        "differs from tuple_count"
                    )
    if (
        sum(row["tuple_count"] for row in crosswalk_rows)
        != EXPECTED_SOURCE_TUPLES
    ):
        raise AmendmentError("amendment source crosswalk rowwise total differs")

    capacity = require_mapping(
        "amendment.capacity_evidence", payload.get("capacity_evidence")
    )
    require_exact_keys(
        "amendment.capacity_evidence",
        capacity,
        CAPACITY_EVIDENCE_OUTPUT_KEYS,
    )
    for field, schema in (
        ("resource_certificate", CAPACITY_SCHEMA),
        ("root_join_certificate", ROOT_JOIN_SCHEMA),
        ("pp_audit", CAPACITY_AUDIT_SCHEMA),
        ("auau_audit", CAPACITY_AUDIT_SCHEMA),
    ):
        validate_pinned_output_artifact(
            f"amendment capacity {field}",
            capacity.get(field),
            expected_schema=schema,
        )
    validation_authority = require_mapping(
        "amendment capacity validation_authority",
        capacity.get("validation_authority"),
    )
    require_exact_keys(
        "amendment capacity validation_authority",
        validation_authority,
        VALIDATION_AUTHORITY_OUTPUT_KEYS,
    )
    for field in VALIDATION_AUTHORITY_OUTPUT_KEYS:
        validate_pinned_output_artifact(
            f"amendment capacity validation {field}",
            validation_authority.get(field),
            extra_keys=("validation_commit",) if field == "controller" else (),
        )


def validate_amendment_payload(payload: Mapping[str, Any]) -> dict[str, Any]:
    expected_keys = {
        "schema",
        "status",
        "authority_state",
        "generated_by",
        "immutable_authority",
        "legacy_preflight",
        "corrected_preflight",
        "source_semantic_crosswalk",
        "capacity_evidence",
        "count_correction",
        "capacity_to_partition_binding",
        "checks",
        "submission_performed",
        "jobs_rerun",
        "science_rerun_required",
        "capacity_reuse_scope",
        "full_training_authority",
        "full_extraction_authority",
        "broad_production_authority",
        "canonical_promotion",
        "amendment_semantic_sha256",
    }
    require_exact_keys("capacity count amendment", payload, expected_keys)
    if payload.get("schema") != AMENDMENT_SCHEMA or payload.get("status") != "PASS":
        raise AmendmentError("amendment schema/status differs")
    if payload.get("authority_state") != AUTHORITY_STATE:
        raise AmendmentError("amendment authority state differs")
    if payload.get("generated_by") != GENERATED_BY:
        raise AmendmentError("amendment generator identity differs")
    validate_nested_amendment_evidence(payload)
    immutable = require_mapping(
        "amendment immutable authority binding",
        payload.get("immutable_authority"),
    )
    require_bool("amendment.submission_performed", payload.get("submission_performed"), False)
    require_exact_int("amendment.jobs_rerun", payload.get("jobs_rerun"), 0)
    require_bool(
        "amendment.science_rerun_required",
        payload.get("science_rerun_required"),
        False,
    )
    if (
        payload.get("capacity_reuse_scope")
        != "group_size_7_resource_runtime_only"
    ):
        raise AmendmentError("amendment capacity reuse scope differs")
    require_exact_int(
        "amendment.full_training_authority",
        payload.get("full_training_authority"),
        0,
    )
    for field in (
        "full_extraction_authority",
        "broad_production_authority",
        "canonical_promotion",
    ):
        require_bool(f"amendment.{field}", payload.get(field), False)
    checks = require_mapping("amendment.checks", payload.get("checks"))
    require_exact_keys("amendment.checks", checks, CHECK_KEYS)
    if any(value is not True for value in checks.values()):
        raise AmendmentError("amendment checks are incomplete or not PASS")
    count = require_mapping(
        "amendment.count_correction", payload.get("count_correction")
    )
    require_exact_keys(
        "amendment.count_correction", count, COUNT_CORRECTION_KEYS
    )
    if count.get("basis") != (
        "ordered_disjoint_group_of_seven_partition_one_execution_and_"
        "source_occurrence_per_output_pair"
    ):
        raise AmendmentError("amendment count-correction basis differs")
    expected_count_values = {
        "source_tuple_count": EXPECTED_SOURCE_TUPLES,
        "corrected_chunk_count": EXPECTED_CHUNKS,
        "corrected_job_count": EXPECTED_CHUNKS,
        "corrected_output_pair_count": EXPECTED_CHUNKS,
        "corrected_analysis_output_count": EXPECTED_CHUNKS,
        "corrected_sidecar_output_count": EXPECTED_CHUNKS,
        "corrected_source_occurrence_count": EXPECTED_CHUNKS,
        "physical_root_artifact_count": EXPECTED_ROOT_ARTIFACTS,
        "source_occurrences_per_output_pair": 1,
        "legacy_occurrence_overstatement": LEGACY_OVERSTATEMENT,
        "group_size": GROUP_SIZE,
    }
    for field, expected in expected_count_values.items():
        require_exact_int(
            f"amendment.count_correction.{field}",
            count.get(field),
            expected,
        )
    require_exact_int(
        "amendment.count_correction.legacy_expected_occurrence_count",
        count.get("legacy_expected_occurrence_count"),
        EXPECTED_SOURCE_TUPLES,
    )
    require_bool(
        "amendment.count_correction.science_semantics_changed",
        count.get("science_semantics_changed"),
        False,
    )
    require_bool(
        "amendment.count_correction.tolerances_changed",
        count.get("tolerances_changed"),
        False,
    )
    capacity_binding = require_mapping(
        "amendment.capacity_to_partition_binding",
        payload.get("capacity_to_partition_binding"),
    )
    require_exact_keys(
        "amendment.capacity_to_partition_binding",
        capacity_binding,
        CAPACITY_BINDING_KEYS,
    )
    if (
        capacity_binding.get("status") != "PASS"
        or capacity_binding.get("reuse_scope")
        != "group_size_7_resource_runtime_only"
    ):
        raise AmendmentError("amendment capacity binding status/scope differs")
    require_sha256(
        "amendment corrected execution-partition SHA-256",
        capacity_binding.get("corrected_execution_partition_sha256"),
    )
    require_sha256(
        "amendment corrected partition-artifact SHA-256",
        capacity_binding.get("corrected_partition_artifact_sha256"),
    )
    selected_rows = require_sequence(
        "amendment capacity selected_rows",
        capacity_binding.get("selected_rows"),
    )
    if (
        len(selected_rows) != len(SELECTED_CAPACITY_ROWS)
        or [
            row.get("row_id") for row in selected_rows
            if isinstance(row, dict)
        ]
        != list(SELECTED_CAPACITY_ROWS)
    ):
        raise AmendmentError("amendment capacity selected-row binding differs")
    corrected = require_mapping(
        "amendment corrected preflight binding",
        payload.get("corrected_preflight"),
    )
    if (
        capacity_binding.get("corrected_execution_partition_sha256")
        != corrected.get("execution_partition_sha256")
        or capacity_binding.get("corrected_partition_artifact_sha256")
        != require_mapping(
            "amendment corrected partition binding",
            corrected.get("partition"),
        ).get("sha256")
    ):
        raise AmendmentError(
            "amendment capacity hashes differ from corrected preflight"
        )
    corrected_semantics_by_id = {
        row["row_id"]: row
        for row in require_sequence(
            "amendment corrected source semantics binding",
            corrected.get("capacity_source_semantics"),
        )
    }
    crosswalk_by_id = {
        row["row_id"]: row
        for row in require_sequence(
            "amendment source crosswalk binding",
            require_mapping(
                "amendment source crosswalk binding",
                payload.get("source_semantic_crosswalk"),
            ).get("rows"),
        )
    }
    for raw in selected_rows:
        row = require_mapping("amendment capacity selected row", raw)
        require_exact_keys(
            f"amendment capacity selected row {row.get('row_id')}",
            row,
            CAPACITY_ROW_OUTPUT_KEYS,
        )
        row_id = str(row["row_id"])
        require_exact_int(f"{row_id}.exit_code", row.get("exit_code"), 0)
        require_exact_int(
            f"{row_id}.num_job_starts", row.get("num_job_starts"), 1
        )
        require_exact_int(f"{row_id}.num_holds", row.get("num_holds"), 0)
        require_exact_int(
            f"{row_id}.request_memory_mb",
            row.get("request_memory_mb"),
            REQUEST_MEMORY_MB,
        )
        memory_usage_mb = row.get("memory_usage_mb")
        resident_set_size_kb = row.get("resident_set_size_kb")
        if (
            isinstance(memory_usage_mb, bool)
            or not isinstance(memory_usage_mb, int)
            or memory_usage_mb <= 0
            or memory_usage_mb > REQUEST_MEMORY_MB
            or isinstance(resident_set_size_kb, bool)
            or not isinstance(resident_set_size_kb, int)
            or resident_set_size_kb <= 0
            or resident_set_size_kb > REQUEST_MEMORY_MB * 1024
        ):
            raise AmendmentError(
                f"{row_id} amendment resource use exceeds its request"
            )
        for field, expected in (
            ("capacity_chunk_index", 1),
            ("corrected_partition_chunk_index", 0),
            ("segment", 1),
            ("tuple_count", GROUP_SIZE),
        ):
            require_exact_int(f"{row_id}.{field}", row.get(field), expected)
        for field in (
            "staged_chunk_sha256",
            "corrected_chunk_fingerprint_sha256",
            "corrected_execution_chunk_sha256",
            "full_source_manifest_sha256",
            "source_identity_canonical_sha256",
            "base_config_sha256",
            "runtime_config_chain_semantic_sha256",
            "capacity_submission_manifest_row_sha256",
        ):
            require_sha256(f"{row_id}.{field}", row.get(field))
        if row.get("base_config_sha256") != corrected_semantics_by_id[
            row_id
        ]["base_config_sha256"]:
            raise AmendmentError(
                f"{row_id} amended base-config binding differs"
            )
        for field in (
            "replay_schema_sha256",
            "training_schema_sha256",
            "semantic_sha256",
            "library_sha256",
            "model_sha256",
        ):
            if field in {
                "replay_schema_sha256",
                "training_schema_sha256",
                "semantic_sha256",
            }:
                expected_authority = immutable[field]
            else:
                artifact_role = (
                    f"{corrected_semantics_by_id[row_id]['system']}_"
                    f"{'library' if field == 'library_sha256' else 'model'}"
                )
                expected_authority = require_mapping(
                    "amendment immutable bundle artifacts",
                    immutable.get("bundle_artifacts"),
                )[artifact_role]["sha256"]
            if corrected_semantics_by_id[row_id][field] != expected_authority:
                raise AmendmentError(
                    f"{row_id} corrected immutable {field} binding differs"
                )
        validate_pinned_output_artifact(
            f"{row_id}.capacity_runtime_config",
            row.get("capacity_runtime_config"),
        )
        validate_pinned_output_artifact(
            f"{row_id}.capacity_materialized_config",
            row.get("capacity_materialized_config"),
        )
        validate_pinned_output_artifact(
            f"{row_id}.capacity_fanout_contract",
            row.get("capacity_fanout_contract"),
        )
        validate_pinned_output_artifact(
            f"{row_id}.snapshot_loader_receipt",
            row.get("snapshot_loader_receipt"),
            expected_schema=SNAPSHOT_LOADER_SCHEMA,
        )
        validate_pinned_output_artifact(
            f"{row_id}.snapshot_manifest",
            row.get("snapshot_manifest"),
            expected_schema=SNAPSHOT_MANIFEST_SCHEMA,
        )
        source_contract = require_mapping(
            f"{row_id}.source_contract", row.get("source_contract")
        )
        require_exact_keys(
            f"{row_id}.source_contract",
            source_contract,
            SOURCE_CONTRACT_KEYS,
        )
        for field in (
            "input_uri_hash",
            "input_file_sha256",
            "source_manifest_sha256",
        ):
            require_sha256(f"{row_id}.source_contract.{field}", source_contract.get(field))
        corrected_semantics = corrected_semantics_by_id[row_id]
        for field in (
            "lane",
            "dataset",
            "sample",
            "period",
            "run",
            "si_di_role",
            "ownership_state",
        ):
            if source_contract.get(field) != corrected_semantics[field]:
                raise AmendmentError(
                    f"{row_id} amended source contract differs from corrected "
                    f"semantics: {field}"
                )
        if (
            row.get("system") != corrected_semantics["system"]
            or row.get("full_source_manifest_sha256")
            != crosswalk_by_id[row_id]["full_source_manifest_sha256"]
            or row.get("staged_chunk_sha256")
            != source_contract["input_uri_hash"]
            or row.get("staged_chunk_sha256")
            != source_contract["input_file_sha256"]
        ):
            raise AmendmentError(
                f"{row_id} amended source/output binding differs"
            )
        expected_identity = canonical_source_identity_sha256(source_contract)
        expected_occurrence = expected_identity[:32]
        if expected_occurrence == "0" * 32:
            expected_occurrence = "0" * 31 + "1"
        if (
            row.get("source_identity_canonical_sha256") != expected_identity
            or row.get("source_occurrence_id_hex") != expected_occurrence
        ):
            raise AmendmentError(
                f"{row_id} amended source identity binding differs"
            )
        for field in ("analysis_output", "sidecar_output"):
            artifact = require_mapping(
                f"{row_id}.{field}", row.get(field)
            )
            require_exact_keys(
                f"{row_id}.{field}", artifact, OUTPUT_ARTIFACT_KEYS
            )
            require_sha256(f"{row_id}.{field}.sha256", artifact.get("sha256"))
            if require_nonnegative_int(
                f"{row_id}.{field}.size_bytes", artifact.get("size_bytes")
            ) <= 0:
                raise AmendmentError(f"{row_id}.{field} must be nonempty")
        validate_pinned_output_artifact(
            f"{row_id}.source_provenance",
            row.get("source_provenance"),
        )
        require_bool(
            f"{row_id}.staged_chunk_equals_corrected_partition_chunk0",
            row.get("staged_chunk_equals_corrected_partition_chunk0"),
            True,
        )
    require_bool(
        "amendment capacity_authority_earned",
        capacity_binding.get("capacity_authority_earned"),
        True,
    )
    require_exact_int(
        "amendment capacity full_training_authority",
        capacity_binding.get("full_training_authority"),
        0,
    )
    require_bool(
        "amendment capacity full_extraction_authority",
        capacity_binding.get("full_extraction_authority"),
        False,
    )
    require_bool(
        "amendment capacity broad_production_authority",
        capacity_binding.get("broad_production_authority"),
        False,
    )
    legacy = require_mapping(
        "amendment legacy preflight binding",
        payload.get("legacy_preflight"),
    )
    if not (
        legacy.get("bundle_manifest_sha256")
        == corrected.get("bundle_manifest_sha256")
        == require_mapping(
            "amendment immutable bundle binding",
            immutable.get("bundle_manifest"),
        ).get("sha256")
        and legacy.get("materialization_receipt_sha256")
        == corrected.get("materialization_receipt_sha256")
        == require_mapping(
            "amendment immutable materialization binding",
            immutable.get("materialization_receipt"),
        ).get("sha256")
    ):
        raise AmendmentError(
            "amendment preflights differ from immutable authority"
        )
    fingerprint_payload = {
        key: value
        for key, value in payload.items()
        if key != "amendment_semantic_sha256"
    }
    expected_semantic_sha = semantic_sha256(fingerprint_payload)
    if require_sha256(
        "amendment semantic SHA-256",
        payload.get("amendment_semantic_sha256"),
    ) != expected_semantic_sha:
        raise AmendmentError("amendment semantic SHA-256 differs")
    return dict(payload)


def build_amendment(spec: Mapping[str, Any]) -> dict[str, Any]:
    require_exact_keys("amendment spec", spec, SPEC_KEYS)
    if spec.get("schema") != SPEC_SCHEMA:
        raise AmendmentError(f"amendment spec schema must be {SPEC_SCHEMA}")
    immutable_spec = require_mapping(
        "spec.immutable_authority", spec["immutable_authority"]
    )
    legacy_spec = require_mapping("spec.legacy", spec["legacy"])
    corrected_spec = require_mapping("spec.corrected", spec["corrected"])
    capacity_spec = require_mapping("spec.capacity", spec["capacity"])
    require_exact_keys(
        "spec.immutable_authority",
        immutable_spec,
        IMMUTABLE_SPEC_KEYS,
    )
    require_exact_keys("spec.legacy", legacy_spec, LEGACY_SPEC_KEYS)
    require_exact_keys("spec.corrected", corrected_spec, CORRECTED_SPEC_KEYS)
    require_exact_keys("spec.capacity", capacity_spec, CAPACITY_SPEC_KEYS)

    immutable_authority = validate_immutable_authority(immutable_spec)

    legacy_plan, legacy_plan_artifact = load_json_artifact(
        "legacy plan",
        legacy_spec["plan"],
        expected_schema=LEGACY_PLAN_SCHEMA,
    )
    legacy_receipt, legacy_receipt_artifact = load_json_artifact(
        "legacy preflight receipt",
        legacy_spec["preflight_receipt"],
        expected_schema=LEGACY_RECEIPT_SCHEMA,
    )
    legacy_source, legacy_source_artifact = load_source_manifest_artifact(
        "legacy source manifest",
        legacy_spec["source_manifest"],
        expected_schema=source_builder.SOURCE_SCHEMA_V1,
    )
    legacy_artifacts = {
        "plan": legacy_plan_artifact,
        "preflight_receipt": legacy_receipt_artifact,
        "source_manifest": legacy_source_artifact,
    }
    legacy_preflight = validate_legacy_preflight(
        legacy_plan,
        legacy_receipt,
        legacy_source,
        legacy_artifacts,
    )
    validate_preflight_immutable_binding(
        "legacy plan",
        legacy_plan,
        legacy_receipt,
        immutable_authority,
    )

    corrected_plan, corrected_plan_artifact = load_json_artifact(
        "corrected plan",
        corrected_spec["plan"],
        expected_schema=CORRECTED_PLAN_SCHEMA,
    )
    corrected_receipt, corrected_receipt_artifact = load_json_artifact(
        "corrected preflight receipt",
        corrected_spec["preflight_receipt"],
        expected_schema=CORRECTED_RECEIPT_SCHEMA,
    )
    corrected_source, corrected_source_artifact = load_source_manifest_artifact(
        "corrected source manifest",
        corrected_spec["source_manifest"],
        expected_schema=CORRECTED_SOURCE_SCHEMA,
    )
    try:
        corrected_source_records = resolver.validate_source_manifest(
            corrected_source
        )
    except resolver.ControllerError as exc:
        raise AmendmentError(
            f"corrected source partition reconstruction failed: {exc}"
        ) from exc
    corrected_rows, corrected_rows_artifact, _rows_bytes = load_jsonl_artifact(
        "corrected rows", corrected_spec["rows"]
    )
    duplicate_text, corrected_duplicate_artifact = load_text_artifact(
        "corrected duplicate fingerprint",
        corrected_spec["duplicate_fingerprint"],
    )
    (
        corrected_partition_records,
        corrected_partition_artifact,
        corrected_partition_bytes,
    ) = load_jsonl_artifact(
        "corrected partition", corrected_spec["partition"]
    )
    corrected_artifacts = {
        "plan": corrected_plan_artifact,
        "preflight_receipt": corrected_receipt_artifact,
        "source_manifest": corrected_source_artifact,
        "rows": corrected_rows_artifact,
        "duplicate_fingerprint": corrected_duplicate_artifact,
        "partition": corrected_partition_artifact,
    }
    corrected_preflight = validate_corrected_preflight(
        corrected_plan,
        corrected_receipt,
        corrected_source,
        corrected_source_records,
        corrected_rows,
        duplicate_text,
        corrected_partition_records,
        corrected_partition_bytes,
        corrected_artifacts,
    )
    validate_preflight_immutable_binding(
        "corrected plan",
        corrected_plan,
        corrected_receipt,
        immutable_authority,
    )
    validate_corrected_bundle_authority(
        corrected_plan,
        immutable_authority,
    )
    source_crosswalk = compare_source_manifests(
        legacy_source, corrected_source
    )

    if (
        legacy_preflight["bundle_manifest_sha256"]
        != corrected_preflight["bundle_manifest_sha256"]
        or legacy_preflight["materialization_receipt_sha256"]
        != corrected_preflight["materialization_receipt_sha256"]
    ):
        raise AmendmentError(
            "V1/V2 immutable bundle or materialization authority differs"
        )
    legacy_rows = require_sequence("legacy plan.rows", legacy_plan["rows"])
    corrected_rows_by_id = {
        row["row_id"]: row for row in corrected_plan["rows"]
    }
    for old in legacy_rows:
        row_id = old["row_id"]
        new = corrected_rows_by_id[row_id]
        for field in (
            "system",
            "lane",
            "dataset",
            "sample",
            "source_role",
            "minimum_bias_gate",
            "photon_id_row_match",
            "source_period",
            "source_si_di_role",
            "training_period_si_contract_sha256",
            "science_contract",
            "bundle_contract",
        ):
            if old.get(field) != new.get(field):
                raise AmendmentError(
                    f"{row_id} V1/V2 scientific or bundle field differs: {field}"
                )

    capacity_resource, capacity_resource_artifact = load_json_artifact(
        "capacity resource certificate",
        capacity_spec["resource_certificate"],
        expected_schema=CAPACITY_SCHEMA,
    )
    capacity_root, capacity_root_artifact = load_json_artifact(
        "capacity root/join certificate",
        capacity_spec["root_join_certificate"],
        expected_schema=ROOT_JOIN_SCHEMA,
    )
    pp_audit, pp_audit_artifact = load_json_artifact(
        "pp capacity audit",
        capacity_spec["pp_audit"],
        expected_schema=CAPACITY_AUDIT_SCHEMA,
    )
    auau_audit, auau_audit_artifact = load_json_artifact(
        "auau capacity audit",
        capacity_spec["auau_audit"],
        expected_schema=CAPACITY_AUDIT_SCHEMA,
    )
    capacity_artifact_inputs = {
        "resource_certificate": capacity_resource_artifact,
        "root_join_certificate": capacity_root_artifact,
        "pp_audit": pp_audit_artifact,
        "auau_audit": auau_audit_artifact,
    }
    capacity_artifacts, capacity_rows = validate_capacity_evidence(
        capacity_resource,
        capacity_root,
        {
            "pp_background_jet8": pp_audit,
            "auau_background_jet12": auau_audit,
        },
        capacity_artifact_inputs,
        corrected_source_records,
        legacy_preflight,
        corrected_preflight,
        immutable_authority,
    )

    count_correction = {
        "basis": (
            "ordered_disjoint_group_of_seven_partition_one_execution_and_"
            "source_occurrence_per_output_pair"
        ),
        "group_size": GROUP_SIZE,
        "source_tuple_count": EXPECTED_SOURCE_TUPLES,
        "legacy_expected_occurrence_count": EXPECTED_SOURCE_TUPLES,
        "corrected_chunk_count": EXPECTED_CHUNKS,
        "corrected_job_count": EXPECTED_CHUNKS,
        "corrected_output_pair_count": EXPECTED_CHUNKS,
        "corrected_analysis_output_count": EXPECTED_CHUNKS,
        "corrected_sidecar_output_count": EXPECTED_CHUNKS,
        "physical_root_artifact_count": EXPECTED_ROOT_ARTIFACTS,
        "corrected_source_occurrence_count": EXPECTED_CHUNKS,
        "source_occurrences_per_output_pair": 1,
        "legacy_occurrence_overstatement": LEGACY_OVERSTATEMENT,
        "science_semantics_changed": False,
        "tolerances_changed": False,
    }
    capacity_to_partition = {
        "status": "PASS",
        "reuse_scope": "group_size_7_resource_runtime_only",
        "corrected_execution_partition_sha256": corrected_preflight[
            "execution_partition_sha256"
        ],
        "corrected_partition_artifact_sha256": corrected_preflight[
            "partition"
        ]["sha256"],
        "selected_rows": capacity_rows,
        "capacity_authority_earned": True,
        "full_training_authority": 0,
        "full_extraction_authority": False,
        "broad_production_authority": False,
    }
    payload: dict[str, Any] = {
        "schema": AMENDMENT_SCHEMA,
        "status": "PASS",
        "authority_state": AUTHORITY_STATE,
        "generated_by": GENERATED_BY,
        "immutable_authority": immutable_authority,
        "legacy_preflight": legacy_preflight,
        "corrected_preflight": corrected_preflight,
        "source_semantic_crosswalk": source_crosswalk,
        "capacity_evidence": capacity_artifacts,
        "count_correction": count_correction,
        "capacity_to_partition_binding": capacity_to_partition,
        "checks": {
            "legacy_v1_evidence_rehashed": True,
            "legacy_tuple_aliases_preserved_as_history_only": True,
            "corrected_v2_evidence_rehashed": True,
            "source_semantics_identical_v1_to_v2": True,
            "rowwise_partition_reconstructed": True,
            "partition_jsonl_byte_exact": True,
            "aggregate_counts_exact": True,
            "capacity_certificate_chain_rehashed": True,
            "capacity_rows_terminal_single_start_no_hold": True,
            "capacity_resource_request_8000_mb": True,
            "valid_empty_and_populated_witnesses_present": True,
            "capacity_chunks_equal_corrected_chunk0": True,
            "resolved_and_materialized_config_chain_revalidated": True,
            "frozen_snapshot_inventory_and_loader_revalidated": True,
            "manifest_runtime_identities_bound_to_immutable_authority": True,
            "audit_validation_authority_cross_bound": True,
            "output_hashes_revalidated": True,
            "no_science_change": True,
            "no_tolerance_change": True,
            "no_submission_or_job_control": True,
        },
        "submission_performed": False,
        "jobs_rerun": 0,
        "science_rerun_required": False,
        "capacity_reuse_scope": "group_size_7_resource_runtime_only",
        "full_training_authority": 0,
        "full_extraction_authority": False,
        "broad_production_authority": False,
        "canonical_promotion": False,
    }
    payload["amendment_semantic_sha256"] = semantic_sha256(payload)
    return validate_amendment_payload(payload)


def write_once_readonly_json(
    path: Path,
    payload: Mapping[str, Any],
    *,
    allowed_output_root: Path | None = None,
) -> str:
    if not path.is_absolute():
        raise AmendmentError(f"output path must be absolute: {path}")
    validate_output_parent(
        "output path",
        path.parent,
        allowed_output_root=allowed_output_root,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    validate_output_parent(
        "output path",
        path.parent,
        allowed_output_root=allowed_output_root,
    )
    content = canonical_json_bytes(payload)
    if path.exists() or path.is_symlink():
        if path.is_symlink() or not path.is_file():
            raise AmendmentError(
                f"existing output is not a regular non-symlink file: {path}"
            )
        if path.lstat().st_nlink != 1:
            raise AmendmentError(
                f"existing output must have exactly one hard link: {path}"
            )
        if path.read_bytes() != content:
            raise AmendmentError(
                f"existing amendment differs; refusing overwrite: {path}"
            )
        os.chmod(path, stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH)
        return "EXISTING_BYTE_IDENTICAL"
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    if temporary.exists() or temporary.is_symlink():
        raise AmendmentError(f"temporary output path already exists: {temporary}")
    try:
        with temporary.open("xb") as stream:
            stream.write(content)
            stream.flush()
            os.fsync(stream.fileno())
        os.chmod(temporary, stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH)
        try:
            os.link(temporary, path, follow_symlinks=False)
        except FileExistsError:
            if path.is_symlink() or not path.is_file():
                raise AmendmentError(
                    "concurrent output is not a regular non-symlink file: "
                    f"{path}"
                )
            if path.lstat().st_nlink != 1:
                raise AmendmentError(
                    "concurrent output must have exactly one hard link: "
                    f"{path}"
                )
            if path.read_bytes() != content:
                raise AmendmentError(
                    "concurrent amendment differs; refusing overwrite: "
                    f"{path}"
                )
            os.chmod(path, stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH)
            return "EXISTING_BYTE_IDENTICAL_RACE"
        directory_fd = os.open(path.parent, os.O_RDONLY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
    finally:
        if temporary.exists():
            temporary.unlink()
    return "CREATED_READONLY"


def load_amendment(path: Path, expected_sha256: str) -> dict[str, Any]:
    payload, _artifact = load_json_artifact(
        "capacity count amendment",
        {"path": str(path), "sha256": expected_sha256},
        expected_schema=AMENDMENT_SCHEMA,
    )
    return validate_amendment_payload(payload)


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="action", required=True)
    build = subparsers.add_parser("build", help="build the immutable amendment")
    build.add_argument("--spec", type=Path, required=True)
    build.add_argument("--expected-spec-sha256", required=True)
    build.add_argument("--output", type=Path, required=True)
    build.add_argument("--allowed-output-root", type=Path, required=True)
    verify = subparsers.add_parser(
        "verify", help="rehash all evidence and byte-compare the amendment"
    )
    verify.add_argument("--spec", type=Path, required=True)
    verify.add_argument("--expected-spec-sha256", required=True)
    verify.add_argument("--amendment", type=Path, required=True)
    verify.add_argument("--expected-amendment-sha256", required=True)
    verify.add_argument("--output-json", type=Path)
    verify.add_argument("--allowed-output-root", type=Path)
    return parser.parse_args(argv)


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        spec = load_spec(args.spec, args.expected_spec_sha256)
        rebuilt = build_amendment(spec)
        if args.action == "build":
            disposition = write_once_readonly_json(
                args.output,
                rebuilt,
                allowed_output_root=args.allowed_output_root,
            )
            summary = {
                "schema": AMENDMENT_SCHEMA,
                "status": "PASS",
                "authority_state": AUTHORITY_STATE,
                "path": str(args.output),
                "sha256": sha256_file(args.output),
                "amendment_semantic_sha256": rebuilt[
                    "amendment_semantic_sha256"
                ],
                "write_disposition": disposition,
                "submission_performed": False,
                "full_training_authority": 0,
                "full_extraction_authority": False,
            }
        else:
            observed = load_amendment(
                args.amendment, args.expected_amendment_sha256
            )
            rebuilt_bytes = canonical_json_bytes(rebuilt)
            try:
                observed_bytes = args.amendment.read_bytes()
            except OSError as exc:
                raise AmendmentError(
                    f"amendment readback cannot be read: {args.amendment}"
                ) from exc
            if (
                observed_bytes != rebuilt_bytes
                or canonical_json_bytes(observed) != rebuilt_bytes
            ):
                raise AmendmentError(
                    "amendment readback differs from current pinned evidence"
                )
            summary = {
                "schema": READBACK_SCHEMA,
                "status": "PASS",
                "authority_state": AUTHORITY_STATE,
                "amendment": str(args.amendment),
                "amendment_sha256": args.expected_amendment_sha256,
                "amendment_semantic_sha256": observed[
                    "amendment_semantic_sha256"
                ],
                "all_evidence_rehashed": True,
                "byte_exact_rebuild": True,
                "submission_performed": False,
                "full_training_authority": 0,
                "full_extraction_authority": False,
            }
            if args.output_json is not None:
                if args.allowed_output_root is None:
                    raise AmendmentError(
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
        AmendmentError,
        source_builder.ManifestError,
        resolver.ControllerError,
        OSError,
        KeyError,
        IndexError,
        TypeError,
        UnicodeDecodeError,
        csv.Error,
        StopIteration,
    ) as exc:
        print(f"THE134_CAPACITY_COUNT_AMENDMENT_ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
