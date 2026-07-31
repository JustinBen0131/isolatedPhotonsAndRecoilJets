#!/usr/bin/env python3
"""Guard the source-complete THE-134 group-of-seven extraction lifecycle.

This controller is intentionally a thin authority layer over the already
certified resolver, immutable-bundle materializer, pre-extraction storage
projector, and RecoilJets submitter.  It does not reconstruct detector or
training semantics, query Condor, control jobs, or physically merge ROOT
files.

The lifecycle is:

``preflight``
    Revalidate every pinned upstream receipt and byte-compare the existing
    dry-materialization directory with the deterministic materializer output.

``seal``
    Consume a separately produced, exact-scope authorization receipt and write
    one immutable execution manifest.  This is the only transition that grants
    this controller permission to call the pinned submitter.

``submit``
    With an additional explicit ``--execute`` switch, invoke the exact pinned
    submitter once for each of the thirteen frozen source rows.  Before the
    first call, the controller reopens every sealed prerequisite and
    revalidates the unexpired authorization.  It inherits only a narrow safe
    operating-system environment allowlist, overlays the sealed row
    environment, requires one exact Condor cluster/job-count report per row,
    and writes a durable partial or complete receipt.  It never releases,
    removes, retries, or resubmits a row.

``terminal``
    Validate an externally captured, hash-pinned terminal snapshot containing
    only this campaign's exact jobs and artifact-health results.  No scheduler
    query occurs here.

``aggregate``
    Produce exact-input-once sidecar lists and source-provenance JSON for the
    downstream THE-134 matrix builders.  This is a logical merge.  Physical
    ``hadd`` is deliberately forbidden because every training sidecar owns one
    source-occurrence identity and per-file metadata.

None of these actions grants model, working-point, science-freeze, THE-121,
THE-122, CANONICAL, or physics-output authority.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import re
import socket
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path, PurePosixPath
from typing import Any, Callable, Iterable, Mapping, Sequence


SHARED_DIRECTORY_MODE = 0o2755
SHARED_FILE_MODE = 0o644


HERE = Path(__file__).resolve().parent
MATERIALIZER_PATH = HERE / "materialize_the134_full_multiview_extraction.py"
CONTRACT_PATH = HERE.parents[2] / "ml" / "contracts" / "the134_h70_contract.py"
SAFE_RESUME_GATE_PATH = (
    HERE.parents[1] / "runtime" / "audit" / "sdcc_safe_resume_gate.py"
)


def _load_local_module(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load local module: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


materializer = _load_local_module(
    "the134_full_materializer_for_execution", MATERIALIZER_PATH
)
h70_contract = _load_local_module("the134_h70_contract_for_execution", CONTRACT_PATH)
safe_resume_gate = _load_local_module(
    "sdcc_safe_resume_gate_for_execution", SAFE_RESUME_GATE_PATH
)


AUTHORIZATION_SCHEMA = "THE134_FULL_EXTRACTION_SUBMISSION_AUTHORIZATION_V1"
EXECUTION_SCHEMA = "THE134_FULL_EXTRACTION_EXECUTION_MANIFEST_V1"
SUBMISSION_SCHEMA = "THE134_FULL_EXTRACTION_SUBMISSION_RECEIPT_V1"
TERMINAL_SNAPSHOT_SCHEMA = "THE134_FULL_EXTRACTION_TERMINAL_SNAPSHOT_V1"
TERMINAL_RECEIPT_SCHEMA = "THE134_FULL_EXTRACTION_TERMINAL_RECEIPT_V1"
AGGREGATE_SCHEMA = "THE134_FULL_EXTRACTION_LOGICAL_AGGREGATE_V1"

READY_STATUS = "READY_TO_SUBMIT_EXACT_GROUP7"
SUBMITTED_STATUS = "SUBMITTED_DURABLE_EXACT_SOURCE_COMPLETE"
TERMINAL_STATUS = "PASS_EXACT_TERMINAL_AND_ARTIFACT_HEALTH"
AGGREGATE_STATUS = "PASS_EXACT_INPUT_ONCE_LOGICAL_AGGREGATE"

EXPECTED_ROW_COUNT = materializer.EXPECTED_ROW_COUNT
EXPECTED_JOB_COUNT = materializer.EXPECTED_JOB_COUNT
EXPECTED_GROUP_SIZE = materializer.EXPECTED_GROUP_SIZE
EXPECTED_REQUEST_MEMORY_MB = materializer.EXPECTED_REQUEST_MEMORY_MB
EXPECTED_OUTPUT_PAIRS = materializer.EXPECTED_OUTPUT_PAIRS
ARTIFACT_PROFILE = materializer.ARTIFACT_PROFILE
SHOWER_VIEWS = tuple(materializer.resolver.SHOWER_VIEWS)

SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
IDENTITY128_RE = re.compile(r"^[0-9a-f]{32}$")
SAFE_LABEL_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.:/@+-]{0,255}$")
SAFE_SDCC_USER_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,63}$")
SUBMIT_REPORT_RE = re.compile(
    r"(?P<count>[0-9]+)\s+job\(s\)\s+submitted\s+to\s+cluster\s+"
    r"(?P<cluster>[0-9]+)",
    re.IGNORECASE,
)
SDCC_USER_ALIAS_ROOT = Path("/sphenix/u")
SDCC_CANONICAL_USER_ROOT = Path("/gpfs/mnt/gpfs02/sphenix/user")

SOURCE_INPUT_RECORD_FIELDS = (
    "path",
    "system",
    "source_sample",
    "source_occurrence_id_hex",
    "input_uri_sha256",
    "input_file_sha256",
    "source_manifest_sha256",
    "config_sha256",
    "code_sha256",
    "training_view_root_sha256",
)

FORBIDDEN_AUTHORITY_TRUE = frozenset(
    {
        "full_training_authority",
        "full_extraction_authority",
        "science_freeze_authority",
        "broad_production_authority",
        "physics_output_authority",
        "canonical_promotion",
        "the121_authority",
        "the122_authority",
    }
)
SAFE_SUBMIT_AMBIENT_KEYS = frozenset(
    {
        "PATH",
        "LD_LIBRARY_PATH",
        "HOME",
        "USER",
        "LOGNAME",
        "SHELL",
        "TMPDIR",
        "TMP",
        "TEMP",
        "LANG",
        "LC_ALL",
        "LC_CTYPE",
        "TZ",
        "TERM",
        "X509_USER_PROXY",
        "KRB5CCNAME",
        "CONDOR_CONFIG",
    }
)
SUBMITTER_TIMEOUT_SECONDS = 600
SUBMITTER_OUTPUT_LIMIT_BYTES = 256 * 1024


class ControllerError(RuntimeError):
    """Fail-closed lifecycle contract violation."""


def canonical_json_bytes(payload: Any) -> bytes:
    return (
        json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
        + "\n"
    ).encode("utf-8")


def canonical_sha256(payload: Any) -> str:
    return hashlib.sha256(canonical_json_bytes(payload).rstrip(b"\n")).hexdigest()


def source_input_records_sha256(records: Sequence[Mapping[str, Any]]) -> str:
    """Hash the exact source-input set, including valid-empty occurrence IDs."""

    normalized = [
        {
            field: str(record.get(field, ""))
            for field in SOURCE_INPUT_RECORD_FIELDS
        }
        for record in records
    ]
    normalized.sort(
        key=lambda record: tuple(
            record[field] for field in SOURCE_INPUT_RECORD_FIELDS
        )
    )
    return canonical_sha256(normalized)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
                digest.update(block)
    except OSError as exc:
        raise ControllerError(f"cannot read artifact: {path}") from exc
    return digest.hexdigest()


def require_sha256(value: object, label: str) -> str:
    if not isinstance(value, str) or SHA256_RE.fullmatch(value) is None:
        raise ControllerError(f"{label} must be a lowercase SHA-256")
    return value


def require_mapping(value: object, label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ControllerError(f"{label} must be an object")
    return value


def require_sequence(value: object, label: str) -> list[Any]:
    if not isinstance(value, list):
        raise ControllerError(f"{label} must be an array")
    return value


def require_nonnegative_int(value: object, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ControllerError(f"{label} must be a non-negative integer")
    return value


def require_positive_int(value: object, label: str) -> int:
    value = require_nonnegative_int(value, label)
    if value == 0:
        raise ControllerError(f"{label} must be positive")
    return value


def require_exact_keys(
    payload: Mapping[str, Any], expected: Iterable[str], label: str
) -> None:
    expected_set = set(expected)
    observed = set(payload)
    if observed != expected_set:
        raise ControllerError(
            f"{label} key inventory differs: "
            f"missing={sorted(expected_set - observed)} "
            f"extra={sorted(observed - expected_set)}"
        )


def load_pinned_json(
    path: Path, expected_sha256: str, label: str
) -> tuple[dict[str, Any], dict[str, Any]]:
    expected = require_sha256(expected_sha256, f"{label} expected SHA-256")
    absolute = path.absolute()
    if not absolute.is_file() or absolute.is_symlink():
        raise ControllerError(f"{label} must be a regular non-symlink file: {absolute}")
    observed = file_sha256(absolute)
    if observed != expected:
        raise ControllerError(
            f"{label} SHA-256 differs: expected={expected} observed={observed}"
        )
    try:
        payload = json.loads(absolute.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ControllerError(f"{label} is not readable JSON") from exc
    return require_mapping(payload, label), {
        "path": str(absolute),
        "sha256": observed,
        "size_bytes": absolute.stat().st_size,
    }


def safe_remote_namespace(value: object, label: str) -> str:
    if (
        not isinstance(value, str)
        or not value.startswith("/")
        or value == "/"
        or value != value.strip()
        or "\x00" in value
        or ".." in PurePosixPath(value).parts
    ):
        raise ControllerError(f"{label} must be a safe absolute namespace")
    return value.rstrip("/")


def path_below(path: str, root: str) -> bool:
    child = PurePosixPath(path)
    parent = PurePosixPath(root)
    return child != parent and parent in child.parents


def safe_fresh_local_output(path: Path, label: str) -> Path:
    if not path.is_absolute() or path == Path("/") or path.name in {"", ".", ".."}:
        raise ControllerError(f"{label} must be a safe absolute path")
    if ".." in path.parts or os.path.lexists(path):
        raise ControllerError(f"{label} already exists or is unsafe: {path}")
    parent = path.parent
    if not parent.is_dir() or parent.is_symlink():
        raise ControllerError(f"{label} parent must be an existing real directory")
    return parent.resolve(strict=True) / path.name


def validate_standard_sdcc_scratch_alias(path: Path, label: str) -> bool:
    """Accept only the exact per-user SDCC scratch alias mapping.

    SDCC exposes ``/sphenix/u/<user>/scratch`` as the stable user-facing alias
    for ``/gpfs/mnt/gpfs02/sphenix/user/<user>``.  Campaign manifests retain
    the stable alias, while this check proves that it still targets the exact
    same user's canonical GPFS root.  All other symlinked parents remain
    forbidden.
    """

    try:
        relative = path.relative_to(SDCC_USER_ALIAS_ROOT)
    except ValueError:
        return False
    parts = relative.parts
    if (
        len(parts) < 3
        or parts[1] != "scratch"
        or SAFE_SDCC_USER_RE.fullmatch(parts[0]) is None
    ):
        raise ControllerError(f"{label} is not a valid SDCC scratch alias path")
    user = parts[0]
    descendant = parts[2:]
    alias_root = SDCC_USER_ALIAS_ROOT / user / "scratch"
    canonical_root = SDCC_CANONICAL_USER_ROOT / user
    expected = canonical_root.joinpath(*descendant)
    try:
        if not alias_root.is_symlink():
            raise ControllerError(
                f"{label} SDCC scratch alias root is not a symlink"
            )
        alias_target = alias_root.resolve(strict=True)
        canonical_target = canonical_root.resolve(strict=True)
    except OSError as exc:
        raise ControllerError(
            f"{label} SDCC scratch alias root cannot be resolved"
        ) from exc
    if alias_target != canonical_root or canonical_target != canonical_root:
        raise ControllerError(
            f"{label} SDCC scratch alias does not target the exact user GPFS root"
        )
    first_resolution = path.resolve(strict=False)
    second_resolution = path.resolve(strict=False)
    if first_resolution != expected or second_resolution != expected:
        raise ControllerError(
            f"{label} SDCC scratch alias resolution is unstable or escapes GPFS"
        )
    ancestor = expected.parent
    while not os.path.lexists(ancestor):
        if ancestor == canonical_root:
            break
        ancestor = ancestor.parent
    if canonical_root not in (ancestor, *ancestor.parents):
        raise ControllerError(f"{label} SDCC scratch target escapes its user root")
    cursor = ancestor
    while cursor != canonical_root:
        if cursor.is_symlink():
            raise ControllerError(
                f"{label} canonical GPFS parent chain contains a symlink"
            )
        cursor = cursor.parent
    return True


def materializer_revalidation_root(path: Path) -> Path:
    """Return a canonical GPFS path only after exact scratch-alias proof."""

    if validate_standard_sdcc_scratch_alias(
        path,
        "materializer execution revalidation root",
    ):
        return path.resolve(strict=False)
    return path


def safe_fresh_local_tree_output(path: Path, label: str) -> Path:
    """Validate a fresh absolute output whose parent chain may not exist yet."""

    if not path.is_absolute() or path == Path("/") or path.name in {"", ".", ".."}:
        raise ControllerError(f"{label} must be a safe absolute path")
    if ".." in path.parts or os.path.lexists(path):
        raise ControllerError(f"{label} already exists or is unsafe: {path}")
    if validate_standard_sdcc_scratch_alias(path, label):
        return path.absolute()
    ancestor = path.parent
    while not os.path.lexists(ancestor):
        if ancestor == ancestor.parent:
            raise ControllerError(f"{label} has no existing parent directory")
        ancestor = ancestor.parent
    if not ancestor.is_dir() or ancestor.is_symlink():
        raise ControllerError(
            f"{label} nearest existing parent must be a real directory"
        )
    cursor = ancestor
    while cursor != cursor.parent:
        if cursor.is_symlink():
            raise ControllerError(f"{label} parent chain contains a symlink")
        cursor = cursor.parent
    return path.absolute()


def write_new_bytes(path: Path, data: bytes) -> None:
    destination = safe_fresh_local_output(path, "output")
    try:
        descriptor = os.open(
            destination,
            os.O_WRONLY | os.O_CREAT | os.O_EXCL,
            SHARED_FILE_MODE,
        )
        os.fchmod(descriptor, SHARED_FILE_MODE)
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
    except OSError as exc:
        raise ControllerError(f"cannot create output: {destination}") from exc
    if destination.read_bytes() != data:
        raise ControllerError(f"output readback differs: {destination}")


def write_new_json(path: Path, payload: Mapping[str, Any]) -> None:
    write_new_bytes(path, canonical_json_bytes(payload))


def write_replace_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(mode=SHARED_DIRECTORY_MODE, parents=False, exist_ok=True)
    path.parent.chmod(SHARED_DIRECTORY_MODE)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    if os.path.lexists(temporary):
        raise ControllerError(f"temporary receipt path already exists: {temporary}")
    data = canonical_json_bytes(payload)
    try:
        descriptor = os.open(
            temporary,
            os.O_WRONLY | os.O_CREAT | os.O_EXCL,
            SHARED_FILE_MODE,
        )
        os.fchmod(descriptor, SHARED_FILE_MODE)
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    except OSError as exc:
        raise ControllerError(f"cannot write durable receipt: {path}") from exc
    if path.read_bytes() != data:
        raise ControllerError(f"receipt readback differs: {path}")


def materializer_context(args: argparse.Namespace) -> dict[str, Any]:
    validation_root = args.validation_root.absolute()
    namespace = argparse.Namespace(
        plan=args.plan,
        plan_sha256=args.plan_sha256,
        preflight_receipt=args.preflight_receipt,
        preflight_receipt_sha256=args.preflight_receipt_sha256,
        storage_certificate=args.storage_certificate,
        storage_certificate_sha256=args.storage_certificate_sha256,
        controller_budget=args.controller_budget,
        controller_budget_sha256=args.controller_budget_sha256,
        artifact_profile=args.artifact_profile,
        staging_root=validation_root,
        overwrite=False,
    )
    try:
        return materializer.validate_plan_and_evidence(namespace)
    except materializer.ControllerError as exc:
        raise ControllerError(
            f"upstream materializer validation failed: {exc}"
        ) from exc


def validate_dry_stage(
    context: Mapping[str, Any],
    dry_staging_root: Path,
    expected_manifest_sha256: str,
) -> dict[str, Any]:
    root = dry_staging_root.absolute()
    if not root.is_dir() or root.is_symlink():
        raise ControllerError("dry staging root must be a real directory")
    observed_names = sorted(path.name for path in root.iterdir())
    expected_names = sorted(materializer.OUTPUT_FILENAMES)
    if observed_names != expected_names:
        raise ControllerError(
            "dry staging inventory differs: "
            f"expected={expected_names} observed={observed_names}"
        )
    expected_artifacts = materializer.build_staged_artifacts(context)
    inventory: dict[str, dict[str, Any]] = {}
    for name in materializer.OUTPUT_FILENAMES:
        path = root / name
        if not path.is_file() or path.is_symlink():
            raise ControllerError(f"dry staged artifact is not a regular file: {path}")
        observed = path.read_bytes()
        expected = expected_artifacts[name]
        if observed != expected:
            raise ControllerError(f"dry staged artifact bytes differ: {name}")
        inventory[name] = {
            "path": str(path),
            "sha256": hashlib.sha256(observed).hexdigest(),
            "size_bytes": len(observed),
        }
    manifest = inventory["materialization_manifest.json"]
    if manifest["sha256"] != require_sha256(
        expected_manifest_sha256, "dry materialization manifest SHA-256"
    ):
        raise ControllerError("dry materialization manifest SHA-256 differs")
    return {
        "root": str(root),
        "manifest": manifest,
        "artifacts": inventory,
    }


def validate_authorization(
    payload: Mapping[str, Any],
    *,
    artifact: Mapping[str, Any],
    context: Mapping[str, Any],
    dry_stage: Mapping[str, Any],
    now_seconds: int | None = None,
) -> dict[str, Any]:
    require_exact_keys(
        payload,
        {
            "schema",
            "status",
            "campaign",
            "bindings",
            "limits",
            "duplicate_guard",
            "capacity_guard",
            "quota_guard",
            "approval",
            "provenance",
            "expires_at_unix",
        },
        "submission authorization",
    )
    if (
        payload.get("schema") != AUTHORIZATION_SCHEMA
        or payload.get("status") != "PASS_EXACT_SUBMISSION_AUTHORIZED"
    ):
        raise ControllerError("submission authorization schema/status differs")

    campaign = require_mapping(payload.get("campaign"), "authorization campaign")
    if campaign != context["campaign"]:
        raise ControllerError("authorization campaign differs")

    bindings = require_mapping(payload.get("bindings"), "authorization bindings")
    require_exact_keys(
        bindings,
        {
            "plan_sha256",
            "dry_materialization_manifest_sha256",
            "storage_certificate_sha256",
            "duplicate_fingerprint_sha256",
            "execution_fingerprint_sha256",
        },
        "authorization bindings",
    )
    expected_bindings = {
        "plan_sha256": context["plan_artifact"]["sha256"],
        "dry_materialization_manifest_sha256": dry_stage["manifest"]["sha256"],
        "storage_certificate_sha256": context["storage_certificate_artifact"]["sha256"],
        "duplicate_fingerprint_sha256": context["duplicate_fingerprint_sha256"],
        "execution_fingerprint_sha256": context["execution_fingerprint_sha256"],
    }
    if bindings != expected_bindings:
        raise ControllerError("authorization bindings differ")

    limits = require_mapping(payload.get("limits"), "authorization limits")
    if limits != {
        "row_count": EXPECTED_ROW_COUNT,
        "job_count": EXPECTED_JOB_COUNT,
        "group_size": EXPECTED_GROUP_SIZE,
        "request_memory_mb": EXPECTED_REQUEST_MEMORY_MB,
        "output_pair_count": EXPECTED_OUTPUT_PAIRS,
    }:
        raise ControllerError("authorization execution limits differ")

    duplicate = require_mapping(
        payload.get("duplicate_guard"), "authorization duplicate guard"
    )
    if duplicate != {
        "status": "PASS_NO_ACTIVE_DUPLICATE",
        "active_matching_jobs": 0,
        "fingerprint_sha256": context["duplicate_fingerprint_sha256"],
    }:
        raise ControllerError("authorization duplicate guard differs")

    capacity = require_mapping(payload.get("capacity_guard"), "capacity guard")
    if capacity != {
        "status": "PASS",
        "group_size": EXPECTED_GROUP_SIZE,
        "request_memory_mb": EXPECTED_REQUEST_MEMORY_MB,
    }:
        raise ControllerError("authorization capacity guard differs")

    quota = require_mapping(payload.get("quota_guard"), "quota guard")
    if quota.get("status") != "PASS" or quota.get("authoritative") is not True:
        raise ControllerError("authoritative quota guard is not PASS")
    require_sha256(quota.get("certificate_sha256"), "quota certificate SHA-256")
    if quota["certificate_sha256"] != context["storage_certificate_artifact"]["sha256"]:
        raise ControllerError("quota guard storage certificate differs")

    approval = require_mapping(payload.get("approval"), "approval")
    if (
        approval.get("explicit") is not True
        or approval.get("scope")
        != "THE134_FULL_SOURCE_COMPLETE_GROUP7_EXTRACTION"
        or not isinstance(approval.get("principal"), str)
        or not approval["principal"].strip()
    ):
        raise ControllerError("explicit exact-scope approval is missing")
    require_sha256(approval.get("evidence_sha256"), "approval evidence SHA-256")

    provenance = require_mapping(payload.get("provenance"), "provenance")
    require_exact_keys(
        provenance, {"codex_chat_name", "codex_thread_id"}, "provenance"
    )
    for field in ("codex_chat_name", "codex_thread_id"):
        value = provenance.get(field)
        if (
            not isinstance(value, str)
            or value in {"", "unknown"}
            or SAFE_LABEL_RE.fullmatch(value) is None
        ):
            raise ControllerError(f"authorization provenance {field} is unsafe")

    expires_at = require_positive_int(
        payload.get("expires_at_unix"), "authorization expiry"
    )
    now = int(time.time()) if now_seconds is None else now_seconds
    if expires_at <= now:
        raise ControllerError("submission authorization is expired")
    return {
        "artifact": dict(artifact),
        "expires_at_unix": expires_at,
        "provenance": dict(provenance),
        "approval": dict(approval),
    }


def public_row(row: Mapping[str, Any], global_job_start: int) -> dict[str, Any]:
    execution = require_mapping(row.get("execution_contract"), "row execution")
    input_contract = require_mapping(row.get("input_contract"), "row input")
    bundle = require_mapping(row.get("bundle_contract"), "row bundle")
    materialization_environment = require_mapping(
        execution.get("materialization_environment"), "materialization environment"
    )
    worker_environment = require_mapping(
        execution.get("worker_environment"), "worker environment"
    )
    expected_jobs = require_positive_int(
        input_contract.get("expected_job_count"), "row expected job count"
    )
    return {
        "row_id": row["row_id"],
        "system": row["system"],
        "lane": row["lane"],
        "dataset": row["dataset"],
        "sample": row["sample"],
        "source_role": row["source_role"],
        "source_period": row["source_period"],
        "source_si_di_role": row["source_si_di_role"],
        "global_job_start": global_job_start,
        "expected_job_count": expected_jobs,
        "row_fingerprint_sha256": row["row_fingerprint_sha256"],
        "row_partition_sha256": input_contract["row_partition_sha256"],
        "full_source_manifest_sha256": input_contract[
            "full_source_manifest_sha256"
        ],
        "config_sha256": bundle["config"]["sha256"],
        "code_sha256": bundle["code_sha256"],
        "submitter": dict(bundle["submitter"]),
        "submitter_argv": list(execution["existing_submitter_argv"]),
        "materialization_environment": dict(
            sorted(materialization_environment.items())
        ),
        "materialization_environment_sha256": materializer.canonical_sha256(
            materialization_environment
        ),
        "worker_environment_sha256": materializer.canonical_sha256(
            worker_environment
        ),
        "analysis_output_namespace": execution["analysis_output_namespace"],
        "training_sidecar_template": execution["multiview_sidecar_template"],
        "submit_namespace": execution["submit_namespace"],
        "evidence_namespace": execution["evidence_namespace"],
    }


def build_execution_manifest(
    context: Mapping[str, Any],
    dry_stage: Mapping[str, Any],
    authorization: Mapping[str, Any],
) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    global_job_start = 0
    for row in context["rows"]:
        record = public_row(row, global_job_start)
        rows.append(record)
        global_job_start += record["expected_job_count"]
    if global_job_start != EXPECTED_JOB_COUNT:
        raise ControllerError("execution row job total differs")
    return {
        "schema": EXECUTION_SCHEMA,
        "status": READY_STATUS,
        "artifact_profile": ARTIFACT_PROFILE,
        "campaign": dict(context["campaign"]),
        "counts": {
            "row_count": EXPECTED_ROW_COUNT,
            "job_count": EXPECTED_JOB_COUNT,
            "group_size": EXPECTED_GROUP_SIZE,
            "request_memory_mb": EXPECTED_REQUEST_MEMORY_MB,
            "output_pair_count": EXPECTED_OUTPUT_PAIRS,
        },
        "bindings": {
            "plan": dict(context["plan_artifact"]),
            "resolver_preflight_receipt": dict(
                context["preflight_receipt_artifact"]
            ),
            "storage_certificate": dict(context["storage_certificate_artifact"]),
            "controller_budget": dict(context["controller_budget_artifact"]),
            "immutable_bundle": dict(context["bundle_artifact"]),
            "immutable_materialization": dict(
                context["materialization_artifact"]
            ),
            "source_authority": dict(context["source_artifact"]),
            "source_partition": dict(context["partition_artifact"]),
            "dry_materialization": dict(dry_stage["manifest"]),
            "authorization": dict(authorization["artifact"]),
            "duplicate_fingerprint_sha256": context[
                "duplicate_fingerprint_sha256"
            ],
            "execution_fingerprint_sha256": context[
                "execution_fingerprint_sha256"
            ],
        },
        "authorization": {
            "expires_at_unix": authorization["expires_at_unix"],
            "provenance": dict(authorization["provenance"]),
            "approval": dict(authorization["approval"]),
        },
        "authority": {
            "submission_performed": False,
            "submission_authority": True,
            "full_training_authority": 0,
            "full_extraction_authority": False,
            "science_freeze_authority": False,
            "broad_production_authority": False,
            "physics_output_authority": False,
            "canonical_promotion": False,
            "the121_authority": False,
            "the122_authority": False,
        },
        "rows": rows,
        "boundaries": [
            "Only the exact pinned existing submitter may be invoked.",
            "No direct Condor query or job control is implemented.",
            "Any partial submission is a hard stop with preserved receipts.",
            "No automatic retry, release, remove, resubmit, merge, or cleanup.",
            (
                "No model, WP, science-freeze, production, physics-output, "
                "or CANONICAL authority."
            ),
        ],
    }


def load_bound_partition_chunks(
    bindings: Mapping[str, Any],
    rows: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    binding = require_mapping(
        bindings.get("source_partition"), "execution source-partition binding"
    )
    path = Path(str(binding.get("path", "")))
    if not path.is_absolute() or not path.is_file() or path.is_symlink():
        raise ControllerError("execution source partition is missing or unsafe")
    expected_sha256 = require_sha256(
        binding.get("sha256"), "execution source-partition SHA-256"
    )
    if file_sha256(path) != expected_sha256:
        raise ControllerError("execution source-partition SHA-256 changed")
    expected_order: list[tuple[str, int]] = []
    for row in rows:
        expected_order.extend(
            (str(row["row_id"]), index)
            for index in range(int(row["expected_job_count"]))
        )
    chunks: list[dict[str, Any]] = []
    try:
        stream = path.open("r", encoding="utf-8")
    except OSError as exc:
        raise ControllerError("cannot read execution source partition") from exc
    with stream:
        for global_index, line in enumerate(stream):
            if global_index >= len(expected_order) or not line.strip():
                raise ControllerError("execution source-partition length differs")
            try:
                chunk = require_mapping(
                    json.loads(line), "execution source-partition chunk"
                )
            except json.JSONDecodeError as exc:
                raise ControllerError(
                    "execution source partition contains invalid JSON"
                ) from exc
            require_exact_keys(
                chunk,
                set(materializer.resolver.EXECUTION_CHUNK_KEYS),
                "execution source-partition chunk",
            )
            expected_row_id, expected_row_index = expected_order[global_index]
            tuple_hashes = require_sequence(
                chunk.get("tuple_input_sha256s"),
                "execution chunk tuple-input SHA-256s",
            )
            tuple_count = require_positive_int(
                chunk.get("tuple_count"), "execution chunk tuple count"
            )
            if (
                chunk.get("schema") != materializer.resolver.CHUNK_SCHEMA
                or chunk.get("global_chunk_index") != global_index
                or chunk.get("row_id") != expected_row_id
                or chunk.get("chunk_index") != expected_row_index
                or chunk.get("group_size") != EXPECTED_GROUP_SIZE
                or tuple_count != len(tuple_hashes)
                or tuple_count > EXPECTED_GROUP_SIZE
            ):
                raise ControllerError(
                    "execution source-partition membership differs"
                )
            for index, value in enumerate(tuple_hashes):
                require_sha256(
                    value,
                    f"execution chunk tuple-input SHA-256 {index}",
                )
            execution_sha256 = require_sha256(
                chunk.get("execution_chunk_sha256"),
                "execution chunk SHA-256",
            )
            if execution_sha256 != materializer.canonical_sha256(
                {
                    key: value
                    for key, value in chunk.items()
                    if key != "execution_chunk_sha256"
                }
            ):
                raise ControllerError("execution chunk SHA-256 differs")
            chunks.append(dict(chunk))
    if len(chunks) != EXPECTED_JOB_COUNT:
        raise ControllerError("execution source-partition job count differs")
    return chunks


def validate_execution_manifest(
    payload: Mapping[str, Any],
    *,
    load_partition_chunks: bool = True,
) -> dict[str, Any]:
    require_exact_keys(
        payload,
        {
            "schema",
            "status",
            "artifact_profile",
            "campaign",
            "counts",
            "bindings",
            "authorization",
            "authority",
            "rows",
            "boundaries",
        },
        "execution manifest",
    )
    if (
        payload.get("schema") != EXECUTION_SCHEMA
        or payload.get("status") != READY_STATUS
        or payload.get("artifact_profile") != ARTIFACT_PROFILE
    ):
        raise ControllerError("execution manifest schema/status/profile differs")
    counts = require_mapping(payload.get("counts"), "execution counts")
    if counts != {
        "row_count": EXPECTED_ROW_COUNT,
        "job_count": EXPECTED_JOB_COUNT,
        "group_size": EXPECTED_GROUP_SIZE,
        "request_memory_mb": EXPECTED_REQUEST_MEMORY_MB,
        "output_pair_count": EXPECTED_OUTPUT_PAIRS,
    }:
        raise ControllerError("execution manifest counts differ")
    campaign = require_mapping(payload.get("campaign"), "execution campaign")
    require_exact_keys(
        campaign,
        {"tag", "output_root", "evidence_root", "submit_root"},
        "execution campaign",
    )
    for field in ("output_root", "evidence_root", "submit_root"):
        safe_remote_namespace(campaign.get(field), f"campaign.{field}")
    tag = campaign.get("tag")
    if not isinstance(tag, str) or not tag or tag != tag.strip():
        raise ControllerError("execution campaign tag is unsafe")
    bindings = require_mapping(payload.get("bindings"), "execution bindings")
    require_exact_keys(
        bindings,
        {
            "plan",
            "resolver_preflight_receipt",
            "storage_certificate",
            "controller_budget",
            "immutable_bundle",
            "immutable_materialization",
            "source_authority",
            "source_partition",
            "dry_materialization",
            "authorization",
            "duplicate_fingerprint_sha256",
            "execution_fingerprint_sha256",
        },
        "execution bindings",
    )
    for name in (
        "plan",
        "resolver_preflight_receipt",
        "storage_certificate",
        "controller_budget",
        "immutable_bundle",
        "immutable_materialization",
        "source_authority",
        "source_partition",
        "dry_materialization",
        "authorization",
    ):
        binding = require_mapping(bindings.get(name), f"execution binding {name}")
        require_sha256(binding.get("sha256"), f"execution binding {name} SHA-256")
    require_sha256(
        bindings.get("duplicate_fingerprint_sha256"),
        "execution duplicate fingerprint SHA-256",
    )
    require_sha256(
        bindings.get("execution_fingerprint_sha256"),
        "execution fingerprint SHA-256",
    )

    authority = require_mapping(payload.get("authority"), "execution authority")
    require_exact_keys(
        authority,
        {
            "submission_performed",
            "submission_authority",
            "full_training_authority",
            "full_extraction_authority",
            "science_freeze_authority",
            "broad_production_authority",
            "physics_output_authority",
            "canonical_promotion",
            "the121_authority",
            "the122_authority",
        },
        "execution authority",
    )
    if (
        authority.get("submission_authority") is not True
        or authority.get("submission_performed") is not False
    ):
        raise ControllerError("sealed execution submission authority differs")
    for field in FORBIDDEN_AUTHORITY_TRUE:
        if authority.get(field) not in (False, 0):
            raise ControllerError(f"execution manifest unlawfully grants {field}")
    authorization = require_mapping(payload.get("authorization"), "authorization")
    require_exact_keys(
        authorization,
        {"expires_at_unix", "provenance", "approval"},
        "execution authorization",
    )
    expires_at = require_positive_int(
        authorization.get("expires_at_unix"), "execution authorization expiry"
    )
    provenance = require_mapping(authorization.get("provenance"), "provenance")
    require_exact_keys(
        provenance, {"codex_chat_name", "codex_thread_id"}, "provenance"
    )
    for field in ("codex_chat_name", "codex_thread_id"):
        value = provenance.get(field)
        if (
            not isinstance(value, str)
            or value in {"", "unknown"}
            or SAFE_LABEL_RE.fullmatch(value) is None
        ):
            raise ControllerError(f"execution provenance {field} is unsafe")
    approval = require_mapping(authorization.get("approval"), "execution approval")
    if (
        approval.get("explicit") is not True
        or approval.get("scope")
        != "THE134_FULL_SOURCE_COMPLETE_GROUP7_EXTRACTION"
        or not isinstance(approval.get("principal"), str)
        or not approval["principal"].strip()
    ):
        raise ControllerError("execution approval differs")
    require_sha256(
        approval.get("evidence_sha256"), "execution approval evidence SHA-256"
    )
    rows = require_sequence(payload.get("rows"), "execution rows")
    if len(rows) != EXPECTED_ROW_COUNT:
        raise ControllerError("execution row count differs")
    expected_inventory = materializer.resolver.inventory_rows()
    expected_ids = [row["row_id"] for row in expected_inventory]
    observed_ids: list[str] = []
    total_jobs = 0
    expected_global_start = 0
    for expected_inventory_row, raw in zip(expected_inventory, rows):
        row = require_mapping(raw, "execution row")
        require_exact_keys(
            row,
            {
                "row_id",
                "system",
                "lane",
                "dataset",
                "sample",
                "source_role",
                "source_period",
                "source_si_di_role",
                "global_job_start",
                "expected_job_count",
                "row_fingerprint_sha256",
                "row_partition_sha256",
                "full_source_manifest_sha256",
                "config_sha256",
                "code_sha256",
                "submitter",
                "submitter_argv",
                "materialization_environment",
                "materialization_environment_sha256",
                "worker_environment_sha256",
                "analysis_output_namespace",
                "training_sidecar_template",
                "submit_namespace",
                "evidence_namespace",
            },
            "execution row",
        )
        row_id = str(row.get("row_id", ""))
        observed_ids.append(row_id)
        for field in ("row_id", "system", "lane", "dataset", "sample", "source_role"):
            if row.get(field) != expected_inventory_row[field]:
                raise ControllerError(f"{row_id}.{field} differs from frozen inventory")
        count = require_positive_int(
            row.get("expected_job_count"), f"{row_id}.expected_job_count"
        )
        if row.get("global_job_start") != expected_global_start:
            raise ControllerError(f"{row_id} global job offset differs")
        expected_global_start += count
        total_jobs += count
        argv = require_sequence(row.get("submitter_argv"), f"{row_id}.argv")
        submitter = require_mapping(row.get("submitter"), f"{row_id}.submitter")
        submitter_path = submitter.get("path")
        if (
            not isinstance(submitter_path, str)
            or not Path(submitter_path).is_absolute()
        ):
            raise ControllerError(f"{row_id} submitter path is not absolute")
        if argv != [
            submitter_path,
            row.get("dataset"),
            "condorDoAll",
            "groupSize",
            str(EXPECTED_GROUP_SIZE),
            f"SAMPLE={row.get('sample')}",
        ]:
            raise ControllerError(f"{row_id} pinned submitter argv differs")
        require_sha256(submitter.get("sha256"), f"{row_id}.submitter SHA-256")
        if not path_below(
            str(row.get("analysis_output_namespace", "")),
            str(campaign["output_root"]),
        ):
            raise ControllerError(f"{row_id} output namespace escapes campaign")
        if not path_below(
            str(row.get("submit_namespace", "")),
            str(campaign["submit_root"]),
        ):
            raise ControllerError(f"{row_id} submit namespace escapes campaign")
        if not path_below(
            str(row.get("evidence_namespace", "")),
            str(campaign["evidence_root"]),
        ):
            raise ControllerError(f"{row_id} evidence namespace escapes campaign")
        template = str(row.get("training_sidecar_template", ""))
        expected_template = (
            f"{row['analysis_output_namespace']}/training_views/"
            "$(Cluster).$(Process).root"
        )
        if template != expected_template:
            raise ControllerError(f"{row_id} training-sidecar template differs")
        environment = require_mapping(
            row.get("materialization_environment"), f"{row_id}.environment"
        )
        if (
            environment.get("RJ_DAG_DRYRUN") != "1"
            or environment.get("RJ_REQUEST_MEMORY")
            != f"{EXPECTED_REQUEST_MEMORY_MB}MB"
            or environment.get("RJ_AUTO_MERGE") != "0"
            or environment.get("RJ_CLEAN_OUTPUT_BASE") != "0"
            or environment.get("RJ_CONDOR_SEALED_ENVIRONMENT") != "1"
        ):
            raise ControllerError(f"{row_id} sealed execution environment differs")
        if materializer.canonical_sha256(environment) != row.get(
            "materialization_environment_sha256"
        ):
            raise ControllerError(f"{row_id} environment SHA-256 differs")
        for field in (
            "row_fingerprint_sha256",
            "row_partition_sha256",
            "full_source_manifest_sha256",
            "config_sha256",
            "code_sha256",
            "worker_environment_sha256",
        ):
            require_sha256(row.get(field), f"{row_id}.{field}")
    if observed_ids != expected_ids or len(set(observed_ids)) != EXPECTED_ROW_COUNT:
        raise ControllerError("execution rows are missing, reordered, or duplicated")
    if total_jobs != EXPECTED_JOB_COUNT:
        raise ControllerError("execution total job count differs")
    normalized = dict(payload)
    normalized["_authorization_expires_at_unix"] = expires_at
    normalized["_rows"] = rows
    normalized["_partition_chunks"] = (
        load_bound_partition_chunks(bindings, rows)
        if load_partition_chunks
        else []
    )
    return normalized


def parse_submit_report(output: str, expected_job_count: int) -> tuple[int, int]:
    matches = [
        (int(match.group("count")), int(match.group("cluster")))
        for match in SUBMIT_REPORT_RE.finditer(output)
    ]
    if len(matches) != 1:
        raise ControllerError(
            f"submitter reported {len(matches)} clusters, expected exactly one"
        )
    submitted_jobs, cluster_id = matches[0]
    if submitted_jobs != expected_job_count:
        raise ControllerError(
            "submitter job count differs: "
            f"expected={expected_job_count} observed={submitted_jobs}"
        )
    return cluster_id, submitted_jobs


def sealed_submit_environment(
    row: Mapping[str, Any],
    provenance: Mapping[str, Any],
    execution_artifact: Mapping[str, Any],
    authorization_artifact: Mapping[str, Any],
    ambient: Mapping[str, str] | None = None,
) -> dict[str, str]:
    source = dict(os.environ if ambient is None else ambient)
    environment = {
        key: value
        for key, value in source.items()
        if key in SAFE_SUBMIT_AMBIENT_KEYS
    }
    library_path = environment.get("LD_LIBRARY_PATH")
    if library_path is not None:
        entries = library_path.split(":")
        if (
            not library_path
            or "\x00" in library_path
            or "\n" in library_path
            or any(
                not entry
                or not Path(entry).is_absolute()
                or ".." in Path(entry).parts
                for entry in entries
            )
        ):
            raise ControllerError("ambient LD_LIBRARY_PATH is unsafe")
    sealed = require_mapping(
        row.get("materialization_environment"), "row materialization environment"
    )
    environment.update({str(key): str(value) for key, value in sealed.items()})
    environment["RJ_DAG_DRYRUN"] = "0"
    environment["RJ_CODEX_CHAT_NAME"] = str(provenance["codex_chat_name"])
    environment["RJ_CODEX_THREAD_ID"] = str(provenance["codex_thread_id"])
    environment["RJ_THE134_EXTRACTION_EXECUTION_MANIFEST"] = str(
        execution_artifact["path"]
    )
    environment["RJ_THE134_EXTRACTION_EXECUTION_MANIFEST_SHA256"] = str(
        execution_artifact["sha256"]
    )
    environment["RJ_THE134_EXTRACTION_AUTHORIZATION_RECEIPT"] = str(
        authorization_artifact["path"]
    )
    environment["RJ_THE134_EXTRACTION_AUTHORIZATION_RECEIPT_SHA256"] = str(
        authorization_artifact["sha256"]
    )
    environment["RJ_THE134_EXTRACTION_ROW_ID"] = str(row["row_id"])
    environment["RJ_THE134_EXTRACTION_ROW_FINGERPRINT_SHA256"] = str(
        row["row_fingerprint_sha256"]
    )
    environment["RJ_THE134_EXTRACTION_SUBMITTER_PATH"] = str(
        row["submitter"]["path"]
    )
    return environment


def revalidate_bound_artifacts(
    execution: Mapping[str, Any],
    *,
    now_seconds: int,
) -> None:
    """Reopen every sealed prerequisite before the first submitter call."""

    bindings = require_mapping(execution.get("bindings"), "execution bindings")
    artifact_names = (
        "plan",
        "resolver_preflight_receipt",
        "storage_certificate",
        "controller_budget",
        "immutable_bundle",
        "immutable_materialization",
        "source_authority",
        "source_partition",
        "dry_materialization",
        "authorization",
    )
    for name in artifact_names:
        record = require_mapping(bindings.get(name), f"execution binding {name}")
        raw_path = record.get("path")
        if not isinstance(raw_path, str) or not Path(raw_path).is_absolute():
            raise ControllerError(f"execution binding {name} path is not absolute")
        path = Path(raw_path)
        if not path.is_file() or path.is_symlink():
            raise ControllerError(
                f"execution binding {name} is missing or symlinked: {path}"
            )
        expected_sha256 = require_sha256(
            record.get("sha256"), f"execution binding {name} SHA-256"
        )
        if file_sha256(path) != expected_sha256:
            raise ControllerError(f"execution binding {name} SHA-256 changed")
        if "size_bytes" in record:
            expected_size = require_nonnegative_int(
                record.get("size_bytes"),
                f"execution binding {name} size_bytes",
            )
            if path.stat().st_size != expected_size:
                raise ControllerError(f"execution binding {name} size changed")

    authorization_binding = require_mapping(
        bindings["authorization"], "execution authorization binding"
    )
    authorization_payload, observed_authorization_artifact = load_pinned_json(
        Path(str(authorization_binding["path"])),
        str(authorization_binding["sha256"]),
        "bound submission authorization",
    )
    for field in ("path", "sha256", "size_bytes"):
        if authorization_binding.get(field) != observed_authorization_artifact[field]:
            raise ControllerError(
                f"bound submission authorization {field} differs"
            )
    context = {
        "campaign": execution["campaign"],
        "plan_artifact": bindings["plan"],
        "storage_certificate_artifact": bindings["storage_certificate"],
        "duplicate_fingerprint_sha256": bindings[
            "duplicate_fingerprint_sha256"
        ],
        "execution_fingerprint_sha256": bindings[
            "execution_fingerprint_sha256"
        ],
    }
    validated = validate_authorization(
        authorization_payload,
        artifact=observed_authorization_artifact,
        context=context,
        dry_stage={"manifest": bindings["dry_materialization"]},
        now_seconds=now_seconds,
    )
    if (
        validated["expires_at_unix"]
        != execution["authorization"]["expires_at_unix"]
        or validated["provenance"] != execution["authorization"]["provenance"]
        or validated["approval"] != execution["authorization"]["approval"]
    ):
        raise ControllerError("bound submission authorization content differs")

    dry_binding = require_mapping(
        bindings["dry_materialization"],
        "execution dry-materialization binding",
    )
    dry_manifest, observed_dry_artifact = load_pinned_json(
        Path(str(dry_binding["path"])),
        str(dry_binding["sha256"]),
        "bound dry-materialization manifest",
    )
    for field in ("path", "sha256", "size_bytes"):
        if dry_binding.get(field) != observed_dry_artifact.get(field):
            raise ControllerError(
                f"bound dry-materialization manifest {field} differs"
            )
    counts = require_mapping(
        dry_manifest.get("counts"), "dry-materialization counts"
    )
    authority = require_mapping(
        dry_manifest.get("authority"), "dry-materialization authority"
    )
    if (
        dry_manifest.get("schema") != materializer.MANIFEST_SCHEMA
        or dry_manifest.get("status")
        != "PASS_DRY_MATERIALIZED_LOCAL_ONLY"
        or dry_manifest.get("artifact_profile") != ARTIFACT_PROFILE
        or counts.get("logical_source_row_count") != EXPECTED_ROW_COUNT
        or counts.get("job_count") != EXPECTED_JOB_COUNT
        or counts.get("group_size") != EXPECTED_GROUP_SIZE
        or counts.get("request_memory_mb") != EXPECTED_REQUEST_MEMORY_MB
        or authority.get("submission_performed") is not False
        or authority.get("submission_authority") is not False
    ):
        raise ControllerError(
            "bound dry-materialization compact authority differs"
        )


def validate_safe_resume_certificate_artifact(
    certificate_path: Path,
    expected_sha256: str,
    execution: Mapping[str, Any],
    execution_artifact: Mapping[str, Any],
    *,
    now_seconds: int | None = None,
) -> dict[str, Any]:
    """Verify site admission and bind it to this exact sealed execution."""

    payload, artifact = load_pinned_json(
        certificate_path,
        expected_sha256,
        "SDCC safe-resume certificate",
    )
    submit_host = socket.gethostname().split(".")[0]
    try:
        verified = safe_resume_gate.verify_certificate(
            certificate_path.absolute(),
            expected_sha256,
            operation_class="full_extraction",
            campaign_tag=str(execution["campaign"]["tag"]),
            submit_host=submit_host,
            now_seconds=now_seconds,
        )
    except safe_resume_gate.SafeResumeError as exc:
        raise ControllerError(
            f"SDCC safe-resume certificate rejected: {exc}"
        ) from exc
    if verified != payload:
        raise ControllerError("SDCC safe-resume certificate readback differs")
    evidence = require_mapping(
        verified.get("evidence"), "SDCC safe-resume evidence"
    )
    execution_binding = require_mapping(
        evidence.get("execution_binding"),
        "SDCC safe-resume execution binding",
    )
    for field in ("path", "sha256", "size_bytes"):
        if execution_binding.get(field) != execution_artifact.get(field):
            raise ControllerError(
                f"SDCC safe-resume execution binding {field} differs"
            )
    operation = require_mapping(
        verified.get("operation"), "SDCC safe-resume operation"
    )
    if {
        "row_count": operation.get("row_count"),
        "logical_job_count": operation.get("logical_job_count"),
        "group_size": operation.get("group_size"),
        "request_memory_mb": operation.get("request_memory_mb"),
        "max_materialize_per_cluster": operation.get(
            "max_materialize_per_cluster"
        ),
        "max_idle_per_cluster": operation.get("max_idle_per_cluster"),
    } != {
        "row_count": EXPECTED_ROW_COUNT,
        "logical_job_count": EXPECTED_JOB_COUNT,
        "group_size": EXPECTED_GROUP_SIZE,
        "request_memory_mb": EXPECTED_REQUEST_MEMORY_MB,
        "max_materialize_per_cluster": 20,
        "max_idle_per_cluster": 5,
    }:
        raise ControllerError("SDCC safe-resume execution limits differ")
    return artifact


def submission_receipt_base(
    execution_artifact: Mapping[str, Any],
    execution: Mapping[str, Any],
    safe_resume_certificate: Mapping[str, Any],
) -> dict[str, Any]:
    return {
        "schema": SUBMISSION_SCHEMA,
        "status": "SUBMISSION_IN_PROGRESS",
        "artifact_profile": ARTIFACT_PROFILE,
        "campaign": dict(execution["campaign"]),
        "execution_manifest": dict(execution_artifact),
        "safe_resume_certificate": dict(safe_resume_certificate),
        "counts": {
            "expected_row_count": EXPECTED_ROW_COUNT,
            "expected_job_count": EXPECTED_JOB_COUNT,
            "submitted_row_count": 0,
            "submitted_job_count": 0,
        },
        "submission_performed": False,
        "attempt_lock": None,
        "full_training_authority": 0,
        "full_extraction_authority": False,
        "science_freeze_authority": False,
        "broad_production_authority": False,
        "canonical_promotion": False,
        "rows": [],
        "first_bad": None,
    }


def acquire_campaign_attempt_lock(
    execution: Mapping[str, Any],
    execution_artifact: Mapping[str, Any],
    safe_resume_certificate: Mapping[str, Any],
    receipt_root: Path,
    *,
    now_seconds: int,
) -> Path:
    evidence_root = Path(str(execution["campaign"]["evidence_root"]))
    lock_path = evidence_root / "the134_full_extraction_submission_attempt.lock"
    if os.path.lexists(lock_path):
        raise ControllerError(
            "campaign submission-attempt lock already exists; preserve it and "
            "do not replay this campaign"
        )
    lock = {
        "schema": "THE134_FULL_EXTRACTION_SUBMISSION_ATTEMPT_LOCK_V1",
        "status": "ATTEMPT_STARTED_NO_AUTOMATIC_RETRY",
        "campaign": dict(execution["campaign"]),
        "execution_manifest": dict(execution_artifact),
        "safe_resume_certificate": dict(safe_resume_certificate),
        "authorization": dict(execution["bindings"]["authorization"]),
        "receipt_root": str(receipt_root),
        "started_at_unix": now_seconds,
        "provenance": dict(execution["authorization"]["provenance"]),
        "submission_performed": False,
        "note": (
            "This O_EXCL lock is intentionally persistent after success or "
            "failure. A bounded correction requires a fresh campaign identity."
        ),
    }
    write_new_json(lock_path, lock)
    return lock_path


def execute_submission(
    execution: Mapping[str, Any],
    execution_artifact: Mapping[str, Any],
    safe_resume_certificate_path: Path,
    safe_resume_certificate_sha256: str,
    receipt_dir: Path,
    *,
    runner: Callable[..., subprocess.CompletedProcess[str]] = subprocess.run,
    now_seconds: int | None = None,
) -> dict[str, Any]:
    now = int(time.time()) if now_seconds is None else now_seconds
    campaign = execution["campaign"]
    evidence_root_requested = Path(str(campaign["evidence_root"]))
    attempt_lock_requested = (
        evidence_root_requested
        / "the134_full_extraction_submission_attempt.lock"
    )
    if os.path.lexists(attempt_lock_requested):
        raise ControllerError(
            "campaign submission-attempt lock already exists; preserve it and "
            "do not replay this campaign"
        )
    sealed_payload, observed_execution_artifact = load_pinned_json(
        Path(str(execution_artifact["path"])),
        str(execution_artifact["sha256"]),
        "sealed execution manifest",
    )
    if observed_execution_artifact != execution_artifact:
        raise ControllerError("sealed execution manifest artifact differs")
    sealed_execution = validate_execution_manifest(
        sealed_payload,
        load_partition_chunks=False,
    )
    observed_execution = {
        key: value for key, value in execution.items() if not key.startswith("_")
    }
    observed_sealed_execution = {
        key: value
        for key, value in sealed_execution.items()
        if not key.startswith("_")
    }
    if canonical_json_bytes(observed_execution) != canonical_json_bytes(
        observed_sealed_execution
    ):
        raise ControllerError("in-memory execution differs from sealed manifest")
    safe_resume_certificate = validate_safe_resume_certificate_artifact(
        safe_resume_certificate_path,
        safe_resume_certificate_sha256,
        execution,
        execution_artifact,
        now_seconds=now,
    )
    evidence_root = safe_fresh_local_tree_output(
        evidence_root_requested, "campaign evidence root"
    )
    receipt_requested = receipt_dir.absolute()
    if (
        receipt_requested.parent != evidence_root
        or receipt_requested.name in {"", ".", ".."}
        or os.path.lexists(receipt_requested)
    ):
        raise ControllerError(
            "submission receipt directory must be one fresh direct child of "
            "the campaign evidence root"
        )
    revalidate_bound_artifacts(execution, now_seconds=now)
    if execution["_authorization_expires_at_unix"] <= now:
        raise ControllerError("execution authorization expired before submission")
    for field in ("output_root", "submit_root", "evidence_root"):
        if os.path.lexists(campaign[field]):
            raise ControllerError(f"campaign {field} already exists")
    try:
        evidence_root.mkdir(
            mode=SHARED_DIRECTORY_MODE,
            parents=True,
            exist_ok=False,
        )
        evidence_root.chmod(SHARED_DIRECTORY_MODE)
    except OSError as exc:
        raise ControllerError(
            "cannot atomically create the fresh campaign evidence root"
        ) from exc
    receipt_root = safe_fresh_local_output(
        receipt_requested, "submission receipt directory"
    )
    attempt_lock_path = acquire_campaign_attempt_lock(
        execution,
        execution_artifact,
        safe_resume_certificate,
        receipt_root,
        now_seconds=now,
    )
    receipt_root.mkdir(mode=SHARED_DIRECTORY_MODE, parents=False, exist_ok=False)
    receipt_root.chmod(SHARED_DIRECTORY_MODE)
    receipt_path = receipt_root / "submission_receipt.json"
    receipt = submission_receipt_base(
        execution_artifact,
        execution,
        safe_resume_certificate,
    )
    receipt["attempt_lock"] = {
        "path": str(attempt_lock_path),
        "sha256": file_sha256(attempt_lock_path),
        "size_bytes": attempt_lock_path.stat().st_size,
    }
    write_replace_json(receipt_path, receipt)

    provenance = execution["authorization"]["provenance"]
    for row in execution["_rows"]:
        row_id = row["row_id"]
        row_output = row["analysis_output_namespace"]
        row_submit = row["submit_namespace"]
        if os.path.lexists(row_output) or os.path.lexists(row_submit):
            receipt["status"] = "PARTIAL_FAILED_HARD_STOP"
            receipt["first_bad"] = {
                "row_id": row_id,
                "failure": "row output or submit namespace already exists",
            }
            write_replace_json(receipt_path, receipt)
            raise ControllerError(
                f"{row_id} output or submit namespace already exists"
            )
        submitter = Path(str(row["submitter"]["path"]))
        if not submitter.is_file() or submitter.is_symlink():
            failure = f"pinned submitter is missing or symlinked: {submitter}"
            receipt["status"] = "PARTIAL_FAILED_HARD_STOP"
            receipt["first_bad"] = {"row_id": row_id, "failure": failure}
            write_replace_json(receipt_path, receipt)
            raise ControllerError(failure)
        if file_sha256(submitter) != row["submitter"]["sha256"]:
            failure = "pinned submitter SHA-256 changed"
            receipt["status"] = "PARTIAL_FAILED_HARD_STOP"
            receipt["first_bad"] = {"row_id": row_id, "failure": failure}
            write_replace_json(receipt_path, receipt)
            raise ControllerError(f"{row_id} {failure}")
        stdout_path = receipt_root / f"{row_id}.stdout.txt"
        stderr_path = receipt_root / f"{row_id}.stderr.txt"
        started = int(time.time())
        submission_environment = sealed_submit_environment(
            row,
            provenance,
            execution_artifact,
            execution["bindings"]["authorization"],
        )
        inherited_environment = {
            key: submission_environment[key]
            for key in sorted(SAFE_SUBMIT_AMBIENT_KEYS)
            if key in submission_environment
        }
        try:
            completed = runner(
                list(row["submitter_argv"]),
                env=submission_environment,
                text=True,
                capture_output=True,
                check=False,
                timeout=SUBMITTER_TIMEOUT_SECONDS,
            )
        except Exception as exc:
            receipt["status"] = "PARTIAL_FAILED_HARD_STOP"
            receipt["submission_performed"] = True
            receipt["first_bad"] = {
                "row_id": row_id,
                "failure": f"submitter execution error:{type(exc).__name__}:{exc}",
            }
            write_replace_json(receipt_path, receipt)
            raise ControllerError(f"{row_id} submitter execution failed") from exc
        receipt["submission_performed"] = True
        stdout_bytes = completed.stdout.encode("utf-8")
        stderr_bytes = completed.stderr.encode("utf-8")
        write_new_bytes(stdout_path, stdout_bytes)
        write_new_bytes(stderr_path, stderr_bytes)
        if (
            len(stdout_bytes) > SUBMITTER_OUTPUT_LIMIT_BYTES
            or len(stderr_bytes) > SUBMITTER_OUTPUT_LIMIT_BYTES
        ):
            receipt["status"] = "PARTIAL_FAILED_HARD_STOP"
            receipt["first_bad"] = {
                "row_id": row_id,
                "failure": "submitter output exceeded bounded capture limit",
                "stdout_size_bytes": len(stdout_bytes),
                "stderr_size_bytes": len(stderr_bytes),
                "limit_bytes": SUBMITTER_OUTPUT_LIMIT_BYTES,
                "stdout_path": str(stdout_path),
                "stderr_path": str(stderr_path),
            }
            write_replace_json(receipt_path, receipt)
            raise ControllerError(
                f"{row_id} submitter output exceeded bounded capture limit"
            )
        row_receipt: dict[str, Any] = {
            "row_id": row_id,
            "returncode": int(completed.returncode),
            "started_at_unix": started,
            "finished_at_unix": int(time.time()),
            "argv": list(row["submitter_argv"]),
            "row_fingerprint_sha256": row["row_fingerprint_sha256"],
            "inherited_environment": inherited_environment,
            "environment_sha256": canonical_sha256(submission_environment),
            "stdout": {
                "path": str(stdout_path),
                "sha256": file_sha256(stdout_path),
                "size_bytes": stdout_path.stat().st_size,
            },
            "stderr": {
                "path": str(stderr_path),
                "sha256": file_sha256(stderr_path),
                "size_bytes": stderr_path.stat().st_size,
            },
            "expected_job_count": row["expected_job_count"],
            "cluster_id": None,
            "submitted_job_count": 0,
        }
        receipt["rows"].append(row_receipt)
        parse_error: ControllerError | None = None
        try:
            cluster_id, submitted_jobs = parse_submit_report(
                completed.stdout + "\n" + completed.stderr,
                row["expected_job_count"],
            )
        except ControllerError as exc:
            parse_error = exc
        else:
            row_receipt["cluster_id"] = cluster_id
            row_receipt["submitted_job_count"] = submitted_jobs
            receipt["counts"]["submitted_row_count"] += 1
            receipt["counts"]["submitted_job_count"] += submitted_jobs
        if completed.returncode != 0:
            receipt["status"] = "PARTIAL_FAILED_HARD_STOP"
            receipt["first_bad"] = {
                "row_id": row_id,
                "failure": f"submitter exit code {completed.returncode}",
                "cluster_identity_preserved": row_receipt["cluster_id"] is not None,
                "parse_failure": str(parse_error) if parse_error else None,
            }
            write_replace_json(receipt_path, receipt)
            raise ControllerError(
                f"{row_id} submitter exited {completed.returncode}; no retry performed"
            )
        if parse_error is not None:
            receipt["status"] = "PARTIAL_FAILED_HARD_STOP"
            receipt["first_bad"] = {
                "row_id": row_id,
                "failure": str(parse_error),
                "cluster_identity_preserved": False,
            }
            write_replace_json(receipt_path, receipt)
            raise parse_error
        write_replace_json(receipt_path, receipt)

    if (
        receipt["counts"]["submitted_row_count"] != EXPECTED_ROW_COUNT
        or receipt["counts"]["submitted_job_count"] != EXPECTED_JOB_COUNT
    ):
        receipt["status"] = "PARTIAL_FAILED_HARD_STOP"
        receipt["first_bad"] = {"failure": "aggregate submitted counts differ"}
        write_replace_json(receipt_path, receipt)
        raise ControllerError("aggregate submitted counts differ")
    receipt["status"] = SUBMITTED_STATUS
    receipt["first_bad"] = None
    write_replace_json(receipt_path, receipt)
    return receipt


def validate_submission_receipt(
    payload: Mapping[str, Any],
    execution: Mapping[str, Any],
    execution_artifact: Mapping[str, Any],
    *,
    verify_attempt_lock: bool = True,
) -> dict[str, Any]:
    require_exact_keys(
        payload,
        {
            "schema",
            "status",
            "artifact_profile",
            "campaign",
            "execution_manifest",
            "safe_resume_certificate",
            "counts",
            "submission_performed",
            "attempt_lock",
            "full_training_authority",
            "full_extraction_authority",
            "science_freeze_authority",
            "broad_production_authority",
            "canonical_promotion",
            "rows",
            "first_bad",
        },
        "submission receipt",
    )
    if (
        payload.get("schema") != SUBMISSION_SCHEMA
        or payload.get("status") != SUBMITTED_STATUS
        or payload.get("artifact_profile") != ARTIFACT_PROFILE
        or payload.get("submission_performed") is not True
        or payload.get("full_training_authority") != 0
        or payload.get("full_extraction_authority") is not False
        or payload.get("science_freeze_authority") is not False
        or payload.get("broad_production_authority") is not False
        or payload.get("canonical_promotion") is not False
        or payload.get("first_bad") is not None
    ):
        raise ControllerError("submission receipt schema/status/authority differs")
    if payload.get("campaign") != execution["campaign"]:
        raise ControllerError("submission receipt campaign differs")
    safe_resume_binding = require_mapping(
        payload.get("safe_resume_certificate"),
        "submission safe-resume certificate",
    )
    require_exact_keys(
        safe_resume_binding,
        {"path", "sha256", "size_bytes"},
        "submission safe-resume certificate",
    )
    safe_resume_payload, observed_safe_resume = load_pinned_json(
        Path(str(safe_resume_binding.get("path"))),
        str(safe_resume_binding.get("sha256")),
        "bound SDCC safe-resume certificate",
    )
    if observed_safe_resume != safe_resume_binding:
        raise ControllerError("submission safe-resume artifact differs")
    safe_operation = require_mapping(
        safe_resume_payload.get("operation"),
        "bound SDCC safe-resume operation",
    )
    safe_evidence = require_mapping(
        safe_resume_payload.get("evidence"),
        "bound SDCC safe-resume evidence",
    )
    safe_execution = require_mapping(
        safe_evidence.get("execution_binding"),
        "bound SDCC safe-resume execution binding",
    )
    if (
        safe_resume_payload.get("schema")
        != safe_resume_gate.CERTIFICATE_SCHEMA
        or safe_resume_payload.get("status") != "SITE_ADMISSION_PASS"
        or safe_resume_payload.get("submission_performed") is not False
        or safe_operation.get("class") != "full_extraction"
        or safe_operation.get("campaign_tag")
        != execution["campaign"]["tag"]
        or any(
            safe_execution.get(field) != execution_artifact.get(field)
            for field in ("path", "sha256", "size_bytes")
        )
    ):
        raise ControllerError(
            "submission safe-resume authority or execution binding differs"
        )
    attempt_lock = require_mapping(
        payload.get("attempt_lock"), "submission attempt lock"
    )
    require_exact_keys(
        attempt_lock,
        {"path", "sha256", "size_bytes"},
        "submission attempt lock",
    )
    expected_lock_path = (
        f"{execution['campaign']['evidence_root']}/"
        "the134_full_extraction_submission_attempt.lock"
    )
    if attempt_lock.get("path") != expected_lock_path:
        raise ControllerError("submission attempt-lock path differs")
    require_sha256(
        attempt_lock.get("sha256"), "submission attempt-lock SHA-256"
    )
    require_positive_int(
        attempt_lock.get("size_bytes"), "submission attempt-lock size"
    )
    if verify_attempt_lock:
        lock_path = Path(str(attempt_lock["path"]))
        if not lock_path.is_file() or lock_path.is_symlink():
            raise ControllerError(
                "submission attempt-lock artifact is missing or symlinked"
            )
        if (
            lock_path.stat().st_size != attempt_lock["size_bytes"]
            or file_sha256(lock_path) != attempt_lock["sha256"]
        ):
            raise ControllerError("submission attempt-lock artifact differs")
        try:
            lock_payload = require_mapping(
                json.loads(lock_path.read_text(encoding="utf-8")),
                "submission attempt-lock payload",
            )
        except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
            raise ControllerError(
                "submission attempt-lock artifact is not readable JSON"
            ) from exc
        require_exact_keys(
            lock_payload,
            {
                "schema",
                "status",
                "campaign",
                "execution_manifest",
                "safe_resume_certificate",
                "authorization",
                "receipt_root",
                "started_at_unix",
                "provenance",
                "submission_performed",
                "note",
            },
            "submission attempt-lock payload",
        )
        receipt_root = safe_remote_namespace(
            lock_payload.get("receipt_root"),
            "submission attempt-lock receipt root",
        )
        if (
            lock_payload.get("schema")
            != "THE134_FULL_EXTRACTION_SUBMISSION_ATTEMPT_LOCK_V1"
            or lock_payload.get("status")
            != "ATTEMPT_STARTED_NO_AUTOMATIC_RETRY"
            or lock_payload.get("campaign") != execution["campaign"]
            or lock_payload.get("execution_manifest") != execution_artifact
            or lock_payload.get("safe_resume_certificate")
            != safe_resume_binding
            or lock_payload.get("authorization")
            != execution["bindings"]["authorization"]
            or lock_payload.get("provenance")
            != execution["authorization"]["provenance"]
            or lock_payload.get("submission_performed") is not False
            or not path_below(
                receipt_root, execution["campaign"]["evidence_root"]
            )
        ):
            raise ControllerError(
                "submission attempt-lock payload identity differs"
            )
        require_positive_int(
            lock_payload.get("started_at_unix"),
            "submission attempt-lock started_at_unix",
        )
        if (
            lock_payload.get("note")
            != "This O_EXCL lock is intentionally persistent after success or "
            "failure. A bounded correction requires a fresh campaign identity."
        ):
            raise ControllerError("submission attempt-lock note differs")
    bound_execution = require_mapping(
        payload.get("execution_manifest"), "submission execution binding"
    )
    if bound_execution != execution_artifact:
        raise ControllerError("submission receipt execution binding differs")
    counts = require_mapping(payload.get("counts"), "submission counts")
    if counts != {
        "expected_row_count": EXPECTED_ROW_COUNT,
        "expected_job_count": EXPECTED_JOB_COUNT,
        "submitted_row_count": EXPECTED_ROW_COUNT,
        "submitted_job_count": EXPECTED_JOB_COUNT,
    }:
        raise ControllerError("submission receipt aggregate counts differ")
    rows = require_sequence(payload.get("rows"), "submission rows")
    expected_rows = execution["_rows"]
    if len(rows) != EXPECTED_ROW_COUNT:
        raise ControllerError("submission receipt row count differs")
    clusters: set[int] = set()
    evidence_paths: set[str] = set()
    total_jobs = 0
    normalized_rows: list[dict[str, Any]] = []
    for expected, raw in zip(expected_rows, rows):
        row = require_mapping(raw, "submission row")
        require_exact_keys(
            row,
            {
                "row_id",
                "returncode",
                "started_at_unix",
                "finished_at_unix",
                "argv",
                "row_fingerprint_sha256",
                "inherited_environment",
                "environment_sha256",
                "stdout",
                "stderr",
                "expected_job_count",
                "cluster_id",
                "submitted_job_count",
            },
            "submission row",
        )
        if (
            row.get("row_id") != expected["row_id"]
            or row.get("returncode") != 0
            or row.get("expected_job_count") != expected["expected_job_count"]
            or row.get("submitted_job_count") != expected["expected_job_count"]
            or row.get("row_fingerprint_sha256")
            != expected["row_fingerprint_sha256"]
        ):
            raise ControllerError(f"{expected['row_id']} submission row differs")
        started = require_positive_int(
            row.get("started_at_unix"), f"{expected['row_id']}.started_at_unix"
        )
        finished = require_positive_int(
            row.get("finished_at_unix"), f"{expected['row_id']}.finished_at_unix"
        )
        if finished < started:
            raise ControllerError(
                f"{expected['row_id']} submission time ordering differs"
            )
        argv = require_sequence(row.get("argv"), f"{expected['row_id']}.argv")
        if argv != expected["submitter_argv"]:
            raise ControllerError(f"{expected['row_id']} submission argv differs")
        inherited_environment = require_mapping(
            row.get("inherited_environment"),
            f"{expected['row_id']}.inherited_environment",
        )
        if (
            not set(inherited_environment).issubset(SAFE_SUBMIT_AMBIENT_KEYS)
            or not all(
                isinstance(key, str) and isinstance(value, str)
                for key, value in inherited_environment.items()
            )
        ):
            raise ControllerError(
                f"{expected['row_id']} inherited environment is unsafe"
            )
        expected_environment_sha256 = canonical_sha256(
            sealed_submit_environment(
                expected,
                execution["authorization"]["provenance"],
                execution_artifact,
                execution["bindings"]["authorization"],
                ambient=inherited_environment,
            )
        )
        if row.get("environment_sha256") != expected_environment_sha256:
            raise ControllerError(
                f"{expected['row_id']} submission environment differs"
            )
        require_sha256(
            row.get("environment_sha256"),
            f"{expected['row_id']}.environment_sha256",
        )
        for stream_name in ("stdout", "stderr"):
            stream = require_mapping(
                row.get(stream_name), f"{expected['row_id']}.{stream_name}"
            )
            require_exact_keys(
                stream,
                {"path", "sha256", "size_bytes"},
                f"{expected['row_id']}.{stream_name}",
            )
            stream_path = safe_remote_namespace(
                stream.get("path"), f"{expected['row_id']}.{stream_name}.path"
            )
            if not path_below(stream_path, execution["campaign"]["evidence_root"]):
                raise ControllerError(
                    f"{expected['row_id']} {stream_name} escapes evidence namespace"
                )
            if stream_path in evidence_paths:
                raise ControllerError("submission receipt duplicates an evidence path")
            evidence_paths.add(stream_path)
            require_sha256(
                stream.get("sha256"),
                f"{expected['row_id']}.{stream_name}.sha256",
            )
            require_nonnegative_int(
                stream.get("size_bytes"),
                f"{expected['row_id']}.{stream_name}.size_bytes",
            )
        cluster = require_positive_int(
            row.get("cluster_id"), f"{expected['row_id']}.cluster_id"
        )
        if cluster in clusters:
            raise ControllerError("submission receipt duplicates a cluster")
        clusters.add(cluster)
        total_jobs += row["submitted_job_count"]
        normalized_rows.append(dict(row))
    if total_jobs != EXPECTED_JOB_COUNT:
        raise ControllerError("submission receipt total job count differs")
    normalized = dict(payload)
    normalized["_rows"] = normalized_rows
    normalized["_row_by_id"] = {row["row_id"]: row for row in normalized_rows}
    return normalized


def validate_health_artifact(
    artifact: object,
    *,
    label: str,
    minimum_bytes: int,
) -> dict[str, Any]:
    record = require_mapping(artifact, label)
    require_exact_keys(
        record,
        {
            "path",
            "size_bytes",
            "sha256",
            "status",
            "readable",
            "zombie",
            "recovered",
        },
        label,
    )
    path = safe_remote_namespace(record.get("path"), f"{label}.path")
    size = require_positive_int(record.get("size_bytes"), f"{label}.size_bytes")
    require_sha256(record.get("sha256"), f"{label}.sha256")
    if (
        size < minimum_bytes
        or record.get("status") != "PASS"
        or record.get("readable") is not True
        or record.get("zombie") is not False
        or record.get("recovered") is not False
    ):
        raise ControllerError(f"{label} artifact health differs")
    return {**record, "path": path}


def validate_staged_chunk_list(
    artifact: object,
    *,
    expected_chunk: Mapping[str, Any],
    row: Mapping[str, Any],
    label: str,
) -> dict[str, Any]:
    record = require_mapping(artifact, label)
    require_exact_keys(
        record,
        {
            "path",
            "size_bytes",
            "sha256",
            "status",
            "readable",
            "execution_chunk_sha256",
            "tuple_input_sha256s",
        },
        label,
    )
    path_text = safe_remote_namespace(record.get("path"), f"{label}.path")
    if not path_below(path_text, row["submit_namespace"]):
        raise ControllerError(f"{label} escapes the row submit namespace")
    path = Path(path_text)
    if not path.is_file() or path.is_symlink():
        raise ControllerError(f"{label} is missing or symlinked")
    size_bytes = require_positive_int(
        record.get("size_bytes"), f"{label}.size_bytes"
    )
    expected_file_sha256 = require_sha256(
        record.get("sha256"), f"{label}.sha256"
    )
    if (
        record.get("status") != "PASS"
        or record.get("readable") is not True
        or path.stat().st_size != size_bytes
        or file_sha256(path) != expected_file_sha256
    ):
        raise ControllerError(f"{label} health/readback differs")
    if (
        record.get("execution_chunk_sha256")
        != expected_chunk["execution_chunk_sha256"]
    ):
        raise ControllerError(f"{label} frozen partition identity differs")
    declared_tuple_hashes = require_sequence(
        record.get("tuple_input_sha256s"), f"{label}.tuple_input_sha256s"
    )
    expected_tuple_hashes = list(expected_chunk["tuple_input_sha256s"])
    if declared_tuple_hashes != expected_tuple_hashes:
        raise ControllerError(f"{label} tuple-input identities differ")

    try:
        data = path.read_bytes()
        text = data.decode("utf-8")
    except (OSError, UnicodeDecodeError) as exc:
        raise ControllerError(f"{label} is not readable UTF-8") from exc
    lines = text.splitlines()
    if len(lines) != expected_chunk["tuple_count"]:
        raise ControllerError(f"{label} tuple count differs")
    list_roles = tuple(materializer.resolver.LIST_ROLES)
    observed_tuple_hashes: list[str] = []
    for line_index, line in enumerate(lines):
        fields = line.split("\t")
        if len(fields) != len(list_roles):
            raise ControllerError(
                f"{label} line {line_index} does not have five inputs"
            )
        inputs: dict[str, str] = {}
        for role, raw in zip(list_roles, fields):
            value = materializer.resolver.executable_line(raw)
            if value is None:
                raise ControllerError(
                    f"{label} line {line_index} contains a non-executable input"
                )
            inputs[role] = value
        observed_tuple_hashes.append(
            materializer.resolver.canonical_sha256(inputs)
        )
    if observed_tuple_hashes != expected_tuple_hashes:
        raise ControllerError(f"{label} bytes differ from frozen tuple inputs")
    return {
        **record,
        "path": path_text,
        "tuple_input_sha256s": declared_tuple_hashes,
    }


def validate_terminal_snapshot(
    payload: Mapping[str, Any],
    *,
    execution: Mapping[str, Any],
    execution_artifact: Mapping[str, Any],
    submission: Mapping[str, Any],
    submission_artifact: Mapping[str, Any],
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    require_exact_keys(
        payload,
        {
            "schema",
            "status",
            "campaign",
            "execution_manifest_sha256",
            "submission_receipt_sha256",
            "captured_at_unix",
            "jobs",
        },
        "terminal snapshot",
    )
    if (
        payload.get("schema") != TERMINAL_SNAPSHOT_SCHEMA
        or payload.get("status") != "COMPLETE"
        or payload.get("campaign") != execution["campaign"]
        or payload.get("execution_manifest_sha256") != execution_artifact["sha256"]
        or payload.get("submission_receipt_sha256") != submission_artifact["sha256"]
    ):
        raise ControllerError("terminal snapshot schema/status/binding differs")
    require_positive_int(payload.get("captured_at_unix"), "snapshot captured_at_unix")
    jobs = require_sequence(payload.get("jobs"), "terminal snapshot jobs")
    if len(jobs) != EXPECTED_JOB_COUNT:
        raise ControllerError("terminal snapshot job count differs")

    expected: list[tuple[dict[str, Any], dict[str, Any], int, int]] = []
    for row in execution["_rows"]:
        submitted = submission["_row_by_id"][row["row_id"]]
        cluster = submitted["cluster_id"]
        for row_job_index in range(row["expected_job_count"]):
            expected.append(
                (
                    row,
                    submitted,
                    row["global_job_start"] + row_job_index,
                    row_job_index,
                )
            )
    seen_scheduler_ids: set[tuple[int, int]] = set()
    seen_paths: set[str] = set()
    seen_occurrences: set[str] = set()
    seen_input_identities: set[tuple[str, str]] = set()
    normalized_jobs: list[dict[str, Any]] = []
    for (row, submitted, global_index, row_index), raw in zip(expected, jobs):
        job = require_mapping(raw, "terminal job")
        require_exact_keys(
            job,
            {
                "row_id",
                "global_job_index",
                "row_job_index",
                "cluster_id",
                "proc_id",
                "job_status",
                "exit_code",
                "num_job_starts",
                "ephemeral_analysis_health",
                "staged_chunk_list",
                "training_sidecar_root",
                "sidecar_health",
            },
            "terminal job",
        )
        row_id = row["row_id"]
        cluster = submitted["cluster_id"]
        if (
            job.get("row_id") != row_id
            or job.get("global_job_index") != global_index
            or job.get("row_job_index") != row_index
            or job.get("cluster_id") != cluster
            or job.get("proc_id") != row_index
            or job.get("job_status") != 4
            or job.get("exit_code") != 0
        ):
            raise ControllerError(f"{row_id}.{row_index} terminal state differs")
        require_positive_int(
            job.get("num_job_starts"), f"{row_id}.{row_index}.num_job_starts"
        )
        scheduler_id = (cluster, row_index)
        if scheduler_id in seen_scheduler_ids:
            raise ControllerError("terminal snapshot duplicates scheduler identity")
        seen_scheduler_ids.add(scheduler_id)
        expected_chunk = execution["_partition_chunks"][global_index]
        staged_chunk = validate_staged_chunk_list(
            job.get("staged_chunk_list"),
            expected_chunk=expected_chunk,
            row=row,
            label=f"{row_id}.{row_index}.staged_chunk_list",
        )

        analysis_health = require_mapping(
            job.get("ephemeral_analysis_health"),
            f"{row_id}.{row_index}.ephemeral_analysis_health",
        )
        require_exact_keys(
            analysis_health,
            {
                "schema",
                "status",
                "mode",
                "analysis_size_bytes",
                "analysis_minimum_bytes",
                "analysis_key_inventory_sha256",
                "analysis_config_present",
                "analysis_directory_present",
                "analysis_histogram_present",
                "analysis_root_non_zombie",
                "analysis_root_non_recovered",
                "analysis_root_retained",
                "sidecar_size_bytes",
                "sidecar_tree_entries",
            },
            "ephemeral analysis health",
        )
        if (
            analysis_health.get("schema")
            != "THE134_EPHEMERAL_ANALYSIS_HEALTH_V1"
            or analysis_health.get("status") != "PASS"
            or analysis_health.get("mode") != "EPHEMERAL_CONDOR_SCRATCH"
            or analysis_health.get("analysis_minimum_bytes") != 50_000
            or analysis_health.get("analysis_config_present") is not True
            or analysis_health.get("analysis_directory_present") is not True
            or analysis_health.get("analysis_histogram_present") is not True
            or analysis_health.get("analysis_root_non_zombie") is not True
            or analysis_health.get("analysis_root_non_recovered") is not True
            or analysis_health.get("analysis_root_retained") is not False
        ):
            raise ControllerError(
                f"{row_id}.{row_index} ephemeral analysis health differs"
            )
        require_positive_int(
            analysis_health.get("analysis_size_bytes"),
            f"{row_id}.{row_index}.analysis_size_bytes",
        )
        if analysis_health["analysis_size_bytes"] < 50_000:
            raise ControllerError(
                f"{row_id}.{row_index} ephemeral analysis ROOT is too small"
            )
        require_sha256(
            analysis_health.get("analysis_key_inventory_sha256"),
            f"{row_id}.{row_index}.analysis_key_inventory_sha256",
        )

        sidecar = validate_health_artifact(
            job.get("training_sidecar_root"),
            label=f"{row_id}.{row_index}.training_sidecar_root",
            minimum_bytes=1,
        )
        expected_sidecar = (
            row["training_sidecar_template"]
            .replace("$(Cluster)", str(cluster))
            .replace("$(Process)", str(row_index))
        )
        if sidecar["path"] != expected_sidecar:
            raise ControllerError(f"{row_id}.{row_index} sidecar path differs")
        for artifact_path in (
            staged_chunk["path"],
            sidecar["path"],
        ):
            if artifact_path in seen_paths:
                raise ControllerError("terminal snapshot duplicates an artifact path")
            seen_paths.add(artifact_path)

        health = require_mapping(
            job.get("sidecar_health"), f"{row_id}.{row_index}.sidecar_health"
        )
        require_exact_keys(
            health,
            {
                "schema",
                "status",
                "profile",
                "certificate_sha256",
                "tree_name",
                "tree_entries",
                "source_occurrence_id_hex",
                "input_uri_sha256",
                "input_file_sha256",
                "source_manifest_sha256",
                "config_sha256",
                "code_sha256",
                "selected_training_rows_by_view",
            },
            "sidecar health",
        )
        if (
            health.get("schema") != "THE134_FULL_EXTRACTION_SIDECAR_HEALTH_V1"
            or health.get("status") != "PASS"
            or health.get("profile") != "photon_training_multiview_v1"
            or health.get("tree_name") != "RJPhotonTrainingViewV1"
            or health.get("source_manifest_sha256")
            != row["full_source_manifest_sha256"]
            or health.get("config_sha256") != row["config_sha256"]
            or health.get("code_sha256") != row["code_sha256"]
        ):
            raise ControllerError(f"{row_id}.{row_index} sidecar identity differs")
        require_sha256(
            health.get("certificate_sha256"), "sidecar health certificate SHA-256"
        )
        input_uri_sha256 = require_sha256(
            health.get("input_uri_sha256"), "sidecar input URI SHA-256"
        )
        input_file_sha256 = require_sha256(
            health.get("input_file_sha256"), "sidecar input file SHA-256"
        )
        input_identity = (input_uri_sha256, input_file_sha256)
        if input_uri_sha256 != input_file_sha256:
            raise ControllerError(
                f"{row_id}.{row_index} staged-chunk input identity differs"
            )
        if input_uri_sha256 != staged_chunk["sha256"]:
            raise ControllerError(
                f"{row_id}.{row_index} sidecar input identity differs from "
                "the frozen staged chunk"
            )
        if input_identity in seen_input_identities:
            raise ControllerError(
                "terminal snapshot duplicates a staged-chunk input identity"
            )
        seen_input_identities.add(input_identity)
        tree_entries = require_nonnegative_int(
            health.get("tree_entries"), "sidecar tree entries"
        )
        if (
            require_positive_int(
                analysis_health.get("sidecar_size_bytes"),
                f"{row_id}.{row_index}.sidecar_size_bytes",
            )
            != sidecar["size_bytes"]
            or require_nonnegative_int(
                analysis_health.get("sidecar_tree_entries"),
                f"{row_id}.{row_index}.sidecar_tree_entries",
            )
            != tree_entries
        ):
            raise ControllerError(
                f"{row_id}.{row_index} worker sidecar health differs"
            )
        occurrence = health.get("source_occurrence_id_hex")
        if (
            not isinstance(occurrence, str)
            or IDENTITY128_RE.fullmatch(occurrence) is None
        ):
            raise ControllerError("sidecar source occurrence identity differs")
        if occurrence in seen_occurrences:
            raise ControllerError("terminal snapshot duplicates source occurrence")
        seen_occurrences.add(occurrence)
        selected = require_mapping(
            health.get("selected_training_rows_by_view"),
            "selected training rows by view",
        )
        if set(selected) != set(SHOWER_VIEWS):
            raise ControllerError("selected training view inventory differs")
        selected_counts: list[int] = []
        for view in SHOWER_VIEWS:
            selected_counts.append(
                require_nonnegative_int(selected[view], f"{view} selected rows")
            )
        if len(set(selected_counts)) != 1:
            raise ControllerError(
                f"{row_id}.{row_index} selected training population differs "
                "across shower views"
            )
        if (
            tree_entries % len(SHOWER_VIEWS) != 0
            or sum(selected_counts) > tree_entries
        ):
            raise ControllerError(
                f"{row_id}.{row_index} seven-view tree population differs"
            )
        normalized_jobs.append(
            {
                **job,
                "staged_chunk_list": staged_chunk,
                "ephemeral_analysis_health": dict(analysis_health),
                "training_sidecar_root": sidecar,
                "sidecar_health": dict(health),
            }
        )
    return dict(payload), normalized_jobs


def terminal_receipt(
    *,
    execution_artifact: Mapping[str, Any],
    submission_artifact: Mapping[str, Any],
    snapshot_artifact: Mapping[str, Any],
    jobs: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    return {
        "schema": TERMINAL_RECEIPT_SCHEMA,
        "status": TERMINAL_STATUS,
        "artifact_profile": ARTIFACT_PROFILE,
        "bindings": {
            "execution_manifest": dict(execution_artifact),
            "submission_receipt": dict(submission_artifact),
            "terminal_snapshot": dict(snapshot_artifact),
        },
        "counts": {
            "row_count": EXPECTED_ROW_COUNT,
            "job_count": len(jobs),
            "ephemeral_analysis_health_count": len(jobs),
            "retained_analysis_root_count": 0,
            "training_sidecar_root_count": len(jobs),
            "output_pair_count": len(jobs),
        },
        "all_exit_code_zero": True,
        "artifact_health_pass": True,
        "exact_input_once_candidate": True,
        "full_training_authority": 0,
        "full_extraction_authority": False,
        "science_freeze_authority": False,
        "broad_production_authority": False,
        "canonical_promotion": False,
        "next_gate": (
            "logical aggregate plus independent full-scope matrix population closure"
        ),
    }


def validate_terminal_receipt(
    payload: Mapping[str, Any],
    *,
    execution_artifact: Mapping[str, Any],
    submission_artifact: Mapping[str, Any],
    snapshot_artifact: Mapping[str, Any],
) -> None:
    require_exact_keys(
        payload,
        {
            "schema",
            "status",
            "artifact_profile",
            "bindings",
            "counts",
            "all_exit_code_zero",
            "artifact_health_pass",
            "exact_input_once_candidate",
            "full_training_authority",
            "full_extraction_authority",
            "science_freeze_authority",
            "broad_production_authority",
            "canonical_promotion",
            "next_gate",
        },
        "terminal receipt",
    )
    if (
        payload.get("schema") != TERMINAL_RECEIPT_SCHEMA
        or payload.get("status") != TERMINAL_STATUS
        or payload.get("artifact_profile") != ARTIFACT_PROFILE
        or payload.get("all_exit_code_zero") is not True
        or payload.get("artifact_health_pass") is not True
        or payload.get("exact_input_once_candidate") is not True
        or payload.get("full_training_authority") != 0
        or payload.get("full_extraction_authority") is not False
        or payload.get("science_freeze_authority") is not False
        or payload.get("broad_production_authority") is not False
        or payload.get("canonical_promotion") is not False
    ):
        raise ControllerError("terminal receipt schema/status/authority differs")
    bindings = require_mapping(payload.get("bindings"), "terminal bindings")
    expected = {
        "execution_manifest": execution_artifact["sha256"],
        "submission_receipt": submission_artifact["sha256"],
        "terminal_snapshot": snapshot_artifact["sha256"],
    }
    for name, sha256 in expected.items():
        record = require_mapping(bindings.get(name), f"terminal binding {name}")
        if record.get("sha256") != sha256:
            raise ControllerError(f"terminal binding {name} differs")
    counts = require_mapping(payload.get("counts"), "terminal counts")
    if (
        counts.get("row_count") != EXPECTED_ROW_COUNT
        or counts.get("job_count") != EXPECTED_JOB_COUNT
        or counts.get("ephemeral_analysis_health_count") != EXPECTED_JOB_COUNT
        or counts.get("retained_analysis_root_count") != 0
        or counts.get("training_sidecar_root_count") != EXPECTED_JOB_COUNT
        or counts.get("output_pair_count") != EXPECTED_JOB_COUNT
    ):
        raise ControllerError("terminal receipt counts differ")


def build_logical_aggregate(
    *,
    execution: Mapping[str, Any],
    jobs: Sequence[Mapping[str, Any]],
    execution_artifact: Mapping[str, Any],
    submission_artifact: Mapping[str, Any],
    snapshot_artifact: Mapping[str, Any],
    terminal_artifact: Mapping[str, Any],
) -> dict[str, bytes]:
    row_by_id = {row["row_id"]: row for row in execution["_rows"]}
    inputs_by_system: dict[str, list[str]] = {"pp": [], "auau": []}
    provenance_by_system: dict[str, list[dict[str, Any]]] = {
        "pp": [],
        "auau": [],
    }
    selected_by_source: dict[str, dict[str, int]] = defaultdict(
        lambda: {view: 0 for view in SHOWER_VIEWS}
    )
    occurrences_by_source: dict[str, set[str]] = defaultdict(set)
    expected_counts_by_source: dict[str, int] = defaultdict(int)

    for job in jobs:
        row = row_by_id[job["row_id"]]
        system = row["system"]
        sample = row["sample"]
        sidecar = job["training_sidecar_root"]
        health = job["sidecar_health"]
        path = sidecar["path"]
        inputs_by_system[system].append(path)
        record = {
            "path": path,
            "row_id": row["row_id"],
            "system": system,
            "source_sample": sample,
            "source_occurrence_id_hex": health["source_occurrence_id_hex"],
            "input_uri_sha256": health["input_uri_sha256"],
            "input_file_sha256": health["input_file_sha256"],
            "source_manifest_sha256": row["full_source_manifest_sha256"],
            "config_sha256": row["config_sha256"],
            "code_sha256": row["code_sha256"],
            "training_view_root_sha256": sidecar["sha256"],
        }
        provenance_by_system[system].append(record)
        expected_counts_by_source[sample] += 1
        occurrences_by_source[sample].add(health["source_occurrence_id_hex"])
        for view in SHOWER_VIEWS:
            selected_by_source[sample][view] += health[
                "selected_training_rows_by_view"
            ][view]

    source_authority_by_system: dict[str, dict[str, Any]] = {"pp": {}, "auau": {}}
    for row in execution["_rows"]:
        system = row["system"]
        source = row["sample"]
        counts = selected_by_source[source]
        if len(set(counts.values())) != 1:
            raise ControllerError(
                f"{source} selected training population differs across shower views"
            )
        records = [
            record
            for record in provenance_by_system[system]
            if record["source_sample"] == source
        ]
        expected_inputs = expected_counts_by_source[source]
        expected_occurrences = len(occurrences_by_source[source])
        if (
            len(records) != expected_inputs
            or expected_inputs <= 0
            or expected_occurrences != expected_inputs
        ):
            raise ControllerError(f"{source} exact input/occurrence closure differs")
        selected_rows = next(iter(counts.values()))
        source_authority_by_system[system][source] = {
            "state": (
                "SOURCE_COMPLETE"
                if selected_rows > 0
                else "SOURCE_COMPLETE_ZERO_IN_DOMAIN"
            ),
            "full_source_manifest_sha256": row[
                "full_source_manifest_sha256"
            ],
            "expected_input_count": expected_inputs,
            "expected_occurrence_count": expected_occurrences,
            "input_records_sha256": source_input_records_sha256(records),
        }

    artifacts: dict[str, bytes] = {}
    for system in ("pp", "auau"):
        paths = inputs_by_system[system]
        if len(paths) != len(set(paths)):
            raise ControllerError(
                f"{system} logical aggregate duplicates an input path"
            )
        artifacts[f"{system}_sidecars.list"] = (
            "\n".join(paths) + "\n"
        ).encode("utf-8")
        provenance = {
            "schema": "THE134_SOURCE_PROVENANCE_V1",
            "inputs": provenance_by_system[system],
            "source_coverage_authority": source_authority_by_system[system],
        }
        artifacts[f"{system}_source_provenance.json"] = canonical_json_bytes(
            provenance
        )

    inventory = {
        name: {
            "sha256": hashlib.sha256(data).hexdigest(),
            "size_bytes": len(data),
        }
        for name, data in artifacts.items()
    }
    aggregate = {
        "schema": AGGREGATE_SCHEMA,
        "status": AGGREGATE_STATUS,
        "artifact_profile": ARTIFACT_PROFILE,
        "bindings": {
            "execution_manifest": dict(execution_artifact),
            "submission_receipt": dict(submission_artifact),
            "terminal_snapshot": dict(snapshot_artifact),
            "terminal_receipt": dict(terminal_artifact),
        },
        "counts": {
            "input_root_count": len(jobs),
            "pp_input_root_count": len(inputs_by_system["pp"]),
            "auau_input_root_count": len(inputs_by_system["auau"]),
            "unique_input_root_count": len(
                set(inputs_by_system["pp"] + inputs_by_system["auau"])
            ),
        },
        "artifacts": inventory,
        "merge_semantics": {
            "mode": "EXACT_INPUT_ONCE_LOGICAL_MANIFEST",
            "physical_hadd_performed": False,
            "reason": (
                "Each sidecar owns one source-occurrence identity and per-file metadata"
            ),
        },
        "source_path_routing": {
            "external_provenance_is_authoritative": True,
            "path_token_inference_required": False,
            "consumer_requirement": (
                "prepare_the134_h70_matrix.py must resolve source from "
                "THE134_SOURCE_PROVENANCE_V1 instead of path tokens"
            ),
        },
        "full_training_authority": 0,
        "full_extraction_authority": False,
        "science_freeze_authority": False,
        "broad_production_authority": False,
        "canonical_promotion": False,
        "next_gate": (
            "independent prepare_the134_h70_matrix.py --scope full validation "
            "for both systems and all seven views"
        ),
    }
    artifacts["logical_aggregate_receipt.json"] = canonical_json_bytes(aggregate)
    return artifacts


def write_aggregate_directory(root: Path, artifacts: Mapping[str, bytes]) -> None:
    destination = safe_fresh_local_output(root, "aggregate output directory")
    destination.mkdir(
        mode=SHARED_DIRECTORY_MODE,
        parents=False,
        exist_ok=False,
    )
    destination.chmod(SHARED_DIRECTORY_MODE)
    try:
        for name in sorted(artifacts):
            write_new_bytes(destination / name, artifacts[name])
    except Exception:
        # Evidence is intentionally preserved on any partial write.  No cleanup
        # or overwrite is attempted by this controller.
        raise


def add_upstream_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--plan-sha256", required=True)
    parser.add_argument("--preflight-receipt", type=Path, required=True)
    parser.add_argument("--preflight-receipt-sha256", required=True)
    parser.add_argument("--storage-certificate", type=Path, required=True)
    parser.add_argument("--storage-certificate-sha256", required=True)
    parser.add_argument("--controller-budget", type=Path, required=True)
    parser.add_argument("--controller-budget-sha256", required=True)
    parser.add_argument(
        "--artifact-profile", choices=(ARTIFACT_PROFILE,), required=True
    )
    parser.add_argument(
        "--validation-root",
        type=Path,
        required=True,
        help="safe absolute path that does not exist and is never created",
    )
    parser.add_argument("--dry-staging-root", type=Path, required=True)
    parser.add_argument("--dry-manifest-sha256", required=True)


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    actions = parser.add_subparsers(dest="action", required=True)

    preflight = actions.add_parser("preflight", help="read-only exact preflight")
    add_upstream_args(preflight)

    seal = actions.add_parser(
        "seal", help="seal an externally authorized exact execution manifest"
    )
    add_upstream_args(seal)
    seal.add_argument("--authorization-receipt", type=Path, required=True)
    seal.add_argument("--authorization-receipt-sha256", required=True)
    seal.add_argument("--output", type=Path, required=True)

    submit = actions.add_parser(
        "submit", help="invoke only the pinned existing submitter for all rows"
    )
    submit.add_argument("--execution-manifest", type=Path, required=True)
    submit.add_argument("--execution-manifest-sha256", required=True)
    submit.add_argument("--safe-resume-certificate", type=Path)
    submit.add_argument("--safe-resume-certificate-sha256")
    submit.add_argument("--receipt-dir", type=Path, required=True)
    submit.add_argument(
        "--execute",
        action="store_true",
        help="required explicit mutation switch; omission always fails closed",
    )

    terminal = actions.add_parser(
        "terminal", help="validate one exact externally captured terminal snapshot"
    )
    terminal.add_argument("--execution-manifest", type=Path, required=True)
    terminal.add_argument("--execution-manifest-sha256", required=True)
    terminal.add_argument("--submission-receipt", type=Path, required=True)
    terminal.add_argument("--submission-receipt-sha256", required=True)
    terminal.add_argument("--terminal-snapshot", type=Path, required=True)
    terminal.add_argument("--terminal-snapshot-sha256", required=True)
    terminal.add_argument("--output", type=Path, required=True)

    aggregate = actions.add_parser(
        "aggregate",
        help="write exact-input-once logical manifests; never physically hadd ROOT",
    )
    aggregate.add_argument("--execution-manifest", type=Path, required=True)
    aggregate.add_argument("--execution-manifest-sha256", required=True)
    aggregate.add_argument("--submission-receipt", type=Path, required=True)
    aggregate.add_argument("--submission-receipt-sha256", required=True)
    aggregate.add_argument("--terminal-snapshot", type=Path, required=True)
    aggregate.add_argument("--terminal-snapshot-sha256", required=True)
    aggregate.add_argument("--terminal-receipt", type=Path, required=True)
    aggregate.add_argument("--terminal-receipt-sha256", required=True)
    aggregate.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args(argv)


def load_execution(
    args: argparse.Namespace,
    *,
    load_partition_chunks: bool = True,
) -> tuple[dict[str, Any], dict[str, Any]]:
    payload, artifact = load_pinned_json(
        args.execution_manifest,
        args.execution_manifest_sha256,
        "execution manifest",
    )
    return (
        validate_execution_manifest(
            payload,
            load_partition_chunks=load_partition_chunks,
        ),
        artifact,
    )


def main(argv: Iterable[str] | None = None) -> int:
    try:
        args = parse_args(argv)
        if args.action in {"preflight", "seal"}:
            context = materializer_context(args)
            dry_stage = validate_dry_stage(
                context, args.dry_staging_root, args.dry_manifest_sha256
            )
            if args.action == "preflight":
                result = {
                    "schema": EXECUTION_SCHEMA,
                    "status": "PASS_PREFLIGHT_NO_WRITE",
                    "artifact_profile": ARTIFACT_PROFILE,
                    "row_count": EXPECTED_ROW_COUNT,
                    "job_count": EXPECTED_JOB_COUNT,
                    "group_size": EXPECTED_GROUP_SIZE,
                    "request_memory_mb": EXPECTED_REQUEST_MEMORY_MB,
                    "dry_materialization_manifest_sha256": dry_stage["manifest"][
                        "sha256"
                    ],
                    "submission_performed": False,
                    "submission_authority": False,
                    "full_training_authority": 0,
                    "full_extraction_authority": False,
                }
                print(json.dumps(result, sort_keys=True))
                return 0
            authorization_payload, authorization_artifact = load_pinned_json(
                args.authorization_receipt,
                args.authorization_receipt_sha256,
                "submission authorization",
            )
            authorization = validate_authorization(
                authorization_payload,
                artifact=authorization_artifact,
                context=context,
                dry_stage=dry_stage,
            )
            execution = build_execution_manifest(context, dry_stage, authorization)
            write_new_json(args.output, execution)
            print(
                json.dumps(
                    {
                        "schema": EXECUTION_SCHEMA,
                        "status": READY_STATUS,
                        "output": str(args.output.absolute()),
                        "output_sha256": file_sha256(args.output.absolute()),
                        "row_count": EXPECTED_ROW_COUNT,
                        "job_count": EXPECTED_JOB_COUNT,
                        "submission_performed": False,
                    },
                    sort_keys=True,
                )
            )
            return 0

        execution, execution_artifact = load_execution(
            args,
            load_partition_chunks=args.action != "submit",
        )
        if args.action == "submit":
            if not args.execute:
                raise ControllerError(
                    "submit requires --execute after exact external authorization"
                )
            if (
                args.safe_resume_certificate is None
                or args.safe_resume_certificate_sha256 is None
            ):
                raise ControllerError(
                    "submit requires one exact SDCC safe-resume certificate"
                )
            receipt = execute_submission(
                execution,
                execution_artifact,
                args.safe_resume_certificate,
                args.safe_resume_certificate_sha256,
                args.receipt_dir,
            )
            print(json.dumps(receipt, sort_keys=True))
            return 0

        submission_payload, submission_artifact = load_pinned_json(
            args.submission_receipt,
            args.submission_receipt_sha256,
            "submission receipt",
        )
        submission = validate_submission_receipt(
            submission_payload, execution, execution_artifact
        )
        snapshot_payload, snapshot_artifact = load_pinned_json(
            args.terminal_snapshot,
            args.terminal_snapshot_sha256,
            "terminal snapshot",
        )
        _, jobs = validate_terminal_snapshot(
            snapshot_payload,
            execution=execution,
            execution_artifact=execution_artifact,
            submission=submission,
            submission_artifact=submission_artifact,
        )
        if args.action == "terminal":
            receipt = terminal_receipt(
                execution_artifact=execution_artifact,
                submission_artifact=submission_artifact,
                snapshot_artifact=snapshot_artifact,
                jobs=jobs,
            )
            write_new_json(args.output, receipt)
            print(json.dumps(receipt, sort_keys=True))
            return 0

        terminal_payload, terminal_artifact = load_pinned_json(
            args.terminal_receipt,
            args.terminal_receipt_sha256,
            "terminal receipt",
        )
        validate_terminal_receipt(
            terminal_payload,
            execution_artifact=execution_artifact,
            submission_artifact=submission_artifact,
            snapshot_artifact=snapshot_artifact,
        )
        artifacts = build_logical_aggregate(
            execution=execution,
            jobs=jobs,
            execution_artifact=execution_artifact,
            submission_artifact=submission_artifact,
            snapshot_artifact=snapshot_artifact,
            terminal_artifact=terminal_artifact,
        )
        write_aggregate_directory(args.output_dir, artifacts)
        print(
            json.dumps(
                {
                    "schema": AGGREGATE_SCHEMA,
                    "status": AGGREGATE_STATUS,
                    "output_dir": str(args.output_dir.absolute()),
                    "artifact_count": len(artifacts),
                    "input_root_count": len(jobs),
                    "physical_hadd_performed": False,
                },
                sort_keys=True,
            )
        )
        return 0
    except (ControllerError, FileNotFoundError, OSError) as exc:
        print(f"[THE134-FULL-EXECUTION][ERROR] {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
