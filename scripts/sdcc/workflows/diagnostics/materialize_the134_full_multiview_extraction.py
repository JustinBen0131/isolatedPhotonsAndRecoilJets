#!/usr/bin/env python3
"""Preflight and locally stage the full THE-134 multiview extraction.

This controller consumes the exact output of
``resolve_the134_full_multiview_extraction.py`` plus the passed pre-extraction
storage certificate and its controller dry-materialization budget.  It has no
submission action.  ``preflight`` is read-only.  ``materialize`` writes only
deterministic controller artifacts below one fresh, caller-provided local
staging root.

The staged submit description is deliberately inert text, not an executable
Condor submit file.  No action invokes Condor, a shell, SSH, a network API, job
control, or any remote namespace creation.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import re
import sys
from pathlib import Path, PurePosixPath
from typing import Any, Iterable, Mapping


HERE = Path(__file__).resolve().parent
RESOLVER_PATH = HERE / "resolve_the134_full_multiview_extraction.py"
PROJECTOR_PATH = HERE / "project_the134_preextraction_storage_quota.py"


def _load_local_module(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load local module: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


resolver = _load_local_module("the134_full_controller_resolver", RESOLVER_PATH)
projector = _load_local_module(
    "the134_full_controller_storage_projector", PROJECTOR_PATH
)


CONTROLLER_SCHEMA = "THE134_FULL_MULTIVIEW_EXTRACTION_CONTROLLER_V1"
ROW_STAGE_SCHEMA = "THE134_FULL_MULTIVIEW_STAGED_ROW_V1"
JOB_STAGE_SCHEMA = "THE134_FULL_MULTIVIEW_STAGED_JOB_V1"
MANIFEST_SCHEMA = "THE134_FULL_MULTIVIEW_DRY_MATERIALIZATION_MANIFEST_V1"
RESOLVED_ARTIFACT_PROFILE = dict(resolver.SIDECAR_ONLY_ARTIFACT_PROFILE)
ARTIFACT_PROFILE = str(RESOLVED_ARTIFACT_PROFILE["artifact_profile"])
AUTHORITY_STATE = "LOCAL_DRY_MATERIALIZATION_ONLY_NO_SUBMISSION"

EXPECTED_ROW_COUNT = 13
EXPECTED_SOURCE_TUPLES = 129_998
EXPECTED_GROUP_SIZE = 7
EXPECTED_JOB_COUNT = 18_577
EXPECTED_OUTPUT_PAIRS = 18_577
EXPECTED_ANALYSIS_OUTPUTS = 18_577
EXPECTED_SIDECAR_OUTPUTS = 18_577
EXPECTED_PHYSICAL_ARTIFACTS = 37_154
EXPECTED_RETAINED_ANALYSIS_OUTPUTS = 0
EXPECTED_DURABLE_ROOT_ARTIFACTS = EXPECTED_SIDECAR_OUTPUTS
EXPECTED_REQUEST_MEMORY_MB = 8_000

OUTPUT_FILENAMES = (
    "the134_full_multiview_rows.jsonl",
    "the134_full_multiview_jobs.jsonl",
    "the134_full_multiview_submit_description.txt",
    "materialization_manifest.json",
)

SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
ENV_KEY_RE = re.compile(r"^[A-Z][A-Z0-9_]*$")
SAFE_SDCC_USER_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,63}$")
SDCC_USER_ALIAS_ROOT = Path("/sphenix/u")
SDCC_CANONICAL_USER_ROOT = Path("/gpfs/mnt/gpfs02/sphenix/user")

FORBIDDEN_POSITIVE_AUTHORITY_FIELDS = frozenset(
    {
        "submission_authority",
        "full_extraction_authority",
        "science_freeze_authority",
        "broad_production_authority",
        "canonical_promotion",
        "physics_output_authority",
        "physics_authority",
        "replay_cache_authority",
        "cache_authority",
        "the121_authority",
        "the122_authority",
    }
)

CONTROLLER_AUTHORITY = {
    "state": AUTHORITY_STATE,
    "artifact_profile": ARTIFACT_PROFILE,
    "submission_performed": False,
    "submission_authority": False,
    "full_training_authority": 0,
    "full_extraction_authority": False,
    "science_freeze_authority": False,
    "broad_production_authority": False,
    "canonical_promotion": False,
    "physics_output_authority": False,
    "replay_cache_authority": False,
    "the121_authority": False,
    "the122_authority": False,
}

PROVENANCE_PLACEHOLDERS = {
    "RJ_CODEX_CHAT_NAME": "__REQUIRED_ONLY_AFTER_EXPLICIT_SUBMISSION_APPROVAL__",
    "RJ_CODEX_THREAD_ID": "__REQUIRED_ONLY_AFTER_EXPLICIT_SUBMISSION_APPROVAL__",
}


class ControllerError(RuntimeError):
    """Fail-closed controller contract violation."""


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


def canonical_sha256(payload: Any) -> str:
    return hashlib.sha256(
        json.dumps(
            payload,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
        ).encode("utf-8")
    ).hexdigest()


RESOLVED_ARTIFACT_PROFILE_SHA256 = canonical_sha256(
    RESOLVED_ARTIFACT_PROFILE
)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
    except OSError as exc:
        raise ControllerError(f"cannot read pinned artifact: {path}") from exc
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


def require_positive_int(value: object, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ControllerError(f"{label} must be a positive integer")
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
    if not absolute.is_file():
        raise ControllerError(f"{label} is missing: {absolute}")
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


def validate_authority_occurrences(payload: Any, label: str) -> None:
    if isinstance(payload, dict):
        for key, value in payload.items():
            path = f"{label}.{key}"
            if key == "submission_performed" and value is not False:
                raise ControllerError(f"{path} must remain false")
            if key == "full_training_authority" and value != 0:
                raise ControllerError(f"{path} must remain zero")
            if key in FORBIDDEN_POSITIVE_AUTHORITY_FIELDS and value not in (
                False,
                0,
                None,
            ):
                raise ControllerError(f"{path} escalates forbidden authority")
            validate_authority_occurrences(value, path)
    elif isinstance(payload, list):
        for index, value in enumerate(payload):
            validate_authority_occurrences(value, f"{label}[{index}]")


def validate_remote_root(value: object, label: str) -> str:
    if (
        not isinstance(value, str)
        or not value.startswith("/")
        or value == "/"
        or "\x00" in value
        or ".." in PurePosixPath(value).parts
        or value != value.strip()
    ):
        raise ControllerError(f"{label} must be a safe absolute namespace")
    return value.rstrip("/")


def resolve_standard_sdcc_scratch_alias(path: Path) -> Path | None:
    """Resolve only SDCC's exact per-user scratch alias to canonical GPFS."""

    try:
        relative = path.relative_to(SDCC_USER_ALIAS_ROOT)
    except ValueError:
        return None
    parts = relative.parts
    if (
        len(parts) < 3
        or parts[1] != "scratch"
        or SAFE_SDCC_USER_RE.fullmatch(parts[0]) is None
    ):
        raise ControllerError("staging root is not a valid SDCC scratch alias path")
    user = parts[0]
    alias_root = SDCC_USER_ALIAS_ROOT / user / "scratch"
    canonical_root = SDCC_CANONICAL_USER_ROOT / user
    expected = canonical_root.joinpath(*parts[2:])
    try:
        if not alias_root.is_symlink():
            raise ControllerError("SDCC scratch alias root is not a symlink")
        alias_target = alias_root.resolve(strict=True)
        canonical_target = canonical_root.resolve(strict=True)
    except OSError as exc:
        raise ControllerError("SDCC scratch alias root cannot be resolved") from exc
    if alias_target != canonical_root or canonical_target != canonical_root:
        raise ControllerError(
            "SDCC scratch alias does not target the exact user GPFS root"
        )
    first_resolution = path.resolve(strict=False)
    second_resolution = path.resolve(strict=False)
    if first_resolution != expected or second_resolution != expected:
        raise ControllerError(
            "SDCC scratch alias resolution is unstable or escapes GPFS"
        )
    ancestor = expected.parent
    while not os.path.lexists(ancestor):
        if ancestor == canonical_root:
            break
        ancestor = ancestor.parent
    if canonical_root not in (ancestor, *ancestor.parents):
        raise ControllerError("SDCC scratch target escapes its user root")
    cursor = ancestor
    while cursor != canonical_root:
        if cursor.is_symlink():
            raise ControllerError(
                "canonical GPFS staging parent chain contains a symlink"
            )
        cursor = cursor.parent
    return expected


def validate_staging_root(path: Path) -> Path:
    if not path.is_absolute():
        raise ControllerError("staging root must be an absolute local path")
    if (
        path == Path("/")
        or path.name in {"", ".", ".."}
        or ".." in path.parts
    ):
        raise ControllerError("staging root is unsafe")
    if os.path.lexists(path):
        raise ControllerError(f"staging root already exists: {path}")
    resolved_alias = resolve_standard_sdcc_scratch_alias(path)
    if resolved_alias is not None:
        if not path.parent.is_dir():
            raise ControllerError(
                f"staging-root parent must already exist locally: {path.parent}"
            )
        if os.path.lexists(resolved_alias):
            raise ControllerError(
                f"resolved staging root already exists: {resolved_alias}"
            )
        return resolved_alias
    parent = path.parent
    if not parent.is_dir():
        raise ControllerError(
            f"staging-root parent must already exist locally: {parent}"
        )
    current = parent
    while True:
        if current.is_symlink():
            raise ControllerError(
                f"staging-root parent/ancestor is a symlink: {current}"
            )
        if current == current.parent:
            break
        current = current.parent
    resolved_parent = parent.resolve(strict=True)
    resolved_leaf = resolved_parent / path.name
    if os.path.lexists(resolved_leaf):
        raise ControllerError(
            f"resolved staging root already exists: {resolved_leaf}"
        )
    return resolved_leaf


def paths_overlap(first: Path, second: Path) -> bool:
    first_absolute = first.resolve(strict=False)
    second_absolute = second.resolve(strict=False)
    return (
        first_absolute == second_absolute
        or first_absolute in second_absolute.parents
        or second_absolute in first_absolute.parents
    )


def validate_bundle(
    payload: dict[str, Any],
) -> tuple[dict[str, Any], dict[str, dict[str, Any]]]:
    if (
        payload.get("schema") != resolver.BUNDLE_SCHEMA
        or payload.get("status") != "PASS"
    ):
        raise ControllerError("immutable bundle schema/status differs")
    for field in (
        "bundle_identity_sha256",
        "semantic_fingerprint_sha256",
        "code_sha256",
        "replay_schema_sha256",
        "training_schema_sha256",
        "semantic_sha256",
    ):
        require_sha256(payload.get(field), f"bundle.{field}")
    runtime = require_mapping(payload.get("runtime"), "bundle.runtime")
    if runtime.get("request_memory_mb") != EXPECTED_REQUEST_MEMORY_MB:
        raise ControllerError("bundle request_memory_mb differs")
    raw_artifacts = require_sequence(payload.get("artifacts"), "bundle.artifacts")
    by_role: dict[str, dict[str, Any]] = {}
    normalized: list[dict[str, Any]] = []
    for raw in raw_artifacts:
        record = require_mapping(raw, "bundle artifact")
        role = record.get("role")
        if not isinstance(role, str) or not role or role in by_role:
            raise ControllerError("bundle artifact roles are invalid or duplicated")
        artifact_path = Path(str(record.get("path", ""))).absolute()
        expected_sha = require_sha256(
            record.get("sha256"), f"bundle artifact {role}.sha256"
        )
        if not artifact_path.is_file():
            raise ControllerError(f"bundle artifact is missing: role={role}")
        observed_sha = file_sha256(artifact_path)
        if observed_sha != expected_sha:
            raise ControllerError(f"bundle artifact hash drift: role={role}")
        size = require_positive_int(
            record.get("size_bytes"), f"bundle artifact {role}.size_bytes"
        )
        if artifact_path.stat().st_size != size:
            raise ControllerError(f"bundle artifact size drift: role={role}")
        normalized_record = {
            "role": role,
            "path": str(artifact_path),
            "resolved_path": str(artifact_path.resolve(strict=True)),
            "sha256": observed_sha,
            "size_bytes": size,
        }
        by_role[role] = normalized_record
        normalized.append(normalized_record)
    required_roles = {
        "submitter",
        "pp_executor",
        "auau_executor",
        "pp_config",
        "auau_config",
        "pp_library",
        "auau_library",
        "pp_model",
        "auau_model",
        "photon_cluster_builder_header",
        "calo_reco_library",
        "release_calo_io",
        "release_clusteriso",
        "release_jetbase",
    }
    missing = sorted(required_roles - set(by_role))
    if missing:
        raise ControllerError(f"bundle lacks controller-required roles: {missing}")
    normalized_bundle = {
        "public_commit": payload.get("public_commit"),
        "bundle_identity_sha256": payload["bundle_identity_sha256"],
        "semantic_fingerprint_sha256": payload[
            "semantic_fingerprint_sha256"
        ],
        "code_sha256": payload["code_sha256"],
        "replay_schema_sha256": payload["replay_schema_sha256"],
        "training_schema_sha256": payload["training_schema_sha256"],
        "semantic_sha256": payload["semantic_sha256"],
        "runtime": dict(runtime),
        "artifact_by_role": by_role,
        "artifacts": normalized,
    }
    return normalized_bundle, by_role


def validate_materialization_binding(
    payload: dict[str, Any],
    *,
    materialization_path: Path,
    bundle_path: Path,
    bundle_sha256: str,
    bundle: Mapping[str, Any],
    plan_record: Mapping[str, Any],
) -> None:
    try:
        validated = (
            resolver.bundle_materializer.validate_materialization_payload(
                payload,
                materialization_receipt_path=materialization_path,
                rehash=True,
            )
        )
    except resolver.bundle_materializer.MaterializationError as exc:
        raise ControllerError(
            f"immutable materialization readback failed: {exc}"
        ) from exc
    bound_bundle = Path(
        str(validated["resolver_bundle_receipt"])
    ).resolve(strict=True)
    expected_bundle = bundle_path.resolve(strict=True)
    if bound_bundle != expected_bundle:
        raise ControllerError("materialization is bound to a different bundle")
    if validated["resolver_bundle_receipt_sha256"] != bundle_sha256:
        raise ControllerError("materialization bundle SHA-256 differs")
    if (
        validated["bundle_identity_sha256"]
        != bundle.get("bundle_identity_sha256")
    ):
        raise ControllerError("materialization bundle identity differs")
    if plan_record.get("readback") != "PASS_SYMLINK_FREE_READONLY_CONTENT_EXACT":
        raise ControllerError("materialization readback profile differs")
    validated_digest_root = Path(
        str(validated["digest_named_bundle_path"])
    ).resolve(strict=True)
    plan_digest_root = Path(
        str(plan_record.get("digest_named_bundle_path", ""))
    ).resolve(strict=True)
    if plan_digest_root != validated_digest_root:
        raise ControllerError("materialization digest-named path differs")
    if plan_record.get("digest_named_bundle_path") != validated.get(
        "digest_named_bundle_path"
    ):
        raise ControllerError("materialization digest-named path differs")
    if materialization_path.absolute() != Path(str(plan_record.get("path", ""))).absolute():
        raise ControllerError("plan materialization path differs")


def artifact_matches(
    observed: object,
    expected: Mapping[str, Any],
    label: str,
) -> None:
    record = require_mapping(observed, label)
    for field in ("path", "sha256"):
        if record.get(field) != expected.get(field):
            raise ControllerError(f"{label}.{field} differs from bundle")


def expected_worker_environment(
    row: Mapping[str, Any],
    *,
    campaign: Mapping[str, str],
    training_period: str,
    bundle: Mapping[str, Any],
) -> dict[str, str]:
    input_contract = require_mapping(row.get("input_contract"), "row input")
    execution = require_mapping(row.get("execution_contract"), "row execution")
    source = {
        "full_source_manifest_sha256": input_contract[
            "full_source_manifest_sha256"
        ]
    }
    row_inventory = {
        key: str(row[key])
        for key in (
            "row_id",
            "system",
            "lane",
            "dataset",
            "sample",
            "source_role",
            "minimum_bias_gate",
            "photon_id_row_match",
        )
    }
    return resolver.runtime_environment(
        row_inventory,
        tag=campaign["tag"],
        pp_period=training_period,
        bundle=dict(bundle),
        source=source,
        sidecar_template=str(execution["multiview_sidecar_template"]),
    )


def expected_materialization_environment(
    row: Mapping[str, Any],
    *,
    campaign: Mapping[str, str],
    bundle: Mapping[str, Any],
    by_role: Mapping[str, Mapping[str, Any]],
    worker_environment: Mapping[str, str],
) -> dict[str, str]:
    system = str(row["system"])
    row_id = str(row["row_id"])
    input_contract = require_mapping(row["input_contract"], "row input contract")
    input_lists = require_sequence(
        input_contract.get("input_lists"), f"{row_id}.input_lists"
    )
    if not input_lists:
        raise ControllerError(f"{row_id} has no input-list records")
    first_list = require_mapping(input_lists[0], f"{row_id}.input_lists[0]")
    sample_root = Path(str(first_list.get("path", ""))).parent
    sim_root = str(sample_root.parent)
    output_namespace = f"{campaign['output_root']}/{row_id}"
    environment = {
        "RJ_DAG_DRYRUN": "1",
        "RJ_CONDOR_SEALED_ENVIRONMENT": "1",
        "RJ_SIM_ROOT_OVERRIDE": sim_root,
        "RJ_CONFIG_YAML": str(by_role[f"{system}_config"]["path"]),
        (
            "RJ_PP_LIBRARY_OVERRIDE"
            if system == "pp"
            else "RJ_AUAU_LIBRARY_OVERRIDE"
        ): str(by_role[f"{system}_library"]["path"]),
        "RJ_PHOTON_CLUSTER_BUILDER_HEADER_OVERRIDE": str(
            by_role["photon_cluster_builder_header"]["path"]
        ),
        "RJ_CALO_RECO_LIBRARY_OVERRIDE": str(
            by_role["calo_reco_library"]["path"]
        ),
        "RJ_PHOTON_CLUSTER_BUILDER_LIBRARY_OVERRIDE": "",
        "RJ_PINNED_CALO_RECO_RELEASE_COMPANIONS": "1",
        "RJ_PINNED_CALO_RECO_SONAME": str(
            bundle["runtime"]["calo_reco_soname"]
        ),
        "RJ_PINNED_RELEASE_NAME": str(bundle["runtime"]["release"]),
        "RJ_PINNED_OFFLINE_MAIN": str(bundle["runtime"]["offline_main"]),
        "RJ_PINNED_RELEASE_CALO_IO_PATH": str(
            by_role["release_calo_io"]["path"]
        ),
        "RJ_PINNED_RELEASE_CALO_IO_SHA256": str(
            by_role["release_calo_io"]["sha256"]
        ),
        "RJ_PINNED_RELEASE_CLUSTERISO_PATH": str(
            by_role["release_clusteriso"]["path"]
        ),
        "RJ_PINNED_RELEASE_CLUSTERISO_SHA256": str(
            by_role["release_clusteriso"]["sha256"]
        ),
        "RJ_PINNED_RELEASE_JETBASE_PATH": str(
            by_role["release_jetbase"]["path"]
        ),
        "RJ_PINNED_RELEASE_JETBASE_SHA256": str(
            by_role["release_jetbase"]["sha256"]
        ),
        "RJ_RELEASE_CORE_LIB_DIR": str(bundle["runtime"]["release_core_lib_dir"]),
        "RJ_RELEASE_CORE_LIB64_DIR": str(
            bundle["runtime"]["release_core_lib64_dir"]
        ),
        "RJ_AUTO_MERGE": "0",
        "RJ_STAGE_EMAIL_MODE": "none",
        "RJ_CLEAN_OUTPUT_BASE": "0",
        "RJ_REQUEST_MEMORY": f"{EXPECTED_REQUEST_MEMORY_MB}MB",
        "RJ_REQUIRE_NON_TINY_OUTPUT": "1",
        "RJ_MIN_OUTPUT_BYTES": "50000",
        "RJ_FAIL_ON_MISSING_CALO_INPUT": "1",
        "RJ_VALIDATE_SIM_INPUT_PATHS": "1",
        "RJ_PROFILE_JOB": "1",
        "RJ_JOB_HEARTBEAT_SECONDS": "120",
        "RJ_DEST_BASE_OVERRIDE": output_namespace,
        "RJ_SUBMISSION_NAMESPACE": row_id,
        "RJ_CONDOR_SUB_DIR": f"{campaign['submit_root']}/{row_id}",
        "RJ_PHOTON_ID_ROW_MATCH": str(row["photon_id_row_match"]),
        "RJ_ID_FANOUT_MAX_ROWS": "1",
        "RJ_SIM_ALLOW_NONE_LISTS": "0" if system == "pp" else "1",
        "RJ_SUBMIT_EXTRA_ENV": resolver.serialized_environment(
            dict(worker_environment)
        ),
    }
    if system == "pp":
        for key, value in worker_environment.items():
            if key.startswith("RJ_PPG12_") or key.startswith(
                "RJ_PP_PHOTONID_"
            ):
                environment[key] = value
    return dict(sorted(environment.items()))


def validate_environment(environment: object, label: str) -> dict[str, str]:
    raw = require_mapping(environment, label)
    normalized: dict[str, str] = {}
    for key, value in raw.items():
        if (
            not isinstance(key, str)
            or ENV_KEY_RE.fullmatch(key) is None
            or not isinstance(value, str)
            or any(character in value for character in ("\n", "\r", "\x00"))
        ):
            raise ControllerError(f"{label} contains an unsafe environment entry")
        normalized[key] = value
    return dict(sorted(normalized.items()))


def validate_rows(
    plan: Mapping[str, Any],
    *,
    campaign: Mapping[str, str],
    bundle: Mapping[str, Any],
    by_role: Mapping[str, Mapping[str, Any]],
    source_manifest: Mapping[str, Any],
) -> tuple[list[dict[str, Any]], dict[str, dict[str, Any]]]:
    rows = require_sequence(plan.get("rows"), "plan.rows")
    expected_inventory = resolver.inventory_rows()
    if len(rows) != EXPECTED_ROW_COUNT:
        raise ControllerError(
            f"logical source row count={len(rows)}, expected={EXPECTED_ROW_COUNT}"
        )
    expected_ids = [row["row_id"] for row in expected_inventory]
    observed_ids = [
        row.get("row_id") if isinstance(row, dict) else None for row in rows
    ]
    if observed_ids != expected_ids or len(set(observed_ids)) != EXPECTED_ROW_COUNT:
        raise ControllerError("logical source rows are missing, reordered, or duplicated")

    source_rows = require_sequence(
        source_manifest.get("rows"), "source manifest rows"
    )
    source_by_id: dict[str, dict[str, Any]] = {}
    for raw in source_rows:
        record = require_mapping(raw, "source manifest row")
        row_id = str(record.get("row_id", ""))
        if not row_id or row_id in source_by_id:
            raise ControllerError("source manifest row identities are duplicated")
        source_by_id[row_id] = record
    if set(source_by_id) != set(expected_ids):
        raise ControllerError("source manifest row closure differs")

    training_contract = require_mapping(
        plan.get("training_period_si_contract"),
        "training_period_si_contract",
    )
    training_period = str(training_contract.get("period", ""))
    if (
        training_contract
        != resolver.training_period_si_contract(training_period)
        or plan.get("training_period_si_contract_sha256")
        != canonical_sha256(training_contract)
    ):
        raise ControllerError("p+p training-period SI contract differs")

    normalized_rows: list[dict[str, Any]] = []
    by_id: dict[str, dict[str, Any]] = {}
    total_tuples = 0
    total_jobs = 0
    row_fingerprints: set[str] = set()
    expected_forbidden = [
        "RJ_PPG12_CROSSING_PERIOD",
        "RJ_PPG12_PHOTON_YIELD_MIX_WEIGHT",
        "RJ_PP_VERTEX_REWEIGHT_FILE",
        "RJ_PPG12_PHOTON_YIELD_TRUTH_VERTEX",
        "RJ_PPG12_PHOTON_YIELD_BUILDER_TRUTH_VERTEX",
        "RJ_PPG12_PHOTON_YIELD_RECO_TRUTH_VERTEX",
        "RJ_THE134_TRUTH_LABEL_DIAGNOSTIC_V1",
        "RJ_THE134_TRUTH_LABEL_COUNTER_DIAGNOSTIC_V1",
        "RJ_THE134_TRUTH_LABEL_PREWEIGHT_DIAGNOSTIC_V1",
    ]
    for expected, raw in zip(expected_inventory, rows):
        row = require_mapping(raw, f"plan row {expected['row_id']}")
        row_id = expected["row_id"]
        if row.get("schema") != resolver.ROW_SCHEMA:
            raise ControllerError(f"{row_id} schema differs")
        for field, value in expected.items():
            if row.get(field) != value:
                raise ControllerError(f"{row_id}.{field} differs")
        validate_authority_occurrences(row, row_id)
        if (
            row.get("requested_scope") != resolver.REQUESTED_SCOPE
            or row.get("full_training_authority") != 0
            or row.get("authority_state") != resolver.PREFLIGHT_AUTHORITY_STATE
        ):
            raise ControllerError(f"{row_id} preflight authority differs")
        source = source_by_id[row_id]
        input_contract = require_mapping(
            row.get("input_contract"), f"{row_id}.input_contract"
        )
        tuple_count = require_positive_int(
            input_contract.get("source_tuple_count"),
            f"{row_id}.source_tuple_count",
        )
        job_count = require_positive_int(
            input_contract.get("expected_job_count"),
            f"{row_id}.expected_job_count",
        )
        expected_jobs = (tuple_count + EXPECTED_GROUP_SIZE - 1) // EXPECTED_GROUP_SIZE
        equal_job_fields = (
            "expected_chunk_count",
            "expected_output_pair_count",
            "expected_analysis_output_count",
            "expected_sidecar_output_count",
            "expected_source_occurrence_count",
        )
        if (
            input_contract.get("group_size") != EXPECTED_GROUP_SIZE
            or input_contract.get("event_limit_per_job") != 0
            or job_count != expected_jobs
            or any(input_contract.get(field) != job_count for field in equal_job_fields)
            or input_contract.get("source_occurrences_per_output_pair") != 1
        ):
            raise ControllerError(f"{row_id} count/partition contract differs")
        for field in (
            "full_source_manifest_sha256",
            "tuple_records_sha256",
            "first_tuple_sha256",
            "last_tuple_sha256",
            "row_partition_sha256",
            "chunk_records_sha256",
            "partition_artifact_sha256",
        ):
            require_sha256(input_contract.get(field), f"{row_id}.{field}")
        for field in (
            "full_source_manifest_sha256",
            "tuple_records_sha256",
            "first_tuple_sha256",
            "last_tuple_sha256",
        ):
            if input_contract.get(field) != source.get(field):
                raise ControllerError(f"{row_id}.{field} differs from source authority")
        if tuple_count != source.get("tuple_count"):
            raise ControllerError(f"{row_id} tuple count differs from source authority")

        bundle_contract = require_mapping(
            row.get("bundle_contract"), f"{row_id}.bundle_contract"
        )
        for field in (
            "public_commit",
            "code_sha256",
            "replay_schema_sha256",
            "training_schema_sha256",
            "semantic_sha256",
        ):
            if bundle_contract.get(field) != bundle.get(field):
                raise ControllerError(f"{row_id} bundle {field} differs")
        system = expected["system"]
        artifact_matches(
            bundle_contract.get("library"),
            by_role[f"{system}_library"],
            f"{row_id}.library",
        )
        artifact_matches(
            bundle_contract.get("model"),
            by_role[f"{system}_model"],
            f"{row_id}.model",
        )
        artifact_matches(
            bundle_contract.get("config"),
            by_role[f"{system}_config"],
            f"{row_id}.config",
        )
        artifact_matches(
            bundle_contract.get("submitter"),
            by_role["submitter"],
            f"{row_id}.submitter",
        )
        artifact_matches(
            bundle_contract.get("executor"),
            by_role[f"{system}_executor"],
            f"{row_id}.executor",
        )

        execution = require_mapping(
            row.get("execution_contract"), f"{row_id}.execution_contract"
        )
        if (
            execution.get("schema") != resolver.EXECUTION_SCHEMA
            or execution.get("state") != "PREFLIGHT_ONLY_NO_CONDOR_MUTATION"
            or execution.get("materialization_requires") != "RJ_DAG_DRYRUN=1"
            or execution.get("required_execution_inputs")
            != ["RJ_CODEX_CHAT_NAME", "RJ_CODEX_THREAD_ID"]
            or execution.get("analysis_output_namespace")
            != f"{campaign['output_root']}/{row_id}"
            or execution.get("submit_namespace")
            != f"{campaign['submit_root']}/{row_id}"
            or execution.get("evidence_namespace")
            != f"{campaign['evidence_root']}/{row_id}"
            or execution.get("multiview_sidecar_template")
            != (
                f"{campaign['output_root']}/{row_id}/training_views/"
                "$(Cluster).$(Process).root"
            )
            or execution.get("artifact_profile")
            != RESOLVED_ARTIFACT_PROFILE
            or execution.get("forbidden_ambient_environment")
            != expected_forbidden
        ):
            raise ControllerError(f"{row_id} execution boundary differs")
        expected_argv = [
            by_role["submitter"]["path"],
            expected["dataset"],
            "condorDoAll",
            "groupSize",
            str(EXPECTED_GROUP_SIZE),
            f"SAMPLE={expected['sample']}",
        ]
        if execution.get("existing_submitter_argv") != expected_argv:
            raise ControllerError(f"{row_id} submitter argv differs")
        worker = validate_environment(
            execution.get("worker_environment"), f"{row_id}.worker_environment"
        )
        expected_worker = expected_worker_environment(
            row,
            campaign=campaign,
            training_period=training_period,
            bundle=bundle,
        )
        if worker != expected_worker:
            raise ControllerError(f"{row_id} exact worker environment differs")
        materialization_environment = validate_environment(
            execution.get("materialization_environment"),
            f"{row_id}.materialization_environment",
        )
        expected_materialization = expected_materialization_environment(
            row,
            campaign=campaign,
            bundle=bundle,
            by_role=by_role,
            worker_environment=worker,
        )
        if materialization_environment != expected_materialization:
            raise ControllerError(
                f"{row_id} exact materialization environment differs"
            )
        if (
            materialization_environment.get("RJ_REQUEST_MEMORY")
            != f"{EXPECTED_REQUEST_MEMORY_MB}MB"
            or worker.get("RJ_REQUEST_MEMORY")
            != f"{EXPECTED_REQUEST_MEMORY_MB}MB"
        ):
            raise ControllerError(f"{row_id} request memory differs")
        fingerprint = require_sha256(
            row.get("row_fingerprint_sha256"), f"{row_id}.row fingerprint"
        )
        fingerprint_payload = dict(row)
        del fingerprint_payload["row_fingerprint_sha256"]
        if fingerprint != canonical_sha256(fingerprint_payload):
            raise ControllerError(f"{row_id} row fingerprint differs")
        if fingerprint in row_fingerprints:
            raise ControllerError("duplicate row fingerprints are forbidden")
        row_fingerprints.add(fingerprint)
        total_tuples += tuple_count
        total_jobs += job_count
        normalized = dict(row)
        normalized["_expected_job_count"] = job_count
        normalized["_expected_tuple_count"] = tuple_count
        normalized["_worker_environment_sha256"] = canonical_sha256(worker)
        normalized["_materialization_environment_sha256"] = canonical_sha256(
            materialization_environment
        )
        normalized_rows.append(normalized)
        by_id[row_id] = normalized

    if total_tuples != EXPECTED_SOURCE_TUPLES:
        raise ControllerError(
            f"source tuple/occurrence count={total_tuples}, "
            f"expected={EXPECTED_SOURCE_TUPLES}"
        )
    if total_jobs != EXPECTED_JOB_COUNT:
        raise ControllerError(
            f"job count={total_jobs}, expected={EXPECTED_JOB_COUNT}"
        )
    return normalized_rows, by_id


def validate_partition(
    partition_path: Path,
    *,
    plan: Mapping[str, Any],
    rows_by_id: Mapping[str, Mapping[str, Any]],
) -> tuple[list[dict[str, Any]], dict[str, list[dict[str, Any]]]]:
    partition_contract = require_mapping(
        plan.get("execution_partition"), "execution_partition"
    )
    expected_counts = {
        "group_size": EXPECTED_GROUP_SIZE,
        "source_tuple_count": EXPECTED_SOURCE_TUPLES,
        "expected_chunk_count": EXPECTED_JOB_COUNT,
        "expected_job_count": EXPECTED_JOB_COUNT,
        "expected_output_pair_count": EXPECTED_OUTPUT_PAIRS,
        "expected_analysis_output_count": EXPECTED_ANALYSIS_OUTPUTS,
        "expected_sidecar_output_count": EXPECTED_SIDECAR_OUTPUTS,
        "expected_physical_root_artifact_count": EXPECTED_PHYSICAL_ARTIFACTS,
        "expected_retained_analysis_output_count": (
            EXPECTED_RETAINED_ANALYSIS_OUTPUTS
        ),
        "expected_durable_root_artifact_count": (
            EXPECTED_DURABLE_ROOT_ARTIFACTS
        ),
        "expected_source_occurrence_count": EXPECTED_JOB_COUNT,
        "source_occurrences_per_output_pair": 1,
    }
    if (
        partition_contract.get("schema") != resolver.PARTITION_SCHEMA
        or any(
            partition_contract.get(field) != expected
            for field, expected in expected_counts.items()
        )
        or partition_contract.get("capacity_canary_required_before_submission")
        is not True
        or partition_contract.get("capacity_authority_earned") is not False
    ):
        raise ControllerError("aggregate execution partition differs")
    artifact = require_mapping(
        partition_contract.get("partition_artifact"),
        "partition artifact",
    )
    if (
        artifact.get("name") != partition_path.name
        or artifact.get("record_count") != EXPECTED_JOB_COUNT
        or artifact.get("sha256") != file_sha256(partition_path)
    ):
        raise ControllerError("partition artifact binding differs")

    chunks: list[dict[str, Any]] = []
    by_row: dict[str, list[dict[str, Any]]] = {
        row_id: [] for row_id in rows_by_id
    }
    tuple_fingerprints: set[str] = set()
    chunk_fingerprints: set[str] = set()
    execution_fingerprints: set[str] = set()
    current_row_index = 0
    row_ids = list(rows_by_id)
    try:
        stream = partition_path.open("r", encoding="utf-8")
    except OSError as exc:
        raise ControllerError(f"cannot read partition artifact: {partition_path}") from exc
    with stream:
        for global_index, line in enumerate(stream):
            if not line.strip():
                raise ControllerError("partition contains an empty record")
            try:
                record = require_mapping(json.loads(line), "partition record")
            except json.JSONDecodeError as exc:
                raise ControllerError("partition contains invalid JSON") from exc
            if set(record) != set(resolver.EXECUTION_CHUNK_KEYS):
                raise ControllerError("partition record field inventory differs")
            row_id = str(record.get("row_id", ""))
            if row_id not in rows_by_id:
                raise ControllerError("partition references an unknown row")
            while current_row_index < len(row_ids) and row_ids[current_row_index] != row_id:
                if by_row[row_ids[current_row_index]]:
                    current_row_index += 1
                else:
                    break
            if current_row_index >= len(row_ids) or row_ids[current_row_index] != row_id:
                raise ControllerError("partition row order differs")
            row = rows_by_id[row_id]
            row_chunks = by_row[row_id]
            chunk_index = len(row_chunks)
            start = sum(int(item["tuple_count"]) for item in row_chunks)
            count = record.get("tuple_count")
            if (
                record.get("schema") != resolver.CHUNK_SCHEMA
                or record.get("group_size") != EXPECTED_GROUP_SIZE
                or record.get("chunk_index") != chunk_index
                or record.get("segment") != chunk_index + 1
                or record.get("global_chunk_index") != global_index
                or record.get("tuple_index_start") != start
                or not isinstance(count, int)
                or isinstance(count, bool)
                or count <= 0
                or count > EXPECTED_GROUP_SIZE
                or record.get("tuple_index_end_exclusive") != start + count
                or record.get("full_source_manifest_sha256")
                != row["input_contract"]["full_source_manifest_sha256"]
            ):
                raise ControllerError(f"{row_id} partition membership differs")
            physical_lines = require_sequence(
                record.get("physical_lines"), f"{row_id} physical lines"
            )
            tuple_hashes = require_sequence(
                record.get("tuple_input_sha256s"), f"{row_id} tuple hashes"
            )
            if (
                len(physical_lines) != count
                or len(tuple_hashes) != count
                or not all(
                    isinstance(value, int) and not isinstance(value, bool) and value > 0
                    for value in physical_lines
                )
                or physical_lines != sorted(physical_lines)
            ):
                raise ControllerError(f"{row_id} source partition members differ")
            for value in tuple_hashes:
                tuple_sha = require_sha256(value, f"{row_id} tuple SHA-256")
                if tuple_sha in tuple_fingerprints:
                    raise ControllerError("duplicate source tuple fingerprints are forbidden")
                tuple_fingerprints.add(tuple_sha)
            require_sha256(
                record.get("tuple_records_sha256"),
                f"{row_id} tuple-record SHA-256",
            )
            chunk_fingerprint = require_sha256(
                record.get("chunk_fingerprint_sha256"),
                f"{row_id} chunk fingerprint",
            )
            chunk_payload = {
                key: value
                for key, value in record.items()
                if key
                not in {"global_chunk_index", "execution_chunk_sha256"}
            }
            if chunk_fingerprint != canonical_sha256(
                {
                    key: value
                    for key, value in chunk_payload.items()
                    if key != "chunk_fingerprint_sha256"
                }
            ):
                raise ControllerError(f"{row_id} chunk fingerprint differs")
            execution_fingerprint = require_sha256(
                record.get("execution_chunk_sha256"),
                f"{row_id} execution fingerprint",
            )
            if execution_fingerprint != canonical_sha256(
                {
                    key: value
                    for key, value in record.items()
                    if key != "execution_chunk_sha256"
                }
            ):
                raise ControllerError(f"{row_id} execution fingerprint differs")
            if (
                chunk_fingerprint in chunk_fingerprints
                or execution_fingerprint in execution_fingerprints
            ):
                raise ControllerError("duplicate partition fingerprints are forbidden")
            chunk_fingerprints.add(chunk_fingerprint)
            execution_fingerprints.add(execution_fingerprint)
            row_chunks.append(record)
            chunks.append(record)
    if len(chunks) != EXPECTED_JOB_COUNT:
        raise ControllerError(
            f"partition job count={len(chunks)}, expected={EXPECTED_JOB_COUNT}"
        )
    if len(tuple_fingerprints) != EXPECTED_SOURCE_TUPLES:
        raise ControllerError(
            f"partition source occurrence count={len(tuple_fingerprints)}, "
            f"expected={EXPECTED_SOURCE_TUPLES}"
        )
    for row_id, row in rows_by_id.items():
        row_chunks = by_row[row_id]
        expected_jobs = row["_expected_job_count"]
        expected_tuples = row["_expected_tuple_count"]
        if (
            len(row_chunks) != expected_jobs
            or sum(record["tuple_count"] for record in row_chunks)
            != expected_tuples
            or any(
                record["tuple_count"] != EXPECTED_GROUP_SIZE
                for record in row_chunks[:-1]
            )
            or canonical_sha256(
                [
                    {
                        key: value
                        for key, value in record.items()
                        if key not in {"global_chunk_index", "execution_chunk_sha256"}
                    }
                    for record in row_chunks
                ]
            )
            != row["input_contract"]["chunk_records_sha256"]
        ):
            raise ControllerError(f"{row_id} partition closure differs")
    return chunks, by_row


def validate_storage_certificate(
    payload: dict[str, Any],
    *,
    certificate_artifact: Mapping[str, Any],
    controller_budget_artifact: Mapping[str, Any],
    plan_artifact: Mapping[str, Any],
    plan: Mapping[str, Any],
    bundle_sha256: str,
    materialization_sha256: str,
) -> dict[str, Any]:
    if (
        payload.get("schema") != projector.STORAGE_CERTIFICATE_SCHEMA
        or payload.get("status") != projector.PASS_STATUS
        or payload.get("gate") != projector.GATE
        or payload.get("blocker_codes") != []
    ):
        raise ControllerError("storage admission certificate is not exact PASS")
    try:
        projector.validate_authority(payload, "storage certificate")
        projector.validate_all_authority_occurrences(payload, "storage certificate")
        projector.verify_semantic_receipt(
            payload, "certificate_semantic_sha256", "storage certificate"
        )
    except projector.ProjectionError as exc:
        raise ControllerError(f"storage certificate authority/hash differs: {exc}") from exc
    manifest_ref = require_mapping(payload.get("manifest"), "storage manifest ref")
    snapshot_ref = require_mapping(
        payload.get("quota_snapshot"), "quota snapshot ref"
    )
    manifest_path = Path(str(manifest_ref.get("path", ""))).absolute()
    snapshot_path = Path(str(snapshot_ref.get("path", ""))).absolute()
    if (
        file_sha256(manifest_path) != manifest_ref.get("sha256")
        or file_sha256(snapshot_path) != snapshot_ref.get("sha256")
    ):
        raise ControllerError("storage certificate evidence hash differs")
    try:
        rebuilt = projector.build_projection(manifest_path, snapshot_path)
    except projector.ProjectionError as exc:
        raise ControllerError(f"storage certificate readback failed: {exc}") from exc
    if rebuilt != payload:
        raise ControllerError("storage certificate live readback/profile differs")
    try:
        manifest_payload = projector.validate_storage_manifest(
            require_mapping(
                json.loads(manifest_path.read_text(encoding="utf-8")),
                "storage manifest",
            )
        )
    except (OSError, UnicodeDecodeError, json.JSONDecodeError, projector.ProjectionError) as exc:
        raise ControllerError(f"storage manifest validation failed: {exc}") from exc
    measurement = require_mapping(
        manifest_payload.get("artifact_measurement"), "artifact measurement"
    )
    bindings = require_mapping(measurement.get("bindings"), "storage bindings")
    for field, expected in (
        ("plan", plan_artifact["sha256"]),
        ("bundle_manifest", bundle_sha256),
        ("materialization_receipt", materialization_sha256),
    ):
        record = require_mapping(bindings.get(field), f"storage binding {field}")
        if record.get("sha256") != expected:
            raise ControllerError(f"storage binding {field} differs")
    controller = require_mapping(
        measurement.get("controller_dry_materialization_budget"),
        "controller dry-materialization budget",
    )
    receipt = require_mapping(controller.get("receipt"), "controller budget receipt")
    if receipt.get("sha256") != controller_budget_artifact["sha256"]:
        raise ControllerError("storage certificate binds a different controller budget")
    count_contract = require_mapping(
        measurement.get("count_contract"), "storage count contract"
    )
    exact_counts = {
        "row_count": EXPECTED_ROW_COUNT,
        "source_tuple_count": EXPECTED_SOURCE_TUPLES,
        "group_size": EXPECTED_GROUP_SIZE,
        "expected_job_count": EXPECTED_JOB_COUNT,
        "expected_output_pair_count": EXPECTED_OUTPUT_PAIRS,
        "expected_analysis_output_count": EXPECTED_ANALYSIS_OUTPUTS,
        "expected_sidecar_output_count": EXPECTED_SIDECAR_OUTPUTS,
        "expected_physical_root_artifact_count": EXPECTED_PHYSICAL_ARTIFACTS,
        "expected_retained_analysis_output_count": (
            EXPECTED_RETAINED_ANALYSIS_OUTPUTS
        ),
        "expected_durable_root_artifact_count": (
            EXPECTED_DURABLE_ROOT_ARTIFACTS
        ),
        "expected_source_occurrence_count": EXPECTED_JOB_COUNT,
        "request_memory_mb": EXPECTED_REQUEST_MEMORY_MB,
    }
    if count_contract != exact_counts:
        raise ControllerError("storage certificate count contract differs")
    roots = require_mapping(
        measurement.get("storage_domain_roots"), "storage domain roots"
    )
    campaign = require_mapping(plan.get("campaign"), "plan campaign")
    if (
        roots.get("bulk_science") != campaign.get("output_root")
        or roots.get("scratch_control") != campaign.get("submit_root")
        or roots.get("scratch_evidence") != campaign.get("evidence_root")
    ):
        raise ControllerError("storage certificate namespace binding differs")
    for projection in require_sequence(
        payload.get("quota_domain_projections"), "quota projections"
    ):
        record = require_mapping(projection, "quota projection")
        if (
            record.get("byte_headroom_at_least_20_percent") is not True
            or record.get("inode_headroom_at_least_20_percent") is not True
        ):
            raise ControllerError("storage certificate lacks required headroom")
    return controller


def validate_budget(
    payload: dict[str, Any],
    artifact: Mapping[str, Any],
) -> dict[str, int]:
    try:
        normalized, blockers = projector.validate_controller_budget(
            {"path": artifact["path"], "sha256": artifact["sha256"]}
        )
    except projector.ProjectionError as exc:
        raise ControllerError(f"controller budget certificate differs: {exc}") from exc
    if normalized is None or blockers:
        raise ControllerError("controller storage budget certificate is missing or blocked")
    validate_authority_occurrences(payload, "controller budget")
    return {
        "projected_bytes": int(normalized["projected_bytes"]),
        "projected_inodes": int(normalized["projected_inodes"]),
    }


def validate_plan_and_evidence(args: argparse.Namespace) -> dict[str, Any]:
    if args.artifact_profile != ARTIFACT_PROFILE:
        raise ControllerError(
            f"artifact profile must be exactly {ARTIFACT_PROFILE}"
        )
    staging_root = validate_staging_root(args.staging_root)
    plan, plan_artifact = load_pinned_json(
        args.plan, args.plan_sha256, "resolved extraction plan"
    )
    require_exact_keys(
        plan,
        {
            "schema",
            "status",
            "execution_state",
            "submission_performed",
            "campaign",
            "training_period_si_contract",
            "training_period_si_contract_sha256",
            "authority",
            "artifact_profile",
            "closure_witness_boundary",
            "source_family_closure",
            "execution_partition",
            "input_manifests",
            "duplicate_contract",
            "duplicate_fingerprint_sha256",
            "execution_fingerprint_sha256",
            "rows",
        },
        "resolved extraction plan",
    )
    if (
        plan.get("schema") != resolver.PLAN_SCHEMA
        or plan.get("status") != "PREFLIGHT_PASS"
        or plan.get("execution_state") != "PREFLIGHT_ONLY_NO_CONDOR_MUTATION"
        or plan.get("submission_performed") is not False
    ):
        raise ControllerError("resolved plan schema/status/state differs")
    if plan.get("artifact_profile") != RESOLVED_ARTIFACT_PROFILE:
        raise ControllerError("resolved sidecar-only artifact profile differs")
    validate_authority_occurrences(plan, "plan")
    try:
        resolver.validate_preflight_authority_payload(plan, label="plan")
    except resolver.ControllerError as exc:
        raise ControllerError(f"resolved plan authority differs: {exc}") from exc
    authority = require_mapping(plan.get("authority"), "plan authority")
    if authority != {
        "requested_scope": "full",
        "full_training_authority": 0,
        "authority_state": resolver.PREFLIGHT_AUTHORITY_STATE,
    }:
        raise ControllerError("resolved plan authority object differs")

    campaign_raw = require_mapping(plan.get("campaign"), "plan campaign")
    require_exact_keys(
        campaign_raw,
        {"tag", "output_root", "evidence_root", "submit_root"},
        "plan campaign",
    )
    tag = campaign_raw.get("tag")
    if not isinstance(tag, str) or resolver.TAG_RE.fullmatch(tag) is None:
        raise ControllerError("campaign tag differs")
    campaign = {
        "tag": tag,
        "output_root": validate_remote_root(
            campaign_raw.get("output_root"), "campaign.output_root"
        ),
        "evidence_root": validate_remote_root(
            campaign_raw.get("evidence_root"), "campaign.evidence_root"
        ),
        "submit_root": validate_remote_root(
            campaign_raw.get("submit_root"), "campaign.submit_root"
        ),
    }
    if len(set(campaign.values())) != 4:
        raise ControllerError("campaign tag/namespaces must be distinct")
    for field in ("output_root", "evidence_root", "submit_root"):
        if PurePosixPath(campaign[field]).name != tag:
            raise ControllerError(f"campaign {field} basename differs from tag")
    for field in ("output_root", "submit_root"):
        if os.path.lexists(campaign[field]):
            raise ControllerError(f"campaign {field} namespace already exists")
    for field in ("output_root", "evidence_root", "submit_root"):
        if paths_overlap(staging_root, Path(campaign[field])):
            raise ControllerError(
                "local staging root overlaps declared campaign namespace: "
                f"{field}"
            )

    receipt, receipt_artifact = load_pinned_json(
        args.preflight_receipt,
        args.preflight_receipt_sha256,
        "resolver preflight receipt",
    )
    if (
        receipt.get("schema") != resolver.RECEIPT_SCHEMA
        or receipt.get("status") != "PASS"
        or receipt.get("submission_performed") is not False
        or receipt.get("row_count") != EXPECTED_ROW_COUNT
        or receipt.get("requested_scope") != "full"
        or receipt.get("full_training_authority") != 0
        or receipt.get("authority_state") != resolver.PREFLIGHT_AUTHORITY_STATE
        or receipt.get("artifact_profile_sha256")
        != RESOLVED_ARTIFACT_PROFILE_SHA256
    ):
        raise ControllerError("resolver preflight receipt differs")
    validate_authority_occurrences(receipt, "preflight receipt")
    receipt_artifacts = require_mapping(
        receipt.get("artifacts"), "preflight receipt artifacts"
    )
    if require_mapping(
        receipt_artifacts.get("plan"), "preflight plan binding"
    ).get("sha256") != plan_artifact["sha256"]:
        raise ControllerError("preflight receipt plan SHA-256 differs")

    inputs = require_mapping(plan.get("input_manifests"), "plan input manifests")
    bundle_ref = require_mapping(inputs.get("bundle"), "plan bundle ref")
    materialization_ref = require_mapping(
        inputs.get("materialization"), "plan materialization ref"
    )
    source_ref = require_mapping(inputs.get("sources"), "plan source ref")
    bundle_path = Path(str(bundle_ref.get("path", ""))).absolute()
    materialization_path = Path(str(materialization_ref.get("path", ""))).absolute()
    source_path = Path(str(source_ref.get("path", ""))).absolute()
    bundle_payload, bundle_artifact = load_pinned_json(
        bundle_path,
        str(bundle_ref.get("sha256", "")),
        "immutable resolver bundle",
    )
    materialization_payload, materialization_artifact = load_pinned_json(
        materialization_path,
        str(materialization_ref.get("sha256", "")),
        "immutable materialization receipt",
    )
    source_payload, source_artifact = load_pinned_json(
        source_path,
        str(source_ref.get("sha256", "")),
        "source authority manifest",
    )
    bundle, by_role = validate_bundle(bundle_payload)
    validate_materialization_binding(
        materialization_payload,
        materialization_path=materialization_path,
        bundle_path=bundle_path,
        bundle_sha256=bundle_artifact["sha256"],
        bundle=bundle,
        plan_record=materialization_ref,
    )
    if (
        source_payload.get("schema") != resolver.SOURCE_SCHEMA
        or source_payload.get("status") != "PASS"
    ):
        raise ControllerError("source authority schema/status differs")
    training_contract = require_mapping(
        plan.get("training_period_si_contract"), "training period"
    )
    training_period = str(training_contract.get("period", ""))
    try:
        source_authority = resolver.validate_source_period_authority(
            source_payload,
            pp_period=training_period,
        )
    except resolver.ControllerError as exc:
        raise ControllerError(str(exc)) from exc
    if source_authority != {
        "pp_period": training_period,
        "pp_si_di_role": training_contract.get("si_di_role"),
    }:
        raise ControllerError("source period/SI authority differs")
    if (
        receipt.get("bundle_manifest_sha256") != bundle_artifact["sha256"]
        or receipt.get("materialization_receipt_sha256")
        != materialization_artifact["sha256"]
        or receipt.get("source_manifest_sha256") != source_artifact["sha256"]
    ):
        raise ControllerError("preflight bundle/materialization/source binding differs")

    rows, rows_by_id = validate_rows(
        plan,
        campaign=campaign,
        bundle=bundle,
        by_role=by_role,
        source_manifest=source_payload,
    )
    partition_ref = require_mapping(
        receipt_artifacts.get("partition"), "preflight partition binding"
    )
    partition_name = str(partition_ref.get("name", ""))
    if (
        not partition_name
        or Path(partition_name).name != partition_name
        or partition_name != plan["execution_partition"]["partition_artifact"]["name"]
    ):
        raise ControllerError("partition artifact name is unsafe or differs")
    partition_path = args.plan.absolute().parent / partition_name
    partition_sha = file_sha256(partition_path)
    if (
        partition_ref.get("sha256") != partition_sha
        or partition_ref.get("record_count") != EXPECTED_JOB_COUNT
        or receipt.get("execution_partition_sha256")
        != canonical_sha256(plan["execution_partition"])
    ):
        raise ControllerError("preflight partition receipt differs")
    chunks, _chunks_by_row = validate_partition(
        partition_path, plan=plan, rows_by_id=rows_by_id
    )

    duplicate = require_mapping(
        plan.get("duplicate_contract"), "duplicate contract"
    )
    duplicate_fingerprint = require_sha256(
        plan.get("duplicate_fingerprint_sha256"), "duplicate fingerprint"
    )
    if (
        duplicate.get("schema") != resolver.DUPLICATE_SCHEMA
        or duplicate.get("requested_scope") != "full"
        or duplicate.get("full_training_authority") != 0
        or duplicate.get("authority_state") != resolver.PREFLIGHT_AUTHORITY_STATE
        or duplicate_fingerprint != canonical_sha256(duplicate)
        or receipt.get("duplicate_fingerprint_sha256") != duplicate_fingerprint
    ):
        raise ControllerError("duplicate fingerprint contract differs")
    duplicate_sources = require_sequence(
        duplicate.get("sources"), "duplicate contract sources"
    )
    duplicate_ids = [
        record.get("row_id") if isinstance(record, dict) else None
        for record in duplicate_sources
    ]
    if duplicate_ids != list(rows_by_id) or len(set(duplicate_ids)) != EXPECTED_ROW_COUNT:
        raise ControllerError("duplicate contract rows are missing or duplicated")

    expected_execution_fingerprint = canonical_sha256(
        {
            "schema": resolver.EXECUTION_SCHEMA,
            "tag": campaign["tag"],
            "output_root": campaign["output_root"],
            "evidence_root": campaign["evidence_root"],
            "submit_root": campaign["submit_root"],
            "materialization_receipt_sha256": materialization_artifact["sha256"],
            "bundle_manifest_sha256": bundle_artifact["sha256"],
            "source_manifest_sha256": source_artifact["sha256"],
            "duplicate_fingerprint_sha256": duplicate_fingerprint,
            "partition_artifact_sha256": partition_sha,
            "execution_partition_sha256": canonical_sha256(
                plan["execution_partition"]
            ),
            "row_fingerprints": [
                row["row_fingerprint_sha256"] for row in rows
            ],
        }
    )
    if (
        plan.get("execution_fingerprint_sha256")
        != expected_execution_fingerprint
        or receipt.get("execution_fingerprint_sha256")
        != expected_execution_fingerprint
    ):
        raise ControllerError("execution fingerprint differs")

    budget_payload, budget_artifact = load_pinned_json(
        args.controller_budget,
        args.controller_budget_sha256,
        "controller storage budget",
    )
    budget = validate_budget(budget_payload, budget_artifact)
    storage_payload, storage_artifact = load_pinned_json(
        args.storage_certificate,
        args.storage_certificate_sha256,
        "pre-extraction storage certificate",
    )
    storage_controller = validate_storage_certificate(
        storage_payload,
        certificate_artifact=storage_artifact,
        controller_budget_artifact=budget_artifact,
        plan_artifact=plan_artifact,
        plan=plan,
        bundle_sha256=bundle_artifact["sha256"],
        materialization_sha256=materialization_artifact["sha256"],
    )
    if (
        storage_controller.get("projected_bytes") != budget["projected_bytes"]
        or storage_controller.get("projected_inodes") != budget["projected_inodes"]
    ):
        raise ControllerError("storage certificate/controller budget projection differs")
    return {
        "staging_root": staging_root,
        "plan": plan,
        "plan_artifact": plan_artifact,
        "preflight_receipt_artifact": receipt_artifact,
        "storage_certificate_artifact": storage_artifact,
        "controller_budget_artifact": budget_artifact,
        "bundle_artifact": bundle_artifact,
        "materialization_artifact": materialization_artifact,
        "source_artifact": source_artifact,
        "partition_artifact": {
            "path": str(partition_path),
            "sha256": partition_sha,
            "size_bytes": partition_path.stat().st_size,
        },
        "campaign": campaign,
        "rows": rows,
        "rows_by_id": rows_by_id,
        "chunks": chunks,
        "budget": budget,
        "duplicate_fingerprint_sha256": duplicate_fingerprint,
        "execution_fingerprint_sha256": expected_execution_fingerprint,
    }


def staged_row(row: Mapping[str, Any]) -> dict[str, Any]:
    public_row = {
        key: value for key, value in row.items() if not key.startswith("_")
    }
    return {
        "schema": ROW_STAGE_SCHEMA,
        "artifact_profile": ARTIFACT_PROFILE,
        "authority": dict(CONTROLLER_AUTHORITY),
        "row_descriptor": public_row,
        "row_fingerprint_sha256": row["row_fingerprint_sha256"],
        "worker_environment_sha256": row["_worker_environment_sha256"],
        "materialization_environment_sha256": row[
            "_materialization_environment_sha256"
        ],
        "provenance_placeholders": dict(PROVENANCE_PLACEHOLDERS),
    }


def staged_job(
    chunk: Mapping[str, Any], row: Mapping[str, Any]
) -> dict[str, Any]:
    execution = require_mapping(row["execution_contract"], "row execution")
    return {
        "schema": JOB_STAGE_SCHEMA,
        "artifact_profile": ARTIFACT_PROFILE,
        "authority": dict(CONTROLLER_AUTHORITY),
        "global_job_index": chunk["global_chunk_index"],
        "row_id": row["row_id"],
        "row_job_index": chunk["chunk_index"],
        "request_memory_mb": EXPECTED_REQUEST_MEMORY_MB,
        "row_fingerprint_sha256": row["row_fingerprint_sha256"],
        "worker_environment_sha256": row["_worker_environment_sha256"],
        "materialization_environment_sha256": row[
            "_materialization_environment_sha256"
        ],
        "source_partition_membership": {
            "full_source_manifest_sha256": chunk[
                "full_source_manifest_sha256"
            ],
            "row_partition_sha256": row["input_contract"][
                "row_partition_sha256"
            ],
            "execution_chunk_sha256": chunk["execution_chunk_sha256"],
            "chunk_fingerprint_sha256": chunk[
                "chunk_fingerprint_sha256"
            ],
            "tuple_index_start": chunk["tuple_index_start"],
            "tuple_index_end_exclusive": chunk[
                "tuple_index_end_exclusive"
            ],
            "tuple_count": chunk["tuple_count"],
            "physical_lines": chunk["physical_lines"],
            "tuple_input_sha256s": chunk["tuple_input_sha256s"],
        },
        "expected_outputs": {
            "output_pair_count": 1,
            "physical_root_artifact_count": 2,
            "analysis_root": {
                "count": 1,
                "namespace": execution["analysis_output_namespace"],
                "authority": "INCIDENTAL_NON_PHYSICS_EXISTING_EXECUTOR_OUTPUT",
                "physics_output_authority": False,
            },
            "training_sidecar_root": {
                "count": 1,
                "path_template": execution["multiview_sidecar_template"],
                "artifact_profile": ARTIFACT_PROFILE,
            },
        },
        "provenance_placeholders": dict(PROVENANCE_PLACEHOLDERS),
    }


def jsonl_bytes(records: Iterable[Mapping[str, Any]]) -> bytes:
    return b"".join(canonical_json_bytes(record) for record in records)


def submit_description_bytes(context: Mapping[str, Any]) -> bytes:
    lines = [
        f"schema={CONTROLLER_SCHEMA}",
        f"authority_state={AUTHORITY_STATE}",
        f"artifact_profile={ARTIFACT_PROFILE}",
        "submission_performed=false",
        "submission_authority=false",
        "condor_submit_argv=FORBIDDEN_AND_NOT_MATERIALIZED",
        "job_control=FORBIDDEN",
        "network_access=FORBIDDEN",
        "remote_namespace_creation=FORBIDDEN",
        f"row_count={EXPECTED_ROW_COUNT}",
        f"source_tuple_occurrence_count={EXPECTED_SOURCE_TUPLES}",
        f"group_size={EXPECTED_GROUP_SIZE}",
        f"job_count={EXPECTED_JOB_COUNT}",
        f"output_pair_count={EXPECTED_OUTPUT_PAIRS}",
        f"physical_root_artifact_count={EXPECTED_PHYSICAL_ARTIFACTS}",
        f"request_memory_mb={EXPECTED_REQUEST_MEMORY_MB}",
        f"plan_sha256={context['plan_artifact']['sha256']}",
        f"partition_sha256={context['partition_artifact']['sha256']}",
        (
            "duplicate_fingerprint_sha256="
            f"{context['duplicate_fingerprint_sha256']}"
        ),
        (
            "execution_fingerprint_sha256="
            f"{context['execution_fingerprint_sha256']}"
        ),
        "row_manifest=the134_full_multiview_rows.jsonl",
        "job_manifest=the134_full_multiview_jobs.jsonl",
        (
            "provenance_RJ_CODEX_CHAT_NAME="
            f"{PROVENANCE_PLACEHOLDERS['RJ_CODEX_CHAT_NAME']}"
        ),
        (
            "provenance_RJ_CODEX_THREAD_ID="
            f"{PROVENANCE_PLACEHOLDERS['RJ_CODEX_THREAD_ID']}"
        ),
        "note=INERT_DESCRIPTION_ONLY_NOT_A_CONDOR_SUBMIT_DESCRIPTION",
    ]
    return ("\n".join(lines) + "\n").encode("utf-8")


def build_staged_artifacts(context: Mapping[str, Any]) -> dict[str, bytes]:
    row_bytes = jsonl_bytes(staged_row(row) for row in context["rows"])
    job_bytes = jsonl_bytes(
        staged_job(chunk, context["rows_by_id"][chunk["row_id"]])
        for chunk in context["chunks"]
    )
    description_bytes = submit_description_bytes(context)
    artifacts = {
        OUTPUT_FILENAMES[0]: row_bytes,
        OUTPUT_FILENAMES[1]: job_bytes,
        OUTPUT_FILENAMES[2]: description_bytes,
    }
    artifact_inventory = {
        name: {
            "sha256": hashlib.sha256(data).hexdigest(),
            "size_bytes": len(data),
        }
        for name, data in artifacts.items()
    }
    manifest = {
        "schema": MANIFEST_SCHEMA,
        "status": "PASS_DRY_MATERIALIZED_LOCAL_ONLY",
        "artifact_profile": ARTIFACT_PROFILE,
        "authority": dict(CONTROLLER_AUTHORITY),
        "counts": {
            "logical_source_row_count": EXPECTED_ROW_COUNT,
            "source_tuple_occurrence_count": EXPECTED_SOURCE_TUPLES,
            "group_size": EXPECTED_GROUP_SIZE,
            "job_count": EXPECTED_JOB_COUNT,
            "output_pair_count": EXPECTED_OUTPUT_PAIRS,
            "analysis_root_count": EXPECTED_ANALYSIS_OUTPUTS,
            "training_sidecar_root_count": EXPECTED_SIDECAR_OUTPUTS,
            "physical_root_artifact_count": EXPECTED_PHYSICAL_ARTIFACTS,
            "request_memory_mb": EXPECTED_REQUEST_MEMORY_MB,
        },
        "bindings": {
            "plan": context["plan_artifact"],
            "resolver_preflight_receipt": context[
                "preflight_receipt_artifact"
            ],
            "storage_certificate": context["storage_certificate_artifact"],
            "controller_budget": context["controller_budget_artifact"],
            "immutable_bundle": context["bundle_artifact"],
            "immutable_materialization": context[
                "materialization_artifact"
            ],
            "source_authority": context["source_artifact"],
            "source_partition": context["partition_artifact"],
            "duplicate_fingerprint_sha256": context[
                "duplicate_fingerprint_sha256"
            ],
            "execution_fingerprint_sha256": context[
                "execution_fingerprint_sha256"
            ],
            "resolved_artifact_profile_sha256": (
                RESOLVED_ARTIFACT_PROFILE_SHA256
            ),
        },
        "staged_artifacts": artifact_inventory,
        "provenance_placeholders": dict(PROVENANCE_PLACEHOLDERS),
        "controller_storage": {
            "actual_bytes": 0,
            "actual_inodes": len(OUTPUT_FILENAMES),
            "certified_bytes": context["budget"]["projected_bytes"],
            "certified_inodes": context["budget"]["projected_inodes"],
        },
        "boundaries": [
            "Local dry materialization only.",
            "No Condor submission or job control.",
            "No remote/output/submit namespace creation.",
            "No THE-121, THE-122, physics-output, replay-cache, science-freeze, or CANONICAL authority.",
        ],
    }
    other_size = sum(len(data) for data in artifacts.values())
    prior = -1
    for _ in range(16):
        manifest_bytes = canonical_json_bytes(manifest)
        total = other_size + len(manifest_bytes)
        if total == prior and manifest["controller_storage"]["actual_bytes"] == total:
            break
        manifest["controller_storage"]["actual_bytes"] = total
        prior = total
    else:
        raise ControllerError("controller manifest byte-size fixed point did not converge")
    manifest_bytes = canonical_json_bytes(manifest)
    final_total = other_size + len(manifest_bytes)
    if manifest["controller_storage"]["actual_bytes"] != final_total:
        raise ControllerError("controller manifest byte accounting differs")
    artifacts[OUTPUT_FILENAMES[3]] = manifest_bytes
    return artifacts


def enforce_controller_budget(
    artifacts: Mapping[str, bytes], budget: Mapping[str, int]
) -> tuple[int, int]:
    actual_bytes = sum(len(data) for data in artifacts.values())
    actual_inodes = len(artifacts)
    if actual_bytes > budget["projected_bytes"]:
        raise ControllerError(
            "controller materialization exceeds certified bytes: "
            f"required={actual_bytes} certified={budget['projected_bytes']}"
        )
    if actual_inodes > budget["projected_inodes"]:
        raise ControllerError(
            "controller materialization exceeds certified inodes: "
            f"required={actual_inodes} certified={budget['projected_inodes']}"
        )
    return actual_bytes, actual_inodes


def materialize_local(
    staging_root: Path, artifacts: Mapping[str, bytes]
) -> None:
    try:
        staging_root.mkdir(mode=0o700, parents=False, exist_ok=False)
    except OSError as exc:
        raise ControllerError(
            f"cannot create fresh local staging root: {staging_root}"
        ) from exc
    try:
        for name in OUTPUT_FILENAMES:
            destination = staging_root / name
            destination.write_bytes(artifacts[name])
        observed = {
            path.name: file_sha256(path)
            for path in sorted(staging_root.iterdir())
            if path.is_file()
        }
        expected = {
            name: hashlib.sha256(data).hexdigest()
            for name, data in artifacts.items()
        }
        if observed != expected:
            raise ControllerError("local staged artifact readback differs")
    except Exception:
        for name in OUTPUT_FILENAMES:
            candidate = staging_root / name
            try:
                if candidate.is_file() and not candidate.is_symlink():
                    candidate.unlink()
            except OSError:
                pass
        try:
            staging_root.rmdir()
        except OSError:
            pass
        raise


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="action", required=True)
    for action, help_text in (
        ("preflight", "validate all inputs and report dry-stage size"),
        (
            "materialize",
            "validate all inputs and write only a fresh local staging root",
        ),
    ):
        command = subparsers.add_parser(action, help=help_text)
        command.add_argument("--plan", type=Path, required=True)
        command.add_argument("--plan-sha256", required=True)
        command.add_argument("--preflight-receipt", type=Path, required=True)
        command.add_argument("--preflight-receipt-sha256", required=True)
        command.add_argument("--storage-certificate", type=Path, required=True)
        command.add_argument("--storage-certificate-sha256", required=True)
        command.add_argument("--controller-budget", type=Path, required=True)
        command.add_argument("--controller-budget-sha256", required=True)
        command.add_argument(
            "--artifact-profile",
            required=True,
            choices=(ARTIFACT_PROFILE,),
        )
        command.add_argument("--staging-root", type=Path, required=True)
        command.add_argument(
            "--overwrite",
            action="store_true",
            help="forbidden; existing staging roots are never overwritten",
        )
    return parser.parse_args(argv)


def main(argv: Iterable[str] | None = None) -> int:
    try:
        args = parse_args(argv)
        if args.overwrite:
            raise ControllerError("--overwrite is explicitly forbidden")
        context = validate_plan_and_evidence(args)
        artifacts = build_staged_artifacts(context)
        actual_bytes, actual_inodes = enforce_controller_budget(
            artifacts, context["budget"]
        )
        if args.action == "materialize":
            materialize_local(context["staging_root"], artifacts)
        result = {
            "schema": CONTROLLER_SCHEMA,
            "status": (
                "PASS_DRY_MATERIALIZED_LOCAL_ONLY"
                if args.action == "materialize"
                else "PASS_PREFLIGHT_NO_WRITE"
            ),
            "artifact_profile": ARTIFACT_PROFILE,
            "authority": dict(CONTROLLER_AUTHORITY),
            "staging_root": str(context["staging_root"]),
            "staging_root_created": args.action == "materialize",
            "required_controller_bytes": actual_bytes,
            "required_controller_inodes": actual_inodes,
            "certified_controller_bytes": context["budget"]["projected_bytes"],
            "certified_controller_inodes": context["budget"]["projected_inodes"],
            "row_count": EXPECTED_ROW_COUNT,
            "source_tuple_occurrence_count": EXPECTED_SOURCE_TUPLES,
            "job_count": EXPECTED_JOB_COUNT,
            "physical_root_artifact_count": EXPECTED_PHYSICAL_ARTIFACTS,
            "submission_performed": False,
        }
        print(json.dumps(result, sort_keys=True))
        return 0
    except (ControllerError, FileNotFoundError) as exc:
        print(f"[THE134-FULL-MATERIALIZE][ERROR] {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
