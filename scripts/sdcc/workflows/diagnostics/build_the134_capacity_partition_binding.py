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
import re
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
FULL_GIT_COMMIT_PATTERN = re.compile(r"[0-9a-f]{40}")

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
CAPACITY_OPERATIONAL_BUNDLE_ROLES = frozenset(
    {"code_manifest", "submitter"}
)

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
        "ephemeral_analysis_health_receipts_revalidated",
        "durable_sidecar_hashes_revalidated",
        "legacy_evidence_consumed",
        "no_science_change",
        "no_tolerance_change",
        "no_submission_or_job_control",
    }
)
CURRENT_ROOT_ROW_KEYS = (
    evidence.ROOT_ROW_KEYS
    - frozenset({"analysis_sha256"})
    | frozenset(
        {
            "analysis_artifact_state",
            "analysis_health_receipt",
            "analysis_key_inventory_sha256",
        }
    )
)
ANALYSIS_HEALTH_RECEIPT_KEYS = frozenset({"path", "payload", "sha256"})
ANALYSIS_HEALTH_PAYLOAD_KEYS = frozenset(
    {
        "schema",
        "status",
        "mode",
        "analysis_size_bytes",
        "analysis_minimum_bytes",
        "analysis_root_non_zombie",
        "analysis_root_non_recovered",
        "analysis_directory_present",
        "analysis_histogram_present",
        "analysis_config_present",
        "analysis_key_inventory_sha256",
        "analysis_root_retained",
        "sidecar_size_bytes",
        "sidecar_tree_entries",
    }
)


def validate_science_validation_commit_binding(
    immutable: Mapping[str, Any],
    validation_authority: Mapping[str, Any],
) -> dict[str, str]:
    """Keep frozen science and validation-only commits explicit and typed."""

    science_commit = immutable.get("public_commit")
    controller = require_mapping(
        validation_authority.get("controller"),
        "capacity resource.validation_authority.controller",
    )
    validation_commit = controller.get("validation_commit")
    for label, value in (
        ("immutable science public_commit", science_commit),
        ("capacity validation controller validation_commit", validation_commit),
    ):
        if not isinstance(value, str) or not FULL_GIT_COMMIT_PATTERN.fullmatch(
            value
        ):
            raise BindingError(f"{label} must be a full 40-character Git commit")
    return {
        "science_commit": science_commit,
        "validation_commit": validation_commit,
        "authority_mode": (
            "SHARED_SCIENCE_AND_VALIDATION_COMMIT"
            if science_commit == validation_commit
            else "FROZEN_SCIENCE_WITH_VALIDATION_ONLY_COMMIT"
        ),
    }


CURRENT_VALIDATION_AUTHORITY_KEYS = frozenset(
    {
        "controller",
        "full_extraction_plan",
        "runtime_authority",
        "submission_journal",
        "submission_manifest",
        "submission_receipt",
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


def _validated_plan_bundle_operational_authority(
    payload: Mapping[str, Any],
) -> dict[str, Any] | None:
    """Validate one plan bundle and separate operational from science roles.

    Capacity measurements remain reusable across a controller/submitter-only
    bundle refresh only when both immutable receipts rehash and every
    non-operational artifact role, hash, and size remains identical.  The
    resulting science fingerprint is retained in the normalized plan, so a
    library, executor, configuration, model, schema, or provider change cannot
    be hidden by bundle-identity normalization.
    """

    raw_manifests = payload.get("input_manifests")
    if raw_manifests is None:
        return None
    manifests = require_mapping(raw_manifests, "capacity plan input_manifests")
    bundle_record = require_mapping(
        manifests.get("bundle"),
        "capacity plan input_manifests.bundle",
    )
    materialization_record = require_mapping(
        manifests.get("materialization"),
        "capacity plan input_manifests.materialization",
    )
    try:
        bundle_path = Path(str(bundle_record.get("path", ""))).resolve(
            strict=True
        )
        materialization_input_path = Path(
            str(materialization_record.get("path", ""))
        )
        materialization_path = materialization_input_path.resolve(strict=True)
    except (FileNotFoundError, OSError, RuntimeError) as exc:
        raise BindingError("capacity plan immutable bundle is missing") from exc
    for label, path, record in (
        ("bundle", bundle_path, bundle_record),
        ("materialization", materialization_path, materialization_record),
    ):
        expected_sha = str(record.get("sha256", ""))
        if (
            not re.fullmatch(r"[0-9a-f]{64}", expected_sha)
            or not path.is_file()
            or file_sha256(path) != expected_sha
        ):
            raise BindingError(
                f"capacity plan immutable {label} receipt hash differs"
            )
    try:
        bundle_payload = require_mapping(
            evidence.strict_json_loads(
                bundle_path.read_text(encoding="utf-8")
            ),
            "capacity plan immutable bundle",
        )
        validated_bundle = resolver.validate_bundle(dict(bundle_payload))
        materialization_payload = require_mapping(
            evidence.strict_json_loads(
                materialization_path.read_text(encoding="utf-8")
            ),
            "capacity plan immutable materialization",
        )
        validated_materialization = resolver.validate_materialization_binding(
            dict(materialization_payload),
            materialization_path=materialization_input_path,
            bundle_path=bundle_path,
            bundle_file_sha256=file_sha256(bundle_path),
            bundle=validated_bundle,
        )
    except (
        OSError,
        UnicodeDecodeError,
        json.JSONDecodeError,
        resolver.ControllerError,
    ) as exc:
        raise BindingError(
            f"capacity plan immutable bundle validation failed: {exc}"
        ) from exc
    by_role = require_mapping(
        validated_bundle.get("artifact_by_role"),
        "capacity plan immutable artifact_by_role",
    )
    submitter = require_mapping(
        by_role.get("submitter"),
        "capacity plan immutable submitter",
    )
    science_inventory = [
        {
            "role": role,
            "sha256": str(record["sha256"]),
            "size_bytes": int(record["size_bytes"]),
        }
        for role, record in sorted(by_role.items())
        if role not in CAPACITY_OPERATIONAL_BUNDLE_ROLES
    ]
    science_fingerprint = semantic_sha256(
        {
            "public_commit": validated_bundle["public_commit"],
            "replay_schema_sha256": validated_bundle[
                "replay_schema_sha256"
            ],
            "training_schema_sha256": validated_bundle[
                "training_schema_sha256"
            ],
            "semantic_sha256": validated_bundle["semantic_sha256"],
            "runtime": {
                key: validated_bundle["runtime"][key]
                for key in (
                    "release",
                    "offline_main",
                    "calo_reco_soname",
                    "request_memory_mb",
                )
            },
            "science_artifacts": science_inventory,
        }
    )
    bundle_root = Path(
        str(validated_materialization["digest_named_bundle_path"])
    )
    return {
        "bundle_roots": sorted(
            {
                str(bundle_root),
                str(bundle_root.resolve(strict=True)),
            },
            key=len,
            reverse=True,
        ),
        "bundle_identity_sha256": validated_bundle[
            "bundle_identity_sha256"
        ],
        "bundle_file_sha256": file_sha256(bundle_path),
        "bundle_semantic_fingerprint_sha256": validated_bundle[
            "semantic_fingerprint_sha256"
        ],
        "materialization_file_sha256": file_sha256(materialization_path),
        "code_sha256": validated_bundle["code_sha256"],
        "submitter_sha256": submitter["sha256"],
        "submitter_size_bytes": int(submitter["size_bytes"]),
        "science_bundle_fingerprint_sha256": science_fingerprint,
    }


def _capacity_plan_namespace_normal_form(
    payload: Mapping[str, Any],
) -> dict[str, Any]:
    """Remove campaign identity and byte-equivalent provider locations.

    A failed zero-cluster submission must use a fresh tag and namespaces.  The
    already-certified capacity pair remains applicable only when the new plan
    is otherwise identical.  A versioned release provider and its immutable
    bundle copy may differ lexically, but only when the exact file at either
    location rehashes to the unchanged pinned SHA-256.  Every source,
    partition, executable, physics, resource, ordering, worker field, and
    provider hash remains authoritative.
    """

    campaign = require_mapping(payload.get("campaign"), "capacity plan campaign")
    tag = campaign.get("tag")
    if not isinstance(tag, str) or not tag:
        raise BindingError("capacity plan campaign tag must be nonempty")

    bundle_authority = _validated_plan_bundle_operational_authority(payload)
    replacements: list[tuple[str, str]] = []
    if bundle_authority is not None:
        replacements.extend(
            (root, "__THE134_IMMUTABLE_BUNDLE__")
            for root in bundle_authority["bundle_roots"]
        )
        replacements.extend(
            (
                (bundle_authority["bundle_identity_sha256"],
                 "__THE134_BUNDLE_IDENTITY__"),
                (bundle_authority["bundle_file_sha256"],
                 "__THE134_BUNDLE_RECEIPT_SHA__"),
                (bundle_authority["bundle_semantic_fingerprint_sha256"],
                 "__THE134_BUNDLE_SEMANTIC_FINGERPRINT__"),
                (bundle_authority["materialization_file_sha256"],
                 "__THE134_MATERIALIZATION_RECEIPT_SHA__"),
                (bundle_authority["code_sha256"],
                 "__THE134_OPERATIONAL_CODE_SHA__"),
                (bundle_authority["submitter_sha256"],
                 "__THE134_OPERATIONAL_SUBMITTER_SHA__"),
            )
        )

    provider_fields = (
        (
            "RJ_PINNED_RELEASE_CALO_IO_PATH",
            "RJ_PINNED_RELEASE_CALO_IO_SHA256",
            "libcalo_io.so",
            "release_calo_io",
        ),
        (
            "RJ_PINNED_RELEASE_CLUSTERISO_PATH",
            "RJ_PINNED_RELEASE_CLUSTERISO_SHA256",
            "libclusteriso.so",
            "release_clusteriso",
        ),
        (
            "RJ_PINNED_RELEASE_JETBASE_PATH",
            "RJ_PINNED_RELEASE_JETBASE_SHA256",
            "libjetbase.so",
            "release_jetbase",
        ),
    )
    provider_contract_keys = {
        key
        for path_key, sha_key, _family, _bundle_role in provider_fields
        for key in (path_key, sha_key)
    }

    def normalize(value: Any) -> Any:
        if isinstance(value, dict):
            result = {
                key: item if key in provider_contract_keys else normalize(item)
                for key, item in value.items()
                if key not in {
                    "duplicate_fingerprint_sha256",
                    "execution_fingerprint_sha256",
                    "row_fingerprint_sha256",
                }
            }
            if (
                bundle_authority is not None
                and result.get("sha256")
                == "__THE134_OPERATIONAL_SUBMITTER_SHA__"
                and "size_bytes" in result
            ):
                result["size_bytes"] = "__THE134_OPERATIONAL_SUBMITTER_SIZE__"
            return result
        if isinstance(value, list):
            return [normalize(item) for item in value]
        if isinstance(value, str):
            result = value.replace(tag, "__THE134_CAMPAIGN_TAG__")
            for old, new in replacements:
                result = result.replace(old, new)
            return result
        return value

    normalized = normalize(dict(payload))
    if not isinstance(normalized, dict):
        raise BindingError("capacity plan normal form must be an object")
    if bundle_authority is not None:
        normalized[
            "_capacity_science_bundle_fingerprint_sha256"
        ] = bundle_authority["science_bundle_fingerprint_sha256"]
    rows = normalized.get("rows")
    if not isinstance(rows, list):
        raise BindingError("capacity plan rows must be a list")

    def beneath(path: Path, root: Path) -> bool:
        try:
            path.relative_to(root)
            return True
        except ValueError:
            return False

    for row_index, raw_row in enumerate(rows):
        row = require_mapping(raw_row, f"capacity plan rows[{row_index}]")
        execution = require_mapping(
            row.get("execution_contract"),
            f"capacity plan rows[{row_index}].execution_contract",
        )
        environment = execution.get("materialization_environment")
        if environment is None:
            continue
        environment = require_mapping(
            environment,
            f"capacity plan rows[{row_index}].materialization_environment",
        )
        release_roots = []
        for key in ("RJ_RELEASE_CORE_LIB64_DIR", "RJ_RELEASE_CORE_LIB_DIR"):
            raw_root = environment.get(key)
            if not isinstance(raw_root, str) or not raw_root.startswith("/"):
                raise BindingError(
                    f"capacity plan provider root is invalid: {key}"
                )
            release_roots.append(Path(raw_root).resolve(strict=True))
        for path_key, sha_key, family, bundle_role in provider_fields:
            raw_path = environment.get(path_key)
            expected_sha = environment.get(sha_key)
            if (
                not isinstance(raw_path, str)
                or not raw_path.startswith("/")
                or not isinstance(expected_sha, str)
                or not re.fullmatch(r"[0-9a-f]{64}", expected_sha)
            ):
                raise BindingError(
                    f"capacity plan pinned provider contract is invalid: {path_key}"
                )
            provider = Path(raw_path)
            try:
                provider_real = provider.resolve(strict=True)
            except (FileNotFoundError, OSError, RuntimeError) as exc:
                raise BindingError(
                    f"capacity plan pinned provider is missing: {path_key}"
                ) from exc
            if not provider_real.is_file() or file_sha256(provider_real) != expected_sha:
                raise BindingError(
                    f"capacity plan pinned provider hash differs: {path_key}"
                )
            in_release = any(
                beneath(provider_real, root)
                for root in release_roots
            )
            immutable_pattern = re.compile(
                r"/immutable_bundles/the134_bundle_sha256_[0-9a-f]{64}/"
                rf"artifacts/{bundle_role}/{re.escape(family)}$"
            )
            in_bundle = immutable_pattern.search(str(provider_real)) is not None
            if not (in_release or in_bundle):
                raise BindingError(
                    f"capacity plan pinned provider escaped immutable authority: "
                    f"{path_key}"
                )
            environment[path_key] = f"__PINNED_RELEASE_PROVIDER__/{family}"
    return normalized


def validate_capacity_plan_binding(
    resource: Mapping[str, Any],
    current: Mapping[str, Any],
) -> dict[str, str]:
    """Require the exact plan or one bounded execution-equivalent fresh plan."""

    resource_path = Path(str(resource.get("full_plan", ""))).absolute()
    resource_sha256 = str(resource.get("full_plan_sha256", ""))
    current_path = Path(str(current["plan"]["path"])).absolute()
    current_sha256 = str(current["plan"]["sha256"])
    if (
        resource_path == current_path
        and resource_sha256 == current_sha256
    ):
        return {
            "mode": "EXACT_PLAN",
            "capacity_plan_sha256": resource_sha256,
            "current_plan_sha256": current_sha256,
            "namespace_normalized_sha256": semantic_sha256(
                _capacity_plan_namespace_normal_form(
                    require_mapping(
                        evidence.strict_json_loads(
                            current_path.read_text(encoding="utf-8")
                        ),
                        "current capacity plan",
                    )
                )
            ),
        }

    if (
        not resource_path.is_file()
        or file_sha256(resource_path) != resource_sha256
        or not current_path.is_file()
        or file_sha256(current_path) != current_sha256
    ):
        raise BindingError("capacity plan equivalence artifacts differ")
    try:
        resource_plan = require_mapping(
            evidence.strict_json_loads(
                resource_path.read_text(encoding="utf-8")
            ),
            "preserved capacity plan",
        )
        current_plan = require_mapping(
            evidence.strict_json_loads(
                current_path.read_text(encoding="utf-8")
            ),
            "current capacity plan",
        )
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise BindingError("capacity plan equivalence JSON is invalid") from exc
    if resource_plan.get("schema") != current_plan.get("schema"):
        raise BindingError("capacity plan equivalence schema differs")
    resource_normal = _capacity_plan_namespace_normal_form(resource_plan)
    current_normal = _capacity_plan_namespace_normal_form(current_plan)
    resource_normal_sha256 = semantic_sha256(resource_normal)
    current_normal_sha256 = semantic_sha256(current_normal)
    if (
        resource_normal != current_normal
        or resource_normal_sha256 != current_normal_sha256
    ):
        raise BindingError(
            "capacity plan differs beyond campaign namespace, fingerprints, "
            "and byte-identical pinned-provider routing"
        )
    return {
        "mode": "FRESH_NAMESPACE_AND_PINNED_PROVIDER_EQUIVALENT_PLAN",
        "capacity_plan_sha256": resource_sha256,
        "current_plan_sha256": current_sha256,
        "namespace_normalized_sha256": current_normal_sha256,
    }


def validate_capacity_preflight_binding(
    resource: Mapping[str, Any],
    current: Mapping[str, Any],
) -> dict[str, str]:
    """Verify the old/current receipts differ only by derived plan identities."""

    resource_plan_path = Path(str(resource.get("full_plan", ""))).absolute()
    resource_path = resource_plan_path.parent / "preflight_receipt.json"
    resource_sha256 = str(resource.get("preflight_receipt_sha256", ""))
    current_path = Path(str(current["preflight_receipt"]["path"])).absolute()
    current_sha256 = str(current["preflight_receipt"]["sha256"])
    if (
        resource_path == current_path
        and resource_sha256 == current_sha256
    ):
        return {
            "mode": "EXACT_PREFLIGHT",
            "capacity_preflight_sha256": resource_sha256,
            "current_preflight_sha256": current_sha256,
        }
    if (
        not resource_path.is_file()
        or file_sha256(resource_path) != resource_sha256
        or not current_path.is_file()
        or file_sha256(current_path) != current_sha256
    ):
        raise BindingError("capacity preflight equivalence artifacts differ")
    try:
        resource_receipt = require_mapping(
            evidence.strict_json_loads(
                resource_path.read_text(encoding="utf-8")
            ),
            "preserved capacity preflight receipt",
        )
        current_receipt = require_mapping(
            evidence.strict_json_loads(
                current_path.read_text(encoding="utf-8")
            ),
            "current capacity preflight receipt",
        )
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise BindingError(
            "capacity preflight equivalence JSON is invalid"
        ) from exc

    def normal_form(payload: Mapping[str, Any]) -> dict[str, Any]:
        normalized = json.loads(json.dumps(payload))
        artifacts = require_mapping(
            normalized.get("artifacts"), "capacity preflight artifacts"
        )
        for field in ("plan", "rows"):
            record = require_mapping(
                artifacts.get(field), f"capacity preflight artifacts.{field}"
            )
            record["sha256"] = "__NAMESPACE_DERIVED_SHA256__"
        duplicate_record = require_mapping(
            artifacts.get("duplicate_fingerprint"),
            "capacity preflight artifacts.duplicate_fingerprint",
        )
        duplicate_record["sha256"] = "__NAMESPACE_DERIVED_SHA256__"
        normalized["execution_fingerprint_sha256"] = (
            "__NAMESPACE_DERIVED_SHA256__"
        )
        normalized["duplicate_fingerprint_sha256"] = (
            "__NAMESPACE_DERIVED_SHA256__"
        )
        normalized["bundle_manifest_sha256"] = (
            "__OPERATIONAL_BUNDLE_DERIVED_SHA256__"
        )
        normalized["materialization_receipt_sha256"] = (
            "__OPERATIONAL_BUNDLE_DERIVED_SHA256__"
        )
        return normalized

    resource_normal = normal_form(resource_receipt)
    current_normal = normal_form(current_receipt)
    if resource_normal != current_normal:
        raise BindingError(
            "capacity preflight differs beyond namespace-derived identities"
        )
    return {
        "mode": "FRESH_NAMESPACE_ONLY_EQUIVALENT_PREFLIGHT",
        "capacity_preflight_sha256": resource_sha256,
        "current_preflight_sha256": current_sha256,
    }


def capacity_resource_bundle_binding_valid(
    resource: Mapping[str, Any],
    immutable: Mapping[str, Any],
    plan_binding: Mapping[str, Any],
) -> bool:
    """Apply direct receipt equality only when the capacity plan is exact."""

    mode = plan_binding.get("mode")
    if mode == "EXACT_PLAN":
        return (
            resource.get("bundle_manifest_sha256")
            == immutable["bundle_manifest"]["sha256"]
            and resource.get("materialization_receipt_sha256")
            == immutable["materialization_receipt"]["sha256"]
        )
    return mode == "FRESH_NAMESPACE_AND_PINNED_PROVIDER_EQUIVALENT_PLAN"


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
    plan_binding: Mapping[str, Any],
    preflight_binding: Mapping[str, Any],
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
    # V18 is intentionally sidecar-only: the analysis ROOT is validated in
    # Condor scratch and then discarded.  Revalidate that current contract
    # directly instead of coercing it into the obsolete durable-ROOT shape
    # expected by the historical count-amendment validator.
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
    require_exact_keys(
        resource, evidence.CAPACITY_KEYS, "current capacity resource certificate"
    )

    if (
        resource.get("schema") != evidence.CAPACITY_SCHEMA
        or resource.get("status") != "PASS"
        or root_join.get("schema") != evidence.ROOT_JOIN_SCHEMA
        or root_join.get("status") != "PASS"
        or root_join.get("scope") != "capacity"
    ):
        raise BindingError("current capacity certificate schema/status differs")
    for label, payload in (
        ("capacity resource", resource),
        ("capacity root/join", root_join),
    ):
        if (
            payload.get("execution_group_size") != EXPECTED_GROUP_SIZE
            or payload.get("source_occurrences_per_output") != 1
            or payload.get("capacity_authority_earned") is not True
            or payload.get("full_training_authority") != 0
        ):
            raise BindingError(f"{label} authority/count contract differs")
    if (
        resource.get("submission_performed") is not True
        or resource.get("selected_rows")
        != list(evidence.SELECTED_CAPACITY_ROWS)
        or not capacity_resource_bundle_binding_valid(
            resource, immutable, plan_binding
        )
        or resource.get("execution_partition_sha256")
        != current["execution_partition_sha256"]
        or resource.get("root_health_identity_join_certificate_sha256")
        != root_join_artifact["sha256"]
    ):
        raise BindingError("capacity resource binding to current plan differs")
    try:
        evidence.require_same_resolved_file(
            "capacity root/join binding",
            resource.get("root_health_identity_join_certificate"),
            root_join_artifact["path"],
        )
    except evidence.AmendmentError as exc:
        raise BindingError(str(exc)) from exc

    validation_authority = require_mapping(
        resource.get("validation_authority"),
        "capacity resource.validation_authority",
    )
    require_exact_keys(
        validation_authority,
        CURRENT_VALIDATION_AUTHORITY_KEYS,
        "capacity resource.validation_authority",
    )
    try:
        for name in CURRENT_VALIDATION_AUTHORITY_KEYS:
            evidence.verify_nested_file_reference(
                f"capacity validation authority {name}",
                validation_authority[name],
                allowed_extra=(
                    ("validation_commit",) if name == "controller" else ()
                ),
            )
        evidence.require_same_resolved_file(
            "capacity validation full plan",
            validation_authority["full_extraction_plan"]["path"],
            resource.get("full_plan"),
        )
    except evidence.AmendmentError as exc:
        raise BindingError(str(exc)) from exc
    if (
        validation_authority["full_extraction_plan"].get("sha256")
        != resource.get("full_plan_sha256")
    ):
        raise BindingError("capacity validation authority differs")
    validation_commit_binding = validate_science_validation_commit_binding(
        immutable, validation_authority
    )
    validation_commit_binding["capacity_plan_binding_mode"] = plan_binding["mode"]
    validation_commit_binding["capacity_plan_sha256"] = plan_binding[
        "capacity_plan_sha256"
    ]
    validation_commit_binding["current_plan_sha256"] = plan_binding[
        "current_plan_sha256"
    ]
    validation_commit_binding["namespace_normalized_sha256"] = plan_binding[
        "namespace_normalized_sha256"
    ]
    validation_commit_binding["capacity_preflight_binding_mode"] = (
        preflight_binding["mode"]
    )
    validation_commit_binding["capacity_preflight_sha256"] = preflight_binding[
        "capacity_preflight_sha256"
    ]
    validation_commit_binding["current_preflight_sha256"] = preflight_binding[
        "current_preflight_sha256"
    ]

    audit_by_row = {
        "pp_background_jet8": pp_audit,
        "auau_background_jet12": auau_audit,
    }
    audit_artifact_by_row = {
        "pp_background_jet8": pp_audit_artifact,
        "auau_background_jet12": auau_audit_artifact,
    }
    for row_id, audit in audit_by_row.items():
        try:
            evidence.require_exact_keys(
                f"{row_id} capacity audit", audit, evidence.AUDIT_KEYS
            )
            evidence.validate_capacity_non_training_audit(
                row_id, audit, audit_artifact_by_row[row_id]
            )
            audit_root_binding = evidence.verify_nested_file_reference(
                f"{row_id} root/join audit binding",
                audit.get("root_health_identity_join_certificate"),
                allowed_extra=("schema", "scope", "status"),
            )
            evidence.require_same_resolved_file(
                f"{row_id} root/join audit binding",
                audit_root_binding["path"],
                root_join_artifact["path"],
            )
        except evidence.AmendmentError as exc:
            raise BindingError(str(exc)) from exc
        if (
            audit_root_binding.get("sha256")
            != root_join_artifact["sha256"]
            or audit_root_binding.get("schema") != evidence.ROOT_JOIN_SCHEMA
            or audit_root_binding.get("scope") != "capacity"
            or audit_root_binding.get("status") != "PASS"
        ):
            raise BindingError(f"{row_id} root/join audit binding differs")

    root_rows = root_join.get("rows")
    if not isinstance(root_rows, list):
        raise BindingError("capacity root/join rows must be a list")
    if [row.get("row_id") for row in root_rows if isinstance(row, dict)] != list(
        evidence.SELECTED_CAPACITY_ROWS
    ):
        raise BindingError("capacity root/join row order/identity differs")
    if (
        root_join.get("row_count") != 2
        or root_join.get("populated_row_count") != 1
        or root_join.get("valid_empty_row_count") != 1
        or root_join.get("populated_rows") != ["auau_background_jet12"]
        or root_join.get("valid_empty_rows") != ["pp_background_jet8"]
        or root_join.get("failures") != []
    ):
        raise BindingError("capacity root/join population contract differs")

    resource_rows = resource.get("rows")
    if not isinstance(resource_rows, list):
        raise BindingError("capacity resource rows must be a list")
    resource_by_row: dict[str, dict[str, Any]] = {}
    try:
        for raw in resource_rows:
            row = require_mapping(raw, "capacity resource row")
            row_id = str(row.get("row_id", ""))
            if row_id in resource_by_row:
                raise BindingError(f"capacity resource duplicates {row_id}")
            resource_by_row[row_id] = evidence.validate_capacity_resource_row(
                row_id, row
            )
    except evidence.AmendmentError as exc:
        raise BindingError(str(exc)) from exc
    if list(resource_by_row) != list(evidence.SELECTED_CAPACITY_ROWS):
        raise BindingError("capacity resource row order/identity differs")

    selected_rows: list[dict[str, Any]] = []
    for raw in root_rows:
        row = require_mapping(raw, "capacity root/join row")
        row_id = str(row.get("row_id", ""))
        require_exact_keys(
            row,
            CURRENT_ROOT_ROW_KEYS,
            f"capacity root/join row {row_id}",
        )
        audit = audit_by_row[row_id]
        expected_population = evidence.EXPECTED_POPULATION_STATE[row_id]
        if (
            row.get("status") != "PASS"
            or row.get("population_state") != expected_population
            or row.get("analysis_artifact_state")
            != "EPHEMERAL_VALIDATED_NOT_RETAINED"
            or row.get("full_training_authority") != 0
            or row.get("source_contract") != audit.get("source_contract")
            or row.get("source_execution_contract")
            != audit.get("source_execution_contract")
            or row.get("source_identity_canonical_sha256")
            != audit.get("source_identity_canonical_sha256")
            or row.get("source_occurrence_id_hex")
            != audit.get("source_occurrence_id_hex")
        ):
            raise BindingError(f"{row_id} capacity identity/state differs")

        analysis_path = Path(str(row.get("analysis_output_root", "")))
        if not analysis_path.is_absolute() or analysis_path.exists():
            raise BindingError(
                f"{row_id} ephemeral analysis output ownership differs"
            )
        analysis_bytes = row.get("analysis_bytes")
        if (
            isinstance(analysis_bytes, bool)
            or not isinstance(analysis_bytes, int)
            or analysis_bytes < 50_000
        ):
            raise BindingError(f"{row_id} analysis size witness differs")
        health = require_mapping(
            row.get("analysis_health_receipt"),
            f"{row_id}.analysis_health_receipt",
        )
        require_exact_keys(
            health,
            ANALYSIS_HEALTH_RECEIPT_KEYS,
            f"{row_id}.analysis_health_receipt",
        )
        health_path = Path(str(health.get("path", "")))
        health_sha256 = str(health.get("sha256", ""))
        if (
            not health_path.is_absolute()
            or not health_path.is_file()
            or file_sha256(health_path) != health_sha256
        ):
            raise BindingError(f"{row_id} analysis health receipt differs")
        payload = require_mapping(
            health.get("payload"), f"{row_id}.analysis_health_receipt.payload"
        )
        require_exact_keys(
            payload,
            ANALYSIS_HEALTH_PAYLOAD_KEYS,
            f"{row_id}.analysis_health_receipt.payload",
        )
        key_inventory_sha256 = str(
            row.get("analysis_key_inventory_sha256", "")
        )
        try:
            evidence.require_sha256(
                f"{row_id}.analysis_key_inventory_sha256",
                key_inventory_sha256,
            )
        except evidence.AmendmentError as exc:
            raise BindingError(str(exc)) from exc
        required_true = (
            "analysis_root_non_zombie",
            "analysis_root_non_recovered",
            "analysis_directory_present",
            "analysis_histogram_present",
            "analysis_config_present",
        )
        if (
            payload.get("schema") != "THE134_EPHEMERAL_ANALYSIS_HEALTH_V1"
            or payload.get("status") != "PASS"
            or payload.get("mode") != "EPHEMERAL_CONDOR_SCRATCH"
            or payload.get("analysis_root_retained") is not False
            or payload.get("analysis_size_bytes") != analysis_bytes
            or payload.get("analysis_minimum_bytes") != 50_000
            or payload.get("analysis_key_inventory_sha256")
            != key_inventory_sha256
            or any(payload.get(field) is not True for field in required_true)
        ):
            raise BindingError(f"{row_id} analysis health payload differs")

        sidecar_path = Path(str(row.get("sidecar", "")))
        sidecar_bytes = row.get("sidecar_bytes")
        sidecar_sha256 = str(row.get("sidecar_sha256", ""))
        if (
            not sidecar_path.is_absolute()
            or not sidecar_path.is_file()
            or not isinstance(sidecar_bytes, int)
            or sidecar_bytes <= 0
            or sidecar_path.stat().st_size != sidecar_bytes
            or file_sha256(sidecar_path) != sidecar_sha256
            or payload.get("sidecar_size_bytes") != sidecar_bytes
            or payload.get("sidecar_tree_entries") != row.get("sidecar_entries")
        ):
            raise BindingError(f"{row_id} durable sidecar hash/size differs")
        resource_row = resource_by_row[row_id]
        selected_rows.append(
            {
                "row_id": row_id,
                "system": evidence.EXPECTED_SYSTEM[row_id],
                "population_state": expected_population,
                "cluster_proc": resource_row["cluster_proc"],
                "exit_code": 0,
                "num_job_starts": 1,
                "num_holds": 0,
                "request_memory_mb": EXPECTED_REQUEST_MEMORY_MB,
                "memory_usage_mb": resource_row["memory_usage_mb"],
                "resident_set_size_kb": resource_row[
                    "resident_set_size_kb"
                ],
                "remote_wall_clock_seconds": resource_row[
                    "remote_wall_clock_seconds"
                ],
                "analysis_health": {
                    "path": str(analysis_path),
                    "artifact_state": row["analysis_artifact_state"],
                    "health_receipt_path": str(health_path),
                    "health_receipt_sha256": health_sha256,
                    "key_inventory_sha256": key_inventory_sha256,
                    "size_bytes": analysis_bytes,
                },
                "sidecar_output": {
                    "path": str(sidecar_path),
                    "sha256": sidecar_sha256,
                    "size_bytes": sidecar_bytes,
                },
            }
        )

    capacity_artifacts = {
        "resource_certificate": evidence.pinned_artifact_record(
            "capacity_resource_certificate",
            resource_artifact,
            schema=evidence.CAPACITY_SCHEMA,
        ),
        "root_join_certificate": evidence.pinned_artifact_record(
            "capacity_root_join_certificate",
            root_join_artifact,
            schema=evidence.ROOT_JOIN_SCHEMA,
        ),
        "pp_audit": evidence.pinned_artifact_record(
            "pp_capacity_audit",
            pp_audit_artifact,
            schema=evidence.CAPACITY_AUDIT_SCHEMA,
        ),
        "auau_audit": evidence.pinned_artifact_record(
            "auau_capacity_audit",
            auau_audit_artifact,
            schema=evidence.CAPACITY_AUDIT_SCHEMA,
        ),
        "validation_authority": validation_authority,
        "validation_commit_binding": validation_commit_binding,
    }
    return capacity_artifacts, selected_rows


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
        capacity_resource_for_binding, _capacity_resource_artifact = (
            evidence.load_json_artifact(
                "capacity resource certificate for plan binding",
                capacity_spec["resource_certificate"],
                expected_schema=evidence.CAPACITY_SCHEMA,
            )
        )
    except evidence.AmendmentError as exc:
        raise BindingError(str(exc)) from exc
    plan_binding = validate_capacity_plan_binding(
        capacity_resource_for_binding, preflight_spec
    )
    preflight_binding = validate_capacity_preflight_binding(
        capacity_resource_for_binding, preflight_spec
    )

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
        plan_binding=plan_binding,
        preflight_binding=preflight_binding,
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
            "ephemeral_analysis_health_receipts_revalidated": True,
            "durable_sidecar_hashes_revalidated": True,
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
