#!/usr/bin/env python3
"""Build and pre-submit revalidate one immutable THE-134 execution bundle.

This tool inventories the exact executors, libraries, models, configuration,
code/schema contracts, CaloReco runtime authorities, and declared dependency
edges consumed by the source-complete THE-134 resolver.  It performs no copy,
install, submission, or job-control action.

``build`` consumes a SHA-pinned JSON spec, rehashes every artifact, validates
all declared hash bindings and dependency providers, and emits the resolver's
``THE134_FULL_EXTRACTION_IMMUTABLE_BUNDLE_V1`` receipt.  The bundle identity is
the SHA-256 of a path-independent semantic inventory; the receipt derives the
full digest-named materialization path but does not create it.

``verify`` requires the exact receipt SHA-256 and repeats path-resolution,
content, size, dependency, hash-binding, and identity checks.  It is the
mandatory readback immediately before submission.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import stat
import sys
from pathlib import Path, PurePosixPath
from typing import Any, Mapping, Sequence


SPEC_SCHEMA = "THE134_IMMUTABLE_BUILD_BUNDLE_SPEC_V1"
BUNDLE_SCHEMA = "THE134_FULL_EXTRACTION_IMMUTABLE_BUNDLE_V1"
READBACK_SCHEMA = "THE134_FULL_EXTRACTION_IMMUTABLE_BUNDLE_READBACK_V1"
AUTHORITY_STATE = "MEASURED_CANDIDATE_NOT_SCIENTIFIC_AUTHORITY"
BUNDLE_NAME_PREFIX = "the134_bundle_sha256_"
REQUEST_MEMORY_MB = 3000
RELEASE = "ana.560"
CALO_RECO_SONAME = "libcalo_reco.so.0"
DECLARED_HASH_NAMES = (
    "code_sha256",
    "replay_schema_sha256",
    "training_schema_sha256",
    "semantic_sha256",
)
REQUIRED_ARTIFACT_ROLES = (
    "submitter",
    "pp_executor",
    "auau_executor",
    "pp_config",
    "auau_config",
    "pp_library",
    "auau_library",
    "pp_model",
    "auau_model",
    "code_manifest",
    "replay_schema_header",
    "training_schema_header",
    "shower_factorial_header",
    "photon_cluster_header",
    "photon_cluster_builder_header",
    "calo_reco_library",
    "calo_reco_build_receipt",
    "calo_reco_source_manifest",
    "runtime_authority_manifest",
    "release_calo_io",
    "release_clusteriso",
    "release_jetbase",
)
EXECUTABLE_ROLES = frozenset({"submitter", "pp_executor", "auau_executor"})
REQUIRED_DEPENDENCY_PROVIDER_ROLES = frozenset(
    {
        "pp_config",
        "auau_config",
        "pp_library",
        "auau_library",
        "pp_model",
        "auau_model",
        "code_manifest",
        "replay_schema_header",
        "training_schema_header",
        "shower_factorial_header",
        "photon_cluster_header",
        "photon_cluster_builder_header",
        "calo_reco_library",
        "calo_reco_build_receipt",
        "calo_reco_source_manifest",
        "runtime_authority_manifest",
        "release_calo_io",
        "release_clusteriso",
        "release_jetbase",
    }
)
HASH_BINDING_MODES = frozenset(
    {"artifact_sha256", "json_pointer", "literal_text"}
)
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
GIT_COMMIT_RE = re.compile(r"^[0-9a-f]{40}$")
ROLE_RE = re.compile(r"^[a-z][a-z0-9_]{1,63}$")
SAFE_KIND_RE = re.compile(r"^[a-z][a-z0-9_.-]{1,63}$")


class BundleError(ValueError):
    """Fail-closed immutable-bundle contract violation."""


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


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require_sha256(label: str, value: object) -> str:
    text = str(value)
    if not SHA256_RE.fullmatch(text):
        raise BundleError(f"{label} must be a lowercase 64-character SHA-256")
    return text


def require_positive_int(label: str, value: object) -> int:
    if isinstance(value, bool):
        raise BundleError(f"{label} must be a positive integer")
    try:
        parsed = int(value)
    except (TypeError, ValueError) as exc:
        raise BundleError(f"{label} must be a positive integer") from exc
    if parsed <= 0 or str(parsed) != str(value):
        raise BundleError(f"{label} must be a canonical positive integer")
    return parsed


def require_safe_text(label: str, value: object) -> str:
    text = str(value)
    if not text or any(character in text for character in "\x00\n\r\t"):
        raise BundleError(f"{label} is empty or contains a control character")
    return text


def require_absolute_file(label: str, value: object) -> Path:
    path = Path(require_safe_text(label, value)).expanduser()
    if not path.is_absolute():
        raise BundleError(f"{label} must be an absolute path: {path}")
    if not path.is_file():
        raise BundleError(f"{label} is not an existing regular file: {path}")
    return path


def require_absolute_dir(label: str, value: object) -> Path:
    path = Path(require_safe_text(label, value)).expanduser()
    if not path.is_absolute():
        raise BundleError(f"{label} must be an absolute path: {path}")
    if not path.is_dir():
        raise BundleError(f"{label} is not an existing directory: {path}")
    return path


def require_normalized_posix_path(label: str, value: object) -> str:
    text = require_safe_text(label, value)
    path = PurePosixPath(text)
    if (
        not path.is_absolute()
        or ".." in path.parts
        or str(path) != text.rstrip("/")
        or str(path) == "/"
    ):
        raise BundleError(
            f"{label} must be a normalized non-root absolute POSIX path"
        )
    return str(path)


def load_pinned_json(
    path: Path,
    expected_sha256: str,
    *,
    expected_schema: str | None = None,
) -> dict[str, Any]:
    expected = require_sha256("expected JSON SHA-256", expected_sha256)
    if not path.is_absolute() or not path.is_file():
        raise BundleError(f"JSON input must be an existing absolute file: {path}")
    observed = sha256_file(path)
    if observed != expected:
        raise BundleError(
            f"JSON file hash drift: expected={expected} observed={observed} "
            f"path={path}"
        )
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise BundleError(f"invalid JSON input: {path}") from exc
    if not isinstance(payload, dict):
        raise BundleError("JSON input must contain one object")
    if expected_schema is not None and payload.get("schema") != expected_schema:
        raise BundleError(
            f"JSON schema must be {expected_schema}, observed "
            f"{payload.get('schema')!r}"
        )
    return payload


def json_pointer_value(payload: Any, pointer: str) -> Any:
    if not pointer.startswith("/"):
        raise BundleError(f"JSON pointer must begin with '/': {pointer!r}")
    current = payload
    for raw in pointer.split("/")[1:]:
        token = raw.replace("~1", "/").replace("~0", "~")
        if isinstance(current, dict):
            if token not in current:
                raise BundleError(f"JSON pointer field is missing: {pointer}")
            current = current[token]
        elif isinstance(current, list):
            if not token.isdigit() or int(token) >= len(current):
                raise BundleError(f"JSON pointer index is invalid: {pointer}")
            current = current[int(token)]
        else:
            raise BundleError(f"JSON pointer traverses a scalar: {pointer}")
    return current


def inspect_artifacts(raw_artifacts: object) -> list[dict[str, Any]]:
    if not isinstance(raw_artifacts, list):
        raise BundleError("artifacts must be a list")
    by_role: dict[str, dict[str, Any]] = {}
    resolved_paths: dict[str, str] = {}
    for raw in raw_artifacts:
        if not isinstance(raw, dict):
            raise BundleError("artifact records must be objects")
        if set(raw) != {"role", "path", "sha256", "size_bytes"}:
            raise BundleError(
                f"artifact record fields differ for role={raw.get('role')!r}"
            )
        role = str(raw.get("role", ""))
        if not ROLE_RE.fullmatch(role) or role in by_role:
            raise BundleError(f"invalid or duplicate artifact role: {role!r}")
        path = require_absolute_file(f"artifact {role}", raw["path"])
        resolved = str(path.resolve(strict=True))
        if resolved in resolved_paths:
            raise BundleError(
                f"artifact roles alias one physical file: "
                f"{resolved_paths[resolved]} and {role}"
            )
        declared_sha = require_sha256(f"artifact {role}", raw["sha256"])
        observed_sha = sha256_file(path)
        if observed_sha != declared_sha:
            raise BundleError(
                f"artifact hash mismatch: role={role} expected={declared_sha} "
                f"observed={observed_sha}"
            )
        declared_size = require_positive_int(
            f"artifact {role} size_bytes", raw["size_bytes"]
        )
        observed_size = path.stat().st_size
        if observed_size != declared_size:
            raise BundleError(
                f"artifact size mismatch: role={role} expected={declared_size} "
                f"observed={observed_size}"
            )
        if role in EXECUTABLE_ROLES and not (
            path.stat().st_mode & stat.S_IXUSR
        ):
            raise BundleError(f"artifact {role} is not owner-executable: {path}")
        by_role[role] = {
            "role": role,
            "path": str(path),
            "resolved_path": resolved,
            "sha256": observed_sha,
            "size_bytes": observed_size,
        }
        resolved_paths[resolved] = role
    expected_roles = set(REQUIRED_ARTIFACT_ROLES)
    if set(by_role) != expected_roles:
        raise BundleError(
            "artifact role closure differs: "
            f"missing={sorted(expected_roles - set(by_role))} "
            f"extra={sorted(set(by_role) - expected_roles)}"
        )
    return [by_role[role] for role in sorted(by_role)]


def index_artifacts(
    artifacts: Sequence[Mapping[str, Any]],
) -> dict[str, Mapping[str, Any]]:
    return {str(artifact["role"]): artifact for artifact in artifacts}


def normalize_dependencies(
    raw_dependencies: object,
    artifacts: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    if not isinstance(raw_dependencies, list):
        raise BundleError("dependencies must be a list")
    artifact_by_role = index_artifacts(artifacts)
    dependencies: dict[str, dict[str, Any]] = {}
    for raw in raw_dependencies:
        if not isinstance(raw, dict):
            raise BundleError("dependency records must be objects")
        required_fields = {
            "provider_role",
            "consumer_roles",
            "kind",
            "required",
            "sha256",
        }
        if set(raw) != required_fields:
            raise BundleError(
                f"dependency record fields differ for "
                f"provider={raw.get('provider_role')!r}"
            )
        provider = str(raw["provider_role"])
        if provider not in artifact_by_role or provider in dependencies:
            raise BundleError(
                f"dependency provider is missing or duplicated: {provider!r}"
            )
        consumers_raw = raw["consumer_roles"]
        if not isinstance(consumers_raw, list) or not consumers_raw:
            raise BundleError(
                f"dependency {provider} consumer_roles must be nonempty"
            )
        consumers = sorted({str(value) for value in consumers_raw})
        if len(consumers) != len(consumers_raw):
            raise BundleError(f"dependency {provider} duplicates consumers")
        if any(
            consumer not in artifact_by_role or consumer == provider
            for consumer in consumers
        ):
            raise BundleError(
                f"dependency {provider} has invalid consumer roles: {consumers}"
            )
        kind = str(raw["kind"])
        if not SAFE_KIND_RE.fullmatch(kind):
            raise BundleError(f"dependency {provider} has invalid kind: {kind!r}")
        if raw["required"] is not True:
            raise BundleError(f"dependency {provider} must be required")
        declared_sha = require_sha256(
            f"dependency {provider}", raw["sha256"]
        )
        artifact_sha = str(artifact_by_role[provider]["sha256"])
        if declared_sha != artifact_sha:
            raise BundleError(
                f"dependency hash mismatch: provider={provider} "
                f"expected_artifact={artifact_sha} declared={declared_sha}"
            )
        dependencies[provider] = {
            "provider_role": provider,
            "consumer_roles": consumers,
            "kind": kind,
            "required": True,
            "sha256": artifact_sha,
        }
    observed = set(dependencies)
    if observed != REQUIRED_DEPENDENCY_PROVIDER_ROLES:
        raise BundleError(
            "dependency provider closure differs: "
            f"missing={sorted(REQUIRED_DEPENDENCY_PROVIDER_ROLES - observed)} "
            f"extra={sorted(observed - REQUIRED_DEPENDENCY_PROVIDER_ROLES)}"
        )
    return [dependencies[role] for role in sorted(dependencies)]


def normalize_declared_hashes(raw_hashes: object) -> dict[str, str]:
    if not isinstance(raw_hashes, dict) or set(raw_hashes) != set(
        DECLARED_HASH_NAMES
    ):
        raise BundleError(
            f"declared_hashes must contain exactly {list(DECLARED_HASH_NAMES)}"
        )
    return {
        name: require_sha256(f"declared hash {name}", raw_hashes[name])
        for name in DECLARED_HASH_NAMES
    }


def binding_observed_value(
    binding: Mapping[str, Any],
    artifact: Mapping[str, Any],
) -> str:
    mode = str(binding["mode"])
    path = Path(str(artifact["path"]))
    if mode == "artifact_sha256":
        if "pointer" in binding:
            raise BundleError("artifact_sha256 binding cannot contain pointer")
        return str(artifact["sha256"])
    if mode == "literal_text":
        if "pointer" in binding:
            raise BundleError("literal_text binding cannot contain pointer")
        return path.read_text(encoding="utf-8").strip()
    if mode == "json_pointer":
        pointer = str(binding.get("pointer", ""))
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise BundleError(
                f"hash-binding artifact is not JSON: {path}"
            ) from exc
        return str(json_pointer_value(payload, pointer))
    raise BundleError(f"unsupported hash-binding mode: {mode!r}")


def normalize_hash_bindings(
    raw_bindings: object,
    declared_hashes: Mapping[str, str],
    artifacts: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    if not isinstance(raw_bindings, list):
        raise BundleError("hash_bindings must be a list")
    artifact_by_role = index_artifacts(artifacts)
    bindings: dict[str, dict[str, Any]] = {}
    for raw in raw_bindings:
        if not isinstance(raw, dict):
            raise BundleError("hash-binding records must be objects")
        allowed_fields = {"name", "artifact_role", "mode", "pointer"}
        if not {"name", "artifact_role", "mode"}.issubset(raw) or not set(
            raw
        ).issubset(allowed_fields):
            raise BundleError(
                f"hash-binding fields differ for name={raw.get('name')!r}"
            )
        name = str(raw["name"])
        if name not in declared_hashes or name in bindings:
            raise BundleError(f"invalid or duplicate hash binding: {name!r}")
        role = str(raw["artifact_role"])
        if role not in artifact_by_role:
            raise BundleError(
                f"hash binding {name} references unknown artifact role {role!r}"
            )
        mode = str(raw["mode"])
        if mode not in HASH_BINDING_MODES:
            raise BundleError(
                f"hash binding {name} has unsupported mode {mode!r}"
            )
        normalized: dict[str, Any] = {
            "name": name,
            "artifact_role": role,
            "mode": mode,
        }
        if "pointer" in raw:
            normalized["pointer"] = str(raw["pointer"])
        observed = binding_observed_value(normalized, artifact_by_role[role])
        expected = declared_hashes[name]
        if observed != expected:
            raise BundleError(
                f"declared hash binding mismatch: name={name} "
                f"expected={expected} observed={observed}"
            )
        normalized["value"] = expected
        normalized["evidence_artifact_sha256"] = artifact_by_role[role][
            "sha256"
        ]
        bindings[name] = normalized
    if set(bindings) != set(DECLARED_HASH_NAMES):
        raise BundleError(
            "hash-binding closure differs: "
            f"missing={sorted(set(DECLARED_HASH_NAMES) - set(bindings))}"
        )
    return [bindings[name] for name in DECLARED_HASH_NAMES]


def normalize_runtime(raw_runtime: object) -> dict[str, Any]:
    if not isinstance(raw_runtime, dict):
        raise BundleError("runtime must be an object")
    required = {
        "release",
        "offline_main",
        "calo_reco_soname",
        "request_memory_mb",
        "release_core_lib_dir",
        "release_core_lib64_dir",
    }
    if set(raw_runtime) != required:
        raise BundleError(
            "runtime field inventory differs: "
            f"missing={sorted(required - set(raw_runtime))} "
            f"extra={sorted(set(raw_runtime) - required)}"
        )
    if raw_runtime["release"] != RELEASE:
        raise BundleError(f"runtime release must remain {RELEASE}")
    if raw_runtime["calo_reco_soname"] != CALO_RECO_SONAME:
        raise BundleError(
            f"runtime calo_reco_soname must remain {CALO_RECO_SONAME}"
        )
    memory = require_positive_int(
        "runtime request_memory_mb", raw_runtime["request_memory_mb"]
    )
    if memory != REQUEST_MEMORY_MB:
        raise BundleError(
            f"runtime request_memory_mb must remain {REQUEST_MEMORY_MB}"
        )
    return {
        "release": RELEASE,
        "offline_main": require_normalized_posix_path(
            "runtime offline_main", raw_runtime["offline_main"]
        ),
        "calo_reco_soname": CALO_RECO_SONAME,
        "request_memory_mb": REQUEST_MEMORY_MB,
        "release_core_lib_dir": str(
            require_absolute_dir(
                "runtime release_core_lib_dir",
                raw_runtime["release_core_lib_dir"],
            )
        ),
        "release_core_lib64_dir": str(
            require_absolute_dir(
                "runtime release_core_lib64_dir",
                raw_runtime["release_core_lib64_dir"],
            )
        ),
    }


def semantic_artifact_inventory(
    artifacts: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    return [
        {
            "role": artifact["role"],
            "sha256": artifact["sha256"],
            "size_bytes": artifact["size_bytes"],
        }
        for artifact in artifacts
    ]


def semantic_binding_inventory(
    bindings: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    return [
        {
            key: binding[key]
            for key in (
                "name",
                "artifact_role",
                "mode",
                "value",
                "evidence_artifact_sha256",
            )
            if key in binding
        }
        | ({"pointer": binding["pointer"]} if "pointer" in binding else {})
        for binding in bindings
    ]


def bundle_identity_payload(
    *,
    public_commit: str,
    declared_hashes: Mapping[str, str],
    runtime: Mapping[str, Any],
    artifacts: Sequence[Mapping[str, Any]],
    dependencies: Sequence[Mapping[str, Any]],
    bindings: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    return {
        "authority_state": AUTHORITY_STATE,
        "public_commit": public_commit,
        "declared_hashes": dict(declared_hashes),
        "runtime": {
            "release": runtime["release"],
            "offline_main": runtime["offline_main"],
            "calo_reco_soname": runtime["calo_reco_soname"],
            "request_memory_mb": runtime["request_memory_mb"],
        },
        "artifacts": semantic_artifact_inventory(artifacts),
        "dependencies": list(dependencies),
        "hash_bindings": semantic_binding_inventory(bindings),
    }


def build_receipt(spec: Mapping[str, Any]) -> dict[str, Any]:
    required_top = {
        "schema",
        "public_commit",
        "bundle_parent",
        "runtime",
        "declared_hashes",
        "artifacts",
        "dependencies",
        "hash_bindings",
    }
    if set(spec) != required_top:
        raise BundleError(
            "bundle spec top-level inventory differs: "
            f"missing={sorted(required_top - set(spec))} "
            f"extra={sorted(set(spec) - required_top)}"
        )
    if spec["schema"] != SPEC_SCHEMA:
        raise BundleError(f"bundle spec schema must be {SPEC_SCHEMA}")
    public_commit = str(spec["public_commit"])
    if not GIT_COMMIT_RE.fullmatch(public_commit):
        raise BundleError("public_commit must be a full lowercase Git SHA")
    bundle_parent = require_normalized_posix_path(
        "bundle_parent", spec["bundle_parent"]
    )
    runtime = normalize_runtime(spec["runtime"])
    declared_hashes = normalize_declared_hashes(spec["declared_hashes"])
    artifacts = inspect_artifacts(spec["artifacts"])
    dependencies = normalize_dependencies(spec["dependencies"], artifacts)
    bindings = normalize_hash_bindings(
        spec["hash_bindings"], declared_hashes, artifacts
    )
    identity_payload = bundle_identity_payload(
        public_commit=public_commit,
        declared_hashes=declared_hashes,
        runtime=runtime,
        artifacts=artifacts,
        dependencies=dependencies,
        bindings=bindings,
    )
    identity = canonical_sha256(identity_payload)
    bundle_name = f"{BUNDLE_NAME_PREFIX}{identity}"
    digest_path = str(PurePosixPath(bundle_parent) / bundle_name)
    return {
        "schema": BUNDLE_SCHEMA,
        "status": "PASS",
        "authority_state": AUTHORITY_STATE,
        "public_commit": public_commit,
        **declared_hashes,
        "runtime": runtime,
        "artifacts": artifacts,
        "dependencies": dependencies,
        "hash_bindings": bindings,
        "bundle_identity_sha256": identity,
        "bundle_name": bundle_name,
        "bundle_parent": bundle_parent,
        "digest_named_bundle_path": digest_path,
        "digest_named_bundle_materialized": False,
        "semantic_fingerprint_sha256": identity,
    }


def receipt_identity_payload(receipt: Mapping[str, Any]) -> dict[str, Any]:
    declared = {name: receipt[name] for name in DECLARED_HASH_NAMES}
    return bundle_identity_payload(
        public_commit=str(receipt["public_commit"]),
        declared_hashes=declared,
        runtime=receipt["runtime"],
        artifacts=receipt["artifacts"],
        dependencies=receipt["dependencies"],
        bindings=receipt["hash_bindings"],
    )


def validate_receipt_payload(
    receipt: Mapping[str, Any],
    *,
    rehash: bool,
) -> dict[str, Any]:
    required_top = {
        "schema",
        "status",
        "authority_state",
        "public_commit",
        *DECLARED_HASH_NAMES,
        "runtime",
        "artifacts",
        "dependencies",
        "hash_bindings",
        "bundle_identity_sha256",
        "bundle_name",
        "bundle_parent",
        "digest_named_bundle_path",
        "digest_named_bundle_materialized",
        "semantic_fingerprint_sha256",
    }
    if set(receipt) != required_top:
        raise BundleError(
            "bundle receipt top-level inventory differs: "
            f"missing={sorted(required_top - set(receipt))} "
            f"extra={sorted(set(receipt) - required_top)}"
        )
    if receipt["schema"] != BUNDLE_SCHEMA or receipt["status"] != "PASS":
        raise BundleError("bundle receipt schema/status is invalid")
    if receipt["authority_state"] != AUTHORITY_STATE:
        raise BundleError("bundle receipt cannot grant scientific authority")
    public_commit = str(receipt["public_commit"])
    if not GIT_COMMIT_RE.fullmatch(public_commit):
        raise BundleError("bundle receipt public_commit is malformed")
    runtime = normalize_runtime(receipt["runtime"])
    declared_hashes = normalize_declared_hashes(
        {name: receipt[name] for name in DECLARED_HASH_NAMES}
    )

    if rehash:
        artifact_specs = [
            {
                "role": artifact["role"],
                "path": artifact["path"],
                "sha256": artifact["sha256"],
                "size_bytes": artifact["size_bytes"],
            }
            for artifact in receipt["artifacts"]
        ]
        artifacts = inspect_artifacts(artifact_specs)
        recorded_resolved = {
            artifact["role"]: artifact["resolved_path"]
            for artifact in receipt["artifacts"]
        }
        for artifact in artifacts:
            if artifact["resolved_path"] != recorded_resolved[artifact["role"]]:
                raise BundleError(
                    f"mutable artifact path drift: role={artifact['role']}"
                )
    else:
        artifacts = [dict(artifact) for artifact in receipt["artifacts"]]
        if {artifact.get("role") for artifact in artifacts} != set(
            REQUIRED_ARTIFACT_ROLES
        ):
            raise BundleError("bundle receipt artifact role closure differs")
        for artifact in artifacts:
            require_sha256(
                f"receipt artifact {artifact.get('role')}",
                artifact.get("sha256", ""),
            )
            require_positive_int(
                f"receipt artifact {artifact.get('role')} size_bytes",
                artifact.get("size_bytes", 0),
            )

    dependencies = normalize_dependencies(receipt["dependencies"], artifacts)
    bindings_spec = [
        {
            key: binding[key]
            for key in ("name", "artifact_role", "mode", "pointer")
            if key in binding
        }
        for binding in receipt["hash_bindings"]
    ]
    bindings = normalize_hash_bindings(
        bindings_spec, declared_hashes, artifacts
    )
    if semantic_binding_inventory(bindings) != semantic_binding_inventory(
        receipt["hash_bindings"]
    ):
        raise BundleError("bundle receipt hash-binding evidence differs")

    identity_payload = bundle_identity_payload(
        public_commit=public_commit,
        declared_hashes=declared_hashes,
        runtime=runtime,
        artifacts=artifacts,
        dependencies=dependencies,
        bindings=bindings,
    )
    identity = canonical_sha256(identity_payload)
    if require_sha256(
        "bundle_identity_sha256", receipt["bundle_identity_sha256"]
    ) != identity:
        raise BundleError("bundle identity digest differs")
    if receipt["semantic_fingerprint_sha256"] != identity:
        raise BundleError("semantic fingerprint differs from bundle identity")
    expected_name = f"{BUNDLE_NAME_PREFIX}{identity}"
    if receipt["bundle_name"] != expected_name:
        raise BundleError("digest-named bundle name differs")
    parent = require_normalized_posix_path(
        "receipt bundle_parent", receipt["bundle_parent"]
    )
    if receipt["digest_named_bundle_path"] != str(
        PurePosixPath(parent) / expected_name
    ):
        raise BundleError("digest-named bundle path differs")
    if receipt["digest_named_bundle_materialized"] is not False:
        raise BundleError(
            "receipt builder must not claim bundle materialization"
        )
    return dict(receipt)


def load_receipt(path: Path, expected_sha256: str) -> dict[str, Any]:
    payload = load_pinned_json(
        path, expected_sha256, expected_schema=BUNDLE_SCHEMA
    )
    return validate_receipt_payload(payload, rehash=True)


def atomic_write_json(path: Path, payload: Mapping[str, Any]) -> None:
    if not path.is_absolute():
        raise BundleError(f"output path must be absolute: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f"{path.name}.tmp.{os.getpid()}")
    temporary.write_bytes(canonical_json_bytes(payload))
    os.replace(temporary, path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    build = subparsers.add_parser("build", help="build one immutable receipt")
    build.add_argument("--spec", type=Path, required=True)
    build.add_argument("--expected-spec-sha256", required=True)
    build.add_argument("--output", type=Path, required=True)

    verify = subparsers.add_parser(
        "verify", help="perform the pinned pre-submit rehash/readback"
    )
    verify.add_argument("--receipt", type=Path, required=True)
    verify.add_argument("--expected-receipt-sha256", required=True)
    verify.add_argument("--output-json", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        if args.command == "build":
            spec = load_pinned_json(
                args.spec,
                args.expected_spec_sha256,
                expected_schema=SPEC_SCHEMA,
            )
            receipt = build_receipt(spec)
            atomic_write_json(args.output, receipt)
            summary = {
                "schema": BUNDLE_SCHEMA,
                "status": "PASS",
                "authority_state": AUTHORITY_STATE,
                "path": str(args.output),
                "sha256": sha256_file(args.output),
                "bundle_identity_sha256": receipt[
                    "bundle_identity_sha256"
                ],
                "bundle_name": receipt["bundle_name"],
                "artifact_count": len(receipt["artifacts"]),
                "dependency_count": len(receipt["dependencies"]),
                "hash_binding_count": len(receipt["hash_bindings"]),
            }
        else:
            receipt = load_receipt(
                args.receipt, args.expected_receipt_sha256
            )
            summary = {
                "schema": READBACK_SCHEMA,
                "status": "PASS",
                "authority_state": AUTHORITY_STATE,
                "receipt": str(args.receipt),
                "receipt_sha256": args.expected_receipt_sha256,
                "bundle_identity_sha256": receipt[
                    "bundle_identity_sha256"
                ],
                "bundle_name": receipt["bundle_name"],
                "artifacts_rehashed": len(receipt["artifacts"]),
                "dependencies_revalidated": len(receipt["dependencies"]),
                "hash_bindings_revalidated": len(receipt["hash_bindings"]),
                "pre_submit_readback": "PASS",
            }
            if args.output_json is not None:
                atomic_write_json(args.output_json, summary)
        print(json.dumps(summary, sort_keys=True))
        return 0
    except BundleError as exc:
        print(f"THE134_IMMUTABLE_BUNDLE_ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
