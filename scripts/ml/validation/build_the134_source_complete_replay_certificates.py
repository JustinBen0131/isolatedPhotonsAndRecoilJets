#!/usr/bin/env python3
"""Build the exact 14 THE-134 system/view replay certificates.

This is a receipt orchestrator, not a scientific validator.  It accepts only
explicit, already-produced per-source validation receipts and registered
artifact paths/hashes.  Every receipt, registry, and artifact is rehashed
before any output is written.  Gate values are copied from exact PASS receipts;
the builder never infers or manufactures scientific closure.

The output directory is published atomically and contains exactly fourteen
``THE134_FACTORIAL_VIEW_REPLAY_CERTIFICATE_V1`` JSON files plus one deterministic
index.  Re-running against an identical existing directory is idempotent;
drift or extra files fail closed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import tempfile
from pathlib import Path
from typing import Any, Iterable

from the134_aggregate_contract import (
    PAIRED_REGISTRY_SCHEMA,
    REPLAY_CERTIFICATE_GATES,
    REPLAY_CERTIFICATE_SCHEMA,
    REQUIRED_SOURCES_BY_SYSTEM,
    REQUIRED_SYSTEM_VIEW_KEYS,
    REQUIRED_VIEW_KEYS,
    SOURCE_WITNESS_GATES,
    ContractError,
    expected_view_semantic,
    read_json,
    require_exact_key_inventory,
    require_exact_true_gates,
    require_git_commit,
    require_sha256,
    sha256_file,
    verify_registered_file,
    verify_registered_json,
)


INPUT_MANIFEST_SCHEMA = "THE134_SOURCE_COMPLETE_REPLAY_INPUT_MANIFEST_V1"
SOURCE_VALIDATION_RECEIPT_SCHEMA = (
    "THE134_SOURCE_VIEW_REPLAY_VALIDATION_RECEIPT_V1"
)
OUTPUT_INDEX_SCHEMA = "THE134_SOURCE_COMPLETE_REPLAY_CERTIFICATE_INDEX_V1"
PROMOTION_STATUS = "CANDIDATE_CURRENT_NOT_CANONICAL"

INPUT_MANIFEST_KEYS = (
    "schema",
    "public_commit",
    "code_sha256",
    "paired_registries",
    "source_validations",
)
REGISTRY_RECORD_KEYS = ("view", "path", "sha256")
SOURCE_RECORD_KEYS = (
    "system",
    "view",
    "source",
    "validation_receipt",
    "artifacts",
)
FILE_RECORD_KEYS = ("path", "sha256")
ARTIFACT_RECORD_KEYS = (
    "source_manifest",
    "direct",
    "writer",
    "cache_receipt",
)
SOURCE_VALIDATION_RECEIPT_KEYS = (
    "schema",
    "status",
    "promotion_status",
    "system",
    "shower_definition",
    "shower_semantic_sha256",
    "source",
    "public_commit",
    "code_sha256",
    "source_manifest_sha256",
    "direct_artifact_sha256",
    "writer_artifact_sha256",
    "cache_receipt_sha256",
    "gates",
)
ARTIFACT_TO_RECEIPT_HASH = {
    "source_manifest": "source_manifest_sha256",
    "direct": "direct_artifact_sha256",
    "writer": "writer_artifact_sha256",
    "cache_receipt": "cache_receipt_sha256",
}
HEALTH_GATES = ("direct_health", "writer_health", "cache_health")
REQUIRED_SOURCE_VIEW_KEYS = tuple(
    (system, view, source)
    for system, view in REQUIRED_SYSTEM_VIEW_KEYS
    for source in REQUIRED_SOURCES_BY_SYSTEM[system]
)
INDEX_NAME = "replay_certificate_index.json"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-manifest", required=True, type=Path)
    parser.add_argument("--output-directory", required=True, type=Path)
    return parser.parse_args()


def require_exact_keys(
    payload: object,
    expected: Iterable[str],
    label: str,
) -> dict[str, Any]:
    if not isinstance(payload, dict):
        raise ContractError(f"{label} must be an object")
    expected_set = set(expected)
    observed_set = set(payload)
    if observed_set != expected_set:
        missing = sorted(expected_set - observed_set)
        extra = sorted(observed_set - expected_set)
        raise ContractError(
            f"{label} field inventory mismatch: missing={missing} extra={extra}"
        )
    return payload


def verify_file_record(record: object, label: str) -> dict[str, str]:
    normalized = require_exact_keys(record, FILE_RECORD_KEYS, label)
    path, observed_sha = verify_registered_file(
        normalized["path"],
        normalized["sha256"],
        label,
    )
    return {"path": str(path.resolve()), "sha256": observed_sha}


def verify_registry(
    view: str,
    record: dict[str, Any],
    public_commit: str,
    code_sha256: str,
) -> dict[str, str]:
    require_exact_keys(record, REGISTRY_RECORD_KEYS, f"{view} registry record")
    if record.get("view") != view:
        raise ContractError(f"{view} registry record view mismatch")
    path, payload, observed_sha = verify_registered_json(
        record, f"{view} paired registry"
    )
    expected = {
        "schema": PAIRED_REGISTRY_SCHEMA,
        "status": "READY",
        "promotion_status": PROMOTION_STATUS,
        "shower_definition": view,
        "shower_semantic_sha256": expected_view_semantic(view),
        "public_commit": public_commit,
        "code_sha256": code_sha256,
    }
    failures = [
        field for field, wanted in expected.items() if payload.get(field) != wanted
    ]
    if failures:
        raise ContractError(
            f"{view} paired registry identity mismatch: {failures}"
        )
    return {"path": str(path.resolve()), "sha256": observed_sha}


def verify_source_validation(
    key: tuple[str, str, str],
    record: dict[str, Any],
    public_commit: str,
    code_sha256: str,
) -> dict[str, Any]:
    system, view, source = key
    label = f"{system}/{view}/{source}"
    require_exact_keys(record, SOURCE_RECORD_KEYS, f"{label} source record")
    if (
        record.get("system"),
        record.get("view"),
        record.get("source"),
    ) != key:
        raise ContractError(f"{label} source-record identity mismatch")

    receipt_record = require_exact_keys(
        record.get("validation_receipt"),
        FILE_RECORD_KEYS,
        f"{label} validation receipt registration",
    )
    receipt_path, receipt, receipt_sha = verify_registered_json(
        receipt_record, f"{label} validation receipt"
    )
    require_exact_keys(
        receipt,
        SOURCE_VALIDATION_RECEIPT_KEYS,
        f"{label} validation receipt",
    )

    expected_identity = {
        "schema": SOURCE_VALIDATION_RECEIPT_SCHEMA,
        "status": "PASS",
        "promotion_status": PROMOTION_STATUS,
        "system": system,
        "shower_definition": view,
        "shower_semantic_sha256": expected_view_semantic(view),
        "source": source,
        "public_commit": public_commit,
        "code_sha256": code_sha256,
    }
    identity_failures = [
        field
        for field, wanted in expected_identity.items()
        if receipt.get(field) != wanted
    ]
    if identity_failures:
        raise ContractError(
            f"{label} validation receipt identity mismatch: "
            f"{identity_failures}"
        )
    require_exact_true_gates(
        receipt.get("gates"),
        SOURCE_WITNESS_GATES,
        f"{label} validation receipt",
    )

    artifact_records = require_exact_keys(
        record.get("artifacts"),
        ARTIFACT_RECORD_KEYS,
        f"{label} artifacts",
    )
    artifacts: dict[str, dict[str, str]] = {}
    for artifact_name in ARTIFACT_RECORD_KEYS:
        artifact = verify_file_record(
            artifact_records[artifact_name],
            f"{label} {artifact_name}",
        )
        receipt_field = ARTIFACT_TO_RECEIPT_HASH[artifact_name]
        receipt_sha_value = require_sha256(
            receipt.get(receipt_field),
            f"{label} validation receipt {receipt_field}",
        )
        if artifact["sha256"] != receipt_sha_value:
            raise ContractError(
                f"{label} {artifact_name} registration disagrees with "
                f"validation receipt"
            )
        artifacts[artifact_name] = artifact

    return {
        "system": system,
        "shower_definition": view,
        "source": source,
        "validation_receipt": str(receipt_path.resolve()),
        "validation_receipt_sha256": receipt_sha,
        "source_manifest": artifacts["source_manifest"]["path"],
        "source_manifest_sha256": artifacts["source_manifest"]["sha256"],
        "direct_artifact": artifacts["direct"]["path"],
        "direct_artifact_sha256": artifacts["direct"]["sha256"],
        "writer_artifact": artifacts["writer"]["path"],
        "writer_artifact_sha256": artifacts["writer"]["sha256"],
        "cache_receipt": artifacts["cache_receipt"]["path"],
        "cache_receipt_sha256": artifacts["cache_receipt"]["sha256"],
        "gates": {
            gate: receipt["gates"][gate] for gate in SOURCE_WITNESS_GATES
        },
    }


def derive_certificate_gates(
    witnesses: list[dict[str, Any]],
    label: str,
) -> dict[str, bool]:
    derived: dict[str, bool] = {}
    for gate in REPLAY_CERTIFICATE_GATES:
        if gate == "artifact_health":
            derived[gate] = all(
                witness["gates"][health_gate] is True
                for witness in witnesses
                for health_gate in HEALTH_GATES
            )
        else:
            derived[gate] = all(
                witness["gates"][gate] is True for witness in witnesses
            )
    require_exact_true_gates(derived, REPLAY_CERTIFICATE_GATES, label)
    return derived


def certificate_name(system: str, view: str) -> str:
    return f"the134_replay_certificate_{system}_{view}.json"


def json_bytes(payload: dict[str, Any]) -> bytes:
    return (
        json.dumps(payload, indent=2, sort_keys=True) + "\n"
    ).encode("utf-8")


def build_outputs(
    manifest_path: Path,
    output_directory: Path,
) -> dict[str, bytes]:
    manifest = read_json(manifest_path, "source-complete replay input manifest")
    require_exact_keys(
        manifest,
        INPUT_MANIFEST_KEYS,
        "source-complete replay input manifest",
    )
    if manifest.get("schema") != INPUT_MANIFEST_SCHEMA:
        raise ContractError(
            f"input manifest schema mismatch: {manifest.get('schema')!r}"
        )
    public_commit = require_git_commit(manifest.get("public_commit"))
    code_sha256 = require_sha256(manifest.get("code_sha256"), "code_sha256")

    registry_records = require_exact_key_inventory(
        manifest.get("paired_registries"),
        REQUIRED_VIEW_KEYS,
        lambda record: record.get("view"),
        "paired registries",
    )
    registries = {
        view: verify_registry(
            view,
            registry_records[view],
            public_commit,
            code_sha256,
        )
        for view in REQUIRED_VIEW_KEYS
    }

    source_records = require_exact_key_inventory(
        manifest.get("source_validations"),
        REQUIRED_SOURCE_VIEW_KEYS,
        lambda record: (
            record.get("system"),
            record.get("view"),
            record.get("source"),
        ),
        "source validation records",
    )
    source_validations = {
        key: verify_source_validation(
            key,
            source_records[key],
            public_commit,
            code_sha256,
        )
        for key in REQUIRED_SOURCE_VIEW_KEYS
    }

    stable_bindings: dict[tuple[str, str], dict[str, str]] = {}
    for system, view, source in REQUIRED_SOURCE_VIEW_KEYS:
        witness = source_validations[(system, view, source)]
        stable = {
            field: witness[field]
            for field in (
                "source_manifest_sha256",
                "direct_artifact_sha256",
                "writer_artifact_sha256",
            )
        }
        key = (system, source)
        previous = stable_bindings.setdefault(key, stable)
        if previous != stable:
            raise ContractError(
                f"{system}/{source} source/direct/writer binding changes "
                "across shower views"
            )

    output_root = output_directory.expanduser().resolve()
    certificates: dict[tuple[str, str], dict[str, Any]] = {}
    outputs: dict[str, bytes] = {}
    for system, view in REQUIRED_SYSTEM_VIEW_KEYS:
        witnesses = [
            source_validations[(system, view, source)]
            for source in REQUIRED_SOURCES_BY_SYSTEM[system]
        ]
        payload = {
            "schema": REPLAY_CERTIFICATE_SCHEMA,
            "status": "PASS",
            "promotion_status": PROMOTION_STATUS,
            "system": system,
            "shower_definition": view,
            "shower_semantic_sha256": expected_view_semantic(view),
            "paired_registry_sha256": registries[view]["sha256"],
            "public_commit": public_commit,
            "code_sha256": code_sha256,
            "gates": derive_certificate_gates(
                witnesses, f"{system}/{view} replay certificate"
            ),
            "source_witnesses": witnesses,
        }
        certificates[(system, view)] = payload
        outputs[certificate_name(system, view)] = json_bytes(payload)

    index_records = []
    for system, view in REQUIRED_SYSTEM_VIEW_KEYS:
        name = certificate_name(system, view)
        index_records.append(
            {
                "system": system,
                "view": view,
                "path": str(output_root / name),
                "sha256": sha256_bytes(outputs[name]),
            }
        )
    index = {
        "schema": OUTPUT_INDEX_SCHEMA,
        "status": "PASS",
        "promotion_status": PROMOTION_STATUS,
        "public_commit": public_commit,
        "code_sha256": code_sha256,
        "input_manifest": str(manifest_path.expanduser().resolve()),
        "input_manifest_sha256": sha256_file(
            manifest_path.expanduser().resolve()
        ),
        "certificate_count": len(certificates),
        "certificates": index_records,
        "boundaries": [
            "Validation receipts are authoritative; this builder does not "
            "recompute scientific gates.",
            "Certificates remain candidate current and are not CANONICAL.",
            "THE-121/THE-122 remain closed pending aggregate P5A-S, P5B, "
            "and P5C certification.",
        ],
    }
    if len(certificates) != 14:
        raise ContractError(
            f"internal certificate cardinality is {len(certificates)}, expected 14"
        )
    outputs[INDEX_NAME] = json_bytes(index)
    return outputs


def sha256_bytes(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def write_file_durable(path: Path, content: bytes) -> None:
    with path.open("xb") as stream:
        stream.write(content)
        stream.flush()
        os.fsync(stream.fileno())


def verify_existing_output(
    output_directory: Path,
    outputs: dict[str, bytes],
) -> None:
    if not output_directory.is_dir() or output_directory.is_symlink():
        raise ContractError(
            f"output path exists but is not a real directory: {output_directory}"
        )
    observed = {path.name for path in output_directory.iterdir()}
    expected = set(outputs)
    if observed != expected:
        raise ContractError(
            "existing output inventory mismatch: "
            f"missing={sorted(expected - observed)} "
            f"extra={sorted(observed - expected)}"
        )
    for name, expected_bytes in outputs.items():
        path = output_directory / name
        if (
            not path.is_file()
            or path.is_symlink()
            or path.read_bytes() != expected_bytes
        ):
            raise ContractError(f"existing output differs: {path}")


def publish_outputs(
    output_directory: Path,
    outputs: dict[str, bytes],
) -> Path:
    raw_output = output_directory.expanduser()
    if raw_output.is_symlink():
        raise ContractError(f"output directory cannot be a symlink: {raw_output}")
    output_root = raw_output.resolve()
    if output_root == output_root.parent:
        raise ContractError("output directory cannot be a filesystem root")
    if output_root.exists():
        verify_existing_output(output_root, outputs)
        return output_root / INDEX_NAME

    output_root.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(
        tempfile.mkdtemp(
            prefix=f".{output_root.name}.",
            suffix=".tmp",
            dir=output_root.parent,
        )
    )
    try:
        for name, content in outputs.items():
            write_file_durable(staging / name, content)
        os.replace(staging, output_root)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    return output_root / INDEX_NAME


def main() -> int:
    args = parse_args()
    try:
        outputs = build_outputs(args.input_manifest, args.output_directory)
        index_path = publish_outputs(args.output_directory, outputs)
    except (ContractError, OSError) as exc:
        print(
            f"THE-134 source-complete replay certificate build failed: {exc}",
            file=os.sys.stderr,
        )
        return 2
    print(index_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
