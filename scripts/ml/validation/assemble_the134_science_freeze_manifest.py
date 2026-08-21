#!/usr/bin/env python3
"""Assemble the exact seven-view THE-134 science-freeze input manifest.

This tool is deliberately a receipt assembler, not a scientific validator.
It accepts one explicit paired-registry path per shower view and one explicit
replay-certificate path per system/view pair, verifies their top-level
identities, computes their byte hashes, and emits the exact manifest consumed
by ``build_the134_science_freeze_certificate.py``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import tempfile
from pathlib import Path

from the134_aggregate_contract import (
    PAIRED_REGISTRY_SCHEMA,
    REPLAY_CERTIFICATE_SCHEMA,
    REQUIRED_SYSTEM_VIEW_KEYS,
    REQUIRED_VIEW_KEYS,
    SCIENCE_MANIFEST_SCHEMA,
    ContractError,
    require_git_commit,
    require_sha256,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--public-commit", required=True)
    parser.add_argument("--code-sha256", required=True)
    parser.add_argument(
        "--registry",
        action="append",
        default=[],
        metavar="VIEW=PATH",
        help="Repeat exactly once for each required shower view.",
    )
    parser.add_argument(
        "--replay-certificate",
        action="append",
        default=[],
        metavar="SYSTEM:VIEW=PATH",
        help="Repeat exactly once for each required system/view pair.",
    )
    parser.add_argument("--json-out", required=True, type=Path)
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path, label: str) -> dict:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ContractError(f"{label} is not readable JSON: {path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise ContractError(f"{label} must be a JSON object: {path}")
    return payload


def parse_binding(value: str, label: str) -> tuple[str, Path]:
    key, separator, raw_path = value.partition("=")
    if not separator or not key or not raw_path:
        raise ContractError(f"{label} must use KEY=PATH syntax: {value!r}")
    path = Path(raw_path).expanduser().resolve(strict=True)
    if not path.is_file():
        raise ContractError(f"{label} is not a regular file: {path}")
    return key, path


def exact_inventory(
    values: list[str],
    expected: tuple[str, ...],
    label: str,
) -> dict[str, Path]:
    observed: dict[str, Path] = {}
    for value in values:
        key, path = parse_binding(value, label)
        if key in observed:
            raise ContractError(f"{label} contains duplicate key={key!r}")
        observed[key] = path
    missing = set(expected) - set(observed)
    extra = set(observed) - set(expected)
    if missing or extra:
        raise ContractError(
            f"{label} inventory mismatch: missing={sorted(missing)} "
            f"extra={sorted(extra)}"
        )
    return observed


def validate_registry(
    view: str,
    path: Path,
    public_commit: str,
    code_sha256: str,
) -> None:
    payload = read_json(path, f"{view} paired registry")
    expected = {
        "schema": PAIRED_REGISTRY_SCHEMA,
        "status": "READY",
        "promotion_status": "CANDIDATE_CURRENT_NOT_CANONICAL",
        "shower_definition": view,
        "public_commit": public_commit,
        "code_sha256": code_sha256,
    }
    failures = [
        field for field, wanted in expected.items() if payload.get(field) != wanted
    ]
    if failures:
        raise ContractError(
            f"{view} paired registry top-level identity mismatch: {failures}"
        )


def validate_replay_certificate(
    system: str,
    view: str,
    path: Path,
    public_commit: str,
    code_sha256: str,
) -> None:
    payload = read_json(path, f"{system}/{view} replay certificate")
    expected = {
        "schema": REPLAY_CERTIFICATE_SCHEMA,
        "status": "PASS",
        "promotion_status": "CANDIDATE_CURRENT_NOT_CANONICAL",
        "system": system,
        "shower_definition": view,
        "public_commit": public_commit,
        "code_sha256": code_sha256,
    }
    failures = [
        field for field, wanted in expected.items() if payload.get(field) != wanted
    ]
    if failures:
        raise ContractError(
            f"{system}/{view} replay certificate top-level identity mismatch: "
            f"{failures}"
        )


def write_json_atomic(path: Path, payload: dict) -> None:
    path = path.expanduser().resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    temporary_path = Path(temporary)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
            json.dump(payload, stream, indent=2, sort_keys=True)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary_path, path)
    finally:
        temporary_path.unlink(missing_ok=True)


def build_manifest(args: argparse.Namespace) -> dict:
    public_commit = require_git_commit(args.public_commit)
    code_sha256 = require_sha256(args.code_sha256, "code_sha256")
    registries = exact_inventory(
        args.registry,
        tuple(REQUIRED_VIEW_KEYS),
        "paired registry",
    )
    replay_expected = tuple(
        f"{system}:{view}" for system, view in REQUIRED_SYSTEM_VIEW_KEYS
    )
    replay_certificates = exact_inventory(
        args.replay_certificate,
        replay_expected,
        "replay certificate",
    )

    registry_records = []
    for view in REQUIRED_VIEW_KEYS:
        path = registries[view]
        validate_registry(view, path, public_commit, code_sha256)
        registry_records.append(
            {"view": view, "path": str(path), "sha256": sha256_file(path)}
        )

    replay_records = []
    for system, view in REQUIRED_SYSTEM_VIEW_KEYS:
        path = replay_certificates[f"{system}:{view}"]
        validate_replay_certificate(
            system,
            view,
            path,
            public_commit,
            code_sha256,
        )
        replay_records.append(
            {
                "system": system,
                "view": view,
                "path": str(path),
                "sha256": sha256_file(path),
            }
        )

    return {
        "schema": SCIENCE_MANIFEST_SCHEMA,
        "public_commit": public_commit,
        "code_sha256": code_sha256,
        "paired_registries": registry_records,
        "replay_certificates": replay_records,
    }


def main() -> int:
    args = parse_args()
    try:
        payload = build_manifest(args)
        write_json_atomic(args.json_out, payload)
    except (ContractError, OSError) as exc:
        print(f"THE-134 science-freeze manifest assembly failed: {exc}", file=sys.stderr)
        return 2
    print(args.json_out.expanduser().resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
