#!/usr/bin/env python3
"""Materialize and revalidate one digest-named THE-134 runtime bundle.

The companion ``build_the134_immutable_bundle_receipt.py`` deliberately
inventories candidate artifacts without copying them or claiming that a
digest-named bundle exists.  This tool crosses that boundary explicitly:

* consume one SHA-pinned, rehashed inventory receipt;
* copy every required artifact into the receipt-derived digest path;
* replace symlinks with their exact regular-file payloads;
* emit a resolver receipt whose artifact paths are inside that digest root;
* make the materialized tree read-only; and
* rehash every copied artifact before reporting success.

The materialization remains candidate execution evidence.  It cannot grant
scientific authority, submit jobs, or overwrite an existing bundle.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import shutil
import stat
import sys
import tempfile
from pathlib import Path
from typing import Any, Mapping


HERE = Path(__file__).resolve().parent
BUILDER_PATH = HERE / "build_the134_immutable_bundle_receipt.py"
SPEC = importlib.util.spec_from_file_location(
    "the134_immutable_bundle_receipt_builder", BUILDER_PATH
)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"cannot load immutable-bundle builder: {BUILDER_PATH}")
builder = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(builder)


MATERIALIZATION_SCHEMA = "THE134_IMMUTABLE_BUILD_BUNDLE_MATERIALIZATION_V1"
READBACK_SCHEMA = (
    "THE134_IMMUTABLE_BUILD_BUNDLE_MATERIALIZATION_READBACK_V1"
)
AUTHORITY_STATE = builder.AUTHORITY_STATE
RESOLVER_RECEIPT_RELATIVE = Path("metadata/resolver_bundle_receipt.json")
MATERIALIZATION_RECEIPT_RELATIVE = Path(
    "metadata/materialization_receipt.json"
)


class MaterializationError(ValueError):
    """Fail-closed materialization or readback violation."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


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


def require_sha256(label: str, value: object) -> str:
    try:
        return builder.require_sha256(label, value)
    except builder.BundleError as exc:
        raise MaterializationError(str(exc)) from exc


def require_absolute_file(label: str, value: object) -> Path:
    try:
        return builder.require_absolute_file(label, value)
    except builder.BundleError as exc:
        raise MaterializationError(str(exc)) from exc


def load_pinned_json(
    path: Path, expected_sha256: str, *, expected_schema: str
) -> dict[str, Any]:
    try:
        return builder.load_pinned_json(
            path,
            expected_sha256,
            expected_schema=expected_schema,
        )
    except builder.BundleError as exc:
        raise MaterializationError(str(exc)) from exc


def target_from_receipt(receipt: Mapping[str, Any]) -> Path:
    target = Path(str(receipt["digest_named_bundle_path"]))
    parent = Path(str(receipt["bundle_parent"]))
    expected_name = str(receipt["bundle_name"])
    if not target.is_absolute() or not parent.is_absolute():
        raise MaterializationError("bundle target and parent must be absolute")
    if target.parent != parent or target.name != expected_name:
        raise MaterializationError(
            "digest bundle target differs from its receipt-derived parent/name"
        )
    if target == Path("/"):
        raise MaterializationError("bundle target cannot be the filesystem root")
    return target


def artifact_relative_path(artifact: Mapping[str, Any]) -> Path:
    role = str(artifact["role"])
    basename = Path(str(artifact["path"])).name
    if not basename or basename in {".", ".."}:
        raise MaterializationError(
            f"artifact {role} has an invalid source basename"
        )
    return Path("artifacts") / role / basename


def copy_artifacts(
    receipt: Mapping[str, Any],
    temporary_root: Path,
    final_root: Path,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    resolver_artifacts: list[dict[str, Any]] = []
    inventory: list[dict[str, Any]] = []
    for artifact in receipt["artifacts"]:
        role = str(artifact["role"])
        relative = artifact_relative_path(artifact)
        source = require_absolute_file(
            f"source artifact {role}", artifact["path"]
        )
        destination = temporary_root / relative
        final_path = final_root / relative
        destination.parent.mkdir(parents=True, exist_ok=False)
        shutil.copyfile(source.resolve(strict=True), destination)

        expected_sha = require_sha256(
            f"source artifact {role}", artifact["sha256"]
        )
        observed_sha = sha256_file(destination)
        expected_size = int(artifact["size_bytes"])
        observed_size = destination.stat().st_size
        if observed_sha != expected_sha or observed_size != expected_size:
            raise MaterializationError(
                f"copied artifact differs: role={role} "
                f"expected_sha={expected_sha} observed_sha={observed_sha} "
                f"expected_size={expected_size} observed_size={observed_size}"
            )
        executable = role in builder.EXECUTABLE_ROLES
        destination.chmod(0o555 if executable else 0o444)
        resolver_artifacts.append(
            {
                "role": role,
                "path": str(final_path),
                "resolved_path": str(final_path.resolve(strict=False)),
                "sha256": observed_sha,
                "size_bytes": observed_size,
            }
        )
        inventory.append(
            {
                "role": role,
                "source_path": str(source),
                "bundle_relative_path": str(relative),
                "bundle_path": str(final_path),
                "sha256": observed_sha,
                "size_bytes": observed_size,
                "executable": executable,
            }
        )
    return (
        sorted(resolver_artifacts, key=lambda record: record["role"]),
        sorted(inventory, key=lambda record: record["role"]),
    )


def make_resolver_receipt(
    source_receipt: Mapping[str, Any],
    resolver_artifacts: list[dict[str, Any]],
) -> dict[str, Any]:
    result = dict(source_receipt)
    result["artifacts"] = resolver_artifacts
    result["digest_named_bundle_materialized"] = False
    return result


def make_materialization_receipt(
    *,
    source_receipt_path: Path,
    source_receipt_sha256: str,
    resolver_receipt_path: Path,
    resolver_receipt_sha256: str,
    source_receipt: Mapping[str, Any],
    artifact_inventory: list[dict[str, Any]],
) -> dict[str, Any]:
    total_bytes = sum(
        int(record["size_bytes"]) for record in artifact_inventory
    )
    semantic_inventory = [
        {
            "role": record["role"],
            "bundle_relative_path": record["bundle_relative_path"],
            "sha256": record["sha256"],
            "size_bytes": record["size_bytes"],
            "executable": record["executable"],
        }
        for record in artifact_inventory
    ]
    return {
        "schema": MATERIALIZATION_SCHEMA,
        "status": "PASS",
        "authority_state": AUTHORITY_STATE,
        "public_commit": source_receipt["public_commit"],
        "bundle_identity_sha256": source_receipt[
            "bundle_identity_sha256"
        ],
        "bundle_name": source_receipt["bundle_name"],
        "bundle_parent": source_receipt["bundle_parent"],
        "digest_named_bundle_path": source_receipt[
            "digest_named_bundle_path"
        ],
        "digest_named_bundle_materialized": True,
        "source_inventory_receipt": str(source_receipt_path),
        "source_inventory_receipt_sha256": source_receipt_sha256,
        "resolver_bundle_receipt": str(resolver_receipt_path),
        "resolver_bundle_receipt_sha256": resolver_receipt_sha256,
        "artifact_count": len(artifact_inventory),
        "dependency_count": len(source_receipt["dependencies"]),
        "hash_binding_count": len(source_receipt["hash_bindings"]),
        "total_artifact_bytes": total_bytes,
        "artifact_inventory_sha256": canonical_sha256(semantic_inventory),
        "artifact_inventory": artifact_inventory,
        "readonly_contract": {
            "regular_artifacts_are_symlink_free": True,
            "executors_mode": "0555",
            "nonexecutors_mode": "0444",
            "directories_mode": "0555",
            "overwrite_existing_bundle": False,
        },
    }


def make_tree_readonly(root: Path) -> None:
    files = sorted(
        (path for path in root.rglob("*") if path.is_file()),
        key=lambda path: len(path.parts),
        reverse=True,
    )
    for path in files:
        relative = path.relative_to(root)
        executable = (
            len(relative.parts) >= 3
            and relative.parts[0] == "artifacts"
            and relative.parts[1] in builder.EXECUTABLE_ROLES
        )
        path.chmod(0o555 if executable else 0o444)
    directories = sorted(
        (path for path in root.rglob("*") if path.is_dir()),
        key=lambda path: len(path.parts),
        reverse=True,
    )
    for path in directories:
        path.chmod(0o555)
    root.chmod(0o555)


def atomic_write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f"{path.name}.tmp.{os.getpid()}")
    temporary.write_bytes(canonical_json_bytes(payload))
    os.replace(temporary, path)


def materialize(
    source_receipt_path: Path,
    expected_source_receipt_sha256: str,
) -> dict[str, Any]:
    expected_source_sha = require_sha256(
        "source inventory receipt", expected_source_receipt_sha256
    )
    try:
        source_receipt = builder.load_receipt(
            source_receipt_path, expected_source_sha
        )
    except builder.BundleError as exc:
        raise MaterializationError(str(exc)) from exc

    target = target_from_receipt(source_receipt)
    if target.exists() or target.is_symlink():
        raise MaterializationError(
            f"refusing to overwrite existing digest bundle: {target}"
        )
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(
        tempfile.mkdtemp(
            prefix=f".{target.name}.tmp.",
            dir=str(target.parent),
        )
    )
    renamed = False
    try:
        resolver_artifacts, artifact_inventory = copy_artifacts(
            source_receipt, temporary, target
        )
        resolver_receipt = make_resolver_receipt(
            source_receipt, resolver_artifacts
        )
        resolver_path_temporary = temporary / RESOLVER_RECEIPT_RELATIVE
        resolver_path_final = target / RESOLVER_RECEIPT_RELATIVE
        atomic_write_json(resolver_path_temporary, resolver_receipt)
        resolver_sha = sha256_file(resolver_path_temporary)

        materialization_path_temporary = (
            temporary / MATERIALIZATION_RECEIPT_RELATIVE
        )
        materialization_path_final = (
            target / MATERIALIZATION_RECEIPT_RELATIVE
        )
        materialization_receipt = make_materialization_receipt(
            source_receipt_path=source_receipt_path,
            source_receipt_sha256=expected_source_sha,
            resolver_receipt_path=resolver_path_final,
            resolver_receipt_sha256=resolver_sha,
            source_receipt=source_receipt,
            artifact_inventory=artifact_inventory,
        )
        atomic_write_json(
            materialization_path_temporary, materialization_receipt
        )
        make_tree_readonly(temporary)
        os.replace(temporary, target)
        renamed = True
        materialization_sha = sha256_file(materialization_path_final)
        validate_materialization_payload(
            materialization_receipt,
            materialization_receipt_path=materialization_path_final,
            rehash=True,
        )
        return {
            "schema": MATERIALIZATION_SCHEMA,
            "status": "PASS",
            "authority_state": AUTHORITY_STATE,
            "bundle_identity_sha256": source_receipt[
                "bundle_identity_sha256"
            ],
            "digest_named_bundle_path": str(target),
            "materialization_receipt": str(materialization_path_final),
            "materialization_receipt_sha256": materialization_sha,
            "resolver_bundle_receipt": str(resolver_path_final),
            "resolver_bundle_receipt_sha256": resolver_sha,
            "artifact_count": len(artifact_inventory),
            "total_artifact_bytes": sum(
                int(record["size_bytes"]) for record in artifact_inventory
            ),
        }
    finally:
        if not renamed and temporary.exists():
            shutil.rmtree(temporary)


def require_readonly_file(path: Path, *, executable: bool) -> None:
    if path.is_symlink() or not path.is_file():
        raise MaterializationError(
            f"materialized artifact is not a symlink-free regular file: {path}"
        )
    mode = stat.S_IMODE(path.stat().st_mode)
    expected = 0o555 if executable else 0o444
    if mode != expected:
        raise MaterializationError(
            f"materialized file mode differs: path={path} "
            f"expected={oct(expected)} observed={oct(mode)}"
        )


def validate_materialization_payload(
    payload: Mapping[str, Any],
    *,
    materialization_receipt_path: Path,
    rehash: bool,
) -> dict[str, Any]:
    required = {
        "schema",
        "status",
        "authority_state",
        "public_commit",
        "bundle_identity_sha256",
        "bundle_name",
        "bundle_parent",
        "digest_named_bundle_path",
        "digest_named_bundle_materialized",
        "source_inventory_receipt",
        "source_inventory_receipt_sha256",
        "resolver_bundle_receipt",
        "resolver_bundle_receipt_sha256",
        "artifact_count",
        "dependency_count",
        "hash_binding_count",
        "total_artifact_bytes",
        "artifact_inventory_sha256",
        "artifact_inventory",
        "readonly_contract",
    }
    if set(payload) != required:
        raise MaterializationError(
            "materialization receipt field inventory differs: "
            f"missing={sorted(required - set(payload))} "
            f"extra={sorted(set(payload) - required)}"
        )
    if (
        payload["schema"] != MATERIALIZATION_SCHEMA
        or payload["status"] != "PASS"
        or payload["authority_state"] != AUTHORITY_STATE
        or payload["digest_named_bundle_materialized"] is not True
    ):
        raise MaterializationError(
            "materialization schema/status/authority contract differs"
        )
    expected_readonly_contract = {
        "regular_artifacts_are_symlink_free": True,
        "executors_mode": "0555",
        "nonexecutors_mode": "0444",
        "directories_mode": "0555",
        "overwrite_existing_bundle": False,
    }
    if payload["readonly_contract"] != expected_readonly_contract:
        raise MaterializationError(
            "materialization readonly/overwrite contract differs"
        )
    identity = require_sha256(
        "bundle identity", payload["bundle_identity_sha256"]
    )
    bundle_name = f"{builder.BUNDLE_NAME_PREFIX}{identity}"
    target = Path(str(payload["digest_named_bundle_path"]))
    if (
        payload["bundle_name"] != bundle_name
        or target.name != bundle_name
        or target.parent != Path(str(payload["bundle_parent"]))
        or not target.is_dir()
        or target.is_symlink()
    ):
        raise MaterializationError(
            "materialized digest path/name/parent contract differs"
        )
    expected_receipt_path = target / MATERIALIZATION_RECEIPT_RELATIVE
    if materialization_receipt_path != expected_receipt_path:
        raise MaterializationError(
            "materialization receipt is outside its digest bundle"
        )
    resolver_path = require_absolute_file(
        "resolver bundle receipt", payload["resolver_bundle_receipt"]
    )
    expected_resolver_path = target / RESOLVER_RECEIPT_RELATIVE
    if resolver_path != expected_resolver_path:
        raise MaterializationError(
            "resolver receipt is outside its digest bundle"
        )
    resolver_sha = require_sha256(
        "resolver bundle receipt", payload["resolver_bundle_receipt_sha256"]
    )
    if sha256_file(resolver_path) != resolver_sha:
        raise MaterializationError("resolver bundle receipt hash differs")
    try:
        resolver = builder.load_receipt(resolver_path, resolver_sha)
    except builder.BundleError as exc:
        raise MaterializationError(str(exc)) from exc
    if (
        resolver["bundle_identity_sha256"] != identity
        or resolver["public_commit"] != payload["public_commit"]
    ):
        raise MaterializationError(
            "resolver receipt identity/public commit differs"
        )
    raw_inventory = payload["artifact_inventory"]
    if not isinstance(raw_inventory, list):
        raise MaterializationError("artifact_inventory must be a list")
    by_role = {
        str(record.get("role", "")): record
        for record in raw_inventory
        if isinstance(record, dict)
    }
    if (
        len(raw_inventory) != len(by_role)
        or set(by_role) != set(builder.REQUIRED_ARTIFACT_ROLES)
        or int(payload["artifact_count"]) != len(by_role)
    ):
        raise MaterializationError(
            "materialized artifact role/count closure differs"
        )
    semantic_inventory: list[dict[str, Any]] = []
    total_bytes = 0
    resolver_by_role = {
        artifact["role"]: artifact for artifact in resolver["artifacts"]
    }
    for role in sorted(by_role):
        record = by_role[role]
        expected_keys = {
            "role",
            "source_path",
            "bundle_relative_path",
            "bundle_path",
            "sha256",
            "size_bytes",
            "executable",
        }
        if set(record) != expected_keys:
            raise MaterializationError(
                f"artifact inventory fields differ: role={role}"
            )
        relative = Path(str(record["bundle_relative_path"]))
        expected_path = target / relative
        path = Path(str(record["bundle_path"]))
        if (
            path != expected_path
            or path != Path(str(resolver_by_role[role]["path"]))
            or relative.is_absolute()
            or ".." in relative.parts
        ):
            raise MaterializationError(
                f"artifact path containment differs: role={role}"
            )
        expected_sha = require_sha256(
            f"materialized artifact {role}", record["sha256"]
        )
        expected_size = int(record["size_bytes"])
        executable = bool(record["executable"])
        if executable != (role in builder.EXECUTABLE_ROLES):
            raise MaterializationError(
                f"artifact executable classification differs: role={role}"
            )
        if rehash:
            require_readonly_file(path, executable=executable)
            if (
                sha256_file(path) != expected_sha
                or path.stat().st_size != expected_size
            ):
                raise MaterializationError(
                    f"materialized artifact content differs: role={role}"
                )
        semantic_inventory.append(
            {
                "role": role,
                "bundle_relative_path": str(relative),
                "sha256": expected_sha,
                "size_bytes": expected_size,
                "executable": executable,
            }
        )
        total_bytes += expected_size
    if (
        require_sha256(
            "artifact inventory", payload["artifact_inventory_sha256"]
        )
        != canonical_sha256(semantic_inventory)
        or int(payload["total_artifact_bytes"]) != total_bytes
        or int(payload["dependency_count"]) != len(resolver["dependencies"])
        or int(payload["hash_binding_count"]) != len(
            resolver["hash_bindings"]
        )
    ):
        raise MaterializationError(
            "materialization aggregate inventory differs"
        )
    if rehash:
        require_readonly_file(
            materialization_receipt_path, executable=False
        )
        require_readonly_file(resolver_path, executable=False)
        for directory in [target, *target.rglob("*")]:
            if directory.is_dir():
                mode = stat.S_IMODE(directory.stat().st_mode)
                if mode != 0o555:
                    raise MaterializationError(
                        f"materialized directory is writable: {directory}"
                    )
    return dict(payload)


def verify(
    materialization_receipt_path: Path,
    expected_materialization_receipt_sha256: str,
) -> dict[str, Any]:
    expected_sha = require_sha256(
        "materialization receipt",
        expected_materialization_receipt_sha256,
    )
    payload = load_pinned_json(
        materialization_receipt_path,
        expected_sha,
        expected_schema=MATERIALIZATION_SCHEMA,
    )
    validated = validate_materialization_payload(
        payload,
        materialization_receipt_path=materialization_receipt_path,
        rehash=True,
    )
    return {
        "schema": READBACK_SCHEMA,
        "status": "PASS",
        "authority_state": AUTHORITY_STATE,
        "materialization_receipt": str(materialization_receipt_path),
        "materialization_receipt_sha256": expected_sha,
        "bundle_identity_sha256": validated["bundle_identity_sha256"],
        "digest_named_bundle_path": validated[
            "digest_named_bundle_path"
        ],
        "resolver_bundle_receipt": validated["resolver_bundle_receipt"],
        "resolver_bundle_receipt_sha256": validated[
            "resolver_bundle_receipt_sha256"
        ],
        "artifacts_rehashed": validated["artifact_count"],
        "total_artifact_bytes": validated["total_artifact_bytes"],
        "readback": "PASS_SYMLINK_FREE_READONLY_CONTENT_EXACT",
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    materialize_parser = subparsers.add_parser(
        "materialize", help="copy one pinned inventory into its digest root"
    )
    materialize_parser.add_argument("--receipt", type=Path, required=True)
    materialize_parser.add_argument(
        "--expected-receipt-sha256", required=True
    )

    verify_parser = subparsers.add_parser(
        "verify", help="rehash and revalidate one materialized bundle"
    )
    verify_parser.add_argument("--receipt", type=Path, required=True)
    verify_parser.add_argument(
        "--expected-receipt-sha256", required=True
    )
    verify_parser.add_argument("--output-json", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        if args.command == "materialize":
            summary = materialize(
                args.receipt, args.expected_receipt_sha256
            )
        else:
            summary = verify(
                args.receipt, args.expected_receipt_sha256
            )
            if args.output_json is not None:
                atomic_write_json(args.output_json, summary)
        print(json.dumps(summary, sort_keys=True))
        return 0
    except (MaterializationError, OSError, ValueError) as exc:
        print(f"THE134_BUNDLE_MATERIALIZATION_ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
