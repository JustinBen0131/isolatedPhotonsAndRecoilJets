#!/usr/bin/env python3
"""Validate the six-lane current RecoilJets continuity projection.

The continuity manifest is a read-only projection of the existing artifact
registry.  It does not register, copy, promote, or delete artifacts.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import sys
from typing import Any


EXPECTED_SAMPLE_KEYS = {
    "pp_data_merged",
    "pp_sim_photonjet_merged",
    "pp_sim_inclusivejet_merged",
    "auau_data_merged",
    "auau_sim_photonjet_merged",
    "auau_sim_inclusivejet_merged",
}
ALLOWED_CLASSIFICATIONS = {
    "CERTIFIED_REFERENCE",
    "QUALIFIED_CURRENT",
    "PROVISIONAL",
    "PARTIAL",
    "INVALID",
}
ALLOWED_ROOT_HEALTH = {"PASS_FRESH_TFILE", "PASS_COMPOSITE"}
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


def read_object(path: Path, label: str) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text())
    except FileNotFoundError as exc:
        raise ValueError(f"{label} is missing: {path}") from exc
    except json.JSONDecodeError as exc:
        raise ValueError(f"{label} is invalid JSON: {path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise ValueError(f"{label} must contain a JSON object: {path}")
    return payload


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require_sha256(value: Any, label: str) -> str:
    if not isinstance(value, str) or not SHA256_RE.fullmatch(value):
        raise ValueError(f"{label} is not a lowercase SHA-256")
    return value


def inventory(path: Path) -> dict[str, Any]:
    try:
        import uproot
    except ImportError as exc:
        raise ValueError(
            "uproot is required for full continuity inventory validation"
        ) from exc
    with uproot.open(path) as root_file:
        classnames = root_file.classnames(recursive=True, cycle=True)
    rows = sorted(f"{key}\t{class_name}" for key, class_name in classnames.items())
    return {
        "objects": len(rows),
        "histograms": sum(
            class_name.startswith(("TH1", "TH2", "TH3", "TProfile"))
            for class_name in classnames.values()
        ),
        "trees": sum(
            class_name.startswith("TTree") for class_name in classnames.values()
        ),
        "directories": sum(
            class_name.startswith("TDirectory")
            for class_name in classnames.values()
        ),
        "sha256": hashlib.sha256(("\n".join(rows) + "\n").encode()).hexdigest(),
    }


def validate(
    registry_path: Path,
    manifest_path: Path,
    *,
    metadata_only: bool = False,
) -> dict[str, Any]:
    registry = read_object(registry_path, "artifact registry")
    manifest = read_object(manifest_path, "continuity manifest")
    failures: list[str] = []

    if manifest.get("schema") != "RecoilJetsCurrentContinuityFoundationV1":
        failures.append("wrong continuity schema")
    if manifest.get("rerun_performed") is not False:
        failures.append("continuity projection records a rerun")
    if manifest.get("copy_performed") is not False:
        failures.append("continuity projection records a copy")
    if manifest.get("promotion_performed") is not False:
        failures.append("continuity projection records a promotion")
    if manifest.get("final_h70_authority") is not False:
        failures.append("current continuity data is mislabeled as final H70")

    registry_ref = manifest.get("registry", {})
    if not isinstance(registry_ref, dict):
        failures.append("registry reference is missing")
    else:
        expected_registry_sha = registry_ref.get("sha256")
        if require_sha256(expected_registry_sha, "registry SHA-256") != sha256_file(
            registry_path
        ):
            failures.append("artifact registry hash drift")

    current = registry.get("current")
    registry_artifacts = registry.get("artifacts")
    if not isinstance(current, dict) or set(current) != EXPECTED_SAMPLE_KEYS:
        failures.append("registry current keys are not the exact six-lane set")
        current = current if isinstance(current, dict) else {}
    if not isinstance(registry_artifacts, list):
        failures.append("registry artifacts is not a list")
        registry_artifacts = []
    artifact_by_id = {
        item.get("id"): item
        for item in registry_artifacts
        if isinstance(item, dict) and isinstance(item.get("id"), str)
    }

    manifest_artifacts = manifest.get("artifacts")
    if not isinstance(manifest_artifacts, list):
        failures.append("manifest artifacts is not a list")
        manifest_artifacts = []
    manifest_by_key = {
        item.get("sample_key"): item
        for item in manifest_artifacts
        if isinstance(item, dict) and isinstance(item.get("sample_key"), str)
    }
    if set(manifest_by_key) != EXPECTED_SAMPLE_KEYS:
        failures.append("manifest sample keys are not the exact six-lane set")

    lane_receipts: list[dict[str, Any]] = []
    for sample_key in sorted(EXPECTED_SAMPLE_KEYS):
        item = manifest_by_key.get(sample_key)
        if item is None:
            continue
        artifact_id = item.get("artifact_id")
        registry_id = current.get(sample_key)
        if artifact_id != registry_id:
            failures.append(
                f"{sample_key}: manifest artifact ID does not match registry current"
            )
        registry_item = artifact_by_id.get(artifact_id)
        if registry_item is None:
            failures.append(f"{sample_key}: artifact ID is absent from registry")
            continue
        if registry_item.get("sample_key") != sample_key:
            failures.append(f"{sample_key}: registry sample key mismatch")

        classification = item.get("classification")
        if classification not in ALLOWED_CLASSIFICATIONS:
            failures.append(f"{sample_key}: invalid continuity classification")
        if not item.get("known_omissions"):
            failures.append(f"{sample_key}: known omissions are empty")
        if not item.get("consumers"):
            failures.append(f"{sample_key}: consumer list is empty")
        if not item.get("source_coverage"):
            failures.append(f"{sample_key}: source coverage is empty")
        if not item.get("semantic_authority"):
            failures.append(f"{sample_key}: semantic authority is empty")

        root_health = item.get("root_health")
        if not isinstance(root_health, dict):
            failures.append(f"{sample_key}: ROOT health receipt is missing")
        elif root_health.get("status") not in ALLOWED_ROOT_HEALTH:
            failures.append(f"{sample_key}: ROOT health status is not accepted")
        elif root_health["status"] == "PASS_FRESH_TFILE":
            if root_health.get("zombie") is not False:
                failures.append(f"{sample_key}: fresh TFile is zombie or unproven")
            if root_health.get("recovered") is not False:
                failures.append(f"{sample_key}: fresh TFile is recovered or unproven")
            if not isinstance(root_health.get("top_level_keys"), int) or (
                root_health["top_level_keys"] <= 0
            ):
                failures.append(f"{sample_key}: fresh TFile key count is invalid")
        else:
            arithmetic = root_health.get("prior_full_arithmetic_audit")
            if root_health.get("fresh_uproot_open") is not True:
                failures.append(f"{sample_key}: composite uproot health is unproven")
            if not isinstance(arithmetic, dict) or arithmetic.get("status") != "PASS":
                failures.append(
                    f"{sample_key}: composite arithmetic certificate is unproven"
                )
            elif (
                arithmetic.get("max_content_delta") != 0.0
                or arithmetic.get("max_sumw2_delta") != 0.0
            ):
                failures.append(
                    f"{sample_key}: composite arithmetic certificate is nonexact"
                )

        local_path = Path(str(item.get("local_path", ""))).expanduser()
        registry_roots = registry_item.get("root_paths", [])
        if str(local_path) not in registry_roots:
            failures.append(f"{sample_key}: local path is not registry-bound")
        remote_path = item.get("remote_path")
        if remote_path not in registry_item.get("remote_paths", []):
            failures.append(f"{sample_key}: remote path is not registry-bound")

        expected_sha = require_sha256(item.get("sha256"), f"{sample_key} SHA-256")
        expected_bytes = item.get("bytes")
        if not isinstance(expected_bytes, int) or expected_bytes <= 0:
            failures.append(f"{sample_key}: invalid byte size")

        observed_inventory: dict[str, Any] | None = None
        if not metadata_only:
            if not local_path.is_file():
                failures.append(f"{sample_key}: local ROOT is missing")
            else:
                observed_bytes = local_path.stat().st_size
                if observed_bytes != expected_bytes:
                    failures.append(f"{sample_key}: byte-size drift")
                if sha256_file(local_path) != expected_sha:
                    failures.append(f"{sample_key}: content hash drift")
                try:
                    observed_inventory = inventory(local_path)
                except Exception as exc:  # preserve exact lane context
                    failures.append(f"{sample_key}: inventory read failed: {exc}")
                expected_inventory = item.get("semantic_inventory")
                if (
                    observed_inventory is not None
                    and isinstance(expected_inventory, dict)
                ):
                    for field in (
                        "objects",
                        "histograms",
                        "trees",
                        "directories",
                        "sha256",
                    ):
                        if observed_inventory.get(field) != expected_inventory.get(
                            field
                        ):
                            failures.append(
                                f"{sample_key}: inventory drift in {field}"
                            )
                elif not isinstance(expected_inventory, dict):
                    failures.append(f"{sample_key}: semantic inventory is missing")

        lane_receipts.append(
            {
                "sample_key": sample_key,
                "artifact_id": artifact_id,
                "classification": classification,
                "bytes": expected_bytes,
                "sha256": expected_sha,
                "inventory": observed_inventory,
            }
        )

    return {
        "schema": "RecoilJetsCurrentContinuityValidationV1",
        "status": "PASS" if not failures else "FAIL",
        "registry": str(registry_path),
        "manifest": str(manifest_path),
        "metadata_only": metadata_only,
        "expected_lanes": len(EXPECTED_SAMPLE_KEYS),
        "validated_lanes": len(lane_receipts),
        "lanes": lane_receipts,
        "failures": failures,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument(
        "--metadata-only",
        action="store_true",
        help="validate bindings without reading or hashing ROOT payloads",
    )
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        receipt = validate(
            args.registry.resolve(),
            args.manifest.resolve(),
            metadata_only=args.metadata_only,
        )
    except ValueError as exc:
        receipt = {
            "schema": "RecoilJetsCurrentContinuityValidationV1",
            "status": "FAIL",
            "failures": [str(exc)],
        }
    encoded = json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded)
    sys.stdout.write(encoded)
    return 0 if receipt["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
