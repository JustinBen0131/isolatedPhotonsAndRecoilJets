#!/usr/bin/env python3
"""Project THE-134 foundation storage with two merge workspaces and headroom."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
from the134_aggregate_contract import (  # noqa: E402
    MINIMUM_FREE_FRACTION,
    REQUIRED_SOURCE_KEYS,
    SCIENCE_CERTIFICATE_GATES,
    SCIENCE_CERTIFICATE_KEYS,
    SCIENCE_CERTIFICATE_SCHEMA,
    SIMULTANEOUS_MERGE_WORKSPACE_COUNT,
    STORAGE_CERTIFICATE_SCHEMA,
    STORAGE_MANIFEST_SCHEMA,
    ContractError,
    require_exact_key_inventory,
    require_exact_true_gates,
    require_sha256,
    read_json,
    semantic_receipt,
    verify_registered_json,
    verify_semantic_receipt,
)


BYTE_FIELDS = (
    "projected_direct_reference_bytes",
    "projected_replay_ttree_bytes",
    "projected_final_merged_bytes",
    "projected_evidence_bytes",
    "merge_workspace_bytes",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--json-out", type=Path, required=True)
    return parser.parse_args()


def require_nonnegative_integer(value: object, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ContractError(f"{label} must be a nonnegative integer")
    return value


def require_positive_integer(value: object, label: str) -> int:
    result = require_nonnegative_integer(value, label)
    if result <= 0:
        raise ContractError(f"{label} must be greater than zero")
    return result


def validate_science_certificate(record: dict) -> dict:
    path, payload, observed_sha = verify_registered_json(
        record, "THE-134 science-freeze certificate"
    )
    observed_keys = set(payload)
    expected_keys = set(SCIENCE_CERTIFICATE_KEYS)
    if observed_keys != expected_keys:
        raise ContractError(
            "science-freeze certificate key inventory mismatch: "
            f"missing={sorted(expected_keys - observed_keys)} "
            f"extra={sorted(observed_keys - expected_keys)}"
        )
    if payload.get("schema") != SCIENCE_CERTIFICATE_SCHEMA:
        raise ContractError("science-freeze certificate schema mismatch")
    if payload.get("status") != "PASS":
        raise ContractError("science-freeze certificate is not PASS")
    if payload.get("gate") != "P5A-S":
        raise ContractError("science-freeze gate identity mismatch")
    if payload.get("promotion_status") != "CANDIDATE_CURRENT_NOT_CANONICAL":
        raise ContractError("science-freeze promotion boundary mismatch")
    require_exact_true_gates(
        payload.get("gates"),
        SCIENCE_CERTIFICATE_GATES,
        "science-freeze certificate",
    )
    aggregate_semantic_sha256 = verify_semantic_receipt(
        payload,
        "aggregate_semantic_sha256",
        "science-freeze certificate",
    )
    if payload.get("training_source_count") != len(REQUIRED_SOURCE_KEYS):
        raise ContractError("science-freeze training-source count mismatch")
    source_bindings = require_exact_key_inventory(
        payload.get("training_source_bindings"),
        REQUIRED_SOURCE_KEYS,
        lambda item: (item.get("system"), item.get("source")),
        "science-freeze training-source bindings",
    )
    for system, source in REQUIRED_SOURCE_KEYS:
        binding = source_bindings[(system, source)]
        if set(binding) != {
            "system",
            "source",
            "source_manifest_sha256",
            "direct_artifact_sha256",
            "writer_artifact_sha256",
        }:
            raise ContractError(
                f"science-freeze {system}/{source} binding fields differ"
            )
        for field in (
            "source_manifest_sha256",
            "direct_artifact_sha256",
            "writer_artifact_sha256",
        ):
            require_sha256(
                binding.get(field),
                f"science-freeze {system}/{source}.{field}",
            )
    return {
        "path": str(path),
        "sha256": observed_sha,
        "aggregate_semantic_sha256": aggregate_semantic_sha256,
    }


def build_projection(manifest_path: Path) -> dict:
    manifest = read_json(manifest_path, "storage projection manifest")
    if manifest.get("schema") != STORAGE_MANIFEST_SCHEMA:
        raise ContractError(
            f"manifest.schema={manifest.get('schema')!r}, "
            f"expected={STORAGE_MANIFEST_SCHEMA!r}"
        )
    science_certificate = validate_science_certificate(
        manifest.get("science_freeze_certificate", {})
    )
    quota = manifest.get("quota")
    if not isinstance(quota, dict):
        raise ContractError("quota must be an object")
    quota_bytes = require_positive_integer(quota.get("quota_bytes"), "quota_bytes")
    used_bytes = require_nonnegative_integer(quota.get("used_bytes"), "used_bytes")
    if used_bytes > quota_bytes:
        raise ContractError("used_bytes exceeds quota_bytes")
    usage_authority_sha256 = require_sha256(
        quota.get("usage_authority_sha256"), "usage_authority_sha256"
    )
    minimum_free_fraction = quota.get(
        "minimum_free_fraction", MINIMUM_FREE_FRACTION
    )
    if (
        isinstance(minimum_free_fraction, bool)
        or not isinstance(minimum_free_fraction, (int, float))
        or float(minimum_free_fraction) < MINIMUM_FREE_FRACTION
        or float(minimum_free_fraction) >= 1.0
    ):
        raise ContractError(
            "minimum_free_fraction must be at least 0.20 and below 1.0"
        )
    minimum_free_fraction = float(minimum_free_fraction)
    merge_workspace_count = manifest.get("simultaneous_merge_workspace_count")
    if merge_workspace_count != SIMULTANEOUS_MERGE_WORKSPACE_COUNT:
        raise ContractError(
            "simultaneous_merge_workspace_count must equal exactly 2"
        )

    projections = require_exact_key_inventory(
        manifest.get("source_projections"),
        REQUIRED_SOURCE_KEYS,
        lambda record: (record.get("system"), record.get("source")),
        "source projections",
    )
    normalized = []
    totals = {field: 0 for field in BYTE_FIELDS}
    for system, source in REQUIRED_SOURCE_KEYS:
        record = projections[(system, source)]
        measurement_receipt_sha256 = require_sha256(
            record.get("measurement_receipt_sha256"),
            f"{system}/{source}.measurement_receipt_sha256",
        )
        values = {
            field: require_positive_integer(
                record.get(field), f"{system}/{source}.{field}"
            )
            for field in BYTE_FIELDS
        }
        for field, value in values.items():
            totals[field] += value
        normalized.append(
            {
                "system": system,
                "source": source,
                "measurement_receipt_sha256": measurement_receipt_sha256,
                **values,
            }
        )

    largest_merge_workspace_bytes = max(
        item["merge_workspace_bytes"] for item in normalized
    )
    merge_workspace_reserve_bytes = (
        SIMULTANEOUS_MERGE_WORKSPACE_COUNT * largest_merge_workspace_bytes
    )
    persistent_increment_bytes = (
        totals["projected_direct_reference_bytes"]
        + totals["projected_replay_ttree_bytes"]
        + totals["projected_final_merged_bytes"]
        + totals["projected_evidence_bytes"]
    )
    projected_peak_increment_bytes = (
        persistent_increment_bytes + merge_workspace_reserve_bytes
    )
    projected_peak_used_bytes = used_bytes + projected_peak_increment_bytes
    projected_peak_free_bytes = quota_bytes - projected_peak_used_bytes
    projected_peak_free_fraction = projected_peak_free_bytes / quota_bytes

    local_hot = manifest.get("local_hot_tier")
    if not isinstance(local_hot, dict):
        raise ContractError("local_hot_tier must be an object")
    local_free_bytes = require_nonnegative_integer(
        local_hot.get("free_bytes"), "local_hot_tier.free_bytes"
    )
    local_protected_free_bytes = require_nonnegative_integer(
        local_hot.get("protected_free_bytes"),
        "local_hot_tier.protected_free_bytes",
    )
    planned_local_pull_bytes = require_nonnegative_integer(
        local_hot.get("planned_pull_bytes"),
        "local_hot_tier.planned_pull_bytes",
    )
    local_free_after_pull_bytes = local_free_bytes - planned_local_pull_bytes

    gates = {
        "science_freeze_certificate_pass": True,
        "all_required_source_projections_present": True,
        "single_dst_pass_budgets_direct_reference_and_replay_ttrees": True,
        "two_largest_merge_workspaces_reserved": merge_workspace_reserve_bytes
        == 2 * largest_merge_workspace_bytes,
        "quota_not_exceeded_at_peak": projected_peak_free_bytes >= 0,
        "minimum_free_headroom_at_least_20_percent": (
            projected_peak_free_fraction >= minimum_free_fraction
        ),
        "local_hot_tier_protected_after_pull": (
            local_free_after_pull_bytes >= local_protected_free_bytes
        ),
    }
    status = "PASS" if all(gates.values()) else "FAIL"
    payload = {
        "schema": STORAGE_CERTIFICATE_SCHEMA,
        "status": status,
        "gate": "P5A-S_STORAGE_REACCEPTANCE",
        "promotion_status": "CANDIDATE_CURRENT_NOT_CANONICAL",
        "manifest": str(manifest_path),
        "manifest_sha256": hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
        "science_freeze_certificate": science_certificate,
        "quota": {
            "quota_bytes": quota_bytes,
            "used_bytes": used_bytes,
            "usage_authority_sha256": usage_authority_sha256,
            "minimum_free_fraction": minimum_free_fraction,
        },
        "source_projection_count": len(normalized),
        "source_projections": normalized,
        "totals": {
            **totals,
            "persistent_increment_bytes": persistent_increment_bytes,
            "simultaneous_merge_workspace_count": (
                SIMULTANEOUS_MERGE_WORKSPACE_COUNT
            ),
            "largest_merge_workspace_bytes": largest_merge_workspace_bytes,
            "merge_workspace_reserve_bytes": merge_workspace_reserve_bytes,
            "projected_peak_increment_bytes": projected_peak_increment_bytes,
            "projected_peak_used_bytes": projected_peak_used_bytes,
            "projected_peak_free_bytes": projected_peak_free_bytes,
            "projected_peak_free_fraction": projected_peak_free_fraction,
        },
        "local_hot_tier": {
            "free_bytes": local_free_bytes,
            "protected_free_bytes": local_protected_free_bytes,
            "planned_pull_bytes": planned_local_pull_bytes,
            "free_after_pull_bytes": local_free_after_pull_bytes,
        },
        "gates": gates,
        "boundaries": [
            "This is a projection receipt, not broad-production authority.",
            "One DST pass per row must emit direct-reference histograms and replay TTrees.",
            "Required replay populations may not be removed to make storage pass.",
            "THE-121/THE-122 remain closed pending THE-131 P5B and joint P5C.",
        ],
    }
    return semantic_receipt(payload, "projection_semantic_sha256")


def main() -> int:
    args = parse_args()
    try:
        payload = build_projection(args.manifest)
    except ContractError as exc:
        print(f"THE-134 storage projection failed: {exc}", file=sys.stderr)
        return 2
    args.json_out.parent.mkdir(parents=True, exist_ok=True)
    args.json_out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(args.json_out)
    return 0 if payload["status"] == "PASS" else 3


if __name__ == "__main__":
    raise SystemExit(main())
