#!/usr/bin/env python3
"""Certify exact V10/V14 fallback versus THE134_FAST_EXTRACTION_V1 parity.

The candidate is allowed to change execution provenance and resource use only.
Every normalized ``RJPhotonTrainingViewV1`` branch and every non-code metadata
identity must remain exactly equal, including row order, stable identities,
features, labels, weights, source fields, and valid-empty state.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import os
import re
import sys
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import uproot


HERE = Path(__file__).resolve()
REPOSITORY = HERE.parents[3]
PREPARER_PATH = (
    REPOSITORY / "scripts" / "ml" / "training" / "prepare_the134_h70_matrix.py"
)
CERTIFICATE_SCHEMA = "THE134_FAST_EXTRACTION_DIFFERENTIAL_CERTIFICATE_V1"
RESOURCE_SCHEMA = "THE134_FAST_EXTRACTION_RESOURCE_PAIR_V1"
TREE_NAME = "RJPhotonTrainingViewV1"
SYSTEMS = ("pp", "auau")
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


class ValidationError(RuntimeError):
    """Fail-closed fast-extraction differential error."""


def load_preparer() -> Any:
    spec = importlib.util.spec_from_file_location(
        "prepare_the134_h70_matrix_fast_differential", PREPARER_PATH
    )
    if spec is None or spec.loader is None:
        raise ValidationError(f"cannot load authoritative reader: {PREPARER_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def exact_value_equal(left: Any, right: Any) -> bool:
    """Compare nested uproot materializations exactly, treating NaN as null."""

    if isinstance(left, np.ndarray) or isinstance(right, np.ndarray):
        left_array = np.asarray(left)
        right_array = np.asarray(right)
        if left_array.shape != right_array.shape:
            return False
        if left_array.dtype == object or right_array.dtype == object:
            return all(
                exact_value_equal(lvalue, rvalue)
                for lvalue, rvalue in zip(left_array.flat, right_array.flat)
            )
        if left_array.dtype.kind in "fc" or right_array.dtype.kind in "fc":
            return bool(
                np.all(
                    (left_array == right_array)
                    | (np.isnan(left_array) & np.isnan(right_array))
                )
            )
        return bool(np.array_equal(left_array, right_array))
    if isinstance(left, (list, tuple)) or isinstance(right, (list, tuple)):
        if not isinstance(left, (list, tuple)) or not isinstance(
            right, (list, tuple)
        ):
            return False
        return len(left) == len(right) and all(
            exact_value_equal(lvalue, rvalue)
            for lvalue, rvalue in zip(left, right)
        )
    if isinstance(left, (float, np.floating)) or isinstance(
        right, (float, np.floating)
    ):
        return bool(
            left == right
            or (
                isinstance(left, (float, np.floating))
                and isinstance(right, (float, np.floating))
                and math.isnan(float(left))
                and math.isnan(float(right))
            )
        )
    return bool(left == right)


def root_inventory(path: Path) -> dict[str, str]:
    with uproot.open(path) as root_file:
        return {
            str(key): str(class_name)
            for key, class_name in root_file.classnames(
                recursive=True, cycle=False
            ).items()
        }


def compare_sidecars(
    baseline: Path,
    candidate: Path,
    *,
    baseline_code_sha256: str,
    candidate_code_sha256: str,
    preparer: Any | None = None,
) -> dict[str, Any]:
    for label, path in (("baseline", baseline), ("candidate", candidate)):
        if not path.is_file():
            raise ValidationError(f"{label} sidecar is missing: {path}")
    for label, value in (
        ("baseline code SHA-256", baseline_code_sha256),
        ("candidate code SHA-256", candidate_code_sha256),
    ):
        if SHA256_RE.fullmatch(value) is None:
            raise ValidationError(f"{label} is malformed")

    reader = preparer or load_preparer()
    baseline_arrays, baseline_branches, baseline_metadata = reader.read_tree(
        baseline, TREE_NAME
    )
    candidate_arrays, candidate_branches, candidate_metadata = reader.read_tree(
        candidate, TREE_NAME
    )
    failures: list[str] = []
    if baseline_branches != candidate_branches:
        failures.append("branch_inventory_or_order")
    if root_inventory(baseline) != root_inventory(candidate):
        failures.append("root_key_or_class_inventory")
    if set(baseline_arrays) != set(candidate_arrays):
        failures.append("materialized_branch_inventory")
    for branch in sorted(set(baseline_arrays) & set(candidate_arrays)):
        if not exact_value_equal(
            baseline_arrays[branch], candidate_arrays[branch]
        ):
            failures.append(f"branch:{branch}")

    if set(baseline_metadata) != set(candidate_metadata):
        failures.append("metadata_inventory")
    if str(baseline_metadata.get("code_sha256", "")) != baseline_code_sha256:
        failures.append("baseline_code_sha256")
    if str(candidate_metadata.get("code_sha256", "")) != candidate_code_sha256:
        failures.append("candidate_code_sha256")
    for name in sorted(set(baseline_metadata) & set(candidate_metadata)):
        if name == "code_sha256":
            continue
        if baseline_metadata[name] != candidate_metadata[name]:
            failures.append(f"metadata:{name}")

    entries = int(baseline_metadata.get("tree_num_entries", -1))
    if entries != int(candidate_metadata.get("tree_num_entries", -2)):
        failures.append("tree_entries")
    return {
        "status": "PASS" if not failures else "FAIL",
        "baseline": {
            "path": str(baseline),
            "sha256": sha256_file(baseline),
            "bytes": baseline.stat().st_size,
            "code_sha256": baseline_code_sha256,
        },
        "candidate": {
            "path": str(candidate),
            "sha256": sha256_file(candidate),
            "bytes": candidate.stat().st_size,
            "code_sha256": candidate_code_sha256,
        },
        "tree_entries": entries,
        "branch_count": len(baseline_branches),
        "failures": failures,
    }


def _positive_number(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValidationError(f"{label} must be numeric")
    result = float(value)
    if not math.isfinite(result) or result <= 0:
        raise ValidationError(f"{label} must be finite and positive")
    return result


def validate_resources(
    payload: Mapping[str, Any], *, minimum_speedup: float
) -> dict[str, Any]:
    if set(payload) != {"schema", "rows"} or payload.get("schema") != RESOURCE_SCHEMA:
        raise ValidationError("resource-pair schema or field inventory differs")
    rows = payload.get("rows")
    if not isinstance(rows, Mapping) or set(rows) != set(SYSTEMS):
        raise ValidationError("resource-pair system inventory differs")
    reports: dict[str, Any] = {}
    failures: list[str] = []
    for system in SYSTEMS:
        row = rows[system]
        if not isinstance(row, Mapping) or set(row) != {"baseline", "candidate"}:
            raise ValidationError(f"{system} resource arm inventory differs")
        arms: dict[str, dict[str, float]] = {}
        for arm in ("baseline", "candidate"):
            record = row[arm]
            expected = {
                "elapsed_seconds",
                "max_rss_kb",
                "analysis_output_bytes",
                "sidecar_bytes",
            }
            if not isinstance(record, Mapping) or set(record) != expected:
                raise ValidationError(f"{system} {arm} resource fields differ")
            arms[arm] = {
                name: _positive_number(
                    record[name], f"{system}.{arm}.{name}"
                )
                for name in expected
            }
        speedup = (
            arms["baseline"]["elapsed_seconds"]
            / arms["candidate"]["elapsed_seconds"]
        )
        if speedup < minimum_speedup:
            failures.append(f"{system}:speedup={speedup:.6g}")
        if arms["candidate"]["max_rss_kb"] > arms["baseline"]["max_rss_kb"]:
            failures.append(f"{system}:rss_regression")
        if (
            arms["candidate"]["analysis_output_bytes"]
            >= arms["baseline"]["analysis_output_bytes"]
        ):
            failures.append(f"{system}:scratch_output_not_reduced")
        reports[system] = {
            "baseline": arms["baseline"],
            "candidate": arms["candidate"],
            "speedup": speedup,
            "rss_ratio": (
                arms["candidate"]["max_rss_kb"]
                / arms["baseline"]["max_rss_kb"]
            ),
            "analysis_output_ratio": (
                arms["candidate"]["analysis_output_bytes"]
                / arms["baseline"]["analysis_output_bytes"]
            ),
        }
    return {
        "status": "PASS" if not failures else "FAIL",
        "minimum_speedup": minimum_speedup,
        "rows": reports,
        "failures": failures,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    for system in SYSTEMS:
        parser.add_argument(f"--baseline-{system}", type=Path, required=True)
        parser.add_argument(f"--candidate-{system}", type=Path, required=True)
    parser.add_argument("--baseline-code-sha256", required=True)
    parser.add_argument("--candidate-code-sha256", required=True)
    parser.add_argument("--resource-pair", type=Path, required=True)
    parser.add_argument("--minimum-speedup", type=float, default=1.5)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def write_new_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    descriptor = os.open(path, flags, 0o644)
    with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
        json.dump(payload, stream, indent=2, sort_keys=True)
        stream.write("\n")


def main() -> int:
    args = parse_args()
    try:
        if not math.isfinite(args.minimum_speedup) or args.minimum_speedup < 1.0:
            raise ValidationError("minimum speedup must be finite and at least one")
        resource_payload = json.loads(args.resource_pair.read_text(encoding="utf-8"))
        preparer = load_preparer()
        sidecars = {
            system: compare_sidecars(
                getattr(args, f"baseline_{system}"),
                getattr(args, f"candidate_{system}"),
                baseline_code_sha256=args.baseline_code_sha256,
                candidate_code_sha256=args.candidate_code_sha256,
                preparer=preparer,
            )
            for system in SYSTEMS
        }
        resources = validate_resources(
            resource_payload, minimum_speedup=args.minimum_speedup
        )
        failures = [
            f"{system}:{failure}"
            for system, report in sidecars.items()
            for failure in report["failures"]
        ]
        failures.extend(f"resource:{failure}" for failure in resources["failures"])
        certificate = {
            "schema": CERTIFICATE_SCHEMA,
            "status": "PASS" if not failures else "FAIL",
            "authority": "BOUNDED_THE134_EXTRACTION_EXECUTOR_ONLY",
            "ordinary_recoiljets_execution_unchanged": True,
            "full_training_authority": 0,
            "broad_production_authority": False,
            "baseline_code_sha256": args.baseline_code_sha256,
            "candidate_code_sha256": args.candidate_code_sha256,
            "sidecars": sidecars,
            "resources": resources,
            "failures": failures,
        }
        write_new_json(args.output, certificate)
        print(
            json.dumps(
                {
                    "status": certificate["status"],
                    "output": str(args.output),
                    "sha256": sha256_file(args.output),
                    "failures": failures,
                },
                sort_keys=True,
            )
        )
        return 0 if not failures else 1
    except (OSError, ValueError, json.JSONDecodeError, ValidationError) as exc:
        print(f"[THE134-FAST-DIFFERENTIAL][ERROR] {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
