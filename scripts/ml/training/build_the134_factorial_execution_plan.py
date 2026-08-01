#!/usr/bin/env python3
"""Build the deterministic 14-lane THE-134 model execution plan.

This tool never trains a model.  It binds the two certified single-read
factorial matrices, all fourteen view-qualified extraction audits, and the
only two frozen model-reuse authorities into exact argv arrays for the
existing per-view runner.  The resulting plan is safe to inspect and route to
a separately authorized bounded executor after complete extraction.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
import tempfile
from pathlib import Path
from typing import Any

_HERE = Path(__file__).resolve()
_CONTRACTS = _HERE.parents[1] / "contracts"
sys.path.insert(0, str(_CONTRACTS))
from the134_h70_contract import (  # noqa: E402
    ALL_SHOWER_VIEWS,
    expected_model_origin,
    shower_semantic_sha256,
)


SCHEMA = "THE134_FACTORIAL_MODEL_EXECUTION_PLAN_V1"
MATRIX_MANIFEST_SCHEMA = "THE134_SINGLE_READ_FACTORIAL_MATRIX_MANIFEST_V1"
MATRIX_AUDIT_SCHEMA = "THE134_FACTORIAL_VIEW_TRAINING_MATRIX_AUDIT_V1"
REUSE_AUDIT_SCHEMA = "THE134_FACTORIAL_VIEW_MODEL_REUSE_AUDIT_V1"
SYSTEMS = ("pp", "auau")
REQUIRED_REUSE = {("pp", "H70"), ("auau", "H0")}
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


class PlanError(ValueError):
    """Fail-closed execution-plan construction error."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path, label: str) -> dict[str, Any]:
    if not path.is_file():
        raise PlanError(f"{label} is missing: {path}")
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise PlanError(f"{label} is not readable JSON: {path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise PlanError(f"{label} must be a JSON object: {path}")
    return payload


def require_exact_keys(payload: dict[str, Any], expected: set[str], label: str) -> None:
    observed = set(payload)
    if observed != expected:
        raise PlanError(
            f"{label} field inventory mismatch: "
            f"missing={sorted(expected - observed)} extra={sorted(observed - expected)}"
        )


def require_sha256(value: object, label: str) -> str:
    text = str(value)
    if not SHA256_RE.fullmatch(text):
        raise PlanError(f"{label} is not a lowercase SHA-256")
    return text


def verify_registered_file(path_value: object, sha_value: object, label: str) -> tuple[Path, str]:
    path = Path(str(path_value)).expanduser().resolve(strict=True)
    if not path.is_file():
        raise PlanError(f"{label} is not a regular file: {path}")
    expected = require_sha256(sha_value, f"{label} SHA-256")
    observed = sha256_file(path)
    if observed != expected:
        raise PlanError(
            f"{label} SHA-256 drift: expected={expected} observed={observed}"
        )
    return path, observed


def require_absolute_output_root(path: Path) -> Path:
    if not path.is_absolute():
        raise PlanError("output root must be an absolute path")
    normalized = Path(os.path.normpath(str(path)))
    if str(normalized) != str(path) or str(path) == "/":
        raise PlanError("output root must be normalized and cannot be filesystem root")
    if path.exists():
        raise PlanError(f"output root already exists; refusing possible duplicate work: {path}")
    return path


def validate_matrix_manifest(path: Path, system: str) -> dict[str, Any]:
    payload = read_json(path, f"{system} factorial matrix manifest")
    require_exact_keys(
        payload,
        {
            "schema",
            "status",
            "system",
            "matrix",
            "matrix_sha256",
            "root_reads",
            "input_count",
            "audits",
        },
        f"{system} factorial matrix manifest",
    )
    if payload.get("schema") != MATRIX_MANIFEST_SCHEMA:
        raise PlanError(f"{system} factorial matrix manifest schema mismatch")
    if payload.get("status") != "PASS" or payload.get("system") != system:
        raise PlanError(f"{system} factorial matrix manifest identity/status mismatch")
    if not isinstance(payload.get("root_reads"), int) or payload["root_reads"] <= 0:
        raise PlanError(f"{system} factorial matrix manifest has invalid root_reads")
    if payload.get("input_count") != payload.get("root_reads"):
        raise PlanError(f"{system} factorial matrix manifest read/input count mismatch")

    matrix, matrix_sha = verify_registered_file(
        payload.get("matrix"), payload.get("matrix_sha256"), f"{system} shared matrix"
    )
    audits = payload.get("audits")
    if not isinstance(audits, dict) or set(audits) != set(ALL_SHOWER_VIEWS):
        observed = sorted(audits) if isinstance(audits, dict) else []
        raise PlanError(
            f"{system} factorial audit inventory mismatch: observed={observed}"
        )

    audit_records: dict[str, dict[str, str]] = {}
    for view in ALL_SHOWER_VIEWS:
        record = audits[view]
        if not isinstance(record, dict):
            raise PlanError(f"{system}/{view} audit record must be an object")
        require_exact_keys(record, {"path", "sha256"}, f"{system}/{view} audit record")
        audit_path, audit_sha = verify_registered_file(
            record.get("path"), record.get("sha256"), f"{system}/{view} extraction audit"
        )
        audit = read_json(audit_path, f"{system}/{view} extraction audit")
        checks = {
            "schema": audit.get("schema") == MATRIX_AUDIT_SCHEMA,
            "status": audit.get("status") == "PASS",
            "system": audit.get("system") == system,
            "view": audit.get("shower_definition") == view,
            "semantic": audit.get("shower_semantic_sha256")
            == shower_semantic_sha256(view),
            "matrix": Path(str(audit.get("matrix", ""))).resolve() == matrix,
            "matrix_sha256": audit.get("matrix_sha256") == matrix_sha,
            "full_training_authority": audit.get("full_training_authority") == 1,
            "source_population_closure": isinstance(
                audit.get("source_population_closure"), dict
            )
            and audit["source_population_closure"].get("status") == "PASS",
            "single_read": isinstance(audit.get("single_read_factorial_projection"), dict)
            and audit["single_read_factorial_projection"].get("status") == "PASS",
        }
        if not all(checks.values()):
            failed = sorted(name for name, passed in checks.items() if not passed)
            raise PlanError(f"{system}/{view} extraction audit failed gates: {failed}")
        audit_records[view] = {"path": str(audit_path), "sha256": audit_sha}

    return {
        "manifest": str(path.resolve()),
        "manifest_sha256": sha256_file(path.resolve()),
        "matrix": str(matrix),
        "matrix_sha256": matrix_sha,
        "root_reads": payload["root_reads"],
        "input_count": payload["input_count"],
        "audits": audit_records,
    }


def parse_reuse_binding(value: str) -> tuple[tuple[str, str], Path]:
    key, separator, raw_path = value.partition("=")
    if not separator:
        raise PlanError(f"reuse audit must use SYSTEM:VIEW=PATH syntax: {value!r}")
    system, colon, view = key.partition(":")
    if not colon or system not in SYSTEMS or view not in ALL_SHOWER_VIEWS:
        raise PlanError(f"invalid reuse audit key: {key!r}")
    return (system, view), Path(raw_path).expanduser().resolve(strict=True)


def validate_reuse_audits(values: list[str]) -> dict[tuple[str, str], dict[str, str]]:
    observed: dict[tuple[str, str], Path] = {}
    for value in values:
        key, path = parse_reuse_binding(value)
        if key in observed:
            raise PlanError(f"duplicate reuse audit binding: {key[0]}:{key[1]}")
        observed[key] = path
    if set(observed) != REQUIRED_REUSE:
        raise PlanError(
            "reuse audit inventory mismatch: "
            f"missing={sorted(REQUIRED_REUSE - set(observed))} "
            f"extra={sorted(set(observed) - REQUIRED_REUSE)}"
        )

    validated: dict[tuple[str, str], dict[str, str]] = {}
    for (system, view), path in observed.items():
        payload = read_json(path, f"{system}/{view} reuse audit")
        checks = {
            "schema": payload.get("schema") == REUSE_AUDIT_SCHEMA,
            "status": payload.get("status") == "PASS",
            "system": payload.get("system") == system,
            "view": payload.get("shower_definition") == view,
            "semantic": payload.get("shower_semantic_sha256")
            == shower_semantic_sha256(view),
            "origin": payload.get("model_origin") == expected_model_origin(system, view),
        }
        if not all(checks.values()):
            failed = sorted(name for name, passed in checks.items() if not passed)
            raise PlanError(f"{system}/{view} reuse audit failed gates: {failed}")
        validated[(system, view)] = {
            "path": str(path),
            "sha256": sha256_file(path),
        }
    return validated


def build_plan(args: argparse.Namespace) -> dict[str, Any]:
    output_root = require_absolute_output_root(args.output_root)
    try:
        python = Path(args.python).expanduser().resolve(strict=True)
    except OSError as exc:
        raise PlanError(f"Python runtime is missing or unreadable: {args.python}") from exc
    if not python.is_file() or not os.access(python, os.X_OK):
        raise PlanError(f"Python runtime is not an executable regular file: {python}")
    runner = args.runner.expanduser().resolve(strict=True)
    if not runner.is_file():
        raise PlanError(f"model runner is not a regular file: {runner}")
    if args.max_parallel_views < 1 or args.max_parallel_views > len(ALL_SHOWER_VIEWS):
        raise PlanError("max parallel views must be between 1 and 7")

    matrices = {
        "pp": validate_matrix_manifest(args.pp_matrix_manifest.resolve(), "pp"),
        "auau": validate_matrix_manifest(args.auau_matrix_manifest.resolve(), "auau"),
    }
    reuse = validate_reuse_audits(args.reuse_audit)
    lanes: list[dict[str, Any]] = []
    for system in SYSTEMS:
        for view in ALL_SHOWER_VIEWS:
            lane_root = output_root / system / view.lower()
            origin = expected_model_origin(system, view)
            command = [
                str(python),
                str(runner),
                "--system",
                system,
                "--view",
                view,
                "--matrix",
                matrices[system]["matrix"],
                "--extraction-audit",
                matrices[system]["audits"][view]["path"],
                "--outdir",
                str(lane_root),
            ]
            reuse_record = reuse.get((system, view))
            if reuse_record is not None:
                command.extend(["--reuse-audit", reuse_record["path"]])
            command.append("--execute")
            lanes.append(
                {
                    "lane_id": f"{system}:{view}",
                    "system": system,
                    "view": view,
                    "model_origin": origin,
                    "status": "READY_NOT_EXECUTED",
                    "matrix_sha256": matrices[system]["matrix_sha256"],
                    "extraction_audit": matrices[system]["audits"][view],
                    "reuse_audit": reuse_record,
                    "output_directory": str(lane_root),
                    "command": command,
                }
            )

    return {
        "schema": SCHEMA,
        "status": "READY_NOT_EXECUTED",
        "promotion_status": "CANDIDATE_CURRENT_NOT_CANONICAL",
        "execution_performed": False,
        "systems": list(SYSTEMS),
        "views": list(ALL_SHOWER_VIEWS),
        "model_count": len(lanes),
        "max_parallel_views": args.max_parallel_views,
        "parallelism_scope": "WITHIN_ONE_CERTIFIED_SYSTEM_AGGREGATE",
        "output_root": str(output_root),
        "python": str(python),
        "runner": {"path": str(runner), "sha256": sha256_file(runner)},
        "matrix_manifests": matrices,
        "reuse_authorities": {
            f"{system}:{view}": record
            for (system, view), record in sorted(reuse.items())
        },
        "lanes": lanes,
        "hard_stops": [
            "NO_AUTHORITY_BEFORE_COMPLETE_SYSTEM_MATRIX_CERTIFICATE",
            "NO_IMPLICIT_RETRAINING_OF_AUTHORIZED_REUSE_LANES",
            "NO_SCIENTIFIC_GATE_OR_TOLERANCE_CHANGE",
            "NO_AUTOMATIC_CANONICAL_PROMOTION",
        ],
    }


def write_json_idempotent(path: Path, payload: dict[str, Any]) -> None:
    content = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if path.exists():
        if path.read_text(encoding="utf-8") == content:
            return
        raise PlanError(f"execution plan exists with different content: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        handle.write(content)
        temporary = Path(handle.name)
    os.replace(temporary, path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pp-matrix-manifest", required=True, type=Path)
    parser.add_argument("--auau-matrix-manifest", required=True, type=Path)
    parser.add_argument(
        "--reuse-audit",
        action="append",
        default=[],
        metavar="SYSTEM:VIEW=PATH",
        help="Repeat exactly for pp:H70 and auau:H0.",
    )
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--json-out", required=True, type=Path)
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument(
        "--runner", type=Path, default=_HERE.with_name("run_the134_h70_model.py")
    )
    parser.add_argument("--max-parallel-views", type=int, default=2)
    return parser.parse_args()


def main() -> int:
    try:
        args = parse_args()
        payload = build_plan(args)
        write_json_idempotent(args.json_out, payload)
    except (OSError, PlanError, RuntimeError) as exc:
        raise SystemExit(str(exc)) from exc
    print(args.json_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
