#!/usr/bin/env python3
"""Check, render, and append ThesisAnalysis artifact provenance records."""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from codex_work_register_common import first_line, load_register


DEFAULT_REGISTRY = Path("agent_context/ARTIFACT_REGISTRY.yaml")

DEFAULT_REQUIRED_FIELDS = [
    "artifact_id",
    "kind",
    "path",
    "thesis_claim",
    "source_inputs",
    "generation_command",
    "qa_state",
    "qa_evidence",
    "limitations",
    "destination",
    "last_verified",
]

DEFAULT_QUALITY_STATES = {"candidate", "qa_passed", "slide_ready", "superseded", "invalid"}


def is_remote_or_virtual(path_text: str) -> bool:
    return path_text.startswith(("http://", "https://", "gs://", "s3://")) or path_text.startswith("remote:")


def validate_registry(data: dict[str, Any]) -> list[str]:
    errors: list[str] = []
    required = data.get("required_artifact_fields") or DEFAULT_REQUIRED_FIELDS
    quality_states = set(data.get("quality_states") or DEFAULT_QUALITY_STATES)
    claim_layers = set(data.get("claim_layers") or [])
    artifacts = data.get("artifacts")
    if not isinstance(artifacts, list):
        return ["artifacts must be a list"]

    seen: set[str] = set()
    for index, item in enumerate(artifacts, start=1):
        prefix = f"artifacts[{index}]"
        if not isinstance(item, dict):
            errors.append(f"{prefix} must be a mapping")
            continue
        artifact_id = item.get("artifact_id")
        if not isinstance(artifact_id, str) or not artifact_id.strip():
            errors.append(f"{prefix}.artifact_id must be a non-empty string")
        elif artifact_id in seen:
            errors.append(f"{prefix}.artifact_id duplicates {artifact_id}")
        else:
            seen.add(artifact_id)
        for field in required:
            if field not in item or item.get(field) in (None, ""):
                errors.append(f"{prefix} missing required field: {field}")
        for list_field in ("source_inputs", "qa_evidence"):
            if list_field in item and not isinstance(item.get(list_field), list):
                errors.append(f"{prefix}.{list_field} must be a list")
        qa_state = item.get("qa_state")
        if qa_state not in quality_states:
            errors.append(f"{prefix}.qa_state has invalid state: {qa_state!r}")
        thesis_claim = item.get("thesis_claim")
        if claim_layers and thesis_claim not in claim_layers:
            errors.append(f"{prefix}.thesis_claim is not in claim_layers: {thesis_claim!r}")
        path_text = first_line(item.get("path"))
        if path_text and not is_remote_or_virtual(path_text) and not Path(path_text).exists():
            errors.append(f"{prefix}.path does not exist locally: {path_text}")
        if qa_state in {"qa_passed", "slide_ready"} and not item.get("qa_evidence"):
            errors.append(f"{prefix}.qa_state={qa_state} requires qa_evidence")
    return errors


def command_check(args: argparse.Namespace) -> int:
    data = load_register(Path(args.registry))
    errors = validate_registry(data)
    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1
    print(f"OK: {args.registry} validates with {len(data.get('artifacts') or [])} artifacts")
    return 0


def command_render(args: argparse.Namespace) -> int:
    data = load_register(Path(args.registry))
    artifacts = data.get("artifacts") or []
    by_claim: dict[str, list[dict[str, Any]]] = {}
    for item in artifacts:
        if isinstance(item, dict):
            by_claim.setdefault(str(item.get("thesis_claim")), []).append(item)

    if args.json:
        print(json.dumps({"artifacts": artifacts, "by_claim": by_claim}, indent=2, sort_keys=True))
        return 0

    print("# Artifact Registry")
    print()
    for claim in sorted(by_claim):
        print(f"## {claim}")
        for item in by_claim[claim]:
            print(
                f"- {item.get('artifact_id')} [{item.get('kind')}, {item.get('qa_state')}]: "
                f"{item.get('path')}"
            )
        print()
    return 0


def command_add(args: argparse.Namespace) -> int:
    try:
        import yaml  # type: ignore
    except ModuleNotFoundError:
        print("ERROR: adding artifacts requires PyYAML in this environment", file=sys.stderr)
        return 2

    path = Path(args.registry)
    data = load_register(path)
    artifacts = data.setdefault("artifacts", [])
    if not isinstance(artifacts, list):
        print("ERROR: artifacts must be a list", file=sys.stderr)
        return 1
    if any(isinstance(item, dict) and item.get("artifact_id") == args.artifact_id for item in artifacts):
        print(f"ERROR: artifact already exists: {args.artifact_id}", file=sys.stderr)
        return 1

    artifact = {
        "artifact_id": args.artifact_id,
        "kind": args.kind,
        "path": args.path,
        "thesis_claim": args.thesis_claim,
        "source_inputs": args.source_input or [],
        "generation_command": args.generation_command,
        "qa_state": args.qa_state,
        "qa_evidence": args.qa_evidence or [],
        "limitations": args.limitations,
        "destination": args.destination,
        "last_verified": args.last_verified or datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }
    artifacts.append(artifact)
    data["last_updated"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    errors = validate_registry(data)
    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1
    path.write_text(yaml.safe_dump(data, sort_keys=False, allow_unicode=False), encoding="utf-8")
    print(f"OK: added artifact {args.artifact_id}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    check = subparsers.add_parser("check", help="validate artifact provenance")
    check.add_argument("registry", nargs="?", default=str(DEFAULT_REGISTRY))
    check.set_defaults(func=command_check)

    render = subparsers.add_parser("render", help="render artifact registry summary")
    render.add_argument("registry", nargs="?", default=str(DEFAULT_REGISTRY))
    render.add_argument("--json", action="store_true")
    render.set_defaults(func=command_render)

    add = subparsers.add_parser("add", help="append one artifact record")
    add.add_argument("--registry", default=str(DEFAULT_REGISTRY))
    add.add_argument("--artifact-id", required=True)
    add.add_argument("--kind", required=True)
    add.add_argument("--path", required=True)
    add.add_argument("--thesis-claim", required=True)
    add.add_argument("--source-input", action="append")
    add.add_argument("--generation-command", required=True)
    add.add_argument("--qa-state", required=True, choices=sorted(DEFAULT_QUALITY_STATES))
    add.add_argument("--qa-evidence", action="append")
    add.add_argument("--limitations", required=True)
    add.add_argument("--destination", required=True)
    add.add_argument("--last-verified")
    add.set_defaults(func=command_add)

    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
