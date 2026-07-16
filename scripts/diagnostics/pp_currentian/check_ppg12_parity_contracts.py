#!/usr/bin/env python3
"""Report PPG12 parity contract state without promoting anything.

The v1 checker is deliberately conservative. It validates the control-plane
registries, resolves current/candidate artifact pointers, checks that expected
local paths exist, and emits structured report rows. It does not turn any
candidate into a blocking gate; canonical gates are enabled only by updating the
contract registry through the canonicalization protocol.
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import subprocess
import sys
from typing import Any


REPO = Path(__file__).resolve().parents[3]
CONTRACT_DIR = REPO / "agent_context" / "analysis_contracts"
LANE_PATH = CONTRACT_DIR / "ppg12_sample_lanes.yaml"
CONTRACT_PATH = CONTRACT_DIR / "ppg12_parity_contracts.yaml"
ARTIFACT_REGISTRY_PATH = REPO / "dataOutput" / "current_recoiljets_artifacts" / "registry.json"
DEFAULT_REPORT_DIR = REPO / "dataOutput" / "ppg12Parity" / "control_plane" / "contract_reports"


def now_stamp() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def now_file_stamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def load_json(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        return json.load(handle)


def load_yaml(path: Path) -> dict[str, Any]:
    try:
        import yaml  # type: ignore

        with path.open() as handle:
            payload = yaml.safe_load(handle)
        return payload or {}
    except ModuleNotFoundError:
        pass

    ruby = subprocess.run(
        [
            "ruby",
            "-ryaml",
            "-rjson",
            "-e",
            "puts JSON.generate(YAML.load_file(ARGV[0]))",
            str(path),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    return json.loads(ruby.stdout)


def duplicate_contract_ids(path: Path) -> list[str]:
    seen: set[str] = set()
    duplicates: list[str] = []
    in_contracts = False
    pattern = re.compile(r"^  ([A-Za-z0-9_]+):\s*$")
    for line in path.read_text().splitlines():
        if line.startswith("contracts:"):
            in_contracts = True
            continue
        if not in_contracts:
            continue
        match = pattern.match(line)
        if not match:
            continue
        contract_id = match.group(1)
        if contract_id in seen:
            duplicates.append(contract_id)
        seen.add(contract_id)
    return duplicates


def rel_or_abs(path_text: str | None) -> Path | None:
    if not path_text:
        return None
    path = Path(path_text)
    if not path.is_absolute():
        path = REPO / path
    return path


def find_entry_by_id(registry: dict[str, Any], entry_id: str | None) -> dict[str, Any] | None:
    if not entry_id:
        return None
    for entry in registry.get("artifacts", []):
        if entry.get("id") == entry_id:
            return entry
    return None


def latest_entry_for_sample(registry: dict[str, Any], sample_key: str) -> dict[str, Any] | None:
    current_id = registry.get("current", {}).get(sample_key)
    current = find_entry_by_id(registry, current_id)
    if current:
        return current
    rows = [entry for entry in registry.get("artifacts", []) if entry.get("sample_key") == sample_key]
    return rows[-1] if rows else None


def validate_registries(lanes: dict[str, Any], contracts_doc: dict[str, Any]) -> list[str]:
    errors: list[str] = []
    duplicates = duplicate_contract_ids(CONTRACT_PATH)
    for contract_id in duplicates:
        errors.append(f"duplicate contract id: {contract_id}")

    lane_ids = set((lanes.get("lanes") or {}).keys())
    valid_lanes = set(contracts_doc.get("valid_lanes") or [])
    valid_rungs = set(contracts_doc.get("valid_rungs") or [])
    valid_statuses = set(contracts_doc.get("valid_statuses") or [])

    missing_lane_decl = lane_ids - valid_lanes
    for lane_id in sorted(missing_lane_decl):
        errors.append(f"lane {lane_id} is not listed in valid_lanes")

    for contract_id, contract in (contracts_doc.get("contracts") or {}).items():
        lane_id = contract.get("sample_lane")
        if lane_id not in lane_ids:
            errors.append(f"{contract_id}: missing lane definition for {lane_id}")
        if lane_id not in valid_lanes:
            errors.append(f"{contract_id}: invalid sample_lane={lane_id}")
        rung = contract.get("rung")
        if rung not in valid_rungs:
            errors.append(f"{contract_id}: invalid rung={rung}")
        status = contract.get("status")
        if status not in valid_statuses:
            errors.append(f"{contract_id}: invalid status={status}")
        checker = (contract.get("checker") or {}).get("command")
        if contract.get("canonical_status") == "canonical" and not checker:
            errors.append(f"{contract_id}: canonical contract lacks checker command")
    return errors


def result_for_contract(
    contract_id: str,
    contract: dict[str, Any],
    lane: dict[str, Any] | None,
    artifact_registry: dict[str, Any],
    validation_errors: list[str],
) -> dict[str, Any]:
    status = contract.get("status", "")
    checker = contract.get("checker") or {}
    checker_command = checker.get("command")
    plot_path = rel_or_abs(contract.get("plot_artifact_path"))
    evidence_paths: list[dict[str, Any]] = []
    missing_paths: list[str] = []
    artifact_statuses: list[str] = []
    artifact_ids: list[str] = []

    if lane:
        for sample_key in lane.get("current_sample_keys") or []:
            entry = latest_entry_for_sample(artifact_registry, sample_key)
            if not entry:
                missing_paths.append(f"artifact_registry:{sample_key}")
                continue
            artifact_statuses.append(f"{sample_key}:{entry.get('status')}")
            artifact_ids.append(str(entry.get("id")))
            for root_path in entry.get("root_paths") or []:
                path = rel_or_abs(root_path)
                exists = bool(path and path.exists())
                evidence_paths.append({"path": str(path), "exists": exists})
                if not exists and path:
                    missing_paths.append(str(path))

    if plot_path:
        exists = plot_path.exists()
        evidence_paths.append({"path": str(plot_path), "exists": exists, "kind": "plot_artifact_path"})
        if not exists:
            missing_paths.append(str(plot_path))

    local_validation = [err for err in validation_errors if err.startswith(f"{contract_id}:")]
    lane_status = lane.get("status", "") if lane else "missing_lane"
    canonical_status = contract.get("canonical_status", "not_canonical")

    if local_validation:
        result = "fail"
        message = "; ".join(local_validation)
    elif not lane:
        result = "fail"
        message = "contract sample lane is missing from lane registry"
    elif status == "blocked_missing_input":
        result = "blocked_missing_input"
        message = "contract is explicitly blocked until required inputs are registered"
    elif missing_paths:
        result = "blocked_missing_input"
        message = "one or more registered artifact or plot paths are missing"
    elif not checker_command:
        result = "draft_not_enforced"
        message = "no checker command is defined; candidate is documented but not enforceable"
    elif status == "paused_unresolved" or "paused" in lane_status:
        result = "warn"
        message = "lane or contract is paused/unresolved and cannot be promoted"
    elif canonical_status != "canonical":
        result = "pass"
        message = "candidate report inputs resolve; result is non-blocking because contract is not canonical"
    else:
        result = "pass"
        message = "canonical contract inputs resolve"

    return {
        "contract_id": contract_id,
        "title": contract.get("title", ""),
        "sample_lane": contract.get("sample_lane", ""),
        "rung": contract.get("rung", ""),
        "status": status,
        "canonical_status": canonical_status,
        "checker_command": checker_command or "",
        "result": result,
        "message": message,
        "lane_status": lane_status,
        "artifact_statuses": artifact_statuses,
        "artifact_ids": artifact_ids,
        "evidence_paths": evidence_paths,
        "missing_paths": missing_paths,
        "gate_behavior": contract.get("gate_behavior", ""),
        "canonicalization_trigger": contract.get("canonicalization_trigger", ""),
    }


def write_reports(payload: dict[str, Any], json_path: Path, csv_path: Path) -> None:
    json_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

    rows = payload["results"]
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "contract_id",
                "result",
                "status",
                "canonical_status",
                "sample_lane",
                "rung",
                "lane_status",
                "artifact_statuses",
                "missing_paths",
                "message",
            ],
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "contract_id": row["contract_id"],
                    "result": row["result"],
                    "status": row["status"],
                    "canonical_status": row["canonical_status"],
                    "sample_lane": row["sample_lane"],
                    "rung": row["rung"],
                    "lane_status": row["lane_status"],
                    "artifact_statuses": ";".join(row["artifact_statuses"]),
                    "missing_paths": ";".join(row["missing_paths"]),
                    "message": row["message"],
                }
            )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", action="append", help="Contract id to report. May be repeated.")
    parser.add_argument("--dry-run", action="store_true", help="Print report JSON and do not write default files.")
    parser.add_argument("--strict", action="store_true", help="Return nonzero on fail or blocked_missing_input rows.")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_REPORT_DIR)
    parser.add_argument("--json-output", type=Path)
    parser.add_argument("--csv-output", type=Path)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    lanes = load_yaml(LANE_PATH)
    contracts_doc = load_yaml(CONTRACT_PATH)
    artifact_registry = load_json(ARTIFACT_REGISTRY_PATH)

    validation_errors = validate_registries(lanes, contracts_doc)
    contracts = contracts_doc.get("contracts") or {}
    selected = args.contract or list(contracts.keys())
    unknown = sorted(set(selected) - set(contracts.keys()))
    if unknown:
        for contract_id in unknown:
            print(f"[ERROR] unknown contract id: {contract_id}", file=sys.stderr)
        return 2

    lane_map = lanes.get("lanes") or {}
    results = [
        result_for_contract(contract_id, contracts[contract_id], lane_map.get(contracts[contract_id].get("sample_lane")), artifact_registry, validation_errors)
        for contract_id in selected
    ]

    payload = {
        "schema_version": 1,
        "generated_at": now_stamp(),
        "dry_run": bool(args.dry_run),
        "contracts_path": str(CONTRACT_PATH),
        "lanes_path": str(LANE_PATH),
        "artifact_registry_path": str(ARTIFACT_REGISTRY_PATH),
        "validation_errors": validation_errors,
        "results": results,
        "summary": {
            "pass": sum(1 for row in results if row["result"] == "pass"),
            "warn": sum(1 for row in results if row["result"] == "warn"),
            "fail": sum(1 for row in results if row["result"] == "fail"),
            "blocked_missing_input": sum(1 for row in results if row["result"] == "blocked_missing_input"),
            "draft_not_enforced": sum(1 for row in results if row["result"] == "draft_not_enforced"),
        },
    }

    write_outputs = not args.dry_run or args.json_output or args.csv_output
    if write_outputs:
        stamp = now_file_stamp()
        json_path = args.json_output or args.output_dir / f"ppg12_contract_report_{stamp}.json"
        csv_path = args.csv_output or args.output_dir / f"ppg12_contract_report_{stamp}.csv"
        write_reports(payload, json_path, csv_path)
        payload["written_outputs"] = {"json": str(json_path), "csv": str(csv_path)}

    print(json.dumps(payload, indent=2, sort_keys=True))

    if validation_errors:
        return 1
    if args.strict and any(row["result"] in {"fail", "blocked_missing_input"} for row in results):
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
