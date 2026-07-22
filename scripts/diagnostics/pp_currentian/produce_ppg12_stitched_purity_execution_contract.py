#!/usr/bin/env python3
"""Produce and verify one canonical THE-97 lane execution contract.

The receipt binds the interaction-specific five-column source graph and the
small environment surface that changes pp-SIM reconstruction.  It is intended
to be written inside the admission or production job after its environment is
frozen; lane extraction re-runs the verifier and carries the semantic digest
into the admission manifest.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[3]
DEFAULT_CONTRACT = (
    REPO / "agent_context/analysis_contracts/ppg12_stitched_purity_closure.yaml"
)
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


class ExecutionContractError(RuntimeError):
    """The executed lane graph cannot be certified."""


def _canonical_bytes(payload: Any) -> bytes:
    return json.dumps(
        payload, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode()


def _payload_sha256(payload: Any) -> str:
    return hashlib.sha256(_canonical_bytes(payload)).hexdigest()


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
    except FileNotFoundError as exc:
        raise ExecutionContractError(f"missing input file: {path}") from exc
    return digest.hexdigest()


def _read_json(path: Path, label: str) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise ExecutionContractError(f"cannot read {label} {path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise ExecutionContractError(f"{label} must contain a JSON object")
    return payload


def _resolve_link(link: Any, parent: Path, label: str) -> dict[str, str]:
    if not isinstance(link, dict):
        raise ExecutionContractError(f"{label} must contain path and sha256")
    raw_path = link.get("path")
    expected = link.get("sha256")
    if not isinstance(raw_path, str) or not raw_path:
        raise ExecutionContractError(f"{label}.path is required")
    if not isinstance(expected, str) or SHA256_RE.fullmatch(expected) is None:
        raise ExecutionContractError(f"{label}.sha256 is invalid")
    path = Path(raw_path).expanduser()
    path = path.resolve() if path.is_absolute() else (parent / path).resolve()
    observed = _file_sha256(path)
    if observed != expected:
        raise ExecutionContractError(
            f"{label} hash mismatch: expected {expected}, observed {observed}"
        )
    return {"path": str(path), "sha256": observed}


def _read_source_rows(path: Path) -> list[list[str]]:
    rows: list[list[str]] = []
    for raw in path.read_text().splitlines():
        text = raw.strip()
        if not text or text.startswith("#"):
            continue
        fields = text.split()
        if len(fields) != 5:
            raise ExecutionContractError(
                f"source graph row {len(rows) + 1} has {len(fields)} rather than 5 columns"
            )
        rows.append(fields)
    if not rows:
        raise ExecutionContractError("source graph is empty")
    return rows


def _identity(lane_id: str, contract: dict[str, Any]) -> dict[str, str]:
    parts = lane_id.split(":")
    if len(parts) != 4:
        raise ExecutionContractError(f"invalid lane_id: {lane_id}")
    family, sample, period, interaction = parts
    expected = {
        f"{candidate_family}:{candidate_sample}:{candidate_period}:{candidate_interaction}"
        for candidate_family, spec in contract["families"].items()
        for candidate_sample in spec["samples"]
        for candidate_period in contract["periods"]
        for candidate_interaction in contract["interactions"]
    }
    if lane_id not in expected:
        raise ExecutionContractError(f"noncanonical lane_id: {lane_id}")
    return {
        "lane_id": lane_id,
        "family": family,
        "sample": sample,
        "period": period,
        "interaction": interaction,
    }


def _expected_environment(
    identity: dict[str, str], contract: dict[str, Any]
) -> dict[str, str]:
    spec = contract.get("lane_execution_contract")
    if not isinstance(spec, dict):
        raise ExecutionContractError("closure contract lacks lane_execution_contract")
    common = spec.get("required_environment_common")
    by_interaction = spec.get("required_environment_by_interaction")
    if not isinstance(common, dict) or not isinstance(by_interaction, dict):
        raise ExecutionContractError("lane execution environment contract is malformed")
    interaction_values = by_interaction.get(identity["interaction"])
    if not isinstance(interaction_values, dict):
        raise ExecutionContractError("lane interaction environment is undefined")
    expected = {str(key): str(value) for key, value in common.items()}
    expected.update({str(key): str(value) for key, value in interaction_values.items()})
    expected["RJ_PPG12_PERIOD"] = identity["period"]
    expected["RJ_PPG12_CROSSING_PERIOD"] = identity["period"]
    return dict(sorted(expected.items()))


def _validate_source_graph(
    rows: list[list[str]], identity: dict[str, str], contract: dict[str, Any]
) -> dict[str, Any]:
    spec = contract["lane_execution_contract"]
    expected_columns = spec.get("source_graph_columns")
    expected_graph = spec.get("source_graph")
    namespaces = spec.get("source_namespace_by_interaction")
    if expected_columns != ["calo", "g4", "truthjet", "global", "mbd"]:
        raise ExecutionContractError("source graph columns are not canonical")
    if expected_graph != ["NONE", "g4", "truthjet", "NONE", "NONE"]:
        raise ExecutionContractError("source graph contract is not canonical")
    if not isinstance(namespaces, dict):
        raise ExecutionContractError("source namespace contract is missing")
    namespace = namespaces.get(identity["interaction"])
    if not isinstance(namespace, str) or not namespace:
        raise ExecutionContractError("source namespace is undefined")
    normalized: list[dict[str, str]] = []
    for index, fields in enumerate(rows, start=1):
        calo, g4, truthjet, global_path, mbd = fields
        if calo != "NONE" or global_path != "NONE" or mbd != "NONE":
            raise ExecutionContractError(
                f"source graph row {index} is not NONE,g4,truthjet,NONE,NONE"
            )
        if f"/{namespace}/g4hits/" not in g4:
            raise ExecutionContractError(
                f"source graph row {index} uses the wrong G4 interaction namespace"
            )
        if f"/{namespace}/nopileup/jets/" not in truthjet:
            raise ExecutionContractError(
                f"source graph row {index} uses the wrong truth-jet interaction namespace"
            )
        normalized.append({"g4": g4, "truthjet": truthjet})
    return {
        "columns": expected_columns,
        "graph": expected_graph,
        "interaction_namespace": namespace,
        "row_count": len(rows),
        "row_identity_sha256": _payload_sha256(normalized),
    }


def _build_receipt(
    lane_id: str,
    source_list_path: Path,
    environment: dict[str, Any],
    contract_path: Path = DEFAULT_CONTRACT,
) -> dict[str, Any]:
    contract_path = contract_path.resolve()
    contract = _read_json(contract_path, "closure contract")
    if contract.get("schema") != "ppg12-stitched-purity-closure-contract/v1":
        raise ExecutionContractError("unsupported closure contract schema")
    identity = _identity(lane_id, contract)
    expected_environment = _expected_environment(identity, contract)
    observed_environment = {str(key): str(value) for key, value in environment.items()}
    if observed_environment != expected_environment:
        missing = sorted(set(expected_environment) - set(observed_environment))
        extra = sorted(set(observed_environment) - set(expected_environment))
        drift = {
            key: {"expected": expected_environment[key], "observed": observed_environment[key]}
            for key in sorted(set(expected_environment) & set(observed_environment))
            if expected_environment[key] != observed_environment[key]
        }
        raise ExecutionContractError(
            f"executed environment differs from canonical lane contract; "
            f"missing={missing} extra={extra} drift={drift}"
        )
    source_list_path = source_list_path.resolve()
    source_graph = _validate_source_graph(
        _read_source_rows(source_list_path), identity, contract
    )
    producer = Path(__file__).resolve()
    semantic = {
        **identity,
        "source_graph": {
            key: source_graph[key]
            for key in ("columns", "graph", "interaction_namespace")
        },
        "environment": expected_environment,
    }
    return {
        "schema": "ppg12-stitched-purity-execution-contract/v1",
        **identity,
        "source_list": {
            "path": str(source_list_path),
            "sha256": _file_sha256(source_list_path),
        },
        "source_graph": source_graph,
        "environment": expected_environment,
        "semantic_contract_sha256": _payload_sha256(semantic),
        "contract": {
            "path": str(contract_path),
            "sha256": _payload_sha256(contract),
        },
        "producer": {"path": str(producer), "sha256": _file_sha256(producer)},
    }


def build_receipt(
    lane_id: str,
    source_list_path: Path,
    environment_path: Path,
    contract_path: Path = DEFAULT_CONTRACT,
) -> dict[str, Any]:
    environment = _read_json(environment_path.resolve(), "executed environment")
    return _build_receipt(lane_id, source_list_path, environment, contract_path)


def validate_receipt(
    receipt_path: Path,
    lane_id: str,
    expected_source_list: dict[str, str],
    contract_path: Path = DEFAULT_CONTRACT,
) -> tuple[dict[str, Any], str]:
    receipt_path = receipt_path.resolve()
    payload = _read_json(receipt_path, "execution contract receipt")
    if payload.get("schema") != "ppg12-stitched-purity-execution-contract/v1":
        raise ExecutionContractError("unsupported execution contract receipt schema")
    producer = _resolve_link(payload.get("producer"), receipt_path.parent, "producer")
    if Path(producer["path"]) != Path(__file__).resolve():
        raise ExecutionContractError("receipt was not made by the canonical producer")
    source = _resolve_link(payload.get("source_list"), receipt_path.parent, "source_list")
    normalized_expected = _resolve_link(
        expected_source_list, receipt_path.parent, "expected source_list"
    )
    if source != normalized_expected:
        raise ExecutionContractError("receipt source list differs from lane evidence")
    environment = payload.get("environment")
    if not isinstance(environment, dict):
        raise ExecutionContractError("execution contract environment is missing")
    reproduced = _build_receipt(
        lane_id,
        Path(source["path"]),
        environment,
        contract_path,
    )
    if payload != reproduced:
        raise ExecutionContractError("execution contract receipt is not canonical producer output")
    digest = payload.get("semantic_contract_sha256")
    if not isinstance(digest, str) or SHA256_RE.fullmatch(digest) is None:
        raise ExecutionContractError("execution contract semantic digest is invalid")
    return payload, digest


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", default=str(DEFAULT_CONTRACT))
    parser.add_argument("--lane-id", required=True)
    parser.add_argument("--source-list", required=True)
    parser.add_argument("--environment-json", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    output = Path(args.output).resolve()
    if output.exists():
        output.unlink()
    try:
        receipt = build_receipt(
            args.lane_id,
            Path(args.source_list),
            Path(args.environment_json),
            Path(args.contract),
        )
        _write_json(output, receipt)
    except (ExecutionContractError, OSError, ValueError) as exc:
        print(json.dumps({"status": "FAIL", "error": str(exc)}, sort_keys=True))
        return 2
    print(
        json.dumps(
            {"status": "PASS", "output": str(output), "sha256": _file_sha256(output)},
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
