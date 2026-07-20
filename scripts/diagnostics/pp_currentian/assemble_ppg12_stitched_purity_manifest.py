#!/usr/bin/env python3
"""Assemble deterministic inputs for the PPG12 stitched-purity closure gate.

This helper is intentionally ROOT-independent.  Upstream extractors write one
compact JSON payload per physics lane; this program verifies those payloads,
their source-file hashes, and their complete bin coverage before assembling
the manifests consumed by ``ppg12_stitched_purity_closure_gate.py``.

Commands
--------
lane-index
    Freeze exactly 32 lane-extract files into a hash-linked index.
manifest
    Assemble a reference, candidate, or production manifest from an exact
    lane index (or directly supplied lane files), provenance, purity, and—when
    applicable—candidate-parity evidence.
merge-audit
    Canonically bind the exact 20 inclusive and 12 photon merge inputs in
    deterministic lane order to the two family arithmetic summaries.  This
    command re-resolves every lane ``merge_input`` and verifies its file hash;
    hand-authored merge summaries are not sufficient admission evidence.
production-wrapper
    Bind a passing admission, full-production manifests, two candidate ROOT
    artifacts, and the historical comparison into the verification wrapper.

The program never submits jobs, merges ROOT files, changes current pointers,
or promotes artifacts.  It removes a stale requested output before validation
and emits no replacement unless every assembly invariant passes.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import re
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any, Iterable


REPO = Path(__file__).resolve().parents[3]
DEFAULT_CONTRACT = (
    REPO / "agent_context/analysis_contracts/ppg12_stitched_purity_closure.yaml"
)
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
_PURITY_REPRODUCTION_CACHE: dict[str, dict[str, Any]] = {}
_LANE_REPRODUCTION_CACHE: dict[str, dict[str, Any]] = {}
REQUIRED_HISTORICAL_SOURCE_ROLES = {
    "historical_purity",
    "historical_abcd",
    "historical_leakage",
    "candidate_purity",
    "candidate_inclusive",
    "candidate_photon",
}
LANE_DIRECT_EVIDENCE_TO_FIELD = {
    "source_list": "source_list_sha256",
    "event_set": "event_set_sha256",
}
LANE_EVIDENCE_KEYS = {"source_list", "event_set", "config"}
RUNTIME_HASH_FIELDS = {
    "config_sha256",
    "reconstruction_sha256",
    "model_set_sha256",
    "ownership_sha256",
    "weight_sha256",
    "estimator_sha256",
}


class AssemblyError(RuntimeError):
    """An input violates the fail-closed assembly contract."""


def _canonical_bytes(payload: Any) -> bytes:
    return json.dumps(
        payload, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode()


def _payload_sha256(payload: Any) -> str:
    return hashlib.sha256(_canonical_bytes(payload)).hexdigest()


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
    except FileNotFoundError as exc:
        raise AssemblyError(f"missing input file: {path}") from exc
    return digest.hexdigest()


def _read_json(path: Path, label: str) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text())
    except FileNotFoundError as exc:
        raise AssemblyError(f"missing {label}: {path}") from exc
    except json.JSONDecodeError as exc:
        raise AssemblyError(f"invalid JSON in {label} {path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise AssemblyError(f"{label} must contain a JSON object: {path}")
    return payload


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    temporary.replace(path)


def _prepare_output(path: Path) -> None:
    if path.exists():
        path.unlink()


def _resolve_path(raw: str, parent: Path) -> Path:
    path = Path(raw).expanduser()
    return path.resolve() if path.is_absolute() else (parent / path).resolve()


def _is_number(value: Any) -> bool:
    return (
        isinstance(value, (int, float))
        and not isinstance(value, bool)
        and math.isfinite(float(value))
    )


def _require_sha256(value: Any, label: str) -> str:
    if not isinstance(value, str) or SHA256_RE.fullmatch(value) is None:
        raise AssemblyError(f"{label} must be a lowercase 64-character SHA-256")
    return value


def _require_schema(payload: dict[str, Any], schema: str, label: str) -> None:
    if payload.get("schema") != schema:
        raise AssemblyError(
            f"unsupported {label} schema: {payload.get('schema')!r}; expected {schema}"
        )


def _read_contract(path: Path) -> dict[str, Any]:
    contract = _read_json(path, "closure contract")
    _require_schema(
        contract, "ppg12-stitched-purity-closure-contract/v1", "closure contract"
    )
    return contract


def _lane_id(family: str, sample: str, period: str, interaction: str) -> str:
    return f"{family}:{sample}:{period}:{interaction}"


def _expected_lanes(contract: dict[str, Any]) -> dict[str, dict[str, str]]:
    expected: dict[str, dict[str, str]] = {}
    for family, spec in contract["families"].items():
        for sample in spec["samples"]:
            for period in contract["periods"]:
                for interaction in contract["interactions"]:
                    identity = {
                        "family": family,
                        "sample": sample,
                        "period": period,
                        "interaction": interaction,
                    }
                    expected[_lane_id(**identity)] = identity
    return expected


def _expected_candidate_parity_lanes(contract: dict[str, Any]) -> set[str]:
    return {
        _lane_id("photon", sample, period, interaction)
        for sample in contract["families"]["photon"]["samples"]
        for period in contract["periods"]
        for interaction in contract["interactions"]
    }


def _finite_vector(value: Any, label: str, *, nonnegative: bool = False) -> list[float]:
    if not isinstance(value, list) or not value:
        raise AssemblyError(f"{label} must be a non-empty list")
    if not all(_is_number(item) for item in value):
        raise AssemblyError(f"{label} contains a non-finite or non-numeric value")
    result = [float(item) for item in value]
    if nonnegative and any(item < 0.0 for item in result):
        raise AssemblyError(f"{label} contains a negative value")
    return result


def _validate_histogram(payload: Any, label: str) -> list[float]:
    if not isinstance(payload, dict):
        raise AssemblyError(f"{label} must be an object")
    edges = _finite_vector(payload.get("bin_edges"), f"{label}.bin_edges")
    sumw = _finite_vector(payload.get("sumw"), f"{label}.sumw", nonnegative=True)
    sumw2 = _finite_vector(payload.get("sumw2"), f"{label}.sumw2", nonnegative=True)
    fills = _finite_vector(payload.get("fills"), f"{label}.fills", nonnegative=True)
    if len(edges) != len(sumw) + 1:
        raise AssemblyError(f"{label} has inconsistent bin-edge and cell counts")
    if len(sumw2) != len(sumw) or len(fills) != len(sumw):
        raise AssemblyError(f"{label} has inconsistent sumw/sumw2/fills lengths")
    if any(right <= left for left, right in zip(edges, edges[1:])):
        raise AssemblyError(f"{label} bin edges are not strictly increasing")
    return edges


def _validate_file_link(link: Any, parent: Path, label: str) -> dict[str, str]:
    if not isinstance(link, dict):
        raise AssemblyError(f"{label} must be an object with path and sha256")
    raw_path = link.get("path")
    if not isinstance(raw_path, str) or not raw_path:
        raise AssemblyError(f"{label}.path is required")
    expected_hash = _require_sha256(link.get("sha256"), f"{label}.sha256")
    path = _resolve_path(raw_path, parent)
    observed_hash = _file_sha256(path)
    if observed_hash != expected_hash:
        raise AssemblyError(
            f"{label} hash mismatch for {path}: expected {expected_hash}, observed {observed_hash}"
        )
    return {"path": str(path), "sha256": observed_hash}


def _canonical_tool_link(
    contract: dict[str, Any], key: str, declared: Any, parent: Path, label: str
) -> dict[str, str]:
    expected_path = (REPO / contract["canonical_tools"][key]).resolve()
    link = _validate_file_link(declared, parent, label)
    if Path(link["path"]) != expected_path:
        raise AssemblyError(
            f"{label} must resolve to canonical {key}: {expected_path}"
        )
    return link


def _load_python_module(path: Path, name: str) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise AssemblyError(f"cannot load canonical evidence checker: {path}")
    module = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(module)
    except Exception as exc:
        raise AssemblyError(f"cannot import canonical evidence checker {path}: {exc}") from exc
    return module


def _text_rows(path: Path, label: str) -> list[str]:
    try:
        return [
            line.strip()
            for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ]
    except UnicodeDecodeError as exc:
        raise AssemblyError(f"{label} is not a text row list: {path}") from exc


def _root_capable_python() -> Path:
    candidates = [
        Path(sys.executable).resolve(),
        (REPO.parent / "analysis/env/bin/python3").resolve(),
    ]
    for candidate in candidates:
        if not candidate.exists():
            continue
        probe = subprocess.run(
            [str(candidate), "-c", "import ROOT, uproot"],
            text=True,
            capture_output=True,
            check=False,
        )
        if probe.returncode == 0:
            return candidate
    raise AssemblyError(
        "ROOT/uproot-capable Python is required to re-run canonical evidence producers"
    )


def _reproduce_lane_extract(
    extractor_path: Path,
    input_sidecar: dict[str, str],
    contract: dict[str, Any],
) -> dict[str, Any]:
    """Re-run the canonical extractor; linked JSON is never accepted on trust."""
    cache_key = _payload_sha256(
        {
            "extractor_sha256": _file_sha256(extractor_path),
            "input_sidecar_sha256": input_sidecar["sha256"],
            "contract": contract,
        }
    )
    cached = _LANE_REPRODUCTION_CACHE.get(cache_key)
    if cached is not None:
        return json.loads(json.dumps(cached))
    with tempfile.TemporaryDirectory(prefix="ppg12-lane-replay-") as temporary:
        temporary_path = Path(temporary)
        contract_path = temporary_path / "contract.json"
        output_path = temporary_path / "lane.json"
        _write_json(contract_path, contract)
        command = [
            str(_root_capable_python()),
            str(extractor_path),
            "--metadata-json",
            input_sidecar["path"],
            "--contract",
            str(contract_path),
            "--output",
            str(output_path),
        ]
        replay = subprocess.run(command, capture_output=True, text=True, check=False)
        if replay.returncode != 0 or not output_path.exists():
            detail = (replay.stderr or replay.stdout).strip()
            raise AssemblyError(f"canonical lane extraction replay failed: {detail}")
        reproduced = _read_json(output_path, "reproduced lane extraction")
    _LANE_REPRODUCTION_CACHE[cache_key] = reproduced
    return json.loads(json.dumps(reproduced))


def _reproduce_purity_evidence(
    contract_path: Path,
    contract: dict[str, Any],
    lanes: list[dict[str, Any]],
    lane_links: list[dict[str, str]],
) -> dict[str, Any]:
    """Replay the canonical estimator for the exact scientific lane payload."""
    normalized_lane_links = sorted(lane_links, key=lambda row: row["lane_id"])

    def materialize(template: dict[str, Any]) -> dict[str, Any]:
        result = json.loads(json.dumps(template))
        result["inputs"] = {
            "contract": {
                "path": str(contract_path.resolve()),
                "sha256": _payload_sha256(contract),
            },
            "lane_index": None,
            "lanes": normalized_lane_links,
            "lane_set_sha256": _payload_sha256(normalized_lane_links),
        }
        return result

    scientific_lanes = [
        {
            "lane_id": lane["lane_id"],
            "family": lane["family"],
            "observables": lane["observables"],
        }
        for lane in sorted(lanes, key=lambda row: row["lane_id"])
    ]
    producer_path = (REPO / contract["canonical_tools"]["purity_producer"]).resolve()
    estimator_path = (
        REPO / contract["canonical_tools"]["ppg12_estimator_source"]
    ).resolve()
    cache_key = _payload_sha256(
        {
            "contract": contract,
            "scientific_lanes": scientific_lanes,
            "producer_sha256": _file_sha256(producer_path),
            "estimator_sha256": _file_sha256(estimator_path),
        }
    )
    cached = _PURITY_REPRODUCTION_CACHE.get(cache_key)
    if cached is not None:
        return materialize(cached)

    with tempfile.TemporaryDirectory(prefix="ppg12-purity-replay-") as temporary:
        output = Path(temporary) / "purity.json"
        command = [
            str(_root_capable_python()),
            str(producer_path),
            "--contract",
            str(contract_path),
            "purity-repeat",
        ]
        for link in sorted(lane_links, key=lambda row: row["lane_id"]):
            command.extend(["--lane-json", link["path"]])
        command.extend(["--output", str(output)])
        replay = subprocess.run(command, capture_output=True, text=True)
        if replay.returncode != 0 or not output.exists():
            detail = (replay.stderr or replay.stdout).strip()
            raise AssemblyError(f"canonical purity replay failed: {detail}")
        payload = _read_json(output, "reproduced purity evidence")
    template = {key: value for key, value in payload.items() if key != "inputs"}
    _PURITY_REPRODUCTION_CACHE[cache_key] = template
    return materialize(template)


def _validate_merge_input(
    payload: Any, lane_path: Path, lane_id: str
) -> dict[str, str]:
    if not isinstance(payload, dict):
        raise AssemblyError(f"{lane_id}.merge_input is required for candidate output")
    input_id = payload.get("input_id")
    if input_id != lane_id:
        raise AssemblyError(
            f"{lane_id}.merge_input.input_id must equal its canonical lane_id"
        )
    link = _validate_file_link(payload, lane_path.parent, f"{lane_id}.merge_input")
    return {"input_id": lane_id, **link}


def _validate_lane_extraction(
    payload: Any,
    lane_path: Path,
    lane_id: str,
    family: str,
    contract: dict[str, Any],
    lane_payload: dict[str, Any],
) -> dict[str, Any]:
    if not isinstance(payload, dict):
        raise AssemblyError(f"{lane_id}.extraction is required")
    _require_schema(
        payload,
        "ppg12-stitched-purity-lane-extraction/v1",
        f"{lane_id}.extraction",
    )
    extractor_link = _canonical_tool_link(
        contract,
        "lane_extractor",
        payload.get("extractor"),
        lane_path.parent,
        f"{lane_id}.extraction.extractor",
    )
    extractor = _load_python_module(
        Path(extractor_link["path"]), "ppg12_stitched_purity_lane_extractor"
    )
    input_sidecar = _validate_file_link(
        payload.get("input_sidecar"), lane_path.parent, f"{lane_id}.input_sidecar"
    )
    root = _validate_file_link(payload.get("root"), lane_path.parent, f"{lane_id}.root")
    fill_evidence = _validate_file_link(
        payload.get("fill_evidence"), lane_path.parent, f"{lane_id}.fill_evidence"
    )
    fill_payload = _read_json(Path(fill_evidence["path"]), f"{lane_id} fill evidence")
    if (
        fill_payload.get("schema") != "ppg12-stitched-purity-fill-evidence/v1"
        or fill_payload.get("lane_id") != lane_id
        or fill_payload.get("root_sha256") != root["sha256"]
    ):
        raise AssemblyError(f"{lane_id} fill evidence is not bound to its exact ROOT")

    evidence = payload.get("evidence")
    if not isinstance(evidence, dict) or set(evidence) != LANE_EVIDENCE_KEYS:
        raise AssemblyError(f"{lane_id}.extraction evidence coverage is not exact")
    normalized_evidence: dict[str, dict[str, str]] = {}
    for name in sorted(LANE_EVIDENCE_KEYS):
        link = _validate_file_link(
            evidence[name], lane_path.parent, f"{lane_id}.extraction.evidence.{name}"
        )
        field = LANE_DIRECT_EVIDENCE_TO_FIELD.get(name)
        if field is not None and lane_payload.get(field) != link["sha256"]:
            raise AssemblyError(f"{lane_id}.{field} is not bound to its evidence file")
        normalized_evidence[name] = link
    runtime_payload = payload.get("runtime_evidence")
    if not isinstance(runtime_payload, dict):
        raise AssemblyError(f"{lane_id}.runtime_evidence is required")
    try:
        runtime_evidence, runtime_source_sets, runtime_hashes = (
            extractor._validate_runtime_manifest(
                runtime_payload.get("manifest"),
                lane_path.parent,
                normalized_evidence["config"],
                contract,
            )
        )
    except Exception as exc:
        raise AssemblyError(f"{lane_id} runtime evidence is invalid: {exc}") from exc
    if runtime_payload != runtime_evidence:
        raise AssemblyError(f"{lane_id} runtime evidence differs from the executed manifest")
    if payload.get("runtime_source_sets") != runtime_source_sets:
        raise AssemblyError(f"{lane_id} runtime source sets are not canonical")
    for field in RUNTIME_HASH_FIELDS:
        if lane_payload.get(field) != runtime_hashes.get(field):
            raise AssemblyError(f"{lane_id}.{field} is not derived from executed runtime roles")

    raw_candidate_links = payload.get("candidate_parity_evidence")
    if not isinstance(raw_candidate_links, list):
        raise AssemblyError(f"{lane_id}.candidate_parity_evidence must be a list")
    candidate_links: list[dict[str, str]] = []
    roles: set[str] = set()
    for index, row in enumerate(raw_candidate_links):
        if not isinstance(row, dict) or not isinstance(row.get("role"), str):
            raise AssemblyError(f"{lane_id} candidate trace {index} is malformed")
        role = row["role"]
        if role in roles:
            raise AssemblyError(f"{lane_id} duplicates candidate trace role {role}")
        roles.add(role)
        candidate_links.append(
            {
                "role": role,
                **_validate_file_link(
                    row, lane_path.parent, f"{lane_id}.candidate_parity_evidence.{role}"
                ),
            }
        )

    group_contract = contract["admission_group_contract"]
    group_count = lane_payload.get("group_count")
    group_size = lane_payload.get("group_size")
    group_index_start = lane_payload.get("group_index_start")
    if (
        not isinstance(group_count, int)
        or isinstance(group_count, bool)
        or group_count < int(group_contract["minimum_group_count"])
        or group_size != int(group_contract["group_size"])
        or group_index_start != 0
    ):
        raise AssemblyError(
            f"{lane_id} must begin with deterministic five-file group zero"
        )
    event_rows = _text_rows(
        Path(normalized_evidence["event_set"]["path"]), f"{lane_id}.event_set"
    )
    expected_event_rows = group_count * group_size
    if (
        lane_payload.get("event_set_row_count") != expected_event_rows
        or len(event_rows) != expected_event_rows
    ):
        raise AssemblyError(
            f"{lane_id} event set does not contain exactly {expected_event_rows} rows"
        )
    try:
        canonical_groups = extractor._canonical_groups(
            lane_id, event_rows, group_size, group_index_start
        )
    except Exception as exc:
        raise AssemblyError(f"{lane_id} group evidence is invalid: {exc}") from exc
    group_set_sha256 = _payload_sha256(canonical_groups)
    if (
        payload.get("groups") != canonical_groups
        or lane_payload.get("groups") != canonical_groups
        or payload.get("group_set_sha256") != group_set_sha256
        or lane_payload.get("group_set_sha256") != group_set_sha256
    ):
        raise AssemblyError(f"{lane_id} deterministic group identities/hashes drifted")
    try:
        extractor._require_exact_source_slices(
            lane_id,
            Path(normalized_evidence["source_list"]["path"]),
            event_rows,
            canonical_groups,
            group_size,
        )
    except Exception as exc:
        raise AssemblyError(
            f"{lane_id} event set is not an exact production-list slice: {exc}"
        ) from exc
    if Path(root["path"]).stat().st_size < int(group_contract["minimum_root_bytes"]):
        raise AssemblyError(f"{lane_id} ROOT is below the terminal canary byte floor")

    raw_group_sidecars = payload.get("group_input_sidecars")
    if not isinstance(raw_group_sidecars, list) or len(raw_group_sidecars) != group_count:
        raise AssemblyError(
            f"{lane_id} requires one canonical input sidecar per five-file group"
        )
    group_sidecars: list[dict[str, Any]] = []
    for index, (row, group) in enumerate(zip(raw_group_sidecars, canonical_groups)):
        if (
            not isinstance(row, dict)
            or row.get("group_index") != group["group_index"]
            or row.get("group_id") != group["group_id"]
        ):
            raise AssemblyError(
                f"{lane_id}.group_input_sidecars[{index}] has a stale group identity"
            )
        group_sidecars.append(
            {
                "group_index": group["group_index"],
                "group_id": group["group_id"],
                **_validate_file_link(
                    row,
                    lane_path.parent,
                    f"{lane_id}.group_input_sidecars[{index}]",
                ),
            }
        )

    reproduced = _reproduce_lane_extract(
        Path(extractor_link["path"]), input_sidecar, contract
    )
    if reproduced != lane_payload:
        raise AssemblyError(
            f"{lane_id} is not the exact output of the canonical lane extractor"
        )

    if family == "photon":
        canary = contract["photon_canary_contract"]
        expected_roles = set(canary["required_candidate_trace_roles"])
        if roles != expected_roles:
            raise AssemblyError(
                f"{lane_id} candidate trace roles are not exact; "
                f"expected={sorted(expected_roles)}, observed={sorted(roles)}"
            )
    elif candidate_links:
        raise AssemblyError(f"{lane_id} inclusive lane carries photon candidate traces")

    top_candidate_links = lane_payload.get("candidate_parity_evidence", [])
    if top_candidate_links != candidate_links:
        raise AssemblyError(f"{lane_id} top-level candidate trace links drifted")
    return {
        "schema": "ppg12-stitched-purity-lane-extraction/v1",
        "extractor": extractor_link,
        "input_sidecar": input_sidecar,
        "root": root,
        "fill_evidence": fill_evidence,
        "evidence": normalized_evidence,
        "runtime_evidence": runtime_evidence,
        "runtime_source_sets": runtime_source_sets,
        "candidate_parity_evidence": candidate_links,
        "groups": canonical_groups,
        "group_set_sha256": group_set_sha256,
        "group_input_sidecars": group_sidecars,
        "object_prefix": payload.get("object_prefix"),
        "fill_count_semantics": payload.get("fill_count_semantics"),
    }


def _validate_lane(
    payload: dict[str, Any],
    path: Path,
    contract: dict[str, Any],
    *,
    require_merge_input: bool,
) -> tuple[dict[str, Any], list[float]]:
    _require_schema(payload, "ppg12-stitched-purity-lane/v1", f"lane extract {path}")
    identity_fields = ("family", "sample", "period", "interaction")
    if any(not isinstance(payload.get(field), str) for field in identity_fields):
        raise AssemblyError(f"lane extract has incomplete identity: {path}")
    identity = {field: payload[field] for field in identity_fields}
    canonical = _lane_id(**identity)
    if payload.get("lane_id") != canonical:
        raise AssemblyError(
            f"lane_id mismatch in {path}: declared {payload.get('lane_id')!r}, expected {canonical}"
        )
    expected = _expected_lanes(contract)
    if canonical not in expected:
        raise AssemblyError(f"unexpected lane identity: {canonical}")
    if identity != expected[canonical]:
        raise AssemblyError(f"lane identity does not match contract: {canonical}")
    if payload.get("abcd_population") != contract["abcd_population"]:
        raise AssemblyError(f"{canonical} must use unsuffixed inclusive A/B/C/D")
    scale = payload.get("external_scale")
    if not _is_number(scale) or abs(float(scale) - 1.0) > 1e-12:
        raise AssemblyError(f"{canonical} uses a forbidden external scale: {scale!r}")
    group_count = payload.get("group_count")
    if not isinstance(group_count, int) or isinstance(group_count, bool) or group_count <= 0:
        raise AssemblyError(f"{canonical}.group_count must be a positive integer")

    lane: dict[str, Any] = {
        "schema": payload["schema"],
        "lane_id": canonical,
        **identity,
        "group_count": group_count,
        "group_size": payload.get("group_size"),
        "group_index_start": payload.get("group_index_start"),
        "event_set_row_count": payload.get("event_set_row_count"),
        "groups": payload.get("groups"),
        "abcd_population": contract["abcd_population"],
        "external_scale": 1.0,
    }
    for field in contract["paired_lane_fields"]:
        lane[field] = _require_sha256(payload.get(field), f"{canonical}.{field}")

    required = set(contract["families"][identity["family"]]["required_observables"])
    observables = payload.get("observables")
    if not isinstance(observables, dict):
        raise AssemblyError(f"{canonical}.observables must be an object")
    observed = set(observables)
    if observed != required:
        missing = sorted(required - observed)
        unexpected = sorted(observed - required)
        raise AssemblyError(
            f"{canonical} observable coverage is not exact; missing={missing}, unexpected={unexpected}"
        )
    common_edges: list[float] | None = None
    lane["observables"] = {}
    for name in sorted(required):
        edges = _validate_histogram(observables[name], f"{canonical}:{name}")
        if common_edges is None:
            common_edges = edges
        elif edges != common_edges:
            raise AssemblyError(f"{canonical} observables do not share one binning")
        lane["observables"][name] = observables[name]
    assert common_edges is not None

    merge_input = payload.get("merge_input")
    if require_merge_input:
        lane["merge_input"] = _validate_merge_input(merge_input, path, canonical)
    elif merge_input is not None:
        lane["merge_input"] = _validate_merge_input(merge_input, path, canonical)
    extraction = _validate_lane_extraction(
        payload.get("extraction"),
        path,
        canonical,
        identity["family"],
        contract,
        payload,
    )
    if require_merge_input and lane.get("merge_input") != {
        "input_id": canonical,
        **extraction["root"],
    }:
        raise AssemblyError(
            f"{canonical}.merge_input is not the exact ROOT extracted for that lane"
        )
    lane["extraction"] = extraction
    if extraction["candidate_parity_evidence"]:
        lane["candidate_parity_evidence"] = extraction["candidate_parity_evidence"]
    return lane, common_edges


def _load_lane_links(
    *,
    index_path: Path | None,
    lane_paths: Iterable[Path],
    contract: dict[str, Any],
) -> list[dict[str, str]]:
    direct = list(lane_paths)
    if index_path is not None and direct:
        raise AssemblyError("use either --lane-index or --lane-json, not both")
    if index_path is None and not direct:
        raise AssemblyError("one --lane-index or exactly 32 --lane-json inputs are required")
    if index_path is not None:
        payload = _read_json(index_path, "lane index")
        _require_schema(payload, "ppg12-stitched-purity-lane-index/v1", "lane index")
        if payload.get("contract_sha256") != _payload_sha256(contract):
            raise AssemblyError("lane index was created from a different closure contract")
        raw_links = payload.get("lanes")
        if not isinstance(raw_links, list):
            raise AssemblyError("lane index .lanes must be a list")
        normalized: list[dict[str, str]] = []
        for index, link in enumerate(raw_links):
            if not isinstance(link, dict) or not isinstance(link.get("lane_id"), str):
                raise AssemblyError(f"lane index entry {index}.lane_id is required")
            normalized.append(
                {
                    "lane_id": link["lane_id"],
                    **_validate_file_link(
                        link, index_path.parent, f"lane index entry {index}"
                    ),
                }
            )
        if payload.get("lane_set_sha256") != _payload_sha256(normalized):
            raise AssemblyError("lane index lane-set hash is stale or malformed")
        return normalized
    links = []
    for path in direct:
        resolved = path.resolve()
        links.append({"path": str(resolved), "sha256": _file_sha256(resolved)})
    return links


def _load_and_validate_lanes(
    links: list[dict[str, str]],
    contract: dict[str, Any],
    *,
    require_merge_input: bool,
) -> tuple[list[dict[str, Any]], list[dict[str, str]], list[float]]:
    expected = _expected_lanes(contract)
    if len(links) != len(expected):
        raise AssemblyError(
            f"exactly {len(expected)} lane files are required; observed {len(links)}"
        )
    lanes: dict[str, dict[str, Any]] = {}
    normalized_links: dict[str, dict[str, str]] = {}
    global_edges: list[float] | None = None
    for index, link in enumerate(links):
        path = Path(link["path"])
        actual_hash = _file_sha256(path)
        expected_hash = _require_sha256(link.get("sha256"), f"lane link {index}.sha256")
        if actual_hash != expected_hash:
            raise AssemblyError(f"lane file changed after indexing: {path}")
        payload = _read_json(path, f"lane extract {index}")
        lane, edges = _validate_lane(
            payload, path, contract, require_merge_input=require_merge_input
        )
        lane_id = lane["lane_id"]
        declared_link_lane = link.get("lane_id")
        if declared_link_lane is not None and declared_link_lane != lane_id:
            raise AssemblyError(
                f"lane index identity mismatch for {path}: declared {declared_link_lane}, observed {lane_id}"
            )
        if lane_id in lanes:
            raise AssemblyError(f"duplicate lane extract: {lane_id}")
        lanes[lane_id] = lane
        normalized_links[lane_id] = {
            "lane_id": lane_id,
            "path": str(path.resolve()),
            "sha256": actual_hash,
        }
        if global_edges is None:
            global_edges = edges
        elif edges != global_edges:
            raise AssemblyError(
                f"lane {lane_id} does not use the shared stitched-purity binning"
            )
    missing = sorted(set(expected) - set(lanes))
    unexpected = sorted(set(lanes) - set(expected))
    if missing or unexpected:
        raise AssemblyError(
            f"lane set is incomplete; missing={missing}, unexpected={unexpected}"
        )
    if require_merge_input:
        merge_rows = [lanes[lane_id]["merge_input"] for lane_id in sorted(lanes)]
        paths = [row["path"] for row in merge_rows]
        hashes = [row["sha256"] for row in merge_rows]
        if len(set(paths)) != len(paths):
            raise AssemblyError(
                "physical lane merge inputs must use 32 unique resolved ROOT paths"
            )
        if len(set(hashes)) != len(hashes):
            raise AssemblyError(
                "physical lane merge inputs must use 32 unique ROOT content hashes"
            )
    assert global_edges is not None
    ordered_ids = sorted(expected)
    return (
        [lanes[lane_id] for lane_id in ordered_ids],
        [normalized_links[lane_id] for lane_id in ordered_ids],
        global_edges,
    )


def _numbers_close(left: float, right: float) -> bool:
    return math.isclose(float(left), float(right), rel_tol=1e-12, abs_tol=1e-12)


def _canonical_group_extracts(
    lanes: list[dict[str, Any]], contract: dict[str, Any]
) -> list[dict[str, Any]]:
    """Replay every independently produced five-file group and close its aggregate."""
    extractor_path = (REPO / contract["canonical_tools"]["lane_extractor"]).resolve()
    stable_fields = (
        "config_sha256",
        "reconstruction_sha256",
        "model_set_sha256",
        "ownership_sha256",
        "weight_sha256",
        "estimator_sha256",
    )
    group_rows: list[dict[str, Any]] = []
    for lane in sorted(lanes, key=lambda row: row["lane_id"]):
        expected_groups = {row["group_index"]: row for row in lane["groups"]}
        components: list[dict[str, Any]] = []
        for link in lane["extraction"]["group_input_sidecars"]:
            reproduced = _reproduce_lane_extract(extractor_path, link, contract)
            group_index = link["group_index"]
            expected_group = expected_groups.get(group_index)
            if (
                expected_group is None
                or reproduced.get("schema") != "ppg12-stitched-purity-lane/v1"
                or reproduced.get("lane_id") != lane["lane_id"]
                or reproduced.get("family") != lane["family"]
                or reproduced.get("sample") != lane["sample"]
                or reproduced.get("period") != lane["period"]
                or reproduced.get("interaction") != lane["interaction"]
                or reproduced.get("group_count") != 1
                or reproduced.get("group_size") != lane["group_size"]
                or reproduced.get("group_index_start") != group_index
                or reproduced.get("event_set_row_count") != lane["group_size"]
                or reproduced.get("groups") != [expected_group]
                or reproduced.get("abcd_population") != contract["abcd_population"]
                or reproduced.get("external_scale") != 1.0
            ):
                raise AssemblyError(
                    f"{lane['lane_id']} group {group_index} is not its canonical five-file extract"
                )
            if any(reproduced.get(field) != lane.get(field) for field in stable_fields):
                raise AssemblyError(
                    f"{lane['lane_id']} group {group_index} changed a frozen execution contract"
                )
            extraction = reproduced.get("extraction")
            if not isinstance(extraction, dict):
                raise AssemblyError(
                    f"{lane['lane_id']} group {group_index} lacks extraction evidence"
                )
            reproduced_sidecars = extraction.get("group_input_sidecars")
            if reproduced_sidecars != [link]:
                raise AssemblyError(
                    f"{lane['lane_id']} group {group_index} sidecar self-link drifted"
                )
            observables = reproduced.get("observables")
            required = set(contract["families"][lane["family"]]["required_observables"])
            if not isinstance(observables, dict) or set(observables) != required:
                raise AssemblyError(
                    f"{lane['lane_id']} group {group_index} observable coverage is not exact"
                )
            for name in sorted(required):
                _validate_histogram(
                    observables[name], f"{lane['lane_id']}:group-{group_index}:{name}"
                )
            component = {
                "lane_id": lane["lane_id"],
                "family": lane["family"],
                "group_index": group_index,
                "group_id": expected_group["group_id"],
                "group_sha256": expected_group["group_sha256"],
                "input_sidecar": {"path": link["path"], "sha256": link["sha256"]},
                "root": extraction["root"],
                "fill_evidence": extraction["fill_evidence"],
                "extraction_payload_sha256": _payload_sha256(reproduced),
                "observables": observables,
            }
            components.append(component)
            group_rows.append(component)

        if sorted(row["group_index"] for row in components) != list(
            range(lane["group_count"])
        ):
            raise AssemblyError(
                f"{lane['lane_id']} component extracts do not cover its exact group prefix"
            )
        for name, aggregate in lane["observables"].items():
            pieces = [row["observables"][name] for row in components]
            if any(piece["bin_edges"] != aggregate["bin_edges"] for piece in pieces):
                raise AssemblyError(
                    f"{lane['lane_id']} group binning differs for {name}"
                )
            for field in ("sumw", "sumw2", "fills"):
                expected = [
                    sum(float(piece[field][index]) for piece in pieces)
                    for index in range(len(aggregate[field]))
                ]
                if any(
                    not _numbers_close(observed, wanted)
                    for observed, wanted in zip(aggregate[field], expected)
                ):
                    raise AssemblyError(
                        f"{lane['lane_id']} aggregate {name}.{field} is not the sum of its groups"
                    )
            aggregate_flow = aggregate.get("flow")
            if not isinstance(aggregate_flow, dict):
                raise AssemblyError(f"{lane['lane_id']} aggregate {name} lacks flow moments")
            for field in ("sumw", "sumw2", "fills"):
                expected_flow = [
                    sum(float(piece["flow"][field][index]) for piece in pieces)
                    for index in range(2)
                ]
                if any(
                    not _numbers_close(observed, wanted)
                    for observed, wanted in zip(aggregate_flow[field], expected_flow)
                ):
                    raise AssemblyError(
                        f"{lane['lane_id']} aggregate {name}.flow.{field} is not the group sum"
                    )
            if not _numbers_close(
                aggregate.get("entries"),
                sum(float(piece["entries"]) for piece in pieces),
            ):
                raise AssemblyError(
                    f"{lane['lane_id']} aggregate {name}.entries is not the group sum"
                )
    return sorted(
        group_rows,
        key=lambda row: (row["group_index"], row["lane_id"], row["group_id"]),
    )


def _required_coverage_cells(
    contract: dict[str, Any], bin_edges: list[float]
) -> list[str]:
    spec = contract["global_minimality_contract"]["required_observables"]
    return [
        f"{family}:{observable}:bin-{index + 1}"
        for family in ("inclusive", "photon")
        for observable in spec[family]
        for index in range(len(bin_edges) - 1)
    ]


def _populated_cells(
    group: dict[str, Any], contract: dict[str, Any], bin_edges: list[float]
) -> set[str]:
    names = contract["global_minimality_contract"]["required_observables"][
        group["family"]
    ]
    populated: set[str] = set()
    for name in names:
        histogram = group["observables"][name]
        if histogram["bin_edges"] != bin_edges:
            raise AssemblyError(
                f"{group['lane_id']} group {group['group_index']} coverage binning drifted"
            )
        for index, (sumw, fills) in enumerate(
            zip(histogram["sumw"], histogram["fills"])
        ):
            if float(sumw) > 0.0 and float(fills) > 0.0:
                populated.add(f"{group['family']}:{name}:bin-{index + 1}")
    return populated


def _build_global_minimality_evidence(
    lanes: list[dict[str, Any]],
    lane_links: list[dict[str, str]],
    bin_edges: list[float],
    contract: dict[str, Any],
) -> dict[str, Any]:
    groups = _canonical_group_extracts(lanes, contract)
    required_cells = _required_coverage_cells(contract, bin_edges)
    required = set(required_cells)
    baseline = [row for row in groups if row["group_index"] == 0]
    extras = [row for row in groups if row["group_index"] > 0]
    if len(baseline) != len(lanes):
        raise AssemblyError("global minimality baseline must contain group zero in all 32 lanes")

    group_sources: list[dict[str, Any]] = []
    populations: dict[tuple[str, int, str], set[str]] = {}
    for group in groups:
        identity = (group["lane_id"], group["group_index"], group["group_id"])
        populated = _populated_cells(group, contract, bin_edges)
        populations[identity] = populated
        group_sources.append(
            {
                "lane_id": group["lane_id"],
                "family": group["family"],
                "group_index": group["group_index"],
                "group_id": group["group_id"],
                "group_sha256": group["group_sha256"],
                "input_sidecar": group["input_sidecar"],
                "root": group["root"],
                "fill_evidence": group["fill_evidence"],
                "extraction_payload_sha256": group["extraction_payload_sha256"],
                "populated_required_cells": sorted(populated),
            }
        )

    covered: set[str] = set()
    baseline_identities: list[dict[str, Any]] = []
    for group in sorted(baseline, key=lambda row: row["lane_id"]):
        identity = (group["lane_id"], group["group_index"], group["group_id"])
        covered.update(populations[identity])
        baseline_identities.append(
            {
                "lane_id": group["lane_id"],
                "group_index": 0,
                "group_id": group["group_id"],
                "group_sha256": group["group_sha256"],
            }
        )
    missing = required - covered
    baseline_missing = sorted(missing)
    additions: list[dict[str, Any]] = []
    for prefix_index, group in enumerate(extras, start=1):
        if not missing:
            raise AssemblyError(
                "global group prefix contains arbitrary extras after coverage became complete"
            )
        identity = (group["lane_id"], group["group_index"], group["group_id"])
        before = sorted(missing)
        newly = missing & populations[identity]
        if not newly:
            raise AssemblyError(
                f"global group prefix addition {group['lane_id']} group {group['group_index']} "
                "does not reduce the missing final-cell set"
            )
        missing -= newly
        additions.append(
            {
                "prefix_index": prefix_index,
                "lane_id": group["lane_id"],
                "group_index": group["group_index"],
                "group_id": group["group_id"],
                "group_sha256": group["group_sha256"],
                "missing_before": before,
                "newly_populated": sorted(newly),
                "missing_after": sorted(missing),
                "covered_cells_sha256": _payload_sha256(sorted(required - missing)),
            }
        )
    if missing:
        raise AssemblyError(
            f"global group prefix is incomplete after all produced groups: {sorted(missing)}"
        )

    payload: dict[str, Any] = {
        "schema": "ppg12-stitched-purity-global-minimality/v1",
        "builder": {
            "path": str(Path(__file__).resolve()),
            "sha256": _file_sha256(Path(__file__).resolve()),
        },
        "contract_sha256": _payload_sha256(contract),
        "lane_links": lane_links,
        "lane_set_sha256": _payload_sha256(lane_links),
        "required_cells": required_cells,
        "group_sources": group_sources,
        "group_source_set_sha256": _payload_sha256(group_sources),
        "baseline": {
            "group_count": len(baseline_identities),
            "groups": baseline_identities,
            "missing_cells": baseline_missing,
            "covered_cells_sha256": _payload_sha256(
                sorted(required - set(baseline_missing))
            ),
        },
        "added_groups": additions,
        "final_missing_cells": [],
        "selected_group_count": len(baseline_identities) + len(additions),
        "stopped_at_first_complete_prefix": True,
    }
    payload["verification_payload_sha256"] = _payload_sha256(payload)
    return payload


def _validate_global_minimality_evidence(
    path: Path,
    lanes: list[dict[str, Any]],
    lane_links: list[dict[str, str]],
    bin_edges: list[float],
    contract: dict[str, Any],
) -> tuple[dict[str, Any], dict[str, str]]:
    payload = _read_json(path, "global minimality evidence")
    unhashed = dict(payload)
    declared_sha = unhashed.pop("verification_payload_sha256", None)
    if declared_sha != _payload_sha256(unhashed):
        raise AssemblyError("global minimality verification payload hash is stale")
    regenerated = _build_global_minimality_evidence(
        lanes, lane_links, bin_edges, contract
    )
    if payload != regenerated:
        raise AssemblyError(
            "global minimality evidence is not canonical replay-derived prefix evidence"
        )
    return payload, {"path": str(path.resolve()), "sha256": _file_sha256(path)}


def _validate_provenance(
    path: Path, contract: dict[str, Any]
) -> tuple[dict[str, str], dict[str, str]]:
    payload = _read_json(path, "provenance")
    _require_schema(
        payload, "ppg12-stitched-purity-provenance/v1", "provenance"
    )
    provenance = payload.get("provenance")
    if not isinstance(provenance, dict):
        raise AssemblyError("provenance.provenance must be an object")
    source_sets = payload.get("source_sets")
    expected_fields = set(contract["frozen_provenance_fields"])
    if not isinstance(source_sets, dict) or set(source_sets) != expected_fields:
        raise AssemblyError("provenance.source_sets must cover every frozen field exactly")
    normalized: dict[str, str] = {}
    for field in contract["frozen_provenance_fields"]:
        source_set = source_sets[field]
        if not isinstance(source_set, dict) or not isinstance(source_set.get("files"), list):
            raise AssemblyError(f"provenance.source_sets.{field}.files is required")
        files = [
            _validate_file_link(
                link, path.parent, f"provenance.source_sets.{field}.files[{index}]"
            )
            for index, link in enumerate(source_set["files"])
        ]
        if not files:
            raise AssemblyError(f"provenance.source_sets.{field}.files must not be empty")
        files.sort(key=lambda row: row["path"])
        if len({row["path"] for row in files}) != len(files):
            raise AssemblyError(f"provenance.source_sets.{field} duplicates a file")
        observed_set_sha = _payload_sha256(files)
        declared_set_sha = _require_sha256(
            source_set.get("set_sha256"), f"provenance.source_sets.{field}.set_sha256"
        )
        declared_provenance = _require_sha256(
            provenance.get(field), f"provenance.{field}"
        )
        if declared_set_sha != observed_set_sha or declared_provenance != observed_set_sha:
            raise AssemblyError(f"provenance.{field} is not the current hash-bound source set")
        normalized[field] = observed_set_sha
    return normalized, {"path": str(path.resolve()), "sha256": _file_sha256(path)}


def _validate_purity(
    path: Path,
    contract: dict[str, Any],
    bin_edges: list[float],
    expected_lane_links: list[dict[str, str]],
    expected_lanes: list[dict[str, Any]],
) -> tuple[dict[str, Any], dict[str, str]]:
    payload = _read_json(path, "purity payload")
    _require_schema(payload, "ppg12-stitched-purity-purity/v1", "purity payload")
    if payload.get("random_seed") != contract["random_seed"]:
        raise AssemblyError("purity payload does not use fixed random seed 42")
    if payload.get("toy_count") != contract["toy_count"]:
        raise AssemblyError("purity payload does not use exactly 20000 toys")
    purity = payload.get("purity")
    if not isinstance(purity, dict):
        raise AssemblyError("purity.purity must be an object")
    edges = _finite_vector(purity.get("bin_edges"), "purity.bin_edges")
    if edges != bin_edges:
        raise AssemblyError("purity and lane-observable binning differ")
    normalized: dict[str, Any] = {"bin_edges": edges}
    for name in ("truth", "raw", "corrected"):
        series = purity.get(name)
        if not isinstance(series, dict):
            raise AssemblyError(f"purity.{name} must be an object")
        values = _finite_vector(series.get("value"), f"purity.{name}.value")
        errors = _finite_vector(
            series.get("error"), f"purity.{name}.error", nonnegative=True
        )
        if len(values) != len(edges) - 1 or len(errors) != len(values):
            raise AssemblyError(f"purity.{name} does not cover every reported bin")
        # The executable estimator fills toys on [-1, 2] and can return a
        # finite fitted mean just outside the physical interval in sparse
        # bins.  Admission tests exact oracle parity; it must not silently
        # clip or reject an otherwise identical executable result.
        normalized[name] = {"value": values, "error": errors}
    repetition = payload.get("fixed_seed_repetition")
    if not isinstance(repetition, dict):
        raise AssemblyError("purity payload lacks fixed_seed_repetition evidence")
    first_sha = _require_sha256(
        repetition.get("first_output_sha256"),
        "fixed_seed_repetition.first_output_sha256",
    )
    repeated_sha = _require_sha256(
        repetition.get("repeated_output_sha256"),
        "fixed_seed_repetition.repeated_output_sha256",
    )
    normalized_sha = _payload_sha256(normalized)
    if first_sha != repeated_sha or first_sha != normalized_sha:
        raise AssemblyError(
            "fixed-seed estimator repetition is non-deterministic or does not bind the purity payload"
        )
    normalized["fixed_seed_repetition"] = {
        "first_output_sha256": first_sha,
        "repeated_output_sha256": repeated_sha,
    }
    diagnostics = payload.get("run_diagnostics")
    if not isinstance(diagnostics, dict):
        raise AssemblyError("purity payload lacks deterministic run diagnostics")
    first_diagnostics = diagnostics.get("first")
    repeated_diagnostics = diagnostics.get("repeated")
    if not isinstance(first_diagnostics, list) or not isinstance(repeated_diagnostics, list):
        raise AssemblyError("purity run diagnostics must contain first/repeated row lists")
    first_diagnostics_sha = _payload_sha256(first_diagnostics)
    repeated_diagnostics_sha = _payload_sha256(repeated_diagnostics)
    if (
        first_diagnostics != repeated_diagnostics
        or diagnostics.get("diagnostics_sha256") != first_diagnostics_sha
        or diagnostics.get("repeated_diagnostics_sha256") != repeated_diagnostics_sha
    ):
        raise AssemblyError("purity run diagnostics are not byte-deterministic")

    algorithm = payload.get("algorithm")
    if not isinstance(algorithm, dict):
        raise AssemblyError("purity payload lacks canonical algorithm provenance")
    _canonical_tool_link(
        contract,
        "ppg12_estimator_source",
        algorithm.get("source"),
        path.parent,
        "purity.algorithm.source",
    )
    _canonical_tool_link(
        contract,
        "purity_producer",
        algorithm.get("producer"),
        path.parent,
        "purity.algorithm.producer",
    )
    if (
        algorithm.get("toys_per_bin") != contract["toy_count"]
        or algorithm.get("random_stream")
        != "one TRandom3(42) stream per full repeated evaluation"
    ):
        raise AssemblyError("purity algorithm seed/toy stream contract drifted")

    inputs = payload.get("inputs")
    if not isinstance(inputs, dict):
        raise AssemblyError("purity payload lacks exact estimator input provenance")
    contract_link = inputs.get("contract")
    if not isinstance(contract_link, dict) or not isinstance(contract_link.get("path"), str):
        raise AssemblyError("purity.inputs.contract link is required")
    linked_contract_path = _resolve_path(contract_link["path"], path.parent)
    linked_contract = _read_contract(linked_contract_path)
    if (
        linked_contract != contract
        or contract_link.get("sha256") != _payload_sha256(contract)
    ):
        raise AssemblyError("purity estimator used a different closure contract")
    raw_lane_links = inputs.get("lanes")
    if not isinstance(raw_lane_links, list):
        raise AssemblyError("purity.inputs.lanes must be the exact 32-lane set")
    normalized_lane_links: list[dict[str, str]] = []
    for index, link in enumerate(raw_lane_links):
        if not isinstance(link, dict) or not isinstance(link.get("lane_id"), str):
            raise AssemblyError(f"purity.inputs.lanes[{index}] is malformed")
        normalized_lane_links.append(
            {
                "lane_id": link["lane_id"],
                **_validate_file_link(
                    link, path.parent, f"purity.inputs.lanes[{index}]"
                ),
            }
        )
    normalized_lane_links.sort(key=lambda row: row["lane_id"])
    expected_links = sorted(expected_lane_links, key=lambda row: row["lane_id"])
    if normalized_lane_links != expected_links:
        raise AssemblyError("purity estimator inputs differ from the manifest lane set")
    if inputs.get("lane_set_sha256") != _payload_sha256(normalized_lane_links):
        raise AssemblyError("purity estimator lane-set hash is stale")
    reproduced = _reproduce_purity_evidence(
        linked_contract_path,
        contract,
        expected_lanes,
        expected_links,
    )
    for field in ("purity", "fixed_seed_repetition", "run_diagnostics", "algorithm"):
        if payload.get(field) != reproduced.get(field):
            raise AssemblyError(
                f"purity {field} does not equal canonical 20000-toy replay"
            )
    return normalized, {"path": str(path.resolve()), "sha256": _file_sha256(path)}


def _rows_by_name(value: Any, label: str) -> dict[str, dict[str, Any]]:
    if not isinstance(value, list):
        raise AssemblyError(f"{label} must be a list")
    result: dict[str, dict[str, Any]] = {}
    for index, row in enumerate(value):
        if not isinstance(row, dict) or not isinstance(row.get("name"), str):
            raise AssemblyError(f"{label}[{index}] has no string name")
        if row["name"] in result:
            raise AssemblyError(f"duplicate {label} entry: {row['name']}")
        result[row["name"]] = row
    return result


def _positive_count(value: Any, label: str) -> int:
    if not isinstance(value, int) or isinstance(value, bool) or value <= 0:
        raise AssemblyError(f"{label} must be a positive integer")
    return value


def _nonnegative_count(value: Any, label: str) -> int:
    if not isinstance(value, int) or isinstance(value, bool) or value < 0:
        raise AssemblyError(f"{label} must be a nonnegative integer")
    return value


def _validate_candidate_parity(
    path: Path, contract: dict[str, Any]
) -> tuple[dict[str, Any], dict[str, str]]:
    payload = _read_json(path, "candidate-parity payload")
    _require_schema(
        payload,
        "ppg12-stitched-purity-candidate-parity/v1",
        "candidate-parity payload",
    )
    parity = payload.get("candidate_parity")
    if not isinstance(parity, dict):
        raise AssemblyError("candidate_parity must be an object")
    trace = payload.get("trace_evidence")
    if not isinstance(trace, dict):
        raise AssemblyError("candidate parity lacks exact raw trace evidence")
    auditor_link = _canonical_tool_link(
        contract,
        "candidate_trace_auditor",
        trace.get("auditor"),
        path.parent,
        "candidate_parity.trace_evidence.auditor",
    )
    builder_link = _canonical_tool_link(
        contract,
        "candidate_summary_builder",
        trace.get("summary_builder"),
        path.parent,
        "candidate_parity.trace_evidence.summary_builder",
    )
    raw_trace_lanes = trace.get("lanes")
    if not isinstance(raw_trace_lanes, list):
        raise AssemblyError("candidate parity trace lanes must be a list")
    trace_rows_by_lane: dict[str, list[dict[str, str]]] = {}
    trace_bundle_by_lane: dict[str, dict[str, Path]] = {}
    normalized_trace_lanes: list[dict[str, Any]] = []
    evidence_paths_by_role: dict[str, set[str]] = {
        role: set()
        for role in (
            "runtime_contract",
            "candidate_csv",
            "executable_aggregate",
            "trace_csv",
            "response_trace_csv",
        )
    }
    evidence_hashes_by_role: dict[str, set[str]] = {
        role: set() for role in evidence_paths_by_role
    }
    auditor = _load_python_module(
        Path(auditor_link["path"]), "ppg12_candidate_trace_auditor"
    )
    builder = _load_python_module(
        Path(builder_link["path"]), "ppg12_candidate_summary_builder"
    )
    for index, row in enumerate(raw_trace_lanes):
        if not isinstance(row, dict) or not isinstance(row.get("lane_id"), str):
            raise AssemblyError(f"candidate parity trace lane {index} is malformed")
        lane_id = row["lane_id"]
        if lane_id in trace_rows_by_lane:
            raise AssemblyError(f"candidate parity trace duplicates lane {lane_id}")
        paired = row.get("paired_evidence")
        if not isinstance(paired, dict) or set(paired) != set(evidence_paths_by_role):
            raise AssemblyError(
                f"candidate parity trace {lane_id} lacks the exact paired-evidence role set"
            )
        normalized_paired: dict[str, dict[str, str]] = {}
        bundle: dict[str, Path] = {}
        for role in sorted(evidence_paths_by_role):
            link = _validate_file_link(
                paired.get(role),
                path.parent,
                f"candidate parity trace {lane_id}.{role}",
            )
            if (
                link["path"] in evidence_paths_by_role[role]
                or link["sha256"] in evidence_hashes_by_role[role]
            ):
                raise AssemblyError(
                    f"candidate parity trace {lane_id} reuses another physical lane's {role}"
                )
            evidence_paths_by_role[role].add(link["path"])
            evidence_hashes_by_role[role].add(link["sha256"])
            normalized_paired[role] = link
            bundle[role] = Path(link["path"])
        try:
            validated_bundle = builder.validate_lane_evidence_bundle(lane_id, bundle)
        except Exception as exc:
            raise AssemblyError(
                f"candidate parity paired evidence failed for {lane_id}: {exc}"
            ) from exc
        csv_path = Path(validated_bundle["paths"]["candidate_csv"])
        rows = validated_bundle["candidate_rows"]
        runtime_contract_sha256 = validated_bundle["contract_sha256"]
        try:
            audit_report = auditor.analyze_candidate_report(
                csv_path,
                expected_lane_id=lane_id,
                expected_runtime_contract_sha256=runtime_contract_sha256,
            )
        except Exception as exc:
            raise AssemblyError(f"candidate parity trace audit failed for {lane_id}: {exc}") from exc
        if audit_report.get("status") != "PASS" or audit_report.get("first_divergence") is not None:
            raise AssemblyError(f"candidate parity raw trace fails executable parity: {lane_id}")
        trace_rows_by_lane[lane_id] = rows
        trace_bundle_by_lane[lane_id] = bundle
        normalized_trace_lanes.append(
            {
                "lane_id": lane_id,
                "paired_evidence": normalized_paired,
                "row_count": len(rows),
                "audit_payload_sha256": _payload_sha256(audit_report),
            }
        )
    normalized_trace_lanes.sort(key=lambda row: row["lane_id"])
    if trace.get("lane_set_sha256") != _payload_sha256(normalized_trace_lanes):
        raise AssemblyError("candidate parity raw-trace set hash is stale")
    tolerances = parity.get("tolerances")
    expected_tolerances = contract["candidate_parity_tolerances"]
    if not isinstance(tolerances, dict) or set(tolerances) != set(expected_tolerances):
        raise AssemblyError("candidate_parity.tolerances does not match the closure contract")
    normalized_tolerances: dict[str, float] = {}
    for name, expected in expected_tolerances.items():
        value = tolerances.get(name)
        if not _is_number(value) or float(value) != float(expected):
            raise AssemblyError(f"candidate_parity tolerance drift: {name}")
        normalized_tolerances[name] = float(value)

    expected_lane_ids = _expected_candidate_parity_lanes(contract)
    if set(trace_rows_by_lane) != expected_lane_ids:
        raise AssemblyError("candidate parity raw trace is not the exact 12 photon lanes")
    lane_rows = parity.get("lane_coverage")
    if not isinstance(lane_rows, list):
        raise AssemblyError("candidate_parity.lane_coverage must be a list")
    indexed_lane_rows: dict[str, dict[str, Any]] = {}
    for index, row in enumerate(lane_rows):
        if not isinstance(row, dict) or not isinstance(row.get("lane_id"), str):
            raise AssemblyError(f"candidate_parity.lane_coverage[{index}] is malformed")
        lane_id = row["lane_id"]
        if lane_id in indexed_lane_rows:
            raise AssemblyError(f"duplicate candidate-parity lane coverage: {lane_id}")
        indexed_lane_rows[lane_id] = row
    if set(indexed_lane_rows) != expected_lane_ids:
        raise AssemblyError("candidate-parity lane coverage is not the exact 12 photon lanes")

    normalized_lane_rows: list[dict[str, Any]] = []
    for lane_id in sorted(expected_lane_ids):
        row = indexed_lane_rows[lane_id]
        lane_population = row.get("population")
        if not isinstance(lane_population, dict):
            raise AssemblyError(f"lane_coverage.{lane_id}.population must be an object")
        features_row = row.get("features")
        scores_row = row.get("scores")
        route_row = row.get("model_route")
        tags_row = row.get("tags")
        if not isinstance(features_row, dict) or not isinstance(scores_row, dict):
            raise AssemblyError(f"lane_coverage.{lane_id} feature/score summaries are required")
        if not isinstance(route_row, dict):
            raise AssemblyError(f"lane_coverage.{lane_id}.model_route must be an object")
        if not isinstance(tags_row, dict) or set(tags_row) != set(contract["required_tag_checks"]):
            raise AssemblyError(f"lane_coverage.{lane_id}.tags coverage is incomplete")
        normalized_lane_rows.append(
            {
                "lane_id": lane_id,
                "population": {
                    "reference_count": _positive_count(
                        lane_population.get("reference_count"),
                        f"lane_coverage.{lane_id}.population.reference_count",
                    ),
                    "candidate_count": _positive_count(
                        lane_population.get("candidate_count"),
                        f"lane_coverage.{lane_id}.population.candidate_count",
                    ),
                    "unmatched_reference": _nonnegative_count(
                        lane_population.get("unmatched_reference"),
                        f"lane_coverage.{lane_id}.population.unmatched_reference",
                    ),
                    "unmatched_candidate": _nonnegative_count(
                        lane_population.get("unmatched_candidate"),
                        f"lane_coverage.{lane_id}.population.unmatched_candidate",
                    ),
                },
                "features": {
                    "comparisons": _positive_count(
                        features_row.get("comparisons"),
                        f"lane_coverage.{lane_id}.features.comparisons",
                    ),
                    "violations": _nonnegative_count(
                        features_row.get("violations"),
                        f"lane_coverage.{lane_id}.features.violations",
                    ),
                    "max_normalized_delta": float(features_row["max_normalized_delta"])
                    if _is_number(features_row.get("max_normalized_delta"))
                    else _raise(
                        f"lane_coverage.{lane_id}.features.max_normalized_delta must be finite"
                    ),
                },
                "scores": {
                    "comparisons": _positive_count(
                        scores_row.get("comparisons"),
                        f"lane_coverage.{lane_id}.scores.comparisons",
                    ),
                    "mismatches": _nonnegative_count(
                        scores_row.get("mismatches"),
                        f"lane_coverage.{lane_id}.scores.mismatches",
                    ),
                    "max_abs_delta": float(scores_row["max_abs_delta"])
                    if _is_number(scores_row.get("max_abs_delta"))
                    else _raise(
                        f"lane_coverage.{lane_id}.scores.max_abs_delta must be finite"
                    ),
                },
                "model_route": {
                    "comparisons": _positive_count(
                        route_row.get("comparisons"),
                        f"lane_coverage.{lane_id}.model_route.comparisons",
                    ),
                    "mismatches": _nonnegative_count(
                        route_row.get("mismatches"),
                        f"lane_coverage.{lane_id}.model_route.mismatches",
                    ),
                },
                "tags": {
                    name: {
                        "comparisons": _positive_count(
                            tags_row[name].get("comparisons")
                            if isinstance(tags_row[name], dict)
                            else None,
                            f"lane_coverage.{lane_id}.tags.{name}.comparisons",
                        ),
                        "mismatches": _nonnegative_count(
                            tags_row[name].get("mismatches")
                            if isinstance(tags_row[name], dict)
                            else None,
                            f"lane_coverage.{lane_id}.tags.{name}.mismatches",
                        ),
                    }
                    for name in contract["required_tag_checks"]
                },
            }
        )
    population = parity.get("population")
    if not isinstance(population, dict):
        raise AssemblyError("candidate_parity.population must be an object")
    normalized: dict[str, Any] = {
        "tolerances": normalized_tolerances,
        "lane_coverage": normalized_lane_rows,
        "population": {
            "reference_count": _positive_count(
                population.get("reference_count"), "population.reference_count"
            ),
            "candidate_count": _positive_count(
                population.get("candidate_count"), "population.candidate_count"
            ),
            "unmatched_reference": _nonnegative_count(
                population.get("unmatched_reference"), "population.unmatched_reference"
            ),
            "unmatched_candidate": _nonnegative_count(
                population.get("unmatched_candidate"), "population.unmatched_candidate"
            ),
        }
    }

    features = _rows_by_name(parity.get("features"), "candidate_parity.features")
    if set(features) != set(contract["required_features"]):
        raise AssemblyError(
            "candidate-parity feature coverage does not match the ordered 11-feature contract"
        )
    normalized["features"] = []
    for name in contract["required_features"]:
        row = features[name]
        normalized["features"].append(
            {
                "name": name,
                "comparisons": _positive_count(
                    row.get("comparisons"), f"features.{name}.comparisons"
                ),
                "violations": _nonnegative_count(
                    row.get("violations"), f"features.{name}.violations"
                ),
                "max_abs_delta": float(row["max_abs_delta"])
                if _is_number(row.get("max_abs_delta"))
                else _raise(f"features.{name}.max_abs_delta must be finite"),
                "max_normalized_delta": float(row["max_normalized_delta"])
                if _is_number(row.get("max_normalized_delta"))
                else _raise(f"features.{name}.max_normalized_delta must be finite"),
            }
        )

    scores = _rows_by_name(parity.get("scores"), "candidate_parity.scores")
    if set(scores) != set(contract["required_score_models"]):
        raise AssemblyError("candidate-parity score coverage must contain baseV3E and base_E")
    normalized["scores"] = []
    for name in contract["required_score_models"]:
        row = scores[name]
        normalized["scores"].append(
            {
                "name": name,
                "comparisons": _positive_count(
                    row.get("comparisons"), f"scores.{name}.comparisons"
                ),
                "mismatches": _nonnegative_count(
                    row.get("mismatches"), f"scores.{name}.mismatches"
                ),
                "max_abs_delta": float(row["max_abs_delta"])
                if _is_number(row.get("max_abs_delta"))
                else _raise(f"scores.{name}.max_abs_delta must be finite"),
            }
        )

    route = parity.get("model_route")
    if not isinstance(route, dict):
        raise AssemblyError("candidate_parity.model_route must be an object")
    normalized["model_route"] = {
        "comparisons": _positive_count(
            route.get("comparisons"), "model_route.comparisons"
        ),
        "mismatches": _nonnegative_count(
            route.get("mismatches"), "model_route.mismatches"
        ),
    }
    tags = parity.get("tags")
    if not isinstance(tags, dict) or set(tags) != set(contract["required_tag_checks"]):
        raise AssemblyError("candidate-parity tag coverage is incomplete")
    normalized["tags"] = {}
    for name in contract["required_tag_checks"]:
        row = tags[name]
        if not isinstance(row, dict):
            raise AssemblyError(f"tags.{name} must be an object")
        normalized["tags"][name] = {
            "comparisons": _positive_count(
                row.get("comparisons"), f"tags.{name}.comparisons"
            ),
            "mismatches": _nonnegative_count(
                row.get("mismatches"), f"tags.{name}.mismatches"
            ),
        }
    try:
        recomputed = builder.build_candidate_parity(
            trace_rows_by_lane,
            created_utc=payload.get("created_utc"),
            trace_bundle_by_lane=trace_bundle_by_lane,
        )["candidate_parity"]
    except Exception as exc:
        raise AssemblyError(f"cannot recompute candidate parity from raw traces: {exc}") from exc
    if recomputed != normalized:
        raise AssemblyError("candidate parity summary does not equal the exact raw traces")
    normalized["trace_evidence"] = {
        "auditor": auditor_link,
        "summary_builder": builder_link,
        "lanes": normalized_trace_lanes,
        "lane_set_sha256": _payload_sha256(normalized_trace_lanes),
    }
    return normalized, {"path": str(path.resolve()), "sha256": _file_sha256(path)}


def _raise(message: str) -> Any:
    raise AssemblyError(message)


def _bind_lane_candidate_traces(
    lanes: list[dict[str, Any]], parity: dict[str, Any]
) -> None:
    trace = parity.get("trace_evidence")
    trace_lanes = trace.get("lanes") if isinstance(trace, dict) else None
    if not isinstance(trace_lanes, list):
        raise AssemblyError("candidate parity normalized trace evidence is missing")
    indexed = {
        row["lane_id"]: row["paired_evidence"]["candidate_csv"]
        for row in trace_lanes
    }
    for lane in lanes:
        if lane["family"] != "photon":
            continue
        links = lane.get("candidate_parity_evidence")
        if not isinstance(links, list) or len(links) != 1:
            raise AssemblyError(f"{lane['lane_id']} lacks its exact candidate-row trace")
        link = links[0]
        expected = indexed.get(lane["lane_id"])
        if link.get("role") != "candidate_rows" or {
            "path": link.get("path"),
            "sha256": link.get("sha256"),
        } != expected:
            raise AssemblyError(
                f"{lane['lane_id']} lane trace differs from candidate-parity evidence"
            )


def _validate_stitched_coverage(
    lanes: list[dict[str, Any]], bin_edges: list[float]
) -> None:
    nbins = len(bin_edges) - 1
    estimator_keys = tuple(
        [f"inclusive:{name}" for name in ("A", "B", "C", "D")]
        + [f"photon:{name}_signal" for name in ("A", "B", "C", "D")]
    )
    sumw_totals = {key: [0.0] * nbins for key in estimator_keys}
    fill_totals = {key: [0.0] * nbins for key in estimator_keys}
    for lane in lanes:
        family = lane["family"]
        names = (
            ("A", "B", "C", "D")
            if family == "inclusive"
            else ("A_signal", "B_signal", "C_signal", "D_signal")
        )
        for name in names:
            key = f"{family}:{name}"
            histogram = lane["observables"][name]
            sumw_totals[key] = [
                left + float(right)
                for left, right in zip(sumw_totals[key], histogram["sumw"])
            ]
            fill_totals[key] = [
                left + float(right)
                for left, right in zip(fill_totals[key], histogram["fills"])
            ]
    for key in estimator_keys:
        missing_sumw = [
            index + 1 for index, value in enumerate(sumw_totals[key]) if value <= 0.0
        ]
        missing_fills = [
            index + 1 for index, value in enumerate(fill_totals[key]) if value <= 0.0
        ]
        if missing_sumw or missing_fills:
            raise AssemblyError(
                f"stitched estimator coverage is empty for {key}; "
                f"sumw_bins={missing_sumw}, fill_bins={missing_fills}"
            )


def _frozen_snapshot(manifest: dict[str, Any], contract: dict[str, Any]) -> dict[str, Any]:
    lane_fields = contract["production_frozen_lane_fields"]
    return {
        "provenance": {
            field: manifest["provenance"].get(field)
            for field in contract["frozen_provenance_fields"]
        },
        "lanes": {
            lane["lane_id"]: {field: lane.get(field) for field in lane_fields}
            for lane in sorted(manifest["lanes"], key=lambda item: item["lane_id"])
        },
        "abcd_population": manifest.get("abcd_population"),
        "external_scale": manifest.get("external_scale"),
        "random_seed": manifest.get("random_seed"),
        "toy_count": manifest.get("toy_count"),
    }


def _command_lane_index(args: argparse.Namespace) -> None:
    contract = _read_contract(Path(args.contract).resolve())
    output = Path(args.output).resolve()
    _prepare_output(output)
    direct = [Path(item).resolve() for item in args.lane_json]
    links = [
        {"path": str(path), "sha256": _file_sha256(path)} for path in direct
    ]
    lanes, normalized_links, _ = _load_and_validate_lanes(
        links, contract, require_merge_input=False
    )
    assert len(lanes) == len(_expected_lanes(contract))
    _write_json(
        output,
        {
            "schema": "ppg12-stitched-purity-lane-index/v1",
            "contract_sha256": _payload_sha256(contract),
            "lanes": [
                {
                    "lane_id": row["lane_id"],
                    "path": row["path"],
                    "sha256": row["sha256"],
                }
                for row in normalized_links
            ],
            "lane_set_sha256": _payload_sha256(normalized_links),
        },
    )


def _command_coverage_evidence(args: argparse.Namespace) -> None:
    contract = _read_contract(Path(args.contract).resolve())
    output = Path(args.output).resolve()
    _prepare_output(output)
    index_path = Path(args.lane_index).resolve() if args.lane_index else None
    links = _load_lane_links(
        index_path=index_path,
        lane_paths=[Path(item).resolve() for item in args.lane_json],
        contract=contract,
    )
    lanes, lane_links, bin_edges = _load_and_validate_lanes(
        links, contract, require_merge_input=False
    )
    _write_json(
        output,
        _build_global_minimality_evidence(
            lanes, lane_links, bin_edges, contract
        ),
    )


def _command_manifest(args: argparse.Namespace) -> None:
    contract_path = Path(args.contract).resolve()
    contract = _read_contract(contract_path)
    output = Path(args.output).resolve()
    _prepare_output(output)
    role = args.role
    index_path = Path(args.lane_index).resolve() if args.lane_index else None
    links = _load_lane_links(
        index_path=index_path,
        lane_paths=[Path(item).resolve() for item in args.lane_json],
        contract=contract,
    )
    lanes, lane_links, bin_edges = _load_and_validate_lanes(
        links, contract, require_merge_input=role in {"candidate", "production"}
    )
    _validate_stitched_coverage(lanes, bin_edges)
    minimality, minimality_link = _validate_global_minimality_evidence(
        Path(args.global_minimality_json).resolve(),
        lanes,
        lane_links,
        bin_edges,
        contract,
    )
    provenance, provenance_link = _validate_provenance(
        Path(args.provenance_json).resolve(), contract
    )
    purity, purity_link = _validate_purity(
        Path(args.purity_json).resolve(), contract, bin_edges, lane_links, lanes
    )
    manifest: dict[str, Any] = {
        "schema": "ppg12-stitched-purity-manifest/v1",
        "role": role,
        "abcd_population": contract["abcd_population"],
        "external_scale": float(contract["external_scale"]),
        "random_seed": contract["random_seed"],
        "toy_count": contract["toy_count"],
        "provenance": provenance,
        "lanes": lanes,
        "global_minimality": minimality,
        "purity": purity,
    }
    evidence: dict[str, Any] = {
        "assembler": {
            "path": str(Path(__file__).resolve()),
            "sha256": _file_sha256(Path(__file__).resolve()),
        },
        "contract": {"path": str(contract_path), "sha256": _payload_sha256(contract)},
        "lanes": lane_links,
        "global_minimality": minimality_link,
        "provenance": provenance_link,
        "purity": purity_link,
    }
    if role in {"candidate", "production"}:
        if not args.candidate_parity_json:
            raise AssemblyError(f"--candidate-parity-json is required for role {role}")
        parity, parity_link = _validate_candidate_parity(
            Path(args.candidate_parity_json).resolve(), contract
        )
        _bind_lane_candidate_traces(lanes, parity)
        manifest["candidate_parity"] = parity
        evidence["candidate_parity"] = parity_link
    elif args.candidate_parity_json:
        raise AssemblyError("reference manifests must not include candidate-parity evidence")
    evidence["input_set_sha256"] = _payload_sha256(evidence)
    manifest["assembly"] = evidence
    _write_json(output, manifest)


def _validate_family_summary(
    path: Path,
    expected_family: str,
    expected_input_links: list[dict[str, str]],
    contract: dict[str, Any],
) -> tuple[dict[str, Any], dict[str, str]]:
    payload = _read_json(path, f"{expected_family} merge summary")
    _require_schema(
        payload,
        "ppg12-stitched-purity-family-merge-summary/v2",
        f"{expected_family} merge summary",
    )
    unhashed = dict(payload)
    declared_payload_sha = unhashed.pop("verification_payload_sha256", None)
    if declared_payload_sha != _payload_sha256(unhashed):
        raise AssemblyError(f"{expected_family} merge verifier payload hash is stale")
    if payload.get("family") != expected_family:
        raise AssemblyError(
            f"merge summary family mismatch: expected {expected_family}, observed {payload.get('family')}"
        )
    status = str(payload.get("status", "")).upper()
    if status not in {"PASS", "FAIL"}:
        raise AssemblyError(f"{expected_family} merge status must be PASS or FAIL")
    failures = payload.get("failures")
    if not isinstance(failures, list):
        raise AssemblyError(f"{expected_family} merge failures must be a list")
    content = payload.get("max_content_delta")
    sumw2 = payload.get("max_sumw2_delta")
    if not _is_number(content) or not _is_number(sumw2):
        raise AssemblyError(f"{expected_family} merge deltas must be finite")
    auditor = _canonical_tool_link(
        contract,
        "merge_auditor",
        payload.get("auditor"),
        path.parent,
        f"{expected_family} merge auditor",
    )
    raw_inputs = payload.get("input_links")
    if not isinstance(raw_inputs, list):
        raise AssemblyError(f"{expected_family} merge summary input_links are required")
    inputs = [
        _validate_file_link(
            link, path.parent, f"{expected_family} merge input_links[{index}]"
        )
        for index, link in enumerate(raw_inputs)
    ]
    expected_files = [
        {"path": row["path"], "sha256": row["sha256"]}
        for row in expected_input_links
    ]
    if inputs != expected_files:
        raise AssemblyError(f"{expected_family} merge verifier used different lane inputs")
    if payload.get("input_set_sha256") != _payload_sha256(inputs):
        raise AssemblyError(f"{expected_family} merge verifier input-set hash is stale")
    output_artifact = _validate_file_link(
        payload.get("output_artifact"), path.parent, f"{expected_family} merge output"
    )
    validation = payload.get("root_validation")
    if not isinstance(validation, list):
        raise AssemblyError(f"{expected_family} merge ROOT validation is required")
    expected_validation = [output_artifact, *inputs]
    if len(validation) != len(expected_validation):
        raise AssemblyError(f"{expected_family} merge ROOT validation coverage is incomplete")
    for index, (row, expected) in enumerate(zip(validation, expected_validation)):
        if not isinstance(row, dict):
            raise AssemblyError(f"{expected_family} root_validation[{index}] is malformed")
        linked = _validate_file_link(
            row, path.parent, f"{expected_family} root_validation[{index}]"
        )
        if linked != expected or row.get("zombie") is not False or row.get("recovered") is not False:
            raise AssemblyError(f"{expected_family} ROOT validation is stale or failed")
    required_tokens = payload.get("required_object_family_tokens")
    if not isinstance(required_tokens, dict):
        raise AssemblyError(f"{expected_family} merge verifier token evidence is missing")
    with tempfile.TemporaryDirectory(prefix="ppg12-merge-reverify-") as temp_dir:
        regenerated_path = Path(temp_dir) / "merge_summary.json"
        command = [
            str(_root_capable_python()),
            auditor["path"],
            "--family",
            expected_family,
            "--output",
            output_artifact["path"],
        ]
        for link in inputs:
            command.extend(["--input", link["path"]])
        command.extend(["--json", str(regenerated_path)])
        for token in required_tokens:
            command.extend(["--required-token", token])
        completed = subprocess.run(command, text=True, capture_output=True, check=False)
        if completed.returncode != 0 or not regenerated_path.exists():
            raise AssemblyError(
                f"{expected_family} ROOT-native merge verifier failed: "
                f"{completed.stderr.strip() or completed.stdout.strip()}"
            )
        regenerated = _read_json(regenerated_path, f"regenerated {expected_family} merge summary")
    if regenerated != payload:
        raise AssemblyError(f"{expected_family} merge summary is not canonical verifier output")
    return (
        {
            "family": expected_family,
            "status": status,
            "failures": failures,
            "max_content_delta": float(content),
            "max_sumw2_delta": float(sumw2),
            "output_artifact": output_artifact,
            "auditor": auditor,
            "verifier_input_set_sha256": _payload_sha256(inputs),
            "verification_payload_sha256": declared_payload_sha,
        },
        {"path": str(path.resolve()), "sha256": _file_sha256(path)},
    )


def verify_manifest_evidence(
    path: Path,
    contract: dict[str, Any],
    allowed_roles: set[str],
) -> dict[str, Any]:
    """Re-resolve and recompute every source used by an assembled manifest."""
    path = path.resolve()
    payload = _read_json(path, "stitched-purity manifest")
    _require_schema(payload, "ppg12-stitched-purity-manifest/v1", "manifest")
    role = payload.get("role")
    if role not in allowed_roles:
        raise AssemblyError(
            f"manifest role {role!r} is not one of {sorted(allowed_roles)}"
        )
    assembly = payload.get("assembly")
    if not isinstance(assembly, dict):
        raise AssemblyError("manifest lacks canonical assembly evidence")
    required_keys = {
        "assembler",
        "contract",
        "lanes",
        "global_minimality",
        "provenance",
        "purity",
        "input_set_sha256",
    }
    if role in {"candidate", "production"}:
        required_keys.add("candidate_parity")
    if set(assembly) != required_keys:
        raise AssemblyError(
            f"manifest assembly evidence coverage is not exact: {sorted(assembly)}"
        )
    _canonical_tool_link(
        contract,
        "assembler",
        assembly.get("assembler"),
        path.parent,
        "manifest.assembly.assembler",
    )
    contract_link = assembly.get("contract")
    if not isinstance(contract_link, dict) or not isinstance(contract_link.get("path"), str):
        raise AssemblyError("manifest assembly contract link is missing")
    contract_path = _resolve_path(contract_link["path"], path.parent)
    if _read_contract(contract_path) != contract or contract_link.get("sha256") != _payload_sha256(contract):
        raise AssemblyError("manifest assembly contract link drifted")

    raw_lane_links = assembly.get("lanes")
    if not isinstance(raw_lane_links, list):
        raise AssemblyError("manifest assembly lane links are missing")
    lanes, lane_links, bin_edges = _load_and_validate_lanes(
        raw_lane_links,
        contract,
        require_merge_input=role in {"candidate", "production"},
    )
    _validate_stitched_coverage(lanes, bin_edges)
    if payload.get("lanes") != lanes:
        raise AssemblyError("manifest lanes differ from their canonical extracts")

    minimality_link = _validate_file_link(
        assembly.get("global_minimality"),
        path.parent,
        "manifest.assembly.global_minimality",
    )
    minimality, normalized_minimality_link = _validate_global_minimality_evidence(
        Path(minimality_link["path"]), lanes, lane_links, bin_edges, contract
    )
    if (
        minimality_link != normalized_minimality_link
        or payload.get("global_minimality") != minimality
    ):
        raise AssemblyError(
            "manifest global minimality differs from canonical prefix evidence"
        )

    provenance_link = _validate_file_link(
        assembly.get("provenance"), path.parent, "manifest.assembly.provenance"
    )
    provenance, normalized_provenance_link = _validate_provenance(
        Path(provenance_link["path"]), contract
    )
    if provenance_link != normalized_provenance_link or payload.get("provenance") != provenance:
        raise AssemblyError("manifest provenance differs from its source sets")

    purity_link = _validate_file_link(
        assembly.get("purity"), path.parent, "manifest.assembly.purity"
    )
    purity, normalized_purity_link = _validate_purity(
        Path(purity_link["path"]), contract, bin_edges, lane_links, lanes
    )
    if purity_link != normalized_purity_link or payload.get("purity") != purity:
        raise AssemblyError("manifest purity differs from canonical producer evidence")

    if role in {"candidate", "production"}:
        parity_link = _validate_file_link(
            assembly.get("candidate_parity"),
            path.parent,
            "manifest.assembly.candidate_parity",
        )
        parity, normalized_parity_link = _validate_candidate_parity(
            Path(parity_link["path"]), contract
        )
        _bind_lane_candidate_traces(lanes, parity)
        if parity_link != normalized_parity_link or payload.get("candidate_parity") != parity:
            raise AssemblyError("manifest candidate parity differs from exact raw traces")
    elif "candidate_parity" in payload:
        raise AssemblyError("reference manifest carries candidate-parity summary")

    unhashed = dict(assembly)
    declared_input_set_sha = unhashed.pop("input_set_sha256", None)
    if declared_input_set_sha != _payload_sha256(unhashed):
        raise AssemblyError("manifest assembly input-set hash is stale")
    return payload


def _load_manifest_for_assembly(
    path: Path, contract: dict[str, Any], allowed_roles: set[str]
) -> dict[str, Any]:
    payload = verify_manifest_evidence(path, contract, allowed_roles)
    if payload.get("abcd_population") != contract["abcd_population"]:
        raise AssemblyError("manifest does not use unsuffixed A/B/C/D")
    if payload.get("external_scale") != 1.0:
        raise AssemblyError("manifest uses a forbidden external scale")
    lanes = payload.get("lanes")
    if not isinstance(lanes, list) or len(lanes) != len(_expected_lanes(contract)):
        raise AssemblyError("manifest does not contain the exact 32-lane set")
    indexed = {lane.get("lane_id"): lane for lane in lanes if isinstance(lane, dict)}
    if set(indexed) != set(_expected_lanes(contract)) or len(indexed) != len(lanes):
        raise AssemblyError("manifest lane identities are incomplete or duplicated")
    return payload


def _command_merge_audit(args: argparse.Namespace) -> None:
    contract = _read_contract(Path(args.contract).resolve())
    output = Path(args.output).resolve()
    _prepare_output(output)
    manifest_path = Path(args.candidate_manifest).resolve()
    manifest = _load_manifest_for_assembly(
        manifest_path, contract, {"candidate", "production"}
    )
    summary_paths = {
        "inclusive": Path(args.inclusive_summary).resolve(),
        "photon": Path(args.photon_summary).resolve(),
    }
    audits: list[dict[str, Any]] = []
    all_input_ids: set[str] = set()
    all_input_paths: set[str] = set()
    all_input_hashes: set[str] = set()
    for family in ("inclusive", "photon"):
        rows = sorted(
            (lane for lane in manifest["lanes"] if lane["family"] == family),
            key=lambda lane: lane["lane_id"],
        )
        input_links: list[dict[str, str]] = []
        for lane in rows:
            merge_input = lane.get("merge_input")
            link = _validate_merge_input(
                merge_input, manifest_path, lane["lane_id"]
            )
            if link["input_id"] in all_input_ids:
                raise AssemblyError(f"duplicate merge input_id: {link['input_id']}")
            if link["path"] in all_input_paths:
                raise AssemblyError(
                    f"duplicate physical merge ROOT path: {link['path']}"
                )
            if link["sha256"] in all_input_hashes:
                raise AssemblyError(
                    f"duplicate physical merge ROOT hash: {link['sha256']}"
                )
            all_input_ids.add(link["input_id"])
            all_input_paths.add(link["path"])
            all_input_hashes.add(link["sha256"])
            input_links.append({"lane_id": lane["lane_id"], **link})
        summary, summary_link = _validate_family_summary(
            summary_paths[family], family, input_links, contract
        )
        audits.append(
            {
                **summary,
                "inputs_fixed_order": [row["input_id"] for row in input_links],
                "inputs_fixed_order_count": len(input_links),
                "input_links": input_links,
                "input_set_sha256": _payload_sha256(input_links),
                "family_summary": summary_link,
            }
        )
    if len(all_input_ids) != len(_expected_lanes(contract)):
        raise AssemblyError("merge audit does not bind exactly 32 unique inputs")
    _write_json(
        output,
        {
            "schema": "ppg12-stitched-purity-merge-audit/v1",
            "candidate_manifest": {
                "path": str(manifest_path),
                "sha256": _file_sha256(manifest_path),
            },
            "audits": audits,
            "input_set_sha256": _payload_sha256(
                [row for audit in audits for row in audit["input_links"]]
            ),
        },
    )


def verify_merge_audit_evidence(
    path: Path,
    candidate_manifest_path: Path,
    contract: dict[str, Any],
) -> tuple[dict[str, Any], list[dict[str, str]]]:
    """Re-resolve the ROOT-native family reports and exact merge outputs."""
    path = path.resolve()
    candidate_manifest_path = candidate_manifest_path.resolve()
    payload = _read_json(path, "merge audit")
    _require_schema(payload, "ppg12-stitched-purity-merge-audit/v1", "merge audit")
    candidate_link = _validate_file_link(
        payload.get("candidate_manifest"), path.parent, "merge audit candidate manifest"
    )
    if Path(candidate_link["path"]) != candidate_manifest_path:
        raise AssemblyError("merge audit is bound to a different candidate manifest")
    manifest = _load_manifest_for_assembly(
        candidate_manifest_path, contract, {"candidate", "production"}
    )
    audits = payload.get("audits")
    if not isinstance(audits, list) or len(audits) != 2:
        raise AssemblyError("merge audit must contain exactly two family audits")
    indexed: dict[str, dict[str, Any]] = {}
    all_links: list[dict[str, str]] = []
    all_paths: set[str] = set()
    all_hashes: set[str] = set()
    outputs: list[dict[str, str]] = []
    for row in audits:
        if not isinstance(row, dict) or row.get("family") not in {"inclusive", "photon"}:
            raise AssemblyError("merge audit contains an invalid family")
        family = row["family"]
        if family in indexed:
            raise AssemblyError(f"merge audit duplicates family {family}")
        indexed[family] = row
        expected_lanes = sorted(
            (lane for lane in manifest["lanes"] if lane["family"] == family),
            key=lambda lane: lane["lane_id"],
        )
        expected_links = [
            {
                "lane_id": lane["lane_id"],
                **_validate_merge_input(
                    lane.get("merge_input"), candidate_manifest_path, lane["lane_id"]
                ),
            }
            for lane in expected_lanes
        ]
        for link in expected_links:
            if link["path"] in all_paths:
                raise AssemblyError(
                    f"merge audit reuses physical ROOT path {link['path']}"
                )
            if link["sha256"] in all_hashes:
                raise AssemblyError(
                    f"merge audit reuses physical ROOT hash {link['sha256']}"
                )
            all_paths.add(link["path"])
            all_hashes.add(link["sha256"])
        if row.get("input_links") != expected_links:
            raise AssemblyError(f"merge audit {family} input links drifted")
        if row.get("inputs_fixed_order") != [link["input_id"] for link in expected_links]:
            raise AssemblyError(f"merge audit {family} fixed order drifted")
        if row.get("inputs_fixed_order_count") != len(expected_links):
            raise AssemblyError(f"merge audit {family} fixed input count drifted")
        if row.get("input_set_sha256") != _payload_sha256(expected_links):
            raise AssemblyError(f"merge audit {family} input-set hash drifted")
        summary_link = _validate_file_link(
            row.get("family_summary"), path.parent, f"merge audit {family} family summary"
        )
        summary, normalized_summary_link = _validate_family_summary(
            Path(summary_link["path"]), family, expected_links, contract
        )
        if summary_link != normalized_summary_link:
            raise AssemblyError(f"merge audit {family} family-summary link drifted")
        for field, value in summary.items():
            if row.get(field) != value:
                raise AssemblyError(f"merge audit {family}.{field} differs from ROOT verifier")
        all_links.extend(expected_links)
        outputs.append({"family": family, **summary["output_artifact"]})
    if set(indexed) != {"inclusive", "photon"}:
        raise AssemblyError("merge audit family coverage is incomplete")
    if payload.get("input_set_sha256") != _payload_sha256(all_links):
        raise AssemblyError("combined merge-audit input-set hash drifted")
    outputs.sort(key=lambda row: row["family"])
    return payload, outputs


def _validate_historical_comparison(path: Path) -> tuple[dict[str, float], dict[str, str]]:
    payload = _read_json(path, "historical comparison")
    _require_schema(
        payload,
        "ppg12-stitched-purity-historical-comparison/v1",
        "historical comparison",
    )
    required = (
        "chi2_ndf",
        "max_abs_pull",
        "weighted_mean_ratio",
        "weighted_mean_ratio_error",
        "abcd_coherent_trend_max_sigma",
        "leakage_coherent_trend_max_sigma",
    )
    result: dict[str, float] = {}
    for field in required:
        if not _is_number(payload.get(field)):
            raise AssemblyError(f"historical comparison {field} must be finite")
        result[field] = float(payload[field])
    if result["weighted_mean_ratio_error"] <= 0.0:
        raise AssemblyError("historical weighted_mean_ratio_error must be positive")
    for field in (
        "chi2_ndf",
        "max_abs_pull",
        "abcd_coherent_trend_max_sigma",
        "leakage_coherent_trend_max_sigma",
    ):
        if result[field] < 0.0:
            raise AssemblyError(f"historical comparison {field} must be nonnegative")
    if payload.get("final_purity_series") != "corrected":
        raise AssemblyError("historical comparison must use the corrected purity series")
    input_evidence = payload.get("input")
    if not isinstance(input_evidence, dict) or input_evidence.get("mode") != "direct_root_objects":
        raise AssemblyError(
            "historical comparison must be derived directly from its ROOT objects"
        )
    for field, expected in (
        ("final_purity", None),
        ("abcd", {"A", "B", "C", "D"}),
        ("leakage", {"cB", "cC", "cD"}),
    ):
        details = payload.get(field)
        if not isinstance(details, dict):
            raise AssemblyError(f"historical comparison lacks detailed {field} evidence")
        if expected is not None and set(details) != expected:
            raise AssemblyError(f"historical comparison {field} coverage is incomplete")
    links = payload.get("source_links")
    if not isinstance(links, list):
        raise AssemblyError("historical comparison source_links must be a list")
    normalized_links: list[dict[str, str]] = []
    observed_roles: set[str] = set()
    for index, link in enumerate(links):
        if not isinstance(link, dict):
            raise AssemblyError(f"historical source_links[{index}] must be an object")
        role = link.get("role")
        if role not in REQUIRED_HISTORICAL_SOURCE_ROLES or role in observed_roles:
            raise AssemblyError(f"invalid or duplicate historical source role: {role!r}")
        observed_roles.add(role)
        raw_path = link.get("path")
        if not isinstance(raw_path, str) or not raw_path:
            raise AssemblyError(f"historical source {role} lacks a path")
        source_path = _resolve_path(raw_path, path.parent)
        expected_sha = _require_sha256(
            link.get("sha256"), f"historical source {role}.sha256"
        )
        observed_sha = _file_sha256(source_path)
        if observed_sha != expected_sha:
            raise AssemblyError(f"historical source hash mismatch for {role}")
        normalized_links.append(
            {"role": role, "path": str(source_path), "sha256": observed_sha}
        )
    if observed_roles != REQUIRED_HISTORICAL_SOURCE_ROLES:
        missing = sorted(REQUIRED_HISTORICAL_SOURCE_ROLES - observed_roles)
        raise AssemblyError(f"historical comparison lacks source roles: {missing}")
    normalized_links.sort(key=lambda row: row["role"])
    if payload.get("source_set_sha256") != _payload_sha256(normalized_links):
        raise AssemblyError("historical comparison source-set hash is stale")
    return result, {"path": str(path.resolve()), "sha256": _file_sha256(path)}


def _execute_canonical_admission_replay(
    gate_path: Path,
    contract_path: Path,
    reference_path: str,
    candidate_path: str,
    merge_path: str,
    outdir: Path,
) -> tuple[int, str]:
    command = [
        str(_root_capable_python()),
        str(gate_path),
        "--contract",
        str(contract_path),
        "admit",
        "--reference-manifest",
        reference_path,
        "--candidate-manifest",
        candidate_path,
        "--merge-audit",
        merge_path,
        "--outdir",
        str(outdir),
    ]
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    return completed.returncode, completed.stderr.strip() or completed.stdout.strip()


def _replay_admission(
    admission_path: Path,
    admission: dict[str, Any],
    contract_path: Path,
    contract: dict[str, Any],
) -> dict[str, str]:
    """Re-run the canonical pair gate; a PASS-shaped JSON is never authority."""
    reference_link = _validate_file_link(
        admission.get("reference_manifest"),
        admission_path.parent,
        "admission reference_manifest",
    )
    candidate_link = _validate_file_link(
        admission.get("candidate_manifest"),
        admission_path.parent,
        "admission candidate_manifest",
    )
    merge_link = _validate_file_link(
        admission.get("merge_audit"),
        admission_path.parent,
        "admission merge_audit",
    )
    gate_path = (REPO / contract["canonical_tools"]["closure_gate"]).resolve()
    if not gate_path.is_file():
        raise AssemblyError(f"canonical closure gate is missing: {gate_path}")
    with tempfile.TemporaryDirectory(prefix="ppg12-admission-replay-") as temporary:
        outdir = Path(temporary) / "gate"
        returncode, detail = _execute_canonical_admission_replay(
            gate_path,
            contract_path,
            reference_link["path"],
            candidate_link["path"],
            merge_link["path"],
            outdir,
        )
        regenerated_path = outdir / "admission_manifest.json"
        if returncode != 0 or not regenerated_path.is_file():
            raise AssemblyError(f"canonical admission replay failed: {detail}")
        regenerated = _read_json(regenerated_path, "replayed admission manifest")
    if regenerated != admission:
        raise AssemblyError(
            "admission manifest is not the exact output of the canonical pair gate"
        )
    return {"path": str(gate_path), "sha256": _file_sha256(gate_path)}


def _command_production_wrapper(args: argparse.Namespace) -> None:
    contract_path = Path(args.contract).resolve()
    contract = _read_contract(contract_path)
    output = Path(args.output).resolve()
    _prepare_output(output)
    admission_path = Path(args.admission_manifest).resolve()
    admission = _read_json(admission_path, "admission manifest")
    _require_schema(
        admission, "ppg12-stitched-purity-admission/v1", "admission manifest"
    )
    if admission.get("status") != "PASS":
        raise AssemblyError("production wrapper requires a passing admission")
    if admission.get("contract_sha256") != _payload_sha256(contract):
        raise AssemblyError("admission was created from a different closure contract")
    admission_gate = _replay_admission(
        admission_path, admission, contract_path, contract
    )

    reference_path = Path(args.reference_manifest).resolve()
    candidate_path = Path(args.candidate_manifest).resolve()
    reference = _load_manifest_for_assembly(reference_path, contract, {"reference"})
    candidate = _load_manifest_for_assembly(candidate_path, contract, {"production"})
    if admission.get("reference_frozen") != _frozen_snapshot(reference, contract):
        raise AssemblyError("production reference contract drifted after admission")
    if admission.get("candidate_frozen") != _frozen_snapshot(candidate, contract):
        raise AssemblyError("production candidate contract drifted after admission")

    artifact_paths = {
        "inclusive": Path(args.inclusive_artifact).resolve(),
        "photon": Path(args.photon_artifact).resolve(),
    }
    if artifact_paths["inclusive"] == artifact_paths["photon"]:
        raise AssemblyError("inclusive and photon candidate artifacts must be distinct")
    artifacts = [
        {
            "family": family,
            "path": str(path),
            "sha256": _file_sha256(path),
        }
        for family, path in sorted(artifact_paths.items())
    ]
    merge_audit_path = Path(args.merge_audit).resolve()
    _, merge_artifacts = verify_merge_audit_evidence(
        merge_audit_path, candidate_path, contract
    )
    if artifacts != merge_artifacts:
        raise AssemblyError(
            "production candidate artifacts are not the exact ROOT-native merge outputs"
        )
    historical, historical_link = _validate_historical_comparison(
        Path(args.historical_comparison).resolve()
    )
    wrapper = {
        "schema": "ppg12-stitched-purity-production/v1",
        "admission_sha256": _file_sha256(admission_path),
        "reference_manifest": {
            "path": str(reference_path),
            "sha256": _file_sha256(reference_path),
        },
        "candidate_manifest": {
            "path": str(candidate_path),
            "sha256": _file_sha256(candidate_path),
        },
        "merge_audit": {
            "path": str(merge_audit_path),
            "sha256": _file_sha256(merge_audit_path),
        },
        "candidate_artifacts": artifacts,
        "historical_archive_comparison": historical,
        "assembly": {
            "contract_sha256": _payload_sha256(contract),
            "admission_replay_gate": admission_gate,
            "historical_comparison": historical_link,
            "candidate_artifact_set_sha256": _payload_sha256(artifacts),
        },
    }
    _write_json(output, wrapper)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--contract", default=str(DEFAULT_CONTRACT), help="JSON-compatible YAML contract"
    )
    commands = parser.add_subparsers(dest="command", required=True)

    lane_index = commands.add_parser("lane-index")
    lane_index.add_argument("--lane-json", action="append", required=True)
    lane_index.add_argument("--output", required=True)

    coverage = commands.add_parser("coverage-evidence")
    coverage_lanes = coverage.add_mutually_exclusive_group(required=True)
    coverage_lanes.add_argument("--lane-index")
    coverage_lanes.add_argument("--lane-json", action="append")
    coverage.add_argument("--output", required=True)

    manifest = commands.add_parser("manifest")
    manifest.add_argument("--role", choices=("reference", "candidate", "production"), required=True)
    lanes = manifest.add_mutually_exclusive_group(required=True)
    lanes.add_argument("--lane-index")
    lanes.add_argument("--lane-json", action="append")
    manifest.add_argument("--provenance-json", required=True)
    manifest.add_argument("--purity-json", required=True)
    manifest.add_argument("--global-minimality-json", required=True)
    manifest.add_argument("--candidate-parity-json")
    manifest.add_argument("--output", required=True)

    merge = commands.add_parser("merge-audit")
    merge.add_argument("--candidate-manifest", required=True)
    merge.add_argument("--inclusive-summary", required=True)
    merge.add_argument("--photon-summary", required=True)
    merge.add_argument("--output", required=True)

    production = commands.add_parser("production-wrapper")
    production.add_argument("--admission-manifest", required=True)
    production.add_argument("--reference-manifest", required=True)
    production.add_argument("--candidate-manifest", required=True)
    production.add_argument("--merge-audit", required=True)
    production.add_argument("--inclusive-artifact", required=True)
    production.add_argument("--photon-artifact", required=True)
    production.add_argument("--historical-comparison", required=True)
    production.add_argument("--output", required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        if args.command == "lane-index":
            _command_lane_index(args)
        elif args.command in {"coverage-evidence", "manifest"}:
            # argparse returns a single value for an optional action=append
            # option only when it was used; normalize it for the loader.
            if args.lane_json is None:
                args.lane_json = []
            if args.command == "coverage-evidence":
                _command_coverage_evidence(args)
            else:
                _command_manifest(args)
        elif args.command == "merge-audit":
            _command_merge_audit(args)
        else:
            _command_production_wrapper(args)
    except (AssemblyError, OSError, ValueError) as exc:
        print(json.dumps({"status": "FAIL", "error": str(exc)}, sort_keys=True))
        return 2
    output = Path(args.output).resolve()
    print(
        json.dumps(
            {"status": "PASS", "output": str(output), "sha256": _file_sha256(output)},
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
