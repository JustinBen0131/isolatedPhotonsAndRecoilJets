#!/usr/bin/env python3
"""Materialize one paired photon oracle as two closure-lane inputs.

The paired oracle writes the information needed by
``extract_ppg12_stitched_purity_lane.py`` in several deliberately independent
artifacts.  This producer joins those artifacts without weakening their
contracts:

* exact fill counts come only from the candidate-level executable trace;
* weighted contents and Sumw2 reconstructed from that trace must reproduce the
  executable aggregate cell by cell;
* the reference side is bound to the uninstrumented preserved-PPG12 ROOT;
* the candidate side is bound to the RecoilJets ROOT from the same event set;
* every file is re-hashed and both generated sidecars are run through the
  canonical lane extractor before the output directory is retained.

This tool is intentionally photon-only.  A photon paired-oracle trace contains
the four signal-leakage cells, but it contains no evidence for the twelve
inclusive-lane observables.  In particular, this code never derives raw fill
counts from ``TH1::GetEntries``, ``sumw``, or ``sumw2``.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import re
import shutil
import sys
import tempfile
from pathlib import Path
from types import ModuleType
from typing import Any, Callable


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
DEFAULT_CLOSURE_CONTRACT = (
    REPO / "agent_context/analysis_contracts/ppg12_stitched_purity_closure.yaml"
)
LANE_EXTRACTOR_PATH = HERE / "extract_ppg12_stitched_purity_lane.py"
AGGREGATE_HELPER_PATH = HERE / "extract_ppg12_recoeff_executable_aggregate.py"
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
SIDES = ("reference", "candidate")
REGIONS = ("A", "B", "C", "D")
PRESERVED_EXECUTABLE_EVIDENCE = "preserved_ppg12_executable"
REQUIRED_CERTIFIED_SCOPES = {
    "tags",
    "isolation_abcd",
    "truth_abcd_fills",
    "weights",
}


class MaterializationError(RuntimeError):
    """The paired artifacts cannot honestly become a closure-lane input."""


def _load_module(path: Path, name: str) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise MaterializationError(f"cannot load canonical helper: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


LANE_EXTRACTOR = _load_module(
    LANE_EXTRACTOR_PATH, "ppg12_stitched_purity_lane_extractor_for_materializer"
)
AGGREGATE_HELPER = _load_module(
    AGGREGATE_HELPER_PATH, "ppg12_recoeff_aggregate_for_materializer"
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
    except FileNotFoundError as exc:
        raise MaterializationError(f"missing evidence file: {path}") from exc
    return digest.hexdigest()


def _payload_sha256(payload: Any) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode()
    return hashlib.sha256(encoded).hexdigest()


def _read_json(path: Path, label: str) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text())
    except FileNotFoundError as exc:
        raise MaterializationError(f"missing {label}: {path}") from exc
    except json.JSONDecodeError as exc:
        raise MaterializationError(f"invalid JSON in {label} {path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise MaterializationError(f"{label} must contain one JSON object: {path}")
    return payload


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    temporary.replace(path)


def _file(path: str | Path, label: str) -> Path:
    resolved = Path(path).expanduser().resolve()
    if not resolved.is_file() or resolved.stat().st_size <= 0:
        raise MaterializationError(f"{label} is missing or empty: {resolved}")
    return resolved


def _link(path: Path) -> dict[str, str]:
    path = _file(path, "linked evidence")
    return {"path": str(path), "sha256": _sha256(path)}


def _require_link(
    value: Any,
    *,
    label: str,
    expected_path: Path | None = None,
) -> dict[str, str]:
    if not isinstance(value, dict):
        raise MaterializationError(f"{label} must contain path and sha256")
    raw_path = value.get("path")
    digest = value.get("sha256")
    if not isinstance(raw_path, str) or not raw_path:
        raise MaterializationError(f"{label}.path is required")
    if not isinstance(digest, str) or SHA256_RE.fullmatch(digest) is None:
        raise MaterializationError(f"{label}.sha256 is malformed")
    path = _file(raw_path, label)
    if expected_path is not None and path != expected_path.resolve():
        raise MaterializationError(
            f"{label} path differs from the paired runtime contract"
        )
    observed = _sha256(path)
    if observed != digest:
        raise MaterializationError(
            f"{label} hash mismatch: expected {digest}, observed {observed}"
        )
    return {"path": str(path), "sha256": observed}


def _canonical_identity(contract: dict[str, Any]) -> dict[str, str]:
    if contract.get("schema_version") != 3:
        raise MaterializationError("paired runtime contract schema_version must be 3")
    lane = contract.get("lane")
    if not isinstance(lane, dict) or lane.get("rows") != 5:
        raise MaterializationError("paired runtime contract must describe five rows")
    lane_id = lane.get("lane_id")
    if not isinstance(lane_id, str):
        raise MaterializationError("paired runtime contract lacks lane_id")
    parts = lane_id.split(":")
    if len(parts) != 4 or parts[0] != "photon":
        raise MaterializationError(
            "paired-oracle lane materialization is photon-only"
        )
    family, sample, period, interaction = parts
    if lane.get("sample", "").lower() != sample:
        raise MaterializationError("paired runtime sample does not match lane_id")
    if lane.get("period") != period or str(lane.get("interaction", "")).lower() != interaction:
        raise MaterializationError("paired runtime period/interaction does not match lane_id")
    return {
        "lane_id": lane_id,
        "family": family,
        "sample": sample,
        "period": period,
        "interaction": interaction,
    }


def _contract_paths(contract: dict[str, Any]) -> dict[str, Path]:
    raw = contract.get("paths")
    if not isinstance(raw, dict):
        raise MaterializationError("paired runtime contract lacks paths")
    required = {
        "g4_full_list",
        "g4_slice",
        "recoil_runtime_manifest",
        "paired_runtime_manifest",
        "recoil_config",
        "ppg_recoeff_period_config",
        "ppg_recoeff_truth_vertex_reweight",
        "ppg_recoeff_yaml_cpp_header_tree_receipt",
        "recoil_root",
        "candidate_csv",
        "ppg_recoeff_baseline_eff_root",
        "ppg12_executable_trace",
        "ppg12_executable_response_trace",
        "ppg12_executable_aggregate",
    }
    missing = sorted(required - set(raw))
    if missing:
        raise MaterializationError(f"paired runtime contract lacks paths: {missing}")
    return {name: _file(raw[name], f"paired path {name}") for name in required}


def _validate_completed_run(contract_path: Path) -> None:
    state = contract_path.parent / "RUN_STATE"
    if not state.is_file() or state.read_text().strip() != "PASS":
        raise MaterializationError(
            f"paired runtime contract is not owned by a completed PASS run: {state}"
        )


def _validate_aggregate(
    aggregate_path: Path,
    aggregate: dict[str, Any],
    *,
    contract_path: Path,
    paths: dict[str, Path],
    lane_id: str,
) -> None:
    if aggregate.get("schema_version") != 1:
        raise MaterializationError("executable aggregate schema_version must be 1")
    if aggregate.get("evidence_source") != "preserved_ppg12_executable_aggregate":
        raise MaterializationError("aggregate is not preserved-executable evidence")
    if aggregate.get("mode") != "full" or aggregate.get("status") != "PASS":
        raise MaterializationError("full executable aggregate did not pass")
    equivalence = aggregate.get("root_equivalence")
    if (
        not isinstance(equivalence, dict)
        or equivalence.get("pass") is not True
        or equivalence.get("efficiency_root_exact") is not True
        or equivalence.get("response_root_exact") is not True
    ):
        raise MaterializationError(
            "instrumented and uninstrumented PPG12 ROOT outputs are not identical"
        )
    comparison = aggregate.get("aggregate_comparison")
    if not isinstance(comparison, dict) or comparison.get("pass") is not True:
        raise MaterializationError(
            "reference and candidate weighted leakage cells do not close"
        )
    if set(aggregate.get("certified_scopes", [])) != REQUIRED_CERTIFIED_SCOPES:
        raise MaterializationError("aggregate does not certify the full fill contract")
    identity = aggregate.get("lane_identity")
    if not isinstance(identity, dict) or identity.get("lane_id") != lane_id:
        raise MaterializationError("aggregate belongs to another lane")
    contract_sha = _sha256(contract_path)
    if identity.get("runtime_contract_sha256") != contract_sha:
        raise MaterializationError("aggregate is bound to a stale runtime contract")
    provenance = aggregate.get("provenance")
    if not isinstance(provenance, dict):
        raise MaterializationError("aggregate lacks provenance")
    expected_links = {
        "runtime_contract": contract_path,
        "candidate_csv": paths["candidate_csv"],
        "runtime_manifest": paths["paired_runtime_manifest"],
        "baseline_root": paths["ppg_recoeff_baseline_eff_root"],
        "trace_csv": paths["ppg12_executable_trace"],
        "response_trace_csv": paths["ppg12_executable_response_trace"],
    }
    for role, expected_path in expected_links.items():
        _require_link(
            provenance.get(role),
            label=f"aggregate provenance {role}",
            expected_path=expected_path,
        )
    if aggregate_path != paths["ppg12_executable_aggregate"]:
        raise MaterializationError("aggregate path differs from paired runtime contract")


def _validate_selected_runtime_manifest(
    manifest_path: Path,
    *,
    period: str,
    ppg_period_config: Path,
    ppg_truth_vertex_reweight: Path,
    ppg_yaml_cpp_header_receipt: Path,
) -> dict[str, dict[str, str]]:
    manifest = _read_json(manifest_path, "paired runtime manifest")
    if manifest.get("selected_period") != period:
        raise MaterializationError("paired runtime selected_period differs from lane")
    files = manifest.get("files")
    if not isinstance(files, list):
        raise MaterializationError("paired runtime manifest lacks files")
    selected_roles = {
        "ppg_recoeff_period_config": ppg_period_config,
        "ppg_recoeff_truth_vertex_reweight": ppg_truth_vertex_reweight,
        "ppg_recoeff_yaml_cpp_header_tree_receipt": ppg_yaml_cpp_header_receipt,
    }
    forbidden = {
        "ppg_recoeff_period_config_0mrad",
        "ppg_recoeff_period_config_1p5mrad",
        "ppg_recoeff_truth_vertex_reweight_0mrad",
        "ppg_recoeff_truth_vertex_reweight_1p5mrad",
    }
    retained = sorted(
        str(entry.get("role"))
        for entry in files
        if isinstance(entry, dict) and entry.get("role") in forbidden
    )
    if retained:
        raise MaterializationError(
            "paired runtime manifest retains unselected period roles: "
            f"{retained}"
        )
    links: dict[str, dict[str, str]] = {}
    for role, expected_path in selected_roles.items():
        selected = [
            entry
            for entry in files
            if isinstance(entry, dict) and entry.get("role") == role
        ]
        if len(selected) != 1:
            raise MaterializationError(
                f"paired runtime manifest must contain exactly one {role} role"
            )
        links[role] = _require_link(
            selected[0],
            label=f"paired runtime {role}",
            expected_path=expected_path,
        )
    return links


def _read_rows(
    path: Path, *, lane_id: str, runtime_contract_sha256: str
) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        rows = list(reader)
        fields = set(reader.fieldnames or [])
    if not rows:
        raise MaterializationError("paired candidate CSV contains no rows")
    required = {
        "lane_id",
        "runtime_contract_sha256",
        "match_status",
        "candidate_identity",
        "rj_cluster_Et",
        "rj_weight_final",
        "rj_signal_fill_multiplicity",
        "ppg12_response_Et",
        "ppg12_weight_final",
        "ppg12_fill_multiplicity",
        "ppg12_tag_evidence_source",
        "ppg12_isolation_abcd_evidence_source",
        "ppg12_truth_response_fill_evidence_source",
        "ppg12_weight_evidence_source",
    }
    required.update(f"rj_signal_fill_{region}" for region in REGIONS)
    required.update(f"ppg12_signal_fill_{region}" for region in REGIONS)
    missing = sorted(required - fields)
    if missing:
        raise MaterializationError(f"paired candidate CSV lacks fields: {missing}")
    identities: set[str] = set()
    for row_number, row in enumerate(rows, start=2):
        if row.get("lane_id") != lane_id:
            raise MaterializationError(f"candidate CSV row {row_number} belongs to another lane")
        if row.get("runtime_contract_sha256") != runtime_contract_sha256:
            raise MaterializationError(
                f"candidate CSV row {row_number} is bound to a stale runtime contract"
            )
        identity = row.get("candidate_identity", "")
        if not identity or identity in identities:
            raise MaterializationError(
                f"missing or duplicate candidate identity at CSV row {row_number}"
            )
        identities.add(identity)
        if row.get("match_status") not in {"matched", "rj_only", "ppg12_only"}:
            raise MaterializationError(f"invalid match_status at CSV row {row_number}")
    return rows


def _integer(value: str, *, label: str) -> int:
    try:
        parsed = int(value)
    except (TypeError, ValueError) as exc:
        raise MaterializationError(f"{label} is not an integer") from exc
    if parsed not in (0, 1):
        raise MaterializationError(f"{label} is not binary")
    return parsed


def _finite_nonnegative(value: str, *, label: str) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise MaterializationError(f"{label} is not numeric") from exc
    if not math.isfinite(parsed) or parsed < 0.0:
        raise MaterializationError(f"{label} must be finite and nonnegative")
    return parsed


def _empty_cells() -> dict[str, list[dict[str, float]]]:
    return {
        region: [
            {"content": 0.0, "sumw2": 0.0, "fills": 0.0}
            for _ in range(len(AGGREGATE_HELPER.RECO_EDGES) + 1)
        ]
        for region in REGIONS
    }


def _aggregate_rows(
    rows: list[dict[str, str]], *, side: str
) -> dict[str, list[dict[str, float]]]:
    if side == "reference":
        prefix = "ppg12"
        multiplicity_field = "ppg12_fill_multiplicity"
        et_field = "ppg12_response_Et"
        weight_field = "ppg12_weight_final"
        present_statuses = {"matched", "ppg12_only"}
    elif side == "candidate":
        prefix = "rj"
        multiplicity_field = "rj_signal_fill_multiplicity"
        et_field = "rj_cluster_Et"
        weight_field = "rj_weight_final"
        present_statuses = {"matched", "rj_only"}
    else:
        raise MaterializationError(f"unsupported paired side: {side}")

    cells = _empty_cells()
    for row_number, row in enumerate(rows, start=2):
        present = row["match_status"] in present_statuses
        flag_fields = [f"{prefix}_signal_fill_{region}" for region in REGIONS]
        raw_values = [row.get(field, "") for field in flag_fields]
        raw_multiplicity = row.get(multiplicity_field, "")
        if not present:
            if any(value != "" for value in [*raw_values, raw_multiplicity]):
                raise MaterializationError(
                    f"absent {side} candidate carries fill evidence at CSV row {row_number}"
                )
            continue
        flags = {
            region: _integer(
                row.get(f"{prefix}_signal_fill_{region}", ""),
                label=f"row {row_number} {prefix}_signal_fill_{region}",
            )
            for region in REGIONS
        }
        try:
            multiplicity = int(raw_multiplicity)
        except (TypeError, ValueError) as exc:
            raise MaterializationError(
                f"row {row_number} {multiplicity_field} is not an integer"
            ) from exc
        if multiplicity not in (0, 1) or multiplicity != sum(flags.values()):
            raise MaterializationError(
                f"row {row_number} {side} fill multiplicity is inconsistent"
            )
        if side == "reference":
            for field in (
                "ppg12_tag_evidence_source",
                "ppg12_isolation_abcd_evidence_source",
                "ppg12_truth_response_fill_evidence_source",
                "ppg12_weight_evidence_source",
            ):
                if row.get(field) != PRESERVED_EXECUTABLE_EVIDENCE:
                    raise MaterializationError(
                        f"row {row_number} lacks preserved-executable {field}"
                    )
        if multiplicity == 0:
            continue
        et = _finite_nonnegative(row.get(et_field, ""), label=f"row {row_number} {et_field}")
        weight = _finite_nonnegative(
            row.get(weight_field, ""), label=f"row {row_number} {weight_field}"
        )
        global_bin = AGGREGATE_HELPER.root_bin_index(
            et, list(AGGREGATE_HELPER.RECO_EDGES)
        )
        region = next(name for name, enabled in flags.items() if enabled)
        cell = cells[region][global_bin]
        cell["content"] += weight
        cell["sumw2"] += weight * weight
        cell["fills"] += 1.0
    return cells


def _aggregate_cells(value: Any, *, label: str) -> dict[str, list[dict[str, float]]]:
    if not isinstance(value, dict) or set(value) != set(REGIONS):
        raise MaterializationError(f"{label} has incomplete region coverage")
    expected_cells = len(AGGREGATE_HELPER.RECO_EDGES) + 1
    normalized: dict[str, list[dict[str, float]]] = {}
    for region in REGIONS:
        raw_cells = value[region]
        if not isinstance(raw_cells, list) or len(raw_cells) != expected_cells:
            raise MaterializationError(f"{label}.{region} has wrong cell count")
        normalized[region] = []
        for index, row in enumerate(raw_cells):
            if not isinstance(row, dict) or set(row) != {"content", "sumw2"}:
                raise MaterializationError(
                    f"{label}.{region}[{index}] is malformed"
                )
            normalized[region].append(
                {
                    "content": _finite_nonnegative(
                        row["content"], label=f"{label}.{region}[{index}].content"
                    ),
                    "sumw2": _finite_nonnegative(
                        row["sumw2"], label=f"{label}.{region}[{index}].sumw2"
                    ),
                }
            )
    return normalized


def _close(left: float, right: float) -> bool:
    return abs(left - right) <= 1.0e-12 + 1.0e-6 * abs(right)


def _require_moment_closure(
    observed: dict[str, list[dict[str, float]]],
    expected: dict[str, list[dict[str, float]]],
    *,
    side: str,
) -> None:
    for region in REGIONS:
        for index, (left, right) in enumerate(zip(observed[region], expected[region])):
            for field in ("content", "sumw2"):
                if not _close(left[field], right[field]):
                    raise MaterializationError(
                        f"{side} candidate trace does not reproduce aggregate "
                        f"{region} global bin {index} {field}: "
                        f"trace={left[field]}, aggregate={right[field]}"
                    )


def _fill_payload(
    *,
    side: str,
    lane_id: str,
    root_link: dict[str, str],
    cells: dict[str, list[dict[str, float]]],
    object_prefix: str,
    candidate_link: dict[str, str],
    aggregate_link: dict[str, str],
    runtime_contract_link: dict[str, str],
) -> dict[str, Any]:
    observables: dict[str, Any] = {}
    edges = [float(value) for value in AGGREGATE_HELPER.RECO_EDGES]
    for region in REGIONS:
        basename = LANE_EXTRACTOR.REGION_OBJECTS[f"{region}_signal"]
        object_path = f"{object_prefix}/{basename}" if object_prefix else basename
        observables[f"{region}_signal"] = {
            "object": object_path,
            "bin_edges": edges,
            "fills": [int(cell["fills"]) for cell in cells[region][1:-1]],
            "flow_fills": [int(cells[region][0]["fills"]), int(cells[region][-1]["fills"])],
        }
    return {
        "schema": "ppg12-stitched-purity-fill-evidence/v1",
        "lane_id": lane_id,
        "root_sha256": root_link["sha256"],
        "side": side,
        "fill_count_semantics": (
            "candidate-level executable fill flags; never inferred from TH1 moments"
        ),
        "candidate_csv": candidate_link,
        "executable_aggregate": aggregate_link,
        "runtime_contract": runtime_contract_link,
        "observables": observables,
    }


def _lane_input(
    *,
    identity: dict[str, str],
    root_link: dict[str, str],
    object_prefix: str,
    source_list_link: dict[str, str],
    event_set_link: dict[str, str],
    config_link: dict[str, str],
    runtime_manifest_link: dict[str, str],
    candidate_link: dict[str, str],
    fill_sha256: str,
) -> dict[str, Any]:
    event_rows = LANE_EXTRACTOR._event_set_rows(Path(event_set_link["path"]))
    groups = LANE_EXTRACTOR._canonical_groups(identity["lane_id"], event_rows, 5)
    if len(groups) != 1:
        raise MaterializationError("paired oracle must contain exactly one five-row group")
    return {
        "schema": "ppg12-stitched-purity-lane-input/v1",
        **identity,
        "group_count": 1,
        "group_size": 5,
        "group_index_start": 0,
        "event_set_row_count": 5,
        "groups": groups,
        "group_set_sha256": _payload_sha256(groups),
        "group_input_sidecars": [],
        "external_scale": 1.0,
        "abcd_population": "unsuffixed",
        "object_prefix": object_prefix,
        "root": root_link,
        "evidence": {
            "source_list": source_list_link,
            "event_set": event_set_link,
            "config": config_link,
        },
        "runtime_manifest": runtime_manifest_link,
        "fill_evidence": {"path": "fill_evidence.json", "sha256": fill_sha256},
        "candidate_parity_evidence": [
            {"role": "candidate_rows", **candidate_link}
        ],
    }


def materialize(
    paired_contract_path: Path,
    output_dir: Path,
    closure_contract_path: Path = DEFAULT_CLOSURE_CONTRACT,
    *,
    reader_factory: Callable[[Path], Any] | None = None,
) -> dict[str, Any]:
    """Create and validate reference/candidate lane inputs atomically."""
    paired_contract_path = paired_contract_path.resolve()
    output_dir = output_dir.resolve()
    closure_contract_path = closure_contract_path.resolve()
    if output_dir.exists():
        raise MaterializationError(f"refusing to overwrite output directory: {output_dir}")
    _file(closure_contract_path, "closure contract")
    paired_contract = _read_json(paired_contract_path, "paired runtime contract")
    identity = _canonical_identity(paired_contract)
    paths = _contract_paths(paired_contract)
    _validate_completed_run(paired_contract_path)
    aggregate_path = paths["ppg12_executable_aggregate"]
    aggregate = _read_json(aggregate_path, "executable aggregate")
    _validate_aggregate(
        aggregate_path,
        aggregate,
        contract_path=paired_contract_path,
        paths=paths,
        lane_id=identity["lane_id"],
    )
    paired_contract_link = _link(paired_contract_path)
    aggregate_link = _link(aggregate_path)
    candidate_link = _link(paths["candidate_csv"])
    selected_runtime_links = _validate_selected_runtime_manifest(
        paths["paired_runtime_manifest"],
        period=identity["period"],
        ppg_period_config=paths["ppg_recoeff_period_config"],
        ppg_truth_vertex_reweight=paths["ppg_recoeff_truth_vertex_reweight"],
        ppg_yaml_cpp_header_receipt=paths[
            "ppg_recoeff_yaml_cpp_header_tree_receipt"
        ],
    )
    runtime_manifest_link = _link(paths["paired_runtime_manifest"])
    rows = _read_rows(
        paths["candidate_csv"],
        lane_id=identity["lane_id"],
        runtime_contract_sha256=paired_contract_link["sha256"],
    )
    cells = {
        "reference": _aggregate_rows(rows, side="reference"),
        "candidate": _aggregate_rows(rows, side="candidate"),
    }
    expected_cells = {
        "reference": _aggregate_cells(aggregate.get("oracle_cells"), label="oracle_cells"),
        "candidate": _aggregate_cells(aggregate.get("recoil_cells"), label="recoil_cells"),
    }
    for side in SIDES:
        _require_moment_closure(cells[side], expected_cells[side], side=side)

    roots = {
        "reference": _link(paths["ppg_recoeff_baseline_eff_root"]),
        "candidate": _link(paths["recoil_root"]),
    }
    prefixes = {"reference": "", "candidate": "SIM"}
    source_list_link = _link(paths["g4_full_list"])
    event_set_link = _link(paths["g4_slice"])
    config_link = _link(paths["recoil_config"])

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(
        tempfile.mkdtemp(prefix=f".{output_dir.name}.staging-", dir=output_dir.parent)
    )
    moved = False
    try:
        for side in SIDES:
            side_dir = staging / side
            fill_payload = _fill_payload(
                side=side,
                lane_id=identity["lane_id"],
                root_link=roots[side],
                cells=cells[side],
                object_prefix=prefixes[side],
                candidate_link=candidate_link,
                aggregate_link=aggregate_link,
                runtime_contract_link=paired_contract_link,
            )
            fill_path = side_dir / "fill_evidence.json"
            _write_json(fill_path, fill_payload)
            lane_input = _lane_input(
                identity=identity,
                root_link=roots[side],
                object_prefix=prefixes[side],
                source_list_link=source_list_link,
                event_set_link=event_set_link,
                config_link=config_link,
                runtime_manifest_link=runtime_manifest_link,
                candidate_link=candidate_link,
                fill_sha256=_sha256(fill_path),
            )
            _write_json(side_dir / "lane_input.json", lane_input)
        staging.rename(output_dir)
        moved = True

        side_links: dict[str, Any] = {}
        for side in SIDES:
            side_dir = output_dir / side
            kwargs: dict[str, Any] = {}
            if reader_factory is not None:
                kwargs["reader_factory"] = reader_factory
            try:
                extracted = LANE_EXTRACTOR.extract_lane(
                    side_dir / "lane_input.json", closure_contract_path, **kwargs
                )
            except Exception as exc:
                raise MaterializationError(
                    f"canonical lane extraction rejected {side} materialization: {exc}"
                ) from exc
            lane_path = side_dir / "lane.json"
            _write_json(lane_path, extracted)
            side_links[side] = {
                "lane_input": _link(side_dir / "lane_input.json"),
                "fill_evidence": _link(side_dir / "fill_evidence.json"),
                "lane": _link(lane_path),
            }

        manifest = {
            "schema": "ppg12-paired-oracle-lane-input-bundle/v1",
            "status": "PASS",
            **identity,
            "paired_runtime_contract": paired_contract_link,
            "executable_aggregate": aggregate_link,
            "candidate_rows": candidate_link,
            "recoil_runtime_manifest": _link(paths["recoil_runtime_manifest"]),
            "paired_runtime_manifest": runtime_manifest_link,
            "ppg_recoeff_period_config": selected_runtime_links[
                "ppg_recoeff_period_config"
            ],
            "ppg_recoeff_truth_vertex_reweight": selected_runtime_links[
                "ppg_recoeff_truth_vertex_reweight"
            ],
            "ppg_recoeff_yaml_cpp_header_tree_receipt": selected_runtime_links[
                "ppg_recoeff_yaml_cpp_header_tree_receipt"
            ],
            "closure_contract": _link(closure_contract_path),
            "sides": side_links,
            "fill_count_semantics": (
                "exact candidate-level executable flags cross-checked against "
                "weighted aggregate cells"
            ),
        }
        manifest_path = output_dir / "materialization_manifest.json"
        _write_json(manifest_path, manifest)
        return {**manifest, "manifest": _link(manifest_path)}
    except Exception:
        if moved:
            shutil.rmtree(output_dir, ignore_errors=True)
        else:
            shutil.rmtree(staging, ignore_errors=True)
        raise


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--paired-contract", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--closure-contract", type=Path, default=DEFAULT_CLOSURE_CONTRACT)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        result = materialize(
            args.paired_contract, args.output_dir, args.closure_contract
        )
    except (MaterializationError, OSError, ValueError) as exc:
        print(json.dumps({"status": "FAIL", "error": str(exc)}, sort_keys=True))
        return 2
    print(
        json.dumps(
            {
                "status": "PASS",
                "manifest": result["manifest"]["path"],
                "sha256": result["manifest"]["sha256"],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
