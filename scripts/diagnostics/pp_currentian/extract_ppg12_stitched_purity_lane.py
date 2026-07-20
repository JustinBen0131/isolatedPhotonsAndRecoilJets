#!/usr/bin/env python3
"""Extract one fail-closed PPG12 stitched-purity lane from a ROOT file.

The 32-lane closure manifest deliberately consumes compact JSON rather than
opening ROOT files itself.  This bridge validates a hash-bound lane sidecar,
reads only the canonical purity histograms, and writes exactly one
``ppg12-stitched-purity-lane/v1`` payload for
``assemble_ppg12_stitched_purity_manifest.py``.

Weighted TH1 objects do not contain per-bin fill counts.  Consequently this
extractor never treats ``GetEntries`` or ``sumw**2/sumw2`` as per-bin fills.
Every lane must provide a hash-bound fill-evidence JSON containing the exact
per-bin and flow-bin fill counts; their total is checked against TH1 entries.

Input sidecar schema (paths may be relative to the sidecar)::

  {
    "schema": "ppg12-stitched-purity-lane-input/v1",
    "lane_id": "inclusive:jet8:0mrad:si",
    "family": "inclusive", "sample": "jet8",
    "period": "0mrad", "interaction": "si",
    "group_count": 1, "group_size": 5, "group_index_start": 0,
    "event_set_row_count": 5,
    "groups": [{
      "group_index": 0,
      "group_id": "group-00000-<hash-prefix>",
      "event_rows": ["file1", "file2", "file3", "file4", "file5"],
      "group_sha256": "..."
    }],
    "group_set_sha256": "...",
    "group_input_sidecars": [],
    "external_scale": 1.0,
    "abcd_population": "unsuffixed", "object_prefix": "SIM",
    "root": {"path": "lane.root", "sha256": "..."},
    "evidence": {
      "source_list": {"path": "...", "sha256": "..."},
      "event_set": {"path": "...", "sha256": "..."},
      "config": {"path": "...", "sha256": "..."}
    },
    "runtime_manifest": {"path": "runtime_manifest.json", "sha256": "..."},
    "fill_evidence": {"path": "fills.json", "sha256": "..."},
    "candidate_parity_evidence": [
      {"role": "candidate_rows", "path": "...", "sha256": "..."}
    ]
  }

The event-set rows are not an independently hashable claim: they must equal
the ordered five-row production-list slices selected by their declared group
indices.  Runtime contract hashes are likewise derived by this extractor from
the actual role-labelled files in the source-locked runtime manifest, rather
than copied from caller-supplied digest strings.

Every lane requires one or more exact deterministic five-file groups.  A
multi-group aggregate additionally links one independently extractable input
sidecar per group.  The assembler re-runs this extractor on those sidecars and
requires their histogram sum to equal the aggregate, so group coverage cannot
be supplied as a hand-authored claim.  Photon lanes additionally require the
exact candidate-row trace; inclusive lanes must not carry photon trace
evidence.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
from pathlib import Path
from typing import Any, Callable, Protocol


REPO = Path(__file__).resolve().parents[3]
DEFAULT_CONTRACT = (
    REPO / "agent_context/analysis_contracts/ppg12_stitched_purity_closure.yaml"
)
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
EVIDENCE_KEYS = ("source_list", "event_set", "config")
DIRECT_HASH_FIELD_BY_EVIDENCE = {
    "source_list": "source_list_sha256",
    "event_set": "event_set_sha256",
}
REGION_OBJECTS = {
    "A": "h_tight_iso_cluster_0",
    "B": "h_tight_noniso_cluster_0",
    "C": "h_nontight_iso_cluster_0",
    "D": "h_nontight_noniso_cluster_0",
    "A_signal": "h_tight_iso_cluster_signal_0",
    "B_signal": "h_tight_noniso_cluster_signal_0",
    "C_signal": "h_nontight_iso_cluster_signal_0",
    "D_signal": "h_nontight_noniso_cluster_signal_0",
    "A_notmatch": "h_tight_iso_cluster_notmatch_0",
    "B_notmatch": "h_tight_noniso_cluster_notmatch_0",
    "C_notmatch": "h_nontight_iso_cluster_notmatch_0",
    "D_notmatch": "h_nontight_noniso_cluster_notmatch_0",
}


class ExtractionError(RuntimeError):
    """The lane cannot be represented without inventing evidence."""


class HistogramReader(Protocol):
    def paths_for_basename(self, basename: str) -> list[str]: ...

    def cycle_count(self, path: str) -> int: ...

    def histogram(self, path: str) -> dict[str, Any]: ...

    def close(self) -> None: ...


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
    except FileNotFoundError as exc:
        raise ExtractionError(f"missing evidence file: {path}") from exc
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
        raise ExtractionError(f"missing {label}: {path}") from exc
    except json.JSONDecodeError as exc:
        raise ExtractionError(f"invalid JSON in {label} {path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise ExtractionError(f"{label} must contain a JSON object: {path}")
    return payload


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    temporary.replace(path)


def _resolve_path(raw: str, parent: Path) -> Path:
    path = Path(raw).expanduser()
    return path.resolve() if path.is_absolute() else (parent / path).resolve()


def _require_sha256(value: Any, label: str) -> str:
    if not isinstance(value, str) or SHA256_RE.fullmatch(value) is None:
        raise ExtractionError(f"{label} must be a lowercase 64-character SHA-256")
    return value


def _validate_link(link: Any, parent: Path, label: str) -> dict[str, str]:
    if not isinstance(link, dict):
        raise ExtractionError(f"{label} must contain path and sha256")
    raw_path = link.get("path")
    if not isinstance(raw_path, str) or not raw_path:
        raise ExtractionError(f"{label}.path is required")
    expected = _require_sha256(link.get("sha256"), f"{label}.sha256")
    path = _resolve_path(raw_path, parent)
    if not path.is_file():
        raise ExtractionError(f"{label} is not a regular file: {path}")
    if path.stat().st_size <= 0:
        raise ExtractionError(f"{label} is empty: {path}")
    observed = _file_sha256(path)
    if observed != expected:
        raise ExtractionError(
            f"{label} hash mismatch for {path}: expected {expected}, observed {observed}"
        )
    return {"path": str(path), "sha256": observed}


def _contract(path: Path) -> dict[str, Any]:
    payload = _read_json(path, "closure contract")
    if payload.get("schema") != "ppg12-stitched-purity-closure-contract/v1":
        raise ExtractionError("unsupported closure contract schema")
    return payload


def _lane_id(family: str, sample: str, period: str, interaction: str) -> str:
    return f"{family}:{sample}:{period}:{interaction}"


def _expected_lanes(contract: dict[str, Any]) -> set[str]:
    return {
        _lane_id(family, sample, period, interaction)
        for family, family_spec in contract["families"].items()
        for sample in family_spec["samples"]
        for period in contract["periods"]
        for interaction in contract["interactions"]
    }


def _finite_vector(value: Any, label: str, *, nonnegative: bool = False) -> list[float]:
    if not isinstance(value, (list, tuple)) or not value:
        raise ExtractionError(f"{label} must be a non-empty vector")
    result: list[float] = []
    for item in value:
        if isinstance(item, bool) or not isinstance(item, (int, float)):
            raise ExtractionError(f"{label} contains a non-numeric value")
        item_float = float(item)
        if not math.isfinite(item_float):
            raise ExtractionError(f"{label} contains a non-finite value")
        if nonnegative and item_float < 0.0:
            raise ExtractionError(f"{label} contains a negative value")
        result.append(item_float)
    return result


def _count_vector(value: Any, label: str, length: int) -> list[float]:
    result = _finite_vector(value, label, nonnegative=True)
    if len(result) != length:
        raise ExtractionError(f"{label} has {len(result)} rather than {length} cells")
    for index, item in enumerate(result):
        if abs(item - round(item)) > 1e-9:
            raise ExtractionError(f"{label}[{index}] is not an exact fill count: {item}")
    return [float(round(item)) for item in result]


def _same_vector(left: list[float], right: list[float], *, tol: float = 1e-12) -> bool:
    return len(left) == len(right) and all(
        abs(a - b) <= tol * max(1.0, abs(a), abs(b)) for a, b in zip(left, right)
    )


class UprootHistogramReader:
    """Small uproot adapter; importing ROOT is intentionally unnecessary."""

    def __init__(self, path: Path):
        try:
            import uproot  # type: ignore
        except ImportError as exc:
            raise ExtractionError(
                "uproot is required; run with the ThesisAnalysis Python runtime"
            ) from exc
        try:
            self._file = uproot.open(str(path))
            cycle_keys = list(self._file.keys(recursive=True, cycle=True))
        except Exception as exc:
            raise ExtractionError(f"unreadable or malformed ROOT file: {path}: {exc}") from exc
        self._paths_by_basename: dict[str, set[str]] = {}
        self._cycles: dict[str, int] = {}
        for raw_key in cycle_keys:
            key = str(raw_key)
            match = re.fullmatch(r"(.+);([0-9]+)", key)
            path_without_cycle = match.group(1) if match else key
            self._cycles[path_without_cycle] = self._cycles.get(path_without_cycle, 0) + 1
            basename = path_without_cycle.rsplit("/", 1)[-1]
            self._paths_by_basename.setdefault(basename, set()).add(path_without_cycle)

    def paths_for_basename(self, basename: str) -> list[str]:
        return sorted(self._paths_by_basename.get(basename, set()))

    def cycle_count(self, path: str) -> int:
        return self._cycles.get(path, 0)

    def histogram(self, path: str) -> dict[str, Any]:
        try:
            obj = self._file[path]
            classname = str(obj.classname)
            if not classname.startswith("TH1") or len(obj.axes) != 1:
                raise ExtractionError(f"{path} is {classname}, not a one-dimensional TH1")
            values = [float(item) for item in obj.values(flow=False)]
            flow_values = [float(item) for item in obj.values(flow=True)]
            edges = [float(item) for item in obj.axis().edges()]
            sumw2_raw = list(obj.member("fSumw2"))
            entries = float(obj.member("fEntries"))
        except ExtractionError:
            raise
        except Exception as exc:
            raise ExtractionError(f"cannot read ROOT histogram {path}: {exc}") from exc
        if len(sumw2_raw) != len(values) + 2:
            raise ExtractionError(
                f"{path} has no stored Sumw2 ({len(sumw2_raw)} cells for {len(values)} bins)"
            )
        sumw2 = [float(item) for item in sumw2_raw[1:-1]]
        flow_sumw2 = [float(sumw2_raw[0]), float(sumw2_raw[-1])]
        if len(edges) != len(values) + 1 or any(
            right <= left for left, right in zip(edges, edges[1:])
        ):
            raise ExtractionError(f"{path} has malformed bin edges")
        for label, vector in (
            ("sumw", values),
            ("sumw2", sumw2),
            ("flow sumw", [flow_values[0], flow_values[-1]]),
            ("flow sumw2", flow_sumw2),
        ):
            if any(not math.isfinite(item) or item < 0.0 for item in vector):
                raise ExtractionError(f"{path} contains invalid {label}")
        if not math.isfinite(entries) or entries < 0.0:
            raise ExtractionError(f"{path} contains invalid TH1 entries")
        return {
            "bin_edges": edges,
            "sumw": values,
            "sumw2": sumw2,
            "entries": entries,
            "flow_sumw": [flow_values[0], flow_values[-1]],
            "flow_sumw2": flow_sumw2,
            "classname": classname,
        }

    def close(self) -> None:
        self._file.close()


def _resolve_object(reader: HistogramReader, basename: str, prefix: str) -> str:
    candidates = reader.paths_for_basename(basename)
    if prefix == "auto":
        if not candidates:
            raise ExtractionError(
                f"missing canonical ROOT object {basename}; class-suffixed substitutes are forbidden"
            )
        if len(candidates) != 1:
            raise ExtractionError(
                f"ambiguous ROOT object {basename}; candidates={candidates}; freeze object_prefix"
            )
        selected = candidates[0]
    else:
        selected = f"{prefix.rstrip('/')}/{basename}" if prefix else basename
        if selected not in candidates:
            raise ExtractionError(
                f"missing canonical ROOT object {selected}; observed basename matches={candidates}; "
                "class-suffixed substitutes are forbidden"
            )
    cycles = reader.cycle_count(selected)
    if cycles != 1:
        raise ExtractionError(
            f"ambiguous ROOT cycles for {selected}: observed {cycles}, expected exactly one"
        )
    return selected


def _fill_payload(
    path: Path,
    lane_id: str,
    root_sha256: str,
    required: set[str],
) -> dict[str, Any]:
    payload = _read_json(path, "fill evidence")
    if payload.get("schema") != "ppg12-stitched-purity-fill-evidence/v1":
        raise ExtractionError("unsupported fill-evidence schema")
    if payload.get("lane_id") != lane_id:
        raise ExtractionError(
            f"fill-evidence lane mismatch: expected {lane_id}, observed {payload.get('lane_id')}"
        )
    if payload.get("root_sha256") != root_sha256:
        raise ExtractionError("fill evidence was produced from a different ROOT file")
    observables = payload.get("observables")
    if not isinstance(observables, dict) or set(observables) != required:
        raise ExtractionError(
            "fill-evidence observable coverage must exactly match the lane contract"
        )
    return payload


def _candidate_links(value: Any, parent: Path, family: str) -> list[dict[str, str]]:
    if value is None:
        value = []
    if not isinstance(value, list):
        raise ExtractionError("candidate_parity_evidence must be a list")
    normalized: list[dict[str, str]] = []
    roles: set[str] = set()
    for index, row in enumerate(value):
        if not isinstance(row, dict) or not isinstance(row.get("role"), str) or not row["role"]:
            raise ExtractionError(f"candidate parity link {index}.role is required")
        role = row["role"]
        if role in roles:
            raise ExtractionError(f"duplicate candidate parity evidence role: {role}")
        roles.add(role)
        normalized.append({"role": role, **_validate_link(row, parent, f"candidate parity {role}")})
    if family == "photon":
        required = {"candidate_rows"}
        if roles != required:
            raise ExtractionError(
                "photon candidate-parity evidence roles must be exact: "
                f"expected={sorted(required)}, observed={sorted(roles)}"
            )
    elif normalized:
        raise ExtractionError("inclusive lanes must not carry photon candidate-parity evidence")
    return normalized


def _event_set_rows(path: Path) -> list[str]:
    """Read the exact nonblank, non-comment rows processed by the canary."""
    try:
        rows = [
            line.strip()
            for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ]
    except UnicodeDecodeError as exc:
        raise ExtractionError(f"event-set evidence is not a text row list: {path}") from exc
    return rows


def _require_exact_source_slices(
    lane_id: str,
    source_list_path: Path,
    event_rows: list[str],
    groups: list[dict[str, Any]],
    group_size: int,
) -> list[str]:
    """Bind selected events to exact ordered production-list group slices."""
    source_rows = _event_set_rows(source_list_path)
    if not source_rows:
        raise ExtractionError(f"{lane_id} production source list is empty")
    expected_rows: list[str] = []
    for group in groups:
        group_index = group["group_index"]
        start = group_index * group_size
        stop = start + group_size
        if stop > len(source_rows):
            raise ExtractionError(
                f"{lane_id} group {group_index} exceeds its bound production source list "
                f"({len(source_rows)} rows)"
            )
        source_slice = source_rows[start:stop]
        if group["event_rows"] != source_slice:
            raise ExtractionError(
                f"{lane_id} group {group_index} is not the exact ordered production-list slice"
            )
        expected_rows.extend(source_slice)
    if event_rows != expected_rows:
        raise ExtractionError(
            f"{lane_id} event set is not the ordered concatenation of its production-list groups"
        )
    return source_rows


def _validate_runtime_manifest(
    runtime_link: Any,
    parent: Path,
    config_link: dict[str, str],
    contract: dict[str, Any],
) -> tuple[dict[str, Any], dict[str, dict[str, Any]], dict[str, str]]:
    """Recompute lane contract hashes from the role-labelled executed files."""
    normalized_manifest_link = _validate_link(
        runtime_link, parent, "runtime_manifest"
    )
    manifest_path = Path(normalized_manifest_link["path"])
    manifest = _read_json(manifest_path, "runtime manifest")
    runtime_contract = contract.get("lane_runtime_contract")
    if not isinstance(runtime_contract, dict):
        raise ExtractionError("closure contract lacks lane_runtime_contract")
    if manifest.get("schema_version") != runtime_contract.get("manifest_schema_version"):
        raise ExtractionError("runtime manifest schema version does not match the closure contract")
    required_values = runtime_contract.get("required_manifest_values")
    if not isinstance(required_values, dict):
        raise ExtractionError("closure contract lacks required runtime-manifest values")
    for field, expected in required_values.items():
        if manifest.get(field) != expected:
            raise ExtractionError(
                f"runtime manifest {field} drift: expected {expected!r}, observed {manifest.get(field)!r}"
            )

    raw_receipt = manifest.get("build_receipt")
    receipt_sha = manifest.get("build_receipt_sha256")
    if runtime_contract.get("build_receipt_required") is not True:
        raise ExtractionError("closure contract must require a runtime build receipt")
    receipt_link = _validate_link(
        {"path": raw_receipt, "sha256": receipt_sha},
        manifest_path.parent,
        "runtime_manifest.build_receipt",
    )

    raw_files = manifest.get("files")
    if not isinstance(raw_files, list) or not raw_files:
        raise ExtractionError("runtime manifest files must be a non-empty list")
    by_role: dict[str, dict[str, str]] = {}
    seen_paths: set[str] = set()
    for index, row in enumerate(raw_files):
        if not isinstance(row, dict) or not isinstance(row.get("role"), str) or not row["role"]:
            raise ExtractionError(f"runtime manifest file {index}.role is required")
        role = row["role"]
        if role in by_role:
            raise ExtractionError(f"runtime manifest duplicates role {role}")
        link = _validate_link(row, manifest_path.parent, f"runtime manifest role {role}")
        if link["path"] in seen_paths:
            raise ExtractionError(
                f"runtime manifest aliases one physical file under multiple roles: {link['path']}"
            )
        seen_paths.add(link["path"])
        by_role[role] = {"role": role, **link}

    if "lane_config" in by_role:
        raise ExtractionError("runtime manifest must not override the lane-local executed config")
    by_role["lane_config"] = {"role": "lane_config", **config_link}
    raw_groups = runtime_contract.get("required_roles_by_lane_field")
    expected_fields = {
        "config_sha256",
        "reconstruction_sha256",
        "model_set_sha256",
        "ownership_sha256",
        "weight_sha256",
        "estimator_sha256",
    }
    if not isinstance(raw_groups, dict) or set(raw_groups) != expected_fields:
        raise ExtractionError("closure runtime role groups are incomplete")
    source_sets: dict[str, dict[str, Any]] = {}
    field_hashes: dict[str, str] = {}
    for field in sorted(expected_fields):
        roles = raw_groups[field]
        if (
            not isinstance(roles, list)
            or not roles
            or any(not isinstance(role, str) or not role for role in roles)
            or len(set(roles)) != len(roles)
        ):
            raise ExtractionError(f"closure runtime role group is malformed: {field}")
        missing = sorted(set(roles) - set(by_role))
        if missing:
            raise ExtractionError(
                f"runtime manifest lacks roles required for {field}: {missing}"
            )
        files = [by_role[role] for role in sorted(roles)]
        portable = [
            {"role": row["role"], "sha256": row["sha256"]}
            for row in files
        ]
        set_sha = _payload_sha256(portable)
        source_sets[field] = {
            "roles": sorted(roles),
            "files": files,
            "portable_role_set_sha256": set_sha,
        }
        field_hashes[field] = set_sha
    runtime_evidence = {
        "manifest": normalized_manifest_link,
        "build_receipt": receipt_link,
        "manifest_payload_sha256": _payload_sha256(manifest),
    }
    return runtime_evidence, source_sets, field_hashes


def _canonical_groups(
    lane_id: str,
    event_rows: list[str],
    group_size: int,
    group_index_start: int = 0,
) -> list[dict[str, Any]]:
    if not event_rows or len(event_rows) % group_size:
        raise ExtractionError(
            f"{lane_id} event set must contain a positive multiple of {group_size} rows"
        )
    if len(set(event_rows)) != len(event_rows):
        raise ExtractionError(f"{lane_id} event set repeats a source row")
    groups: list[dict[str, Any]] = []
    if (
        not isinstance(group_index_start, int)
        or isinstance(group_index_start, bool)
        or group_index_start < 0
    ):
        raise ExtractionError(f"{lane_id} group_index_start must be a nonnegative integer")
    for relative_index, offset in enumerate(range(0, len(event_rows), group_size)):
        group_index = group_index_start + relative_index
        rows = event_rows[offset : offset + group_size]
        group_sha = _payload_sha256(
            {
                "lane_id": lane_id,
                "group_index": group_index,
                "event_rows": rows,
            }
        )
        groups.append(
            {
                "group_index": group_index,
                "group_id": f"group-{group_index:05d}-{group_sha[:16]}",
                "event_rows": rows,
                "group_sha256": group_sha,
            }
        )
    return groups


def extract_lane(
    metadata_path: Path,
    contract_path: Path = DEFAULT_CONTRACT,
    *,
    reader_factory: Callable[[Path], HistogramReader] = UprootHistogramReader,
) -> dict[str, Any]:
    metadata_path = metadata_path.resolve()
    metadata = _read_json(metadata_path, "lane input sidecar")
    if metadata.get("schema") != "ppg12-stitched-purity-lane-input/v1":
        raise ExtractionError("unsupported lane-input schema")
    contract = _contract(contract_path.resolve())
    identity_fields = ("family", "sample", "period", "interaction")
    if any(not isinstance(metadata.get(field), str) for field in identity_fields):
        raise ExtractionError("lane input has incomplete identity")
    identity = {field: metadata[field] for field in identity_fields}
    lane_id = _lane_id(**identity)
    if metadata.get("lane_id") != lane_id or lane_id not in _expected_lanes(contract):
        raise ExtractionError(
            f"lane identity mismatch or noncanonical lane: declared={metadata.get('lane_id')!r}, canonical={lane_id}"
        )
    group_count = metadata.get("group_count")
    group_contract = contract["admission_group_contract"]
    if (
        not isinstance(group_count, int)
        or isinstance(group_count, bool)
        or group_count < int(group_contract["minimum_group_count"])
    ):
        raise ExtractionError("group_count must be a positive integer")
    scale = metadata.get("external_scale")
    if isinstance(scale, bool) or not isinstance(scale, (int, float)) or abs(float(scale) - 1.0) > 1e-12:
        raise ExtractionError(f"forbidden non-unit external scale: {scale!r}")
    if metadata.get("abcd_population") != contract["abcd_population"]:
        raise ExtractionError("lane input must use unsuffixed A/B/C/D")
    prefix = metadata.get("object_prefix", "auto")
    if not isinstance(prefix, str) or prefix.startswith("/") or ".." in prefix.split("/"):
        raise ExtractionError("object_prefix must be a relative ROOT directory or 'auto'")

    root_link = _validate_link(metadata.get("root"), metadata_path.parent, "root")
    evidence = metadata.get("evidence")
    if not isinstance(evidence, dict) or set(evidence) != set(EVIDENCE_KEYS):
        raise ExtractionError(
            f"evidence keys must be exact: expected={list(EVIDENCE_KEYS)}, observed={sorted(evidence) if isinstance(evidence, dict) else None}"
        )
    normalized_evidence = {
        name: _validate_link(evidence[name], metadata_path.parent, f"evidence.{name}")
        for name in EVIDENCE_KEYS
    }
    runtime_evidence, runtime_source_sets, runtime_field_hashes = (
        _validate_runtime_manifest(
            metadata.get("runtime_manifest"),
            metadata_path.parent,
            normalized_evidence["config"],
            contract,
        )
    )
    group_size = metadata.get("group_size")
    if group_size != int(group_contract["group_size"]):
        raise ExtractionError("every admission lane must use deterministic five-file groups")
    event_set_path = Path(normalized_evidence["event_set"]["path"])
    event_rows = _event_set_rows(event_set_path)
    event_set_row_count = len(event_rows)
    declared_rows = metadata.get("event_set_row_count")
    expected_rows = group_count * group_size
    if declared_rows != expected_rows or event_set_row_count != expected_rows:
        raise ExtractionError(
            f"lane event-set evidence must contain exactly {expected_rows} source rows"
        )
    group_index_start = metadata.get("group_index_start", 0)
    if (
        not isinstance(group_index_start, int)
        or isinstance(group_index_start, bool)
        or group_index_start < 0
        or (group_count > 1 and group_index_start != 0)
    ):
        raise ExtractionError(
            "aggregate group_index_start must be zero; component indices must be nonnegative"
        )
    groups = _canonical_groups(
        lane_id, event_rows, group_size, group_index_start
    )
    if len(groups) != group_count or metadata.get("groups") != groups:
        raise ExtractionError("lane group identities/hashes do not match the exact event set")
    group_set_sha256 = _payload_sha256(groups)
    if metadata.get("group_set_sha256") != group_set_sha256:
        raise ExtractionError("lane group-set hash is stale")
    _require_exact_source_slices(
        lane_id,
        Path(normalized_evidence["source_list"]["path"]),
        event_rows,
        groups,
        group_size,
    )
    raw_group_sidecars = metadata.get("group_input_sidecars")
    group_sidecars: list[dict[str, Any]] = []
    if group_count == 1:
        if raw_group_sidecars not in (None, []):
            raise ExtractionError(
                "single-group component sidecars must not recursively declare group inputs"
            )
        group_sidecars = [
            {
                "group_index": groups[0]["group_index"],
                "group_id": groups[0]["group_id"],
                "path": str(metadata_path),
                "sha256": _file_sha256(metadata_path),
            }
        ]
    else:
        if not isinstance(raw_group_sidecars, list) or len(raw_group_sidecars) != group_count:
            raise ExtractionError(
                "multi-group aggregates require one exact input sidecar per group"
            )
        for index, (row, group) in enumerate(zip(raw_group_sidecars, groups)):
            if (
                not isinstance(row, dict)
                or row.get("group_index") != group["group_index"]
                or row.get("group_id") != group["group_id"]
            ):
                raise ExtractionError(
                    f"group_input_sidecars[{index}] does not match its canonical group identity"
                )
            group_sidecars.append(
                {
                    "group_index": group["group_index"],
                    "group_id": group["group_id"],
                    **_validate_link(
                        row, metadata_path.parent, f"group_input_sidecars[{index}]"
                    ),
                }
            )
    root_bytes = Path(root_link["path"]).stat().st_size
    if root_bytes < int(group_contract["minimum_root_bytes"]):
        raise ExtractionError(
            f"lane ROOT is not terminal/non-tiny: {root_bytes} bytes"
        )
    fill_link = _validate_link(
        metadata.get("fill_evidence"), metadata_path.parent, "fill_evidence"
    )
    candidate_links = _candidate_links(
        metadata.get("candidate_parity_evidence"), metadata_path.parent, identity["family"]
    )
    required = set(contract["families"][identity["family"]]["required_observables"])
    fills_payload = _fill_payload(
        Path(fill_link["path"]), lane_id, root_link["sha256"], required
    )

    reader: HistogramReader | None = None
    observables: dict[str, Any] = {}
    try:
        reader = reader_factory(Path(root_link["path"]))
        for name in sorted(required):
            basename = REGION_OBJECTS[name]
            object_path = _resolve_object(reader, basename, prefix)
            histogram = reader.histogram(object_path)
            edges = _finite_vector(histogram.get("bin_edges"), f"{object_path}.bin_edges")
            sumw = _finite_vector(histogram.get("sumw"), f"{object_path}.sumw", nonnegative=True)
            sumw2 = _finite_vector(histogram.get("sumw2"), f"{object_path}.sumw2", nonnegative=True)
            flow_sumw = _finite_vector(
                histogram.get("flow_sumw"), f"{object_path}.flow_sumw", nonnegative=True
            )
            flow_sumw2 = _finite_vector(
                histogram.get("flow_sumw2"), f"{object_path}.flow_sumw2", nonnegative=True
            )
            if len(edges) != len(sumw) + 1 or len(sumw2) != len(sumw):
                raise ExtractionError(f"malformed histogram vectors for {object_path}")
            if len(flow_sumw) != 2 or len(flow_sumw2) != 2:
                raise ExtractionError(f"{object_path} lacks exact flow-bin moments")
            fill_row = fills_payload["observables"][name]
            if not isinstance(fill_row, dict) or fill_row.get("object") != object_path:
                raise ExtractionError(
                    f"fill evidence {name} does not bind exact ROOT object {object_path}"
                )
            fill_edges = _finite_vector(
                fill_row.get("bin_edges"), f"fill_evidence.{name}.bin_edges"
            )
            if not _same_vector(edges, fill_edges):
                raise ExtractionError(f"fill-evidence binning mismatch for {name}")
            fills = _count_vector(fill_row.get("fills"), f"fill_evidence.{name}.fills", len(sumw))
            flow_fills = _count_vector(
                fill_row.get("flow_fills"), f"fill_evidence.{name}.flow_fills", 2
            )
            entries = histogram.get("entries")
            if isinstance(entries, bool) or not isinstance(entries, (int, float)) or not math.isfinite(float(entries)):
                raise ExtractionError(f"{object_path} has invalid TH1 entries")
            expected_entries = sum(fills) + sum(flow_fills)
            if abs(float(entries) - expected_entries) > 1e-7 * max(1.0, expected_entries):
                raise ExtractionError(
                    f"fill evidence does not close TH1 entries for {object_path}: fills={expected_entries}, entries={entries}"
                )
            observables[name] = {
                "bin_edges": edges,
                "sumw": sumw,
                "sumw2": sumw2,
                "fills": fills,
                "object_path": object_path,
                "entries": float(entries),
                "flow": {
                    "sumw": flow_sumw,
                    "sumw2": flow_sumw2,
                    "fills": flow_fills,
                },
            }
    except ExtractionError:
        raise
    except Exception as exc:
        raise ExtractionError(f"bad ROOT input {root_link['path']}: {exc}") from exc
    finally:
        if reader is not None:
            reader.close()

    common_edges: list[float] | None = None
    for name in sorted(observables):
        edges = observables[name]["bin_edges"]
        if common_edges is None:
            common_edges = edges
        elif not _same_vector(common_edges, edges):
            raise ExtractionError(f"lane observables do not share one binning: {name}")

    output: dict[str, Any] = {
        "schema": "ppg12-stitched-purity-lane/v1",
        "lane_id": lane_id,
        **identity,
        "group_count": group_count,
        "group_size": group_size,
        "group_index_start": group_index_start,
        "event_set_row_count": event_set_row_count,
        "groups": groups,
        "group_set_sha256": group_set_sha256,
        "abcd_population": contract["abcd_population"],
        "external_scale": 1.0,
        "merge_input": {"input_id": lane_id, **root_link},
        "observables": observables,
        "extraction": {
            "schema": "ppg12-stitched-purity-lane-extraction/v1",
            "extractor": {
                "path": str(Path(__file__).resolve()),
                "sha256": _file_sha256(Path(__file__).resolve()),
            },
            "input_sidecar": {
                "path": str(metadata_path),
                "sha256": _file_sha256(metadata_path),
            },
            "root": root_link,
            "fill_evidence": fill_link,
            "evidence": normalized_evidence,
            "runtime_evidence": runtime_evidence,
            "runtime_source_sets": runtime_source_sets,
            "candidate_parity_evidence": candidate_links,
            "groups": groups,
            "group_set_sha256": group_set_sha256,
            "group_input_sidecars": group_sidecars,
            "object_prefix": prefix,
            "fill_count_semantics": "exact external per-bin counts checked against TH1 entries",
        },
    }
    for evidence_name, field_name in DIRECT_HASH_FIELD_BY_EVIDENCE.items():
        output[field_name] = normalized_evidence[evidence_name]["sha256"]
    output.update(runtime_field_hashes)
    if candidate_links:
        output["candidate_parity_evidence"] = candidate_links
    return output


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata-json", required=True, help="hash-bound lane input sidecar")
    parser.add_argument("--contract", default=str(DEFAULT_CONTRACT))
    parser.add_argument("--output", required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    output = Path(args.output).resolve()
    if output.exists():
        output.unlink()
    try:
        payload = extract_lane(
            Path(args.metadata_json), Path(args.contract)
        )
        _write_json(output, payload)
    except (ExtractionError, OSError, ValueError) as exc:
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
