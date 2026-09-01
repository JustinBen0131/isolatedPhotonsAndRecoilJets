#!/usr/bin/env python3
"""Reduce accepted schema-10 simulation TTrees to additive stitching histograms.

Only the event and source-stage truth tables are read.  The reducer does not
open DSTs, run reconstruction, regenerate TTrees, or apply ownership cuts.  It
records the ownership decision as an audit while retaining every finite,
positive event maximum in the histogram or its flow counters.
"""

from __future__ import annotations

import argparse
import copy
import json
import math
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Iterator, Mapping, Sequence


SCHEMA = "THE248Schema10StitchingHistogramV1"
STATUS = "PASS"
DIRECTORY_NAME = "ReplayFoundationV1"
EVENT_TREE_NAME = "RJEventV1"
PHOTON_TREE_NAME = "RJTruthPhotonV1"
JET_TREE_NAME = "RJTruthJetV1"
BIN_MIN_GEV = 0.0
BIN_MAX_GEV = 100.0
BIN_WIDTH_GEV = 0.25
BIN_COUNT = int(round((BIN_MAX_GEV - BIN_MIN_GEV) / BIN_WIDTH_GEV))
BIN_EDGES_GEV = tuple(BIN_MIN_GEV + index * BIN_WIDTH_GEV for index in range(BIN_COUNT + 1))


@dataclass(frozen=True)
class SampleContract:
    sample_id: str
    system: str
    source_class: str
    observable: str
    ownership_low_gev: float
    ownership_high_gev: float


@dataclass(frozen=True)
class Event:
    key: str
    weight: float


@dataclass(frozen=True)
class TruthPhoton:
    event_key: str
    pt: float
    source_role: int


@dataclass(frozen=True)
class TruthJet:
    event_key: str
    algorithm: str
    radius: float
    pt: float


def _contract(
    sample_id: str,
    system: str,
    source_class: str,
    observable: str,
    low: float,
    high: float,
) -> SampleContract:
    return SampleContract(sample_id, system, source_class, observable, low, high)


SAMPLE_CONTRACTS: dict[str, SampleContract] = {
    "pp_photon5": _contract("pp_photon5", "pp", "photon_signal", "max_truth_photon_pt", 0.0, 14.0),
    "pp_photon10": _contract("pp_photon10", "pp", "photon_signal", "max_truth_photon_pt", 14.0, 22.0),
    "pp_photon20": _contract("pp_photon20", "pp", "photon_signal", "max_truth_photon_pt", 22.0, math.inf),
    "pp_jet8": _contract("pp_jet8", "pp", "inclusive_background", "max_r04_truth_jet_pt", 9.0, 14.0),
    "pp_jet12": _contract("pp_jet12", "pp", "inclusive_background", "max_r04_truth_jet_pt", 14.0, 21.0),
    "pp_jet20": _contract("pp_jet20", "pp", "inclusive_background", "max_r04_truth_jet_pt", 21.0, 32.0),
    "pp_jet30": _contract("pp_jet30", "pp", "inclusive_background", "max_r04_truth_jet_pt", 32.0, 42.0),
    "pp_jet40": _contract("pp_jet40", "pp", "inclusive_background", "max_r04_truth_jet_pt", 42.0, 100.0),
    "auau_photon12": _contract("auau_photon12", "auau", "photon_signal", "max_truth_photon_pt", 12.0, 21.0),
    "auau_photon20": _contract("auau_photon20", "auau", "photon_signal", "max_truth_photon_pt", 21.0, math.inf),
    "auau_jet12": _contract("auau_jet12", "auau", "inclusive_background", "max_r04_truth_jet_pt", 12.0, 21.0),
    "auau_jet20": _contract("auau_jet20", "auau", "inclusive_background", "max_r04_truth_jet_pt", 21.0, 31.0),
    "auau_jet30": _contract("auau_jet30", "auau", "inclusive_background", "max_r04_truth_jet_pt", 31.0, 41.0),
    "auau_jet40": _contract("auau_jet40", "auau", "inclusive_background", "max_r04_truth_jet_pt", 41.0, math.inf),
}


def identity_key(high: Any, low: Any) -> str:
    """Retain both 64-bit halves without routing them through a float."""
    return f"{int(high):016x}:{int(low):016x}"


def sample_contract(sample_id: str) -> SampleContract:
    try:
        return SAMPLE_CONTRACTS[sample_id]
    except KeyError as error:
        raise ValueError(f"unsupported stitching sample {sample_id}") from error


def _ownership_payload(contract: SampleContract) -> dict[str, Any]:
    finite_high = math.isfinite(contract.ownership_high_gev)
    return {
        "low_gev": contract.ownership_low_gev,
        "high_gev": contract.ownership_high_gev if finite_high else None,
        "high_is_unbounded": not finite_high,
        "lower_edge_inclusive": True,
        "upper_edge_inclusive": True,
    }


def ownership_contains(value: float, contract: SampleContract) -> bool:
    """Match the current implementation's inclusive finite upper edge."""
    if not math.isfinite(value) or value < contract.ownership_low_gev:
        return False
    return math.isinf(contract.ownership_high_gev) or value <= contract.ownership_high_gev


def _flow() -> dict[str, float | int]:
    return {"raw_count": 0, "sumw": 0.0, "sumw2": 0.0}


def _reader_contract(contract: SampleContract) -> dict[str, Any]:
    observable_tree = PHOTON_TREE_NAME if contract.observable == "max_truth_photon_pt" else JET_TREE_NAME
    observable_branches = (
        ["event_id_hi", "event_id_lo", "pt", "source_role"]
        if observable_tree == PHOTON_TREE_NAME
        else ["event_id_hi", "event_id_lo", "algorithm", "radius", "pt"]
    )
    return {
        "directory": DIRECTORY_NAME,
        "event_tree": EVENT_TREE_NAME,
        "event_branches": ["event_id_hi", "event_id_lo", "event_weight"],
        "observable_tree": observable_tree,
        "observable_branches": observable_branches,
        "dst_reads": 0,
        "ttree_regeneration": False,
    }


def _empty_payload(contract: SampleContract) -> dict[str, Any]:
    return {
        "schema": SCHEMA,
        "status": STATUS,
        "sample_id": contract.sample_id,
        "system": contract.system,
        "source_class": contract.source_class,
        "observable": contract.observable,
        "ownership_window": _ownership_payload(contract),
        "ownership_audit": {
            "application": "audit_only_production_ownership_preserved",
            "entries_dropped": 0,
            "events_in_window": 0,
            "events_outside_window": 0,
            "events_on_lower_edge": 0,
            "events_on_upper_edge": 0,
        },
        "binning": {
            "minimum_gev": BIN_MIN_GEV,
            "maximum_gev": BIN_MAX_GEV,
            "width_gev": BIN_WIDTH_GEV,
            "edges_gev": list(BIN_EDGES_GEV),
        },
        "histogram": {
            "raw_counts": [0] * BIN_COUNT,
            "sumw": [0.0] * BIN_COUNT,
            "sumw2": [0.0] * BIN_COUNT,
            "underflow": _flow(),
            "overflow": _flow(),
            "missing": _flow(),
        },
        "event_totals": {
            "generated_events": 0,
            "events_with_observable": 0,
            "events_in_histogram_range": 0,
            "events_underflow": 0,
            "events_overflow": 0,
            "events_missing_observable": 0,
            "events_with_finite_weight": 0,
            "events_with_weight_fallback": 0,
            "sum_event_weights": 0.0,
            "sum_event_weights2": 0.0,
        },
        "truth_record_totals": {
            "records_read": 0,
            "eligible_records": 0,
        },
        "reader_contract": _reader_contract(contract),
        "inputs": [],
        "input_count": 0,
        "dst_reads": 0,
        "ttree_regeneration": False,
    }


def _accumulate_flow(flow: dict[str, float | int], weight: float) -> None:
    flow["raw_count"] = int(flow["raw_count"]) + 1
    flow["sumw"] = float(flow["sumw"]) + weight
    flow["sumw2"] = float(flow["sumw2"]) + weight * weight


def _eligible_pt(record: TruthPhoton | TruthJet, contract: SampleContract) -> float | None:
    if contract.observable == "max_truth_photon_pt":
        if not isinstance(record, TruthPhoton):
            raise TypeError("photon sample requires TruthPhoton records")
        return record.pt if record.source_role == 1 and math.isfinite(record.pt) and record.pt > 0.0 else None
    if not isinstance(record, TruthJet):
        raise TypeError("jet sample requires TruthJet records")
    valid = (
        record.algorithm == "antikt"
        and math.isfinite(record.radius)
        and math.isclose(record.radius, 0.4, rel_tol=0.0, abs_tol=1.0e-9)
        and math.isfinite(record.pt)
        and record.pt > 0.0
    )
    return record.pt if valid else None


def reduce_records(
    sample_id: str,
    events: Iterable[Event],
    truth_records: Iterable[TruthPhoton | TruthJet],
) -> dict[str, Any]:
    """Pure event reduction used by both the ROOT reader and unit tests."""
    contract = sample_contract(sample_id)
    event_rows: dict[str, Event] = {}
    for event in events:
        if event.key in event_rows:
            raise ValueError(f"duplicate event identity {event.key}")
        event_rows[event.key] = event

    maxima: dict[str, float] = {}
    records_read = 0
    eligible_records = 0
    for record in truth_records:
        records_read += 1
        if record.event_key not in event_rows:
            raise ValueError(f"truth record references unknown event {record.event_key}")
        pt = _eligible_pt(record, contract)
        if pt is None:
            continue
        eligible_records += 1
        previous = maxima.get(record.event_key)
        if previous is None or pt > previous:
            maxima[record.event_key] = pt

    payload = _empty_payload(contract)
    payload["truth_record_totals"] = {
        "records_read": records_read,
        "eligible_records": eligible_records,
    }
    totals = payload["event_totals"]
    audit = payload["ownership_audit"]
    histogram = payload["histogram"]
    totals["generated_events"] = len(event_rows)

    for event in event_rows.values():
        finite_weight = math.isfinite(event.weight)
        weight = event.weight if finite_weight else 1.0
        totals["events_with_finite_weight" if finite_weight else "events_with_weight_fallback"] += 1
        totals["sum_event_weights"] += weight
        totals["sum_event_weights2"] += weight * weight

        pt = maxima.get(event.key)
        if pt is None:
            totals["events_missing_observable"] += 1
            _accumulate_flow(histogram["missing"], weight)
            continue

        totals["events_with_observable"] += 1
        in_window = ownership_contains(pt, contract)
        audit["events_in_window" if in_window else "events_outside_window"] += 1
        if math.isclose(pt, contract.ownership_low_gev, rel_tol=0.0, abs_tol=1.0e-12):
            audit["events_on_lower_edge"] += 1
        if math.isfinite(contract.ownership_high_gev) and math.isclose(
            pt, contract.ownership_high_gev, rel_tol=0.0, abs_tol=1.0e-12
        ):
            audit["events_on_upper_edge"] += 1

        if pt < BIN_MIN_GEV:
            totals["events_underflow"] += 1
            _accumulate_flow(histogram["underflow"], weight)
            continue
        if pt >= BIN_MAX_GEV:
            totals["events_overflow"] += 1
            _accumulate_flow(histogram["overflow"], weight)
            continue

        index = int((pt - BIN_MIN_GEV) / BIN_WIDTH_GEV)
        histogram["raw_counts"][index] += 1
        histogram["sumw"][index] += weight
        histogram["sumw2"][index] += weight * weight
        totals["events_in_histogram_range"] += 1

    return payload


def _add_numeric_mapping(destination: dict[str, Any], source: Mapping[str, Any], keys: Iterable[str]) -> None:
    for key in keys:
        destination[key] += source[key]


def merge_histograms(payloads: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """Merge equal-contract histogram payloads and reject repeated inputs."""
    if not payloads:
        raise ValueError("at least one stitching histogram is required")
    merged = copy.deepcopy(payloads[0])
    reference = payloads[0]
    contract_fields = (
        "schema",
        "status",
        "sample_id",
        "system",
        "source_class",
        "observable",
        "ownership_window",
        "binning",
        "reader_contract",
        "dst_reads",
        "ttree_regeneration",
    )
    seen_paths: set[str] = set()
    merged_inputs: list[dict[str, Any]] = []

    for payload in payloads:
        for field in contract_fields:
            if payload.get(field) != reference.get(field):
                raise ValueError(f"incompatible stitching histogram field {field}")
        for item in payload.get("inputs", []):
            if not isinstance(item, Mapping) or not isinstance(item.get("path"), str):
                raise ValueError("invalid stitching input receipt")
            path = str(Path(item["path"]).resolve())
            if path in seen_paths:
                raise ValueError(f"duplicate input path {path}")
            seen_paths.add(path)
            merged_inputs.append({"path": path, "size_bytes": int(item["size_bytes"])})

    merged["inputs"] = merged_inputs
    merged["input_count"] = len(merged_inputs)
    merged["histogram"] = copy.deepcopy(reference["histogram"])
    merged["event_totals"] = copy.deepcopy(reference["event_totals"])
    merged["truth_record_totals"] = copy.deepcopy(reference["truth_record_totals"])
    merged["ownership_audit"] = copy.deepcopy(reference["ownership_audit"])

    for payload in payloads[1:]:
        for field in ("raw_counts", "sumw", "sumw2"):
            source_values = payload["histogram"][field]
            destination_values = merged["histogram"][field]
            if len(source_values) != len(destination_values):
                raise ValueError(f"incompatible stitching histogram shape {field}")
            merged["histogram"][field] = [left + right for left, right in zip(destination_values, source_values)]
        for flow_name in ("underflow", "overflow", "missing"):
            _add_numeric_mapping(
                merged["histogram"][flow_name],
                payload["histogram"][flow_name],
                ("raw_count", "sumw", "sumw2"),
            )
        _add_numeric_mapping(merged["event_totals"], payload["event_totals"], merged["event_totals"])
        _add_numeric_mapping(
            merged["truth_record_totals"], payload["truth_record_totals"], merged["truth_record_totals"]
        )
        for key in (
            "events_in_window",
            "events_outside_window",
            "events_on_lower_edge",
            "events_on_upper_edge",
            "entries_dropped",
        ):
            merged["ownership_audit"][key] += payload["ownership_audit"][key]
    return merged


def _required_tree(directory: Any, name: str) -> Any:
    tree = directory.Get(name)
    if tree is None or not tree:
        raise RuntimeError(f"missing required tree {DIRECTORY_NAME}/{name}")
    return tree


def _leaf_rows(tree: Any, fields: Sequence[str]) -> Iterator[dict[str, Any]]:
    leaves = {field: tree.GetLeaf(field) for field in fields}
    missing = [field for field, leaf in leaves.items() if leaf is None]
    if missing:
        raise RuntimeError(f"{tree.GetName()} missing branches: {', '.join(missing)}")
    for index in range(int(tree.GetEntries())):
        if tree.GetEntry(index) <= 0:
            raise RuntimeError(f"failed to read {tree.GetName()} entry {index}")
        yield {field: getattr(tree, field) for field in fields}


def _string_leaf_rows(
    tree: Any,
    scalar_fields: Sequence[str],
    string_fields: Sequence[str],
) -> Iterator[dict[str, Any]]:
    leaves = {field: tree.GetLeaf(field) for field in scalar_fields}
    branches = {field: tree.GetBranch(field) for field in string_fields}
    missing = [field for field, leaf in leaves.items() if leaf is None]
    missing.extend(field for field, branch in branches.items() if branch is None)
    if missing:
        raise RuntimeError(f"{tree.GetName()} missing branches: {', '.join(missing)}")
    for index in range(int(tree.GetEntries())):
        if tree.GetEntry(index) <= 0:
            raise RuntimeError(f"failed to read {tree.GetName()} entry {index}")
        row = {field: getattr(tree, field) for field in scalar_fields}
        row.update({field: str(getattr(tree, field)) for field in string_fields})
        yield row


def read_schema10_root(path: Path, sample_id: str) -> dict[str, Any]:
    """Read the two required ReplayFoundationV1 trees from one ROOT file."""
    contract = sample_contract(sample_id)
    try:
        import ROOT  # type: ignore[import-not-found]
    except ImportError as error:  # pragma: no cover - exercised in the analysis runtime
        raise RuntimeError("PyROOT is required for schema-10 stitching reduction") from error

    root_file = ROOT.TFile.Open(str(path), "READ")
    if root_file is None or not root_file or root_file.IsZombie():
        raise RuntimeError(f"cannot open ROOT input {path}")
    try:
        directory = root_file.Get(DIRECTORY_NAME)
        if directory is None or not directory:
            raise RuntimeError(f"{path} has no {DIRECTORY_NAME} directory")
        event_tree = _required_tree(directory, EVENT_TREE_NAME)
        events = [
            Event(
                identity_key(row["event_id_hi"], row["event_id_lo"]),
                float(row["event_weight"]),
            )
            for row in _leaf_rows(event_tree, ("event_id_hi", "event_id_lo", "event_weight"))
        ]
        if contract.observable == "max_truth_photon_pt":
            truth_tree = _required_tree(directory, PHOTON_TREE_NAME)
            truth_records: list[TruthPhoton | TruthJet] = [
                TruthPhoton(
                    identity_key(row["event_id_hi"], row["event_id_lo"]),
                    float(row["pt"]),
                    int(row["source_role"]),
                )
                for row in _leaf_rows(truth_tree, ("event_id_hi", "event_id_lo", "pt", "source_role"))
            ]
        else:
            truth_tree = _required_tree(directory, JET_TREE_NAME)
            truth_records = [
                TruthJet(
                    identity_key(row["event_id_hi"], row["event_id_lo"]),
                    row["algorithm"],
                    float(row["radius"]),
                    float(row["pt"]),
                )
                for row in _string_leaf_rows(
                    truth_tree,
                    ("event_id_hi", "event_id_lo", "radius", "pt"),
                    ("algorithm",),
                )
            ]
        return reduce_records(sample_id, events, truth_records)
    finally:
        root_file.Close()


def reduce_root_files(sample_id: str, inputs: Sequence[Path]) -> dict[str, Any]:
    if not inputs:
        raise ValueError("at least one accepted schema-10 input is required")
    canonical_inputs: list[Path] = []
    seen_paths: set[str] = set()
    for supplied in inputs:
        path = supplied.expanduser().resolve()
        canonical = str(path)
        if canonical in seen_paths:
            raise ValueError(f"duplicate input path {canonical}")
        seen_paths.add(canonical)
        if not path.is_file():
            raise ValueError(f"accepted schema-10 input is not a regular file: {path}")
        canonical_inputs.append(path)

    payloads: list[dict[str, Any]] = []
    for path in canonical_inputs:
        size_before = path.stat().st_size
        payload = read_schema10_root(path, sample_id)
        size_after = path.stat().st_size
        if size_after != size_before:
            raise RuntimeError(f"schema-10 input size changed while reading: {path}")
        payload["inputs"] = [{"path": str(path), "size_bytes": size_after}]
        payload["input_count"] = 1
        payloads.append(payload)
    return merge_histograms(payloads)


def canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"


def write_output(payload: Mapping[str, Any], output: str) -> None:
    serialized = canonical_json(payload)
    if output == "-":
        sys.stdout.write(serialized)
        return
    path = Path(output).expanduser().resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o640)
    with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
        stream.write(serialized)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-id", required=True, choices=sorted(SAMPLE_CONTRACTS))
    parser.add_argument("--input", action="append", type=Path, required=True)
    parser.add_argument("--output", required=True, help="output JSON path, or - for canonical JSON on stdout")
    args = parser.parse_args(argv)
    payload = reduce_root_files(args.sample_id, args.input)
    write_output(payload, args.output)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
