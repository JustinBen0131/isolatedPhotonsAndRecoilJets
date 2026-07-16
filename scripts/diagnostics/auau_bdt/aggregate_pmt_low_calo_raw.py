#!/usr/bin/env python3
"""Aggregate THE-95 PMT diagnostics without ROOT labeled-axis merging."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any

import boost_histogram as bh
import numpy as np
import uproot


PREFIX = "SIM/"
REQUIRED = {
    "SIM/h_pmtDiag_audit",
    "SIM/h3_pmtDiag_mbdPmtOccupancyByChannel",
    "SIM/h3_pmtDiag_mbdPmtChargeByChannel",
    "SIM/h3_pmtDiag_mbdNFiredTotalByClass",
    "SIM/h3_pmtDiag_mbdChargeTotalByClass",
    "SIM/h3_pmtDiag_mbdChargeAsymmetryByClass",
    "SIM/h3_pmtDiag_mbdNFiredSouthVsNorth",
    "SIM/h3_pmtDiag_mbdChargeSouthVsNorth",
    "SIM/h3_pmtDiag_logTotalCaloVsMbdCharge",
    "SIM/h3_pmtDiag_emcalEnergyVsMbdCharge",
    "SIM/h3_pmtDiag_ihcalEnergyVsMbdCharge",
    "SIM/h3_pmtDiag_ohcalEnergyVsMbdCharge",
    "SIM/h3_pmtDiag_totalCaloEnergyVsMbdCharge",
}
SPLIT_BASES = (
    "h3_pmtDiag_mbdPmtOccupancyByChannel",
    "h3_pmtDiag_mbdPmtChargeByChannel",
    "h3_pmtDiag_mbdNFiredTotalByClass",
    "h3_pmtDiag_mbdChargeTotalByClass",
    "h3_pmtDiag_mbdChargeAsymmetryByClass",
    "h3_pmtDiag_logTotalCaloVsMbdCharge",
    "h3_pmtDiag_totalCaloEnergyVsMbdCharge",
)
for _base in SPLIT_BASES:
    REQUIRED.add(f"SIM/{_base}_mbPass")
    REQUIRED.add(f"SIM/{_base}_mbFail")

CLASS_AXIS_BASES = {
    "h3_pmtDiag_mbdPmtOccupancyByChannel",
    "h3_pmtDiag_mbdPmtChargeByChannel",
    "h3_pmtDiag_mbdNFiredTotalByClass",
    "h3_pmtDiag_mbdChargeTotalByClass",
    "h3_pmtDiag_mbdChargeAsymmetryByClass",
}


@dataclass(frozen=True)
class Component:
    name: str
    weight: float
    directory: Path


def parse_component(text: str) -> Component:
    fields = text.split("=", 2)
    if len(fields) != 3:
        raise argparse.ArgumentTypeError("component must be NAME=WEIGHT=DIRECTORY")
    name, weight_text, directory = fields
    weight = float(weight_text)
    if not name or not np.isfinite(weight) or weight <= 0.0:
        raise argparse.ArgumentTypeError(f"invalid component: {text}")
    return Component(name, weight, Path(directory))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--component", action="append", required=True, type=parse_component)
    parser.add_argument("--expected-files", type=int, default=1000)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=6)
    return parser.parse_args()


def selected_key(key: str) -> bool:
    return (
        key == "SIM/h_pmtDiag_audit"
        or key.startswith("SIM/h3_pmtDiag_")
    )


def read_object(obj: Any) -> dict[str, Any]:
    values = np.asarray(obj.values(flow=False), dtype=np.float64)
    raw_variances = obj.variances(flow=False)
    variances = np.abs(values) if raw_variances is None else np.asarray(raw_variances, dtype=np.float64)
    edges = tuple(np.asarray(axis.edges(flow=False), dtype=np.float64) for axis in obj.axes)
    return {"values": values, "variances": variances, "edges": edges}


def check_compatible(key: str, left: dict[str, Any], right: dict[str, Any], path: Path) -> None:
    if left["values"].shape != right["values"].shape:
        raise RuntimeError(f"{path}: shape drift for {key}: {right['values'].shape} != {left['values'].shape}")
    if len(left["edges"]) != len(right["edges"]):
        raise RuntimeError(f"{path}: dimension drift for {key}")
    for expected, observed in zip(left["edges"], right["edges"]):
        if not np.array_equal(expected, observed):
            raise RuntimeError(f"{path}: axis drift for {key}")


def process_component(component: Component, expected_files: int) -> dict[str, Any]:
    files = sorted(component.directory.glob("*.root"))
    if len(files) != expected_files:
        raise RuntimeError(f"{component.name}: expected {expected_files} ROOTs, found {len(files)}")
    accumulated: dict[str, dict[str, Any]] = {}
    total_bytes = 0
    for index, path in enumerate(files, start=1):
        if path.stat().st_size < 1024:
            raise RuntimeError(f"{component.name}: tiny ROOT: {path}")
        total_bytes += path.stat().st_size
        with uproot.open(path) as root_file:
            keys = [key for key in root_file.keys(recursive=True, cycle=False) if selected_key(key)]
            for key in keys:
                current = read_object(root_file[key])
                short_name = key.removeprefix(PREFIX).split("_mb", 1)[0]
                if short_name in CLASS_AXIS_BASES:
                    expected_edges = np.asarray([0.5, 1.5, 2.5, 3.5])
                    if len(current["edges"]) != 3 or not np.array_equal(current["edges"][2], expected_edges):
                        raise RuntimeError(f"{path}: corrupted all/low/main axis for {key}: {current['edges'][-1]}")
                if key not in accumulated:
                    accumulated[key] = current
                else:
                    check_compatible(key, accumulated[key], current, path)
                    accumulated[key]["values"] += current["values"]
                    accumulated[key]["variances"] += current["variances"]
        if index % 100 == 0:
            print(f"[{component.name}] {index}/{len(files)}", flush=True)
    missing = sorted(REQUIRED - set(accumulated))
    if missing:
        raise RuntimeError(f"{component.name}: missing required objects: {missing}")
    return {
        "name": component.name,
        "weight": component.weight,
        "directory": str(component.directory),
        "files": len(files),
        "bytes": total_bytes,
        "histograms": accumulated,
    }


def combine(results: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    combined: dict[str, dict[str, Any]] = {}
    for result in results:
        weight = float(result["weight"])
        for key, hist in result["histograms"].items():
            if key not in combined:
                combined[key] = {
                    "values": weight * hist["values"],
                    "variances": weight * weight * hist["variances"],
                    "edges": hist["edges"],
                }
            else:
                check_compatible(key, combined[key], hist, Path(result["directory"]))
                combined[key]["values"] += weight * hist["values"]
                combined[key]["variances"] += weight * weight * hist["variances"]
    missing = sorted(REQUIRED - set(combined))
    if missing:
        raise RuntimeError(f"combined output missing required objects: {missing}")
    return combined


def write_root(path: Path, histograms: dict[str, dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with uproot.recreate(path) as output:
        for key, payload in sorted(histograms.items()):
            axes = [bh.axis.Variable(edges) for edges in payload["edges"]]
            hist = bh.Histogram(*axes, storage=bh.storage.Weight())
            view = hist.view(flow=False)
            view.value[...] = payload["values"]
            view.variance[...] = payload["variances"]
            output[key] = hist


def main() -> int:
    args = parse_args()
    components = list(args.component)
    if len({component.name for component in components}) != len(components):
        raise RuntimeError("component names must be unique")
    results: list[dict[str, Any]] = []
    with ProcessPoolExecutor(max_workers=min(args.workers, len(components))) as executor:
        futures = {
            executor.submit(process_component, component, args.expected_files): component.name
            for component in components
        }
        for future in as_completed(futures):
            result = future.result()
            print(f"[{result['name']}] complete", flush=True)
            results.append(result)
    results.sort(key=lambda row: row["name"])
    combined = combine(results)
    write_root(args.output_root, combined)
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    manifest = {
        "schema": "AUAU_PMT_RAW_DIAGNOSTIC_AGGREGATE_V1",
        "output_root": str(args.output_root),
        "expected_files_per_component": args.expected_files,
        "weighting": "producer event weight retained; canonical finalStitch sample weight applied once",
        "categorical_axis_contract": "all/low/main axis must have exactly three bins with edges 0.5,1.5,2.5,3.5",
        "components": [
            {key: value for key, value in result.items() if key != "histograms"}
            for result in results
        ],
        "histogram_count": len(combined),
    }
    args.manifest.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
