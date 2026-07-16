#!/usr/bin/env python3
"""Extract the exact 11-bin signal ABCD family from two audited ROOT files."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import ROOT


OBJECT_NAMES = {
    "A": "h_tight_iso_cluster_signal_0",
    "B": "h_tight_noniso_cluster_signal_0",
    "C": "h_nontight_iso_cluster_signal_0",
    "D": "h_nontight_noniso_cluster_signal_0",
}
PPG12_OBJECTS = {region: name for region, name in OBJECT_NAMES.items()}
CURRENT_OBJECTS = {region: f"SIM/{name}" for region, name in OBJECT_NAMES.items()}
EDGES = [10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def open_root(path: Path) -> ROOT.TFile:
    handle = ROOT.TFile.Open(str(path), "READ")
    if not handle or handle.IsZombie() or handle.TestBit(ROOT.TFile.kRecovered):
        raise RuntimeError(f"unreadable, zombie, or recovered ROOT: {path}")
    return handle


def extract(handle: ROOT.TFile, label: str, objects: dict[str, str]) -> dict[str, ROOT.TH1]:
    result: dict[str, ROOT.TH1] = {}
    for region, path in objects.items():
        hist = handle.Get(path)
        if not hist:
            raise RuntimeError(f"missing {path} in {label}")
        if hist.GetNbinsX() != 11:
            raise RuntimeError(f"{path} has {hist.GetNbinsX()} rather than 11 bins")
        got = [hist.GetXaxis().GetBinLowEdge(1)]
        got.extend(hist.GetXaxis().GetBinUpEdge(index) for index in range(1, 12))
        if any(abs(left - right) > 1.0e-12 for left, right in zip(got, EDGES)):
            raise RuntimeError(f"{path} bin edges do not match {EDGES}: {got}")
        result[region] = hist
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ppg12-root", type=Path, required=True)
    parser.add_argument("--current-root", type=Path, required=True)
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--expected-ppg12-sha256", required=True)
    args = parser.parse_args()

    ppg_sha = sha256(args.ppg12_root)
    if ppg_sha != args.expected_ppg12_sha256:
        raise RuntimeError(f"PPG12 SHA mismatch: {ppg_sha}")
    current_sha = sha256(args.current_root)
    ppg_file = open_root(args.ppg12_root)
    current_file = open_root(args.current_root)
    ppg = extract(ppg_file, "PPG12", PPG12_OBJECTS)
    current = extract(current_file, "current", CURRENT_OBJECTS)

    rows: list[dict[str, float | int]] = []
    for index in range(1, 12):
        row: dict[str, float | int] = {
            "bin_index": index,
            "pt_lo": EDGES[index - 1],
            "pt_hi": EDGES[index],
        }
        for source, histograms in (("ppg12", ppg), ("current", current)):
            for region in OBJECT_NAMES:
                row[f"{source}_{region}"] = float(histograms[region].GetBinContent(index))
                row[f"{source}_{region}_err"] = float(histograms[region].GetBinError(index))
        rows.append(row)

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    args.manifest.write_text(json.dumps({
        "schema": "THE97_SIGNAL_ABCD_EXTRACT_V1",
        "status": "neutral noncanonical candidate extraction",
        "ppg12_root": str(args.ppg12_root),
        "ppg12_root_sha256": ppg_sha,
        "current_root": str(args.current_root),
        "current_root_sha256": current_sha,
        "objects": {
            "ppg12": PPG12_OBJECTS,
            "current": CURRENT_OBJECTS,
        },
        "bin_edges_gev": EDGES,
        "rows": len(rows),
    }, indent=2) + "\n")
    print(json.dumps({
        "output_csv": str(args.output_csv),
        "manifest": str(args.manifest),
        "ppg12_root_sha256": ppg_sha,
        "current_root_sha256": current_sha,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
