#!/usr/bin/env python3
"""Extract and compare preserved PPG12 RecoEff leakage aggregates.

This is the only bridge from the preserved executable to certification of the
weighted photon-leakage cells.  A Python selection shadow is intentionally not
accepted as a substitute.  The executable is run twice: an uninstrumented
staged copy and a trace-only instrumented copy.  Their physics ROOT outputs
must be exactly equivalent before the trace or aggregate can certify anything.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Any, Iterable


SCHEMA_VERSION = 1
EVIDENCE_SOURCE = "preserved_ppg12_executable_aggregate"
RECO_EDGES = [10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 32.0, 36.0]
REGION_HISTS = {
    "A": "h_tight_iso_cluster_signal_0",
    "B": "h_tight_noniso_cluster_signal_0",
    "C": "h_nontight_iso_cluster_signal_0",
    "D": "h_nontight_noniso_cluster_signal_0",
}
CERTIFIED_SCOPES = ["tags", "isolation_abcd", "truth_abcd_fills", "weights"]


class ExtractFailure(RuntimeError):
    pass


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def require_file(path: Path, label: str) -> Path:
    if not path.is_file() or path.stat().st_size <= 0:
        raise ExtractFailure(f"{label} is missing or empty: {path}")
    return path


def root_bin_index(value: float, edges: list[float] = RECO_EDGES) -> int:
    """Return ROOT-style global bin index for an explicit variable axis."""
    if not math.isfinite(value):
        raise ExtractFailure("non-finite cluster ET in candidate aggregate")
    if value < edges[0]:
        return 0
    if value >= edges[-1]:
        return len(edges)
    return bisect.bisect_right(edges, value)


def empty_cells(edges: list[float] = RECO_EDGES) -> dict[str, list[dict[str, float]]]:
    return {
        region: [
            {"content": 0.0, "sumw2": 0.0} for _ in range(len(edges) + 1)
        ]
        for region in REGION_HISTS
    }


def aggregate_recoil_rows(
    rows: Iterable[dict[str, str]], edges: list[float] = RECO_EDGES
) -> dict[str, list[dict[str, float]]]:
    cells = empty_cells(edges)
    identities: set[str] = set()
    for row_number, row in enumerate(rows, start=2):
        identity = row.get("candidate_identity", "")
        if not identity:
            raise ExtractFailure(f"candidate row {row_number} lacks identity")
        if identity in identities:
            raise ExtractFailure(f"duplicate candidate identity: {identity}")
        identities.add(identity)
        flags: dict[str, int] = {}
        try:
            for region in REGION_HISTS:
                raw = row.get(f"rj_signal_fill_{region}", "")
                flags[region] = 0 if raw == "" else int(raw)
        except ValueError as exc:
            raise ExtractFailure(f"invalid fill flag at candidate row {row_number}") from exc
        if any(value not in (0, 1) for value in flags.values()):
            raise ExtractFailure(f"non-binary fill flag at candidate row {row_number}")
        multiplicity = sum(flags.values())
        raw_multiplicity = row.get("rj_signal_fill_multiplicity", "")
        if raw_multiplicity != "" and int(raw_multiplicity) != multiplicity:
            raise ExtractFailure(f"fill multiplicity mismatch at candidate row {row_number}")
        if multiplicity > 1:
            raise ExtractFailure(f"candidate fills multiple ABCD regions at row {row_number}")
        if multiplicity == 0:
            continue
        try:
            et = float(row["rj_cluster_Et"])
            weight = float(row["rj_weight_final"])
        except (KeyError, ValueError) as exc:
            raise ExtractFailure(f"invalid ET/weight at candidate row {row_number}") from exc
        if not math.isfinite(weight):
            raise ExtractFailure(f"non-finite weight at candidate row {row_number}")
        ibin = root_bin_index(et, edges)
        for region, enabled in flags.items():
            if enabled:
                cells[region][ibin]["content"] += weight
                cells[region][ibin]["sumw2"] += weight * weight
    return cells


def leakage_ratios(
    cells: dict[str, list[dict[str, float]]]
) -> dict[str, list[float | None]]:
    output: dict[str, list[float | None]] = {}
    for numerator in ("B", "C", "D"):
        values: list[float | None] = []
        for num, den in zip(cells[numerator], cells["A"]):
            denominator = den["content"]
            values.append(num["content"] / denominator if denominator != 0.0 else None)
        output[f"{numerator}_over_A"] = values
    return output


def close(left: float, right: float, rtol: float, atol: float) -> bool:
    return abs(left - right) <= atol + rtol * abs(right)


def compare_cells(
    oracle: dict[str, list[dict[str, float]]],
    recoil: dict[str, list[dict[str, float]]],
    *,
    rtol: float = 1.0e-6,
    atol: float = 1.0e-9,
) -> dict[str, Any]:
    records: list[dict[str, Any]] = []
    passed = True
    for region in REGION_HISTS:
        if len(oracle[region]) != len(recoil[region]):
            raise ExtractFailure(f"cell count differs for region {region}")
        for ibin, (left, right) in enumerate(zip(oracle[region], recoil[region])):
            content_pass = close(right["content"], left["content"], rtol, atol)
            sumw2_pass = close(right["sumw2"], left["sumw2"], rtol, atol)
            passed &= content_pass and sumw2_pass
            records.append(
                {
                    "region": region,
                    "global_bin": ibin,
                    "oracle_content": left["content"],
                    "recoil_content": right["content"],
                    "oracle_sumw2": left["sumw2"],
                    "recoil_sumw2": right["sumw2"],
                    "content_pass": content_pass,
                    "sumw2_pass": sumw2_pass,
                }
            )
    oracle_leak = leakage_ratios(oracle)
    recoil_leak = leakage_ratios(recoil)
    leakage_records: list[dict[str, Any]] = []
    for name in oracle_leak:
        for ibin, (left, right) in enumerate(zip(oracle_leak[name], recoil_leak[name])):
            ratio_pass = (left is None and right is None) or (
                left is not None and right is not None and close(right, left, rtol, atol)
            )
            passed &= ratio_pass
            leakage_records.append(
                {
                    "ratio": name,
                    "global_bin": ibin,
                    "oracle": left,
                    "recoil": right,
                    "pass": ratio_pass,
                }
            )
    return {"pass": passed, "cells": records, "leakage": leakage_records}


def _hist_payload(hist: Any) -> dict[str, Any]:
    dimension = int(hist.GetDimension())
    axes = []
    for axis in (hist.GetXaxis(), hist.GetYaxis(), hist.GetZaxis())[:dimension]:
        axes.append(
            [float(axis.GetBinLowEdge(i)) for i in range(1, axis.GetNbins() + 2)]
        )
    n_cells = int(hist.GetNcells())
    sumw2 = hist.GetSumw2()
    return {
        "class": hist.ClassName(),
        "dimension": dimension,
        "axes": axes,
        "content": [float(hist.GetBinContent(i)) for i in range(n_cells)],
        "error2": [float(hist.GetBinError(i)) ** 2 for i in range(n_cells)],
        "sumw2": [float(sumw2.At(i)) for i in range(int(sumw2.GetSize()))],
        "entries": float(hist.GetEntries()),
    }


def _object_payload(obj: Any) -> dict[str, Any]:
    if obj.InheritsFrom("TH1"):
        return _hist_payload(obj)
    if obj.InheritsFrom("TEfficiency"):
        return {
            "class": obj.ClassName(),
            "passed": _hist_payload(obj.GetPassedHistogram()),
            "total": _hist_payload(obj.GetTotalHistogram()),
        }
    if obj.InheritsFrom("RooUnfoldResponse"):
        payload: dict[str, Any] = {"class": obj.ClassName()}
        for name in ("Hresponse", "Hmeasured", "Htruth", "Hfakes"):
            value = getattr(obj, name)()
            payload[name] = _hist_payload(value) if value else None
        return payload
    if obj.InheritsFrom("TObjString"):
        return {"class": obj.ClassName(), "string": str(obj.GetString())}
    raise ExtractFailure(f"unsupported ROOT object class in equivalence gate: {obj.ClassName()}")


def root_payload(path: Path) -> dict[str, Any]:
    try:
        import ROOT  # type: ignore
    except ImportError as exc:
        raise ExtractFailure("PyROOT is required for executable ROOT extraction") from exc
    root_file = ROOT.TFile.Open(str(path))
    if not root_file or root_file.IsZombie():
        raise ExtractFailure(f"unreadable ROOT file: {path}")
    if root_file.TestBit(ROOT.TFile.kRecovered):
        raise ExtractFailure(f"recovered ROOT file is inadmissible: {path}")
    payload: dict[str, Any] = {}
    for key in root_file.GetListOfKeys():
        name = str(key.GetName())
        cycle = int(key.GetCycle())
        identity = f"{name};{cycle}"
        if identity in payload:
            raise ExtractFailure(f"duplicate ROOT key identity is inadmissible: {identity}")
        payload[identity] = _object_payload(key.ReadObj())
    root_file.Close()
    return payload


def extract_oracle_cells(path: Path) -> dict[str, list[dict[str, float]]]:
    try:
        import ROOT  # type: ignore
    except ImportError as exc:
        raise ExtractFailure("PyROOT is required for executable aggregate extraction") from exc
    root_file = ROOT.TFile.Open(str(path))
    if not root_file or root_file.IsZombie() or root_file.TestBit(ROOT.TFile.kRecovered):
        raise ExtractFailure(f"inadmissible RecoEff ROOT file: {path}")
    cells = empty_cells()
    for region, name in REGION_HISTS.items():
        hist = root_file.Get(name)
        if not hist or not hist.InheritsFrom("TH1"):
            raise ExtractFailure(f"missing executable leakage histogram: {name}")
        observed_edges = [
            float(hist.GetXaxis().GetBinLowEdge(i))
            for i in range(1, hist.GetNbinsX() + 2)
        ]
        if observed_edges != RECO_EDGES:
            raise ExtractFailure(f"unexpected binning for {name}: {observed_edges}")
        for ibin in range(hist.GetNbinsX() + 2):
            cells[region][ibin] = {
                "content": float(hist.GetBinContent(ibin)),
                "sumw2": float(hist.GetBinError(ibin)) ** 2,
            }
    root_file.Close()
    return cells


def parse_assets(values: list[str]) -> dict[str, Path]:
    assets: dict[str, Path] = {}
    for value in values:
        if "=" not in value:
            raise ExtractFailure(f"--asset must be ROLE=PATH: {value}")
        role, raw_path = value.split("=", 1)
        if not role or role in assets:
            raise ExtractFailure(f"duplicate/empty asset role: {role}")
        assets[role] = require_file(Path(raw_path), f"asset {role}")
    return assets


def validate_candidate_identity_binding(
    path: Path, *, lane_id: str, runtime_contract_sha256: str
) -> list[dict[str, str]]:
    """Load the paired rows only after proving their physical-lane binding."""
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ExtractFailure("paired candidate CSV contains no rows")
    for row_number, row in enumerate(rows, start=2):
        if row.get("lane_id") != lane_id:
            raise ExtractFailure(
                f"candidate row {row_number} belongs to another physical lane"
            )
        if row.get("runtime_contract_sha256") != runtime_contract_sha256:
            raise ExtractFailure(
                f"candidate row {row_number} is bound to another runtime contract"
            )
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-root", type=Path, required=True)
    parser.add_argument("--instrumented-root", type=Path, required=True)
    parser.add_argument("--baseline-response-root", type=Path, required=True)
    parser.add_argument("--instrumented-response-root", type=Path, required=True)
    parser.add_argument("--trace-csv", type=Path)
    parser.add_argument("--response-trace-csv", type=Path)
    parser.add_argument("--candidate-csv", type=Path)
    parser.add_argument("--lane-id")
    parser.add_argument("--runtime-contract", type=Path)
    parser.add_argument("--runtime-manifest", type=Path, required=True)
    parser.add_argument(
        "--root-equivalence-only",
        action="store_true",
        help="Gate the dual executable ROOT payloads before consuming trace evidence.",
    )
    parser.add_argument("--asset", action="append", default=[])
    parser.add_argument("--out-json", type=Path, required=True)
    args = parser.parse_args()

    try:
        required = {
            "baseline_root": require_file(args.baseline_root, "baseline RecoEff ROOT"),
            "instrumented_root": require_file(args.instrumented_root, "instrumented RecoEff ROOT"),
            "baseline_response_root": require_file(args.baseline_response_root, "baseline response ROOT"),
            "instrumented_response_root": require_file(args.instrumented_response_root, "instrumented response ROOT"),
            "runtime_manifest": require_file(args.runtime_manifest, "runtime manifest"),
        }
        if not args.root_equivalence_only:
            if (
                args.trace_csv is None
                or args.response_trace_csv is None
                or args.candidate_csv is None
                or not args.lane_id
                or args.runtime_contract is None
            ):
                raise ExtractFailure(
                    "full aggregate extraction requires --trace-csv, "
                    "--response-trace-csv, --candidate-csv, --lane-id, "
                    "and --runtime-contract"
                )
            required.update(
                {
                    "trace_csv": require_file(args.trace_csv, "executable candidate trace"),
                    "response_trace_csv": require_file(args.response_trace_csv, "executable response trace"),
                    "candidate_csv": require_file(args.candidate_csv, "paired candidate CSV"),
                    "runtime_contract": require_file(
                        args.runtime_contract, "paired runtime contract"
                    ),
                }
            )
        assets = parse_assets(args.asset)
        baseline_payload = root_payload(required["baseline_root"])
        instrumented_payload = root_payload(required["instrumented_root"])
        baseline_response_payload = root_payload(required["baseline_response_root"])
        instrumented_response_payload = root_payload(required["instrumented_response_root"])
        equivalence = {
            "efficiency_root_exact": baseline_payload == instrumented_payload,
            "response_root_exact": baseline_response_payload == instrumented_response_payload,
        }
        equivalence["pass"] = all(equivalence.values())
        oracle_cells = extract_oracle_cells(required["baseline_root"])
        comparison: dict[str, Any] | None = None
        recoil_cells: dict[str, list[dict[str, float]]] | None = None
        lane_identity: dict[str, str] | None = None
        if not args.root_equivalence_only:
            contract_digest = sha256(required["runtime_contract"])
            candidate_rows = validate_candidate_identity_binding(
                required["candidate_csv"],
                lane_id=args.lane_id,
                runtime_contract_sha256=contract_digest,
            )
            recoil_cells = aggregate_recoil_rows(candidate_rows)
            comparison = compare_cells(oracle_cells, recoil_cells)
            lane_identity = {
                "lane_id": args.lane_id,
                "runtime_contract_sha256": contract_digest,
            }
        status = "PASS" if equivalence["pass"] and (
            comparison is None or comparison["pass"]
        ) else "FAIL"
        provenance = {
            role: {"path": str(path.resolve()), "sha256": sha256(path)}
            for role, path in {**required, **assets}.items()
        }
        data = {
            "schema_version": SCHEMA_VERSION,
            "evidence_source": EVIDENCE_SOURCE,
            "status": status,
            "mode": "root_equivalence_only" if args.root_equivalence_only else "full",
            "certified_scopes": (
                ["instrumentation_side_channel"]
                if status == "PASS" and args.root_equivalence_only
                else CERTIFIED_SCOPES if status == "PASS" else []
            ),
            "uncertified_scopes": ["response_candidate_identity"],
            "root_equivalence": equivalence,
            "aggregate_comparison": comparison,
            "oracle_cells": oracle_cells,
            "recoil_cells": recoil_cells,
            "lane_identity": lane_identity,
            "provenance": provenance,
        }
        args.out_json.parent.mkdir(parents=True, exist_ok=True)
        args.out_json.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
        print(f"PPG12_RECOEFF_AGGREGATE_{status} output={args.out_json}")
        if status != "PASS":
            return 2
    except (OSError, ValueError, ExtractFailure) as exc:
        print(f"PPG12_RECOEFF_AGGREGATE_FAIL: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
