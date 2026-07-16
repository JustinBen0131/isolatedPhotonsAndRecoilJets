#!/usr/bin/env python3
"""Audit a fixed-order additive RecoilJets SIM merge with ROOT-native arithmetic."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ROOT


def open_root(path: str):
    handle = ROOT.TFile.Open(path, "READ")
    if not handle or handle.IsZombie() or handle.TestBit(ROOT.TFile.kRecovered):
        raise RuntimeError(f"unreadable, zombie, or recovered ROOT: {path}")
    return handle


def collect_histograms(directory, prefix: str = "") -> dict[str, object]:
    result: dict[str, object] = {}
    for key in directory.GetListOfKeys():
        name = key.GetName()
        path = f"{prefix}/{name}" if prefix else name
        obj = directory.Get(name)
        if obj.InheritsFrom("TDirectory"):
            result.update(collect_histograms(obj, path))
        elif obj.InheritsFrom("TH1"):
            result[path] = obj
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--input", action="append", required=True, dest="inputs")
    parser.add_argument("--json", dest="json_path")
    parser.add_argument("--required-token", action="append", default=[])
    args = parser.parse_args()

    if len(args.inputs) < 2:
        raise SystemExit("at least two --input files are required")

    files = [open_root(path) for path in [args.output, *args.inputs]]
    output_hists = collect_histograms(files[0])
    input_hists = [collect_histograms(handle) for handle in files[1:]]
    failures: list[str] = []
    max_content_delta = 0.0
    max_sumw2_delta = 0.0
    compared_cells = 0
    compared_sumw2 = 0

    union_paths = set().union(*(set(items) for items in input_hists))
    if set(output_hists) != union_paths:
        missing = sorted(union_paths - set(output_hists))
        extra = sorted(set(output_hists) - union_paths)
        failures.append(f"histogram key union mismatch: missing={missing[:10]} extra={extra[:10]}")

    for path in sorted(union_paths & set(output_hists)):
        output_hist = output_hists[path]
        present = [(index, items[path]) for index, items in enumerate(input_hists) if path in items]
        components = [hist for _, hist in present]
        if any(hist.GetNcells() != output_hist.GetNcells() for hist in components):
            failures.append(f"cell-count mismatch: {path}")
            continue
        first_index, first_hist = present[0]
        expected_hist = first_hist.Clone(f"audit_expected_{abs(hash(path))}")
        expected_hist.SetDirectory(0)
        for index, items in enumerate(input_hists):
            if index <= first_index or path not in items:
                continue
            component = items[path]
            expected_hist.Add(component)
        for index in range(output_hist.GetNcells()):
            expected = float(expected_hist.GetBinContent(index))
            actual = float(output_hist.GetBinContent(index))
            delta = abs(actual - expected)
            max_content_delta = max(max_content_delta, delta)
            tolerance = max(1.0e-12, 1.0e-12 * abs(expected))
            if not math.isfinite(actual) or delta > tolerance:
                failures.append(
                    f"content mismatch {path} cell={index} actual={actual} expected={expected} delta={delta}"
                )
                if len(failures) >= 100:
                    break
            compared_cells += 1

        output_sumw2 = output_hist.GetSumw2()
        expected_sumw2 = expected_hist.GetSumw2()
        if output_sumw2.GetSize() or expected_sumw2.GetSize():
            if output_sumw2.GetSize() != expected_sumw2.GetSize():
                failures.append(
                    f"Sumw2-storage mismatch: {path} "
                    f"actual_size={output_sumw2.GetSize()} expected_size={expected_sumw2.GetSize()}"
                )
                continue
            if output_sumw2.GetSize() != output_hist.GetNcells():
                failures.append(f"Sumw2-size mismatch: {path}")
                continue
            for index in range(output_hist.GetNcells()):
                expected = float(expected_sumw2.At(index))
                actual = float(output_sumw2.At(index))
                delta = abs(actual - expected)
                max_sumw2_delta = max(max_sumw2_delta, delta)
                tolerance = max(1.0e-12, 1.0e-12 * abs(expected))
                if not math.isfinite(actual) or delta > tolerance:
                    failures.append(
                        f"Sumw2 mismatch {path} cell={index} actual={actual} expected={expected} delta={delta}"
                    )
                    if len(failures) >= 100:
                        break
                compared_sumw2 += 1
        if len(failures) >= 100:
            break

    key_text = "\n".join(sorted(output_hists))
    required_tokens = args.required_token
    token_presence = {token: token.lower() in key_text.lower() for token in required_tokens}
    for token, present in token_presence.items():
        if not present:
            failures.append(f"required object-family token absent: {token}")

    config = files[0].Get("analysis_config_yaml")
    config_present = bool(config)
    if not config_present:
        failures.append("analysis_config_yaml absent")

    report = {
        "status": "PASS" if not failures else "FAIL",
        "output": args.output,
        "inputs_fixed_order": args.inputs,
        "output_bytes": Path(args.output).stat().st_size,
        "histograms": len(output_hists),
        "compared_cells": compared_cells,
        "compared_sumw2_cells": compared_sumw2,
        "max_content_delta": max_content_delta,
        "max_sumw2_delta": max_sumw2_delta,
        "analysis_config_yaml_present": config_present,
        "required_object_family_tokens": token_presence,
        "failures": failures,
    }
    rendered = json.dumps(report, indent=2, sort_keys=True)
    print(rendered)
    if args.json_path:
        Path(args.json_path).write_text(rendered + "\n", encoding="utf-8")
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
