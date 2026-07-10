#!/usr/bin/env python3
"""Validate RecoilJets AuAu embedded truth-isolation ROOT outputs."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import ROOT


CENTRALITIES = ("cent_0_20", "cent_20_50", "cent_50_80")
PT_INTERVALS = ("pT_10_15", "pT_15_20", "pT_25_30")
PHOTON_CLASSES = ("direct", "fragmentation")
SPECTRUM_CLASSES = ("total", "direct", "fragmentation")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("roots", nargs="+", help="ROOT file(s) to validate")
    parser.add_argument("--json-out", type=Path)
    return parser.parse_args()


def get_object(root_file: Any, name: str) -> Any:
    obj = root_file.Get(f"SIM/{name}")
    if obj:
        return obj
    obj = root_file.Get(name)
    if obj:
        return obj
    return None


def expected_names() -> list[str]:
    names: list[str] = []
    for centrality in CENTRALITIES:
        for photon_class in PHOTON_CLASSES:
            for pt_interval in PT_INTERVALS:
                names.append(f"h_auauTruthIso_{photon_class}_{pt_interval}_{centrality}")
        for photon_class in SPECTRUM_CLASSES:
            names.append(f"h_auauTruthPt_{photon_class}_iso4_{centrality}")
    names.append("h_auauTruthIsoDiag_audit")
    return names


def validate_file(path: str) -> dict[str, Any]:
    root_path = Path(path)
    result: dict[str, Any] = {
        "path": str(root_path),
        "exists": root_path.is_file(),
        "bytes": root_path.stat().st_size if root_path.is_file() else 0,
        "errors": [],
        "warnings": [],
    }
    if not root_path.is_file():
        result["errors"].append("file missing")
        return result

    root_file = ROOT.TFile.Open(str(root_path), "READ")
    if not root_file or root_file.IsZombie():
        result["errors"].append("ROOT file is zombie or failed to open")
        return result
    result["recovered"] = bool(root_file.TestBit(ROOT.TFile.kRecovered))
    if result["recovered"]:
        result["errors"].append("ROOT file has the recovered bit set")

    missing = [name for name in expected_names() if not get_object(root_file, name)]
    result["expected_object_count"] = len(expected_names())
    result["missing_objects"] = missing
    if missing:
        result["errors"].append(f"missing {len(missing)} expected diagnostic objects")

    audit = get_object(root_file, "h_auauTruthIsoDiag_audit")
    if audit:
        audit_values = [float(audit.GetBinContent(index)) for index in range(1, 17)]
        result["audit_values"] = audit_values
        result["audit"] = {
            "accepted_events": audit_values[0],
            "direct_class": audit_values[2],
            "fragmentation_class": audit_values[3],
            "invalid_truth_match": audit_values[5],
            "vertex_rejection": audit_values[7],
            "invalid_centrality": audit_values[8],
            "events_by_centrality": audit_values[13:16],
        }
        if audit_values[2] <= 0.0 or audit_values[3] <= 0.0:
            result["errors"].append("direct or fragmentation audit population is zero")
        if not math.isclose(audit_values[0], sum(audit_values[13:16]), rel_tol=1.0e-9, abs_tol=1.0e-9):
            result["errors"].append("accepted event count does not equal centrality event sum")

    closure_failures: list[dict[str, Any]] = []
    for centrality in CENTRALITIES:
        total = get_object(root_file, f"h_auauTruthPt_total_iso4_{centrality}")
        direct = get_object(root_file, f"h_auauTruthPt_direct_iso4_{centrality}")
        fragmentation = get_object(root_file, f"h_auauTruthPt_fragmentation_iso4_{centrality}")
        if not total or not direct or not fragmentation:
            continue
        if total.GetNbinsX() != 25 or total.GetXaxis().GetXmin() != 10.0 or total.GetXaxis().GetXmax() != 35.0:
            result["errors"].append(f"wrong truth-pT binning for {centrality}")
        for bin_index in range(0, total.GetNbinsX() + 2):
            difference = float(
                total.GetBinContent(bin_index)
                - direct.GetBinContent(bin_index)
                - fragmentation.GetBinContent(bin_index)
            )
            scale = max(abs(float(total.GetBinContent(bin_index))), 1.0)
            if abs(difference) > 1.0e-9 * scale:
                closure_failures.append(
                    {"centrality": centrality, "bin": bin_index, "difference": difference}
                )
    result["truth_pt_closure_failures"] = closure_failures
    if closure_failures:
        result["errors"].append("truth-pT total != direct + fragmentation")

    fraction_checks = 0
    empty_fraction_histograms: list[str] = []
    for centrality in CENTRALITIES:
        for photon_class in PHOTON_CLASSES:
            for pt_interval in PT_INTERVALS:
                name = f"h_auauTruthIso_{photon_class}_{pt_interval}_{centrality}"
                hist = get_object(root_file, name)
                if not hist:
                    continue
                denominator = float(hist.Integral(0, hist.GetNbinsX() + 1))
                if denominator <= 0.0:
                    empty_fraction_histograms.append(name)
                    continue
                previous = -1.0
                for cutoff in range(1, 21):
                    upper_bin = hist.GetXaxis().FindFixBin(float(cutoff) - 1.0e-9)
                    numerator = float(hist.Integral(0, upper_bin))
                    fraction = numerator / denominator
                    if not math.isfinite(fraction) or fraction < -1.0e-10 or fraction > 1.0 + 1.0e-10:
                        result["errors"].append(f"invalid cumulative fraction for {name} at {cutoff} GeV")
                    if fraction + 1.0e-10 < previous:
                        result["errors"].append(f"non-monotonic cumulative fraction for {name}")
                    previous = fraction
                    fraction_checks += 1
    result["fraction_cutoff_checks"] = fraction_checks
    result["empty_fraction_histograms"] = empty_fraction_histograms

    legacy_truth_iso = get_object(root_file, "h_EisoTruth")
    legacy_truth_decision = get_object(root_file, "h_EisoTruthDecision")
    result["legacy_truth_objects"] = {
        "h_EisoTruth": bool(legacy_truth_iso),
        "h_EisoTruthDecision": bool(legacy_truth_decision),
    }
    if not legacy_truth_iso or not legacy_truth_decision:
        result["errors"].append("legacy truth-isolation QA objects are missing")

    result["ok"] = not result["errors"]
    root_file.Close()
    return result


def main() -> int:
    args = parse_args()
    results = [validate_file(path) for path in args.roots]
    payload = {
        "schema_version": 1,
        "ok": all(result.get("ok", False) for result in results),
        "files": results,
    }
    text = json.dumps(payload, indent=2, sort_keys=True)
    print(text)
    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(text + "\n", encoding="utf-8")
    return 0 if payload["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
