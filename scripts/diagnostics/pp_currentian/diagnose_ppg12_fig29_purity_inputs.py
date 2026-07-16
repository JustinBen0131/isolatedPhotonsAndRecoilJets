#!/usr/bin/env python3
"""Diagnose whether current pp RecoilJets outputs can reproduce PPG12 Fig. 29.

This is intentionally an audit, not a plotting helper.  PPG12 Fig. 29 is made
from the efficiencytool/CalculatePhotonYield.C contract, while the RecoilJets
table-QA output contains a separate ABCD diagnostic family.  This script records
the exact binning/key mismatch so slide generators do not silently present the
diagnostic as an apples-to-apples IAN reproduction.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import ROOT


ROOT.gROOT.SetBatch(True)

DEFAULT_BASE = Path(
    "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/ppg12TableQA/"
    "THE42_ppg12_tableqa_v1_basev3e_20260611"
)
DEFAULT_DATA_ROOT = (
    DEFAULT_BASE
    / "merged_roots/RecoilJets_pp_ALL_preselectionNewPPG12_tightNewPPG12_nonTightNewPPG12.root"
)
DEFAULT_SIGNAL_ROOT = DEFAULT_BASE / "merged_roots/RecoilJets_photonjet5plus10plus20_MERGED.root"
DEFAULT_OUTDIR = DEFAULT_BASE / "purity_current_pp"

TRIGGER_DIR = "PPG12_scaledtrigger30"
SIM_DIR = "SIM"
ISO_TOKEN = "isoR40_fixedIso2GeV"

PPG12_IAN_BINS = [(10, 12), (12, 14), (14, 16), (16, 18), (18, 20), (20, 22), (22, 24), (24, 26), (26, 28), (28, 32), (32, 36)]
RECOILJETS_JES3_BINS = [(5, 8), (8, 10), (10, 12), (12, 14), (14, 16), (16, 18), (18, 20), (20, 22), (22, 24), (24, 26), (26, 35)]


def open_file(path: Path) -> ROOT.TFile:
    f = ROOT.TFile.Open(str(path))
    if not f or f.IsZombie():
        raise RuntimeError(f"Could not open ROOT file: {path}")
    return f


def dir_keys(directory: ROOT.TDirectory) -> set[str]:
    return {key.GetName() for key in directory.GetListOfKeys()}


def bin1(directory: ROOT.TDirectory, name: str) -> float | None:
    h = directory.Get(name)
    if not h:
        return None
    return float(h.GetBinContent(1))


def raw_purity(a: float | None, b: float | None, c: float | None, d: float | None) -> float | None:
    if a is None or b is None or c is None or d is None or a <= 0.0 or d <= 0.0:
        return None
    return max(a - b * c / d, 0.0) / a


def continuous_counts(directory: ROOT.TDirectory, lo: float, hi: float) -> dict[str, float | None]:
    out: dict[str, float | None] = {}
    for region in "ABCD":
        h = directory.Get(f"h_pTgamma_ABCD_{region}_{ISO_TOKEN}")
        if not h:
            out[region] = None
            continue
        axis = h.GetXaxis()
        ilo = axis.FindBin(lo + 1.0e-6)
        ihi = axis.FindBin(hi - 1.0e-6)
        out[region] = float(h.Integral(ilo, ihi))
    return out


def hist_edges(directory: ROOT.TDirectory, hist_name: str) -> list[float] | None:
    h = directory.Get(hist_name)
    if not h:
        return None
    axis = h.GetXaxis()
    return [float(axis.GetBinLowEdge(i)) for i in range(1, h.GetNbinsX() + 1)] + [float(axis.GetBinUpEdge(h.GetNbinsX()))]


def finite_or_none(value: float | None) -> float | None:
    if value is None or not math.isfinite(value):
        return None
    return value


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", type=Path, default=DEFAULT_DATA_ROOT)
    ap.add_argument("--signal-root", type=Path, default=DEFAULT_SIGNAL_ROOT)
    ap.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    f_data = open_file(args.data_root)
    f_sig = open_file(args.signal_root)
    d_data = f_data.Get(TRIGGER_DIR)
    d_sig = f_sig.Get(SIM_DIR)
    if not d_data:
        raise RuntimeError(f"Missing trigger directory {TRIGGER_DIR} in {args.data_root}")
    if not d_sig:
        raise RuntimeError(f"Missing SIM directory {SIM_DIR} in {args.signal_root}")

    keys = dir_keys(d_data)
    one_bin_prefixes = [
        "h_isIsolated_isTight",
        "h_notIsolated_isTight",
        "h_isIsolated_notTight",
        "h_notIsolated_notTight",
    ]

    ian_rows: list[dict[str, Any]] = []
    for lo, hi in PPG12_IAN_BINS:
        suffix = f"_{ISO_TOKEN}_pT_{lo}_{hi}"
        one_bin_counts = {prefix: bin1(d_data, prefix + suffix) for prefix in one_bin_prefixes}
        continuous = continuous_counts(d_data, lo, hi)
        ian_rows.append(
            {
                "bin": [lo, hi],
                "one_bin_keys_present": all((prefix + suffix) in keys for prefix in one_bin_prefixes),
                "one_bin_counts": one_bin_counts,
                "continuous_counts": continuous,
                "continuous_raw_purity": finite_or_none(
                    raw_purity(continuous["A"], continuous["B"], continuous["C"], continuous["D"])
                ),
            }
        )

    recoiljets_rows: list[dict[str, Any]] = []
    for lo, hi in RECOILJETS_JES3_BINS:
        suffix = f"_{ISO_TOKEN}_pT_{lo}_{hi}"
        counts = {prefix: bin1(d_data, prefix + suffix) for prefix in one_bin_prefixes}
        recoiljets_rows.append(
            {
                "bin": [lo, hi],
                "one_bin_keys_present": all((prefix + suffix) in keys for prefix in one_bin_prefixes),
                "A": counts["h_isIsolated_isTight"],
                "B": counts["h_notIsolated_isTight"],
                "C": counts["h_isIsolated_notTight"],
                "D": counts["h_notIsolated_notTight"],
                "raw_purity": finite_or_none(raw_purity(counts["h_isIsolated_isTight"], counts["h_notIsolated_isTight"], counts["h_isIsolated_notTight"], counts["h_notIsolated_notTight"])),
            }
        )

    report: dict[str, Any] = {
        "data_root": str(args.data_root),
        "signal_root": str(args.signal_root),
        "trigger_dir": TRIGGER_DIR,
        "sim_dir": SIM_DIR,
        "diagnosis": "Current RecoilJets table-QA output is not an apples-to-apples PPG12 Fig. 29 purity input.",
        "primary_reasons": [
            "PPG12 Fig. 29 uses efficiencytool/CalculatePhotonYield.C graph outputs gpurity, gpurity_leak, and grFineConf_leak.",
            "The current ROOT has RecoilJets ABCD diagnostic histograms h_isIsolated_* and h_pTgamma_ABCD_* instead.",
            "The RecoilJets photon pT binning is 5,8,10,12,14,16,18,20,22,24,26,35, while the current IAN config uses 10,12,14,16,18,20,22,24,26,28,32,36.",
            "The current output has no one-bin ABCD keys for 26-28, 28-32, or 32-36; it has one coarse 26-35 bin.",
        ],
        "ppg12_ian_bins": PPG12_IAN_BINS,
        "recoiljets_jes3_bins": RECOILJETS_JES3_BINS,
        "continuous_abcd_edges": hist_edges(d_data, f"h_pTgamma_ABCD_A_{ISO_TOKEN}"),
        "ian_bin_audit": ian_rows,
        "recoiljets_bin_audit": recoiljets_rows,
    }

    json_path = args.outdir / "ppg12_fig29_purity_input_mismatch_diagnostic.json"
    json_path.write_text(json.dumps(report, indent=2))

    md_path = args.outdir / "ppg12_fig29_purity_input_mismatch_diagnostic.md"
    lines = [
        "# PPG12 Fig. 29 purity input mismatch diagnostic",
        "",
        f"Data ROOT: `{args.data_root}`",
        f"Signal ROOT: `{args.signal_root}`",
        "",
        "## Diagnosis",
        "",
        report["diagnosis"],
        "",
        "## Evidence",
        "",
    ]
    for reason in report["primary_reasons"]:
        lines.append(f"- {reason}")
    lines.extend(["", "## RecoilJets ABCD bins", "", "| bin | A | B | C | D | raw purity |", "| --- | ---: | ---: | ---: | ---: | ---: |"])
    for row in recoiljets_rows:
        lines.append(
            f"| {row['bin'][0]}-{row['bin'][1]} | {row['A']} | {row['B']} | {row['C']} | {row['D']} | {row['raw_purity']} |"
        )
    lines.extend(["", "## PPG12 IAN bin key presence", "", "| bin | one-bin ABCD keys present | continuous A/B/C/D |", "| --- | --- | --- |"])
    for row in ian_rows:
        cc = row["continuous_counts"]
        lines.append(
            f"| {row['bin'][0]}-{row['bin'][1]} | {row['one_bin_keys_present']} | "
            f"{cc['A']}, {cc['B']}, {cc['C']}, {cc['D']} |"
        )
    lines.append("")
    md_path.write_text("\n".join(lines))

    print(f"wrote {json_path}")
    print(f"wrote {md_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
