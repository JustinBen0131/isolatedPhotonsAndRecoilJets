#!/usr/bin/env python3
"""Classify scaled-trigger run-by-run QA shapes from per-run ROOT histograms.

This script is intended to run in the SDCC analysis environment where PyROOT is
available. It reads the already-produced per-run ROOT files and prints a CSV to
stdout. The rules are deliberately simple and ordered so each run lands in one
maintainable, slide-facing category.
"""

from __future__ import annotations

import csv
import glob
import math
import os
import re
import sys

import ROOT  # type: ignore


DEFAULT_INPUT_DIR = (
    "/sphenix/u/patsfan753/scratch/thesisAnalysis/output/auau/perRun/"
    "jetMinPt5_7pi_8_vz60_isoR40_isSliding_baseVariant_preselectionReference_"
    "tightReference_nonTightReference_scaledTriggerStudy"
)

HISTS = {
    "mbd": (
        "MBD_NS_geq_2_vtx_lt_150/"
        "h_maxEnergyClus_NewTriggerFilling_perRunCorrected_MBD_NS_geq_2_vtx_lt_150"
    ),
    "p10": "Photon_10/h_maxEnergyClus_NewTriggerFilling_perRunCorrected_Photon_10",
    "p12": "Photon_12/h_maxEnergyClus_NewTriggerFilling_perRunCorrected_Photon_12",
}

RANGES = {
    "total": (1.0, 20.0),
    "low": (1.0, 3.0),
    "pre": (2.0, 4.0),
    "turn3": (3.0, 5.0),
    "turn4": (4.0, 6.0),
    "mid": (6.0, 9.0),
    "plateau": (9.0, 13.0),
    "tail": (15.0, 20.0),
}


def integral(hist, lo: float, hi: float) -> float:
    total = 0.0
    xaxis = hist.GetXaxis()
    for ibin in range(1, hist.GetNbinsX() + 1):
        center = xaxis.GetBinCenter(ibin)
        if lo <= center < hi:
            total += float(hist.GetBinContent(ibin))
    return total


def ratio(num: float, den: float) -> float:
    return num / den if den > 0 else math.nan


def classify(v: dict[str, float]) -> tuple[str, str]:
    """Return an ordered slide-facing category and a compact reason string."""
    near_tail = (
        v["mbd_tail"] >= 100.0
        and 0.97 <= v["r10_tail"] <= 1.08
        and 0.97 <= v["r12_tail"] <= 1.08
    )
    clean_shape = (
        v["r10_low"] < 0.05
        and v["r12_low"] < 0.02
        and v["r10_pre"] < 0.10
        and v["r12_pre"] < 0.04
        and v["r10_turn3"] < 0.25
        and v["r12_turn3"] < 0.10
        and v["r10_turn4"] < 0.50
        and v["r12_turn4"] < 0.22
    )

    max_tail_ratio = max(v["r10_tail"], v["r12_tail"])
    min_tail_ratio = min(v["r10_tail"], v["r12_tail"])
    max_total_ratio = max(v["r10_total"], v["r12_total"])

    # These are visually empty or triggerless: the photon-trigger spectra are
    # absent or nearly absent relative to the MBD reference.
    if v["mbd_tail"] >= 100.0 and max_tail_ratio < 0.05:
        return "missing_photon_trigger_data", "tail ratios <0.05"
    if v["mbd_total"] >= 1.0e6 and max_total_ratio < 0.01:
        return "missing_photon_trigger_data", "total trigger/MBD <0.01"

    # These are not missing data; they are effectively identical to MBD from
    # threshold, which makes the turn-on overlay flat near unity and unphysical
    # for a photon trigger.
    if near_tail and v["r10_low"] > 0.90 and v["r12_low"] > 0.90:
        return "flat_unity_from_threshold", "low ratios near 1"

    # Tail unity alone is not enough: 71260-like runs start too high before the
    # physical turn-on, even though the high-E tail agrees with MBD.
    if near_tail and not clean_shape:
        return "early_turn_on_outlier", "tail unity, low-E already on"

    if near_tail and clean_shape:
        return "clean_full_turn_on", "tail unity and low-E near zero"

    if v["mbd_tail"] < 100.0:
        return "low_stat_tail", "MBD tail <100"

    if min_tail_ratio < 0.80:
        return "tail_inefficient", "tail ratio <0.80"

    if max_tail_ratio > 1.25:
        return "tail_overscaled", "tail ratio >1.25"

    if min_tail_ratio < 0.92 or max_tail_ratio > 1.15:
        return "tail_intermediate_deviation", "tail outside [0.92,1.15]"

    return "mild_tail_deviation", "tail near but not clean-unity"


def main() -> int:
    ROOT.gROOT.SetBatch(True)
    input_dir = sys.argv[1] if len(sys.argv) > 1 else os.environ.get("SCALED_TRIGGER_INPUT_DIR", DEFAULT_INPUT_DIR)
    paths = sorted(glob.glob(os.path.join(input_dir, "chunkMerge_run_*.root")))

    writer = csv.writer(sys.stdout)
    header = [
        "run",
        "category",
        "reason",
        "mbd_total",
        "p10_total",
        "p12_total",
        "r10_total",
        "r12_total",
        "mbd_low",
        "r10_low",
        "r12_low",
        "mbd_pre",
        "r10_pre",
        "r12_pre",
        "r10_turn3",
        "r12_turn3",
        "r10_turn4",
        "r12_turn4",
        "r10_mid",
        "r12_mid",
        "mbd_tail",
        "p10_tail",
        "p12_tail",
        "r10_tail",
        "r12_tail",
        "png",
        "root_path",
    ]
    writer.writerow(header)

    for path in paths:
        match = re.search(r"run_(\d+)\.root$", path)
        if not match:
            continue
        run = int(match.group(1))
        root_file = ROOT.TFile.Open(path)
        if not root_file or root_file.IsZombie():
            continue
        hists = {name: root_file.Get(hist_path) for name, hist_path in HISTS.items()}
        if any(not hist for hist in hists.values()):
            root_file.Close()
            continue

        vals: dict[str, float] = {}
        for label, (lo, hi) in RANGES.items():
            mbd = integral(hists["mbd"], lo, hi)
            p10 = integral(hists["p10"], lo, hi)
            p12 = integral(hists["p12"], lo, hi)
            vals[f"mbd_{label}"] = mbd
            vals[f"p10_{label}"] = p10
            vals[f"p12_{label}"] = p12
            vals[f"r10_{label}"] = ratio(p10, mbd)
            vals[f"r12_{label}"] = ratio(p12, mbd)

        category, reason = classify(vals)
        writer.writerow(
            [
                run,
                category,
                reason,
                vals["mbd_total"],
                vals["p10_total"],
                vals["p12_total"],
                vals["r10_total"],
                vals["r12_total"],
                vals["mbd_low"],
                vals["r10_low"],
                vals["r12_low"],
                vals["mbd_pre"],
                vals["r10_pre"],
                vals["r12_pre"],
                vals["r10_turn3"],
                vals["r12_turn3"],
                vals["r10_turn4"],
                vals["r12_turn4"],
                vals["r10_mid"],
                vals["r12_mid"],
                vals["mbd_tail"],
                vals["p10_tail"],
                vals["p12_tail"],
                vals["r10_tail"],
                vals["r12_tail"],
                f"png_all/run_{run:08d}_scaledTriggerQA_1x2.png",
                path,
            ]
        )
        root_file.Close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
