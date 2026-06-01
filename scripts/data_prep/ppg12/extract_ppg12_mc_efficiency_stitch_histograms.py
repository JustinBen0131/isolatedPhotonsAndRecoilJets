#!/usr/bin/env python3
"""Dump PPG12 MC_efficiency stitch histograms as CSV.

This script is intended to be streamed to SDCC with:
  ssh ... "ssh ... 'python3 -'" < scripts/extract_ppg12_mc_efficiency_stitch_histograms.py

It reads Shuhang's existing MC_efficiency ROOT outputs and prints compact CSV
to stdout. It does not write remote files.
"""

from __future__ import annotations
# Keep purpose-folder helpers runnable when invoked directly.
import sys as _codex_sys
from pathlib import Path as _CodexPath
_CODEX_THIS_FILE = _CodexPath(__file__).resolve()
_CODEX_SCRIPTS_DIR = next((p for p in _CODEX_THIS_FILE.parents if p.name == "scripts"), _CODEX_THIS_FILE.parent)
_CODEX_SCRIPTS_DIR_STR = str(_CODEX_SCRIPTS_DIR)
if _CODEX_SCRIPTS_DIR_STR not in _codex_sys.path:
    _codex_sys.path.append(_CODEX_SCRIPTS_DIR_STR)
del _CODEX_THIS_FILE, _CODEX_SCRIPTS_DIR, _CODEX_SCRIPTS_DIR_STR

import csv
import sys

import ROOT


ROOT.gROOT.SetBatch(True)

BASE = "/sphenix/user/shuhangli/ppg12/efficiencytool/results"
SPECS = [
    ("photon", "photon5", "MC_efficiency_photon5_bdt_nom.root", "h_max_photon_pT", 0.0, 14.0, "#E7298A"),
    ("photon", "photon10", "MC_efficiency_photon10_bdt_nom.root", "h_max_photon_pT", 14.0, 22.0, "#33A02C"),
    ("photon", "photon20", "MC_efficiency_photon20_bdt_nom.root", "h_max_photon_pT", 22.0, 200.0, "#1F78B4"),
    ("jet", "jet8", "MC_efficiency_jet8_bdt_nom.root", "h_max_truth_jet_pT", 9.0, 14.0, "#E7298A"),
    ("jet", "jet12", "MC_efficiency_jet12_bdt_nom.root", "h_max_truth_jet_pT", 14.0, 21.0, "#33A02C"),
    ("jet", "jet20", "MC_efficiency_jet20_bdt_nom.root", "h_max_truth_jet_pT", 21.0, 32.0, "#1F78B4"),
    ("jet", "jet30", "MC_efficiency_jet30_bdt_nom.root", "h_max_truth_jet_pT", 32.0, 42.0, "#FF7F00"),
    ("jet", "jet40", "MC_efficiency_jet40_bdt_nom.root", "h_max_truth_jet_pT", 42.0, 200.0, "#E7298A"),
]


def main() -> None:
    print("BEGIN_CSV")
    writer = csv.writer(sys.stdout)
    writer.writerow(
        [
            "group",
            "sample",
            "bin_low",
            "bin_high",
            "bin_center",
            "value",
            "error",
            "stitch_window_low",
            "stitch_window_high",
            "used_in_stitch",
            "color",
            "source_root",
            "hist_name",
            "bin_width_GeV",
            "raw_integral",
        ]
    )
    for group, sample, fname, hist_name, lo, hi, color in SPECS:
        path = f"{BASE}/{fname}"
        root_file = ROOT.TFile.Open(path, "READ")
        if not root_file or root_file.IsZombie():
            raise SystemExit(f"missing ROOT file: {path}")
        hist = root_file.Get(hist_name)
        if not hist:
            raise SystemExit(f"missing histogram {hist_name}: {path}")
        rebinned = hist.Clone(f"{sample}_{hist_name}_rebin10")
        rebinned.SetDirectory(0)
        rebinned.Rebin(10)
        bin_width = rebinned.GetBinWidth(1)
        raw_integral = hist.Integral()
        for bin_idx in range(1, rebinned.GetNbinsX() + 1):
            bin_low = rebinned.GetBinLowEdge(bin_idx)
            bin_high = bin_low + bin_width
            center = rebinned.GetBinCenter(bin_idx)
            used = int(center >= lo and center < hi)
            writer.writerow(
                [
                    group,
                    sample,
                    f"{bin_low:.8g}",
                    f"{bin_high:.8g}",
                    f"{center:.8g}",
                    f"{rebinned.GetBinContent(bin_idx):.12g}",
                    f"{rebinned.GetBinError(bin_idx):.12g}",
                    f"{lo:.8g}",
                    f"{hi:.8g}",
                    used,
                    color,
                    path,
                    hist_name,
                    f"{bin_width:.8g}",
                    f"{raw_integral:.12g}",
                ]
            )
        root_file.Close()
    print("END_CSV")


if __name__ == "__main__":
    main()
