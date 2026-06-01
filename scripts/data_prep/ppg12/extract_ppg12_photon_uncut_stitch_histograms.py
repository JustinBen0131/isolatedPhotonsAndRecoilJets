#!/usr/bin/env python3
"""Dump PPG12 photon_max_pT_uncut histograms as CSV.

This reads the exact histogram source used by PPG12's `plot_combine_uncut.C`
for the analysis-note photon stitching figure. It prints CSV to stdout and
does not write remote files.
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

SOURCE = "/sphenix/user/shuhangli/ppg12/plotting/photon_max_pT_uncut.root"
SPECS = [
    ("photon5", "h_max_photon_pT_photon5", "h_max_photon_pT_photon5_sumw2", 0.0, 14.0, "#E7298A"),
    ("photon10", "h_max_photon_pT_photon10", "h_max_photon_pT_photon10_sumw2", 14.0, 22.0, "#33A02C"),
    ("photon20", "h_max_photon_pT_photon20", "h_max_photon_pT_photon20_sumw2", 22.0, 200.0, "#1F78B4"),
]


def main() -> None:
    root_file = ROOT.TFile.Open(SOURCE, "READ")
    if not root_file or root_file.IsZombie():
        raise SystemExit(f"missing ROOT file: {SOURCE}")

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
            "sumw2_hist_name",
            "bin_width_GeV",
            "raw_integral",
        ]
    )
    for sample, hist_name, sumw2_name, lo, hi, color in SPECS:
        hist = root_file.Get(hist_name)
        sumw2 = root_file.Get(sumw2_name)
        if not hist or not sumw2:
            raise SystemExit(f"missing {hist_name} or {sumw2_name}: {SOURCE}")
        bin_width = hist.GetBinWidth(1)
        raw_integral = hist.Integral()
        for bin_idx in range(1, hist.GetNbinsX() + 1):
            bin_low = hist.GetBinLowEdge(bin_idx)
            bin_high = bin_low + bin_width
            center = hist.GetBinCenter(bin_idx)
            used = int(center >= lo and center < hi)
            err2 = sumw2.GetBinContent(bin_idx)
            writer.writerow(
                [
                    "photon",
                    sample,
                    f"{bin_low:.8g}",
                    f"{bin_high:.8g}",
                    f"{center:.8g}",
                    f"{hist.GetBinContent(bin_idx):.12g}",
                    f"{err2 ** 0.5 if err2 > 0 else 0.0:.12g}",
                    f"{lo:.8g}",
                    f"{hi:.8g}",
                    used,
                    color,
                    SOURCE,
                    hist_name,
                    sumw2_name,
                    f"{bin_width:.8g}",
                    f"{raw_integral:.12g}",
                ]
            )
    print("END_CSV")
    root_file.Close()


if __name__ == "__main__":
    main()
