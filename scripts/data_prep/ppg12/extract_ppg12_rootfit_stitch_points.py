#!/usr/bin/env python3
"""Dump PPG12-style stitched spectra after ROOT TF1 fits.

This script is intended to be streamed to SDCC with:
  ssh ... "ssh ... 'python3 -'" < scripts/extract_ppg12_rootfit_stitch_points.py

It reads Shuhang's PPG12 ROOT histogram sources, applies the same ownership
windows used in the IAN stitching figures, fits the stitched histograms with
the PPG12 ROOT TF1 form, and prints compact CSV to stdout. It does not write
remote files.
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
import math
import sys

import ROOT


ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

PHOTON_SOURCE = "/sphenix/user/shuhangli/ppg12/plotting/photon_max_pT_uncut.root"
EFFICIENCY_DIR = "/sphenix/user/shuhangli/ppg12/efficiencytool/results"

PHOTON_SPECS = [
    ("photon5", "h_max_photon_pT_photon5", "h_max_photon_pT_photon5_sumw2", 0.0, 14.0, "#E7298A"),
    ("photon10", "h_max_photon_pT_photon10", "h_max_photon_pT_photon10_sumw2", 14.0, 22.0, "#33A02C"),
    ("photon20", "h_max_photon_pT_photon20", "h_max_photon_pT_photon20_sumw2", 22.0, 200.0, "#1F78B4"),
]

JET_SPECS = [
    ("jet8", "MC_efficiency_jet8_bdt_nom.root", "h_max_truth_jet_pT", 9.0, 14.0, "#E7298A"),
    ("jet12", "MC_efficiency_jet12_bdt_nom.root", "h_max_truth_jet_pT", 14.0, 21.0, "#33A02C"),
    ("jet20", "MC_efficiency_jet20_bdt_nom.root", "h_max_truth_jet_pT", 21.0, 32.0, "#1F78B4"),
    ("jet30", "MC_efficiency_jet30_bdt_nom.root", "h_max_truth_jet_pT", 32.0, 42.0, "#FF7F00"),
    ("jet40", "MC_efficiency_jet40_bdt_nom.root", "h_max_truth_jet_pT", 42.0, 200.0, "#E7298A"),
]


def set_sumw2_errors(hist, sumw2_hist) -> None:
    hist.Sumw2()
    for idx in range(1, hist.GetNbinsX() + 1):
        err2 = sumw2_hist.GetBinContent(idx)
        hist.SetBinError(idx, math.sqrt(err2) if err2 > 0 else 0.0)


def zero_outside(hist, lo: float, hi: float) -> None:
    for idx in range(1, hist.GetNbinsX() + 1):
        center = hist.GetBinCenter(idx)
        if center < lo or center >= hi:
            hist.SetBinContent(idx, 0.0)
            hist.SetBinError(idx, 0.0)


def fetch_photon_hists() -> tuple[list[tuple[str, object, float, float, str]], object]:
    root_file = ROOT.TFile.Open(PHOTON_SOURCE, "READ")
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"cannot open {PHOTON_SOURCE}")

    hists = []
    stitched = None
    for sample, hist_name, sumw2_name, lo, hi, color in PHOTON_SPECS:
        hist = root_file.Get(hist_name)
        sumw2 = root_file.Get(sumw2_name)
        if not hist or not sumw2:
            raise RuntimeError(f"missing {hist_name} or {sumw2_name} in {PHOTON_SOURCE}")
        clone = hist.Clone(f"{sample}_windowed")
        clone.SetDirectory(0)
        set_sumw2_errors(clone, sumw2)
        zero_outside(clone, lo, hi)
        hists.append((sample, clone, lo, hi, color))
        if stitched is None:
            stitched = clone.Clone("photon_stitched")
            stitched.SetDirectory(0)
        else:
            stitched.Add(clone)
    root_file.Close()
    return hists, stitched


def fetch_jet_hists() -> tuple[list[tuple[str, object, float, float, str]], object]:
    hists = []
    stitched = None
    for sample, filename, hist_name, lo, hi, color in JET_SPECS:
        path = f"{EFFICIENCY_DIR}/{filename}"
        root_file = ROOT.TFile.Open(path, "READ")
        if not root_file or root_file.IsZombie():
            raise RuntimeError(f"cannot open {path}")
        hist = root_file.Get(hist_name)
        if not hist:
            raise RuntimeError(f"missing {hist_name} in {path}")
        clone = hist.Clone(f"{sample}_windowed")
        clone.SetDirectory(0)
        root_file.Close()
        clone.Rebin(10)
        zero_outside(clone, lo, hi)
        hists.append((sample, clone, lo, hi, color))
        if stitched is None:
            stitched = clone.Clone("jet_stitched")
            stitched.SetDirectory(0)
        else:
            stitched.Add(clone)
    return hists, stitched


def fit_hist(hist, name: str, xmin: float, xmax: float, fit_xmax: float):
    func = ROOT.TF1(
        f"f1_{name}",
        "[0]*pow([1]/x,[2]+[3]*log(x/[1])+ [4]*x)",
        xmin,
        xmax,
    )
    func.SetParameters(2.09375e9, 1.0, 1.0, 2.0, 0.01)
    hist.Fit(func, "REMNQ", "", xmin, fit_xmax)
    hist.Fit(func, "REMNQ", "", xmin, fit_xmax)
    return func


def write_group(writer, group: str, hists, stitched, ratio_label: str, fit_xmax: float) -> None:
    fit = fit_hist(stitched, group, 10.0, 50.0, fit_xmax)
    params = [fit.GetParameter(i) for i in range(5)]
    param_str = ";".join(f"{p:.12g}" for p in params)

    for sample, hist, lo, hi, color in hists:
        for idx in range(1, hist.GetNbinsX() + 1):
            x = hist.GetBinCenter(idx)
            y = hist.GetBinContent(idx)
            if x < 0 or x > 100:
                continue
            writer.writerow(
                [
                    group,
                    sample,
                    f"{hist.GetBinLowEdge(idx):.8g}",
                    f"{hist.GetBinLowEdge(idx) + hist.GetBinWidth(idx):.8g}",
                    f"{x:.8g}",
                    f"{y:.12g}",
                    f"{hist.GetBinError(idx):.12g}",
                    f"{lo:.8g}",
                    f"{hi:.8g}",
                    int(y > 0),
                    color,
                    "",
                    "",
                    "",
                    "",
                    ratio_label,
                ]
            )

    for idx in range(1, stitched.GetNbinsX() + 1):
        x = stitched.GetBinCenter(idx)
        y = stitched.GetBinContent(idx)
        ey = stitched.GetBinError(idx)
        fy = fit.Eval(x) if x > 0 else 0.0
        ratio = y / fy if fy > 0 and y > 0 else 0.0
        ratio_err = ey / fy if fy > 0 and y > 0 else 0.0
        writer.writerow(
            [
                group,
                "stitched",
                f"{stitched.GetBinLowEdge(idx):.8g}",
                f"{stitched.GetBinLowEdge(idx) + stitched.GetBinWidth(idx):.8g}",
                f"{x:.8g}",
                f"{y:.12g}",
                f"{ey:.12g}",
                "",
                "",
                int(y > 0),
                "#000000",
                f"{fy:.12g}",
                f"{ratio:.12g}",
                f"{ratio_err:.12g}",
                param_str,
                ratio_label,
            ]
        )


def main() -> None:
    photon_hists, photon_stitched = fetch_photon_hists()
    jet_hists, jet_stitched = fetch_jet_hists()

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
            "root_fit_value",
            "root_fit_ratio",
            "root_fit_ratio_error",
            "root_fit_params",
            "ratio_label",
        ]
    )
    write_group(writer, "photon", photon_hists, photon_stitched, "Data / Fit", 36.0)
    write_group(writer, "jet", jet_hists, jet_stitched, "MC / Fit", 50.0)
    print("END_CSV")


if __name__ == "__main__":
    main()
