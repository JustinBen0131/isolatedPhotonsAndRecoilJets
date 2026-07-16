#!/usr/bin/env python3
"""Reproduce PPG12 IAN Fig. 3 from SDCC-extracted source points.

Source convention:
  ppg12codeGit/efficiencytool/FindTruthETCut.C
  /sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_noiso.root
  h_direct_pT_truth_isoET_0, h_frag_pT_truth_isoET_0

The source CSV is expected to contain the PPG12 ROOT projection after summing
normal E_T^iso bins with upper edge <= 4 GeV and applying the same Rebin(10)
contract as the PPG12 macro.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
import ROOT


REPO = Path(__file__).resolve().parents[3]
DEFAULT_OUTDIR = REPO / "dataOutput/ppg12Parity/the76_ppg12_fig3_sdcc_reference"
DEFAULT_POINTS_CSV = DEFAULT_OUTDIR / "ppg12_fig3_sdcc_root_points.csv"

SERIES = {
    "total": {
        "label": "Total",
        "color": "#d4148e",
        "marker": "o",
        "zorder": 4,
    },
    "direct": {
        "label": "Direct",
        "color": "#2f680d",
        "marker": "s",
        "zorder": 3,
    },
    "frag": {
        "label": "Fragmentation",
        "color": "#1f78ff",
        "marker": "^",
        "zorder": 2,
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--points-csv", type=Path, default=DEFAULT_POINTS_CSV)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--tag", default="ppg12_fig3_sdcc_root_reproduction")
    return parser.parse_args()


def read_points(path: Path) -> dict[str, dict[str, np.ndarray]]:
    rows_by_series: dict[str, list[dict[str, str]]] = {key: [] for key in SERIES}
    with path.open() as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if row.get("series") in rows_by_series:
                rows_by_series[row["series"]].append(row)
    out: dict[str, dict[str, np.ndarray]] = {}
    for key, rows in rows_by_series.items():
        if not rows:
            raise RuntimeError(f"{path} has no rows for series={key}")
        rows = sorted(rows, key=lambda row: float(row["xcenter_gev"]))
        out[key] = {
            "x": np.array([float(row["xcenter_gev"]) for row in rows], dtype=float),
            "y": np.array([float(row["count"]) for row in rows], dtype=float),
            "err": np.array([float(row.get("stat_err") or 0.0) for row in rows], dtype=float),
        }
    return out


def draw(points: dict[str, dict[str, np.ndarray]], out_png: Path) -> None:
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)

    root_colors = {
        "total": ROOT.kPink + 8,
        "direct": ROOT.kSpring - 7,
        "frag": ROOT.kAzure - 3,
    }
    root_markers = {"total": 20, "direct": 21, "frag": 22}
    root_labels = {"total": "Total", "direct": "Direct", "frag": "Fragmentation"}

    hists = {}
    for key in SERIES:
        x = points[key]["x"]
        y = points[key]["y"]
        err = points[key]["err"]
        if len(x) < 2:
            raise RuntimeError(f"Need at least two x points for {key}")
        width = float(np.median(np.diff(x)))
        low = float(x[0] - 0.5 * width)
        high = float(x[-1] + 0.5 * width)
        hist = ROOT.TH1F(f"h_{key}", f"h_{key}", len(x), low, high)
        hist.SetDirectory(0)
        hist.SetTitle("")
        for i, (yy, ee) in enumerate(zip(y, err), start=1):
            hist.SetBinContent(i, float(yy))
            hist.SetBinError(i, float(ee))
        hist.SetMarkerStyle(root_markers[key])
        hist.SetMarkerColor(root_colors[key])
        hist.SetLineColor(root_colors[key])
        hist.SetMarkerSize(1.0)
        hists[key] = hist

    can = ROOT.TCanvas("can_ppg12_fig3", "", 800, 889)
    can.Divide(1, 2)

    pad1 = can.cd(1)
    pad1.SetPad(0, 0.4, 1, 1)
    pad1.SetTopMargin(0.12)
    pad1.SetLeftMargin(0.13)
    pad1.SetBottomMargin(0.035)
    pad1.SetRightMargin(0.08)
    pad1.SetLogy()
    pad1.SetTicks(1, 1)

    hists["total"].GetXaxis().SetRangeUser(10, 35)
    hists["total"].GetYaxis().SetRangeUser(7.5e3, 1.35e8)
    hists["total"].GetXaxis().SetTitle("")
    hists["total"].GetYaxis().SetTitle("Counts")
    hists["total"].GetXaxis().SetTitleOffset(0.98)
    hists["total"].GetYaxis().SetTitleOffset(1.15)
    hists["total"].GetXaxis().SetLabelSize(0)
    hists["total"].GetYaxis().SetLabelSize(0.045)
    hists["total"].GetXaxis().SetTitleSize(0.045)
    hists["total"].GetYaxis().SetTitleSize(0.045)
    hists["total"].GetXaxis().SetLabelOffset(2)
    hists["total"].GetXaxis().CenterTitle()
    hists["total"].GetYaxis().CenterTitle()
    hists["total"].Draw("P")
    hists["direct"].Draw("P SAME")
    hists["frag"].Draw("P SAME")

    def ndc_text(x: float, y: float, text: str, size: float = 0.04, font: int = 42) -> ROOT.TLatex:
        latex = ROOT.TLatex()
        latex.SetNDC(True)
        latex.SetTextFont(font)
        latex.SetTextSize(size)
        latex.DrawLatex(x, y, text)
        return latex

    held_text = [
        ndc_text(0.065, 0.97, "#bf{#it{sPHENIX}} Simulation", 0.04),
        ndc_text(0.065, 0.92, "Photon Jet Samples", 0.04),
        ndc_text(0.16, 0.34, "Pythia, #sqrt{s}=200 GeV", 0.035),
        ndc_text(0.16, 0.29, "|#eta^{#gamma}| < 0.7", 0.035),
        ndc_text(0.16, 0.24, "R = 0.3, E_{T}^{iso} < 4 GeV", 0.035),
    ]

    legend = ROOT.TLegend(0.64, 0.68, 0.90, 0.83)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.035)
    for key in ["total", "direct", "frag"]:
        legend.AddEntry(hists[key], root_labels[key], "lep")
    legend.Draw()

    pad2 = can.cd(2)
    pad2.SetPad(0, 0, 1, 0.4)
    pad2.SetTopMargin(0.02)
    pad2.SetLeftMargin(0.13)
    pad2.SetBottomMargin(0.25)
    pad2.SetRightMargin(0.08)
    pad2.SetTicks(1, 1)

    ratio = hists["direct"].Clone("h_ratio_direct_over_total")
    ratio.SetDirectory(0)
    ratio.SetTitle("")
    ratio.Divide(hists["total"])
    ratio.SetMarkerStyle(root_markers["direct"])
    ratio.SetMarkerColor(root_colors["direct"])
    ratio.SetLineColor(root_colors["direct"])
    ratio.GetXaxis().SetTitle("p_{T} [GeV]")
    ratio.GetYaxis().SetTitle("Direct/Total")
    ratio.GetXaxis().SetRangeUser(10, 35)
    ratio.GetYaxis().SetRangeUser(0.5, 1.0)
    ratio.GetXaxis().CenterTitle()
    ratio.GetYaxis().CenterTitle()
    ratio.GetXaxis().SetTitleOffset(0.98)
    ratio.GetYaxis().SetTitleOffset(hists["total"].GetYaxis().GetTitleOffset() * 4 / 6.0)
    ratio.GetYaxis().SetLabelOffset(hists["total"].GetYaxis().GetLabelOffset() * 4 / 6.0)
    ratio.GetXaxis().SetLabelSize(0.0675)
    ratio.GetYaxis().SetLabelSize(0.0675)
    ratio.GetXaxis().SetTitleSize(0.0675)
    ratio.GetYaxis().SetTitleSize(0.0675)
    ratio.Draw("HIST")

    line = ROOT.TLine(10, 1, 35, 1)
    line.SetLineColor(ROOT.kBlack)
    line.SetLineStyle(3)
    line.Draw("SAME")

    out_png.parent.mkdir(parents=True, exist_ok=True)
    can.SaveAs(str(out_png))
    # Keep PyROOT objects alive until after SaveAs.
    _ = held_text, legend, line, ratio, hists


def main() -> int:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    points = read_points(args.points_csv)
    out_png = args.outdir / f"{args.tag}.png"
    draw(points, out_png)
    manifest = {
        "artifact": "PPG12 IAN Fig.3 source reproduction from SDCC ROOT points",
        "ppg12_remote_root": "/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_noiso.root",
        "ppg12_objects": ["h_direct_pT_truth_isoET_0", "h_frag_pT_truth_isoET_0"],
        "ppg12_macro": "ppg12codeGit/efficiencytool/FindTruthETCut.C",
        "source_points_csv": str(args.points_csv),
        "png": str(out_png),
        "projection": "sum normal E_T^iso bins with upper edge <= 4 GeV, then Rebin(10) to 1 GeV pT bins; total=direct+fragmentation; bottom=direct/total",
        "style_target": "PPG12 IAN Fig.3 qualitative reproduction; no DataThief input",
    }
    manifest_path = args.outdir / f"{args.tag}_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"png": str(out_png), "manifest": str(manifest_path)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
