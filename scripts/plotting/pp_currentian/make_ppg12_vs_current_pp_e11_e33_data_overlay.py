#!/usr/bin/env python3
"""Overlay PPG12 Fig. 13 data with current RecoilJets pp-data E11/E33.

PPG12 is read from the validated SDCC ROOT-extraction JSON.  Current pp data is
read from a RecoilJets ROOT output and combined over 22-28 GeV using the
inclusive/no-NPB shower-shape histograms.
"""

from __future__ import annotations

import argparse
import json
import math
from array import array
from pathlib import Path

import ROOT


BASE = Path(
    "dataOutput/ppg12PhotonYield/"
    "ppg12_photon_yield_v1_data_20260620"
)
DEFAULT_PPG12_JSON = (
    BASE
    / "shower_shape_reference_validation/fig13_e11_e33/"
    / "ppg12_sdcc_fig13_e11_to_e33_histograms.json"
)
DEFAULT_CURRENT_ROOT = (
    BASE
    / "purity_fig29_comparison/globalmbd_mbddigi_componentmix_20260629_current_pp/"
    / "canary_ppdata_currentcode_20260629/merged_root/"
    / "RecoilJets_ppdata_currentcode_canary_MERGED.root"
)
DEFAULT_OUT = (
    BASE
    / "shower_shape_reference_validation/fig13_e11_e33/"
    / "ppg12_sdcc_vs_current_pp_canary_e11_e33_data_overlay.png"
)

HIST_DIR = "Photon_4_GeV_plus_MBD_NS_geq_1"
CURRENT_HISTS = [
    "h_ss_e11e33_inclusive_pT_22_24",
    "h_ss_e11e33_inclusive_pT_24_26",
    "h_ss_e11e33_inclusive_pT_26_28",
]


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ppg12-json", type=Path, default=DEFAULT_PPG12_JSON)
    ap.add_argument("--current-root", type=Path, default=DEFAULT_CURRENT_ROOT)
    ap.add_argument("--current-label", default="Current pp data canary")
    ap.add_argument("--output", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--canvas-width", type=int, default=1050)
    ap.add_argument("--canvas-height", type=int, default=900)
    return ap.parse_args()


def ppg12_graph(payload: dict) -> ROOT.TGraphErrors:
    data = payload["data"]
    xs = array("d", [float(x) for x in data["centers"]])
    ys = array("d", [float(y) for y in data["values"]])
    ex = array("d", [0.0 for _ in xs])
    ey = array("d", [float(e) for e in data["errors"]])
    g = ROOT.TGraphErrors(len(xs), xs, ys, ex, ey)
    g.SetName("g_ppg12_data")
    g.SetMarkerStyle(20)
    g.SetMarkerSize(1.0)
    g.SetMarkerColor(ROOT.kBlack)
    g.SetLineColor(ROOT.kBlack)
    g.SetLineWidth(2)
    return g


def current_graph(root_path: Path) -> tuple[ROOT.TGraphErrors, ROOT.TGraphErrors, float]:
    f = ROOT.TFile.Open(str(root_path))
    if not f or f.IsZombie():
        raise SystemExit(f"Could not open current ROOT: {root_path}")

    combined = None
    raw_entries = 0.0
    missing: list[str] = []
    for hist_name in CURRENT_HISTS:
        full = f"{HIST_DIR}/{hist_name}"
        h = f.Get(full)
        if not h:
            missing.append(full)
            continue
        raw_entries += float(h.Integral())
        if combined is None:
            combined = h.Clone("h_current_combined_raw")
            combined.SetDirectory(0)
        else:
            combined.Add(h)
    f.Close()
    if missing:
        raise SystemExit(f"Missing current histograms: {missing}")
    if combined is None or combined.Integral() <= 0:
        raise SystemExit("Current pp combined histogram is empty")

    # Current histograms are 0.01 wide over 0-1.2; PPG12 is 0.04 wide over 0-1.
    rebinned = combined.Rebin(4, "h_current_rebinned_004")
    rebinned.SetDirectory(0)
    first = rebinned.GetXaxis().FindBin(0.000001)
    last = rebinned.GetXaxis().FindBin(0.999999)
    norm = float(rebinned.Integral(first, last))
    if norm <= 0:
        raise SystemExit("Current pp has no 0-1 entries after rebin")

    xs: list[float] = []
    ys: list[float] = []
    ex: list[float] = []
    ey: list[float] = []
    for ibin in range(first, last + 1):
        xlo = rebinned.GetXaxis().GetBinLowEdge(ibin)
        xhi = rebinned.GetXaxis().GetBinUpEdge(ibin)
        xs.append(0.5 * (xlo + xhi))
        content = float(rebinned.GetBinContent(ibin))
        err = float(rebinned.GetBinError(ibin))
        ys.append(content / norm)
        ex.append(0.0)
        ey.append(err / norm)

    g = ROOT.TGraphErrors(
        len(xs),
        array("d", xs),
        array("d", ys),
        array("d", ex),
        array("d", ey),
    )
    g.SetName("g_current_pp_data")
    g.SetMarkerStyle(24)
    g.SetMarkerSize(1.15)
    g.SetMarkerColor(ROOT.kBlue + 1)
    g.SetLineColor(ROOT.kBlue + 1)
    g.SetLineWidth(2)

    return g, rebinned, raw_entries


def ratio_graph(current: ROOT.TGraphErrors, ppg12: ROOT.TGraphErrors) -> ROOT.TGraphErrors:
    xs: list[float] = []
    ys: list[float] = []
    ex: list[float] = []
    ey: list[float] = []
    for i in range(current.GetN()):
        x = current.GetPointX(i)
        y = current.GetPointY(i)
        ey_cur = current.GetErrorY(i)
        y_ref = ppg12.GetPointY(i)
        ey_ref = ppg12.GetErrorY(i)
        if y_ref <= 0:
            continue
        ratio = y / y_ref
        # Standard independent relative-error propagation for visual QA.
        rel2 = 0.0
        if y > 0:
            rel2 += (ey_cur / y) ** 2
        if y_ref > 0:
            rel2 += (ey_ref / y_ref) ** 2
        xs.append(x)
        ys.append(ratio)
        ex.append(0.0)
        ey.append(abs(ratio) * math.sqrt(rel2))

    g = ROOT.TGraphErrors(
        len(xs),
        array("d", xs),
        array("d", ys),
        array("d", ex),
        array("d", ey),
    )
    g.SetName("g_current_over_ppg12")
    g.SetMarkerStyle(20)
    g.SetMarkerSize(0.9)
    g.SetMarkerColor(ROOT.kBlue + 1)
    g.SetLineColor(ROOT.kBlue + 1)
    g.SetLineWidth(2)
    return g


def main() -> None:
    args = parse_args()
    with args.ppg12_json.open() as f:
        payload = json.load(f)

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetTextFont(42)
    ROOT.gStyle.SetTitleFont(42, "XYZ")
    ROOT.gStyle.SetLabelFont(42, "XYZ")
    ROOT.gStyle.SetEndErrorSize(0)

    ppg12 = ppg12_graph(payload)
    current, _, raw_entries = current_graph(args.current_root)
    ratio = ratio_graph(current, ppg12)

    canvas = ROOT.TCanvas(
        "c_ppg12_vs_current_pp_e11e33_data",
        "",
        int(args.canvas_width),
        int(args.canvas_height),
    )
    canvas.SetFillColor(ROOT.kWhite)
    top = ROOT.TPad("top", "", 0.0, 0.30, 1.0, 1.0)
    bottom = ROOT.TPad("bottom", "", 0.0, 0.0, 1.0, 0.30)
    for pad in (top, bottom):
        pad.SetFillColor(ROOT.kWhite)
        pad.SetTicks(1, 1)
        pad.SetLeftMargin(0.13)
        pad.SetRightMargin(0.04)
    top.SetBottomMargin(0.02)
    top.SetTopMargin(0.09)
    bottom.SetTopMargin(0.04)
    bottom.SetBottomMargin(0.33)
    top.Draw()
    bottom.Draw()

    top.cd()
    frame = top.DrawFrame(0.0, 0.0, 1.0, 0.18)
    frame.SetTitle("")
    frame.GetYaxis().SetTitle("unit-normalized counts")
    frame.GetYaxis().SetTitleSize(0.055)
    frame.GetYaxis().SetLabelSize(0.045)
    frame.GetYaxis().SetTitleOffset(1.08)
    frame.GetXaxis().SetLabelSize(0.0)
    ppg12.Draw("P SAME")
    current.Draw("P SAME")

    latex = ROOT.TLatex()
    latex.SetNDC(True)
    latex.SetTextFont(42)
    latex.SetTextSize(0.034)
    latex.DrawLatex(0.16, 0.875, "#it{#bf{sPHENIX}} Internal")
    latex.SetTextSize(0.028)
    latex.DrawLatex(0.16, 0.818, "p+p #sqrt{s}=200 GeV, |#eta^{#gamma}| < 0.7")
    latex.DrawLatex(0.16, 0.768, f"E11/E33 data only, 22 < E_{{T}} < 28 GeV, no NPB cut; current N={raw_entries:.0f}")

    leg = ROOT.TLegend(0.56, 0.73, 0.94, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.030)
    leg.AddEntry(ppg12, "PPG12 SDCC data", "pe")
    leg.AddEntry(current, args.current_label, "pe")
    leg.Draw()
    top.RedrawAxis()

    bottom.cd()
    rframe = bottom.DrawFrame(0.0, 0.0, 1.0, 2.6)
    rframe.SetTitle("")
    rframe.GetXaxis().SetTitle("e11_to_e33")
    rframe.GetXaxis().SetTitleSize(0.105)
    rframe.GetXaxis().SetLabelSize(0.080)
    rframe.GetXaxis().SetTitleOffset(1.10)
    rframe.GetYaxis().SetTitle("Current / PPG12")
    rframe.GetYaxis().SetTitleSize(0.085)
    rframe.GetYaxis().SetLabelSize(0.070)
    rframe.GetYaxis().SetTitleOffset(0.63)
    rframe.GetYaxis().SetNdivisions(505)
    one = ROOT.TLine(0.0, 1.0, 1.0, 1.0)
    one.SetLineColor(ROOT.kBlack)
    one.SetLineStyle(2)
    one.SetLineWidth(2)
    one.Draw("SAME")
    ratio.Draw("P SAME")
    bottom.RedrawAxis()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    canvas.SaveAs(str(args.output))
    print(args.output)


if __name__ == "__main__":
    main()
