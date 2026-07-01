#!/usr/bin/env python3
"""Replot the PPG12 Fig. 13 E11/E33 panel from extracted SDCC ROOT histograms.

The input JSON is produced from Shuhang's SDCC ROOT files.  This helper does
not read pixels from the IAN PDF or screenshot; it draws the ROOT-extracted bin
contents with a layout that keeps labels and legend off the data.
"""

from __future__ import annotations

import argparse
import json
from array import array
from pathlib import Path

import ROOT


DEFAULT_DIR = Path(
    "dataOutput/ppg12PhotonYield/"
    "ppg12_photon_yield_v1_data_20260620/"
    "shower_shape_reference_validation/fig13_e11_e33"
)
DEFAULT_JSON = DEFAULT_DIR / "ppg12_sdcc_fig13_e11_to_e33_histograms.json"
DEFAULT_OUT = DEFAULT_DIR / "ppg12_sdcc_root_replot_fig13_e11_to_e33_clean_no_overlap.png"


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", type=Path, default=DEFAULT_JSON)
    ap.add_argument("--output", type=Path, default=DEFAULT_OUT)
    return ap.parse_args()


def make_hist(name: str, payload: dict, color: int, width: int = 3) -> ROOT.TH1D:
    edges = array("d", [float(x) for x in payload["edges"]])
    h = ROOT.TH1D(name, "", len(edges) - 1, edges)
    h.SetDirectory(0)
    for i, value in enumerate(payload["values"], start=1):
        h.SetBinContent(i, float(value))
        h.SetBinError(i, float(payload["errors"][i - 1]))
    h.SetLineColor(color)
    h.SetLineWidth(width)
    h.SetFillStyle(0)
    h.SetStats(False)
    return h


def make_graph(name: str, payload: dict) -> ROOT.TGraphErrors:
    xs = array("d", [float(x) for x in payload["centers"]])
    ys = array("d", [float(y) for y in payload["values"]])
    ex = array("d", [0.0 for _ in xs])
    ey = array("d", [float(e) for e in payload["errors"]])
    g = ROOT.TGraphErrors(len(xs), xs, ys, ex, ey)
    g.SetName(name)
    g.SetMarkerStyle(20)
    g.SetMarkerSize(0.95)
    g.SetMarkerColor(ROOT.kBlack)
    g.SetLineColor(ROOT.kBlack)
    g.SetLineWidth(2)
    return g


def style_axis(axis: ROOT.TAxis, title: str, title_size: float, label_size: float) -> None:
    axis.SetTitle(title)
    axis.SetTitleSize(title_size)
    axis.SetLabelSize(label_size)
    axis.SetTitleFont(42)
    axis.SetLabelFont(42)
    axis.CenterTitle(False)


def main() -> None:
    args = parse_args()
    with args.input.open() as f:
        payload = json.load(f)

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetTextFont(42)
    ROOT.gStyle.SetTitleFont(42, "XYZ")
    ROOT.gStyle.SetLabelFont(42, "XYZ")
    ROOT.gStyle.SetEndErrorSize(0)

    signal = make_hist("h_signal_sdcc", payload["signal"], ROOT.kRed + 1)
    inclusive = make_hist("h_inclusive_sdcc", payload["inclusive"], ROOT.kBlue + 1)
    npb = make_hist("h_npb_sdcc", payload["npb"], ROOT.kGreen + 2)
    data = make_graph("g_data_sdcc", payload["data"])
    residual = make_graph("g_residual_sdcc", payload["residual"])

    canvas = ROOT.TCanvas("c_ppg12_fig13_e11e33_clean", "", 1150, 1000)
    canvas.SetFillColor(ROOT.kWhite)
    canvas.SetFrameFillColor(ROOT.kWhite)

    top = ROOT.TPad("top", "", 0.0, 0.28, 1.0, 1.0)
    bottom = ROOT.TPad("bottom", "", 0.0, 0.0, 1.0, 0.28)
    for pad in (top, bottom):
        pad.SetFillColor(ROOT.kWhite)
        pad.SetTicks(1, 1)
        pad.SetLeftMargin(0.14)
        pad.SetRightMargin(0.34)
    top.SetTopMargin(0.08)
    top.SetBottomMargin(0.02)
    bottom.SetTopMargin(0.04)
    bottom.SetBottomMargin(0.34)
    top.Draw()
    bottom.Draw()

    top.cd()
    frame = top.DrawFrame(0.0, 0.0, 1.0, 0.235)
    frame.SetTitle("")
    frame.SetStats(False)
    style_axis(frame.GetYaxis(), "normalized counts", 0.055, 0.045)
    frame.GetYaxis().SetTitleOffset(1.18)
    frame.GetYaxis().SetNdivisions(506)
    frame.GetXaxis().SetLabelSize(0.0)
    frame.GetXaxis().SetTickLength(0.03)

    signal.Draw("HIST SAME")
    inclusive.Draw("HIST SAME")
    npb.Draw("HIST SAME")
    data.Draw("P SAME")

    latex = ROOT.TLatex()
    latex.SetNDC(False)
    latex.SetTextFont(42)
    latex.SetTextSize(0.034)
    latex.DrawLatex(0.050, 0.223, "#it{#bf{sPHENIX}} Internal")
    latex.SetTextSize(0.0275)
    latex.DrawLatex(0.050, 0.207, "p+p #sqrt{s}=200 GeV")
    latex.DrawLatex(0.050, 0.194, "|#eta^{#gamma}| < 0.7")
    latex.DrawLatex(0.050, 0.181, "22 < p_{T} < 28 GeV, w/o nbkg cut")
    latex.DrawLatex(
        0.050,
        0.168,
        f"#chi^{{2}}/ndf = {payload['chi2']:.1f}/{int(payload['ndf'])} = {payload['chi2_ndf']:.2f}",
    )
    latex.DrawLatex(0.050, 0.155, "p-value = 0.0000")

    leg = ROOT.TLegend(0.690, 0.710, 0.985, 0.925)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.040)
    leg.AddEntry(data, "Data", "pe")
    leg.AddEntry(signal, "Signal MC", "l")
    leg.AddEntry(inclusive, "Inclusive MC", "l")
    leg.AddEntry(npb, "NPB-tagged data", "l")
    leg.Draw()
    top.RedrawAxis()

    bottom.cd()
    rframe = bottom.DrawFrame(0.0, -0.075, 1.0, 0.125)
    rframe.SetTitle("")
    rframe.SetStats(False)
    style_axis(rframe.GetXaxis(), "e11_to_e33", 0.11, 0.095)
    style_axis(rframe.GetYaxis(), "Data - Incl. MC", 0.10, 0.075)
    rframe.GetXaxis().SetTitleOffset(1.13)
    rframe.GetYaxis().SetTitleOffset(0.58)
    rframe.GetYaxis().SetNdivisions(405)
    zero = ROOT.TLine(0.0, 0.0, 1.0, 0.0)
    zero.SetLineColor(ROOT.kBlack)
    zero.SetLineStyle(2)
    zero.SetLineWidth(2)
    zero.Draw("SAME")
    residual.Draw("P SAME")
    bottom.RedrawAxis()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    canvas.SaveAs(str(args.output))
    print(args.output)


if __name__ == "__main__":
    main()
