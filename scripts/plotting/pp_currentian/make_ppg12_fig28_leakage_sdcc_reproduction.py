#!/usr/bin/env python3
"""Reproduce PPG12 IAN Fig. 28 signal-leakage panel from SDCC ROOT data.

The target PPG12 macro is ppg12codeGit/plotting/plot_sideband_selection.C:
it draws the signal-MC leakage fractions

  B/A = h_tight_noniso_cluster_0 / h_tight_iso_cluster_0
  C/A = h_nontight_iso_cluster_0 / h_tight_iso_cluster_0
  D/A = h_nontight_noniso_cluster_0 / h_tight_iso_cluster_0

from Shuhang's MC_efficiency_<suffix>.root files.
"""

from __future__ import annotations

import argparse
import csv
import json
from array import array
from pathlib import Path

import ROOT


REPO = Path(__file__).resolve().parents[3]
DEFAULT_IN = (
    REPO
    / "dataOutput/ppg12Parity/the76_ppg12_fig28_leakage_reference"
    / "ppg12_fig28_leakage_sdcc_points.csv"
)
DEFAULT_OUTDIR = REPO / "dataOutput/ppg12Parity/the76_ppg12_fig28_leakage_reference"


SERIES = [
    (
        "B",
        "tight_noniso",
        ROOT.kBlack,
        "#it{N}^{sig}_{B}/#it{N}^{sig}_{A} tight noniso",
    ),
    (
        "C",
        "nontight_iso",
        ROOT.kRed,
        "#it{N}^{sig}_{C}/#it{N}^{sig}_{A} nontight iso",
    ),
    (
        "D",
        "nontight_noniso",
        ROOT.kBlue,
        "#it{N}^{sig}_{D}/#it{N}^{sig}_{A} nontight noniso",
    ),
]


def read_rows(path: Path, source_kind: str) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        rows = [row for row in csv.DictReader(handle) if row["source_kind"] == source_kind]
    if not rows:
        raise RuntimeError(f"no rows with source_kind={source_kind!r} in {path}")
    return rows


def make_hist(rows: list[dict[str, str]], region: str, name: str) -> ROOT.TH1D:
    region_rows = [row for row in rows if row["region"] == region]
    region_rows.sort(key=lambda row: int(row["bin_index"]))
    if not region_rows:
        raise RuntimeError(f"no rows for region {region}")
    edges = [float(region_rows[0]["x_low_gev"])]
    edges.extend(float(row["x_high_gev"]) for row in region_rows)
    hist = ROOT.TH1D(name, "", len(edges) - 1, array("d", edges))
    hist.SetDirectory(0)
    for ibin, row in enumerate(region_rows, start=1):
        hist.SetBinContent(ibin, float(row["leakage"]))
        hist.SetBinError(ibin, float(row["error"]))
    return hist


def draw_line_legend(text_x: float, y: float, color: int, text: str, size: float) -> None:
    """Mirror PPG12 BlairUtils::myMarkerLineText with marker size zero."""
    line = ROOT.TLine()
    line.SetLineColor(color)
    line.SetLineWidth(2)
    line.DrawLineNDC(text_x - 0.95 * size, y, text_x - 0.15 * size, y)
    latex = ROOT.TLatex()
    latex.SetNDC(True)
    latex.SetTextFont(42)
    latex.SetTextSize(size)
    latex.SetTextAlign(12)
    latex.DrawLatex(text_x, y, text)


def render(rows: list[dict[str, str]], png: Path, source_kind: str, style: str) -> dict[str, object]:
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetTextFont(42)
    ROOT.gStyle.SetLabelFont(42, "XYZ")
    ROOT.gStyle.SetTitleFont(42, "XYZ")
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)

    hists = {}
    for region, _, color, _ in SERIES:
        hist = make_hist(rows, region, f"h_fig28_{region}")
        hist.SetLineColor(color)
        hist.SetLineWidth(2)
        hist.SetMarkerSize(0)
        hists[region] = hist

    canvas = ROOT.TCanvas("c_ppg12_fig28_leakage", "", 435, 439)
    canvas.SetLeftMargin(0.16)
    canvas.SetRightMargin(0.025)
    canvas.SetTopMargin(0.08)
    canvas.SetBottomMargin(0.15)
    canvas.SetTicks(1, 1)

    frame = ROOT.TH1F("frame_fig28", "", 43, 7.0, 50.0)
    frame.SetStats(False)
    if style == "sdcc_pdf":
        # Match the IAN screenshot frame Justin is comparing by eye.  This
        # exposes the 32-36 GeV C/A source bin present in the SDCC ROOT data.
        frame.GetXaxis().SetRangeUser(10.0, 35.0)
        frame.GetYaxis().SetRangeUser(0.0, 1.3)
        style_target = "PPG12 IAN screenshot axes: x=10-35 GeV, y=0-1.3"
    else:
        # Diagnostic exposure of the full high bin.  This is useful to prove
        # the source ROOT contents, but it is not the visual match to the IAN
        # screenshot/note panel because it reveals the 32-36 GeV C/A bin.
        frame.GetXaxis().SetRangeUser(10.0, 35.0)
        frame.GetYaxis().SetRangeUser(0.0, 1.3)
        style_target = "Extended diagnostic axes: x=10-35 GeV, y=0-1.3"
    frame.GetXaxis().SetTitle("#it{E}_{T}^{#gamma,rec} [GeV]")
    frame.GetYaxis().SetTitle("Signal leakage")
    frame.GetXaxis().SetTitleSize(0.052)
    frame.GetYaxis().SetTitleSize(0.052)
    frame.GetXaxis().SetLabelSize(0.047)
    frame.GetYaxis().SetLabelSize(0.047)
    frame.GetXaxis().SetTitleOffset(1.02)
    frame.GetYaxis().SetTitleOffset(1.15)
    frame.GetXaxis().SetNdivisions(505)
    frame.GetYaxis().SetNdivisions(507)
    frame.Draw("axis")

    for region, _, _, _ in SERIES:
        hists[region].Draw("same hist")

    latex = ROOT.TLatex()
    latex.SetNDC(True)
    latex.SetTextFont(42)
    latex.SetTextSize(0.040)
    latex.SetTextAlign(12)
    latex.DrawLatex(0.17, 0.885, "#bf{#it{sPHENIX}} Internal")
    latex.DrawLatex(0.17, 0.830, "#it{p}+#it{p} #kern[-0.1]{#sqrt{#it{s}} = 200 GeV}")
    latex.DrawLatex(0.17, 0.775, "|#it{#eta^{#gamma}}| < 0.7")
    latex.DrawLatex(0.17, 0.720, "PYTHIA Signal")

    for y, (_, _, color, label) in zip([0.875, 0.795, 0.715], SERIES):
        draw_line_legend(0.50, y, color, label, 0.042)

    canvas.RedrawAxis()
    png.parent.mkdir(parents=True, exist_ok=True)
    canvas.SaveAs(str(png))

    sources = sorted({row["source_root"] for row in rows})
    return {
        "png": str(png),
        "input_csv": str(args.input),
        "source_kind": source_kind,
        "source_roots": sources,
        "ppg12_macro_reference": "ppg12codeGit/plotting/plot_sideband_selection.C lines 114-134",
        "style_target": style_target,
        "definition": {
            "B": "h_tight_noniso_cluster_0 / h_tight_iso_cluster_0",
            "C": "h_nontight_iso_cluster_0 / h_tight_iso_cluster_0",
            "D": "h_nontight_noniso_cluster_0 / h_tight_iso_cluster_0",
        },
    }


def main() -> None:
    global args
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_IN)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument(
        "--source-kind",
        choices=["macro_recomputed", "persisted_final"],
        default="macro_recomputed",
        help="macro_recomputed follows plot_sideband_selection.C directly; "
        "persisted_final reads h_leak_B/C/D from Photon_final_bdt_nom_mc.root.",
    )
    parser.add_argument(
        "--style",
        choices=["sdcc_pdf", "extended"],
        default="sdcc_pdf",
        help="sdcc_pdf matches the PPG12 analysis-note PDF render; extended "
        "exposes the 32-36 GeV source bin.",
    )
    args = parser.parse_args()

    rows = read_rows(args.input, args.source_kind)
    tag = "macro_source" if args.source_kind == "macro_recomputed" else "persisted_final"
    style_tag = "note_style" if args.style == "sdcc_pdf" else "extended_axes"
    png = args.outdir / f"ppg12_fig28_signal_leakage_{tag}_{style_tag}_reproduction.png"
    manifest = render(rows, png, args.source_kind, args.style)
    manifest_path = png.with_suffix(".manifest.json")
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(png)
    print(manifest_path)


if __name__ == "__main__":
    main()
