#!/usr/bin/env python3
"""Reproduce the PPG12 paper Fig. 3 isolation template with THE97 outputs.

The histogram summation, variable rebinning, tail normalization, signal-MC
normalization, stack order, colors, and displayed range mirror
``ppg12codeGit/plotting/CONF_plots.C`` and
``ppg12codeGit/plotting/paper/plot_paper_iso_template.C``.
"""

from __future__ import annotations

import argparse
from array import array
import csv
import hashlib
import json
from pathlib import Path

import ROOT


REPO = Path(__file__).resolve().parents[3]
PT_INDEXES = (3, 4, 5)
TIGHT_BASE = "h_tight_isoET_0_"
NONTIGHT_BASE = "h_nontight_isoET_0_"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def open_root(path: Path) -> ROOT.TFile:
    root_file = ROOT.TFile.Open(str(path))
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"could not open ROOT file: {path}")
    if root_file.TestBit(ROOT.TFile.kRecovered):
        raise RuntimeError(f"refusing recovered ROOT file: {path}")
    return root_file


def sum_histograms(paths: list[Path], directory: str, base: str, name: str) -> ROOT.TH1:
    result = None
    files = []
    try:
        for path in paths:
            root_file = open_root(path)
            files.append(root_file)
            for index in PT_INDEXES:
                key = f"{directory}/{base}{index}"
                hist = root_file.Get(key)
                if not hist:
                    raise KeyError(f"missing {key} in {path}")
                if result is None:
                    result = hist.Clone(name)
                    result.SetDirectory(0)
                    if result.GetSumw2N() == 0:
                        result.Sumw2()
                else:
                    if hist.GetNbinsX() != result.GetNbinsX():
                        raise ValueError(f"incompatible bin count for {key} in {path}")
                    result.Add(hist)
    finally:
        for root_file in files:
            root_file.Close()
    if result is None:
        raise RuntimeError(f"no histograms found for {directory}/{base}{{3,4,5}}")
    return result


def variable_rebin(hist: ROOT.TH1, name: str) -> ROOT.TH1:
    edges = [float(hist.GetXaxis().GetBinLowEdge(1))]
    index = 1
    while index <= hist.GetNbinsX():
        group_size = 1 if hist.GetXaxis().GetBinLowEdge(index) < 2.5 else 5
        last = min(index + group_size - 1, hist.GetNbinsX())
        edges.append(
            float(
                hist.GetXaxis().GetBinLowEdge(last)
                + hist.GetXaxis().GetBinWidth(last)
            )
        )
        index += group_size
    rebinned = hist.Rebin(len(edges) - 1, name, array("d", edges))
    rebinned.SetDirectory(0)
    rebinned.Scale(1.0, "width")
    return rebinned


def make_template(args: argparse.Namespace) -> tuple[ROOT.TH1, ROOT.TH1, ROOT.TH1, dict[str, float]]:
    tight = sum_histograms(args.data_root, args.data_dir, TIGHT_BASE, "h_the97_tight_sum")
    nontight = sum_histograms(args.data_root, args.data_dir, NONTIGHT_BASE, "h_the97_nontight_sum")
    signal_mc = sum_histograms([args.mc_root], args.mc_dir, TIGHT_BASE, "h_the97_signal_mc_sum")

    tight = variable_rebin(tight, "h_the97_tight_rebinned")
    nontight = variable_rebin(nontight, "h_the97_nontight_rebinned")
    signal_mc = variable_rebin(signal_mc, "h_the97_signal_mc_rebinned")

    if not (
        tight.GetNbinsX() == nontight.GetNbinsX() == signal_mc.GetNbinsX()
    ):
        raise ValueError("rebinned data and MC bin counts do not match")
    for bin_index in range(1, tight.GetNbinsX() + 1):
        x = tight.GetXaxis().GetBinCenter(bin_index)
        if abs(x - nontight.GetXaxis().GetBinCenter(bin_index)) > 1.0e-9:
            raise ValueError("tight and non-tight binning do not match")
        if abs(x - signal_mc.GetXaxis().GetBinCenter(bin_index)) > 1.0e-9:
            raise ValueError("data and signal-MC binning do not match")

    tail_first_bin = tight.FindBin(args.tail_threshold)
    tail_last_bin = tight.GetNbinsX()
    tight_tail = float(tight.Integral(tail_first_bin, tail_last_bin))
    nontight_tail_before = float(nontight.Integral(tail_first_bin, tail_last_bin))
    if nontight_tail_before <= 0:
        raise RuntimeError("non-tight tail integral is not positive")
    background_scale = tight_tail / nontight_tail_before
    nontight.Scale(background_scale)

    tight_integral = float(tight.Integral())
    background_integral = float(nontight.Integral())
    data_signal = tight_integral - background_integral
    mc_integral_before = float(signal_mc.Integral())
    if data_signal <= 0 or mc_integral_before <= 0:
        raise RuntimeError(
            f"invalid PPG12 normalization: data_signal={data_signal}, mc={mc_integral_before}"
        )
    mc_scale = data_signal / mc_integral_before
    signal_mc.Scale(mc_scale)

    metrics = {
        "tail_threshold_gev": float(args.tail_threshold),
        "tail_first_bin": int(tail_first_bin),
        "tight_tail_integral": tight_tail,
        "nontight_tail_integral_before_scale": nontight_tail_before,
        "background_scale": background_scale,
        "tight_integral": tight_integral,
        "background_integral_after_scale": background_integral,
        "data_signal_difference": data_signal,
        "mc_integral_before_scale": mc_integral_before,
        "mc_scale": mc_scale,
        "mc_integral_after_scale": float(signal_mc.Integral()),
    }
    return tight, nontight, signal_mc, metrics


def draw(
    out_png: Path,
    tight: ROOT.TH1,
    nontight: ROOT.TH1,
    signal_mc: ROOT.TH1,
    args: argparse.Namespace,
) -> None:
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetTextFont(42)
    ROOT.gStyle.SetLabelFont(42, "XYZ")
    ROOT.gStyle.SetTitleFont(42, "XYZ")

    tight.SetStats(False)
    tight.SetTitle("")
    tight.SetMarkerStyle(20)
    tight.SetMarkerSize(0.90)
    tight.SetMarkerColor(ROOT.kBlack)
    tight.SetLineColor(ROOT.kBlack)
    tight.SetLineWidth(1)

    nontight.SetStats(False)
    nontight.SetFillColorAlpha(ROOT.kRed, 0.30)
    nontight.SetFillStyle(1001)
    nontight.SetLineColor(ROOT.kRed - 1)
    nontight.SetLineWidth(2)

    signal_mc.SetStats(False)
    signal_mc.SetFillColorAlpha(ROOT.kBlue, 0.30)
    signal_mc.SetFillStyle(1001)
    signal_mc.SetLineColor(ROOT.kBlue - 1)
    signal_mc.SetLineWidth(2)

    canvas = ROOT.TCanvas("c_the97_ppg12_fig3_iso", "", 900, 840)
    canvas.SetLeftMargin(0.16)
    canvas.SetRightMargin(0.045)
    canvas.SetTopMargin(0.055)
    canvas.SetBottomMargin(0.14)
    canvas.SetTicks(1, 1)

    tight.GetXaxis().SetRangeUser(-1.0, 15.0)
    tight.GetXaxis().SetTitle("#it{E}_{T}^{iso,reco} [GeV]")
    tight.GetYaxis().SetTitle("Counts / Bin Width")
    tight.GetXaxis().SetTitleSize(0.052)
    tight.GetYaxis().SetTitleSize(0.052)
    tight.GetXaxis().SetLabelSize(0.043)
    tight.GetYaxis().SetLabelSize(0.043)
    tight.GetXaxis().SetTitleOffset(1.12)
    tight.GetYaxis().SetTitleOffset(1.34)
    tight.GetXaxis().CenterTitle(False)
    tight.GetYaxis().CenterTitle(False)
    ymax = max(
        float(tight.GetMaximum()),
        max(
            float(nontight.GetBinContent(i) + signal_mc.GetBinContent(i))
            for i in range(1, tight.GetNbinsX() + 1)
        ),
    )
    tight.SetMinimum(0.0)
    tight.SetMaximum(ymax * 1.16)
    tight.Draw("E X0")

    stack = ROOT.THStack("hs_the97_ppg12_fig3_iso", "")
    stack.Add(nontight)
    stack.Add(signal_mc)
    stack.Draw("HIST SAME")
    tight.Draw("E X0 SAME")
    tight.Draw("AXIS SAME")

    latex = ROOT.TLatex()
    latex.SetNDC(True)
    latex.SetTextFont(42)
    latex.SetTextAlign(31)
    latex.SetTextSize(0.046)
    latex.DrawLatex(0.91, 0.89, "#bf{#it{sPHENIX}} Internal")
    latex.SetTextSize(0.043)
    latex.DrawLatex(0.91, 0.835, "#it{p}+#it{p} #sqrt{#it{s}} = 200 GeV")
    latex.DrawLatex(0.91, 0.780, "16 < #it{E}_{T}^{#gamma} < 22 GeV")
    latex.DrawLatex(0.91, 0.725, "|#eta^{#gamma}| < 0.7")

    legend = ROOT.TLegend(0.55, 0.515, 0.91, 0.685)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.042)
    legend.AddEntry(tight, "Data (Signal)", "pe")
    legend.AddEntry(nontight, "Data (Background)", "f")
    legend.AddEntry(signal_mc, "Signal MC", "f")
    legend.Draw()

    coverage = 100.0 * args.partial_completed / args.partial_expected
    coverage_note = ROOT.TLatex()
    coverage_note.SetNDC(True)
    coverage_note.SetTextFont(42)
    coverage_note.SetTextSize(0.027)
    coverage_note.SetTextAlign(11)
    if args.partial_completed != args.partial_expected:
        coverage_note.SetTextColor(ROOT.kRed + 1)
        coverage_note.DrawLatex(0.56, 0.452, "PRELIMINARY")
        coverage_note.DrawLatex(0.56, 0.416, f"PARTIAL-COVERAGE: {coverage:.3f}%")
        coverage_note.DrawLatex(0.56, 0.380, f"{args.partial_completed:,} / {args.partial_expected:,} files")

    out_png.parent.mkdir(parents=True, exist_ok=True)
    canvas.SaveAs(str(out_png))
    _ = stack, legend, latex, coverage_note


def write_csv(path: Path, tight: ROOT.TH1, nontight: ROOT.TH1, signal_mc: ROOT.TH1) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "bin_low",
                "bin_high",
                "bin_center",
                "data_signal_selection_counts_per_width",
                "data_signal_selection_error",
                "data_background_counts_per_width",
                "data_background_error",
                "signal_mc_counts_per_width",
                "signal_mc_error",
                "stack_total_counts_per_width",
            ]
        )
        for i in range(1, tight.GetNbinsX() + 1):
            writer.writerow(
                [
                    tight.GetXaxis().GetBinLowEdge(i),
                    tight.GetXaxis().GetBinUpEdge(i),
                    tight.GetXaxis().GetBinCenter(i),
                    tight.GetBinContent(i),
                    tight.GetBinError(i),
                    nontight.GetBinContent(i),
                    nontight.GetBinError(i),
                    signal_mc.GetBinContent(i),
                    signal_mc.GetBinError(i),
                    nontight.GetBinContent(i) + signal_mc.GetBinContent(i),
                ]
            )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, nargs="+", required=True)
    parser.add_argument("--mc-root", type=Path, required=True)
    parser.add_argument("--data-dir", default="PPG12_scaledtrigger30")
    parser.add_argument("--mc-dir", default="SIM")
    parser.add_argument("--tail-threshold", type=float, default=6.0)
    parser.add_argument("--partial-completed", type=int, required=True)
    parser.add_argument("--partial-expected", type=int, required=True)
    parser.add_argument("--campaign-tag", required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--tag", default="the97_ppg12_fig3_iso_template_PRELIMINARY_PARTIAL_COVERAGE")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.data_root = [path.resolve() for path in args.data_root]
    args.mc_root = args.mc_root.resolve()
    if not (0 < args.partial_completed <= args.partial_expected):
        raise ValueError("partial coverage must satisfy 0 < completed <= expected")

    tight, nontight, signal_mc, metrics = make_template(args)
    args.outdir.mkdir(parents=True, exist_ok=True)
    png_path = args.outdir / f"{args.tag}.png"
    csv_path = args.outdir / f"{args.tag}.csv"
    manifest_path = args.outdir / f"{args.tag}_manifest.json"
    draw(png_path, tight, nontight, signal_mc, args)
    write_csv(csv_path, tight, nontight, signal_mc)

    coverage = 100.0 * args.partial_completed / args.partial_expected
    manifest = {
        "artifact": "THE97 current-output reproduction of PPG12 paper Fig. 3 isolation template",
        "status": (
            "FULL_STAT_CURRENT"
            if args.partial_completed == args.partial_expected
            else "PRELIMINARY PARTIAL-COVERAGE"
        ),
        "campaign": args.campaign_tag,
        "coverage": {
            "completed_files": args.partial_completed,
            "expected_files": args.partial_expected,
            "percent": coverage,
        },
        "selection": {
            "photon_et_gev": "16 < E_T^gamma < 22",
            "eta": "|eta^gamma| < 0.7",
            "pt_indexes_summed": list(PT_INDEXES),
            "data_directory": args.data_dir,
            "mc_directory": args.mc_dir,
        },
        "ppg12_contract": {
            "macro": "ppg12codeGit/plotting/CONF_plots.C",
            "paper_macro": "ppg12codeGit/plotting/paper/plot_paper_iso_template.C",
            "variable_rebin": "group 1 original bin below 2.5 GeV, group 5 otherwise; divide by rebinned bin width",
            "background": "scale non-tight data so its width-scaled integral from FindBin(6 GeV) through the last normal bin equals tight data",
            "signal_mc": "scale tight photon+jet MC width-scaled integral to tight-data integral minus scaled non-tight-data integral",
            "stack_order": ["Data (Background)", "Signal MC"],
            "x_range_gev": [-1.0, 15.0],
        },
        "inputs": {
            "data_roots": [str(path) for path in args.data_root],
            "data_sha256": {str(path): sha256_file(path) for path in args.data_root},
            "mc_root": str(args.mc_root),
            "mc_sha256": sha256_file(args.mc_root),
            "data_keys": [
                f"{args.data_dir}/{base}{index}"
                for base in (TIGHT_BASE, NONTIGHT_BASE)
                for index in PT_INDEXES
            ],
            "mc_keys": [f"{args.mc_dir}/{TIGHT_BASE}{index}" for index in PT_INDEXES],
        },
        "normalization_metrics": metrics,
        "outputs": {
            "png": str(png_path.resolve()),
            "csv": str(csv_path.resolve()),
            "manifest": str(manifest_path.resolve()),
        },
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(manifest["outputs"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
