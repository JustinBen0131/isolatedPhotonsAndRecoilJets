#!/usr/bin/env python3
"""Make Blair-requested background spectra from stitched RecoilJets ABCD hists."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import ROOT


ET_EDGES = [15.0, 17.0, 19.0, 21.0, 23.0, 25.0, 27.0, 30.0, 35.0]
CENT_SUFFIXES = [
    ("0-20%", "_cent_0_20"),
    ("20-50%", "_cent_20_50"),
    ("50-80%", "_cent_50_80"),
]
REGIONS = ["A", "B", "C", "D"]


def root_hist(root_file: ROOT.TFile, name: str) -> ROOT.TH1:
    obj = root_file.Get(f"SIM/{name}")
    if not obj:
        raise KeyError(f"missing SIM/{name}")
    hist = obj.Clone(f"{name}_clone")
    hist.SetDirectory(0)
    return hist


def abcd_sum(root_file: ROOT.TFile, iso_tag: str, cent_suffix: str) -> ROOT.TH1:
    total = None
    for region in REGIONS:
        hist = root_hist(root_file, f"h_pTgamma_ABCD_{region}_{iso_tag}{cent_suffix}")
        if total is None:
            total = hist.Clone(f"abcdsum_{iso_tag}{cent_suffix}")
            total.SetDirectory(0)
        else:
            total.Add(hist)
    if total is None:
        raise RuntimeError("no ABCD histograms found")
    return total


def integrate(hist: ROOT.TH1, lo: float, hi: float) -> tuple[float, float]:
    value = 0.0
    err2 = 0.0
    axis = hist.GetXaxis()
    for ibin in range(1, hist.GetNbinsX() + 1):
        center = axis.GetBinCenter(ibin)
        if center < lo or center >= hi:
            continue
        value += hist.GetBinContent(ibin)
        err = hist.GetBinError(ibin)
        err2 += err * err
    return value, math.sqrt(err2)


def spectra(root_file: ROOT.TFile, iso_tag: str) -> dict[str, list[dict[str, float | str]]]:
    out: dict[str, list[dict[str, float | str]]] = {}
    for label, suffix in CENT_SUFFIXES:
        hist = abcd_sum(root_file, iso_tag, suffix)
        rows = []
        for lo, hi in zip(ET_EDGES[:-1], ET_EDGES[1:]):
            val, err = integrate(hist, lo, hi)
            rows.append(
                {
                    "centrality": label,
                    "et_lo": lo,
                    "et_hi": hi,
                    "et_mid": 0.5 * (lo + hi),
                    "et_half_width": 0.5 * (hi - lo),
                    "yield": val,
                    "stat_err": err,
                }
            )
        out[label] = rows
    return out


def write_csv(path: Path, background: dict[str, list[dict[str, float | str]]], signal=None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["sample", "centrality", "et_lo", "et_hi", "yield", "stat_err"]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for sample_name, spectra_dict in (("background", background), ("signal", signal)):
            if spectra_dict is None:
                continue
            for rows in spectra_dict.values():
                for row in rows:
                    writer.writerow(
                        {
                            "sample": sample_name,
                            "centrality": row["centrality"],
                            "et_lo": row["et_lo"],
                            "et_hi": row["et_hi"],
                            "yield": f"{row['yield']:.9g}",
                            "stat_err": f"{row['stat_err']:.9g}",
                        }
                    )


def graph_from_rows(rows: list[dict[str, float | str]], color: int, marker: int) -> ROOT.TGraphErrors:
    graph = ROOT.TGraphErrors(len(rows))
    for idx, row in enumerate(rows):
        graph.SetPoint(idx, float(row["et_mid"]), float(row["yield"]))
        graph.SetPointError(idx, 0.0, float(row["stat_err"]))
    graph.SetMarkerStyle(marker)
    graph.SetMarkerSize(1.15)
    graph.SetMarkerColor(color)
    graph.SetLineColor(color)
    graph.SetLineWidth(2)
    return graph


def max_y(*spectra_dicts) -> float:
    ymax = 0.0
    for spectra_dict in spectra_dicts:
        if spectra_dict is None:
            continue
        for rows in spectra_dict.values():
            for row in rows:
                ymax = max(ymax, float(row["yield"]) + float(row["stat_err"]))
    return ymax


def draw_label(sample_text: str, note: str, x: float = 0.15) -> None:
    tex = ROOT.TLatex()
    tex.SetNDC(True)
    tex.SetTextFont(42)
    tex.SetTextSize(0.040)
    tex.DrawLatex(x, 0.87, "#it{#bf{sPHENIX}} Internal")
    tex.SetTextSize(0.034)
    tex.DrawLatex(x, 0.82, sample_text)
    tex.SetTextSize(0.030)
    tex.DrawLatex(x, 0.775, note)


def draw_background_only(out: Path, background: dict[str, list[dict[str, float | str]]], iso_tag: str) -> None:
    colors = [ROOT.kAzure + 1, ROOT.kOrange + 7, ROOT.kGreen + 2]
    markers = [20, 21, 22]
    canvas = ROOT.TCanvas("c_background", "", 1050, 780)
    canvas.SetLeftMargin(0.13)
    canvas.SetRightMargin(0.05)
    canvas.SetBottomMargin(0.13)
    canvas.SetTopMargin(0.08)
    frame = canvas.DrawFrame(14.0, 0.0, 36.0, max_y(background) * 1.38)
    frame.GetXaxis().SetTitle("Candidate E_{T} [GeV]")
    frame.GetYaxis().SetTitle("Weighted background candidates")
    frame.GetXaxis().SetTitleSize(0.045)
    frame.GetYaxis().SetTitleSize(0.045)
    frame.GetXaxis().SetLabelSize(0.038)
    frame.GetYaxis().SetLabelSize(0.038)
    frame.SetStats(False)

    leg = ROOT.TLegend(0.67, 0.68, 0.92, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.034)
    keepalive = []
    for idx, (label, rows) in enumerate(background.items()):
        graph = graph_from_rows(rows, colors[idx], markers[idx])
        keepalive.append(graph)
        graph.Draw("P SAME")
        leg.AddEntry(graph, label, "pe")
    leg.Draw()
    draw_label("Au+Au embedded inclusive Jet12+20+30", "ABCD-summed reco-photon candidates; fixed isolation")
    canvas.SetLogy(False)
    canvas.SaveAs(str(out / f"background_abcdsum_spectrum_cent3_{iso_tag}.png"))


def draw_background_fractional_uncertainty(
    out: Path, background: dict[str, list[dict[str, float | str]]], iso_tag: str
) -> None:
    colors = [ROOT.kAzure + 1, ROOT.kOrange + 7, ROOT.kGreen + 2]
    markers = [20, 21, 22]
    canvas = ROOT.TCanvas("c_background_frac", "", 1050, 780)
    canvas.SetLeftMargin(0.13)
    canvas.SetRightMargin(0.05)
    canvas.SetBottomMargin(0.13)
    canvas.SetTopMargin(0.08)

    frac_rows: dict[str, list[dict[str, float | str]]] = {}
    ymax = 0.0
    for label, rows in background.items():
        frac_rows[label] = []
        for row in rows:
            val = float(row["yield"])
            if val <= 0.0:
                continue
            err = float(row["stat_err"])
            frac = err / val if val > 0.0 else 0.0
            ymax = max(ymax, frac)
            new_row = dict(row)
            new_row["yield"] = frac
            new_row["stat_err"] = 0.0
            frac_rows[label].append(new_row)

    frame = canvas.DrawFrame(14.0, 0.0, 36.0, max(0.08, ymax * 1.45))
    frame.GetXaxis().SetTitle("Candidate E_{T} [GeV]")
    frame.GetYaxis().SetTitle("Fractional statistical uncertainty")
    frame.GetXaxis().SetTitleSize(0.045)
    frame.GetYaxis().SetTitleSize(0.045)
    frame.GetXaxis().SetLabelSize(0.038)
    frame.GetYaxis().SetLabelSize(0.038)
    frame.SetStats(False)

    leg = ROOT.TLegend(0.67, 0.68, 0.92, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.034)
    keepalive = []
    for idx, (label, rows) in enumerate(frac_rows.items()):
        graph = graph_from_rows(rows, colors[idx], markers[idx])
        keepalive.append(graph)
        graph.Draw("P SAME")
        leg.AddEntry(graph, label, "p")
    leg.Draw()
    draw_label("Au+Au embedded inclusive Jet12+20+30", "ABCD-summed reco-photon candidates; fixed isolation")
    canvas.SaveAs(str(out / f"background_abcdsum_fractional_stat_uncertainty_cent3_{iso_tag}.png"))


def draw_overlay(out: Path, background, signal, iso_tag: str) -> None:
    b_colors = [ROOT.kAzure + 1, ROOT.kOrange + 7, ROOT.kGreen + 2]
    s_colors = [ROOT.kAzure - 6, ROOT.kOrange - 3, ROOT.kGreen + 3]
    canvas = ROOT.TCanvas("c_overlay", "", 1120, 820)
    canvas.SetLeftMargin(0.13)
    canvas.SetRightMargin(0.05)
    canvas.SetBottomMargin(0.13)
    canvas.SetTopMargin(0.08)
    frame = canvas.DrawFrame(14.0, 0.0, 36.0, max_y(background, signal) * 1.45)
    frame.GetXaxis().SetTitle("Candidate E_{T} [GeV]")
    frame.GetYaxis().SetTitle("Weighted candidates")
    frame.GetXaxis().SetTitleSize(0.045)
    frame.GetYaxis().SetTitleSize(0.045)
    frame.GetXaxis().SetLabelSize(0.038)
    frame.GetYaxis().SetLabelSize(0.038)
    frame.SetStats(False)

    leg = ROOT.TLegend(0.68, 0.56, 0.94, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.029)
    keepalive = []
    for idx, label in enumerate(background.keys()):
        gb = graph_from_rows(background[label], b_colors[idx], 20 + idx)
        gs = graph_from_rows(signal[label], s_colors[idx], 24 + idx)
        keepalive.extend([gb, gs])
        gb.Draw("P SAME")
        gs.Draw("P SAME")
        leg.AddEntry(gb, f"bkg {label}", "pe")
        leg.AddEntry(gs, f"signal {label}", "pe")
    leg.Draw()
    draw_label("Photon12+20 signal vs Jet12+20+30 background", "ABCD-summed; fixed isolation; stitched weights")
    canvas.SaveAs(str(out / f"signal_background_abcdsum_overlay_cent3_{iso_tag}.png"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--background-root", required=True, type=Path)
    parser.add_argument("--signal-root", type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--iso-tag", default="isoR30_fixedIso4GeV")
    args = parser.parse_args()

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetTitleFont(42, "XYZ")
    ROOT.gStyle.SetLabelFont(42, "XYZ")
    ROOT.gStyle.SetTextFont(42)
    ROOT.gStyle.SetEndErrorSize(6)
    ROOT.TH1.SetDefaultSumw2(True)

    args.outdir.mkdir(parents=True, exist_ok=True)
    bfile = ROOT.TFile.Open(str(args.background_root))
    if not bfile or bfile.IsZombie():
        raise OSError(args.background_root)
    background = spectra(bfile, args.iso_tag)

    signal = None
    if args.signal_root:
        sfile = ROOT.TFile.Open(str(args.signal_root))
        if not sfile or sfile.IsZombie():
            raise OSError(args.signal_root)
        signal = spectra(sfile, args.iso_tag)

    write_csv(args.outdir / f"abcdsum_spectrum_cent3_{args.iso_tag}.csv", background, signal)
    draw_background_only(args.outdir, background, args.iso_tag)
    draw_background_fractional_uncertainty(args.outdir, background, args.iso_tag)
    if signal is not None:
        draw_overlay(args.outdir, background, signal, args.iso_tag)


if __name__ == "__main__":
    main()
