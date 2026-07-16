#!/usr/bin/env python3
"""Render a Fig. 27-style sideband-ratio overlay from the existing partial merge.

This consumes the audited raw-ABCD points already rendered from the combined
0mrad+1p5mrad partial pp-data roots.  It deliberately performs no new merge,
normalization, weighting, or remote action.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import ROOT


ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

REPO = Path(__file__).resolve().parents[3]
BASE = REPO / "dataOutput/ppg12Parity/the97_ppg12_final_parity_full_20260709_2230/partial_pp_data_20260710T160658Z"
DEFAULT_INPUT = BASE / "data_abcd_fig27_ppg12_scaledtrigger30/the97_pp_data_abcd_yield_sdcc_vs_current_PRELIMINARY_PARTIAL_COVERAGE_points.csv"
DEFAULT_OUT = BASE / "data_abcd_fig27_ppg12_scaledtrigger30"
COVERAGE = "43.826% (9,420 / 21,494 files)"

SPECS = (
    ("BoverA", "B", "A", "B/A: tight noniso", ROOT.kBlack, 20, -0.18),
    ("CoverA", "C", "A", "C/A: nontight iso", ROOT.kRed + 1, 21, 0.00),
    ("DoverA", "D", "A", "D/A: nontight noniso", ROOT.kBlue + 1, 22, 0.18),
)


def load_rows(path: Path) -> dict[str, list[dict[str, float]]]:
    regions = {"A": [], "B": [], "C": [], "D": []}
    with path.open(newline="") as handle:
        for raw in csv.DictReader(handle):
            region = raw["region"]
            if region not in regions:
                raise RuntimeError(f"unexpected region {region!r}")
            regions[region].append({key: float(raw[key]) for key in (
                "pt_lo", "pt_hi", "center", "half_width", "current", "current_error", "ppg12", "ppg12_error",
            )})
    for region, rows in regions.items():
        if len(rows) != 11:
            raise RuntimeError(f"{region}: expected 11 bins, got {len(rows)}")
        rows.sort(key=lambda item: item["pt_lo"])
    return regions


def ratio(num: dict[str, float], den: dict[str, float], prefix: str) -> dict[str, float]:
    value = num[prefix] / den[prefix] if den[prefix] else math.nan
    error_key = f"{prefix}_error"
    error = abs(value) * math.hypot(
        num[error_key] / num[prefix] if num[prefix] else 0.0,
        den[error_key] / den[prefix] if den[prefix] else 0.0,
    )
    return {"value": value, "error": error}


def build_points(regions: dict[str, list[dict[str, float]]]) -> dict[str, list[dict[str, float]]]:
    result: dict[str, list[dict[str, float]]] = {}
    for key, numerator, denominator, *_ in SPECS:
        points: list[dict[str, float]] = []
        for num, den in zip(regions[numerator], regions[denominator]):
            if (num["pt_lo"], num["pt_hi"]) != (den["pt_lo"], den["pt_hi"]):
                raise RuntimeError(f"{key}: numerator/denominator binning mismatch")
            current = ratio(num, den, "current")
            ppg12 = ratio(num, den, "ppg12")
            current_over_ppg12 = current["value"] / ppg12["value"] if ppg12["value"] else math.nan
            current_over_ppg12_error = abs(current_over_ppg12) * math.hypot(
                current["error"] / current["value"] if current["value"] else 0.0,
                ppg12["error"] / ppg12["value"] if ppg12["value"] else 0.0,
            )
            points.append({
                "pt_lo": num["pt_lo"], "pt_hi": num["pt_hi"], "center": num["center"], "half_width": num["half_width"],
                "current": current["value"], "current_error": current["error"],
                "ppg12": ppg12["value"], "ppg12_error": ppg12["error"],
                "current_over_ppg12": current_over_ppg12,
                "current_over_ppg12_error": current_over_ppg12_error,
            })
        result[key] = points
    return result


def graph(name: str, rows: list[dict[str, float]], value: str, error: str, color: int, marker: int, offset: float) -> ROOT.TGraphErrors:
    result = ROOT.TGraphErrors(len(rows))
    result.SetName(name)
    result.SetLineColor(color); result.SetMarkerColor(color)
    result.SetMarkerStyle(marker); result.SetMarkerSize(1.05); result.SetLineWidth(2)
    for index, row in enumerate(rows):
        result.SetPoint(index, row["center"] + offset * (2.0 * row["half_width"]), row[value])
        result.SetPointError(index, 0.0, row[error])
    return result


def render(points: dict[str, list[dict[str, float]]], output: Path) -> None:
    canvas = ROOT.TCanvas("the97_partial_fig27_sideband", "", 1120, 980)
    top = ROOT.TPad("fig27_top", "", 0.0, 0.31, 1.0, 1.0)
    bottom = ROOT.TPad("fig27_bottom", "", 0.0, 0.0, 1.0, 0.31)
    for pad in (canvas, top, bottom):
        pad.SetFillColor(ROOT.kWhite)
    top.SetBottomMargin(0.02); bottom.SetTopMargin(0.02); bottom.SetBottomMargin(0.29)
    top.SetLeftMargin(0.12); bottom.SetLeftMargin(0.12); top.SetRightMargin(0.28); bottom.SetRightMargin(0.28)
    top.Draw(); bottom.Draw()

    current_graphs: dict[str, ROOT.TGraphErrors] = {}
    ppg12_graphs: dict[str, ROOT.TGraphErrors] = {}
    ratio_graphs: dict[str, ROOT.TGraphErrors] = {}
    for key, _, _, _, color, marker, offset in SPECS:
        current_graphs[key] = graph(f"g_current_{key}", points[key], "current", "current_error", color, marker, offset)
        ppg12_graphs[key] = graph(f"g_ppg12_{key}", points[key], "ppg12", "ppg12_error", color, 24, offset)
        ratio_graphs[key] = graph(f"g_current_over_ppg12_{key}", points[key], "current_over_ppg12", "current_over_ppg12_error", color, marker, offset)

    top.cd()
    frame = ROOT.TH1F("fig27_sideband_top_frame", "", 26, 10.0, 36.0)
    frame.SetStats(False); frame.SetFillColor(ROOT.kWhite)
    frame.GetYaxis().SetRangeUser(0.0, 2.0)
    frame.GetYaxis().SetTitle("ABCD yield ratio")
    frame.GetYaxis().SetTitleSize(0.058); frame.GetYaxis().SetTitleOffset(0.91); frame.GetYaxis().SetLabelSize(0.047)
    frame.GetXaxis().SetLabelSize(0.0); frame.Draw("axis")
    for key, *_ in SPECS:
        ppg12_graphs[key].Draw("P SAME")
        current_graphs[key].Draw("P SAME")
    label = ROOT.TLatex(); label.SetNDC(True); label.SetTextFont(42)
    label.SetTextSize(0.040); label.DrawLatex(0.15, 0.90, "#bf{#it{sPHENIX}} Internal")
    label.SetTextSize(0.031); label.DrawLatex(0.15, 0.845, "p+p  #sqrt{s} = 200 GeV")
    label.DrawLatex(0.15, 0.790, "|#eta^{#gamma}| < 0.7")
    label.DrawLatex(0.15, 0.735, "Data: PPG12_scaledtrigger30")
    label.SetTextColor(ROOT.kRed + 1); label.SetTextSize(0.027)
    # The right margin is intentionally reserved by the pad geometry.  Keeping
    # the long preliminary marker there avoids all B/A, C/A, and D/A points.
    label.DrawLatex(0.70, 0.435, "PRELIMINARY PARTIAL-COVERAGE")
    label.DrawLatex(0.70, 0.390, COVERAGE)
    label.SetTextColor(ROOT.kBlack)
    source_leg = ROOT.TLegend(0.68, 0.78, 0.98, 0.94)
    source_leg.SetBorderSize(0); source_leg.SetFillStyle(0); source_leg.SetTextFont(42); source_leg.SetTextSize(0.031)
    source_leg.SetHeader("source", "C")
    source_leg.AddEntry(ppg12_graphs["BoverA"], "PPG12 SDCC", "p")
    source_leg.AddEntry(current_graphs["BoverA"], "Current output", "p")
    source_leg.Draw()
    ratio_leg = ROOT.TLegend(0.68, 0.52, 0.985, 0.75)
    ratio_leg.SetBorderSize(0); ratio_leg.SetFillStyle(0); ratio_leg.SetTextFont(42); ratio_leg.SetTextSize(0.030)
    ratio_leg.SetHeader("ratio definition", "C")
    for key, _, _, text, _, _, _ in SPECS:
        ratio_leg.AddEntry(current_graphs[key], text, "p")
    ratio_leg.Draw(); top.RedrawAxis()

    bottom.cd()
    minimum, maximum = 0.45, 1.75
    for graph_item in ratio_graphs.values():
        for index in range(graph_item.GetN()):
            y = float(graph_item.GetPointY(index))
            if math.isfinite(y):
                minimum = min(minimum, y - graph_item.GetErrorY(index))
                maximum = max(maximum, y + graph_item.GetErrorY(index))
    span = max(0.1, maximum - minimum)
    minimum = max(0.0, math.floor((minimum - 0.08 * span) * 20.0) / 20.0)
    maximum = math.ceil((maximum + 0.08 * span) * 20.0) / 20.0
    ratio_frame = ROOT.TH1F("fig27_sideband_bottom_frame", "", 26, 10.0, 36.0)
    ratio_frame.SetStats(False); ratio_frame.SetFillColor(ROOT.kWhite)
    ratio_frame.GetYaxis().SetRangeUser(minimum, maximum)
    ratio_frame.GetYaxis().SetTitle("Current / PPG12")
    ratio_frame.GetXaxis().SetTitle("E_{T}^{#gamma,rec} [GeV]")
    ratio_frame.GetYaxis().SetTitleSize(0.079); ratio_frame.GetYaxis().SetTitleOffset(0.68); ratio_frame.GetYaxis().SetLabelSize(0.064)
    ratio_frame.GetXaxis().SetTitleSize(0.106); ratio_frame.GetXaxis().SetTitleOffset(1.05); ratio_frame.GetXaxis().SetLabelSize(0.082)
    ratio_frame.Draw("axis")
    line = ROOT.TLine(10.0, 1.0, 36.0, 1.0); line.SetLineColor(ROOT.kGray + 2); line.SetLineStyle(7); line.Draw()
    for key, *_ in SPECS:
        ratio_graphs[key].Draw("P SAME")
    bottom.RedrawAxis(); canvas.SaveAs(str(output))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    regions = load_rows(args.input)
    points = build_points(regions)
    args.outdir.mkdir(parents=True, exist_ok=True)
    png = args.outdir / "the97_pp_data_fig27_sideband_ratios_sdcc_vs_current_PRELIMINARY_PARTIAL_COVERAGE_overlay_ratio.png"
    csv_path = args.outdir / "the97_pp_data_fig27_sideband_ratios_sdcc_vs_current_PRELIMINARY_PARTIAL_COVERAGE_points.csv"
    manifest = args.outdir / "the97_pp_data_fig27_sideband_ratios_sdcc_vs_current_PRELIMINARY_PARTIAL_COVERAGE_manifest.json"
    render(points, png)
    rows = []
    for key, numerator, denominator, label, *_ in SPECS:
        for point in points[key]:
            rows.append({"ratio": key, "definition": f"{numerator}/{denominator}", "label": label, **point})
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    manifest.write_text(json.dumps({
        "schema": "THE97_PARTIAL_DATA_FIG27_SIDEBAND_RATIO_OVERLAY_V1",
        "status": "PRELIMINARY PARTIAL-COVERAGE",
        "coverage": COVERAGE,
        "input_raw_abcd_points": str(args.input),
        "input_provenance": "Existing raw-ABCD CSV from the audited combined 0mrad+1p5mrad partial merged ROOTs; no fresh merge or transfer.",
        "data_namespace": "PPG12_scaledtrigger30",
        "ratios": ["B/A: tight noniso / tight iso", "C/A: nontight iso / tight iso", "D/A: nontight noniso / tight iso"],
        "normalization": "none; numerator and denominator are raw yields from the same partial data set in each E_T bin",
        "bottom_panel": "Current output sideband ratio / PPG12 SDCC sideband ratio, with propagated independent statistical errors",
        "png": str(png), "points_csv": str(csv_path),
    }, indent=2) + "\n")
    print(png); print(csv_path); print(manifest)


if __name__ == "__main__":
    main()
