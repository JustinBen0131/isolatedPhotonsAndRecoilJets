#!/usr/bin/env python3
"""Render the frozen-cutoff THE97 partial-data raw-ABCD yield comparison.

The result is intentionally a preliminary diagnostic.  It uses only the
summed 0 mrad + 1.5 mrad partial data roots and the SHA-recorded PPG12 Fig. 27
SDCC yield points; neither set of data is normalized, reweighted, or fitted.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import ROOT

from make_the97_pp_data_partial_purity_contract_v1 import (
    CFG,
    DATA_DIR,
    DATA_NAMES,
    DEFAULT_0MRAD,
    DEFAULT_1P5MRAD,
    EXPECTED_EDGES,
    assert_binning,
    open_root,
    summed_data_hist,
)


ROOT.gROOT.SetBatch(True)

REPO = Path(__file__).resolve().parents[3]
TAG = "the97_pp_data_partial_diagnostic_20260710T160658Z"
COVERAGE = {
    "0mrad_files": "4112/10747",
    "1p5mrad_files": "5308/10747",
    "combined_files": "9420/21494",
    "combined_percent": 43.826,
}
BASE = REPO / "dataOutput/ppg12Parity/the97_ppg12_final_parity_full_20260709_2230" / "partial_pp_data_20260710T160658Z"
DEFAULT_REFERENCE = REPO / "dataOutput/ppg12Parity/the93_ppg12_canonical_full_20260706_2145/data_abcd_fig27_ppg12_scaledtrigger30/fig27_abcd_yield_current_vs_ppg12_overlay_ratio.csv"
DEFAULT_OUT = BASE / "data_abcd_fig27_ppg12_scaledtrigger30"

REGIONS = (
    ("A", "A: tight iso", ROOT.kBlack, 20, 24),
    ("B", "B: tight noniso", ROOT.kRed + 1, 21, 25),
    ("C", "C: nontight iso", ROOT.kBlue + 1, 22, 26),
    ("D", "D: nontight noniso", ROOT.kMagenta + 1, 23, 32),
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def graph(name: str, rows: list[dict[str, float]], value: str, error: str, color: int, marker: int) -> ROOT.TGraphErrors:
    result = ROOT.TGraphErrors(len(rows))
    result.SetName(name)
    result.SetMarkerColor(color)
    result.SetLineColor(color)
    result.SetMarkerStyle(marker)
    result.SetMarkerSize(1.05)
    result.SetLineWidth(2)
    for index, row in enumerate(rows):
        result.SetPoint(index, row["center"], row[value])
        result.SetPointError(index, row["half_width"], row[error])
    return result


def load_reference(path: Path) -> dict[str, list[dict[str, float]]]:
    result: dict[str, list[dict[str, float]]] = {key: [] for key, *_ in REGIONS}
    with path.open(newline="") as handle:
        for raw in csv.DictReader(handle):
            region = raw["region"]
            if region not in result:
                continue
            lo, hi = float(raw["pt_lo"]), float(raw["pt_hi"])
            result[region].append({
                "center": 0.5 * (lo + hi),
                "half_width": 0.5 * (hi - lo),
                "ppg12": float(raw["ppg12_yield"]),
                "ppg12_error": float(raw["ppg12_error"]),
            })
    for region, values in result.items():
        if len(values) != 11:
            raise RuntimeError(f"PPG12 reference has {len(values)} rather than 11 {region} bins")
        got_edges = tuple(values[index]["center"] - values[index]["half_width"] for index in range(11)) + (values[-1]["center"] + values[-1]["half_width"],)
        if got_edges != EXPECTED_EDGES:
            raise RuntimeError(f"PPG12 reference binning mismatch for {region}: {got_edges}")
    return result


def data_rows(data: dict[str, ROOT.TH1]) -> dict[str, list[dict[str, float]]]:
    rows: dict[str, list[dict[str, float]]] = {}
    for region, hist in data.items():
        rows[region] = []
        for index in range(1, hist.GetNbinsX() + 1):
            lo, hi = EXPECTED_EDGES[index - 1], EXPECTED_EDGES[index]
            rows[region].append({
                "center": 0.5 * (lo + hi),
                "half_width": 0.5 * (hi - lo),
                "current": float(hist.GetBinContent(index)),
                "current_error": float(hist.GetBinError(index)),
            })
    return rows


def combined_rows(current: dict[str, list[dict[str, float]]], reference: dict[str, list[dict[str, float]]]) -> list[dict[str, float | str]]:
    result: list[dict[str, float | str]] = []
    for region, label, *_ in REGIONS:
        for index, (cur, ref) in enumerate(zip(current[region], reference[region])):
            ratio = cur["current"] / ref["ppg12"] if ref["ppg12"] else math.nan
            ratio_error = abs(ratio) * math.hypot(
                cur["current_error"] / cur["current"] if cur["current"] else 0.0,
                ref["ppg12_error"] / ref["ppg12"] if ref["ppg12"] else 0.0,
            )
            result.append({
                "region": region,
                "label": label,
                "pt_lo": EXPECTED_EDGES[index],
                "pt_hi": EXPECTED_EDGES[index + 1],
                **cur,
                **ref,
                "current_over_ppg12": ratio,
                "current_over_ppg12_error": ratio_error,
            })
    return result


def render(current: dict[str, list[dict[str, float]]], reference: dict[str, list[dict[str, float]]], rows: list[dict[str, float | str]], output: Path) -> None:
    canvas = ROOT.TCanvas("the97_partial_data_abcd_yields", "", 1120, 980)
    top = ROOT.TPad("top", "", 0.0, 0.30, 1.0, 1.0)
    bottom = ROOT.TPad("bottom", "", 0.0, 0.0, 1.0, 0.30)
    for pad in (canvas, top, bottom):
        pad.SetFillColor(ROOT.kWhite)
    top.SetLogy(True)
    top.SetBottomMargin(0.02); bottom.SetTopMargin(0.02); bottom.SetBottomMargin(0.28)
    top.SetLeftMargin(0.12); bottom.SetLeftMargin(0.12); top.SetRightMargin(0.31); bottom.SetRightMargin(0.31)
    top.Draw(); bottom.Draw()

    top.cd()
    frame = ROOT.TH1F("abcd_yield_frame", "", 26, 10.0, 36.0)
    frame.SetStats(False); frame.SetFillColor(ROOT.kWhite)
    frame.GetYaxis().SetRangeUser(1.0, 2.5e5)
    frame.GetYaxis().SetTitle("Raw ABCD yield")
    frame.GetYaxis().SetTitleSize(0.058); frame.GetYaxis().SetLabelSize(0.046); frame.GetYaxis().SetTitleOffset(0.90)
    frame.GetXaxis().SetLabelSize(0.0); frame.Draw("axis")

    ppg_graphs: dict[str, ROOT.TGraphErrors] = {}
    cur_graphs: dict[str, ROOT.TGraphErrors] = {}
    for region, _, color, filled, open_marker in REGIONS:
        ppg_graphs[region] = graph(f"g_ppg12_{region}", reference[region], "ppg12", "ppg12_error", color, open_marker)
        cur_graphs[region] = graph(f"g_current_{region}", current[region], "current", "current_error", color, filled)
        ppg_graphs[region].Draw("P SAME")
        cur_graphs[region].Draw("P SAME")

    label = ROOT.TLatex(); label.SetNDC(True); label.SetTextFont(42)
    label.SetTextSize(0.041); label.DrawLatex(0.46, 0.89, "#bf{#it{sPHENIX}} Internal")
    label.SetTextSize(0.032); label.DrawLatex(0.46, 0.83, "p+p  #sqrt{s} = 200 GeV")
    label.DrawLatex(0.46, 0.78, "|#eta^{#gamma}| < 0.7")
    label.DrawLatex(0.46, 0.73, "Data: PPG12_scaledtrigger30")
    # Deliberately drawn in the lower-left low-y white space, away from the
    # low-ET high-y data points and high-ET low-y data points.
    label.SetTextColor(ROOT.kRed + 1); label.SetTextSize(0.029)
    label.DrawLatex(0.16, 0.145, "PRELIMINARY PARTIAL-COVERAGE: 43.826%")
    label.SetTextColor(ROOT.kBlack)

    legend = ROOT.TLegend(0.715, 0.48, 0.985, 0.86)
    legend.SetBorderSize(0); legend.SetFillStyle(0); legend.SetTextFont(42); legend.SetTextSize(0.030)
    legend.AddEntry(ppg_graphs["A"], "PPG12 SDCC (open)", "p")
    legend.AddEntry(cur_graphs["A"], "Current output (filled)", "p")
    for region, region_label, _, _, _ in REGIONS:
        legend.AddEntry(cur_graphs[region], region_label, "p")
    legend.Draw(); top.RedrawAxis()

    bottom.cd()
    ratio_frame = ROOT.TH1F("abcd_yield_ratio_frame", "", 26, 10.0, 36.0)
    ratio_frame.SetStats(False); ratio_frame.SetFillColor(ROOT.kWhite)
    ratio_frame.GetYaxis().SetRangeUser(0.0, 2.2)
    ratio_frame.GetXaxis().SetTitle("E_{T}^{#gamma,rec} [GeV]")
    ratio_frame.GetYaxis().SetTitle("Current / PPG12")
    ratio_frame.GetXaxis().SetTitleSize(0.105); ratio_frame.GetXaxis().SetLabelSize(0.082)
    ratio_frame.GetYaxis().SetTitleSize(0.080); ratio_frame.GetYaxis().SetLabelSize(0.064)
    ratio_frame.GetYaxis().SetTitleOffset(0.67); ratio_frame.GetXaxis().SetTitleOffset(1.06)
    ratio_frame.Draw("axis")
    unity = ROOT.TLine(10.0, 1.0, 36.0, 1.0)
    unity.SetLineColor(ROOT.kGray + 2); unity.SetLineStyle(7); unity.Draw()
    ratio_graphs: list[ROOT.TGraphErrors] = []
    for region, _, color, filled, _ in REGIONS:
        region_rows = [row for row in rows if row["region"] == region]
        ratio_graph = graph(f"g_ratio_{region}", region_rows, "current_over_ppg12", "current_over_ppg12_error", color, filled)
        ratio_graphs.append(ratio_graph)
        ratio_graph.Draw("P SAME")
    bottom.RedrawAxis()
    canvas.SaveAs(str(output))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-0mrad", type=Path, default=DEFAULT_0MRAD)
    parser.add_argument("--data-1p5mrad", type=Path, default=DEFAULT_1P5MRAD)
    parser.add_argument("--reference-csv", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    first, second = open_root(args.data_0mrad), open_root(args.data_1p5mrad)
    data = {key: summed_data_hist(first, second, name) for key, name in DATA_NAMES.items()}
    assert_binning(list(data.values()))
    current = data_rows(data)
    reference = load_reference(args.reference_csv)
    rows = combined_rows(current, reference)
    args.outdir.mkdir(parents=True, exist_ok=True)
    png = args.outdir / "the97_pp_data_abcd_yield_sdcc_vs_current_PRELIMINARY_PARTIAL_COVERAGE_overlay_ratio.png"
    csv_path = args.outdir / "the97_pp_data_abcd_yield_sdcc_vs_current_PRELIMINARY_PARTIAL_COVERAGE_points.csv"
    manifest = args.outdir / "the97_pp_data_abcd_yield_sdcc_vs_current_PRELIMINARY_PARTIAL_COVERAGE_manifest.json"
    render(current, reference, rows, png)
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    manifest.write_text(json.dumps({
        "schema": "THE97_PARTIAL_DATA_ABCD_YIELD_SDCC_OVERLAY_V1",
        "tag": TAG,
        "status": "PRELIMINARY PARTIAL-COVERAGE",
        "coverage": COVERAGE,
        "data_roots": {"0mrad": str(args.data_0mrad), "1p5mrad": str(args.data_1p5mrad), "0mrad_sha256": sha256(args.data_0mrad), "1p5mrad_sha256": sha256(args.data_1p5mrad)},
        "data_namespace": DATA_DIR,
        "histograms": DATA_NAMES,
        "ppg12_sdcc_reference_csv": str(args.reference_csv),
        "ppg12_sdcc_reference_csv_sha256": sha256(args.reference_csv),
        "ppg12_sdcc_reference_root": "/sphenix/user/shuhangli/ppg12/efficiencytool/results/data_histo_bdt_nom.root",
        "comparison": "raw ABCD bin yields; PPG12 open markers and THE97 partial-current filled markers; no normalization or scale factor",
        "ratio_panel": "THE97 current / PPG12 SDCC, propagated independent statistical errors",
        "bin_edges_gev": EXPECTED_EDGES,
        "png": str(png), "points_csv": str(csv_path),
    }, indent=2) + "\n")
    print(png); print(csv_path); print(manifest)


if __name__ == "__main__":
    main()
