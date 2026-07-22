#!/usr/bin/env python3
"""Render the contract-locked THE97 partial pp-data ABCD purity preview.

This intentionally accepts only the canonical combined photon+jet signal
leakage template and the two audited partial-data period ROOTs.  It is a
non-final diagnostic: the final product must be regenerated after raw data
production drains.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import ROOT

from make_ppg12_fig3_purity_sim_sdcc_vs_current_overlay import (
    effective_count,
    leak_ratio,
    ppg12_toy_estimate,
)


ROOT.gROOT.SetBatch(True)

REPO = Path(__file__).resolve().parents[3]
CFG = "jetMinPtScan_dphiScan_vz60_isoR40_isSliding_preselectionNewPPG12_tightNewPPG12_nonTightNewPPG12"
TAG = "the97_pp_data_partial_diagnostic_20260710T160658Z"
CONTRACT_ID = "the97_ppg12_data_purity_leakage_v1"
SOURCE_SHA256 = "c352e5c70a05777a34167f7793d3dce7727813793a9dcf1cdc4ba5c8145f13ad"
EXPECTED_EDGES = (10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36)
DATA_DIR = "PPG12_scaledtrigger30"
SIM_DIR = "SIM"
DATA_NAMES = {
    "A": "h_tight_iso_cluster_0",
    "B": "h_tight_noniso_cluster_0",
    "C": "h_nontight_iso_cluster_0",
    "D": "h_nontight_noniso_cluster_0",
}
SIM_NAMES = {
    "A": "h_tight_iso_cluster_signal_0",
    "B": "h_tight_noniso_cluster_signal_0",
    "C": "h_nontight_iso_cluster_signal_0",
    "D": "h_nontight_noniso_cluster_signal_0",
}
BASE = REPO / "dataOutput/ppg12Parity/the97_ppg12_final_parity_full_20260709_2230" / "partial_pp_data_20260710T160658Z"
DEFAULT_0MRAD = BASE / "remote_roots_0mrad_tar/0mrad/pp" / f"RecoilJets_pp_ALL_{CFG}.root"
DEFAULT_1P5MRAD = BASE / "remote_roots_1p5mrad_tar/1p5mrad/pp" / f"RecoilJets_pp_ALL_{CFG}.root"
DEFAULT_SIM = REPO / "InputFiles/the97_ppg12_final_parity_full_20260709_2230/final_merged_roots/final_combined_canonical_20260710/RecoilJets_photonjet5plus10plus20_si_di_period_combined_MERGED.root"
DEFAULT_REFERENCE = REPO / "dataOutput/ppg12Parity/the76_ppg12_parity_full_20260701_003024/fig29_purity_datathief_audit/ppg12_fig29_purity_sdcc_points.csv"
DEFAULT_OUT = BASE / "purity_contract_v1"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def open_root(path: Path) -> ROOT.TFile:
    root_file = ROOT.TFile.Open(str(path), "READ")
    if not root_file or root_file.IsZombie() or root_file.TestBit(ROOT.TFile.kRecovered):
        raise RuntimeError(f"invalid ROOT input: {path}")
    return root_file


def require_hist(root_file: ROOT.TFile, path: str) -> ROOT.TH1:
    hist = root_file.Get(path)
    if not hist or not hist.InheritsFrom("TH1"):
        raise RuntimeError(f"missing TH1: {path}")
    return hist


def edges(hist: ROOT.TH1) -> tuple[float, ...]:
    axis = hist.GetXaxis()
    return tuple(float(axis.GetBinLowEdge(index)) for index in range(1, hist.GetNbinsX() + 1)) + (
        float(axis.GetBinUpEdge(hist.GetNbinsX())),
    )


def assert_binning(histograms: list[ROOT.TH1]) -> None:
    for hist in histograms:
        if hist.GetNbinsX() != 11 or edges(hist) != EXPECTED_EDGES:
            raise RuntimeError(f"unexpected purity binning in {hist.GetName()}: {edges(hist)}")


def summed_data_hist(first: ROOT.TFile, second: ROOT.TFile, name: str) -> ROOT.TH1:
    initial = require_hist(first, f"{DATA_DIR}/{name}")
    other = require_hist(second, f"{DATA_DIR}/{name}")
    result = initial.Clone(f"partial_combined_{name}")
    result.SetDirectory(0)
    result.Add(other)
    return result


def graph(name: str, xs: list[float], exs: list[float], ys: list[float], eys: list[float], color: int, marker: int) -> ROOT.TGraphErrors:
    result = ROOT.TGraphErrors(len(xs))
    result.SetName(name)
    result.SetMarkerColor(color); result.SetLineColor(color)
    result.SetMarkerStyle(marker); result.SetMarkerSize(1.25); result.SetLineWidth(2)
    for index, (x, ex, y, ey) in enumerate(zip(xs, exs, ys, eys)):
        result.SetPoint(index, x, y); result.SetPointError(index, ex, ey)
    return result


def render(rows: list[dict[str, float]], output: Path) -> None:
    canvas = ROOT.TCanvas("the97_partial_data_purity", "", 930, 880)
    top = ROOT.TPad("top", "", 0.0, 0.32, 1.0, 1.0)
    bottom = ROOT.TPad("bottom", "", 0.0, 0.0, 1.0, 0.32)
    canvas.SetFillColor(ROOT.kWhite); top.SetFillColor(ROOT.kWhite); bottom.SetFillColor(ROOT.kWhite)
    top.SetBottomMargin(0.02); bottom.SetTopMargin(0.02); bottom.SetBottomMargin(0.27)
    top.SetLeftMargin(0.13); bottom.SetLeftMargin(0.13); top.SetRightMargin(0.04); bottom.SetRightMargin(0.04)
    top.Draw(); bottom.Draw()
    xs = [row["center"] for row in rows]; exs = [row["half_width"] for row in rows]
    raw = graph("g_partial_raw", xs, exs, [row["raw"] for row in rows], [row["raw_error"] for row in rows], ROOT.kBlack, 20)
    corrected = graph("g_partial_corrected", xs, exs, [row["corrected"] for row in rows], [row["corrected_error"] for row in rows], ROOT.kBlue + 1, 21)
    ratio = graph("g_partial_correction_ratio", xs, exs, [row["corrected_over_raw"] for row in rows], [row["corrected_over_raw_error"] for row in rows], ROOT.kBlue + 1, 21)
    top.cd()
    frame = ROOT.TH1F("frame_partial_purity", "", 26, 10.0, 36.0)
    frame.SetStats(False); frame.SetFillColor(ROOT.kWhite); frame.GetYaxis().SetRangeUser(0.0, 1.2)
    frame.GetYaxis().SetTitle("Photon purity"); frame.GetYaxis().SetTitleSize(0.060); frame.GetYaxis().SetLabelSize(0.048); frame.GetYaxis().SetTitleOffset(0.90)
    frame.GetXaxis().SetLabelSize(0.0); frame.Draw("axis")
    raw.Draw("P SAME"); corrected.Draw("P SAME")
    legend = ROOT.TLegend(0.52, 0.14, 0.94, 0.31)
    legend.SetBorderSize(0); legend.SetFillStyle(0); legend.SetTextFont(42); legend.SetTextSize(0.037)
    legend.AddEntry(raw, "Raw ABCD", "p"); legend.AddEntry(corrected, "Signal-leakage corrected", "p"); legend.Draw()
    label = ROOT.TLatex(); label.SetNDC(True); label.SetTextFont(42)
    label.SetTextSize(0.043); label.DrawLatex(0.15, 0.88, "#bf{#it{sPHENIX}} Internal")
    label.SetTextSize(0.033); label.DrawLatex(0.15, 0.82, "p+p  #sqrt{s} = 200 GeV   PPG12 scaled trigger 30")
    label.SetTextColor(ROOT.kRed + 1); label.SetTextSize(0.035); label.DrawLatex(0.15, 0.75, "PRELIMINARY PARTIAL-COVERAGE: 43.826% (9,420 / 21,494 files)")
    label.SetTextColor(ROOT.kBlack); label.SetTextSize(0.029); label.DrawLatex(0.15, 0.69, "Both crossing-angle periods; canonical combined photon+jet signal leakage")
    offscale = [row for row in rows if row["raw"] > 1.2]
    if offscale:
        row = offscale[0]
        label.SetTextColor(ROOT.kRed + 1); label.SetTextSize(0.027)
        label.DrawLatex(0.15, 0.63, f"Raw {row['pt_lo']:.0f}-{row['pt_hi']:.0f} GeV toy fit = {row['raw']:.2f} #pm {row['raw_error']:.2f} (off scale)")
        label.SetTextColor(ROOT.kBlack)
    bottom.cd()
    ratio_frame = ROOT.TH1F("frame_partial_ratio", "", 26, 10.0, 36.0)
    ratio_frame.SetStats(False); ratio_frame.SetFillColor(ROOT.kWhite); ratio_frame.GetYaxis().SetRangeUser(0.0, 1.4)
    ratio_frame.GetXaxis().SetTitle("E_{T}^{#gamma,rec} [GeV]"); ratio_frame.GetYaxis().SetTitle("Corrected / raw")
    ratio_frame.GetXaxis().SetTitleSize(0.105); ratio_frame.GetXaxis().SetLabelSize(0.085); ratio_frame.GetYaxis().SetTitleSize(0.085); ratio_frame.GetYaxis().SetLabelSize(0.070)
    ratio_frame.GetYaxis().SetTitleOffset(0.64); ratio_frame.GetXaxis().SetTitleOffset(1.05); ratio_frame.Draw("axis")
    unity = ROOT.TLine(10.0, 1.0, 36.0, 1.0); unity.SetLineStyle(7); unity.SetLineColor(ROOT.kGray + 2); unity.Draw()
    ratio.Draw("P SAME")
    canvas.SaveAs(str(output))


def load_reference(path: Path) -> dict[str, list[dict[str, float]]]:
    reference = {"raw": [], "leakage_corrected": []}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            series = row["series"]
            if series not in reference:
                continue
            reference[series].append({key: float(row[key]) for key in ("x", "purity", "ex_low", "ex_high", "ey_low", "ey_high")})
    if any(len(reference[key]) != 11 for key in reference):
        raise RuntimeError(f"PPG12 reference does not contain 11 points per series: {path}")
    return reference


def overlay_rows(rows: list[dict[str, float]], reference: dict[str, list[dict[str, float]]]) -> list[dict[str, float | str]]:
    result: list[dict[str, float | str]] = []
    for series, current_key, error_key in (("raw", "raw", "raw_error"), ("leakage_corrected", "corrected", "corrected_error")):
        for current, ppg12 in zip(rows, reference[series]):
            value = float(current[current_key]); error = float(current[error_key])
            pvalue = ppg12["purity"]; perror = 0.5 * (ppg12["ey_low"] + ppg12["ey_high"])
            ratio = value / pvalue if pvalue else math.nan
            ratio_error = abs(ratio) * math.hypot(error / value, perror / pvalue) if value and pvalue else math.nan
            result.append({"series": series, "x": current["center"], "current": value, "current_error": error, "ppg12_sdcc": pvalue, "ppg12_error": perror, "current_over_ppg12": ratio, "ratio_error": ratio_error})
    return result


def render_overlay(rows: list[dict[str, float]], reference: dict[str, list[dict[str, float]]], output: Path) -> None:
    canvas = ROOT.TCanvas("the97_partial_data_purity_sdcc_overlay", "", 1120, 980)
    top = ROOT.TPad("overlay_top", "", 0.0, 0.31, 1.0, 1.0)
    bottom = ROOT.TPad("overlay_bottom", "", 0.0, 0.0, 1.0, 0.31)
    canvas.SetFillColor(ROOT.kWhite); top.SetFillColor(ROOT.kWhite); bottom.SetFillColor(ROOT.kWhite)
    top.SetBottomMargin(0.02); bottom.SetTopMargin(0.02); bottom.SetBottomMargin(0.28)
    top.SetLeftMargin(0.11); bottom.SetLeftMargin(0.11); top.SetRightMargin(0.27); bottom.SetRightMargin(0.27)
    top.Draw(); bottom.Draw()
    xs = [row["center"] for row in rows]; exs = [row["half_width"] for row in rows]
    ppg_raw = graph("g_ppg12_raw", xs, exs, [row["purity"] for row in reference["raw"]], [0.5 * (row["ey_low"] + row["ey_high"]) for row in reference["raw"]], ROOT.kGray + 2, 24)
    cur_raw = graph("g_current_raw", xs, exs, [row["raw"] for row in rows], [row["raw_error"] for row in rows], ROOT.kBlack, 20)
    ppg_corr = graph("g_ppg12_corrected", xs, exs, [row["purity"] for row in reference["leakage_corrected"]], [0.5 * (row["ey_low"] + row["ey_high"]) for row in reference["leakage_corrected"]], ROOT.kAzure + 6, 25)
    cur_corr = graph("g_current_corrected", xs, exs, [row["corrected"] for row in rows], [row["corrected_error"] for row in rows], ROOT.kBlue + 1, 21)
    overlay = overlay_rows(rows, reference)
    raw_ratio_rows = [row for row in overlay if row["series"] == "raw"]
    corr_ratio_rows = [row for row in overlay if row["series"] == "leakage_corrected"]
    raw_ratio = graph("g_current_over_ppg12_raw", xs, exs, [float(row["current_over_ppg12"]) for row in raw_ratio_rows], [float(row["ratio_error"]) for row in raw_ratio_rows], ROOT.kBlack, 20)
    corr_ratio = graph("g_current_over_ppg12_corrected", xs, exs, [float(row["current_over_ppg12"]) for row in corr_ratio_rows], [float(row["ratio_error"]) for row in corr_ratio_rows], ROOT.kBlue + 1, 21)
    top.cd()
    frame = ROOT.TH1F("frame_partial_purity_sdcc", "", 26, 10.0, 36.0)
    frame.SetStats(False); frame.SetFillColor(ROOT.kWhite); frame.GetYaxis().SetRangeUser(0.25, 1.18)
    frame.GetYaxis().SetTitle("Purity"); frame.GetYaxis().SetTitleSize(0.060); frame.GetYaxis().SetLabelSize(0.047); frame.GetYaxis().SetTitleOffset(0.88)
    frame.GetXaxis().SetLabelSize(0.0); frame.Draw("axis")
    ppg_corr.Draw("P SAME"); cur_corr.Draw("P SAME"); ppg_raw.Draw("P SAME"); cur_raw.Draw("P SAME")
    label = ROOT.TLatex(); label.SetNDC(True); label.SetTextFont(42)
    label.SetTextSize(0.042); label.DrawLatex(0.15, 0.89, "#bf{#it{sPHENIX}} Internal")
    label.SetTextSize(0.032); label.DrawLatex(0.15, 0.83, "p+p  #sqrt{s} = 200 GeV   |#eta| < 0.7")
    label.DrawLatex(0.15, 0.78, "PPG12 SDCC vs registered July 16 output")
    label.SetTextColor(ROOT.kRed + 1); label.SetTextSize(0.033); label.DrawLatex(0.15, 0.72, "PRELIMINARY PARTIAL-COVERAGE: 43.826% (9,420 / 21,494 files)")
    label.SetTextColor(ROOT.kBlack)
    legend = ROOT.TLegend(0.745, 0.48, 0.985, 0.70)
    legend.SetBorderSize(0); legend.SetFillStyle(0); legend.SetTextFont(42); legend.SetTextSize(0.031)
    legend.AddEntry(ppg_corr, "PPG12 leakage corrected", "p"); legend.AddEntry(cur_corr, "Registered output corrected", "p")
    legend.AddEntry(ppg_raw, "PPG12 raw", "p"); legend.AddEntry(cur_raw, "Registered output raw", "p"); legend.Draw()
    top.RedrawAxis()
    bottom.cd()
    ratio_frame = ROOT.TH1F("frame_partial_sdcc_ratio", "", 26, 10.0, 36.0)
    ratio_frame.SetStats(False); ratio_frame.SetFillColor(ROOT.kWhite); ratio_frame.GetYaxis().SetRangeUser(0.50, 1.25)
    ratio_frame.GetXaxis().SetTitle("E_{T}^{#gamma,rec} [GeV]"); ratio_frame.GetYaxis().SetTitle("Current / PPG12")
    ratio_frame.GetXaxis().SetTitleSize(0.105); ratio_frame.GetXaxis().SetLabelSize(0.082); ratio_frame.GetYaxis().SetTitleSize(0.082); ratio_frame.GetYaxis().SetLabelSize(0.066)
    ratio_frame.GetYaxis().SetTitleOffset(0.67); ratio_frame.GetXaxis().SetTitleOffset(1.05); ratio_frame.Draw("axis")
    unity = ROOT.TLine(10.0, 1.0, 36.0, 1.0); unity.SetLineStyle(7); unity.SetLineColor(ROOT.kGray + 2); unity.Draw()
    raw_ratio.Draw("P SAME"); corr_ratio.Draw("P SAME")
    ratio_legend = ROOT.TLegend(0.745, 0.70, 0.985, 0.94)
    ratio_legend.SetBorderSize(0); ratio_legend.SetFillStyle(0); ratio_legend.SetTextFont(42); ratio_legend.SetTextSize(0.062)
    ratio_legend.AddEntry(corr_ratio, "leakage corrected", "p"); ratio_legend.AddEntry(raw_ratio, "raw", "p"); ratio_legend.Draw()
    bottom.RedrawAxis(); canvas.SaveAs(str(output))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-0mrad", type=Path, default=DEFAULT_0MRAD)
    parser.add_argument("--data-1p5mrad", type=Path, default=DEFAULT_1P5MRAD)
    parser.add_argument("--signal-root", type=Path, default=DEFAULT_SIM)
    parser.add_argument("--reference-csv", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    if sha256(args.signal_root) != SOURCE_SHA256:
        raise RuntimeError("canonical photon+jet leakage ROOT SHA-256 mismatch")
    data_zero, data_one, signal = open_root(args.data_0mrad), open_root(args.data_1p5mrad), open_root(args.signal_root)
    data = {key: summed_data_hist(data_zero, data_one, name) for key, name in DATA_NAMES.items()}
    signal_hists = {key: require_hist(signal, f"{SIM_DIR}/{name}") for key, name in SIM_NAMES.items()}
    assert_binning(list(data.values()) + list(signal_hists.values()))
    leakage = {key: leak_ratio(signal_hists[key], signal_hists["A"], f"leak_{key}") for key in ("B", "C", "D")}
    rng = ROOT.TRandom3(20260710)
    rows: list[dict[str, float]] = []
    for index in range(1, 12):
        values = tuple(float(data[key].GetBinContent(index)) for key in ("A", "B", "C", "D"))
        n_eff = tuple(effective_count(data[key], index) for key in ("A", "B", "C", "D"))
        leak = tuple(float(leakage[key].GetBinContent(index)) for key in ("B", "C", "D"))
        leak_error = tuple(float(leakage[key].GetBinError(index)) for key in ("B", "C", "D"))
        raw, raw_error, corrected, corrected_error, toys = ppg12_toy_estimate(rng, values, n_eff, leak, leak_error, f"partial_{index}")
        center = 0.5 * (EXPECTED_EDGES[index - 1] + EXPECTED_EDGES[index])
        ratio = corrected / raw if raw else math.nan
        ratio_error = abs(ratio) * math.hypot(corrected_error / corrected, raw_error / raw) if corrected and raw else math.nan
        rows.append({"pt_lo": EXPECTED_EDGES[index - 1], "pt_hi": EXPECTED_EDGES[index], "center": center, "half_width": 0.5 * (EXPECTED_EDGES[index] - EXPECTED_EDGES[index - 1]), "A": values[0], "B": values[1], "C": values[2], "D": values[3], "n_eff_A": n_eff[0], "n_eff_B": n_eff[1], "n_eff_C": n_eff[2], "n_eff_D": n_eff[3], "cB": leak[0], "cC": leak[1], "cD": leak[2], "cB_error": leak_error[0], "cC_error": leak_error[1], "cD_error": leak_error[2], "raw": raw, "raw_error": raw_error, "corrected": corrected, "corrected_error": corrected_error, "corrected_over_raw": ratio, "corrected_over_raw_error": ratio_error, **toys})
    args.outdir.mkdir(parents=True, exist_ok=True)
    png = args.outdir / "the97_pp_data_purity_PRELIMINARY_PARTIAL_COVERAGE.png"
    overlay_png = args.outdir / "the97_pp_data_purity_sdcc_vs_current_PRELIMINARY_PARTIAL_COVERAGE_overlay_ratio.png"
    csv_path = args.outdir / "the97_pp_data_purity_PRELIMINARY_PARTIAL_COVERAGE_points.csv"
    overlay_csv = args.outdir / "the97_pp_data_purity_sdcc_vs_current_PRELIMINARY_PARTIAL_COVERAGE_points.csv"
    manifest = args.outdir / "the97_pp_data_purity_PRELIMINARY_PARTIAL_COVERAGE_manifest.json"
    render(rows, png)
    reference = load_reference(args.reference_csv)
    render_overlay(rows, reference, overlay_png)
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    overlay_points = overlay_rows(rows, reference)
    with overlay_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(overlay_points[0])); writer.writeheader(); writer.writerows(overlay_points)
    manifest.write_text(json.dumps({"contract_id": CONTRACT_ID, "status": "PRELIMINARY PARTIAL-COVERAGE", "coverage": {"0mrad_files": "4112/10747", "1p5mrad_files": "5308/10747", "combined_files": "9420/21494", "combined_percent": 43.826}, "data_roots": {"0mrad": str(args.data_0mrad), "1p5mrad": str(args.data_1p5mrad), "0mrad_sha256": sha256(args.data_0mrad), "1p5mrad_sha256": sha256(args.data_1p5mrad)}, "signal_leakage_root": str(args.signal_root), "signal_leakage_sha256": SOURCE_SHA256, "data_namespace": DATA_DIR, "data_histograms": DATA_NAMES, "signal_histograms": SIM_NAMES, "bin_edges_gev": EXPECTED_EDGES, "solver": "PPG12 quadratic CalculatePhotonYield closure; 20000 effective-Poisson data toys plus Gaussian cB/cC/cD throws and PPG12 Gaussian-fit estimator", "ppg12_sdcc_reference_csv": str(args.reference_csv), "ppg12_sdcc_reference_source": "/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_bdt_nom.root -> gpurity/gpurity_leak", "ratio_panel": "THE97 current / PPG12 SDCC; uncertainty propagates both current toy-fit and PPG12 graph errors", "offscale_raw_bin_note": "The 32-36 GeV raw partial-data toy fit and current/PPG12 ratio are retained in CSV and explicitly marked off scale in the overlay PNG.", "forbidden_substitutes_used": False, "png": str(png), "overlay_png": str(overlay_png), "points_csv": str(csv_path), "overlay_points_csv": str(overlay_csv)}, indent=2) + "\n")
    print(png); print(overlay_png); print(csv_path); print(overlay_csv); print(manifest)


if __name__ == "__main__":
    main()
