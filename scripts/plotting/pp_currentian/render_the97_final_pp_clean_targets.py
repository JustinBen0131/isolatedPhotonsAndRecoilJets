#!/usr/bin/env python3
"""Render clean full-coverage PPG12 purity and raw-ABCD pp-data comparisons.

This is a local rendering consumer of the two audited period ROOTs.  It keeps
the existing partial diagnostic intact and reuses the locked PPG12 quadratic
20,000-toy purity estimator and SHA-pinned photon+jet leakage template.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
from pathlib import Path

import ROOT

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import make_the97_pp_data_partial_purity_contract_v1 as purity
import make_the97_pp_data_partial_abcd_yield_overlay as abcd
import make_the97_pp_data_partial_fig27_sideband_ratio_overlay as sideband

ROOT.gROOT.SetBatch(True)

REPO = Path(__file__).resolve().parents[3]
TAG = "the97_pp_data_full_coverage_20260713T1635"
OUT = REPO / "dataOutput/ppg12Parity/the97_ppg12_final_parity_full_20260709_2230/final_pp_data_full_20260713T1635"
ROOT_NAME = f"RecoilJets_pp_ALL_{purity.CFG}.root"
DATA_0 = OUT / "remote_roots_0mrad" / ROOT_NAME
DATA_1 = OUT / "remote_roots_1p5mrad" / ROOT_NAME
COVERAGE = {"0mrad_files": "10747/10747", "1p5mrad_files": "10747/10747", "combined_files": "21494/21494", "combined_percent": 100.0}
PPG12_TOY_SEED = 42


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def graph(name: str, xs, exs, ys, eys, color: int, marker: int, size: float = 1.15) -> ROOT.TGraphErrors:
    result = ROOT.TGraphErrors(len(xs))
    result.SetName(name)
    result.SetMarkerColor(color); result.SetLineColor(color)
    result.SetMarkerStyle(marker); result.SetMarkerSize(size); result.SetLineWidth(2)
    for i, (x, ex, y, ey) in enumerate(zip(xs, exs, ys, eys)):
        result.SetPoint(i, x, y); result.SetPointError(i, ex, ey)
    return result


def shared_purity_rows(data_zero: ROOT.TFile, data_one: ROOT.TFile, signal: ROOT.TFile) -> list[dict[str, float]]:
    data = {key: purity.summed_data_hist(data_zero, data_one, name) for key, name in purity.DATA_NAMES.items()}
    signal_hists = {key: purity.require_hist(signal, f"{purity.SIM_DIR}/{name}") for key, name in purity.SIM_NAMES.items()}
    purity.assert_binning(list(data.values()) + list(signal_hists.values()))
    leakage = {key: purity.leak_ratio(signal_hists[key], signal_hists["A"], f"final_leak_{key}") for key in ("B", "C", "D")}
    # Match efficiencytool/CalculatePhotonYield.C exactly.  The 20,000-toy
    # estimator is sequential across ET bins, so a date-derived seed changes
    # the audience-facing central values and errors even for identical inputs.
    rng = ROOT.TRandom3(PPG12_TOY_SEED)
    rows = []
    for index in range(1, 12):
        values = tuple(float(data[key].GetBinContent(index)) for key in ("A", "B", "C", "D"))
        neff = tuple(purity.effective_count(data[key], index) for key in ("A", "B", "C", "D"))
        leak = tuple(float(leakage[key].GetBinContent(index)) for key in ("B", "C", "D"))
        leak_err = tuple(float(leakage[key].GetBinError(index)) for key in ("B", "C", "D"))
        raw, raw_err, corrected, corrected_err, toys = purity.ppg12_toy_estimate(rng, values, neff, leak, leak_err, f"final_{index}")
        lo, hi = purity.EXPECTED_EDGES[index - 1], purity.EXPECTED_EDGES[index]
        rows.append({"pt_lo": lo, "pt_hi": hi, "center": 0.5 * (lo + hi), "half_width": 0.5 * (hi - lo), "A": values[0], "B": values[1], "C": values[2], "D": values[3], "n_eff_A": neff[0], "n_eff_B": neff[1], "n_eff_C": neff[2], "n_eff_D": neff[3], "cB": leak[0], "cC": leak[1], "cD": leak[2], "cB_error": leak_err[0], "cC_error": leak_err[1], "cD_error": leak_err[2], "raw": raw, "raw_error": raw_err, "corrected": corrected, "corrected_error": corrected_err, **toys})
    return rows


def draw_purity(rows: list[dict[str, float]], reference: dict, output: Path) -> list[dict]:
    canvas = ROOT.TCanvas("the97_final_purity", "", 1120, 980)
    top, bottom = ROOT.TPad("purity_top", "", 0, .31, 1, 1), ROOT.TPad("purity_bottom", "", 0, 0, 1, .31)
    for pad in (canvas, top, bottom): pad.SetFillColor(ROOT.kWhite)
    top.SetBottomMargin(.02); bottom.SetTopMargin(.02); bottom.SetBottomMargin(.28)
    top.SetLeftMargin(.11); bottom.SetLeftMargin(.11); top.SetRightMargin(.27); bottom.SetRightMargin(.27)
    top.Draw(); bottom.Draw()
    xs, exs = [r["center"] for r in rows], [r["half_width"] for r in rows]
    p_raw = graph("ppg12_raw", xs, exs, [r["purity"] for r in reference["raw"]], [.5 * (r["ey_low"] + r["ey_high"]) for r in reference["raw"]], ROOT.kGray + 2, 24)
    c_raw = graph("current_raw", xs, exs, [r["raw"] for r in rows], [r["raw_error"] for r in rows], ROOT.kBlack, 20)
    p_cor = graph("ppg12_corrected", xs, exs, [r["purity"] for r in reference["leakage_corrected"]], [.5 * (r["ey_low"] + r["ey_high"]) for r in reference["leakage_corrected"]], ROOT.kAzure + 6, 25)
    c_cor = graph("current_corrected", xs, exs, [r["corrected"] for r in rows], [r["corrected_error"] for r in rows], ROOT.kBlue + 1, 21)
    top.cd()
    frame = ROOT.TH1F("final_purity_frame", "", 26, 10, 36); frame.SetStats(False); frame.SetFillColor(ROOT.kWhite)
    frame.GetYaxis().SetRangeUser(.25, 1.25); frame.GetYaxis().SetTitle("Purity"); frame.GetYaxis().SetTitleSize(.060); frame.GetYaxis().SetLabelSize(.047); frame.GetYaxis().SetTitleOffset(.88); frame.GetXaxis().SetLabelSize(0); frame.Draw("axis")
    for item in (p_cor, c_cor, p_raw, c_raw): item.Draw("P SAME")
    text = ROOT.TLatex(); text.SetNDC(True); text.SetTextFont(42); text.SetTextSize(.042); text.DrawLatex(.15, .84, "#bf{#it{sPHENIX}} Internal")
    text.SetTextSize(.032); text.DrawLatex(.15, .78, "p+p  #sqrt{s} = 200 GeV   |#eta| < 0.7")
    text.DrawLatex(.15, .73, "PPG12 SDCC vs Current output   data: PPG12_scaledtrigger30")
    leg = ROOT.TLegend(.745, .48, .985, .70); leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextFont(42); leg.SetTextSize(.031)
    leg.AddEntry(p_cor, "PPG12 leakage corrected", "p"); leg.AddEntry(c_cor, "Current output corrected", "p"); leg.AddEntry(p_raw, "PPG12 raw", "p"); leg.AddEntry(c_raw, "Current output raw", "p"); leg.Draw()
    top.RedrawAxis()
    overlay = purity.overlay_rows(rows, reference)
    raw_rows, cor_rows = [r for r in overlay if r["series"] == "raw"], [r for r in overlay if r["series"] == "leakage_corrected"]
    bottom.cd(); rframe = ROOT.TH1F("final_purity_ratio_frame", "", 26, 10, 36); rframe.SetStats(False); rframe.SetFillColor(ROOT.kWhite); rframe.GetYaxis().SetRangeUser(.50, 1.35)
    rframe.GetXaxis().SetTitle("E_{T}^{#gamma,rec} [GeV]"); rframe.GetYaxis().SetTitle("Current / PPG12"); rframe.GetXaxis().SetTitleSize(.105); rframe.GetXaxis().SetLabelSize(.082); rframe.GetYaxis().SetTitleSize(.082); rframe.GetYaxis().SetLabelSize(.066); rframe.GetYaxis().SetTitleOffset(.67); rframe.GetXaxis().SetTitleOffset(1.05); rframe.Draw("axis")
    line = ROOT.TLine(10, 1, 36, 1); line.SetLineStyle(7); line.SetLineColor(ROOT.kGray + 2); line.Draw()
    raw_ratio = graph("final_raw_ratio", xs, exs, [float(r["current_over_ppg12"]) for r in raw_rows], [float(r["ratio_error"]) for r in raw_rows], ROOT.kBlack, 20)
    corrected_ratio = graph("final_cor_ratio", xs, exs, [float(r["current_over_ppg12"]) for r in cor_rows], [float(r["ratio_error"]) for r in cor_rows], ROOT.kBlue + 1, 21)
    raw_ratio.Draw("P SAME"); corrected_ratio.Draw("P SAME")
    bottom.RedrawAxis(); canvas.SaveAs(str(output))
    return overlay


def draw_abcd(data_zero: ROOT.TFile, data_one: ROOT.TFile, output: Path) -> list[dict]:
    data = {key: purity.summed_data_hist(data_zero, data_one, name) for key, name in purity.DATA_NAMES.items()}
    purity.assert_binning(list(data.values()))
    current, reference = abcd.data_rows(data), abcd.load_reference(abcd.DEFAULT_REFERENCE)
    rows = abcd.combined_rows(current, reference)
    canvas = ROOT.TCanvas("the97_final_abcd", "", 1120, 980); top, bottom = ROOT.TPad("abcd_top", "", 0, .30, 1, 1), ROOT.TPad("abcd_bottom", "", 0, 0, 1, .30)
    for pad in (canvas, top, bottom): pad.SetFillColor(ROOT.kWhite)
    top.SetLogy(True); top.SetBottomMargin(.02); bottom.SetTopMargin(.02); bottom.SetBottomMargin(.28); top.SetLeftMargin(.12); bottom.SetLeftMargin(.12); top.SetRightMargin(.31); bottom.SetRightMargin(.31); top.Draw(); bottom.Draw()
    top.cd(); frame = ROOT.TH1F("final_abcd_frame", "", 26, 10, 36); frame.SetStats(False); frame.SetFillColor(ROOT.kWhite); frame.GetYaxis().SetRangeUser(1, 2.5e5); frame.GetYaxis().SetTitle("Raw ABCD yield"); frame.GetYaxis().SetTitleSize(.058); frame.GetYaxis().SetLabelSize(.046); frame.GetYaxis().SetTitleOffset(.90); frame.GetXaxis().SetLabelSize(0); frame.Draw("axis")
    pgraphs, cgraphs = {}, {}
    for region, _, color, filled, open_marker in abcd.REGIONS:
        pgraphs[region] = abcd.graph(f"p_{region}", reference[region], "ppg12", "ppg12_error", color, open_marker)
        cgraphs[region] = abcd.graph(f"c_{region}", current[region], "current", "current_error", color, filled)
        pgraphs[region].Draw("P SAME"); cgraphs[region].Draw("P SAME")
    text = ROOT.TLatex(); text.SetNDC(True); text.SetTextFont(42); text.SetTextSize(.041); text.DrawLatex(.46, .84, "#bf{#it{sPHENIX}} Internal"); text.SetTextSize(.032); text.DrawLatex(.46, .78, "p+p  #sqrt{s} = 200 GeV"); text.DrawLatex(.46, .73, "|#eta^{#gamma}| < 0.7"); text.DrawLatex(.46, .68, "Data: PPG12_scaledtrigger30")
    leg = ROOT.TLegend(.715, .48, .985, .86); leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextFont(42); leg.SetTextSize(.030); leg.AddEntry(pgraphs["A"], "PPG12 SDCC (open)", "p"); leg.AddEntry(cgraphs["A"], "Current output (filled)", "p")
    for region, label, _, _, _ in abcd.REGIONS: leg.AddEntry(cgraphs[region], label, "p")
    leg.Draw(); top.RedrawAxis()
    bottom.cd(); rframe = ROOT.TH1F("final_abcd_ratio_frame", "", 26, 10, 36); rframe.SetStats(False); rframe.SetFillColor(ROOT.kWhite); rframe.GetYaxis().SetRangeUser(0, 2.2); rframe.GetXaxis().SetTitle("E_{T}^{#gamma,rec} [GeV]"); rframe.GetYaxis().SetTitle("Current / PPG12"); rframe.GetXaxis().SetTitleSize(.105); rframe.GetXaxis().SetLabelSize(.082); rframe.GetYaxis().SetTitleSize(.080); rframe.GetYaxis().SetLabelSize(.064); rframe.GetYaxis().SetTitleOffset(.67); rframe.GetXaxis().SetTitleOffset(1.06); rframe.Draw("axis")
    line = ROOT.TLine(10, 1, 36, 1); line.SetLineColor(ROOT.kGray + 2); line.SetLineStyle(7); line.Draw()
    ratio_graphs = []
    for region, _, color, filled, _ in abcd.REGIONS:
        rr = [r for r in rows if r["region"] == region]
        ratio_graph = abcd.graph(f"r_{region}", rr, "current_over_ppg12", "current_over_ppg12_error", color, filled)
        ratio_graphs.append(ratio_graph)
        ratio_graph.Draw("P SAME")
    bottom.RedrawAxis(); canvas.SaveAs(str(output))
    return rows


def draw_sideband(points: dict[str, list[dict[str, float]]], output: Path) -> None:
    canvas = ROOT.TCanvas("the97_final_sideband", "", 1120, 980)
    top, bottom = ROOT.TPad("sideband_top", "", 0, .31, 1, 1), ROOT.TPad("sideband_bottom", "", 0, 0, 1, .31)
    for pad in (canvas, top, bottom): pad.SetFillColor(ROOT.kWhite)
    top.SetBottomMargin(.02); bottom.SetTopMargin(.02); bottom.SetBottomMargin(.29); top.SetLeftMargin(.12); bottom.SetLeftMargin(.12); top.SetRightMargin(.28); bottom.SetRightMargin(.28); top.Draw(); bottom.Draw()
    current, ppg12, ratios = {}, {}, {}
    for key, _, _, _, color, marker, offset in sideband.SPECS:
        current[key] = sideband.graph(f"final_current_{key}", points[key], "current", "current_error", color, marker, offset)
        ppg12[key] = sideband.graph(f"final_ppg12_{key}", points[key], "ppg12", "ppg12_error", color, 24, offset)
        ratios[key] = sideband.graph(f"final_ratio_{key}", points[key], "current_over_ppg12", "current_over_ppg12_error", color, marker, offset)
    top.cd(); frame = ROOT.TH1F("final_sideband_frame", "", 26, 10, 36); frame.SetStats(False); frame.SetFillColor(ROOT.kWhite); frame.GetYaxis().SetRangeUser(0, 2); frame.GetYaxis().SetTitle("ABCD yield ratio"); frame.GetYaxis().SetTitleSize(.058); frame.GetYaxis().SetTitleOffset(.91); frame.GetYaxis().SetLabelSize(.047); frame.GetXaxis().SetLabelSize(0); frame.Draw("axis")
    for key, *_ in sideband.SPECS: ppg12[key].Draw("P SAME"); current[key].Draw("P SAME")
    text = ROOT.TLatex(); text.SetNDC(True); text.SetTextFont(42); text.SetTextSize(.040); text.DrawLatex(.15, .84, "#bf{#it{sPHENIX}} Internal"); text.SetTextSize(.031); text.DrawLatex(.15, .79, "p+p  #sqrt{s} = 200 GeV"); text.DrawLatex(.15, .74, "|#eta^{#gamma}| < 0.7"); text.DrawLatex(.15, .69, "Data: PPG12_scaledtrigger30")
    source = ROOT.TLegend(.72, .79, .985, .95); source.SetBorderSize(0); source.SetFillStyle(0); source.SetTextFont(42); source.SetTextSize(.031); source.SetHeader("source", "C"); source.AddEntry(ppg12["BoverA"], "PPG12 SDCC", "p"); source.AddEntry(current["BoverA"], "Current output", "p"); source.Draw()
    definitions = ROOT.TLegend(.72, .53, .985, .75); definitions.SetBorderSize(0); definitions.SetFillStyle(0); definitions.SetTextFont(42); definitions.SetTextSize(.030); definitions.SetHeader("ratio definition", "C")
    for key, _, _, label, _, _, _ in sideband.SPECS: definitions.AddEntry(current[key], label, "p")
    definitions.Draw(); top.RedrawAxis()
    minimum, maximum = .45, 1.75
    for item in ratios.values():
        for i in range(item.GetN()):
            value = float(item.GetPointY(i))
            if math.isfinite(value): minimum, maximum = min(minimum, value - item.GetErrorY(i)), max(maximum, value + item.GetErrorY(i))
    span = max(.1, maximum - minimum); minimum = max(0., math.floor((minimum - .08 * span) * 20) / 20); maximum = math.ceil((maximum + .08 * span) * 20) / 20
    bottom.cd(); rframe = ROOT.TH1F("final_sideband_ratio_frame", "", 26, 10, 36); rframe.SetStats(False); rframe.SetFillColor(ROOT.kWhite); rframe.GetYaxis().SetRangeUser(minimum, maximum); rframe.GetYaxis().SetTitle("Current / PPG12"); rframe.GetXaxis().SetTitle("E_{T}^{#gamma,rec} [GeV]"); rframe.GetYaxis().SetTitleSize(.079); rframe.GetYaxis().SetTitleOffset(.68); rframe.GetYaxis().SetLabelSize(.064); rframe.GetXaxis().SetTitleSize(.106); rframe.GetXaxis().SetTitleOffset(1.05); rframe.GetXaxis().SetLabelSize(.082); rframe.Draw("axis")
    line = ROOT.TLine(10, 1, 36, 1); line.SetLineColor(ROOT.kGray + 2); line.SetLineStyle(7); line.Draw()
    for key, *_ in sideband.SPECS: ratios[key].Draw("P SAME")
    bottom.RedrawAxis(); canvas.SaveAs(str(output))


def main() -> None:
    if sha256(purity.DEFAULT_SIM) != purity.SOURCE_SHA256:
        raise RuntimeError("canonical combined photon+jet signal leakage SHA-256 mismatch")
    zero, one, signal = purity.open_root(DATA_0), purity.open_root(DATA_1), purity.open_root(purity.DEFAULT_SIM)
    outdir = OUT / "full_coverage_clean_targets"; outdir.mkdir(parents=True, exist_ok=True)
    rows = shared_purity_rows(zero, one, signal); ref = purity.load_reference(purity.DEFAULT_REFERENCE)
    purity_png = outdir / "the97_pp_data_purity_sdcc_vs_current_FULL_COVERAGE_overlay_ratio.png"
    purity_overlay = draw_purity(rows, ref, purity_png)
    abcd_png = outdir / "the97_pp_data_abcd_yield_sdcc_vs_current_FULL_COVERAGE_overlay_ratio.png"
    abcd_rows = draw_abcd(zero, one, abcd_png)
    with (outdir / "the97_pp_data_purity_FULL_COVERAGE_points.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    with (outdir / "the97_pp_data_purity_FULL_COVERAGE_overlay_points.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(purity_overlay[0])); writer.writeheader(); writer.writerows(purity_overlay)
    with (outdir / "the97_pp_data_abcd_yield_FULL_COVERAGE_points.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(abcd_rows[0])); writer.writeheader(); writer.writerows(abcd_rows)
    sideband_points = sideband.build_points(sideband.load_rows(outdir / "the97_pp_data_abcd_yield_FULL_COVERAGE_points.csv"))
    sideband_png = outdir / "the97_pp_data_fig27_sideband_ratios_sdcc_vs_current_FULL_COVERAGE_overlay_ratio.png"
    draw_sideband(sideband_points, sideband_png)
    (outdir / "manifest.json").write_text(json.dumps({"tag": TAG, "status": "FULL RAW COVERAGE", "coverage": COVERAGE, "contract_id": purity.CONTRACT_ID, "data_roots": {"0mrad": str(DATA_0), "1p5mrad": str(DATA_1), "0mrad_sha256": sha256(DATA_0), "1p5mrad_sha256": sha256(DATA_1)}, "signal_leakage_root": str(purity.DEFAULT_SIM), "signal_leakage_sha256": purity.SOURCE_SHA256, "purity_solver": "PPG12 quadratic CalculatePhotonYield closure; 20000 effective-Poisson data toys and Gaussian leakage throws", "purity_toy_seed": PPG12_TOY_SEED, "abcd_comparison": "raw unnormalized ABCD yields; current/PPG12 lower ratio has propagated independent statistical errors", "pngs": [str(purity_png), str(abcd_png), str(sideband_png)]}, indent=2) + "\n")
    print(purity_png); print(abcd_png); print(sideband_png)


if __name__ == "__main__":
    main()
