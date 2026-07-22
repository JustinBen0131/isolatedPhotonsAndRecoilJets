#!/usr/bin/env python3
"""Render the matched-size THE-110 purity and Fig. 6 efficiency slide pair."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import struct
import sys
from pathlib import Path

import ROOT


REPO = Path(__file__).resolve().parents[3]
PLOT_DIR = REPO / "scripts/plotting/pp_currentian"
sys.path.insert(0, str(PLOT_DIR))

import make_the76_fig6_efficiency_fixed_rerun_overlay as efficiency  # noqa: E402
import make_the97_current_pp_data_purity_overlay as purity  # noqa: E402


TRUTH = "h_photonEffPpg12Fig6TruthDen_pTgamma_0"
RECO = "h_photonEffPpg12Fig6Reco_pTgamma_0"
RECO_ISO = "h_photonEffPpg12Fig6RecoIso_pTgamma_0"
TIGHT_ISO = "h_photonEffPpg12Fig6RecoTightIso_pTgamma_0"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def png_size(path: Path) -> tuple[int, int]:
    with path.open("rb") as handle:
        header = handle.read(24)
    if header[:8] != b"\x89PNG\r\n\x1a\n":
        raise RuntimeError(f"not a PNG: {path}")
    return struct.unpack(">II", header[16:24])


def weighted_subset(num: float, den: float, den_error2: float) -> tuple[float, float]:
    return efficiency.weighted_subset_efficiency(num, den, math.sqrt(max(0.0, den_error2)))


def build_current_efficiency(
    sums_path: Path,
    ppg12: dict[str, list[dict[str, float]]],
) -> dict[str, list[dict[str, float]]]:
    payload = json.loads(sums_path.read_text())
    if payload.get("input_count") != 6000:
        raise RuntimeError(f"expected 6000 photon ROOTs, got {payload.get('input_count')}")
    mids = [float(value) for value in payload["pt_mid"]]
    expected = [row["pt_mid"] for row in ppg12["reco"]]
    if mids != expected:
        raise RuntimeError(f"efficiency pT centers differ: {mids} != {expected}")
    sums = payload["sums"]
    out: dict[str, list[dict[str, float]]] = {stage: [] for stage, *_ in efficiency.CURVES}
    for index, ref in enumerate(ppg12["reco"]):
        truth, truth_error2 = map(float, sums[TRUTH][index])
        reco, _reco_error2 = map(float, sums[RECO][index])
        reco_iso, reco_iso_error2 = map(float, sums[RECO_ISO][index])
        tight_iso, _tight_iso_error2 = map(float, sums[TIGHT_ISO][index])

        reco_eff, reco_err = weighted_subset(reco, truth, truth_error2)
        id_eff, id_err = weighted_subset(tight_iso, reco_iso, reco_iso_error2)
        all_eff, all_err = weighted_subset(tight_iso, truth, truth_error2)
        product_eff = reco_eff * id_eff
        product_err = math.sqrt((id_eff * reco_err) ** 2 + (reco_eff * id_err) ** 2)

        common = {
            "pt_low": ref["pt_low"],
            "pt_high": ref["pt_high"],
            "pt_mid": ref["pt_mid"],
        }
        out["reco"].append(
            {**common, "eff": reco_eff, "err_low": reco_err, "err_high": reco_err, "num": reco, "den": truth}
        )
        out["reco_id"].append(
            {**common, "eff": product_eff, "err_low": product_err, "err_high": product_err, "num": tight_iso, "den": reco_iso}
        )
        out["reco_id_iso"].append(
            {**common, "eff": all_eff, "err_low": all_err, "err_high": all_err, "num": tight_iso, "den": truth}
        )
    return out


def write_current_csv(path: Path, rows: dict[str, list[dict[str, float]]]) -> None:
    fields = ["stage", "pt_low", "pt_high", "pt_mid", "current_eff", "current_err", "current_num", "current_den"]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for stage, *_ in efficiency.CURVES:
            for row in rows[stage]:
                writer.writerow(
                    {
                        "stage": stage,
                        "pt_low": row["pt_low"],
                        "pt_high": row["pt_high"],
                        "pt_mid": row["pt_mid"],
                        "current_eff": row["eff"],
                        "current_err": row["err_low"],
                        "current_num": row["num"],
                        "current_den": row["den"],
                    }
                )


def make_efficiency_graph(
    name: str,
    rows: list[dict[str, float]],
    color: int,
    marker: int,
    *,
    reference: bool,
) -> ROOT.TGraphAsymmErrors:
    graph = ROOT.TGraphAsymmErrors(len(rows))
    graph.SetName(name)
    graph.SetMarkerStyle(marker)
    graph.SetMarkerSize(1.12 if reference else 1.05)
    graph.SetMarkerColor(color)
    graph.SetLineColor(color)
    graph.SetLineWidth(2)
    for index, row in enumerate(rows):
        graph.SetPoint(index, row["pt_mid"], row["eff"])
        graph.SetPointError(
            index,
            row["pt_mid"] - row["pt_low"],
            row["pt_high"] - row["pt_mid"],
            row["err_low"],
            row["err_high"],
        )
    return graph


def make_efficiency_ratio_graph(
    name: str,
    reference: list[dict[str, float]],
    current: list[dict[str, float]],
    color: int,
) -> ROOT.TGraphAsymmErrors:
    graph = ROOT.TGraphAsymmErrors(len(reference))
    graph.SetName(name)
    graph.SetMarkerStyle(20)
    graph.SetMarkerSize(1.05)
    graph.SetMarkerColor(color)
    graph.SetLineColor(color)
    graph.SetLineWidth(2)
    for index, (ref, cur) in enumerate(zip(reference, current)):
        ratio = cur["eff"] / ref["eff"]
        low = ratio * math.hypot(cur["err_low"] / cur["eff"], ref["err_high"] / ref["eff"])
        high = ratio * math.hypot(cur["err_high"] / cur["eff"], ref["err_low"] / ref["eff"])
        graph.SetPoint(index, ref["pt_mid"], ratio)
        graph.SetPointError(
            index,
            ref["pt_mid"] - ref["pt_low"],
            ref["pt_high"] - ref["pt_mid"],
            low,
            high,
        )
    return graph


def render_efficiency_root(
    output: Path,
    ppg12: dict[str, list[dict[str, float]]],
    current: dict[str, list[dict[str, float]]],
) -> dict[str, dict[str, float]]:
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetEndErrorSize(3)
    colors = {"reco": ROOT.kBlack, "reco_id": ROOT.kMagenta + 1, "reco_id_iso": ROOT.kGreen + 2}
    labels = {
        "reco": "#varepsilon_{reco}",
        "reco_id": "#varepsilon_{reco} #times #varepsilon_{ID}",
        "reco_id_iso": "#varepsilon_{reco} #times #varepsilon_{ID} #times #varepsilon_{iso}",
    }

    canvas = ROOT.TCanvas("the110_bin36_efficiency", "", 860, 900)
    canvas.SetFillColor(ROOT.kWhite)
    top = ROOT.TPad("efficiency_top", "", 0.0, 0.34, 1.0, 1.0)
    bottom = ROOT.TPad("efficiency_ratio", "", 0.0, 0.0, 1.0, 0.34)
    for pad in (top, bottom):
        pad.SetFillColor(ROOT.kWhite)
        pad.SetLeftMargin(0.14)
        pad.SetRightMargin(0.045)
        pad.SetTicks(1, 1)
    top.SetTopMargin(0.045)
    top.SetBottomMargin(0.018)
    bottom.SetTopMargin(0.025)
    bottom.SetBottomMargin(0.27)
    top.Draw()
    bottom.Draw()

    reference_graphs = {
        stage: make_efficiency_graph(f"g_ppg12_{stage}", ppg12[stage], colors[stage], 24, reference=True)
        for stage, *_ in efficiency.CURVES
    }
    current_graphs = {
        stage: make_efficiency_graph(f"g_current_{stage}", current[stage], colors[stage], 20, reference=False)
        for stage, *_ in efficiency.CURVES
    }
    ratio_graphs = {
        stage: make_efficiency_ratio_graph(
            f"g_current_over_ppg12_{stage}", ppg12[stage], current[stage], colors[stage]
        )
        for stage, *_ in efficiency.CURVES
    }

    top.cd()
    frame = ROOT.TH1F("efficiency_frame", "", 26, 10.0, 36.0)
    frame.SetDirectory(0)
    frame.SetStats(False)
    frame.GetYaxis().SetRangeUser(0.0, 1.15)
    frame.GetYaxis().SetTitle("Efficiency")
    frame.GetYaxis().SetTitleSize(0.064)
    frame.GetYaxis().SetTitleOffset(0.88)
    frame.GetYaxis().SetLabelSize(0.048)
    frame.GetYaxis().SetNdivisions(508)
    frame.GetXaxis().SetLabelSize(0.0)
    frame.Draw("AXIS")
    for stage, *_ in efficiency.CURVES:
        reference_graphs[stage].Draw("PZ SAME")
        current_graphs[stage].Draw("PZ SAME")

    label = ROOT.TLatex()
    label.SetNDC(True)
    label.SetTextFont(42)
    label.SetTextColor(ROOT.kBlack)
    # Held off the top frame edge, at larger type.  "Pythia" replaces the
    # collision-system label, so the separate PYTHIA8 line is redundant.
    label.SetTextSize(0.052)
    label.DrawLatex(0.195, 0.880, "#it{#bf{sPHENIX}} Internal")
    label.SetTextSize(0.040)
    label.DrawLatex(0.195, 0.815, "Pythia, #sqrt{s} = 200 GeV, |#eta^{#gamma}| < 0.7")

    # Two columns filling the empty band under the lowest curve.  The fill is
    # opaque, so the top edge must clear the lowest drawn point: reco_id_iso
    # bottoms out at 0.2934, which is NDC y=0.257 on this 0-1.15 frame.
    legend = ROOT.TLegend(0.165, 0.045, 0.945, 0.237)
    legend.SetBorderSize(0)
    legend.SetFillColor(ROOT.kWhite)
    legend.SetFillStyle(1001)
    legend.SetTextFont(42)
    legend.SetTextSize(0.040)
    legend.SetNColumns(2)
    legend.SetColumnSeparation(0.01)
    legend.SetMargin(0.16)
    legend.AddEntry(reference_graphs["reco"], "PPG12 SDCC source", "p")
    legend.AddEntry(current_graphs["reco"], "This analysis output", "p")
    for stage, *_ in efficiency.CURVES:
        legend.AddEntry(current_graphs[stage], labels[stage], "p")
    legend.Draw()
    top.RedrawAxis()

    bottom.cd()
    ratio_frame = ROOT.TH1F("efficiency_ratio_frame", "", 26, 10.0, 36.0)
    ratio_frame.SetDirectory(0)
    ratio_frame.SetStats(False)
    ratio_frame.GetYaxis().SetRangeUser(0.55, 1.45)
    ratio_frame.GetXaxis().SetTitle("E_{T}^{#gamma,truth} [GeV]")
    ratio_frame.GetYaxis().SetTitle("This analysis / PPG12")
    ratio_frame.GetXaxis().SetTitleSize(0.100)
    ratio_frame.GetXaxis().SetTitleOffset(1.10)
    ratio_frame.GetXaxis().SetLabelSize(0.076)
    ratio_frame.GetYaxis().SetTitleSize(0.078)
    ratio_frame.GetYaxis().SetTitleOffset(0.72)
    ratio_frame.GetYaxis().SetLabelSize(0.061)
    ratio_frame.GetYaxis().SetNdivisions(505)
    ratio_frame.Draw("AXIS")
    unity = ROOT.TLine(10.0, 1.0, 36.0, 1.0)
    unity.SetLineColor(ROOT.kGray + 1)
    unity.SetLineStyle(7)
    unity.SetLineWidth(2)
    unity.Draw()
    for stage, *_ in efficiency.CURVES:
        ratio_graphs[stage].Draw("PZ SAME")
    bottom.RedrawAxis()
    canvas.SaveAs(str(output))

    summary: dict[str, dict[str, float]] = {}
    for stage, *_ in efficiency.CURVES:
        ratios = [cur["eff"] / ref["eff"] for ref, cur in zip(ppg12[stage], current[stage])]
        summary[stage] = {
            "min_ratio": min(ratios),
            "max_ratio": max(ratios),
            "mean_ratio": sum(ratios) / len(ratios),
        }
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--photon-root", type=Path, required=True)
    parser.add_argument("--efficiency-sums", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument(
        "--artifact-status",
        default="DIAGNOSTIC_PARTIAL_PP_21492_OF_21494_COMPLETE_PHOTON_SIM",
        help="Truthful status recorded in the output manifest.",
    )
    parser.add_argument(
        "--data-coverage",
        default="21492/21494 pp groups; excludes two pathological parents and their active recovery children",
    )
    parser.add_argument(
        "--rerender-requirement",
        default=(
            "Replace the pp data input and rerender purity after all eight recovery children pass "
            "the exact source-replacement audit; the efficiency panel is already full-stat for this campaign."
        ),
    )
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    purity_png = args.outdir / "the110_bin36_data_purity_ppg12_vs_current.png"
    purity_csv = args.outdir / "the110_bin36_data_purity_ppg12_vs_current_points.csv"
    current_purity = purity.calculate_current(args.data_root, args.photon_root)
    purity_rows = purity.merge_reference(current_purity, purity.load_reference(purity.REFERENCE_CSV))
    purity.write_csv(purity_csv, purity_rows)
    purity_ratio_range = purity.render(purity_rows, purity_png)
    width, height = png_size(purity_png)

    ppg12_eff = efficiency.read_ppg12(efficiency.DEFAULT_PPG12_CSV)
    current_eff = build_current_efficiency(args.efficiency_sums, ppg12_eff)
    current_eff_csv = args.outdir / "the110_bin36_fig6_efficiency_current_compact.csv"
    comparison_csv = args.outdir / "the110_bin36_fig6_efficiency_ppg12_vs_current_points.csv"
    efficiency_png = args.outdir / "the110_bin36_fig6_efficiency_ppg12_vs_current.png"
    write_current_csv(current_eff_csv, current_eff)
    efficiency.write_csv(comparison_csv, ppg12_eff, current_eff)
    efficiency_summary = render_efficiency_root(efficiency_png, ppg12_eff, current_eff)
    if png_size(efficiency_png) != (width, height):
        raise RuntimeError("matched-pair PNG dimensions differ")

    output_files = [purity_png, purity_csv, current_eff_csv, comparison_csv, efficiency_png]
    manifest = {
        "schema": "THE110_BIN36_SLIDE4_MATCHED_PAIR_V1",
        "artifact_status": args.artifact_status,
        "canonical": False,
        "slides_mutated": False,
        "pixel_dimensions": {"width": width, "height": height},
        "purity": {
            "data_root": str(args.data_root),
            "data_root_sha256": sha256(args.data_root),
            "data_coverage": args.data_coverage,
            "signal_leakage_root": str(args.photon_root),
            "signal_leakage_root_sha256": sha256(args.photon_root),
            "signal_leakage_campaign": "the97_ppg12_bin36_final_triple_full_20260717_1322",
            "signal_leakage_coverage": "6000/6000 photonjet raw groups, plain additive target-histogram sum",
            "estimator": "PPG12 quadratic physical branch; ROOT.TRandom3(42); one sequential stream; 20,000 toys per bin",
            "ratio_ylim": purity_ratio_range,
        },
        "efficiency": {
            "raw_sums": str(args.efficiency_sums),
            "raw_sums_sha256": sha256(args.efficiency_sums),
            "campaign": "the97_ppg12_bin36_final_triple_full_20260717_1322",
            "coverage": "6000/6000 photonjet raw ROOTs",
            "input_list_sha256": json.loads(args.efficiency_sums.read_text())["input_list_sha256"],
            "definition": "Exact PPG12 Fig. 6 construction: reco=Reco/TruthDen; reco*ID=(Reco/TruthDen)*(TightIso/RecoIso); reco*ID*iso=TightIso/TruthDen",
            "ratio_summary": efficiency_summary,
        },
        "ppg12_references": {
            "purity_csv": str(purity.REFERENCE_CSV),
            "purity_csv_sha256": sha256(purity.REFERENCE_CSV),
            "efficiency_csv": str(efficiency.DEFAULT_PPG12_CSV),
            "efficiency_csv_sha256": sha256(efficiency.DEFAULT_PPG12_CSV),
        },
        "outputs": {str(path.name): {"path": str(path), "sha256": sha256(path)} for path in output_files},
        "rerender_requirement": args.rerender_requirement,
        "renderer": str(Path(__file__).resolve()),
        "renderer_sha256": sha256(Path(__file__).resolve()),
        "code_dependencies": {
            "purity_renderer": str(Path(purity.__file__).resolve()),
            "purity_renderer_sha256": sha256(Path(purity.__file__).resolve()),
            "efficiency_renderer": str(Path(efficiency.__file__).resolve()),
            "efficiency_renderer_sha256": sha256(Path(efficiency.__file__).resolve()),
        },
        "visual_qa": {"status": "pending_manual_inspection"},
    }
    manifest_path = args.outdir / "the110_bin36_slide4_matched_pair_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"purity_png": str(purity_png), "efficiency_png": str(efficiency_png), "manifest": str(manifest_path)}, indent=2))


if __name__ == "__main__":
    main()
