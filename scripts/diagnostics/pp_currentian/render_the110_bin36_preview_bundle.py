#!/usr/bin/env python3
"""Render the three clean THE-110 bin36 preview PNGs from diagnostic sums."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

import ROOT


ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetEndErrorSize(3)

REPO = Path(__file__).resolve().parents[3]
PLOT_DIR = REPO / "scripts/plotting/pp_currentian"
sys.path.insert(0, str(PLOT_DIR))

import make_ppg12_fig3_purity_sim_sdcc_vs_current_overlay as sim_plot  # noqa: E402
import make_the97_final_pp_data_abcd_three_panel as data_plot  # noqa: E402


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(path: Path, rows: list[dict[str, float]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def raw_graphs(rows: list[dict[str, float]]) -> tuple[ROOT.TGraphAsymmErrors, ROOT.TGraphErrors]:
    ppg12 = ROOT.TGraphAsymmErrors(len(rows))
    current = ROOT.TGraphErrors(len(rows))
    ppg12.SetName("g_ppg12_sdcc_pulled_raw_purity")
    current.SetName("g_current_ppg19_raw_purity")
    for graph, marker, size in ((ppg12, 24, 1.28), (current, 20, 1.18)):
        graph.SetMarkerStyle(marker)
        graph.SetMarkerSize(size)
        graph.SetMarkerColor(ROOT.kBlack)
        graph.SetLineColor(ROOT.kBlack)
        graph.SetLineWidth(2)
    for index, row in enumerate(rows):
        ppg12.SetPoint(index, row["center"], row["raw_purity_ppg12"])
        ppg12.SetPointError(
            index,
            row["half_width"], row["half_width"],
            row["raw_purity_ppg12_error"], row["raw_purity_ppg12_error"],
        )
        current.SetPoint(index, row["center"], row["raw_purity_current"])
        current.SetPointError(index, row["half_width"], row["raw_purity_current_error"])
    return ppg12, current


def render_square(rows: list[dict[str, float]], output: Path) -> None:
    canvas = ROOT.TCanvas("the110_bin36_raw_purity_square", "", 800, 800)
    canvas.SetCanvasSize(800, 800)
    canvas.SetFillColor(ROOT.kWhite)
    canvas.SetLeftMargin(0.16)
    canvas.SetRightMargin(0.045)
    canvas.SetTopMargin(0.045)
    canvas.SetBottomMargin(0.14)
    canvas.SetTicks(1, 1)
    frame = ROOT.TH1F("the110_bin36_raw_purity_frame", "", 26, 10.0, 36.0)
    frame.SetDirectory(0)
    frame.SetStats(False)
    frame.GetYaxis().SetRangeUser(0.0, 1.16)
    frame.GetYaxis().SetTitle("Purity")
    frame.GetXaxis().SetTitle("E_{T}^{#gamma,rec} [GeV]")
    frame.GetXaxis().SetTitleSize(0.052)
    frame.GetYaxis().SetTitleSize(0.052)
    frame.GetXaxis().SetTitleOffset(1.12)
    frame.GetYaxis().SetTitleOffset(1.22)
    frame.GetXaxis().SetLabelSize(0.043)
    frame.GetYaxis().SetLabelSize(0.043)
    frame.GetYaxis().SetNdivisions(506)
    frame.Draw("AXIS")
    ppg12, current = raw_graphs(rows)
    current.Draw("PZ SAME")
    ppg12.Draw("PZ SAME")
    label = ROOT.TLatex()
    label.SetNDC(True)
    label.SetTextFont(42)
    # Header block dropped away from the top frame edge; the purity points sit
    # near y~0.4-0.75 of a 0-1.16 axis, so there is room below the frame top.
    label.SetTextSize(0.049)
    label.DrawLatex(0.19, 0.878, "#it{#bf{sPHENIX}} Internal")
    label.SetTextSize(0.039)
    label.DrawLatex(0.19, 0.822, "p+p  #sqrt{s} = 200 GeV")
    label.DrawLatex(0.19, 0.770, "|#eta^{#gamma}| < 0.7")
    # The lower third of the frame is empty, so the legend can be much larger
    # without touching any point. Anchored where it already sat.
    legend = ROOT.TLegend(0.185, 0.150, 0.930, 0.330)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.046)
    legend.SetMargin(0.16)
    legend.AddEntry(ppg12, "PPG12 SDCC pulled raw purity", "p")
    legend.AddEntry(current, "Current PPG19 analysis code", "p")
    legend.Draw()
    canvas.RedrawAxis()
    canvas.SaveAs(str(output))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--inclusive-root", type=Path, required=True)
    parser.add_argument("--photon-root", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument(
        "--artifact-status",
        default="DIAGNOSTIC_PARTIAL_PP_NO_CANVAS_CAVEAT",
        help="Truthful status recorded in the output manifest.",
    )
    parser.add_argument(
        "--pp-input-note",
        default=(
            "Excludes only the two original pathological pp parent groups and all eight "
            "recovery children; rerender after accepted recovery substitution."
        ),
    )
    parser.add_argument(
        "--sim-input-note",
        default="Complete 6000 photonjet plus 40000 inclusivejet raw populations; independent of pp recovery.",
    )
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    current = data_plot.load_current(args.data_root)
    reference = data_plot.load_reference(data_plot.REFERENCE)
    data_rows = data_plot.combine(current, reference)
    data_png = args.outdir / "the110_bin36_data_raw_abcd_yield_and_purity.png"
    data_csv = args.outdir / "the110_bin36_data_raw_abcd_yield_and_purity_points.csv"
    data_plot.render(data_rows, data_png)
    write_csv(data_csv, data_rows)

    square_png = args.outdir / "the110_bin36_raw_purity_square_overlay.png"
    square_csv = args.outdir / "the110_bin36_raw_purity_square_overlay_points.csv"
    render_square(data_rows, square_png)
    write_csv(square_csv, data_rows)

    reference_sim = sim_plot.reference_points(sim_plot.DEFAULT_REFERENCE_ROOT)
    current_sim, counts = sim_plot.current_points(
        args.inclusive_root, args.photon_root, "unsuffixed"
    )
    sim_png = args.outdir / "the110_bin36_sim_purity_ppg12_vs_current.png"
    sim_csv = args.outdir / "the110_bin36_sim_purity_ppg12_vs_current_points.csv"
    ratios = sim_plot.render(
        reference_sim,
        current_sim,
        sim_png,
        "unsuffixed",
        "Current PPG19 analysis code",
    )
    sim_plot.write_table(sim_csv, reference_sim, current_sim, ratios)
    counts_path = args.outdir / "the110_bin36_sim_purity_current_counts.json"
    counts_path.write_text(json.dumps(counts, indent=2) + "\n")

    manifest = {
        "schema": "THE110_BIN36_THREE_PNG_PREVIEW_BUNDLE_V1",
        "artifact_status": args.artifact_status,
        "canonical": False,
        "slides_mutated": False,
        "pp_input_note": args.pp_input_note,
        "sim_input_note": args.sim_input_note,
        "data_root": str(args.data_root),
        "data_root_sha256": sha256(args.data_root),
        "inclusive_root": str(args.inclusive_root),
        "inclusive_root_sha256": sha256(args.inclusive_root),
        "photon_root": str(args.photon_root),
        "photon_root_sha256": sha256(args.photon_root),
        "outputs": {
            "raw_purity_square": {"path": str(square_png), "sha256": sha256(square_png)},
            "data_raw_abcd_three_panel": {"path": str(data_png), "sha256": sha256(data_png)},
            "sim_purity_overlay": {"path": str(sim_png), "sha256": sha256(sim_png)},
            "data_points": {"path": str(data_csv), "sha256": sha256(data_csv)},
            "square_points": {"path": str(square_csv), "sha256": sha256(square_csv)},
            "sim_points": {"path": str(sim_csv), "sha256": sha256(sim_csv)},
            "sim_counts": {"path": str(counts_path), "sha256": sha256(counts_path)},
        },
        "renderer": str(Path(__file__).resolve()),
        "renderer_sha256": sha256(Path(__file__).resolve()),
    }
    manifest_path = args.outdir / "the110_bin36_three_png_preview_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"pngs": [str(square_png), str(data_png), str(sim_png)], "manifest": str(manifest_path)}, indent=2))


if __name__ == "__main__":
    main()
