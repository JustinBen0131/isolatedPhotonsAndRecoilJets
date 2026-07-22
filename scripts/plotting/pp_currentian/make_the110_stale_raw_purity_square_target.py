#!/usr/bin/env python3
"""Render the pre-bin36 square raw-purity overlay target for THE-110.

This deliberately reads the immutable July-17 full-stat point table rather
than resolving a ``current`` pointer.  It is a presentation/layout target only
until the bin36 pp-data merge passes every frozen campaign gate.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import ROOT


ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetEndErrorSize(3)

REPO = Path(__file__).resolve().parents[3]
SOURCE_CAMPAIGN = "the97_ppg12_final_accepted_triple_full_20260714_1550"
TARGET_CAMPAIGN = "the97_ppg12_bin36_final_triple_full_20260717_1322"
SOURCE_CSV = (
    REPO
    / "dataOutput/ppg12Parity"
    / SOURCE_CAMPAIGN
    / "final_pp_data_canonical_20260717/purity_overlay_current"
    / "the97_pp_data_purity_ppg12_vs_current_fullstat_points.csv"
)
OUT_DIR = (
    REPO
    / "dataOutput/ppg12Parity"
    / TARGET_CAMPAIGN
    / "preplacement_stale_raw_purity_square"
)
PNG = OUT_DIR / "the110_ppg12_raw_purity_square_overlay_STALE_prebin36.png"
POINTS_CSV = OUT_DIR / "the110_ppg12_raw_purity_square_overlay_STALE_prebin36_points.csv"
MANIFEST = OUT_DIR / "the110_ppg12_raw_purity_square_overlay_STALE_prebin36_manifest.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_rows() -> list[dict[str, float]]:
    required = {
        "pt_lo",
        "pt_hi",
        "center",
        "half_width",
        "ppg12_raw",
        "ppg12_raw_error_low",
        "ppg12_raw_error_high",
        "current_raw",
        "current_raw_error",
    }
    with SOURCE_CSV.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise RuntimeError(f"stale source lacks required raw-purity fields: {SOURCE_CSV}")
        rows = [{key: float(record[key]) for key in required} for record in reader]
    if len(rows) != 11:
        raise RuntimeError(f"expected 11 raw-purity points, found {len(rows)}: {SOURCE_CSV}")
    return rows


def make_ppg12_graph(rows: list[dict[str, float]]) -> ROOT.TGraphAsymmErrors:
    graph = ROOT.TGraphAsymmErrors(len(rows))
    graph.SetName("g_ppg12_sdcc_pulled_raw_purity_stale")
    graph.SetMarkerStyle(24)
    graph.SetMarkerSize(1.28)
    graph.SetMarkerColor(ROOT.kBlack)
    graph.SetLineColor(ROOT.kBlack)
    graph.SetLineWidth(2)
    for index, row in enumerate(rows):
        graph.SetPoint(index, row["center"], row["ppg12_raw"])
        graph.SetPointError(
            index,
            row["half_width"],
            row["half_width"],
            row["ppg12_raw_error_low"],
            row["ppg12_raw_error_high"],
        )
    return graph


def make_current_graph(rows: list[dict[str, float]]) -> ROOT.TGraphErrors:
    graph = ROOT.TGraphErrors(len(rows))
    graph.SetName("g_current_ppg19_raw_purity_stale")
    graph.SetMarkerStyle(20)
    graph.SetMarkerSize(1.18)
    graph.SetMarkerColor(ROOT.kBlack)
    graph.SetLineColor(ROOT.kBlack)
    graph.SetLineWidth(2)
    for index, row in enumerate(rows):
        graph.SetPoint(index, row["center"], row["current_raw"])
        graph.SetPointError(index, row["half_width"], row["current_raw_error"])
    return graph


def render(rows: list[dict[str, float]]) -> None:
    canvas = ROOT.TCanvas("the110_stale_raw_purity_square", "", 800, 800)
    canvas.SetFillColor(ROOT.kWhite)
    canvas.SetLeftMargin(0.16)
    canvas.SetRightMargin(0.045)
    canvas.SetTopMargin(0.045)
    canvas.SetBottomMargin(0.14)
    canvas.SetTicks(1, 1)

    frame = ROOT.TH1F("the110_stale_raw_purity_frame", "", 26, 10.0, 36.0)
    frame.SetDirectory(0)
    frame.SetStats(False)
    frame.GetYaxis().SetRangeUser(0.0, 1.16)
    # Retain the Figure-5 axis wording; the legend fixes the plotted series as
    # raw ABCD purity rather than the leakage-corrected series.
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

    ppg12 = make_ppg12_graph(rows)
    current = make_current_graph(rows)
    current.Draw("PZ SAME")
    # Draw the open PPG12 circles last: several pre-bin36 raw points coincide
    # numerically, and this preserves both requested marker conventions without
    # introducing an artificial x-offset.
    ppg12.Draw("PZ SAME")

    label = ROOT.TLatex()
    label.SetNDC(True)
    label.SetTextFont(42)
    label.SetTextColor(ROOT.kBlack)
    label.SetTextSize(0.049)
    label.DrawLatex(0.19, 0.915, "#it{#bf{sPHENIX}}")
    label.SetTextSize(0.039)
    label.DrawLatex(0.19, 0.862, "p+p  #sqrt{s} = 200 GeV")
    label.DrawLatex(0.19, 0.813, "|#eta^{#gamma}| < 0.7")
    label.SetTextColor(ROOT.kRed + 1)
    label.SetTextFont(62)
    label.SetTextSize(0.027)
    label.DrawLatex(0.51, 0.915, "STALE PRE-BIN36 PREVIEW")

    legend = ROOT.TLegend(0.190, 0.182, 0.748, 0.287)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.031)
    legend.SetMargin(0.22)
    legend.AddEntry(ppg12, "PPG12 SDCC pulled raw purity", "p")
    legend.AddEntry(current, "Current PPG19 analysis code", "p")
    legend.Draw()
    canvas.RedrawAxis()
    canvas.SaveAs(str(PNG))


def write_points(rows: list[dict[str, float]]) -> None:
    fields = [
        "pt_lo",
        "pt_hi",
        "center",
        "half_width",
        "ppg12_raw",
        "ppg12_raw_error_low",
        "ppg12_raw_error_high",
        "current_raw",
        "current_raw_error",
    ]
    with POINTS_CSV.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows([{field: row[field] for field in fields} for row in rows])


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = read_rows()
    write_points(rows)
    render(rows)
    manifest = {
        "schema_version": 1,
        "artifact_status": "STALE_PRE_BIN36_TARGET_ONLY",
        "target_campaign": TARGET_CAMPAIGN,
        "source_campaign": SOURCE_CAMPAIGN,
        "purpose": "Square raw-purity slide layout target; replace only after accepted bin36 pp-data merged ROOT is available.",
        "ian_mutated": False,
        "slides_mutated": False,
        "source_points_csv": str(SOURCE_CSV),
        "source_points_csv_sha256": sha256(SOURCE_CSV),
        "raw_purity_definition": "PPG12 quadratic physical branch; existing accepted pre-bin36 20,000-toy TRandom3(42) point table.",
        "marker_contract": {
            "ppg12": "open black circle with SDCC asymmetric statistical errors and bin-width horizontal errors",
            "current": "closed black circle with pre-bin36 toy statistical errors and bin-width horizontal errors",
        },
        "legend": ["PPG12 SDCC pulled raw purity", "Current PPG19 analysis code"],
        "canvas_pixels": [800, 800],
        "replacement_gate": "accepted bin36 pp-data merge, exact input-once audit, ROOT/provenance/content validation, and regenerated source CSV/manifest",
        "outputs": {
            "png": str(PNG),
            "png_sha256": sha256(PNG),
            "points_csv": str(POINTS_CSV),
            "points_csv_sha256": sha256(POINTS_CSV),
        },
        "renderer": {
            "path": str(Path(__file__).resolve()),
            "sha256": sha256(Path(__file__).resolve()),
        },
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"png": str(PNG), "points_csv": str(POINTS_CSV), "manifest": str(MANIFEST)}, indent=2))


if __name__ == "__main__":
    main()
