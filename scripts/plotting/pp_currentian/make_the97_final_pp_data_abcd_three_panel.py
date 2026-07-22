#!/usr/bin/env python3
"""Render the final full-stat THE97 pp-data raw-ABCD parity diagnostic."""

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
    ppg12_toy_estimate,
)
from make_the97_pp_data_partial_purity_contract_v1 import (
    DATA_DIR,
    DATA_NAMES,
    EXPECTED_EDGES,
    assert_binning,
    open_root,
)


ROOT.gROOT.SetBatch(True)

REPO = Path(__file__).resolve().parents[3]
CAMPAIGN = "the97_ppg12_final_accepted_triple_full_20260714_1550"
SCHEMA = "THE97_FINAL_PP_DATA_RAW_ABCD_THREE_PANEL_V1"
POINTER = REPO / "dataOutput/current_recoiljets_artifacts/current/pp_data_merged/current.json"
REFERENCE = REPO / "dataOutput/ppg12Parity/the93_ppg12_canonical_full_20260706_2145/data_abcd_fig27_ppg12_scaledtrigger30/fig27_abcd_yield_current_vs_ppg12_overlay_ratio.csv"
OUTDIR = REPO / f"dataOutput/ppg12Parity/{CAMPAIGN}/final_pp_data_canonical_20260717/raw_abcd_three_panel"

REGIONS = (
    # Colour encodes the region; marker fill encodes the source (filled circle
    # = this analysis, open circle = PPG12 SDCC).  All four use circles so the
    # two encodings stay independent.
    ("A", "A: tight iso", ROOT.kGreen + 3, 20, 24),
    ("B", "B: tight noniso", ROOT.kRed + 1, 20, 24),
    ("C", "C: nontight iso", ROOT.kBlue + 1, 20, 24),
    ("D", "D: nontight noniso", ROOT.kMagenta + 1, 20, 24),
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def resolve_current(pointer: Path) -> tuple[Path, dict[str, object]]:
    payload = json.loads(pointer.read_text())
    if payload.get("sample_key") != "pp_data_merged":
        raise RuntimeError(f"wrong current sample key in {pointer}")
    if not payload.get("current_entry_id"):
        raise RuntimeError(f"pp_data_merged pointer has no current_entry_id: {pointer}")
    if payload.get("campaign_tag") != CAMPAIGN:
        raise RuntimeError(
            f"pp_data_merged points to {payload.get('campaign_tag')}, not {CAMPAIGN}"
        )
    roots = payload.get("root_paths") or []
    if len(roots) != 1:
        raise RuntimeError(f"expected one current pp-data ROOT, found {roots}")
    root = Path(str(roots[0]))
    if not root.is_file():
        raise RuntimeError(f"missing current pp-data ROOT: {root}")
    return root, payload


def load_current(root_path: Path) -> dict[str, list[dict[str, float]]]:
    root_file = open_root(root_path)
    histograms: dict[str, ROOT.TH1] = {}
    for region, name in DATA_NAMES.items():
        hist = root_file.Get(f"{DATA_DIR}/{name}")
        if not hist or not hist.InheritsFrom("TH1"):
            raise RuntimeError(f"missing TH1 {DATA_DIR}/{name} in {root_path}")
        clone = hist.Clone(f"final_current_{name}")
        clone.SetDirectory(0)
        histograms[region] = clone
    root_file.Close()
    assert_binning(list(histograms.values()))
    rows: dict[str, list[dict[str, float]]] = {key: [] for key, *_ in REGIONS}
    for region, *_ in REGIONS:
        hist = histograms[region]
        for index in range(1, 12):
            lo, hi = EXPECTED_EDGES[index - 1], EXPECTED_EDGES[index]
            rows[region].append(
                {
                    "pt_lo": float(lo),
                    "pt_hi": float(hi),
                    "center": 0.5 * (lo + hi),
                    "half_width": 0.5 * (hi - lo),
                    "yield": float(hist.GetBinContent(index)),
                    "error": float(hist.GetBinError(index)),
                    "effective_count": effective_count(hist, index),
                }
            )
    return rows


def load_reference(path: Path) -> dict[str, list[dict[str, float]]]:
    rows: dict[str, list[dict[str, float]]] = {key: [] for key, *_ in REGIONS}
    with path.open(newline="") as handle:
        for raw in csv.DictReader(handle):
            region = raw["region"]
            if region not in rows:
                continue
            value = float(raw["ppg12_yield"])
            error = float(raw["ppg12_error"])
            if value <= 0.0 or error <= 0.0:
                raise RuntimeError(f"nonpositive PPG12 yield/error in {region}: {raw}")
            lo, hi = float(raw["pt_lo"]), float(raw["pt_hi"])
            rows[region].append(
                {
                    "pt_lo": lo,
                    "pt_hi": hi,
                    "center": 0.5 * (lo + hi),
                    "half_width": 0.5 * (hi - lo),
                    "yield": value,
                    "error": error,
                    "effective_count": value * value / (error * error),
                }
            )
    for region, values in rows.items():
        if len(values) != 11:
            raise RuntimeError(f"reference contains {len(values)} rather than 11 {region} bins")
        edges = tuple(row["pt_lo"] for row in values) + (values[-1]["pt_hi"],)
        if edges != EXPECTED_EDGES:
            raise RuntimeError(f"reference binning mismatch for {region}: {edges}")
    return rows


def raw_purity(
    samples: dict[str, list[dict[str, float]]], label: str
) -> list[dict[str, float]]:
    # Resetting both samples to the identical seed is intentional: each side
    # receives the same PPG12 20k-toy contract rather than a seed-order artifact.
    rng = ROOT.TRandom3(42)
    output: list[dict[str, float]] = []
    for index in range(11):
        values = tuple(samples[key][index]["yield"] for key in "ABCD")
        counts = tuple(samples[key][index]["effective_count"] for key in "ABCD")
        raw, error, _, _, diagnostics = ppg12_toy_estimate(
            rng,
            values,
            counts,
            (0.0, 0.0, 0.0),
            (0.0, 0.0, 0.0),
            f"{label}_{index}",
        )
        output.append(
            {
                "purity": raw,
                "error": error,
                "toy_entries": diagnostics["raw_toy_entries"],
                "toy_underflow": diagnostics["raw_toy_underflow"],
                "toy_overflow": diagnostics["raw_toy_overflow"],
            }
        )
    return output


def ratio(value: float, error: float, reference: float, reference_error: float) -> tuple[float, float]:
    if value <= 0.0 or reference <= 0.0:
        return math.nan, math.nan
    result = value / reference
    uncertainty = abs(result) * math.hypot(error / value, reference_error / reference)
    return result, uncertainty


def combine(
    current: dict[str, list[dict[str, float]]],
    reference: dict[str, list[dict[str, float]]],
) -> list[dict[str, float]]:
    current_purity = raw_purity(current, "current_raw")
    reference_purity = raw_purity(reference, "ppg12_raw")
    rows: list[dict[str, float]] = []
    for index in range(11):
        lo, hi = EXPECTED_EDGES[index], EXPECTED_EDGES[index + 1]
        row: dict[str, float] = {
            "pt_lo": float(lo),
            "pt_hi": float(hi),
            "center": 0.5 * (lo + hi),
            "half_width": 0.5 * (hi - lo),
        }
        for region in "ABCD":
            cur, ref = current[region][index], reference[region][index]
            val, err = ratio(cur["yield"], cur["error"], ref["yield"], ref["error"])
            row.update(
                {
                    f"{region}_current": cur["yield"],
                    f"{region}_current_error": cur["error"],
                    f"{region}_current_effective_count": cur["effective_count"],
                    f"{region}_ppg12": ref["yield"],
                    f"{region}_ppg12_error": ref["error"],
                    f"{region}_ppg12_effective_count": ref["effective_count"],
                    f"{region}_current_over_ppg12": val,
                    f"{region}_current_over_ppg12_error": err,
                }
            )
        pval, perr = ratio(
            current_purity[index]["purity"],
            current_purity[index]["error"],
            reference_purity[index]["purity"],
            reference_purity[index]["error"],
        )
        row.update(
            {
                "raw_purity_current": current_purity[index]["purity"],
                "raw_purity_current_error": current_purity[index]["error"],
                "raw_purity_current_toy_entries": current_purity[index]["toy_entries"],
                "raw_purity_ppg12": reference_purity[index]["purity"],
                "raw_purity_ppg12_error": reference_purity[index]["error"],
                "raw_purity_ppg12_toy_entries": reference_purity[index]["toy_entries"],
                "raw_purity_current_over_ppg12": pval,
                "raw_purity_current_over_ppg12_error": perr,
            }
        )
        rows.append(row)
    return rows


def graph(
    name: str,
    rows: list[dict[str, float]],
    value_key: str,
    error_key: str,
    color: int,
    marker: int,
) -> ROOT.TGraphErrors:
    result = ROOT.TGraphErrors(len(rows))
    result.SetName(name)
    result.SetMarkerColor(color)
    result.SetLineColor(color)
    result.SetMarkerStyle(marker)
    result.SetMarkerSize(1.05)
    result.SetLineWidth(2)
    for index, row in enumerate(rows):
        result.SetPoint(index, row["center"], row[value_key])
        result.SetPointError(index, row["half_width"], row[error_key])
    return result


def render(rows: list[dict[str, float]], output: Path) -> None:
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetTextFont(42)
    ROOT.gStyle.SetLabelFont(42, "XYZ")
    ROOT.gStyle.SetTitleFont(42, "XYZ")

    canvas = ROOT.TCanvas("c_final_raw_abcd_three_panel", "", 940, 1120)
    canvas.SetFillColor(ROOT.kWhite)
    top = ROOT.TPad("top", "", 0.0, 0.54, 1.0, 1.0)
    middle = ROOT.TPad("middle", "", 0.0, 0.27, 1.0, 0.54)
    bottom = ROOT.TPad("bottom", "", 0.0, 0.0, 1.0, 0.27)
    for pad in (top, middle, bottom):
        pad.SetFillColor(ROOT.kWhite)
        pad.SetLeftMargin(0.135)
        pad.SetRightMargin(0.035)
    top.SetTopMargin(0.045); top.SetBottomMargin(0.015); top.SetLogy(True)
    middle.SetTopMargin(0.018); middle.SetBottomMargin(0.015)
    bottom.SetTopMargin(0.018); bottom.SetBottomMargin(0.29)
    top.Draw(); middle.Draw(); bottom.Draw()

    top.cd()
    ymax = max(row[f"{region}_current"] + row[f"{region}_current_error"] for row in rows for region in "ABCD")
    frame = ROOT.TH1F("frame_abcd_yield", "", 26, 10.0, 36.0)
    frame.SetStats(False)
    frame.GetYaxis().SetRangeUser(1.0, max(2.5e5, 3.5 * ymax))
    frame.GetYaxis().SetTitle("Raw ABCD yield")
    frame.GetYaxis().SetTitleSize(0.062); frame.GetYaxis().SetLabelSize(0.048)
    frame.GetYaxis().SetTitleOffset(0.95); frame.GetXaxis().SetLabelSize(0.0)
    frame.Draw("AXIS")
    reference_graphs = {}
    current_graphs = {}
    for region, _, color, filled, opened in REGIONS:
        reference_graphs[region] = graph(
            f"g_ppg12_{region}", rows, f"{region}_ppg12", f"{region}_ppg12_error", color, opened
        )
        current_graphs[region] = graph(
            f"g_current_{region}", rows, f"{region}_current", f"{region}_current_error", color, filled
        )
        reference_graphs[region].Draw("P SAME")
        current_graphs[region].Draw("P SAME")
    # The yield spectrum falls left-to-right, so the top-right corner is the
    # empty one for the header block and the bottom-left is free for a large
    # legend.  Both were previously the other way round.
    label = ROOT.TLatex(); label.SetNDC(True); label.SetTextFont(42)
    # Right-aligned but pulled in from the frame edge and down from the top, so
    # the block sits inside the empty upper-right corner rather than against it.
    label.SetTextAlign(31)
    label.SetTextSize(0.055); label.DrawLatex(0.900, 0.870, "#it{#bf{sPHENIX}} Internal")
    label.SetTextSize(0.038); label.DrawLatex(0.900, 0.800, "p+p  #sqrt{s} = 200 GeV, |#eta^{#gamma}| < 0.7")
    # Two columns, filled row-wise: source pair heads the columns, then the four
    # ABCD regions pair beneath.  Bounded to x < 0.575 and y < 0.34, which is the
    # region the falling spectrum leaves empty (points at E_T >= 24 GeV sit at
    # x >= 0.582, y = 0.14-0.33 in pad NDC).
    legend = ROOT.TLegend(0.145, 0.055, 0.575, 0.335)
    legend.SetBorderSize(0); legend.SetFillStyle(0); legend.SetTextFont(42); legend.SetTextSize(0.036)
    legend.SetNColumns(2); legend.SetMargin(0.20); legend.SetColumnSeparation(0.01)
    # Neutral black keys: these two entries describe the marker-fill convention,
    # not a region, so they must not borrow region A's colour.
    source_keys = []
    for style in (24, 20):
        key = ROOT.TGraphErrors(1)
        key.SetPoint(0, -1.0e6, -1.0e6)  # off-frame; legend key only
        key.SetMarkerColor(ROOT.kBlack); key.SetLineColor(ROOT.kBlack)
        key.SetMarkerStyle(style); key.SetMarkerSize(1.05); key.SetLineWidth(2)
        source_keys.append(key)
    legend.AddEntry(source_keys[0], "PPG12 SDCC", "p")
    legend.AddEntry(source_keys[1], "This analysis", "p")
    for region, region_label, *_ in REGIONS:
        legend.AddEntry(current_graphs[region], region_label, "p")
    legend.Draw(); top.RedrawAxis()

    middle.cd()
    mid = ROOT.TH1F("frame_abcd_ratio", "", 26, 10.0, 36.0)
    mid.SetStats(False); mid.GetYaxis().SetRangeUser(0.45, 2.15)
    # Per-panel y titles are replaced by one shared label spanning both ratio
    # pads (drawn on the canvas after the pads); this also stops the long
    # per-pad titles from being clipped at the left canvas edge.
    mid.GetYaxis().SetTitle("")
    mid.GetYaxis().SetTitleSize(0.090); mid.GetYaxis().SetLabelSize(0.068)
    mid.GetYaxis().SetTitleOffset(0.64); mid.GetXaxis().SetLabelSize(0.0)
    mid.GetYaxis().SetNdivisions(506); mid.Draw("AXIS")
    unity_mid = ROOT.TLine(10.0, 1.0, 36.0, 1.0)
    unity_mid.SetLineColor(ROOT.kGray + 2); unity_mid.SetLineStyle(7); unity_mid.Draw()
    ratio_graphs = []
    for region, _, color, filled, _ in REGIONS:
        ratio_graphs.append(graph(
            f"g_yield_ratio_{region}", rows,
            f"{region}_current_over_ppg12", f"{region}_current_over_ppg12_error",
            color, filled,
        ))
        ratio_graphs[-1].Draw("P SAME")
    mid_caption = ROOT.TLatex(); mid_caption.SetNDC(True); mid_caption.SetTextFont(42)
    mid_caption.SetTextSize(0.078)
    mid_caption.DrawLatex(0.155, 0.075, "Ratio of Raw Counts in Regions ABCD")
    middle.RedrawAxis()

    bottom.cd()
    finite = [
        (row["raw_purity_current_over_ppg12"], row["raw_purity_current_over_ppg12_error"])
        for row in rows
        if math.isfinite(row["raw_purity_current_over_ppg12"])
        and math.isfinite(row["raw_purity_current_over_ppg12_error"])
    ]
    low = min(0.72, min(value - error for value, error in finite) - 0.06) if finite else 0.65
    high = max(1.18, max(value + error for value, error in finite) + 0.06) if finite else 1.25
    bot = ROOT.TH1F("frame_raw_purity_ratio", "", 26, 10.0, 36.0)
    bot.SetStats(False); bot.GetYaxis().SetRangeUser(max(0.0, low), high)
    bot.GetXaxis().SetTitle("E_{T}^{#gamma,rec} [GeV]")
    bot.GetYaxis().SetTitle("")  # shared label drawn across both ratio pads
    bot.GetXaxis().SetTitleSize(0.102); bot.GetXaxis().SetLabelSize(0.077)
    bot.GetYaxis().SetTitleSize(0.082); bot.GetYaxis().SetLabelSize(0.067)
    bot.GetYaxis().SetTitleOffset(0.72); bot.GetXaxis().SetTitleOffset(1.08)
    bot.GetYaxis().SetNdivisions(505); bot.Draw("AXIS")
    unity_bot = ROOT.TLine(10.0, 1.0, 36.0, 1.0)
    unity_bot.SetLineColor(ROOT.kGray + 2); unity_bot.SetLineStyle(7); unity_bot.Draw()
    purity_ratio_graph = graph(
        "g_raw_purity_ratio", rows,
        "raw_purity_current_over_ppg12", "raw_purity_current_over_ppg12_error",
        ROOT.kBlack, 20,
    )
    purity_ratio_graph.Draw("P SAME")
    bot_caption = ROOT.TLatex(); bot_caption.SetNDC(True); bot_caption.SetTextFont(42)
    bot_caption.SetTextSize(0.072)
    bot_caption.DrawLatex(0.155, 0.345, "Raw ABCD purity ratio")
    bottom.RedrawAxis()

    # One y-axis label for both ratio pads, centred on their combined span
    # (bottom pad 0.00-0.27, middle pad 0.27-0.54 in canvas NDC).
    canvas.cd()
    shared_y = ROOT.TLatex(); shared_y.SetNDC(True); shared_y.SetTextFont(42)
    shared_y.SetTextAngle(90); shared_y.SetTextAlign(22); shared_y.SetTextSize(0.026)
    shared_y.DrawLatex(0.038, 0.27, "This analysis / PPG12 SDCC")
    canvas.SaveAs(str(output))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pointer", type=Path, default=POINTER)
    parser.add_argument("--reference-csv", type=Path, default=REFERENCE)
    parser.add_argument("--outdir", type=Path, default=OUTDIR)
    args = parser.parse_args()

    root_path, pointer = resolve_current(args.pointer)
    current = load_current(root_path)
    reference = load_reference(args.reference_csv)
    rows = combine(current, reference)
    args.outdir.mkdir(parents=True, exist_ok=True)
    png = args.outdir / "the97_final_pp_data_raw_abcd_yield_and_purity_current_over_ppg12.png"
    csv_path = args.outdir / "the97_final_pp_data_raw_abcd_yield_and_purity_current_over_ppg12.csv"
    manifest = args.outdir / "the97_final_pp_data_raw_abcd_yield_and_purity_current_over_ppg12_manifest.json"
    render(rows, png)
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    manifest.write_text(
        json.dumps(
            {
                "schema": SCHEMA,
                "status": "FULL_STAT_CURRENT",
                "campaign_tag": CAMPAIGN,
                "current_artifact_id": pointer["current_entry_id"],
                "current_pointer": str(args.pointer),
                "current_pointer_sha256": sha256(args.pointer),
                "current_root": str(root_path),
                "current_root_sha256": sha256(root_path),
                "current_root_bytes": root_path.stat().st_size,
                "data_namespace": DATA_DIR,
                "data_histograms": DATA_NAMES,
                "ppg12_reference_csv": str(args.reference_csv),
                "ppg12_reference_csv_sha256": sha256(args.reference_csv),
                "ppg12_reference_root": "/sphenix/user/shuhangli/ppg12/efficiencytool/results/data_histo_bdt_nom.root",
                "bin_edges_gev": EXPECTED_EDGES,
                "yield_comparison": "raw unsuffixed A/B/C/D bin yields; PPG12 open and current full-stat filled; no normalization or fitted scale",
                "raw_purity_method": "each sample's own raw A/B/C/D; PPG12 quadratic physical branch; TRandom3(42) reset per sample; 20000 effective-Poisson toys per bin; PPG12 Gaussian-fit estimator",
                "sim_leakage_used": False,
                "partial_coverage_label_used": False,
                "png": str(png),
                "png_sha256": sha256(png),
                "source_csv": str(csv_path),
                "source_csv_sha256": sha256(csv_path),
                "plot_script": str(Path(__file__).resolve()),
                "plot_script_sha256": sha256(Path(__file__).resolve()),
            },
            indent=2,
            sort_keys=True,
        ) + "\n"
    )
    print(png)
    print(csv_path)
    print(manifest)


if __name__ == "__main__":
    main()
