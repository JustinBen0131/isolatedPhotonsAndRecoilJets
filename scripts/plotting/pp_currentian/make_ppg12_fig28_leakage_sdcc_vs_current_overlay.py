#!/usr/bin/env python3
"""Overlay PPG12 Fig. 28 leakage curves with current RecoilJets output.

Top panel:
  PPG12 SDCC source as step lines.
  Current RecoilJets photon+jet output as filled markers.

Bottom panel:
  PPG12 SDCC / Current output.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from array import array
from pathlib import Path

import ROOT


REPO = Path(__file__).resolve().parents[3]
DEFAULT_PPG12_CSV = (
    REPO
    / "dataOutput/ppg12Parity/the76_ppg12_fig28_leakage_reference"
    / "ppg12_fig28_leakage_sdcc_points.csv"
)
DEFAULT_CURRENT_JSON = (
    REPO
    / "dataOutput/current_recoiljets_artifacts/current/pp_sim_photonjet_merged/current.json"
)
DEFAULT_OUTDIR = REPO / "dataOutput/ppg12Parity/the76_ppg12_fig28_leakage_reference"

BASE_SERIES = [
    {
        "region": "B",
        "color": ROOT.kBlack,
        "marker": 20,
        "label": "#it{N}^{sig}_{B}/#it{N}^{sig}_{A} tight noniso",
        "base_num": "h_tight_noniso_cluster",
    },
    {
        "region": "C",
        "color": ROOT.kRed,
        "marker": 20,
        "label": "#it{N}^{sig}_{C}/#it{N}^{sig}_{A} nontight iso",
        "base_num": "h_nontight_iso_cluster",
    },
    {
        "region": "D",
        "color": ROOT.kBlue,
        "marker": 20,
        "label": "#it{N}^{sig}_{D}/#it{N}^{sig}_{A} nontight noniso",
        "base_num": "h_nontight_noniso_cluster",
    },
]

CURRENT_OBJECT_FAMILIES = {
    "signal": {
        "suffix": "_signal_0",
        "description": "truth-matched signal histograms used by the final photon-yield leakage correction",
    },
    "macro_raw": {
        "suffix": "_0",
        "description": (
            "non-suffixed all-reconstructed-candidate signal-MC ABCD histograms "
            "matching the PPG12 Fig.28 display-macro object names; these are not "
            "the truth-matched leakage factors used by the purity correction"
        ),
    },
}


PPG12_SOURCE_ROOTS = {
    "analysis_note_fig28": "preserved PPG12 IAN Fig.28 analysis-note backup values",
    "macro_recomputed": "/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_bdt_nom.root",
    "persisted_final": "/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_bdt_nom_mc.root",
    "persisted_nom_combined": "/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_bdt_nom_mc.root",
    "persisted_all": "/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_bdt_all_mc.root",
}


def build_series(current_object_family: str) -> tuple[str, list[dict[str, object]]]:
    if current_object_family not in CURRENT_OBJECT_FAMILIES:
        raise RuntimeError(f"unknown current object family {current_object_family}")
    suffix = CURRENT_OBJECT_FAMILIES[current_object_family]["suffix"]
    current_den = f"SIM/h_tight_iso_cluster{suffix}"
    series = []
    for spec in BASE_SERIES:
        item = dict(spec)
        item["current_num"] = f"SIM/{spec['base_num']}{suffix}"
        series.append(item)
    return current_den, series


def resolve_current_root(path: Path) -> Path:
    payload = json.loads(path.read_text())
    roots = payload.get("root_paths") or []
    if not roots:
        raise RuntimeError(f"no root_paths in {path}")
    root = Path(roots[0])
    if not root.is_absolute():
        root = REPO / root
    if not root.exists():
        raise FileNotFoundError(root)
    return root


def read_ppg12_rows(
    path: Path, source_kind: str, series: list[dict[str, object]]
) -> dict[str, list[dict[str, float]]]:
    out = {spec["region"]: [] for spec in series}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            if row["source_kind"] != source_kind:
                continue
            region = row["region"]
            if region not in out:
                continue
            out[region].append(
                {
                    "bin_index": int(row["bin_index"]),
                    "x_low": float(row["x_low_gev"]),
                    "x_high": float(row["x_high_gev"]),
                    "x_center": float(row["x_center_gev"]),
                    "y": float(row["leakage"]),
                    "err": float(row["error"]),
                }
            )
    for rows in out.values():
        rows.sort(key=lambda r: r["bin_index"])
    return out


def read_ppg12_root(
    path: Path, series: list[dict[str, object]]
) -> dict[str, list[dict[str, float]]]:
    f = ROOT.TFile.Open(str(path))
    if not f or f.IsZombie():
        raise RuntimeError(f"failed to open PPG12 ROOT: {path}")
    out: dict[str, list[dict[str, float]]] = {}
    for spec in series:
        region = spec["region"]
        hist_name = f"h_leak_{region}"
        hist = f.Get(hist_name)
        if not hist:
            raise RuntimeError(f"missing {hist_name} in {path}")
        rows: list[dict[str, float]] = []
        for i in range(1, hist.GetNbinsX() + 1):
            rows.append(
                {
                    "bin_index": i,
                    "x_low": hist.GetXaxis().GetBinLowEdge(i),
                    "x_high": hist.GetXaxis().GetBinUpEdge(i),
                    "x_center": hist.GetXaxis().GetBinCenter(i),
                    "y": hist.GetBinContent(i),
                    "err": hist.GetBinError(i),
                }
            )
        out[region] = rows
    f.Close()
    return out


def read_analysis_note_rows(
    path: Path, series: list[dict[str, object]]
) -> dict[str, list[dict[str, float]]]:
    fields = {
        "B": "ppg12_fig28_f_b",
        "C": "ppg12_fig28_f_c",
        "D": "ppg12_fig28_f_d",
    }
    source_rows = list(csv.DictReader(path.open(newline="")))
    out: dict[str, list[dict[str, float]]] = {}
    for spec in series:
        region = spec["region"]
        out[region] = [
            {
                "bin_index": i,
                "x_low": float(row["pt_lo"]),
                "x_high": float(row["pt_hi"]),
                "x_center": float(row["pt_center"]),
                "y": float(row[fields[region]]),
                "err": 0.0,
            }
            for i, row in enumerate(source_rows, start=1)
        ]
    return out


def make_ppg12_hist(rows: list[dict[str, float]], name: str, color: int) -> ROOT.TH1D:
    edges = [rows[0]["x_low"]]
    edges.extend(row["x_high"] for row in rows)
    h = ROOT.TH1D(name, "", len(edges) - 1, array("d", edges))
    h.SetDirectory(0)
    h.SetLineColor(color)
    h.SetLineWidth(2)
    h.SetMarkerSize(0)
    for i, row in enumerate(rows, start=1):
        h.SetBinContent(i, row["y"])
        h.SetBinError(i, row["err"])
    return h


def current_points(
    root_path: Path, current_den: str, series: list[dict[str, object]]
) -> dict[str, list[dict[str, float]]]:
    f = ROOT.TFile.Open(str(root_path))
    if not f or f.IsZombie():
        raise RuntimeError(f"failed to open current ROOT: {root_path}")
    den = f.Get(current_den)
    if not den:
        raise RuntimeError(f"missing denominator {current_den} in {root_path}")

    out: dict[str, list[dict[str, float]]] = {}
    for spec in series:
        num = f.Get(spec["current_num"])
        if not num:
            raise RuntimeError(f"missing numerator {spec['current_num']} in {root_path}")
        rows: list[dict[str, float]] = []
        for i in range(1, num.GetNbinsX() + 1):
            n = num.GetBinContent(i)
            d = den.GetBinContent(i)
            ne = num.GetBinError(i)
            de = den.GetBinError(i)
            y = n / d if d else 0.0
            rel2 = 0.0
            if n > 0:
                rel2 += (ne / n) ** 2
            if d > 0:
                rel2 += (de / d) ** 2
            err = abs(y) * math.sqrt(rel2) if y else 0.0
            rows.append(
                {
                    "bin_index": i,
                    "x_low": num.GetXaxis().GetBinLowEdge(i),
                    "x_high": num.GetXaxis().GetBinUpEdge(i),
                    "x_center": num.GetXaxis().GetBinCenter(i),
                    "x_err": 0.0,
                    "y": y,
                    "err": err,
                    "num": n,
                    "den": d,
                    "num_err": ne,
                    "den_err": de,
                }
            )
        out[spec["region"]] = rows
    f.Close()
    return out


def make_graph(rows: list[dict[str, float]], name: str, color: int, marker: int) -> ROOT.TGraphErrors:
    g = ROOT.TGraphErrors(len(rows))
    g.SetName(name)
    g.SetMarkerStyle(marker)
    g.SetMarkerSize(0.9)
    g.SetMarkerColor(color)
    g.SetLineColor(color)
    for i, row in enumerate(rows):
        g.SetPoint(i, row["x_center"], row["y"])
        g.SetPointError(i, row.get("x_err", 0.0), row.get("err", 0.0))
    return g


def ratio_rows(
    ppg12: list[dict[str, float]], current: list[dict[str, float]]
) -> list[dict[str, float]]:
    by_edge = {(row["x_low"], row["x_high"]): row for row in current}
    rows: list[dict[str, float]] = []
    for row in ppg12:
        cur = by_edge.get((row["x_low"], row["x_high"]))
        if not cur:
            continue
        y = row["y"] / cur["y"] if cur["y"] else 0.0
        rel2 = 0.0
        if row["y"] > 0:
            rel2 += (row["err"] / row["y"]) ** 2
        if cur["y"] > 0:
            rel2 += (cur["err"] / cur["y"]) ** 2
        rows.append(
            {
                "bin_index": row["bin_index"],
                "x_low": row["x_low"],
                "x_high": row["x_high"],
                "x_center": row["x_center"],
                "x_err": 0.0,
                "y": y,
                "err": abs(y) * math.sqrt(rel2) if y else 0.0,
            }
        )
    return rows


def draw_line_legend(text_x: float, y: float, color: int, text: str, size: float) -> None:
    line = ROOT.TLine()
    line.SetLineColor(color)
    line.SetLineWidth(2)
    line.DrawLineNDC(text_x - 0.95 * size, y, text_x - 0.15 * size, y)
    latex = ROOT.TLatex()
    latex.SetNDC(True)
    latex.SetTextFont(42)
    latex.SetTextSize(size)
    latex.SetTextAlign(12)
    latex.DrawLatex(text_x, y, text)


def ratio_axis_bounds(ratio_rows_by_region: dict[str, list[dict[str, float]]]) -> tuple[float, float]:
    lower_values: list[float] = []
    upper_values: list[float] = []
    for rows in ratio_rows_by_region.values():
        for row in rows:
            y = row.get("y")
            err = row.get("err", 0.0)
            if y is None or err is None:
                continue
            if math.isfinite(y) and math.isfinite(err) and y > 0.0 and err >= 0.0:
                lower_values.append(y - err)
                upper_values.append(y + err)
    if not upper_values:
        return 0.5, 3.0

    # Bound the full statistical error bars, not only their central values.
    low = min(min(lower_values), 1.0)
    high = max(max(upper_values), 1.0)
    span = high - low
    pad = max(0.025, 0.07 * span)
    low = max(0.0, low - pad)
    high = high + pad

    # Round outward to simple tick-friendly boundaries while preserving a tight frame.
    step = 0.05
    low = math.floor(low / step) * step
    high = math.ceil(high / step) * step
    if high - low < 0.2:
        mid = 0.5 * (high + low)
        low = max(0.0, mid - 0.1)
        high = mid + 0.1
    return low, high


def render(
    ppg12: dict[str, list[dict[str, float]]],
    current: dict[str, list[dict[str, float]]],
    output_png: Path,
    series: list[dict[str, object]],
    source_kind: str,
    current_object_family: str,
) -> dict[str, object]:
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetTextFont(42)
    ROOT.gStyle.SetLabelFont(42, "XYZ")
    ROOT.gStyle.SetTitleFont(42, "XYZ")
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)

    canvas = ROOT.TCanvas("c_fig28_overlay", "", 640, 720)
    top = ROOT.TPad("top", "", 0.0, 0.31, 1.0, 1.0)
    bot = ROOT.TPad("bot", "", 0.0, 0.0, 1.0, 0.31)
    top.SetLeftMargin(0.15)
    top.SetRightMargin(0.03)
    top.SetTopMargin(0.07)
    top.SetBottomMargin(0.02)
    bot.SetLeftMargin(0.15)
    bot.SetRightMargin(0.03)
    bot.SetTopMargin(0.03)
    bot.SetBottomMargin(0.32)
    for pad in (top, bot):
        pad.SetTicks(1, 1)
        pad.Draw()

    ppg_hists = {}
    ppg_graphs = {}
    cur_graphs = {}
    ratio_graphs = {}
    ratio_rows_by_region: dict[str, list[dict[str, float]]] = {}
    ratio_table: list[dict[str, object]] = []

    for spec in series:
        region = spec["region"]
        ppg_hists[region] = make_ppg12_hist(
            ppg12[region], f"h_ppg12_{region}", spec["color"]
        )
        ppg_graphs[region] = make_graph(
            ppg12[region], f"g_ppg12_{region}", spec["color"], 24
        )
        cur_graphs[region] = make_graph(
            current[region], f"g_current_{region}", spec["color"], spec["marker"]
        )
        rr = ratio_rows(ppg12[region], current[region])
        ratio_rows_by_region[region] = rr
        ratio_graphs[region] = make_graph(
            rr, f"g_ratio_{region}", spec["color"], spec["marker"]
        )
        by_edge = {(row["x_low"], row["x_high"]): row for row in current[region]}
        for row in ppg12[region]:
            cur = by_edge.get((row["x_low"], row["x_high"]))
            if not cur:
                continue
            ratio = row["y"] / cur["y"] if cur["y"] else None
            ratio_table.append(
                {
                    "region": region,
                    "x_low_gev": row["x_low"],
                    "x_high_gev": row["x_high"],
                    "x_center_gev": row["x_center"],
                    "ppg12_sdcc": row["y"],
                    "ppg12_sdcc_error": row["err"],
                    "current_output": cur["y"],
                    "current_output_error": cur["err"],
                    "sdcc_over_current": ratio,
                }
            )

    top.cd()
    frame_top = ROOT.TH1F("frame_fig28_overlay_top", "", 43, 7.0, 50.0)
    frame_top.SetStats(False)
    frame_top.GetXaxis().SetRangeUser(10.0, 35.0)
    frame_top.GetYaxis().SetRangeUser(0.0, 1.3)
    frame_top.GetYaxis().SetTitle("Signal leakage")
    frame_top.GetYaxis().SetTitleSize(0.060)
    frame_top.GetYaxis().SetLabelSize(0.048)
    frame_top.GetYaxis().SetTitleOffset(0.95)
    frame_top.GetXaxis().SetLabelSize(0.0)
    frame_top.GetXaxis().SetTitleSize(0.0)
    frame_top.GetXaxis().SetNdivisions(505)
    frame_top.GetYaxis().SetNdivisions(507)
    frame_top.Draw("axis")

    ppg_as_markers = source_kind == "persisted_nom_combined"
    for spec in series:
        if ppg_as_markers:
            ppg_graphs[spec["region"]].Draw("same p")
        else:
            ppg_hists[spec["region"]].Draw("same hist")
    for spec in series:
        cur_graphs[spec["region"]].Draw("same p")

    latex = ROOT.TLatex()
    latex.SetNDC(True)
    latex.SetTextFont(42)
    latex.SetTextSize(0.037)
    latex.SetTextAlign(12)
    latex.DrawLatex(0.17, 0.89, "#bf{#it{sPHENIX}} Internal")
    latex.DrawLatex(0.17, 0.84, "#it{p}+#it{p} #kern[-0.1]{#sqrt{#it{s}} = 200 GeV}")
    latex.DrawLatex(0.17, 0.79, "|#it{#eta^{#gamma}}| < 0.7")
    latex.DrawLatex(0.17, 0.74, "PYTHIA Signal")
    for y, spec in zip([0.89, 0.81, 0.73], series):
        draw_line_legend(0.47, y, spec["color"], spec["label"], 0.038)

    src_leg = ROOT.TLegend(0.17, 0.48, 0.61, 0.61)
    src_leg.SetBorderSize(0)
    src_leg.SetFillStyle(0)
    src_leg.SetTextFont(42)
    src_leg.SetTextSize(0.034)
    dummy_line = ROOT.TLine()
    dummy_line.SetLineColor(ROOT.kBlack)
    dummy_line.SetLineWidth(2)
    dummy_marker = ROOT.TGraph()
    dummy_marker.SetMarkerColor(ROOT.kBlack)
    dummy_marker.SetMarkerStyle(20)
    dummy_marker.SetMarkerSize(0.9)
    dummy_ppg_marker = ROOT.TGraph()
    dummy_ppg_marker.SetMarkerColor(ROOT.kBlack)
    dummy_ppg_marker.SetMarkerStyle(24)
    dummy_ppg_marker.SetMarkerSize(0.9)
    if source_kind == "analysis_note_fig28" and current_object_family == "macro_raw":
        src_leg.AddEntry(dummy_line, "PPG12 IAN Fig.28", "l")
        src_leg.AddEntry(dummy_marker, "Current all-reco Fig.28 family", "p")
    elif source_kind == "persisted_nom_combined" and current_object_family == "signal":
        src_leg.AddEntry(dummy_ppg_marker, "PPG12 correction input (SI+DI)", "p")
        src_leg.AddEntry(dummy_marker, "Current correction input (SI+DI)", "p")
    elif source_kind == "persisted_all" and current_object_family == "signal":
        src_leg.AddEntry(dummy_line, "PPG12 combined 0+1.5 mrad", "l")
        src_leg.AddEntry(dummy_marker, "Current combined 0+1.5 mrad SI+DI", "p")
    elif source_kind == "persisted_final" and current_object_family == "signal":
        src_leg.AddEntry(dummy_line, "PPG12 nominal correction input", "l")
        src_leg.AddEntry(dummy_marker, "Current combined 0+1.5 mrad SI+DI", "p")
    else:
        src_leg.AddEntry(dummy_line, "PPG12 SDCC", "l")
        src_leg.AddEntry(dummy_marker, "Current output", "p")
    src_leg.Draw()
    top.RedrawAxis()

    bot.cd()
    frame_bot = ROOT.TH1F("frame_fig28_overlay_ratio", "", 43, 7.0, 50.0)
    frame_bot.SetStats(False)
    frame_bot.GetXaxis().SetRangeUser(10.0, 35.0)
    ratio_y_min, ratio_y_max = ratio_axis_bounds(ratio_rows_by_region)
    frame_bot.GetYaxis().SetRangeUser(ratio_y_min, ratio_y_max)
    frame_bot.GetXaxis().SetTitle("#it{E}_{T}^{#gamma,rec} [GeV]")
    frame_bot.GetYaxis().SetTitle(
        "IAN Fig.28 / Current"
        if source_kind == "analysis_note_fig28"
        else "SDCC / Current"
    )
    frame_bot.GetXaxis().SetTitleSize(0.095)
    frame_bot.GetYaxis().SetTitleSize(0.080)
    frame_bot.GetXaxis().SetLabelSize(0.080)
    frame_bot.GetYaxis().SetLabelSize(0.070)
    frame_bot.GetXaxis().SetTitleOffset(1.0)
    frame_bot.GetYaxis().SetTitleOffset(0.70)
    frame_bot.GetXaxis().SetNdivisions(505)
    frame_bot.GetYaxis().SetNdivisions(505)
    frame_bot.Draw("axis")
    line = ROOT.TLine(10.0, 1.0, 35.0, 1.0)
    line.SetLineColor(ROOT.kGray + 2)
    line.SetLineStyle(7)
    line.SetLineWidth(2)
    line.Draw()
    for spec in series:
        ratio_graphs[spec["region"]].Draw("same p")
    bot.RedrawAxis()

    output_png.parent.mkdir(parents=True, exist_ok=True)
    canvas.SaveAs(str(output_png))
    return {
        "ratio_table": ratio_table,
        "ratio_y_min": ratio_y_min,
        "ratio_y_max": ratio_y_max,
    }


def write_points_csv(rows: list[dict[str, object]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "region",
        "x_low_gev",
        "x_high_gev",
        "x_center_gev",
        "ppg12_sdcc",
        "ppg12_sdcc_error",
        "current_output",
        "current_output_error",
        "sdcc_over_current",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ppg12-csv", type=Path, default=DEFAULT_PPG12_CSV)
    parser.add_argument(
        "--ppg12-root",
        type=Path,
        help="read canonical h_leak_B/C/D directly from an audited PPG12 ROOT",
    )
    parser.add_argument(
        "--ppg12-analysis-note-csv",
        type=Path,
        help="read preserved historical PPG12 IAN Fig.28 B/C/D values",
    )
    parser.add_argument("--current-json", type=Path, default=DEFAULT_CURRENT_JSON)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--source-kind", default="macro_recomputed")
    parser.add_argument(
        "--current-object-family",
        choices=sorted(CURRENT_OBJECT_FAMILIES),
        default="macro_raw",
        help="current ROOT object family to compare against the selected PPG12 source",
    )
    args = parser.parse_args()

    current_root = resolve_current_root(args.current_json)
    current_den, series = build_series(args.current_object_family)
    if args.ppg12_analysis_note_csv:
        ppg12 = read_analysis_note_rows(args.ppg12_analysis_note_csv, series)
    elif args.ppg12_root:
        ppg12 = read_ppg12_root(args.ppg12_root, series)
    else:
        ppg12 = read_ppg12_rows(args.ppg12_csv, args.source_kind, series)
    current = current_points(current_root, current_den, series)

    stem = (
        "ppg12_fig28_signal_leakage_sdcc_vs_current"
        f"_{args.source_kind}_{args.current_object_family}"
        "_overlay_ratio"
    )
    png = args.outdir / f"{stem}.png"
    payload = render(
        ppg12,
        current,
        png,
        series,
        args.source_kind,
        args.current_object_family,
    )
    csv_path = png.with_name(png.stem + "_points.csv")
    json_path = png.with_name(png.stem + "_manifest.json")
    write_points_csv(payload["ratio_table"], csv_path)

    manifest = {
        "png": str(png),
        "points_csv": str(csv_path),
        "ppg12_csv": str(args.ppg12_csv),
        "ppg12_input_root": str(args.ppg12_root) if args.ppg12_root else None,
        "ppg12_analysis_note_csv": (
            str(args.ppg12_analysis_note_csv)
            if args.ppg12_analysis_note_csv
            else None
        ),
        "ppg12_source_kind": args.source_kind,
        "ppg12_source_root": PPG12_SOURCE_ROOTS.get(args.source_kind, "unknown"),
        "current_json": str(args.current_json),
        "current_root": str(current_root),
        "current_object_family": args.current_object_family,
        "current_object_family_description": CURRENT_OBJECT_FAMILIES[args.current_object_family]["description"],
        "current_denominator": current_den,
        "current_numerators": {spec["region"]: spec["current_num"] for spec in series},
        "ratio_definition": (
            "PPG12 IAN Fig.28 / Current output"
            if args.source_kind == "analysis_note_fig28"
            else "PPG12 SDCC / Current output"
        ),
        "ratio_y_range": [payload["ratio_y_min"], payload["ratio_y_max"]],
        "top_panel": (
            "PPG12 IAN Fig.28 as solid step curves; current all-reconstructed-candidate "
            "Fig.28 display-macro family as filled markers; not the purity-correction inputs"
            if args.source_kind == "analysis_note_fig28"
            else (
                "PPG12 SDCC as open markers; Current output as filled markers"
                if args.source_kind == "persisted_nom_combined"
                else "PPG12 SDCC as solid step curves; Current output as filled markers"
            )
        ),
        "comparison_scope": (
            "Published PPG12 IAN Fig.28 leakage values versus the current "
            "unsuffixed all-reconstructed-candidate photon-MC ABCD ratios used for "
            "Fig.28 display-macro parity; this is not a comparison of the truth-matched "
            "leakage factors used by the purity correction"
            if args.source_kind == "analysis_note_fig28"
            and args.current_object_family == "macro_raw"
            else (
                "PPG12 persisted purity-correction leakage factors versus the current "
                "fully combined, preweighted 0+1.5 mrad SI+DI photon sample"
                if args.source_kind
                in {"persisted_final", "persisted_nom_combined", "persisted_all"}
                and args.current_object_family == "signal"
                else "PPG12 and current leakage comparison"
            )
        ),
        "canonical_for_combined_nominal_parity": (
            args.source_kind == "persisted_nom_combined"
            and args.current_object_family == "signal"
        ),
        "canonical_for_published_fig28_parity": (
            args.source_kind == "analysis_note_fig28"
            and args.current_object_family == "macro_raw"
        ),
        "current_values_used_by_purity_correction": (
            args.current_object_family == "signal"
        ),
        "reference_status": (
            "preserved historical PPG12 IAN Fig.28 reference"
            if args.source_kind == "analysis_note_fig28"
            else (
                "current June combined-nominal PPG12 correction product"
                if args.source_kind == "persisted_nom_combined"
                else (
                    "legacy April bdt_all product; not canonical for current combined-nominal parity"
                    if args.source_kind == "persisted_all"
                    else "diagnostic or alternate PPG12 reference"
                )
            )
        ),
    }
    json_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(png)
    print(csv_path)
    print(json_path)


if __name__ == "__main__":
    main()
