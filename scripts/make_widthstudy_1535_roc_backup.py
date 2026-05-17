#!/usr/bin/env python3
"""Make backup ROC-curve table for the 15-35 GeV width-input BDT study."""

from __future__ import annotations

import csv
import math
from array import array
from pathlib import Path

import numpy as np
import ROOT


REPORT_DIR = Path("dataOutput/auauTightBDTValidation/model_validation_condor_20260511_171943")
OUT_DIR = Path("dataOutput/jstg_slide_candidates/slide10_widthstudy_1535")

PRODUCTS = [
    ("centAsFeat_pt15to35", "Base widths", ROOT.kGray + 2, 1),
    ("centAsFeat3x3_pt15to35", "3x3 widths", ROOT.kOrange + 7, 2),
    ("centAsFeatBase3x3_pt15to35", "Base + 3x3", ROOT.kGreen + 3, 1),
]

CENT_BINS = [
    ("0_20", 0.0, 20.0, "0-20% central"),
    ("20_50", 20.0, 50.0, "20-50% mid-central"),
    ("50_80", 50.0, 80.0, "50-80% peripheral"),
]

PT_BINS = [
    ("15_20", 15.0, 20.0, "15-20 GeV"),
    ("20_25", 20.0, 25.0, "20-25 GeV"),
    ("25_35", 25.0, 35.0, "25-35 GeV"),
]


def load_score_caches(cache_dir: Path) -> dict[str, np.ndarray]:
    parts: dict[str, list[np.ndarray]] = {
        "is_signal": [],
        "cluster_Et": [],
        "centrality": [],
    }
    for product, _, _, _ in PRODUCTS:
        parts[f"score_{product}"] = []

    for path in sorted(cache_dir.glob("score_cache_*.npz")):
        with np.load(path, allow_pickle=True) as zf:
            for key in parts:
                if key in zf.files:
                    parts[key].append(np.asarray(zf[key]))

    merged = {}
    for key, arrays in parts.items():
        if not arrays:
            raise RuntimeError(f"No arrays found for {key} in {cache_dir}")
        merged[key] = np.concatenate(arrays)
    return merged


def read_auc_table(path: Path) -> dict[tuple[str, str, str], float]:
    out: dict[tuple[str, str, str], float] = {}
    with path.open() as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            out[(row["product"], row["centrality_bin"], row["pt_bin"])] = float(row["auc"])
    return out


def roc_curve_unweighted(y_true: np.ndarray, score: np.ndarray):
    mask = np.isfinite(score) & np.isin(y_true, [0, 1])
    y = np.asarray(y_true[mask], dtype=np.int8)
    s = np.asarray(score[mask], dtype=np.float64)
    if y.size < 2 or np.unique(y).size < 2:
        return math.nan, np.array([]), np.array([]), np.array([])

    order = np.argsort(s, kind="mergesort")[::-1]
    y = y[order]
    s = s[order]

    distinct = np.where(np.diff(s))[0]
    threshold_idxs = np.r_[distinct, y.size - 1]
    tps = np.cumsum(y, dtype=np.float64)[threshold_idxs]
    fps = 1.0 + threshold_idxs.astype(np.float64) - tps
    thresholds = s[threshold_idxs]

    tps = np.r_[0.0, tps]
    fps = np.r_[0.0, fps]
    thresholds = np.r_[np.inf, thresholds]

    if tps[-1] <= 0 or fps[-1] <= 0:
        return math.nan, np.array([]), np.array([]), np.array([])
    tpr = tps / tps[-1]
    fpr = fps / fps[-1]
    auc = float(np.trapz(tpr, fpr))
    return auc, fpr, tpr, thresholds


def graph_from_curve(fpr: np.ndarray, tpr: np.ndarray, color: int, style: int) -> ROOT.TGraph:
    graph = ROOT.TGraph(
        len(fpr),
        array("d", np.asarray(fpr, dtype=np.float64)),
        array("d", np.asarray(tpr, dtype=np.float64)),
    )
    graph.SetLineColor(color)
    graph.SetLineWidth(2)
    graph.SetLineStyle(style)
    return graph


def draw_sphenix_label(x: float, y: float):
    text = ROOT.TLatex()
    text.SetNDC()
    text.SetTextAlign(11)
    text.SetTextFont(42)
    text.SetTextSize(0.024)
    text.DrawLatex(x, y, "#it{#bf{sPHENIX}} Internal")
    text.SetTextSize(0.018)
    text.DrawLatex(x, y - 0.030, "Pythia overlay, #sqrt{s_{NN}} = 200 GeV")


def main() -> None:
    ROOT.gROOT.SetBatch(True)
    ROOT.gROOT.LoadMacro("macros/sPhenixStyle.C")
    ROOT.SetsPhenixStyle()
    ROOT.gStyle.SetOptStat(0)

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    frame = load_score_caches(REPORT_DIR / "score_caches")
    y = frame["is_signal"].astype(np.int8)
    et = frame["cluster_Et"].astype(np.float64)
    cent = frame["centrality"].astype(np.float64)
    auc_table = read_auc_table(REPORT_DIR / "validation_auc_table.csv")

    curves = {}
    max_auc_diff = 0.0
    raw_csv = OUT_DIR / "width_inputs_1535_cell_roc_points.csv"
    auc_csv = OUT_DIR / "width_inputs_1535_cell_roc_auc_summary.csv"

    with raw_csv.open("w", newline="") as raw_handle, auc_csv.open("w", newline="") as auc_handle:
        raw_writer = csv.writer(raw_handle)
        raw_writer.writerow(
            ["product", "model_label", "centrality_bin", "pt_bin", "point_index", "background_fake_rate", "signal_efficiency", "threshold"]
        )
        auc_writer = csv.writer(auc_handle)
        auc_writer.writerow(["product", "model_label", "centrality_bin", "pt_bin", "auc_from_curve", "auc_from_validation_table", "entries", "signal_entries", "background_entries"])

        for product, label, color, style in PRODUCTS:
            score = frame[f"score_{product}"].astype(np.float64)
            for cent_key, clo, chi, cent_label in CENT_BINS:
                for pt_key, plo, phi, pt_label in PT_BINS:
                    mask = (
                        np.isfinite(score)
                        & np.isin(y, [0, 1])
                        & np.isfinite(cent)
                        & (cent >= clo)
                        & (cent < chi)
                        & np.isfinite(et)
                        & (et >= plo)
                        & (et < phi)
                    )
                    auc, fpr, tpr, thresholds = roc_curve_unweighted(y[mask], score[mask])
                    ref_auc = auc_table[(product, cent_key, pt_key)]
                    max_auc_diff = max(max_auc_diff, abs(auc - ref_auc))
                    curves[(product, cent_key, pt_key)] = {
                        "auc": auc,
                        "fpr": fpr,
                        "tpr": tpr,
                        "thresholds": thresholds,
                        "label": label,
                        "color": color,
                        "style": style,
                    }
                    auc_writer.writerow(
                        [
                            product,
                            label,
                            cent_key,
                            pt_key,
                            f"{auc:.10f}",
                            f"{ref_auc:.10f}",
                            int(mask.sum()),
                            int((mask & (y == 1)).sum()),
                            int((mask & (y == 0)).sum()),
                        ]
                    )
                    for idx, (x, yy, thr) in enumerate(zip(fpr, tpr, thresholds)):
                        raw_writer.writerow(
                            [
                                product,
                                label,
                                cent_key,
                                pt_key,
                                idx,
                                f"{float(x):.10g}",
                                f"{float(yy):.10g}",
                                "" if not np.isfinite(thr) else f"{float(thr):.10g}",
                            ]
                        )

    if max_auc_diff > 5.0e-7:
        raise RuntimeError(f"Regenerated AUC differs from validation table by {max_auc_diff:.3g}")

    canvas = ROOT.TCanvas("c_width_inputs_1535_roc_table", "", 2400, 1350)
    canvas.SetMargin(0, 0, 0, 0)

    title = ROOT.TLatex()
    title.SetNDC()
    title.SetTextFont(42)
    title.SetTextAlign(22)
    title.SetTextSize(0.030)
    title.DrawLatex(0.53, 0.973, "Backup: ROC curves behind the 15-35 GeV width-input AUC table")
    title.SetTextSize(0.019)
    title.DrawLatex(0.53, 0.940, "Same validation cells as the winner map; curves are raw ROC points from score caches.")
    draw_sphenix_label(0.055, 0.972)

    legend = ROOT.TLegend(0.510, 0.882, 0.965, 0.912)
    legend.SetNColumns(3)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.017)
    dummy_graphs = []
    for _, label, color, style in PRODUCTS:
        g = ROOT.TGraph()
        g.SetLineColor(color)
        g.SetLineWidth(3)
        g.SetLineStyle(style)
        dummy_graphs.append(g)
        legend.AddEntry(g, label, "l")
    legend.Draw()

    pads = []
    graphs = []
    frames = []
    diagonals = []

    x0, x1 = 0.060, 0.975
    y0, y1 = 0.075, 0.855
    gap_x, gap_y = 0.018, 0.028
    pad_w = (x1 - x0 - 2.0 * gap_x) / 3.0
    pad_h = (y1 - y0 - 2.0 * gap_y) / 3.0

    for ic, (cent_key, _, _, cent_label) in enumerate(CENT_BINS):
        for ip, (pt_key, _, _, pt_label) in enumerate(PT_BINS):
            px0 = x0 + ip * (pad_w + gap_x)
            px1 = px0 + pad_w
            py1 = y1 - ic * (pad_h + gap_y)
            py0 = py1 - pad_h
            pad = ROOT.TPad(f"pad_{ic}_{ip}", "", px0, py0, px1, py1)
            pad.SetLeftMargin(0.13)
            pad.SetRightMargin(0.035)
            pad.SetTopMargin(0.105)
            pad.SetBottomMargin(0.145)
            pad.Draw()
            pads.append(pad)
            pad.cd()
            pad.SetTicks(1, 1)

            frame_hist = ROOT.TH2D(
                f"frame_{ic}_{ip}",
                ";Background fake rate;Signal efficiency",
                10,
                0.0,
                1.0,
                10,
                0.0,
                1.0,
            )
            frame_hist.SetDirectory(0)
            frame_hist.SetStats(False)
            frame_hist.GetXaxis().SetNdivisions(505)
            frame_hist.GetYaxis().SetNdivisions(505)
            frame_hist.GetXaxis().SetTitleSize(0.050)
            frame_hist.GetXaxis().SetLabelSize(0.042)
            frame_hist.GetXaxis().SetTitleOffset(1.05)
            frame_hist.GetYaxis().SetTitleSize(0.050)
            frame_hist.GetYaxis().SetLabelSize(0.042)
            frame_hist.GetYaxis().SetTitleOffset(1.18)
            frame_hist.Draw("AXIS")
            frames.append(frame_hist)

            diag = ROOT.TLine(0.0, 0.0, 1.0, 1.0)
            diag.SetLineColor(ROOT.kGray + 1)
            diag.SetLineStyle(2)
            diag.SetLineWidth(1)
            diag.Draw()
            diagonals.append(diag)

            for product, _, _, _ in PRODUCTS:
                item = curves[(product, cent_key, pt_key)]
                graph = graph_from_curve(item["fpr"], item["tpr"], item["color"], item["style"])
                graph.Draw("L SAME")
                graphs.append(graph)

            text = ROOT.TLatex()
            text.SetNDC()
            text.SetTextFont(42)
            text.SetTextAlign(11)
            text.SetTextSize(0.045)
            text.DrawLatex(0.17, 0.915, f"{cent_label}, {pt_label}")
            text.SetTextSize(0.037)
            yy = 0.265
            for product, label, color, _ in PRODUCTS:
                auc = curves[(product, cent_key, pt_key)]["auc"]
                text.SetTextColor(color)
                text.DrawLatex(0.54, yy, f"{label}: AUC {auc:.3f}")
                yy -= 0.055
            text.SetTextColor(ROOT.kBlack)
            frame_hist.Draw("AXIS SAME")

            canvas.cd()

    png = OUT_DIR / "width_inputs_1535_cell_roc_curve_table.png"
    canvas.SaveAs(str(png))
    print(f"Wrote {png}")
    print(f"Wrote {auc_csv}")
    print(f"Wrote {raw_csv}")


if __name__ == "__main__":
    main()
