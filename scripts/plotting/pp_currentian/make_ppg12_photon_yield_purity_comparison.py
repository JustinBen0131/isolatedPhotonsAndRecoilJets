#!/usr/bin/env python3
"""Make a PPG12 Fig. 29 style purity comparison from photon-yield outputs."""

from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image, ImageDraw, ImageFont
import ROOT


ROOT.gROOT.SetBatch(True)

REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
CFG_TAG = "jetMinPtScan_dphiScan_vz60_isoR40_isSliding_preselectionNewPPG12_tightNewPPG12_nonTightNewPPG12"
DEFAULT_DATA_ROOT = (
    REPO
    / "InputFiles/pp24/ppg12_photon_yield_v1_data_20260620/pp"
    / f"RecoilJets_pp_ALL_{CFG_TAG}.root"
)
DEFAULT_SIGNAL_ROOT = (
    REPO
    / "dataOutput/ppg12PhotonYield/THE76_ppg12_photon_yield_v1_sim_20260616/merged_roots"
    / CFG_TAG
    / "photonJet5and10and20merged_SIM/RecoilJets_photonjet5plus10plus20_MERGED.root"
)
DEFAULT_OUTDIR = REPO / "dataOutput/ppg12PhotonYield/ppg12_photon_yield_v1_data_20260620/purity_fig29_comparison"
IAN_PDF = REPO / "usefulDocs/PPG12_analysis_note_2026-05-21_v4_current_IAN.pdf"
PDFTOPPM = Path("/Users/patsfan753/.cache/codex-runtimes/codex-primary-runtime/dependencies/bin/pdftoppm")
TRIGGER_DIR = "Photon_4_GeV_plus_MBD_NS_geq_1"
SIM_DIR = "SIM"
PT_EDGES = [10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36]


@dataclass(frozen=True)
class PurityPoint:
    pt_lo: float
    pt_hi: float
    a: float
    b: float
    c: float
    d: float
    sig_a: float
    sig_b: float
    sig_c: float
    sig_d: float
    f_b: float
    f_c: float
    f_d: float
    raw: float
    raw_err: float
    corrected: float
    corrected_err: float
    correction_ok: bool

    @property
    def x(self) -> float:
        return 0.5 * (self.pt_lo + self.pt_hi)

    @property
    def ex(self) -> float:
        return 0.5 * (self.pt_hi - self.pt_lo)


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    candidates = [
        "/System/Library/Fonts/Supplemental/Times New Roman Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Times New Roman.ttf",
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/Library/Fonts/Arial Bold.ttf" if bold else "/Library/Fonts/Arial.ttf",
    ]
    for candidate in candidates:
        if Path(candidate).exists():
            return ImageFont.truetype(candidate, size=size)
    return ImageFont.load_default()


def open_root(path: Path) -> ROOT.TFile:
    f = ROOT.TFile.Open(str(path))
    if not f or f.IsZombie():
        raise RuntimeError(f"could not open ROOT file: {path}")
    return f


def require_hist(directory: ROOT.TDirectory, name: str) -> ROOT.TH1:
    hist = directory.Get(name)
    if not hist:
        raise RuntimeError(f"missing histogram {name} in {directory.GetName()}")
    return hist


def axis_edges(hist: ROOT.TH1) -> list[float]:
    axis = hist.GetXaxis()
    return [float(axis.GetBinLowEdge(i)) for i in range(1, hist.GetNbinsX() + 1)] + [
        float(axis.GetBinUpEdge(hist.GetNbinsX()))
    ]


def check_edges(hist: ROOT.TH1, expected: list[int], name: str) -> None:
    got = axis_edges(hist)
    if len(got) != len(expected) or any(abs(a - b) > 1.0e-6 for a, b in zip(got, expected)):
        raise RuntimeError(f"{name} has edges {got}, expected {expected}")


def raw_purity(a: float, b: float, c: float, d: float) -> float:
    if a <= 0.0 or d <= 0.0:
        return float("nan")
    return max(a - b * c / d, 0.0) / a


def raw_purity_error(a: float, b: float, c: float, d: float) -> float:
    if a <= 0.0 or d <= 0.0:
        return 0.0
    d_p_da = (b * c) / (a * a * d)
    d_p_db = -c / (a * d)
    d_p_dc = -b / (a * d)
    d_p_dd = (b * c) / (a * d * d)
    var = 0.0
    for derivative, count in ((d_p_da, a), (d_p_db, b), (d_p_dc, c), (d_p_dd, d)):
        if count > 0.0:
            var += derivative * derivative * count
    return math.sqrt(max(var, 0.0))


def solve_leakage_corrected_sa(a: float, b: float, c: float, d: float, f_b: float, f_c: float, f_d: float) -> tuple[float, bool]:
    if a <= 0.0:
        return 0.0, True
    s = min(max(a - b * c / d, 0.0), a) if d != 0.0 else a

    def fixed_point(value: float) -> float:
        denom = d - f_d * value
        if denom == 0.0:
            return float("nan")
        return a - (b - f_b * value) * (c - f_c * value) / denom

    damping = 0.25
    for iteration in range(200):
        if f_d > 0.0:
            s_max = d / f_d * 0.999
            if math.isfinite(s_max):
                s = min(s, max(0.0, s_max))
        next_value = fixed_point(s)
        if not math.isfinite(next_value):
            return s, False
        s_new = (1.0 - damping) * s + damping * next_value
        if not math.isfinite(s_new):
            return s, False
        s_new = min(max(s_new, 0.0), a)
        delta = abs(s_new - s)
        s = s_new
        if delta < 1.0e-6:
            return s, True
        if iteration > 10 and delta > 0.5 * a and damping > 0.05:
            damping *= 0.5
    return s, False


def corrected_purity(a: float, b: float, c: float, d: float, f_b: float, f_c: float, f_d: float) -> tuple[float, bool]:
    s, ok = solve_leakage_corrected_sa(a, b, c, d, f_b, f_c, f_d)
    if ok and a > 0.0:
        return s / a, True
    return raw_purity(a, b, c, d), False


def corrected_purity_error(a: float, b: float, c: float, d: float, f_b: float, f_c: float, f_d: float, ok: bool) -> float:
    if not ok:
        return raw_purity_error(a, b, c, d)
    if a <= 0.0:
        return 0.0
    base = [a, b, c, d]
    widths = [math.sqrt(max(value, 1.0)) for value in base]

    def value(vals: list[float]) -> float:
        out, _ = corrected_purity(vals[0], vals[1], vals[2], vals[3], f_b, f_c, f_d)
        return out

    var = 0.0
    for idx, count in enumerate(base):
        lo = base.copy()
        hi = base.copy()
        lo[idx] = max(0.0, lo[idx] - widths[idx])
        hi[idx] += widths[idx]
        denom = hi[idx] - lo[idx]
        derivative = (value(hi) - value(lo)) / denom if denom > 0.0 else 0.0
        if count > 0.0:
            var += derivative * derivative * count
    return math.sqrt(max(var, 0.0))


def read_points(data_root: Path, signal_root: Path) -> list[PurityPoint]:
    data_file = open_root(data_root)
    signal_file = open_root(signal_root)
    data_dir = data_file.Get(TRIGGER_DIR)
    sig_dir = signal_file.Get(SIM_DIR)
    if not data_dir:
        raise RuntimeError(f"missing {TRIGGER_DIR} in {data_root}")
    if not sig_dir:
        raise RuntimeError(f"missing {SIM_DIR} in {signal_root}")

    h_a = require_hist(data_dir, "h_tight_iso_cluster_0")
    h_b = require_hist(data_dir, "h_tight_noniso_cluster_0")
    h_c = require_hist(data_dir, "h_nontight_iso_cluster_0")
    h_d = require_hist(data_dir, "h_nontight_noniso_cluster_0")
    h_common = require_hist(data_dir, "h_common_cluster_0")
    h_sig_a = require_hist(sig_dir, "h_tight_iso_cluster_signal_0")
    h_sig_b = require_hist(sig_dir, "h_tight_noniso_cluster_signal_0")
    h_sig_c = require_hist(sig_dir, "h_nontight_iso_cluster_signal_0")
    h_sig_d = require_hist(sig_dir, "h_nontight_noniso_cluster_signal_0")

    for name, hist in (
        ("h_tight_iso_cluster_0", h_a),
        ("h_tight_noniso_cluster_0", h_b),
        ("h_nontight_iso_cluster_0", h_c),
        ("h_nontight_noniso_cluster_0", h_d),
        ("h_common_cluster_0", h_common),
        ("h_tight_iso_cluster_signal_0", h_sig_a),
        ("h_tight_noniso_cluster_signal_0", h_sig_b),
        ("h_nontight_iso_cluster_signal_0", h_sig_c),
        ("h_nontight_noniso_cluster_signal_0", h_sig_d),
    ):
        check_edges(hist, PT_EDGES, name)

    points: list[PurityPoint] = []
    for idx, (lo, hi) in enumerate(zip(PT_EDGES[:-1], PT_EDGES[1:]), start=1):
        a = float(h_a.GetBinContent(idx))
        b = float(h_b.GetBinContent(idx))
        c = float(h_c.GetBinContent(idx))
        d = float(h_d.GetBinContent(idx))
        sig_a = float(h_sig_a.GetBinContent(idx))
        sig_b = float(h_sig_b.GetBinContent(idx))
        sig_c = float(h_sig_c.GetBinContent(idx))
        sig_d = float(h_sig_d.GetBinContent(idx))
        f_b = sig_b / sig_a if sig_a > 0.0 else 0.0
        f_c = sig_c / sig_a if sig_a > 0.0 else 0.0
        f_d = sig_d / sig_a if sig_a > 0.0 else 0.0
        raw = raw_purity(a, b, c, d)
        raw_err = raw_purity_error(a, b, c, d)
        corr, ok = corrected_purity(a, b, c, d, f_b, f_c, f_d)
        corr_err = corrected_purity_error(a, b, c, d, f_b, f_c, f_d, ok)
        points.append(
            PurityPoint(
                pt_lo=float(lo),
                pt_hi=float(hi),
                a=a,
                b=b,
                c=c,
                d=d,
                sig_a=sig_a,
                sig_b=sig_b,
                sig_c=sig_c,
                sig_d=sig_d,
                f_b=f_b,
                f_c=f_c,
                f_d=f_d,
                raw=raw,
                raw_err=raw_err,
                corrected=corr,
                corrected_err=corr_err,
                correction_ok=ok,
            )
        )
    data_file.Close()
    signal_file.Close()
    return points


def write_points(points: list[PurityPoint], out: Path) -> None:
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(PurityPoint.__dataclass_fields__.keys()))
        writer.writeheader()
        for point in points:
            writer.writerow(point.__dict__)


def render_current_panel(points: list[PurityPoint], out: Path) -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "mathtext.fontset": "dejavuserif",
            "axes.linewidth": 1.25,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )
    x = np.array([p.x for p in points])
    ex = np.array([p.ex for p in points])
    y_raw = np.array([p.raw for p in points])
    y_corr = np.array([p.corrected for p in points])
    ey_raw = np.array([p.raw_err for p in points])
    ey_corr = np.array([p.corrected_err for p in points])

    fig, ax = plt.subplots(figsize=(6.2, 5.6), dpi=240)
    fig.patch.set_facecolor("white")
    ax.errorbar(
        x,
        y_corr,
        xerr=ex,
        yerr=ey_corr,
        linestyle="None",
        marker="o",
        markersize=7.6,
        markerfacecolor="white",
        markeredgecolor="#2255ff",
        markeredgewidth=1.5,
        ecolor="#2255ff",
        elinewidth=1.1,
        capsize=2.5,
        label="w/ sig. leak. corr.",
    )
    ax.errorbar(
        x,
        y_raw,
        xerr=ex,
        yerr=ey_raw,
        linestyle="None",
        marker="o",
        markersize=7.2,
        markerfacecolor="black",
        markeredgecolor="black",
        ecolor="black",
        elinewidth=1.1,
        capsize=2.5,
        label="w/o sig. leak. corr.",
    )
    ax.set_xlim(10, 36)
    ax.set_ylim(0, 1.2)
    ax.set_xlabel(r"$E_T^{\gamma,\mathrm{rec}}$ [GeV]", fontsize=16)
    ax.set_ylabel("Purity", fontsize=17)
    ax.tick_params(which="both", labelsize=13, length=5)
    ax.minorticks_on()
    ax.text(0.08, 0.95, r"$\bf{\it{sPHENIX}}$ Internal", transform=ax.transAxes, ha="left", va="top", fontsize=14)
    ax.text(0.08, 0.885, r"$p{+}p$ $\sqrt{s}=200$ GeV", transform=ax.transAxes, ha="left", va="top", fontsize=13)
    ax.legend(loc="lower left", bbox_to_anchor=(0.07, 0.07), frameon=False, fontsize=12.2, handlelength=1.5)
    fig.subplots_adjust(left=0.16, right=0.98, top=0.98, bottom=0.15)
    fig.savefig(out)
    plt.close(fig)


def ensure_ian_crop(outdir: Path) -> Path:
    outdir.mkdir(parents=True, exist_ok=True)
    page = outdir / "ppg12_current_ian_purity_page-039.png"
    crop = outdir / "ppg12_current_ian_fig29_left_purity.png"
    if crop.exists() and crop.stat().st_size > 50_000:
        return crop
    if not page.exists() or page.stat().st_size < 100_000:
        prefix = outdir / "ppg12_current_ian_purity_page"
        subprocess.run(
            [str(PDFTOPPM), "-png", "-f", "39", "-l", "39", "-r", "220", str(IAN_PDF), str(prefix)],
            check=True,
        )
    Image.open(page).convert("RGB").crop((285, 245, 920, 885)).save(crop)
    return crop


def trim_white(img: Image.Image, tolerance: int = 248, pad: int = 8) -> Image.Image:
    rgb = img.convert("RGB")
    arr = np.asarray(rgb)
    mask = np.any(arr < tolerance, axis=2)
    if not mask.any():
        return rgb
    ys, xs = np.where(mask)
    return rgb.crop(
        (
            max(int(xs.min()) - pad, 0),
            max(int(ys.min()) - pad, 0),
            min(int(xs.max()) + pad + 1, rgb.width),
            min(int(ys.max()) + pad + 1, rgb.height),
        )
    )


def fit_image(img: Image.Image, size: tuple[int, int]) -> Image.Image:
    canvas = Image.new("RGB", size, "white")
    scale = min(size[0] / img.width, size[1] / img.height)
    fitted = img.resize((max(1, int(img.width * scale)), max(1, int(img.height * scale))), Image.Resampling.LANCZOS)
    canvas.paste(fitted, ((size[0] - fitted.width) // 2, (size[1] - fitted.height) // 2))
    return canvas


def compose_comparison(reference_path: Path, current_path: Path, out: Path) -> None:
    slide = Image.new("RGB", (2560, 1440), "white")
    draw = ImageDraw.Draw(slide)
    title_font = font(78, bold=True)
    subtitle_font = font(48)
    label_font = font(42, bold=True)
    note_font = font(38)
    title = "Photon Purity - pp Baseline Cross-Check"
    subtitle = "Signal-leakage correction with PPG12 photon-yield binning"
    draw.text((70, 42), title, font=title_font, fill=(0, 0, 0))
    draw.text((70, 150), subtitle, font=subtitle_font, fill=(25, 35, 52))

    card_size = (940, 800)
    left_xy = (220, 292)
    right_xy = (1400, 292)
    ref = fit_image(trim_white(Image.open(reference_path).convert("RGB")), card_size)
    cur = fit_image(trim_white(Image.open(current_path).convert("RGB")), card_size)
    slide.paste(ref, left_xy)
    slide.paste(cur, right_xy)

    labels = [("PPG12 IAN Fig. 29", left_xy), ("Current pp data", right_xy)]
    for label, xy in labels:
        x = xy[0] + (card_size[0] - int(draw.textlength(label, font=label_font))) // 2
        draw.text((x, xy[1] + card_size[1] + 14), label, font=label_font, fill=(0, 0, 0))

    notes = [
        "Data ROOT: full pp GRL merge, 1564 runs / 29216 jobs, exact 10-36 GeV reco bins.",
        "Leakage correction uses the matching pp signal MC photon-yield contract output.",
    ]
    y = 1212
    for text in notes:
        draw.polygon([(88, y + 9), (118, y + 24), (88, y + 39)], fill=(35, 95, 165))
        draw.text((142, y), text, font=note_font, fill=(25, 35, 52))
        y += 62
    slide.save(out)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", type=Path, default=DEFAULT_DATA_ROOT)
    ap.add_argument("--signal-root", type=Path, default=DEFAULT_SIGNAL_ROOT)
    ap.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    points = read_points(args.data_root, args.signal_root)
    csv_path = args.outdir / "current_pp_photon_yield_purity_points.csv"
    panel_path = args.outdir / "current_pp_photon_yield_purity_panel.png"
    comparison_path = args.outdir / "current_pp_vs_ppg12_ian_fig29_purity_comparison.png"
    write_points(points, csv_path)
    render_current_panel(points, panel_path)
    reference_path = ensure_ian_crop(args.outdir)
    compose_comparison(reference_path, panel_path, comparison_path)
    manifest = {
        "data_root": str(args.data_root),
        "signal_root": str(args.signal_root),
        "points_csv": str(csv_path),
        "current_panel_png": str(panel_path),
        "comparison_png": str(comparison_path),
        "ian_reference_pdf": str(IAN_PDF),
        "ian_reference_page": 39,
        "ian_reference_crop": str(reference_path),
        "pt_edges": PT_EDGES,
        "method": "Raw ABCD purity = max(A - B*C/D, 0)/A; corrected purity solves the same signal-leakage fixed point used by AnalyzeRecoilJets with signal MC fractions Bsig/Asig, Csig/Asig, Dsig/Asig per pT bin.",
        "audience_label_note": "No THE-76 label is drawn on the PNG.",
    }
    manifest_path = args.outdir / "current_pp_vs_ppg12_ian_fig29_purity_comparison_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(comparison_path)
    print(panel_path)
    print(csv_path)
    print(manifest_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
