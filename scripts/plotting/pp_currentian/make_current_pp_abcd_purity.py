#!/usr/bin/env python3
"""Plot current pp raw and leakage-corrected ABCD photon purity.

This is a local plotting helper for the recovered THE-74 pp table-QA root.
It mirrors the RecoilJets ABCD convention:
  A = isolated & tight
  B = non-isolated & tight
  C = isolated & non-tight
  D = non-isolated & non-tight
and uses the same fixed-point leakage solver implemented in
macros/AnalyzeRecoilJets.h.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import ROOT


DEFAULT_BASE = Path(
    "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/ppg12TableQA/"
    "THE42_ppg12_tableqa_v1_basev3e_20260611"
)
DEFAULT_DATA_ROOT = (
    DEFAULT_BASE
    / "merged_roots/RecoilJets_pp_ALL_preselectionNewPPG12_tightNewPPG12_nonTightNewPPG12.root"
)
DEFAULT_SIGNAL_ROOT = DEFAULT_BASE / "merged_roots/RecoilJets_photonjet5plus10plus20_MERGED.root"
DEFAULT_OUTDIR = DEFAULT_BASE / "purity_current_pp"

TRIGGER_DIR = "Photon_4_GeV_plus_MBD_NS_geq_1"
SIM_DIR = "SIM"
ISO_TOKEN = "isoR40_fixedIso2GeV"
PT_BINS = [
    (5, 8),
    (8, 10),
    (10, 12),
    (12, 14),
    (14, 16),
    (16, 18),
    (18, 20),
    (20, 22),
    (22, 24),
    (24, 26),
    (26, 35),
]


@dataclass(frozen=True)
class Point:
    pt_lo: float
    pt_hi: float
    a: float
    b: float
    c: float
    d: float
    raw: float
    raw_err: float
    corrected: float
    corrected_err: float
    lead_corrected: float
    f_b: float
    f_c: float
    f_d: float
    lead_f_b: float
    lead_f_c: float
    lead_f_d: float
    correction_ok: bool

    @property
    def x(self) -> float:
        return 0.5 * (self.pt_lo + self.pt_hi)

    @property
    def ex(self) -> float:
        return 0.5 * (self.pt_hi - self.pt_lo)


def read_bin1(directory: ROOT.TDirectory, name: str) -> float:
    hist = directory.Get(name)
    if not hist:
        return float("nan")
    return float(hist.GetBinContent(1))


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
    if a > 0.0:
        var += d_p_da * d_p_da * a
    if b > 0.0:
        var += d_p_db * d_p_db * b
    if c > 0.0:
        var += d_p_dc * d_p_dc * c
    if d > 0.0:
        var += d_p_dd * d_p_dd * d
    return math.sqrt(var) if var > 0.0 else 0.0


def solve_leakage_corrected_sa(
    a: float, b: float, c: float, d: float, f_b: float, f_c: float, f_d: float
) -> tuple[float, bool]:
    """Python copy of ARJ::SolveLeakageCorrectedSA."""
    if a <= 0.0:
        return 0.0, True

    s = a
    if d != 0.0:
        s = min(max(a - b * (c / d), 0.0), a)

    def fixed_point(value: float) -> float:
        denom = d - f_d * value
        if denom == 0.0:
            return float("nan")
        return a - (b - f_b * value) * (c - f_c * value) / denom

    lam = 0.25
    for iteration in range(200):
        if f_d > 0.0:
            s_max = (d / f_d) * 0.999
            if math.isfinite(s_max):
                s = min(s, max(0.0, s_max))

        f_val = fixed_point(s)
        if not math.isfinite(f_val):
            return s, False

        s_new = (1.0 - lam) * s + lam * f_val
        if not math.isfinite(s_new):
            return s, False

        s_new = min(max(s_new, 0.0), a)
        delta = abs(s_new - s)
        s = s_new
        if delta < 1.0e-6:
            return s, True
        if iteration > 10 and delta > 0.5 * a and lam > 0.05:
            lam *= 0.5

    return s, False


def corrected_purity_value(
    a: float, b: float, c: float, d: float, f_b: float, f_c: float, f_d: float
) -> tuple[float, bool]:
    sa, ok = solve_leakage_corrected_sa(a, b, c, d, f_b, f_c, f_d)
    if ok and a > 0.0:
        return sa / a, True
    return raw_purity(a, b, c, d), False


def corrected_purity_error(
    a: float, b: float, c: float, d: float, f_b: float, f_c: float, f_d: float, has_correction: bool
) -> float:
    if not has_correction:
        return raw_purity_error(a, b, c, d)
    if a <= 0.0:
        return 0.0

    widths = [math.sqrt(max(v, 1.0)) for v in (a, b, c, d)]
    lows = [max(0.0, v - w) for v, w in zip((a, b, c, d), widths)]
    highs = [v + w for v, w in zip((a, b, c, d), widths)]

    def value(vals: Iterable[float]) -> float:
        aa, bb, cc, dd = vals
        pur, _ = corrected_purity_value(aa, bb, cc, dd, f_b, f_c, f_d)
        return pur

    derivs = []
    for idx in range(4):
        up = [a, b, c, d]
        down = [a, b, c, d]
        up[idx] = highs[idx]
        down[idx] = lows[idx]
        denom = up[idx] - down[idx]
        derivs.append((value(up) - value(down)) / denom if denom > 0.0 else 0.0)

    var = 0.0
    for derivative, count in zip(derivs, (a, b, c, d)):
        if count > 0.0:
            var += derivative * derivative * count
    return math.sqrt(var) if var > 0.0 else 0.0


def leakage_fractions(sim_dir: ROOT.TDirectory, prefix: str, suffix: str) -> tuple[float, float, float]:
    hist = sim_dir.Get(prefix + suffix)
    if not hist:
        return 0.0, 0.0, 0.0
    a_sig = float(hist.GetBinContent(1))
    if a_sig <= 0.0:
        return 0.0, 0.0, 0.0
    return (
        float(hist.GetBinContent(2)) / a_sig,
        float(hist.GetBinContent(3)) / a_sig,
        float(hist.GetBinContent(4)) / a_sig,
    )


def build_points(data_root: Path, signal_root: Path) -> list[Point]:
    data_file = ROOT.TFile.Open(str(data_root))
    if not data_file or data_file.IsZombie():
        raise RuntimeError(f"Could not open data ROOT: {data_root}")
    signal_file = ROOT.TFile.Open(str(signal_root))
    if not signal_file or signal_file.IsZombie():
        raise RuntimeError(f"Could not open signal ROOT: {signal_root}")

    data_dir = data_file.Get(TRIGGER_DIR)
    sim_dir = signal_file.Get(SIM_DIR)
    if not data_dir:
        raise RuntimeError(f"Missing data directory {TRIGGER_DIR} in {data_root}")
    if not sim_dir:
        raise RuntimeError(f"Missing SIM directory in {signal_root}")

    points: list[Point] = []
    for pt_lo, pt_hi in PT_BINS:
        suffix = f"_{ISO_TOKEN}_pT_{pt_lo}_{pt_hi}"
        a = read_bin1(data_dir, "h_isIsolated_isTight" + suffix)
        b = read_bin1(data_dir, "h_notIsolated_isTight" + suffix)
        c = read_bin1(data_dir, "h_isIsolated_notTight" + suffix)
        d = read_bin1(data_dir, "h_notIsolated_notTight" + suffix)
        raw = raw_purity(a, b, c, d)
        raw_err = raw_purity_error(a, b, c, d)

        f_b, f_c, f_d = leakage_fractions(sim_dir, "h_sigABCD_MC", suffix)
        corrected, ok = corrected_purity_value(a, b, c, d, f_b, f_c, f_d)
        corrected_err = corrected_purity_error(a, b, c, d, f_b, f_c, f_d, ok)

        lead_f_b, lead_f_c, lead_f_d = leakage_fractions(sim_dir, "h_xJpurityLead_sigABCD_MC", suffix)
        lead_corrected, _ = corrected_purity_value(a, b, c, d, lead_f_b, lead_f_c, lead_f_d)

        points.append(
            Point(
                pt_lo=pt_lo,
                pt_hi=pt_hi,
                a=a,
                b=b,
                c=c,
                d=d,
                raw=raw,
                raw_err=raw_err,
                corrected=corrected,
                corrected_err=corrected_err,
                lead_corrected=lead_corrected,
                f_b=f_b,
                f_c=f_c,
                f_d=f_d,
                lead_f_b=lead_f_b,
                lead_f_c=lead_f_c,
                lead_f_d=lead_f_d,
                correction_ok=ok,
            )
        )

    data_file.Close()
    signal_file.Close()
    return points


def write_csv(points: list[Point], path: Path) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(Point.__dataclass_fields__.keys()))
        writer.writeheader()
        for point in points:
            writer.writerow(point.__dict__)


def make_plot(points: list[Point], output: Path, max_plot_hi: float = 26.0) -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.linewidth": 1.4,
            "xtick.major.width": 1.2,
            "ytick.major.width": 1.2,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )

    visible_points = [p for p in points if p.pt_hi <= max_plot_hi]
    hidden_points = [p for p in points if p.pt_hi > max_plot_hi]

    x = [p.x for p in visible_points]
    ex = [p.ex for p in visible_points]
    y_raw = [p.raw for p in visible_points]
    y_corr = [p.corrected for p in visible_points]
    err_raw = [p.raw_err for p in visible_points]
    err_corr = [p.corrected_err for p in visible_points]
    fig, axes = plt.subplots(1, 2, figsize=(13.6, 5.8), sharex=True, sharey=True)
    panel_specs = [
        (axes[0], y_raw, err_raw, "Raw ABCD purity", "black", "black", "black", "full"),
        (axes[1], y_corr, err_corr, "Leakage-corrected purity", "#d62728", "white", "#d62728", "none"),
    ]

    for ax, yvals, yerrs, title, color, face, edge, _fill in panel_specs:
        ax.errorbar(
            x,
            yvals,
            xerr=ex,
            yerr=yerrs,
            linestyle="None",
            marker="o",
            markersize=9.0,
            markerfacecolor=face,
            markeredgecolor=edge,
            markeredgewidth=2.0 if face == "white" else 1.4,
            ecolor=color,
            elinewidth=1.25,
            capsize=2.5,
        )
        ax.set_title(title, fontsize=18, fontweight="bold", pad=10)
        ax.set_ylim(0.0, 1.05)
        ax.set_xlim(4.6, max_plot_hi + 0.6)
        ax.set_xlabel(r"$p_T^\gamma$ [GeV]", fontsize=17)
        ax.tick_params(labelsize=13, length=6)
        ax.grid(axis="y", color="0.88", linewidth=0.9)
    axes[0].set_ylabel("ABCD photon purity", fontsize=17)

    header = (
        r"$\bf{\it{sPHENIX}}$ Internal   "
        + r"$p{+}p$ $\sqrt{s}=200$ GeV   "
        + r"Photon 4 GeV trigger, $|\eta|<0.7$   "
        + r"iso $R=0.4$, $E_T^{iso}<2$ GeV"
    )
    fig.text(0.5, 0.965, header, ha="center", va="top", fontsize=16.5)
    fig.text(
        0.5,
        0.913,
        "Current complete pp table-QA data; same binning and axis scale in both panels",
        ha="center",
        va="top",
        fontsize=13.8,
        color="0.25",
    )
    fig.subplots_adjust(left=0.08, right=0.985, top=0.82, bottom=0.17, wspace=0.08)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-root", type=Path, default=DEFAULT_DATA_ROOT)
    parser.add_argument("--signal-root", type=Path, default=DEFAULT_SIGNAL_ROOT)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    points = build_points(args.data_root, args.signal_root)

    png_path = args.outdir / "current_pp_abcd_purity_raw_vs_leakage_corrected.png"
    csv_path = args.outdir / "current_pp_abcd_purity_raw_vs_leakage_corrected_points.csv"
    manifest_path = args.outdir / "current_pp_abcd_purity_raw_vs_leakage_corrected_manifest.json"

    write_csv(points, csv_path)
    make_plot(points, png_path)
    manifest = {
        "plot": str(png_path),
        "points_csv": str(csv_path),
        "data_root": str(args.data_root),
        "signal_root_for_leakage": str(args.signal_root),
        "data_directory": TRIGGER_DIR,
        "sim_directory": SIM_DIR,
        "iso_token": ISO_TOKEN,
        "formula": {
            "raw": "max(A - B*C/D, 0) / A",
            "corrected": "SolveLeakageCorrectedSA(A,B,C,D,fB,fC,fD) / A",
            "leakage_fractions": "candidate-level h_sigABCD_MC bins B,C,D divided by bin A",
        },
        "cross_check": "event-leading h_xJpurityLead_sigABCD_MC corrected purity is written to CSV as lead_corrected",
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"[WROTE] {png_path}")
    print(f"[WROTE] {csv_path}")
    print(f"[WROTE] {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
