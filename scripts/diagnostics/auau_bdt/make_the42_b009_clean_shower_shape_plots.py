#!/usr/bin/env python3
"""Render clean THE-42 b009 shower-shape and NPB-tagged data diagnostics."""

from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


THIS_FILE = Path(__file__).resolve()
REPO = next((p for p in THIS_FILE.parents if (p / "AGENTS.md").exists()), THIS_FILE.parents[3])
BASE = REPO / "dataOutput/auauTightBDTValidation/THE42_wp80_centlinear_ss_20260604"
CFG = "preselectionNewPPG12_tightAuAuCentInputBase3x3BDT_nonTightAuAuBDTComplement_baseVariant"
DEFAULT_DATA = BASE / "b009_merged_data_roots_20260609" / f"RecoilJets_auau_ALL_{CFG}.root"
DEFAULT_SIGNAL = BASE / "sim_roots" / CFG / "photonJet12and20merged_SIM" / "RecoilJets_embeddedPhoton12plus20_MERGED.root"
DEFAULT_INCLUSIVE = (
    BASE
    / "sim_roots"
    / CFG
    / "embeddedJet12and20and30and40merged_SIM"
    / "RecoilJets_embeddedJet12plus20plus30plus40_MERGED.root"
)
DEFAULT_OUTDIR = BASE / "clean_shower_shape_b009_20260610"
PATH_INDEX: dict[Path, dict[str, list[str]]] = {}

PT_BINS = ((22, 24), (24, 26), (26, 28))
CENT_BINS = ((0, 10), (10, 20), (20, 30), (30, 40), (40, 50), (50, 60), (60, 80))
VAR_LABELS = {
    "e11e33": r"$E_{11}/E_{33}$",
    "e32e35": r"$E_{32}/E_{35}$",
    "et1": r"$E_{1}/E_{\mathrm{cluster}}$",
    "weta33": r"$w_{\eta,3\times3}$",
    "wphi33": r"$w_{\phi,3\times3}$",
    "npbScore": "NPB score",
}
XRANGES = {
    "e11e33": (0.0, 1.08),
    "e32e35": (0.0, 1.08),
    "et1": (0.0, 1.08),
    "weta33": (0.0, 0.35),
    "wphi33": (0.0, 0.35),
    "npbScore": (0.0, 1.0),
}


@dataclass
class Curve:
    label: str
    tag: str
    edges: np.ndarray
    centers: np.ndarray
    values: np.ndarray
    errors: np.ndarray
    raw_integral: float
    matches: list[str]


def parse_bins(spec: str) -> list[tuple[int, int]]:
    bins: list[tuple[int, int]] = []
    for chunk in spec.split(","):
        if not chunk.strip():
            continue
        lo_s, hi_s = chunk.split(":", 1)
        lo, hi = int(float(lo_s)), int(float(hi_s))
        if hi <= lo:
            raise ValueError(f"invalid bin {chunk!r}")
        bins.append((lo, hi))
    if not bins:
        raise ValueError(f"no bins parsed from {spec!r}")
    return bins


def open_root(path: Path):
    import ROOT

    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.AddDirectory(False)
    f = ROOT.TFile.Open(str(path), "READ")
    if not f or f.IsZombie():
        raise OSError(f"could not open ROOT file: {path}")
    return f


def walk(directory, prefix: str = ""):
    for key in directory.GetListOfKeys():
        obj = key.ReadObj()
        path = f"{prefix}/{key.GetName()}" if prefix else key.GetName()
        yield path, obj
        if obj.InheritsFrom("TDirectory"):
            yield from walk(obj, path)


def hist_name(var: str, tag: str, pt_bin: tuple[int, int], cent_bin: tuple[int, int]) -> str:
    pt_lo, pt_hi = pt_bin
    c_lo, c_hi = cent_bin
    return f"h_ss_{var}_{tag}_pT_{pt_lo}_{pt_hi}_cent_{c_lo}_{c_hi}"


def build_path_index(root_path: Path) -> dict[str, list[str]]:
    cached = PATH_INDEX.get(root_path)
    if cached is not None:
        return cached
    root_file = open_root(root_path)
    index: dict[str, list[str]] = {}
    try:
        for path, obj in walk(root_file):
            if not obj.InheritsFrom("TH1"):
                continue
            index.setdefault(Path(path).name, []).append(path)
    finally:
        root_file.Close()
    PATH_INDEX[root_path] = index
    return index


def find_exact_hist(root_file, target: str, paths: list[str]):
    for path in paths:
        obj = root_file.Get(path)
        if not obj.InheritsFrom("TH1"):
            continue
        h = obj.Clone(f"{target}_{abs(hash(path))}")
        h.SetDirectory(0)
        return path, h
    return None, None


def sum_hist(root_path: Path, var: str, tag: str, pt_bins, cent_bins) -> tuple[object, list[str]]:
    path_index = build_path_index(root_path)
    root_file = open_root(root_path)
    acc = None
    matches: list[str] = []
    missing = 0
    try:
        for pt_bin in pt_bins:
            for cent_bin in cent_bins:
                target = hist_name(var, tag, pt_bin, cent_bin)
                obj_path, hist = find_exact_hist(root_file, target, path_index.get(target, []))
                if hist is None:
                    missing += 1
                    continue
                matches.append(obj_path)
                if acc is None:
                    acc = hist.Clone(f"{var}_{tag}_sum")
                    acc.SetDirectory(0)
                else:
                    acc.Add(hist)
    finally:
        root_file.Close()
    if acc is None:
        raise RuntimeError(f"no histograms found for {root_path.name} var={var} tag={tag}")
    if missing:
        print(f"warning: {root_path.name} var={var} tag={tag} missing {missing} requested bins", file=sys.stderr)
    return acc, matches


def rebin(hist, factor: int):
    if factor <= 1:
        return hist
    nb = hist.GetNbinsX()
    if nb % factor != 0:
        return hist
    out = hist.Rebin(factor, f"{hist.GetName()}_rebin{factor}")
    out.SetDirectory(0)
    return out


def curve_from_hist(hist, matches: list[str], label: str, tag: str, x_min: float, x_max: float, rebin_factor: int) -> Curve:
    hist = rebin(hist, rebin_factor)
    nb = hist.GetNbinsX()
    edges = np.array([hist.GetXaxis().GetBinLowEdge(i) for i in range(1, nb + 2)], dtype=float)
    centers = np.array([hist.GetXaxis().GetBinCenter(i) for i in range(1, nb + 1)], dtype=float)
    counts = np.array([hist.GetBinContent(i) for i in range(1, nb + 1)], dtype=float)
    errs = np.array([hist.GetBinError(i) for i in range(1, nb + 1)], dtype=float)
    mask = (centers >= x_min) & (centers < x_max)
    idx = np.flatnonzero(mask)
    if len(idx):
        edges = np.concatenate(([edges[idx[0]]], edges[idx + 1]))
    centers = centers[mask]
    counts = counts[mask]
    errs = errs[mask]
    integral = float(np.sum(counts))
    values = np.zeros_like(counts)
    errors = np.zeros_like(errs)
    if integral > 0:
        values = counts / integral
        errors = errs / integral
    return Curve(label, tag, edges, centers, values, errors, integral, matches)


def load_curve(root_path: Path, var: str, tag: str, label: str, pt_bins, cent_bins, rebin_factor: int) -> Curve:
    x_min, x_max = XRANGES.get(var, (0.0, 1.0))
    hist, matches = sum_hist(root_path, var, tag, pt_bins, cent_bins)
    return curve_from_hist(hist, matches, label, tag, x_min, x_max, rebin_factor)


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "axes.edgecolor": "#111827",
            "axes.linewidth": 1.0,
            "xtick.color": "#111827",
            "ytick.color": "#111827",
        }
    )


def bin_range_label(pt_bins) -> str:
    lo = min(b[0] for b in pt_bins)
    hi = max(b[1] for b in pt_bins)
    return rf"${lo:g} \leq E_T^\gamma < {hi:g}$ GeV"


def draw_step(ax, curve: Curve, color: str, *, lw: float = 2.2, ls: str = "-", marker: str | None = None, alpha: float = 1.0):
    ax.step(curve.edges[:-1], curve.values, where="post", color=color, lw=lw, ls=ls, alpha=alpha, label=f"{curve.label} ({curve.raw_integral:.0f})")
    if marker:
        ax.errorbar(curve.centers, curve.values, yerr=curve.errors, fmt=marker, ms=3.2, color=color, lw=0.9, alpha=alpha)


def decorate(ax, var: str, title: str) -> None:
    ax.set_title(title, loc="left", fontsize=13.5, fontweight="bold")
    ax.set_xlabel(VAR_LABELS[var], fontsize=12.5)
    ax.set_ylabel("Unit-normalized candidates", fontsize=12.5)
    ax.grid(True, color="#E5E7EB", lw=0.7, alpha=0.85)
    ax.tick_params(labelsize=10.5)
    ax.text(0.035, 0.955, r"$\it{\bf{sPHENIX}}$ Internal", transform=ax.transAxes, ha="left", va="top", fontsize=10.5)


def render_ppg12_like_npb(data_root: Path, signal_root: Path, inclusive_root: Path, outdir: Path, pt_bins, cent_bins) -> tuple[Path, dict]:
    var = "e11e33"
    curves = [
        load_curve(data_root, var, "inclusive", "Au+Au data, no NPB cut", pt_bins, cent_bins, 2),
        load_curve(signal_root, var, "inclusive_sig", "Embedded direct photons", pt_bins, cent_bins, 2),
        load_curve(inclusive_root, var, "inclusive_bkg", "Embedded inclusive jets", pt_bins, cent_bins, 2),
        load_curve(data_root, var, "npbPass", "Au+Au data, NPB-tagged", pt_bins, cent_bins, 2),
    ]
    colors = ["#111827", "#C0262D", "#2563EB", "#16A34A"]

    fig, ax = plt.subplots(figsize=(9.2, 6.2), dpi=180)
    for curve, color in zip(curves, colors):
        draw_step(ax, curve, color, marker="o" if "data" in curve.label else None, lw=2.3 if "NPB" not in curve.label else 2.8)
    decorate(ax, var, "PPG12-style NPB-tagged data check")
    ax.text(
        0.035,
        0.875,
        rf"b009 v008 Au+Au, {bin_range_label(pt_bins)}, 0-80% centrality",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=10.5,
    )
    ax.text(
        0.035,
        0.815,
        "Green curve is explicit h_ss_e11e33_npbPass from the merged data ROOT.",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=10.2,
        color="#14532D",
    )
    ax.set_xlim(*XRANGES[var])
    ymax = max(float(np.max(c.values)) for c in curves if len(c.values))
    ax.set_ylim(0, ymax * 1.22 if ymax > 0 else 1)
    ax.legend(loc="upper center", bbox_to_anchor=(0.54, 1.01), ncols=2, frameon=False, fontsize=9.6, handlelength=2.3)
    fig.tight_layout(pad=1.0)
    out = outdir / "the42_b009_ppg12_style_npb_tagged_e11e33.png"
    fig.savefig(out)
    plt.close(fig)
    return out, {c.tag: {"label": c.label, "integral": c.raw_integral, "matches": len(c.matches)} for c in curves}


def render_tight_grid(data_root: Path, signal_root: Path, inclusive_root: Path, outdir: Path, pt_bins, cent_bins) -> tuple[Path, dict]:
    specs = [
        ("e11e33", "E11/E33: compact core fraction"),
        ("e32e35", "E32/E35: 3x2 over 3x5 energy"),
        ("et1", "E1/Ecluster: hottest-tower share"),
        ("weta33", "3x3 eta width"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(12.4, 8.6), dpi=180)
    summary: dict[str, dict] = {}
    for ax, (var, title) in zip(axes.ravel(), specs):
        curves = [
            load_curve(data_root, var, "tight", "Au+Au data, tight WP80", pt_bins, cent_bins, 2),
            load_curve(signal_root, var, "tight_sig", "Embedded direct photons, tight", pt_bins, cent_bins, 2),
            load_curve(inclusive_root, var, "tight_bkg", "Embedded inclusive jets, tight", pt_bins, cent_bins, 2),
            load_curve(data_root, var, "npbPass", "Au+Au data, NPB-tagged", pt_bins, cent_bins, 2),
        ]
        for curve, color, style in zip(curves, ["#111827", "#C0262D", "#2563EB", "#16A34A"], ["-", "-", "-", "--"]):
            draw_step(ax, curve, color, lw=2.2 if style == "-" else 1.9, ls=style, marker="o" if curve.label.startswith("Au+Au data, tight") else None, alpha=0.92)
        decorate(ax, var, title)
        ax.set_xlim(*XRANGES[var])
        ymax = max(float(np.max(c.values)) for c in curves if len(c.values))
        ax.set_ylim(0, ymax * 1.2 if ymax > 0 else 1)
        summary[var] = {c.tag: {"label": c.label, "integral": c.raw_integral, "matches": len(c.matches)} for c in curves}
    handles, labels = axes.ravel()[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.58, 0.945), ncols=2, frameon=False, fontsize=9.7)
    fig.suptitle("THE-42 b009 clean shower-shape readout", x=0.07, y=0.985, ha="left", fontsize=17.5, fontweight="bold")
    fig.text(
        0.07,
        0.935,
        rf"Latest-production v008 Au+Au, {bin_range_label(pt_bins)}, 0-80% centrality; all curves area-normalized.",
        fontsize=11.2,
        color="#374151",
    )
    fig.tight_layout(rect=[0.02, 0.02, 0.98, 0.875])
    out = outdir / "the42_b009_clean_tight_wp80_shower_shapes.png"
    fig.savefig(out)
    plt.close(fig)
    return out, summary


def render_npb_score(data_root: Path, outdir: Path, pt_bins, cent_bins) -> tuple[Path, dict]:
    var = "npbScore"
    curves = [
        load_curve(data_root, var, "inclusive", "Au+Au data, all NPB-scored candidates", pt_bins, cent_bins, 1),
        load_curve(data_root, var, "npbPass", "Au+Au data, NPB-tagged", pt_bins, cent_bins, 1),
    ]
    fig, ax = plt.subplots(figsize=(8.7, 5.6), dpi=180)
    draw_step(ax, curves[0], "#111827", marker="o", lw=2.4)
    draw_step(ax, curves[1], "#16A34A", marker="o", lw=2.8)
    decorate(ax, var, "NPB score distribution in b009 data")
    ax.text(0.035, 0.875, rf"{bin_range_label(pt_bins)}, 0-80% centrality", transform=ax.transAxes, ha="left", va="top", fontsize=10.5)
    ax.axvline(0.5, color="#16A34A", ls="--", lw=1.5, alpha=0.75)
    ax.text(0.515, 0.80, "NPB > 0.5", transform=ax.transAxes, fontsize=10.2, color="#14532D")
    ax.set_xlim(*XRANGES[var])
    ymax = max(float(np.max(c.values)) for c in curves if len(c.values))
    ax.set_ylim(0, ymax * 1.22 if ymax > 0 else 1)
    ax.legend(loc="upper left", bbox_to_anchor=(0.33, 0.99), frameon=False, fontsize=9.7)
    fig.tight_layout()
    out = outdir / "the42_b009_data_npb_score_distribution.png"
    fig.savefig(out)
    plt.close(fig)
    return out, {c.tag: {"label": c.label, "integral": c.raw_integral, "matches": len(c.matches)} for c in curves}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", type=Path, default=DEFAULT_DATA)
    ap.add_argument("--signal-root", type=Path, default=DEFAULT_SIGNAL)
    ap.add_argument("--inclusive-root", type=Path, default=DEFAULT_INCLUSIVE)
    ap.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    ap.add_argument("--pt-bins", default="22:24,24:26,26:28")
    ap.add_argument("--cent-bins", default="0:10,10:20,20:30,30:40,40:50,50:60,60:80")
    args = ap.parse_args()

    setup_style()
    args.outdir.mkdir(parents=True, exist_ok=True)
    pt_bins = parse_bins(args.pt_bins)
    cent_bins = parse_bins(args.cent_bins)
    for path in (args.data_root, args.signal_root, args.inclusive_root):
        if not path.exists():
            raise FileNotFoundError(path)

    npb_plot, npb_summary = render_ppg12_like_npb(args.data_root, args.signal_root, args.inclusive_root, args.outdir, pt_bins, cent_bins)
    grid_plot, grid_summary = render_tight_grid(args.data_root, args.signal_root, args.inclusive_root, args.outdir, pt_bins, cent_bins)
    score_plot, score_summary = render_npb_score(args.data_root, args.outdir, pt_bins, cent_bins)
    manifest = {
        "inputs": {
            "data_root": str(args.data_root),
            "signal_root": str(args.signal_root),
            "inclusive_root": str(args.inclusive_root),
        },
        "pt_bins": pt_bins,
        "cent_bins": cent_bins,
        "plots": {
            "ppg12_style_npb_e11e33": str(npb_plot),
            "tight_wp80_grid": str(grid_plot),
            "npb_score": str(score_plot),
        },
        "summaries": {
            "ppg12_style_npb_e11e33": npb_summary,
            "tight_wp80_grid": grid_summary,
            "npb_score": score_summary,
        },
        "note": "NPB-tagged data uses explicit h_ss_*_npbPass histograms from the b009 merged AuAu data ROOT.",
    }
    manifest_path = args.outdir / "the42_b009_clean_shower_shape_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
