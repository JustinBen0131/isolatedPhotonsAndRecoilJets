#!/usr/bin/env python3
"""Make PPG12-IAN-style BDT score PNGs for THE-42 b009 0-20% AuAu."""

from __future__ import annotations

import json
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


THIS_FILE = Path(__file__).resolve()
REPO = next((p for p in THIS_FILE.parents if (p / "AGENTS.md").exists()), THIS_FILE.parents[5])
SCRIPT_DIR = THIS_FILE.parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.append(str(SCRIPT_DIR))

from make_the42_b009_pp_overlay_three_shape_slide import (  # noqa: E402
    DATA_ROOT,
    INCLUSIVE_ROOT,
    SIGNAL_ROOT,
    OUTDIR,
    build_index,
    open_root,
    setup_style,
)


PT_BINS = [(15, 18), (18, 20), (20, 22), (22, 24), (24, 26), (26, 28), (28, 30), (30, 35)]
CENT_BINS = [(0, 10), (10, 20)]


@dataclass
class Curve:
    label: str
    centers: np.ndarray
    edges: np.ndarray
    values: np.ndarray
    errors: np.ndarray
    integral: float
    matches: list[str]
    missing: list[str]


def target(tag: str, pt: tuple[int, int], cent: tuple[int, int]) -> str:
    return f"h_ss_npbScore_{tag}_pT_{pt[0]}_{pt[1]}_cent_{cent[0]}_{cent[1]}"


def hist_to_curve(hist, label: str, matches: list[str], missing: list[str]) -> Curve:
    nb = hist.GetNbinsX()
    edges = np.array([hist.GetXaxis().GetBinLowEdge(i) for i in range(1, nb + 2)], dtype=float)
    centers = np.array([hist.GetXaxis().GetBinCenter(i) for i in range(1, nb + 1)], dtype=float)
    counts = np.array([hist.GetBinContent(i) for i in range(1, nb + 1)], dtype=float)
    errors = np.array([hist.GetBinError(i) for i in range(1, nb + 1)], dtype=float)
    mask = (centers >= 0.0) & (centers <= 1.0)
    keep = np.flatnonzero(mask)
    if len(keep):
        edges = np.concatenate(([edges[keep[0]]], edges[keep + 1]))
    centers = centers[mask]
    counts = counts[mask]
    errors = errors[mask]
    integral = float(np.sum(counts))
    if integral > 0:
        return Curve(label, centers, edges, counts / integral, errors / integral, integral, matches, missing)
    return Curve(label, centers, edges, counts, errors, integral, matches, missing)


def load_curve(root_path: Path, index: dict[str, list[str]], tag: str, label: str) -> Curve:
    f = open_root(root_path)
    acc = None
    matches: list[str] = []
    missing: list[str] = []
    try:
        for pt in PT_BINS:
            for cent in CENT_BINS:
                name = target(tag, pt, cent)
                paths = index.get(name, [])
                if not paths:
                    missing.append(name)
                    continue
                h0 = f.Get(paths[0])
                if not h0 or not h0.InheritsFrom("TH1"):
                    missing.append(paths[0])
                    continue
                h = h0.Clone(f"{name}_{label}_{tag}")
                h.SetDirectory(0)
                matches.append(paths[0])
                if acc is None:
                    acc = h.Clone(f"sum_{label}_{tag}")
                    acc.SetDirectory(0)
                else:
                    acc.Add(h)
    finally:
        f.Close()
    if acc is None:
        raise RuntimeError(f"no histograms for {root_path} tag={tag}")
    return hist_to_curve(acc, label, matches, missing)


def draw_sphenix(ax) -> None:
    ax.text(0.04, 0.955, "sPHENIX", transform=ax.transAxes, ha="left", va="top", fontsize=13.0, fontstyle="italic", fontweight="bold")
    ax.text(0.205, 0.955, "Internal", transform=ax.transAxes, ha="left", va="top", fontsize=13.0)


def draw_curve(ax, c: Curve, color: str, *, marker: bool = False) -> None:
    if marker:
        ax.errorbar(c.centers, c.values, yerr=c.errors, fmt="o", color=color, mfc=color, mec=color, ms=4.2, lw=1.0, capsize=1.8)
    else:
        ax.stairs(c.values, c.edges, color=color, lw=1.7)


def chi2_ndf(data: Curve, ref: Curve) -> tuple[float, int]:
    n = min(len(data.values), len(ref.values))
    chi2 = 0.0
    ndof = 0
    for i in range(n):
        err = data.errors[i] if data.errors[i] > 0 else np.sqrt(max(data.values[i], 0.0)) if data.values[i] > 0 else 0.0
        if err <= 0:
            continue
        chi2 += float(((data.values[i] - ref.values[i]) / err) ** 2)
        ndof += 1
    return chi2, max(0, ndof - 1)


def plot_ppg12_style(
    out: Path,
    title_cut: str,
    data: Curve,
    signal: Curve,
    inclusive: Curve,
    npb_template: Curve | None = None,
) -> dict:
    fig = plt.figure(figsize=(5.4, 7.0), dpi=180)
    gs = fig.add_gridspec(2, 1, height_ratios=[3.4, 1.15], hspace=0.035)
    ax = fig.add_subplot(gs[0])
    rax = fig.add_subplot(gs[1], sharex=ax)
    draw_curve(ax, data, "black", marker=True)
    draw_curve(ax, signal, "#F0442E")
    draw_curve(ax, inclusive, "#4367FF")
    if npb_template is not None:
        draw_curve(ax, npb_template, "#2CA02C")
    draw_sphenix(ax)
    chi2, ndof = chi2_ndf(data, inclusive)
    lines = [
        r"AuAu $\sqrt{s_{NN}}=200$ GeV",
        "0-20%, 15 < $p_T^\\gamma$ < 35 GeV",
        title_cut,
        rf"$\chi^2$/ndf = {chi2:.1f}/{ndof:d} = {chi2 / ndof:.2f}" if ndof > 0 else r"$\chi^2$/ndf unavailable",
    ]
    ax.text(0.04, 0.875, "\n".join(lines), transform=ax.transAxes, ha="left", va="top", fontsize=10.6)
    labels = ["Data", "Signal MC", "Inclusive MC"]
    handles = [
        plt.Line2D([], [], color="black", marker="o", lw=1.0, ms=4.5, label="Data"),
        plt.Line2D([], [], color="#F0442E", lw=1.7, label="Signal MC"),
        plt.Line2D([], [], color="#4367FF", lw=1.7, label="Inclusive MC"),
    ]
    if npb_template is not None:
        labels.append("NPB-tagged data")
        handles.append(plt.Line2D([], [], color="#2CA02C", lw=1.7, label="NPB-tagged data"))
    ax.legend(handles=handles, loc="upper right", frameon=False, fontsize=10.5, handlelength=1.8)
    ax.set_ylabel("normalized counts", fontsize=13.0)
    ax.set_xlim(0.0, 1.0)
    ymax = max(np.max(data.values), np.max(signal.values), np.max(inclusive.values), np.max(npb_template.values) if npb_template is not None else 0.0)
    ax.set_ylim(0.0, ymax * 1.35 if ymax > 0 else 1.0)
    ax.tick_params(direction="in", which="both", top=True, right=True, labelbottom=False)
    ax.minorticks_on()

    residual = data.values - inclusive.values
    rax.axhline(0.0, color="#555555", lw=0.8, ls=":")
    rax.errorbar(data.centers, residual, yerr=data.errors, fmt="o", color="black", mfc="black", ms=3.2, lw=0.9, capsize=1.5)
    rax.set_ylabel("Data - Incl. MC", fontsize=10.5)
    rax.set_xlabel("bdt", fontsize=12.5, loc="right")
    max_abs = max(0.03, float(np.nanmax(np.abs(residual))) * 1.25 if residual.size else 0.03)
    rax.set_ylim(-max_abs, max_abs)
    rax.tick_params(direction="in", which="both", top=True, right=True)
    rax.minorticks_on()
    fig.savefig(out)
    plt.close(fig)
    return {
        "output": str(out),
        "data_integral": data.integral,
        "signal_integral": signal.integral,
        "inclusive_integral": inclusive.integral,
        "npb_template_integral": npb_template.integral if npb_template else None,
        "chi2": chi2,
        "ndof": ndof,
        "data_matches": len(data.matches),
        "signal_matches": len(signal.matches),
        "inclusive_matches": len(inclusive.matches),
        "npb_template_matches": len(npb_template.matches) if npb_template else None,
        "missing": {
            "data": data.missing,
            "signal": signal.missing,
            "inclusive": inclusive.missing,
            "npb_template": npb_template.missing if npb_template else [],
        },
    }


def plot_unavailable(out: Path) -> dict:
    fig, ax = plt.subplots(figsize=(5.4, 7.0), dpi=180)
    ax.axis("off")
    draw_sphenix(ax)
    ax.text(0.04, 0.86, "AuAu $\\sqrt{s_{NN}}=200$ GeV\n0-20%, 15 < $p_T^\\gamma$ < 35 GeV", transform=ax.transAxes, fontsize=12.0, va="top")
    ax.text(
        0.04,
        0.68,
        "Fig. 20-style tight-BDT score\ncannot be drawn from b009",
        transform=ax.transAxes,
        fontsize=17.0,
        fontweight="bold",
        va="top",
    )
    ax.text(
        0.04,
        0.50,
        "\n".join(
            [
                "Evidence:",
                "• No TTree is stored in data/signal/inclusive merged ROOTs.",
                "• No h_tightBDTScore_preselected histograms are present.",
                "• h_tightFail_bdt is a 1-bin count histogram, not a score distribution.",
                "",
                "Needed for the true Fig. 20 analog:",
                "fill and merge final photon-ID tight_bdt_score histograms",
                "before and after the tight BDT selection.",
            ]
        ),
        transform=ax.transAxes,
        fontsize=12.0,
        va="top",
    )
    fig.savefig(out)
    plt.close(fig)
    return {"output": str(out), "status": "unavailable", "reason": "final tight-BDT score distribution is not stored in b009 merged outputs"}


def main() -> int:
    setup_style()
    OUTDIR.mkdir(parents=True, exist_ok=True)
    indexes = {
        "data": build_index(DATA_ROOT),
        "signal": build_index(SIGNAL_ROOT),
        "inclusive": build_index(INCLUSIVE_ROOT),
    }
    rows = []
    no_cut = {
        "data": load_curve(DATA_ROOT, indexes["data"], "inclusive", "Data"),
        "signal": load_curve(SIGNAL_ROOT, indexes["signal"], "inclusive", "Signal MC"),
        "inclusive": load_curve(INCLUSIVE_ROOT, indexes["inclusive"], "inclusive", "Inclusive MC"),
        "npb": load_curve(DATA_ROOT, indexes["data"], "npbFail", "NPB-tagged data"),
    }
    rows.append(
        plot_ppg12_style(
            OUTDIR / "the42_b009_0_20_ppg12_fig13_style_bdt_no_npb_cut.png",
            "w/o NPB cut",
            no_cut["data"],
            no_cut["signal"],
            no_cut["inclusive"],
            no_cut["npb"],
        )
    )
    npb_cut = {
        "data": load_curve(DATA_ROOT, indexes["data"], "npbPass", "Data"),
        "signal": load_curve(SIGNAL_ROOT, indexes["signal"], "npbPass", "Signal MC"),
        "inclusive": load_curve(INCLUSIVE_ROOT, indexes["inclusive"], "npbPass", "Inclusive MC"),
    }
    rows.append(
        plot_ppg12_style(
            OUTDIR / "the42_b009_0_20_ppg12_fig19_style_bdt_with_npb_cut.png",
            "w/ NPB cut",
            npb_cut["data"],
            npb_cut["signal"],
            npb_cut["inclusive"],
            None,
        )
    )
    rows.append(plot_unavailable(OUTDIR / "the42_b009_0_20_ppg12_fig20_style_tight_bdt_score_unavailable.png"))
    manifest = {
        "schema": "THE42_B009_PPG12_BDT_SCORE_PNGS_V1",
        "ppg12_reference": "usefulDocs/PPG12_analysis_note_2026-05-21_v4_current_IAN.pdf figures 13, 19, 20",
        "data_root": str(DATA_ROOT),
        "signal_root": str(SIGNAL_ROOT),
        "inclusive_root": str(INCLUSIVE_ROOT),
        "centrality": "0-20%",
        "pt_bins": PT_BINS,
        "outputs": rows,
    }
    manifest_path = OUTDIR / "the42_b009_0_20_ppg12_bdt_score_pngs_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    for row in rows:
        print(row["output"])
    print(manifest_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
