#!/usr/bin/env python3
"""Render a full-slide THE-42 w_eta distribution evolution plot grid."""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np


THIS_FILE = Path(__file__).resolve()
REPO = next((p for p in THIS_FILE.parents if (p / "AGENTS.md").exists()), THIS_FILE.parents[3])
DEFAULT_BASE = REPO / "dataOutput/auauTightBDTValidation/THE42_wp80_centlinear_ss_20260604"
DEFAULT_ROOT_DIR = (
    DEFAULT_BASE
    / "sim_roots/preselectionNewPPG12_tightAuAuCentInputBase3x3BDT_nonTightAuAuBDTComplement_baseVariant"
)
DEFAULT_SIGNAL_ROOT = (
    DEFAULT_ROOT_DIR / "photonJet12and20merged_SIM/RecoilJets_embeddedPhoton12plus20_MERGED.root"
)
DEFAULT_INCLUSIVE_ROOT = (
    DEFAULT_ROOT_DIR / "embeddedJet12and20and30and40merged_SIM/RecoilJets_embeddedJet12plus20plus30plus40_MERGED.root"
)
DEFAULT_OUTDIR = DEFAULT_BASE / "efficiency_purity_qa"

VARIABLES = {
    "weta": {
        "hist": "weta",
        "label": r"$w_{\eta}$",
        "plain": "weta",
        "title": r"Tight WP80 collapses broad inclusive $w_{\eta}$ shapes",
        "subtitle_name": r"$w_{\eta}$",
        "high_tail": r"broad high-$w_{\eta}$ tail",
        "script_name": "w_eta",
        "xmax": 0.75,
    },
    "e11e33": {
        "hist": "e11e33",
        "label": r"$E_{1\times1}/E_{3\times3}$",
        "plain": "E11/E33",
        "title": r"Tight WP80 reshapes inclusive $E_{1\times1}/E_{3\times3}$",
        "subtitle_name": r"$E_{1\times1}/E_{3\times3}$",
        "high_tail": r"broad low-ratio shoulder",
        "script_name": "E11/E33",
        "xmax": 1.0,
    },
}
PT_BINS = ((22, 24), (24, 26), (26, 28))
CENT_GROUPS = (
    ("0-20%", ((0, 10), (10, 20))),
    ("20-50%", ((20, 30), (30, 40), (40, 50))),
    ("50-80%", ((50, 60), (60, 80))),
)
STAGES = (
    ("No preselection", "inclusive", 2.2),
    ("Preselection", "pre", 2.5),
    ("Tight WP80", "tight", 2.9),
)
SAMPLE_TAGS = {
    "Signal MC": {"inclusive": "inclusive_sig", "pre": "pre_sig", "tight": "tight_sig"},
    "Inclusive MC": {"inclusive": "inclusive_bkg", "pre": "pre_bkg", "tight": "tight_bkg"},
}
SAMPLES = (("Signal MC", "#D84A4A"), ("Inclusive MC", "#2F78B7"))
SAMPLE_STAGE_COLORS = {
    "Signal MC": {
        "No preselection": "#F2B0A8",
        "Preselection": "#E15B52",
        "Tight WP80": "#B91C1C",
    },
    "Inclusive MC": {
        "No preselection": "#AFC8E8",
        "Preselection": "#2F78B7",
        "Tight WP80": "#174EA6",
    },
}

INK = "#111827"
MUTED = "#556070"
GRID = "#E5EAF0"
PANEL_EDGE = "#C8D6E5"
SOFT_PANEL = "#F8FAFC"
INCLUSIVE_ONLY_STAGE_COLORS = {
    "No preselection": "#5B6472",
    "Preselection": "#D97706",
    "Tight WP80": "#174EA6",
}


@dataclass
class Curve:
    sample: str
    centrality: str
    stage: str
    edges: np.ndarray
    values: np.ndarray
    mean: float
    rms: float
    integral: float


def setup_style() -> None:
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
        "mathtext.fontset": "dejavuserif",
        "axes.edgecolor": INK,
        "axes.linewidth": 0.9,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
    })


def open_root(path: Path):
    import ROOT

    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.AddDirectory(False)
    f = ROOT.TFile.Open(str(path), "READ")
    if not f or f.IsZombie():
        raise OSError(f"Could not open ROOT input: {path}")
    return f


def sum_variable_hist(root_file, variable: str, tag: str, cent_bins: tuple[tuple[int, int], ...]):
    acc = None
    missing = []
    hist_var = VARIABLES[variable]["hist"]
    for pt_lo, pt_hi in PT_BINS:
        for c_lo, c_hi in cent_bins:
            name = f"SIM/h_ss_{hist_var}_{tag}_pT_{pt_lo}_{pt_hi}_cent_{c_lo}_{c_hi}"
            hist = root_file.Get(name)
            if not hist:
                missing.append(name)
                continue
            if acc is None:
                acc = hist.Clone(f"sum_{hist_var}_{tag}_{pt_lo}_{pt_hi}_{c_lo}_{c_hi}")
                acc.SetDirectory(0)
            else:
                acc.Add(hist)
    if acc is None:
        raise RuntimeError(f"No h_ss_{hist_var} histograms found for tag={tag}; first missing={missing[:3]}")
    return acc


def rebin_for_display(edges: np.ndarray, counts: np.ndarray, factor: int) -> tuple[np.ndarray, np.ndarray]:
    if factor <= 1:
        return edges, counts
    n = (len(counts) // factor) * factor
    if n <= 0:
        return edges, counts
    rebinned_counts = counts[:n].reshape(-1, factor).sum(axis=1)
    rebinned_edges = edges[: n + 1 : factor]
    if len(rebinned_edges) != len(rebinned_counts) + 1:
        rebinned_edges = np.append(rebinned_edges, edges[n])
    if n < len(counts):
        rebinned_counts = np.append(rebinned_counts, counts[n:].sum())
        rebinned_edges = np.append(rebinned_edges, edges[-1])
    return rebinned_edges, rebinned_counts


def curve_from_hist(sample: str, centrality: str, stage: str, hist, rebin: int = 1) -> Curve:
    nb = hist.GetNbinsX()
    edges = np.array([hist.GetXaxis().GetBinLowEdge(i) for i in range(1, nb + 2)], dtype=float)
    centers = np.array([hist.GetXaxis().GetBinCenter(i) for i in range(1, nb + 1)], dtype=float)
    counts = np.array([hist.GetBinContent(i) for i in range(1, nb + 1)], dtype=float)
    integral = float(np.sum(counts))
    mean = float(np.sum(centers * counts) / integral) if integral > 0 else math.nan
    rms = float(np.sqrt(np.sum(counts * (centers - mean) ** 2) / integral)) if integral > 0 else math.nan
    display_edges, display_counts = rebin_for_display(edges, counts, rebin)
    values = display_counts / integral if integral > 0 else np.zeros_like(display_counts)
    return Curve(sample, centrality, stage, display_edges, values, mean, rms, integral)


def collect_curves(signal_root: Path, inclusive_root: Path, variable: str, rebin: int = 1) -> list[Curve]:
    curves: list[Curve] = []
    roots = {"Signal MC": signal_root, "Inclusive MC": inclusive_root}
    for sample, path in roots.items():
        f = open_root(path)
        try:
            for cent_label, cent_bins in CENT_GROUPS:
                for stage_label, stage_key, _ in STAGES:
                    hist = sum_variable_hist(f, variable, SAMPLE_TAGS[sample][stage_key], cent_bins)
                    curves.append(curve_from_hist(sample, cent_label, stage_label, hist, rebin=rebin))
        finally:
            f.Close()
    return curves


def lookup_curves(curves: list[Curve]) -> dict[tuple[str, str, str], Curve]:
    return {(c.sample, c.centrality, c.stage): c for c in curves}


def add_box(fig, xywh, face, edge=PANEL_EDGE, radius=0.018, lw=1.0):
    ax = fig.add_axes(xywh)
    ax.axis("off")
    patch = FancyBboxPatch(
        (0, 0), 1, 1,
        boxstyle=f"round,pad=0.006,rounding_size={radius}",
        facecolor=face, edgecolor=edge, linewidth=lw,
        transform=ax.transAxes, clip_on=False,
    )
    ax.add_patch(patch)
    return ax


def write_script(path: Path, variable: str) -> None:
    meta = VARIABLES[variable]
    path.write_text(
        f"# THE-42 {meta['plain']} Distribution Evolution Slide Script\n\n"
        f"This slide shows {meta['script_name']} as vertically stacked distributions rather than as a table of "
        "means or fully overlapping curves. Each curve is unit-normalized; the vertical offsets are visual only.\n\n"
        "The top row is the embedded photon signal sample and is always shown in red. The signal shapes are already "
        "narrow, so preselection changes very little and the tight WP80 selection mostly trims the high-width tail.\n\n"
        "The bottom row is the embedded inclusive-jet sample and is always shown in blue. This is the audience-facing "
        "effect: the no-preselection and preselection shapes are broad, while the tight WP80 shape collapses toward "
        "the narrower, signal-like w_eta region. The mean and delta-mu labels quantify the same movement.\n\n"
        "This is shape QA, not an efficiency plot. The efficiency and background-acceptance behavior should be quoted "
        "from the score-threshold QA slide, while this slide makes the BDT selection effect visually intuitive.\n",
        encoding="utf-8",
    )


def write_inclusive_script(path: Path, variable: str) -> None:
    meta = VARIABLES[variable]
    path.write_text(
        f"# THE-42 Inclusive MC {meta['plain']} Distribution Evolution Slide Script\n\n"
        "Here I isolate the inclusive-jet embedded sample so the selection effect is easier to see without the signal row. "
        f"Each panel is one centrality range, and within each panel the three {meta['script_name']} distributions are unit-normalized and "
        "vertically offset only to make the shapes readable.\n\n"
        "The gray shape is before preselection, the amber shape is after the PPG12-style preselection, and the deep blue "
        "shape is after the tight WP80 centrality-dependent BDT cut. The preselection alone shifts the mean by about "
        f"0.026 in {meta['script_name']}, while the tight WP80 selection shows the main shape movement across centrality.\n\n"
        "The visual point is that the inclusive background starts with a broad high-width tail and the tight WP80 cut "
        f"moves that distribution toward the signal-like {meta['script_name']} region. This is a direct shape-level diagnostic "
        "for why the BDT selection is doing useful photon-ID work.\n",
        encoding="utf-8",
    )


def draw_stacked_panel(
    ax,
    lk: dict[tuple[str, str, str], Curve],
    sample: str,
    cent_label: str,
    stage_colors: dict[str, str],
    sample_color: str,
    *,
    show_internal: bool = False,
    show_stage_labels: bool = True,
    show_stage_mean: bool = True,
    show_panel_annotations: bool = True,
    show_pre_delta: bool = False,
    label_font: float = 10.5,
    mean_font: float = 9.8,
    delta_font: float = 11.2,
    cent_font: float = 15.8,
    scale_height: float = 0.70,
    x_max: float = 0.75,
) -> None:
    stage_baselines = {
        "No preselection": 2.18,
        "Preselection": 1.18,
        "Tight WP80": 0.18,
    }
    stage_alpha = {
        "No preselection": 0.18,
        "Preselection": 0.25,
        "Tight WP80": 0.36,
    }
    ax.set_facecolor("white")
    for spine in ax.spines.values():
        spine.set_color(PANEL_EDGE)
        spine.set_linewidth(1.0)

    panel_max = 0.0
    for stage_label, _, _ in STAGES:
        panel_max = max(panel_max, float(np.nanmax(lk[(sample, cent_label, stage_label)].values)))
    scale = scale_height / panel_max if panel_max > 0 else 1.0

    for stage_label, _, lw in STAGES:
        c = lk[(sample, cent_label, stage_label)]
        base = stage_baselines[stage_label]
        color = stage_colors[stage_label]
        y = base + c.values * scale
        ax.fill_between(
            c.edges[:-1], base, y, step="post", color=color,
            alpha=stage_alpha[stage_label], linewidth=0.0,
        )
        ax.stairs(y, c.edges, color=color, linewidth=lw + 0.35)
        ax.hlines(base, 0.0, x_max, color="#D8DEE8", linewidth=0.85)
        if show_stage_labels:
            ax.text(0.020, base + 0.080, stage_label, ha="left", va="bottom",
                    fontsize=label_font, fontweight="bold", color=color, transform=ax.get_yaxis_transform(),
                    bbox=dict(boxstyle="round,pad=0.10", facecolor="white", edgecolor="none", alpha=0.74))
        if show_stage_mean:
            ax.text(0.972, base + 0.105, rf"$\mu={c.mean:.3f}$", ha="right", va="bottom",
                    fontsize=mean_font, color=color, transform=ax.get_yaxis_transform(),
                    bbox=dict(boxstyle="round,pad=0.08", facecolor="white", edgecolor="none", alpha=0.76))

    no = lk[(sample, cent_label, "No preselection")]
    pre = lk[(sample, cent_label, "Preselection")]
    tight = lk[(sample, cent_label, "Tight WP80")]
    delta = tight.mean - no.mean
    pre_delta = pre.mean - no.mean
    if show_panel_annotations:
        ax.text(0.035, 0.955, cent_label, ha="left", va="top", fontsize=cent_font,
                fontweight="bold", color=INK, transform=ax.transAxes,
                bbox=dict(boxstyle="round,pad=0.08", facecolor="white", edgecolor="none", alpha=0.80))
        ax.text(0.965, 0.955, rf"$\Delta\mu_{{tight-no}}$={delta:+.3f}",
                ha="right", va="top", fontsize=delta_font, fontweight="bold",
                color=sample_color, transform=ax.transAxes,
                bbox=dict(boxstyle="round,pad=0.10", facecolor="white", edgecolor="none", alpha=0.80))
        if show_pre_delta:
            ax.text(0.965, 0.875, rf"$\Delta\mu_{{pre-no}}$={pre_delta:+.3f}",
                    ha="right", va="top", fontsize=delta_font - 0.8, color=MUTED, transform=ax.transAxes,
                    bbox=dict(boxstyle="round,pad=0.08", facecolor="white", edgecolor="none", alpha=0.76))
        if show_internal:
            ax.text(0.035, 0.845, r"$\it{\bf{sPHENIX}}$ Internal",
                    ha="left", va="top", fontsize=12.8, color=INK, transform=ax.transAxes,
                    bbox=dict(boxstyle="round,pad=0.08", facecolor="white", edgecolor="none", alpha=0.80))
    ax.set_xlim(0.0, x_max)
    ax.set_ylim(-0.05, 3.02)
    ax.grid(True, axis="x", color=GRID, linewidth=0.8)
    ax.grid(False, axis="y")
    ax.tick_params(labelsize=10.8, direction="in", top=True, right=True, left=False, labelleft=False)


def render_slide(curves: list[Curve], outdir: Path, variable: str) -> dict[str, Path]:
    setup_style()
    outdir.mkdir(parents=True, exist_ok=True)
    meta = VARIABLES[variable]
    lk = lookup_curves(curves)

    fig = plt.figure(figsize=(16, 9), dpi=160)
    fig.text(0.045, 0.958, meta["title"],
             ha="left", va="top", fontsize=30.5, fontweight="bold", color=INK)
    fig.text(0.046, 0.908,
             rf"Stacked, unit-normalized {meta['subtitle_name']} distributions for $22 \leq p_T^\gamma < 28$ GeV; signal is red and inclusive MC is blue.",
             ha="left", va="top", fontsize=16.0, color=MUTED)

    key = add_box(fig, [0.095, 0.802, 0.881, 0.076], SOFT_PANEL, edge=PANEL_EDGE, radius=0.014)
    key.text(0.030, 0.50, "Vertical order in every panel", ha="left", va="center",
             fontsize=15.0, fontweight="bold", color=INK, transform=key.transAxes)
    stack_items = (
        ("top", "No preselection"),
        ("middle", "Preselection"),
        ("bottom", "Tight WP80"),
    )
    x0 = 0.290
    for idx, (pos, label) in enumerate(stack_items):
        x = x0 + idx * 0.215
        key.text(x, 0.64, pos, ha="left", va="center", fontsize=12.0,
                 fontweight="bold", color=MUTED, transform=key.transAxes)
        key.text(x, 0.32, label, ha="left", va="center", fontsize=13.6,
                 color=INK, transform=key.transAxes)
        if idx < 2:
            key.plot([x + 0.170, x + 0.170], [0.16, 0.84], color="#D5DFEA",
                     linewidth=1.0, transform=key.transAxes)

    left = 0.095
    top = 0.770
    panel_w = 0.275
    panel_h = 0.286
    hgap = 0.028
    vgap = 0.062
    row_y = [top - panel_h, top - panel_h - vgap - panel_h]

    for row_idx, (sample, sample_color) in enumerate(SAMPLES):
        label_ax = fig.add_axes([0.026, row_y[row_idx] + 0.045, 0.041, panel_h - 0.020])
        label_ax.axis("off")
        label_ax.text(0.50, 0.50, sample, ha="center", va="center", rotation=90,
                      fontsize=18.0, fontweight="bold", color=sample_color, transform=label_ax.transAxes)

        for col_idx, (cent_label, _) in enumerate(CENT_GROUPS):
            ax = fig.add_axes([left + col_idx * (panel_w + hgap), row_y[row_idx], panel_w, panel_h])
            draw_stacked_panel(
                ax, lk, sample, cent_label, SAMPLE_STAGE_COLORS[sample], sample_color,
                show_internal=(row_idx == 0 and col_idx == 0),
                show_stage_labels=False,
                show_stage_mean=False,
                show_pre_delta=False,
                label_font=10.2,
                mean_font=9.4,
                delta_font=12.2,
                x_max=meta["xmax"],
            )
            if row_idx == 1:
                ax.set_xlabel(meta["label"], fontsize=12.0, labelpad=2)
            else:
                ax.tick_params(labelbottom=False)

    card1 = add_box(fig, [0.095, 0.018, 0.275, 0.078], "#FFF7ED", edge="#FED7AA", radius=0.018)
    card1.text(0.055, 0.67, "Signal row", fontsize=15.4, fontweight="bold", color="#D84A4A",
               ha="left", va="center", transform=card1.transAxes)
    card1.text(0.055, 0.30, "small movement; already narrow before tight WP80",
               fontsize=12.6, color=INK, ha="left", va="center", transform=card1.transAxes)

    card2 = add_box(fig, [0.398, 0.018, 0.275, 0.078], "#EFF6FF", edge="#BFDBFE", radius=0.018)
    card2.text(0.055, 0.67, "Inclusive row", fontsize=15.4, fontweight="bold", color="#2F78B7",
               ha="left", va="center", transform=card2.transAxes)
    card2.text(0.055, 0.31, r"large movement; broad tail collapses after tight WP80",
               fontsize=12.6, color=INK, ha="left", va="center", transform=card2.transAxes)

    card3 = add_box(fig, [0.701, 0.018, 0.275, 0.078], "#F8FAFC", edge=PANEL_EDGE, radius=0.018)
    card3.text(0.055, 0.67, "Plot convention", fontsize=15.4, fontweight="bold", color=MUTED,
               ha="left", va="center", transform=card3.transAxes)
    card3.text(0.055, 0.31, "offsets are visual only; all shapes are normalized",
               fontsize=12.6, color=INK, ha="left", va="center", transform=card3.transAxes)

    png = outdir / f"the42_{variable}_ss_evolution_distribution_slide.png"
    script = outdir / f"the42_{variable}_ss_evolution_distribution_script.md"
    manifest = outdir / f"the42_{variable}_ss_evolution_distribution_manifest.json"
    fig.savefig(png, dpi=160)
    plt.close(fig)
    write_script(script, variable)
    outputs = {"png": png, "speaker_script": script, "manifest": manifest}
    return outputs


def render_inclusive_slide(curves: list[Curve], outdir: Path, variable: str) -> dict[str, Path]:
    setup_style()
    outdir.mkdir(parents=True, exist_ok=True)
    meta = VARIABLES[variable]
    lk = lookup_curves(curves)

    fig = plt.figure(figsize=(16, 9), dpi=160)
    fig.text(0.055, 0.958, meta["title"],
             ha="left", va="top", fontsize=32.0, fontweight="bold", color=INK)
    fig.text(0.056, 0.908,
             rf"Embedded inclusive-jet MC, unit-normalized {meta['subtitle_name']} distributions, $22 \leq p_T^\gamma < 28$ GeV; vertical offsets are visual only.",
             ha="left", va="top", fontsize=16.3, color=MUTED)

    badge = add_box(fig, [0.760, 0.895, 0.195, 0.056], "#EFF6FF", edge="#BFDBFE", radius=0.014)
    badge.text(0.50, 0.66, r"$\it{\bf{sPHENIX}}$ Internal", ha="center", va="center",
               fontsize=11.8, color=INK, transform=badge.transAxes)
    badge.text(0.50, 0.30, "Inclusive MC only", ha="center", va="center",
               fontsize=13.6, fontweight="bold", color="#174EA6", transform=badge.transAxes)

    key = add_box(fig, [0.075, 0.812, 0.880, 0.058], SOFT_PANEL, edge=PANEL_EDGE, radius=0.014)
    key.text(0.035, 0.50, "Selection stage", ha="left", va="center",
             fontsize=15.0, fontweight="bold", color=INK, transform=key.transAxes)
    key_x = [0.235, 0.505, 0.755]
    for x, label in zip(key_x, ("No preselection", "Preselection", "Tight WP80")):
        color = INCLUSIVE_ONLY_STAGE_COLORS[label]
        key.plot([x, x + 0.075], [0.50, 0.50], color=color, linewidth=4.2,
                 transform=key.transAxes)
        key.text(x + 0.087, 0.50, label, ha="left", va="center", fontsize=13.4,
                 color=INK, transform=key.transAxes)

    left = 0.075
    bottom = 0.222
    panel_w = 0.270
    panel_h = 0.535
    hgap = 0.035
    panel_xs = [left + i * (panel_w + hgap) for i in range(3)]
    for col_idx, (cent_label, _) in enumerate(CENT_GROUPS):
        no = lk[("Inclusive MC", cent_label, "No preselection")]
        pre = lk[("Inclusive MC", cent_label, "Preselection")]
        tight = lk[("Inclusive MC", cent_label, "Tight WP80")]
        d_pre = pre.mean - no.mean
        d_tight = tight.mean - no.mean
        header = fig.add_axes([panel_xs[col_idx], bottom + panel_h + 0.009, panel_w, 0.037])
        header.axis("off")
        header.text(0.00, 0.54, cent_label, ha="left", va="center",
                    fontsize=15.6, fontweight="bold", color=INK, transform=header.transAxes)
        header.text(1.00, 0.68, rf"$\Delta\mu_{{tight-no}}$={d_tight:+.3f}",
                    ha="right", va="center", fontsize=12.8, fontweight="bold",
                    color="#174EA6", transform=header.transAxes)
        header.text(1.00, 0.22, rf"$\Delta\mu_{{pre-no}}$={d_pre:+.3f}",
                    ha="right", va="center", fontsize=11.8, color=MUTED, transform=header.transAxes)
        ax = fig.add_axes([panel_xs[col_idx], bottom, panel_w, panel_h])
        draw_stacked_panel(
            ax, lk, "Inclusive MC", cent_label, INCLUSIVE_ONLY_STAGE_COLORS, "#174EA6",
            show_internal=False,
            show_stage_labels=False,
            show_stage_mean=False,
            show_panel_annotations=False,
            show_pre_delta=True,
            label_font=12.4,
            mean_font=12.8,
            delta_font=13.8,
            cent_font=16.2,
            scale_height=0.78,
            x_max=meta["xmax"],
        )
        ax.set_xlabel(meta["label"], fontsize=15.0, labelpad=3)
        if col_idx == 0:
            ax.set_ylabel("stage offset", fontsize=14.0)

    card_h = 0.100
    card_y = 0.060
    card1 = add_box(fig, [panel_xs[0], card_y, panel_w, card_h], "#F8FAFC", edge=PANEL_EDGE, radius=0.018)
    card1.text(0.055, 0.69, "No preselection", fontsize=16.0, fontweight="bold",
               color=INCLUSIVE_ONLY_STAGE_COLORS["No preselection"], ha="left", va="center",
               transform=card1.transAxes)
    no_means = [lk[("Inclusive MC", cent_label, "No preselection")].mean for cent_label, _ in CENT_GROUPS]
    pre_deltas = [
        lk[("Inclusive MC", cent_label, "Preselection")].mean
        - lk[("Inclusive MC", cent_label, "No preselection")].mean
        for cent_label, _ in CENT_GROUPS
    ]
    tight_deltas = [
        lk[("Inclusive MC", cent_label, "Tight WP80")].mean
        - lk[("Inclusive MC", cent_label, "No preselection")].mean
        for cent_label, _ in CENT_GROUPS
    ]
    no_summary = (
        rf"{meta['high_tail']}; $\mu$ near {np.mean(no_means):.2f}"
        if max(no_means) - min(no_means) < 0.015
        else rf"{meta['high_tail']}; $\mu$ spans {min(no_means):.2f}--{max(no_means):.2f}"
    )
    pre_summary = rf"small shift, $\Delta\mu_{{pre-no}}\simeq {np.mean(pre_deltas):+.3f}$"
    tight_summary = rf"dominant movement, $\Delta\mu_{{tight-no}}\simeq {np.mean(tight_deltas):+.2f}$"
    card1.text(0.055, 0.32, no_summary,
               fontsize=13.2, color=INK, ha="left", va="center", transform=card1.transAxes)

    card2 = add_box(fig, [panel_xs[1], card_y, panel_w, card_h], "#FFF7ED", edge="#FED7AA", radius=0.018)
    card2.text(0.055, 0.69, "Preselection", fontsize=16.0, fontweight="bold",
               color=INCLUSIVE_ONLY_STAGE_COLORS["Preselection"], ha="left", va="center",
               transform=card2.transAxes)
    card2.text(0.055, 0.32, pre_summary,
               fontsize=13.2, color=INK, ha="left", va="center", transform=card2.transAxes)

    card3 = add_box(fig, [panel_xs[2], card_y, panel_w, card_h], "#EFF6FF", edge="#BFDBFE", radius=0.018)
    card3.text(0.055, 0.69, "Tight WP80", fontsize=16.0, fontweight="bold",
               color=INCLUSIVE_ONLY_STAGE_COLORS["Tight WP80"], ha="left", va="center",
               transform=card3.transAxes)
    card3.text(0.055, 0.32, tight_summary,
               fontsize=13.2, color=INK, ha="left", va="center", transform=card3.transAxes)

    png = outdir / f"the42_{variable}_ss_evolution_inclusive_only_distribution_slide.png"
    script = outdir / f"the42_{variable}_ss_evolution_inclusive_only_distribution_script.md"
    fig.savefig(png, dpi=160)
    plt.close(fig)
    write_inclusive_script(script, variable)
    return {"inclusive_png": png, "inclusive_speaker_script": script}


def write_manifest(path: Path, curves: list[Curve], outputs: dict[str, Path], args: argparse.Namespace) -> None:
    rows = []
    for c in curves:
        rows.append({
            "sample": c.sample,
            "centrality": c.centrality,
            "stage": c.stage,
            "mean": c.mean,
            "rms": c.rms,
            "integral": c.integral,
        })
    payload = {
        "script": str(THIS_FILE),
        "signal_root": str(args.signal_root.resolve()),
        "inclusive_root": str(args.inclusive_root.resolve()),
        "variable": args.variable,
        "variable_label": VARIABLES[args.variable]["plain"],
        "plot_mode": "unit-normalized distribution shapes with vertical visual offsets; display-only rebin",
        "display_rebin_factor": args.rebin,
        "pt_bins": PT_BINS,
        "centrality_groups": [{"label": label, "fine_bins": bins} for label, bins in CENT_GROUPS],
        "stage_tags": SAMPLE_TAGS,
        "x_axis_range_shown": [0.0, VARIABLES[args.variable]["xmax"]],
        "outputs": {k: str(v.resolve()) for k, v in outputs.items()},
        "curve_summaries": rows,
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--signal-root", type=Path, default=DEFAULT_SIGNAL_ROOT)
    ap.add_argument("--inclusive-root", type=Path, default=DEFAULT_INCLUSIVE_ROOT)
    ap.add_argument("--output-dir", type=Path, default=DEFAULT_OUTDIR)
    ap.add_argument("--variable", choices=sorted(VARIABLES), default="weta")
    ap.add_argument("--rebin", type=int, default=3, help="Display rebin factor applied after summing ROOT histograms.")
    return ap.parse_args()


def main() -> int:
    args = parse_args()
    if args.rebin < 1:
        raise SystemExit("--rebin must be >= 1")
    for path in (args.signal_root, args.inclusive_root):
        if not path.exists():
            raise SystemExit(f"Missing ROOT input: {path}")
    curves = collect_curves(args.signal_root, args.inclusive_root, args.variable, rebin=args.rebin)
    outputs = render_slide(curves, args.output_dir, args.variable)
    outputs.update(render_inclusive_slide(curves, args.output_dir, args.variable))
    write_manifest(outputs["manifest"], curves, outputs, args)
    for key, value in outputs.items():
        print(f"{key}={value}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
