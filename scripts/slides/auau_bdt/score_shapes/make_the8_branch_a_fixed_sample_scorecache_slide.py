#!/usr/bin/env python3
"""Build a full-slide PNG comparing THE-8 Branch A fixed-sample validations."""

from __future__ import annotations

# Keep purpose-folder helpers runnable when invoked directly.
import sys as _codex_sys
from pathlib import Path as _CodexPath

_CODEX_THIS_FILE = _CodexPath(__file__).resolve()
_CODEX_SCRIPTS_DIR = next((p for p in _CODEX_THIS_FILE.parents if p.name == "scripts"), _CODEX_THIS_FILE.parent)
_CODEX_SCRIPTS_DIR_STR = str(_CODEX_SCRIPTS_DIR)
if _CODEX_SCRIPTS_DIR_STR not in _codex_sys.path:
    _codex_sys.path.append(_CODEX_SCRIPTS_DIR_STR)
del _CODEX_THIS_FILE, _CODEX_SCRIPTS_DIR, _CODEX_SCRIPTS_DIR_STR

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch, Rectangle

from slides.common.slide_defaults import SLIDE_DPI, slide_figsize


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
DEFAULT_JSON = (
    REPO
    / "dataOutput/auauTightBDTValidation/THE8_branchA_ladder_scorecache_fullstat_20260527"
    / "fixed_sample_controls/the8_branch_a_ladder_fixed_sample_scorecache_controls_truthSignal_inclusiveJet_v1.json"
)
DEFAULT_OUTDIR = DEFAULT_JSON.parent
DEFAULT_OUTPNG = DEFAULT_OUTDIR / "the8_branchA_fixed_sample_scorecache_controls_truthSignal_inclusiveJet_slide_v1.png"
DEFAULT_PP_BASELINE_JSON = (
    REPO
    / "dataOutput/ppPhotonMLPipeline/ppg12_basev3E_currentIAN_rawOverlayEnvFix_20260527_1630"
    / "validation/fullsim_shuhang_overlay_raw_inclusive_pt1535"
    / "pp_currentian_basev3e_bdt_score_overlay_vs_ppg12_pp_noCent_bdt_15_35_noNPB_eta0-pt3-cut0_summary.json"
)
PRODUCT = "globalEtCent1535_bdt_noIso"

NAVY = "#101828"
INK = "#1D2939"
MUTED = "#475467"
BORDER = "#C9D2DE"
SOFT_GRAY = "#F6F7F9"
SOFT_LAVENDER = "#F3EEF8"
SIGNAL = "#C5392F"
INCLUSIVE = "#1F77B4"
GRID = "#D9DEE7"
SAMPLE_ACCENTS = (
    {"fill": "#ECF7EE", "edge": "#6EA77B", "strip": "#3F8F5A"},
    {"fill": "#FFF4CC", "edge": "#D6A84A", "strip": "#B88400"},
    {"fill": "#F3ECFA", "edge": "#9A7CC3", "strip": "#7651A6"},
)


def add_box(
    fig: plt.Figure,
    xy: tuple[float, float],
    wh: tuple[float, float],
    face: str,
    *,
    edge: str = BORDER,
    lw: float = 1.0,
    radius: float = 0.010,
    zorder: int = -2,
) -> None:
    fig.patches.append(
        FancyBboxPatch(
            xy,
            wh[0],
            wh[1],
            boxstyle=f"round,pad=0.010,rounding_size={radius}",
            transform=fig.transFigure,
            facecolor=face,
            edgecolor=edge,
            linewidth=lw,
            zorder=zorder,
        )
    )


def add_text(
    fig: plt.Figure,
    x: float,
    y: float,
    text: str,
    *,
    size: float,
    weight: str = "normal",
    color: str = INK,
    ha: str = "left",
    va: str = "top",
    linespacing: float = 1.22,
) -> None:
    fig.text(
        x,
        y,
        text,
        fontsize=size,
        fontweight=weight,
        color=color,
        ha=ha,
        va=va,
        linespacing=linespacing,
    )


def step_xy(edges: np.ndarray, density: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return np.repeat(edges, 2)[1:-1], np.repeat(density, 2)


def auc_of(branch: dict, centrality_key: str | None = None) -> float:
    auc_key = f"{PRODUCT}_auc_{centrality_key}" if centrality_key else f"{PRODUCT}_auc"
    value = branch["summary"].get(auc_key, branch["summary"].get(f"{PRODUCT}_auc", "nan"))
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def entry_counts(branch: dict, centrality_key: str | None = None) -> tuple[int, int]:
    if centrality_key:
        block = branch["by_centrality"][centrality_key]
        return int(block["signal"]["entries"]), int(block["background"]["entries"])
    summary = branch["summary"]
    return int(summary.get("signal_entries", 0)), int(summary.get("background_entries", 0))


def wp_fake_rate_at_signal_eff(
    branch: dict,
    *,
    centrality_key: str | None = None,
    target_efficiency: float = 0.80,
) -> dict[str, float]:
    block = hist_block(branch, centrality_key)
    signal_counts = np.asarray(block["signal"]["counts"], dtype=float)
    background_counts = np.asarray(block["background"]["counts"], dtype=float)
    edges = np.asarray(branch["bin_edges"], dtype=float)
    sig_total = float(signal_counts.sum())
    bkg_total = float(background_counts.sum())
    if sig_total <= 0 or bkg_total <= 0:
        return {
            "threshold": float("nan"),
            "signal_efficiency": float("nan"),
            "background_fake_rate": float("nan"),
        }
    signal_above = np.cumsum(signal_counts[::-1])[::-1] / sig_total
    background_above = np.cumsum(background_counts[::-1])[::-1] / bkg_total
    idx = int(np.nanargmin(np.abs(signal_above - target_efficiency)))
    return {
        "threshold": float(edges[idx]),
        "signal_efficiency": float(signal_above[idx]),
        "background_fake_rate": float(background_above[idx]),
    }


def quantile_from_hist(counts: np.ndarray, edges: np.ndarray, quantile: float) -> float:
    total = float(np.sum(counts))
    if total <= 0:
        return float("nan")
    idx = int(np.searchsorted(np.cumsum(counts), quantile * total, side="left"))
    idx = max(0, min(idx, len(counts) - 1))
    return float(0.5 * (edges[idx] + edges[idx + 1]))


def median_bdt_score_gap(branch: dict, *, centrality_key: str | None = None) -> dict[str, float]:
    block = hist_block(branch, centrality_key)
    signal_counts = np.asarray(block["signal"]["counts"], dtype=float)
    background_counts = np.asarray(block["background"]["counts"], dtype=float)
    edges = np.asarray(branch["bin_edges"], dtype=float)
    signal_median = quantile_from_hist(signal_counts, edges, 0.50)
    background_median = quantile_from_hist(background_counts, edges, 0.50)
    return {
        "signal_median": signal_median,
        "background_median": background_median,
        "median_score_gap": signal_median - background_median,
    }


def pp_baseline_median_gaps(path: Path = DEFAULT_PP_BASELINE_JSON) -> dict[str, float | str]:
    payload = json.loads(path.read_text())
    edges = np.asarray(payload["bins"], dtype=float)

    def gap(signal_key: str, background_key: str) -> float:
        signal_hist = np.asarray(payload[signal_key], dtype=float)
        background_hist = np.asarray(payload[background_key], dtype=float)
        signal_hist = signal_hist / signal_hist.sum()
        background_hist = background_hist / background_hist.sum()
        return quantile_from_hist(signal_hist, edges, 0.50) - quantile_from_hist(background_hist, edges, 0.50)

    return {
        "this_analysis_gap": gap("this_analysis_signal_hist", "this_analysis_inclusive_hist"),
        "shuhang_reference_gap": gap("shuhang_signal_hist", "shuhang_inclusive_hist"),
        "source": str(path),
        "cuts": str(payload.get("cuts", {})),
    }


def hist_block(branch: dict, centrality_key: str | None = None) -> dict:
    if centrality_key:
        return branch["by_centrality"][centrality_key]
    return branch["inclusive"]


def formatted_sample_label(label: str) -> str:
    if label == "Jet12+20+30":
        return "Jet12+20\n+30"
    if label == "Jet12+20+30+40":
        return "Jet12+20\n+30+40"
    return label


def density_max(sample: dict, centrality_key: str | None = None) -> float:
    vals: list[float] = []
    for branch in sample["branches"]:
        block = hist_block(branch, centrality_key)
        vals.extend(block["signal"]["density"])
        vals.extend(block["background"]["density"])
    return max(vals) if vals else 1.0


def draw_cell(
    ax: plt.Axes,
    branch: dict,
    *,
    ymax: float,
    show_y: bool,
    centrality_key: str | None,
    annotate_wp80_fake: bool,
    annotate_median_gap: bool,
) -> dict[str, float]:
    edges = np.asarray(branch["bin_edges"], dtype=float)
    block = hist_block(branch, centrality_key)
    sig = np.asarray(block["signal"]["density"], dtype=float)
    inc = np.asarray(block["background"]["density"], dtype=float)
    xs, ys = step_xy(edges, sig)
    xb, yb = step_xy(edges, inc)
    ax.fill_between(xb, 0.0, yb, step="pre", color=INCLUSIVE, alpha=0.16, linewidth=0)
    ax.fill_between(xs, 0.0, ys, step="pre", color=SIGNAL, alpha=0.13, linewidth=0)
    ax.step(edges[:-1], inc, where="post", color=INCLUSIVE, lw=2.0)
    ax.step(edges[:-1], sig, where="post", color=SIGNAL, lw=2.2)
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, ymax)
    ax.set_xticks(np.linspace(0.0, 1.0, 6))
    ax.tick_params(direction="in", top=True, right=True, length=4.8, labelsize=10.8)
    if not show_y:
        ax.set_yticklabels([])
    else:
        ax.set_ylabel("Unit-area\ndensity", fontsize=11.6, color=INK, labelpad=1)
    ax.grid(True, color=GRID, alpha=0.72, lw=0.8)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_color("#344054")
        spine.set_linewidth(0.9)

    sig_n, inc_n = entry_counts(branch, centrality_key)
    auc = auc_of(branch, centrality_key)
    wp80 = wp_fake_rate_at_signal_eff(branch, centrality_key=centrality_key)
    med_gap = median_bdt_score_gap(branch, centrality_key=centrality_key)
    if annotate_median_gap:
        ax.text(
            0.045,
            0.925,
            f"Median BDT gap\n{med_gap['median_score_gap']:.2f}",
            transform=ax.transAxes,
            fontsize=13.2,
            fontweight="bold",
            ha="left",
            va="top",
            color="#111827",
            linespacing=1.04,
            bbox=dict(facecolor="#FFFFFF", edgecolor="#94A3B8", boxstyle="round,pad=0.22", alpha=0.96),
        )
    elif annotate_wp80_fake:
        ax.text(
            0.045,
            0.945,
            f"AUC {auc:.3f}",
            transform=ax.transAxes,
            fontsize=12.3,
            fontweight="bold",
            ha="left",
            va="top",
            color=INK,
            bbox=dict(facecolor="white", edgecolor="#D0D5DD", boxstyle="round,pad=0.18", alpha=0.94),
        )
        ax.text(
            0.045,
            0.735,
            f"WP80 fake {100.0 * wp80['background_fake_rate']:.1f}%",
            transform=ax.transAxes,
            fontsize=12.3,
            fontweight="bold",
            ha="left",
            va="top",
            color="#C2410C",
            bbox=dict(facecolor="#FFF7ED", edgecolor="#FDBA74", boxstyle="round,pad=0.18", alpha=0.96),
        )
        ax.text(
            0.045,
            0.535,
            f"S {sig_n:,}\nIncl {inc_n:,}",
            transform=ax.transAxes,
            fontsize=10.4,
            ha="left",
            va="top",
            color=MUTED,
            linespacing=1.12,
        )
    else:
        ax.text(
            0.045,
            0.935,
            f"AUC {auc:.3f}\nS {sig_n:,}\nIncl {inc_n:,}",
            transform=ax.transAxes,
            fontsize=11.0,
            ha="left",
            va="top",
            color=INK,
            linespacing=1.18,
            bbox=dict(facecolor="white", edgecolor="#D0D5DD", boxstyle="round,pad=0.28", alpha=0.94),
        )
    return {
        "auc": auc,
        "wp80_threshold": wp80["threshold"],
        "wp80_signal_efficiency": wp80["signal_efficiency"],
        "wp80_background_fake_rate": wp80["background_fake_rate"],
        "signal_median_bdt": med_gap["signal_median"],
        "background_median_bdt": med_gap["background_median"],
        "median_bdt_score_gap": med_gap["median_score_gap"],
        "signal_entries": sig_n,
        "background_entries": inc_n,
    }


def add_row_header_card(
    fig: plt.Figure,
    *,
    xy: tuple[float, float],
    wh: tuple[float, float],
    role: str,
    label: str,
    accent: dict[str, str],
) -> None:
    x, y = xy
    w, h = wh
    add_box(fig, xy, wh, accent["fill"], edge=accent["edge"], lw=1.1, radius=0.008)
    fig.patches.append(
        Rectangle(
            (x + 0.006, y + 0.010),
            0.005,
            max(0.001, h - 0.020),
            transform=fig.transFigure,
            facecolor=accent["strip"],
            edgecolor="none",
            zorder=-1,
        )
    )
    add_text(
        fig,
        x + 0.5 * w + 0.009,
        y + 0.5 * h + 0.028,
        role,
        size=11.6,
        weight="bold",
        ha="center",
        va="center",
        color=NAVY,
    )
    sample_label = formatted_sample_label(label)
    add_text(
        fig,
        x + 0.5 * w + 0.009,
        y + 0.5 * h - 0.008,
        sample_label,
        size=13.2 if "\n" in sample_label else 13.8,
        weight="bold",
        ha="center",
        va="center",
        color=NAVY,
        linespacing=1.18,
    )


def validation_scope_label(payload: dict) -> str:
    first = payload["common_samples"][0]["branches"][0]["summary"]
    row_scope = first.get("row_scope", "")
    if "direct_model_scoring" in row_scope:
        return "10% held-out validation rows, directly scored by each trained model"
    if "holdout" in row_scope:
        return "10% held-out validation rows from the source score caches"
    return "full score-cache validation rows"


def build_slide(
    payload: dict,
    outpng: Path,
    *,
    centrality_key: str | None = None,
    centrality_label: str | None = None,
    annotate_wp80_fake: bool = False,
    annotate_median_gap: bool = False,
) -> Path:
    outpng.parent.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
            "axes.unicode_minus": False,
        }
    )
    fig = plt.figure(figsize=slide_figsize(SLIDE_DPI), dpi=SLIDE_DPI, facecolor="white")
    fig.patches.append(Rectangle((0, 0), 1, 1, transform=fig.transFigure, facecolor="white", zorder=-10))

    add_text(
        fig,
        0.040,
        0.958,
        f"{centrality_label or '0-80%'} centrality: fixed validation sample vs trained BDT",
        size=28.5,
        weight="bold",
        color=NAVY,
    )
    add_box(fig, (0.050, 0.812), (0.425, 0.078), SOFT_GRAY)
    add_text(fig, 0.070, 0.875, "Validation rows", size=17.6, weight="bold", color=NAVY)
    add_text(
        fig,
        0.070,
        0.837,
        "Compare models on the same holdout sample.",
        size=14.2,
        color=INK,
        linespacing=1.24,
    )
    add_box(fig, (0.525, 0.812), (0.425, 0.078), SOFT_LAVENDER, edge="#CAB8DA")
    add_text(fig, 0.545, 0.875, "Training columns", size=17.6, weight="bold", color=NAVY)
    add_text(
        fig,
        0.545,
        0.837,
        "Move left to right to change only training input.",
        size=14.2,
        color=INK,
        linespacing=1.34,
    )

    add_box(fig, (0.222, 0.742), (0.710, 0.041), "white", edge="#A9B8CC", lw=1.0, radius=0.008)
    fig.lines.append(plt.Line2D([0.262, 0.304], [0.762, 0.762], transform=fig.transFigure, color=SIGNAL, lw=3.0))
    add_text(
        fig,
        0.314,
        0.772,
        "Signal MC (truth-isolated prompt)",
        size=12.7,
        weight="bold",
        color=INK,
        va="top",
    )
    fig.lines.append(plt.Line2D([0.618, 0.660], [0.762, 0.762], transform=fig.transFigure, color=INCLUSIVE, lw=3.0))
    add_text(
        fig,
        0.670,
        0.772,
        "Inclusive MC (embedded jet; no truth filter)",
        size=12.7,
        weight="bold",
        color=INK,
        va="top",
    )

    samples = payload["common_samples"]
    cols = [branch["label"] for branch in samples[0]["branches"]]
    summary_rows: list[dict[str, object]] = []
    left, right = 0.185, 0.968
    bottom, top = 0.198, 0.647
    wspace, hspace = 0.026, 0.072
    ncols, nrows = 3, len(samples)
    if nrows >= 3:
        hspace = 0.038
    ax_w = (right - left - wspace * (ncols - 1)) / ncols
    ax_h = (top - bottom - hspace * (nrows - 1)) / nrows

    for c, col in enumerate(cols):
        accent = SAMPLE_ACCENTS[c % len(SAMPLE_ACCENTS)]
        x_header = left + c * (ax_w + wspace)
        add_box(fig, (x_header + 0.004, 0.672), (ax_w - 0.008, 0.035), accent["fill"], edge=accent["edge"], lw=1.1, radius=0.007)
        add_text(
            fig,
            left + c * (ax_w + wspace) + ax_w / 2,
            0.699,
            f"TRAINING: {col}",
            size=13.1,
            weight="bold",
            ha="center",
            color=NAVY,
        )

    for r, sample in enumerate(samples):
        y0 = top - (r + 1) * ax_h - r * hspace
        row_accent = SAMPLE_ACCENTS[r % len(SAMPLE_ACCENTS)]
        add_row_header_card(
            fig,
            xy=(0.026, y0 + 0.012),
            wh=(0.096, ax_h - 0.024),
            role="VALIDATION",
            label=sample["validation_sample"],
            accent=row_accent,
        )
        ymax = max(1.0, density_max(sample, centrality_key) * 1.34)
        for c, branch in enumerate(sample["branches"]):
            x0 = left + c * (ax_w + wspace)
            ax = fig.add_axes([x0, y0, ax_w, ax_h])
            metrics = draw_cell(
                ax,
                branch,
                ymax=ymax,
                show_y=(c == 0),
                centrality_key=centrality_key,
                annotate_wp80_fake=annotate_wp80_fake,
                annotate_median_gap=annotate_median_gap,
            )
            summary_rows.append(
                {
                    "validation_sample": sample["validation_sample"],
                    "training_sample": branch["label"],
                    "centrality": centrality_label or centrality_key or "0-80%",
                    "auc": metrics["auc"],
                    "wp80_threshold_binned": metrics["wp80_threshold"],
                    "wp80_signal_efficiency_binned": metrics["wp80_signal_efficiency"],
                    "wp80_background_fake_rate_binned": metrics["wp80_background_fake_rate"],
                    "signal_median_bdt_binned": metrics["signal_median_bdt"],
                    "background_median_bdt_binned": metrics["background_median_bdt"],
                    "median_bdt_score_gap_binned": metrics["median_bdt_score_gap"],
                    "signal_entries": metrics["signal_entries"],
                    "background_entries": metrics["background_entries"],
                }
            )
            if r == nrows - 1:
                ax.set_xlabel("BDT score", fontsize=12.2, labelpad=4)
            else:
                ax.set_xticklabels([])

    row_gains = []
    jet40_steps = []
    for sample in samples:
        aucs = [auc_of(branch, centrality_key) for branch in sample["branches"]]
        if len(aucs) >= 3:
            row_gains.append(aucs[-1] - aucs[0])
            jet40_steps.append(aucs[-1] - aucs[-2])
    gains_text = " / ".join(f"+{gain:.3f}" for gain in row_gains)
    max_jet40_step = max(jet40_steps) if jet40_steps else float("nan")
    median_gap_steps = []
    for sample in samples:
        gaps = [median_bdt_score_gap(branch, centrality_key=centrality_key)["median_score_gap"] for branch in sample["branches"]]
        if len(gaps) >= 3:
            median_gap_steps.append((gaps[0], gaps[-1]))
    median_gap_text = " / ".join(f"{left_gap:.2f}->{right_gap:.2f}" for left_gap, right_gap in median_gap_steps)
    pp_baseline = pp_baseline_median_gaps() if annotate_median_gap else None

    add_box(fig, (0.050, 0.022), (0.900, 0.118), "#FFF7ED", edge="#FDBA74", lw=1.1, radius=0.010)
    if annotate_median_gap:
        add_text(fig, 0.070, 0.130, "What median BDT gap means", size=15.0, weight="bold", color="#9A3412")
        add_text(
            fig,
            0.070,
            0.103,
            "Definition: median(signal BDT) - median(inclusive-jet BDT).",
            size=12.8,
            color=INK,
            linespacing=1.22,
        )
        add_text(
            fig,
            0.070,
            0.079,
            "Physical read: larger gap = photon-like showers sit farther from jet-like candidates.",
            size=12.6,
            color=MUTED,
            linespacing=1.22,
        )
        if pp_baseline is not None:
            add_text(
                fig,
                0.070,
                0.048,
                "Right callout: pp trained baseV3E split in the same 15-35 GeV pT range.",
                size=11.9,
                color="#7C2D12",
                linespacing=1.20,
            )
        add_box(fig, (0.700, 0.032), (0.225, 0.096), "#FFFFFF", edge="#CBD5E1", lw=1.0, radius=0.008, zorder=-1)
        add_text(fig, 0.718, 0.119, "pp trained baseV3E", size=11.6, weight="bold", color=NAVY)
        add_text(fig, 0.718, 0.100, "same 15-35 GeV pT range", size=9.8, color=MUTED)
        add_text(fig, 0.718, 0.083, "raw-inclusive pp validation", size=9.8, color=MUTED)
        add_text(fig, 0.718, 0.063, "median BDT gap", size=10.2, weight="bold", color=INK)
        add_text(fig, 0.866, 0.068, f"{pp_baseline['this_analysis_gap']:.2f}", size=25.0, weight="bold", color="#9A3412", ha="center")
    else:
        add_text(
            fig,
            0.070,
            0.124,
            f"What the {centrality_label or '0-80%'} matrix shows",
            size=15.0,
            weight="bold",
            color="#9A3412",
        )
        add_text(
            fig,
            0.070,
            0.095,
            f"AUC gain from Jet12+20 to Jet12+20+30+40 training: {gains_text} across the validation rows.",
            size=12.8,
            color=INK,
            linespacing=1.26,
        )
        add_text(
            fig,
            0.070,
            0.068,
            (
                "WP80 fake rate is mostly flat within a fixed validation row; "
                f"Jet40 remains a saturation check (max +Jet40 step {max_jet40_step:.3f} AUC)."
                if annotate_wp80_fake
                else f"Conceptually: Jet30 provides the visible central-bin separation gain; Jet40 is a saturation check (max +Jet40 step {max_jet40_step:.3f} AUC)."
            ),
            size=12.9,
            color=MUTED,
            linespacing=1.32,
        )

    fig.savefig(outpng, dpi=SLIDE_DPI)
    plt.close(fig)
    if summary_rows:
        if annotate_median_gap:
            csv_path = outpng.with_suffix(".median_gap_summary.csv")
        elif annotate_wp80_fake:
            csv_path = outpng.with_suffix(".wp80_summary.csv")
        else:
            csv_path = outpng.with_suffix(".summary.csv")
        with csv_path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0]))
            writer.writeheader()
            writer.writerows(summary_rows)
    return outpng


def write_speaker_script(
    payload: dict,
    outpng: Path,
    input_json: Path,
    *,
    centrality_key: str | None = None,
    centrality_label: str | None = None,
    annotate_median_gap: bool = False,
) -> Path:
    outmd = outpng.with_suffix(".speaker_script.md")
    sample_lines = []
    for sample in payload["common_samples"]:
        if annotate_median_gap:
            metrics = ", ".join(
                f"{branch['label']} median BDT gap {median_bdt_score_gap(branch, centrality_key=centrality_key)['median_score_gap']:.2f}"
                for branch in sample["branches"]
            )
        else:
            metrics = ", ".join(
                (
                    f"{branch['label']} AUC {auc_of(branch, centrality_key):.3f}, "
                    f"WP80 fake {100.0 * wp_fake_rate_at_signal_eff(branch, centrality_key=centrality_key)['background_fake_rate']:.1f}%"
                )
                for branch in sample["branches"]
            )
        sample_lines.append(f"- {sample['validation_sample']}: {metrics}")
    interpretation = (
        "The important check is that the same-sample score curves do visibly polarize as the training sample expands. "
        "The median BDT score gap captures that: it is larger left-to-right because the signal median moves farther above the inclusive-jet median. "
        "Physically, that means the BDT is putting photon-like shower-shape patterns farther to the high-score side than the inclusive jet-like population. "
        "This is a useful sign that the model is using a more decisive score scale, but it is not by itself a tighter working-point rejection claim."
        if annotate_median_gap
        else "The important check is that the top row does not show a large WP80 fake-rate drop when we keep the Jet12+20 validation sample fixed. The score distributions shift, but at a re-derived 80% signal working point the background acceptance is roughly flat."
    )
    summary_title = "Binned median BDT score-gap summary:" if annotate_median_gap else "Weighted source-defined AUC and binned WP80 fake-rate summary:"
    pp_lines: list[str] = []
    if annotate_median_gap:
        pp_baseline = pp_baseline_median_gaps()
        pp_lines = [
            "",
            "Same-range pp trained-BDT reference:",
            f"- Our trained pp baseV3E raw-inclusive 15-35 GeV median BDT gap: {pp_baseline['this_analysis_gap']:.2f}",
            "- Caveat: the same-pT pp artifact is the available raw-inclusive 15-35 GeV contract; the later truth-window-fixed reference is available at 22-28 GeV, not this exact range.",
        ]
    outmd.write_text(
        "\n".join(
            [
                "# THE-8 Branch A fixed-sample validation slide",
                "",
                "Speaker framing:",
                f"This slide shows only the {centrality_label or '0-80%'} centrality bin.",
                "We are holding the validation population fixed by row and changing only which Branch A BDT is applied by column.",
                f"Validation scope: {validation_scope_label(payload)}.",
                "The red curve is truth-isolated prompt signal MC from embedded-photon rows. The blue curve is the inclusive embedded-jet MC candidate population with no truth-background filter.",
                "That means the legend terms are literal: Signal MC is the clean truth-isolated prompt-photon template, and Inclusive MC is the generic inclusive-jet candidate population after selections.",
                interpretation,
                "",
                summary_title,
                *sample_lines,
                *pp_lines,
                "",
                f"PNG: {outpng}",
                f"JSON: {input_json}",
            ]
        )
        + "\n"
    )
    return outmd


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-json", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--outpng", type=Path, default=DEFAULT_OUTPNG)
    parser.add_argument("--centrality-key", default=None, choices=[None, "0_20", "20_50", "50_80"])
    parser.add_argument("--centrality-label", default=None)
    parser.add_argument(
        "--annotate-wp80-fake",
        action="store_true",
        help="Print binned WP80 background fake rate below each cell AUC.",
    )
    parser.add_argument(
        "--annotate-median-gap",
        action="store_true",
        help="Replace cell AUC labels with binned median BDT score gap.",
    )
    args = parser.parse_args()
    payload = json.loads(args.input_json.read_text())
    build_slide(
        payload,
        args.outpng,
        centrality_key=args.centrality_key,
        centrality_label=args.centrality_label,
        annotate_wp80_fake=args.annotate_wp80_fake,
        annotate_median_gap=args.annotate_median_gap,
    )
    write_speaker_script(
        payload,
        args.outpng,
        args.input_json,
        centrality_key=args.centrality_key,
        centrality_label=args.centrality_label,
        annotate_median_gap=args.annotate_median_gap,
    )
    print(args.outpng)
    print(args.outpng.with_suffix(".speaker_script.md"))
    if args.annotate_median_gap:
        print(args.outpng.with_suffix(".median_gap_summary.csv"))
    elif args.annotate_wp80_fake:
        print(args.outpng.with_suffix(".wp80_summary.csv"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
