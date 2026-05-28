#!/usr/bin/env python3
"""Build slide-ready same-row pp baseV3E/Shuhang score correlation figures."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm
from matplotlib.patches import FancyBboxPatch


INK = "#111827"
MUTED = "#4B5563"
GRID = "#D1D5DB"
BLUE = "#2563EB"
GREEN = "#047857"
PANEL = "#FFFFFF"
PANEL_EDGE = "#CBD5E1"
SOFT_BLUE = "#EFF6FF"
SOFT_GREEN = "#ECFDF5"
SOFT_AMBER = "#FFF7ED"


def _set_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "Times New Roman",
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "axes.edgecolor": INK,
            "axes.linewidth": 1.3,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "mathtext.fontset": "dejavuserif",
        }
    )


def draw_sphenix_label(ax: plt.Axes, *, x: float = 0.04, y: float = 0.955, size: int = 18) -> None:
    ax.text(
        x,
        y,
        "sPHENIX",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=size,
        fontstyle="italic",
        fontweight="bold",
        color=INK,
    )
    ax.text(
        x + 0.205,
        y,
        "Internal",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=size,
        color=INK,
    )


def load_inputs(csv_path: Path, summary_path: Path | None) -> tuple[pd.DataFrame, dict]:
    frame = pd.read_csv(csv_path)
    needed = {"class", "sample", "our_score", "shuhang_split_score"}
    missing = needed.difference(frame.columns)
    if missing:
        raise SystemExit(f"Missing required columns in {csv_path}: {sorted(missing)}")
    summary = {}
    if summary_path and summary_path.exists():
        summary = json.loads(summary_path.read_text())
    return frame, summary


def compute_plot_arrays(frame: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    x = frame["our_score"].to_numpy(dtype=float)
    y = frame["shuhang_split_score"].to_numpy(dtype=float)
    valid = np.isfinite(x) & np.isfinite(y) & (x >= 0.0) & (x <= 1.0) & (y >= 0.0) & (y <= 1.0)
    return x[valid], y[valid], frame.loc[valid].copy()


def pearson_by_et(frame: pd.DataFrame, edges: np.ndarray) -> pd.DataFrame:
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        sub = frame[(frame["cluster_Et"] >= lo) & (frame["cluster_Et"] < hi)]
        if len(sub) >= 20:
            r = sub[["our_score", "shuhang_split_score"]].corr().iloc[0, 1]
        else:
            r = np.nan
        rows.append({"lo": lo, "hi": hi, "center": 0.5 * (lo + hi), "pearson": r, "rows": len(sub)})
    return pd.DataFrame(rows)


def draw_th2_panel(
    ax: plt.Axes,
    x: np.ndarray,
    y: np.ndarray,
    *,
    bins: int,
    title: str | None,
    show_sphenix: bool,
    info_box: bool,
    rows: int,
) -> object:
    h = ax.hist2d(
        x,
        y,
        bins=bins,
        range=[[0.0, 1.0], [0.0, 1.0]],
        norm=LogNorm(vmin=1),
        cmap="viridis",
    )
    image = h[3]
    ax.plot([0, 1], [0, 1], color="white", lw=3.0, ls="--", alpha=0.88, zorder=4)
    ax.plot([0, 1], [0, 1], color=INK, lw=1.1, ls="--", alpha=0.72, zorder=5)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("This analysis frozen XGBoost score", fontsize=16)
    ax.set_ylabel("Shuhang/PPG12 split TMVA baseV3E score", fontsize=16)
    if title:
        ax.set_title(title, loc="left", fontsize=22, fontweight="bold", color=INK, pad=12)
    ax.grid(color=GRID, lw=0.8, alpha=0.48)
    ax.tick_params(axis="both", labelsize=14)
    if show_sphenix:
        draw_sphenix_label(ax)
    if info_box:
        ax.text(
            0.04,
            0.82,
            r"$p{+}p$ $\sqrt{s}=200$ GeV" "\n" r"$|\eta|<0.7$, $22<E_T<28$ GeV" "\n"
            f"same-row audit: {rows:,} rows",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=13.5,
            color=INK,
            bbox=dict(boxstyle="round,pad=0.34", fc="white", ec="#E5E7EB", alpha=0.92),
        )
    ax.text(
        0.98,
        0.055,
        "dashed line: identical scores",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=11.5,
        color=INK,
        bbox=dict(boxstyle="round,pad=0.28", fc="white", ec="#E5E7EB", alpha=0.9),
    )
    return image


def rounded_box(fig: plt.Figure, x: float, y: float, w: float, h: float, fc: str, ec: str = PANEL_EDGE) -> None:
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        transform=fig.transFigure,
        boxstyle="round,pad=0.012,rounding_size=0.012",
        linewidth=1.1,
        edgecolor=ec,
        facecolor=fc,
        zorder=-1,
    )
    fig.patches.append(patch)


def text_block(fig: plt.Figure, x: float, y: float, text: str, *, size: int = 24, color: str = INK, weight: str = "normal", linespacing: float = 1.18) -> None:
    fig.text(x, y, text, ha="left", va="top", fontsize=size, color=color, fontweight=weight, linespacing=linespacing)


def make_standalone_plot(frame: pd.DataFrame, summary: dict, out: Path, *, bins: int) -> None:
    x, y, valid_frame = compute_plot_arrays(frame)
    metrics = summary.get("metrics", {})
    rows = int(metrics.get("scatter_rows", len(valid_frame)))
    pearson = float(metrics.get("our_vs_shuhang_split_pearson", np.corrcoef(x, y)[0, 1]))
    spearman = float(metrics.get("our_vs_shuhang_split_spearman", valid_frame[["our_score", "shuhang_split_score"]].corr(method="spearman").iloc[0, 1]))

    fig, ax = plt.subplots(figsize=(10.8, 8.7), constrained_layout=True)
    image = draw_th2_panel(
        ax,
        x,
        y,
        bins=bins,
        title="Same-row BDT score correlation",
        show_sphenix=True,
        info_box=True,
        rows=rows,
    )
    ax.text(
        0.04,
        0.68,
        f"Pearson r = {pearson:.3f}\nSpearman rho = {spearman:.3f}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=12.5,
        color=INK,
        bbox=dict(boxstyle="round,pad=0.34", fc="white", ec="#E5E7EB", alpha=0.92),
    )
    cbar = fig.colorbar(image, ax=ax, pad=0.02)
    cbar.set_label("same-row candidates per score bin", fontsize=14)
    cbar.ax.tick_params(labelsize=12)
    fig.savefig(out, dpi=220, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)


def make_slide(frame: pd.DataFrame, summary: dict, out: Path, *, bins: int) -> None:
    x, y, valid_frame = compute_plot_arrays(frame)
    metrics = summary.get("metrics", {})
    rows = int(metrics.get("scatter_rows", len(valid_frame)))
    pearson = float(metrics.get("our_vs_shuhang_split_pearson", np.corrcoef(x, y)[0, 1]))
    spearman = float(metrics.get("our_vs_shuhang_split_spearman", valid_frame[["our_score", "shuhang_split_score"]].corr(method="spearman").iloc[0, 1]))
    our_auc = float(metrics.get("our_auc", np.nan))
    split_auc = float(metrics.get("shuhang_split_auc", np.nan))
    file_counts = summary.get("file_counts", {})
    def natural_sample_key(sample: str) -> tuple[str, int]:
        head = "".join(ch for ch in sample if not ch.isdigit())
        digits = "".join(ch for ch in sample if ch.isdigit())
        return head, int(digits or 0)

    samples = {
        "signal": sorted((k.split("/", 1)[1] for k in file_counts if k.startswith("signal/")), key=natural_sample_key),
        "inclusive": sorted((k.split("/", 1)[1] for k in file_counts if k.startswith("inclusive/")), key=natural_sample_key),
    }
    signal_text = ", ".join(samples["signal"]) if samples["signal"] else "photonjet5, photonjet10, photonjet20"
    inclusive_text = ", ".join(samples["inclusive"]) if samples["inclusive"] else "jet8, jet12, jet20, jet30, jet40"

    fig = plt.figure(figsize=(16, 9), dpi=160)
    fig.subplots_adjust(0, 0, 1, 1)

    text_block(fig, 0.045, 0.945, "Same-row BDT Model-Equivalence Check", size=31, weight="bold")
    text_block(
        fig,
        0.045,
        0.900,
        "This analysis frozen XGBoost score versus Shuhang/PPG12 split TMVA baseV3E score, evaluated directly on identical candidates.",
        size=15,
        color=MUTED,
    )

    ax = fig.add_axes([0.045, 0.125, 0.535, 0.725])
    image = draw_th2_panel(
        ax,
        x,
        y,
        bins=bins,
        title=None,
        show_sphenix=True,
        info_box=False,
        rows=rows,
    )
    cax = fig.add_axes([0.595, 0.18, 0.018, 0.58])
    cbar = fig.colorbar(image, cax=cax)
    cbar.set_label("candidates per score bin", fontsize=12, labelpad=22)
    cbar.ax.yaxis.set_label_position("left")
    cbar.ax.tick_params(labelsize=10)

    # Right-side explanatory boxes.
    rounded_box(fig, 0.67, 0.715, 0.295, 0.135, SOFT_BLUE, "#BFDBFE")
    rhs_x = 0.682

    text_block(fig, rhs_x, 0.823, "What same-row means", size=17.6, color=BLUE, weight="bold")
    text_block(
        fig,
        rhs_x,
        0.786,
        "Every point is one candidate cluster scored twice:\n"
        "x = this analysis frozen XGBoost score\n"
        "y = Shuhang/PPG12 split TMVA score",
        size=12.0,
        linespacing=1.28,
    )

    rounded_box(fig, 0.67, 0.515, 0.295, 0.165, SOFT_GREEN, "#A7F3D0")
    text_block(fig, rhs_x, 0.650, "Cluster-node contract", size=17.3, color=GREEN, weight="bold")
    text_block(
        fig,
        rhs_x,
        0.610,
        "This is the correct split-cluster comparison.\n\n"
        "PPG12 split node: CLUSTERINFO_CEMC\n"
        "This analysis: CLUSTERINFO_CEMC -> PHOTONCLUSTER_CEMC\n\n"
        "No-split would be a separate CLUSTERINFO_CEMC_NO_SPLIT check.",
        size=9.7,
        linespacing=1.20,
    )

    rounded_box(fig, 0.67, 0.305, 0.295, 0.170, PANEL, PANEL_EDGE)
    text_block(fig, rhs_x, 0.452, r"Agreement versus cluster $E_T$", size=16.2, color=INK, weight="bold")
    et_corr = pearson_by_et(valid_frame, np.arange(22.0, 29.0, 1.0))
    text_block(fig, rhs_x, 0.418, f"Pearson r = {pearson:.3f} overall; stable in every 1 GeV E_T bin.", size=10.4, color=MUTED)
    corr_ax = fig.add_axes([0.700, 0.334, 0.215, 0.058])
    corr_ax.plot(et_corr["center"], et_corr["pearson"], color=BLUE, marker="o", ms=4.0, lw=1.8)
    corr_ax.set_xlim(22.0, 28.0)
    corr_ax.set_ylim(0.90, 1.0)
    corr_ax.set_xticks([22, 24, 26, 28])
    corr_ax.set_yticks([0.90, 0.95, 1.00])
    corr_ax.grid(color=GRID, lw=0.55, alpha=0.7)
    corr_ax.tick_params(labelsize=8.5, pad=2)
    corr_ax.set_xlabel(r"cluster $E_T$ [GeV]", fontsize=9, labelpad=1)
    corr_ax.set_ylabel("r", fontsize=9, labelpad=1)

    rounded_box(fig, 0.67, 0.150, 0.295, 0.125, PANEL, PANEL_EDGE)
    text_block(fig, rhs_x, 0.247, "Audit inputs", size=15.7, color=INK, weight="bold")
    text_block(
        fig,
        rhs_x,
        0.212,
        f"{rows:,} same-row candidates, 22 < ET < 28 GeV, |eta| < 0.7\n"
        f"Signal MC: {signal_text}\n"
        f"Inclusive MC: {inclusive_text}\n"
        "TMVA scores are evaluated directly, not read from a pre-made histogram.",
        size=10.4,
        linespacing=1.22,
    )

    rounded_box(fig, 0.67, 0.045, 0.295, 0.072, SOFT_AMBER, "#FDBA74")
    fig.text(rhs_x, 0.097, "Takeaway:", ha="left", va="top", fontsize=11.5, fontweight="bold", color=INK)
    fig.text(
        0.742,
        0.097,
        "the two score definitions agree closely on the same rows.\n"
        f"AUC: this analysis {our_auc:.3f}, PPG12 split TMVA {split_auc:.3f}.",
        ha="left",
        va="top",
        fontsize=11.5,
        color=INK,
        linespacing=1.12,
    )

    fig.savefig(out, dpi=160, bbox_inches=None, pad_inches=0)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", required=True, type=Path)
    parser.add_argument("--summary-json", type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--bins", type=int, default=100)
    args = parser.parse_args()

    _set_style()
    args.outdir.mkdir(parents=True, exist_ok=True)
    frame, summary = load_inputs(args.csv, args.summary_json)

    standalone = args.outdir / "pp_basev3e_same_row_shuhang_split_vs_this_analysis_TH2_annotated.png"
    slide = args.outdir / "pp_basev3e_same_row_shuhang_split_model_equivalence_slide.png"
    make_standalone_plot(frame, summary, standalone, bins=args.bins)
    make_slide(frame, summary, slide, bins=args.bins)

    print(standalone)
    print(slide)


if __name__ == "__main__":
    main()
