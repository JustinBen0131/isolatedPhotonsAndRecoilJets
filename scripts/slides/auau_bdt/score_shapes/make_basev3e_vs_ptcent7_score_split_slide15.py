#!/usr/bin/env python3
"""Build a full-slide replacement PNG for the slide-15 routed-BDT score split."""

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

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import font_manager


ROOT = Path(
    "dataOutput/auauMLDiagnosticRuns/"
    "global_etcent_inclusive3_sixpack_20260516_135439"
)
SOURCE_DIR = ROOT / "slideReady/basev3e_vs_best_bdt_score_split"
SOURCE_CSV = SOURCE_DIR / "basev3e_withCentrality_vs_best_ptCent7_score_split_coarse3_2x3.csv"
OUT_PNG = SOURCE_DIR / "basev3e_withCentrality_vs_best_ptCent7_score_split_slide15_clean.png"

TIMES_FONT_FILES = [
    Path("/System/Library/Fonts/Supplemental/Times New Roman.ttf"),
    Path("/System/Library/Fonts/Supplemental/Times New Roman Bold.ttf"),
    Path("/System/Library/Fonts/Supplemental/Times New Roman Italic.ttf"),
    Path("/System/Library/Fonts/Supplemental/Times New Roman Bold Italic.ttf"),
]
TIMES_FAMILY = "Times New Roman"

INK = "#111827"
MUTED = "#374151"
LIGHT_MUTED = "#64748B"
GRID = "#E5E7EB"
SIGNAL = "#1F77B4"
BACKGROUND = "#D95F02"
ROW_SHADE = "#F8FAFC"

MODEL_ORDER = [
    ("baseBDT_v3E_withCentrality", "Base v3E\n+ centrality", "1 BDT", "12 inputs"),
    ("globalEtCent1535_bdt_noIso_ptCent7", r"Routed BDT" "\n" r"$E_T \times$ centrality", "56 BDTs", "32 inputs"),
]
CENT_ORDER = ["0-20%", "20-50%", "50-80%"]


def register_times_new_roman() -> None:
    missing = [path for path in TIMES_FONT_FILES if not path.exists()]
    if missing:
        raise RuntimeError(f"Missing Times New Roman font files: {missing}")
    for path in TIMES_FONT_FILES:
        font_manager.fontManager.addfont(str(path))


def setup_matplotlib() -> None:
    register_times_new_roman()
    plt.rcParams.update(
        {
            "font.family": TIMES_FAMILY,
            "font.serif": [TIMES_FAMILY],
            "mathtext.fontset": "custom",
            "mathtext.rm": TIMES_FAMILY,
            "mathtext.it": f"{TIMES_FAMILY}:italic",
            "mathtext.bf": f"{TIMES_FAMILY}:bold",
            "axes.unicode_minus": False,
        }
    )


def add_header(fig: plt.Figure) -> None:
    fig.text(
        0.025,
        0.962,
        r"$E_T \times$ centrality routing sharpens BDT score separation",
        ha="left",
        va="top",
        fontsize=28.5,
        fontweight="bold",
        color=INK,
    )
    fig.text(
        0.027,
        0.902,
        r"Same validation sample and score plot in each panel: Photon12+20 signal vs Jet12+20+30 inclusive background, $15<E_T<35$ GeV.",
        ha="left",
        va="top",
        fontsize=15.2,
        color=MUTED,
    )
    fig.text(
        0.027,
        0.862,
        "Area-normalized score density: blue = signal, orange = background. More signal at high BDT score means cleaner photon ranking.",
        ha="left",
        va="top",
        fontsize=15.2,
        color=MUTED,
    )
    fig.text(
        0.865,
        0.956,
        "sPHENIX",
        ha="right",
        va="top",
        fontsize=18.5,
        fontstyle="italic",
        fontweight="bold",
        color=INK,
    )
    fig.text(
        0.868,
        0.956,
        " Internal",
        ha="left",
        va="top",
        fontsize=18.5,
        color=INK,
    )


def add_row_label(fig: plt.Figure, y: float, model_label: str, bdt_count: str, input_count: str, color: str) -> None:
    fig.text(
        0.050,
        y + 0.070,
        model_label,
        ha="left",
        va="center",
        fontsize=17.8,
        fontweight="bold",
        color=INK,
        linespacing=0.95,
    )
    fig.text(
        0.050,
        y - 0.020,
        f"{bdt_count}\n{input_count}",
        ha="left",
        va="center",
        fontsize=15.6,
        color=MUTED,
        linespacing=1.12,
    )
    fig.add_artist(
        plt.Rectangle(
            (0.036, y - 0.080),
            0.006,
            0.160,
            transform=fig.transFigure,
            facecolor=color,
            edgecolor=color,
            linewidth=0,
        )
    )


def style_axis(ax: plt.Axes) -> None:
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 20.0)
    ax.grid(True, color=GRID, linewidth=0.75)
    ax.tick_params(direction="in", top=True, right=True, labelsize=10.8, pad=2)
    for spine in ax.spines.values():
        spine.set_color("#1F2937")
        spine.set_linewidth(1.0)


def plot_panel(ax: plt.Axes, panel: pd.DataFrame, *, show_ylabel: bool, show_xlabel: bool) -> None:
    panel = panel.sort_values("score_center")
    ax.step(
        panel["score_center"],
        panel["signal_density"],
        where="mid",
        color=SIGNAL,
        linewidth=2.4,
        label="Signal",
    )
    ax.step(
        panel["score_center"],
        panel["background_density"],
        where="mid",
        color=BACKGROUND,
        linewidth=2.4,
        label="Background",
    )
    style_axis(ax)
    if show_ylabel:
        ax.set_ylabel("Density", fontsize=12.5, labelpad=5)
    else:
        ax.set_yticklabels([])
    if show_xlabel:
        ax.set_xticks([0.2, 0.4, 0.6, 0.8, 1.0])
        ax.set_xlabel("BDT score", fontsize=12.5, labelpad=5)
    else:
        ax.set_xticklabels([])

    auc = float(panel["auc"].iloc[0])
    ax.text(
        0.045,
        0.885,
        f"AUC {auc:.3f}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=13.5,
        fontweight="bold",
        color=INK,
        bbox={
            "facecolor": "white",
            "edgecolor": "#D1D5DB",
            "boxstyle": "round,pad=0.22",
            "linewidth": 0.9,
            "alpha": 0.96,
        },
    )


def add_footer(fig: plt.Figure, df: pd.DataFrame) -> None:
    aucs = (
        df.groupby(["model_key", "cent_plot"], sort=False)["auc"]
        .first()
        .unstack("cent_plot")
        .reindex(index=[row[0] for row in MODEL_ORDER], columns=CENT_ORDER)
    )
    gains = 100.0 * (aucs.iloc[1] / aucs.iloc[0] - 1.0)
    gain_text = ", ".join(f"{cent}: +{gains[cent]:.1f}% AUC" for cent in CENT_ORDER)

    fig.text(
        0.027,
        0.065,
        f"Routing gain by centrality: {gain_text}.",
        ha="left",
        va="bottom",
        fontsize=15.7,
        fontweight="bold",
        color=INK,
    )
    fig.text(
        0.027,
        0.034,
        r"Interpretation: the global BDT learns one average response; the routed BDT lets each $E_T$ and centrality region learn its local signal/background shape.",
        ha="left",
        va="bottom",
        fontsize=13.8,
        color=LIGHT_MUTED,
    )


def main() -> None:
    if not SOURCE_CSV.exists():
        raise FileNotFoundError(SOURCE_CSV)
    setup_matplotlib()
    df = pd.read_csv(SOURCE_CSV)

    fig = plt.figure(figsize=(16.0, 9.0), dpi=200)
    fig.patch.set_facecolor("white")
    add_header(fig)

    left = 0.225
    width = 0.225
    gap = 0.030
    top_y = 0.515
    bottom_y = 0.205
    height = 0.235
    axes = []

    fig.add_artist(
        plt.Rectangle(
            (0.018, top_y - 0.035),
            0.957,
            height + 0.088,
            transform=fig.transFigure,
            facecolor=ROW_SHADE,
            edgecolor="none",
            zorder=-1,
        )
    )
    fig.add_artist(
        plt.Rectangle(
            (0.018, bottom_y - 0.035),
            0.957,
            height + 0.088,
            transform=fig.transFigure,
            facecolor="white",
            edgecolor="none",
            zorder=-1,
        )
    )

    add_row_label(fig, top_y + 0.118, MODEL_ORDER[0][1], MODEL_ORDER[0][2], MODEL_ORDER[0][3], "#94A3B8")
    add_row_label(fig, bottom_y + 0.118, MODEL_ORDER[1][1], MODEL_ORDER[1][2], MODEL_ORDER[1][3], "#2563EB")

    for col, cent in enumerate(CENT_ORDER):
        x = left + col * (width + gap)
        fig.text(
            x + width / 2,
            top_y + height + 0.038,
            cent,
            ha="center",
            va="bottom",
            fontsize=19,
            fontweight="bold",
            color=INK,
        )
        for row, (model_key, _, _, _) in enumerate(MODEL_ORDER):
            y = top_y if row == 0 else bottom_y
            ax = fig.add_axes([x, y, width, height])
            panel = df[df["model_key"].eq(model_key) & df["cent_plot"].eq(cent)]
            if panel.empty:
                raise RuntimeError(f"Missing panel for {model_key}, {cent}")
            plot_panel(ax, panel, show_ylabel=col == 0, show_xlabel=row == 1)
            axes.append(ax)

    add_footer(fig, df)

    fig.savefig(OUT_PNG, facecolor="white")
    plt.close(fig)
    print(OUT_PNG)


if __name__ == "__main__":
    main()
