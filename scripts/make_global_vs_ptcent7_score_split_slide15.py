#!/usr/bin/env python3
"""Build a full-slide PNG comparing global and routed 32-input noIso BDTs."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import font_manager


ROOT = Path(
    "dataOutput/auauMLDiagnosticRuns/"
    "global_etcent_inclusive3_sixpack_20260516_135439"
)
GLOBAL_VALIDATION = ROOT / "validation/bdt_finished_only_20260516_180817"
ROUTED_VALIDATION = ROOT / "validation/bdt_binned_sidecars_fullstat_20260517_2152"
OUT_DIR = ROOT / "slideReady/binned_bdt_comparison/cent7"
OUT_PNG = OUT_DIR / "global_noiso_bdt_vs_ptCent7_routed_score_split_slide15_clean.png"
OUT_CSV = OUT_DIR / "global_noiso_bdt_vs_ptCent7_routed_score_split_slide15_clean.csv"

GLOBAL_PRODUCT = "globalEtCent1535_bdt_noIso"
ROUTED_PRODUCT = "globalEtCent1535_bdt_noIso_ptCent7"
CENT_BINS = [("0_20", "0-20%"), ("20_50", "20-50%"), ("50_80", "50-80%")]

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
    (GLOBAL_PRODUCT, "Global BDT\nexpanded inputs", "1 BDT", "32 inputs"),
    (ROUTED_PRODUCT, r"Routed BDT" "\n" r"$E_T \times$ centrality", "56 BDTs", "32 inputs"),
]


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


def load_table() -> pd.DataFrame:
    global_hist = json.loads((GLOBAL_VALIDATION / "score_histograms.json").read_text())
    global_metrics = json.loads((GLOBAL_VALIDATION / "validation_metrics.json").read_text())
    edges = [float(x) for x in global_hist["bin_edges"]]
    rows = []

    global_payload = global_hist["products"][GLOBAL_PRODUCT]["by_centrality"]
    global_auc = global_metrics["products"][GLOBAL_PRODUCT]["auc_by_centrality"]
    for cent_key, cent_label in CENT_BINS:
        payload = global_payload[cent_key]
        for idx, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
            rows.append(
                {
                    "model_key": GLOBAL_PRODUCT,
                    "row_model": "Global 32-input BDT",
                    "cent_bin": cent_key,
                    "cent_plot": cent_label,
                    "bin_lo": lo,
                    "bin_hi": hi,
                    "score_center": 0.5 * (lo + hi),
                    "signal_density": float(payload["signal"]["density"][idx]),
                    "background_density": float(payload["background"]["density"][idx]),
                    "signal_entries": int(payload["signal"]["entries"]),
                    "background_entries": int(payload["background"]["entries"]),
                    "auc": float(global_auc[cent_key]),
                }
            )

    routed = pd.read_csv(ROUTED_VALIDATION / "coarse_centrality_score_histograms_binned_noIso.csv")
    routed = routed[routed["product"].eq(ROUTED_PRODUCT)].copy()
    if routed.empty:
        raise RuntimeError(f"No rows found for {ROUTED_PRODUCT}")
    for _, row in routed.iterrows():
        rows.append(
            {
                "model_key": ROUTED_PRODUCT,
                "row_model": "Routed 32-input BDT",
                "cent_bin": str(row["cent_bin"]),
                "cent_plot": str(row["cent_label"]),
                "bin_lo": float(row["bin_lo"]),
                "bin_hi": float(row["bin_hi"]),
                "score_center": 0.5 * (float(row["bin_lo"]) + float(row["bin_hi"])),
                "signal_density": float(row["signal_density"]),
                "background_density": float(row["background_density"]),
                "signal_entries": int(row["signal_entries"]),
                "background_entries": int(row["background_entries"]),
                "auc": float(row["auc"]),
            }
        )

    df = pd.DataFrame(rows)
    df["model_key"] = pd.Categorical(df["model_key"], [row[0] for row in MODEL_ORDER], ordered=True)
    df["cent_plot"] = pd.Categorical(df["cent_plot"], [row[1] for row in CENT_BINS], ordered=True)
    return df.sort_values(["model_key", "cent_plot", "score_center"]).reset_index(drop=True)


def add_header(fig: plt.Figure) -> None:
    fig.text(
        0.025,
        0.962,
        r"Routing alone improves the same 32-input BDT",
        ha="left",
        va="top",
        fontsize=30.0,
        fontweight="bold",
        color=INK,
    )
    fig.text(
        0.027,
        0.902,
        r"Only the model routing changes: one global no-isolation BDT versus separate $E_T \times$ centrality BDTs using the same 32-input family.",
        ha="left",
        va="top",
        fontsize=15.2,
        color=MUTED,
    )
    fig.text(
        0.027,
        0.862,
        r"Area-normalized BDT score density for Photon12+20 signal vs Jet12+20+30 inclusive background, $15<E_T<35$ GeV; blue = signal, orange = background.",
        ha="left",
        va="top",
        fontsize=15.2,
        color=MUTED,
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
    )
    ax.step(
        panel["score_center"],
        panel["background_density"],
        where="mid",
        color=BACKGROUND,
        linewidth=2.4,
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
        df.groupby(["model_key", "cent_plot"], observed=True)["auc"]
        .first()
        .unstack("cent_plot")
        .reindex(index=[row[0] for row in MODEL_ORDER], columns=[row[1] for row in CENT_BINS])
    )
    gains = 100.0 * (aucs.iloc[1] / aucs.iloc[0] - 1.0)
    gain_text = ", ".join(f"{cent}: +{gains[cent]:.1f}% AUC" for _, cent in CENT_BINS)

    fig.text(
        0.027,
        0.065,
        f"Routing-only gain by centrality: {gain_text}.",
        ha="left",
        va="bottom",
        fontsize=15.7,
        fontweight="bold",
        color=INK,
    )
    fig.text(
        0.027,
        0.034,
        r"Interpretation: the same input family ranks photons better when the decision boundary is allowed to adapt across $E_T$ and Au+Au occupancy.",
        ha="left",
        va="bottom",
        fontsize=13.8,
        color=LIGHT_MUTED,
    )


def plot(df: pd.DataFrame) -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_CSV, index=False)

    fig = plt.figure(figsize=(16.0, 9.0), dpi=200)
    fig.patch.set_facecolor("white")
    add_header(fig)

    left = 0.225
    width = 0.225
    gap = 0.030
    top_y = 0.515
    bottom_y = 0.205
    height = 0.235

    for y, fill in ((top_y, ROW_SHADE), (bottom_y, "white")):
        fig.add_artist(
            plt.Rectangle(
                (0.018, y - 0.035),
                0.957,
                height + 0.088,
                transform=fig.transFigure,
                facecolor=fill,
                edgecolor="none",
                zorder=-1,
            )
        )

    add_row_label(fig, top_y + 0.118, MODEL_ORDER[0][1], MODEL_ORDER[0][2], MODEL_ORDER[0][3], "#94A3B8")
    add_row_label(fig, bottom_y + 0.118, MODEL_ORDER[1][1], MODEL_ORDER[1][2], MODEL_ORDER[1][3], "#2563EB")

    for col, (_, cent_label) in enumerate(CENT_BINS):
        x = left + col * (width + gap)
        fig.text(
            x + width / 2,
            top_y + height + 0.038,
            cent_label,
            ha="center",
            va="bottom",
            fontsize=19,
            fontweight="bold",
            color=INK,
        )
        for row, (model_key, _, _, _) in enumerate(MODEL_ORDER):
            y = top_y if row == 0 else bottom_y
            ax = fig.add_axes([x, y, width, height])
            panel = df[df["model_key"].astype(str).eq(model_key) & df["cent_plot"].astype(str).eq(cent_label)]
            if panel.empty:
                raise RuntimeError(f"Missing panel for {model_key}, {cent_label}")
            plot_panel(ax, panel, show_ylabel=col == 0, show_xlabel=row == 1)

    add_footer(fig, df)
    fig.savefig(OUT_PNG, facecolor="white")
    plt.close(fig)
    return OUT_PNG


def main() -> None:
    setup_matplotlib()
    out_png = plot(load_table())
    print(out_png)
    print(OUT_CSV)


if __name__ == "__main__":
    main()
