#!/usr/bin/env python3
"""Make slide-ready score-separation panels from compact route validation CSVs."""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager  # noqa: E402


ROOT = Path("dataOutput/auauMLDiagnosticRuns/ppg12_weighted_route_slide_recovery_20260522")
OUT_DIR = ROOT / "slideReady"

TIMES_FILES = [
    Path("/System/Library/Fonts/Supplemental/Times New Roman.ttf"),
    Path("/System/Library/Fonts/Supplemental/Times New Roman Bold.ttf"),
    Path("/System/Library/Fonts/Supplemental/Times New Roman Italic.ttf"),
    Path("/System/Library/Fonts/Supplemental/Times New Roman Bold Italic.ttf"),
]
TIMES = "Times New Roman"

SIGNAL = "#1F77B4"
BACKGROUND = "#D95F02"
INK = "#111827"
MUTED = "#334155"
GRID = "#E5E7EB"

CENT_ORDER = ["0_20", "20_50", "50_80"]


def setup_style() -> None:
    for path in TIMES_FILES:
        if path.exists():
            font_manager.fontManager.addfont(str(path))
    plt.rcParams.update(
        {
            "font.family": TIMES,
            "font.serif": [TIMES],
            "mathtext.fontset": "custom",
            "mathtext.rm": TIMES,
            "mathtext.it": f"{TIMES}:italic",
            "mathtext.bf": f"{TIMES}:bold",
            "axes.unicode_minus": False,
        }
    )


def read_rows(path: Path) -> list[dict]:
    with path.open() as handle:
        return list(csv.DictReader(handle))


def group_rows(rows: list[dict]) -> dict[tuple[str, str], list[dict]]:
    groups: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for row in rows:
        groups[(row["product"], row["cent_bin"])].append(row)
    for key in groups:
        groups[key].sort(key=lambda row: float(row["bin_lo"]))
    return groups


def read_summary(path: Path) -> dict[tuple[str, str], dict]:
    with path.open() as handle:
        return {(row["product"], row["cent_bin"]): row for row in csv.DictReader(handle)}


def step(ax, rows: list[dict], key: str, color: str, label: str) -> None:
    xs = [float(r["bin_lo"]) for r in rows] + [float(rows[-1]["bin_hi"])]
    ys = [float(r[key]) for r in rows]
    ys = ys + [ys[-1]]
    ax.step(xs, ys, where="post", color=color, linewidth=1.9, label=label)


def style_axis(ax, *, show_xlabel: bool, show_ylabel: bool) -> None:
    ax.set_xlim(0, 1)
    ax.grid(True, color=GRID, linewidth=0.55)
    ax.tick_params(direction="in", top=True, right=True, labelsize=8.8, pad=2)
    for spine in ax.spines.values():
        spine.set_linewidth(0.9)
        spine.set_color("#1F2937")
    if show_xlabel:
        ax.set_xlabel("BDT score", fontsize=12.0, labelpad=4)
    else:
        ax.set_xticklabels([])
    if show_ylabel:
        ax.set_ylabel("Unit area", fontsize=11.3, labelpad=5)


def plot_route_slide(
    *,
    csv_path: Path,
    summary_path: Path,
    products: list[str],
    title: str,
    subtitle: str,
    note: str,
    out_png: Path,
    row_colors: list[str],
) -> None:
    rows = read_rows(csv_path)
    groups = group_rows(rows)
    summary = read_summary(summary_path)
    label_by_product = {row["product"]: row["model_label"] for row in rows}
    route_by_product = {row["product"]: row["route_label"] for row in rows}
    cent_label = {row["cent_bin"]: row["cent_label"] for row in rows}

    ymax = 0.0
    for product in products:
        for cent in CENT_ORDER:
            panel = groups[(product, cent)]
            ymax = max(
                ymax,
                max(float(r["signal_density"]) for r in panel),
                max(float(r["background_density"]) for r in panel),
            )
    ymax *= 1.18

    fig, axes = plt.subplots(len(products), 3, figsize=(16, 9), sharex=True, sharey=True)
    if len(products) == 1:
        axes = [axes]

    fig.patch.set_facecolor("white")
    fig.text(0.035, 0.960, title, ha="left", va="top", fontsize=25.5, fontweight="bold", color=INK)
    fig.text(0.035, 0.912, subtitle, ha="left", va="top", fontsize=14.0, color=MUTED)
    fig.text(0.035, 0.875, note, ha="left", va="top", fontsize=12.4, color=MUTED)

    fig.text(0.775, 0.956, "Signal", ha="left", va="center", fontsize=13.0, color=SIGNAL, fontweight="bold")
    fig.text(0.840, 0.956, "Background", ha="left", va="center", fontsize=13.0, color=BACKGROUND, fontweight="bold")
    fig.add_artist(plt.Line2D([0.735, 0.770], [0.956, 0.956], transform=fig.transFigure, color=SIGNAL, lw=2.2))
    fig.add_artist(plt.Line2D([0.800, 0.835], [0.956, 0.956], transform=fig.transFigure, color=BACKGROUND, lw=2.2))

    if len(products) == 3:
        row_y_centers = [0.705, 0.455, 0.205]
        strip_height = 0.170
    elif len(products) == 2:
        row_y_centers = [0.640, 0.300]
        strip_height = 0.220
    else:
        row_y_centers = [0.500]
        strip_height = 0.240

    for i, product in enumerate(products):
        row_y0 = row_y_centers[i]
        fig.add_artist(
            plt.Rectangle(
                (0.026, max(0.08, row_y0 - strip_height / 2)),
                0.010,
                strip_height,
                transform=fig.transFigure,
                facecolor=row_colors[i],
                edgecolor=row_colors[i],
                linewidth=0,
            )
        )
        fig.text(
            0.045,
            row_y0 + 0.030,
            label_by_product[product],
            ha="left",
            va="center",
            fontsize=14.5,
            fontweight="bold",
            color=INK,
        )
        fig.text(
            0.045,
            row_y0 - 0.020,
            route_by_product[product],
            ha="left",
            va="center",
            fontsize=12.2,
            color=MUTED,
        )

        for j, cent in enumerate(CENT_ORDER):
            ax = axes[i][j]
            panel = groups[(product, cent)]
            step(ax, panel, "signal_density", SIGNAL, "Signal")
            step(ax, panel, "background_density", BACKGROUND, "Background")
            ax.set_ylim(0, ymax)
            style_axis(ax, show_xlabel=(i == len(products) - 1), show_ylabel=(j == 0))
            if i == 0:
                ax.set_title(cent_label[cent], fontsize=14.2, fontweight="bold", pad=8, color=INK)
            s = summary[(product, cent)]
            auc = float(s["auc"])
            sig_n = int(float(s["signal_entries"]))
            bkg_n = int(float(s["background_entries"]))
            ax.text(
                0.045,
                0.900,
                f"AUC {auc:.3f}",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=11.0,
                fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="#CBD5E1", alpha=0.94),
            )
            ax.text(
                0.045,
                0.785,
                f"S {sig_n/1e6:.2f}M  B {bkg_n/1e3:.0f}k",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=9.4,
                color=MUTED,
                bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="#E5E7EB", alpha=0.88),
            )

    left = 0.270 if len(products) == 2 else 0.225
    fig.subplots_adjust(left=left, right=0.985, top=0.815, bottom=0.090, wspace=0.150, hspace=0.210)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=220)
    plt.close(fig)


def main() -> None:
    setup_style()
    plot_route_slide(
        csv_path=ROOT / "global32_route_comparison/global32_route_comparison_score_histograms_compact.csv",
        summary_path=ROOT / "global32_route_comparison/global32_route_comparison_summary.csv",
        products=[
            "globalEtCent1535_bdt_noIso",
            "globalEtCent1535_bdt_noIso_ptCent3",
            "globalEtCent1535_bdt_noIso_ptCent7",
        ],
        title="Routing improves the same 32-input BDT",
        subtitle=r"Held-out embedded validation, $15<E_{T}<35$ GeV; each row uses the same input feature family and PPG12-style training weights.",
        note="The 8x7 routed BDT gives the cleanest separation in every centrality slice; the trend is visible without changing inputs.",
        out_png=OUT_DIR / "slide27_global32_vs_routed_8x3_8x7_fullstat.png",
        row_colors=["#F472B6", "#FB923C", "#22C55E"],
    )
    plot_route_slide(
        csv_path=ROOT / "centinput_global_vs_routed/centinput_global_vs_routed_score_histograms_compact.csv",
        summary_path=ROOT / "centinput_global_vs_routed/centinput_global_vs_routed_summary.csv",
        products=[
            "centInput_pt1535",
            "ptFine_cent7",
        ],
        title="Routing the Stack-input BDT sharpens separation",
        subtitle=r"Same baseV3E + $w_{\eta,3\times3}/w_{\phi,3\times3}$ feature family; top row is the global centInput BDT, bottom row is the new 8x7 routed analogue.",
        note="This isolates the routing change for the BDT that feeds the stack, using full-stat held-out embedded validation.",
        out_png=OUT_DIR / "slide27_centinput_global_vs_ptfine_cent7_fullstat.png",
        row_colors=["#F472B6", "#FB923C"],
    )


if __name__ == "__main__":
    main()
