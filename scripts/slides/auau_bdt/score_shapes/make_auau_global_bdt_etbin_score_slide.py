#!/usr/bin/env python3
"""Make a polished full-slide 4x2 E_T-bin BDT score-separation PNG."""

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
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


REPO = Path(__file__).resolve().parents[1]
DEFAULT_INDIR = (
    REPO
    / "dataOutput/auauMLDiagnosticRuns/ppg12_weighted_basev3e_stack_20260521_231600/"
    "slideReady/global_bdt_etbin_4x2_fullstat"
)
DEFAULT_CSV = DEFAULT_INDIR / "centInput_pt1535_score_separation_etbins_fullstat_test.csv"
DEFAULT_META = DEFAULT_INDIR / "centInput_pt1535_score_separation_etbins_fullstat_test_metadata.json"
DEFAULT_OUT = DEFAULT_INDIR / "slide28_global_centInput_pt1535_score_separation_by_et_4x2.png"

SIGNAL = "#CC334E"
BACKGROUND = "#0072B2"
INK = "#111827"
MUTED = "#4B5563"
GRID = "#E5E7EB"
BLUE_TINT = "#EAF2FF"
BLUE_EDGE = "#BFD7FF"
PINK_TINT = "#FDE7F1"
PINK_EDGE = "#F4B9D4"
YELLOW_TINT = "#FFF4C7"
YELLOW_EDGE = "#F4D35E"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_META)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    return parser.parse_args()


def load_histograms(path: Path):
    rows_by_et: dict[str, dict[str, list[dict[str, float]]]] = defaultdict(lambda: defaultdict(list))
    with path.open() as handle:
        for row in csv.DictReader(handle):
            item = {
                "et_lo": float(row["et_lo"]),
                "et_hi": float(row["et_hi"]),
                "bin_low": float(row["bin_low"]),
                "bin_high": float(row["bin_high"]),
                "density": float(row["density"]),
                "count": int(row["count"]),
                "entries": int(row["entries"]),
                "auc": float(row["auc"]),
                "signal_entries": int(row["signal_entries"]),
                "background_entries": int(row["background_entries"]),
            }
            rows_by_et[row["et_label"]][row["class"]].append(item)
    for cls_map in rows_by_et.values():
        for values in cls_map.values():
            values.sort(key=lambda item: item["bin_low"])
    return rows_by_et


def plot_step(ax, rows, color, label, *, lw=2.2):
    x_vals: list[float] = []
    y_vals: list[float] = []
    for idx, row in enumerate(rows):
        lo = row["bin_low"]
        hi = row["bin_high"]
        val = row["density"]
        if idx == 0:
            x_vals.extend([lo, lo])
            y_vals.extend([0.0, val])
        else:
            x_vals.extend([lo, lo])
            y_vals.extend([y_vals[-1], val])
        x_vals.append(hi)
        y_vals.append(val)
    if rows:
        x_vals.append(rows[-1]["bin_high"])
        y_vals.append(0.0)
    ax.plot(x_vals, y_vals, color=color, lw=lw, label=label, solid_joinstyle="miter")


def add_card(fig, xywh, title, body, *, face, edge, title_color=INK, body_color=MUTED):
    ax = fig.add_axes(xywh)
    ax.axis("off")
    ax.add_patch(
        plt.Rectangle(
            (0, 0),
            1,
            1,
            transform=ax.transAxes,
            facecolor=face,
            edgecolor=edge,
            linewidth=1.0,
        )
    )
    ax.text(0.035, 0.69, title, fontsize=13.5, fontweight="bold", color=title_color, ha="left", va="center")
    ax.text(0.035, 0.28, body, fontsize=11.6, color=body_color, ha="left", va="center", linespacing=1.15)
    return ax


def draw_slide(csv_path: Path, metadata_path: Path, out_path: Path) -> None:
    rows_by_et = load_histograms(csv_path)
    metadata = json.loads(metadata_path.read_text())
    et_labels = [f"{int(lo)}-{int(hi)}" for lo, hi in metadata["et_bins"]]
    slide_title = metadata.get(
        "slide_title",
        r"Global BDT score separation weakens at mid-$E_T$, then stabilizes",
    )
    slide_subtitle = metadata.get(
        "slide_subtitle",
        r"One $15 < E_T < 35$ GeV Au+Au BDT, evaluated on full-stat held-out validation; each panel is unit-area signal vs background.",
    )
    internal_subtitle = metadata.get("internal_subtitle", "PYTHIA8 Au+Au embedded validation")
    model_body = metadata.get("model_label", r"global $\bf{centInput\_pt1535}$ BDT")
    inputs_title = metadata.get("inputs_title", "Inputs")
    inputs_body = metadata.get("inputs_label", r"baseV3E + centrality + $w_{\eta,33}/w_{\phi,33}$")
    training_body = metadata.get("training_label", r"PPG12-style $E_T$ and $\eta$ weights")
    training_title = metadata.get("training_title", "Training prior control")
    readout_note = metadata.get(
        "readout_note",
        rf"Background shifts toward signal-like scores around 21-25 GeV; full-stat held-out test has {metadata.get('rows_loaded', 0):,} scored rows, $15 \leq E_T < 35$ GeV, 0-80% centrality.",
    )

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "axes.linewidth": 1.0,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "mathtext.fontset": "dejavuserif",
        }
    )

    fig = plt.figure(figsize=(16, 9), dpi=220, facecolor="white")

    fig.text(
        0.045,
        0.965,
        slide_title,
        fontsize=26.2,
        fontweight="bold",
        color=INK,
        ha="left",
        va="top",
    )
    fig.text(
        0.046,
        0.915,
        slide_subtitle,
        fontsize=13.2,
        color=MUTED,
        ha="left",
        va="top",
    )
    fig.text(0.794, 0.955, r"$\bf{\it{sPHENIX}}$ Internal", fontsize=13.6, ha="left", va="top", color=INK)
    fig.text(0.794, 0.925, internal_subtitle, fontsize=10.8, ha="left", va="top", color=INK)

    add_card(
        fig,
        [0.045, 0.828, 0.285, 0.060],
        "Model shown",
        model_body,
        face=PINK_TINT,
        edge=PINK_EDGE,
    )
    add_card(
        fig,
        [0.344, 0.828, 0.300, 0.060],
        inputs_title,
        inputs_body,
        face="#F8FAFC",
        edge="#CBD5E1",
    )
    add_card(
        fig,
        [0.658, 0.828, 0.215, 0.060],
        training_title,
        training_body,
        face=BLUE_TINT,
        edge=BLUE_EDGE,
        title_color="#1E3A8A",
        body_color="#1E3A8A",
    )

    legend_ax = fig.add_axes([0.885, 0.828, 0.110, 0.060])
    legend_ax.axis("off")
    legend_ax.set_xlim(0, 1)
    legend_ax.set_ylim(0, 1)
    legend_ax.plot([0.02, 0.28], [0.68, 0.68], color=SIGNAL, lw=4.0, solid_capstyle="butt")
    legend_ax.text(0.34, 0.68, "Signal", fontsize=11.0, color=INK, ha="left", va="center", clip_on=False)
    legend_ax.plot([0.02, 0.28], [0.30, 0.30], color=BACKGROUND, lw=4.0, solid_capstyle="butt")
    legend_ax.text(0.34, 0.30, "Background", fontsize=11.0, color=INK, ha="left", va="center", clip_on=False)

    axes = []
    left = 0.055
    bottom = 0.178
    width = 0.914
    height = 0.607
    cols = 4
    rows = 2
    x_gap = 0.027
    y_gap = 0.070
    panel_w = (width - x_gap * (cols - 1)) / cols
    panel_h = (height - y_gap * (rows - 1)) / rows
    for i in range(rows):
        for j in range(cols):
            axes.append(
                fig.add_axes(
                    [
                        left + j * (panel_w + x_gap),
                        bottom + (rows - 1 - i) * (panel_h + y_gap),
                        panel_w,
                        panel_h,
                    ]
                )
            )

    ymax = 0.0
    for label in et_labels:
        for cls in ("signal", "background"):
            ymax = max(ymax, max(row["density"] for row in rows_by_et[label][cls]))
    ymax = max(1.0, ymax * 1.14)

    for idx, label in enumerate(et_labels):
        ax = axes[idx]
        sig = rows_by_et[label]["signal"]
        bkg = rows_by_et[label]["background"]
        plot_step(ax, bkg, BACKGROUND, "Background")
        plot_step(ax, sig, SIGNAL, "Signal")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, ymax)
        ax.grid(True, color=GRID, lw=0.55)
        ax.set_axisbelow(True)
        lo, hi = label.split("-")
        ax.set_title(rf"${lo}<E_T<{hi}$ GeV", fontsize=14.0, fontweight="bold", pad=5)
        ax.tick_params(labelsize=10.4, length=4.5, width=0.9)
        ax.set_xticks(np.linspace(0, 1, 6))
        if idx % 4 == 0:
            ax.set_ylabel("Area-normalized density", fontsize=11.6)
        else:
            ax.set_yticklabels([])
        if idx >= 4:
            ax.set_xlabel("BDT score", fontsize=12.0, labelpad=3)
        else:
            ax.set_xticklabels([])

        auc = sig[0]["auc"]
        n_sig = sig[0]["signal_entries"]
        n_bkg = sig[0]["background_entries"]
        ax.text(
            0.046,
            0.925,
            f"AUC {auc:.3f}",
            transform=ax.transAxes,
            fontsize=12.0,
            fontweight="bold",
            ha="left",
            va="top",
            bbox=dict(
                boxstyle="round,pad=0.25",
                facecolor="white",
                edgecolor="#D1D5DB",
                linewidth=0.9,
                alpha=0.96,
            ),
        )
        ax.text(
            0.046,
            0.745,
            f"S {n_sig:,}\nB {n_bkg:,}",
            transform=ax.transAxes,
            fontsize=9.7,
            color=MUTED,
            ha="left",
            va="top",
            linespacing=1.18,
            bbox=dict(boxstyle="round,pad=0.18", facecolor="white", edgecolor="none", alpha=0.76),
        )

    aucs = [rows_by_et[label]["signal"][0]["auc"] for label in et_labels]
    best_idx = int(np.argmax(aucs))
    worst_idx = int(np.argmin(aucs))
    total_rows = metadata.get("rows_loaded", 0)

    band = fig.add_axes([0.045, 0.043, 0.918, 0.093])
    band.axis("off")
    band.add_patch(
        plt.Rectangle((0, 0), 1, 1, transform=band.transAxes, facecolor=YELLOW_TINT, edgecolor=YELLOW_EDGE, linewidth=1.0)
    )
    band.text(0.025, 0.62, "Readout", fontsize=14.2, fontweight="bold", color=INK, ha="left", va="center")
    band.text(
        0.130,
        0.62,
        rf"Best separation: {et_labels[best_idx]} GeV, AUC {aucs[best_idx]:.3f}; weakest bin: {et_labels[worst_idx]} GeV, AUC {aucs[worst_idx]:.3f}.",
        fontsize=13.0,
        color=INK,
        ha="left",
        va="center",
    )
    band.text(
        0.130,
        0.24,
        readout_note,
        fontsize=11.2,
        color=MUTED,
        ha="left",
        va="center",
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    draw_slide(args.csv, args.metadata, args.out)
    print(args.out)


if __name__ == "__main__":
    main()
