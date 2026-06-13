#!/usr/bin/env python3
"""Build a full-slide PNG explaining single-value BDT performance metrics."""

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

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch, Rectangle

from slides.common.slide_defaults import SLIDE_DPI, slide_figsize

OUTPUT_DPI = 300


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
CONTROL_DIR = (
    REPO
    / "dataOutput/auauTightBDTValidation/THE8_branchA_ladder_scorecache_fullstat_20260527"
    / "fixed_sample_controls"
)
EXACT_METRICS_CSV = CONTROL_DIR / "the8_branchA_fixed_sample_holdout_3x3_direct_exact_metrics_from_sdcc_20260611.csv"
TABLE_METRICS_CSV = (
    CONTROL_DIR
    / "the8_branchA_fixed_sample_holdout_3x3_direct_0to20_truthSignal_inclusiveJet_slide_v35_logloss_metric_table.metric_table_summary.csv"
)
OUTDIR = CONTROL_DIR / "single_value_metric_explainer_20260611"
OUTPNG = OUTDIR / "the8_single_value_metric_explainer_slide.png"
OUTSCRIPT = OUTDIR / "the8_single_value_metric_explainer_slide_script.md"
OUTMANIFEST = OUTDIR / "the8_single_value_metric_explainer_manifest.json"
OUTCSV = OUTDIR / "the8_single_value_metric_explainer_signed_changes.csv"

NAVY = "#101828"
INK = "#1D2939"
MUTED = "#475467"
BORDER = "#CBD5E1"
GRID = "#E4E9F2"
SOFT_BLUE = "#EFF6FF"
SOFT_GOLD = "#FFF7D6"
SOFT_RED = "#FEF3F2"
SOFT_PURPLE = "#F5F0FF"
SOFT_GREEN = "#ECFDF3"
SIGNAL = "#C5392F"
BLUE = "#1F77B4"
GOLD = "#C28400"
RED = "#C24135"
PURPLE = "#7651A6"
GREEN = "#228A4D"


def add_box(
    fig: plt.Figure,
    x: float,
    y: float,
    w: float,
    h: float,
    face: str,
    *,
    edge: str = BORDER,
    lw: float = 1.2,
    radius: float = 0.012,
    zorder: int = -2,
) -> None:
    fig.patches.append(
        FancyBboxPatch(
            (x, y),
            w,
            h,
            transform=fig.transFigure,
            boxstyle=f"round,pad=0.010,rounding_size={radius}",
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
    linespacing: float = 1.15,
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
        family="Times New Roman",
    )


def read_exact_metrics() -> dict[tuple[str, str], dict[str, float]]:
    rows: dict[tuple[str, str], dict[str, float]] = {}
    with EXACT_METRICS_CSV.open(newline="") as handle:
        for row in csv.DictReader(handle):
            if row["scope"] != "0_20":
                continue
            rows[(row["validation_sample"], row["training_sample"])] = {
                "auc": float(row["weighted_auc_source_classes"]),
                "logloss": float(row["weighted_logloss_source_classes"]),
                "brier": float(row["weighted_brier_source_classes"]),
                "gap": float(row["median_gap"]),
                "signal_entries": float(row["signal_entries"]),
                "background_entries": float(row["background_entries"]),
                "entries": float(row["entries"]),
            }
    return rows


def read_wp80_fake_rates() -> dict[tuple[str, str], float]:
    rows: dict[tuple[str, str], float] = {}
    with TABLE_METRICS_CSV.open(newline="") as handle:
        for row in csv.DictReader(handle):
            rows[(row["validation_sample"], row["training_sample"])] = 100.0 * float(
                row["wp80_background_fake_rate_binned"]
            )
    return rows


def combined_metrics() -> dict[tuple[str, str], dict[str, float]]:
    exact = read_exact_metrics()
    fake = read_wp80_fake_rates()
    for key, rate in fake.items():
        exact[key]["fake"] = rate
    return exact


def relative_better_change(left: float, right: float, *, higher_is_better: bool) -> float:
    if not np.isfinite(left) or left == 0:
        return float("nan")
    if higher_is_better:
        return 100.0 * (right - left) / abs(left)
    return 100.0 * (left - right) / abs(left)


def row_change(metrics: dict[tuple[str, str], dict[str, float]], validation: str) -> dict[str, float]:
    left = metrics[(validation, "Jet12+20")]
    right = metrics[(validation, "Jet12+20+30+40")]
    return {
        "AUC": relative_better_change(left["auc"], right["auc"], higher_is_better=True),
        "Logloss": relative_better_change(left["logloss"], right["logloss"], higher_is_better=False),
        "Median gap": relative_better_change(left["gap"], right["gap"], higher_is_better=True),
        "WP80 fake": relative_better_change(left["fake"], right["fake"], higher_is_better=False),
    }


def metric_definition_card(
    fig: plt.Figure,
    *,
    x: float,
    y: float,
    w: float,
    h: float,
    title: str,
    tag: str,
    accent: str,
    answers: str,
    read_as: str,
    watch: str,
) -> None:
    """Draw one open metric strip with bounded audience-facing text."""
    fig.patches.append(
        Rectangle(
            (x, y + h),
            w,
            0.0011,
            transform=fig.transFigure,
            facecolor="#E6EBF2",
            edgecolor="none",
            zorder=-1,
        )
    )
    fig.patches.append(
        Rectangle(
            (x + 0.004, y + 0.017),
            0.0065,
            h - 0.034,
            transform=fig.transFigure,
            facecolor=accent,
            edgecolor="none",
            zorder=-1,
        )
    )
    label_panel_w = 0.145
    title_size = 18.6 if len(title) <= 12 else 16.5
    add_text(
        fig,
        x + 0.020 + label_panel_w / 2.0,
        y + h - 0.034,
        title,
        size=title_size,
        weight="bold",
        color=NAVY,
        ha="center",
    )
    add_text(
        fig,
        x + 0.020 + label_panel_w / 2.0,
        y + h - 0.088,
        tag,
        size=11.5,
        weight="bold",
        color=accent,
        ha="center",
        linespacing=1.12,
    )

    fig.patches.append(
        Rectangle(
            (x + 0.185, y + 0.017),
            0.0012,
            h - 0.034,
            transform=fig.transFigure,
            facecolor=accent,
            edgecolor="none",
            alpha=0.30,
            zorder=-1,
        )
    )

    label_x = x + 0.203
    text_x = x + 0.285
    rows = [
        ("Question", answers, 12.8, "bold"),
        ("Use", read_as, 12.0, "normal"),
        ("Limit", watch, 12.0, "normal"),
    ]
    for row_i, (label, text, size, weight) in enumerate(rows):
        yy = y + h - (0.031 + 0.041 * row_i)
        add_text(fig, label_x, yy, label, size=10.8, weight="bold", color=accent, ha="left")
        add_text(fig, text_x, yy, text, size=size, weight=weight, color=INK, ha="left", linespacing=1.04)


def draw_overlay_plot(fig: plt.Figure, metrics: dict[tuple[str, str], dict[str, float]]) -> None:
    ax = fig.add_axes([0.625, 0.382, 0.318, 0.335])
    labels = ["AUC", "Logloss", "Median\ngap", "WP80\nfake"]
    x = np.arange(len(labels))
    top = row_change(metrics, "Jet12+20")
    bottom = row_change(metrics, "Jet12+20+30+40")
    y_top = np.array([top["AUC"], top["Logloss"], top["Median gap"], top["WP80 fake"]])
    y_bottom = np.array([bottom["AUC"], bottom["Logloss"], bottom["Median gap"], bottom["WP80 fake"]])

    ax.axhline(0.0, color="#667085", lw=1.0)
    value_label_box = {"facecolor": "white", "edgecolor": "none", "alpha": 0.84, "pad": 0.8}
    ax.plot(x, y_top, "-o", color=GREEN, lw=3.1, markersize=8.6, label="Jet12+20 validation")
    ax.plot(
        x,
        y_bottom,
        "-s",
        color=PURPLE,
        lw=3.1,
        markersize=8.2,
        label="Jet12+20+30+40 validation",
    )
    for xx, yy in zip(x, y_top):
        ax.text(
            xx,
            yy + (2.0 if yy >= 0 else -3.2),
            f"{yy:+.1f}%",
            ha="center",
            va="bottom" if yy >= 0 else "top",
            fontsize=11.2,
            color=GREEN,
            fontweight="bold",
            bbox=value_label_box,
        )
    for idx, (xx, yy) in enumerate(zip(x, y_bottom)):
        if idx == 2:
            ax.text(
                xx + 0.22,
                yy - 0.8,
                f"{yy:+.1f}%",
                ha="left",
                va="center",
                fontsize=11.5,
                color=PURPLE,
                fontweight="bold",
                bbox=value_label_box,
            )
            continue
        ax.text(
            xx + 0.050,
            yy - 3.2 if yy >= 0 else yy + 2.4,
            f"{yy:+.1f}%",
            ha="center",
            va="top" if yy >= 0 else "bottom",
            fontsize=11.2,
            color=PURPLE,
            fontweight="bold",
            bbox=value_label_box,
        )

    ax.set_ylabel("Gain relative to Jet12+20 training (%)", fontsize=11.7, color=INK, labelpad=4)
    ax.set_xticks(x, labels, fontsize=13.5)
    ax.set_ylim(-12, 70)
    ax.set_xlim(-0.35, len(labels) - 0.65)
    ax.grid(True, axis="y", color=GRID, lw=1.0)
    ax.set_axisbelow(True)
    ax.legend(
        loc="upper left",
        bbox_to_anchor=(0.00, 1.00),
        frameon=True,
        facecolor="white",
        edgecolor="none",
        framealpha=0.88,
        fontsize=11.7,
        handlelength=1.75,
        borderpad=0.25,
        labelspacing=0.25,
    )
    ax.tick_params(direction="in", top=False, right=False, labelsize=12.2, width=0.9)
    for spine in ax.spines.values():
        spine.set_color("#98A2B3")
        spine.set_linewidth(0.9)


def draw_under_plot_bullets(fig: plt.Figure, metrics: dict[tuple[str, str], dict[str, float]]) -> None:
    top_l = metrics[("Jet12+20", "Jet12+20")]
    top_r = metrics[("Jet12+20", "Jet12+20+30+40")]
    bot_l = metrics[("Jet12+20+30+40", "Jet12+20")]
    bot_r = metrics[("Jet12+20+30+40", "Jet12+20+30+40")]

    bullet_x = 0.608
    bullet_w = 0.335
    row_h = 0.054
    row_gap = 0.012
    add_text(fig, bullet_x, 0.323, "Readout from the metric comparison", size=17.5, weight="bold", color=NAVY)
    bullets = [
        (
            f"AUC: top-row rank changes only {top_l['auc']:.3f} to {top_r['auc']:.3f}.",
            BLUE,
        ),
        (
            "Logloss + median gap: score scale visibly cleans up.",
            PURPLE,
        ),
        (
            f"WP80 fake: hard-row leakage falls {bot_l['fake']:.1f}% to {bot_r['fake']:.1f}%.",
            RED,
        ),
    ]
    y_top = 0.226
    for text, color in bullets:
        row_y = y_top
        add_box(fig, bullet_x, row_y, bullet_w, row_h, "#F8FAFC", edge="#DCE4EE", lw=0.9, radius=0.006)
        bar_h = row_h - 0.022
        fig.patches.append(
            Rectangle(
                (bullet_x + 0.014, row_y + (row_h - bar_h) / 2.0),
                0.0052,
                bar_h,
                transform=fig.transFigure,
                facecolor=color,
                edgecolor="none",
                zorder=-1,
            )
        )
        add_text(
            fig,
            bullet_x + 0.030,
            row_y + row_h / 2.0,
            text,
            size=12.7,
            color=INK,
            linespacing=1.08,
            va="center",
        )
        y_top -= row_h + row_gap


def make_script(metrics: dict[tuple[str, str], dict[str, float]]) -> str:
    top_l = metrics[("Jet12+20", "Jet12+20")]
    top_r = metrics[("Jet12+20", "Jet12+20+30+40")]
    bot_l = metrics[("Jet12+20+30+40", "Jet12+20")]
    bot_r = metrics[("Jet12+20+30+40", "Jet12+20+30+40")]
    return f"""# Speaker Script: Single-Value BDT Metrics

This slide is the definition slide for the numbers printed on slide 12.

The key point is that the four numbers are not interchangeable. AUC asks a ranking question: if I choose one truth-isolated prompt photon candidate and one inclusive-jet background candidate from the same validation cell, how often does the signal candidate get the higher BDT score? That is why AUC can stay almost flat when the visible red and blue distributions pull apart but their pairwise ordering changes only slightly.

Logloss is different. It is the weighted binary cross-entropy on the candidate rows: minus the weighted average of y log p plus one minus y times log one minus p, where p is the BDT signal probability. These are exact SDCC matrix values for each validation-row and training-column cell, not a histogram read-off. Lower logloss means the model is assigning more decisive and more correct probabilities.

Median BDT gap is the direct score-axis split: median signal score minus median inclusive-background score in the same validation cell. That is why it tracks the visual separation on slide 12 better than AUC does.

WP80 fake rate is the analysis-cut readout. For each cell I set the BDT threshold so the binned signal efficiency is about 80 percent, then I report the weighted fraction of inclusive-jet background above that threshold. Lower is better.

The right plot puts all four metrics on one signed-improvement axis, comparing Jet12+20 training to Jet12+20+30+40 training. In the top row, AUC moves only {100.0 * (top_r['auc'] - top_l['auc']) / top_l['auc']:.1f} percent, but logloss improves by {100.0 * (top_l['logloss'] - top_r['logloss']) / top_l['logloss']:.1f} percent and the median gap grows by {100.0 * (top_r['gap'] - top_l['gap']) / top_l['gap']:.1f} percent. That is the clean explanation for why the histograms look more separated even though AUC barely changes.

In the bottom row, the validation sample itself includes the higher-reach Jet30 and Jet40 background support. There AUC and WP80 fake rate also respond more: AUC goes from {bot_l['auc']:.3f} to {bot_r['auc']:.3f}, and WP80 fake rate goes from {bot_l['fake']:.1f} percent to {bot_r['fake']:.1f} percent.
"""


def write_change_csv(metrics: dict[tuple[str, str], dict[str, float]]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for validation in ("Jet12+20", "Jet12+20+30+40"):
        left = metrics[(validation, "Jet12+20")]
        right = metrics[(validation, "Jet12+20+30+40")]
        for metric, key, higher in (
            ("AUC", "auc", True),
            ("Logloss", "logloss", False),
            ("Median gap", "gap", True),
            ("WP80 fake", "fake", False),
        ):
            rows.append(
                {
                    "validation_sample": validation,
                    "metric": metric,
                    "left_training": "Jet12+20",
                    "right_training": "Jet12+20+30+40",
                    "left_value": f"{left[key]:.12g}",
                    "right_value": f"{right[key]:.12g}",
                    "signed_change_percent_positive_better": f"{relative_better_change(left[key], right[key], higher_is_better=higher):.6g}",
                }
            )
    with OUTCSV.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return rows


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    metrics = combined_metrics()
    rows = write_change_csv(metrics)

    plt.rcParams.update(
        {
            "font.family": "Times New Roman",
            "axes.titleweight": "bold",
            "figure.facecolor": "#FFFFFF",
            "savefig.facecolor": "#FFFFFF",
        }
    )

    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI)
    fig.subplots_adjust(0, 0, 1, 1)

    add_text(fig, 0.055, 0.948, "The metrics answer different BDT-performance questions", size=25.0, weight="bold", color=NAVY)

    add_box(fig, 0.050, 0.085, 0.505, 0.770, "#FFFFFF", edge="#D5DCE8", lw=1.15, radius=0.014)
    add_text(fig, 0.070, 0.826, "Different BDT ranking metrics", size=18.2, weight="bold", color=NAVY)

    card_w = 0.470
    card_h = 0.147
    metric_definition_card(
        fig,
        x=0.066,
        y=0.625,
        w=card_w,
        h=card_h,
        title="AUC",
        tag="dataset-specific\nnon-additive rank",
        accent=BLUE,
        answers="Does signal outrank background?",
        read_as="Pair ordering in one validation row.",
        watch="Misses score distance or confidence.",
    )
    metric_definition_card(
        fig,
        x=0.066,
        y=0.450,
        w=card_w,
        h=card_h,
        title="Logloss",
        tag="candidate-level\nadditive loss",
        accent=PURPLE,
        answers="Are probabilities correct?",
        read_as="Confidence penalty; lower is better.",
        watch="Not a direct peak-split metric.",
    )
    metric_definition_card(
        fig,
        x=0.066,
        y=0.275,
        w=card_w,
        h=card_h,
        title="Median BDT gap",
        tag="dataset-local\nscore-shape summary",
        accent=GOLD,
        answers="How far apart are centers?",
        read_as="Visible red-blue median separation.",
        watch="Not a working-point leakage metric.",
    )
    metric_definition_card(
        fig,
        x=0.066,
        y=0.100,
        w=card_w,
        h=card_h,
        title="WP80 fake rate",
        tag="cut-specific\ndataset leakage",
        accent=RED,
        answers="How much bkg passes WP80?",
        read_as="Leakage at 80% signal efficiency.",
        watch="Ignores ranking away from the cut.",
    )

    add_box(fig, 0.595, 0.085, 0.360, 0.770, "#FFFFFF", edge="#D5DCE8", lw=1.15, radius=0.014)
    rhs_center_x = 0.595 + 0.360 / 2.0
    add_text(
        fig,
        rhs_center_x,
        0.833,
        "Gain from adding Jet30+40 to training",
        size=17.8,
        weight="bold",
        color=NAVY,
        ha="center",
    )
    add_text(
        fig,
        rhs_center_x,
        0.791,
        "Each point compares two trainings on the same validation row.\nBaseline = Jet12+20-trained BDT; plotted point = Jet12+20+30+40.\nPositive means better separation or less leakage.",
        size=11.2,
        color=MUTED,
        ha="center",
        linespacing=1.20,
    )
    draw_overlay_plot(fig, metrics)
    draw_under_plot_bullets(fig, metrics)

    fig.savefig(OUTPNG, dpi=OUTPUT_DPI)
    plt.close(fig)

    OUTSCRIPT.write_text(make_script(metrics))
    OUTMANIFEST.write_text(
        json.dumps(
            {
                "schema": "THE8_SINGLE_VALUE_METRIC_EXPLAINER_V2",
                "status": "READY",
                "slide": str(OUTPNG),
                "script": str(OUTSCRIPT),
                "signed_changes_csv": str(OUTCSV),
                "exact_metrics_csv": str(EXACT_METRICS_CSV),
                "wp80_metrics_csv": str(TABLE_METRICS_CSV),
                "source_definitions": {
                    "auc": "weighted roc_auc_score on source-class labels in the fixed validation cell",
                    "logloss": "weighted sklearn log_loss on source-class labels and BDT signal probability from exact SDCC matrices",
                    "median_bdt_gap": "median signal BDT score minus median inclusive-jet BDT score in the same cell",
                    "wp80_fake_rate": "binned inclusive-jet background fraction above the threshold closest to 80% binned signal efficiency",
                },
                "signed_change_rows": rows,
                "claim": "AUC is a rank-order health check; logloss and median BDT gap explain the visible score-axis separation; WP80 fake rate is the working-point leakage metric.",
                "google_slides_mutated": False,
            },
            indent=2,
        )
        + "\n"
    )
    print(OUTPNG)
    print(OUTSCRIPT)
    print(OUTMANIFEST)
    print(OUTCSV)


if __name__ == "__main__":
    main()
