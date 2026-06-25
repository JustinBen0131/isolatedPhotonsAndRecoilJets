#!/usr/bin/env python3
"""Render THE-37 fixed-validation sample-grid score-shape slides."""

from __future__ import annotations

import argparse
import json
import sys as _codex_sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Rectangle
from matplotlib.ticker import MaxNLocator

_CODEX_THIS_FILE = Path(__file__).resolve()
_CODEX_SCRIPTS_DIR = next((p for p in _CODEX_THIS_FILE.parents if p.name == "scripts"), _CODEX_THIS_FILE.parent)
if str(_CODEX_SCRIPTS_DIR) not in _codex_sys.path:
    _codex_sys.path.append(str(_CODEX_SCRIPTS_DIR))

from slides.common.slide_defaults import SLIDE_DPI, slide_figsize


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
OUT_DIR = REPO / "dataOutput/auauTightBDTValidation/THE37_sample_grid_slides_20260621"
DEFAULT_PAYLOAD = OUT_DIR / "the37_fixed_validation_score_shape_payload_20260621.json"
DEFAULT_SIGNAL_TRAINVAL_PAYLOAD = OUT_DIR / "the37_auau_fixed_inclusive_signal_trainval_matrix_payload_20260621.json"
DEFAULT_PP_SIGNAL_TRAINVAL_PAYLOAD = OUT_DIR / "the37_pp_fixed_inclusive_signal_trainval_matrix_payload_20260621.json"
DEFAULT_PP_INCLUSIVE_TRAINVAL_PAYLOAD = OUT_DIR / "the37_pp_fixed_signal_inclusive_trainval_matrix_payload_20260621.json"
TIMES_DIR = Path("/System/Library/Fonts/Supplemental")
TIMES_FONTS = [
    TIMES_DIR / "Times New Roman.ttf",
    TIMES_DIR / "Times New Roman Bold.ttf",
    TIMES_DIR / "Times New Roman Italic.ttf",
    TIMES_DIR / "Times New Roman Bold Italic.ttf",
]
FONT_FAMILY = "Times New Roman"

SIGNAL = "#d13b35"
BACKGROUND = "#1f78b4"
INK = "#111827"
MUTED = "#667085"
GRID = "#dde5f0"
BOX_EDGE = "#9aaec9"

PALETTE = [
    ("#4c9f68", "#eef8f0"),
    ("#d89b00", "#fff4cf"),
    ("#8a62c4", "#f1eafa"),
    ("#2e90fa", "#eef6ff"),
]


def configure_fonts() -> None:
    for font_path in TIMES_FONTS:
        if font_path.exists():
            font_manager.fontManager.addfont(str(font_path))
    plt.rcParams.update(
        {
            "font.family": FONT_FAMILY,
            "font.serif": [FONT_FAMILY],
            "mathtext.fontset": "custom",
            "mathtext.rm": FONT_FAMILY,
            "mathtext.it": f"{FONT_FAMILY}:italic",
            "mathtext.bf": f"{FONT_FAMILY}:bold",
        }
    )


def font_px(points: float) -> float:
    return points * SLIDE_DPI / 72.0


def bbox_px(fig: plt.Figure, artist) -> list[float]:
    fig.canvas.draw()
    bbox = artist.get_window_extent(renderer=fig.canvas.get_renderer())
    scale_x = fig.bbox.width / 2560.0
    scale_y = fig.bbox.height / 1440.0
    height = fig.bbox.height
    return [
        float(bbox.x0 / scale_x),
        float((height - bbox.y1) / scale_y),
        float(bbox.x1 / scale_x),
        float((height - bbox.y0) / scale_y),
    ]


def step_xy(edges: list[float], density: list[float]) -> tuple[list[float], list[float]]:
    x: list[float] = []
    y: list[float] = []
    for idx, value in enumerate(density):
        x.extend([edges[idx], edges[idx + 1]])
        y.extend([value, value])
    return x, y


def metric_rows(metrics: dict) -> list[tuple[str, str]]:
    return [
        ("AUC", f"{metrics['auc']:.3f}"),
        ("Logloss", f"{metrics['logloss']:.3f}"),
        ("Median gap", f"{metrics['median_gap']:.2f}"),
        ("WP80 incl.", f"{100.0 * metrics['wp80_background_fake_rate']:.1f}%"),
    ]


def metric_pair_rows(metrics: dict) -> list[tuple[str, str, str, str]]:
    return [
        ("AUC", f"{metrics['auc']:.3f}", "Gap", f"{metrics['median_gap']:.2f}"),
        ("Logloss", f"{metrics['logloss']:.3f}", "WP80", f"{100.0 * metrics['wp80_background_fake_rate']:.1f}%"),
    ]


def draw_metric_box(
    ax,
    metrics: dict,
    *,
    x: float = 0.02,
    y: float = 0.64,
    w: float = 0.72,
    h: float = 0.32,
    fontsize: float = 9.4,
    paired: bool = True,
) -> list:
    box = Rectangle(
        (x, y),
        w,
        h,
        transform=ax.transAxes,
        facecolor="white",
        edgecolor=BOX_EDGE,
        linewidth=0.8,
        zorder=5,
    )
    ax.add_patch(box)
    artists = []
    if paired:
        label_fontsize = max(7.8, fontsize - 1.8)
        for idx, (left_name, left_value, right_name, right_value) in enumerate(metric_pair_rows(metrics)):
            yy = y + h - 0.098 - idx * (h - 0.17)
            for text, xx, ha, bold in [
                (left_name, x + 0.026, "left", False),
                (left_value, x + 0.395, "right", True),
                (right_name, x + 0.535, "left", False),
                (right_value, x + w - 0.026, "right", True),
            ]:
                artists.append(
                    ax.text(
                        xx,
                        yy,
                        text,
                        transform=ax.transAxes,
                        ha=ha,
                        va="center",
                        fontsize=fontsize if bold else label_fontsize,
                        family=FONT_FAMILY,
                        fontweight="bold" if bold else "normal",
                        color=INK,
                        zorder=6,
                    )
                )
        return artists
    for idx, (name, value) in enumerate(metric_rows(metrics)):
        yy = y + h - 0.075 - idx * (h - 0.12) / 3.0
        artists.append(
            ax.text(
                x + 0.015,
                yy,
                name,
                transform=ax.transAxes,
                ha="left",
                va="center",
                fontsize=fontsize,
                family=FONT_FAMILY,
                color=INK,
                zorder=6,
            )
        )
        artists.append(
            ax.text(
                x + w - 0.015,
                yy,
                value,
                transform=ax.transAxes,
                ha="right",
                va="center",
                fontsize=fontsize,
                family=FONT_FAMILY,
                fontweight="bold",
                color=INK,
                zorder=6,
            )
        )
    return artists


def draw_axis(ax, cell: dict, *, ymax: float, show_y: bool, show_x: bool, compact: bool = False) -> list:
    edges = cell["bin_edges"]
    xs, ys = step_xy(edges, cell["signal"]["density"])
    xb, yb = step_xy(edges, cell["background"]["density"])
    ax.fill_between(xs, ys, step="pre", color=SIGNAL, alpha=0.10, linewidth=0)
    ax.fill_between(xb, yb, step="pre", color=BACKGROUND, alpha=0.10, linewidth=0)
    ax.plot(xs, ys, color=SIGNAL, linewidth=2.1)
    ax.plot(xb, yb, color=BACKGROUND, linewidth=2.1)
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, ymax)
    ax.grid(True, color=GRID, linewidth=0.7)
    ax.tick_params(axis="both", which="both", direction="in", top=True, right=True, labelsize=9.6, width=1.1)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=3, integer=True))
    if show_y:
        ax.set_ylabel("Unit-area\ndensity", fontsize=13, family=FONT_FAMILY, labelpad=0)
    else:
        ax.set_yticklabels([])
    if show_x:
        ax.set_xlabel("BDT score", fontsize=15, family=FONT_FAMILY, labelpad=0)
    else:
        ax.set_xticklabels([])
    for spine in ax.spines.values():
        spine.set_color("#1f2937")
        spine.set_linewidth(1.1)
    for item in ax.get_xticklabels() + ax.get_yticklabels():
        item.set_fontfamily(FONT_FAMILY)
    return draw_metric_box(ax, cell["metrics"], w=0.94 if compact else 0.82, fontsize=11.1 if compact else 10.8)


def draw_legend(fig: plt.Figure, y: float = 0.846) -> list:
    artists = []
    legend_x, legend_w, legend_h = 0.215, 0.74, 0.061
    fig.patches.append(
        Rectangle(
            (legend_x, y),
            legend_w,
            legend_h,
            transform=fig.transFigure,
            facecolor="white",
            edgecolor=BOX_EDGE,
            linewidth=0.9,
            joinstyle="round",
        )
    )
    fig.lines.append(
        plt.Line2D([legend_x + 0.05, legend_x + 0.095], [y + 0.031, y + 0.031], transform=fig.transFigure, color=SIGNAL, linewidth=3.0)
    )
    artists.append(
        fig.text(
            legend_x + 0.106,
            y + 0.031,
            "Signal MC",
            ha="left",
            va="center",
            fontsize=14,
            family=FONT_FAMILY,
            fontweight="bold",
            color=INK,
        )
    )
    fig.lines.append(
        plt.Line2D([legend_x + 0.41, legend_x + 0.455], [y + 0.031, y + 0.031], transform=fig.transFigure, color=BACKGROUND, linewidth=3.0)
    )
    artists.append(
        fig.text(
            legend_x + 0.465,
            y + 0.031,
            "Inclusive MC",
            ha="left",
            va="center",
            fontsize=14,
            family=FONT_FAMILY,
            fontweight="bold",
            color=INK,
        )
    )
    return artists


def label_text(label: str, max_parts: int = 3) -> str:
    parts = label.split("+")
    if len(parts) <= max_parts:
        return label
    return "+".join(parts[:max_parts]) + "\n+" + "+".join(parts[max_parts:])


def compact_training_label(label: str) -> str:
    parts = label.split("+")
    if len(parts) <= 2:
        return label
    return "+".join(parts[:2]) + "\n+" + "+".join(parts[2:])


def draw_row_label(fig: plt.Figure, label: str, y0: float, height: float, idx: int, *, prefix: str = "SIGNAL") -> list:
    accent, fill = PALETTE[idx % len(PALETTE)]
    x0, width = 0.015, 0.123
    fig.patches.append(
        Rectangle((x0, y0 + 0.005), width, height - 0.010, transform=fig.transFigure, facecolor=fill, edgecolor=accent, linewidth=0.9, joinstyle="round")
    )
    fig.patches.append(
        Rectangle((x0 + 0.018, y0 + 0.035), 0.006, height - 0.070, transform=fig.transFigure, facecolor=accent, edgecolor="none")
    )
    artists = [
        fig.text(
            x0 + 0.071,
            y0 + height * 0.67,
            f"TRAINING\n{prefix}",
            ha="center",
            va="center",
            fontsize=11.8,
            family=FONT_FAMILY,
            fontweight="bold",
            linespacing=0.86,
            color=INK,
        ),
        fig.text(
            x0 + 0.071,
            y0 + height * 0.30,
            label_text(label, max_parts=2),
            ha="center",
            va="center",
            fontsize=12.6,
            family=FONT_FAMILY,
            fontweight="bold",
            linespacing=0.82,
            color=INK,
        ),
    ]
    return artists


def draw_col_label(fig: plt.Figure, label: str, x0: float, width: float, idx: int, *, compact: bool = False) -> list:
    accent, fill = PALETTE[idx % len(PALETTE)]
    y, h = (0.748, 0.070) if compact else (0.760, 0.055)
    text = f"INCLUSIVE TRAINING\n{compact_training_label(label)}" if compact else f"INCLUSIVE TRAINING: {label}"
    fig.patches.append(
        Rectangle((x0, y), width, h, transform=fig.transFigure, facecolor=fill, edgecolor=accent, linewidth=0.9, joinstyle="round")
    )
    return [
        fig.text(
            x0 + width / 2,
            y + h / 2,
            text,
            ha="center",
            va="center",
            fontsize=11.5 if compact else 14.2,
            family=FONT_FAMILY,
            fontweight="bold",
            linespacing=0.82 if compact else 1.0,
            color=INK,
        )
    ]


def render_grid(payload: dict, system_key: str, out_path: Path) -> Path:
    config = payload[system_key]
    rows = config["row_order"]
    cols = config["col_order"]
    cells = {(c["row_key"], c["col_key"]): c for c in config["cells"]}
    n_rows, n_cols = len(rows), len(cols)
    compact = n_cols >= 4

    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI, facecolor="white")
    text_nodes = []
    title = fig.text(
        0.04,
        0.94,
        config["title"],
        fontsize=28 if compact else 30,
        family=FONT_FAMILY,
        fontweight="bold",
        color=INK,
    )
    text_nodes.append(("title", title, 28 if compact else 30, "title"))
    for idx, artist in enumerate(draw_legend(fig)):
        text_nodes.append((f"legend {idx}", artist, 14, "audience"))

    grid_left = 0.188
    gap_x = 0.018 if compact else 0.026
    plot_w = (0.965 - grid_left - gap_x * (n_cols - 1)) / n_cols
    row_top = 0.742
    gap_y = 0.046 if n_rows <= 2 else 0.042
    bottom_margin = 0.105 if n_rows <= 2 else 0.110
    plot_h = (row_top - bottom_margin - gap_y * (n_rows - 1)) / n_rows
    ymax = float(config.get("ymax") or max(max(c["signal"]["density"] + c["background"]["density"]) for c in config["cells"]) * 1.08)

    for col_idx, col in enumerate(cols):
        x0 = grid_left + col_idx * (plot_w + gap_x)
        for artist in draw_col_label(fig, col["label"], x0, plot_w, col_idx, compact=compact):
            text_nodes.append((f"column {col['key']}", artist, 11.5 if compact else 14.2, "audience"))

    for row_idx, row in enumerate(rows):
        y0 = row_top - row_idx * (plot_h + gap_y) - plot_h
        for idx, artist in enumerate(draw_row_label(fig, row["label"], y0, plot_h, row_idx, prefix=config.get("row_label_prefix", "SIGNAL"))):
            text_nodes.append((f"row {row['key']} {idx}", artist, 11.8 if idx == 0 else 12.6, "audience"))
        for col_idx, col in enumerate(cols):
            x0 = grid_left + col_idx * (plot_w + gap_x)
            ax = fig.add_axes([x0, y0, plot_w, plot_h])
            cell = cells[(row["key"], col["key"])]
            for idx, artist in enumerate(draw_axis(ax, cell, ymax=ymax, show_y=col_idx == 0, show_x=row_idx == n_rows - 1, compact=compact)):
                text_nodes.append((f"metric {row['key']} {col['key']} {idx}", artist, 9.0 if compact else 9.4, "plot_annotation"))

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=SLIDE_DPI)

    layout_path = out_path.with_suffix(".layout_nodes.json")
    nodes = []
    for name, artist, points, role in text_nodes:
        nodes.append(
            {
                "name": name,
                "kind": "text",
                "role": role,
                "text": artist.get_text(),
                "font_px": font_px(points),
                "bbox": bbox_px(fig, artist),
                "title_anchor": role == "title",
                "colon_style_exception": name.startswith("column"),
            }
        )
    layout_path.write_text(
        json.dumps(
            {
                "schema": "slide_layout_nodes_v1",
                "title_axis_x": 0.04 * 2560,
                "minimum_audience_font_px": font_px(11.5),
                "minimum_title_font_px": font_px(24),
                "minimum_plot_annotation_font_px": font_px(9.0),
                "nodes": nodes,
            },
            indent=2,
        )
        + "\n"
    )
    plt.close(fig)
    return out_path


def draw_signal_col_label(fig: plt.Figure, label: str, x0: float, width: float, idx: int) -> list:
    accent, fill = PALETTE[idx % len(PALETTE)]
    y, h = 0.690, 0.064
    fig.patches.append(
        Rectangle((x0, y), width, h, transform=fig.transFigure, facecolor=fill, edgecolor=accent, linewidth=0.9, joinstyle="round")
    )
    return [
        fig.text(
            x0 + width / 2,
            y + h / 2,
            f"SIGNAL TRAINING: {label}",
            ha="center",
            va="center",
            fontsize=15.0,
            family=FONT_FAMILY,
            fontweight="bold",
            color=INK,
        )
    ]


def draw_validation_signal_col_label(
    fig: plt.Figure,
    label: str,
    x0: float,
    width: float,
    idx: int,
    *,
    compact: bool = False,
) -> list:
    accent, fill = PALETTE[idx % len(PALETTE)]
    y, h = (0.690, 0.070) if compact else (0.700, 0.060)
    fig.patches.append(
        Rectangle((x0, y), width, h, transform=fig.transFigure, facecolor=fill, edgecolor=accent, linewidth=0.9, joinstyle="round")
    )
    text = f"VALIDATION SIGNAL\n{compact_training_label(label)}" if compact else f"VALIDATION SIGNAL: {label}"
    return [
        fig.text(
            x0 + width / 2,
            y + h / 2,
            text,
            ha="center",
            va="center",
            fontsize=12.4 if compact else 14.8,
            family=FONT_FAMILY,
            fontweight="bold",
            linespacing=0.86 if compact else 1.0,
            color=INK,
        )
    ]


def render_auau_full_bkg_signal_focus(payload: dict, out_path: Path) -> Path:
    config = payload["auau"]
    fixed_col_key = "12_20_30_40"
    rows = config["row_order"]
    cells = {(c["row_key"], c["col_key"]): c for c in config["cells"]}

    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI, facecolor="white")
    text_nodes = []
    title = fig.text(
        0.04,
        0.94,
        "Au+Au 0-20%: fixed full embedded inclusive MC vs signal-training sample",
        fontsize=27,
        family=FONT_FAMILY,
        fontweight="bold",
        color=INK,
    )
    text_nodes.append(("title", title, 27, "title"))
    for idx, artist in enumerate(draw_legend(fig, y=0.845)):
        text_nodes.append((f"legend {idx}", artist, 14, "audience"))

    fixed_x, fixed_y, fixed_w, fixed_h = 0.185, 0.772, 0.78, 0.060
    fig.patches.append(
        Rectangle(
            (fixed_x, fixed_y),
            fixed_w,
            fixed_h,
            transform=fig.transFigure,
            facecolor="#f1eafa",
            edgecolor="#8a62c4",
            linewidth=0.9,
            joinstyle="round",
        )
    )
    fixed = fig.text(
        fixed_x + fixed_w / 2,
        fixed_y + fixed_h / 2,
        "FIXED INCLUSIVE TRAINING: Jet12+20+30+40",
        ha="center",
        va="center",
        fontsize=16.0,
        family=FONT_FAMILY,
        fontweight="bold",
        color=INK,
    )
    text_nodes.append(("fixed background", fixed, 16.0, "audience"))

    grid_left = 0.185
    gap_x = 0.045
    plot_w = (0.965 - grid_left - gap_x) / 2.0
    plot_y = 0.195
    plot_h = 0.465
    ymax = float(config.get("ymax") or max(max(c["signal"]["density"] + c["background"]["density"]) for c in config["cells"]) * 1.08)

    for col_idx, row in enumerate(rows):
        x0 = grid_left + col_idx * (plot_w + gap_x)
        for artist in draw_signal_col_label(fig, row["label"], x0, plot_w, col_idx):
            text_nodes.append((f"signal column {row['key']}", artist, 15.0, "audience"))
        ax = fig.add_axes([x0, plot_y, plot_w, plot_h])
        cell = cells[(row["key"], fixed_col_key)]
        for idx, artist in enumerate(draw_axis(ax, cell, ymax=ymax, show_y=col_idx == 0, show_x=True, compact=False)):
            text_nodes.append((f"metric {row['key']} fixed_bkg {idx}", artist, 9.4, "plot_annotation"))

    left_metrics = cells[(rows[0]["key"], fixed_col_key)]["metrics"]
    right_metrics = cells[(rows[1]["key"], fixed_col_key)]["metrics"]
    readout = fig.text(
        0.185,
        0.095,
        (
            "Readout: adding Photon20 signal coverage raises AUC "
            f"{left_metrics['auc']:.3f} -> {right_metrics['auc']:.3f}, "
            f"median gap {left_metrics['median_gap']:.2f} -> {right_metrics['median_gap']:.2f},\n"
            f"and lowers WP80 inclusive pass {100.0 * left_metrics['wp80_background_fake_rate']:.1f}% -> "
            f"{100.0 * right_metrics['wp80_background_fake_rate']:.1f}%."
        ),
        ha="left",
        va="center",
        fontsize=15.0,
        family=FONT_FAMILY,
        linespacing=1.0,
        color=INK,
    )
    text_nodes.append(("readout", readout, 15.0, "audience"))

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=SLIDE_DPI)

    layout_path = out_path.with_suffix(".layout_nodes.json")
    nodes = []
    for name, artist, points, role in text_nodes:
        nodes.append(
            {
                "name": name,
                "kind": "text",
                "role": role,
                "text": artist.get_text(),
                "font_px": font_px(points),
                "bbox": bbox_px(fig, artist),
                "title_anchor": role == "title",
                "colon_style_exception": name in {"fixed background", "readout"} or name.startswith("signal column"),
            }
        )
    layout_path.write_text(
        json.dumps(
            {
                "schema": "slide_layout_nodes_v1",
                "title_axis_x": 0.04 * 2560,
                "minimum_audience_font_px": font_px(11.5),
                "minimum_title_font_px": font_px(24),
                "minimum_plot_annotation_font_px": font_px(9.0),
                "nodes": nodes,
            },
            indent=2,
        )
        + "\n"
    )
    plt.close(fig)
    return out_path


def render_auau_signal_trainval_matrix(payload: dict, out_path: Path) -> Path:
    config = payload["auau_signal_trainval"]
    rows = config["row_order"]
    cols = config["col_order"]
    cells = {(c["row_key"], c["col_key"]): c for c in config["cells"]}

    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI, facecolor="white")
    text_nodes = []
    title = fig.text(
        0.04,
        0.94,
        "Au+Au 0-20%: signal train/validation matrix with fixed inclusive MC",
        fontsize=28.0,
        family=FONT_FAMILY,
        fontweight="bold",
        color=INK,
    )
    text_nodes.append(("title", title, 28.0, "title"))
    for idx, artist in enumerate(draw_legend(fig, y=0.845)):
        text_nodes.append((f"legend {idx}", artist, 14, "audience"))

    fixed_x, fixed_y, fixed_w, fixed_h = 0.185, 0.775, 0.78, 0.055
    fig.patches.append(
        Rectangle(
            (fixed_x, fixed_y),
            fixed_w,
            fixed_h,
            transform=fig.transFigure,
            facecolor="#f1eafa",
            edgecolor="#8a62c4",
            linewidth=0.9,
            joinstyle="round",
        )
    )
    fixed = fig.text(
        fixed_x + fixed_w / 2,
        fixed_y + fixed_h / 2,
        config.get("fixed_inclusive_label", "FIXED INCLUSIVE MC: Jet12+20+30+40"),
        ha="center",
        va="center",
        fontsize=15.8,
        family=FONT_FAMILY,
        fontweight="bold",
        color=INK,
    )
    text_nodes.append(("fixed inclusive mc", fixed, 15.8, "audience"))

    grid_left = 0.185
    gap_x = 0.045
    plot_w = (0.965 - grid_left - gap_x) / 2.0
    row_top = 0.675
    gap_y = 0.058
    bottom_margin = 0.125
    plot_h = (row_top - bottom_margin - gap_y) / 2.0
    ymax = float(config.get("ymax") or max(max(c["signal"]["density"] + c["background"]["density"]) for c in config["cells"]) * 1.08)

    for col_idx, col in enumerate(cols):
        x0 = grid_left + col_idx * (plot_w + gap_x)
        for artist in draw_validation_signal_col_label(fig, col["label"], x0, plot_w, col_idx):
            text_nodes.append((f"validation signal column {col['key']}", artist, 14.8, "audience"))

    for row_idx, row in enumerate(rows):
        y0 = row_top - row_idx * (plot_h + gap_y) - plot_h
        for idx, artist in enumerate(draw_row_label(fig, row["label"], y0, plot_h, row_idx, prefix="SIGNAL")):
            text_nodes.append((f"train signal row {row['key']} {idx}", artist, 11.8 if idx == 0 else 12.6, "audience"))
        for col_idx, col in enumerate(cols):
            x0 = grid_left + col_idx * (plot_w + gap_x)
            ax = fig.add_axes([x0, y0, plot_w, plot_h])
            cell = cells[(row["key"], col["key"])]
            for idx, artist in enumerate(draw_axis(ax, cell, ymax=ymax, show_y=col_idx == 0, show_x=row_idx == len(rows) - 1, compact=False)):
                text_nodes.append((f"metric {row['key']} {col['key']} {idx}", artist, 9.4, "plot_annotation"))

    bottom_left = cells[(rows[0]["key"], cols[0]["key"])]["metrics"]
    bottom_right = cells[(rows[1]["key"], cols[1]["key"])]["metrics"]
    readout = fig.text(
        0.185,
        0.062,
        (
            "Readout: Photon12-only validation is the transfer stress test with inclusive MC fixed.\n"
            f"Matched Photon12+20 baseline: AUC {bottom_right['auc']:.3f}, logloss {bottom_right['logloss']:.3f}, "
            f"median gap {bottom_right['median_gap']:.2f}, WP80 inclusive pass {100.0 * bottom_right['wp80_background_fake_rate']:.1f}%."
        ),
        ha="left",
        va="center",
        fontsize=13.4,
        family=FONT_FAMILY,
        linespacing=0.95,
        color=INK,
    )
    text_nodes.append(("readout", readout, 13.4, "audience"))

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=SLIDE_DPI)

    layout_path = out_path.with_suffix(".layout_nodes.json")
    nodes = []
    for name, artist, points, role in text_nodes:
        nodes.append(
            {
                "name": name,
                "kind": "text",
                "role": role,
                "text": artist.get_text(),
                "font_px": font_px(points),
                "bbox": bbox_px(fig, artist),
                "title_anchor": role == "title",
                "colon_style_exception": name in {"fixed inclusive mc", "readout"} or name.startswith("validation signal column"),
            }
        )
    layout_path.write_text(
        json.dumps(
            {
                "schema": "slide_layout_nodes_v1",
                "title_axis_x": 0.04 * 2560,
                "minimum_audience_font_px": font_px(11.5),
                "minimum_title_font_px": font_px(24),
                "minimum_plot_annotation_font_px": font_px(9.0),
                "nodes": nodes,
            },
            indent=2,
        )
        + "\n"
    )
    plt.close(fig)
    return out_path


def render_pp_signal_trainval_matrix(payload: dict, out_path: Path) -> Path:
    config = payload["pp_signal_trainval"]
    rows = config["row_order"]
    cols = config["col_order"]
    cells = {(c["row_key"], c["col_key"]): c for c in config["cells"]}

    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI, facecolor="white")
    text_nodes = []
    title = fig.text(
        0.04,
        0.94,
        "pp 15-35 GeV: signal train/validation matrix with fixed inclusive MC",
        fontsize=27.0,
        family=FONT_FAMILY,
        fontweight="bold",
        color=INK,
    )
    text_nodes.append(("title", title, 27.0, "title"))
    for idx, artist in enumerate(draw_legend(fig, y=0.845)):
        text_nodes.append((f"legend {idx}", artist, 14, "audience"))

    fixed_x, fixed_y, fixed_w, fixed_h = 0.185, 0.775, 0.78, 0.055
    fig.patches.append(
        Rectangle(
            (fixed_x, fixed_y),
            fixed_w,
            fixed_h,
            transform=fig.transFigure,
            facecolor="#f1eafa",
            edgecolor="#8a62c4",
            linewidth=0.9,
            joinstyle="round",
        )
    )
    fixed = fig.text(
        fixed_x + fixed_w / 2,
        fixed_y + fixed_h / 2,
        config.get("fixed_inclusive_label", "FIXED INCLUSIVE MC: Jet8+12+20+30+40"),
        ha="center",
        va="center",
        fontsize=15.8,
        family=FONT_FAMILY,
        fontweight="bold",
        color=INK,
    )
    text_nodes.append(("fixed inclusive mc", fixed, 15.8, "audience"))

    grid_left = 0.185
    gap_x = 0.026
    plot_w = (0.965 - grid_left - gap_x * (len(cols) - 1)) / len(cols)
    row_top = 0.675
    gap_y = 0.037
    bottom_margin = 0.155
    plot_h = (row_top - bottom_margin - gap_y * (len(rows) - 1)) / len(rows)
    ymax = float(config.get("ymax") or max(max(c["signal"]["density"] + c["background"]["density"]) for c in config["cells"]) * 1.08)

    for col_idx, col in enumerate(cols):
        x0 = grid_left + col_idx * (plot_w + gap_x)
        for artist in draw_validation_signal_col_label(fig, col["label"], x0, plot_w, col_idx, compact=True):
            text_nodes.append((f"validation signal column {col['key']}", artist, 12.4, "audience"))

    for row_idx, row in enumerate(rows):
        y0 = row_top - row_idx * (plot_h + gap_y) - plot_h
        for idx, artist in enumerate(draw_row_label(fig, row["label"], y0, plot_h, row_idx, prefix="SIGNAL")):
            text_nodes.append((f"train signal row {row['key']} {idx}", artist, 11.8 if idx == 0 else 12.6, "audience"))
        for col_idx, col in enumerate(cols):
            x0 = grid_left + col_idx * (plot_w + gap_x)
            ax = fig.add_axes([x0, y0, plot_w, plot_h])
            cell = cells[(row["key"], col["key"])]
            for idx, artist in enumerate(
                draw_axis(ax, cell, ymax=ymax, show_y=col_idx == 0, show_x=row_idx == len(rows) - 1, compact=True)
            ):
                text_nodes.append((f"metric {row['key']} {col['key']} {idx}", artist, 9.0, "plot_annotation"))

    matched = cells[(rows[-1]["key"], cols[-1]["key"])]["metrics"]
    transfer = cells[(rows[0]["key"], cols[-1]["key"])]["metrics"]
    readout = fig.text(
        0.185,
        0.047,
        (
            "Readout: fixed blue = all Jet8+12+20+30+40 inclusive MC; columns change only signal validation.\n"
            f"Matched Photon5+10+20 endpoint: AUC {matched['auc']:.3f}, logloss {matched['logloss']:.3f}, "
            f"gap {matched['median_gap']:.2f}, WP80 incl. pass {100.0 * matched['wp80_background_fake_rate']:.1f}%; "
            f"Photon5-trained full-validation AUC {transfer['auc']:.3f}."
        ),
        ha="left",
        va="center",
        fontsize=11.6,
        family=FONT_FAMILY,
        linespacing=0.94,
        color=INK,
    )
    text_nodes.append(("readout", readout, 11.6, "audience"))

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=SLIDE_DPI)

    layout_path = out_path.with_suffix(".layout_nodes.json")
    nodes = []
    for name, artist, points, role in text_nodes:
        nodes.append(
            {
                "name": name,
                "kind": "text",
                "role": role,
                "text": artist.get_text(),
                "font_px": font_px(points),
                "bbox": bbox_px(fig, artist),
                "title_anchor": role == "title",
                "colon_style_exception": name in {"fixed inclusive mc", "readout"} or name.startswith("validation signal column"),
            }
        )
    layout_path.write_text(
        json.dumps(
            {
                "schema": "slide_layout_nodes_v1",
                "title_axis_x": 0.04 * 2560,
                "minimum_audience_font_px": font_px(11.5),
                "minimum_title_font_px": font_px(24),
                "minimum_plot_annotation_font_px": font_px(9.0),
                "nodes": nodes,
            },
            indent=2,
        )
        + "\n"
    )
    plt.close(fig)
    return out_path


def draw_validation_inclusive_row_label(fig: plt.Figure, label: str, y0: float, height: float, idx: int) -> list:
    accent, fill = PALETTE[idx % len(PALETTE)]
    x0, width = 0.015, 0.132
    fig.patches.append(
        Rectangle((x0, y0 + 0.004), width, height - 0.008, transform=fig.transFigure, facecolor=fill, edgecolor=accent, linewidth=0.9, joinstyle="round")
    )
    fig.patches.append(
        Rectangle((x0 + 0.019, y0 + 0.026), 0.006, height - 0.052, transform=fig.transFigure, facecolor=accent, edgecolor="none")
    )
    return [
        fig.text(
            x0 + 0.076,
            y0 + height * 0.69,
            "VALIDATION\nINCLUSIVE",
            ha="center",
            va="center",
            fontsize=10.2,
            family=FONT_FAMILY,
            fontweight="bold",
            linespacing=0.82,
            color=INK,
        ),
        fig.text(
            x0 + 0.076,
            y0 + height * 0.30,
            label_text(label, max_parts=2),
            ha="center",
            va="center",
            fontsize=10.8,
            family=FONT_FAMILY,
            fontweight="bold",
            linespacing=0.78,
            color=INK,
        ),
    ]


def draw_inclusive_training_col_label(fig: plt.Figure, label: str, x0: float, width: float, idx: int) -> list:
    accent, fill = PALETTE[idx % len(PALETTE)]
    y, h = 0.680, 0.062
    fig.patches.append(
        Rectangle((x0, y), width, h, transform=fig.transFigure, facecolor=fill, edgecolor=accent, linewidth=0.9, joinstyle="round")
    )
    return [
        fig.text(
            x0 + width / 2,
            y + h / 2,
            f"INCLUSIVE TRAINING\n{label_text(label, max_parts=2)}",
            ha="center",
            va="center",
            fontsize=10.4,
            family=FONT_FAMILY,
            fontweight="bold",
            linespacing=0.78,
            color=INK,
        )
    ]


def render_pp_inclusive_trainval_matrix(payload: dict, out_path: Path) -> Path:
    config = payload["pp_inclusive_trainval"]
    rows = config["row_order"]
    cols = config["col_order"]
    cells = {(c["row_key"], c["col_key"]): c for c in config["cells"]}

    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI, facecolor="white")
    text_nodes = []
    title = fig.text(
        0.04,
        0.94,
        "pp 15-35 GeV: fixed signal vs inclusive-MC train/validation sample",
        fontsize=25.6,
        family=FONT_FAMILY,
        fontweight="bold",
        color=INK,
    )
    text_nodes.append(("title", title, 25.6, "title"))
    for idx, artist in enumerate(draw_legend(fig, y=0.848)):
        text_nodes.append((f"legend {idx}", artist, 14, "audience"))

    fixed_x, fixed_y, fixed_w, fixed_h = 0.188, 0.758, 0.777, 0.046
    fig.patches.append(
        Rectangle(
            (fixed_x, fixed_y),
            fixed_w,
            fixed_h,
            transform=fig.transFigure,
            facecolor="#eef6ff",
            edgecolor="#2e90fa",
            linewidth=0.9,
            joinstyle="round",
        )
    )
    fixed = fig.text(
        fixed_x + fixed_w / 2,
        fixed_y + fixed_h / 2,
        config.get("fixed_signal_label", "FIXED SIGNAL MC: Photon5+10+20"),
        ha="center",
        va="center",
        fontsize=14.2,
        family=FONT_FAMILY,
        fontweight="bold",
        color=INK,
    )
    text_nodes.append(("fixed signal mc", fixed, 15.2, "audience"))

    grid_left = 0.188
    gap_x = 0.018
    plot_w = (0.965 - grid_left - gap_x * (len(cols) - 1)) / len(cols)
    row_top = 0.650
    gap_y = 0.030
    bottom_margin = 0.083
    plot_h = (row_top - bottom_margin - gap_y * (len(rows) - 1)) / len(rows)
    ymax = float(config.get("ymax") or max(max(c["signal"]["density"] + c["background"]["density"]) for c in config["cells"]) * 1.08)

    for col_idx, col in enumerate(cols):
        x0 = grid_left + col_idx * (plot_w + gap_x)
        for artist in draw_inclusive_training_col_label(fig, col["label"], x0, plot_w, col_idx):
            text_nodes.append((f"inclusive training column {col['key']}", artist, 10.4, "audience"))

    for row_idx, row in enumerate(rows):
        y0 = row_top - row_idx * (plot_h + gap_y) - plot_h
        for idx, artist in enumerate(draw_validation_inclusive_row_label(fig, row["label"], y0, plot_h, row_idx)):
            text_nodes.append((f"validation inclusive row {row['key']} {idx}", artist, 10.2 if idx == 0 else 10.8, "audience"))
        for col_idx, col in enumerate(cols):
            x0 = grid_left + col_idx * (plot_w + gap_x)
            ax = fig.add_axes([x0, y0, plot_w, plot_h])
            cell = cells[(row["key"], col["key"])]
            for idx, artist in enumerate(
                draw_axis(ax, cell, ymax=ymax, show_y=col_idx == 0, show_x=row_idx == len(rows) - 1, compact=True)
            ):
                text_nodes.append((f"metric {row['key']} {col['key']} {idx}", artist, 9.0, "plot_annotation"))

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=SLIDE_DPI)

    layout_path = out_path.with_suffix(".layout_nodes.json")
    nodes = []
    for name, artist, points, role in text_nodes:
        nodes.append(
            {
                "name": name,
                "kind": "text",
                "role": role,
                "text": artist.get_text(),
                "font_px": font_px(points),
                "bbox": bbox_px(fig, artist),
                "title_anchor": role == "title",
                "colon_style_exception": name in {"fixed signal mc"} or name.startswith("inclusive training column"),
            }
        )
    layout_path.write_text(
        json.dumps(
            {
                "schema": "slide_layout_nodes_v1",
                "title_axis_x": 0.04 * 2560,
                "minimum_audience_font_px": font_px(10.0),
                "minimum_title_font_px": font_px(24),
                "minimum_plot_annotation_font_px": font_px(9.0),
                "nodes": nodes,
            },
            indent=2,
        )
        + "\n"
    )
    plt.close(fig)
    return out_path


def cell_metrics(config: dict, row_key: str, col_key: str) -> dict:
    for cell in config["cells"]:
        if cell["row_key"] == row_key and cell["col_key"] == col_key:
            return cell["metrics"]
    raise KeyError(f"missing cell {row_key=} {col_key=}")


def fmt_metrics(metrics: dict) -> str:
    return (
        f"AUC {metrics['auc']:.3f}, logloss {metrics['logloss']:.3f}, "
        f"median gap {metrics['median_gap']:.2f}, "
        f"WP80 inclusive pass {100.0 * metrics['wp80_background_fake_rate']:.1f}%"
    )


def write_speaker_scripts(
    payload: dict,
    outdir: Path,
    signal_trainval_payload: dict | None = None,
    pp_signal_trainval_payload: dict | None = None,
    pp_inclusive_trainval_payload: dict | None = None,
) -> dict[str, Path]:
    auau_best = cell_metrics(payload["auau"], "20", "12_20_30_40")
    pp_full = cell_metrics(payload["pp"], "5_10_20", "8_12_20_30_40")
    pp_photon5 = cell_metrics(payload["pp"], "5", "8_12_20_30_40")
    auau_photon12 = cell_metrics(payload["auau"], "12", "12_20_30_40")
    scripts = {
        "auau": outdir / "the37_auau_fixed_validation_score_shape_grid_speaker_script.md",
        "pp": outdir / "the37_pp_fixed_validation_score_shape_grid_speaker_script.md",
        "auau_signal_focus": outdir / "the37_auau_full_inclusive_signal_variation_score_shape_speaker_script.md",
    }
    if signal_trainval_payload is not None:
        scripts["auau_signal_trainval"] = outdir / "the37_auau_fixed_inclusive_signal_trainval_matrix_score_shape_speaker_script.md"
    if pp_signal_trainval_payload is not None:
        scripts["pp_signal_trainval"] = outdir / "the37_pp_fixed_inclusive_signal_trainval_matrix_score_shape_speaker_script.md"
    if pp_inclusive_trainval_payload is not None:
        scripts["pp_inclusive_trainval"] = outdir / "the37_pp_fixed_signal_inclusive_trainval_matrix_score_shape_speaker_script.md"
    scripts["auau"].write_text(
        "\n".join(
            [
                "# AuAu sample-grid score-shape slide script",
                "",
                "This slide holds the AuAu validation sample fixed: 0-20% centrality, 15-35 GeV, Photon12+20 signal, and Jet12+20+30+40 embedded inclusive MC with no truth filter.",
                "Each pad changes only the trained BDT composition. Rows vary the photon-jet signal training sample, and columns vary the embedded inclusive MC training sample.",
                "The box in each pad reports AUC, logloss, the median signal-inclusive BDT-score gap, and the inclusive MC pass rate at the WP80 cut.",
                f"The cleanest full-composition AuAu point is Photon12+20 signal with Jet12+20+30+40 inclusive MC: {fmt_metrics(auau_best)}.",
                "The main readout is that adding harder embedded inclusive MC keeps pushing the inclusive score distribution lower and improves the WP80 inclusive-pass/logloss behavior, while adding Photon20 helps the full-composition row.",
                "",
            ]
        )
    )
    scripts["pp"].write_text(
        "\n".join(
            [
                "# pp sample-grid score-shape slide script",
                "",
                "This slide holds the pp validation sample fixed: 15-35 GeV with the full Photon5+10+20 signal and Jet8+12+20+30+40 inclusive MC mixture.",
                "Rows vary photon-jet signal training coverage, and columns vary inclusive MC training coverage.",
                "The full signal plus full inclusive MC pp model is the clean endpoint: " + fmt_metrics(pp_full) + ".",
                "The Photon5-only row is a useful failure mode, not a plotting mistake. On the full 15-35 GeV validation sample it has near-zero median gap and very large logloss; even with full inclusive MC its metrics are " + fmt_metrics(pp_photon5) + ".",
                "So the pp grid isolates the signal-coverage effect: high-threshold signal coverage is required before the BDT score split and calibration become trustworthy across the full validation window.",
                "",
            ]
        )
    )
    scripts["auau_signal_focus"].write_text(
        "\n".join(
            [
                "# AuAu full-inclusive-MC signal-variation score-shape script",
                "",
                "This focused slide fixes the blue comparison class to embedded inclusive MC: Jet12+20+30+40, with no truth filter. It is not a truth-tagged background sample.",
                "The two pads vary only the photon-jet signal training sample: Photon12 only versus Photon12+20.",
                "Both pads are scored on the same fixed AuAu validation sample: 0-20% centrality, 15-35 GeV, Photon12+20 signal, and Jet12+20+30+40 embedded inclusive MC.",
                "The Photon12-only model gives " + fmt_metrics(auau_photon12) + ".",
                "Adding Photon20 gives " + fmt_metrics(auau_best) + ".",
                "The takeaway is that holding inclusive MC fixed, adding Photon20 signal coverage raises AUC and slightly improves the median score gap and WP80 inclusive pass rate.",
                "",
            ]
        )
    )
    if signal_trainval_payload is not None:
        cfg = signal_trainval_payload["auau_signal_trainval"]
        cells = {(c["row_key"], c["col_key"]): c for c in cfg["cells"]}
        matched = cells[("photon12_20", "photon12_20")]["metrics"]
        transfer = cells[("photon12_20", "photon12")]["metrics"]
        scripts["auau_signal_trainval"].write_text(
            "\n".join(
                [
                    "# AuAu fixed-inclusive-MC signal train/validation matrix script",
                    "",
                    "This slide fixes the blue comparison class everywhere: embedded inclusive MC Jet12+20+30+40, with no truth filter, in both training and validation.",
                    "Rows vary the signal sample used to train the BDT. Columns vary the signal sample shown in validation.",
                    "The upper-left pad is Photon12 trained and Photon12 validated; the upper-right asks whether the Photon12-trained model transfers onto Photon12+20 validation.",
                    "The lower-left is the stress test in the other direction: the full Photon12+20 baseline model scored on Photon12-only signal validation.",
                    "The lower-right is the matched full-signal/full-inclusive baseline point: " + fmt_metrics(matched) + ".",
                    "The transfer-stress lower-left point is " + fmt_metrics(transfer) + "; it shows why the validation signal definition must be stated explicitly.",
                    "",
                ]
            )
        )
    if pp_signal_trainval_payload is not None:
        cfg = pp_signal_trainval_payload["pp_signal_trainval"]
        cells = {(c["row_key"], c["col_key"]): c for c in cfg["cells"]}
        matched = cells[("photon5_10_20", "photon5_10_20")]["metrics"]
        photon5_transfer = cells[("photon5", "photon5_10_20")]["metrics"]
        scripts["pp_signal_trainval"].write_text(
            "\n".join(
                [
                    "# pp fixed-inclusive-MC signal train/validation matrix script",
                    "",
                    "This slide fixes the blue comparison class everywhere: pp inclusive MC Jet8+12+20+30+40, with no truth filter, in both training and validation.",
                    "Rows vary the photon-jet signal samples used to train the BDT. Columns vary the photon-jet signal samples shown in validation.",
                    "The upper-right pad asks whether a Photon5-trained model transfers to the full Photon5+10+20 validation signal. That transfer point is " + fmt_metrics(photon5_transfer) + ".",
                    "The lower-right matched endpoint is Photon5+10+20 trained and Photon5+10+20 validated against the same fixed inclusive MC: " + fmt_metrics(matched) + ".",
                    "The point of the slide is source coverage, not truth-background tagging: the blue distribution is all selected inclusive-jet rows from the fixed Jet8+12+20+30+40 sample.",
                    "",
                ]
            )
        )
    if pp_inclusive_trainval_payload is not None:
        cfg = pp_inclusive_trainval_payload["pp_inclusive_trainval"]
        cells = {(c["row_key"], c["col_key"]): c for c in cfg["cells"]}
        matched = cells[("jet8_12_20_30_40", "jet8_12_20_30_40")]["metrics"]
        light_to_full = cells[("jet8_12_20_30_40", "jet8_12")]["metrics"]
        scripts["pp_inclusive_trainval"].write_text(
            "\n".join(
                [
                    "# pp fixed-signal inclusive-MC train/validation matrix script",
                    "",
                    "This slide fixes the red signal class everywhere: pp Photon5+10+20 signal in the 15-35 GeV window.",
                    "Rows vary the source-defined inclusive MC sample shown in validation. Columns vary the source-defined inclusive MC sample used to train the BDT.",
                    "The blue class is inclusive MC from the named Jet8/12/20/30/40 source rows with no truth-background filter.",
                    "The lower-right matched endpoint is full inclusive MC trained and full inclusive MC validated: " + fmt_metrics(matched) + ".",
                    "The lower-left transfer point asks whether a Jet8+12-trained BDT handles the full Jet8+12+20+30+40 validation mix: " + fmt_metrics(light_to_full) + ".",
                    "",
                ]
            )
        )
    return scripts


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--payload", type=Path, default=DEFAULT_PAYLOAD)
    ap.add_argument("--signal-trainval-payload", type=Path, default=DEFAULT_SIGNAL_TRAINVAL_PAYLOAD)
    ap.add_argument("--pp-signal-trainval-payload", type=Path, default=DEFAULT_PP_SIGNAL_TRAINVAL_PAYLOAD)
    ap.add_argument("--pp-inclusive-trainval-payload", type=Path, default=DEFAULT_PP_INCLUSIVE_TRAINVAL_PAYLOAD)
    ap.add_argument("--outdir", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    configure_fonts()
    payload = json.loads(args.payload.read_text())
    signal_trainval_payload = json.loads(args.signal_trainval_payload.read_text()) if args.signal_trainval_payload.exists() else None
    pp_signal_trainval_payload = json.loads(args.pp_signal_trainval_payload.read_text()) if args.pp_signal_trainval_payload.exists() else None
    pp_inclusive_trainval_payload = (
        json.loads(args.pp_inclusive_trainval_payload.read_text()) if args.pp_inclusive_trainval_payload.exists() else None
    )
    outputs = {
        "auau": render_grid(payload, "auau", args.outdir / "the37_auau_fixed_validation_score_shape_grid.png"),
        "pp": render_grid(payload, "pp", args.outdir / "the37_pp_fixed_validation_score_shape_grid.png"),
        "auau_signal_focus": render_auau_full_bkg_signal_focus(
            payload,
            args.outdir / "the37_auau_full_inclusive_signal_variation_score_shape_grid.png",
        ),
    }
    if signal_trainval_payload is not None:
        outputs["auau_signal_trainval"] = render_auau_signal_trainval_matrix(
            signal_trainval_payload,
            args.outdir / "the37_auau_fixed_inclusive_signal_trainval_matrix_score_shape_grid.png",
        )
    if pp_signal_trainval_payload is not None:
        outputs["pp_signal_trainval"] = render_pp_signal_trainval_matrix(
            pp_signal_trainval_payload,
            args.outdir / "the37_pp_fixed_inclusive_signal_trainval_matrix_score_shape_grid.png",
        )
    if pp_inclusive_trainval_payload is not None:
        outputs["pp_inclusive_trainval"] = render_pp_inclusive_trainval_matrix(
            pp_inclusive_trainval_payload,
            args.outdir / "the37_pp_fixed_signal_inclusive_trainval_matrix_score_shape_grid.png",
        )
    scripts = write_speaker_scripts(
        payload,
        args.outdir,
        signal_trainval_payload,
        pp_signal_trainval_payload,
        pp_inclusive_trainval_payload,
    )
    manifest = args.outdir / "the37_fixed_validation_score_shape_slides_manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "schema": "THE37_FIXED_VALIDATION_SCORE_SHAPE_SLIDES_V1",
                "payload": str(args.payload),
                "signal_trainval_payload": str(args.signal_trainval_payload) if signal_trainval_payload is not None else None,
                "pp_signal_trainval_payload": str(args.pp_signal_trainval_payload) if pp_signal_trainval_payload is not None else None,
                "pp_inclusive_trainval_payload": str(args.pp_inclusive_trainval_payload)
                if pp_inclusive_trainval_payload is not None
                else None,
                "outputs": {key: str(path) for key, path in outputs.items()},
                "speaker_scripts": {key: str(path) for key, path in scripts.items()},
                "font_family": FONT_FAMILY,
                "font_files": [str(p) for p in TIMES_FONTS if p.exists()],
                "size_px": [2560, 1440],
                "dpi": SLIDE_DPI,
                "slide_numbers_baked_in": False,
                "deck_mutated": False,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    for path in outputs.values():
        print(path)
    for path in scripts.values():
        print(path)
    print(manifest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
