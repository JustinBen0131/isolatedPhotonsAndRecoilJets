#!/usr/bin/env python3
"""Render bar-chart summary slides for corrected THE8 routed-BDT metrics."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from matplotlib.patches import FancyBboxPatch, Patch, Rectangle

THIS_FILE = Path(__file__).resolve()
SCRIPTS_DIR = next((p for p in THIS_FILE.parents if p.name == "scripts"), THIS_FILE.parent)
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.append(str(SCRIPTS_DIR))

from slides.common.slide_defaults import SLIDE_DPI, slide_figsize


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
BASE_OUT = REPO / "dataOutput/auauTightBDTValidation/THE59_model_comparison_remakes_20260620"
SOURCE_DIR = BASE_OUT / "routing_comprehensive_20260622"
PAYLOAD = SOURCE_DIR / "the8_routed_bdt_comprehensive_baseline_holdout_metrics_weighted.json"
OUT_DIR = BASE_OUT / "routing_bar_summary_20260622"

TIMES_DIR = Path("/System/Library/Fonts/Supplemental")
TIMES_FONTS = [
    TIMES_DIR / "Times New Roman.ttf",
    TIMES_DIR / "Times New Roman Bold.ttf",
    TIMES_DIR / "Times New Roman Italic.ttf",
    TIMES_DIR / "Times New Roman Bold Italic.ttf",
]
FONT = "Times New Roman"

INK = "#111827"
MUTED = "#64748b"
GRID = "#dbe4ef"
EDGE = "#c8d3e2"
READOUT_FILL = "#fff7e6"
READOUT_EDGE = "#e6b85d"
PANEL_FILL = "#fbfdff"
BASELINE_FILL = "#ffffff"
BASELINE_EDGE = "#ef4444"
BASELINE_TEXT = "#991b1b"
HIGHLIGHT_YELLOW = "#FDF287"

ROUTES = [
    ("cent3", "C3", "#7c3aed"),
    ("cent7", "C7", "#a855f7"),
    ("et", r"$E_T$", "#0284c7"),
    ("etcent3", r"$E_T$×C3", "#f59e0b"),
    ("etcent7", r"$E_T$×C7", "#15803d"),
]

ROUTE_BDT_COUNTS = {
    "cent3": 3,
    "cent7": 7,
    "et": 8,
    "etcent3": 24,
    "etcent7": 56,
}

METRICS = {
    "auc": {
        "title": "AUC improvement",
        "axis": "AUC improvement ×10³",
        "scale": 1000.0,
        "good_positive": True,
        "fmt": "{:+.1f}",
    },
    "logloss": {
        "title": "Logloss improvement",
        "axis": "logloss reduction ×10³",
        "scale": 1000.0,
        "good_positive": False,
        "fmt": "{:+.1f}",
    },
    "wp80_inclusive_fake": {
        "title": "Fake-rate improvement",
        "axis": "inclusive WP80 pass-rate reduction (percentage points)",
        "scale": 100.0,
        "good_positive": False,
        "fmt": "{:+.1f}",
    },
    "median_gap": {
        "title": "Score-gap improvement",
        "axis": "median-gap gain ×10²",
        "scale": 100.0,
        "good_positive": True,
        "fmt": "{:+.1f}",
    },
}


def configure_fonts() -> None:
    for font in TIMES_FONTS:
        if font.exists():
            font_manager.fontManager.addfont(str(font))
    plt.rcParams.update(
        {
            "font.family": FONT,
            "font.serif": [FONT, "Times", "DejaVu Serif"],
            "mathtext.fontset": "custom",
            "mathtext.rm": FONT,
            "mathtext.it": f"{FONT}:italic",
            "mathtext.bf": f"{FONT}:bold",
            "axes.unicode_minus": False,
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


def add_text(
    fig: plt.Figure,
    nodes: list[tuple[str, object, float, str, dict]],
    name: str,
    x: float,
    y: float,
    text: str,
    *,
    size: float,
    role: str = "audience",
    ha: str = "left",
    va: str = "center",
    weight: str = "normal",
    color: str = INK,
    linespacing: float = 1.0,
    **flags,
) -> object:
    artist = fig.text(
        x,
        y,
        text,
        ha=ha,
        va=va,
        fontsize=size,
        family=FONT,
        fontweight=weight,
        color=color,
        linespacing=linespacing,
        zorder=20,
    )
    nodes.append((name, artist, size, role, flags))
    return artist


def rounded_box(fig: plt.Figure, x: float, y: float, w: float, h: float, *, face: str, edge: str, lw: float = 1.0, radius: float = 0.010, zorder: int = -10) -> None:
    fig.patches.append(
        FancyBboxPatch(
            (x, y),
            w,
            h,
            transform=fig.transFigure,
            boxstyle=f"round,pad=0.006,rounding_size={radius}",
            facecolor=face,
            edgecolor=edge,
            linewidth=lw,
            zorder=zorder,
        )
    )


def add_title(
    fig: plt.Figure,
    nodes: list[tuple[str, object, float, str, dict]],
    title: str,
    tag: str | None = None,
    *,
    show_rule: bool = True,
) -> None:
    add_text(fig, nodes, "title", 0.040, 0.935, title, size=35.5, role="title", weight="bold", title_anchor=True)
    if not show_rule:
        return
    rule = "Bars show improvement versus the global 14-feature baseline on the same holdout; positive is better."
    add_text(fig, nodes, "reading rule", 0.043, 0.850, rule, size=13.0, color=MUTED, colon_style_exception=True)


def add_baseline_callout(fig: plt.Figure, nodes: list[tuple[str, object, float, str, dict]]) -> None:
    rounded_box(fig, 0.040, 0.826, 0.545, 0.060, face=BASELINE_FILL, edge=BASELINE_EDGE, lw=1.35, radius=0.004, zorder=5)
    add_text(
        fig,
        nodes,
        "baseline callout label",
        0.056,
        0.866,
        "RELATIVE TO BASELINE",
        size=12.8,
        weight="bold",
        color=BASELINE_TEXT,
        colon_style_exception=True,
    )
    add_text(
        fig,
        nodes,
        "baseline callout value",
        0.056,
        0.842,
        "14-feature baseV3E + cent + weta33/wphi33 baseline; positive = better",
        size=15.0,
        weight="bold",
        color=BASELINE_TEXT,
        colon_style_exception=True,
    )


def add_route_key(fig: plt.Figure, nodes: list[tuple[str, object, float, str, dict]]) -> None:
    legend_positions = [
        (0.600, 0.860, ROUTES[0]),
        (0.745, 0.860, ROUTES[1]),
        (0.870, 0.860, ROUTES[2]),
        (0.615, 0.829, ROUTES[3]),
        (0.790, 0.829, ROUTES[4]),
    ]
    for x, y, (route, route_label, color) in legend_positions:
        display_label = f"{route_label} ({ROUTE_BDT_COUNTS[route]} BDTs)"
        fig.patches.append(
            Rectangle(
                (x, y - 0.010),
                0.030,
                0.018,
                transform=fig.transFigure,
                facecolor=color,
                edgecolor="none",
                zorder=20,
            )
        )
        add_text(
            fig,
            nodes,
            f"route legend {display_label}",
            x + 0.037,
            y,
            display_label,
            size=11.9,
            weight="bold",
            color=INK,
        )
    add_text(
        fig,
        nodes,
        "route key line 1",
        0.055,
        0.776,
        "C3 = 0-20, 20-50, 50-80%",
        size=15.8,
        weight="bold",
        color=INK,
        colon_style_exception=True,
    )
    add_text(
        fig,
        nodes,
        "route key line 2",
        0.330,
        0.776,
        "C7 = 0-10, 10-20, 20-30, 30-40, 40-50, 50-60, 60-80%",
        size=15.8,
        weight="bold",
        color=INK,
        colon_style_exception=True,
    )
    add_text(
        fig,
        nodes,
        "route key line 3",
        0.055,
        0.735,
        r"$E_T$ = x-axis cluster-$E_T$ bin; $E_T$×C = $E_T$ bin × centrality bin",
        size=16.4,
        weight="bold",
        color=INK,
        colon_style_exception=True,
    )


def add_readout(fig: plt.Figure, nodes: list[tuple[str, object, float, str, dict]], lines: list[str]) -> None:
    rounded_box(fig, 0.050, 0.055, 0.910, 0.098, face=READOUT_FILL, edge=READOUT_EDGE, lw=1.0)
    add_text(fig, nodes, "readout title", 0.066, 0.109, "Readout", size=17.8, weight="bold")
    y0 = 0.122 if len(lines) > 1 else 0.104
    for i, line in enumerate(lines[:2]):
        add_text(
            fig,
            nodes,
            f"readout body {i + 1}",
            0.158,
            y0 - 0.033 * i,
            line,
            size=12.3,
            color=INK,
            colon_style_exception=True,
        )


def add_color_legend(
    fig: plt.Figure,
    routes: list[tuple[str, str, str]],
    *,
    x: float = 0.060,
    y: float = 0.810,
    ncol: int | None = None,
    fontsize: float = 11.5,
    columnspacing: float = 1.3,
) -> None:
    handles = [Patch(facecolor=color, edgecolor="none", label=label) for _key, label, color in routes]
    fig.legend(
        handles=handles,
        loc="upper left",
        bbox_to_anchor=(x, y),
        ncol=ncol or len(routes),
        frameon=False,
        fontsize=fontsize,
        columnspacing=columnspacing,
        handlelength=1.8,
        handletextpad=0.45,
    )


def style_axis(
    ax: plt.Axes,
    title: str,
    ylabel: str | None = None,
    *,
    title_size: float = 14.6,
    tick_size: float = 11.0,
) -> None:
    ax.set_facecolor(PANEL_FILL)
    ax.set_title(title, loc="left", fontsize=title_size, fontweight="bold", pad=7, color=INK)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=12.8, color=INK)
    ax.axhline(0.0, color="#334155", lw=1.0, zorder=2)
    ax.grid(axis="y", color=GRID, linewidth=1.0, zorder=1)
    ax.tick_params(axis="both", labelsize=tick_size, colors=INK)
    for spine in ax.spines.values():
        spine.set_edgecolor(EDGE)
        spine.set_linewidth(1.0)


def cell_lookup(payload: dict) -> dict[tuple[str, str, str], dict]:
    return {(c["group_kind"], c["group_label"], c["product"]): c for c in payload["cells"]}


def group_labels(payload: dict, kind: str) -> list[str]:
    labels: list[str] = []
    for cell in payload["cells"]:
        if cell["group_kind"] == kind and cell["product"] == "global" and cell["group_label"] not in labels:
            labels.append(cell["group_label"])
    return labels


def benefit(cells: dict[tuple[str, str, str], dict], kind: str, label: str, route: str, metric: str) -> float:
    spec = METRICS[metric]
    base = cells[(kind, label, "global")]["metrics"][metric]
    cur = cells[(kind, label, route)]["metrics"][metric]
    delta = float(cur - base)
    signed = delta if spec["good_positive"] else -delta
    return signed * spec["scale"]


def values_for(payload: dict, kind: str, labels: list[str], metric: str, routes: list[tuple[str, str, str]] = ROUTES) -> np.ndarray:
    cells = cell_lookup(payload)
    out = np.zeros((len(routes), len(labels)), dtype="float64")
    for i, (route, _label, _color) in enumerate(routes):
        for j, label in enumerate(labels):
            out[i, j] = benefit(cells, kind, label, route, metric)
    return out


def label_for_group(label: str) -> str:
    return label.replace(" GeV", "").replace("%", "")


def draw_grouped_bars(
    ax: plt.Axes,
    labels: list[str],
    values: np.ndarray,
    *,
    routes: list[tuple[str, str, str]],
    title: str,
    ylabel: str,
    value_labels: bool = False,
    legend: bool = False,
    rotate: int = 0,
    show_xticklabels: bool = True,
    show_ylabel: bool = True,
    force_zero_floor: bool = False,
    emphasize_bins: bool = False,
    xtick_size: float = 11.0,
    xtick_weight: str = "normal",
    xlabel: str | None = None,
    title_size: float = 14.6,
    tick_size: float = 11.0,
    xlabel_size: float = 13.5,
) -> None:
    style_axis(ax, title, ylabel if show_ylabel else None, title_size=title_size, tick_size=tick_size)
    x = np.arange(len(labels))
    n = len(routes)
    width = min(0.78 / n, 0.16)
    offsets = (np.arange(n) - (n - 1) / 2.0) * width
    ymax = 0.0
    ymin = 0.0
    if emphasize_bins:
        for j in range(len(labels)):
            if j % 2 == 0:
                ax.axvspan(j - 0.5, j + 0.5, color="#eef4fb", alpha=0.55, zorder=0)
    for i, (_route, route_label, color) in enumerate(routes):
        bars = ax.bar(x + offsets[i], values[i], width=width * 0.92, color=color, label=route_label, zorder=3)
        ymax = max(ymax, float(np.nanmax(values[i])))
        ymin = min(ymin, float(np.nanmin(values[i])))
        if value_labels:
            for bar, val in zip(bars, values[i]):
                va = "bottom" if val >= 0 else "top"
                y = bar.get_height() + (0.04 * max(abs(ymax), abs(ymin), 1.0) if val >= 0 else -0.04 * max(abs(ymax), abs(ymin), 1.0))
                ax.text(
                    bar.get_x() + bar.get_width() / 2.0,
                    y,
                    f"{val:+.1f}",
                    ha="center",
                    va=va,
                    fontsize=8.7,
                    fontweight="bold",
                    color=INK,
                    rotation=0,
                    zorder=4,
                )
    if force_zero_floor:
        lim = max(ymax, 0.8)
        ax.set_ylim(0, lim * 1.18)
    elif ymin >= 0:
        lim = max(ymax, 0.8)
        ax.set_ylim(0, lim * 1.22)
    else:
        upper = max(ymax, 0.8)
        lower = min(ymin * 1.45, -0.08 * upper)
        ax.set_ylim(lower, upper * 1.22)
    ax.set_xticks(x)
    if show_xticklabels:
        ax.set_xticklabels([label_for_group(l) for l in labels], rotation=rotate, ha="right" if rotate else "center")
        ax.tick_params(axis="x", labelsize=xtick_size, pad=5)
        for tick in ax.get_xticklabels():
            tick.set_fontweight(xtick_weight)
    else:
        ax.set_xticklabels([])
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=xlabel_size, fontweight="bold", color=INK, labelpad=7)
    if legend:
        ax.legend(
            loc="upper left",
            bbox_to_anchor=(0.0, 1.02),
            ncol=len(routes),
            frameon=False,
            fontsize=10.7,
            columnspacing=1.0,
            handlelength=1.6,
        )


def draw_horizontal_route_bars(ax: plt.Axes, payload: dict, metric: str, *, title: str) -> None:
    cells = cell_lookup(payload)
    label = group_labels(payload, "all")[0]
    route_labels = [r[1] for r in ROUTES]
    colors = [r[2] for r in ROUTES]
    values = np.asarray([benefit(cells, "all", label, r[0], metric) for r in ROUTES], dtype="float64")
    ax.set_facecolor(PANEL_FILL)
    ax.set_title(METRICS[metric]["axis"], loc="left", fontsize=14.8, fontweight="bold", pad=7, color=INK)
    y = np.arange(len(ROUTES))
    bars = ax.barh(y, values, color=colors, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(route_labels, fontsize=12.3, fontweight="bold")
    ax.invert_yaxis()
    vmax = max(float(np.nanmax(np.abs(values))), 0.8)
    ax.set_xlim(min(0, float(np.nanmin(values))) - 0.14 * vmax, float(np.nanmax(values)) + 0.24 * vmax)
    ax.axvline(0.0, color="#334155", lw=1.0, zorder=2)
    ax.grid(axis="x", color=GRID, linewidth=1.0, zorder=1)
    ax.grid(axis="y", visible=False)
    ax.tick_params(axis="x", labelsize=10.8, colors=INK)
    ax.tick_params(axis="y", colors=INK)
    for spine in ax.spines.values():
        spine.set_edgecolor(EDGE)
        spine.set_linewidth(1.0)
    for bar, val in zip(bars, values):
        x = val + (0.035 * vmax if val >= 0 else -0.035 * vmax)
        ha = "left" if val >= 0 else "right"
        ax.text(x, bar.get_y() + bar.get_height() / 2.0, f"{val:+.1f}", va="center", ha=ha, fontsize=11.0, fontweight="bold", color=INK)


def save_slide(fig: plt.Figure, nodes: list[tuple[str, object, float, str, dict]], output: Path, payload: dict, *, extra: dict) -> Path:
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=SLIDE_DPI)
    layout_path = output.with_suffix(".layout_nodes.json")
    layout_nodes = []
    for name, artist, points, role, flags in nodes:
        layout_nodes.append(
            {
                "name": name,
                "kind": "text",
                "role": role,
                "text": artist.get_text(),
                "font_px": font_px(points),
                "bbox": bbox_px(fig, artist),
                "title_anchor": role == "title",
                **flags,
            }
        )
    layout_path.write_text(
        json.dumps(
            {
                "schema": "slide_layout_nodes_v1",
                "title_axis_x": 0.040 * 2560,
                "minimum_audience_font_px": font_px(10.5),
                "minimum_title_font_px": font_px(35.0),
                "minimum_plot_annotation_font_px": font_px(8.5),
                "nodes": layout_nodes,
            },
            indent=2,
        )
        + "\n"
    )
    manifest_path = output.with_suffix(".manifest.json")
    manifest = {
        "schema": "the8_routing_bar_summary_slide_v1",
        "png": str(output),
        "payload": str(PAYLOAD),
        "layout_nodes": str(layout_path),
        "metric_convention": "All bars are signed as improvements relative to the global 14-feature baseline; positive is better.",
        "source": payload.get("source"),
        "selection": payload.get("selection"),
        "labels": payload.get("labels"),
        "metric_weighting": payload.get("metric_weighting"),
        **extra,
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    plt.close(fig)
    return output


def write_bar_datapoints(payload: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    labels_by_kind = {
        "all": group_labels(payload, "all"),
        "coarse_cent": group_labels(payload, "coarse_cent"),
        "fine_cent": group_labels(payload, "fine_cent"),
        "et": group_labels(payload, "et"),
        "et_coarse_cent": group_labels(payload, "et_coarse_cent"),
        "et_fine_cent": group_labels(payload, "et_fine_cent"),
    }
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["group_kind", "group_label", "route", "metric", "improvement_value", "axis_units"])
        cells = cell_lookup(payload)
        for kind, labels in labels_by_kind.items():
            for label in labels:
                for route, _route_label, _color in ROUTES:
                    if (kind, label, route) not in cells:
                        continue
                    for metric, spec in METRICS.items():
                        writer.writerow([kind, label, route, metric, benefit(cells, kind, label, route, metric), spec["axis"]])


def overall_slide(payload: dict) -> Path:
    configure_fonts()
    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI, facecolor="white")
    nodes: list[tuple[str, object, float, str, dict]] = []
    add_title(fig, nodes, "Routed BDT route ranking", "same holdout | positive bars mean better")
    positions = [
        [0.075, 0.505, 0.395, 0.245],
        [0.555, 0.505, 0.385, 0.245],
        [0.075, 0.205, 0.395, 0.245],
        [0.555, 0.205, 0.385, 0.245],
    ]
    for metric, pos in zip(["auc", "logloss", "wp80_inclusive_fake", "median_gap"], positions):
        ax = fig.add_axes(pos)
        draw_horizontal_route_bars(ax, payload, metric, title=METRICS[metric]["title"])
    add_readout(
        fig,
        nodes,
        [
            "ET×fine centrality is the strongest overall route on every main diagnostic shown here.",
            "Use this as the executive slide; the next slides show where the gains sit in centrality and ET.",
        ],
    )
    return save_slide(fig, nodes, OUT_DIR / "the8_routing_bar_overall_route_ranking_v1.png", payload, extra={"slide_role": "overall_route_ranking"})


def centrality_slide(payload: dict, *, metrics: tuple[str, str], suffix: str, title: str, readout: list[str]) -> Path:
    configure_fonts()
    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI, facecolor="white")
    nodes: list[tuple[str, object, float, str, dict]] = []
    add_title(fig, nodes, title, "coarse + fine centrality bins | same holdout")
    coarse = group_labels(payload, "coarse_cent")
    fine = group_labels(payload, "fine_cent")
    positions = [
        [0.065, 0.500, 0.415, 0.235],
        [0.545, 0.500, 0.405, 0.235],
        [0.065, 0.205, 0.415, 0.235],
        [0.545, 0.205, 0.405, 0.235],
    ]
    add_color_legend(fig, ROUTES, y=0.810)
    panel_specs = [
        ("coarse_cent", coarse, metrics[0], f"Coarse centrality: {METRICS[metrics[0]]['title']}", True),
        ("fine_cent", fine, metrics[0], f"Fine centrality: {METRICS[metrics[0]]['title']}", False),
        ("coarse_cent", coarse, metrics[1], f"Coarse centrality: {METRICS[metrics[1]]['title']}", True),
        ("fine_cent", fine, metrics[1], f"Fine centrality: {METRICS[metrics[1]]['title']}", False),
    ]
    for pos, (kind, labels, metric, panel_title, labels_on) in zip(positions, panel_specs):
        ax = fig.add_axes(pos)
        vals = values_for(payload, kind, labels, metric)
        draw_grouped_bars(
            ax,
            labels,
            vals,
            routes=ROUTES,
            title=panel_title,
            ylabel=METRICS[metric]["axis"],
            value_labels=labels_on,
            legend=False,
            rotate=0,
        )
    add_readout(fig, nodes, readout)
    return save_slide(fig, nodes, OUT_DIR / f"the8_routing_bar_centrality_{suffix}_v1.png", payload, extra={"slide_role": f"centrality_{suffix}", "metrics": metrics})


def et_slide(payload: dict, *, metrics: tuple[str, str], suffix: str, title: str, readout: list[str]) -> Path:
    configure_fonts()
    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI, facecolor="white")
    nodes: list[tuple[str, object, float, str, dict]] = []
    add_title(fig, nodes, title, "ET bins | same holdout | positive bars mean better")
    labels = group_labels(payload, "et")
    positions = [
        [0.065, 0.470, 0.885, 0.235],
        [0.065, 0.195, 0.885, 0.205],
    ]
    add_color_legend(fig, ROUTES, y=0.810)
    for idx, (pos, metric) in enumerate(zip(positions, metrics)):
        ax = fig.add_axes(pos)
        vals = values_for(payload, "et", labels, metric)
        draw_grouped_bars(
            ax,
            labels,
            vals,
            routes=ROUTES,
            title=METRICS[metric]["title"],
            ylabel=METRICS[metric]["axis"],
            value_labels=False,
            legend=False,
            rotate=0,
            show_xticklabels=idx == 1,
        )
    add_readout(fig, nodes, readout)
    return save_slide(fig, nodes, OUT_DIR / f"the8_routing_bar_et_{suffix}_v1.png", payload, extra={"slide_role": f"et_{suffix}", "metrics": metrics})


def et_cent_slide(payload: dict, *, kind: str, route: tuple[str, str, str], suffix: str, title: str, readout: list[str]) -> Path:
    configure_fonts()
    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI, facecolor="white")
    nodes: list[tuple[str, object, float, str, dict]] = []
    add_title(fig, nodes, title, f"{route[1]} only | ET × centrality cells")
    et_labels = group_labels(payload, "et")
    if kind == "et_coarse_cent":
        cent_labels = group_labels(payload, "coarse_cent")
        colors = ["#2563eb", "#f59e0b", "#16a34a"]
    else:
        cent_labels = group_labels(payload, "fine_cent")
        colors = ["#2563eb", "#0ea5e9", "#14b8a6", "#22c55e", "#a3e635", "#f59e0b", "#ef4444"]
    routes = [(route[0], label_for_group(c), colors[i % len(colors)]) for i, c in enumerate(cent_labels)]
    cells = cell_lookup(payload)

    def matrix(metric: str) -> np.ndarray:
        out = np.zeros((len(cent_labels), len(et_labels)), dtype="float64")
        for i, cent in enumerate(cent_labels):
            for j, et in enumerate(et_labels):
                group_label = f"{et} × {cent}"
                out[i, j] = benefit(cells, kind, group_label, route[0], metric)
        return out

    positions = [[0.065, 0.470, 0.885, 0.235], [0.065, 0.195, 0.885, 0.205]]
    add_color_legend(fig, routes, y=0.810, ncol=min(len(routes), 7))
    for idx, (pos, metric) in enumerate(zip(positions, ["auc", "wp80_inclusive_fake"])):
        ax = fig.add_axes(pos)
        draw_grouped_bars(
            ax,
            et_labels,
            matrix(metric),
            routes=routes,
            title=METRICS[metric]["title"],
            ylabel=METRICS[metric]["axis"],
            value_labels=False,
            legend=False,
            rotate=0,
            show_xticklabels=idx == 1,
        )
    add_readout(fig, nodes, readout)
    return save_slide(fig, nodes, OUT_DIR / f"the8_routing_bar_et_cent_{suffix}_v1.png", payload, extra={"slide_role": f"et_cent_{suffix}", "route": route[0]})


def et_0_20_three_metric_slide(payload: dict) -> Path:
    configure_fonts()
    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI, facecolor="white")
    nodes: list[tuple[str, object, float, str, dict]] = []
    add_title(fig, nodes, r"0-20% centrality: routed-BDT gains versus $E_T$", show_rule=False)
    add_baseline_callout(fig, nodes)
    add_route_key(fig, nodes)
    cells = cell_lookup(payload)
    labels = [label for label in group_labels(payload, "et_coarse_cent") if label.endswith("0-20%")]
    et_labels = [label.split(" × ", 1)[0] for label in labels]
    positions = [
        [0.065, 0.500, 0.885, 0.140],
        [0.065, 0.305, 0.885, 0.130],
        [0.065, 0.120, 0.885, 0.130],
    ]
    panel_titles = {
        "auc": "AUC improvement ×10³",
        "wp80_inclusive_fake": "Inclusive WP80 pass-rate reduction (percentage points)",
        "logloss": "Logloss reduction ×10³",
    }
    for idx, (pos, metric) in enumerate(zip(positions, ["auc", "wp80_inclusive_fake", "logloss"])):
        ax = fig.add_axes(pos)
        vals = np.zeros((len(ROUTES), len(labels)), dtype="float64")
        for i, (route, _route_label, _color) in enumerate(ROUTES):
            for j, label in enumerate(labels):
                vals[i, j] = benefit(cells, "et_coarse_cent", label, route, metric)
        draw_grouped_bars(
            ax,
            et_labels,
            vals,
            routes=ROUTES,
            title=panel_titles[metric],
            ylabel=METRICS[metric]["axis"],
            value_labels=False,
            legend=False,
            rotate=0,
            show_xticklabels=idx == 2,
            show_ylabel=False,
            force_zero_floor=True,
            emphasize_bins=True,
            xtick_size=13.2,
            xtick_weight="bold",
            xlabel=r"cluster $E_T$ bin (GeV)" if idx == 2 else None,
        )
    return save_slide(
        fig,
        nodes,
        OUT_DIR / "the8_routing_bar_0_20_et_auc_fake_logloss_v1.png",
        payload,
        extra={"slide_role": "0_20_et_auc_fake_logloss", "group_kind": "et_coarse_cent", "centrality": "0-20%"},
    )


def _text_width_frac(fig: plt.Figure, text: str, size: float, weight: str = "normal") -> float:
    """Measure rendered text width as a fraction of figure width."""
    probe = fig.text(0, -1, text, fontsize=size, family=FONT, fontweight=weight)
    fig.canvas.draw()
    width = probe.get_window_extent(renderer=fig.canvas.get_renderer()).width / fig.bbox.width
    probe.remove()
    return float(width)


def add_audience_route_key(fig: plt.Figure, nodes: list, y: float, *, size: float = 15.5) -> None:
    """Centred single-row legend using the compact C3/C7/ET shorthand.

    Widths are measured from rendered glyphs, so swatches never collide with
    the neighbouring label.  Route order and colours are untouched.
    """
    entries = [(r, f"{short} ({ROUTE_BDT_COUNTS[r]} BDTs)", c) for r, short, c in ROUTES]
    swatch_w, text_pad, gap = 0.020, 0.007, 0.030
    widths = [_text_width_frac(fig, lab, size) for _r, lab, _c in entries]
    total = sum(swatch_w + text_pad + w for w in widths) + gap * (len(entries) - 1)
    x = (1.0 - total) / 2.0
    for (route, label, color), w in zip(entries, widths):
        fig.patches.append(
            Rectangle((x, y - 0.010), swatch_w, 0.019, transform=fig.transFigure,
                      facecolor=color, edgecolor="none", zorder=20)
        )
        add_text(fig, nodes, f"legend {route}", x + swatch_w + text_pad, y, label, size=size)
        x += swatch_w + text_pad + w + gap


def et_0_20_three_metric_audience_slide(payload: dict) -> Path:
    """Audience restyle of the 0-20% routed-BDT slide. Presentation only."""
    configure_fonts()
    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI, facecolor="white")
    nodes: list[tuple[str, object, float, str, dict]] = []

    title_text = "Au+Au BDT design choice: global vs routed models in 0–20% centrality"
    title_size = 33.0
    while title_size > 18.0 and _text_width_frac(fig, title_text, title_size, "bold") > 0.920:
        title_size -= 0.5
    add_text(fig, nodes, "title", 0.040, 0.945, title_text, size=title_size,
             role="title", weight="bold", title_anchor=True)

    fig.text(0.040, 0.888, "▶", fontsize=15.0, color="#2468A8", ha="left",
             va="center", fontfamily="DejaVu Sans")
    add_text(fig, nodes, "subtitle", 0.061, 0.888,
             r"Finer $E_T\times$centrality routing improves all three metrics; one global BDT remains nominal for simplicity",
             size=17.5, color=MUTED, colon_style_exception=True)

    # Reference line: centred, red outline, no fill; "Reference:" bold.
    ref_size = 15.0
    ref_bold, ref_rest = "Reference:", r" global BDT — 0–80%, $15 \leq E_T < 35$ GeV"
    w_bold = _text_width_frac(fig, ref_bold, ref_size, "bold")
    w_rest = _text_width_frac(fig, ref_rest, ref_size)
    ref_y = 0.827
    ref_x0 = 0.5 - (w_bold + w_rest) / 2.0
    rounded_box(fig, ref_x0 - 0.012, ref_y - 0.021, w_bold + w_rest + 0.024, 0.042,
                face="none", edge=BASELINE_EDGE, lw=1.35, radius=0.004, zorder=5)
    add_text(fig, nodes, "reference label", ref_x0, ref_y, ref_bold, size=ref_size,
             weight="bold", colon_style_exception=True)
    add_text(fig, nodes, "reference value", ref_x0 + w_bold, ref_y, ref_rest,
             size=ref_size, color=MUTED, colon_style_exception=True)

    add_audience_route_key(fig, nodes, 0.768)
    add_text(fig, nodes, "routing notation", 0.5, 0.726,
             "C3 = {0–20, 20–50, 50–80%}    C7 = {0–10, 10–20, …, 60–80%}    $E_T$ = x-axis cluster-$E_T$ bin",
             size=13.5, color=MUTED, ha="center", colon_style_exception=True)

    cells = cell_lookup(payload)
    labels = [l for l in group_labels(payload, "et_coarse_cent") if l.endswith("0-20%")]
    et_labels = [l.split(" × ", 1)[0] for l in labels]
    positions = [
        [0.108, 0.500, 0.842, 0.163],
        [0.108, 0.292, 0.842, 0.163],
        [0.108, 0.084, 0.842, 0.163],
    ]
    panel_titles = {
        "auc": r"AUC improvement $\times10^{3}$",
        "wp80_inclusive_fake": "Background acceptance reduction at WP80 [percentage points]",
        "logloss": r"Log-loss reduction $\times10^{3}$",
    }
    for idx, (pos, metric) in enumerate(zip(positions, ["auc", "wp80_inclusive_fake", "logloss"])):
        ax = fig.add_axes(pos)
        vals = np.zeros((len(ROUTES), len(labels)), dtype="float64")
        for i, (route, _rl, _c) in enumerate(ROUTES):
            for j, label in enumerate(labels):
                vals[i, j] = benefit(cells, "et_coarse_cent", label, route, metric)
        draw_grouped_bars(
            ax, et_labels, vals, routes=ROUTES, title=panel_titles[metric],
            ylabel=METRICS[metric]["axis"], value_labels=False, legend=False, rotate=0,
            show_xticklabels=idx == 2, show_ylabel=False, force_zero_floor=True,
            emphasize_bins=True, xtick_size=15.0, xtick_weight="bold",
            xlabel=r"Cluster $E_T$ bin [GeV]" if idx == 2 else None,
            title_size=17.0, tick_size=13.2, xlabel_size=16.0,
        )

    # Rotated spanning label on the left, reading like a shared y-axis title.
    span_lo, span_hi = positions[2][1], positions[0][1] + positions[0][3]
    span_mid = (span_lo + span_hi) / 2.0
    label_size = 13.5
    label_text = "Positive values = improvement wrt reference"
    # Rotated 90 deg: the measured *width* fraction spans figure *height*, so
    # rescale by the 2560x1440 aspect before sizing the highlight.
    label_len = _text_width_frac(fig, label_text, label_size, "bold") * (2560.0 / 1440.0)
    # Yellow highlighter behind bold black text, matching the deck's callout style.
    highlight_x, highlight_w = 0.029, 0.012
    rounded_box(fig, highlight_x, span_mid - label_len / 2.0 - 0.002, highlight_w, label_len + 0.004,
                face=HIGHLIGHT_YELLOW, edge=HIGHLIGHT_YELLOW, lw=0.8, radius=0.002, zorder=5)
    artist = fig.text(highlight_x + highlight_w / 2.0, span_mid, label_text,
                      fontsize=label_size, family=FONT, color=INK, fontweight="bold",
                      ha="center", va="center", rotation=90, zorder=20)
    nodes.append(("spanning improvement label", artist, label_size, "audience",
                  {"colon_style_exception": True}))

    return save_slide(
        fig, nodes,
        OUT_DIR / "the8_routing_bar_0_20_et_auc_fake_logloss_v2_audience.png",
        payload,
        extra={"slide_role": "0_20_et_auc_fake_logloss_audience",
               "group_kind": "et_coarse_cent", "centrality": "0-20%",
               "presentation_only": True},
    )


def et_coarse_cent_auc_rows_slide(payload: dict) -> Path:
    configure_fonts()
    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI, facecolor="white")
    nodes: list[tuple[str, object, float, str, dict]] = []
    add_title(fig, nodes, r"AUC improvement versus $E_T$ by centrality", show_rule=False)
    add_baseline_callout(fig, nodes)
    add_route_key(fig, nodes)
    cells = cell_lookup(payload)
    cent_labels = group_labels(payload, "coarse_cent")
    et_labels = group_labels(payload, "et")
    matrices: list[np.ndarray] = []
    for cent in cent_labels:
        vals = np.zeros((len(ROUTES), len(et_labels)), dtype="float64")
        for i, (route, _route_label, _color) in enumerate(ROUTES):
            for j, et_label in enumerate(et_labels):
                vals[i, j] = benefit(cells, "et_coarse_cent", f"{et_label} × {cent}", route, "auc")
        matrices.append(vals)
    common_top = max(float(np.nanmax(vals)) for vals in matrices)
    common_top = max(common_top, 0.8) * 1.18
    positions = [
        [0.065, 0.505, 0.885, 0.160],
        [0.065, 0.300, 0.885, 0.160],
        [0.065, 0.095, 0.885, 0.160],
    ]
    for idx, (pos, cent, vals) in enumerate(zip(positions, cent_labels, matrices)):
        ax = fig.add_axes(pos)
        draw_grouped_bars(
            ax,
            et_labels,
            vals,
            routes=ROUTES,
            title=f"{cent} centrality | AUC improvement ×10³",
            ylabel=METRICS["auc"]["axis"],
            value_labels=False,
            legend=False,
            rotate=0,
            show_xticklabels=idx == 2,
            show_ylabel=False,
            force_zero_floor=True,
            emphasize_bins=True,
            xtick_size=13.2,
            xtick_weight="bold",
            xlabel=r"cluster $E_T$ bin (GeV)" if idx == 2 else None,
        )
        ax.set_ylim(0, common_top)
    return save_slide(
        fig,
        nodes,
        OUT_DIR / "the8_routing_bar_auc_et_by_coarse_centrality_v1.png",
        payload,
        extra={"slide_role": "auc_et_by_coarse_centrality", "group_kind": "et_coarse_cent"},
    )


def write_script(outputs: list[Path]) -> None:
    script = OUT_DIR / "the8_routing_bar_summary_speaker_script.md"
    script.write_text(
        "\n".join(
            [
                "# Corrected THE8 routed-BDT bar-summary slides",
                "",
                "These slides convert the comprehensive fixed-holdout routed-BDT metrics into bar charts.",
                "Every bar is an improvement relative to the global 14-feature baseline, so positive is better.",
                "The class definition is source-class Signal MC from embeddedPhoton rows with is_signal==1 versus inclusive embeddedJet MC with no truth-background filter.",
                "The display weights use the PPG12-style source-class diagnostic weights recorded in the payload.",
                "",
                "Suggested flow:",
                "1. Start with the route-ranking slide to show the global ordering.",
                "2. Use the centrality and ET slides to show that the improvement is not confined to one bin.",
                "3. Use the ET×centrality slides only if the audience asks where the two-dimensional routing is doing the work.",
                "",
                "Generated PNGs:",
                *[f"- {path}" for path in outputs],
                "",
            ]
        )
        + "\n"
    )


def render_all() -> list[Path]:
    payload = json.loads(PAYLOAD.read_text())
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    write_bar_datapoints(payload, OUT_DIR / "the8_routing_bar_summary_datapoints.csv")
    outputs = [
        overall_slide(payload),
        centrality_slide(
            payload,
            metrics=("auc", "wp80_inclusive_fake"),
            suffix="primary",
            title="Centrality view: separation and fake-rate gains",
            readout=[
                "ET×fine centrality gives the largest centrality-bin gains; fine centrality alone also beats coarse-only routing.",
                "This slide is the centrality-facing answer: where AUC rises and inclusive-MC leakage falls.",
            ],
        ),
        centrality_slide(
            payload,
            metrics=("logloss", "median_gap"),
            suffix="calibration_gap",
            title="Centrality view: calibration and score-gap checks",
            readout=[
                "Logloss reduction tracks the AUC/fake-rate story, so the stronger route is not only stretching the score scale.",
                "Median-gap bars show how much the signal-inclusive score split opens in each centrality bin.",
            ],
        ),
        et_slide(
            payload,
            metrics=("auc", "wp80_inclusive_fake"),
            suffix="primary",
            title="ET view: separation and fake-rate gains",
            readout=[
                "ET-aware routing improves every ET bin, with ET×fine centrality remaining the highest route across the full window.",
                "Use this when the question is whether routing only helps one cluster-ET region.",
            ],
        ),
        et_slide(
            payload,
            metrics=("logloss", "median_gap"),
            suffix="calibration_gap",
            title="ET view: calibration and score-gap checks",
            readout=[
                "The calibration diagnostic improves coherently with ET-aware routing, especially once centrality is included.",
                "The score gap stays positive in the ET view, which supports the same conclusion without relying only on AUC.",
            ],
        ),
        et_cent_slide(
            payload,
            kind="et_coarse_cent",
            route=("etcent3", "ET×C3", "#f59e0b"),
            suffix="coarse_primary",
            title="ET × coarse centrality: where the 24 routed BDTs help",
            readout=[
                "Each ET bin is split by the three coarse centrality regions; bars show improvement versus global.",
                "This is the compact diagnostic for the ET×coarse route before going to the fine 56-model view.",
            ],
        ),
        et_cent_slide(
            payload,
            kind="et_fine_cent",
            route=("etcent7", "ET×C7", "#15803d"),
            suffix="fine_primary",
            title="ET × fine centrality: where the 56 routed BDTs help",
            readout=[
                "The fine route is positive across the two-dimensional grid, not just in the all-rows aggregate.",
                "Use this as the detailed backup when explaining why ET×fine is the strongest candidate route.",
            ],
        ),
        et_0_20_three_metric_slide(payload),
        et_coarse_cent_auc_rows_slide(payload),
    ]
    write_script(outputs)
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--only-audience-0-20", action="store_true",
                        help="Render only the audience restyle; leaves all other outputs untouched.")
    args = parser.parse_args()
    if args.only_audience_0_20:
        payload = json.loads(PAYLOAD.read_text())
        OUT_DIR.mkdir(parents=True, exist_ok=True)
        print(et_0_20_three_metric_audience_slide(payload))
        return
    for output in render_all():
        print(output)


if __name__ == "__main__":
    main()
