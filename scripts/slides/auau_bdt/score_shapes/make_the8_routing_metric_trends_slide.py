#!/usr/bin/env python3
"""Render bin-wise routing-gain trends for corrected THE8 BDT diagnostics."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Rectangle
from matplotlib.ticker import MaxNLocator

THIS_FILE = Path(__file__).resolve()
SCRIPTS_DIR = next((p for p in THIS_FILE.parents if p.name == "scripts"), THIS_FILE.parent)
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.append(str(SCRIPTS_DIR))

from slides.common.slide_defaults import SLIDE_DPI, slide_figsize


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
BASE_OUT = REPO / "dataOutput/auauTightBDTValidation/THE59_model_comparison_remakes_20260620"
ET_PAYLOAD = BASE_OUT / "et_routing_20260622/the8_et_routing_score_payload_baseline_holdout.json"
CENT_PAYLOAD = BASE_OUT / "centrality_routing_20260622/the8_centrality_routing_fine_rows_payload_baseline_holdout.json"
OUT_DIR = BASE_OUT / "routing_trends_20260622"
DEFAULT_OUTPUT = OUT_DIR / "the8_routing_metric_trends_baseline_holdout_v1.png"

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
GRID = "#dfe7f1"
EDGE = "#c8d3e2"
ZERO = "#94a3b8"
ET_COLOR = "#0f9f8e"
COARSE_COLOR = "#d89b00"
FINE_COLOR = "#8a62c4"
GOOD_FILL = "#eaf7ed"
READOUT_FILL = "#fff7e6"
READOUT_EDGE = "#e7b85d"


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
    zorder: int = 6,
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
        zorder=zorder,
    )
    nodes.append((name, artist, size, role, flags))
    return artist


def rounded_box(fig: plt.Figure, x: float, y: float, w: float, h: float, *, face: str, edge: str, lw: float = 1.0) -> None:
    fig.patches.append(
        FancyBboxPatch(
            (x, y),
            w,
            h,
            transform=fig.transFigure,
            boxstyle="round,pad=0.006,rounding_size=0.006",
            facecolor=face,
            edgecolor=edge,
            linewidth=lw,
            zorder=1,
        )
    )


def cell_map(payload: dict) -> dict[tuple[str, str], dict]:
    return {(cell["row"], cell["col"]): cell for cell in payload["cells"]}


def bin_centers(payload: dict, axis: str) -> list[float]:
    cells = cell_map(payload)
    centers: list[float] = []
    for row in payload["row_order"]:
        cell = cells[(row, "global")]
        lo, hi = cell[axis]
        centers.append((lo + hi) / 2.0)
    return centers


def deltas(payload: dict, col: str, metric: str, *, percent_points: bool = False) -> list[float]:
    cells = cell_map(payload)
    values = []
    for row in payload["row_order"]:
        routed = cells[(row, col)]["metrics"][metric]
        global_value = cells[(row, "global")]["metrics"][metric]
        delta = routed - global_value
        values.append(100.0 * delta if percent_points else delta)
    return values


def draw_top(fig: plt.Figure, nodes: list[tuple[str, object, float, str, dict]]) -> None:
    add_text(
        fig,
        nodes,
        "title",
        0.040,
        0.946,
        "Routing gains are modest and coherent",
        size=34.5,
        role="title",
        va="top",
        weight="bold",
        title_anchor=True,
    )
    add_text(
        fig,
        nodes,
        "contract",
        0.042,
        0.864,
        "Each point is routed model minus global 14-feature BDT in the same baseline holdout bin.",
        size=15.6,
        color=MUTED,
        title_band_exception=True,
    )
    rounded_box(fig, 0.604, 0.815, 0.344, 0.058, face="white", edge=EDGE, lw=0.95)
    legend_items = [
        (0.622, ET_COLOR, "ET routed"),
        (0.728, COARSE_COLOR, "Cent3 routed"),
        (0.854, FINE_COLOR, "Cent7 routed"),
    ]
    for idx, (x, color, label) in enumerate(legend_items):
        fig.patches.append(Rectangle((x, 0.841), 0.036, 0.006, transform=fig.transFigure, facecolor=color, edgecolor=color, zorder=3))
        add_text(fig, nodes, f"legend {idx}", x + 0.044, 0.845, label, size=12.8, weight="bold", title_band_exception=True)


def draw_group_header(
    fig: plt.Figure,
    nodes: list[tuple[str, object, float, str, dict]],
    x: float,
    y: float,
    w: float,
    title: str,
    subtitle: str,
    edge: str,
    fill: str,
) -> None:
    rounded_box(fig, x, y, w, 0.052, face=fill, edge=edge, lw=1.15)
    add_text(fig, nodes, f"{title} header", x + w / 2.0, y + 0.032, title, size=17.0, ha="center", weight="bold")
    add_text(fig, nodes, f"{title} subheader", x + w / 2.0, y + 0.014, subtitle, size=10.8, ha="center", color=MUTED)


def style_axis(ax, *, bottom: bool) -> None:
    ax.set_facecolor("white")
    for side in ["top", "right"]:
        ax.spines[side].set_visible(False)
    for side in ["left", "bottom"]:
        ax.spines[side].set_color("#9fb1c6")
        ax.spines[side].set_linewidth(0.8)
    ax.grid(True, axis="y", color=GRID, linewidth=0.8)
    ax.grid(True, axis="x", color="#edf2f7", linewidth=0.6)
    ax.tick_params(axis="both", labelsize=9.2, colors=MUTED, length=3.0, width=0.8)
    if not bottom:
        ax.tick_params(labelbottom=False)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4))


def plot_metric(
    ax,
    x: list[float],
    series: list[tuple[str, list[float], str, str]],
    *,
    ylim: tuple[float, float],
    good_positive: bool,
    bottom: bool,
    xlabels: list[str] | None = None,
) -> None:
    ax.set_ylim(*ylim)
    if good_positive:
        ax.axhspan(0.0, ylim[1], color=GOOD_FILL, alpha=0.70, zorder=0)
    else:
        ax.axhspan(ylim[0], 0.0, color=GOOD_FILL, alpha=0.70, zorder=0)
    ax.axhline(0.0, color=ZERO, linewidth=1.1, zorder=1)
    for label, values, color, marker in series:
        ax.plot(
            x,
            values,
            color=color,
            linewidth=1.7,
            marker=marker,
            markersize=5.2,
            markerfacecolor=color,
            markeredgecolor="white",
            markeredgewidth=0.7,
            label=label,
            zorder=3,
        )
    ax.set_xlim(min(x) - 1.0, max(x) + 1.0)
    if xlabels:
        ax.set_xticks(x)
        ax.set_xticklabels(xlabels, fontsize=8.8, family=FONT, color=MUTED)
    style_axis(ax, bottom=bottom)


def build_summary(et_payload: dict, cent_payload: dict) -> tuple[str, str]:
    et_auc = deltas(et_payload, "et", "auc")
    et_ll = deltas(et_payload, "et", "logloss")
    et_fake = deltas(et_payload, "et", "wp80_inclusive_fake", percent_points=True)
    cent_cols = ["coarse", "fine"]
    cent_auc = [value for col in cent_cols for value in deltas(cent_payload, col, "auc")]
    cent_ll = [value for col in cent_cols for value in deltas(cent_payload, col, "logloss")]
    cent_fake = [value for col in cent_cols for value in deltas(cent_payload, col, "wp80_inclusive_fake", percent_points=True)]
    et_wins = sum(1 for a, ll, fake in zip(et_auc, et_ll, et_fake) if a > 0.0 and ll < 0.0 and fake < 0.0)
    cent_auc_ll_wins = sum(1 for a, ll in zip(cent_auc, cent_ll) if a > 0.0 and ll < 0.0)
    cent_fake_wins = sum(1 for fake in cent_fake if fake < 0.0)
    line1 = (
        f"ET routing improves AUC, logloss, and fake rate in {et_wins}/8 ET bins; "
        f"fake-rate change spans {min(et_fake):+.2f} to {max(et_fake):+.2f} pp."
    )
    line2 = (
        f"Centrality routing improves AUC and logloss in {cent_auc_ll_wins}/14 routed curves, "
        f"but fake-rate change is mixed with {cent_fake_wins}/14 below global."
    )
    return line1, line2


def render(output_path: Path = DEFAULT_OUTPUT) -> Path:
    configure_fonts()
    et_payload = json.loads(ET_PAYLOAD.read_text())
    cent_payload = json.loads(CENT_PAYLOAD.read_text())

    et_x = bin_centers(et_payload, "row_et_range")
    cent_x = bin_centers(cent_payload, "row_cent_range")
    et_labels = [row.replace(" GeV", "") for row in et_payload["row_order"]]
    cent_labels = [row.replace("%", "") for row in cent_payload["row_order"]]

    metrics = [
        ("AUC change", "auc", False, (-0.0004, 0.0042), True),
        ("Logloss change", "logloss", False, (-0.0086, 0.0006), False),
        ("Fake-rate change pp", "wp80_inclusive_fake", True, (-0.78, 0.62), False),
        ("Median-gap change", "median_gap", False, (-0.022, 0.036), True),
    ]

    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI, facecolor="white")
    nodes: list[tuple[str, object, float, str, dict]] = []
    draw_top(fig, nodes)

    left_x, right_x, panel_w = 0.150, 0.595, 0.350
    header_y = 0.750
    draw_group_header(fig, nodes, left_x, header_y, panel_w, "ET-bin routing", "cluster ET bins; base14_perEt minus global", ET_COLOR, "#ecfeff")
    draw_group_header(fig, nodes, right_x, header_y, panel_w, "Centrality-bin routing", "centrality bins; base14_perCent3/7 minus global", FINE_COLOR, "#f5f0ff")

    row_h = 0.123
    row_gap = 0.021
    top_y = 0.607
    axes = []
    for idx, (ylabel, metric, percent_points, ylim, good_positive) in enumerate(metrics):
        y = top_y - idx * (row_h + row_gap)
        ax_l = fig.add_axes([left_x, y, panel_w, row_h])
        ax_r = fig.add_axes([right_x, y, panel_w, row_h])
        bottom = idx == len(metrics) - 1
        add_text(fig, nodes, f"metric label {idx}", 0.042, y + row_h * 0.52, ylabel, size=12.3, weight="bold", color=INK, linespacing=0.92)
        et_series = [("ET routed", deltas(et_payload, "et", metric, percent_points=percent_points), ET_COLOR, "o")]
        cent_series = [
            ("Cent3 routed", deltas(cent_payload, "coarse", metric, percent_points=percent_points), COARSE_COLOR, "s"),
            ("Cent7 routed", deltas(cent_payload, "fine", metric, percent_points=percent_points), FINE_COLOR, "o"),
        ]
        plot_metric(ax_l, et_x, et_series, ylim=ylim, good_positive=good_positive, bottom=bottom, xlabels=et_labels if bottom else None)
        plot_metric(ax_r, cent_x, cent_series, ylim=ylim, good_positive=good_positive, bottom=bottom, xlabels=cent_labels if bottom else None)
        axes.extend([ax_l, ax_r])
    line1, line2 = build_summary(et_payload, cent_payload)
    rounded_box(fig, 0.042, 0.046, 0.918, 0.082, face=READOUT_FILL, edge=READOUT_EDGE, lw=1.0)
    add_text(fig, nodes, "readout title", 0.058, 0.091, "Readout", size=16.8, weight="bold")
    add_text(fig, nodes, "readout body one", 0.150, 0.100, line1, size=11.7, color=INK)
    add_text(fig, nodes, "readout body two", 0.150, 0.072, line2, size=11.7, color=INK)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=SLIDE_DPI)

    layout_path = output_path.with_suffix(".layout_nodes.json")
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
                "minimum_audience_font_px": font_px(10.4),
                "minimum_title_font_px": font_px(30.0),
                "minimum_plot_annotation_font_px": font_px(8.8),
                "nodes": layout_nodes,
            },
            indent=2,
        )
        + "\n"
    )
    manifest_path = output_path.with_suffix(".manifest.json")
    manifest_path.write_text(
        json.dumps(
            {
                "schema": "the8_routing_metric_trends_slide_v1",
                "png": str(output_path),
                "et_payload": str(ET_PAYLOAD),
                "centrality_payload": str(CENT_PAYLOAD),
                "layout_nodes": str(layout_path),
                "selection": et_payload.get("selection"),
                "labels": et_payload.get("labels"),
                "model_columns": {
                    "et": et_payload.get("model_columns"),
                    "centrality": cent_payload.get("model_columns"),
                },
                "readout": [line1, line2],
            },
            indent=2,
        )
        + "\n"
    )
    script_path = output_path.with_name(output_path.stem + "_speaker_script.md")
    script_path.write_text(
        "\n".join(
            [
                "# Routing gains are modest and coherent",
                "",
                "This slide shows bin-wise changes relative to the global 14-feature BDT, using the same reconstructed baseline holdout sample.",
                "",
                "On the left, the ET-routed model is compared with the global BDT in each cluster ET bin. It moves in the right direction for AUC, logloss, and fake rate in every ET bin, but the changes are modest.",
                "",
                "On the right, the centrality-routed models are compared with the same global BDT in fine centrality rows. Coarse and fine routing both improve AUC and logloss consistently, but the fake-rate movement is mixed, so centrality routing is not a clean fake-rate win by itself.",
                "",
                "The green half of each panel marks the favorable direction for that metric.",
                "",
            ]
        )
    )
    plt.close(fig)
    return output_path


def main() -> None:
    out = render()
    print(out)
    print(out.with_suffix(".manifest.json"))
    print(out.with_suffix(".layout_nodes.json"))
    print(out.with_name(out.stem + "_speaker_script.md"))


if __name__ == "__main__":
    main()
