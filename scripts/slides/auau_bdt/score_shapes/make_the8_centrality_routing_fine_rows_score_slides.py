#!/usr/bin/env python3
"""Render fine-centrality score-shape slides for centrality-routed THE8 BDTs."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import FancyBboxPatch, Rectangle
from matplotlib.ticker import MaxNLocator

THIS_FILE = Path(__file__).resolve()
SCRIPTS_DIR = next((p for p in THIS_FILE.parents if p.name == "scripts"), THIS_FILE.parent)
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.append(str(SCRIPTS_DIR))

from slides.common.slide_defaults import SLIDE_DPI, slide_figsize


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
OUT_DIR = (
    REPO
    / "dataOutput/auauTightBDTValidation/THE59_model_comparison_remakes_20260620"
    / "centrality_routing_20260622"
)
DEFAULT_PAYLOAD = OUT_DIR / "the8_centrality_routing_fine_rows_payload_baseline_holdout.json"
OUTPUT_STEM = "the8_centrality_routing_fine_rows_score_shapes_baseline_holdout"

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
RED = "#d13b35"
BLUE = "#1f78b4"

COLUMN_PALETTE = [
    ("#5265d6", "#eef2ff"),
    ("#0f9f8e", "#ecfeff"),
    ("#8a62c4", "#f5f0ff"),
]
ROW_ACCENTS = ["#4c9f68", "#64b67a", "#d89b00", "#e0aa2b", "#e7b85d", "#8a62c4", "#a985d9"]

PAGES = [
    ("part1", ["0-10%", "10-20%", "20-30%", "30-40%"], "fine bins 0-40%"),
    ("part2", ["40-50%", "50-60%", "60-80%"], "fine bins 40-80%"),
]


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


def step_xy(edges: list[float], density: list[float]) -> tuple[list[float], list[float]]:
    x: list[float] = []
    y: list[float] = []
    for idx, value in enumerate(density):
        x.extend([edges[idx], edges[idx + 1]])
        y.extend([value, value])
    return x, y


def fmt_count(value: int) -> str:
    if value >= 1_000_000:
        return f"{value / 1_000_000:.2f}M"
    if value >= 1000:
        return f"{value / 1000:.0f}k"
    return str(value)


def cell_map(payload: dict) -> dict[tuple[str, str], dict]:
    return {(cell["row"], cell["col"]): cell for cell in payload["cells"]}


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
    zorder: int = 5,
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


def draw_top_band(fig: plt.Figure, nodes: list[tuple[str, object, float, str, dict]], page_label: str) -> None:
    add_text(
        fig,
        nodes,
        "title",
        0.040,
        0.940,
        f"Centrality routing score separation, {page_label}",
        size=35.0,
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
        0.872,
        "Same baseline holdout sample in every pad; only evaluator routing changes.",
        size=16.5,
        color=MUTED,
        title_band_exception=True,
    )
    legend_items = [
        (0.042, 0.232, RED, "Signal MC", "embedded Photon12+20"),
        (0.298, 0.350, BLUE, "Inclusive MC", "Jet12+20+30+40, no truth filter"),
    ]
    for idx, (x, w, color, label, detail) in enumerate(legend_items):
        rounded_box(fig, x, 0.804, w, 0.047, face="white", edge=EDGE, lw=0.9)
        fig.patches.append(Rectangle((x + 0.016, 0.823), 0.046, 0.0065, transform=fig.transFigure, facecolor=color, edgecolor=color, zorder=3))
        add_text(fig, nodes, f"legend label {idx}", x + 0.072, 0.831, label, size=13.8, weight="bold")
        add_text(fig, nodes, f"legend detail {idx}", x + 0.072, 0.814, detail, size=11.4, color=MUTED)


def draw_column_header(
    fig: plt.Figure,
    nodes: list[tuple[str, object, float, str, dict]],
    x: float,
    y: float,
    w: float,
    idx: int,
    title: str,
    subtitle: str,
) -> None:
    accent, fill = COLUMN_PALETTE[idx]
    rounded_box(fig, x, y, w, 0.054, face=fill, edge=accent, lw=1.1)
    add_text(fig, nodes, f"column title {idx}", x + w / 2, y + 0.034, title, size=15.8, ha="center", weight="bold")
    add_text(fig, nodes, f"column subtitle {idx}", x + w / 2, y + 0.014, subtitle, size=11.2, ha="center", color=MUTED)


def draw_row_label(
    fig: plt.Figure,
    nodes: list[tuple[str, object, float, str, dict]],
    x: float,
    y: float,
    w: float,
    h: float,
    row: str,
    counts: dict,
    idx: int,
) -> None:
    accent = ROW_ACCENTS[idx % len(ROW_ACCENTS)]
    fig.patches.append(Rectangle((x, y), w, h, transform=fig.transFigure, facecolor="white", edgecolor=EDGE, linewidth=0.9, zorder=1))
    fig.patches.append(Rectangle((x, y), 0.0065, h, transform=fig.transFigure, facecolor=accent, edgecolor=accent, linewidth=0, zorder=2))
    add_text(fig, nodes, f"row label {row}", x + 0.018, y + h * 0.61, row, size=18.5, weight="bold")
    add_text(
        fig,
        nodes,
        f"row counts {row}",
        x + 0.018,
        y + h * 0.34,
        f"S {fmt_count(counts[row]['signal'])}\nI {fmt_count(counts[row]['inclusive'])}",
        size=10.8,
        role="plot_annotation",
        color=MUTED,
        linespacing=0.92,
    )


def draw_metrics_box(ax, cell: dict, *, compact: bool) -> list[object]:
    metrics = cell["metrics"]
    artists: list[object] = []
    box_h = 0.305 if compact else 0.330
    box = Rectangle((0.023, 0.666), 0.615, box_h, transform=ax.transAxes, facecolor="white", edgecolor="#b8c4d6", linewidth=0.85, zorder=5, alpha=0.94)
    ax.add_patch(box)
    label_size = 7.8 if compact else 8.7
    value_size = 8.7 if compact else 9.8
    entries = [
        ("AUC", f"{metrics['auc']:.3f}", "Gap", f"{metrics['median_gap']:.2f}"),
        ("Logloss", f"{metrics['logloss']:.3f}", "Fake rate", f"{100.0 * metrics['wp80_inclusive_fake']:.1f}%"),
    ]
    for ridx, row in enumerate(entries):
        y = 0.886 - 0.145 * ridx if compact else 0.888 - 0.155 * ridx
        for text, x, ha, bold, fs in [
            (row[0], 0.045, "left", False, label_size),
            (row[1], 0.265, "right", True, value_size),
            (row[2], 0.318, "left", False, label_size),
            (row[3], 0.610, "right", True, value_size),
        ]:
            artists.append(
                ax.text(
                    x,
                    y,
                    text,
                    transform=ax.transAxes,
                    ha=ha,
                    va="center",
                    fontsize=fs,
                    family=FONT,
                    fontweight="bold" if bold else "normal",
                    color=INK if bold else MUTED,
                    zorder=6,
                )
            )
    return artists


def draw_axis(ax, cell: dict, edges: list[float], y_max: float, *, show_y: bool, show_x: bool, compact: bool) -> list[object]:
    signal_x, signal_y = step_xy(edges, cell["signal"]["density"])
    inclusive_x, inclusive_y = step_xy(edges, cell["inclusive"]["density"])
    ax.plot(signal_x, signal_y, color=RED, linewidth=1.65)
    ax.fill_between(signal_x, signal_y, step="pre", color=RED, alpha=0.075)
    ax.plot(inclusive_x, inclusive_y, color=BLUE, linewidth=1.65)
    ax.fill_between(inclusive_x, inclusive_y, step="pre", color=BLUE, alpha=0.075)
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, y_max)
    ax.grid(True, color=GRID, linewidth=0.70, alpha=0.92)
    ax.tick_params(axis="both", direction="in", top=True, right=True, labelsize=8.6 if compact else 9.5, length=3.5, width=0.8, colors="#374151")
    ax.yaxis.set_major_locator(MaxNLocator(nbins=3, integer=True))
    if not show_y:
        ax.set_yticklabels([])
    if not show_x:
        ax.set_xticklabels([])
    for spine in ax.spines.values():
        spine.set_linewidth(0.95)
        spine.set_color("#4b5563")
    artists = draw_metrics_box(ax, cell, compact=compact)
    artists.append(
        ax.text(
            0.965,
            0.935,
            f"S {fmt_count(cell['signal']['entries'])}\nI {fmt_count(cell['inclusive']['entries'])}",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=8.4 if compact else 9.4,
            family=FONT,
            fontweight="bold",
            color=MUTED,
            linespacing=0.92,
            zorder=6,
        )
    )
    return artists


def page_readout(payload: dict, rows: list[str]) -> str:
    cells = cell_map(payload)
    fine_minus_global = [cells[(row, "fine")]["metrics"]["auc"] - cells[(row, "global")]["metrics"]["auc"] for row in rows]
    fake_delta = [
        100.0 * (cells[(row, "fine")]["metrics"]["wp80_inclusive_fake"] - cells[(row, "global")]["metrics"]["wp80_inclusive_fake"])
        for row in rows
    ]
    return (
        f"Fine routing AUC shift versus global across these rows: {min(fine_minus_global):+.3f} to {max(fine_minus_global):+.3f}; "
        f"fake-rate shift: {min(fake_delta):+.2f} to {max(fake_delta):+.2f} percentage points."
    )


def render_page(payload: dict, page_key: str, rows: list[str], page_label: str, output_path: Path) -> Path:
    cells = cell_map(payload)
    cols = [
        ("global", "GLOBAL", "1 BDT"),
        ("coarse", "COARSE ROUTING", "3 BDTs"),
        ("fine", "FINE ROUTING", "7 BDTs"),
    ]
    compact = len(rows) >= 4
    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI, facecolor="white")
    nodes: list[tuple[str, object, float, str, dict]] = []

    draw_top_band(fig, nodes, page_label)
    grid_left = 0.178
    grid_right = 0.966
    row_label_x = 0.042
    row_label_w = 0.116
    gap_x = 0.020
    plot_w = (grid_right - grid_left - 2 * gap_x) / 3.0
    header_y = 0.724 if compact else 0.712
    row_top = 0.700 if compact else 0.685
    row_bottom = 0.145
    gap_y = 0.025 if compact else 0.038
    plot_h = (row_top - row_bottom - (len(rows) - 1) * gap_y) / len(rows)

    for c_idx, (_, title, subtitle) in enumerate(cols):
        draw_column_header(fig, nodes, grid_left + c_idx * (plot_w + gap_x), header_y, plot_w, c_idx, title, subtitle)
    add_text(
        fig,
        nodes,
        f"shared y label {page_key}",
        grid_left - 0.022,
        (row_top + row_bottom) / 2,
        "Unit-area density",
        size=12.0,
        role="plot_annotation",
        ha="center",
        va="center",
    )
    nodes[-1][1].set_rotation(90)

    edges = payload["bin_edges"]
    for r_idx, row in enumerate(rows):
        y = row_top - (r_idx + 1) * plot_h - r_idx * gap_y
        global_row_idx = payload["row_order"].index(row)
        draw_row_label(fig, nodes, row_label_x, y, row_label_w, plot_h, row, payload["raw_available_counts_by_row"], global_row_idx)
        row_y_max = max(max(cells[(row, col)]["signal"]["density"] + cells[(row, col)]["inclusive"]["density"]) for col, _, _ in cols)
        row_y_max = max(4.0, row_y_max * 1.28)
        for c_idx, (col, _, _) in enumerate(cols):
            x = grid_left + c_idx * (plot_w + gap_x)
            ax = fig.add_axes([x, y, plot_w, plot_h])
            for a_idx, artist in enumerate(draw_axis(ax, cells[(row, col)], edges, row_y_max, show_y=c_idx == 0, show_x=r_idx == len(rows) - 1, compact=compact)):
                nodes.append((f"plot text {page_key} {row} {col} {a_idx}", artist, 8.2 if compact else 9.2, "plot_annotation", {}))

    rounded_box(fig, 0.042, 0.052, 0.918, 0.066, face="#fff7e6", edge="#e7b85d", lw=1.0)
    add_text(fig, nodes, f"readout title {page_key}", 0.058, 0.086, "Readout", size=16.0, weight="bold")
    add_text(fig, nodes, f"readout body {page_key}", 0.142, 0.086, page_readout(payload, rows), size=12.3, color=INK)

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
                "minimum_audience_font_px": font_px(11.0),
                "minimum_title_font_px": font_px(30.0),
                "minimum_plot_annotation_font_px": font_px(7.6),
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
                "schema": "the8_centrality_routing_fine_rows_score_shapes_slide_v1",
                "png": str(output_path),
                "payload": str(DEFAULT_PAYLOAD),
                "layout_nodes": str(layout_path),
                "page": page_key,
                "rows": rows,
                "campaign_tag": payload.get("campaign_tag"),
                "baseline_model": payload.get("baseline_model"),
                "model_columns": payload.get("model_columns"),
                "selection": payload.get("selection"),
                "labels": payload.get("labels"),
                "holdout_reconstruction": payload.get("holdout_reconstruction"),
                "readout": page_readout(payload, rows),
            },
            indent=2,
        )
        + "\n"
    )
    plt.close(fig)
    return output_path


def render(payload_path: Path = DEFAULT_PAYLOAD) -> list[Path]:
    configure_fonts()
    payload = json.loads(payload_path.read_text())
    outputs = []
    for page_key, rows, page_label in PAGES:
        output = OUT_DIR / f"{OUTPUT_STEM}_{page_key}.png"
        outputs.append(render_page(payload, page_key, rows, page_label, output))
    script_path = OUT_DIR / f"{OUTPUT_STEM}_speaker_script.md"
    script_path.write_text(
        "\n".join(
            [
                "# Centrality routing score-shape slides",
                "",
                "- These candidates use the reconstructed baseline 10% row holdout, not a capped preview sample.",
                "- Every pad uses the same source-class signal and inclusive holdout rows. The columns change only the evaluator routing.",
                "- Signal MC is embedded Photon12+20. Inclusive MC is embedded Jet12+20+30+40 with no truth-background filter.",
                "- Page 1 shows 0-40% fine bins. Page 2 shows 40-80% fine bins.",
                "",
            ]
        )
    )
    return outputs


def main() -> None:
    for output in render():
        print(output)
        print(output.with_suffix(".manifest.json"))
        print(output.with_suffix(".layout_nodes.json"))


if __name__ == "__main__":
    main()
