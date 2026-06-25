#!/usr/bin/env python3
"""Render a centrality-routing BDT score-shape comparison slide."""

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
DEFAULT_PAYLOAD = OUT_DIR / "the8_centrality_routing_score_payload_baseline_holdout.json"
DEFAULT_OUTPUT = OUT_DIR / "the8_centrality_routing_ladder_baseline_holdout_v1.png"

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
PANEL = "#f8fafc"
EDGE = "#c7d2e1"
RED = "#d13b35"
BLUE = "#1f78b4"
GOLD = "#d89b00"

ROW_PALETTE = [
    ("#4c9f68", "#eef8f0"),
    ("#d89b00", "#fff4cf"),
    ("#8a62c4", "#f1eafa"),
]
COLUMN_PALETTE = [
    ("#5265d6", "#eef2ff"),
    ("#0f9f8e", "#ecfeff"),
    ("#8a62c4", "#f5f0ff"),
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


def draw_legend(fig: plt.Figure, text_nodes: list[tuple[str, object, float, str]]) -> None:
    y0, h = 0.802, 0.058
    cards = [
        (0.290, 0.312, RED, "Signal MC\nembedded Photon12+20", 14.6),
        (0.625, 0.335, BLUE, "Inclusive MC\nJet12+20+30+40, no truth filter", 13.9),
    ]
    for i, (x0, w, color, label, size) in enumerate(cards):
        fig.patches.append(
            FancyBboxPatch(
                (x0, y0),
                w,
                h,
                transform=fig.transFigure,
                boxstyle="round,pad=0.006,rounding_size=0.006",
                facecolor="white",
                edgecolor="#c9d5e6",
                linewidth=1.0,
                zorder=1,
            )
        )
        fig.patches.append(
            Rectangle(
                (x0 + 0.020, y0 + 0.019),
                0.050,
                0.007,
                transform=fig.transFigure,
                facecolor=color,
                edgecolor=color,
                linewidth=1.0,
                zorder=2,
            )
        )
        text = fig.text(
            x0 + 0.085,
            y0 + h / 2,
            label,
            ha="left",
            va="center",
            fontsize=size,
            family=FONT,
            fontweight="bold",
            color=INK,
            linespacing=0.82,
            zorder=2,
        )
        text_nodes.append((f"legend {i}", text, size, "audience"))


def draw_column_label(
    fig: plt.Figure,
    x: float,
    w: float,
    idx: int,
    title: str,
    subtitle: str,
    text_nodes: list[tuple[str, object, float, str]],
) -> None:
    accent, fill = COLUMN_PALETTE[idx]
    y, h = 0.726, 0.068
    fig.patches.append(
        FancyBboxPatch(
            (x, y),
            w,
            h,
            transform=fig.transFigure,
            boxstyle="round,pad=0.006,rounding_size=0.006",
            facecolor=fill,
            edgecolor=accent,
            linewidth=1.1,
            zorder=1,
        )
    )
    t1 = fig.text(
        x + w / 2,
        y + h * 0.60,
        title,
        ha="center",
        va="center",
        fontsize=17.2,
        family=FONT,
        fontweight="bold",
        color=INK,
        zorder=2,
    )
    t2 = fig.text(
        x + w / 2,
        y + h * 0.25,
        subtitle,
        ha="center",
        va="center",
        fontsize=12.5,
        family=FONT,
        color=MUTED,
        zorder=2,
    )
    text_nodes.extend([(f"column title {idx}", t1, 17.2, "audience"), (f"column subtitle {idx}", t2, 12.5, "audience")])


def draw_row_label(
    fig: plt.Figure,
    y: float,
    h: float,
    idx: int,
    label: str,
    counts: dict,
    text_nodes: list[tuple[str, object, float, str]],
) -> None:
    accent, fill = ROW_PALETTE[idx]
    x, w = 0.036, 0.119
    fig.patches.append(
        Rectangle(
            (x, y),
            w,
            h,
            transform=fig.transFigure,
            facecolor=fill,
            edgecolor=accent,
            linewidth=1.1,
            zorder=1,
        )
    )
    fig.patches.append(Rectangle((x + 0.017, y + 0.18 * h), 0.006, 0.64 * h, transform=fig.transFigure, facecolor=accent, edgecolor=accent, zorder=2))
    t1 = fig.text(
        x + 0.067,
        y + 0.63 * h,
        label,
        ha="center",
        va="center",
        fontsize=24.0,
        family=FONT,
        fontweight="bold",
        color=INK,
        zorder=2,
    )
    t2 = fig.text(
        x + 0.067,
        y + 0.35 * h,
        "baseline\nholdout",
        ha="center",
        va="center",
        fontsize=12.0,
        family=FONT,
        color=INK,
        linespacing=0.9,
        zorder=2,
    )
    t3 = fig.text(
        x + 0.067,
        y + 0.14 * h,
        f"S {fmt_count(counts[label]['signal'])}  I {fmt_count(counts[label]['inclusive'])}",
        ha="center",
        va="center",
        fontsize=9.2,
        family=FONT,
        color=MUTED,
        zorder=2,
    )
    text_nodes.extend([(f"row label {idx}", t1, 24.0, "audience"), (f"row sublabel {idx}", t2, 12.0, "audience"), (f"row count {idx}", t3, 9.2, "plot_annotation")])


def draw_metrics(ax, metrics: dict, *, fontsize: float = 11.0) -> list:
    artists = []
    box = Rectangle((0.025, 0.640), 0.785, 0.325, transform=ax.transAxes, facecolor="white", edgecolor="#b8c4d6", linewidth=0.9, zorder=5)
    ax.add_patch(box)
    rows = [
        ("AUC", f"{metrics['auc']:.3f}", "Gap", f"{metrics['median_gap']:.2f}"),
        ("Logloss", f"{metrics['logloss']:.3f}", "Fake rate", f"{100.0 * metrics['wp80_inclusive_fake']:.1f}%"),
    ]
    for idx, (l1, v1, l2, v2) in enumerate(rows):
        yy = 0.875 - idx * 0.165
        for text, xx, ha, bold, fs in [
            (l1, 0.045, "left", False, fontsize - 1.5),
            (v1, 0.330, "right", True, fontsize),
            (l2, 0.430, "left", False, fontsize - 2.0),
            (v2, 0.785, "right", True, fontsize),
        ]:
            artists.append(
                ax.text(
                    xx,
                    yy,
                    text,
                    transform=ax.transAxes,
                    ha=ha,
                    va="center",
                    fontsize=fs,
                    family=FONT,
                    fontweight="bold" if bold else "normal",
                    color=INK,
                    zorder=6,
                )
            )
    return artists


def draw_axis(ax, cell: dict, edges: list[float], y_max: float, show_y: bool, show_x: bool) -> list:
    signal_x, signal_y = step_xy(edges, cell["signal"]["density"])
    incl_x, incl_y = step_xy(edges, cell["inclusive"]["density"])
    ax.plot(signal_x, signal_y, color=RED, linewidth=1.9)
    ax.fill_between(signal_x, signal_y, step="pre", color=RED, alpha=0.08)
    ax.plot(incl_x, incl_y, color=BLUE, linewidth=1.9)
    ax.fill_between(incl_x, incl_y, step="pre", color=BLUE, alpha=0.08)

    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, y_max)
    ax.grid(True, color=GRID, linewidth=0.75, alpha=0.9)
    ax.tick_params(axis="both", direction="in", top=True, right=True, labelsize=10.2, length=4, width=0.8, colors="#374151")
    ax.yaxis.set_major_locator(MaxNLocator(nbins=3, integer=True))
    if not show_y:
        ax.set_yticklabels([])
    else:
        ax.set_ylabel("Unit-area\ndensity", fontsize=13.2, family=FONT, color=INK, labelpad=2)
    if show_x:
        ax.set_xlabel("BDT score", fontsize=13.4, family=FONT, color=INK, labelpad=2)
    else:
        ax.set_xticklabels([])
    for spine in ax.spines.values():
        spine.set_linewidth(1.0)
        spine.set_color("#4b5563")

    artists = draw_metrics(ax, cell["metrics"])
    artists.append(
        ax.text(
            0.965,
            0.940,
            f"S {fmt_count(cell['signal']['entries'])}\nI {fmt_count(cell['inclusive']['entries'])}",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=10.7,
            family=FONT,
            fontweight="bold",
            color=MUTED,
            linespacing=0.95,
            zorder=6,
        )
    )
    return artists


def metric_delta_text(payload: dict) -> str:
    rows = ["0-20%", "20-50%", "50-80%"]
    cells = cell_map(payload)
    fine_auc = []
    fine_wp = []
    coarse_auc = []
    coarse_wp = []
    for row in rows:
        g = cells[(row, "global")]["metrics"]
        c = cells[(row, "coarse")]["metrics"]
        f = cells[(row, "fine")]["metrics"]
        coarse_auc.append(c["auc"] - g["auc"])
        fine_auc.append(f["auc"] - g["auc"])
        coarse_wp.append(100.0 * (c["wp80_inclusive_fake"] - g["wp80_inclusive_fake"]))
        fine_wp.append(100.0 * (f["wp80_inclusive_fake"] - g["wp80_inclusive_fake"]))
    all_auc = coarse_auc + fine_auc
    all_wp = coarse_wp + fine_wp
    return (
        "Same baseline holdout sample in every pad. Only the evaluator changes: 1 global BDT, 3 coarse-centrality BDTs, or 7 fine-centrality BDTs.\n"
        f"Routing changes AUC by {min(all_auc):+.3f} to {max(all_auc):+.3f}; "
        f"fake rate moves by {min(all_wp):+.2f} to {max(all_wp):+.2f} percentage points."
    )


def render(payload_path: Path, output_path: Path) -> Path:
    configure_fonts()
    payload = json.loads(payload_path.read_text())
    cells = cell_map(payload)
    rows = ["0-20%", "20-50%", "50-80%"]
    cols = [
        ("global", "GLOBAL", "one 14-feature BDT"),
        ("coarse", "COARSE ROUTING", "3 centrality BDTs"),
        ("fine", "FINE ROUTING", "7 centrality BDTs"),
    ]
    edges = payload["bin_edges"]
    y_max = max(max(cell["signal"]["density"] + cell["inclusive"]["density"]) for cell in payload["cells"]) * 1.20
    y_max = max(8.0, y_max)

    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI, facecolor="white")
    text_nodes: list[tuple[str, object, float, str]] = []
    title = fig.text(
        0.040,
        0.936,
        "Centrality routing: is one 14-feature BDT enough?",
        ha="left",
        va="top",
        fontsize=35.0,
        family=FONT,
        fontweight="bold",
        color=INK,
    )
    text_nodes.append(("title", title, 35.0, "title"))
    draw_legend(fig, text_nodes)

    grid_left = 0.190
    grid_right = 0.970
    gap_x = 0.026
    plot_w = (grid_right - grid_left - 2 * gap_x) / 3.0
    row_top = 0.704
    row_bottom = 0.170
    gap_y = 0.050
    plot_h = (row_top - row_bottom - 2 * gap_y) / 3.0

    for idx, (_, title_text, subtitle) in enumerate(cols):
        draw_column_label(fig, grid_left + idx * (plot_w + gap_x), plot_w, idx, title_text, subtitle, text_nodes)

    for r_idx, row in enumerate(rows):
        y = row_top - (r_idx + 1) * plot_h - r_idx * gap_y
        draw_row_label(fig, y, plot_h, r_idx, row, payload["raw_available_counts_by_row"], text_nodes)
        for c_idx, (col, _, _) in enumerate(cols):
            x = grid_left + c_idx * (plot_w + gap_x)
            ax = fig.add_axes([x, y, plot_w, plot_h])
            for a_idx, artist in enumerate(draw_axis(ax, cells[(row, col)], edges, y_max, show_y=c_idx == 0, show_x=r_idx == 2)):
                text_nodes.append((f"plot text {row} {col} {a_idx}", artist, 10.5, "plot_annotation"))

    readout_box = FancyBboxPatch(
        (0.040, 0.052),
        0.930,
        0.078,
        transform=fig.transFigure,
        boxstyle="round,pad=0.010,rounding_size=0.006",
        facecolor="#fff7e6",
        edgecolor="#e7b85d",
        linewidth=1.0,
        zorder=1,
    )
    fig.patches.append(readout_box)
    ro_title = fig.text(
        0.055,
        0.104,
        "Readout",
        ha="left",
        va="center",
        fontsize=17.0,
        family=FONT,
        fontweight="bold",
        color=INK,
        zorder=2,
    )
    ro = fig.text(
        0.145,
        0.092,
        metric_delta_text(payload),
        ha="left",
        va="center",
        fontsize=12.6,
        family=FONT,
        color=INK,
        zorder=2,
        linespacing=1.08,
    )
    text_nodes.extend([("readout title", ro_title, 17.0, "audience"), ("readout body", ro, 13.0, "audience")])

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=SLIDE_DPI)

    layout_path = output_path.with_suffix(".layout_nodes.json")
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
                "colon_style_exception": name.startswith("legend"),
            }
        )
    layout_path.write_text(
        json.dumps(
            {
                "schema": "slide_layout_nodes_v1",
                "title_axis_x": 0.040 * 2560,
                "minimum_audience_font_px": font_px(12.0),
                "minimum_title_font_px": font_px(26.0),
                "minimum_plot_annotation_font_px": font_px(9.0),
                "nodes": nodes,
            },
            indent=2,
        )
        + "\n"
    )

    manifest_path = output_path.with_suffix(".manifest.json")
    manifest_path.write_text(
        json.dumps(
            {
                "schema": "the8_centrality_routing_ladder_slide_v1",
                "png": str(output_path),
                "payload": str(payload_path),
                "layout_nodes": str(layout_path),
                "campaign_tag": payload.get("campaign_tag"),
                "model_columns": payload.get("model_columns"),
                "selection": payload.get("selection"),
                "labels": payload.get("labels"),
                "holdout_reconstruction": payload.get("holdout_reconstruction"),
                "raw_available_counts_by_row": payload.get("raw_available_counts_by_row"),
                "kept_counts_by_row": payload.get("kept_counts_by_row"),
                "readout": metric_delta_text(payload),
            },
            indent=2,
        )
        + "\n"
    )

    script_path = output_path.with_name(output_path.stem + "_speaker_script.md")
    script_path.write_text(
        "\n".join(
            [
                "# Centrality-routed BDT score-shape comparison",
                "",
                "- This candidate uses the reconstructed baseline 10% row holdout, not a capped scored-cache preview.",
                "- Every column uses the same corrected 14-feature input family; only the model routing changes: one global model, three coarse centrality models, or seven fine centrality models.",
                "- The rows are broad centrality slices within the fixed baseline holdout. Signal is embedded Photon12+20 source-class MC, and inclusive is embedded Jet12+20+30+40 source-class MC with no truth-background filter.",
                f"- {metric_delta_text(payload)}",
                "- The decision point for the final slide is whether coarse or fine centrality routing gives a meaningful enough fake-rate or median-gap gain over the global model to justify the extra routing complexity.",
                "",
            ]
        )
    )
    plt.close(fig)
    return output_path


def main() -> None:
    out = render(DEFAULT_PAYLOAD, DEFAULT_OUTPUT)
    print(out)
    print(out.with_suffix(".manifest.json"))
    print(out.with_suffix(".layout_nodes.json"))
    print(out.with_name(out.stem + "_speaker_script.md"))


if __name__ == "__main__":
    main()
