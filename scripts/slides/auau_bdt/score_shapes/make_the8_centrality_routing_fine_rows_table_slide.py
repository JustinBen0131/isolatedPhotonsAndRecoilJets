#!/usr/bin/env python3
"""Render a fine-centrality table for centrality-routed THE8 BDT diagnostics."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import FancyBboxPatch, Rectangle

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
DEFAULT_OUTPUT = OUT_DIR / "the8_centrality_routing_fine_rows_table_baseline_holdout_v1.png"

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
BLUE = "#1f78b4"
RED = "#d13b35"
EDGE = "#c8d3e2"
GRID = "#e5ebf3"
PAPER = "#ffffff"
READOUT_FILL = "#fff7e6"
READOUT_EDGE = "#e7b85d"

COLUMN_PALETTE = {
    "global": ("#5265d6", "#eef2ff"),
    "coarse": ("#0f9f8e", "#ecfeff"),
    "fine": ("#8a62c4", "#f5f0ff"),
}
ROW_ACCENTS = ["#4c9f68", "#64b67a", "#d89b00", "#e0aa2b", "#e7b85d", "#8a62c4", "#a985d9"]


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
    **node_flags,
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
    nodes.append((name, artist, size, role, node_flags))
    return artist


def rounded_box(
    fig: plt.Figure,
    x: float,
    y: float,
    w: float,
    h: float,
    *,
    face: str,
    edge: str,
    lw: float = 1.0,
    radius: float = 0.006,
    zorder: int = 1,
) -> None:
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


def draw_contract_and_legend(fig: plt.Figure, nodes: list[tuple[str, object, float, str, dict]]) -> None:
    add_text(
        fig,
        nodes,
        "contract",
        0.042,
        0.872,
        "Same baseline holdout sample in every cell; only evaluator routing changes.",
        size=16.4,
        color=MUTED,
        title_band_exception=True,
    )
    legend_y = 0.797
    legend_items = [
        (0.042, 0.230, RED, "Signal MC", "embedded Photon12+20"),
        (0.292, 0.350, BLUE, "Inclusive MC", "Jet12+20+30+40, no truth filter"),
    ]
    for idx, (x, w, color, label, detail) in enumerate(legend_items):
        rounded_box(fig, x, legend_y, w, 0.045, face=PAPER, edge=EDGE, lw=0.9)
        fig.patches.append(
            Rectangle(
                (x + 0.015, legend_y + 0.018),
                0.045,
                0.006,
                transform=fig.transFigure,
                facecolor=color,
                edgecolor=color,
                linewidth=1.0,
                zorder=3,
            )
        )
        add_text(
            fig,
            nodes,
            f"legend label {idx}",
            x + 0.071,
            legend_y + 0.027,
            label,
            size=13.6,
            weight="bold",
            va="center",
        )
        add_text(
            fig,
            nodes,
            f"legend detail {idx}",
            x + 0.071,
            legend_y + 0.011,
            detail,
            size=11.2,
            color=MUTED,
            va="center",
        )


def draw_column_header(
    fig: plt.Figure,
    nodes: list[tuple[str, object, float, str, dict]],
    x: float,
    y: float,
    w: float,
    h: float,
    col: str,
    title: str,
    subtitle: str,
) -> None:
    accent, fill = COLUMN_PALETTE[col]
    rounded_box(fig, x, y, w, h, face=fill, edge=accent, lw=1.15)
    add_text(fig, nodes, f"column title {col}", x + w / 2, y + h * 0.61, title, size=17.4, ha="center", weight="bold")
    add_text(fig, nodes, f"column subtitle {col}", x + w / 2, y + h * 0.27, subtitle, size=12.4, ha="center", color=MUTED)


def metric_pairs(cell: dict) -> list[tuple[str, str]]:
    metrics = cell["metrics"]
    return [
        ("AUC", f"{metrics['auc']:.3f}"),
        ("Logloss", f"{metrics['logloss']:.3f}"),
        ("Gap", f"{metrics['median_gap']:.2f}"),
        ("Fake rate", f"{100.0 * metrics['wp80_inclusive_fake']:.1f}%"),
    ]


def cell_fill(cell: dict, global_cell: dict | None) -> str:
    if global_cell is None:
        return "#f8fafc"
    delta_fake_pp = 100.0 * (cell["metrics"]["wp80_inclusive_fake"] - global_cell["metrics"]["wp80_inclusive_fake"])
    if delta_fake_pp <= -0.35:
        return "#eef8f0"
    if delta_fake_pp >= 0.35:
        return "#fff1f2"
    return "#ffffff"


def draw_metric_cell(
    fig: plt.Figure,
    nodes: list[tuple[str, object, float, str, dict]],
    x: float,
    y: float,
    w: float,
    h: float,
    *,
    row: str,
    col: str,
    cell: dict,
    global_cell: dict | None,
) -> None:
    face = cell_fill(cell, global_cell)
    fig.patches.append(
        Rectangle(
            (x, y),
            w,
            h,
            transform=fig.transFigure,
            facecolor=face,
            edgecolor=EDGE,
            linewidth=0.85,
            zorder=1,
        )
    )
    left_x = x + 0.018
    mid_x = x + w * 0.505
    value_left = x + w * 0.385
    value_right = x + w * 0.885
    ys = [y + h * 0.66, y + h * 0.35]
    metrics = metric_pairs(cell)
    entries = [
        (metrics[0], left_x, value_left, ys[0]),
        (metrics[1], mid_x, value_right, ys[0]),
        (metrics[2], left_x, value_left, ys[1]),
        (metrics[3], mid_x, value_right, ys[1]),
    ]
    for idx, ((label, value), lx, vx, yy) in enumerate(entries):
        add_text(
            fig,
            nodes,
            f"metric label {row} {col} {idx}",
            lx,
            yy,
            label,
            size=11.6,
            role="plot_annotation",
            color=MUTED,
        )
        add_text(
            fig,
            nodes,
            f"metric value {row} {col} {idx}",
            vx,
            yy,
            value,
            size=13.6,
            role="plot_annotation",
            ha="right",
            weight="bold",
            color=INK,
        )


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
    fig.patches.append(Rectangle((x, y), w, h, transform=fig.transFigure, facecolor="#ffffff", edgecolor=EDGE, linewidth=0.85, zorder=1))
    fig.patches.append(Rectangle((x, y), 0.007, h, transform=fig.transFigure, facecolor=accent, edgecolor=accent, linewidth=0, zorder=2))
    add_text(fig, nodes, f"row centrality {row}", x + 0.020, y + h * 0.63, row, size=17.8, weight="bold")
    add_text(
        fig,
        nodes,
        f"row counts {row}",
        x + 0.020,
        y + h * 0.34,
        f"S {fmt_count(counts[row]['signal'])}    I {fmt_count(counts[row]['inclusive'])}",
        size=11.8,
        role="plot_annotation",
        color=MUTED,
    )


def readout_text(payload: dict) -> str:
    cells = cell_map(payload)
    auc_deltas: list[float] = []
    fake_deltas: list[float] = []
    best_rows: list[str] = []
    for row in payload["row_order"]:
        global_metrics = cells[(row, "global")]["metrics"]
        best = max(("global", "coarse", "fine"), key=lambda col: cells[(row, col)]["metrics"]["auc"])
        best_rows.append(best)
        for col in ("coarse", "fine"):
            metrics = cells[(row, col)]["metrics"]
            auc_deltas.append(metrics["auc"] - global_metrics["auc"])
            fake_deltas.append(100.0 * (metrics["wp80_inclusive_fake"] - global_metrics["wp80_inclusive_fake"]))
    fine_wins = best_rows.count("fine")
    coarse_wins = best_rows.count("coarse")
    global_wins = best_rows.count("global")
    return (
        f"Routing changes AUC by {min(auc_deltas):+.3f} to {max(auc_deltas):+.3f} and fake rate by "
        f"{min(fake_deltas):+.2f} to {max(fake_deltas):+.2f} percentage points versus the global BDT.\n"
        f"Best AUC by fine bin, global {global_wins}, coarse {coarse_wins}, fine {fine_wins}."
    )


def render(payload_path: Path, output_path: Path) -> Path:
    configure_fonts()
    payload = json.loads(payload_path.read_text())
    cells = cell_map(payload)
    rows = payload["row_order"]
    cols = [
        ("global", "GLOBAL", "1 BDT"),
        ("coarse", "COARSE ROUTING", "3 BDTs"),
        ("fine", "FINE ROUTING", "7 BDTs"),
    ]

    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI, facecolor="white")
    nodes: list[tuple[str, object, float, str, dict]] = []

    add_text(
        fig,
        nodes,
        "title",
        0.040,
        0.940,
        "Centrality routing across fine holdout bins",
        size=39.0,
        role="title",
        va="top",
        weight="bold",
        title_anchor=True,
    )
    draw_contract_and_legend(fig, nodes)

    grid_x = 0.042
    row_label_w = 0.142
    gap_x = 0.012
    col_w = (0.960 - grid_x - row_label_w - 3 * gap_x) / 3.0
    header_y = 0.708
    header_h = 0.056
    row_top = 0.690
    row_bottom = 0.147
    gap_y = 0.006
    row_h = (row_top - row_bottom - (len(rows) - 1) * gap_y) / len(rows)

    add_text(fig, nodes, "row header", grid_x + 0.020, header_y + header_h * 0.52, "Fine bin", size=14.4, weight="bold", color=MUTED)
    for c_idx, (col, title, subtitle) in enumerate(cols):
        x = grid_x + row_label_w + gap_x + c_idx * (col_w + gap_x)
        draw_column_header(fig, nodes, x, header_y, col_w, header_h, col, title, subtitle)

    for r_idx, row in enumerate(rows):
        y = row_top - (r_idx + 1) * row_h - r_idx * gap_y
        draw_row_label(fig, nodes, grid_x, y, row_label_w, row_h, row, payload["raw_available_counts_by_row"], r_idx)
        global_cell = cells[(row, "global")]
        for c_idx, (col, _, _) in enumerate(cols):
            x = grid_x + row_label_w + gap_x + c_idx * (col_w + gap_x)
            draw_metric_cell(
                fig,
                nodes,
                x,
                y,
                col_w,
                row_h,
                row=row,
                col=col,
                cell=cells[(row, col)],
                global_cell=None if col == "global" else global_cell,
            )

    rounded_box(fig, 0.042, 0.055, 0.918, 0.070, face=READOUT_FILL, edge=READOUT_EDGE, lw=1.0)
    add_text(fig, nodes, "readout title", 0.057, 0.091, "Readout", size=17.5, weight="bold")
    add_text(fig, nodes, "readout body", 0.143, 0.090, readout_text(payload), size=12.8, color=INK, linespacing=1.06)

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
                "minimum_plot_annotation_font_px": font_px(10.5),
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
                "schema": "the8_centrality_routing_fine_rows_table_slide_v1",
                "png": str(output_path),
                "payload": str(payload_path),
                "layout_nodes": str(layout_path),
                "campaign_tag": payload.get("campaign_tag"),
                "baseline_model": payload.get("baseline_model"),
                "model_columns": payload.get("model_columns"),
                "selection": payload.get("selection"),
                "labels": payload.get("labels"),
                "holdout_reconstruction": payload.get("holdout_reconstruction"),
                "raw_available_counts_by_row": payload.get("raw_available_counts_by_row"),
                "readout": readout_text(payload),
            },
            indent=2,
        )
        + "\n"
    )

    script_path = output_path.with_name(output_path.stem + "_speaker_script.md")
    script_path.write_text(
        "\n".join(
            [
                "# Centrality routing across fine holdout bins",
                "",
                "- This is the table companion to the three-row score-shape slide.",
                "- Every cell uses the same reconstructed baseline 10% row holdout, with no capped preview sample.",
                "- The three columns change only the evaluator: one global 14-feature BDT, three coarse centrality BDTs, or seven fine centrality BDTs.",
                "- Signal is source-class embedded Photon12+20. Inclusive is source-class embedded Jet12+20+30+40 with no truth-background filter.",
                f"- {readout_text(payload)}",
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
