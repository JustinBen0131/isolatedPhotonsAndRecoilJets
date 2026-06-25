#!/usr/bin/env python3
"""Render an ET-routing table slide for corrected THE8 BDT diagnostics."""

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
    / "et_routing_20260622"
)
DEFAULT_PAYLOAD = OUT_DIR / "the8_et_routing_score_payload_baseline_holdout.json"
DEFAULT_OUTPUT = OUT_DIR / "the8_et_routing_table_baseline_holdout_v1.png"

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
EDGE = "#c8d3e2"
RED = "#d13b35"
BLUE = "#1f78b4"
GREEN = "#15803d"
GOLD = "#d89b00"
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


def metric_lines(cell: dict) -> list[tuple[str, str]]:
    m = cell["metrics"]
    return [
        ("AUC", f"{m['auc']:.3f}"),
        ("Logloss", f"{m['logloss']:.3f}"),
        ("Gap", f"{m['median_gap']:.2f}"),
        ("Fake rate", f"{100.0 * m['wp80_inclusive_fake']:.1f}%"),
    ]


def delta_payload(global_cell: dict, et_cell: dict) -> dict[str, float]:
    g = global_cell["metrics"]
    e = et_cell["metrics"]
    return {
        "auc": e["auc"] - g["auc"],
        "logloss": e["logloss"] - g["logloss"],
        "gap": e["median_gap"] - g["median_gap"],
        "fake_pp": 100.0 * (e["wp80_inclusive_fake"] - g["wp80_inclusive_fake"]),
    }


def delta_face(delta: dict[str, float]) -> str:
    if delta["auc"] > 0.0 and delta["logloss"] < 0.0 and delta["fake_pp"] < 0.0:
        if delta["fake_pp"] <= -0.50:
            return "#eaf7ed"
        return "#f3fbf5"
    return "#fff7ed"


def draw_top(fig: plt.Figure, nodes: list[tuple[str, object, float, str, dict]]) -> None:
    add_text(
        fig,
        nodes,
        "title",
        0.040,
        0.940,
        "ET-dependent BDT routing across the fixed holdout",
        size=38.0,
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
        0.873,
        "Fixed baseline holdout; rows split only by cluster ET.",
        size=16.2,
        color=MUTED,
        title_band_exception=True,
    )
    add_text(
        fig,
        nodes,
        "metric note",
        0.042,
        0.846,
        "Fake rate is inclusive MC above the row WP80 threshold.",
        size=13.2,
        color=MUTED,
        title_band_exception=True,
    )
    legend_items = [
        (0.502, 0.210, RED, "Signal MC", "embedded Photon12+20"),
        (0.732, 0.246, BLUE, "Inclusive MC", "Jet12+20+30+40, no truth filter"),
    ]
    for idx, (x, w, color, label, detail) in enumerate(legend_items):
        rounded_box(fig, x, 0.842, w, 0.050, face="white", edge=EDGE, lw=0.9)
        fig.patches.append(Rectangle((x + 0.014, 0.862), 0.040, 0.006, transform=fig.transFigure, facecolor=color, edgecolor=color, zorder=3))
        add_text(fig, nodes, f"legend label {idx}", x + 0.064, 0.870, label, size=12.8, weight="bold", title_band_exception=True)
        add_text(fig, nodes, f"legend detail {idx}", x + 0.064, 0.852, detail, size=10.4, color=MUTED, title_band_exception=True)


def draw_header(
    fig: plt.Figure,
    nodes: list[tuple[str, object, float, str, dict]],
    x: float,
    y: float,
    w: float,
    title: str,
    subtitle: str,
    *,
    edge: str,
    fill: str,
) -> None:
    rounded_box(fig, x, y, w, 0.058, face=fill, edge=edge, lw=1.15)
    add_text(fig, nodes, f"header {title}", x + w / 2, y + 0.036, title, size=16.0, ha="center", weight="bold")
    add_text(fig, nodes, f"header sub {title}", x + w / 2, y + 0.015, subtitle, size=10.8, ha="center", color=MUTED)


def draw_metric_cell(
    fig: plt.Figure,
    nodes: list[tuple[str, object, float, str, dict]],
    x: float,
    y: float,
    w: float,
    h: float,
    cell: dict,
    name: str,
    *,
    face: str = "#ffffff",
) -> None:
    fig.patches.append(Rectangle((x, y), w, h, transform=fig.transFigure, facecolor=face, edgecolor=EDGE, linewidth=0.85, zorder=1))
    pairs = metric_lines(cell)
    positions = [
        (pairs[0], x + 0.015, x + w * 0.390, y + h * 0.67),
        (pairs[1], x + w * 0.505, x + w * 0.885, y + h * 0.67),
        (pairs[2], x + 0.015, x + w * 0.390, y + h * 0.35),
        (pairs[3], x + w * 0.505, x + w * 0.885, y + h * 0.35),
    ]
    for idx, ((label, value), lx, vx, yy) in enumerate(positions):
        add_text(fig, nodes, f"{name} label {idx}", lx, yy, label, size=10.6, role="plot_annotation", color=MUTED)
        add_text(fig, nodes, f"{name} value {idx}", vx, yy, value, size=12.8, role="plot_annotation", ha="right", weight="bold")


def draw_delta_cell(
    fig: plt.Figure,
    nodes: list[tuple[str, object, float, str, dict]],
    x: float,
    y: float,
    w: float,
    h: float,
    delta: dict[str, float],
    name: str,
) -> None:
    face = delta_face(delta)
    fig.patches.append(Rectangle((x, y), w, h, transform=fig.transFigure, facecolor=face, edgecolor=EDGE, linewidth=0.85, zorder=1))
    win = delta["auc"] > 0 and delta["logloss"] < 0 and delta["fake_pp"] < 0
    verdict = "ET better" if win else "mixed"
    verdict_color = GREEN if win else GOLD
    add_text(fig, nodes, f"{name} verdict", x + 0.015, y + h * 0.72, verdict, size=12.8, role="plot_annotation", weight="bold", color=verdict_color)
    add_text(fig, nodes, f"{name} dauc", x + 0.015, y + h * 0.44, f"AUC {delta['auc']:+.3f}", size=10.9, role="plot_annotation", color=INK)
    add_text(fig, nodes, f"{name} dfake", x + w * 0.520, y + h * 0.44, f"Fake {delta['fake_pp']:+.2f} pp", size=10.9, role="plot_annotation", color=INK)
    add_text(fig, nodes, f"{name} dll", x + 0.015, y + h * 0.22, f"Logloss {delta['logloss']:+.3f}", size=10.1, role="plot_annotation", color=MUTED)
    add_text(fig, nodes, f"{name} dgap", x + w * 0.520, y + h * 0.22, f"Gap {delta['gap']:+.2f}", size=10.1, role="plot_annotation", color=MUTED)


def draw_et_label(
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
    accent = "#4c9f68" if idx < 3 else "#d89b00" if idx < 6 else "#8a62c4"
    fig.patches.append(Rectangle((x, y), w, h, transform=fig.transFigure, facecolor="#ffffff", edgecolor=EDGE, linewidth=0.85, zorder=1))
    fig.patches.append(Rectangle((x, y), 0.0065, h, transform=fig.transFigure, facecolor=accent, edgecolor=accent, linewidth=0, zorder=2))
    add_text(fig, nodes, f"row label {row}", x + 0.018, y + h * 0.64, row, size=15.5, weight="bold")
    add_text(
        fig,
        nodes,
        f"row counts {row}",
        x + 0.018,
        y + h * 0.34,
        f"S {fmt_count(counts[row]['signal'])}    I {fmt_count(counts[row]['inclusive'])}",
        size=10.5,
        role="plot_annotation",
        color=MUTED,
    )


def readout_text(payload: dict) -> str:
    cells = cell_map(payload)
    deltas = [delta_payload(cells[(row, "global")], cells[(row, "et")]) for row in payload["row_order"]]
    auc = [d["auc"] for d in deltas]
    fake = [d["fake_pp"] for d in deltas]
    wins = sum(1 for d in deltas if d["auc"] > 0 and d["logloss"] < 0 and d["fake_pp"] < 0)
    best_idx = min(range(len(fake)), key=lambda i: fake[i])
    return (
        f"ET routing improves AUC, logloss, and fake rate in {wins}/8 ET bins. "
        f"AUC gain spans {min(auc):+.3f} to {max(auc):+.3f};\n"
        f"fake-rate change spans {min(fake):+.2f} to {max(fake):+.2f} percentage points. "
        f"Largest fake-rate drop is {payload['row_order'][best_idx]}."
    )


def render(payload_path: Path, output_path: Path) -> Path:
    configure_fonts()
    payload = json.loads(payload_path.read_text())
    cells = cell_map(payload)
    rows = payload["row_order"]
    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI, facecolor="white")
    nodes: list[tuple[str, object, float, str, dict]] = []
    draw_top(fig, nodes)

    grid_x = 0.042
    et_w = 0.132
    gap_x = 0.012
    global_w = 0.250
    etroute_w = 0.250
    delta_w = 0.284
    header_y = 0.760
    row_top = 0.742
    row_bottom = 0.148
    gap_y = 0.006
    row_h = (row_top - row_bottom - (len(rows) - 1) * gap_y) / len(rows)

    add_text(fig, nodes, "row header", grid_x + 0.018, header_y + 0.032, "ET bin", size=13.6, weight="bold", color=MUTED)
    x_global = grid_x + et_w + gap_x
    x_et = x_global + global_w + gap_x
    x_delta = x_et + etroute_w + gap_x
    draw_header(fig, nodes, x_global, header_y, global_w, "GLOBAL", "1 BDT over 15-35 GeV", edge="#5265d6", fill="#eef2ff")
    draw_header(fig, nodes, x_et, header_y, etroute_w, "ET ROUTING", "8 BDTs by ET bin", edge="#0f9f8e", fill="#ecfeff")
    draw_header(fig, nodes, x_delta, header_y, delta_w, "ET ROUTING CHANGE", "ET-routed minus global", edge="#15803d", fill="#eef8f0")

    counts = payload["raw_available_counts_by_row"]
    for idx, row in enumerate(rows):
        y = row_top - (idx + 1) * row_h - idx * gap_y
        draw_et_label(fig, nodes, grid_x, y, et_w, row_h, row, counts, idx)
        global_cell = cells[(row, "global")]
        et_cell = cells[(row, "et")]
        draw_metric_cell(fig, nodes, x_global, y, global_w, row_h, global_cell, f"{row} global", face="#f8fafc")
        draw_metric_cell(fig, nodes, x_et, y, etroute_w, row_h, et_cell, f"{row} et")
        draw_delta_cell(fig, nodes, x_delta, y, delta_w, row_h, delta_payload(global_cell, et_cell), f"{row} delta")

    rounded_box(fig, 0.042, 0.046, 0.918, 0.074, face=READOUT_FILL, edge=READOUT_EDGE, lw=1.0)
    add_text(fig, nodes, "readout title", 0.058, 0.085, "Readout", size=16.2, weight="bold")
    add_text(fig, nodes, "readout body", 0.142, 0.084, readout_text(payload), size=10.8, color=INK, linespacing=1.05)

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
                "minimum_plot_annotation_font_px": font_px(9.6),
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
                "schema": "the8_et_routing_table_slide_v1",
                "png": str(output_path),
                "payload": str(payload_path),
                "layout_nodes": str(layout_path),
                "campaign_tag": payload.get("campaign_tag"),
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
                "# ET-dependent BDT routing across the fixed holdout",
                "",
                "- This slide compares the global 14-feature BDT with the ET-only routed BDT family, not centrality routing.",
                "- Every row uses the same reconstructed baseline 10% row holdout, split by cluster ET.",
                "- Signal is source-class embedded Photon12+20. Inclusive is source-class embedded Jet12+20+30+40 with no truth-background filter.",
                "- Fake rate is the inclusive MC fraction above the row threshold set at 80% signal efficiency.",
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
