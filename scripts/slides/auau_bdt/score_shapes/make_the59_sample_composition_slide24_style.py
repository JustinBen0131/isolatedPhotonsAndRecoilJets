#!/usr/bin/env python3
"""Render the THE-59 sample-composition controls slide-24-style PNG."""

from __future__ import annotations

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


OUT_DIR = Path(
    "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/"
    "auauTightBDTValidation/THE59_model_comparison_remakes_20260620"
)
METRICS_JSON = OUT_DIR / "the59_sample_composition_sourceclass_metrics.json"
TIMES_DIR = Path("/System/Library/Fonts/Supplemental")
TIMES_FONTS = [
    TIMES_DIR / "Times New Roman.ttf",
    TIMES_DIR / "Times New Roman Bold.ttf",
    TIMES_DIR / "Times New Roman Italic.ttf",
    TIMES_DIR / "Times New Roman Bold Italic.ttf",
]
FONT_FAMILY = "Times New Roman"

ROW_ORDER = ["sample_jet12_20", "sample_jet12_20_30", "baseline"]
COL_ORDER = ["sample_jet12_20", "sample_jet12_20_30", "baseline"]
LABELS = {
    "sample_jet12_20": "Jet12+20",
    "sample_jet12_20_30": "Jet12+20\n+30",
    "baseline": "Jet12+20\n+30+40",
}
COL_LABELS = {
    "sample_jet12_20": "Jet12+20",
    "sample_jet12_20_30": "Jet12+20+30",
    "baseline": "Jet12+20+30+40",
}
ROW_COLORS = {
    "sample_jet12_20": ("#4c9f68", "#eef8f0"),
    "sample_jet12_20_30": ("#d89b00", "#fff4cf"),
    "baseline": ("#8a62c4", "#f1eafa"),
}
COL_COLORS = {
    "sample_jet12_20": ("#4c9f68", "#eef8f0"),
    "sample_jet12_20_30": ("#d89b00", "#fff4cf"),
    "baseline": ("#8a62c4", "#f1eafa"),
}


def _configure_fonts() -> None:
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


def _font_px(points: float) -> float:
    return points * SLIDE_DPI / 72.0


def _bbox_px(fig: plt.Figure, artist) -> list[float]:
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


def _load_cells() -> dict[tuple[str, str], dict]:
    payload = json.loads(METRICS_JSON.read_text())
    cells: dict[tuple[str, str], dict] = {}
    for sample in payload["common_samples"]:
        row_key = sample["key"]
        for branch in sample["branches"]:
            cells[(row_key, branch["key"])] = branch
    return cells


def _step_xy(branch: dict, class_name: str) -> tuple[list[float], list[float]]:
    edges = branch["bin_edges"]
    density = branch[class_name]["density"]
    x = []
    y = []
    for idx, value in enumerate(density):
        x.extend([edges[idx], edges[idx + 1]])
        y.extend([value, value])
    return x, y


def _format_metrics(metrics: dict) -> list[tuple[str, str]]:
    return [
        ("AUC", f"{metrics['auc']:.3f}"),
        ("Logloss", f"{metrics['logloss']:.3f}"),
        ("Median gap", f"{metrics['median_gap']:.2f}"),
        ("WP80 fake", f"{100.0 * metrics['wp80_background_fake_rate']:.1f}%"),
    ]


def _draw_metric_box(ax, metrics: dict) -> list:
    rows = _format_metrics(metrics)
    box = Rectangle(
        (0.02, 0.56),
        0.36,
        0.40,
        transform=ax.transAxes,
        facecolor="white",
        edgecolor="#9aaec9",
        linewidth=0.8,
        zorder=5,
    )
    ax.add_patch(box)
    artists = []
    for i, (name, value) in enumerate(rows):
        y = 0.88 - i * 0.105
        artists.append(ax.text(
            0.035,
            y,
            name,
            transform=ax.transAxes,
            ha="left",
            va="center",
            fontsize=10,
            family=FONT_FAMILY,
            color="#111827",
            zorder=6,
        ))
        artists.append(ax.text(
            0.365,
            y,
            value,
            transform=ax.transAxes,
            ha="right",
            va="center",
            fontsize=10,
            family=FONT_FAMILY,
            fontweight="bold",
            color="#111827",
            zorder=6,
        ))
    return artists


def _draw_row_label(fig, row_key: str, y0: float, height: float) -> list:
    accent, fill = ROW_COLORS[row_key]
    artists = []
    x0 = 0.015
    width = 0.122
    card = Rectangle(
        (x0, y0 + 0.005),
        width,
        height - 0.010,
        transform=fig.transFigure,
        facecolor=fill,
        edgecolor=accent,
        linewidth=0.9,
        joinstyle="round",
    )
    fig.patches.append(card)
    fig.patches.append(
        Rectangle(
            (x0 + 0.018, y0 + 0.035),
            0.006,
            height - 0.070,
            transform=fig.transFigure,
            facecolor=accent,
            edgecolor="none",
        )
    )
    artists.append(fig.text(
        x0 + 0.070,
        y0 + height * 0.630,
        "VALIDATION",
        ha="center",
        va="center",
        fontsize=13,
        family=FONT_FAMILY,
        fontweight="bold",
        color="#111827",
    ))
    artists.append(fig.text(
        x0 + 0.070,
        y0 + height * 0.340,
        LABELS[row_key],
        ha="center",
        va="center",
        fontsize=13,
        family=FONT_FAMILY,
        fontweight="bold",
        linespacing=0.82,
        color="#111827",
    ))
    return artists


def _draw_col_label(fig, col_key: str, x0: float, width: float) -> list:
    accent, fill = COL_COLORS[col_key]
    card = Rectangle(
        (x0, 0.760),
        width,
        0.055,
        transform=fig.transFigure,
        facecolor=fill,
        edgecolor=accent,
        linewidth=0.9,
        joinstyle="round",
    )
    fig.patches.append(card)
    artist = fig.text(
        x0 + width / 2,
        0.787,
        f"TRAINING: {COL_LABELS[col_key]}",
        ha="center",
        va="center",
        fontsize=15,
        family=FONT_FAMILY,
        fontweight="bold",
        color="#111827",
    )
    return [artist]


def render() -> Path:
    _configure_fonts()
    cells = _load_cells()
    output = OUT_DIR / "the59_sample_composition_controls_candidate_v10_sourceclass_slide24_style_times_checker.png"
    manifest = output.with_name(output.stem + "_manifest.json")
    layout_path = output.with_suffix(".layout_nodes.json")

    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI, facecolor="white")
    text_nodes = []
    title_artist = fig.text(
        0.04,
        0.94,
        "0-20% centrality: fixed validation sample vs trained BDT",
        fontsize=30,
        family=FONT_FAMILY,
        fontweight="bold",
        color="#111827",
    )
    text_nodes.append(("title", title_artist, 30, "title"))

    legend_x, legend_y, legend_w, legend_h = 0.215, 0.846, 0.74, 0.061
    fig.patches.append(
        Rectangle(
            (legend_x, legend_y),
            legend_w,
            legend_h,
            transform=fig.transFigure,
            facecolor="white",
            edgecolor="#9aaec9",
            linewidth=0.9,
            joinstyle="round",
        )
    )
    fig.lines.append(
        plt.Line2D(
            [legend_x + 0.05, legend_x + 0.095],
            [legend_y + 0.031, legend_y + 0.031],
            transform=fig.transFigure,
            color="#d13b35",
            linewidth=3.0,
        )
    )
    legend_signal = fig.text(
        legend_x + 0.106,
        legend_y + 0.031,
        "Signal MC (truth-isolated prompt)",
        ha="left",
        va="center",
        fontsize=13,
        family=FONT_FAMILY,
        fontweight="bold",
        color="#111827",
    )
    text_nodes.append(("legend signal label", legend_signal, 13, "audience"))
    fig.lines.append(
        plt.Line2D(
            [legend_x + 0.41, legend_x + 0.455],
            [legend_y + 0.031, legend_y + 0.031],
            transform=fig.transFigure,
            color="#1f78b4",
            linewidth=3.0,
        )
    )
    legend_background = fig.text(
        legend_x + 0.465,
        legend_y + 0.031,
        "Inclusive MC (embedded jet; no truth filter)",
        ha="left",
        va="center",
        fontsize=13,
        family=FONT_FAMILY,
        fontweight="bold",
        color="#111827",
    )
    text_nodes.append(("legend inclusive label", legend_background, 13, "audience"))

    grid_left = 0.188
    plot_w = 0.246
    gap_x = 0.026
    row_top = 0.742
    plot_h = 0.192
    gap_y = 0.042
    axes = []
    for col_idx, col_key in enumerate(COL_ORDER):
        x0 = grid_left + col_idx * (plot_w + gap_x)
        for artist in _draw_col_label(fig, col_key, x0, plot_w):
            text_nodes.append((f"column label {col_key}", artist, 15, "audience"))
    for row_idx, row_key in enumerate(ROW_ORDER):
        y0 = row_top - row_idx * (plot_h + gap_y) - plot_h
        for idx, artist in enumerate(_draw_row_label(fig, row_key, y0, plot_h)):
            points = 13
            text_nodes.append((f"row label {row_key} {idx}", artist, points, "audience"))
        for col_idx, col_key in enumerate(COL_ORDER):
            x0 = grid_left + col_idx * (plot_w + gap_x)
            ax = fig.add_axes([x0, y0, plot_w, plot_h])
            branch = cells[(row_key, col_key)]
            xs, ys = _step_xy(branch, "signal")
            xb, yb = _step_xy(branch, "background")
            ax.fill_between(xs, ys, step="pre", color="#d13b35", alpha=0.10, linewidth=0)
            ax.fill_between(xb, yb, step="pre", color="#1f78b4", alpha=0.10, linewidth=0)
            ax.plot(xs, ys, color="#d13b35", linewidth=2.1)
            ax.plot(xb, yb, color="#1f78b4", linewidth=2.1)
            ax.set_xlim(0.0, 1.0)
            ymax = 12 if row_key == "baseline" else 10
            ax.set_ylim(0, ymax)
            ax.grid(True, color="#dde5f0", linewidth=0.7)
            ax.tick_params(
                axis="both",
                which="both",
                direction="in",
                top=True,
                right=True,
                labelsize=10,
                width=1.1,
            )
            ax.yaxis.set_major_locator(MaxNLocator(nbins=3, integer=True))
            if col_idx == 0:
                ax.set_ylabel("Unit-area\ndensity", fontsize=13, family=FONT_FAMILY, labelpad=0)
            else:
                ax.set_yticklabels([])
            if row_idx == len(ROW_ORDER) - 1:
                ax.set_xlabel("BDT score", fontsize=15, family=FONT_FAMILY, labelpad=0)
            else:
                ax.set_xticklabels([])
            for spine in ax.spines.values():
                spine.set_color("#1f2937")
                spine.set_linewidth(1.1)
            for metric_idx, artist in enumerate(_draw_metric_box(ax, branch["metrics"])):
                text_nodes.append((f"metric box {row_key} {col_key} {metric_idx}", artist, 10, "plot_annotation"))
            axes.append(ax)

    for ax in axes:
        for item in ax.get_xticklabels() + ax.get_yticklabels():
            item.set_fontfamily(FONT_FAMILY)
        ax.xaxis.label.set_fontfamily(FONT_FAMILY)
        ax.yaxis.label.set_fontfamily(FONT_FAMILY)

    nodes = []
    for name, artist, points, role in text_nodes:
        node = {
            "name": name,
            "kind": "text",
            "role": role,
            "text": artist.get_text(),
            "font_px": _font_px(points),
            "bbox": _bbox_px(fig, artist),
            "title_anchor": role == "title",
        }
        if name.startswith("column label "):
            node["colon_style_exception"] = True
        nodes.append(node)
    layout_path.write_text(
        json.dumps(
            {
                "schema": "slide_layout_nodes_v1",
                "title_axis_x": 0.04 * 2560,
                "minimum_audience_font_px": _font_px(13),
                "minimum_title_font_px": _font_px(24),
                "minimum_plot_annotation_font_px": _font_px(10),
                "nodes": nodes,
            },
            indent=2,
        )
        + "\n"
    )

    fig.savefig(output, dpi=SLIDE_DPI)
    plt.close(fig)

    summary = []
    for row_key in ROW_ORDER:
        for col_key in COL_ORDER:
            branch = cells[(row_key, col_key)]
            metrics = branch["metrics"]
            entries = branch["entries"]
            summary.append(
                {
                    "validation_sample": COL_LABELS[row_key].replace("\n", ""),
                    "training_sample": COL_LABELS[col_key],
                    "signal_entries": entries["signal"],
                    "background_entries": entries["background"],
                    "auc": metrics["auc"],
                    "logloss": metrics["logloss"],
                    "median_gap": metrics["median_gap"],
                    "wp80_fake": metrics["wp80_background_fake_rate"],
                }
            )

    manifest.write_text(
        json.dumps(
            {
                "campaign_tag": "THE8_corrected_truthdefault_20260615",
                "metrics_json": str(METRICS_JSON),
                "png": str(output),
                "layout_nodes": str(layout_path),
                "change": "row label text contained; canonical 2560x1440 render; Times New Roman registered from macOS font files",
                "font_family": FONT_FAMILY,
                "font_files": [str(p) for p in TIMES_FONTS if p.exists()],
                "size_px": [2560, 1440],
                "dpi": SLIDE_DPI,
                "slide_number_baked_in": False,
                "deck_mutated": False,
                "tiny_provenance_footer_baked_in": False,
                "summary": summary,
            },
            indent=2,
        )
        + "\n"
    )
    return output


def main() -> None:
    print(render())


if __name__ == "__main__":
    main()
