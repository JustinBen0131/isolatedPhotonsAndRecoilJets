#!/usr/bin/env python3
"""Render a plot-first score-distribution summary for PPG19 slides 25-28."""

from __future__ import annotations

import json
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from matplotlib.patches import FancyBboxPatch, Rectangle

THIS_FILE = Path(__file__).resolve()
SCRIPTS_DIR = next((p for p in THIS_FILE.parents if p.name == "scripts"), THIS_FILE.parent)
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.append(str(SCRIPTS_DIR))

from slides.common.slide_defaults import SLIDE_DPI, slide_figsize


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
THE37_DIR = REPO / "dataOutput/auauTightBDTValidation/THE37_sample_grid_slides_20260621"
THE8_CONTROL_DIR = (
    REPO
    / "dataOutput/auauTightBDTValidation/THE8_branchA_ladder_scorecache_fullstat_20260527"
    / "fixed_sample_controls"
)
THE8_SCORE_JSON = THE8_CONTROL_DIR / "the8_branch_a_ladder_fixed_sample_holdout_3x3_direct_truthSignal_inclusiveJet_v3_withCentAuc.json"
PP_INCLUSIVE_TRAINVAL_PAYLOAD = THE37_DIR / "the37_pp_fixed_signal_inclusive_trainval_matrix_payload_20260621.json"
AUAU_SIGNAL_TRAINVAL_PAYLOAD = THE37_DIR / "the37_auau_fixed_inclusive_signal_trainval_matrix_payload_20260621.json"
PP_SIGNAL_TRAINVAL_PAYLOAD = THE37_DIR / "the37_pp_fixed_inclusive_signal_trainval_matrix_payload_20260621.json"

OUT_DIR = THE37_DIR / "score_distribution_story_20260629"
OUT_PNG = OUT_DIR / "the37_slides25_28_score_distribution_story_slide.png"

TIMES_DIR = Path("/System/Library/Fonts/Supplemental")
TIMES_FONTS = [
    TIMES_DIR / "Times New Roman.ttf",
    TIMES_DIR / "Times New Roman Bold.ttf",
    TIMES_DIR / "Times New Roman Italic.ttf",
    TIMES_DIR / "Times New Roman Bold Italic.ttf",
]
FONT = "Times New Roman"

INK = "#111827"
MUTED = "#586273"
GRID = "#dbe5f0"
SIGNAL = "#d92828"
INCL = "#0876bd"
EDGE = "#bbc8d8"
HEADER = "#eef6ff"
NOTE = "#fff5df"
NOTE_EDGE = "#d99b1e"


@dataclass(frozen=True)
class HistPair:
    edges: np.ndarray
    signal_density: np.ndarray
    inclusive_density: np.ndarray
    signal_weights: np.ndarray
    inclusive_weights: np.ndarray


@dataclass(frozen=True)
class Panel:
    key: str
    title: str
    limited_label: str
    matched_label: str
    claim: str
    source_slide: int
    limited: HistPair
    matched: HistPair


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
    zorder: int = 10,
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


def rounded_box(
    fig: plt.Figure,
    x: float,
    y: float,
    w: float,
    h: float,
    *,
    face: str,
    edge: str = EDGE,
    lw: float = 1.0,
    radius: float = 0.008,
    zorder: int = 1,
) -> None:
    fig.patches.append(
        FancyBboxPatch(
            (x, y),
            w,
            h,
            transform=fig.transFigure,
            boxstyle=f"round,pad=0.005,rounding_size={radius}",
            facecolor=face,
            edgecolor=edge,
            linewidth=lw,
            zorder=zorder,
        )
    )


def weights_from_hist(hist: dict) -> np.ndarray:
    values = hist.get("weighted_counts") or hist.get("counts")
    return np.asarray(values, dtype=float)


def hist_pair_from_cell(cell: dict) -> HistPair:
    return HistPair(
        edges=np.asarray(cell["bin_edges"], dtype=float),
        signal_density=np.asarray(cell["signal"]["density"], dtype=float),
        inclusive_density=np.asarray(cell["background"]["density"], dtype=float),
        signal_weights=weights_from_hist(cell["signal"]),
        inclusive_weights=weights_from_hist(cell["background"]),
    )


def hist_pair_from_the8(validation_sample: str, training_sample: str) -> HistPair:
    payload = json.loads(THE8_SCORE_JSON.read_text())
    sample = next(s for s in payload["common_samples"] if s["validation_sample"] == validation_sample)
    branch = next(b for b in sample["branches"] if b["label"] == training_sample)
    cell = branch["by_centrality"]["0_20"]
    return HistPair(
        edges=np.asarray(branch["bin_edges"], dtype=float),
        signal_density=np.asarray(cell["signal"]["density"], dtype=float),
        inclusive_density=np.asarray(cell["background"]["density"], dtype=float),
        signal_weights=weights_from_hist(cell["signal"]),
        inclusive_weights=weights_from_hist(cell["background"]),
    )


def payload_cell(path: Path, key: str, row_label: str, col_label: str) -> dict:
    payload = json.loads(path.read_text())[key]
    return next(c for c in payload["cells"] if c["row_label"] == row_label and c["col_label"] == col_label)


def build_panels() -> list[Panel]:
    return [
        Panel(
            key="auau_inclusive",
            title="Au+Au inclusive MC ladder",
            limited_label="limited: Jet12+20 BDT",
            matched_label="matched: Jet12+20+30+40 BDT",
            claim="Jet30/40 coverage separates red from blue without changing the physics class definition.",
            source_slide=25,
            limited=hist_pair_from_the8("Jet12+20+30+40", "Jet12+20"),
            matched=hist_pair_from_the8("Jet12+20+30+40", "Jet12+20+30+40"),
        ),
        Panel(
            key="pp_inclusive",
            title="pp inclusive MC ladder",
            limited_label="limited: Jet8+12 BDT",
            matched_label="matched: Jet8+12+20+30+40 BDT",
            claim="Hard-jet coverage moves inclusive MC away from the signal-like high-score tail.",
            source_slide=26,
            limited=hist_pair_from_cell(payload_cell(PP_INCLUSIVE_TRAINVAL_PAYLOAD, "pp_inclusive_trainval", "Jet8+12+20+30+40", "Jet8+12")),
            matched=hist_pair_from_cell(payload_cell(PP_INCLUSIVE_TRAINVAL_PAYLOAD, "pp_inclusive_trainval", "Jet8+12+20+30+40", "Jet8+12+20+30+40")),
        ),
        Panel(
            key="auau_signal",
            title="Au+Au signal MC ladder",
            limited_label="limited: Photon12 BDT",
            matched_label="matched: Photon12+20 BDT",
            claim="Full signal coverage is stable; the off-diagonal cells are the diagnostic stress test.",
            source_slide=27,
            limited=hist_pair_from_cell(payload_cell(AUAU_SIGNAL_TRAINVAL_PAYLOAD, "auau_signal_trainval", "Photon12 only", "Photon12+20")),
            matched=hist_pair_from_cell(payload_cell(AUAU_SIGNAL_TRAINVAL_PAYLOAD, "auau_signal_trainval", "Photon12+20", "Photon12+20")),
        ),
        Panel(
            key="pp_signal",
            title="pp signal MC ladder",
            limited_label="limited: Photon5 BDT",
            matched_label="matched: Photon5+10+20 BDT",
            claim="Photon5-only training does not describe the full signal mixture; matched coverage restores separation.",
            source_slide=28,
            limited=hist_pair_from_cell(payload_cell(PP_SIGNAL_TRAINVAL_PAYLOAD, "pp_signal_trainval", "Photon5", "Photon5+10+20")),
            matched=hist_pair_from_cell(payload_cell(PP_SIGNAL_TRAINVAL_PAYLOAD, "pp_signal_trainval", "Photon5+10+20", "Photon5+10+20")),
        ),
    ]


def weighted_quantile(edges: np.ndarray, weights: np.ndarray, probs: list[float]) -> list[float]:
    centers = 0.5 * (edges[:-1] + edges[1:])
    weights = np.asarray(weights, dtype=float)
    if float(weights.sum()) <= 0:
        return [float("nan")] * len(probs)
    order = np.argsort(centers)
    centers = centers[order]
    weights = weights[order]
    cdf = np.cumsum(weights) / float(weights.sum())
    return [float(np.interp(p, cdf, centers)) for p in probs]


def draw_ridge(ax: plt.Axes, hist: HistPair, density: np.ndarray, weights: np.ndarray, y0: float, color: str, label: str, max_density: float) -> None:
    centers = 0.5 * (hist.edges[:-1] + hist.edges[1:])
    scale = 0.54 / max_density if max_density > 0 else 0.0
    clipped_density = np.minimum(density, max_density)
    y = y0 + clipped_density * scale
    ax.fill_between(centers, y0, y, step="mid", color=color, alpha=0.18, linewidth=0)
    ax.step(centers, y, where="mid", color=color, linewidth=1.8)
    q10, q25, q50, q75, q90 = weighted_quantile(hist.edges, weights, [0.10, 0.25, 0.50, 0.75, 0.90])
    ax.plot([q10, q90], [y0 - 0.070, y0 - 0.070], color=color, linewidth=2.0, solid_capstyle="round")
    ax.add_patch(Rectangle((q25, y0 - 0.135), q75 - q25, 0.13, facecolor=color, edgecolor="white", alpha=0.72, linewidth=0.8))
    ax.plot([q50, q50], [y0 - 0.155, y0 + 0.005], color="white", linewidth=1.6)
    ax.text(
        0.020,
        y0 + 0.24,
        label,
        ha="left",
        va="center",
        fontsize=8.7,
        color=color,
        fontweight="bold",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 0.6},
        clip_on=False,
    )


def draw_panel(ax: plt.Axes, panel: Panel) -> None:
    all_densities = [
        panel.limited.signal_density,
        panel.limited.inclusive_density,
        panel.matched.signal_density,
        panel.matched.inclusive_density,
    ]
    positive = np.concatenate([d[d > 0] for d in all_densities if len(d)])
    max_density = float(np.quantile(positive, 0.97)) if len(positive) else 1.0
    y_positions = [3.15, 2.28, 1.13, 0.26]
    draw_ridge(ax, panel.limited, panel.limited.signal_density, panel.limited.signal_weights, y_positions[0], SIGNAL, "limited signal", max_density)
    draw_ridge(ax, panel.limited, panel.limited.inclusive_density, panel.limited.inclusive_weights, y_positions[1], INCL, "limited inclusive", max_density)
    draw_ridge(ax, panel.matched, panel.matched.signal_density, panel.matched.signal_weights, y_positions[2], SIGNAL, "matched signal", max_density)
    draw_ridge(ax, panel.matched, panel.matched.inclusive_density, panel.matched.inclusive_weights, y_positions[3], INCL, "matched inclusive", max_density)
    ax.axvspan(0.80, 1.00, color="#f6e4e4", alpha=0.22, zorder=-5)
    ax.axvline(0.80, color="#9b2f2f", linewidth=0.9, linestyle=(0, (4, 3)), alpha=0.70)
    ax.text(0.020, 4.10, panel.title, ha="left", va="top", fontsize=12.5, fontweight="bold", color=INK)
    ax.text(0.020, 3.78, panel.limited_label, ha="left", va="top", fontsize=9.0, color=INK)
    ax.text(0.020, 1.82, panel.matched_label, ha="left", va="top", fontsize=9.0, color=INK)
    ax.text(0.805, -0.34, "signal-like scores", ha="left", va="center", fontsize=9.0, color="#8a3333")
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.45, 4.22)
    ax.set_xticks([0.0, 0.5, 1.0])
    ax.set_xticklabels(["0", "0.5", "1"], fontsize=9.5)
    ax.set_yticks([])
    ax.grid(axis="x", color=GRID, linewidth=0.8)
    ax.tick_params(axis="x", length=3, colors=INK, pad=2)
    for spine in ax.spines.values():
        spine.set_color(EDGE)
        spine.set_linewidth(0.9)


def write_script(path: Path) -> None:
    path.write_text(
        "\n".join(
            [
                "# PPG19 Slides 25-28 Score-Distribution Story Script",
                "",
                "This version keeps the actual BDT score distributions on the slide. Red is signal MC and blue is inclusive MC. Within each panel, the top two ridges are the limited training choice and the bottom two ridges are the matched full-coverage endpoint. The small boxes show the middle of each distribution, so the audience can see the separation without reading off a table.",
                "",
                "The physics reading is simple. A better BDT pushes the red signal distribution toward high score while keeping the blue inclusive-MC distribution toward low score. In the inclusive-MC rows, adding the harder jet samples mostly removes the high-score inclusive tail. In the signal-MC rows, the Au+Au endpoint is stable, while the pp Photon5-only training clearly does not describe the full Photon5+10+20 validation mixture.",
                "",
                "So the clean takeaway is not a number. It is the shape: matched sample coverage makes the BDT score axis physically interpretable, while mismatched coverage creates score-shape failures that show up directly in the red and blue distributions.",
                "",
            ]
        )
        + "\n"
    )


def render(output_path: Path = OUT_PNG) -> Path:
    configure_fonts()
    panels = build_panels()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI)
    fig.patch.set_facecolor("white")
    nodes: list[tuple[str, object, float, str, dict]] = []

    add_text(
        fig,
        nodes,
        "title",
        0.040,
        0.948,
        "The sample ladder is visible in the BDT score shapes",
        size=30.0,
        role="title",
        va="top",
        weight="bold",
        title_axis_align="left",
    )
    rounded_box(fig, 0.040, 0.816, 0.920, 0.062, face=HEADER, edge="#a7c7ea", lw=1.0, radius=0.006, zorder=1)
    add_text(
        fig,
        nodes,
        "reading_rule",
        0.060,
        0.852,
        "Plot over scorecard - each panel shows actual binned BDT score distributions from slides 25-28.",
        size=14.2,
        weight="bold",
        color=INK,
    )
    add_text(
        fig,
        nodes,
        "reading_key",
        0.060,
        0.827,
        "Good separation means red signal lives at high BDT score while blue inclusive MC stays low; boxes mark score quantiles.",
        size=12.1,
        color=MUTED,
    )
    fig.patches.append(Rectangle((0.760, 0.846), 0.018, 0.014, transform=fig.transFigure, facecolor=SIGNAL, alpha=0.55, edgecolor=SIGNAL, zorder=5))
    add_text(fig, nodes, "legend_signal", 0.782, 0.853, "Signal MC", size=11.2, color=INK)
    fig.patches.append(Rectangle((0.858, 0.846), 0.018, 0.014, transform=fig.transFigure, facecolor=INCL, alpha=0.55, edgecolor=INCL, zorder=5))
    add_text(fig, nodes, "legend_inclusive", 0.880, 0.853, "Inclusive MC", size=11.2, color=INK)

    axes_positions = [
        [0.075, 0.527, 0.405, 0.270],
        [0.555, 0.527, 0.405, 0.270],
        [0.075, 0.186, 0.405, 0.270],
        [0.555, 0.186, 0.405, 0.270],
    ]
    for panel, pos in zip(panels, axes_positions):
        ax = fig.add_axes(pos)
        draw_panel(ax, panel)

    rounded_box(fig, 0.040, 0.046, 0.920, 0.080, face=NOTE, edge=NOTE_EDGE, lw=1.05, radius=0.006, zorder=1)
    add_text(
        fig,
        nodes,
        "takeaway",
        0.060,
        0.094,
        "Physics readout - matched coverage makes the score axis interpretable.",
        size=13.6,
        weight="bold",
        color=INK,
    )
    add_text(
        fig,
        nodes,
        "takeaway_sub",
        0.060,
        0.066,
        "Mismatched coverage shows up as red/blue shape failure, not just as a metric change.",
        size=11.6,
        color=MUTED,
    )

    script_path = output_path.with_name(output_path.stem + "_speaker_script.md")
    write_script(script_path)
    fig.savefig(output_path, dpi=SLIDE_DPI)

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
    layout_nodes.append(
        {
            "name": "takeaway_box",
            "kind": "box",
            "role": "layout",
            "bbox": [0.040 * 2560, (1.0 - (0.046 + 0.080)) * 1440, (0.040 + 0.920) * 2560, (1.0 - 0.046) * 1440],
        }
    )
    layout_path = output_path.with_suffix(".layout_nodes.json")
    layout_path.write_text(
        json.dumps(
            {
                "schema": "slide_layout_nodes_v1",
                "title_axis_x": 0.040 * 2560,
                "minimum_audience_font_px": font_px(9.0),
                "minimum_title_font_px": font_px(28.0),
                "minimum_plot_annotation_font_px": font_px(8.8),
                "vertical_margin_balance": {
                    "top_node": "title",
                    "bottom_node": "takeaway_box",
                    "target_gap_px": 74.0,
                    "tolerance_px": 32.0,
                },
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
                "schema": "the37_slides25_28_score_distribution_story_v1",
                "png": str(output_path),
                "layout_nodes": str(layout_path),
                "speaker_script": str(script_path),
                "deck_mutated": False,
                "source_deck": {
                    "title": "PPG19_6_29_26",
                    "presentation_id": "13G56p6lZl05W3GUtUdKHzh5zcSDRdFjs0UDi0sLiLMM",
                    "slides_summarized": [25, 26, 27, 28],
                },
                "inputs": {
                    "slide25_the8_score_json": str(THE8_SCORE_JSON),
                    "slide26_pp_inclusive_trainval_payload": str(PP_INCLUSIVE_TRAINVAL_PAYLOAD),
                    "slide27_auau_signal_trainval_payload": str(AUAU_SIGNAL_TRAINVAL_PAYLOAD),
                    "slide28_pp_signal_trainval_payload": str(PP_SIGNAL_TRAINVAL_PAYLOAD),
                },
                "panels": [
                    {
                        "key": panel.key,
                        "source_slide": panel.source_slide,
                        "title": panel.title,
                        "limited_label": panel.limited_label,
                        "matched_label": panel.matched_label,
                        "claim": panel.claim,
                    }
                    for panel in panels
                ],
                "plot_contract": "Ridge/box plots use binned BDT score distributions; boxes are weighted quantiles from bin contents, not scalar metric bars.",
            },
            indent=2,
        )
        + "\n"
    )
    plt.close(fig)
    return output_path


def main() -> None:
    out = render()
    print(out)
    print(out.with_suffix(".layout_nodes.json"))


if __name__ == "__main__":
    main()
