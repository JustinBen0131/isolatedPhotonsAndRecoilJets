#!/usr/bin/env python3
"""Build a spacious AUC ordering-versus-median-gap teaching slide."""

from __future__ import annotations

import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.text import Text


THIS_FILE = Path(__file__).resolve()
REPO_ROOT = next(parent for parent in THIS_FILE.parents if (parent / "agent_context").exists())
SCRIPTS_DIR = REPO_ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.append(str(SCRIPTS_DIR))

from slides.common.slide_defaults import (  # noqa: E402
    SLIDE_DPI,
    SLIDE_HEIGHT_PX,
    SLIDE_WIDTH_PX,
    slide_figsize,
)


OUT_DIR = (
    REPO_ROOT
    / "dataOutput/auauTightBDTValidation/THE8_branchA_ladder_scorecache_fullstat_20260527/"
    / "fixed_sample_controls/auc_pair_ordering_teaching_20260611"
)
PNG_PATH = OUT_DIR / "auc_pair_ordering_teaching_slide.png"
NOTE_PATH = OUT_DIR / "auc_pair_ordering_teaching_note.md"
SCRIPT_PATH = OUT_DIR / "auc_pair_ordering_teaching_speaker_script.md"
MANIFEST_PATH = OUT_DIR / "auc_pair_ordering_teaching_manifest.json"
LAYOUT_PATH = OUT_DIR / "auc_pair_ordering_teaching_layout_nodes.json"


COLORS = {
    "ink": "#172033",
    "muted": "#606A78",
    "rule": "#D8DEE8",
    "soft_rule": "#EDF1F6",
    "readout_fill": "#FBFCFE",
    "readout_header": "#F0F4FA",
    "readout_alt": "#F7F9FD",
    "readout_delta": "#F4F7FC",
    "signal": "#C53A32",
    "background": "#1E78B8",
}


Role = Literal["title", "audience"]


@dataclass
class TextRecord:
    name: str
    role: Role
    artist: Text
    font_px: float
    text: str
    title_anchor: bool = False
    title_alignment: str | None = None
    title_band_exception: bool = False
    colon_style_exception: bool = False


TEXT_RECORDS: list[TextRecord] = []


def add_text(
    ax: plt.Axes,
    x: float,
    y: float,
    text: str,
    *,
    name: str,
    size: float,
    role: Role = "audience",
    color: str = COLORS["ink"],
    weight: str = "normal",
    ha: str = "left",
    va: str = "center",
    linespacing: float = 1.15,
    title_alignment: str | None = None,
    title_band_exception: bool = False,
    colon_style_exception: bool = False,
) -> Text:
    artist = ax.text(
        x,
        y,
        text,
        ha=ha,
        va=va,
        fontsize=size,
        color=color,
        fontweight=weight,
        family="Times New Roman",
        linespacing=linespacing,
        transform=ax.transAxes,
    )
    TEXT_RECORDS.append(
        TextRecord(
            name=name,
            role=role,
            artist=artist,
            font_px=size * SLIDE_DPI / 72.0,
            text=text,
            title_anchor=(role == "title"),
            title_alignment=title_alignment,
            title_band_exception=title_band_exception,
            colon_style_exception=colon_style_exception,
        )
    )
    return artist


def add_centered_segments(
    ax: plt.Axes,
    center_x: float,
    y: float,
    segments: list[tuple[str, str, str]],
    *,
    name: str,
    size: float,
    gap: float = 0.012,
) -> None:
    """Draw three short segments centered around a common visual midpoint."""
    if len(segments) != 3:
        raise ValueError("This helper expects left, middle, right segments")
    left_text, left_color, left_weight = segments[0]
    middle_text, middle_color, middle_weight = segments[1]
    right_text, right_color, right_weight = segments[2]

    add_text(
        ax,
        center_x - gap,
        y,
        left_text,
        name=f"{name} left",
        size=size,
        color=left_color,
        weight=left_weight,
        ha="right",
    )
    add_text(
        ax,
        center_x,
        y,
        middle_text,
        name=f"{name} operator",
        size=size,
        color=middle_color,
        weight=middle_weight,
        ha="center",
    )
    add_text(
        ax,
        center_x + gap,
        y,
        right_text,
        name=f"{name} right",
        size=size,
        color=right_color,
        weight=right_weight,
        ha="left",
    )


def add_score_row(ax: plt.Axes, x_label: float, x_value: float, y: float, label: str, value: str, color: str, name: str) -> None:
    add_text(ax, x_label, y, label, name=f"{name} label", size=26, color=color, weight="bold")
    add_text(ax, x_value, y, value, name=f"{name} value", size=34, color=color, weight="bold", ha="left")


def add_gap(ax: plt.Axes, center_x: float, y: float, value: str, name: str) -> None:
    add_text(
        ax,
        center_x - 0.02,
        y,
        "gap =",
        name=f"{name} gap label",
        size=28,
        color=COLORS["ink"],
        weight="bold",
        ha="right",
    )
    add_text(
        ax,
        center_x + 0.02,
        y,
        value,
        name=f"{name} gap value",
        size=34,
        color=COLORS["ink"],
        weight="bold",
        ha="left",
    )


def text_nodes_from_renderer(fig: plt.Figure) -> list[dict[str, object]]:
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    canvas_width = float(renderer.width)
    canvas_height = float(renderer.height)
    scale_x = SLIDE_WIDTH_PX / float(canvas_width)
    scale_y = SLIDE_HEIGHT_PX / float(canvas_height)
    nodes: list[dict[str, object]] = []
    for record in TEXT_RECORDS:
        bbox = record.artist.get_window_extent(renderer=renderer)
        nodes.append(
            {
                "name": record.name,
                "kind": "text",
                "role": record.role,
                "text": record.text,
                "font_px": record.font_px,
                "bbox": [
                    round(float(bbox.x0) * scale_x, 3),
                    round(float(canvas_height - bbox.y1) * scale_y, 3),
                    round(float(bbox.x1) * scale_x, 3),
                    round(float(canvas_height - bbox.y0) * scale_y, 3),
                ],
                "title_anchor": record.title_anchor,
                "title_alignment": record.title_alignment,
                "title_band_exception": record.title_band_exception,
                "colon_style_exception": record.colon_style_exception,
            }
        )
    return nodes


def write_layout_nodes(fig: plt.Figure) -> None:
    nodes = [
        {
            "name": "old model table column",
            "kind": "panel",
            "bbox": [905, 390, 1430, 930],
            "symmetry_group": "model table columns",
        },
        {
            "name": "new model table column",
            "kind": "panel",
            "bbox": [1615, 390, 2140, 930],
            "symmetry_group": "model table columns",
        },
    ]
    nodes.extend(text_nodes_from_renderer(fig))
    payload = {
        "slide_size_px": [SLIDE_WIDTH_PX, SLIDE_HEIGHT_PX],
        "minimum_audience_font_px": 13 * SLIDE_DPI / 72,
        "minimum_title_font_px": 24 * SLIDE_DPI / 72,
        "nodes": nodes,
    }
    LAYOUT_PATH.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def draw_slide() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    TEXT_RECORDS.clear()
    plt.rcParams.update(
        {
            "font.family": ["Times New Roman", "Times", "DejaVu Serif"],
            "figure.dpi": SLIDE_DPI,
            "savefig.dpi": SLIDE_DPI,
        }
    )

    fig = plt.figure(figsize=slide_figsize(SLIDE_DPI), dpi=SLIDE_DPI)
    fig.patch.set_facecolor("white")
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_facecolor("white")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    add_text(
        ax,
        0.060,
        0.920,
        "AUC checks order; median gap checks distance",
        name="slide title",
        role="title",
        size=34,
        weight="bold",
        ha="left",
    )
    add_text(
        ax,
        0.500,
        0.824,
        "Toy pair, not data:",
        name="slide subtitle lead",
        size=20.8,
        color=COLORS["ink"],
        weight="bold",
        ha="center",
        colon_style_exception=True,
    )
    add_text(
        ax,
        0.500,
        0.772,
        "Correct ordering is unchanged; median-gap separation grows.",
        name="slide subtitle second line",
        size=20.8,
        color=COLORS["muted"],
        ha="center",
        colon_style_exception=True,
    )

    table_left = 0.075
    table_right = 0.925
    table_top = 0.680
    table_bottom = 0.340
    label_x = 0.120
    old_x = 0.505
    new_x = 0.780
    col1 = 0.380
    col2 = 0.635

    ax.plot([table_left, table_right], [table_top, table_top], color=COLORS["rule"], lw=2.0, transform=ax.transAxes)
    ax.plot([table_left, table_right], [table_bottom, table_bottom], color=COLORS["rule"], lw=2.0, transform=ax.transAxes)
    ax.plot([col1, col1], [table_bottom, table_top], color=COLORS["soft_rule"], lw=2.0, transform=ax.transAxes)
    ax.plot([col2, col2], [table_bottom, table_top], color=COLORS["soft_rule"], lw=2.0, transform=ax.transAxes)
    for y in [0.640, 0.565, 0.485, 0.410]:
        ax.plot([table_left, table_right], [y, y], color=COLORS["soft_rule"], lw=1.4, transform=ax.transAxes)

    add_text(ax, label_x, 0.656, "TOY PAIR EXAMPLE", name="table label header", size=14.8, color=COLORS["ink"], weight="bold")
    add_text(ax, old_x, 0.656, "OLD BDT", name="old model header", size=20.5, weight="bold", ha="center")
    add_text(ax, new_x, 0.656, "NEW BDT", name="new model header", size=20.5, weight="bold", ha="center")

    add_text(ax, label_x, 0.603, "signal score", name="signal row label", size=18, color=COLORS["signal"], weight="bold")
    add_text(ax, old_x, 0.603, "0.62", name="old signal value", size=26, color=COLORS["signal"], weight="bold", ha="center")
    add_text(ax, new_x, 0.603, "0.86", name="new signal value", size=26, color=COLORS["signal"], weight="bold", ha="center")

    add_text(ax, label_x, 0.525, "background score", name="background row label", size=18, color=COLORS["background"], weight="bold")
    add_text(ax, old_x, 0.525, "0.35", name="old background value", size=26, color=COLORS["background"], weight="bold", ha="center")
    add_text(ax, new_x, 0.525, "0.10", name="new background value", size=26, color=COLORS["background"], weight="bold", ha="center")

    add_text(ax, label_x, 0.452, "AUC question", name="auc row label", size=18, color=COLORS["ink"], weight="bold")
    add_centered_segments(
        ax,
        old_x,
        0.468,
        [
            ("0.62", COLORS["signal"], "bold"),
            (">", COLORS["ink"], "bold"),
            ("0.35", COLORS["background"], "bold"),
        ],
        name="old model inequality",
        size=21,
        gap=0.016,
    )
    add_text(ax, old_x, 0.430, "YES", name="old auc answer", size=16.5, color=COLORS["muted"], weight="bold", ha="center")
    add_centered_segments(
        ax,
        new_x,
        0.468,
        [
            ("0.86", COLORS["signal"], "bold"),
            (">", COLORS["ink"], "bold"),
            ("0.10", COLORS["background"], "bold"),
        ],
        name="new model inequality",
        size=21,
        gap=0.016,
    )
    add_text(ax, new_x, 0.430, "same YES", name="new auc answer", size=16.5, color=COLORS["muted"], weight="bold", ha="center")

    add_text(ax, label_x, 0.371, "score gap", name="gap row label", size=18.5, color=COLORS["ink"], weight="bold")
    add_text(ax, old_x, 0.371, "0.27", name="old model gap value", size=28, color=COLORS["ink"], weight="bold", ha="center")
    add_text(ax, new_x, 0.371, "0.76", name="new model gap value", size=28, color=COLORS["ink"], weight="bold", ha="center")

    add_text(
        ax,
        0.500,
        0.286,
        "Same ordering gives little AUC gain; larger score distance gives the cleaner split.",
        name="central takeaway line",
        size=21,
        color=COLORS["ink"],
        weight="bold",
        ha="center",
    )
    ax.plot([0.130, 0.870], [0.254, 0.254], color=COLORS["rule"], lw=1.4, transform=ax.transAxes)

    read_left = 0.075
    read_right = 0.900
    label_x = 0.115
    auc_x = 0.650
    gap_x = 0.815
    row_size = 17.4
    delta_size = row_size * 0.96
    read_width = read_right - read_left
    for y, height, fill in [
        (0.064, 0.174, COLORS["readout_fill"]),
        (0.207, 0.031, COLORS["readout_header"]),
        (0.158, 0.049, COLORS["readout_alt"]),
        (0.064, 0.049, COLORS["readout_delta"]),
    ]:
        ax.add_patch(
            Rectangle(
                (read_left, y),
                read_width,
                height,
                transform=ax.transAxes,
                facecolor=fill,
                edgecolor="none",
                zorder=0.1,
            )
        )
    add_text(ax, label_x, 0.225, "REAL READOUT", name="bottom label column header", size=16.5, color=COLORS["ink"], weight="bold")
    add_text(ax, auc_x, 0.222, "AUC", name="bottom auc column header", size=13.5, color=COLORS["muted"], weight="bold", ha="center")
    add_text(ax, gap_x, 0.222, "median gap", name="bottom gap column header", size=13.5, color=COLORS["muted"], weight="bold", ha="center")
    ax.plot([read_left, read_right], [0.207, 0.207], color=COLORS["soft_rule"], lw=1.4, transform=ax.transAxes)
    ax.plot([0.565, 0.565], [0.073, 0.225], color=COLORS["soft_rule"], lw=1.4, transform=ax.transAxes)
    ax.plot([0.735, 0.735], [0.073, 0.225], color=COLORS["soft_rule"], lw=1.4, transform=ax.transAxes)

    add_text(ax, label_x, 0.180, "pp reference, 15-35 GeV", name="bottom pp label", size=row_size, color=COLORS["ink"], weight="bold")
    add_text(ax, auc_x, 0.180, "0.934", name="bottom pp auc", size=row_size, color=COLORS["ink"], ha="center")
    add_text(ax, gap_x, 0.180, "0.84", name="bottom pp gap", size=row_size, color=COLORS["ink"], ha="center")

    add_text(ax, label_x, 0.135, "best 0-20% AuAu fixed-validation cell", name="bottom auau label", size=row_size, color=COLORS["ink"], weight="bold")
    add_text(ax, auc_x, 0.135, "0.870", name="bottom auau auc", size=row_size, color=COLORS["ink"], ha="center")
    add_text(ax, gap_x, 0.135, "0.60", name="bottom auau gap", size=row_size, color=COLORS["ink"], ha="center")

    add_text(ax, label_x, 0.095, "pp vs AuAu difference", name="bottom pp vs auau label", size=delta_size, color=COLORS["ink"])
    add_text(ax, auc_x, 0.095, "+7.4%", name="bottom pp vs auau auc", size=delta_size, color=COLORS["ink"], ha="center")
    add_text(ax, gap_x, 0.095, "+40%", name="bottom pp vs auau gap", size=delta_size, color=COLORS["ink"], ha="center")

    ax.plot([read_left, read_right], [0.064, 0.064], color=COLORS["soft_rule"], lw=1.4, transform=ax.transAxes)
    add_text(
        ax,
        0.500,
        0.033,
        "Readout: AUC is the ranking check; median gap is the visible separation check.",
        name="bottom interpretation bullet",
        size=17.0,
        color=COLORS["ink"],
        ha="center",
        colon_style_exception=True,
    )

    write_layout_nodes(fig)
    fig.savefig(PNG_PATH, bbox_inches=None, pad_inches=0)
    plt.close(fig)


def write_notes() -> None:
    NOTE_PATH.write_text(
        """# AUC versus median gap teaching slide

This slide is a toy example, not a data example. The point is to isolate what
AUC sees versus what the median BDT score gap sees. The subtitle intentionally
separates the toy warning from the metric read: the ordering can stay unchanged
while the median-gap separation grows.

## Toy values

| Quantity | Old BDT | New BDT |
| --- | ---: | ---: |
| Signal score | 0.62 | 0.86 |
| Background score | 0.35 | 0.10 |
| Ordering | 0.62 > 0.35 | 0.86 > 0.10 |
| Median-gap analogue | 0.27 | 0.76 |

## Exact spoken interpretation

This is a toy pair, not a data point. The reason it is useful is that it
separates the two questions. AUC asks whether the signal score is above the
background score. In both cases the answer is yes, so this pair gives almost no
new AUC information. Median gap asks a different question: how far apart did the
BDT push the signal-like and background-like scores? That grows from 0.27 to
0.76. So for the visual question of why the peaks look more separated, median
gap is the more direct number. AUC is still important, but it is the ranking
check, not the distance-on-the-score-axis check.

The bottom bullets connect that concept to the score-shape comparisons. The pp
reference in 15-35 GeV has AUC 0.934 and median gap 0.84, computed from the
same binned pp baseV3E score-shape histograms. The best 0-20% AuAu fixed-sample
cell shown here has AUC 0.870 and median gap 0.60. On the displayed values, pp
is 7.4% higher in AUC and 40% higher in median gap than the AuAu cell, which is
why the median-gap contrast is the more visually obvious separation readout.
""",
        encoding="utf-8",
    )

    SCRIPT_PATH.write_text(
        """# THE-51 Calibration Slide Script - AUC versus median gap

This is a toy pair, not a data point. It is here to isolate the question we are
trying to answer.

AUC asks a yes-or-no ordering question: is the signal score above the background
score? In the old BDT, the answer is yes. In the new BDT, the answer is still
yes. So AUC does not get much new information from this pair.

But the median-gap idea asks a different question: how far apart are the two
classes on the score axis? In the toy example, that gap grows from 0.27 to 0.76.
That is the kind of change our eye sees when the signal peak moves closer to one
and the background peak moves closer to zero.

So if the question is whether the BDT ranks signal above background, AUC is the
right number. If the question is why the score distributions look more cleanly
split, the median gap is the more direct number.

That is why the bottom readout compares both. The pp reference has a larger
median gap, 0.84, with AUC about 0.934. The best 0-20 percent AuAu cell here has
AUC about 0.870 and median gap 0.60. On those displayed values, pp is only
about 7.4 percent higher in AUC, but 40 percent higher in median gap. So the
median gap is the number that more directly tracks the visual peak separation.
""",
        encoding="utf-8",
    )


def write_manifest() -> None:
    manifest = {
        "artifact": "auc_pair_ordering_teaching_slide",
        "purpose": "pedagogical toy slide explaining that AUC checks ranking while median gap checks score-axis separation",
        "slide_size_px": [SLIDE_WIDTH_PX, SLIDE_HEIGHT_PX],
        "source_script": str(THIS_FILE.relative_to(REPO_ROOT)),
        "outputs": {
            "png": str(PNG_PATH.relative_to(REPO_ROOT)),
            "note": str(NOTE_PATH.relative_to(REPO_ROOT)),
            "speaker_script": str(SCRIPT_PATH.relative_to(REPO_ROOT)),
            "layout_nodes": str(LAYOUT_PATH.relative_to(REPO_ROOT)),
            "audit_json": str(PNG_PATH.with_suffix(".slide_audit.json").relative_to(REPO_ROOT)),
        },
        "toy_values": {
            "old_model": {"signal": 0.62, "background": 0.35, "ordering": "0.62 > 0.35", "gap": 0.27},
            "new_model": {"signal": 0.86, "background": 0.10, "ordering": "0.86 > 0.10", "gap": 0.76},
        },
        "bottom_metrics": {
            "pp_reference_15_35": {
                "auc_binned_score_hist": 0.9342455454026949,
                "median_gap_binned_score_hist": 0.84,
                "source": "dataOutput/ppPhotonMLPipeline/ppg12_basev3E_currentIAN_rawOverlayEnvFix_20260527_1630/validation/fullsim_shuhang_overlay_raw_inclusive_pt1535/pp_currentian_basev3e_bdt_score_overlay_vs_ppg12_pp_noCent_bdt_15_35_noNPB_eta0-pt3-cut0_summary.json",
            },
            "auau_best_0_20_fixed_sample_cell": {
                "validation_sample": "Jet12+20+30+40",
                "training_sample": "Jet12+20+30+40",
                "auc": 0.8698867691691528,
                "median_gap_binned": 0.6000000000000001,
                "source": "dataOutput/auauTightBDTValidation/THE8_branchA_ladder_scorecache_fullstat_20260527/fixed_sample_controls/the8_branchA_fixed_sample_holdout_3x3_direct_0to20_truthSignal_inclusiveJet_slide_v35_logloss_metric_table.metric_table_summary.csv",
            },
            "displayed_pp_vs_auau_percent_difference": {
                "auc_percent_higher": 7.4,
                "median_gap_percent_higher": 40.0,
                "basis": "percent higher than the displayed AuAu values: (pp - AuAu) / AuAu",
            },
        },
        "conceptual_boundary": "This toy pair is not a data point. AUC is the ranking check; median gap is the score-separation readout.",
        "deck_mutation": "none",
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    draw_slide()
    write_notes()
    write_manifest()
    print(PNG_PATH)


if __name__ == "__main__":
    main()
