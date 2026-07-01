#!/usr/bin/env python3
"""Render a scalar metric-delta summary for THE37 sample-grid slides."""

from __future__ import annotations

import csv
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
THE8_METRIC_TABLE = (
    THE8_CONTROL_DIR
    / "the8_branchA_fixed_sample_holdout_3x3_direct_0to20_truthSignal_inclusiveJet_slide_v35_logloss_metric_table.metric_table_summary.csv"
)
FIXED_VALIDATION_PAYLOAD = THE37_DIR / "the37_fixed_validation_score_shape_payload_20260621.json"
AUAU_SIGNAL_TRAINVAL_PAYLOAD = THE37_DIR / "the37_auau_fixed_inclusive_signal_trainval_matrix_payload_20260621.json"
PP_INCLUSIVE_TRAINVAL_PAYLOAD = THE37_DIR / "the37_pp_fixed_signal_inclusive_trainval_matrix_payload_20260621.json"
PP_SIGNAL_TRAINVAL_PAYLOAD = THE37_DIR / "the37_pp_fixed_inclusive_signal_trainval_matrix_payload_20260621.json"

OUT_DIR = THE37_DIR / "metric_delta_summary_20260628"
OUT_PNG = OUT_DIR / "the37_sample_grid_metric_delta_summary_slide.png"

TIMES_DIR = Path("/System/Library/Fonts/Supplemental")
TIMES_FONTS = [
    TIMES_DIR / "Times New Roman.ttf",
    TIMES_DIR / "Times New Roman Bold.ttf",
    TIMES_DIR / "Times New Roman Italic.ttf",
    TIMES_DIR / "Times New Roman Bold Italic.ttf",
]
FONT = "Times New Roman"

INK = "#111827"
MUTED = "#5b6472"
GRID = "#dfe7f1"
EDGE = "#c9d5e7"
SOFT_BLUE = "#eef6ff"
SOFT_YELLOW = "#fff6d8"
SOFT_RED = "#fff1ef"
READOUT_EDGE = "#d8a924"
COLORS = {
    "auau_inclusive": "#17844a",
    "pp_inclusive": "#1585c7",
    "auau_signal": "#d99500",
    "pp_signal": "#7e55c7",
}


@dataclass(frozen=True)
class Metrics:
    auc: float
    logloss: float
    median_gap: float
    wp80_pass_percent: float


@dataclass(frozen=True)
class Comparison:
    key: str
    label: str
    short_label: str
    color: str
    start_label: str
    end_label: str
    source_slide: str
    start: Metrics
    end: Metrics

    def deltas(self) -> dict[str, float]:
        return {
            "auc": 1000.0 * (self.end.auc - self.start.auc),
            "logloss": 1000.0 * (self.start.logloss - self.end.logloss),
            "median_gap": 1000.0 * (self.end.median_gap - self.start.median_gap),
            "wp80": self.start.wp80_pass_percent - self.end.wp80_pass_percent,
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
    zorder: int = 8,
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
    radius: float = 0.009,
    zorder: int = -10,
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


def load_the8_auau_inclusive_metrics() -> dict[tuple[str, str], Metrics]:
    rows: dict[tuple[str, str], Metrics] = {}
    with THE8_METRIC_TABLE.open(newline="") as handle:
        for row in csv.DictReader(handle):
            rows[(row["validation_sample"], row["training_sample"])] = Metrics(
                auc=float(row["exact_auc_source_classes"]),
                logloss=float(row["exact_logloss_source_classes"]),
                median_gap=float(row["median_bdt_score_gap_binned"]),
                wp80_pass_percent=100.0 * float(row["wp80_background_fake_rate_binned"]),
            )
    return rows


def payload_cells(path: Path, key: str) -> dict[tuple[str, str], Metrics]:
    payload = json.loads(path.read_text())[key]
    out: dict[tuple[str, str], Metrics] = {}
    for cell in payload["cells"]:
        metrics = cell["metrics"]
        out[(cell["row_label"], cell["col_label"])] = Metrics(
            auc=float(metrics["auc"]),
            logloss=float(metrics["logloss"]),
            median_gap=float(metrics["median_gap"]),
            wp80_pass_percent=100.0 * float(metrics["wp80_background_fake_rate"]),
        )
    return out


def build_comparisons() -> list[Comparison]:
    auau_inclusive = load_the8_auau_inclusive_metrics()
    pp_inclusive = payload_cells(PP_INCLUSIVE_TRAINVAL_PAYLOAD, "pp_inclusive_trainval")
    auau_signal = payload_cells(AUAU_SIGNAL_TRAINVAL_PAYLOAD, "auau_signal_trainval")
    pp_signal = payload_cells(PP_SIGNAL_TRAINVAL_PAYLOAD, "pp_signal_trainval")
    return [
        Comparison(
            key="auau_inclusive",
            label="Au+Au inclusive-MC ladder",
            short_label="AuAu incl.",
            color=COLORS["auau_inclusive"],
            start_label="Jet12+20 trained BDT on full validation",
            end_label="Jet12+20+30+40 trained BDT on full validation",
            source_slide="slide 30",
            start=auau_inclusive[("Jet12+20+30+40", "Jet12+20")],
            end=auau_inclusive[("Jet12+20+30+40", "Jet12+20+30+40")],
        ),
        Comparison(
            key="pp_inclusive",
            label="pp inclusive-MC ladder",
            short_label="pp incl.",
            color=COLORS["pp_inclusive"],
            start_label="Jet8+12 trained BDT on full validation",
            end_label="Jet8+12+20+30+40 trained BDT on full validation",
            source_slide="slide 31",
            start=pp_inclusive[("Jet8+12+20+30+40", "Jet8+12")],
            end=pp_inclusive[("Jet8+12+20+30+40", "Jet8+12+20+30+40")],
        ),
        Comparison(
            key="auau_signal",
            label="Au+Au signal ladder",
            short_label="AuAu sig.",
            color=COLORS["auau_signal"],
            start_label="Photon12-only trained BDT on full validation",
            end_label="Photon12+20 trained BDT on full validation",
            source_slide="slide 32",
            start=auau_signal[("Photon12 only", "Photon12+20")],
            end=auau_signal[("Photon12+20", "Photon12+20")],
        ),
        Comparison(
            key="pp_signal",
            label="pp signal ladder",
            short_label="pp sig.",
            color=COLORS["pp_signal"],
            start_label="Photon5 trained BDT on full validation",
            end_label="Photon5+10+20 trained BDT on full validation",
            source_slide="slide 33",
            start=pp_signal[("Photon5", "Photon5+10+20")],
            end=pp_signal[("Photon5+10+20", "Photon5+10+20")],
        ),
    ]


def stress_cells() -> list[dict[str, object]]:
    auau_signal = payload_cells(AUAU_SIGNAL_TRAINVAL_PAYLOAD, "auau_signal_trainval")
    pp_signal = payload_cells(PP_SIGNAL_TRAINVAL_PAYLOAD, "pp_signal_trainval")
    return [
        {
            "label": "Au+Au off-diagonal",
            "text": "Photon12+20 train -> Photon12-only validation",
            "metrics": auau_signal[("Photon12+20", "Photon12 only")],
            "source_slide": "slide 32",
        },
        {
            "label": "pp off-diagonal",
            "text": "Photon5+10 train -> Photon5+10+20 validation",
            "metrics": pp_signal[("Photon5+10", "Photon5+10+20")],
            "source_slide": "slide 33",
        },
    ]


def draw_metric_axis(
    ax: plt.Axes,
    comps: list[Comparison],
    metric_key: str,
    *,
    title: str,
    ylim: tuple[float, float],
    label_fmt,
    clip_at: float | None = None,
    clipped_label: str | None = None,
) -> None:
    x = np.arange(len(comps))
    values = np.array([comp.deltas()[metric_key] for comp in comps], dtype=float)
    plotted = np.minimum(values, clip_at) if clip_at is not None else values
    colors = [comp.color for comp in comps]
    bars = ax.bar(x, plotted, width=0.58, color=colors, edgecolor="white", linewidth=1.0, zorder=3)
    ax.axhline(0, color="#94a3b8", linewidth=1.1, zorder=2)
    ax.set_ylim(*ylim)
    ax.set_xlim(-0.55, len(comps) - 0.45)
    ax.grid(axis="y", color=GRID, linewidth=0.9, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(axis="both", labelsize=12.5, width=1.0, colors=INK)
    ax.set_xticks(x)
    ax.set_xticklabels([comp.short_label for comp in comps], fontsize=12.5, fontweight="bold")
    ax.text(
        -0.52,
        ylim[1] * 0.90,
        title,
        ha="left",
        va="top",
        fontsize=15.7,
        fontweight="bold",
        color=INK,
    )
    for spine in ax.spines.values():
        spine.set_color("#b7c4d8")
        spine.set_linewidth(0.9)
    for bar, real_value, shown_value in zip(bars, values, plotted):
        label = label_fmt(real_value)
        if clip_at is not None and real_value > clip_at and clipped_label:
            label = clipped_label.format(value=real_value)
        y = shown_value + (ylim[1] - ylim[0]) * 0.035
        ax.text(
            bar.get_x() + bar.get_width() / 2.0,
            y,
            label,
            ha="center",
            va="bottom",
            fontsize=12.0,
            fontweight="bold",
            color=INK,
            zorder=4,
        )
    if clip_at is not None:
        idxs = np.where(values > clip_at)[0]
        for idx in idxs:
            ax.plot([idx], [clip_at], marker="^", markersize=10, color=colors[idx], markeredgecolor="white", zorder=5)


def write_datapoints(comps: list[Comparison], path: Path) -> None:
    rows = []
    for comp in comps:
        deltas = comp.deltas()
        rows.append(
            {
                "key": comp.key,
                "label": comp.label,
                "source_slide": comp.source_slide,
                "start_label": comp.start_label,
                "end_label": comp.end_label,
                "start_auc": comp.start.auc,
                "end_auc": comp.end.auc,
                "delta_auc_x1e3": deltas["auc"],
                "start_logloss": comp.start.logloss,
                "end_logloss": comp.end.logloss,
                "delta_logloss_reduction_x1e3": deltas["logloss"],
                "start_median_gap": comp.start.median_gap,
                "end_median_gap": comp.end.median_gap,
                "delta_median_gap_x1e3": deltas["median_gap"],
                "start_wp80_pass_percent": comp.start.wp80_pass_percent,
                "end_wp80_pass_percent": comp.end.wp80_pass_percent,
                "delta_wp80_pass_reduction_pp": deltas["wp80"],
            }
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_script(path: Path) -> None:
    path.write_text(
        "\n".join(
            [
                "# WP GammaJets Sample-Ladder Metric Summary Script",
                "",
                "This slide compresses the four score-shape matrices into scalar changes. Every bar is read the same way: positive means better performance when I move from a minimal or transferred sample choice to the matched full-coverage endpoint on the same validation row.",
                "",
                "The first row is AUC, so it asks whether signal is ranked above inclusive background more often. The second row is logloss, which is sensitive to confidence and probability calibration. The third row is the BDT median gap, so it shows the score-axis separation that we visually saw in the histograms. The last row is the WP80 inclusive or fake pass-rate, which is the cut-level leakage metric.",
                "",
                "The important pattern is that the signal-ladder stress tests are not identical between Au+Au and pp. The Au+Au signal ladder only nudges the matched endpoint, while the pp signal ladder has a huge failure mode when the low-threshold photon training does not cover the full validation signal mixture. That is why I keep the off-diagonal stress cells in the small warning band rather than hiding them in the matrix.",
                "",
                "So the takeaway is: matched sample coverage helps, but the metric tells us which part of the BDT changed. Median gap explains visible score separation, and WP80 tells us whether that separation actually removes leakage at the working point.",
                "",
            ]
        )
    )


def render(output_path: Path = OUT_PNG) -> Path:
    configure_fonts()
    comps = build_comparisons()
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
        "Matched sample coverage reduces BDT leakage",
        size=34.0,
        role="title",
        va="top",
        weight="bold",
        title_axis_align="left",
    )
    rounded_box(fig, 0.040, 0.795, 0.920, 0.061, face=SOFT_BLUE, edge="#a7c7ea", lw=1.0, radius=0.006)
    add_text(
        fig,
        nodes,
        "reading rule",
        0.058,
        0.834,
        "Deltas compare minimal or transferred coverage against the matched full-coverage endpoint on the same validation row; positive = better.",
        size=14.0,
        color=INK,
        weight="bold",
    )
    add_text(
        fig,
        nodes,
        "reading caveat",
        0.058,
        0.807,
        "AUC/logloss/median-gap summarize ranking, calibration, and score-axis separation; WP80 is the cut-level leakage check.",
        size=12.2,
        color=MUTED,
    )

    axes = [
        fig.add_axes([0.078, 0.665, 0.855, 0.108]),
        fig.add_axes([0.078, 0.505, 0.855, 0.108]),
        fig.add_axes([0.078, 0.345, 0.855, 0.108]),
        fig.add_axes([0.078, 0.185, 0.855, 0.108]),
    ]
    draw_metric_axis(
        axes[0],
        comps,
        "auc",
        title=r"AUC improvement $\times 10^3$",
        ylim=(0.0, 135.0),
        label_fmt=lambda v: f"+{v:.0f}",
    )
    draw_metric_axis(
        axes[1],
        comps,
        "logloss",
        title=r"Logloss reduction $\times 10^3$",
        ylim=(0.0, 430.0),
        label_fmt=lambda v: f"+{v:.0f}",
        clip_at=410.0,
        clipped_label="+{value:.0f}\nclipped",
    )
    draw_metric_axis(
        axes[2],
        comps,
        "median_gap",
        title=r"BDT median-gap increase $\times 10^3$",
        ylim=(0.0, 970.0),
        label_fmt=lambda v: f"+{v:.0f}" if v >= 0 else f"{v:.0f}",
    )
    draw_metric_axis(
        axes[3],
        comps,
        "wp80",
        title="WP80 pass reduction (pp)",
        ylim=(0.0, 31.0),
        label_fmt=lambda v: f"+{v:.1f}",
    )
    for ax in axes[:-1]:
        ax.set_xticklabels([])
        ax.tick_params(axis="x", length=0)

    stress = stress_cells()
    rounded_box(fig, 0.040, 0.040, 0.920, 0.074, face=SOFT_YELLOW, edge=READOUT_EDGE, lw=1.1, radius=0.007)
    add_text(fig, nodes, "stress label", 0.058, 0.087, "Stress cells to keep in backup", size=13.6, weight="bold", color=INK)
    for idx, item in enumerate(stress):
        x0 = 0.330 + idx * 0.315
        metrics = item["metrics"]
        assert isinstance(metrics, Metrics)
        fig.patches.append(Rectangle((x0, 0.058), 0.0065, 0.034, transform=fig.transFigure, facecolor="#c24135", edgecolor="none", zorder=6))
        add_text(
            fig,
            nodes,
            f"stress {idx}",
            x0 + 0.014,
            0.087,
            f"{item['label']}",
            size=10.8,
            weight="bold",
            color=INK,
            va="center",
        )
        add_text(
            fig,
            nodes,
            f"stress metrics {idx}",
            x0 + 0.014,
            0.061,
            f"AUC {metrics.auc:.3f}, gap {metrics.median_gap:.2f}, WP80 {metrics.wp80_pass_percent:.1f}%",
            size=10.7,
            color=INK,
            va="center",
        )

    datapoints_path = output_path.with_name(output_path.stem + "_datapoints.csv")
    write_datapoints(comps, datapoints_path)
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
    layout_path = output_path.with_suffix(".layout_nodes.json")
    layout_path.write_text(
        json.dumps(
            {
                "schema": "slide_layout_nodes_v1",
                "title_axis_x": 0.040 * 2560,
                "minimum_audience_font_px": font_px(10.2),
                "minimum_title_font_px": font_px(28.0),
                "minimum_plot_annotation_font_px": font_px(9.0),
                "vertical_margin_balance": {
                    "top_node": "title",
                    "bottom_node": "stress metrics 1",
                    "target_gap_px": 58.0,
                    "tolerance_px": 42.0,
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
                "schema": "the37_sample_grid_metric_delta_summary_v1",
                "png": str(output_path),
                "datapoints_csv": str(datapoints_path),
                "layout_nodes": str(layout_path),
                "speaker_script": str(script_path),
                "deck_mutated": False,
                "slide_numbers_baked_in": False,
                "inputs": {
                    "the8_metric_table": str(THE8_METRIC_TABLE),
                    "fixed_validation_payload": str(FIXED_VALIDATION_PAYLOAD),
                    "auau_signal_trainval_payload": str(AUAU_SIGNAL_TRAINVAL_PAYLOAD),
                    "pp_inclusive_trainval_payload": str(PP_INCLUSIVE_TRAINVAL_PAYLOAD),
                    "pp_signal_trainval_payload": str(PP_SIGNAL_TRAINVAL_PAYLOAD),
                },
                "comparison_contract": [
                    {
                        "key": comp.key,
                        "label": comp.label,
                        "source_slide": comp.source_slide,
                        "start": comp.start_label,
                        "end": comp.end_label,
                        "deltas": comp.deltas(),
                    }
                    for comp in comps
                ],
                "outlier_note": "pp signal-ladder logloss reduction is clipped visually at 410 x10^3; exact value is preserved in label, CSV, and manifest.",
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
