#!/usr/bin/env python3
"""Render a box-chart story summary for PPG19 slides 25-28."""

from __future__ import annotations

import csv
import json
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import to_rgb
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
PP_INCLUSIVE_TRAINVAL_PAYLOAD = THE37_DIR / "the37_pp_fixed_signal_inclusive_trainval_matrix_payload_20260621.json"
AUAU_SIGNAL_TRAINVAL_PAYLOAD = THE37_DIR / "the37_auau_fixed_inclusive_signal_trainval_matrix_payload_20260621.json"
PP_SIGNAL_TRAINVAL_PAYLOAD = THE37_DIR / "the37_pp_fixed_inclusive_signal_trainval_matrix_payload_20260621.json"

OUT_DIR = THE37_DIR / "box_story_summary_20260628"
OUT_PNG = OUT_DIR / "the37_slides25_28_boxchart_story_summary_slide.png"

TIMES_DIR = Path("/System/Library/Fonts/Supplemental")
TIMES_FONTS = [
    TIMES_DIR / "Times New Roman.ttf",
    TIMES_DIR / "Times New Roman Bold.ttf",
    TIMES_DIR / "Times New Roman Italic.ttf",
    TIMES_DIR / "Times New Roman Bold Italic.ttf",
]
FONT = "Times New Roman"
ARROW = "\u2192"

INK = "#111827"
MUTED = "#536070"
GRID = "#d9e2ee"
EDGE = "#c5d0df"
HEADER_BLUE = "#eef6ff"
WARNING_FILL = "#fff5e5"
WARNING_EDGE = "#d89a1e"
RED = "#c23a35"
GREEN = "#16864a"

COLORS = {
    "auau_inclusive": "#17844a",
    "pp_inclusive": "#1585c7",
    "auau_signal": "#d99500",
    "pp_signal": "#7e55c7",
}

METRIC_SCALES = {
    "auc": 120.0,
    "logloss": 420.0,
    "gap": 900.0,
    "wp80": 20.0,
}


@dataclass(frozen=True)
class Metrics:
    auc: float
    logloss: float
    median_gap: float
    wp80_pass_percent: float


@dataclass(frozen=True)
class StoryRow:
    key: str
    group: str
    title: str
    message: str
    comparison: str
    validation: str
    color: str
    source_slide: int
    start_label: str
    end_label: str
    start: Metrics
    end: Metrics

    def deltas(self) -> dict[str, float]:
        return {
            "auc": 1000.0 * (self.end.auc - self.start.auc),
            "logloss": 1000.0 * (self.start.logloss - self.end.logloss),
            "gap": 1000.0 * (self.end.median_gap - self.start.median_gap),
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
    alpha: float = 1.0,
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
            alpha=alpha,
            zorder=zorder,
        )
    )


def blend(color: str, amount: float) -> tuple[float, float, float]:
    base = to_rgb(color)
    amount = max(0.0, min(1.0, amount))
    return tuple((1.0 - amount) + amount * c for c in base)


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


def build_story_rows() -> list[StoryRow]:
    auau_inclusive = load_the8_auau_inclusive_metrics()
    pp_inclusive = payload_cells(PP_INCLUSIVE_TRAINVAL_PAYLOAD, "pp_inclusive_trainval")
    auau_signal = payload_cells(AUAU_SIGNAL_TRAINVAL_PAYLOAD, "auau_signal_trainval")
    pp_signal = payload_cells(PP_SIGNAL_TRAINVAL_PAYLOAD, "pp_signal_trainval")

    return [
        StoryRow(
            key="auau_inclusive",
            group="Inclusive MC coverage",
            title="Au+Au inclusive MC",
            message="Modest but consistent gain",
            comparison="Jet12+20 train to Jet12+20+30+40 train",
            validation="Full validation Jet12+20+30+40",
            color=COLORS["auau_inclusive"],
            source_slide=25,
            start_label="Jet12+20",
            end_label="Jet12+20+30+40",
            start=auau_inclusive[("Jet12+20+30+40", "Jet12+20")],
            end=auau_inclusive[("Jet12+20+30+40", "Jet12+20+30+40")],
        ),
        StoryRow(
            key="pp_inclusive",
            group="Inclusive MC coverage",
            title="pp inclusive MC",
            message="Hard-jet coverage removes leakage",
            comparison="Jet8+12 train to Jet8+12+20+30+40 train",
            validation="Full validation Jet8+12+20+30+40",
            color=COLORS["pp_inclusive"],
            source_slide=26,
            start_label="Jet8+12",
            end_label="Jet8+12+20+30+40",
            start=pp_inclusive[("Jet8+12+20+30+40", "Jet8+12")],
            end=pp_inclusive[("Jet8+12+20+30+40", "Jet8+12+20+30+40")],
        ),
        StoryRow(
            key="auau_signal",
            group="Signal MC coverage",
            title="Au+Au signal MC",
            message="Endpoint stable; mismatch is diagnostic",
            comparison="Photon12 train to Photon12+20 train",
            validation="Full validation Photon12+20",
            color=COLORS["auau_signal"],
            source_slide=27,
            start_label="Photon12",
            end_label="Photon12+20",
            start=auau_signal[("Photon12 only", "Photon12+20")],
            end=auau_signal[("Photon12+20", "Photon12+20")],
        ),
        StoryRow(
            key="pp_signal",
            group="Signal MC coverage",
            title="pp signal MC",
            message="Unmatched low-threshold training fails",
            comparison="Photon5 train to Photon5+10+20 train",
            validation="Full validation Photon5+10+20",
            color=COLORS["pp_signal"],
            source_slide=28,
            start_label="Photon5",
            end_label="Photon5+10+20",
            start=pp_signal[("Photon5", "Photon5+10+20")],
            end=pp_signal[("Photon5+10+20", "Photon5+10+20")],
        ),
    ]


def stress_cells() -> list[dict[str, object]]:
    auau_signal = payload_cells(AUAU_SIGNAL_TRAINVAL_PAYLOAD, "auau_signal_trainval")
    pp_signal = payload_cells(PP_SIGNAL_TRAINVAL_PAYLOAD, "pp_signal_trainval")
    return [
        {
            "label": "Au+Au signal off-diagonal",
            "text": "Photon12+20 train on Photon12-only validation",
            "metrics": auau_signal[("Photon12+20", "Photon12 only")],
            "source_slide": 27,
        },
        {
            "label": "pp signal off-diagonal",
            "text": "Photon5+10 train on Photon5+10+20 validation",
            "metrics": pp_signal[("Photon5+10", "Photon5+10+20")],
            "source_slide": 28,
        },
    ]


def metric_specs() -> list[dict[str, object]]:
    return [
        {
            "key": "auc",
            "header": "AUC",
            "subheader": r"improvement $\times 10^3$",
            "delta_label": lambda v: f"+{v:.0f}",
            "value_label": lambda a, b: f"{a.auc:.3f} {ARROW} {b.auc:.3f}",
        },
        {
            "key": "logloss",
            "header": "Logloss",
            "subheader": r"reduction $\times 10^3$",
            "delta_label": lambda v: f"+{v:.0f}",
            "value_label": lambda a, b: f"{a.logloss:.3f} {ARROW} {b.logloss:.3f}",
        },
        {
            "key": "gap",
            "header": "BDT median gap",
            "subheader": r"increase $\times 10^3$",
            "delta_label": lambda v: f"+{v:.0f}",
            "value_label": lambda a, b: f"{a.median_gap:.2f} {ARROW} {b.median_gap:.2f}",
        },
        {
            "key": "wp80",
            "header": "WP80 incl. pass",
            "subheader": "reduction, pp",
            "delta_label": lambda v: f"+{v:.1f}",
            "value_label": lambda a, b: f"{a.wp80_pass_percent:.1f}% {ARROW} {b.wp80_pass_percent:.1f}%",
        },
    ]


def draw_metric_cell(
    fig: plt.Figure,
    nodes: list[tuple[str, object, float, str, dict]],
    row: StoryRow,
    spec: dict[str, object],
    x: float,
    y: float,
    w: float,
    h: float,
) -> None:
    deltas = row.deltas()
    key = str(spec["key"])
    value = float(deltas[key])
    norm = min(abs(value) / METRIC_SCALES[key], 1.0)
    face = blend(GREEN if value >= 0 else RED, 0.08 + 0.22 * norm)
    edge = row.color
    rounded_box(fig, x, y, w, h, face=face, edge=edge, lw=1.15, radius=0.006, zorder=2)

    delta_text = spec["delta_label"](value)  # type: ignore[index,operator]
    if key == "logloss" and value > METRIC_SCALES["logloss"]:
        delta_text = f"+{value:.0f}"
    add_text(
        fig,
        nodes,
        f"{row.key}_{key}_delta",
        x + w / 2.0,
        y + h * 0.66,
        delta_text,
        size=19.0,
        ha="center",
        weight="bold",
        color=GREEN if value >= 0 else RED,
    )
    add_text(
        fig,
        nodes,
        f"{row.key}_{key}_values",
        x + w / 2.0,
        y + h * 0.37,
        spec["value_label"](row.start, row.end),  # type: ignore[index,operator]
        size=11.9,
        ha="center",
        color=INK,
    )
    if key == "logloss" and value > METRIC_SCALES["logloss"]:
        add_text(
            fig,
            nodes,
            f"{row.key}_{key}_clip",
            x + w / 2.0,
            y + h * 0.17,
            "scale-clipped",
            size=10.8,
            ha="center",
            color=MUTED,
        )


def write_datapoints(rows: list[StoryRow], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "key",
        "source_slide",
        "group",
        "title",
        "comparison",
        "validation",
        "start_label",
        "end_label",
        "start_auc",
        "end_auc",
        "delta_auc_x1e3",
        "start_logloss",
        "end_logloss",
        "delta_logloss_reduction_x1e3",
        "start_median_gap",
        "end_median_gap",
        "delta_median_gap_x1e3",
        "start_wp80_pass_percent",
        "end_wp80_pass_percent",
        "delta_wp80_pass_reduction_pp",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            deltas = row.deltas()
            writer.writerow(
                {
                    "key": row.key,
                    "source_slide": row.source_slide,
                    "group": row.group,
                    "title": row.title,
                    "comparison": row.comparison,
                    "validation": row.validation,
                    "start_label": row.start_label,
                    "end_label": row.end_label,
                    "start_auc": row.start.auc,
                    "end_auc": row.end.auc,
                    "delta_auc_x1e3": deltas["auc"],
                    "start_logloss": row.start.logloss,
                    "end_logloss": row.end.logloss,
                    "delta_logloss_reduction_x1e3": deltas["logloss"],
                    "start_median_gap": row.start.median_gap,
                    "end_median_gap": row.end.median_gap,
                    "delta_median_gap_x1e3": deltas["gap"],
                    "start_wp80_pass_percent": row.start.wp80_pass_percent,
                    "end_wp80_pass_percent": row.end.wp80_pass_percent,
                    "delta_wp80_pass_reduction_pp": deltas["wp80"],
                }
            )


def write_script(path: Path) -> None:
    path.write_text(
        "\n".join(
            [
                "# PPG19 Slides 25-28 Box-Chart Summary Script",
                "",
                "This slide is meant to replace four detailed matrix slides with one readout. Each row takes the same idea from the source slide: compare a limited training sample to the matched full-coverage endpoint on the same full validation row.",
                "",
                "The first two rows are the inclusive-MC ladder. In Au+Au, adding Jet30 and Jet40 gives a modest but consistent improvement: AUC moves up, logloss comes down, the BDT median gap grows, and the WP80 inclusive pass rate drops by about one percentage point. In pp, the same coverage question is much sharper. The low-coverage Jet8+12 training leaves much more inclusive leakage on the full validation sample, and the matched full inclusive endpoint cuts that leakage by nearly twenty percentage points.",
                "",
                "The next two rows are the signal-MC ladder. Au+Au is fairly stable when I move from Photon12 to Photon12+20 on the Photon12+20 validation row. The pp signal ladder is not stable under that mismatch: a Photon5-only training model is badly calibrated on the full Photon5+10+20 validation mixture, and the matched full signal endpoint restores both the median gap and the WP80 leakage behavior.",
                "",
                "The warning band is important. The off-diagonal cells are diagnostic stress tests, not production candidates. They show that unmatched signal coverage can create very large failures, especially in pp. The take-away is that the matched full-coverage endpoint is the defensible baseline, and the off-diagonal cells explain why the sample ladder matters.",
                "",
            ]
        )
        + "\n"
    )


def render(output_path: Path = OUT_PNG) -> Path:
    configure_fonts()
    rows = build_story_rows()
    specs = metric_specs()
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
        "Full sample coverage stabilizes the BDT",
        size=34.0,
        role="title",
        va="top",
        weight="bold",
        title_axis_align="left",
    )
    rounded_box(fig, 0.040, 0.820, 0.920, 0.064, face=HEADER_BLUE, edge="#a8c5e8", lw=1.0, radius=0.006, zorder=1)
    add_text(
        fig,
        nodes,
        "reading_rule",
        0.058,
        0.856,
        "Each box compares limited training coverage with the matched full-coverage endpoint on the same full validation row.",
        size=14.7,
        weight="bold",
        color=INK,
    )
    add_text(
        fig,
        nodes,
        "reading_caveat",
        0.058,
        0.829,
        "Green boxes are improvements; BDT median gap is included because it explains the visible score separation in the source matrices.",
        size=12.7,
        color=MUTED,
    )

    label_x, label_w = 0.048, 0.252
    metric_x, metric_w, metric_gap = 0.322, 0.145, 0.018
    row_h, row_gap = 0.101, 0.021
    row_tops = [0.735, 0.613, 0.491, 0.369]
    header_y = 0.792

    for idx, spec in enumerate(specs):
        x = metric_x + idx * (metric_w + metric_gap)
        add_text(
            fig,
            nodes,
            f"metric_header_{idx}",
            x + metric_w / 2.0,
            header_y,
            str(spec["header"]),
            size=14.0,
            ha="center",
            weight="bold",
            color=INK,
        )
        add_text(
            fig,
            nodes,
            f"metric_subheader_{idx}",
            x + metric_w / 2.0,
            header_y - 0.025,
            str(spec["subheader"]),
            size=11.1,
            ha="center",
            color=MUTED,
        )

    for row_idx, row in enumerate(rows):
        y = row_tops[row_idx] - row_h
        rounded_box(fig, label_x, y, label_w, row_h, face=blend(row.color, 0.10), edge=row.color, lw=1.2, radius=0.006, zorder=2)
        fig.patches.append(Rectangle((label_x + 0.010, y + 0.018), 0.006, row_h - 0.036, transform=fig.transFigure, facecolor=row.color, edgecolor="none", zorder=4))
        add_text(
            fig,
            nodes,
            f"{row.key}_title",
            label_x + 0.024,
            y + row_h * 0.73,
            row.title,
            size=15.4,
            weight="bold",
            color=INK,
        )
        add_text(
            fig,
            nodes,
            f"{row.key}_message",
            label_x + 0.024,
            y + row_h * 0.51,
            row.message,
            size=12.3,
            weight="bold",
            color=row.color,
        )
        add_text(
            fig,
            nodes,
            f"{row.key}_comparison",
            label_x + 0.024,
            y + row_h * 0.30,
            row.comparison,
            size=10.9,
            color=INK,
        )
        add_text(
            fig,
            nodes,
            f"{row.key}_validation",
            label_x + 0.024,
            y + row_h * 0.13,
            row.validation,
            size=10.8,
            color=MUTED,
        )

        for metric_idx, spec in enumerate(specs):
            x = metric_x + metric_idx * (metric_w + metric_gap)
            draw_metric_cell(fig, nodes, row, spec, x, y, metric_w, row_h)

    warning_y, warning_h = 0.058, 0.142
    rounded_box(fig, 0.040, warning_y, 0.920, warning_h, face=WARNING_FILL, edge=WARNING_EDGE, lw=1.2, radius=0.007, zorder=1)
    add_text(
        fig,
        nodes,
        "warning_title",
        0.058,
        warning_y + warning_h * 0.73,
        "Off-diagonal boxes are stress tests, not production endpoints",
        size=15.2,
        weight="bold",
        color=INK,
    )
    add_text(
        fig,
        nodes,
        "warning_summary",
        0.058,
        warning_y + warning_h * 0.41,
        "Matched full-coverage endpoint = baseline.\nMismatched cells explain the ladder.",
        size=12.4,
        color=MUTED,
        linespacing=1.12,
    )
    for idx, item in enumerate(stress_cells()):
        x0 = 0.492 + idx * 0.236
        metrics = item["metrics"]
        assert isinstance(metrics, Metrics)
        rounded_box(fig, x0, warning_y + 0.026, 0.216, 0.086, face="#fffafa", edge=RED, lw=1.05, radius=0.006, zorder=2)
        fig.patches.append(Rectangle((x0 + 0.012, warning_y + 0.043), 0.005, 0.052, transform=fig.transFigure, facecolor=RED, edgecolor="none", zorder=5))
        add_text(
            fig,
            nodes,
            f"stress_{idx}_label",
            x0 + 0.024,
            warning_y + 0.091,
            str(item["label"]),
            size=11.4,
            weight="bold",
            color=INK,
        )
        add_text(
            fig,
            nodes,
            f"stress_{idx}_text",
            x0 + 0.024,
            warning_y + 0.068,
            "Photon12+20 train on Photon12 val" if idx == 0 else "Photon5+10 train on full signal val",
            size=10.1,
            color=MUTED,
        )
        add_text(
            fig,
            nodes,
            f"stress_{idx}_metrics",
            x0 + 0.024,
            warning_y + 0.044,
            f"AUC {metrics.auc:.3f}, gap {metrics.median_gap:.2f}, WP80 {metrics.wp80_pass_percent:.1f}%",
            size=10.9,
            weight="bold",
            color=RED,
        )

    datapoints_path = output_path.with_name(output_path.stem + "_datapoints.csv")
    script_path = output_path.with_name(output_path.stem + "_speaker_script.md")
    write_datapoints(rows, datapoints_path)
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
            "name": "warning_box",
            "kind": "box",
            "role": "layout",
            "bbox": [
                0.040 * 2560,
                (1.0 - (warning_y + warning_h)) * 1440,
                (0.040 + 0.920) * 2560,
                (1.0 - warning_y) * 1440,
            ],
        }
    )

    layout_path = output_path.with_suffix(".layout_nodes.json")
    layout_path.write_text(
        json.dumps(
            {
                "schema": "slide_layout_nodes_v1",
                "title_axis_x": 0.040 * 2560,
                "minimum_audience_font_px": font_px(9.8),
                "minimum_title_font_px": font_px(28.0),
                "minimum_plot_annotation_font_px": font_px(9.8),
                "vertical_margin_balance": {
                    "top_node": "title",
                    "bottom_node": "warning_box",
                    "target_gap_px": 80.0,
                    "tolerance_px": 30.0,
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
                "schema": "the37_slides25_28_boxchart_story_summary_v1",
                "png": str(output_path),
                "datapoints_csv": str(datapoints_path),
                "layout_nodes": str(layout_path),
                "speaker_script": str(script_path),
                "deck_mutated": False,
                "source_deck": {
                    "title": "PPG19_6_29_26",
                    "presentation_id": "13G56p6lZl05W3GUtUdKHzh5zcSDRdFjs0UDi0sLiLMM",
                    "slides_summarized": [25, 26, 27, 28],
                    "slide_object_ids": {
                        "25": "g3efaeb2c01d_0_529",
                        "26": "g3efaeb2c01d_0_534",
                        "27": "g3efaeb2c01d_0_539",
                        "28": "g3efaeb2c01d_0_544",
                    },
                },
                "inputs": {
                    "slide25_the8_metric_table": str(THE8_METRIC_TABLE),
                    "slide26_pp_inclusive_trainval_payload": str(PP_INCLUSIVE_TRAINVAL_PAYLOAD),
                    "slide27_auau_signal_trainval_payload": str(AUAU_SIGNAL_TRAINVAL_PAYLOAD),
                    "slide28_pp_signal_trainval_payload": str(PP_SIGNAL_TRAINVAL_PAYLOAD),
                },
                "comparison_rows": [
                    {
                        "key": row.key,
                        "source_slide": row.source_slide,
                        "title": row.title,
                        "comparison": row.comparison,
                        "validation": row.validation,
                        "start_label": row.start_label,
                        "end_label": row.end_label,
                        "deltas": row.deltas(),
                    }
                    for row in rows
                ],
                "off_diagonal_note": "Stress cells are intentionally shown as diagnostics, not as recommended endpoint models.",
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
