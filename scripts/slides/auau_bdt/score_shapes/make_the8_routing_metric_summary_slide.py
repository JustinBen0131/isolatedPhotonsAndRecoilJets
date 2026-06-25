#!/usr/bin/env python3
"""Render a compact metric-only summary for corrected THE8 routing variants."""

from __future__ import annotations

import csv
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
    / "routing_metric_summary_20260622"
)
CENT_PAYLOAD = (
    REPO
    / "dataOutput/auauTightBDTValidation/THE59_model_comparison_remakes_20260620"
    / "centrality_routing_20260622/the8_centrality_routing_fine_rows_payload_baseline_holdout.json"
)
ET_PAYLOAD = (
    REPO
    / "dataOutput/auauTightBDTValidation/THE59_model_comparison_remakes_20260620"
    / "et_routing_20260622/the8_et_routing_score_payload_baseline_holdout.json"
)
DEFAULT_OUTPUT = OUT_DIR / "the8_routing_metric_summary_baseline_holdout_v1.png"

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
PANEL = "#f8fafc"
GLOBAL = "#5265d6"
COARSE = "#0f9f8e"
FINE = "#8a62c4"
ET = "#d97706"
READOUT_FILL = "#fff7e6"
READOUT_EDGE = "#e7b85d"

VARIANTS = {
    "global": ("Global", GLOBAL, "o"),
    "coarse": ("3 cent BDTs", COARSE, "s"),
    "fine": ("7 cent BDTs", FINE, "^"),
    "et": ("8 ET BDTs", ET, "D"),
}
METRICS = [
    ("auc", "AUC", "higher is better", "{:.3f}", "higher"),
    ("wp80_inclusive_fake", "Fake rate", "inclusive MC above 80% signal-efficiency threshold", "{:.1%}", "lower"),
    ("logloss", "Logloss", "lower is better", "{:.3f}", "lower"),
    ("median_gap", "Median score gap", "signal median minus inclusive median", "{:.2f}", "higher"),
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
    edge: str,
    lw: float = 1.0,
    radius: float = 0.010,
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


def cell_map(payload: dict) -> dict[tuple[str, str], dict]:
    return {(cell["row"], cell["col"]): cell for cell in payload["cells"]}


def values(payload: dict, cols: list[str], metric: str) -> dict[str, list[float]]:
    cells = cell_map(payload)
    return {col: [float(cells[(row, col)]["metrics"][metric]) for row in payload["row_order"]] for col in cols}


def write_datapoints(cent: dict, et_payload: dict, output_csv: Path) -> None:
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    for axis, payload, cols in [
        ("centrality", cent, ["global", "coarse", "fine"]),
        ("et", et_payload, ["global", "et"]),
    ]:
        cells = cell_map(payload)
        for row in payload["row_order"]:
            for col in cols:
                metrics = cells[(row, col)]["metrics"]
                rows.append(
                    {
                        "axis": axis,
                        "bin": row,
                        "variant_key": col,
                        "variant_label": VARIANTS[col][0],
                        "auc": metrics["auc"],
                        "logloss": metrics["logloss"],
                        "median_gap": metrics["median_gap"],
                        "fake_rate": metrics["wp80_inclusive_fake"],
                    }
                )
    with output_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def metric_limits(series: dict[str, list[float]], metric: str) -> tuple[float, float]:
    vals = [v for arr in series.values() for v in arr]
    lo, hi = min(vals), max(vals)
    if metric == "wp80_inclusive_fake":
        pad = max(0.004, (hi - lo) * 0.22)
    elif metric == "auc":
        pad = max(0.0015, (hi - lo) * 0.28)
    else:
        pad = max(0.002, (hi - lo) * 0.24)
    return lo - pad, hi + pad


def setup_axis(ax: plt.Axes, title: str, note: str, xlabels: list[str], *, show_xlabel: bool) -> None:
    ax.set_facecolor("white")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.spines["left"].set_color("#94a3b8")
    ax.spines["bottom"].set_color("#94a3b8")
    ax.grid(True, color=GRID, linewidth=0.8, alpha=0.95)
    ax.tick_params(labelsize=8.8, colors=MUTED, length=3.2, width=0.8, pad=2)
    ax.set_title(title, fontsize=12.2, fontweight="bold", loc="left", pad=3, color=INK)
    ax.text(0.99, 1.03, note, transform=ax.transAxes, ha="right", va="bottom", fontsize=8.6, color=MUTED)
    ax.set_xticks(range(len(xlabels)))
    ax.set_xticklabels(xlabels if show_xlabel else ["" for _ in xlabels], rotation=0)


def plot_metric_block(
    fig: plt.Figure,
    payload: dict,
    cols: list[str],
    *,
    x: float,
    y: float,
    w: float,
    h: float,
    title: str,
    x_label: str,
    header_color: str,
    nodes: list[tuple[str, object, float, str, dict]],
) -> None:
    rounded_box(fig, x, y, w, h, face=PANEL, edge=EDGE, lw=1.0, radius=0.012, zorder=0)
    fig.patches.append(Rectangle((x, y + h - 0.041), w, 0.041, transform=fig.transFigure, facecolor="white", edgecolor="none", zorder=1))
    fig.patches.append(Rectangle((x, y + h - 0.041), 0.008, 0.041, transform=fig.transFigure, facecolor=header_color, edgecolor=header_color, zorder=2))
    add_text(fig, nodes, f"{title} header", x + 0.019, y + h - 0.020, title, size=15.8, weight="bold")
    add_text(fig, nodes, f"{title} x label", x + w - 0.018, y + h - 0.020, x_label, size=10.7, color=MUTED, ha="right")

    left = x + 0.040
    right = x + w - 0.018
    top = y + h - 0.070
    bottom = y + 0.034
    row_gap = 0.024
    ax_h = (top - bottom - 3 * row_gap) / 4.0
    xlabels = payload["row_order"]

    for idx, (metric, label, note, fmt, direction) in enumerate(METRICS):
        ay = top - (idx + 1) * ax_h - idx * row_gap
        ax = fig.add_axes([left, ay, right - left, ax_h], zorder=3)
        series = values(payload, cols, metric)
        setup_axis(ax, label, note, xlabels, show_xlabel=(idx == len(METRICS) - 1))
        ax.set_ylim(*metric_limits(series, metric))
        for col in cols:
            _, color, marker = VARIANTS[col]
            ax.plot(
                range(len(xlabels)),
                series[col],
                color=color,
                lw=2.1,
                marker=marker,
                markersize=4.8,
                label=VARIANTS[col][0],
            )
        if metric == "wp80_inclusive_fake":
            ax.yaxis.set_major_formatter(lambda val, pos: f"{100.0 * val:.0f}%")
        elif metric == "median_gap":
            ax.yaxis.set_major_formatter(lambda val, pos: f"{val:.2f}")
        else:
            ax.yaxis.set_major_formatter(lambda val, pos: fmt.format(val))
        if idx == 0:
            leg = ax.legend(
                loc="upper left",
                bbox_to_anchor=(0.01, 0.98),
                frameon=True,
                facecolor="white",
                edgecolor=EDGE,
                fontsize=8.8,
                ncol=len(cols),
                handlelength=1.7,
                borderpad=0.35,
                columnspacing=1.0,
            )
            leg.get_frame().set_alpha(0.94)


def centrality_readout(payload: dict) -> str:
    cells = cell_map(payload)
    auc_deltas: list[float] = []
    fake_deltas: list[float] = []
    for row in payload["row_order"]:
        global_metrics = cells[(row, "global")]["metrics"]
        for col in ("coarse", "fine"):
            metrics = cells[(row, col)]["metrics"]
            auc_deltas.append(metrics["auc"] - global_metrics["auc"])
            fake_deltas.append(100.0 * (metrics["wp80_inclusive_fake"] - global_metrics["wp80_inclusive_fake"]))
    return (
        f"Centrality-only routing gives small AUC gains ({min(auc_deltas):+.3f} to {max(auc_deltas):+.3f}) "
        f"but the fake-rate change is mixed ({min(fake_deltas):+.2f} to {max(fake_deltas):+.2f} pp)."
    )


def et_readout(payload: dict) -> str:
    cells = cell_map(payload)
    auc_deltas: list[float] = []
    fake_deltas: list[float] = []
    logloss_deltas: list[float] = []
    wins = 0
    for row in payload["row_order"]:
        global_metrics = cells[(row, "global")]["metrics"]
        metrics = cells[(row, "et")]["metrics"]
        d_auc = metrics["auc"] - global_metrics["auc"]
        d_fake = 100.0 * (metrics["wp80_inclusive_fake"] - global_metrics["wp80_inclusive_fake"])
        d_logloss = metrics["logloss"] - global_metrics["logloss"]
        auc_deltas.append(d_auc)
        fake_deltas.append(d_fake)
        logloss_deltas.append(d_logloss)
        if d_auc > 0 and d_fake < 0 and d_logloss < 0:
            wins += 1
    best_fake_row = payload["row_order"][min(range(len(fake_deltas)), key=lambda i: fake_deltas[i])]
    return (
        f"ET routing is cleaner. AUC/logloss/fake-rate all improve in {wins}/8 ET bins; "
        f"fake-rate drops {min(fake_deltas):+.2f} to {max(fake_deltas):+.2f} pp, largest at {best_fake_row}."
    )


def render(output_path: Path = DEFAULT_OUTPUT) -> Path:
    configure_fonts()
    cent = json.loads(CENT_PAYLOAD.read_text())
    et_payload = json.loads(ET_PAYLOAD.read_text())

    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI, facecolor="white")
    nodes: list[tuple[str, object, float, str, dict]] = []

    add_text(
        fig,
        nodes,
        "title",
        0.040,
        0.943,
        "Routed BDT metrics on one fixed holdout",
        size=40.0,
        role="title",
        va="top",
        weight="bold",
        title_anchor=True,
    )
    rounded_box(fig, 0.040, 0.835, 0.922, 0.052, face="white", edge=EDGE, lw=0.95, radius=0.006)
    add_text(
        fig,
        nodes,
        "holdout contract",
        0.058,
        0.863,
        "Fixed holdout = 10% row holdout, 15-35 GeV, 0-80%",
        size=12.2,
        color=INK,
        weight="bold",
        title_band_exception=True,
    )
    add_text(
        fig,
        nodes,
        "sample contract",
        0.502,
        0.863,
        "Signal Photon12+20; inclusive Jet12+20+30+40, no truth filter",
        size=11.4,
        color=MUTED,
        title_band_exception=True,
    )

    plot_metric_block(
        fig,
        cent,
        ["global", "coarse", "fine"],
        x=0.040,
        y=0.170,
        w=0.452,
        h=0.632,
        title="Centrality routing",
        x_label="fine centrality bins",
        header_color=FINE,
        nodes=nodes,
    )
    plot_metric_block(
        fig,
        et_payload,
        ["global", "et"],
        x=0.515,
        y=0.170,
        w=0.447,
        h=0.632,
        title="ET routing",
        x_label="cluster ET bins",
        header_color=ET,
        nodes=nodes,
    )

    rounded_box(fig, 0.040, 0.050, 0.922, 0.094, face=READOUT_FILL, edge=READOUT_EDGE, lw=1.05, radius=0.008)
    add_text(fig, nodes, "readout label", 0.058, 0.103, "Readout", size=17.2, weight="bold")
    add_text(fig, nodes, "readout cent", 0.145, 0.113, centrality_readout(cent), size=12.2, color=INK)
    add_text(fig, nodes, "readout et", 0.145, 0.081, et_readout(et_payload), size=12.2, color=INK)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    datapoints_path = output_path.with_name(output_path.stem + "_datapoints.csv")
    write_datapoints(cent, et_payload, datapoints_path)
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
                "minimum_title_font_px": font_px(30.0),
                "minimum_plot_annotation_font_px": font_px(8.5),
                "nodes": layout_nodes,
            },
            indent=2,
        )
        + "\n"
    )
    manifest_path = output_path.with_suffix(".manifest.json")
    manifest = {
        "schema": "the8_routing_metric_summary_slide_v1",
        "png": str(output_path),
        "centrality_payload": str(CENT_PAYLOAD),
        "et_payload": str(ET_PAYLOAD),
        "datapoints_csv": str(datapoints_path),
        "layout_nodes": str(layout_path),
        "campaign_tag": cent.get("campaign_tag"),
        "selection": cent.get("selection"),
        "labels": cent.get("labels"),
        "centrality_model_columns": cent.get("model_columns"),
        "et_model_columns": et_payload.get("model_columns"),
        "centrality_readout": centrality_readout(cent),
        "et_readout": et_readout(et_payload),
        "scope_note": "Current corrected same-holdout payloads cover global, base14_perCent3, base14_perCent7, and base14_perEt. ET x centrality routed payloads are not mixed in here until scored on the same holdout contract.",
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    script_path = output_path.with_name(output_path.stem + "_speaker_script.md")
    script_path.write_text(
        "\n".join(
            [
                "# Metric view: routed BDTs over the same holdout",
                "",
                "This slide removes the score-shape grid and plots the four scalar diagnostics directly over the physical bins.",
                "",
                "The key contract is that every point is evaluated on the same reconstructed baseline 10% row holdout. Signal is embedded Photon12+20. Inclusive MC is embedded Jet12+20+30+40 with no truth-background filter.",
                "",
                f"{centrality_readout(cent)}",
                "",
                f"{et_readout(et_payload)}",
                "",
                "The missing comparison for a later slide is ET x centrality routing on this same corrected holdout contract. I would not mix the older ptCent7 validation tables into this slide because that would change both the model family and the validation contract.",
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
    print(out.with_name(out.stem + "_datapoints.csv"))
    print(out.with_name(out.stem + "_speaker_script.md"))


if __name__ == "__main__":
    main()
