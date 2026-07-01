#!/usr/bin/env python3
"""Build a THE-79 expanded-BDT pT/centrality AUC slide candidate."""

from __future__ import annotations

import csv
import json
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager
from matplotlib.patches import FancyBboxPatch, Patch, Rectangle


THIS_FILE = Path(__file__).resolve()
REPO = next((p for p in THIS_FILE.parents if (p / "AGENTS.md").exists()), THIS_FILE.parents[4])
SCRIPTS = REPO / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.append(str(SCRIPTS))

from slides.common.slide_defaults import SLIDE_DPI, SLIDE_HEIGHT_PX, SLIDE_WIDTH_PX, slide_figsize  # noqa: E402


REPORT = REPO / "dataOutput/auauTightBDTValidation/model_validation_condor_20260509_192942"
SOURCE_TABLE = REPORT / "validation_auc_table.csv"
OUT_DIR = REPO / "dataOutput/auauBDTModelSelection/THE79_expanded_bdt_cent3_pt_auc_20260627"
OUT_PNG = OUT_DIR / "the79_expanded_bdt_cent3_pt_auc_slide.png"
OUT_CSV = OUT_DIR / "the79_expanded_bdt_cent3_pt_auc_values.csv"
OUT_MANIFEST = OUT_DIR / "the79_expanded_bdt_cent3_pt_auc_manifest.json"
OUT_SCRIPT = OUT_DIR / "the79_expanded_bdt_cent3_pt_auc_speaker_script.md"
OUT_LAYOUT = OUT_DIR / "the79_expanded_bdt_cent3_pt_auc_layout_nodes.json"
OUT_AUDIT = OUT_DIR / "the79_expanded_bdt_cent3_pt_auc_slide_audit.json"

TIMES_DIR = Path("/System/Library/Fonts/Supplemental")
TIMES_FONTS = [
    TIMES_DIR / "Times New Roman.ttf",
    TIMES_DIR / "Times New Roman Bold.ttf",
    TIMES_DIR / "Times New Roman Italic.ttf",
    TIMES_DIR / "Times New Roman Bold Italic.ttf",
]

INK = "#111827"
MUTED = "#475569"
GRID = "#dbe4ef"
EDGE = "#c8d3e2"
PANEL_FILL = "#fbfdff"
NOTE_FILL = "#fff7ed"
NOTE_EDGE = "#f59e0b"
GREEN = "#15803d"


@dataclass(frozen=True)
class ProductSpec:
    product: str
    label: str
    color: str
    count: str


PRODUCTS = [
    ProductSpec("centAsFeat_pt5to40", "Global C", "#6b7280", "1 BDT"),
    ProductSpec("centDepBDTs_pt5to40", "C3", "#7c3aed", "3 BDTs"),
    ProductSpec("centDepFineBDTs_pt5to40", "C7", "#a855f7", "7 BDTs"),
    ProductSpec("ptBinCentAsFeat", r"$E_T$", "#0284c7", "5 pT bins"),
    ProductSpec("ptCentDep3", r"$E_T$×C3", "#f59e0b", "pT × 3 cent"),
    ProductSpec("ptCentDepFine", r"$E_T$×C7", "#15803d", "pT × 7 cent"),
]

CENT_KEYS = ["0_20", "20_50", "50_80"]
CENT_LABELS = {"0_20": "0-20%", "20_50": "20-50%", "50_80": "50-80%"}
PT_KEYS = ["6_10", "10_15", "15_20", "20_25", "25_35"]
PT_LABELS = {"6_10": "6-10", "10_15": "10-15", "15_20": "15-20", "20_25": "20-25", "25_35": "25-35"}


def setup_style() -> None:
    for font in TIMES_FONTS:
        if font.exists():
            font_manager.fontManager.addfont(str(font))
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "custom",
            "mathtext.rm": "Times New Roman",
            "mathtext.it": "Times New Roman:italic",
            "mathtext.bf": "Times New Roman:bold",
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "axes.facecolor": PANEL_FILL,
            "axes.edgecolor": EDGE,
            "axes.linewidth": 1.0,
            "xtick.color": INK,
            "ytick.color": INK,
            "axes.unicode_minus": False,
        }
    )


def font_px(points: float) -> float:
    return points * SLIDE_DPI / 72.0


def text_bbox_px(fig: plt.Figure, artist) -> list[float]:
    fig.canvas.draw()
    bbox = artist.get_window_extent(renderer=fig.canvas.get_renderer())
    scale_x = fig.bbox.width / SLIDE_WIDTH_PX
    scale_y = fig.bbox.height / SLIDE_HEIGHT_PX
    height = fig.bbox.height
    return [
        float(bbox.x0 / scale_x),
        float((height - bbox.y1) / scale_y),
        float(bbox.x1 / scale_x),
        float((height - bbox.y0) / scale_y),
    ]


def add_text(
    fig: plt.Figure,
    nodes: list[dict],
    name: str,
    x: float,
    y: float,
    text: str,
    *,
    size: float,
    role: str = "audience",
    weight: str = "normal",
    color: str = INK,
    ha: str = "left",
    va: str = "center",
    linespacing: float = 1.0,
    **flags,
) -> None:
    artist = fig.text(
        x,
        y,
        text,
        fontsize=size,
        fontweight=weight,
        color=color,
        ha=ha,
        va=va,
        linespacing=linespacing,
    )
    nodes.append(
        {
            "name": name,
            "kind": "text",
            "role": role,
            "text": text,
            "font_px": font_px(size),
            "bbox": text_bbox_px(fig, artist),
            "title_anchor": role == "title",
            **flags,
        }
    )


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
            zorder=4,
        )
    )


def load_values() -> pd.DataFrame:
    df = pd.read_csv(SOURCE_TABLE)
    want = {p.product for p in PRODUCTS}
    sub = df[df["product"].isin(want) & df["centrality_bin"].isin(CENT_KEYS) & df["pt_bin"].isin(PT_KEYS)].copy()
    missing = []
    for product in want:
        for cent in CENT_KEYS:
            for pt in PT_KEYS:
                if sub[(sub["product"] == product) & (sub["centrality_bin"] == cent) & (sub["pt_bin"] == pt)].empty:
                    missing.append((product, cent, pt))
    if missing:
        raise RuntimeError(f"missing AUC cells: {missing[:8]} ... total={len(missing)}")
    sub["product_label"] = sub["product"].map({p.product: p.label for p in PRODUCTS})
    sub["centrality_label"] = sub["centrality_bin"].map(CENT_LABELS)
    sub["pt_label"] = sub["pt_bin"].map(PT_LABELS)
    return sub


def matrix_for(df: pd.DataFrame, cent: str) -> np.ndarray:
    out = np.full((len(PRODUCTS), len(PT_KEYS)), np.nan)
    for i, spec in enumerate(PRODUCTS):
        for j, pt in enumerate(PT_KEYS):
            row = df[(df["product"] == spec.product) & (df["centrality_bin"] == cent) & (df["pt_bin"] == pt)].iloc[0]
            out[i, j] = float(row["auc"])
    return out


def draw_panel(ax: plt.Axes, values: np.ndarray, cent_label: str, *, show_xticks: bool) -> None:
    x = np.arange(len(PT_KEYS))
    n = len(PRODUCTS)
    width = 0.78 / n
    offsets = (np.arange(n) - (n - 1) / 2.0) * width
    for j in range(len(PT_KEYS)):
        if j % 2 == 0:
            ax.axvspan(j - 0.5, j + 0.5, color="#eef4fb", alpha=0.55, zorder=0)
    for i, spec in enumerate(PRODUCTS):
        bars = ax.bar(x + offsets[i], values[i], width=width * 0.90, color=spec.color, label=spec.label, zorder=3)
        for bar, val in zip(bars, values[i]):
            if val >= 0.790 or (cent_label == "0-20%" and val >= 0.750):
                ax.text(
                    bar.get_x() + bar.get_width() / 2.0,
                    val + 0.004,
                    f"{val:.2f}",
                    ha="center",
                    va="bottom",
                    fontsize=8.4,
                    color=INK,
                    fontweight="bold",
                    rotation=90,
                    zorder=4,
                )
    ax.set_title(f"{cent_label} centrality | ROC AUC", loc="left", fontsize=15.6, fontweight="bold", pad=8, color=INK)
    ax.set_ylim(0.60, 0.91)
    ax.set_yticks([0.60, 0.70, 0.80, 0.90])
    ax.grid(axis="y", color=GRID, linewidth=1.0, zorder=1)
    ax.tick_params(axis="y", labelsize=12.5)
    ax.tick_params(axis="x", labelsize=14.0, pad=5)
    for spine in ax.spines.values():
        spine.set_edgecolor(EDGE)
        spine.set_linewidth(1.0)
    ax.set_xticks(x)
    if show_xticks:
        ax.set_xticklabels([PT_LABELS[k] for k in PT_KEYS], fontsize=14.0, fontweight="bold")
        ax.set_xlabel(r"cluster $E_T$ bin [GeV]", fontsize=15.2, fontweight="bold", labelpad=8, color=INK)
    else:
        ax.set_xticklabels([])


def write_outputs(df: pd.DataFrame, summary: dict) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    ordered = df.sort_values(
        by=["centrality_bin", "pt_bin", "product"],
        key=lambda s: s.map({**{c: i for i, c in enumerate(CENT_KEYS)}, **{p: i for i, p in enumerate(PT_KEYS)}, **{spec.product: i for i, spec in enumerate(PRODUCTS)}}).fillna(99)
        if s.name in {"centrality_bin", "pt_bin", "product"}
        else s,
    )
    with OUT_CSV.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["product", "product_label", "centrality_bin", "pt_bin", "auc", "entries", "signal_entries", "background_entries"])
        for row in ordered.itertuples():
            writer.writerow([row.product, row.product_label, row.centrality_label, row.pt_label, row.auc, row.entries, row.signal_entries, row.background_entries])
    OUT_MANIFEST.write_text(
        json.dumps(
            {
                "schema": "THE79_EXPANDED_BDT_CENT3_PT_AUC_SLIDE_V1",
                "png": str(OUT_PNG),
                "source_table": str(SOURCE_TABLE),
                "source_report": str(REPORT),
                "slide_40_reference": {
                    "deck": "WP_GammaJets_6_10_26",
                    "presentation_id": "167x-He2rOOBO2i4nNS6Pdcqu7Wv03GeFMuWH9tRRx-8",
                    "slide_object_id": "g3ee13a1139c_4_77",
                    "source_generator": "scripts/slides/auau_bdt/score_shapes/make_the8_routing_bar_summary_slides.py",
                    "source_png": "dataOutput/auauTightBDTValidation/THE59_model_comparison_remakes_20260620/routing_bar_summary_20260622/the8_routing_bar_auc_et_by_coarse_centrality_v1.png",
                },
                "metric": "absolute ROC AUC; higher is better",
                "pt_bins": [PT_LABELS[k] for k in PT_KEYS],
                "centrality_bins": [CENT_LABELS[k] for k in CENT_KEYS],
                "products": [spec.__dict__ for spec in PRODUCTS],
                "summary": summary,
            },
            indent=2,
        )
        + "\n"
    )
    OUT_SCRIPT.write_text(
        "\n".join(
            [
                "# Speaker Script: Expanded BDT pT and centrality performance",
                "",
                "This is the first low-pT view of the expanded AuAu BDT family, using the same visual grammar as slide 40 but plotting absolute AUC instead of improvement.",
                "",
                "Each row is one of the three coarse centrality bins. The x axis is the cluster E_T validation bin. The first available local bin is 6 to 10 GeV, so this is the low-pT check we can make from the existing pulled report without pretending we have a separate 5 to 6 bin.",
                "",
                "The main pattern is that centrality-only routing is already competitive in the lowest E_T bin. Once we are at 15 GeV and above, the pT by centrality routed model, especially E_T times C7, becomes the cleanest performer in nearly every cell.",
                "",
                f"Quantitatively, E_T times C7 is the best model in {summary['et_c7_win_count']} of {summary['cell_count']} centrality-by-E_T cells. C7 wins {summary['c7_win_count']} low-E_T cells, and E_T-only wins {summary['et_win_count']} cell.",
                "",
                "So yes, it is plausible and useful to check a higher-threshold trained model below 15 GeV, but we should label that as extrapolation. The existing expanded-BDT evidence says the low-E_T region should not be assumed to behave like the 15 to 35 GeV region.",
                "",
                "Generated artifacts:",
                f"- {OUT_PNG}",
                f"- {OUT_CSV}",
                f"- {OUT_MANIFEST}",
                "",
            ]
        )
    )


def render() -> None:
    setup_style()
    df = load_values()
    win_counts = {"ptCentDepFine": 0, "centDepFineBDTs_pt5to40": 0, "ptBinCentAsFeat": 0}
    cell_count = 0
    for cent in CENT_KEYS:
        for pt in PT_KEYS:
            cell = df[(df["centrality_bin"] == cent) & (df["pt_bin"] == pt)]
            best = cell.loc[cell["auc"].idxmax(), "product"]
            if best in win_counts:
                win_counts[best] += 1
            cell_count += 1
    summary = {
        "cell_count": cell_count,
        "et_c7_win_count": win_counts["ptCentDepFine"],
        "c7_win_count": win_counts["centDepFineBDTs_pt5to40"],
        "et_win_count": win_counts["ptBinCentAsFeat"],
    }

    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI, facecolor="white")
    nodes: list[dict] = []
    add_text(
        fig,
        nodes,
        "title",
        0.040,
        0.944,
        r"Expanded BDT performance down to low $E_T$",
        size=34.5,
        role="title",
        weight="bold",
        va="top",
    )
    add_text(
        fig,
        nodes,
        "reading line",
        0.043,
        0.842,
        r"Absolute ROC AUC in each centrality and cluster-$E_T$ bin; higher bars mean cleaner truth-photon/background ranking.",
        size=16.5,
        weight="bold",
        color=INK,
        title_band_exception=True,
    )

    rounded_box(fig, 0.040, 0.772, 0.555, 0.064, face="white", edge="#ef4444", lw=1.4)
    add_text(fig, nodes, "scope label", 0.056, 0.816, "EXPANDED 5-40 / ALL-RANGE VALIDATION", size=12.4, weight="bold", color="#991b1b")
    add_text(
        fig,
        nodes,
        "scope body",
        0.056,
        0.791,
        r"first local low-$E_T$ bin is 6-10 GeV; this slide does not invent a separate 5-6 bin",
        size=13.9,
        weight="bold",
        color="#991b1b",
    )

    legend_handles = [Patch(facecolor=spec.color, edgecolor="none", label=f"{spec.label} ({spec.count})") for spec in PRODUCTS]
    fig.legend(
        handles=legend_handles,
        loc="upper left",
        bbox_to_anchor=(0.610, 0.838),
        ncol=3,
        frameon=False,
        fontsize=12.8,
        columnspacing=1.4,
        handlelength=1.7,
        handletextpad=0.5,
    )

    positions = [
        [0.065, 0.555, 0.885, 0.155],
        [0.065, 0.355, 0.885, 0.155],
        [0.065, 0.155, 0.885, 0.155],
    ]
    for idx, (cent, pos) in enumerate(zip(CENT_KEYS, positions)):
        ax = fig.add_axes(pos)
        draw_panel(ax, matrix_for(df, cent), CENT_LABELS[cent], show_xticks=idx == 2)

    rounded_box(fig, 0.040, 0.035, 0.920, 0.090, face=NOTE_FILL, edge=NOTE_EDGE, lw=1.1)
    fig.patches.append(Rectangle((0.057, 0.052), 0.017, 0.056, transform=fig.transFigure, facecolor=GREEN, edgecolor="none", zorder=8))
    add_text(fig, nodes, "readout lead", 0.092, 0.095, r"Readout", size=18.3, weight="bold")
    add_text(
        fig,
        nodes,
        "readout body",
        0.207,
        0.095,
        r"$E_T$×C7 wins 11/15 cells; from 15-35 GeV it is the consistent best route.",
        size=17.0,
        color=INK,
    )
    add_text(
        fig,
        nodes,
        "readout body second",
        0.207,
        0.061,
        r"In 6-10 GeV, C7/global-centrality are competitive, so low-$E_T$ extrapolation must be labeled.",
        size=15.6,
        color=MUTED,
    )

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=SLIDE_DPI)
    OUT_LAYOUT.write_text(
        json.dumps(
            {
                "schema": "slide_layout_nodes_v1",
                "title_axis_x": 0.040 * SLIDE_WIDTH_PX,
                "minimum_audience_font_px": font_px(12.0),
                "minimum_title_font_px": font_px(34.0),
                "minimum_plot_annotation_font_px": font_px(8.0),
                "nodes": nodes,
            },
            indent=2,
        )
        + "\n"
    )
    plt.close(fig)
    write_outputs(df, summary)

    audit = SCRIPTS / "slides/common/post_render_slide_audit.py"
    subprocess.run(
        [
            sys.executable,
            str(audit),
            "--png",
            str(OUT_PNG),
            "--layout-nodes",
            str(OUT_LAYOUT),
            "--output",
            str(OUT_AUDIT),
        ],
        check=True,
    )
    print(f"Wrote {OUT_PNG}")
    print(f"Wrote {OUT_SCRIPT}")
    print(f"Wrote {OUT_MANIFEST}")
    print(f"Wrote {OUT_AUDIT}")


if __name__ == "__main__":
    render()
