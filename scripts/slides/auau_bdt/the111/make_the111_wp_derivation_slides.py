#!/usr/bin/env python3
"""Render two JSTG-ready THE-111 candidate working-point slides.

The slides are regenerated only from the frozen combined-correction
simulation-validation payload.  THE-112 consumes this model/WP surface in its
sideband discovery scan; it does not rederive these thresholds.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Rectangle
import numpy as np


def find_repo() -> Path:
    for parent in Path(__file__).resolve().parents:
        if (parent / "AGENTS.md").exists() and (parent / "agent_context").exists():
            return parent
    raise RuntimeError("Could not locate ThesisAnalysis root")


REPO = find_repo()
SOURCE = Path(
    "/Users/patsfan753/Desktop/PPG19_IAN_WP/Plots/photon_id/"
    "combined_candidate_wp_derivation/source/the111_combined_candidate_wp.json"
)
OUTDIR = REPO / "dataOutput/slides/the45_jstg_20260720/the111_wp_derivation"
SLIDE11 = OUTDIR / "the111_wp_flat_fits_by_centrality_slide.png"
SLIDE12 = OUTDIR / "the111_wp80_centrality_fit_slide.png"
SLIDE12_CLEANED = OUTDIR / "the111_wp80_centrality_fit_cleaned_slide.png"
SLIDE11_NOTES = OUTDIR / "the111_wp_flat_fits_by_centrality_speaker_notes.md"
SLIDE12_NOTES = OUTDIR / "the111_wp80_centrality_fit_speaker_notes.md"
SLIDE12_CLEANED_NOTES = OUTDIR / "the111_wp80_centrality_fit_cleaned_speaker_notes.md"
CSV_PATH = OUTDIR / "the111_wp_derivation_values.csv"
MANIFEST = OUTDIR / "the111_wp_derivation_slide_manifest.json"

WIDTH, HEIGHT, DPI = 2560, 1440, 160
INK = "#152238"
MUTED = "#56657a"
GRID = "#d8e0ea"
BLUE_BAND = "#f5f9ff"
BLUE_EDGE = "#bad1ec"
# Match the established JSTG bullet treatment used by the BDT-score overlay
# slide: literal DejaVu Sans arrowhead, #2468A8, and a 0.021-wide hanging
# indent.  Keep this as a shared deck convention rather than a near-match.
BLUE_BULLET = "#2468A8"
WP_COLORS = {"WP90": "#d84a9b", "WP80": "#15966b", "WP70": "#3978c5"}
WP_MARKERS = {"WP90": "o", "WP80": "s", "WP70": "^"}
WP_ORDER = ("WP90", "WP80", "WP70")

# Slide 11 audience layout: drop the shaded header panel and the flat-fit
# threshold table, then spend the reclaimed height on larger symmetric panels.
CLEAN_LAYOUT = False


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def configure() -> None:
    plt.rcParams.update(
        {
            "font.family": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
            "axes.linewidth": 1.05,
            "axes.labelsize": 14.0,
            "xtick.labelsize": 11.0,
            "ytick.labelsize": 11.0,
        }
    )


def add_box(fig: plt.Figure, xywh: tuple[float, float, float, float], *, face: str, edge: str, lw: float = 1.1, radius: float = 0.006) -> None:
    fig.add_artist(
        FancyBboxPatch(
            xywh[:2],
            xywh[2],
            xywh[3],
            transform=fig.transFigure,
            boxstyle=f"round,pad=0.005,rounding_size={radius}",
            facecolor=face,
            edgecolor=edge,
            linewidth=lw,
            zorder=-10,
        )
    )


def add_sphenix(ax: plt.Axes) -> None:
    ax.text(
        0.975,
        0.965,
        r"$\bf{\it{sPHENIX}}$ Internal",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=9.0,
        color=INK,
    )


def decorate(ax: plt.Axes) -> None:
    ax.grid(axis="y", color=GRID, alpha=0.76, lw=0.65)
    ax.tick_params(direction="in", top=True, right=True, length=4.2, width=0.9)
    for spine in ax.spines.values():
        spine.set_color(INK)
        spine.set_linewidth(0.9)


def centrality_labels(payload: dict) -> list[str]:
    # Preserve the payload's ASCII-hyphen key for lookups; translate only at
    # draw time so audience-facing slide labels use a proper en dash.
    return [f"{int(lo)}-{int(hi)}%" for lo, hi in zip(payload["cent_edges"][:-1], payload["cent_edges"][1:])]


def display_centrality(label: str) -> str:
    return label.replace("-", "–")


def flat_row(payload: dict, label: str, wp: str) -> dict:
    return next(row for row in payload["flat_rows"] if row["centrality_label"] == label and row["wp_label"] == wp)


def et_rows(payload: dict, label: str, wp: str) -> list[dict]:
    return sorted(
        (row for row in payload["et_cells"] if row["centrality_label"] == label and row["wp_label"] == wp),
        key=lambda row: row["et_center"],
    )


def figure() -> plt.Figure:
    return plt.figure(figsize=(WIDTH / DPI, HEIGHT / DPI), dpi=DPI, facecolor="white")


def title_and_band(
    fig: plt.Figure,
    title: str,
    detail: str,
    *,
    detail_fontsize: float = 14.6,
    title_fontsize: float = 27.0,
    boxed_detail: bool = True,
    plain_detail: bool = False,
    detail_y: float = 0.870,
) -> None:
    fig.text(0.050, 0.940, title, ha="left", va="top", fontsize=title_fontsize, weight="bold", color=INK)
    if plain_detail:
        # Same wording and baseline as the boxed form, with the panel removed
        # and the text pulled back into line with the title.
        detail_x = 0.052
    elif boxed_detail:
        add_box(fig, (0.050, 0.834, 0.900, 0.066), face=BLUE_BAND, edge=BLUE_EDGE)
        detail_x = 0.069
    else:
        # Exact JSTG BDT-score-overlay arrowhead convention.
        detail_x = 0.058
        fig.text(
            detail_x,
            detail_y,
            "▶",
            ha="left",
            # The slide-13 glyph/colour/indent are retained, but this larger
            # single-line subtitle is centred on the arrow tip rather than
            # sharing its top edge.
            va="center",
            fontsize=14.5,
            color=BLUE_BULLET,
            fontfamily="DejaVu Sans",
        )
        detail_x += 0.021
    fig.text(detail_x, detail_y, detail, ha="left", va="center", fontsize=detail_fontsize, color=INK)


def draw_wp_legend(
    fig: plt.Figure,
    x: float,
    y: float,
    *,
    fontsize: float = 12.6,
    markersize: float = 5.1,
    linewidth: float = 1.6,
    columnspacing: float = 1.55,
) -> None:
    handles = [
        Line2D([0], [0], marker=WP_MARKERS[wp], color=WP_COLORS[wp], lw=linewidth, markersize=markersize, label=f"{wp}  {int(100 * {'WP90': .9, 'WP80': .8, 'WP70': .7}[wp])}%")
        for wp in WP_ORDER
    ]
    fig.legend(
        handles=handles,
        loc="center",
        bbox_to_anchor=(x, y),
        ncol=3,
        frameon=False,
        fontsize=fontsize,
        handletextpad=0.45,
        columnspacing=columnspacing,
    )


def render_slide11(payload: dict) -> list[dict]:
    fig = figure()
    title_and_band(
        fig,
        "Flat fits define one BDT threshold per centrality interval",
        r"Corrected 14-feature Au+Au candidate: eight $E_T^{\gamma}$ points in $15\leq E_T^{\gamma}<35$ GeV are compressed to one threshold per centrality interval.",
        plain_detail=CLEAN_LAYOUT,
        detail_fontsize=16.4 if CLEAN_LAYOUT else 14.6,
    )
    if CLEAN_LAYOUT:
        # Larger panel type to match the enlarged axes.
        plt.rcParams.update({"axes.labelsize": 18.5, "xtick.labelsize": 15.0, "ytick.labelsize": 15.0})
        draw_wp_legend(fig, 0.760, 0.795, fontsize=18.0, markersize=9.0, linewidth=2.7, columnspacing=2.1)
        # The scatter about each constant fit runs 12-89x the per-point
        # statistical error (median 29x), so the residual E_T structure is real
        # and the flat fit is stated as interim rather than as a description.
        fig.text(0.052, 0.795, r"$E_T^{\gamma}$ dependence is still under study; flat fits interpolate one threshold per centrality for now.", fontsize=14.6, color=MUTED, va="center")
        # Four over three, bottom row centred on the same column pitch.
        panel_w, panel_h = 0.207, 0.285
        top_y, bottom_y = 0.445, 0.085
        top_x = [0.050, 0.281, 0.512, 0.743]
        bottom_x = [0.1655, 0.3965, 0.6275]
        positions = [(x, top_y) for x in top_x] + [(x, bottom_y) for x in bottom_x]
        title_fs, marker_s = 17.5, 5.6
    else:
        draw_wp_legend(fig, 0.758, 0.789)
        fig.text(0.052, 0.788, "Points show fixed signal-efficiency thresholds; lines are inverse-variance weighted constant fits.", fontsize=12.1, color=MUTED)
        panel_w, panel_h = 0.190, 0.177
        positions = [(0.060, 0.525), (0.282, 0.525), (0.504, 0.525), (0.726, 0.525), (0.171, 0.300), (0.393, 0.300), (0.615, 0.300)]
        title_fs, marker_s = 12.8, 3.3

    labels = centrality_labels(payload)
    axes: list[plt.Axes] = []
    for idx, (label, (x, y)) in enumerate(zip(labels, positions, strict=True)):
        ax = fig.add_axes([x, y, panel_w, panel_h])
        axes.append(ax)
        for wp in WP_ORDER:
            points = et_rows(payload, label, wp)
            ex = np.asarray([row["et_center"] for row in points], dtype=float)
            ey = np.asarray([row["threshold"] for row in points], dtype=float)
            elo = np.asarray([row["threshold_stat_err_low"] for row in points], dtype=float)
            ehi = np.asarray([row["threshold_stat_err_high"] for row in points], dtype=float)
            ax.errorbar(
                ex,
                ey,
                yerr=np.vstack([elo, ehi]),
                fmt=WP_MARKERS[wp],
                color=WP_COLORS[wp],
                markersize=marker_s,
                markeredgecolor="white",
                markeredgewidth=0.45 if not CLEAN_LAYOUT else 0.7,
                capsize=1.5 if not CLEAN_LAYOUT else 2.4,
                elinewidth=0.75 if not CLEAN_LAYOUT else 1.1,
                zorder=3,
            )
            ax.axhline(flat_row(payload, label, wp)["threshold"], color=WP_COLORS[wp], lw=1.35 if not CLEAN_LAYOUT else 2.0, alpha=0.92, zorder=2)
        ax.set_title(display_centrality(label), fontsize=title_fs, weight="bold", color=INK, pad=4 if not CLEAN_LAYOUT else 7)
        ax.set_xlim(14.5, 35.5)
        ax.set_ylim(0.35, 0.78)
        ax.set_xticks((15, 25, 35))
        ax.set_yticks((0.4, 0.5, 0.6, 0.7))
        if idx in (0, 4):
            ax.set_ylabel("BDT score cut", labelpad=3)
        else:
            ax.set_yticklabels([])
        if idx >= 4:
            ax.set_xlabel(r"$E_T^{\gamma}$ [GeV]", labelpad=2)
        else:
            ax.set_xticklabels([])
        decorate(ax)

    if CLEAN_LAYOUT:
        fig.savefig(SLIDE11, facecolor="white")
        plt.close(fig)
        return [
            {"name": "slide11 title", "kind": "text", "role": "title", "title_anchor": True, "title_axis_align": "left", "bbox": [128, 55, 2290, 130], "font_px": 68, "text": "Flat fits define one BDT threshold per centrality interval"},
            {"name": "slide11 explanatory line", "kind": "text", "role": "audience", "title_axis_align": "left", "bbox": [133, 165, 2400, 215], "font_px": 36, "text": "Corrected 14-feature AuAu candidate uses 15 to 35 GeV photon energy points."},
        ]

    add_box(fig, (0.100, 0.072, 0.800, 0.160), face="#f8fafc", edge="#d5dee9", lw=0.95, radius=0.003)
    table_ax = fig.add_axes([0.116, 0.088, 0.768, 0.127])
    table_ax.set_xlim(0, 1)
    table_ax.set_ylim(0, 1)
    table_ax.axis("off")
    table_ax.text(0.0, 0.87, "Flat-fit threshold", fontsize=11.1, weight="bold", color=INK, va="center")
    columns = [0.25 + i * (0.75 / 7.0) for i in range(7)]
    for x, label in zip(columns, labels, strict=True):
        table_ax.text(x, 0.87, display_centrality(label), fontsize=10.6, color=MUTED, ha="center", va="center")
    rows_y = {"WP90": 0.61, "WP80": 0.38, "WP70": 0.15}
    for wp in WP_ORDER:
        table_ax.text(0.0, rows_y[wp], wp, fontsize=10.8, color=WP_COLORS[wp], weight="bold", va="center")
        for x, label in zip(columns, labels, strict=True):
            value = flat_row(payload, label, wp)["threshold"]
            table_ax.text(x, rows_y[wp], f"{value:.3f}", fontsize=10.7, color=INK, ha="center", va="center")
    table_ax.text(1.0, 0.02, "Embedded simulation  •  corrected shower inputs  •  PPG12 source-role labels", fontsize=9.2, color=MUTED, ha="right", va="bottom")
    fig.savefig(SLIDE11, facecolor="white")
    plt.close(fig)

    return [
        {"name": "slide11 title", "kind": "text", "role": "title", "title_anchor": True, "title_axis_align": "left", "bbox": [128, 55, 2290, 130], "font_px": 68, "text": "Flat fits define one BDT threshold per centrality interval"},
        {"name": "slide11 explanatory band", "kind": "text", "role": "audience", "title_axis_align": "left", "bbox": [128, 160, 2380, 224], "font_px": 31, "text": "Corrected 14-feature AuAu candidate uses 15 to 35 GeV photon energy points."},
        {"name": "slide11 threshold table", "kind": "box", "bbox": [256, 1094, 2304, 1290]},
    ]


def render_slide12(payload: dict) -> list[dict]:
    fig = figure()
    title_and_band(
        fig,
        "Centrality-dependent WP80 candidate threshold",
        r"Seven WP80 thresholds from $15\leq E_T^{\gamma}<35$ GeV are interpolated versus centrality to define one frozen candidate surface.",
        title_fontsize=30.5,
        detail_fontsize=18.8,
        boxed_detail=False,
        detail_y=0.842,
    )
    fits = {row["wp_label"]: row for row in payload["continuous_fits"]}
    x_line = np.linspace(0.0, 80.0, 401)
    # Preserve a clear gutter between the x-axis label and the lower contract
    # band; the previous placement allowed those two text regions to touch.
    ax = fig.add_axes([0.075, 0.240, 0.585, 0.525])
    for wp in WP_ORDER:
        selected = sorted((row for row in payload["flat_rows"] if row["wp_label"] == wp), key=lambda row: row["centrality_center"])
        x = np.asarray([row["centrality_center"] for row in selected], dtype=float)
        y = np.asarray([row["threshold"] for row in selected], dtype=float)
        yerr = np.asarray([row["threshold_stat_err"] for row in selected], dtype=float)
        fit = fits[wp]
        center = fit["intercept"] + fit["slope_per_centrality_percentile"] * x_line
        ax.errorbar(
            x,
            y,
            yerr=yerr,
            fmt=WP_MARKERS[wp],
            color=WP_COLORS[wp],
            markersize=6.4,
            markeredgecolor="white",
            markeredgewidth=0.6,
            capsize=2.1,
            elinewidth=0.95,
            label=f"{wp} constants",
            zorder=3,
        )
        ax.plot(x_line, center, color=WP_COLORS[wp], lw=2.0, label=f"{wp} fit", zorder=2)
        if wp == "WP80":
            ax.fill_between(x_line, center - fit["rms_residual"], center + fit["rms_residual"], color=WP_COLORS[wp], alpha=0.14, lw=0, label="WP80 RMS fit scatter")
    ax.set_xlim(0, 80)
    ax.set_ylim(0.36, 0.78)
    ax.set_xlabel("Centrality percentile $c$")
    ax.set_ylabel("BDT score cut")
    ax.set_xticks((0, 10, 20, 30, 40, 50, 60, 70, 80))
    add_sphenix(ax)
    decorate(ax)
    ax.legend(loc="upper left", ncol=2, fontsize=9.5, frameon=True, framealpha=0.95, edgecolor="#d3dce8", borderpad=0.55, labelspacing=0.45)

    add_box(fig, (0.695, 0.462, 0.255, 0.303), face="#eefbf4", edge="#65b98c", lw=1.25)
    fig.text(0.716, 0.726, "Candidate tight selection", fontsize=16.2, color=INK, weight="bold")
    fig.text(0.716, 0.677, "WP80 centrality surface", fontsize=12.7, color=MUTED)
    # The source/notes retain the full coefficients; the audience-facing card
    # uses a readable rounded display without crowding the right boundary.
    fig.text(0.716, 0.627, r"$T_{80}^{\rm cand}(c)=0.5544+0.001550\,c$", fontsize=15.1, color="#16885d", weight="bold")
    fig.text(0.716, 0.571, r"Training range:  $15\leq E_T^{\gamma}<35$ GeV", fontsize=12.6, color=INK)
    fig.text(0.716, 0.529, r"Tight selection:  score $>T_{80}^{\rm cand}(c)$", fontsize=12.6, color=INK)
    fig.text(0.716, 0.484, "Defined from simulation; data application is separate.", fontsize=10.8, color=MUTED)

    wp80 = fits["WP80"]
    add_box(fig, (0.695, 0.258, 0.255, 0.160), face="#fbfcfe", edge="#d5dee9", lw=1.0)
    fig.text(0.716, 0.378, "WP80 interpolation check", fontsize=16.0, color=INK, weight="bold")
    fig.text(0.716, 0.326, "RMS residual", fontsize=13.6, color=INK)
    fig.text(0.928, 0.326, f"{wp80['rms_residual']:.5f}", fontsize=13.6, color=INK, ha="right")
    fig.text(0.716, 0.280, "Maximum residual", fontsize=13.6, color=INK)
    fig.text(0.928, 0.280, f"{wp80['max_abs_residual']:.5f}", fontsize=13.6, color=INK, ha="right")

    add_box(fig, (0.075, 0.040, 0.875, 0.085), face="#f8fafc", edge="#d5dee9", lw=0.9, radius=0.003)
    fig.text(
        0.512,
        0.094,
        "This centrality-dependent WP80 surface defines the candidate BDT cut.",
        fontsize=16.4,
        color=INK,
        weight="bold",
        ha="center",
    )
    fig.text(
        0.512,
        0.060,
        "The non-tight sideband is being tuned using BDT–isolation interdependence; work in progress.",
        fontsize=15.8,
        color=MUTED,
        ha="center",
    )
    fig.savefig(SLIDE12_CLEANED, facecolor="white")
    plt.close(fig)

    return [
        {"name": "slide12 title", "kind": "text", "role": "title", "title_anchor": True, "title_axis_align": "left", "bbox": [128, 45, 2380, 135], "font_px": 76, "text": "Centrality-dependent WP80 candidate threshold"},
        {"name": "slide12 explanatory band", "kind": "text", "role": "audience", "title_axis_align": "left", "bbox": [128, 185, 2420, 250], "font_px": 42, "text": "Seven 15 to 35 GeV WP80 thresholds are interpolated versus centrality."},
        {"name": "slide12 bottom contract", "kind": "text", "role": "audience", "bbox": [192, 1260, 2432, 1380], "font_px": 36.5, "text": "The centrality-dependent WP80 surface defines the candidate BDT cut."},
    ]


def write_csv(payload: dict) -> None:
    rows: list[dict] = []
    for row in payload["et_cells"]:
        rows.append({"record": "et_point", **row})
    for row in payload["flat_rows"]:
        rows.append({"record": "flat_constant", **row})
    for row in payload["continuous_fits"]:
        rows.append({"record": "centrality_fit", **row})
    fields = sorted({key for row in rows for key in row})
    with CSV_PATH.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_notes(payload: dict) -> None:
    fit = next(row for row in payload["continuous_fits"] if row["wp_label"] == "WP80")
    SLIDE11_NOTES.write_text(
        "\n".join(
            [
                "# Flat fits define one BDT threshold per centrality interval",
                "",
                "These are the fixed-efficiency thresholds from the combined corrected-shower and PPG12-source-role candidate. In each centrality interval, the eight 15–35 GeV points have only modest residual energy dependence, so an inverse-variance weighted constant is sufficient for each efficiency target.",
                "",
                "The green WP80 row is the working-point input. The other two rows are context curves showing the same construction at 90% and 70% signal efficiency.",
            ]
        )
        + "\n"
    )
    SLIDE12_CLEANED_NOTES.write_text(
        "\n".join(
            [
                "# Centrality-dependent WP80 candidate threshold",
                "",
                f"The green points are the seven centrality constants from the previous slide. Their weighted linear interpolation gives T80(c) = {fit['intercept']:.7f} + {fit['slope_per_centrality_percentile']:.8f} c.",
                "",
                f"The RMS scatter is {fit['rms_residual']:.5f} and the largest residual is {fit['max_abs_residual']:.5f}. This band describes interpolation scatter, not an uncertainty on the physics result.",
                "",
                "THE-112 uses this frozen candidate surface when it evaluates non-tight sideband choices. It does not retrain or retune the BDT.",
            ]
        )
        + "\n"
    )


def write_manifest(payload: dict, layout11: list[dict], layout12: list[dict]) -> None:
    wp80 = next(row for row in payload["continuous_fits"] if row["wp_label"] == "WP80")
    manifest = {
        "schema": "THE45_THE111_WORKING_POINT_SLIDES_V1",
        "campaign": "THE-45 JSTG July 22 practice deck",
        "reference_deck": {
            "presentation_id": "1eTv92LXUw9Q6B7yR1R1dbc91BV55fdn8Teab2laf4cQ",
            "slide11_object_id": "g3f56ef1138e_0_407",
            "slide12_object_id": "g3f56ef1138e_0_412",
            "replaced_image_titles": [
                "the79_binned5to40_wp80_8to40_flat_points_slide.png",
                "the79_binned5to40_wp80_8to40_linear_fits_slide.png",
            ],
        },
        "source": {
            "path": str(SOURCE),
            "sha256": sha256(SOURCE),
            "campaign": payload["campaign"],
            "model_label": payload["model_label"],
            "selection": "15 <= reconstructed photon ET < 35 GeV; seven centrality intervals; embedded simulation",
            "label_contract": payload["label_contract"],
            "shower_contract": "corrected calibrated TowerInfo full good 7x7; zero AuAu shower-cell floor",
        },
        "working_point": {
            "wp80_formula": wp80["formula"],
            "wp80_rms_residual": wp80["rms_residual"],
            "wp80_max_abs_residual": wp80["max_abs_residual"],
        },
        "outputs": {
            "slide11_png": str(SLIDE11),
            "slide12_png": str(SLIDE12_CLEANED),
            "slide11_speaker_notes": str(SLIDE11_NOTES),
            "slide12_speaker_notes": str(SLIDE12_CLEANED_NOTES),
            "csv": str(CSV_PATH),
        },
        "boundaries": [
            "THE-111 simulation-validation working-point payload; not a promoted canonical model",
            "THE-112 consumes the frozen candidate and does not rederive or retune it",
            "no data scoring, production rerun, IAN mutation, or Google Slides mutation",
        ],
        "layouts": {"slide11": layout11, "slide12": layout12},
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")
    SLIDE11.with_suffix(".layout_nodes.json").write_text(json.dumps({"schema": "slide_layout_nodes_v1", "title_axis_x": 128, "nodes": layout11}, indent=2) + "\n")
    SLIDE12_CLEANED.with_suffix(".layout_nodes.json").write_text(json.dumps({"schema": "slide_layout_nodes_v1", "title_axis_x": 128, "nodes": layout12}, indent=2) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--slide11-only",
        action="store_true",
        help="render only slide 11; leaves slide 12, CSV, notes, and manifest untouched",
    )
    parser.add_argument(
        "--clean-layout",
        action="store_true",
        help="slide 11 audience layout: no header panel, no threshold table, larger symmetric panels, larger legend",
    )
    parser.add_argument(
        "--variant-suffix",
        default="",
        help="suffix appended to the slide 11 filename, e.g. '_clean'",
    )
    args = parser.parse_args()

    global CLEAN_LAYOUT, SLIDE11
    CLEAN_LAYOUT = args.clean_layout
    if args.variant_suffix:
        SLIDE11 = SLIDE11.with_name(f"{SLIDE11.stem}{args.variant_suffix}{SLIDE11.suffix}")

    if not SOURCE.exists():
        raise SystemExit(f"Missing frozen THE-111 working-point source: {SOURCE}")
    payload = json.loads(SOURCE.read_text(encoding="utf-8"))
    if payload.get("campaign") != "THE-111" or payload.get("status") != "SIMULATION_VALIDATION_ONLY_NO_PROMOTION":
        raise SystemExit("Unexpected working-point payload identity or status")
    if len(payload["et_cells"]) != 168 or len(payload["flat_rows"]) != 21 or len(payload["continuous_fits"]) != 3:
        raise SystemExit("THE-111 payload does not have expected 8×7×3 working-point coverage")
    OUTDIR.mkdir(parents=True, exist_ok=True)
    configure()
    layout11 = render_slide11(payload)
    if args.slide11_only:
        print(SLIDE11)
        return
    layout12 = render_slide12(payload)
    write_csv(payload)
    write_notes(payload)
    write_manifest(payload, layout11, layout12)
    for path in (SLIDE11, SLIDE12_CLEANED, SLIDE11_NOTES, SLIDE12_CLEANED_NOTES, CSV_PATH, MANIFEST):
        print(path)


if __name__ == "__main__":
    main()
