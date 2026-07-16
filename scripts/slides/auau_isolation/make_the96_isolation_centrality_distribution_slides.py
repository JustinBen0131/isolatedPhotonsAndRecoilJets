#!/usr/bin/env python3
"""Build THE-96 AuAu isolation-distribution slides for R=0.3 and R=0.4."""

from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
import sys
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import uproot
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
from scipy.optimize import curve_fit

THIS_FILE = Path(__file__).resolve()
REPO = next((p for p in THIS_FILE.parents if (p / "AGENTS.md").exists()), THIS_FILE.parents[4])
SCRIPTS = REPO / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.append(str(SCRIPTS))

from slides.common.slide_defaults import (  # noqa: E402
    SLIDE_DPI,
    SLIDE_HEIGHT_PX,
    SLIDE_WIDTH_PX,
    slide_figsize,
)

DEFAULT_INPUT = (
    REPO
    / "InputFiles/the96_auau_iso_baseline5pct_20260709_1335"
    / "expanded_5to40_baseline_simembedded"
    / "RecoilJets_embeddedPhoton12plus20_MERGED.root"
)
DEFAULT_OUTPUT = (
    REPO
    / "dataOutput/auau/isolation_baseline5pct_the96_20260709"
    / "baseline_nominal_15to35_centrality_overlays"
)
PT_BINS = [(15, 17), (17, 19), (19, 21), (21, 23), (23, 26), (26, 30), (30, 35)]
CENTRALITY_BINS = [(lo, lo + 5) for lo in range(0, 80, 5)]
REBIN_FACTOR = 2

TIMES_DIR = Path("/System/Library/Fonts/Supplemental")
TIMES_FONTS = [
    TIMES_DIR / "Times New Roman.ttf",
    TIMES_DIR / "Times New Roman Bold.ttf",
    TIMES_DIR / "Times New Roman Italic.ttf",
    TIMES_DIR / "Times New Roman Bold Italic.ttf",
]
INK = "#111827"
MUTED = "#475569"
GRID = "#D7DEE8"
ZERO = "#6B7280"


@dataclass
class GaussianFit:
    cone: str
    cent_min: int
    cent_max: int
    amplitude_counts: float
    location_gev: float
    location_error_gev: float
    sigma_left_gev: float
    sigma_left_error_gev: float
    sigma_right_gev: float
    sigma_right_error_gev: float
    distribution_mean_gev: float
    distribution_mean_error_gev: float
    fit_min_gev: float
    fit_max_gev: float
    chi2: float
    ndf: int
    sumw: float
    effective_entries: float


@dataclass
class AggregatedHistogram:
    cone: str
    cent_min: int
    cent_max: int
    values: np.ndarray
    variances: np.ndarray
    edges: np.ndarray
    sumw: float
    effective_entries: float


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
            "axes.facecolor": "white",
            "axes.edgecolor": INK,
            "axes.linewidth": 1.0,
            "xtick.color": INK,
            "ytick.color": INK,
            "axes.unicode_minus": False,
        }
    )


def bifurcated_gaussian(
    x: np.ndarray,
    amplitude: float,
    location: float,
    sigma_left: float,
    sigma_right: float,
) -> np.ndarray:
    sigma = np.where(x < location, sigma_left, sigma_right)
    return amplitude * np.exp(-0.5 * ((x - location) / sigma) ** 2)


def rebin_histogram(values, variances, edges, factor: int):
    if factor <= 0 or len(values) % factor:
        raise ValueError(f"cannot rebin {len(values)} bins by factor {factor}")
    groups = len(values) // factor
    return (
        values.reshape(groups, factor).sum(axis=1),
        variances.reshape(groups, factor).sum(axis=1),
        edges[::factor],
    )


def _fit_pass(
    x,
    values,
    variances,
    *,
    mode,
    peak,
    location,
    sigma_left,
    sigma_right,
    fit_min,
    fit_max,
):
    mask = (
        (x >= fit_min)
        & (x <= fit_max)
        & (values > 0.0)
        & (variances > 0.0)
        & np.isfinite(values)
        & np.isfinite(variances)
    )
    if int(mask.sum()) < 4:
        raise RuntimeError(f"Gaussian fit has only {int(mask.sum())} usable bins")
    return curve_fit(
        bifurcated_gaussian,
        x[mask],
        values[mask],
        p0=[
            max(float(values[np.argmin(np.abs(x - location))]), 1.0),
            location,
            sigma_left,
            sigma_right,
        ],
        sigma=np.sqrt(variances[mask]),
        absolute_sigma=True,
        bounds=(
            [0.0, mode - 3.0, 0.10, 0.10],
            [10.0 * max(peak, 1.0), mode + 3.0, 20.0, 20.0],
        ),
        maxfev=50000,
    )


def fit_gaussian_iterative(hist: AggregatedHistogram) -> GaussianFit:
    """Fit a true bifurcated Gaussian using the macro's mode/HWHM seed strategy."""
    values, variances, edges = rebin_histogram(
        hist.values, hist.variances, hist.edges, REBIN_FACTOR
    )
    centers = 0.5 * (edges[:-1] + edges[1:])
    maximum_bin = int(np.argmax(values))
    peak = float(values[maximum_bin])
    if not peak > 0.0:
        raise RuntimeError(f"empty peak for {hist.cone} {hist.cent_min}-{hist.cent_max}%")
    mode = float(centers[maximum_bin])
    half_maximum = 0.5 * peak
    hwhm_to_sigma = 1.0 / math.sqrt(2.0 * math.log(2.0))
    hwhm_left = next(
        (mode - float(centers[i]) for i in range(maximum_bin - 1, -1, -1) if values[i] <= half_maximum),
        -1.0,
    )
    hwhm_right = next(
        (float(centers[i]) - mode for i in range(maximum_bin + 1, len(values)) if values[i] <= half_maximum),
        -1.0,
    )
    if hwhm_left > 0.0:
        sigma_left = hwhm_left * hwhm_to_sigma
    else:
        positive = np.clip(values, 0.0, None)
        weighted_mean = float(np.average(centers, weights=positive))
        weighted_rms = math.sqrt(float(np.average((centers - weighted_mean) ** 2, weights=positive)))
        sigma_left = max(0.5, 0.35 * weighted_rms)
    if hwhm_right > 0.0:
        sigma_right = hwhm_right * hwhm_to_sigma
    else:
        sigma_right = sigma_left
    sigma_left = sigma_left if sigma_left > 0.0 and math.isfinite(sigma_left) else 0.5
    sigma_right = sigma_right if sigma_right > 0.0 and math.isfinite(sigma_right) else 0.5

    location = mode
    covariance = np.zeros((4, 4), dtype=float)
    for _ in range(4):
        fitted, covariance = _fit_pass(
            centers,
            values,
            variances,
            mode=mode,
            peak=peak,
            location=location,
            sigma_left=sigma_left,
            sigma_right=sigma_right,
            fit_min=location - 1.75 * sigma_left,
            fit_max=location + 1.75 * sigma_right,
        )
        location = float(fitted[1])
        sigma_left = abs(float(fitted[2]))
        sigma_right = abs(float(fitted[3]))

    fit_min = location - 1.75 * sigma_left
    fit_max = location + 1.75 * sigma_right
    fitted, covariance = _fit_pass(
        centers,
        values,
        variances,
        mode=mode,
        peak=peak,
        location=location,
        sigma_left=sigma_left,
        sigma_right=sigma_right,
        fit_min=fit_min,
        fit_max=fit_max,
    )
    errors = np.sqrt(np.diag(covariance))
    fit_mask = (
        (centers >= fit_min)
        & (centers <= fit_max)
        & (values > 0.0)
        & (variances > 0.0)
    )
    residual = values[fit_mask] - bifurcated_gaussian(centers[fit_mask], *fitted)
    chi2 = float(np.sum((residual**2) / variances[fit_mask]))
    mean_gradient = np.array([1.0, -math.sqrt(2.0 / math.pi), math.sqrt(2.0 / math.pi)])
    mean_covariance = covariance[1:4, 1:4]
    distribution_mean = float(fitted[1] + math.sqrt(2.0 / math.pi) * (fitted[3] - fitted[2]))
    distribution_mean_error = math.sqrt(max(float(mean_gradient @ mean_covariance @ mean_gradient), 0.0))
    return GaussianFit(
        cone=hist.cone,
        cent_min=hist.cent_min,
        cent_max=hist.cent_max,
        amplitude_counts=float(fitted[0]),
        location_gev=float(fitted[1]),
        location_error_gev=float(errors[1]),
        sigma_left_gev=abs(float(fitted[2])),
        sigma_left_error_gev=float(errors[2]),
        sigma_right_gev=abs(float(fitted[3])),
        sigma_right_error_gev=float(errors[3]),
        distribution_mean_gev=distribution_mean,
        distribution_mean_error_gev=distribution_mean_error,
        fit_min_gev=float(fit_min),
        fit_max_gev=float(fit_max),
        chi2=chi2,
        ndf=max(int(fit_mask.sum()) - 4, 0),
        sumw=hist.sumw,
        effective_entries=hist.effective_entries,
    )


def load_aggregated_histograms(input_root: Path, cone: str) -> list[AggregatedHistogram]:
    output = []
    with uproot.open(input_root) as root_file:
        if "SIM" not in root_file:
            raise RuntimeError(f"SIM directory missing from {input_root}")
        directory = root_file["SIM"]
        for cent_min, cent_max in CENTRALITY_BINS:
            total_values = total_variances = common_edges = None
            for pt_min, pt_max in PT_BINS:
                name = (
                    f"h_EisoReco_truthSigMatched_iso{cone}_pT_{pt_min}_{pt_max}"
                    f"_cent_{cent_min}_{cent_max}"
                )
                if name not in directory:
                    raise RuntimeError(f"missing histogram SIM/{name}")
                histogram = directory[name]
                values, edges = histogram.to_numpy(flow=False)
                variances = histogram.variances(flow=False)
                if variances is None:
                    raise RuntimeError(f"missing Sumw2 variances for SIM/{name}")
                values = np.asarray(values, dtype=float)
                variances = np.asarray(variances, dtype=float)
                edges = np.asarray(edges, dtype=float)
                if common_edges is None:
                    common_edges = edges
                    total_values = np.zeros_like(values)
                    total_variances = np.zeros_like(variances)
                elif not np.array_equal(common_edges, edges):
                    raise RuntimeError(f"axis mismatch for SIM/{name}")
                total_values += values
                total_variances += variances
            assert total_values is not None and total_variances is not None and common_edges is not None
            if float(np.min(total_values)) < -1.0e-9:
                raise RuntimeError(f"negative aggregate bin for {cone} {cent_min}-{cent_max}%")
            sumw = float(total_values.sum())
            sumw2 = float(total_variances.sum())
            if not (sumw > 0.0 and sumw2 > 0.0):
                raise RuntimeError(f"empty aggregate for {cone} {cent_min}-{cent_max}%")
            output.append(
                AggregatedHistogram(
                    cone,
                    cent_min,
                    cent_max,
                    total_values,
                    total_variances,
                    common_edges,
                    sumw,
                    (sumw * sumw) / sumw2,
                )
            )
    return output


def font_px(points: float) -> float:
    return points * SLIDE_DPI / 72.0


def text_bbox_px(fig: plt.Figure, artist) -> list[float]:
    fig.canvas.draw()
    bbox = artist.get_window_extent(renderer=fig.canvas.get_renderer())
    scale_x = fig.bbox.width / SLIDE_WIDTH_PX
    scale_y = fig.bbox.height / SLIDE_HEIGHT_PX
    return [
        float(bbox.x0 / scale_x),
        float((fig.bbox.height - bbox.y1) / scale_y),
        float(bbox.x1 / scale_x),
        float((fig.bbox.height - bbox.y0) / scale_y),
    ]


def add_text(fig, nodes, name, x, y, text, *, size, role="audience", weight="normal", color=INK, ha="left", va="center", **flags):
    artist = fig.text(
        x, y, text, fontsize=size, fontweight=weight, color=color, ha=ha, va=va, zorder=10
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
    return artist


def add_artist_node(fig, nodes, name, artist, *, size, text):
    nodes.append(
        {
            "name": name,
            "kind": "text",
            "role": "plot_annotation",
            "text": text,
            "font_px": font_px(size),
            "bbox": text_bbox_px(fig, artist),
            "title_anchor": False,
        }
    )


def distribution_x_limits(histograms):
    lows, highs = [], []
    for hist in histograms:
        values, _, edges = rebin_histogram(hist.values, hist.variances, hist.edges, REBIN_FACTOR)
        cumulative = np.cumsum(np.clip(values, 0.0, None))
        cumulative /= cumulative[-1]
        low_index = int(np.searchsorted(cumulative, 0.001, side="left"))
        high_index = int(np.searchsorted(cumulative, 0.999, side="left"))
        lows.append(float(edges[max(low_index, 0)]))
        highs.append(float(edges[min(high_index + 1, len(edges) - 1)]))
    return math.floor(min(lows) - 1.0), math.ceil(max(highs) + 1.0)


def padded_limits(values, errors):
    low = float(np.min(values - errors))
    high = float(np.max(values + errors))
    span = max(high - low, 0.2)
    low -= 0.16 * span
    high += 0.16 * span
    if low < 0.0 < high:
        pad = max(abs(low), abs(high))
        low, high = -1.05 * pad, 1.05 * pad
    return low, high


def render_slide(cone, histograms, fits, output_dir, *, variant_key, title_subject, subtitle, algorithm):
    radius = "0.3" if cone == "R30" else "0.4"
    radius_tag = "r03" if cone == "R30" else "r04"
    stem = output_dir / f"auau_{radius_tag}_{variant_key}_isolation_centrality_overlay"
    png = Path(f"{stem}_slide.png")
    layout = Path(f"{stem}_layout_nodes.json")
    audit = Path(f"{stem}_slide_audit.json")
    script = Path(f"{stem}_speaker_script.md")
    colors = [plt.colormaps["viridis"](v) for v in np.linspace(0.03, 0.95, len(histograms))]

    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI, facecolor="white")
    nodes = []
    add_text(
        fig, nodes, "title", 0.045, 0.950,
        f"R = {radius} {title_subject} by 5% centrality, 15-35 GeV",
        size=31.0, role="title", weight="bold", va="top",
    )
    subtitle_artist = add_text(
        fig, nodes, "subtitle", 0.050, 0.842,
        subtitle,
        size=15.0, color=MUTED,
    )
    main_ax = fig.add_axes([0.065, 0.130, 0.590, 0.545])
    mean_ax = fig.add_axes([0.715, 0.130, 0.245, 0.545])
    legend_handles = []
    counts_maximum = 0.0
    x_min, x_max = distribution_x_limits(histograms)
    for index, (hist, fit, color) in enumerate(zip(histograms, fits, colors)):
        values, variances, edges = rebin_histogram(hist.values, hist.variances, hist.edges, REBIN_FACTOR)
        centers = 0.5 * (edges[:-1] + edges[1:])
        visible = (centers >= x_min) & (centers <= x_max)
        counts_maximum = max(counts_maximum, float(values[visible].max()))
        main_ax.errorbar(
            centers[visible], values[visible], yerr=np.sqrt(variances[visible]),
            color=color, marker="o", markersize=1.9, linewidth=0.9,
            elinewidth=0.55, capsize=0.0, alpha=0.92, zorder=2 + index,
        )
        fit_x = np.linspace(fit.fit_min_gev, fit.fit_max_gev, 220)
        fit_y = bifurcated_gaussian(
            fit_x,
            fit.amplitude_counts,
            fit.location_gev,
            fit.sigma_left_gev,
            fit.sigma_right_gev,
        )
        main_ax.plot(fit_x, fit_y, color=color, linestyle="--", linewidth=1.3, alpha=0.95, zorder=20 + index)
        legend_handles.append(
            Line2D([0], [0], color=color, linewidth=2.2, label=f"{hist.cent_min}-{hist.cent_max}%")
        )
    main_ax.set_xlim(x_min, x_max)
    main_ax.set_ylim(0.0, 1.16 * counts_maximum)
    main_ax.set_xlabel(r"Reconstructed $E_T^{\mathrm{iso}}$ [GeV]", fontsize=16.0, labelpad=7)
    main_title = main_ax.set_title("Rebinned isolation distributions", loc="left", fontsize=17.0, fontweight="bold", pad=9)
    main_ax.tick_params(labelsize=13.0, length=4)
    main_ax.yaxis.set_major_locator(MaxNLocator(nbins=6, min_n_ticks=4))
    main_ax.grid(color=GRID, linewidth=0.8, alpha=0.75)
    main_ax.set_axisbelow(True)
    for spine in main_ax.spines.values():
        spine.set_color(INK)
        spine.set_linewidth(1.0)

    internal_word = main_ax.text(
        0.975, 0.955, "Internal", transform=main_ax.transAxes,
        ha="right", va="top", fontsize=13.0, color=INK,
    )
    fig.canvas.draw()
    internal_bbox = internal_word.get_window_extent(renderer=fig.canvas.get_renderer())
    sphenix_x = main_ax.transAxes.inverted().transform((internal_bbox.x0 - 5.0, internal_bbox.y0))[0]
    sphenix_word = main_ax.text(
        sphenix_x, 0.955, "sPHENIX", transform=main_ax.transAxes,
        ha="right", va="top", fontstyle="italic", fontweight="bold",
        fontsize=13.0, color=INK,
    )
    sample_word = main_ax.text(
        0.975, 0.905, "Embedded Photon+Jet 12+20", transform=main_ax.transAxes,
        ha="right", va="top", fontsize=12.2, color=INK,
    )

    centrality = np.array([0.5 * (fit.cent_min + fit.cent_max) for fit in fits])
    means = np.array([fit.distribution_mean_gev for fit in fits])
    mean_errors = np.array([fit.distribution_mean_error_gev for fit in fits])
    mean_ax.axhline(0.0, color=ZERO, linestyle="--", linewidth=1.3, zorder=1)
    mean_ax.plot(centrality, means, color=INK, linewidth=1.35, alpha=0.75, zorder=2)
    for index, (x_value, mean, error, color) in enumerate(zip(centrality, means, mean_errors, colors)):
        mean_ax.errorbar(
            [x_value], [mean], yerr=[error], color=color, marker="o",
            markeredgecolor=INK, markeredgewidth=0.45, markersize=6.4,
            elinewidth=1.2, capsize=2.0, zorder=3 + index,
        )
    y_min, y_max = padded_limits(means, mean_errors)
    mean_ax.set_xlim(0.0, 80.0)
    mean_ax.set_ylim(y_min, y_max)
    mean_ax.set_xticks([0, 20, 40, 60, 80])
    mean_ax.yaxis.set_major_locator(MaxNLocator(nbins=6, min_n_ticks=4))
    mean_ax.set_xlabel("Centrality [%]", fontsize=15.0, labelpad=7)
    mean_ax.set_ylabel("Asymmetric-Gaussian mean [GeV]", fontsize=15.0, labelpad=8)
    mean_title = mean_ax.set_title("Asymmetric Gaussian mean", loc="left", fontsize=17.0, fontweight="bold", pad=9)
    mean_ax.tick_params(labelsize=12.5, length=4)
    mean_ax.grid(color=GRID, linewidth=0.8, alpha=0.75)
    mean_ax.set_axisbelow(True)
    for spine in mean_ax.spines.values():
        spine.set_color(INK)
        spine.set_linewidth(1.0)
    mean_ax.text(
        0.96, 0.945,
        "Asymmetric Gaussian fit\n"
        r"$f(x)=A\exp\!\left[-\frac{(x-\mu)^2}{2\sigma_{L/R}^{2}}\right]$"
        "\n"
        r"$\sigma_L\ (x<\mu),\quad \sigma_R\ (x\geq\mu)$",
        transform=mean_ax.transAxes, ha="right", va="top",
        fontsize=10.8, color=MUTED,
    )

    legend_order = [index for pair in zip(range(8), range(8, 16)) for index in pair]
    legend = fig.legend(
        handles=[legend_handles[index] for index in legend_order],
        loc="upper center", bbox_to_anchor=(0.510, 0.785),
        ncol=8, frameon=False, fontsize=11.8, handlelength=1.35,
        handletextpad=0.38, columnspacing=1.05,
    )
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    subtitle_bbox = subtitle_artist.get_window_extent(renderer=renderer)
    main_title_bbox = main_title.get_window_extent(renderer=renderer)
    legend_bbox = legend.get_window_extent(renderer=renderer)
    # Center the legend in the actual rendered gap between the subtitle and
    # the LHS panel title; this keeps both slide variants geometrically matched.
    gap_center_display = 0.5 * (subtitle_bbox.y0 + main_title_bbox.y1)
    legend_top_display = gap_center_display + 0.5 * legend_bbox.height
    legend.set_bbox_to_anchor(
        (0.510, legend_top_display / fig.bbox.height), transform=fig.transFigure
    )
    fig.canvas.draw()
    add_artist_node(fig, nodes, "centrality legend", legend.get_texts()[0], size=11.8, text="0-5% through 75-80%")
    add_artist_node(fig, nodes, "main axis label", main_ax.xaxis.label, size=16.0, text="Reconstructed isolation")
    add_artist_node(fig, nodes, "mean panel title", mean_title, size=17.0, text="Asymmetric Gaussian mean")
    add_artist_node(fig, nodes, "sPHENIX label", sphenix_word, size=13.0, text="sPHENIX")
    add_artist_node(fig, nodes, "Internal label", internal_word, size=13.0, text="Internal")
    add_artist_node(fig, nodes, "sample label", sample_word, size=12.2, text="Embedded Photon+Jet 12+20")

    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(png, dpi=SLIDE_DPI)
    layout.write_text(
        json.dumps(
            {
                "schema": "slide_layout_nodes_v1",
                "title_axis_x": 0.045 * SLIDE_WIDTH_PX,
                "minimum_audience_font_px": font_px(14.0),
                "minimum_title_font_px": font_px(30.0),
                "minimum_plot_annotation_font_px": font_px(11.5),
                "nodes": nodes,
            },
            indent=2,
        ) + "\n",
        encoding="utf-8",
    )
    plt.close(fig)
    script.write_text(
        "\n".join(
            [
                f"# Speaker Script: R={radius} {algorithm} isolation distributions by centrality",
                "",
                f"This slide overlays the truth-matched embedded-photon isolation distributions for R equals {radius} in sixteen five-percent centrality bins.",
                "The photon transverse-momentum range is fifteen to thirty-five GeV, built by summing the seven native pT histograms in each centrality interval.",
                "The distributions are shown as weighted counts after rebinning from 0.1 to 0.2 GeV, with statistical error bars.",
                "Each core is fit with a true bifurcated Gaussian that has independent left and right widths, seeded from the mode and half-maximum crossings.",
                "The right panel shows the mathematical mean of that bifurcated Gaussian, including the shift from unequal left and right widths.",
                "The centrality dependence and sign of the fitted means are determined by the reconstructed distributions and are not imposed by the fit.",
                "",
            ]
        ),
        encoding="utf-8",
    )
    return {
        "cone": cone,
        "radius": radius,
        "png": str(png),
        "layout_nodes": str(layout),
        "audit": str(audit),
        "speaker_script": str(script),
        "x_range_gev": [x_min, x_max],
        "y_range_counts": [0.0, 1.16 * counts_maximum],
        "mean_y_range_gev": [y_min, y_max],
    }


def write_fit_table(output_dir, cone, fits, *, variant_key):
    radius_tag = "r03" if cone == "R30" else "r04"
    path = output_dir / f"auau_{radius_tag}_{variant_key}_isolation_gaussian_fits_15to35.csv"
    rows = [asdict(fit) for fit in fits]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def run_audits(slides):
    audit_script = SCRIPTS / "slides/common/post_render_slide_audit.py"
    for slide in slides:
        subprocess.run(
            [
                sys.executable, str(audit_script), "--png", slide["png"],
                "--layout-nodes", slide["layout_nodes"], "--output", slide["audit"],
                "--require-pass",
            ],
            check=True,
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-root", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--variant-key", default="baseline")
    parser.add_argument("--title-subject", default="isolation counts")
    parser.add_argument(
        "--subtitle",
        default="Embedded Photon+Jet 12+20 | truth-matched photons | baseline tower-cone | 0.2 GeV rebinned counts",
    )
    parser.add_argument("--campaign-tag", default="the96_auau_iso_baseline5pct_20260709_1335")
    parser.add_argument("--algorithm", default="baseline tower-cone without topocluster isolation")
    parser.add_argument(
        "--manifest-schema",
        default="THE96_AUAU_BASELINE_ISOLATION_CENTRALITY_OVERLAY_SLIDES_V2",
    )
    parser.add_argument(
        "--production-note",
        default=(
            "The dedicated nominal baseline signal lane was intentionally closed at 7/8 after Justin accepted the slide result. "
            "The completed expanded baseline Photon12+20 ROOT is used because its generated YAML differs from the nominal YAML "
            "only in photon-pT bin edges; this plot sums only the native 15-35 GeV histogram bins."
        ),
    )
    args = parser.parse_args()
    input_root = args.input_root.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    variant_key = args.variant_key.strip()
    if not variant_key or any(character not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-" for character in variant_key):
        raise ValueError(f"invalid --variant-key: {args.variant_key!r}")
    if not input_root.is_file():
        raise FileNotFoundError(input_root)
    setup_style()
    output_dir.mkdir(parents=True, exist_ok=True)
    slides, fit_tables, fit_summary = [], {}, {}
    for cone in ("R30", "R40"):
        histograms = load_aggregated_histograms(input_root, cone)
        fits = [fit_gaussian_iterative(histogram) for histogram in histograms]
        slides.append(
            render_slide(
                cone,
                histograms,
                fits,
                output_dir,
                variant_key=variant_key,
                title_subject=args.title_subject,
                subtitle=args.subtitle,
                algorithm=args.algorithm,
            )
        )
        fit_tables[cone] = str(write_fit_table(output_dir, cone, fits, variant_key=variant_key))
        fit_summary[cone] = [asdict(fit) for fit in fits]
    run_audits(slides)
    manifest = output_dir / f"auau_{variant_key}_isolation_centrality_overlay_slides_manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "schema": args.manifest_schema,
                "generated_at": datetime.now(timezone.utc).isoformat(),
                "campaign_tag": args.campaign_tag,
                "input_root": str(input_root),
                "root_directory": "SIM",
                "sample": "Embedded Photon+Jet 12+20 truth-matched signal",
                "algorithm": args.algorithm,
                "pt_range_gev": [15, 35],
                "pt_bins_summed": PT_BINS,
                "centrality_bins": CENTRALITY_BINS,
                "histogram_family": "h_EisoReco_truthSigMatched_isoR{30,40}_pT_*_cent_*",
                "normalization": "none; weighted counts are shown after rebinning from 0.1 to 0.2 GeV",
                "fit_contract": {
                    "source_reference": "mode/HWHM seeding adapted from macros/AnalyzeRecoilJets_RunIsoQA_UEComparisons_AuAu.cpp::FitGaussianIterative",
                    "rebin_factor": REBIN_FACTOR,
                    "bin_width_after_rebin_gev": 0.2,
                    "method": "true bifurcated Gaussian with independent left/right widths; four iterative fits over [location-1.75 sigma_left, location+1.75 sigma_right]",
                    "distribution_mean_definition": "location + sqrt(2/pi) * (sigma_right - sigma_left)",
                },
                "production_artifact_exception": args.production_note,
                "fit_tables": fit_tables,
                "fit_results": fit_summary,
                "slides": slides,
                "mutation_boundary": "local PNG candidates only; no Google Slides mutation",
            },
            indent=2,
        ) + "\n",
        encoding="utf-8",
    )
    print(f"Wrote {manifest}")
    for slide in slides:
        print(f"Wrote {slide['png']}")
        print(f"Wrote {slide['audit']}")


if __name__ == "__main__":
    main()
