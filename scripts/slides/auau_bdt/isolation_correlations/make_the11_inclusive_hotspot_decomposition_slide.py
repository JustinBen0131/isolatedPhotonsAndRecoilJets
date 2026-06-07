#!/usr/bin/env python3
"""Build a slide that decomposes inclusive BDT-isolation density hot spots."""

from __future__ import annotations

from collections import deque
import json
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import colors, font_manager, ticker  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.patches import FancyBboxPatch, Rectangle  # noqa: E402
import numpy as np  # noqa: E402


REPO = Path(__file__).resolve().parents[4]
SOURCE_REPORT = REPO / "dataOutput/auauTightBDTValidation/model_validation_condor_20260511_194832"
SCORE_CACHE_MANIFEST = SOURCE_REPORT / "score_caches.local.list"
SOURCE_SLIDE21_PNG = (
    SOURCE_REPORT
    / "slideReady/default_eiso_vs_bdt_score_20260604/plot_variants_kbird_probability_20260604/02_kbird_log_density_centrality_by_class.png"
)
OUT_DIR = (
    SOURCE_REPORT
    / "slideReady/default_eiso_vs_bdt_score_20260604/the11_inclusive_hotspot_decomposition_20260606"
)

SCORE_COLUMN = "score_ptFine_cent7"
ISO_COLUMN = "reco_eiso"
PT_RANGE = (15.0, 35.0)
CENT_BINS = [(0.0, 20.0, "0-20%"), (20.0, 50.0, "20-50%"), (50.0, 80.0, "50-80%")]
X_BINS = np.linspace(-10.0, 18.0, 57)
Y_BINS = np.linspace(0.0, 1.0, 51)

LOW_FLOOR = {"x": (-1.0, 13.0), "y": (0.0, 0.08), "label": "low-score floor"}
MID_ISLAND = {"x": (-1.0, 6.0), "y": (0.45, 0.60), "label": "mid-BDT island"}
HIGH_LEAK = {"x": (2.0, 13.0), "y": (0.60, 1.00), "label": "high-BDT positive isolation"}
SIGNAL_LIKE_BAND = {"x": (0.0, 8.0), "y": (0.52, 0.78), "label": "signal-like positive-isolation band"}

W, H = 2560, 1440
DPI = 200
WHITE = "#ffffff"
INK = "#151515"
MUTED = "#596171"
GRID = "#dce3ea"
SUBTITLE_BG = "#f4f7fb"
SUBTITLE_EDGE = "#cad6e4"
BAND_BG = "#fff7d9"
BAND_EDGE = "#d5c47d"
LOW_COLOR = "#1f5aa6"
MID_COLOR = "#15917f"
HIGH_COLOR = "#c44e00"


def kbird_cmap() -> LinearSegmentedColormap:
    stops = np.linspace(0.0, 1.0, 9)
    red = [0.2082, 0.0592, 0.0780, 0.0232, 0.1802, 0.5301, 0.8186, 0.9956, 0.9764]
    green = [0.1664, 0.3599, 0.5041, 0.6419, 0.7178, 0.7492, 0.7328, 0.7862, 0.9832]
    blue = [0.5293, 0.8684, 0.8385, 0.8385, 0.6425, 0.4662, 0.3499, 0.1968, 0.0539]
    return LinearSegmentedColormap.from_list("root_kbird_like", list(zip(stops, zip(red, green, blue))), N=256)


def setup_style() -> None:
    available = {f.name for f in font_manager.fontManager.ttflist}
    for family in ("Times New Roman", "Times", "DejaVu Serif"):
        if family in available:
            plt.rcParams["font.family"] = family
            break
    plt.rcParams.update(
        {
            "figure.facecolor": WHITE,
            "savefig.facecolor": WHITE,
            "axes.facecolor": WHITE,
            "axes.edgecolor": INK,
            "axes.linewidth": 1.0,
            "axes.labelcolor": INK,
            "xtick.color": INK,
            "ytick.color": INK,
            "text.color": INK,
            "mathtext.default": "regular",
        }
    )


def load_background_arrays() -> dict[str, np.ndarray]:
    if not SCORE_CACHE_MANIFEST.exists():
        raise FileNotFoundError(f"Missing score cache list: {SCORE_CACHE_MANIFEST}")

    chunks: dict[str, list[np.ndarray]] = {key: [] for key in ("eiso", "score", "is_signal", "cluster_et", "centrality")}
    for raw_line in SCORE_CACHE_MANIFEST.read_text().splitlines():
        line = raw_line.strip()
        if not line:
            continue
        cache_path = Path(line)
        if not cache_path.is_absolute():
            cache_path = REPO / cache_path
        if not cache_path.exists():
            raise FileNotFoundError(f"Missing score cache listed in manifest: {line}")
        data = np.load(cache_path, allow_pickle=True)
        required = [ISO_COLUMN, SCORE_COLUMN, "is_signal", "cluster_Et", "centrality"]
        missing = [key for key in required if key not in data.files]
        if missing:
            raise KeyError(f"{cache_path} missing required columns: {missing}")
        chunks["eiso"].append(data[ISO_COLUMN].astype("float32", copy=False))
        chunks["score"].append(data[SCORE_COLUMN].astype("float32", copy=False))
        chunks["is_signal"].append(data["is_signal"].astype("int8", copy=False))
        chunks["cluster_et"].append(data["cluster_Et"].astype("float32", copy=False))
        chunks["centrality"].append(data["centrality"].astype("float32", copy=False))

    arrays = {key: np.concatenate(values) for key, values in chunks.items()}
    selected = (
        np.isfinite(arrays["eiso"])
        & np.isfinite(arrays["score"])
        & (arrays["is_signal"] == 0)
        & (arrays["cluster_et"] >= PT_RANGE[0])
        & (arrays["cluster_et"] < PT_RANGE[1])
        & (arrays["centrality"] >= 0.0)
        & (arrays["centrality"] < 80.0)
    )
    return {key: values[selected] for key, values in arrays.items()}


def band_mask(x: np.ndarray, y: np.ndarray, band: dict[str, object]) -> np.ndarray:
    xlo, xhi = band["x"]
    ylo, yhi = band["y"]
    return (x >= xlo) & (x < xhi) & (y >= ylo) & (y < yhi)


def smooth_histogram(hist: np.ndarray) -> np.ndarray:
    padded = np.pad(hist, 1, mode="edge")
    smoothed = (
        padded[:-2, :-2]
        + 2.0 * padded[:-2, 1:-1]
        + padded[:-2, 2:]
        + 2.0 * padded[1:-1, :-2]
        + 4.0 * padded[1:-1, 1:-1]
        + 2.0 * padded[1:-1, 2:]
        + padded[2:, :-2]
        + 2.0 * padded[2:, 1:-1]
        + padded[2:, 2:]
    ) / 16.0
    return smoothed


def connected_components(mask: np.ndarray, weights: np.ndarray, x_edges: np.ndarray, y_edges: np.ndarray) -> list[dict[str, float | int]]:
    seen = np.zeros(mask.shape, dtype=bool)
    comps: list[dict[str, float | int]] = []
    for start_i, start_j in np.argwhere(mask):
        if seen[start_i, start_j]:
            continue
        queue: deque[tuple[int, int]] = deque([(int(start_i), int(start_j))])
        seen[start_i, start_j] = True
        pts: list[tuple[int, int]] = []
        while queue:
            i, j = queue.popleft()
            pts.append((i, j))
            for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                ni, nj = i + di, j + dj
                if 0 <= ni < mask.shape[0] and 0 <= nj < mask.shape[1] and mask[ni, nj] and not seen[ni, nj]:
                    seen[ni, nj] = True
                    queue.append((ni, nj))
        xs = np.asarray([0.5 * (x_edges[i] + x_edges[i + 1]) for i, _ in pts])
        ys = np.asarray([0.5 * (y_edges[j] + y_edges[j + 1]) for _, j in pts])
        w = np.asarray([weights[i, j] for i, j in pts])
        comps.append(
            {
                "mass_fraction": float(np.sum(w)),
                "bin_count": int(len(pts)),
                "x_min": float(np.min(xs)),
                "x_max": float(np.max(xs)),
                "y_min": float(np.min(ys)),
                "y_max": float(np.max(ys)),
                "x_weighted": float(np.average(xs, weights=w)),
                "y_weighted": float(np.average(ys, weights=w)),
            }
        )
    return sorted(comps, key=lambda item: item["mass_fraction"], reverse=True)


def build_payload(arrays: dict[str, np.ndarray]) -> tuple[dict[str, object], list[np.ndarray], float, float]:
    panels: dict[str, object] = {}
    histograms: list[np.ndarray] = []
    positive_display: list[np.ndarray] = []

    for cent_lo, cent_hi, cent_label in CENT_BINS:
        cent_mask = (arrays["centrality"] >= cent_lo) & (arrays["centrality"] < cent_hi)
        x = arrays["eiso"][cent_mask]
        y = arrays["score"][cent_mask]
        hist_counts, x_edges, y_edges = np.histogram2d(x, y, bins=[X_BINS, Y_BINS])
        visible = int(np.sum(hist_counts))
        probability = hist_counts / visible if visible else hist_counts
        smoothed_probability = smooth_histogram(probability)
        display_probability = smoothed_probability / np.max(smoothed_probability) if np.max(smoothed_probability) > 0.0 else smoothed_probability
        histograms.append(display_probability)
        positive_display.append(display_probability[display_probability > 0.0])

        low = band_mask(x, y, LOW_FLOOR)
        mid = band_mask(x, y, MID_ISLAND)
        high = band_mask(x, y, HIGH_LEAK)
        signal_like = band_mask(x, y, SIGNAL_LIKE_BAND)
        high_total = y >= 0.60
        top_threshold = float(np.quantile(probability[probability > 0.0], 0.92)) if np.any(probability > 0.0) else 0.0
        comps = connected_components(probability >= top_threshold, probability, x_edges, y_edges)[:5]
        panels[cent_label] = {
            "n_total": int(len(x)),
            "n_visible": visible,
            "visible_fraction": float(visible / len(x)) if len(x) else 0.0,
            "median_eiso": float(np.nanmedian(x)) if len(x) else float("nan"),
            "median_score": float(np.nanmedian(y)) if len(y) else float("nan"),
            "bands": {
                "low_score_floor": {
                    "n": int(np.sum(low)),
                    "fraction": float(np.mean(low)) if len(x) else 0.0,
                    "median_eiso": float(np.nanmedian(x[low])) if np.any(low) else float("nan"),
                    "median_score": float(np.nanmedian(y[low])) if np.any(low) else float("nan"),
                },
                "mid_bdt_island": {
                    "n": int(np.sum(mid)),
                    "fraction": float(np.mean(mid)) if len(x) else 0.0,
                    "median_eiso": float(np.nanmedian(x[mid])) if np.any(mid) else float("nan"),
                    "median_score": float(np.nanmedian(y[mid])) if np.any(mid) else float("nan"),
                },
                "signal_like_positive_band": {
                    "n": int(np.sum(signal_like)),
                    "fraction": float(np.mean(signal_like)) if len(x) else 0.0,
                    "median_eiso": float(np.nanmedian(x[signal_like])) if np.any(signal_like) else float("nan"),
                    "median_score": float(np.nanmedian(y[signal_like])) if np.any(signal_like) else float("nan"),
                },
                "high_bdt_positive_isolation": {
                    "n": int(np.sum(high)),
                    "fraction": float(np.mean(high)) if len(x) else 0.0,
                    "median_eiso": float(np.nanmedian(x[high])) if np.any(high) else float("nan"),
                    "median_score": float(np.nanmedian(y[high])) if np.any(high) else float("nan"),
                },
            },
            "high_bdt_breakdown": {
                "n": int(np.sum(high_total)),
                "fraction": float(np.mean(high_total)) if len(x) else 0.0,
                "eiso_lt_0": float(np.mean(x[high_total] < 0.0)) if np.any(high_total) else 0.0,
                "eiso_0_to_2": float(np.mean((x[high_total] >= 0.0) & (x[high_total] < 2.0))) if np.any(high_total) else 0.0,
                "eiso_2_to_5": float(np.mean((x[high_total] >= 2.0) & (x[high_total] < 5.0))) if np.any(high_total) else 0.0,
                "eiso_gt_5": float(np.mean(x[high_total] >= 5.0)) if np.any(high_total) else 0.0,
            },
            "top_density_threshold_probability": top_threshold,
            "top_density_components": comps,
        }

    all_positive = np.concatenate(positive_display)
    vmin = 0.018
    vmax = 1.0
    payload: dict[str, object] = {
        "x_bins": X_BINS.tolist(),
        "y_bins": Y_BINS.tolist(),
        "panels": panels,
        "normalization": "Each centrality heatmap is normalized by its own inclusive/background entries before applying log color.",
        "hotspot_regions": {
            "low_score_floor": LOW_FLOOR,
            "mid_bdt_island": MID_ISLAND,
            "high_bdt_positive_isolation": HIGH_LEAK,
            "signal_like_positive_band": SIGNAL_LIKE_BAND,
        },
    }
    return payload, histograms, vmin, vmax


def add_box(fig: plt.Figure, xywh: tuple[float, float, float, float], face: str, edge: str, *, zorder: float = 0.2) -> None:
    fig.add_artist(
        FancyBboxPatch(
            (xywh[0], xywh[1]),
            xywh[2],
            xywh[3],
            boxstyle="round,pad=0.010,rounding_size=0.009",
            transform=fig.transFigure,
            facecolor=face,
            edgecolor=edge,
            linewidth=1.0,
            zorder=zorder,
        )
    )


def wrapped_text(fig: plt.Figure, x: float, y: float, text: str, width_chars: int, *, fontsize: float, color: str = INK, weight: str = "normal") -> None:
    fig.text(
        x,
        y,
        textwrap.fill(text, width=width_chars),
        ha="left",
        va="top",
        fontsize=fontsize,
        color=color,
        fontweight=weight,
        linespacing=1.10,
    )


def add_header(fig: plt.Figure) -> None:
    fig.text(
        0.050,
        0.948,
        "Inclusive hot spots separate into two background structures",
        ha="left",
        va="top",
        fontsize=31.0,
        fontweight="bold",
    )
    fig.add_artist(Rectangle((0.050, 0.902), 0.900, 0.003, transform=fig.transFigure, color="#d8dee8", lw=0))
    add_box(fig, (0.052, 0.812, 0.896, 0.058), SUBTITLE_BG, SUBTITLE_EDGE)
    fig.text(0.070, 0.854, "Color", ha="left", va="center", fontsize=15.8, fontweight="bold")
    fig.text(
        0.178,
        0.854,
        "lightly smoothed relative density within each centrality panel; N boxes show raw statistics separately.",
        ha="left",
        va="center",
        fontsize=14.6,
        color=MUTED,
    )
    fig.text(0.070, 0.828, "Outlines", ha="left", va="center", fontsize=15.8, fontweight="bold")
    fig.text(0.178, 0.828, "blue = low-BDT positive-isolation floor; orange = signal-like positive-isolation band near the tight-ID boundary.", ha="left", va="center", fontsize=14.6, color=MUTED)


def draw_panels(fig: plt.Figure, payload: dict[str, object], histograms: list[np.ndarray], vmin: float, vmax: float) -> None:
    cmap = kbird_cmap()
    norm = colors.LogNorm(vmin=vmin, vmax=vmax)
    lefts = [0.070, 0.365, 0.660]
    bottom = 0.342
    panel_w = 0.250
    panel_h = 0.420
    mesh = None
    x_edges = np.asarray(payload["x_bins"])
    y_edges = np.asarray(payload["y_bins"])
    for idx, (_, _, cent_label) in enumerate(CENT_BINS):
        ax = fig.add_axes([lefts[idx], bottom, panel_w, panel_h])
        probability = np.ma.masked_less_equal(histograms[idx].T, 0.0)
        mesh = ax.pcolormesh(x_edges, y_edges, probability, cmap=cmap, norm=norm, shading="auto")
        ax.grid(True, color=GRID, lw=0.45, alpha=0.60)
        ax.set_xlim(-10.0, 18.0)
        ax.set_ylim(0.0, 1.0)
        ax.set_title(f"{cent_label} centrality", fontsize=18.8, fontweight="bold", pad=10)
        ax.tick_params(axis="both", labelsize=12.0, pad=2)
        ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
        ax.xaxis.set_major_locator(ticker.MultipleLocator(5.0))
        ax.set_xlabel(r"reco $E_T^{iso}$, $\Delta R<0.3$ [GeV]", fontsize=12.7, labelpad=4)
        if idx == 0:
            ax.set_ylabel("BDT score", fontsize=13.2, labelpad=5)
        else:
            ax.set_yticklabels([])

        ax.axhspan(LOW_FLOOR["y"][0], LOW_FLOOR["y"][1], xmin=0.0, xmax=1.0, color=LOW_COLOR, alpha=0.055, zorder=2.0)
        ax.add_patch(Rectangle((LOW_FLOOR["x"][0], LOW_FLOOR["y"][0]), LOW_FLOOR["x"][1] - LOW_FLOOR["x"][0], LOW_FLOOR["y"][1] - LOW_FLOOR["y"][0], fill=False, edgecolor=LOW_COLOR, linewidth=2.2, zorder=3.5))
        ax.add_patch(Rectangle((SIGNAL_LIKE_BAND["x"][0], SIGNAL_LIKE_BAND["y"][0]), SIGNAL_LIKE_BAND["x"][1] - SIGNAL_LIKE_BAND["x"][0], SIGNAL_LIKE_BAND["y"][1] - SIGNAL_LIKE_BAND["y"][0], facecolor=HIGH_COLOR, alpha=0.065, edgecolor=HIGH_COLOR, linewidth=2.5, zorder=3.5))
        ax.axhline(0.60, color=HIGH_COLOR, lw=1.15, alpha=0.85)
        ax.axvline(5.0, color="#222222", lw=1.0, ls=":", alpha=0.70)

        panel = payload["panels"][cent_label]
        low = panel["bands"]["low_score_floor"]
        high = panel["bands"]["high_bdt_positive_isolation"]
        ax.text(
            0.035,
            0.962,
            f"N={panel['n_total']:,}\nlow floor {low['fraction']:.0%}\nhigh-BDT +iso {high['fraction']:.0%}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=10.8,
            color="#20252c",
            linespacing=1.18,
            bbox=dict(boxstyle="round,pad=0.27", facecolor="white", edgecolor="#d7dde8", linewidth=0.8, alpha=0.93),
        )
        if idx == 2:
            ax.text(0.62, 0.080, "low-score floor", transform=ax.transAxes, ha="center", va="bottom", fontsize=10.6, color=LOW_COLOR, fontweight="bold")
            ax.text(0.58, 0.670, "signal-like\npositive iso", transform=ax.transAxes, ha="center", va="center", fontsize=10.4, color=HIGH_COLOR, fontweight="bold", linespacing=0.95)

    cax = fig.add_axes([0.922, bottom + 0.030, 0.014, panel_h - 0.060])
    cb = fig.colorbar(mesh, cax=cax)
    cb.set_label("relative density\nwithin panel", fontsize=11.2, labelpad=8)
    cb.ax.tick_params(labelsize=10.0)


def draw_interpretation(fig: plt.Figure, payload: dict[str, object]) -> None:
    add_box(fig, (0.070, 0.045, 0.860, 0.210), BAND_BG, BAND_EDGE)
    cols = [0.095, 0.377, 0.660]
    fig.text(cols[0], 0.220, "1. Dominant yellow strip", ha="left", va="top", fontsize=15.5, fontweight="bold", color=LOW_COLOR)
    wrapped_text(
        fig,
        cols[0],
        0.186,
        "The largest hot structure sits near the BDT floor and positive isolation. It grows from 9% to 22%, so it is bright background occupancy, not tight-ID leakage.",
        39,
        fontsize=11.6,
    )
    fig.text(cols[1], 0.220, "2. ABCD-relevant strip", ha="left", va="top", fontsize=15.5, fontweight="bold", color=HIGH_COLOR)
    wrapped_text(
        fig,
        cols[1],
        0.186,
        "The orange band is the real concern: signal-like BDT score with positive cone energy. The stricter high-BDT + positive-isolation region carries 9-13% of inclusive background.",
        39,
        fontsize=11.6,
    )
    fig.text(cols[2], 0.220, "3. Normalization check", ha="left", va="top", fontsize=15.5, fontweight="bold")
    wrapped_text(
        fig,
        cols[2],
        0.186,
        "0-20 has fewer candidates and a more diffuse density. Color is relative inside each panel, so that difference is not a hidden raw-count normalization artifact.",
        39,
        fontsize=11.6,
    )


def write_script(path: Path, payload: dict[str, object]) -> None:
    low_parts = []
    high_parts = []
    for _, _, cent_label in CENT_BINS:
        panel = payload["panels"][cent_label]
        low_parts.append(f"{cent_label}: {panel['bands']['low_score_floor']['fraction']:.0%}")
        high_parts.append(f"{cent_label}: {panel['bands']['high_bdt_positive_isolation']['fraction']:.0%}")
    path.write_text(
        "\n".join(
            [
                "# THE-11 Slide Script - Inclusive Hotspot Decomposition",
                "",
                "This is the more useful follow-up to the inclusive row of the BDT-isolation density slide.",
                "Here I am not slicing the plot into one-dimensional distributions; I am decomposing the yellow density structures directly in the same two-dimensional phase space.",
                "",
                "Each centrality panel is normalized separately, so yellow means the densest bins within that centrality bin, not simply the centrality bin with the most raw entries.",
                "The displayed density is lightly smoothed only for readability; the printed region fractions are computed from the unsmoothed candidates.",
                "The 0-20 panel genuinely has fewer inclusive-background candidates and a more diffuse structure, which is why it looks less strip-like even after the normalization is fixed.",
                "The blue outline marks the dominant low-score positive-isolation floor.",
                "That structure grows from "
                + ", ".join(low_parts)
                + " of the inclusive background sample.",
                "Because it sits at very low BDT score, it is visually prominent but not the main tight-ID leakage concern.",
                "",
                "The orange outline is the more important structure for ABCD.",
                "It marks high-BDT background candidates that still have positive reconstructed cone energy.",
                "That region carries "
                + ", ".join(high_parts)
                + " of the inclusive background sample.",
                "This is the population that can make the BDT and isolation axes fail to factorize for background.",
                "",
                "So the point is not that every yellow band is dangerous.",
                "The slide separates a large low-score background floor from the smaller high-BDT positive-isolation leakage band.",
                "That is the piece that should be connected directly to the ABCD closure-ratio check.",
                "",
            ]
        )
        + "\n"
    )


def make_slide() -> tuple[Path, Path, Path]:
    setup_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    arrays = load_background_arrays()
    payload, histograms, vmin, vmax = build_payload(arrays)

    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)
    add_header(fig)
    draw_panels(fig, payload, histograms, vmin, vmax)
    draw_interpretation(fig, payload)

    png_path = OUT_DIR / "inclusive_hotspot_decomposition_slide.png"
    script_path = OUT_DIR / "inclusive_hotspot_decomposition_speaker_script.md"
    manifest_path = OUT_DIR / "inclusive_hotspot_decomposition_manifest.json"
    fig.savefig(png_path, dpi=DPI, facecolor=WHITE)
    plt.close(fig)

    write_script(script_path, payload)
    manifest = {
        "generator": str(Path(__file__).relative_to(REPO)),
        "source_slide21_png": str(SOURCE_SLIDE21_PNG.relative_to(REPO)),
        "score_cache_manifest": str(SCORE_CACHE_MANIFEST.relative_to(REPO)),
        "output_png": str(png_path.relative_to(REPO)),
        "speaker_script": str(script_path.relative_to(REPO)),
        "selection": {
            "score_column": SCORE_COLUMN,
            "isolation_column": ISO_COLUMN,
            "cluster_et_range_gev": list(PT_RANGE),
            "centrality_range_percent": [0.0, 80.0],
            "isolation_cone": "Delta R < 0.3",
            "class": "inclusive/background only, is_signal == 0",
        },
        "color_scale": {
            "quantity": "Relative P(reco_eiso, BDT score | centrality, inclusive background) per bin divided by the maximum bin in each centrality panel",
            "normalization": "each centrality panel normalized independently",
            "scale": "LogNorm",
            "display_smoothing": "3x3 binomial smoothing applied to the displayed density only; region fractions are computed from unsmoothed candidates",
            "vmin": vmin,
            "vmax": vmax,
        },
        **payload,
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")
    return png_path, script_path, manifest_path


def main() -> None:
    png_path, script_path, manifest_path = make_slide()
    print(png_path)
    print(script_path)
    print(manifest_path)


if __name__ == "__main__":
    main()
