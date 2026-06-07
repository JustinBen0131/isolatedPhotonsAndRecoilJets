#!/usr/bin/env python3
"""Make kBird/log-z density and signal-probability variants for THE-11."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402
from matplotlib.cm import ScalarMappable  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, LogNorm, Normalize  # noqa: E402
import numpy as np  # noqa: E402
from PIL import Image, ImageDraw, ImageFont  # noqa: E402


REPO = Path(__file__).resolve().parents[3]
REPORT = REPO / "dataOutput/auauTightBDTValidation/model_validation_condor_20260511_194832"
MANIFEST = REPORT / "score_caches.local.list"
OUT_DIR = REPORT / "slideReady/default_eiso_vs_bdt_score_20260604/plot_variants_kbird_probability_20260604"

SCORE = "score_ptFine_cent7"
EISO = "reco_eiso"
PT_RANGE = (15.0, 35.0)
CENT_BINS = [(0.0, 20.0, "0-20%"), (20.0, 50.0, "20-50%"), (50.0, 80.0, "50-80%")]

WIDE = (12.8, 7.2)
DPI = 200
INK = "#161616"
MUTED = "#59616f"
GRID = "#dfe4ea"
WHITE = "#ffffff"


def root_kbird() -> LinearSegmentedColormap:
    """ROOT kBird-style palette, close to TColor::kBird."""
    colors = [
        (0.2082, 0.1664, 0.5293),
        (0.0592, 0.3599, 0.8683),
        (0.0200, 0.5000, 0.9000),
        (0.0280, 0.6800, 0.7900),
        (0.1500, 0.7800, 0.6000),
        (0.4000, 0.8600, 0.3500),
        (0.7200, 0.9000, 0.2000),
        (0.9500, 0.9000, 0.2800),
        (0.9963, 0.9303, 0.5083),
    ]
    cmap = LinearSegmentedColormap.from_list("root_kbird_like", colors, N=256)
    cmap.set_bad("#f2f2f2")
    cmap.set_under("#f8f8f8")
    return cmap


KBIRD = root_kbird()


def setup_style() -> None:
    plt.rcParams.update(
        {
            "figure.facecolor": WHITE,
            "savefig.facecolor": WHITE,
            "axes.facecolor": WHITE,
            "axes.edgecolor": INK,
            "axes.linewidth": 1.05,
            "axes.labelcolor": INK,
            "xtick.color": INK,
            "ytick.color": INK,
            "text.color": INK,
            "font.family": ["Times New Roman", "DejaVu Serif"],
            "mathtext.default": "regular",
        }
    )


def load_arrays() -> dict[str, np.ndarray]:
    chunks: dict[str, list[np.ndarray]] = {k: [] for k in ["eiso", "score", "is_signal", "et", "cent"]}
    for line in MANIFEST.read_text().splitlines():
        if not line.strip():
            continue
        path = REPO / line.strip()
        data = np.load(path, allow_pickle=True)
        for key in [EISO, SCORE, "is_signal", "cluster_Et", "centrality"]:
            if key not in data.files:
                raise KeyError(f"{path} missing {key}")
        chunks["eiso"].append(data[EISO].astype("float32", copy=False))
        chunks["score"].append(data[SCORE].astype("float32", copy=False))
        chunks["is_signal"].append(data["is_signal"].astype("int8", copy=False))
        chunks["et"].append(data["cluster_Et"].astype("float32", copy=False))
        chunks["cent"].append(data["centrality"].astype("float32", copy=False))
    arrays = {key: np.concatenate(vals) for key, vals in chunks.items()}
    selected = (
        np.isfinite(arrays["eiso"])
        & np.isfinite(arrays["score"])
        & (arrays["et"] >= PT_RANGE[0])
        & (arrays["et"] < PT_RANGE[1])
        & (arrays["cent"] >= 0.0)
        & (arrays["cent"] < 80.0)
    )
    return {key: val[selected] for key, val in arrays.items()}


def add_title(fig, title: str, subtitle: str) -> None:
    fig.text(0.045, 0.945, title, ha="left", va="top", fontsize=22, fontweight="bold")
    fig.text(0.045, 0.905, subtitle, ha="left", va="top", fontsize=12.5, color=MUTED)


def style_ax(ax, title: str, xlim: tuple[float, float], xlabel: bool = True, ylabel: bool = True) -> None:
    ax.set_title(title, fontsize=12.5, fontweight="bold", pad=6)
    ax.set_xlim(*xlim)
    ax.set_ylim(0.0, 1.0)
    if xlabel:
        ax.set_xlabel(r"reco $E_T^{iso}$, $\Delta R<0.3$ [GeV]", fontsize=10.8)
    if ylabel:
        ax.set_ylabel("BDT score", fontsize=10.8)
    ax.grid(True, color=GRID, lw=0.65)
    ax.tick_params(labelsize=9.8)


def hist2d(x: np.ndarray, y: np.ndarray, xedges: np.ndarray, yedges: np.ndarray) -> np.ndarray:
    h, _, _ = np.histogram2d(x, y, bins=[xedges, yedges])
    return h.astype("float64")


def percentile_rank(values: np.ndarray) -> np.ndarray:
    if len(values) == 0:
        return np.empty(0, dtype="float32")
    if len(values) == 1:
        return np.full(1, 50.0, dtype="float32")
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype="float32")
    ranks[order] = np.linspace(0.0, 100.0, len(values), dtype="float32")
    return ranks


def probability_grid(data: dict[str, np.ndarray], mask: np.ndarray, xedges: np.ndarray, yedges: np.ndarray, min_count: int = 10):
    signal = data["is_signal"].astype(bool)
    counts = hist2d(data["eiso"][mask], data["score"][mask], xedges, yedges)
    sig_counts = hist2d(data["eiso"][mask & signal], data["score"][mask & signal], xedges, yedges)
    prob = np.divide(sig_counts, counts, out=np.full_like(counts, np.nan), where=counts >= min_count)
    return counts, prob


def add_log_count_contours(ax, counts: np.ndarray, xedges: np.ndarray, yedges: np.ndarray) -> None:
    positive = counts[counts > 0]
    if positive.size < 6:
        return
    levels = np.unique(np.nanpercentile(positive, [55, 75, 90, 97]))
    levels = levels[levels > 0]
    if len(levels) < 2:
        return
    xx = 0.5 * (xedges[:-1] + xedges[1:])
    yy = 0.5 * (yedges[:-1] + yedges[1:])
    ax.contour(xx, yy, counts.T, levels=levels, colors="black", linewidths=0.7, alpha=0.55)


def save_all_density_facets(data: dict[str, np.ndarray], xlim: tuple[float, float]) -> Path:
    xedges = np.linspace(xlim[0], xlim[1], 76)
    yedges = np.linspace(0.0, 1.0, 62)
    hists = []
    for lo, hi, _ in CENT_BINS:
        mask = (data["cent"] >= lo) & (data["cent"] < hi)
        hists.append(hist2d(data["eiso"][mask], data["score"][mask], xedges, yedges))
    vmax = max(float(np.nanmax(h)) for h in hists)
    norm = LogNorm(vmin=1.0, vmax=max(2.0, vmax))

    fig, axes = plt.subplots(1, 3, figsize=WIDE, dpi=DPI, sharex=True, sharey=True)
    fig.subplots_adjust(left=0.065, right=0.905, top=0.790, bottom=0.145, wspace=0.18)
    add_title(
        fig,
        "kBird log-density view highlights the populated BDT-isolation hot spots",
        r"All candidates, $15 \leq E_T^{cluster}<35$ GeV; shared log-z scale across centrality bins.",
    )
    mesh = None
    for ax, h, (_, _, label) in zip(axes, hists, CENT_BINS):
        mesh = ax.pcolormesh(xedges, yedges, h.T, cmap=KBIRD, norm=norm, shading="auto")
        style_ax(ax, label, xlim, ylabel=ax is axes[0])
    cax = fig.add_axes([0.925, 0.190, 0.018, 0.540])
    fig.colorbar(mesh, cax=cax, label="candidates / bin, log scale")
    out = OUT_DIR / "01_kbird_log_density_centrality_all.png"
    fig.savefig(out)
    plt.close(fig)
    return out


def save_class_density_facets(data: dict[str, np.ndarray], xlim: tuple[float, float]) -> Path:
    xedges = np.linspace(xlim[0], xlim[1], 68)
    yedges = np.linspace(0.0, 1.0, 56)
    signal = data["is_signal"].astype(bool)
    panel_hists: dict[str, list[np.ndarray]] = {"truth photons": [], "inclusive jets": []}
    panel_counts: dict[str, list[int]] = {"truth photons": [], "inclusive jets": []}
    panel_stats: dict[str, dict[str, list[float]]] = {
        "truth photons": {"score_median": [], "eiso_median": []},
        "inclusive jets": {"score_median": [], "eiso_median": []},
    }
    for class_label, class_mask in [("truth photons", signal), ("inclusive jets", ~signal)]:
        for lo, hi, _ in CENT_BINS:
            mask = class_mask & (data["cent"] >= lo) & (data["cent"] < hi)
            h = hist2d(data["eiso"][mask], data["score"][mask], xedges, yedges)
            panel_hists[class_label].append(h)
            panel_counts[class_label].append(int(np.count_nonzero(mask)))
            panel_stats[class_label]["score_median"].append(float(np.nanmedian(data["score"][mask])))
            panel_stats[class_label]["eiso_median"].append(float(np.nanmedian(data["eiso"][mask])))

    fig, axes = plt.subplots(2, 3, figsize=WIDE, dpi=DPI, sharex=True, sharey=True)
    fig.subplots_adjust(left=0.115, right=0.875, top=0.790, bottom=0.300, wspace=0.105, hspace=0.150)
    fig.text(
        0.050,
        0.950,
        "Truth photons and inclusive jets occupy distinct BDT-isolation regions",
        ha="left",
        va="top",
        fontsize=23.0,
        fontweight="bold",
    )
    fig.text(
        0.050,
        0.895,
        r"Au+Au simulation, $15 \leq E_T^{cluster}<35$ GeV. "
        "Color = relative density within each panel on a log scale; raw N is printed in each panel.",
        ha="left",
        va="top",
        fontsize=11.4,
        color=MUTED,
    )
    meshes = []
    class_colors = {"truth photons": "#9f2337", "inclusive jets": "#245a92"}
    norm = LogNorm(vmin=0.015, vmax=1.0)
    for row, class_label in enumerate(["truth photons", "inclusive jets"]):
        for col, ((_, _, cent_label), h, panel_n) in enumerate(
            zip(CENT_BINS, panel_hists[class_label], panel_counts[class_label])
        ):
            ax = axes[row, col]
            max_count = float(np.nanmax(h)) if np.any(h > 0.0) else 0.0
            relative_density = h / max_count if max_count > 0.0 else h
            relative_density = np.ma.masked_less_equal(relative_density, 0.0)
            mesh = ax.pcolormesh(xedges, yedges, relative_density.T, cmap=KBIRD, norm=norm, shading="auto")
            meshes.append(mesh)
            style_ax(ax, "", xlim, xlabel=False, ylabel=col == 0)
            if row == 0:
                ax.set_title(f"{cent_label} centrality", fontsize=12.6, fontweight="bold", pad=7)
            ax.set_facecolor("#fbfbfb")
            ax.text(
                0.035,
                0.925,
                f"N={panel_n:,}",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=8.8,
                color="#2f3540",
                bbox=dict(boxstyle="round,pad=0.18", facecolor="white", edgecolor="#d7dde8", linewidth=0.6, alpha=0.88),
            )
    fig.text(
        0.060,
        0.650,
        "truth photons",
        ha="center",
        va="center",
        rotation=90,
        fontsize=14.0,
        fontweight="bold",
        color=class_colors["truth photons"],
    )
    fig.text(
        0.060,
        0.400,
        "inclusive jets",
        ha="center",
        va="center",
        rotation=90,
        fontsize=14.0,
        fontweight="bold",
        color=class_colors["inclusive jets"],
    )
    fig.text(
        0.495,
        0.225,
        r"reconstructed $E_T^{iso}$, $\Delta R<0.3$ [GeV]",
        ha="center",
        va="center",
        fontsize=12.1,
    )
    cax0 = fig.add_axes([0.900, 0.575, 0.014, 0.175])
    cbar0 = fig.colorbar(meshes[0], cax=cax0)
    cbar0.set_label("truth relative density\nlog scale", fontsize=9.4)
    cbar0.ax.tick_params(labelsize=8.8)
    cax1 = fig.add_axes([0.900, 0.325, 0.014, 0.175])
    cbar1 = fig.colorbar(meshes[-1], cax=cax1)
    cbar1.set_label("inclusive relative density\nlog scale", fontsize=9.4)
    cbar1.ax.tick_params(labelsize=8.8)

    sig_score = panel_stats["truth photons"]["score_median"]
    bkg_score = panel_stats["inclusive jets"]["score_median"]
    sig_iso = panel_stats["truth photons"]["eiso_median"]
    bkg_iso = panel_stats["inclusive jets"]["eiso_median"]
    band = FancyBboxPatch(
        (0.085, 0.040),
        0.790,
        0.140,
        boxstyle="round,pad=0.010,rounding_size=0.008",
        transform=fig.transFigure,
        facecolor="#f8f8f5",
        edgecolor="#d7d0ad",
        linewidth=0.9,
    )
    fig.add_artist(band)
    fig.text(
        0.105,
        0.155,
        "Physical meaning of the BDT-isolation correlation",
        ha="left",
        va="center",
        fontsize=12.8,
        fontweight="bold",
    )
    fig.text(
        0.105,
        0.119,
        "Truth photons",
        ha="left",
        va="center",
        fontsize=10.7,
        fontweight="bold",
        color=class_colors["truth photons"],
    )
    fig.text(
        0.195,
        0.119,
        f"compact EM showers stay near zero cone energy "
        f"({min(sig_iso):.1f} to {max(sig_iso):.1f} GeV median) and high BDT score "
        f"({min(sig_score):.2f}-{max(sig_score):.2f} median).",
        ha="left",
        va="center",
        fontsize=10.3,
    )
    fig.text(
        0.105,
        0.088,
        "Inclusive jets",
        ha="left",
        va="center",
        fontsize=10.7,
        fontweight="bold",
        color=class_colors["inclusive jets"],
    )
    fig.text(
        0.195,
        0.088,
        f"nearby fragments add cone energy and make the shower less photon-like, "
        f"filling a lower-score positive-isolation tail ({min(bkg_score):.2f}-{max(bkg_score):.2f} median).",
        ha="left",
        va="center",
        fontsize=10.3,
    )
    fig.text(
        0.105,
        0.057,
        "Physical meaning",
        ha="left",
        va="center",
        fontsize=10.7,
        fontweight="bold",
    )
    fig.text(
        0.205,
        0.057,
        "the anti-correlation is a class-composition effect: both axes respond to hadronic activity around the EM cluster.",
        ha="left",
        va="center",
        fontsize=10.3,
    )
    out = OUT_DIR / "02_kbird_log_density_centrality_by_class.png"
    fig.savefig(out)
    plt.close(fig)
    return out


def save_isolation_rank_class_density_facets(data: dict[str, np.ndarray]) -> tuple[Path, dict[str, object]]:
    xedges = np.linspace(0.0, 100.0, 69)
    yedges = np.linspace(0.0, 1.0, 56)
    signal = data["is_signal"].astype(bool)
    panel_hists: dict[str, list[np.ndarray]] = {"truth photons": [], "inclusive jets": []}
    panel_stats: dict[str, dict[str, list[float]]] = {
        "truth photons": {"rank_median": [], "top20_fraction": [], "top20_score_median": []},
        "inclusive jets": {"rank_median": [], "top20_fraction": [], "top20_score_median": []},
    }
    eff80_cut_rank: list[float] = []
    eff80_inclusive_pass: list[float] = []

    for lo, hi, _ in CENT_BINS:
        cent_mask = (data["cent"] >= lo) & (data["cent"] < hi)
        cent_indices = np.flatnonzero(cent_mask)
        cent_eiso = data["eiso"][cent_mask]
        ranks = percentile_rank(data["eiso"][cent_mask])
        cent_signal = signal[cent_mask]
        iso_eff80 = float(np.nanpercentile(cent_eiso[cent_signal], 80.0))
        eff80_cut_rank.append(float(100.0 * np.mean(cent_eiso <= iso_eff80)))
        eff80_inclusive_pass.append(float(np.mean(cent_eiso[~cent_signal] <= iso_eff80)))
        for class_label, class_mask in [("truth photons", signal[cent_mask]), ("inclusive jets", ~signal[cent_mask])]:
            class_ranks = ranks[class_mask]
            class_scores = data["score"][cent_indices[class_mask]]
            panel_hists[class_label].append(hist2d(class_ranks, class_scores, xedges, yedges))
            high = class_ranks >= 80.0
            panel_stats[class_label]["rank_median"].append(float(np.nanmedian(class_ranks)))
            panel_stats[class_label]["top20_fraction"].append(float(np.mean(high)))
            panel_stats[class_label]["top20_score_median"].append(float(np.nanmedian(class_scores[high])) if np.any(high) else float("nan"))

    fig, axes = plt.subplots(2, 3, figsize=WIDE, dpi=DPI, sharex=True, sharey=True)
    fig.subplots_adjust(left=0.132, right=0.872, top=0.752, bottom=0.360, wspace=0.115, hspace=0.155)
    fig.text(
        0.050,
        0.958,
        "Isolation rank shows the BDT separates photons from busy jet cones",
        ha="left",
        va="top",
        fontsize=24.7,
        fontweight="bold",
    )
    subtitle_box = FancyBboxPatch(
        (0.047, 0.806),
        0.825,
        0.096,
        boxstyle="round,pad=0.010,rounding_size=0.007",
        transform=fig.transFigure,
        facecolor="#f4f6fb",
        edgecolor="#cbd5e1",
        linewidth=0.9,
    )
    fig.add_artist(subtitle_box)
    fig.text(0.065, 0.879, "Yellow band", ha="left", va="center", fontsize=12.4, fontweight="bold", color="#9a7400")
    fig.text(
        0.153,
        0.879,
        "= density hot spot: most candidates in that class/centrality panel, shown with log counts per bin.",
        ha="left",
        va="center",
        fontsize=11.9,
        color=MUTED,
    )
    fig.text(0.065, 0.854, "Rank axis", ha="left", va="center", fontsize=12.4, fontweight="bold", color=INK)
    fig.text(
        0.153,
        0.854,
        "= 0 cleanest cone to 100 busiest cone, computed separately within each centrality bin.",
        ha="left",
        va="center",
        fontsize=11.9,
        color=MUTED,
    )
    fig.text(0.065, 0.829, "Orange marker", ha="left", va="center", fontsize=12.4, fontweight="bold", color="#a45305")
    fig.text(
        0.153,
        0.829,
        "= approximate 80% photon-efficiency isolation point; the kept side is left/lower rank.",
        ha="left",
        va="center",
        fontsize=11.9,
        color=MUTED,
    )

    meshes = []
    class_colors = {"truth photons": "#9f2337", "inclusive jets": "#245a92"}
    for row, class_label in enumerate(["truth photons", "inclusive jets"]):
        row_vmax = max(float(np.nanmax(h)) for h in panel_hists[class_label])
        norm = LogNorm(vmin=1.0, vmax=max(2.0, row_vmax))
        for col, ((_, _, cent_label), h) in enumerate(zip(CENT_BINS, panel_hists[class_label])):
            ax = axes[row, col]
            mesh = ax.pcolormesh(xedges, yedges, h.T, cmap=KBIRD, norm=norm, shading="auto")
            meshes.append(mesh)
            ax.set_xlim(0.0, 100.0)
            ax.set_ylim(0.0, 1.0)
            ax.grid(True, color=GRID, lw=0.65)
            ax.tick_params(labelsize=10.7)
            if row == 0:
                ax.set_title(f"{cent_label} centrality", fontsize=14.4, fontweight="bold", pad=7)
            cut_rank = eff80_cut_rank[col]
            ax.plot(
                [cut_rank],
                [0.018],
                marker="v",
                markersize=7.2,
                color="#d87306",
                markeredgecolor="white",
                markeredgewidth=0.7,
                zorder=5,
            )
            if row == 1:
                ax.annotate(
                    "",
                    xy=(6.0, -0.145),
                    xytext=(cut_rank, -0.145),
                    xycoords=ax.get_xaxis_transform(),
                    textcoords=ax.get_xaxis_transform(),
                    arrowprops=dict(arrowstyle="->", color="#b96005", lw=1.15, shrinkA=0, shrinkB=0),
                    clip_on=False,
                )
                ax.text(
                    0.5 * (6.0 + cut_rank),
                    -0.200,
                    "kept by isolation cut",
                    transform=ax.get_xaxis_transform(),
                    ha="center",
                    va="top",
                    fontsize=9.4,
                    color="#8a4a05",
                    clip_on=False,
                )
            ax.set_facecolor("#fbfbfb")

    fig.text(
        0.044,
        0.650,
        "truth photons",
        ha="center",
        va="center",
        rotation=90,
        fontsize=16.3,
        fontweight="bold",
        color=class_colors["truth photons"],
    )
    fig.text(
        0.044,
        0.400,
        "inclusive jets",
        ha="center",
        va="center",
        rotation=90,
        fontsize=16.3,
        fontweight="bold",
        color=class_colors["inclusive jets"],
    )
    fig.text(
        0.091,
        0.565,
        "BDT photon-likeness score",
        ha="center",
        va="center",
        rotation=90,
        fontsize=13.0,
    )
    fig.text(
        0.495,
        0.274,
        r"within-centrality isolation percentile rank  (0 = cleanest cone, 100 = busiest cone)",
        ha="center",
        va="center",
        fontsize=13.2,
    )
    cax0 = fig.add_axes([0.900, 0.580, 0.014, 0.165])
    cbar0 = fig.colorbar(meshes[0], cax=cax0)
    cbar0.set_label("truth photons per bin\nlog scale", fontsize=10.3)
    cbar0.ax.tick_params(labelsize=9.6)
    cax1 = fig.add_axes([0.900, 0.375, 0.014, 0.165])
    cbar1 = fig.colorbar(meshes[-1], cax=cax1)
    cbar1.set_label("inclusive jets per bin\nlog scale", fontsize=10.3)
    cbar1.ax.tick_params(labelsize=9.6)

    band = FancyBboxPatch(
        (0.065, 0.038),
        0.855,
        0.188,
        boxstyle="round,pad=0.010,rounding_size=0.008",
        transform=fig.transFigure,
        facecolor="#faf9f1",
        edgecolor="#d7d0ad",
        linewidth=0.9,
    )
    fig.add_artist(band)
    fig.text(0.090, 0.198, "What the rank is", ha="left", va="center", fontsize=15.0, fontweight="bold")
    fig.text(
        0.090,
        0.116,
        r"Sort by reco $E_T^{iso}$ inside" "\n"
        "each centrality bin.\n"
        "Rank 0 is cleanest; 100 is busiest.\n"
        r"Formula: 100 $\times$ sorted position / (N-1).",
        ha="left",
        va="center",
        fontsize=12.0,
        linespacing=1.26,
    )
    fig.text(
        0.355,
        0.198,
        "Interpretation",
        ha="left",
        va="center",
        fontsize=15.0,
        fontweight="bold",
        color="#a45305",
    )
    fig.text(
        0.355,
        0.116,
        "Low rank means little nearby energy.\n"
        "High rank means a busier cone.\n"
        "Yellow photon band stays high-score.\n"
        "Yellow jet density piles up busy-side.",
        ha="left",
        va="center",
        fontsize=12.0,
        linespacing=1.26,
    )
    fig.text(
        0.635,
        0.198,
        "BDT interpretation",
        ha="left",
        va="center",
        fontsize=15.0,
        fontweight="bold",
        color=class_colors["inclusive jets"],
    )
    fig.text(
        0.635,
        0.116,
        "Isolation removes the busy tail.\n"
        "Truth photons remain high-score.\n"
        "Jets spread lower as cones get busy.\n"
        "BDT and isolation are complementary.",
        ha="left",
        va="center",
        fontsize=12.0,
        linespacing=1.26,
    )

    out = OUT_DIR / "08_kbird_isolation_percentile_rank_by_class.png"
    fig.savefig(out)
    plt.close(fig)
    panel_stats["eff80_isolation_cut"] = {
        "rank_in_all_candidates": eff80_cut_rank,
        "inclusive_pass_fraction": eff80_inclusive_pass,
        "note": "Computed from the 80th percentile of truth-photon reco_eiso in each centrality bin for this R<0.3 diagnostic.",
    }
    return out, panel_stats


def save_probability_facets(data: dict[str, np.ndarray], xlim: tuple[float, float]) -> Path:
    xedges = np.linspace(xlim[0], xlim[1], 68)
    yedges = np.linspace(0.0, 1.0, 56)
    fig, axes = plt.subplots(1, 3, figsize=WIDE, dpi=DPI, sharex=True, sharey=True)
    fig.subplots_adjust(left=0.065, right=0.905, top=0.790, bottom=0.145, wspace=0.18)
    add_title(
        fig,
        "Signal-probability map shows class composition at each BDT-isolation location",
        "Color is truth-photon fraction; black contours are log-density occupancy hot spots.",
    )
    norm = Normalize(vmin=0.0, vmax=1.0)
    mesh = None
    medians = []
    for ax, (lo, hi, label) in zip(axes, CENT_BINS):
        mask = (data["cent"] >= lo) & (data["cent"] < hi)
        counts, prob = probability_grid(data, mask, xedges, yedges, min_count=10)
        mesh = ax.pcolormesh(xedges, yedges, np.ma.masked_invalid(prob.T), cmap=KBIRD, norm=norm, shading="auto")
        add_log_count_contours(ax, counts, xedges, yedges)
        style_ax(ax, label, xlim, ylabel=ax is axes[0])
        medians.append(float(np.nanmedian(prob)))
    cax = fig.add_axes([0.925, 0.190, 0.018, 0.540])
    fig.colorbar(mesh, cax=cax, label="truth-photon fraction in populated bin")
    out = OUT_DIR / "03_kbird_signal_probability_centrality_with_log_density_contours.png"
    fig.savefig(out)
    plt.close(fig)
    return out


def save_probability_alpha_facets(data: dict[str, np.ndarray], xlim: tuple[float, float]) -> Path:
    xedges = np.linspace(xlim[0], xlim[1], 68)
    yedges = np.linspace(0.0, 1.0, 56)
    fig, axes = plt.subplots(1, 3, figsize=WIDE, dpi=DPI, sharex=True, sharey=True)
    fig.subplots_adjust(left=0.065, right=0.905, top=0.790, bottom=0.145, wspace=0.18)
    add_title(
        fig,
        "Probability map with log-occupancy masking suppresses empty-bin noise",
        "Color is truth-photon fraction; only populated bins are shown, with stronger opacity in denser bins.",
    )
    norm = Normalize(vmin=0.0, vmax=1.0)
    mesh = ScalarMappable(norm=norm, cmap=KBIRD)
    mesh.set_array([])
    for ax, (lo, hi, label) in zip(axes, CENT_BINS):
        mask = (data["cent"] >= lo) & (data["cent"] < hi)
        counts, prob = probability_grid(data, mask, xedges, yedges, min_count=8)
        alpha = np.zeros_like(counts)
        if np.nanmax(counts) > 0:
            alpha = 0.22 + 0.78 * np.log1p(counts) / np.log1p(np.nanmax(counts))
        alpha = np.where(np.isfinite(prob), alpha, 0.0)
        rgba = KBIRD(norm(np.nan_to_num(prob.T, nan=0.0)))
        rgba[..., -1] = alpha.T
        ax.imshow(
            rgba,
            origin="lower",
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
            aspect="auto",
            interpolation="nearest",
        )
        add_log_count_contours(ax, counts, xedges, yedges)
        style_ax(ax, label, xlim, ylabel=ax is axes[0])
    cax = fig.add_axes([0.925, 0.190, 0.018, 0.540])
    fig.colorbar(mesh, cax=cax, label="truth-photon fraction")
    out = OUT_DIR / "04_kbird_signal_probability_centrality_log_occupancy_alpha.png"
    fig.savefig(out)
    plt.close(fig)
    return out


def save_probability_contrast_facets(data: dict[str, np.ndarray], xlim: tuple[float, float]) -> Path:
    xedges = np.linspace(xlim[0], xlim[1], 68)
    yedges = np.linspace(0.0, 1.0, 56)
    grids = []
    counts_list = []
    for lo, hi, _ in CENT_BINS:
        mask = (data["cent"] >= lo) & (data["cent"] < hi)
        counts, prob = probability_grid(data, mask, xedges, yedges, min_count=10)
        grids.append(prob)
        counts_list.append(counts)
    all_probs = np.concatenate([g[np.isfinite(g)] for g in grids])
    vmin = max(0.0, float(np.nanpercentile(all_probs, 5.0)))
    vmax = min(1.0, float(np.nanpercentile(all_probs, 98.0)))
    if vmax - vmin < 0.08:
        vmin, vmax = 0.55, 1.0

    fig, axes = plt.subplots(1, 3, figsize=WIDE, dpi=DPI, sharex=True, sharey=True)
    fig.subplots_adjust(left=0.065, right=0.905, top=0.790, bottom=0.145, wspace=0.18)
    add_title(
        fig,
        "Contrast-stretched probability view separates the high-purity band",
        f"Color is truth-photon fraction, stretched over populated-bin range {vmin:.2f}-{vmax:.2f}; black contours are log-density hot spots.",
    )
    norm = Normalize(vmin=vmin, vmax=vmax)
    mesh = None
    for ax, prob, counts, (_, _, label) in zip(axes, grids, counts_list, CENT_BINS):
        mesh = ax.pcolormesh(xedges, yedges, np.ma.masked_invalid(prob.T), cmap=KBIRD, norm=norm, shading="auto")
        add_log_count_contours(ax, counts, xedges, yedges)
        style_ax(ax, label, xlim, ylabel=ax is axes[0])
    cax = fig.add_axes([0.925, 0.190, 0.018, 0.540])
    fig.colorbar(mesh, cax=cax, label="truth-photon fraction, contrast-stretched")
    out = OUT_DIR / "06_kbird_signal_probability_centrality_contrast_stretched.png"
    fig.savefig(out)
    plt.close(fig)
    return out


def save_inclusive_probability_log_facets(data: dict[str, np.ndarray], xlim: tuple[float, float]) -> Path:
    xedges = np.linspace(xlim[0], xlim[1], 68)
    yedges = np.linspace(0.0, 1.0, 56)
    fig, axes = plt.subplots(1, 3, figsize=WIDE, dpi=DPI, sharex=True, sharey=True)
    fig.subplots_adjust(left=0.065, right=0.905, top=0.790, bottom=0.145, wspace=0.18)
    add_title(
        fig,
        "Log-scaled inclusive-jet probability highlights the background tail",
        "Color is inclusive-jet fraction in populated bins, log-scaled so small but structured background components remain visible.",
    )
    norm = LogNorm(vmin=0.01, vmax=1.0)
    mesh = None
    for ax, (lo, hi, label) in zip(axes, CENT_BINS):
        mask = (data["cent"] >= lo) & (data["cent"] < hi)
        counts, prob_signal = probability_grid(data, mask, xedges, yedges, min_count=10)
        prob_inclusive = 1.0 - prob_signal
        prob_inclusive = np.where(prob_inclusive >= 0.01, prob_inclusive, np.nan)
        mesh = ax.pcolormesh(
            xedges,
            yedges,
            np.ma.masked_invalid(prob_inclusive.T),
            cmap=KBIRD,
            norm=norm,
            shading="auto",
        )
        add_log_count_contours(ax, counts, xedges, yedges)
        style_ax(ax, label, xlim, ylabel=ax is axes[0])
    cax = fig.add_axes([0.925, 0.190, 0.018, 0.540])
    fig.colorbar(mesh, cax=cax, label="inclusive-jet fraction, log scale")
    out = OUT_DIR / "07_kbird_inclusive_probability_centrality_log_scale.png"
    fig.savefig(out)
    plt.close(fig)
    return out


def save_central_bin_four_panel(data: dict[str, np.ndarray], xlim: tuple[float, float]) -> Path:
    xedges = np.linspace(xlim[0], xlim[1], 72)
    yedges = np.linspace(0.0, 1.0, 58)
    cent = (data["cent"] >= 0.0) & (data["cent"] < 20.0)
    signal = data["is_signal"].astype(bool)
    panels = [
        ("all candidates log density", cent),
        ("truth photons log density", cent & signal),
        ("inclusive jets log density", cent & ~signal),
    ]
    density_hists = [hist2d(data["eiso"][m], data["score"][m], xedges, yedges) for _, m in panels]
    vmax = max(float(np.nanmax(h)) for h in density_hists)
    density_norm = LogNorm(vmin=1.0, vmax=max(2.0, vmax))
    counts, prob = probability_grid(data, cent, xedges, yedges, min_count=8)

    fig, axes = plt.subplots(2, 2, figsize=WIDE, dpi=DPI, sharex=True, sharey=True)
    fig.subplots_adjust(left=0.070, right=0.900, top=0.805, bottom=0.135, wspace=0.16, hspace=0.28)
    add_title(
        fig,
        "0-20% centrality: density and probability views of the same structure",
        "Central-bin focus: log-density hot spots plus the local truth-photon probability map.",
    )
    meshes = []
    for ax, (title, _), h in zip(axes.ravel()[:3], panels, density_hists):
        mesh = ax.pcolormesh(xedges, yedges, h.T, cmap=KBIRD, norm=density_norm, shading="auto")
        meshes.append(mesh)
        style_ax(ax, title, xlim, xlabel=ax in axes[1], ylabel=ax in axes[:, 0])
    prob_mesh = axes[1, 1].pcolormesh(
        xedges,
        yedges,
        np.ma.masked_invalid(prob.T),
        cmap=KBIRD,
        norm=Normalize(vmin=0.0, vmax=1.0),
        shading="auto",
    )
    add_log_count_contours(axes[1, 1], counts, xedges, yedges)
    style_ax(axes[1, 1], "truth-photon probability + log-density contours", xlim, xlabel=True, ylabel=False)
    cax0 = fig.add_axes([0.922, 0.475, 0.016, 0.280])
    fig.colorbar(meshes[0], cax=cax0, label="candidates / bin, log")
    cax1 = fig.add_axes([0.922, 0.165, 0.016, 0.280])
    fig.colorbar(prob_mesh, cax=cax1, label="truth-photon fraction")
    out = OUT_DIR / "05_kbird_020_density_probability_four_panel.png"
    fig.savefig(out)
    plt.close(fig)
    return out


def make_contact_sheet(paths: list[Path]) -> Path:
    thumbs: list[Image.Image] = []
    for path in paths:
        img = Image.open(path).convert("RGB")
        img.thumbnail((760, 430), Image.Resampling.LANCZOS)
        canvas = Image.new("RGB", (800, 500), WHITE)
        canvas.paste(img, ((800 - img.width) // 2, 20))
        draw = ImageDraw.Draw(canvas)
        try:
            font = ImageFont.truetype("Arial.ttf", 20)
        except OSError:
            font = ImageFont.load_default()
        draw.text((24, 462), path.name, fill=(20, 20, 20), font=font)
        thumbs.append(canvas)
    rows = (len(thumbs) + 1) // 2
    sheet = Image.new("RGB", (1600, 500 * rows), WHITE)
    for idx, img in enumerate(thumbs):
        sheet.paste(img, ((idx % 2) * 800, (idx // 2) * 500))
    out = OUT_DIR / "00_contact_sheet_kbird_probability_variants.png"
    sheet.save(out)
    return out


def main() -> None:
    setup_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    data = load_arrays()
    xlim = tuple(float(v) for v in np.nanpercentile(data["eiso"], [0.5, 99.5]))
    rank_path, rank_stats = save_isolation_rank_class_density_facets(data)
    paths = [
        save_all_density_facets(data, xlim),
        save_class_density_facets(data, xlim),
        save_probability_facets(data, xlim),
        save_probability_alpha_facets(data, xlim),
        save_central_bin_four_panel(data, xlim),
        save_probability_contrast_facets(data, xlim),
        save_inclusive_probability_log_facets(data, xlim),
        rank_path,
    ]
    sheet = make_contact_sheet(paths)
    manifest = OUT_DIR / "kbird_probability_variants_manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "schema": "THE11_KBIRD_PROBABILITY_VARIANTS_V1",
                "source_report": str(REPORT),
                "score_cache_manifest": str(MANIFEST),
                "score_column": SCORE,
                "isolation_column": EISO,
                "selection": {
                    "cluster_Et_min_inclusive": PT_RANGE[0],
                    "cluster_Et_max_exclusive": PT_RANGE[1],
                    "centrality_min_inclusive": 0.0,
                    "centrality_max_exclusive": 80.0,
                    "selected_entries": int(len(data["score"])),
                    "signal_entries": int(data["is_signal"].sum()),
                    "background_entries": int((~data["is_signal"].astype(bool)).sum()),
                },
                "visual_encoding": {
                    "density": "ROOT-like kBird palette with LogNorm. 01 uses raw log occupancy; 02 class-density panels use per-panel relative density so color shows shape/hot spots while raw N is printed in each panel.",
                    "probability": "truth-photon fraction in populated bins using kBird; log-count structure shown by contours or occupancy alpha",
                    "isolation_rank": "reco_eiso percentile rank computed separately within each centrality bin before class splitting",
                    "canvas": "white",
                },
                "isolation_rank_panel_stats": rank_stats,
                "xlim_percentile_0p5_99p5": list(xlim),
                "outputs": [str(p) for p in [sheet, *paths]],
                "google_slides_mutated": False,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    print(f"wrote {sheet}")
    for path in paths:
        print(f"wrote {path}")
    print(f"wrote {manifest}")


if __name__ == "__main__":
    main()
