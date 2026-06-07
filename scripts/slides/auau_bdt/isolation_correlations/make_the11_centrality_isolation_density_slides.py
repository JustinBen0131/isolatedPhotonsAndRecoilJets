#!/usr/bin/env python3
"""Build THE-11 centrality-aware BDT/isolation slide candidates."""

from __future__ import annotations

import csv
import json
import math
import textwrap
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402
from matplotlib.patches import FancyBboxPatch, Rectangle  # noqa: E402
import numpy as np  # noqa: E402


REPO = Path(__file__).resolve().parents[4]
REPORT = REPO / "dataOutput/auauTightBDTValidation/model_validation_condor_20260511_194832"
MANIFEST = REPORT / "score_caches.local.list"
CORR_CSV = REPO / (
    "dataOutput/auauMLDiagnosticRuns/global_etcent_inclusive3_sixpack_20260516_135439/"
    "slideReady/isolation_feature_correlations/noiso_score_diagnostic/"
    "isolation_feature_score_correlations.csv"
)
AUC_CSV = REPO / (
    "dataOutput/auauMLDiagnosticRuns/global_etcent_inclusive3_sixpack_20260516_135439/"
    "slideReady/eiso_cone_raw_comparison/eiso_raw_auc_summary.csv"
)
OUT_DIR = REPORT / "slideReady/default_eiso_vs_bdt_score_20260604/the11_centrality_isolation_story_20260604"

SCORE = "score_ptFine_cent7"
EISO = "reco_eiso"
PT_RANGE = (15.0, 35.0)
CENT_BINS = [(0.0, 20.0, "0-20%"), (20.0, 50.0, "20-50%"), (50.0, 80.0, "50-80%")]
ET_BINS_020 = [(15.0, 20.0, "15-20 GeV"), (20.0, 25.0, "20-25 GeV"), (25.0, 35.0, "25-35 GeV")]

W, H = 2560, 1440
DPI = 200

BG = "#ffffff"
INK = "#161616"
MUTED = "#59616f"
GRID = "#dfe4ea"
RED = "#c84c4c"
BLUE = "#315f9c"
PURPLE = "#6750a4"
TEAL = "#1f8a83"
GOLD = "#c9932e"
YELLOW = "#fff3c7"


@dataclass
class Peak:
    eiso: float
    score: float
    count: int


def setup_style() -> None:
    available = {f.name for f in font_manager.fontManager.ttflist}
    for name in ["Times New Roman", "Times", "DejaVu Serif"]:
        if name in available:
            plt.rcParams["font.family"] = name
            break
    plt.rcParams.update(
        {
            "figure.facecolor": BG,
            "axes.facecolor": "white",
            "savefig.facecolor": BG,
            "axes.edgecolor": INK,
            "axes.linewidth": 1.1,
            "axes.labelcolor": INK,
            "xtick.color": INK,
            "ytick.color": INK,
            "text.color": INK,
            "mathtext.default": "regular",
        }
    )


def corr(x: np.ndarray, y: np.ndarray) -> float:
    finite = np.isfinite(x) & np.isfinite(y)
    if finite.sum() < 3:
        return math.nan
    return float(np.corrcoef(x[finite].astype("float64"), y[finite].astype("float64"))[0, 1])


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
    arrays = {k: np.concatenate(v) for k, v in chunks.items()}
    mask = (
        np.isfinite(arrays["eiso"])
        & np.isfinite(arrays["score"])
        & (arrays["et"] >= PT_RANGE[0])
        & (arrays["et"] < PT_RANGE[1])
        & (arrays["cent"] >= 0.0)
        & (arrays["cent"] < 80.0)
    )
    return {k: v[mask] for k, v in arrays.items()}


def mask_cent(data: dict[str, np.ndarray], lo: float, hi: float) -> np.ndarray:
    return (data["cent"] >= lo) & (data["cent"] < hi)


def mask_et(data: dict[str, np.ndarray], lo: float, hi: float) -> np.ndarray:
    return (data["et"] >= lo) & (data["et"] < hi)


def percentile_rank(x: np.ndarray) -> np.ndarray:
    order = np.argsort(x)
    ranks = np.empty(len(x), dtype="float32")
    ranks[order] = np.linspace(0.0, 100.0, len(x), dtype="float32")
    return ranks


def peak_density(x: np.ndarray, y: np.ndarray, xlim: tuple[float, float]) -> Peak:
    h, xe, ye = np.histogram2d(x, y, bins=[70, 55], range=[xlim, [0.0, 1.0]])
    i, j = np.unravel_index(int(np.nanargmax(h)), h.shape)
    return Peak(float(0.5 * (xe[i] + xe[i + 1])), float(0.5 * (ye[j] + ye[j + 1])), int(h[i, j]))


def rounded_box(fig, xywh, facecolor="white", edgecolor="#d7dde7", lw=1.4, radius=0.016):
    patch = FancyBboxPatch(
        (xywh[0], xywh[1]),
        xywh[2],
        xywh[3],
        boxstyle=f"round,pad=0.010,rounding_size={radius}",
        transform=fig.transFigure,
        facecolor=facecolor,
        edgecolor=edgecolor,
        linewidth=lw,
        zorder=0.2,
    )
    fig.add_artist(patch)
    return patch


def add_header(fig, title: str, subtitle: str) -> None:
    fig.text(0.045, 0.935, title, ha="left", va="top", fontsize=27, fontweight="bold")
    fig.text(0.045, 0.895, subtitle, ha="left", va="top", fontsize=15.5, color=MUTED)
    fig.add_artist(Rectangle((0.045, 0.865), 0.91, 0.003, transform=fig.transFigure, color="#d7dce3", lw=0))


def add_bottom_cards(fig, cards: list[tuple[str, str, str]], y=0.050, h=0.140) -> None:
    x0, gap = 0.055, 0.020
    w = (0.89 - gap * (len(cards) - 1)) / len(cards)
    for i, (label, body, color) in enumerate(cards):
        x = x0 + i * (w + gap)
        rounded_box(fig, (x, y, w, h), facecolor="white", edgecolor="#d9dfe8")
        fig.text(x + 0.018, y + h - 0.030, label, fontsize=15.0, fontweight="bold", color=color, va="top")
        fig.text(
            x + 0.018,
            y + h - 0.066,
            textwrap.fill(body, width=48),
            fontsize=10.8,
            color=INK,
            va="top",
            linespacing=1.10,
        )


def style_panel(ax, title: str, xlim: tuple[float, float], xlabel: str | None = None, ylabel: str | None = None) -> None:
    ax.set_title(title, fontsize=14.0, fontweight="bold", pad=6)
    ax.set_xlim(*xlim)
    ax.set_ylim(0.0, 1.0)
    ax.grid(True, color=GRID, lw=0.7)
    ax.tick_params(labelsize=10.5)
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=12.5)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=12.5)


def hist_density(ax, x, y, xlim, title, cmap="magma", show_ylabel=False):
    h = ax.hist2d(x, y, bins=(72, 56), range=[xlim, [0.0, 1.0]], norm=LogNorm(vmin=1), cmap=cmap)
    style_panel(
        ax,
        title,
        xlim,
        xlabel=r"reco $E_T^{iso}$, $\Delta R < 0.3$ [GeV]",
        ylabel="BDT score" if show_ylabel else None,
    )
    return h


def hist_rank(ax, rank, score, title, show_ylabel=False):
    h = ax.hist2d(rank, score, bins=(72, 56), range=[[0.0, 100.0], [0.0, 1.0]], norm=LogNorm(vmin=1), cmap="cividis")
    style_panel(
        ax,
        title,
        (0.0, 100.0),
        xlabel="within-bin isolation percentile",
        ylabel="BDT score" if show_ylabel else None,
    )
    return h


def save_script(stem: str, title: str, paragraphs: list[str]) -> Path:
    path = OUT_DIR / f"{stem}_script.md"
    path.write_text("# THE-11 Script - " + title + "\n\n" + "\n\n".join(paragraphs) + "\n")
    return path


def slide_log_density_centrality(data: dict[str, np.ndarray], xlim: tuple[float, float]) -> tuple[Path, Path, dict]:
    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)
    add_header(
        fig,
        "Default isolation and BDT score have a centrality-dependent geometry",
        r"Row-level validation view, $15 \leq E_T^{cluster}<35$ GeV; default Au+Au isolation uses $\Delta R<0.3$.",
    )
    axes = [fig.add_axes([0.060 + i * 0.300, 0.285, 0.255, 0.510]) for i in range(3)]
    peaks, rhos = [], []
    last = None
    for i, (ax, (lo, hi, label)) in enumerate(zip(axes, CENT_BINS)):
        m = mask_cent(data, lo, hi)
        last = hist_density(ax, data["eiso"][m], data["score"][m], xlim, label, show_ylabel=(i == 0))
        peak = peak_density(data["eiso"][m], data["score"][m], xlim)
        rho = corr(data["eiso"][m], data["score"][m])
        peaks.append(peak)
        rhos.append(rho)
        ax.text(0.04, 0.94, rf"$\rho={rho:+.2f}$", transform=ax.transAxes, fontsize=12.5, color=INK, va="top", bbox={"facecolor": "white", "edgecolor": "#d9dfe8", "boxstyle": "round,pad=0.25"})
    cax = fig.add_axes([0.930, 0.315, 0.015, 0.450])
    fig.colorbar(last[3], cax=cax, label="candidates / bin")
    add_bottom_cards(
        fig,
        [
            ("What the color shows", "Log color shows where candidates live, avoiding the raw scatter-plot blob.", TEAL),
            ("Densest region", f"All bins peak near BDT score ~0.6; 0-20% peaks at reco Eiso {peaks[0].eiso:+.1f} GeV.", GOLD),
            ("Physical read", "Higher nearby activity mainly adds a broader low-score tail.", PURPLE),
        ],
    )
    stem = "slide01_r30_log_density_centrality"
    png = OUT_DIR / f"{stem}.png"
    fig.savefig(png)
    plt.close(fig)
    script = save_script(
        stem,
        "Centrality Log Density",
        [
            "This slide replaces the raw scatter with a log-density view, split into the three centrality bins.",
            "The important feature is not an individual point cloud. It is the geometry of the population: most candidates sit near the same BDT-score band, while the higher-isolation tail pulls the full-sample correlation negative.",
            f"In the 0-20 percent bin, the Pearson correlation is {rhos[0]:+.3f}; that means isolation and score are related, but not interchangeable.",
        ],
    )
    return png, script, {"peaks": [p.__dict__ for p in peaks], "rhos": rhos}


def slide_percentile_centrality(data: dict[str, np.ndarray]) -> tuple[Path, Path, dict]:
    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)
    add_header(
        fig,
        "Isolation percentile exposes the same score ordering without the scale shift",
        r"Each panel ranks reco $E_T^{iso}$ within its centrality bin, so the x-axis compares relative nearby activity.",
    )
    axes = [fig.add_axes([0.060 + i * 0.300, 0.285, 0.255, 0.510]) for i in range(3)]
    high_tail_scores = []
    last = None
    for i, (ax, (lo, hi, label)) in enumerate(zip(axes, CENT_BINS)):
        m = mask_cent(data, lo, hi)
        rank = percentile_rank(data["eiso"][m])
        last = hist_rank(ax, rank, data["score"][m], label, show_ylabel=(i == 0))
        high = data["score"][m][rank >= 80.0]
        high_tail_scores.append(float(np.nanmedian(high)))
        ax.axvline(80, color="white", lw=1.5, ls="--", alpha=0.85)
        ax.text(
            0.965,
            0.94,
            "top 20% iso",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=9.8,
            color=INK,
            bbox={"facecolor": "white", "edgecolor": "#d9dfe8", "boxstyle": "round,pad=0.20", "alpha": 0.88},
        )
    cax = fig.add_axes([0.930, 0.315, 0.015, 0.450])
    fig.colorbar(last[3], cax=cax, label="candidates / bin")
    add_bottom_cards(
        fig,
        [
            ("Why percentile helps", "Ranks isolation within each centrality bin, removing the raw scale shift.", TEAL),
            ("0-20% read", f"Top-isolation quintile median score is {high_tail_scores[0]:.2f}; it still overlaps the main band.", GOLD),
            ("Model implication", "Nearby activity is context for the BDT, not a replacement axis.", PURPLE),
        ],
    )
    stem = "slide02_r30_percentile_centrality"
    png = OUT_DIR / f"{stem}.png"
    fig.savefig(png)
    plt.close(fig)
    script = save_script(
        stem,
        "Isolation Percentile Centrality",
        [
            "Here I convert reconstructed isolation into a percentile within each centrality bin.",
            "That makes the comparison more robust to the fact that the raw isolation scale changes with event activity.",
            "The visible point is that even the high-isolation percentile region still overlaps the main BDT-score band, so isolation and the BDT score are correlated but not the same axis.",
        ],
    )
    return png, script, {"high_tail_median_score": high_tail_scores}


def slide_class_density_centrality(data: dict[str, np.ndarray], xlim: tuple[float, float]) -> tuple[Path, Path, dict]:
    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)
    add_header(
        fig,
        "Truth photons and inclusive jets occupy different score-isolation regions",
        r"Class densities are split before plotting so the inclusive-jet tail does not hide the photon band.",
    )
    axes = []
    for row in range(2):
        for col in range(3):
            axes.append(fig.add_axes([0.060 + col * 0.300, 0.555 - row * 0.275, 0.255, 0.190]))
    sig = data["is_signal"].astype(bool)
    medians = {"signal": [], "background": []}
    for col, (lo, hi, label) in enumerate(CENT_BINS):
        cm = mask_cent(data, lo, hi)
        for row, (class_mask, color_name, class_label, key) in enumerate(
            [(sig, "Reds", "Truth photons", "signal"), (~sig, "Blues", "Inclusive jets", "background")]
        ):
            ax = axes[row * 3 + col]
            m = cm & class_mask
            ax.hist2d(data["eiso"][m], data["score"][m], bins=(58, 44), range=[xlim, [0.0, 1.0]], norm=LogNorm(vmin=1), cmap=color_name)
            style_panel(
                ax,
                f"{label} | {class_label}",
                xlim,
                xlabel=r"reco $E_T^{iso}$, $\Delta R<0.3$ [GeV]" if row == 1 else None,
                ylabel="BDT score" if col == 0 else None,
            )
            medians[key].append(float(np.nanmedian(data["score"][m])))
    add_bottom_cards(
        fig,
        [
            ("Signal band", f"Truth photons stay high: median score {medians['signal'][0]:.2f} in 0-20%.", RED),
            ("Inclusive tail", f"Inclusive jets are broader/lower: median score {medians['background'][0]:.2f} in 0-20%.", BLUE),
            ("Physical read", "Isolation exposes where the jet-like tail enters the score plane.", PURPLE),
        ],
    )
    stem = "slide03_r30_class_density_centrality"
    png = OUT_DIR / f"{stem}.png"
    fig.savefig(png)
    plt.close(fig)
    script = save_script(
        stem,
        "Class Density Centrality",
        [
            "This slide splits the same two-dimensional space by truth label before plotting.",
            "The red truth-photon density stays concentrated at higher BDT score, while the blue inclusive-jet density is lower and broader, especially into the higher-isolation tail.",
            "That is the physical reason the full-sample correlation is visible without being a complete one-variable story.",
        ],
    )
    return png, script, medians


def slide_et_log_density_020(data: dict[str, np.ndarray], xlim: tuple[float, float]) -> tuple[Path, Path, dict]:
    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)
    add_header(
        fig,
        "ET dependence should be read first inside 0-20% centrality",
        r"Central-bin-only view, split by cluster $E_T$; this follows the thesis-priority plotting rule.",
    )
    base = mask_cent(data, 0.0, 20.0)
    axes = [fig.add_axes([0.060 + i * 0.300, 0.285, 0.255, 0.510]) for i in range(3)]
    peaks = []
    last = None
    for i, (ax, (lo, hi, label)) in enumerate(zip(axes, ET_BINS_020)):
        m = base & mask_et(data, lo, hi)
        last = hist_density(ax, data["eiso"][m], data["score"][m], xlim, label, show_ylabel=(i == 0))
        peak = peak_density(data["eiso"][m], data["score"][m], xlim)
        peaks.append(peak)
        ax.text(0.04, 0.94, f"n={int(m.sum()):,}", transform=ax.transAxes, fontsize=11.5, color=INK, va="top", bbox={"facecolor": "white", "edgecolor": "#d9dfe8", "boxstyle": "round,pad=0.25"})
    cax = fig.add_axes([0.930, 0.315, 0.015, 0.450])
    fig.colorbar(last[3], cax=cax, label="candidates / bin")
    add_bottom_cards(
        fig,
        [
            ("Why 0-20%", "ET trends should be read first in the thesis-priority central bin.", TEAL),
            ("Density peak", f"15-20 GeV peak: score {peaks[0].score:.2f}, reco Eiso {peaks[0].eiso:+.1f} GeV.", GOLD),
            ("Interpretation", "The score band persists; the low-score isolation tail changes.", PURPLE),
        ],
    )
    stem = "slide04_r30_et_slices_020_log_density"
    png = OUT_DIR / f"{stem}.png"
    fig.savefig(png)
    plt.close(fig)
    script = save_script(
        stem,
        "Central Bin ET Dependence",
        [
            "This slide applies the centrality-first rule directly: when we ask about E_T dependence, the first view is only 0-20 percent centrality.",
            "Across these three E_T intervals, the densest BDT-score band remains visible, while the isolation tail changes how much low-score inclusive structure appears.",
            "This is the safer way to talk about E_T dependence without mixing together different underlying-event regimes.",
        ],
    )
    return png, script, {"peaks": [p.__dict__ for p in peaks]}


def slide_et_class_density_020(data: dict[str, np.ndarray], xlim: tuple[float, float]) -> tuple[Path, Path, dict]:
    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)
    add_header(
        fig,
        "In the 0-20% bin, the class split explains the apparent isolation tail",
        r"Truth photons and inclusive jets are shown separately across $E_T$ slices.",
    )
    base = mask_cent(data, 0.0, 20.0)
    sig = data["is_signal"].astype(bool)
    axes = []
    for row in range(2):
        for col in range(3):
            axes.append(fig.add_axes([0.060 + col * 0.300, 0.555 - row * 0.275, 0.255, 0.190]))
    class_medians = {"signal": [], "background": []}
    for col, (lo, hi, label) in enumerate(ET_BINS_020):
        em = mask_et(data, lo, hi)
        for row, (class_mask, cmap, class_label, key) in enumerate(
            [(sig, "Reds", "Truth photons", "signal"), (~sig, "Blues", "Inclusive jets", "background")]
        ):
            m = base & em & class_mask
            ax = axes[row * 3 + col]
            ax.hist2d(data["eiso"][m], data["score"][m], bins=(58, 44), range=[xlim, [0.0, 1.0]], norm=LogNorm(vmin=1), cmap=cmap)
            style_panel(
                ax,
                f"{label} | {class_label}",
                xlim,
                xlabel=r"reco $E_T^{iso}$, $\Delta R<0.3$ [GeV]" if row == 1 else None,
                ylabel="BDT score" if col == 0 else None,
            )
            class_medians[key].append(float(np.nanmedian(data["score"][m])) if m.any() else math.nan)
    add_bottom_cards(
        fig,
        [
            ("Photon band", "Truth photons remain compact across ET slices in central events.", RED),
            ("Inclusive band", "Inclusive jets are lower/wider and drive the low-score structure.", BLUE),
            ("Use in the talk", "Best backup if asked whether this is only a mixture effect.", PURPLE),
        ],
    )
    stem = "slide05_r30_et_slices_020_class_density"
    png = OUT_DIR / f"{stem}.png"
    fig.savefig(png)
    plt.close(fig)
    script = save_script(
        stem,
        "Central Bin Class Density ET Dependence",
        [
            "This is the class-split version of the central-bin E_T check.",
            "The truth photons remain concentrated at high score, while the inclusive sample carries more of the broad low-score and high-isolation population.",
            "That distinction is important because it prevents us from overinterpreting the inclusive two-dimensional density as if it were a single physical class.",
        ],
    )
    return png, script, class_medians


def load_correlation_rows() -> list[dict[str, str]]:
    with CORR_CSV.open() as f:
        return list(csv.DictReader(f))


def corr_value(rows: list[dict[str, str]], scope: str, klass: str, iso: str, target: str) -> float:
    matches = [r for r in rows if r["scope"] == scope and r["class"] == klass and r["iso_variable"] == iso and r["target_variable"] == target]
    if not matches:
        return math.nan
    return float(matches[0]["pearson"])


def slide_synthesis(data: dict[str, np.ndarray], row_stats: dict, corr_rows: list[dict[str, str]]) -> tuple[Path, Path, dict]:
    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)
    add_header(
        fig,
        "The BDT-isolation relationship is real, but it is a population effect",
        r"Summary of the row-level $\Delta R<0.3$ density view plus exact R=0.3/R=0.4 aggregate correlation checks.",
    )
    ax1 = fig.add_axes([0.075, 0.300, 0.385, 0.470])
    labels = [c[2] for c in CENT_BINS]
    x = np.arange(3)
    all_rho = [corr(data["eiso"][mask_cent(data, lo, hi)], data["score"][mask_cent(data, lo, hi)]) for lo, hi, _ in CENT_BINS]
    sig = data["is_signal"].astype(bool)
    sig_rho = [corr(data["eiso"][mask_cent(data, lo, hi) & sig], data["score"][mask_cent(data, lo, hi) & sig]) for lo, hi, _ in CENT_BINS]
    bkg_rho = [corr(data["eiso"][mask_cent(data, lo, hi) & ~sig], data["score"][mask_cent(data, lo, hi) & ~sig]) for lo, hi, _ in CENT_BINS]
    ax1.axhline(0, color="#777", lw=1)
    ax1.bar(x - 0.22, all_rho, width=0.22, color=PURPLE, label="all")
    ax1.bar(x, sig_rho, width=0.22, color=RED, label="truth photons")
    ax1.bar(x + 0.22, bkg_rho, width=0.22, color=BLUE, label="inclusive jets")
    ax1.set_xticks(x)
    ax1.set_xticklabels(labels, fontsize=12)
    ax1.set_ylim(-0.35, 0.12)
    ax1.set_ylabel(r"Pearson $\rho$ with BDT score", fontsize=13)
    ax1.set_title(r"Row-level default $\Delta R<0.3$ isolation", fontsize=15, fontweight="bold")
    ax1.grid(axis="y", color=GRID)
    ax1.legend(frameon=False, fontsize=10)

    ax2 = fig.add_axes([0.550, 0.300, 0.345, 0.470])
    scopes = ["cent_0_20", "cent_20_50", "cent_50_80"]
    r30 = [corr_value(corr_rows, s, "all", "reco_eiso_r30", "score_globalEtCent1535_bdt_noIso_ptCent7") for s in scopes]
    r40 = [corr_value(corr_rows, s, "all", "reco_eiso_r40", "score_globalEtCent1535_bdt_noIso_ptCent7") for s in scopes]
    ax2.axhline(0, color="#777", lw=1)
    ax2.plot(x, r30, marker="o", ms=7, color=TEAL, lw=2.5, label=r"$\Delta R<0.3$")
    ax2.plot(x, r40, marker="s", ms=7, color=GOLD, lw=2.5, label=r"$\Delta R<0.4$")
    ax2.set_xticks(x)
    ax2.set_xticklabels(labels, fontsize=12)
    ax2.set_ylim(-0.70, 0.05)
    ax2.set_ylabel(r"Aggregate Pearson $\rho$", fontsize=13)
    ax2.set_title("Exact no-isolation-score aggregate check", fontsize=15, fontweight="bold")
    ax2.grid(axis="y", color=GRID)
    ax2.legend(frameon=False, fontsize=11)

    add_bottom_cards(
        fig,
        [
            ("Meaning of rho", "Negative rho means higher nearby energy tends to come with lower BDT score.", PURPLE),
            ("Where it is dense", "The photon-like band is near score ~0.6; the tail is lower and broader.", GOLD),
            ("Physics interpretation", "The BDT is not isolation, but nearby activity is entangled with its inputs.", TEAL),
        ],
    )
    stem = "slide06_r30_synthesis_and_cone_context"
    png = OUT_DIR / f"{stem}.png"
    fig.savefig(png)
    plt.close(fig)
    script = save_script(
        stem,
        "Synthesis",
        [
            "This summary slide puts the row-level density result next to the exact aggregate cone-size check.",
            "The left panel is the row-level default isolation view. The right panel uses the exact no-isolation-score aggregate diagnostic for R=0.3 and R=0.4.",
            "The conclusion is careful: the BDT score is not simply isolation, but it is not independent of nearby activity either. That relationship is strongest as a population effect involving the inclusive background tail.",
        ],
    )
    return png, script, {"row_level_rho": {"all": all_rho, "signal": sig_rho, "background": bkg_rho}, "aggregate_r30": r30, "aggregate_r40": r40}


def slide_backup_cone_correlations(corr_rows: list[dict[str, str]]) -> tuple[Path, Path, dict]:
    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)
    add_header(
        fig,
        "Backup: R=0.3 and R=0.4 give the same correlation story",
        r"Exact aggregate diagnostic: reco $E_T^{iso}$ versus the no-isolation BDT score.",
    )
    axes = [fig.add_axes([0.090 + i * 0.285, 0.315, 0.235, 0.455]) for i in range(3)]
    scopes = ["cent_0_20", "cent_20_50", "cent_50_80"]
    classes = ["all", "signal", "background"]
    colors = [PURPLE, RED, BLUE]
    values = {}
    for ax, scope, (_, _, label) in zip(axes, scopes, CENT_BINS):
        x = np.arange(len(classes))
        r30 = [corr_value(corr_rows, scope, c, "reco_eiso_r30", "score_globalEtCent1535_bdt_noIso_ptCent7") for c in classes]
        r40 = [corr_value(corr_rows, scope, c, "reco_eiso_r40", "score_globalEtCent1535_bdt_noIso_ptCent7") for c in classes]
        values[scope] = {"r30": r30, "r40": r40}
        ax.axhline(0, color="#777", lw=1)
        ax.bar(x - 0.17, r30, width=0.32, color=TEAL, label=r"$\Delta R<0.3$")
        ax.bar(x + 0.17, r40, width=0.32, color=GOLD, label=r"$\Delta R<0.4$")
        ax.set_xticks(x)
        ax.set_xticklabels(["all", "truth", "incl. jet"], fontsize=11)
        ax.set_ylim(-0.72, 0.10)
        ax.set_title(label, fontsize=15, fontweight="bold")
        ax.grid(axis="y", color=GRID)
        if ax is axes[0]:
            ax.set_ylabel(r"Pearson $\rho$", fontsize=13)
            ax.legend(frameon=False, fontsize=10)
    add_bottom_cards(
        fig,
        [
            ("Main result", "All candidates anticorrelate with score in every bin for both cone sizes.", TEAL),
            ("Signal check", "Truth photons alone are nearly flat; this is not a pure signal effect.", RED),
            ("Why backup", "Supports R=0.3 as mainline while keeping R=0.4 as a check.", GOLD),
        ],
    )
    stem = "backup01_r30_r40_correlation_by_centrality"
    png = OUT_DIR / f"{stem}.png"
    fig.savefig(png)
    plt.close(fig)
    script = save_script(
        stem,
        "Backup Cone Correlations",
        [
            "This backup slide uses the exact aggregate no-isolation-score diagnostic for both cone sizes.",
            "The important point is that R=0.3 and R=0.4 tell the same qualitative story: the full sample anticorrelates with BDT score, while the truth-photon-only correlation is much weaker.",
        ],
    )
    return png, script, values


def slide_backup_auc() -> tuple[Path, Path, dict]:
    rows = list(csv.DictReader(AUC_CSV.open()))
    keep = [
        ("Baseline routed BDT 32 inputs", "Baseline"),
        ("Baseline + raw R=0.3 $E_T^{iso}$", "+R=0.3"),
        ("Baseline + raw R=0.4 $E_T^{iso}$", "+R=0.4"),
        ("Baseline + raw R=0.3 and R=0.4 $E_T^{iso}$", "+both"),
    ]
    rowmap = {r["label"]: r for r in rows}
    vals = [[float(rowmap[label][f"auc_{cent}"]) for cent in ["0_20", "20_50", "50_80"]] for label, _ in keep]
    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)
    add_header(
        fig,
        "Backup: raw isolation is powerful when allowed as a BDT input",
        r"This is not the production-safe no-isolation BDT; it is the ceiling/check showing how much isolation carries.",
    )
    ax = fig.add_axes([0.090, 0.275, 0.800, 0.500])
    x = np.arange(3)
    width = 0.18
    colors = ["#6b7280", TEAL, GOLD, PURPLE]
    for i, ((_, label), v) in enumerate(zip(keep, vals)):
        ax.bar(x + (i - 1.5) * width, v, width=width, color=colors[i], label=label)
    ax.set_xticks(x)
    ax.set_xticklabels(["0-20%", "20-50%", "50-80%"], fontsize=13)
    ax.set_ylim(0.80, 0.97)
    ax.set_ylabel("AUC", fontsize=14)
    ax.set_title("AUC by centrality when raw isolation is included as an input", fontsize=16, fontweight="bold")
    ax.grid(axis="y", color=GRID)
    ax.legend(frameon=False, ncol=4, fontsize=11, loc="upper left")
    gain_020 = vals[1][0] - vals[0][0]
    add_bottom_cards(
        fig,
        [
            ("0-20% ceiling", f"Raw R=0.3 raises central-bin AUC by {gain_020:+.3f}.", TEAL),
            ("Cone-size read", "R=0.3 and R=0.4 are close; using both adds little beyond R=0.3.", GOLD),
            ("Why not mainline", "Diagnostic ceiling only; not the ABCD-safe no-isolation score.", PURPLE),
        ],
    )
    stem = "backup02_raw_cone_auc_ceiling"
    png = OUT_DIR / f"{stem}.png"
    fig.savefig(png)
    plt.close(fig)
    script = save_script(
        stem,
        "Backup AUC Ceiling",
        [
            "This backup shows what happens when raw isolation is allowed into the BDT input list.",
            "It is useful as a ceiling test: isolation contains strong separation information, especially in the 0-20 percent bin.",
            "It is not the production-safe photon-ID strategy because putting isolation into the BDT would complicate the ABCD purity logic.",
        ],
    )
    return png, script, {"auc": vals, "gain_0_20_r30": gain_020}


def slide_backup_cone_redundancy(corr_rows: list[dict[str, str]]) -> tuple[Path, Path, dict]:
    scopes = ["cent_0_20", "cent_20_50", "cent_50_80"]
    vals = [corr_value(corr_rows, s, "all", "reco_eiso_r30", "reco_eiso_r40") for s in scopes]
    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)
    add_header(
        fig,
        "Backup: R=0.3 and R=0.4 are highly redundant nearby-activity axes",
        r"Exact aggregate correlation between the two reconstructed isolation radii.",
    )
    ax = fig.add_axes([0.145, 0.310, 0.710, 0.455])
    x = np.arange(3)
    bars = ax.bar(x, vals, color=[TEAL, GOLD, PURPLE], width=0.55)
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v - 0.045, f"{v:.2f}", ha="center", va="top", fontsize=18, fontweight="bold", color="white")
    ax.set_xticks(x)
    ax.set_xticklabels(["0-20%", "20-50%", "50-80%"], fontsize=14)
    ax.set_ylim(0.0, 1.05)
    ax.set_ylabel(r"Pearson $\rho(R=0.3, R=0.4)$", fontsize=14)
    ax.set_title("Cone-size redundancy by centrality", fontsize=17, fontweight="bold")
    ax.grid(axis="y", color=GRID)
    add_bottom_cards(
        fig,
        [
            ("Redundancy", "The cone radii move together strongly; R=0.4 is a backup axis.", TEAL),
            ("Centrality stability", "The correlation stays high across all three centrality bins.", GOLD),
            ("Practical read", "Lead with R=0.3; use R=0.4 as a robustness check.", PURPLE),
        ],
    )
    stem = "backup03_r30_r40_redundancy"
    png = OUT_DIR / f"{stem}.png"
    fig.savefig(png)
    plt.close(fig)
    script = save_script(
        stem,
        "Backup Cone Redundancy",
        [
            "This backup answers why the main story can focus on R=0.3.",
            "R=0.3 and R=0.4 reconstructed isolation are very highly correlated in the existing aggregate diagnostic, so the larger cone is mainly a robustness and systematics check.",
        ],
    )
    return png, script, {"r30_r40_corr": vals}


def main() -> None:
    setup_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    data = load_arrays()
    xlim = tuple(float(v) for v in np.nanpercentile(data["eiso"], [0.5, 99.5]))
    corr_rows = load_correlation_rows()
    outputs: list[dict[str, object]] = []
    for maker in [
        lambda: slide_log_density_centrality(data, xlim),
        lambda: slide_percentile_centrality(data),
        lambda: slide_class_density_centrality(data, xlim),
        lambda: slide_et_log_density_020(data, xlim),
        lambda: slide_et_class_density_020(data, xlim),
        lambda: slide_synthesis(data, {}, corr_rows),
        lambda: slide_backup_cone_correlations(corr_rows),
        slide_backup_auc,
        lambda: slide_backup_cone_redundancy(corr_rows),
    ]:
        png, script, stats = maker()
        outputs.append({"png": str(png), "script": str(script), "stats": stats})
        print(f"wrote {png}")
        print(f"wrote {script}")

    manifest = OUT_DIR / "the11_centrality_isolation_story_manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "schema": "THE11_CENTRALITY_ISOLATION_STORY_V1",
                "source_report": str(REPORT),
                "score_cache_manifest": str(MANIFEST),
                "score_column": SCORE,
                "isolation_column": EISO,
                "row_level_selection": {
                    "cluster_Et_min_inclusive": PT_RANGE[0],
                    "cluster_Et_max_exclusive": PT_RANGE[1],
                    "centrality_min_inclusive": 0.0,
                    "centrality_max_exclusive": 80.0,
                    "selected_entries": int(len(data["score"])),
                    "signal_entries": int(data["is_signal"].sum()),
                    "background_entries": int((~data["is_signal"].astype(bool)).sum()),
                },
                "cone_label_basis": [
                    "Row-level density slides use local default reco_eiso caches.",
                    "RecoilJets_AuAu initializes m_isoConeR to 0.3 and eiso() evaluates eisoForCone(clus, m_isoConeR).",
                    "No local row-level reco_eiso_r40 score-cache shards were found, so R=0.4 is represented only by exact aggregate correlation/AUC backup evidence.",
                ],
                "aggregate_inputs": {
                    "correlation_csv": str(CORR_CSV),
                    "auc_csv": str(AUC_CSV),
                },
                "outputs": outputs,
                "google_slides_mutated": False,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    print(f"wrote {manifest}")


if __name__ == "__main__":
    main()
