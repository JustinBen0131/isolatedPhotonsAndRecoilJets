#!/usr/bin/env python3
"""Build the THE-11 truth-photon BDT-score isolation-slice follow-up slide."""

from __future__ import annotations

import json
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager, ticker  # noqa: E402
from matplotlib.patches import FancyBboxPatch, Rectangle  # noqa: E402
import numpy as np  # noqa: E402


REPO = Path(__file__).resolve().parents[4]
SOURCE_REPORT = REPO / "dataOutput/auauTightBDTValidation/model_validation_condor_20260511_194832"
SCORE_CACHE_MANIFEST = SOURCE_REPORT / "score_caches.local.list"
SOURCE_SLIDE21_GENERATOR = REPO / "scripts/plotting/auau_bdt/make_default_eiso_bdt_kbird_probability_variants.py"
SOURCE_SLIDE21_PNG = (
    SOURCE_REPORT
    / "slideReady/default_eiso_vs_bdt_score_20260604/plot_variants_kbird_probability_20260604/02_kbird_log_density_centrality_by_class.png"
)
OUT_DIR = (
    SOURCE_REPORT
    / "slideReady/default_eiso_vs_bdt_score_20260604/the11_truth_photon_bdt_isolation_slices_20260605"
)

SCORE_COLUMN = "score_ptFine_cent7"
ISO_COLUMN = "reco_eiso"
PT_RANGE = (15.0, 35.0)
CENT_BINS = [(0.0, 20.0, "0-20%"), (20.0, 50.0, "20-50%"), (50.0, 80.0, "50-80%")]
ISO_SLICES = [
    ("negative", r"$E_T^{iso}<0$", lambda eiso: eiso < 0.0, "#2166ac"),
    ("near_neutral", r"$0 \leq E_T^{iso}<2$", lambda eiso: (eiso >= 0.0) & (eiso < 2.0), "#1b9e77"),
    ("moderate", r"$2 \leq E_T^{iso}<5$", lambda eiso: (eiso >= 2.0) & (eiso < 5.0), "#e69f00"),
    ("busy_tail", r"$E_T^{iso}\geq5$", lambda eiso: eiso >= 5.0, "#7b3294"),
]
BDT_BINS = np.linspace(0.0, 1.0, 42)

W, H = 2560, 1440
DPI = 200
WHITE = "#ffffff"
INK = "#171717"
MUTED = "#59616f"
GRID = "#dfe4ea"
PANEL_BG = "#fbfbfb"
SUBTITLE_BG = "#f4f7fb"
SUBTITLE_EDGE = "#cad6e4"
BAND_BG = "#fff8dc"
BAND_EDGE = "#d8c983"
PHOTON_BAND = "#fff1a8"


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
            "axes.facecolor": PANEL_BG,
            "axes.edgecolor": INK,
            "axes.linewidth": 1.2,
            "axes.labelcolor": INK,
            "xtick.color": INK,
            "ytick.color": INK,
            "text.color": INK,
            "mathtext.default": "regular",
        }
    )


def load_arrays() -> dict[str, np.ndarray]:
    if not SCORE_CACHE_MANIFEST.exists():
        raise FileNotFoundError(f"Missing score cache list: {SCORE_CACHE_MANIFEST}")

    chunks: dict[str, list[np.ndarray]] = {
        "eiso": [],
        "score": [],
        "is_signal": [],
        "cluster_et": [],
        "centrality": [],
    }
    for raw_line in SCORE_CACHE_MANIFEST.read_text().splitlines():
        line = raw_line.strip()
        if not line:
            continue
        path = Path(line)
        cache_path = path if path.is_absolute() else REPO / path
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
        & (arrays["cluster_et"] >= PT_RANGE[0])
        & (arrays["cluster_et"] < PT_RANGE[1])
        & (arrays["centrality"] >= 0.0)
        & (arrays["centrality"] < 80.0)
    )
    return {key: values[selected] for key, values in arrays.items()}


def describe_slice(score: np.ndarray, total_cent: int, hist: np.ndarray, edges: np.ndarray) -> dict[str, float | int]:
    if len(score) == 0:
        return {
            "n": 0,
            "fraction_within_centrality": 0.0,
            "median_bdt": float("nan"),
            "q25_bdt": float("nan"),
            "q75_bdt": float("nan"),
            "fraction_bdt_gt_0p5": float("nan"),
            "fraction_bdt_gt_0p6": float("nan"),
            "peak_bin_center": float("nan"),
        }
    peak = int(np.nanargmax(hist))
    return {
        "n": int(len(score)),
        "fraction_within_centrality": float(len(score) / total_cent) if total_cent else 0.0,
        "median_bdt": float(np.nanmedian(score)),
        "q25_bdt": float(np.nanquantile(score, 0.25)),
        "q75_bdt": float(np.nanquantile(score, 0.75)),
        "fraction_bdt_gt_0p5": float(np.mean(score > 0.5)),
        "fraction_bdt_gt_0p6": float(np.mean(score > 0.6)),
        "peak_bin_center": float(0.5 * (edges[peak] + edges[peak + 1])),
    }


def build_histograms(arrays: dict[str, np.ndarray]) -> tuple[dict[str, dict[str, object]], float]:
    signal = arrays["is_signal"].astype(bool)
    hists: dict[str, dict[str, object]] = {}
    ymax = 0.0
    for cent_lo, cent_hi, cent_label in CENT_BINS:
        cent_mask = signal & (arrays["centrality"] >= cent_lo) & (arrays["centrality"] < cent_hi)
        cent_total = int(np.sum(cent_mask))
        cent_payload: dict[str, object] = {
            "centrality": cent_label,
            "truth_photon_n": cent_total,
            "slices": {},
        }
        for slice_key, _, selector, _ in ISO_SLICES:
            mask = cent_mask & selector(arrays["eiso"])
            scores = arrays["score"][mask]
            counts, edges = np.histogram(scores, bins=BDT_BINS)
            frac = counts.astype("float64")
            if frac.sum() > 0:
                frac /= frac.sum()
            ymax = max(ymax, float(np.nanmax(frac)) if len(frac) else 0.0)
            cent_payload["slices"][slice_key] = {
                "hist_fraction": frac.tolist(),
                "bin_edges": edges.tolist(),
                "summary": describe_slice(scores, cent_total, frac, edges),
            }
        hists[cent_label] = cent_payload
    return hists, ymax


def add_box(fig: plt.Figure, xywh: tuple[float, float, float, float], face: str, edge: str, zorder: float = 0.2) -> None:
    fig.add_artist(
        FancyBboxPatch(
            (xywh[0], xywh[1]),
            xywh[2],
            xywh[3],
            boxstyle="round,pad=0.010,rounding_size=0.009",
            transform=fig.transFigure,
            facecolor=face,
            edgecolor=edge,
            linewidth=1.05,
            zorder=zorder,
        )
    )


def fig_textbox(
    fig: plt.Figure,
    x: float,
    y: float,
    text: str,
    width_chars: int,
    *,
    fontsize: float,
    color: str = INK,
    weight: str = "normal",
    linespacing: float = 1.18,
) -> None:
    fig.text(
        x,
        y,
        textwrap.fill(text, width=width_chars),
        ha="left",
        va="top",
        fontsize=fontsize,
        color=color,
        fontweight=weight,
        linespacing=linespacing,
    )


def median_range(panel: dict[str, object]) -> tuple[float, float]:
    medians = [
        payload["summary"]["median_bdt"]
        for payload in panel["slices"].values()
        if payload["summary"]["n"] > 0
    ]
    return float(min(medians)), float(max(medians))


def minmax_fraction_gt(panel: dict[str, object], key: str) -> tuple[float, float]:
    vals = [
        payload["summary"][key]
        for payload in panel["slices"].values()
        if payload["summary"]["n"] > 0
    ]
    return float(min(vals)), float(max(vals))


def add_title_and_subtitle(fig: plt.Figure) -> None:
    fig.text(
        0.050,
        0.948,
        "Truth-photon BDT scores stay high across isolation slices",
        ha="left",
        va="top",
        fontsize=30.5,
        fontweight="bold",
    )
    fig.add_artist(Rectangle((0.050, 0.897), 0.900, 0.0030, transform=fig.transFigure, color="#d8dee8", lw=0))
    add_box(fig, (0.050, 0.805, 0.900, 0.068), SUBTITLE_BG, SUBTITLE_EDGE)
    fig.text(0.068, 0.852, "What is plotted", ha="left", va="center", fontsize=15.4, fontweight="bold")
    fig.text(
        0.190,
        0.852,
        r"same Slide-21 truth photons, $15 \leq E_T^{cluster}<35$ GeV, split by raw reco $E_T^{iso}$, $\Delta R<0.3$.",
        ha="left",
        va="center",
        fontsize=15.0,
        color=MUTED,
    )
    fig.text(0.068, 0.824, "Question", ha="left", va="center", fontsize=15.4, fontweight="bold")
    fig.text(
        0.190,
        0.824,
        "each curve sums to 100%; this tests whether the high-score band persists across clean, neutral, moderate, and busy cones.",
        ha="left",
        va="center",
        fontsize=15.0,
        color=MUTED,
    )


def plot_panels(fig: plt.Figure, hists: dict[str, dict[str, object]], ymax: float) -> None:
    lefts = [0.070, 0.372, 0.674]
    axes = []
    for left, (_, _, cent_label) in zip(lefts, CENT_BINS):
        ax = fig.add_axes([left, 0.365, 0.258, 0.350])
        axes.append(ax)
        ax.axvspan(0.50, 0.76, color=PHOTON_BAND, alpha=0.50, zorder=0)
        panel = hists[cent_label]
        for slice_key, _, _, color in ISO_SLICES:
            payload = panel["slices"][slice_key]
            edges = np.asarray(payload["bin_edges"])
            values = np.asarray(payload["hist_fraction"])
            ax.stairs(values, edges, color=color, linewidth=2.9, zorder=3)
        lo, hi = median_range(panel)
        gt50_lo, gt50_hi = minmax_fraction_gt(panel, "fraction_bdt_gt_0p5")
        ax.set_title(f"{cent_label} centrality", fontsize=18.0, fontweight="bold", pad=9)
        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(0.0, max(0.105, ymax * 1.22))
        ax.grid(True, color=GRID, lw=0.75)
        ax.tick_params(labelsize=13.3)
        ax.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1.0, decimals=0))
        ax.xaxis.set_major_locator(ticker.MultipleLocator(0.2))
        ax.text(
            0.035,
            0.952,
            f"median BDT {lo:.2f}-{hi:.2f}\nscore>0.5: {gt50_lo:.0%}-{gt50_hi:.0%}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=12.7,
            color="#28313f",
            bbox=dict(boxstyle="round,pad=0.28", facecolor="white", edgecolor="#d8dee8", linewidth=0.85, alpha=0.92),
        )
        if ax is axes[0]:
            ax.set_ylabel("Share of each slice / BDT bin", fontsize=15.0, labelpad=8)
            ax.text(
                0.545,
                max(0.105, ymax * 1.22) * 0.885,
                "Slide-21\nhigh-score band",
                ha="left",
                va="center",
                fontsize=12.5,
                color="#806000",
            )
        else:
            ax.set_yticklabels([])
    fig.text(0.500, 0.305, "BDT photon-likeness score", ha="center", va="center", fontsize=17.0)

    fig.text(0.072, 0.764, "raw isolation slices:", ha="left", va="center", fontsize=15.4, fontweight="bold")
    legend_items = [
        (0.235, "< 0 GeV", ISO_SLICES[0][3]),
        (0.382, "0-2 GeV", ISO_SLICES[1][3]),
        (0.542, "2-5 GeV", ISO_SLICES[2][3]),
        (0.700, "> 5 GeV", ISO_SLICES[3][3]),
    ]
    for x, label, color in legend_items:
        fig.add_artist(Rectangle((x, 0.759), 0.034, 0.006, transform=fig.transFigure, color=color, lw=0))
        fig.text(x + 0.041, 0.764, label, ha="left", va="center", fontsize=14.5)


def add_interpretation(fig: plt.Figure, hists: dict[str, dict[str, object]]) -> None:
    ranges = {label: median_range(hists[label]) for _, _, label in CENT_BINS}
    all_gt50 = [
        payload["summary"]["fraction_bdt_gt_0p5"]
        for panel in hists.values()
        for payload in panel["slices"].values()
        if payload["summary"]["n"] > 0
    ]
    add_box(fig, (0.070, 0.052, 0.860, 0.200), BAND_BG, BAND_EDGE)
    col_x = [0.095, 0.382, 0.666]
    fig.text(col_x[0], 0.226, "How to read it", ha="left", va="top", fontsize=17.0, fontweight="bold")
    fig_textbox(
        fig,
        col_x[0],
        0.190,
        "Each curve is truth photons in one raw isolation slice. Curves all sum to 100%, so line height is shape, not slice size.",
        34,
        fontsize=13.8,
        color="#2f333a",
    )
    fig.text(col_x[1], 0.226, "Numerical read", ha="left", va="top", fontsize=17.0, fontweight="bold")
    fig_textbox(
        fig,
        col_x[1],
        0.190,
        f"Median BDT shifts little: 0-20% {ranges['0-20%'][0]:.2f}-{ranges['0-20%'][1]:.2f}, "
        f"20-50% {ranges['20-50%'][0]:.2f}-{ranges['20-50%'][1]:.2f}, "
        f"50-80% {ranges['50-80%'][0]:.2f}-{ranges['50-80%'][1]:.2f}. "
        f"{min(all_gt50):.0%}-{max(all_gt50):.0%} stay above score 0.5.",
        34,
        fontsize=13.8,
        color="#2f333a",
    )
    fig.text(col_x[2], 0.226, "Interpretation", ha="left", va="top", fontsize=17.0, fontweight="bold")
    fig_textbox(
        fig,
        col_x[2],
        0.190,
        "The Slide-21 truth-photon band is not only negative or near-zero isolation. Genuine photons remain BDT-photon-like even when extra cone energy is present.",
        34,
        fontsize=13.8,
        color="#2f333a",
    )


def write_speaker_script(path: Path, hists: dict[str, dict[str, object]]) -> None:
    ranges = {label: median_range(hists[label]) for _, _, label in CENT_BINS}
    all_gt50 = [
        payload["summary"]["fraction_bdt_gt_0p5"]
        for panel in hists.values()
        for payload in panel["slices"].values()
        if payload["summary"]["n"] > 0
    ]
    path.write_text(
        "\n".join(
            [
                "# Truth-photon BDT score in isolation slices",
                "",
                "This is the direct follow-up to the Slide 21 BDT-isolation density view.",
                "It uses the same local score-cache sample and the same 15 <= cluster E_T < 35 GeV selection, but keeps only truth photons.",
                "",
                "The point of the slide is to test whether the bright truth-photon band on Slide 21 is caused by one narrow raw-isolation region.",
                "To remove the population-size effect, every colored curve is normalized to unit area within its own isolation slice.",
                "",
                "The result is that the BDT-score shape is broadly stable across raw reco E_T^iso slices.",
                f"The median ranges are 0-20%: {ranges['0-20%'][0]:.3f}-{ranges['0-20%'][1]:.3f}, "
                f"20-50%: {ranges['20-50%'][0]:.3f}-{ranges['20-50%'][1]:.3f}, "
                f"50-80%: {ranges['50-80%'][0]:.3f}-{ranges['50-80%'][1]:.3f}.",
                f"Across all shown slices, {min(all_gt50):.1%}-{max(all_gt50):.1%} of truth photons have BDT score greater than 0.5.",
                "",
                "The physics read is that raw isolation is mostly describing the cone environment around the EM cluster.",
                "For genuine photons, additional cone energy does not remove the photon-like BDT response.",
                "So the Slide 21 yellow photon band should be read as a stable truth-photon population, not as a feature produced only by the cleanest isolation slice.",
                "",
            ]
        )
        + "\n"
    )


def make_slide() -> tuple[Path, Path, Path]:
    setup_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    arrays = load_arrays()
    hists, ymax = build_histograms(arrays)

    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)
    add_title_and_subtitle(fig)
    plot_panels(fig, hists, ymax)
    add_interpretation(fig, hists)

    png_path = OUT_DIR / "truth_photon_bdt_isolation_slices_slide.png"
    script_path = OUT_DIR / "truth_photon_bdt_isolation_slices_speaker_script.md"
    manifest_path = OUT_DIR / "truth_photon_bdt_isolation_slices_manifest.json"
    fig.savefig(png_path, dpi=DPI, facecolor=WHITE)
    plt.close(fig)

    write_speaker_script(script_path, hists)
    manifest = {
        "generator": str(Path(__file__).relative_to(REPO)),
        "source_slide21_generator": str(SOURCE_SLIDE21_GENERATOR.relative_to(REPO)),
        "source_slide21_png": str(SOURCE_SLIDE21_PNG.relative_to(REPO)),
        "score_cache_manifest": str(SCORE_CACHE_MANIFEST.relative_to(REPO)),
        "output_png": str(png_path.relative_to(REPO)),
        "speaker_script": str(script_path.relative_to(REPO)),
        "selection": {
            "truth_photons_only": True,
            "score_column": SCORE_COLUMN,
            "isolation_column": ISO_COLUMN,
            "cluster_et_range_gev": list(PT_RANGE),
            "centrality_range_percent": [0.0, 80.0],
            "isolation_cone": "Delta R < 0.3",
        },
        "isolation_slices_gev": [
            {"key": "negative", "label": "E_T^iso < 0"},
            {"key": "near_neutral", "label": "0 <= E_T^iso < 2"},
            {"key": "moderate", "label": "2 <= E_T^iso < 5"},
            {"key": "busy_tail", "label": "E_T^iso >= 5"},
        ],
        "centrality_panels": hists,
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
