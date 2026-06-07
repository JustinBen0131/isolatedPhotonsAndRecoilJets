#!/usr/bin/env python3
"""Build a background-only BDT-slice isolation normalization slide."""

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
SOURCE_SLIDE21_PNG = (
    SOURCE_REPORT
    / "slideReady/default_eiso_vs_bdt_score_20260604/plot_variants_kbird_probability_20260604/02_kbird_log_density_centrality_by_class.png"
)
OUT_DIR = (
    SOURCE_REPORT
    / "slideReady/default_eiso_vs_bdt_score_20260604/the11_background_bdt_slice_normalization_20260606"
)

SCORE_COLUMN = "score_ptFine_cent7"
ISO_COLUMN = "reco_eiso"
PT_RANGE = (15.0, 35.0)
CENT_BINS = [(0.0, 20.0, "0-20%"), (20.0, 50.0, "20-50%"), (50.0, 80.0, "50-80%")]
BDT_SLICES = [
    (0.00, 0.45, "low", "low BDT", "score < 0.45", "#56616f"),
    (0.45, 0.60, "mid", "middle BDT", "0.45-0.60", "#15807a"),
    (0.60, 1.01, "high", "high BDT", r"score $\geq$ 0.60", "#c85a00"),
]

W, H = 2560, 1440
DPI = 200
WHITE = "#ffffff"
INK = "#151515"
MUTED = "#586171"
GRID = "#dde3ea"
PANEL_BG = "#fbfbfb"
SUBTITLE_BG = "#f4f7fb"
SUBTITLE_EDGE = "#cad6e4"
BAND_BG = "#fff7d9"
BAND_EDGE = "#d5c47d"
CARD_BG = "#ffffff"
CARD_EDGE = "#d6dde8"


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
            "axes.linewidth": 1.05,
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


def add_box(
    fig: plt.Figure,
    xywh: tuple[float, float, float, float],
    face: str,
    edge: str,
    *,
    radius: float = 0.010,
    lw: float = 1.0,
    zorder: float = 0.2,
) -> None:
    fig.add_artist(
        FancyBboxPatch(
            (xywh[0], xywh[1]),
            xywh[2],
            xywh[3],
            boxstyle=f"round,pad=0.010,rounding_size={radius}",
            transform=fig.transFigure,
            facecolor=face,
            edgecolor=edge,
            linewidth=lw,
            zorder=zorder,
        )
    )


def wrapped_text(
    fig: plt.Figure,
    x: float,
    y: float,
    text: str,
    width_chars: int,
    *,
    fontsize: float,
    color: str = INK,
    weight: str = "normal",
    linespacing: float = 1.12,
    ha: str = "left",
) -> None:
    fig.text(
        x,
        y,
        textwrap.fill(text, width=width_chars),
        ha=ha,
        va="top",
        fontsize=fontsize,
        color=color,
        fontweight=weight,
        linespacing=linespacing,
    )


def compute_payload(arrays: dict[str, np.ndarray]) -> tuple[dict[str, object], np.ndarray, np.ndarray, float]:
    eiso = arrays["eiso"]
    xlo, xhi = np.nanpercentile(eiso, [0.7, 99.3])
    xlo = min(float(xlo), -9.0)
    xhi = max(float(xhi), 17.0)
    edges = np.linspace(xlo, xhi, 49)
    centers = 0.5 * (edges[:-1] + edges[1:])

    panels: dict[str, object] = {}
    ymax = 0.0
    for cent_lo, cent_hi, cent_label in CENT_BINS:
        panel: dict[str, object] = {"centrality": cent_label, "slices": {}}
        cent_mask = (arrays["centrality"] >= cent_lo) & (arrays["centrality"] < cent_hi)
        for score_lo, score_hi, key, label, score_label, color in BDT_SLICES:
            mask = cent_mask & (arrays["score"] >= score_lo) & (arrays["score"] < score_hi)
            values_all = arrays["eiso"][mask]
            visible = (values_all >= edges[0]) & (values_all <= edges[-1])
            values = values_all[visible]
            counts, _ = np.histogram(values, bins=edges)
            frac = counts.astype("float64")
            if frac.sum() > 0:
                frac /= frac.sum()
            ymax = max(ymax, float(np.max(frac)) if len(frac) else 0.0)
            panel["slices"][key] = {
                "label": label,
                "score_label": score_label,
                "score_range": [score_lo, score_hi],
                "color": color,
                "hist_fraction": frac.tolist(),
                "summary": {
                    "n_total": int(len(values_all)),
                    "n_visible": int(np.sum(visible)),
                    "visible_fraction": float(np.sum(visible) / len(values_all)) if len(values_all) else 0.0,
                    "median_eiso": float(np.nanmedian(values_all)) if len(values_all) else float("nan"),
                    "q25_eiso": float(np.nanquantile(values_all, 0.25)) if len(values_all) else float("nan"),
                    "q75_eiso": float(np.nanquantile(values_all, 0.75)) if len(values_all) else float("nan"),
                    "fraction_eiso_lt_2": float(np.mean(values_all < 2.0)) if len(values_all) else float("nan"),
                    "fraction_eiso_gt_5": float(np.mean(values_all > 5.0)) if len(values_all) else float("nan"),
                },
            }
        panels[cent_label] = panel

    payload: dict[str, object] = {
        "x_range_gev": [float(edges[0]), float(edges[-1])],
        "bin_edges": edges.tolist(),
        "panels": panels,
        "normalization": (
            "For every centrality panel and BDT-score slice, the inclusive/background "
            "isolation histogram is divided by its own total entries."
        ),
    }
    return payload, edges, centers, ymax


def add_title_and_subtitle(fig: plt.Figure) -> None:
    fig.text(
        0.050,
        0.948,
        "Normalized background slices test BDT-isolation structure",
        ha="left",
        va="top",
        fontsize=31.0,
        fontweight="bold",
    )
    fig.add_artist(Rectangle((0.050, 0.902), 0.900, 0.003, transform=fig.transFigure, color="#d8dee8", lw=0))
    add_box(fig, (0.052, 0.812, 0.896, 0.063), SUBTITLE_BG, SUBTITLE_EDGE)
    fig.text(0.070, 0.856, "Sample", ha="left", va="center", fontsize=16.0, fontweight="bold")
    fig.text(
        0.225,
        0.857,
        r"background only; same Slide-21 sample, 15-35 GeV clusters, reco $E_T^{iso}$ with $\Delta R<0.3$.",
        ha="left",
        va="center",
        fontsize=14.6,
        color=MUTED,
    )
    fig.text(0.070, 0.831, "Normalization", ha="left", va="center", fontsize=16.0, fontweight="bold")
    fig.text(
        0.225,
        0.831,
        "each curve is area-normalized inside one centrality bin; N values keep raw occupancy separate.",
        ha="left",
        va="center",
        fontsize=14.6,
        color=MUTED,
    )


def draw_curve_panels(fig: plt.Figure, payload: dict[str, object], centers: np.ndarray, ymax: float) -> None:
    lefts = [0.080, 0.375, 0.670]
    panel_w = 0.250
    panel_h = 0.310
    bottom = 0.445
    y_upper = max(0.110, ymax * 1.22)

    for col, (_, _, cent_label) in enumerate(CENT_BINS):
        ax = fig.add_axes([lefts[col], bottom, panel_w, panel_h])
        panel = payload["panels"][cent_label]
        for idx, (_, _, key, label, score_label, color) in enumerate(BDT_SLICES):
            item = panel["slices"][key]
            hist = np.asarray(item["hist_fraction"])
            ax.fill_between(centers, hist, step="mid", color=color, alpha=0.12 if idx < 2 else 0.16, linewidth=0)
            ax.step(centers, hist, where="mid", color=color, linewidth=2.55, label=f"{label} ({score_label})")
            median = item["summary"]["median_eiso"]
            ax.axvline(median, color=color, linewidth=1.2, alpha=0.68)

        ax.axvline(2.0, color="#202020", lw=1.0, ls=":", alpha=0.65)
        ax.axvline(5.0, color="#202020", lw=1.0, ls="--", alpha=0.65)
        ax.text(
            0.97,
            0.93,
            "2 GeV\n5 GeV",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=10.0,
            color="#343b45",
            linespacing=1.15,
        )
        ax.set_title(f"{cent_label} centrality", fontsize=18.5, pad=10, fontweight="bold")
        ax.set_xlim(payload["x_range_gev"])
        ax.set_ylim(0.0, y_upper)
        ax.grid(True, color=GRID, lw=0.60)
        ax.tick_params(axis="both", labelsize=12.4, pad=2)
        ax.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1.0, decimals=0))
        if col == 0:
            ax.set_ylabel("fraction of background candidates / bin", fontsize=13.0, labelpad=5)
            ax.legend(loc="upper left", bbox_to_anchor=(0.01, 0.995), fontsize=11.0, frameon=True, facecolor="white", edgecolor="#d7dde8")
        else:
            ax.set_yticklabels([])
    fig.text(
        0.500,
        0.399,
        r"reco $E_T^{iso}$, $\Delta R<0.3$ [GeV]",
        ha="center",
        va="center",
        fontsize=14.2,
    )


def draw_summary_cards(fig: plt.Figure, payload: dict[str, object]) -> None:
    lefts = [0.080, 0.375, 0.670]
    card_w = 0.250
    card_h = 0.145
    y = 0.214
    for col, (_, _, cent_label) in enumerate(CENT_BINS):
        add_box(fig, (lefts[col], y, card_w, card_h), CARD_BG, CARD_EDGE, radius=0.008, lw=0.9)
        fig.text(lefts[col] + 0.017, y + card_h - 0.024, "raw N and shape read", ha="left", va="top", fontsize=13.4, fontweight="bold")
        panel = payload["panels"][cent_label]
        row_y = [y + 0.090, y + 0.057, y + 0.024]
        for idx, (_, _, key, label, _, color) in enumerate(BDT_SLICES):
            s = panel["slices"][key]["summary"]
            fig.add_artist(Rectangle((lefts[col] + 0.018, row_y[idx] - 0.006), 0.018, 0.007, transform=fig.transFigure, color=color, lw=0))
            fig.text(lefts[col] + 0.042, row_y[idx], label, ha="left", va="center", fontsize=11.1, color=color, fontweight="bold")
            fig.text(
                lefts[col] + 0.122,
                row_y[idx],
                f"N={s['n_total']:,}   med={s['median_eiso']:.1f}   >5={s['fraction_eiso_gt_5']:.0%}",
                ha="left",
                va="center",
                fontsize=10.9,
                color="#2f3540",
            )


def draw_interpretation(fig: plt.Figure, payload: dict[str, object]) -> None:
    add_box(fig, (0.070, 0.038, 0.860, 0.150), BAND_BG, BAND_EDGE, radius=0.009, lw=1.0)
    cols = [0.095, 0.385, 0.662]
    fig.text(cols[0], 0.160, "What changed from Slide 21", ha="left", va="top", fontsize=15.4, fontweight="bold")
    wrapped_text(
        fig,
        cols[0],
        0.129,
        "The hot spots are no longer raw counts. Area-normalized curves prevent the 0-20% statistics from visually dominating the peripheral panels.",
        38,
        fontsize=11.8,
        color="#2f333a",
    )
    fig.text(cols[1], 0.160, "What the normalized test says", ha="left", va="top", fontsize=15.4, fontweight="bold")
    wrapped_text(
        fig,
        cols[1],
        0.129,
        "Across low, middle, and high BDT slices, the background shape stays near 4 GeV with a large >5 GeV busy tail.",
        36,
        fontsize=11.8,
        color="#2f333a",
    )
    fig.text(cols[2], 0.160, "ABCD implication", ha="left", va="top", fontsize=15.4, fontweight="bold")
    wrapped_text(
        fig,
        cols[2],
        0.129,
        "This is the background-only shape check: isolation remains useful after BDT selection; final validity comes from the ABCD closure map.",
        41,
        fontsize=11.8,
        color="#2f333a",
    )


def write_script(path: Path, payload: dict[str, object]) -> None:
    high_bits = []
    for _, _, cent_label in CENT_BINS:
        s = payload["panels"][cent_label]["slices"]["high"]["summary"]
        high_bits.append(f"{cent_label}: median {s['median_eiso']:.1f} GeV, >5 GeV tail {s['fraction_eiso_gt_5']:.0%}")
    path.write_text(
        "\n".join(
            [
                "# THE-11 Slide Script - Background Isolation Shapes In BDT Slices",
                "",
                "Here I am taking the inclusive-background row from the previous BDT-isolation density slide and removing the most dangerous visual ambiguity.",
                "Instead of comparing raw occupancies, each curve is normalized to unit area inside one centrality bin and one BDT-score slice.",
                "",
                "The three columns are the same centrality bins as before.",
                "Within a column, the gray, teal, and orange curves show the background isolation distribution in low, middle, and high BDT-score regions.",
                "The printed N values are the original counts, so we can see the statistics without letting the central 0-20% sample visually dominate the plot.",
                "",
                "The important qualitative result is that the background isolation shape does not collapse to a clean photon-like isolation distribution when the BDT score becomes high.",
                "For the high-BDT slice, the inclusive-background medians and busy tails are "
                + "; ".join(high_bits)
                + ".",
                "So the background candidates that score photon-like still carry substantial positive cone energy.",
                "",
                "That is the ABCD-relevant reading.",
                "The BDT and isolation axes are related, but isolation remains a separate handle on nearby hadronic activity.",
                "This normalized shape check helps interpret the hot spots in the density plot, while the actual ABCD validity statement should still come from the closure-ratio map.",
                "",
            ]
        )
        + "\n"
    )


def make_slide() -> tuple[Path, Path, Path]:
    setup_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    arrays = load_background_arrays()
    payload, _edges, centers, ymax = compute_payload(arrays)

    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)
    add_title_and_subtitle(fig)
    draw_curve_panels(fig, payload, centers, ymax)
    draw_summary_cards(fig, payload)
    draw_interpretation(fig, payload)

    png_path = OUT_DIR / "background_bdt_slice_normalization_slide.png"
    script_path = OUT_DIR / "background_bdt_slice_normalization_speaker_script.md"
    manifest_path = OUT_DIR / "background_bdt_slice_normalization_manifest.json"
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
        "bdt_slices": [
            {"key": key, "score_range": [lo, hi], "label": label, "score_label": score_label}
            for lo, hi, key, label, score_label, _ in BDT_SLICES
        ],
        "interpretation_limit": (
            "This slide diagnoses normalized background isolation-shape stability across BDT score slices. "
            "It does not replace the explicit ABCD closure-ratio test."
        ),
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
