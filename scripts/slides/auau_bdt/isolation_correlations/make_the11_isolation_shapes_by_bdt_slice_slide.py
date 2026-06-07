#!/usr/bin/env python3
"""Build the THE-11 flipped BDT-slice isolation-shape slide."""

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
    / "slideReady/default_eiso_vs_bdt_score_20260604/the11_isolation_shapes_by_bdt_slice_20260605"
)

SCORE_COLUMN = "score_ptFine_cent7"
ISO_COLUMN = "reco_eiso"
PT_RANGE = (15.0, 35.0)
CENT_BINS = [(0.0, 20.0, "0-20%"), (20.0, 50.0, "20-50%"), (50.0, 80.0, "50-80%")]
BDT_SLICES = [
    (0.00, 0.45, "background-like", "low BDT\nscore < 0.45"),
    (0.45, 0.60, "transition", "middle BDT\n0.45-0.60"),
    (0.60, 1.01, "photon-like", "high BDT\n" + r"score $\geq 0.60$"),
]

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
SIGNAL = "#b0182d"
BACKGROUND = "#1f5aa6"
HIGH_ROW = "#fff1a8"


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
            "axes.linewidth": 1.0,
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


def add_box(fig: plt.Figure, xywh: tuple[float, float, float, float], face: str, edge: str) -> None:
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
            zorder=0.2,
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
    linespacing: float = 1.16,
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


def summary(values: np.ndarray, visible_mask: np.ndarray, total_before_visible: int) -> dict[str, float | int]:
    if len(values) == 0:
        return {
            "n": 0,
            "median_eiso": float("nan"),
            "q25_eiso": float("nan"),
            "q75_eiso": float("nan"),
            "fraction_eiso_lt_2": float("nan"),
            "fraction_eiso_gt_5": float("nan"),
            "visible_fraction": 0.0,
        }
    return {
        "n": int(len(values)),
        "median_eiso": float(np.nanmedian(values)),
        "q25_eiso": float(np.nanquantile(values, 0.25)),
        "q75_eiso": float(np.nanquantile(values, 0.75)),
        "fraction_eiso_lt_2": float(np.mean(values < 2.0)),
        "fraction_eiso_gt_5": float(np.mean(values > 5.0)),
        "visible_fraction": float(np.sum(visible_mask) / total_before_visible) if total_before_visible else 0.0,
    }


def build_payload(arrays: dict[str, np.ndarray]) -> tuple[dict[str, object], np.ndarray, float]:
    eiso = arrays["eiso"]
    xlo, xhi = np.nanpercentile(eiso, [0.5, 99.5])
    pad = 0.05 * (xhi - xlo)
    xlo, xhi = float(xlo - pad), float(xhi + pad)
    edges = np.linspace(xlo, xhi, 55)

    signal = arrays["is_signal"].astype(bool)
    payload: dict[str, object] = {
        "x_range_gev": [xlo, xhi],
        "bin_edges": edges.tolist(),
        "panels": {},
    }
    ymax = 0.0
    for cent_lo, cent_hi, cent_label in CENT_BINS:
        for score_lo, score_hi, slice_key, slice_label in BDT_SLICES:
            panel_key = f"{cent_label}_{slice_key}".replace("%", "pct")
            mask_base = (
                (arrays["centrality"] >= cent_lo)
                & (arrays["centrality"] < cent_hi)
                & (arrays["score"] >= score_lo)
                & (arrays["score"] < score_hi)
            )
            panel: dict[str, object] = {
                "centrality": cent_label,
                "bdt_slice_key": slice_key,
                "bdt_slice_label": slice_label,
                "score_range": [score_lo, score_hi],
                "classes": {},
            }
            for class_key, class_mask in [("truth_photons", signal), ("inclusive_jets", ~signal)]:
                mask = mask_base & class_mask
                class_values_all = arrays["eiso"][mask]
                visible = (class_values_all >= xlo) & (class_values_all <= xhi)
                class_values = class_values_all[visible]
                counts, _ = np.histogram(class_values, bins=edges)
                frac = counts.astype("float64")
                if frac.sum() > 0:
                    frac /= frac.sum()
                ymax = max(ymax, float(np.nanmax(frac)) if len(frac) else 0.0)
                panel["classes"][class_key] = {
                    "hist_fraction": frac.tolist(),
                    "summary": summary(class_values_all, visible, len(class_values_all)),
                }
            payload["panels"][panel_key] = panel
    return payload, edges, ymax


def add_title_and_subtitle(fig: plt.Figure) -> None:
    fig.text(
        0.050,
        0.950,
        "BDT slices show isolation remains a separate background handle",
        ha="left",
        va="top",
        fontsize=30.0,
        fontweight="bold",
    )
    fig.add_artist(Rectangle((0.050, 0.899), 0.900, 0.0030, transform=fig.transFigure, color="#d8dee8", lw=0))
    add_box(fig, (0.050, 0.809, 0.900, 0.066), SUBTITLE_BG, SUBTITLE_EDGE)
    fig.text(0.068, 0.855, "What is plotted", ha="left", va="center", fontsize=15.4, fontweight="bold")
    fig.text(
        0.190,
        0.855,
        r"same Slide-21 sample; columns are centrality, rows are BDT-score phase-space slices.",
        ha="left",
        va="center",
        fontsize=15.0,
        color=MUTED,
    )
    fig.text(0.068, 0.828, "Curve meaning", ha="left", va="center", fontsize=15.4, fontweight="bold")
    fig.text(
        0.190,
        0.828,
        r"red truth-photon and blue inclusive-jet isolation shapes are each normalized to 100% within that panel.",
        ha="left",
        va="center",
        fontsize=15.0,
        color=MUTED,
    )


def panel_key(cent_label: str, slice_key: str) -> str:
    return f"{cent_label}_{slice_key}".replace("%", "pct")


def plot_panels(fig: plt.Figure, payload: dict[str, object], edges: np.ndarray, ymax: float) -> None:
    lefts = [0.140, 0.415, 0.690]
    bottoms = [0.590, 0.413, 0.236]
    row_height = 0.132
    panel_width = 0.220
    y_upper = max(0.125, ymax * 1.16)
    centers = 0.5 * (edges[:-1] + edges[1:])

    fig.text(0.120, 0.769, "truth photons", ha="left", va="center", fontsize=15.2, fontweight="bold", color=SIGNAL)
    fig.add_artist(Rectangle((0.068, 0.766), 0.043, 0.006, transform=fig.transFigure, color=SIGNAL, lw=0))
    fig.text(0.300, 0.769, "inclusive jets", ha="left", va="center", fontsize=15.2, fontweight="bold", color=BACKGROUND)
    fig.add_artist(Rectangle((0.248, 0.766), 0.043, 0.006, transform=fig.transFigure, color=BACKGROUND, lw=0))
    fig.text(0.486, 0.769, "all curves normalized inside panel", ha="left", va="center", fontsize=14.3, color=MUTED)

    for col, (_, _, cent_label) in enumerate(CENT_BINS):
        fig.text(
            lefts[col] + panel_width / 2,
            0.724,
            f"{cent_label} centrality",
            ha="center",
            va="bottom",
            fontsize=17.5,
            fontweight="bold",
        )

    panels = payload["panels"]
    for row, (_, _, slice_key, slice_label) in enumerate(BDT_SLICES):
        fig.text(
            0.075,
            bottoms[row] + row_height / 2,
            slice_label,
            ha="center",
            va="center",
            rotation=0,
            fontsize=13.5,
            fontweight="bold" if row == 2 else "normal",
            color="#6b5200" if row == 2 else INK,
            linespacing=1.05,
        )
        for col, (_, _, cent_label) in enumerate(CENT_BINS):
            ax = fig.add_axes([lefts[col], bottoms[row], panel_width, row_height])
            if row == 2:
                ax.set_facecolor("#fffdf1")
                ax.axhspan(0.0, y_upper, color=HIGH_ROW, alpha=0.23, zorder=0)
            p = panels[panel_key(cent_label, slice_key)]
            truth = np.asarray(p["classes"]["truth_photons"]["hist_fraction"])
            jets = np.asarray(p["classes"]["inclusive_jets"]["hist_fraction"])
            ax.fill_between(centers, truth, step="mid", color=SIGNAL, alpha=0.22, linewidth=0)
            ax.step(centers, truth, where="mid", color=SIGNAL, linewidth=2.2)
            ax.fill_between(centers, jets, step="mid", color=BACKGROUND, alpha=0.19, linewidth=0)
            ax.step(centers, jets, where="mid", color=BACKGROUND, linewidth=2.2)

            t_median = p["classes"]["truth_photons"]["summary"]["median_eiso"]
            b_median = p["classes"]["inclusive_jets"]["summary"]["median_eiso"]
            ax.axvline(t_median, color=SIGNAL, linewidth=1.2, alpha=0.82)
            ax.axvline(b_median, color=BACKGROUND, linewidth=1.2, alpha=0.82)
            if row == 2:
                ax.text(
                    0.03,
                    0.88,
                    rf"med: {t_median:.1f} vs {b_median:.1f} GeV",
                    transform=ax.transAxes,
                    ha="left",
                    va="top",
                    fontsize=10.1,
                    color="#293241",
                    bbox=dict(boxstyle="round,pad=0.22", facecolor="white", edgecolor="#d8dee8", linewidth=0.75, alpha=0.92),
                )
            ax.set_xlim(edges[0], edges[-1])
            ax.set_ylim(0.0, y_upper)
            ax.grid(True, color=GRID, lw=0.55)
            ax.tick_params(labelsize=10.5, pad=2)
            ax.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1.0, decimals=0))
            if col != 0:
                ax.set_yticklabels([])
            else:
                ax.set_ylabel("")
            if row != 2:
                ax.set_xticklabels([])
            else:
                ax.set_xlabel(r"reco $E_T^{iso}$, $\Delta R<0.3$ [GeV]", fontsize=11.2, labelpad=3)


def high_bdt_numbers(payload: dict[str, object]) -> tuple[list[str], list[str]]:
    med_lines: list[str] = []
    tail_lines: list[str] = []
    panels = payload["panels"]
    for _, _, cent_label in CENT_BINS:
        p = panels[panel_key(cent_label, "photon-like")]
        truth = p["classes"]["truth_photons"]["summary"]
        jets = p["classes"]["inclusive_jets"]["summary"]
        med_lines.append(f"{cent_label}: {truth['median_eiso']:.1f} vs {jets['median_eiso']:.1f} GeV")
        tail_lines.append(
            f"{cent_label}: {truth['fraction_eiso_gt_5']:.0%} vs {jets['fraction_eiso_gt_5']:.0%}"
        )
    return med_lines, tail_lines


def add_interpretation(fig: plt.Figure, payload: dict[str, object]) -> None:
    med_lines, tail_lines = high_bdt_numbers(payload)
    add_box(fig, (0.070, 0.040, 0.860, 0.154), BAND_BG, BAND_EDGE)
    col_x = [0.095, 0.382, 0.666]
    fig.text(col_x[0], 0.168, "How to read it", ha="left", va="top", fontsize=15.8, fontweight="bold")
    fig_textbox(
        fig,
        col_x[0],
        0.137,
        "Rows move from low to high BDT score. Within every panel, compare red and blue isolation shapes, not class yields.",
        34,
        fontsize=12.2,
        color="#2f333a",
    )
    fig.text(col_x[1], 0.168, "High-BDT numerical read", ha="left", va="top", fontsize=15.8, fontweight="bold")
    fig_textbox(
        fig,
        col_x[1],
        0.137,
        "Median isolation, truth vs jets: "
        + "; ".join(med_lines)
        + ". Jet busy tails above 5 GeV remain much larger.",
        36,
        fontsize=12.2,
        color="#2f333a",
    )
    fig.text(col_x[2], 0.168, "Interpretation", ha="left", va="top", fontsize=15.8, fontweight="bold")
    fig_textbox(
        fig,
        col_x[2],
        0.137,
        "Even inside photon-like BDT space, inclusive jets retain extra cone energy. BDT and isolation are correlated, but not redundant.",
        34,
        fontsize=12.2,
        color="#2f333a",
    )


def write_speaker_script(path: Path, payload: dict[str, object]) -> None:
    med_lines, tail_lines = high_bdt_numbers(payload)
    path.write_text(
        "\n".join(
            [
                "# THE-11 Slide Script - Isolation Shapes In BDT Slices",
                "",
                "Now I am flipping the previous diagnostic around.",
                "Instead of asking how the BDT score changes in isolation slices, I am taking slices of BDT score and looking at the raw isolation distribution inside each slice.",
                "",
                "The columns are the same three centrality bins as the previous BDT-isolation view.",
                "The rows move from background-like BDT scores to the photon-like BDT region.",
                "In every panel, the red curve is truth photons and the blue curve is inclusive jets, and each curve is normalized within that panel.",
                "So the comparison is the shape of the isolation distribution, not how many objects are in each class.",
                "",
                "The important thing is the bottom row.",
                "Even after selecting high, photon-like BDT scores, the inclusive-jet isolation distribution remains shifted to positive cone energy.",
                "The high-BDT median isolation values for truth photons versus inclusive jets are "
                + "; ".join(med_lines)
                + ".",
                "The same pattern shows up in the busy tail above 5 GeV: "
                + "; ".join(tail_lines)
                + ".",
                "",
                "That is the clean physics interpretation.",
                "The BDT is learning photon-like shower structure, while isolation still carries independent information about nearby hadronic activity in the cone.",
                "They are related, but they are not the same handle, which is why the isolation cut remains useful even after the BDT selection.",
                "",
            ]
        )
        + "\n"
    )


def make_slide() -> tuple[Path, Path, Path]:
    setup_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    arrays = load_arrays()
    payload, edges, ymax = build_payload(arrays)

    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)
    add_title_and_subtitle(fig)
    plot_panels(fig, payload, edges, ymax)
    add_interpretation(fig, payload)

    png_path = OUT_DIR / "isolation_shapes_by_bdt_slice_slide.png"
    script_path = OUT_DIR / "isolation_shapes_by_bdt_slice_speaker_script.md"
    manifest_path = OUT_DIR / "isolation_shapes_by_bdt_slice_manifest.json"
    fig.savefig(png_path, dpi=DPI, facecolor=WHITE)
    plt.close(fig)

    write_speaker_script(script_path, payload)
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
            "classes": {
                "truth_photons": "is_signal == 1",
                "inclusive_jets": "is_signal == 0",
            },
        },
        "normalization": "Each class curve is normalized to unit area within each centrality and BDT-score slice.",
        "bdt_slices": [
            {"key": key, "score_range": [lo, hi], "label": label}
            for lo, hi, key, label in BDT_SLICES
        ],
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
