#!/usr/bin/env python3
"""Build a THE-11 slide showing tight-BDT isolation-shape closure."""

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
OUT_DIR = SOURCE_REPORT / "slideReady/default_eiso_vs_bdt_score_20260604/the11_bdt_tight_id_isolation_closure_20260605"

SCORE_COLUMN = "score_ptFine_cent7"
ISO_COLUMN = "reco_eiso"
PT_RANGE = (15.0, 35.0)
CENT_BINS = [(0.0, 20.0, "0-20%"), (20.0, 50.0, "20-50%"), (50.0, 80.0, "50-80%")]
TARGET_SIGNAL_EFFICIENCY = 0.80

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
TRUTH = "#b0182d"
POOL_BEFORE = "#6b7280"
POOL_TIGHT = "#1b9e77"
POOL_REJECTED = "#5b6472"
TAIL_SHADE = "#fff1a8"


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


def fig_textbox(fig: plt.Figure, x: float, y: float, text: str, width_chars: int, *, fontsize: float) -> None:
    fig.text(
        x,
        y,
        textwrap.fill(text, width=width_chars),
        ha="left",
        va="top",
        fontsize=fontsize,
        color="#2f333a",
        linespacing=1.14,
    )


def hist_fraction(values: np.ndarray, edges: np.ndarray) -> np.ndarray:
    counts, _ = np.histogram(values, bins=edges)
    frac = counts.astype("float64")
    if frac.sum() > 0:
        frac /= frac.sum()
    return frac


def distribution_summary(values: np.ndarray, visible: np.ndarray | None = None) -> dict[str, float | int]:
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
        "visible_fraction": float(np.mean(visible)) if visible is not None and len(visible) else 1.0,
    }


def total_variation(a: np.ndarray, b: np.ndarray) -> float:
    return float(0.5 * np.sum(np.abs(a - b)))


def build_payload(arrays: dict[str, np.ndarray]) -> tuple[dict[str, object], np.ndarray, float]:
    eiso = arrays["eiso"]
    xlo, xhi = np.nanpercentile(eiso, [0.5, 99.5])
    pad = 0.05 * (xhi - xlo)
    edges = np.linspace(float(xlo - pad), float(xhi + pad), 64)
    centers = 0.5 * (edges[:-1] + edges[1:])

    signal = arrays["is_signal"].astype(bool)
    panels: dict[str, object] = {}
    ymax = 0.0
    for cent_lo, cent_hi, cent_label in CENT_BINS:
        cent_mask = (arrays["centrality"] >= cent_lo) & (arrays["centrality"] < cent_hi)
        truth_mask = cent_mask & signal
        if not np.any(truth_mask):
            raise ValueError(f"No truth photons in centrality bin {cent_label}")
        threshold = float(np.nanquantile(arrays["score"][truth_mask], 1.0 - TARGET_SIGNAL_EFFICIENCY))
        tight_mask = cent_mask & (arrays["score"] > threshold)
        rejected_mask = cent_mask & ~tight_mask
        truth_tight_mask = truth_mask & tight_mask
        background_mask = cent_mask & ~signal
        background_tight_mask = background_mask & tight_mask

        samples = {
            "truth_tight": arrays["eiso"][truth_tight_mask],
            "pool_before": arrays["eiso"][cent_mask],
            "pool_tight": arrays["eiso"][tight_mask],
            "pool_rejected": arrays["eiso"][rejected_mask],
            "background_before": arrays["eiso"][background_mask],
            "background_tight": arrays["eiso"][background_tight_mask],
            "truth_before": arrays["eiso"][truth_mask],
        }
        visible_samples: dict[str, np.ndarray] = {}
        hists: dict[str, list[float]] = {}
        for key, values in samples.items():
            visible = (values >= edges[0]) & (values <= edges[-1])
            visible_samples[key] = visible
            hist = hist_fraction(values[visible], edges)
            hists[key] = hist.tolist()
            ymax = max(ymax, float(np.max(hist)) if hist.size else 0.0)

        truth_hist = np.asarray(hists["truth_tight"])
        before_hist = np.asarray(hists["pool_before"])
        tight_hist = np.asarray(hists["pool_tight"])
        rejected_hist = np.asarray(hists["pool_rejected"])
        panel = {
            "centrality": cent_label,
            "score_threshold": threshold,
            "truth_fraction_before": float(np.mean(signal[cent_mask])),
            "truth_fraction_after_tight": float(np.mean(signal[tight_mask])) if np.any(tight_mask) else float("nan"),
            "truth_efficiency": float(np.sum(truth_tight_mask) / np.sum(truth_mask)),
            "background_pass_fraction": float(np.sum(background_tight_mask) / np.sum(background_mask)) if np.any(background_mask) else float("nan"),
            "shape_distance_before_to_truth_tight": total_variation(before_hist, truth_hist),
            "shape_distance_tight_to_truth_tight": total_variation(tight_hist, truth_hist),
            "shape_distance_rejected_to_truth_tight": total_variation(rejected_hist, truth_hist),
            "classes": {
                key: {
                    "hist_fraction": hists[key],
                    "summary": distribution_summary(samples[key], visible_samples[key]),
                }
                for key in samples
            },
        }
        panels[cent_label] = panel
    return {"x_range_gev": [float(edges[0]), float(edges[-1])], "bin_edges": edges.tolist(), "panels": panels}, centers, ymax


def add_title_and_subtitle(fig: plt.Figure) -> None:
    fig.text(
        0.050,
        0.950,
        "Tight BDT ID makes reco-candidate isolation photon-like",
        ha="left",
        va="top",
        fontsize=30.0,
        fontweight="bold",
    )
    fig.add_artist(Rectangle((0.050, 0.900), 0.900, 0.0030, transform=fig.transFigure, color="#d8dee8", lw=0))
    add_box(fig, (0.050, 0.812, 0.900, 0.064), SUBTITLE_BG, SUBTITLE_EDGE)
    fig.text(0.068, 0.858, "What changed", ha="left", va="center", fontsize=15.3, fontweight="bold")
    fig.text(
        0.190,
        0.858,
        r"compare reco candidates that pass a same-cache WP80 BDT cut with candidates removed by that same cut.",
        ha="left",
        va="center",
        fontsize=14.7,
        color=MUTED,
    )
    fig.text(0.068, 0.831, "Target shape", ha="left", va="center", fontsize=15.3, fontweight="bold")
    fig.text(
        0.190,
        0.831,
        r"red is truth photons after the cut; green tracking red means the selected reco pool keeps photon-like cone activity.",
        ha="left",
        va="center",
        fontsize=14.7,
        color=MUTED,
    )


def add_legend(fig: plt.Figure) -> None:
    y = 0.770
    items = [
        (0.090, TRUTH, "truth photons after tight ID", "solid"),
        (0.352, POOL_REJECTED, "reco candidates rejected by BDT", "dashed"),
        (0.660, POOL_TIGHT, "selected reco candidates after tight ID", "solid"),
    ]
    for x, color, label, style in items:
        if style == "dashed":
            fig.add_artist(Rectangle((x, y - 0.003), 0.044, 0.006, transform=fig.transFigure, color=color, alpha=0.40, lw=0))
            fig.text(x + 0.052, y, label, ha="left", va="center", fontsize=14.2, color=color)
        else:
            fig.add_artist(Rectangle((x, y - 0.003), 0.044, 0.006, transform=fig.transFigure, color=color, lw=0))
            fig.text(x + 0.052, y, label, ha="left", va="center", fontsize=14.2, color=color, fontweight="bold")


def plot_panels(fig: plt.Figure, payload: dict[str, object], centers: np.ndarray, ymax: float) -> None:
    edges = np.asarray(payload["bin_edges"])
    lefts = [0.070, 0.372, 0.674]
    panel_width = 0.257
    y_upper = max(0.112, ymax * 1.18)
    for left, (_, _, cent_label) in zip(lefts, CENT_BINS):
        panel = payload["panels"][cent_label]
        ax = fig.add_axes([left, 0.345, panel_width, 0.365])
        ax.axvspan(5.0, edges[-1], color=TAIL_SHADE, alpha=0.44, zorder=0)
        truth = np.asarray(panel["classes"]["truth_tight"]["hist_fraction"])
        before = np.asarray(panel["classes"]["pool_before"]["hist_fraction"])
        tight = np.asarray(panel["classes"]["pool_tight"]["hist_fraction"])
        rejected = np.asarray(panel["classes"]["pool_rejected"]["hist_fraction"])
        ax.fill_between(centers, rejected, step="mid", color=POOL_REJECTED, alpha=0.10, linewidth=0)
        ax.step(centers, rejected, where="mid", color=POOL_REJECTED, linewidth=2.0, alpha=0.76, linestyle=(0, (4, 3)))
        ax.fill_between(centers, tight, step="mid", color=POOL_TIGHT, alpha=0.22, linewidth=0)
        ax.step(centers, tight, where="mid", color=POOL_TIGHT, linewidth=2.7)
        ax.step(centers, truth, where="mid", color=TRUTH, linewidth=2.7)
        ax.step(centers, before, where="mid", color=POOL_BEFORE, linewidth=1.4, alpha=0.35, linestyle=(0, (1, 2)))

        truth_median = panel["classes"]["truth_tight"]["summary"]["median_eiso"]
        tight_median = panel["classes"]["pool_tight"]["summary"]["median_eiso"]
        ax.axvline(truth_median, color=TRUTH, linewidth=1.15, alpha=0.82)
        ax.axvline(tight_median, color=POOL_TIGHT, linewidth=1.15, alpha=0.82)
        ax.set_title(f"{cent_label} centrality", fontsize=17.8, fontweight="bold", pad=9)
        ax.set_xlim(edges[0], edges[-1])
        ax.set_ylim(0.0, y_upper)
        ax.grid(True, color=GRID, lw=0.65)
        ax.tick_params(labelsize=12.3)
        ax.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1.0, decimals=0))
        ax.set_xlabel(r"reco $E_T^{iso}$, $\Delta R<0.3$ [GeV]", fontsize=13.4, labelpad=5)
        if left == lefts[0]:
            ax.set_ylabel("share per isolation bin", fontsize=13.2)
        else:
            ax.set_yticklabels([])
        ax.text(
            0.030,
            0.960,
            f"WP80 score cut: {panel['score_threshold']:.3f}\n"
            f"val truth frac.: {panel['truth_fraction_before']:.0%} -> {panel['truth_fraction_after_tight']:.0%}\n"
            f"distance to red: pass {panel['shape_distance_tight_to_truth_tight']:.3f}, reject {panel['shape_distance_rejected_to_truth_tight']:.3f}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=11.4,
            color="#28313f",
            bbox=dict(boxstyle="round,pad=0.28", facecolor="white", edgecolor="#d8dee8", linewidth=0.85, alpha=0.94),
        )
        tail_truth = panel["classes"]["truth_tight"]["summary"]["fraction_eiso_gt_5"]
        tail_before = panel["classes"]["pool_before"]["summary"]["fraction_eiso_gt_5"]
        tail_tight = panel["classes"]["pool_tight"]["summary"]["fraction_eiso_gt_5"]
        tail_rejected = panel["classes"]["pool_rejected"]["summary"]["fraction_eiso_gt_5"]
        ax.text(
            0.975,
            0.065,
            f"tail >5 GeV\npass {tail_tight:.0%}  truth {tail_truth:.0%}\nreject {tail_rejected:.0%}",
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=11.2,
            color="#6b5200",
            bbox=dict(boxstyle="round,pad=0.26", facecolor="#fffdf1", edgecolor="#e8d78b", linewidth=0.75, alpha=0.92),
        )


def add_interpretation(fig: plt.Figure, payload: dict[str, object]) -> None:
    panels = payload["panels"]
    purities = [(p["truth_fraction_before"], p["truth_fraction_after_tight"]) for p in panels.values()]
    distances = [(p["shape_distance_before_to_truth_tight"], p["shape_distance_tight_to_truth_tight"]) for p in panels.values()]
    reject_distances = [p["shape_distance_rejected_to_truth_tight"] for p in panels.values()]
    bkg_pass = [p["background_pass_fraction"] for p in panels.values()]
    add_box(fig, (0.070, 0.040, 0.860, 0.167), BAND_BG, BAND_EDGE)
    col_x = [0.095, 0.382, 0.666]
    fig.text(col_x[0], 0.181, "Selection effect", ha="left", va="top", fontsize=15.6, fontweight="bold")
    fig_textbox(
        fig,
        col_x[0],
        0.150,
        "Ask what the selected reco pool looks like. The rejected pool is the contrast, so this tests the BDT selection rather than repeating a fixed label split.",
        34,
        fontsize=11.6,
    )
    fig.text(col_x[1], 0.181, "Numerical read", ha="left", va="top", fontsize=15.6, fontweight="bold")
    fig_textbox(
        fig,
        col_x[1],
        0.150,
        f"Selected-pool distance to the red target is {min(b for _, b in distances):.3f}-{max(b for _, b in distances):.3f}; "
        f"rejected candidates sit farther away at {min(reject_distances):.3f}-{max(reject_distances):.3f}. "
        f"Validation truth fraction rises {min(a for a, _ in purities):.0%}-{max(a for a, _ in purities):.0%} to {min(b for _, b in purities):.0%}-{max(b for _, b in purities):.0%}.",
        34,
        fontsize=11.6,
    )
    fig.text(col_x[2], 0.181, "Interpretation", ha="left", va="top", fontsize=15.6, fontweight="bold")
    fig_textbox(
        fig,
        col_x[2],
        0.150,
        f"The BDT selects photon-like shower structure. Isolation remains the residual busy-cone handle, not a duplicate of the score; background-label pass: {min(bkg_pass):.0%}-{max(bkg_pass):.0%}.",
        34,
        fontsize=11.6,
    )


def write_speaker_script(path: Path, payload: dict[str, object]) -> None:
    panels = payload["panels"]
    purity_lines = [
        f"{cent}: {p['truth_fraction_before']:.1%} to {p['truth_fraction_after_tight']:.1%}"
        for cent, p in panels.items()
    ]
    tv_lines = [
        f"{cent}: {p['shape_distance_before_to_truth_tight']:.3f} to {p['shape_distance_tight_to_truth_tight']:.3f}"
        for cent, p in panels.items()
    ]
    tail_lines = [
        f"{cent}: {p['classes']['pool_before']['summary']['fraction_eiso_gt_5']:.1%} to "
        f"{p['classes']['pool_tight']['summary']['fraction_eiso_gt_5']:.1%}, compared with "
        f"{p['classes']['truth_tight']['summary']['fraction_eiso_gt_5']:.1%} for tight truth photons and "
        f"{p['classes']['pool_rejected']['summary']['fraction_eiso_gt_5']:.1%} for rejected candidates"
        for cent, p in panels.items()
    ]
    path.write_text(
        "\n".join(
            [
                "# THE-11 Slide Script - Tight BDT Isolation-Shape Closure",
                "",
                "This version is the more useful way to ask the BDT-isolation question.",
                "Instead of comparing fixed truth and background labels, I am asking what happens to the reco-candidate pool after applying the tight BDT selection.",
                "",
                "The red curve is the target: truth photons that pass the same tight BDT cut.",
                "The gray dashed curve is the reco-candidate pool rejected by the BDT cut; the faint gray dotted line is the full pool before the cut.",
                "The green curve is the reco-candidate pool after the BDT cut.",
                "All curves are normalized within each centrality bin, so the visual comparison is shape.",
                "",
                "The tight BDT cut here is derived inside the same Slide-21 score cache as an 80 percent truth-photon efficiency cut in each broad centrality bin.",
                "After applying that cut, the candidate-pool isolation shape moves onto the truth-photon target shape.",
                "The truth fraction in this validation mixture changes from "
                + "; ".join(purity_lines)
                + ".",
                "The total-variation shape distance to the tight-truth target changes from "
                + "; ".join(tv_lines)
                + ".",
                "The busy isolation tail above 5 GeV changes from "
                + "; ".join(tail_lines)
                + ".",
                "",
                "The interpretation is that the BDT is doing what we want from a photon-ID selector.",
                "It does not make every background-like reco cluster physically signal-like.",
                "It enriches the selected candidate pool in photon-like shower structure, and the remaining pool has an isolation distribution close to the truth-photon target.",
                "That is why the BDT and isolation should be treated as related but not redundant handles.",
                "",
            ]
        )
        + "\n"
    )


def make_slide() -> tuple[Path, Path, Path]:
    setup_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    arrays = load_arrays()
    payload, centers, ymax = build_payload(arrays)

    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)
    add_title_and_subtitle(fig)
    add_legend(fig)
    plot_panels(fig, payload, centers, ymax)
    add_interpretation(fig, payload)

    png_path = OUT_DIR / "bdt_tight_id_isolation_closure_slide.png"
    script_path = OUT_DIR / "bdt_tight_id_isolation_closure_speaker_script.md"
    manifest_path = OUT_DIR / "bdt_tight_id_isolation_closure_manifest.json"
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
            "tight_id_definition": "same-cache centrality-bin WP80: score threshold is the 20th percentile of truth-photon scores in each broad centrality bin; tight if score > threshold",
            "target_signal_efficiency": TARGET_SIGNAL_EFFICIENCY,
            "classes": {
                "truth_photons": "is_signal == 1",
                "inclusive_jet_labeled_background": "is_signal == 0",
                "reco_candidate_pool": "all selected rows in the validation mixture",
            },
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
