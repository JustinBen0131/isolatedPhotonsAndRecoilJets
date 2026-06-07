#!/usr/bin/env python3
"""Build a single-slide default-isolation vs BDT-score scatter candidate."""

from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager  # noqa: E402
from matplotlib.colors import to_rgba  # noqa: E402
from matplotlib.patches import FancyBboxPatch, Rectangle  # noqa: E402
import numpy as np  # noqa: E402


REPO = Path(__file__).resolve().parents[4]
SOURCE_REPORT = REPO / "dataOutput/auauTightBDTValidation/model_validation_condor_20260511_194832"
MANIFEST = SOURCE_REPORT / "score_caches.local.list"
OUT_DIR = SOURCE_REPORT / "slideReady/default_eiso_vs_bdt_score_20260604"

SCORE_COLUMN = "score_ptFine_cent7"
ISO_COLUMN = "reco_eiso"
PT_RANGE = (15.0, 35.0)
CENT_RANGE = (0.0, 80.0)
MAX_POINTS_PER_CLASS = 150_000
RNG_SEED = 20260604

W, H = 2560, 1440
DPI = 200

BG = "#fbfaf7"
INK = "#161616"
MUTED = "#59616f"
GRID = "#dfe4ea"
RED = "#c84c4c"
BLUE = "#315f9c"
PURPLE = "#6750a4"
SOFT_YELLOW = "#fff3c7"


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
            "axes.edgecolor": "#222222",
            "axes.linewidth": 1.2,
            "xtick.color": INK,
            "ytick.color": INK,
            "text.color": INK,
            "axes.labelcolor": INK,
            "savefig.facecolor": BG,
        }
    )


def corr(x: np.ndarray, y: np.ndarray) -> float:
    finite = np.isfinite(x) & np.isfinite(y)
    if finite.sum() < 3:
        return math.nan
    xf = x[finite].astype("float64", copy=False)
    yf = y[finite].astype("float64", copy=False)
    return float(np.corrcoef(xf, yf)[0, 1])


def load_arrays() -> dict[str, np.ndarray]:
    if not MANIFEST.exists():
        raise FileNotFoundError(f"Missing score-cache manifest: {MANIFEST}")

    chunks: dict[str, list[np.ndarray]] = {
        "eiso": [],
        "score": [],
        "is_signal": [],
        "cluster_et": [],
        "centrality": [],
    }
    paths = [Path(line.strip()) for line in MANIFEST.read_text().splitlines() if line.strip()]
    missing: list[str] = []
    for path in paths:
        p = path if path.is_absolute() else REPO / path
        if not p.exists():
            missing.append(str(path))
            continue
        data = np.load(p, allow_pickle=True)
        for required in [ISO_COLUMN, SCORE_COLUMN, "is_signal", "cluster_Et", "centrality"]:
            if required not in data.files:
                raise KeyError(f"{p} missing required column {required}")
        chunks["eiso"].append(data[ISO_COLUMN].astype("float32", copy=False))
        chunks["score"].append(data[SCORE_COLUMN].astype("float32", copy=False))
        chunks["is_signal"].append(data["is_signal"].astype("int8", copy=False))
        chunks["cluster_et"].append(data["cluster_Et"].astype("float32", copy=False))
        chunks["centrality"].append(data["centrality"].astype("float32", copy=False))

    if missing:
        raise FileNotFoundError(f"Missing {len(missing)} cache files; first missing: {missing[0]}")

    arrays = {key: np.concatenate(vals) for key, vals in chunks.items()}
    selected = (
        np.isfinite(arrays["eiso"])
        & np.isfinite(arrays["score"])
        & (arrays["cluster_et"] >= PT_RANGE[0])
        & (arrays["cluster_et"] < PT_RANGE[1])
        & (arrays["centrality"] >= CENT_RANGE[0])
        & (arrays["centrality"] < CENT_RANGE[1])
    )
    return {key: val[selected] for key, val in arrays.items()}


def sample_class(mask: np.ndarray, max_points: int, rng: np.random.Generator) -> np.ndarray:
    idx = np.flatnonzero(mask)
    if len(idx) <= max_points:
        return idx
    return np.sort(rng.choice(idx, size=max_points, replace=False))


def binned_median_line(x: np.ndarray, y: np.ndarray, bins: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    centers: list[float] = []
    medians: list[float] = []
    lows: list[float] = []
    highs: list[float] = []
    for lo, hi in zip(bins[:-1], bins[1:]):
        mask = (x >= lo) & (x < hi)
        if mask.sum() < 50:
            continue
        vals = y[mask]
        centers.append(float(0.5 * (lo + hi)))
        medians.append(float(np.nanmedian(vals)))
        lows.append(float(np.nanpercentile(vals, 25)))
        highs.append(float(np.nanpercentile(vals, 75)))
    return np.asarray(centers), np.asarray(medians), np.asarray(lows), np.asarray(highs)


def rounded_box(fig, xywh, facecolor="white", edgecolor="#d4dae3", lw=1.4, radius=0.016):
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


def add_title(fig) -> None:
    fig.text(0.045, 0.94, "Default AuAu isolation vs BDT score", ha="left", va="top", fontsize=30, fontweight="bold")
    fig.text(
        0.045,
        0.900,
        r"Baseline AuAu tight-BDT validation: default reconstructed $E_T^{iso}$ ($reco\_eiso$) against $score\_ptFine\_cent7$.",
        ha="left",
        va="top",
        fontsize=15.5,
        color=MUTED,
    )
    fig.add_artist(Rectangle((0.045, 0.867), 0.91, 0.003, transform=fig.transFigure, color="#d7dce3", lw=0))


def add_metric_card(fig, x, y, w, h, label, value, detail, color) -> None:
    rounded_box(fig, (x, y, w, h), facecolor="white", edgecolor="#d8dfe8")
    fig.text(x + 0.018, y + h - 0.030, label, fontsize=15.5, fontweight="bold", color=color, va="top")
    if detail:
        fig.text(x + w - 0.018, y + h - 0.031, detail, fontsize=12.2, color=MUTED, ha="right", va="top")
    fig.text(x + 0.018, y + 0.033, value, fontsize=26.5, fontweight="bold", color=color, va="bottom")


def make_slide() -> tuple[Path, Path, Path]:
    setup_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    arrays = load_arrays()

    eiso = arrays["eiso"]
    score = arrays["score"]
    is_signal = arrays["is_signal"].astype(bool)
    is_background = ~is_signal
    rng = np.random.default_rng(RNG_SEED)
    sig_idx = sample_class(is_signal, MAX_POINTS_PER_CLASS, rng)
    bkg_idx = sample_class(is_background, MAX_POINTS_PER_CLASS, rng)

    xlo, xhi = np.nanpercentile(eiso, [0.5, 99.5])
    pad = 0.06 * (xhi - xlo)
    xlo = float(xlo - pad)
    xhi = float(xhi + pad)
    ylo, yhi = 0.0, min(1.0, float(np.nanpercentile(score, 99.8) + 0.06))

    rho_all = corr(eiso, score)
    rho_sig = corr(eiso[is_signal], score[is_signal])
    rho_bkg = corr(eiso[is_background], score[is_background])

    bins = np.linspace(xlo, xhi, 24)
    xmid, ymed, yq25, yq75 = binned_median_line(eiso, score, bins)

    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)
    add_title(fig)

    ax = fig.add_axes([0.075, 0.285, 0.665, 0.535])
    ax.scatter(
        eiso[bkg_idx],
        score[bkg_idx],
        s=6,
        color=to_rgba(BLUE, 0.075),
        edgecolors="none",
        label="Inclusive jet candidates",
    )
    ax.scatter(
        eiso[sig_idx],
        score[sig_idx],
        s=6,
        color=to_rgba(RED, 0.085),
        edgecolors="none",
        label="Truth photons",
    )
    if len(xmid):
        ax.fill_between(xmid, yq25, yq75, color=to_rgba(INK, 0.12), lw=0, label="middle 50% by isolation")
        ax.plot(xmid, ymed, color=INK, lw=3.0, label="median BDT score")

    ax.set_xlim(xlo, xhi)
    ax.set_ylim(ylo, yhi)
    ax.set_xlabel(r"default reconstructed $E_T^{iso}$ ($reco\_eiso$) [GeV]", fontsize=19)
    ax.set_ylabel("BDT score", fontsize=19)
    ax.tick_params(axis="both", labelsize=16)
    ax.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    leg = ax.legend(loc="upper right", frameon=True, fontsize=13.5, handlelength=1.5)
    leg.get_frame().set_facecolor("white")
    leg.get_frame().set_edgecolor("#d8dfe8")
    leg.get_frame().set_linewidth(1.0)
    ax.text(
        0.018,
        0.965,
        r"$15 \leq E_T^{cluster} < 35$ GeV, $0 \leq centrality < 80$",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=14.5,
        color=MUTED,
        bbox={"facecolor": "white", "edgecolor": "#d8dfe8", "boxstyle": "round,pad=0.32"},
    )

    add_metric_card(
        fig,
        0.775,
        0.675,
        0.18,
        0.145,
        "All candidates",
        rf"$\rho$ = {rho_all:+.2f}",
        "",
        PURPLE,
    )
    add_metric_card(
        fig,
        0.775,
        0.505,
        0.18,
        0.145,
        "Truth photons",
        rf"$\rho$ = {rho_sig:+.2f}",
        "",
        RED,
    )
    add_metric_card(
        fig,
        0.775,
        0.335,
        0.18,
        0.145,
        "Inclusive jets",
        rf"$\rho$ = {rho_bkg:+.2f}",
        "",
        BLUE,
    )

    rounded_box(fig, (0.075, 0.050, 0.88, 0.092), facecolor=SOFT_YELLOW, edgecolor="#e6c76d")
    fig.text(0.100, 0.117, "Clean read", fontsize=17.5, fontweight="bold", color="#6f4d00", va="top")
    fig.text(
        0.100,
        0.083,
        "The direct 2D view shows a modest full-sample anticorrelation; signal and inclusive-jet candidates occupy different score bands.",
        fontsize=14.2,
        color=INK,
        va="top",
    )

    png = OUT_DIR / "the11_default_eiso_vs_bdt_score_scatter_slide.png"
    fig.savefig(png, dpi=DPI)
    plt.close(fig)

    script = OUT_DIR / "the11_default_eiso_vs_bdt_score_scatter_script.md"
    script.write_text(
        "\n".join(
            [
                "# THE-11 Default Isolation vs BDT Score",
                "",
                "This is the simplest visual diagnostic: put the baseline AuAu reconstructed isolation on the x-axis and the BDT score on the y-axis.",
                "",
                f"The full selected validation sample gives Pearson rho = {rho_all:+.3f}. Truth photons are nearly flat by themselves, rho = {rho_sig:+.3f}, while the inclusive jet candidates have a weak tail-driven anticorrelation, rho = {rho_bkg:+.3f}.",
                "",
                "The audience-facing point is measured, not overstated: isolation and BDT score are not dramatically locked together in this default validation view, but the full-sample trend is visible and the truth classes occupy different regions.",
                "",
            ]
        )
    )

    manifest = OUT_DIR / "the11_default_eiso_vs_bdt_score_scatter_manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "schema": "THE11_DEFAULT_EISO_BDT_SCATTER_V1",
                "source_report": str(SOURCE_REPORT),
                "score_cache_manifest": str(MANIFEST),
                "score_column": SCORE_COLUMN,
                "isolation_column": ISO_COLUMN,
                "selection": {
                    "cluster_Et_min_inclusive": PT_RANGE[0],
                    "cluster_Et_max_exclusive": PT_RANGE[1],
                    "centrality_min_inclusive": CENT_RANGE[0],
                    "centrality_max_exclusive": CENT_RANGE[1],
                },
                "entries": {
                    "selected": int(len(eiso)),
                    "signal": int(is_signal.sum()),
                    "background": int(is_background.sum()),
                    "plotted_signal": int(len(sig_idx)),
                    "plotted_background": int(len(bkg_idx)),
                },
                "correlations": {
                    "pearson_all": rho_all,
                    "pearson_signal": rho_sig,
                    "pearson_background": rho_bkg,
                },
                "outputs": {
                    "png": str(png),
                    "speaker_script": str(script),
                },
                "notes": [
                    "This slide uses local row-level AuAu tight-BDT validation caches for the default reco_eiso variable.",
                    "The separate full-stat no-isolation score diagnostic currently has local aggregate correlations for reco_eiso_r30/r40, not local row-level reco_eiso scatter inputs.",
                    "Google Slides was not mutated.",
                ],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    return png, script, manifest


def main() -> None:
    png, script, manifest = make_slide()
    print(f"wrote {png}")
    print(f"wrote {script}")
    print(f"wrote {manifest}")


if __name__ == "__main__":
    main()
