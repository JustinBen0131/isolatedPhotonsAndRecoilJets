#!/usr/bin/env python3
"""Build focused THE-11 class-split isolation/BDT correlation slides."""

from __future__ import annotations

import csv
import json
import math
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402
from matplotlib.patches import FancyBboxPatch, Rectangle  # noqa: E402
import numpy as np  # noqa: E402


REPO = Path(__file__).resolve().parents[4]
CORR_CSV = REPO / (
    "dataOutput/auauMLDiagnosticRuns/global_etcent_inclusive3_sixpack_20260516_135439/"
    "slideReady/isolation_feature_correlations/noiso_score_diagnostic/"
    "isolation_feature_score_correlations.csv"
)
ROW_REPORT = REPO / "dataOutput/auauTightBDTValidation/model_validation_condor_20260511_194832"
ROW_MANIFEST = ROW_REPORT / "score_caches.local.list"
OUT_DIR = REPO / (
    "dataOutput/auauMLDiagnosticRuns/global_etcent_inclusive3_sixpack_20260516_135439/"
    "slideReady/isolation_feature_correlations/noiso_score_diagnostic/"
    "the11_class_split_correlation_story_20260604"
)

TARGET = "score_globalEtCent1535_bdt_noIso_ptCent7"
ROW_SCORE = "score_ptFine_cent7"
ROW_EISO = "reco_eiso"
PT_RANGE = (15.0, 35.0)
CENTRALITY = [
    ("cent_0_20", "0-20%"),
    ("cent_20_50", "20-50%"),
    ("cent_50_80", "50-80%"),
]
CENT_BINS = [
    (0.0, 20.0, "0-20%"),
    (20.0, 50.0, "20-50%"),
    (50.0, 80.0, "50-80%"),
]
CLASSES = [
    ("all", "all candidates", "#6a51a3"),
    ("signal", "truth photons", "#c94b4b"),
    ("background", "inclusive jets", "#3264a8"),
]
CONES = [
    ("reco_eiso_r30", "0.3", "slide01_r30_class_split_correlation"),
    ("reco_eiso_r40", "0.4", "backup01_r40_class_split_correlation"),
]

W, H = 2560, 1440
DPI = 200
BG = "#ffffff"
INK = "#161616"
MUTED = "#59616f"
GRID = "#dfe4ea"
TEAL = "#1f8a83"
GOLD = "#c9932e"
PURPLE = "#6750a4"
RED = "#c94b4b"
BLUE = "#3264a8"


def setup_style() -> None:
    available = {f.name for f in font_manager.fontManager.ttflist}
    for name in ["Times New Roman", "Times", "DejaVu Serif"]:
        if name in available:
            plt.rcParams["font.family"] = name
            break
    plt.rcParams.update(
        {
            "figure.facecolor": BG,
            "savefig.facecolor": BG,
            "axes.facecolor": "white",
            "axes.edgecolor": INK,
            "axes.linewidth": 1.15,
            "axes.labelcolor": INK,
            "xtick.color": INK,
            "ytick.color": INK,
            "text.color": INK,
            "mathtext.default": "regular",
        }
    )


def load_rows() -> list[dict[str, str]]:
    with CORR_CSV.open() as handle:
        return list(csv.DictReader(handle))


def load_density_arrays() -> dict[str, np.ndarray]:
    chunks: dict[str, list[np.ndarray]] = {k: [] for k in ["eiso", "score", "is_signal", "et", "cent"]}
    for line in ROW_MANIFEST.read_text().splitlines():
        if not line.strip():
            continue
        path = REPO / line.strip()
        data = np.load(path, allow_pickle=True)
        for key in [ROW_EISO, ROW_SCORE, "is_signal", "cluster_Et", "centrality"]:
            if key not in data.files:
                raise KeyError(f"{path} missing {key}")
        chunks["eiso"].append(data[ROW_EISO].astype("float32", copy=False))
        chunks["score"].append(data[ROW_SCORE].astype("float32", copy=False))
        chunks["is_signal"].append(data["is_signal"].astype("int8", copy=False))
        chunks["et"].append(data["cluster_Et"].astype("float32", copy=False))
        chunks["cent"].append(data["centrality"].astype("float32", copy=False))
    arrays = {key: np.concatenate(vals) for key, vals in chunks.items()}
    mask = (
        np.isfinite(arrays["eiso"])
        & np.isfinite(arrays["score"])
        & (arrays["et"] >= PT_RANGE[0])
        & (arrays["et"] < PT_RANGE[1])
        & (arrays["cent"] >= 0.0)
        & (arrays["cent"] < 80.0)
    )
    return {key: val[mask] for key, val in arrays.items()}


def value(rows: list[dict[str, str]], scope: str, klass: str, iso: str, metric: str) -> float:
    matches = [
        row
        for row in rows
        if row["scope"] == scope
        and row["class"] == klass
        and row["iso_variable"] == iso
        and row["target_variable"] == TARGET
    ]
    if not matches:
        return math.nan
    return float(matches[0][metric])


def add_header(fig, title: str, subtitle: str) -> None:
    fig.text(0.045, 0.935, title, ha="left", va="top", fontsize=27, fontweight="bold")
    fig.text(0.045, 0.895, subtitle, ha="left", va="top", fontsize=15.5, color=MUTED)
    fig.add_artist(Rectangle((0.045, 0.865), 0.91, 0.003, transform=fig.transFigure, color="#d7dce3", lw=0))


def rounded_box(fig, xywh, facecolor="white", edgecolor="#d9dfe8", lw=1.4, radius=0.016):
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


def add_bottom_cards(fig, cards: list[tuple[str, str, str]]) -> None:
    y, h = 0.045, 0.125
    x0, gap = 0.055, 0.020
    w = (0.89 - gap * 2) / 3
    for i, (label, body, color) in enumerate(cards):
        x = x0 + i * (w + gap)
        rounded_box(fig, (x, y, w, h))
        fig.text(x + 0.018, y + h - 0.030, label, fontsize=15.2, fontweight="bold", color=color, va="top")
        fig.text(
            x + 0.018,
            y + h - 0.067,
            textwrap.fill(body, width=48),
            fontsize=10.8,
            linespacing=1.10,
            va="top",
        )


def style_density_panel(ax, title: str, xlim: tuple[float, float], xlabel: str | None, ylabel: str | None) -> None:
    ax.set_title(title, fontsize=13.8, fontweight="bold", pad=6)
    ax.set_xlim(*xlim)
    ax.set_ylim(0.0, 1.0)
    ax.grid(True, color=GRID, lw=0.7)
    ax.tick_params(labelsize=10.2)
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=12.0)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=12.4)


def make_density_slide(data: dict[str, np.ndarray]) -> tuple[Path, Path, dict]:
    xlim = tuple(float(v) for v in np.nanpercentile(data["eiso"], [0.5, 99.5]))
    signal = data["is_signal"].astype(bool)
    medians = {"signal": [], "background": []}

    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)
    add_header(
        fig,
        "Class densities show the BDT-isolation geometry",
        r"Row-level visual companion, $15 \leq E_T^{cluster}<35$ GeV; default Au+Au isolation uses $\Delta R<0.3$.",
    )
    axes = []
    for row in range(2):
        for col in range(3):
            axes.append(fig.add_axes([0.060 + col * 0.300, 0.575 - row * 0.265, 0.255, 0.180]))

    for col, (lo, hi, label) in enumerate(CENT_BINS):
        cent_mask = (data["cent"] >= lo) & (data["cent"] < hi)
        for row, (class_mask, cmap, class_label, key) in enumerate(
            [
                (signal, "Reds", "Truth photons", "signal"),
                (~signal, "Blues", "Inclusive jets", "background"),
            ]
        ):
            mask = cent_mask & class_mask
            ax = axes[row * 3 + col]
            ax.hist2d(
                data["eiso"][mask],
                data["score"][mask],
                bins=(58, 44),
                range=[xlim, [0.0, 1.0]],
                norm=LogNorm(vmin=1),
                cmap=cmap,
            )
            style_density_panel(
                ax,
                f"{label} | {class_label}",
                xlim,
                xlabel=r"reco $E_T^{iso}$, $\Delta R<0.3$ [GeV]" if row == 1 else None,
                ylabel="BDT score" if col == 0 else None,
            )
            medians[key].append(float(np.nanmedian(data["score"][mask])))

    add_bottom_cards(
        fig,
        [
            ("Visual anatomy", "Split classes first: the photon band and inclusive tail occupy different score-isolation regions.", TEAL),
            ("0-20% read", f"Central-bin median score: photons {medians['signal'][0]:.2f}, inclusive jets {medians['background'][0]:.2f}.", PURPLE),
            ("Use with rho slide", "This shows the population geometry; the next slide quantifies the correlation split.", GOLD),
        ],
    )
    png = OUT_DIR / "slide02_r30_centrality_class_density.png"
    fig.savefig(png)
    plt.close(fig)

    script = OUT_DIR / "slide02_r30_centrality_class_density_script.md"
    script.write_text(
        "\n\n".join(
            [
                "# THE-11 Script - centrality class density R=0.3",
                "This is the visual companion to the class-split correlation slide.",
                "Truth photons and inclusive jets are plotted separately so the audience can see why the aggregate BDT-isolation relationship should not be read as one single population.",
                f"In the 0-20 percent centrality bin, the median score is {medians['signal'][0]:.3f} for truth photons and {medians['background'][0]:.3f} for inclusive jets.",
            ]
        )
        + "\n"
    )
    return png, script, {"score_medians": medians, "selected_entries": int(len(data["score"]))}


def plot_class_bars(ax, vals: dict[str, list[float]], cone_label: str) -> None:
    x = np.arange(len(CENTRALITY))
    width = 0.22
    ax.axhspan(-0.07, 0.03, color="#f5e8e8", zorder=0, alpha=0.75)
    ax.axhline(0.0, color="#777777", lw=1.1)
    for offset, (klass, label, color) in zip([-width, 0.0, width], CLASSES):
        bars = ax.bar(x + offset, vals[klass], width=width, color=color, label=label, zorder=3)
        for bar, rho in zip(bars, vals[klass]):
            y = rho + 0.025 if rho < -0.12 else rho + 0.014
            va = "bottom"
            color_text = "white" if rho < -0.16 else INK
            ax.text(
                bar.get_x() + bar.get_width() / 2.0,
                y,
                f"{rho:+.2f}",
                ha="center",
                va=va,
                fontsize=11,
                color=color_text,
                fontweight="bold",
                zorder=4,
            )
    ax.text(
        2.54,
        -0.020,
        "truth-photon\nnear-flat band",
        ha="right",
        va="center",
        fontsize=10.5,
        color="#8f3f45",
    )
    ax.set_xticks(x)
    ax.set_xticklabels([label for _, label in CENTRALITY], fontsize=13)
    ax.set_ylim(-0.60, 0.08)
    ax.set_xlim(-0.55, 2.70)
    ax.set_ylabel("Pearson rho with no-isolation BDT score", fontsize=14)
    ax.set_title(rf"reco $E_T^{{iso}}$, $\Delta R < {cone_label}$", fontsize=17, fontweight="bold")
    ax.grid(axis="y", color=GRID, lw=0.8, zorder=1)
    ax.legend(frameon=False, fontsize=12, loc="lower left", ncol=3)


def make_slide(rows: list[dict[str, str]], iso: str, cone_label: str, stem: str) -> tuple[Path, Path, dict]:
    vals = {
        klass: [value(rows, scope, klass, iso, "pearson") for scope, _ in CENTRALITY]
        for klass, _, _ in CLASSES
    }
    spearman = {
        klass: [value(rows, scope, klass, iso, "spearman") for scope, _ in CENTRALITY]
        for klass, _, _ in CLASSES
    }
    is_backup = stem.startswith("backup")
    title = (
        f"Backup: class split persists for R={cone_label}"
        if is_backup
        else "Class split explains the BDT-isolation correlation"
    )
    subtitle = (
        "Exact aggregate diagnostic: isolation versus the no-isolation BDT score, split by truth class and centrality."
    )
    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)
    add_header(fig, title, subtitle)
    ax = fig.add_axes([0.085, 0.280, 0.830, 0.500])
    plot_class_bars(ax, vals, cone_label)

    all_central = vals["all"][0]
    signal_max = max(abs(v) for v in vals["signal"])
    background_span = (min(vals["background"]), max(vals["background"]))
    add_bottom_cards(
        fig,
        [
            (
                "Central-bin read",
                f"In 0-20%, all candidates have rho={all_central:+.2f}, while truth photons are only {vals['signal'][0]:+.2f}.",
                TEAL,
            ),
            (
                "What splits the story",
                f"Truth photons stay near flat (max |rho|={signal_max:.2f}); inclusive jets carry the visible anticorrelation.",
                PURPLE,
            ),
            (
                "Physical meaning",
                f"Higher nearby activity is mainly a background/context tail, not the BDT simply replaying isolation.",
                GOLD,
            ),
        ],
    )
    png = OUT_DIR / f"{stem}.png"
    fig.savefig(png)
    plt.close(fig)

    script = OUT_DIR / f"{stem}_script.md"
    script.write_text(
        "\n\n".join(
            [
                f"# THE-11 Script - class split correlation R={cone_label}",
                "This slide is the tightened version of the isolation story.",
                f"For R={cone_label}, the all-candidate correlation with the no-isolation BDT score is negative in each centrality bin, but the truth-photon-only correlation stays close to zero.",
                f"The inclusive-jet correlation ranges from {background_span[0]:+.3f} to {background_span[1]:+.3f}, so the full-sample trend is best read as a background/context-tail effect rather than a signal-photon law.",
            ]
        )
        + "\n"
    )
    return png, script, {"pearson": vals, "spearman": spearman}


def main() -> None:
    setup_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = load_rows()
    density_data = load_density_arrays()
    outputs = []
    png, script, stats = make_density_slide(density_data)
    outputs.append({"png": str(png), "script": str(script), "stats": stats, "source": "row_level_default_score"})
    print(f"wrote {png}")
    print(f"wrote {script}")
    for iso, cone_label, stem in CONES:
        png, script, stats = make_slide(rows, iso, cone_label, stem)
        outputs.append({"png": str(png), "script": str(script), "stats": stats, "source": "exact_aggregate_noiso_correlation"})
        print(f"wrote {png}")
        print(f"wrote {script}")
    manifest = OUT_DIR / "the11_class_split_correlation_story_manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "schema": "THE11_CLASS_SPLIT_CORRELATION_STORY_V1",
                "correlation_csv": str(CORR_CSV),
                "row_level_density_source_report": str(ROW_REPORT),
                "row_level_density_manifest": str(ROW_MANIFEST),
                "row_level_density_score_column": ROW_SCORE,
                "row_level_density_isolation_column": ROW_EISO,
                "row_level_density_selection": {
                    "cluster_Et_min_inclusive": PT_RANGE[0],
                    "cluster_Et_max_exclusive": PT_RANGE[1],
                    "centrality_min_inclusive": 0.0,
                    "centrality_max_exclusive": 80.0,
                    "selected_entries": int(len(density_data["score"])),
                    "signal_entries": int(density_data["is_signal"].sum()),
                    "background_entries": int((~density_data["is_signal"].astype(bool)).sum()),
                },
                "row_level_density_caveat": "The class-density companion uses local default validation score caches. The correlation bars use the exact aggregate no-isolation diagnostic.",
                "target_variable": TARGET,
                "centrality_bins": CENTRALITY,
                "class_order": [label for _, label, _ in CLASSES],
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
