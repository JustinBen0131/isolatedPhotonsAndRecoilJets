#!/usr/bin/env python3
"""Make THE-11 shower-width versus isolation correlation example plots."""

from __future__ import annotations

import json
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.patheffects as pe  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, LogNorm, TwoSlopeNorm  # noqa: E402
import numpy as np  # noqa: E402
from PIL import Image, ImageDraw, ImageFont  # noqa: E402


REPO = Path(__file__).resolve().parents[3]
REPORT = REPO / "dataOutput/auauTightBDTValidation/model_validation_condor_20260511_194832"
MANIFEST = REPORT / "score_caches.local.list"
OUT_DIR = REPORT / "slideReady/default_eiso_vs_bdt_score_20260604/plot_variants_width_isolation_20260605"
BDT_WP_CSV = REPO / (
    "dataOutput/auauTightBDTValidation/THE8_branchA_ladder_scorecache_fullstat_20260527/"
    "the41_centdep_wp_slides_20260604/the41_centdep_bdt_wp_summary.csv"
)
ISO_FINE_CSV = REPO / (
    "dataOutput/auauMLDiagnosticRuns/global_etcent_inclusive3_sixpack_20260516_135439/"
    "slideReady/isolation_wp_checks/"
    "auau_sliding_isolation_wp_r03_r04_eff70_80_90_fine7_no_numbers_coefficients.csv"
)

EISO = "reco_eiso"
SCORE = "score_ptFine_cent7"
FEATURES = {
    "weta": {
        "column": "cluster_weta_cogx",
        "label": r"cluster $w_{\eta}^{COGx}$",
        "short": r"$w_{\eta}$",
        "ylim": (0.0, 0.60),
        "title": r"$w_{\eta}$ versus reconstructed isolation",
    },
    "wphi": {
        "column": "cluster_wphi_cogx",
        "label": r"cluster $w_{\phi}^{COGx}$",
        "short": r"$w_{\phi}$",
        "ylim": (0.0, 0.80),
        "title": r"$w_{\phi}$ versus reconstructed isolation",
    },
}
PT_RANGE = (15.0, 35.0)
CENT_BINS = [(0.0, 20.0, "0-20%"), (20.0, 50.0, "20-50%"), (50.0, 80.0, "50-80%")]
CENT7 = [
    (0.0, 10.0, "0-10%"),
    (10.0, 20.0, "10-20%"),
    (20.0, 30.0, "20-30%"),
    (30.0, 40.0, "30-40%"),
    (40.0, 50.0, "40-50%"),
    (50.0, 60.0, "50-60%"),
    (60.0, 80.0, "60-80%"),
]
ET8 = [
    (15.0, 17.0, "15-17"),
    (17.0, 19.0, "17-19"),
    (19.0, 21.0, "19-21"),
    (21.0, 23.0, "21-23"),
    (23.0, 25.0, "23-25"),
    (25.0, 27.0, "25-27"),
    (27.0, 30.0, "27-30"),
    (30.0, 35.0, "30-35"),
]
EISO_RANGE = (-12.0, 18.0)
EISO_DIAG_RANGE = (-10.0, 14.0)
EISO_TREND_EDGES = np.array([-10.0, -6.0, -3.0, -1.0, 1.0, 3.0, 6.0, 10.0, 16.0])
EISO_DIAG_MEDIAN_EDGES = np.array([-10.0, -7.0, -5.0, -3.0, -1.0, 1.0, 3.0, 5.0, 8.0, 11.0, 14.0])

WIDE = (12.8, 7.2)
DPI = 200
WHITE = "#ffffff"
INK = "#171717"
MUTED = "#59616f"
GRID = "#dfe4ea"
TRUTH = "#b0182d"
ALL = "#34495e"
TEAL = "#168f7a"
ORANGE = "#cf6a1b"


def root_kbird() -> LinearSegmentedColormap:
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
    cmap.set_bad("#f7f7f7")
    cmap.set_under("#ffffff")
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
    required = [EISO, SCORE, "cluster_weta_cogx", "cluster_wphi_cogx", "cluster_Et", "centrality", "is_signal"]
    chunks: dict[str, list[np.ndarray]] = {key: [] for key in required}
    for raw_line in MANIFEST.read_text().splitlines():
        line = raw_line.strip()
        if not line:
            continue
        path = Path(line)
        cache_path = path if path.is_absolute() else REPO / path
        data = np.load(cache_path, allow_pickle=True)
        missing = [key for key in required if key not in data.files]
        if missing:
            raise KeyError(f"{cache_path} missing required columns: {missing}")
        for key in required:
            dtype = "int8" if key == "is_signal" else "float32"
            chunks[key].append(data[key].astype(dtype, copy=False))

    arrays = {key: np.concatenate(values) for key, values in chunks.items()}
    selected = (
        np.isfinite(arrays[EISO])
        & np.isfinite(arrays[SCORE])
        & np.isfinite(arrays["cluster_weta_cogx"])
        & np.isfinite(arrays["cluster_wphi_cogx"])
        & (arrays["cluster_Et"] >= PT_RANGE[0])
        & (arrays["cluster_Et"] < PT_RANGE[1])
        & (arrays["centrality"] >= 0.0)
        & (arrays["centrality"] < 80.0)
    )
    return {key: values[selected] for key, values in arrays.items()}


def pearson(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 3:
        return float("nan")
    x = x.astype("float64", copy=False)
    y = y.astype("float64", copy=False)
    x = x - np.nanmean(x)
    y = y - np.nanmean(y)
    denom = np.sqrt(np.nansum(x * x) * np.nansum(y * y))
    if denom <= 0:
        return float("nan")
    return float(np.nansum(x * y) / denom)


def rankdata(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype="float64")
    ranks[order] = np.arange(len(values), dtype="float64")
    return ranks


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 3:
        return float("nan")
    return pearson(rankdata(x), rankdata(y))


def hist2d_density(x: np.ndarray, y: np.ndarray, xedges: np.ndarray, yedges: np.ndarray) -> np.ndarray:
    h, _, _ = np.histogram2d(x, y, bins=[xedges, yedges])
    h = h.astype("float64")
    total = h.sum()
    if total > 0:
        h /= total
    h[h <= 0] = np.nan
    return h


def sample_mask(arrays: dict[str, np.ndarray], sample: str) -> np.ndarray:
    if sample == "all":
        return np.ones(len(arrays[EISO]), dtype=bool)
    if sample == "truth":
        return arrays["is_signal"].astype(bool)
    if sample == "inclusive":
        return ~arrays["is_signal"].astype(bool)
    raise ValueError(sample)


def panel_stats(arrays: dict[str, np.ndarray], mask: np.ndarray, feature_col: str) -> dict[str, float | int]:
    if mask.sum() == 0:
        return {
            "entries": 0,
            "pearson": float("nan"),
            "spearman": float("nan"),
            "median_feature": float("nan"),
            "median_eiso": float("nan"),
            "tail_feature_gt_95_all": float("nan"),
        }
    return {
        "entries": int(mask.sum()),
        "pearson": pearson(arrays[EISO][mask], arrays[feature_col][mask]),
        "spearman": spearman(arrays[EISO][mask], arrays[feature_col][mask]),
        "median_feature": float(np.nanmedian(arrays[feature_col][mask])),
        "median_eiso": float(np.nanmedian(arrays[EISO][mask])),
    }


def add_title(fig: plt.Figure, title: str, subtitle: str) -> None:
    fig.text(0.050, 0.955, title, ha="left", va="top", fontsize=22.2, fontweight="bold")
    fig.text(0.050, 0.907, subtitle, ha="left", va="top", fontsize=12.4, color=MUTED)


def save_density_facets(arrays: dict[str, np.ndarray], feature_key: str) -> tuple[Path, dict[str, object]]:
    info = FEATURES[feature_key]
    feature_col = str(info["column"])
    xedges = np.linspace(EISO_RANGE[0], EISO_RANGE[1], 76)
    yedges = np.linspace(info["ylim"][0], info["ylim"][1], 66)
    panels: dict[str, dict[str, object]] = {}
    hists: list[np.ndarray] = []
    for sample in ["all", "truth"]:
        for lo, hi, label in CENT_BINS:
            mask = sample_mask(arrays, sample) & (arrays["centrality"] >= lo) & (arrays["centrality"] < hi)
            h = hist2d_density(arrays[EISO][mask], arrays[feature_col][mask], xedges, yedges)
            hists.append(h)
            panels[f"{sample}_{label}"] = panel_stats(arrays, mask, feature_col)
            panels[f"{sample}_{label}"]["visible_fraction"] = float(
                np.mean(
                    (arrays[EISO][mask] >= EISO_RANGE[0])
                    & (arrays[EISO][mask] <= EISO_RANGE[1])
                    & (arrays[feature_col][mask] >= info["ylim"][0])
                    & (arrays[feature_col][mask] <= info["ylim"][1])
                )
            )

    positive = np.concatenate([h[np.isfinite(h)] for h in hists if np.isfinite(h).any()])
    vmax = float(np.nanpercentile(positive, 99.6)) if positive.size else 1.0
    norm = LogNorm(vmin=2.0e-5, vmax=max(vmax, 2.0e-4))

    fig, axes = plt.subplots(2, 3, figsize=WIDE, dpi=DPI, sharex=True, sharey=True)
    fig.subplots_adjust(left=0.106, right=0.880, top=0.800, bottom=0.140, wspace=0.105, hspace=0.160)
    add_title(
        fig,
        f"{info['title']} shows whether isolation is tied to shower width",
        r"Rows compare all scored candidates to truth-tagged photons; columns split centrality. "
        r"Color is panel-normalized density on a log scale, $15 \leq E_T^{cluster}<35$ GeV.",
    )
    mesh = None
    for row, sample in enumerate(["all", "truth"]):
        for col, (lo, hi, cent_label) in enumerate(CENT_BINS):
            ax = axes[row, col]
            mask = sample_mask(arrays, sample) & (arrays["centrality"] >= lo) & (arrays["centrality"] < hi)
            h = hist2d_density(arrays[EISO][mask], arrays[feature_col][mask], xedges, yedges)
            mesh = ax.pcolormesh(xedges, yedges, h.T, cmap=KBIRD, norm=norm, shading="auto")
            rho = pearson(arrays[EISO][mask], arrays[feature_col][mask])
            ax.text(
                0.035,
                0.935,
                rf"$\rho$ = {rho:+.2f}",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=10.4,
                bbox=dict(boxstyle="round,pad=0.22", facecolor="white", edgecolor="#cfd7e3", alpha=0.92),
            )
            ax.set_xlim(*EISO_RANGE)
            ax.set_ylim(*info["ylim"])
            ax.grid(True, color=GRID, lw=0.58)
            ax.tick_params(labelsize=9.6)
            if row == 0:
                ax.set_title(f"{cent_label} centrality", fontsize=12.8, fontweight="bold", pad=8)
            if row == 1:
                ax.set_xlabel(r"reco $E_T^{iso}$, $\Delta R<0.3$ [GeV]", fontsize=10.8)
            if col == 0:
                ax.set_ylabel(info["label"], fontsize=10.8)
            else:
                ax.set_ylabel("")
    fig.text(0.055, 0.600, "all candidates", ha="center", va="center", rotation=90, fontsize=13.8, fontweight="bold", color=ALL)
    fig.text(0.055, 0.335, "truth-tagged photons", ha="center", va="center", rotation=90, fontsize=13.8, fontweight="bold", color=TRUTH)
    cax = fig.add_axes([0.900, 0.235, 0.018, 0.470])
    fig.colorbar(mesh, cax=cax, label="share per bin, log scale")
    out = OUT_DIR / f"01_{feature_key}_vs_eiso_density_by_centrality.png"
    fig.savefig(out)
    plt.close(fig)
    return out, panels


def median_by_eiso(
    arrays: dict[str, np.ndarray],
    mask: np.ndarray,
    feature_col: str,
    edges: np.ndarray = EISO_TREND_EDGES,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    centers = 0.5 * (edges[:-1] + edges[1:])
    med = np.full(len(centers), np.nan, dtype="float64")
    qlo = np.full(len(centers), np.nan, dtype="float64")
    qhi = np.full(len(centers), np.nan, dtype="float64")
    for idx, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
        bin_mask = mask & (arrays[EISO] >= lo) & (arrays[EISO] < hi)
        if bin_mask.sum() < 40:
            continue
        vals = arrays[feature_col][bin_mask]
        med[idx] = np.nanmedian(vals)
        qlo[idx], qhi[idx] = np.nanpercentile(vals, [25, 75])
    return centers, med, np.vstack([med - qlo, qhi - med])


def isolation_wp90_bands() -> dict[str, tuple[float, float, float]]:
    fine_rows: list[tuple[float, float, float]] = []
    with ISO_FINE_CSV.open() as handle:
        for row in csv.DictReader(handle):
            if row["cone"] != "R0.3" or abs(float(row["efficiency"]) - 0.90) > 1e-9:
                continue
            fine_rows.append((float(row["cent_min"]), float(row["cent_max"]), float(row["threshold_gev"])))

    out: dict[str, tuple[float, float, float]] = {}
    for broad_lo, broad_hi, label in CENT_BINS:
        vals = [thr for lo, hi, thr in fine_rows if lo >= broad_lo and hi <= broad_hi]
        if not vals:
            continue
        out[label] = (float(min(vals)), float(max(vals)), float(np.mean(vals)))
    return out


def configured_bdt_wp80() -> dict[tuple[float, float, float, float], float]:
    out: dict[tuple[float, float, float, float], float] = {}
    with BDT_WP_CSV.open() as handle:
        for row in csv.DictReader(handle):
            if row["row_type"] != "et_bin_point" or row["wp_label"] != "WP80":
                continue
            key = (
                float(row["centrality_min"]),
                float(row["centrality_max"]),
                float(row["et_min"]),
                float(row["et_max"]),
            )
            out[key] = float(row["threshold"])
    return out


def bdt_wp80_pass_mask(arrays: dict[str, np.ndarray]) -> np.ndarray:
    thresholds = configured_bdt_wp80()
    passed = np.zeros(len(arrays[SCORE]), dtype=bool)
    for cent_lo, cent_hi, _ in CENT7:
        cent_mask = (arrays["centrality"] >= cent_lo) & (arrays["centrality"] < cent_hi)
        for et_lo, et_hi, _ in ET8:
            key = (cent_lo, cent_hi, et_lo, et_hi)
            if key not in thresholds:
                raise KeyError(f"Missing WP80 threshold for {key}")
            mask = cent_mask & (arrays["cluster_Et"] >= et_lo) & (arrays["cluster_Et"] < et_hi)
            passed[mask] = arrays[SCORE][mask] >= thresholds[key]
    return passed


def median_width_on_side(
    arrays: dict[str, np.ndarray],
    mask: np.ndarray,
    feature_col: str,
    cut: float,
    clean_side: bool,
) -> float:
    side = arrays[EISO] < cut if clean_side else arrays[EISO] >= cut
    vals = arrays[feature_col][mask & side]
    if len(vals) == 0:
        return float("nan")
    return float(np.nanmedian(vals))


def save_020_before_after_wp80(arrays: dict[str, np.ndarray], feature_key: str) -> tuple[Path, dict[str, object]]:
    info = FEATURES[feature_key]
    feature_col = str(info["column"])
    cut_min, cut_max, cut_mean = isolation_wp90_bands()["0-20%"]
    wp80_pass = bdt_wp80_pass_mask(arrays)
    cent020 = (arrays["centrality"] >= 0.0) & (arrays["centrality"] < 20.0)
    columns = [("before WP80", cent020), ("after WP80", cent020 & wp80_pass)]
    rows = [("truth", "truth-tagged photons", TRUTH), ("inclusive", "inclusive candidates", "#1f5aa6")]

    xedges = np.linspace(EISO_DIAG_RANGE[0], EISO_DIAG_RANGE[1], 82)
    yedges = np.linspace(info["ylim"][0], info["ylim"][1], 74)
    hists: list[np.ndarray] = []
    for sample, _, _ in rows:
        for _, col_mask in columns:
            mask = col_mask & sample_mask(arrays, sample)
            hists.append(hist2d_density(arrays[EISO][mask], arrays[feature_col][mask], xedges, yedges))
    positive = np.concatenate([h[np.isfinite(h)] for h in hists if np.isfinite(h).any()])
    vmax = float(np.nanpercentile(positive, 99.7)) if positive.size else 1.0
    norm = LogNorm(vmin=1.4e-5, vmax=max(vmax, 1.5e-4))

    fig, axes = plt.subplots(2, 2, figsize=WIDE, dpi=DPI, sharex=True, sharey=True)
    fig.subplots_adjust(left=0.115, right=0.880, top=0.775, bottom=0.175, wspace=0.135, hspace=0.170)
    add_title(
        fig,
        f"0-20% WP80 check: {info['short']} versus $R<0.3$ isolation",
        r"Rows compare truth-tagged photons and inclusive candidates. Columns show the same sample before and after the configured BDT WP80 tight selection.",
    )
    stats: dict[str, object] = {}
    mesh = None
    for row_idx, (sample, sample_label, label_color) in enumerate(rows):
        for col_idx, (col_label, col_mask) in enumerate(columns):
            ax = axes[row_idx, col_idx]
            mask = col_mask & sample_mask(arrays, sample)
            h = hist2d_density(arrays[EISO][mask], arrays[feature_col][mask], xedges, yedges)
            mesh = ax.pcolormesh(xedges, yedges, h.T, cmap=KBIRD, norm=norm, shading="auto")
            ax.axvspan(cut_min, cut_max, color=ORANGE, alpha=0.16, zorder=2)
            ax.axvline(cut_mean, color=ORANGE, lw=1.5, ls="--", zorder=3)
            centers, med, _ = median_by_eiso(arrays, mask, feature_col, EISO_DIAG_MEDIAN_EDGES)
            line = ax.plot(centers, med, color=INK, lw=2.35, zorder=5)[0]
            line.set_path_effects([pe.Stroke(linewidth=4.4, foreground="white"), pe.Normal()])
            rho = pearson(arrays[EISO][mask], arrays[feature_col][mask])
            clean_med = median_width_on_side(arrays, mask, feature_col, cut_mean, True)
            tail_med = median_width_on_side(arrays, mask, feature_col, cut_mean, False)
            stats[f"{feature_key}_{sample}_{col_label.replace(' ', '_').lower()}"] = {
                "entries": int(mask.sum()),
                "pearson": rho,
                "spearman": spearman(arrays[EISO][mask], arrays[feature_col][mask]),
                "median_width_clean_side": clean_med,
                "median_width_busy_side": tail_med,
                "r03_wp90_cut_band_gev": [cut_min, cut_max],
                "r03_wp90_cut_mean_gev": cut_mean,
            }
            ax.text(
                0.035,
                0.940,
                rf"$\rho={rho:+.2f}$" + "\n" + f"N = {int(mask.sum()):,}" + "\n" + rf"median: {clean_med:.3f} $\rightarrow$ {tail_med:.3f}",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=9.4,
                linespacing=1.13,
                bbox=dict(boxstyle="round,pad=0.24", facecolor="white", edgecolor="#cfd7e3", alpha=0.93),
            )
            if row_idx == 0:
                ax.set_title(col_label, fontsize=13.5, fontweight="bold", pad=8)
            if row_idx == 1:
                ax.set_xlabel(r"reco $E_T^{iso}$, $\Delta R<0.3$ [GeV]", fontsize=11.0)
            if col_idx == 0:
                ax.set_ylabel(info["label"], fontsize=11.0)
            ax.set_xlim(*EISO_DIAG_RANGE)
            ax.set_ylim(*info["ylim"])
            ax.grid(True, color=GRID, lw=0.55)
            ax.tick_params(labelsize=10.0)
        fig.text(0.064, 0.572 - row_idx * 0.310, sample_label, ha="center", va="center", rotation=90, fontsize=13.5, fontweight="bold", color=label_color)

    cax = fig.add_axes([0.902, 0.288, 0.018, 0.405])
    fig.colorbar(mesh, cax=cax, label="share per bin, log scale")
    fig.text(
        0.115,
        0.070,
        r"Orange band is the $R<0.3$ WP90 isolation cut range in 0-20%; the BDT WP80 cut is applied event-by-event using fine centrality and $E_T$ bins.",
        ha="left",
        va="center",
        fontsize=11.8,
        color=MUTED,
    )
    out = OUT_DIR / f"07_{feature_key}_vs_eiso_020_truth_inclusive_before_after_wp80.png"
    fig.savefig(out)
    plt.close(fig)
    return out, stats


def save_2d_scatter_density_diagnostic(
    arrays: dict[str, np.ndarray],
    feature_key: str,
    comparison: str = "all_truth",
) -> tuple[Path, dict[str, object]]:
    info = FEATURES[feature_key]
    feature_col = str(info["column"])
    cut_bands = isolation_wp90_bands()
    if comparison == "all_truth":
        sample_rows = [("all", "all candidates", ALL), ("truth", "truth-tagged photons", TRUTH)]
        title = f"2D width-isolation check: {info['short']} versus $R<0.3$ isolation"
        subtitle_tail = ""
        out_name = f"05_{feature_key}_vs_eiso_2d_scatter_density_clean_cut.png"
    elif comparison == "truth_inclusive":
        sample_rows = [("truth", "truth-tagged photons", TRUTH), ("inclusive", "inclusive candidates", "#1f5aa6")]
        title = f"Truth vs inclusive width-isolation phase space: {info['short']}"
        subtitle_tail = ""
        out_name = f"06_{feature_key}_vs_eiso_truth_inclusive_2d_scatter_density.png"
    else:
        raise ValueError(comparison)
    xedges = np.linspace(EISO_DIAG_RANGE[0], EISO_DIAG_RANGE[1], 82)
    yedges = np.linspace(info["ylim"][0], info["ylim"][1], 74)
    hists: list[np.ndarray] = []
    stats: dict[str, object] = {}
    for sample, _, _ in sample_rows:
        for lo, hi, cent_label in CENT_BINS:
            mask = sample_mask(arrays, sample) & (arrays["centrality"] >= lo) & (arrays["centrality"] < hi)
            hists.append(hist2d_density(arrays[EISO][mask], arrays[feature_col][mask], xedges, yedges))
    positive = np.concatenate([h[np.isfinite(h)] for h in hists if np.isfinite(h).any()])
    vmax = float(np.nanpercentile(positive, 99.75)) if positive.size else 1.0
    norm = LogNorm(vmin=1.4e-5, vmax=max(vmax, 1.5e-4))

    fig, axes = plt.subplots(2, 3, figsize=WIDE, dpi=DPI, sharex=True, sharey=True)
    fig.subplots_adjust(left=0.104, right=0.882, top=0.778, bottom=0.195, wspace=0.112, hspace=0.175)
    add_title(
        fig,
        title,
        r"Color = log density; black curve = median width. Orange band = current 90% photon-efficiency isolation cut; accepted side is left."
        + subtitle_tail,
    )
    mesh = None
    for row, (sample, sample_label, label_color) in enumerate(sample_rows):
        for col, (lo, hi, cent_label) in enumerate(CENT_BINS):
            ax = axes[row, col]
            mask = sample_mask(arrays, sample) & (arrays["centrality"] >= lo) & (arrays["centrality"] < hi)
            h = hist2d_density(arrays[EISO][mask], arrays[feature_col][mask], xedges, yedges)
            mesh = ax.pcolormesh(xedges, yedges, h.T, cmap=KBIRD, norm=norm, shading="auto")

            cut_min, cut_max, cut_mean = cut_bands[cent_label]
            ax.axvspan(cut_min, cut_max, color=ORANGE, alpha=0.16, zorder=2)
            ax.axvline(cut_mean, color=ORANGE, lw=1.5, ls="--", zorder=3)
            centers, med, _ = median_by_eiso(arrays, mask, feature_col, EISO_DIAG_MEDIAN_EDGES)
            line = ax.plot(centers, med, color=INK, lw=2.2, zorder=5)[0]
            line.set_path_effects([pe.Stroke(linewidth=4.2, foreground="white"), pe.Normal()])

            rho = pearson(arrays[EISO][mask], arrays[feature_col][mask])
            clean_med = median_width_on_side(arrays, mask, feature_col, cut_mean, True)
            tail_med = median_width_on_side(arrays, mask, feature_col, cut_mean, False)
            stats[f"{feature_key}_{sample}_{cent_label}"] = {
                "entries": int(mask.sum()),
                "pearson": rho,
                "spearman": spearman(arrays[EISO][mask], arrays[feature_col][mask]),
                "r03_wp90_cut_band_gev": [cut_min, cut_max],
                "r03_wp90_cut_mean_gev": cut_mean,
                "median_width_clean_side": clean_med,
                "median_width_busy_side": tail_med,
            }
            ax.text(
                0.035,
                0.940,
                rf"$\rho={rho:+.2f}$" + "\n" + rf"median: {clean_med:.3f} $\rightarrow$ {tail_med:.3f}",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=9.3,
                linespacing=1.18,
                bbox=dict(boxstyle="round,pad=0.24", facecolor="white", edgecolor="#cfd7e3", alpha=0.93),
            )
            if row == 0 and col == 0:
                ax.annotate(
                    "kept by isolation",
                    xy=(cut_min - 0.35, info["ylim"][1] * 0.77),
                    xytext=(cut_min - 4.6, info["ylim"][1] * 0.77),
                    arrowprops=dict(arrowstyle="<-", color=ORANGE, lw=1.4),
                    color=ORANGE,
                    fontsize=9.4,
                    fontweight="bold",
                    ha="center",
                    va="center",
                )
            ax.set_xlim(*EISO_DIAG_RANGE)
            ax.set_ylim(*info["ylim"])
            ax.grid(True, color=GRID, lw=0.55)
            ax.tick_params(labelsize=9.7)
            if row == 0:
                ax.set_title(f"{cent_label} centrality", fontsize=12.8, fontweight="bold", pad=8)
            if row == 1:
                ax.set_xlabel(r"reco $E_T^{iso}$, $\Delta R<0.3$ [GeV]", fontsize=10.8)
            if col == 0:
                ax.set_ylabel(info["label"], fontsize=10.8)
            else:
                ax.set_ylabel("")
        fig.text(0.055, 0.597 - row * 0.271, sample_label, ha="center", va="center", rotation=90, fontsize=13.3, fontweight="bold", color=label_color)
    cax = fig.add_axes([0.902, 0.300, 0.018, 0.380])
    fig.colorbar(mesh, cax=cax, label="share per bin, log scale")
    fig.text(
        0.103,
        0.060,
        (
            r"Fast read: truth photons are nearly flat; inclusive non-truth-tagged candidates sit wider "
            r"and show the stronger positive isolation-width trend."
            if comparison == "truth_inclusive"
            else r"Fast read: a strong width-isolation correlation would appear as a steep black median curve; here the truth-photon curve is nearly flat."
        ),
        ha="left",
        va="center",
        fontsize=11.9,
        color=MUTED,
    )
    out = OUT_DIR / out_name
    fig.savefig(out)
    plt.close(fig)
    return out, stats


def save_median_trends(arrays: dict[str, np.ndarray]) -> tuple[Path, dict[str, object]]:
    fig, axes = plt.subplots(2, 3, figsize=WIDE, dpi=DPI, sharex=True)
    fig.subplots_adjust(left=0.080, right=0.965, top=0.795, bottom=0.135, wspace=0.160, hspace=0.215)
    add_title(
        fig,
        "The width-isolation trend is visible in all candidates, but nearly flat for truth photons",
        r"Points are median shower width in reconstructed-isolation bins; vertical bars show the central 50% interval.",
    )
    out_stats: dict[str, object] = {}
    for row, feature_key in enumerate(["weta", "wphi"]):
        info = FEATURES[feature_key]
        feature_col = str(info["column"])
        for col, (lo, hi, cent_label) in enumerate(CENT_BINS):
            ax = axes[row, col]
            cent_mask = (arrays["centrality"] >= lo) & (arrays["centrality"] < hi)
            for sample, color, marker, label in [
                ("all", ALL, "o", "all candidates"),
                ("truth", TRUTH, "s", "truth-tagged photons"),
            ]:
                mask = cent_mask & sample_mask(arrays, sample)
                centers, med, yerr = median_by_eiso(arrays, mask, feature_col)
                ax.errorbar(
                    centers,
                    med,
                    yerr=yerr,
                    color=color,
                    marker=marker,
                    markersize=3.9,
                    linewidth=1.7,
                    capsize=2.0,
                    alpha=0.94,
                    label=label if row == 0 and col == 0 else None,
                )
                out_stats[f"{feature_key}_{sample}_{cent_label}"] = {
                    "pearson": pearson(arrays[EISO][mask], arrays[feature_col][mask]),
                    "spearman": spearman(arrays[EISO][mask], arrays[feature_col][mask]),
                    "entries": int(mask.sum()),
                    "median_by_eiso": [
                        {"center": float(c), "median": float(m) if np.isfinite(m) else None}
                        for c, m in zip(centers, med)
                    ],
                }
            ax.set_xlim(EISO_TREND_EDGES[0], EISO_TREND_EDGES[-1])
            ax.set_ylim(*info["ylim"])
            ax.grid(True, color=GRID, lw=0.60)
            ax.tick_params(labelsize=9.5)
            if row == 0:
                ax.set_title(f"{cent_label} centrality", fontsize=12.8, fontweight="bold", pad=8)
            if row == 1:
                ax.set_xlabel(r"reco $E_T^{iso}$, $\Delta R<0.3$ [GeV]", fontsize=10.8)
            if col == 0:
                ax.set_ylabel(info["label"], fontsize=10.8)
            if row == 0 and col == 0:
                ax.legend(frameon=False, fontsize=10.4, loc="upper left")
            if col == 2:
                ax.text(
                    1.035,
                    0.50,
                    info["short"],
                    transform=ax.transAxes,
                    ha="left",
                    va="center",
                    fontsize=15.0,
                    fontweight="bold",
                    color=TEAL,
                )
    out = OUT_DIR / "03_width_vs_eiso_median_trends_by_centrality.png"
    fig.savefig(out)
    plt.close(fig)
    return out, out_stats


def save_correlation_summary(arrays: dict[str, np.ndarray]) -> tuple[Path, dict[str, object]]:
    stats: dict[str, object] = {}
    fig, axes = plt.subplots(1, 2, figsize=WIDE, dpi=DPI)
    fig.subplots_adjust(left=0.135, right=0.875, top=0.780, bottom=0.235, wspace=0.320)
    add_title(
        fig,
        "All-candidate width-isolation correlation is mostly not a truth-photon effect",
        r"Pearson correlation between reconstructed isolation and shower width, split by centrality and sample definition.",
    )
    images = []
    for ax_idx, (ax, feature_key) in enumerate(zip(axes, ["weta", "wphi"])):
        info = FEATURES[feature_key]
        feature_col = str(info["column"])
        matrix = np.zeros((2, len(CENT_BINS)), dtype="float64")
        for row, sample in enumerate(["all", "truth"]):
            for col, (lo, hi, cent_label) in enumerate(CENT_BINS):
                mask = sample_mask(arrays, sample) & (arrays["centrality"] >= lo) & (arrays["centrality"] < hi)
                matrix[row, col] = pearson(arrays[EISO][mask], arrays[feature_col][mask])
                stats[f"{feature_key}_{sample}_{cent_label}"] = panel_stats(arrays, mask, feature_col)
        im = ax.imshow(matrix, cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-0.16, vcenter=0.0, vmax=0.16), aspect="auto")
        images.append(im)
        ax.set_title(info["title"], fontsize=14.4, fontweight="bold", pad=10)
        ax.set_xticks(np.arange(len(CENT_BINS)), [label for _, _, label in CENT_BINS], fontsize=11.4)
        if ax_idx == 0:
            ax.set_yticks([0, 1], ["all candidates", "truth-tagged photons"], fontsize=11.4)
        else:
            ax.set_yticks([0, 1], ["", ""], fontsize=11.4)
        for row in range(2):
            for col in range(len(CENT_BINS)):
                val = matrix[row, col]
                ax.text(
                    col,
                    row,
                    f"{val:+.2f}",
                    ha="center",
                    va="center",
                    fontsize=17.0,
                    fontweight="bold",
                    color="white" if abs(val) > 0.08 else INK,
                )
        ax.tick_params(length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
    cax = fig.add_axes([0.905, 0.315, 0.019, 0.340])
    fig.colorbar(images[-1], cax=cax, label=r"Pearson $\rho(E_T^{iso}, width)$")
    fig.text(
        0.135,
        0.122,
        r"Reading: positive values mean wider showers at larger cone energy. "
        r"All candidates show +0.07 to +0.12; truth-tagged photons stay near zero.",
        ha="left",
        va="center",
        fontsize=12.4,
        color=MUTED,
    )
    out = OUT_DIR / "04_width_isolation_correlation_summary.png"
    fig.savefig(out)
    plt.close(fig)
    return out, stats


def make_contact_sheet(paths: list[Path]) -> Path:
    thumbs: list[Image.Image] = []
    for path in paths:
        img = Image.open(path).convert("RGB")
        img.thumbnail((640, 360), Image.Resampling.LANCZOS)
        canvas = Image.new("RGB", (680, 430), "white")
        canvas.paste(img, ((680 - img.width) // 2, 18))
        draw = ImageDraw.Draw(canvas)
        try:
            font = ImageFont.truetype("/System/Library/Fonts/Supplemental/Times New Roman.ttf", 24)
        except OSError:
            font = ImageFont.load_default()
        draw.text((24, 385), path.name, fill=INK, font=font)
        thumbs.append(canvas)
    cols = 2
    rows = int(np.ceil(len(thumbs) / cols))
    sheet = Image.new("RGB", (680 * cols, 430 * rows), "white")
    for idx, thumb in enumerate(thumbs):
        sheet.paste(thumb, ((idx % cols) * 680, (idx // cols) * 430))
    out = OUT_DIR / "00_contact_sheet_width_isolation_examples.png"
    sheet.save(out)
    return out


def write_manifest(
    paths: list[Path],
    density_stats: dict[str, object],
    diag_stats: dict[str, object],
    wp80_stats: dict[str, object],
    trend_stats: dict[str, object],
    corr_stats: dict[str, object],
    arrays: dict[str, np.ndarray],
) -> Path:
    signal = arrays["is_signal"].astype(bool)
    summary = {
        "generator": str(Path(__file__).relative_to(REPO)),
        "score_cache_manifest": str(MANIFEST.relative_to(REPO)),
        "source_definition": "Local row-level score caches from model_validation_condor_20260511_194832; compares all scored candidates to truth-tagged photon rows using is_signal == 1.",
        "selection": {
            "cluster_et_range_gev": list(PT_RANGE),
            "centrality_range_percent": [0.0, 80.0],
            "isolation_column": EISO,
            "isolation_label": "reco E_T^iso, Delta R < 0.3, as used in existing Slide-21 local cache plots",
            "isolation_wp90_source": str(ISO_FINE_CSV.relative_to(REPO)),
            "score_column": SCORE,
            "bdt_wp80_source": str(BDT_WP_CSV.relative_to(REPO)),
            "features": {key: val["column"] for key, val in FEATURES.items()},
        },
        "counts": {
            "selected_entries": int(len(arrays[EISO])),
            "truth_tagged_entries": int(signal.sum()),
            "truth_tagged_fraction": float(signal.mean()),
        },
        "outputs": [str(path.relative_to(REPO)) for path in paths],
        "density_panel_stats": density_stats,
        "scatter_density_cut_stats": diag_stats,
        "wp80_before_after_stats": wp80_stats,
        "trend_stats": trend_stats,
        "correlation_stats": corr_stats,
    }
    out = OUT_DIR / "width_isolation_correlation_examples_manifest.json"
    out.write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    return out


def main() -> None:
    setup_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    arrays = load_arrays()
    paths: list[Path] = []
    density_stats: dict[str, object] = {}
    for feature_key in ["weta", "wphi"]:
        path, stats = save_density_facets(arrays, feature_key)
        paths.append(path)
        density_stats[feature_key] = stats
    diag_stats: dict[str, object] = {}
    for feature_key in ["weta", "wphi"]:
        path, stats = save_2d_scatter_density_diagnostic(arrays, feature_key, comparison="all_truth")
        paths.append(path)
        diag_stats[f"{feature_key}_all_truth"] = stats
    for feature_key in ["weta", "wphi"]:
        path, stats = save_2d_scatter_density_diagnostic(arrays, feature_key, comparison="truth_inclusive")
        paths.append(path)
        diag_stats[f"{feature_key}_truth_inclusive"] = stats
    wp80_stats: dict[str, object] = {}
    for feature_key in ["weta", "wphi"]:
        path, stats = save_020_before_after_wp80(arrays, feature_key)
        paths.append(path)
        wp80_stats[feature_key] = stats
    trend_path, trend_stats = save_median_trends(arrays)
    corr_path, corr_stats = save_correlation_summary(arrays)
    paths.extend([trend_path, corr_path])
    contact = make_contact_sheet(paths)
    manifest = write_manifest([contact, *paths], density_stats, diag_stats, wp80_stats, trend_stats, corr_stats, arrays)
    for path in [contact, *paths, manifest]:
        print(path)


if __name__ == "__main__":
    main()
