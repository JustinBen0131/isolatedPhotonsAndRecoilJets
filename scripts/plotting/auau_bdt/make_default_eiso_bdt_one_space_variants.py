#!/usr/bin/env python3
"""Make one-space BDT/isolation relationship variants for THE-11 iteration."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, LogNorm, Normalize  # noqa: E402
from PIL import Image, ImageDraw, ImageFont  # noqa: E402


REPO = Path(__file__).resolve().parents[3]
REPORT = REPO / "dataOutput/auauTightBDTValidation/model_validation_condor_20260511_194832"
MANIFEST = REPORT / "score_caches.local.list"
ISO_FINE_CSV = REPO / (
    "dataOutput/auauMLDiagnosticRuns/global_etcent_inclusive3_sixpack_20260516_135439/"
    "slideReady/isolation_wp_checks/"
    "auau_sliding_isolation_wp_r03_r04_eff70_80_90_fine7_no_numbers_coefficients.csv"
)
OUT_DIR = REPORT / "slideReady/default_eiso_vs_bdt_score_20260604/plot_variants_one_space_20260604"

SCORE = "score_ptFine_cent7"
EISO = "reco_eiso"
PT_RANGE = (15.0, 35.0)
CENT_BINS_7 = [
    (0.0, 10.0, "0-10%"),
    (10.0, 20.0, "10-20%"),
    (20.0, 30.0, "20-30%"),
    (30.0, 40.0, "30-40%"),
    (40.0, 50.0, "40-50%"),
    (50.0, 60.0, "50-60%"),
    (60.0, 80.0, "60-80%"),
]
CENT_BINS_3 = [(0.0, 20.0, "0-20%"), (20.0, 50.0, "20-50%"), (50.0, 80.0, "50-80%")]
ET_BINS = [(15.0, 20.0, "15-20"), (20.0, 25.0, "20-25"), (25.0, 30.0, "25-30"), (30.0, 35.0, "30-35")]

WIDE = (12.8, 7.2)
DPI = 200
INK = "#121826"
MUTED = "#5c6678"
GRID = "#dfe5ee"
RED = "#c94747"
BLUE = "#2869a6"
TEAL = "#148f86"
ORANGE = "#d97706"
PURPLE = "#7057a8"
GREEN = "#188a55"


def kbird() -> LinearSegmentedColormap:
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
    cmap.set_bad("#f8fafc")
    cmap.set_under("#f8fafc")
    return cmap


def setup() -> None:
    plt.rcParams.update(
        {
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "axes.facecolor": "white",
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
        data = np.load(REPO / line.strip(), allow_pickle=True)
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


def iso_thresholds_r03() -> dict[tuple[float, float], float]:
    out: dict[tuple[float, float], float] = {}
    with ISO_FINE_CSV.open() as handle:
        for row in csv.DictReader(handle):
            if row["cone"] == "R0.3" and abs(float(row["efficiency"]) - 0.90) < 1e-9:
                out[(float(row["cent_min"]), float(row["cent_max"]))] = float(row["threshold_gev"])
    return out


def add_working_margins(data: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    thresholds = iso_thresholds_r03()
    iso_cut = np.full(len(data["score"]), np.nan, dtype="float32")
    bdt_cut = np.full(len(data["score"]), np.nan, dtype="float32")
    signal = data["is_signal"].astype(bool)
    for lo, hi, _ in CENT_BINS_7:
        mask = (data["cent"] >= lo) & (data["cent"] < hi)
        iso_cut[mask] = thresholds[(lo, hi)]
        sig_mask = mask & signal
        bdt_cut[mask] = np.nanpercentile(data["score"][sig_mask], 20.0)
    valid = np.isfinite(iso_cut) & np.isfinite(bdt_cut)
    out = {key: val[valid] for key, val in data.items()}
    out["iso_cut"] = iso_cut[valid]
    out["bdt_cut"] = bdt_cut[valid]
    out["iso_margin"] = out["iso_cut"] - out["eiso"]
    out["bdt_margin"] = out["score"] - out["bdt_cut"]
    return out


def hist2d(x: np.ndarray, y: np.ndarray, xedges: np.ndarray, yedges: np.ndarray) -> np.ndarray:
    hist, _, _ = np.histogram2d(x, y, bins=[xedges, yedges])
    return hist.astype("float64")


def smooth(hist: np.ndarray, passes: int = 2) -> np.ndarray:
    kernel = np.array([1.0, 2.0, 1.0], dtype="float64")
    kernel /= kernel.sum()
    out = hist.astype("float64", copy=True)
    for _ in range(passes):
        out = np.apply_along_axis(lambda v: np.convolve(v, kernel, mode="same"), 0, out)
        out = np.apply_along_axis(lambda v: np.convolve(v, kernel, mode="same"), 1, out)
    return out


def corr(x: np.ndarray, y: np.ndarray) -> float:
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3:
        return float("nan")
    return float(np.corrcoef(x[ok].astype("float64"), y[ok].astype("float64"))[0, 1])


def binned_line(x: np.ndarray, y: np.ndarray, bins: np.ndarray, min_n: int = 80) -> tuple[np.ndarray, np.ndarray]:
    centers: list[float] = []
    meds: list[float] = []
    for lo, hi in zip(bins[:-1], bins[1:]):
        mask = (x >= lo) & (x < hi)
        if mask.sum() >= min_n:
            centers.append(float(0.5 * (lo + hi)))
            meds.append(float(np.nanmedian(y[mask])))
    return np.asarray(centers), np.asarray(meds)


def draw_class_contours(ax, data: dict[str, np.ndarray], xkey: str, ykey: str, xedges, yedges) -> None:
    signal = data["is_signal"].astype(bool)
    for mask, color, label in [(signal, RED, "truth photons"), (~signal, BLUE, "inclusive jets")]:
        hist = smooth(hist2d(data[xkey][mask], data[ykey][mask], xedges, yedges), passes=2)
        positive = hist[hist > 0]
        if positive.size < 10:
            continue
        levels = np.unique(np.nanpercentile(positive, [70, 86, 95]))
        levels = levels[levels > 0]
        xx = 0.5 * (xedges[:-1] + xedges[1:])
        yy = 0.5 * (yedges[:-1] + yedges[1:])
        ax.contour(xx, yy, hist.T, levels=levels, colors=color, linewidths=1.9, alpha=0.92)
        ax.plot([], [], color=color, lw=2.5, label=label)


def fit_label(ax, data: dict[str, np.ndarray], xkey: str, ykey: str, xlim: tuple[float, float], loc="raw") -> dict[str, float]:
    signal = data["is_signal"].astype(bool)
    stats = {}
    for mask, color, name, yoff in [
        (signal, RED, "truth photons", 0.0),
        (~signal, BLUE, "inclusive jets", -0.045),
    ]:
        x = data[xkey][mask]
        y = data[ykey][mask]
        ok = np.isfinite(x) & np.isfinite(y) & (x >= xlim[0]) & (x <= xlim[1])
        if ok.sum() < 20:
            continue
        coef = np.polyfit(x[ok].astype("float64"), y[ok].astype("float64"), 1)
        xs = np.linspace(xlim[0], xlim[1], 80)
        ax.plot(xs, coef[0] * xs + coef[1], color=color, ls="--", lw=2.0, alpha=0.95)
        rho = corr(x[ok], y[ok])
        stats[name] = rho
        if loc == "raw":
            ax.text(
                0.03,
                0.94 + yoff,
                rf"{name}: fit slope {coef[0]:+.3f}, $\rho={rho:+.2f}$",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=11.2,
                color=color,
                bbox={"facecolor": "white", "edgecolor": "#d7dee8", "boxstyle": "round,pad=0.24", "alpha": 0.92},
            )
    return stats


def save_probability_map(data: dict[str, np.ndarray]) -> tuple[Path, dict]:
    signal = data["is_signal"].astype(bool)
    mask = (data["cent"] >= 0.0) & (data["cent"] < 20.0)
    x = data["eiso"][mask]
    y = data["score"][mask]
    sig = signal[mask]
    xlim = tuple(np.nanpercentile(x, [0.5, 99.2]))
    xpad = 0.04 * (xlim[1] - xlim[0])
    xlim = (float(xlim[0] - xpad), float(xlim[1] + xpad))
    xedges = np.linspace(xlim[0], xlim[1], 76)
    yedges = np.linspace(0.0, 1.0, 66)
    counts = hist2d(x, y, xedges, yedges)
    sig_counts = hist2d(x[sig], y[sig], xedges, yedges)
    prob = np.divide(sig_counts, counts, out=np.full_like(counts, np.nan), where=counts >= 8)

    fig, ax = plt.subplots(figsize=WIDE, dpi=DPI)
    fig.subplots_adjust(left=0.085, right=0.875, top=0.835, bottom=0.140)
    fig.text(0.055, 0.955, r"One-look BDT-isolation map in 0-20% centrality", fontsize=23.5, fontweight="bold", ha="left", va="top")
    fig.text(
        0.055,
        0.905,
        r"Color = truth-photon fraction; black contours = total density; orange/green lines = isolation and BDT working cuts.",
        fontsize=12.8,
        color=MUTED,
        ha="left",
        va="top",
    )
    mesh = ax.pcolormesh(xedges, yedges, prob.T, cmap="RdYlBu_r", norm=Normalize(vmin=0, vmax=1), shading="auto")
    positive = counts[counts > 0]
    if positive.size:
        levels = np.unique(np.nanpercentile(positive, [55, 75, 90, 97]))
        xx = 0.5 * (xedges[:-1] + xedges[1:])
        yy = 0.5 * (yedges[:-1] + yedges[1:])
        ax.contour(xx, yy, counts.T, levels=levels, colors="black", linewidths=0.8, alpha=0.62)
    sub = {key: val[mask] for key, val in data.items()}
    draw_class_contours(ax, sub, "eiso", "score", xedges, yedges)
    stats = {
        "pearson_all": corr(sub["eiso"], sub["score"]),
        "pearson_signal": corr(sub["eiso"][sub["is_signal"].astype(bool)], sub["score"][sub["is_signal"].astype(bool)]),
        "pearson_background": corr(sub["eiso"][~sub["is_signal"].astype(bool)], sub["score"][~sub["is_signal"].astype(bool)]),
    }

    thresholds = iso_thresholds_r03()
    iso_band = [thresholds[(0.0, 10.0)], thresholds[(10.0, 20.0)]]
    bdt_wp80 = float(np.nanpercentile(data["score"][mask & signal], 20.0))
    ax.axvspan(min(iso_band), max(iso_band), color=ORANGE, alpha=0.16, lw=0)
    ax.axvline(float(np.mean(iso_band)), color=ORANGE, lw=2.2, ls="-")
    ax.axhline(bdt_wp80, color=GREEN, lw=2.2, ls="-")
    ax.text(float(np.mean(iso_band)), 0.035, "iso WP90", color=ORANGE, fontsize=11.5, ha="center", va="bottom", fontweight="bold")
    ax.text(xlim[1] - 0.05 * (xlim[1] - xlim[0]), bdt_wp80 + 0.020, "BDT WP80", color=GREEN, fontsize=11.5, ha="right", va="bottom", fontweight="bold")
    ax.text(
        0.700,
        0.910,
        "Read in one pass:\n"
        "red upper band = truth-rich\n"
        "blue contour = inclusive tail\n"
        "ET medians stay near the same ridge",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=12.4,
        color=INK,
        linespacing=1.25,
        bbox={"facecolor": "white", "edgecolor": "#d7dee8", "boxstyle": "round,pad=0.35", "alpha": 0.94},
    )

    et_colors = ["#7c3aed", "#0f766e", "#d97706", "#be123c"]
    for (lo, hi, label), color in zip(ET_BINS, et_colors):
        em = mask & (data["et"] >= lo) & (data["et"] < hi)
        if em.sum() < 20:
            continue
        ax.scatter(
            [np.nanmedian(data["eiso"][em])],
            [np.nanmedian(data["score"][em])],
            s=105,
            marker="o",
            facecolor="white",
            edgecolor=color,
            linewidth=2.5,
            zorder=5,
            label=f"{label} GeV median",
        )
    ax.set_xlim(*xlim)
    ax.set_ylim(0, 1)
    ax.set_xlabel(r"default reco $E_T^{iso}$ [GeV]", fontsize=14.5)
    ax.set_ylabel("BDT score", fontsize=14.5)
    ax.tick_params(labelsize=12.2)
    ax.grid(True, color=GRID, lw=0.7)
    ax.legend(loc="lower left", ncol=2, frameon=True, facecolor="white", edgecolor="#d7dee8", fontsize=9.8)
    cax = fig.add_axes([0.900, 0.200, 0.018, 0.590])
    fig.colorbar(mesh, cax=cax, label="truth-photon fraction")
    out = OUT_DIR / "01_probability_density_020_with_cuts_and_fits.png"
    fig.savefig(out)
    plt.close(fig)
    return out, {
        "centrality": "0-20%",
        "bdt_wp80_local": bdt_wp80,
        "iso_wp90_band_r03": iso_band,
        "fit_correlations": stats,
    }


def save_margin_map(data: dict[str, np.ndarray]) -> tuple[Path, dict]:
    x = data["iso_margin"]
    y = data["bdt_margin"]
    xlim = (-7.5, 9.5)
    ylim = (-0.42, 0.36)
    xedges = np.linspace(*xlim, 86)
    yedges = np.linspace(*ylim, 70)
    counts = hist2d(x, y, xedges, yedges)
    fig, ax = plt.subplots(figsize=WIDE, dpi=DPI)
    fig.subplots_adjust(left=0.090, right=0.880, top=0.835, bottom=0.145)
    fig.text(0.055, 0.955, "Cut-margin density: both decisions in one coordinate system", fontsize=23.0, fontweight="bold", ha="left", va="top")
    fig.text(
        0.055,
        0.905,
        r"$x=I_{90}(c)-E_T^{iso}$, $y=score-T_{80}(c)$. Zero lines are the isolation and BDT pass boundaries.",
        fontsize=12.8,
        color=MUTED,
        ha="left",
        va="top",
    )
    mesh = ax.pcolormesh(xedges, yedges, counts.T, cmap=kbird(), norm=LogNorm(vmin=1, vmax=max(2, float(np.nanmax(counts)))), shading="auto")
    draw_class_contours(ax, data, "iso_margin", "bdt_margin", xedges, yedges)
    stats = fit_label(ax, data, "iso_margin", "bdt_margin", xlim, loc="none")
    colors = [PURPLE, TEAL, ORANGE]
    bins = np.linspace(xlim[0], xlim[1], 24)
    for (lo, hi, label), color in zip(CENT_BINS_3, colors):
        cmask = (data["cent"] >= lo) & (data["cent"] < hi)
        xs, ys = binned_line(data["iso_margin"][cmask], data["bdt_margin"][cmask], bins, min_n=60)
        if len(xs):
            ax.plot(xs, ys, color=color, lw=3.0, marker="o", ms=4, label=f"{label} median path")
    ax.axvline(0, color=ORANGE, lw=2.4)
    ax.axhline(0, color=GREEN, lw=2.4)
    ax.text(0.02, 0.96, "fails isolation", transform=ax.transAxes, fontsize=11.5, color=ORANGE, ha="left", va="top", fontweight="bold")
    ax.text(0.98, 0.96, "passes isolation", transform=ax.transAxes, fontsize=11.5, color=ORANGE, ha="right", va="top", fontweight="bold")
    ax.text(0.985, 0.525, "BDT pass", transform=ax.transAxes, fontsize=11.5, color=GREEN, ha="right", va="bottom", fontweight="bold")
    ax.text(0.985, 0.485, "BDT fail", transform=ax.transAxes, fontsize=11.5, color=GREEN, ha="right", va="top", fontweight="bold")
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xlabel(r"isolation margin: $I_{90}(c)-$ default reco $E_T^{iso}$ [GeV]", fontsize=14.2)
    ax.set_ylabel(r"BDT margin: score $-T_{80}(c)$", fontsize=14.2)
    ax.tick_params(labelsize=12)
    ax.grid(True, color=GRID, lw=0.65)
    ax.legend(loc="lower right", frameon=True, facecolor="white", edgecolor="#d7dee8", fontsize=10.4)
    cax = fig.add_axes([0.905, 0.210, 0.018, 0.560])
    fig.colorbar(mesh, cax=cax, label="candidates / bin, log")
    out = OUT_DIR / "02_cut_margin_density_all_centrality.png"
    fig.savefig(out)
    plt.close(fig)
    return out, {
        "pearson_margin_all": corr(x, y),
        "pearson_margin_signal": corr(x[data["is_signal"].astype(bool)], y[data["is_signal"].astype(bool)]),
        "pearson_margin_background": corr(x[~data["is_signal"].astype(bool)], y[~data["is_signal"].astype(bool)]),
        "fit_correlations": stats,
    }


def save_margin_et_map(data: dict[str, np.ndarray]) -> tuple[Path, dict]:
    cmask = (data["cent"] >= 0.0) & (data["cent"] < 20.0)
    sub = {key: val[cmask] for key, val in data.items()}
    xlim = (-7.5, 9.5)
    ylim = (-0.42, 0.36)
    xedges = np.linspace(*xlim, 86)
    yedges = np.linspace(*ylim, 70)
    counts = hist2d(sub["iso_margin"], sub["bdt_margin"], xedges, yedges)
    fig, ax = plt.subplots(figsize=WIDE, dpi=DPI)
    fig.subplots_adjust(left=0.090, right=0.880, top=0.835, bottom=0.145)
    fig.text(0.055, 0.955, r"0-20% centrality: $E_T$ medians move inside the same cut-margin space", fontsize=23.0, fontweight="bold", ha="left", va="top")
    fig.text(
        0.055,
        0.905,
        r"Density is all candidates. Colored points are $E_T$-bin medians; contours separate truth photons from inclusive jets.",
        fontsize=12.8,
        color=MUTED,
        ha="left",
        va="top",
    )
    mesh = ax.pcolormesh(xedges, yedges, counts.T, cmap=kbird(), norm=LogNorm(vmin=1, vmax=max(2, float(np.nanmax(counts)))), shading="auto")
    draw_class_contours(ax, sub, "iso_margin", "bdt_margin", xedges, yedges)
    et_colors = ["#7c3aed", "#0f766e", "#d97706", "#be123c"]
    medians = []
    for (lo, hi, label), color in zip(ET_BINS, et_colors):
        mask = (sub["et"] >= lo) & (sub["et"] < hi)
        if mask.sum() < 30:
            continue
        xm = float(np.nanmedian(sub["iso_margin"][mask]))
        ym = float(np.nanmedian(sub["bdt_margin"][mask]))
        medians.append({"et_bin": label, "iso_margin_median": xm, "bdt_margin_median": ym, "entries": int(mask.sum())})
        ax.scatter([xm], [ym], s=170, facecolor="white", edgecolor=color, linewidth=3.0, zorder=6)
        ax.text(xm, ym, label.replace("-", "\n"), color=color, ha="center", va="center", fontsize=8.8, fontweight="bold", zorder=7)
    if len(medians) >= 2:
        ax.plot([m["iso_margin_median"] for m in medians], [m["bdt_margin_median"] for m in medians], color=INK, lw=2.0, alpha=0.70, zorder=5)
    ax.axvline(0, color=ORANGE, lw=2.4)
    ax.axhline(0, color=GREEN, lw=2.4)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xlabel(r"isolation margin: $I_{90}(c)-$ default reco $E_T^{iso}$ [GeV]", fontsize=14.2)
    ax.set_ylabel(r"BDT margin: score $-T_{80}(c)$", fontsize=14.2)
    ax.tick_params(labelsize=12)
    ax.grid(True, color=GRID, lw=0.65)
    ax.legend(loc="lower right", frameon=True, facecolor="white", edgecolor="#d7dee8", fontsize=10.6)
    cax = fig.add_axes([0.905, 0.210, 0.018, 0.560])
    fig.colorbar(mesh, cax=cax, label="candidates / bin, log")
    out = OUT_DIR / "03_cut_margin_density_020_et_medians.png"
    fig.savefig(out)
    plt.close(fig)
    return out, {"et_medians": medians}


def make_contact_sheet(paths: list[Path]) -> Path:
    thumbs = []
    for path in paths:
        img = Image.open(path).convert("RGB")
        img.thumbnail((780, 440), Image.Resampling.LANCZOS)
        thumbs.append((path, img.copy()))
    w, h = 820, 500
    sheet = Image.new("RGB", (w * len(thumbs), h), "white")
    draw = ImageDraw.Draw(sheet)
    try:
        font = ImageFont.truetype("Arial.ttf", 24)
    except OSError:
        font = ImageFont.load_default()
    for idx, (path, img) in enumerate(thumbs):
        x = idx * w + 20
        sheet.paste(img, (x, 42))
        draw.text((x, 12), path.name, fill=(20, 24, 34), font=font)
    out = OUT_DIR / "00_contact_sheet_one_space_variants.png"
    sheet.save(out)
    return out


def main() -> None:
    setup()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    data = add_working_margins(load_arrays())
    outputs = []
    p1, s1 = save_probability_map(data)
    outputs.append(p1)
    p2, s2 = save_margin_map(data)
    outputs.append(p2)
    p3, s3 = save_margin_et_map(data)
    outputs.append(p3)
    contact = make_contact_sheet(outputs)
    manifest = {
        "schema": "THE11_ONE_SPACE_ISO_BDT_VARIANTS_V1",
        "source": str(MANIFEST),
        "score_column": SCORE,
        "isolation_column": EISO,
        "notes": [
            "These variants use local row-level default reco_eiso caches; local row-level reco_eiso_r40 density is not available.",
            "The cut-margin plots use R=0.3 WP90 isolation constants because reco_eiso is the available default row-level isolation variable.",
            "BDT WP80 margins are computed from the same local score cache by centrality bin at 80% truth-photon efficiency.",
            "Google Slides was not mutated.",
        ],
        "outputs": {
            "contact_sheet": str(contact),
            "probability_density_020": str(p1),
            "cut_margin_all_centrality": str(p2),
            "cut_margin_020_et_medians": str(p3),
        },
        "stats": {
            "probability_density_020": s1,
            "cut_margin_all_centrality": s2,
            "cut_margin_020_et_medians": s3,
        },
    }
    manifest_path = OUT_DIR / "one_space_variants_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
