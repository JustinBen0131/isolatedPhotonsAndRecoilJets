#!/usr/bin/env python3
"""Make a clean/busy isolation grid for BDT score versus cluster ET."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import Normalize  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from PIL import Image, ImageDraw, ImageFont  # noqa: E402


REPO = Path(__file__).resolve().parents[3]
REPORT = REPO / "dataOutput/auauTightBDTValidation/model_validation_condor_20260511_194832"
MANIFEST = REPORT / "score_caches.local.list"
ISO_FINE_CSV = REPO / (
    "dataOutput/auauMLDiagnosticRuns/global_etcent_inclusive3_sixpack_20260516_135439/"
    "slideReady/isolation_wp_checks/"
    "auau_sliding_isolation_wp_r03_r04_eff70_80_90_fine7_no_numbers_coefficients.csv"
)
OUT_DIR = REPORT / "slideReady/default_eiso_vs_bdt_score_20260604/plot_variants_clean_busy_grid_20260604"

SCORE = "score_ptFine_cent7"
EISO = "reco_eiso"
CENT_ROWS = [(0.0, 10.0, "0-10%"), (10.0, 20.0, "10-20%")]
ET_BINS = [
    (15.0, 17.0),
    (17.0, 19.0),
    (19.0, 21.0),
    (21.0, 23.0),
    (23.0, 25.0),
    (25.0, 27.0),
    (27.0, 30.0),
    (30.0, 35.0),
]
PT_RANGE = (15.0, 35.0)

W, H, DPI = 2560, 1440, 200
INK = "#121826"
MUTED = "#5b6472"
GRID = "#dfe5ee"
GREEN = "#14804a"
ORANGE = "#d97706"
RED = "#c94747"
BLUE = "#2869a6"


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
        & (arrays["cent"] < 20.0)
    )
    return {key: val[selected] for key, val in arrays.items()}


def iso_thresholds_r03() -> dict[tuple[float, float], float]:
    out: dict[tuple[float, float], float] = {}
    with ISO_FINE_CSV.open() as handle:
        for row in csv.DictReader(handle):
            if row["cone"] == "R0.3" and abs(float(row["efficiency"]) - 0.90) < 1e-9:
                out[(float(row["cent_min"]), float(row["cent_max"]))] = float(row["threshold_gev"])
    return out


def wp80_curve(data: dict[str, np.ndarray], cent_range: tuple[float, float]) -> list[tuple[float, float]]:
    signal = data["is_signal"].astype(bool)
    out = []
    clo, chi = cent_range
    for elo, ehi in ET_BINS:
        mask = signal & (data["cent"] >= clo) & (data["cent"] < chi) & (data["et"] >= elo) & (data["et"] < ehi)
        if mask.sum() < 50:
            continue
        out.append((0.5 * (elo + ehi), float(np.nanpercentile(data["score"][mask], 20.0)), int(mask.sum())))
    return out


def histograms(data: dict[str, np.ndarray], mask: np.ndarray, xedges, yedges):
    signal = data["is_signal"].astype(bool)
    total, _, _ = np.histogram2d(data["et"][mask], data["score"][mask], bins=[xedges, yedges])
    sig, _, _ = np.histogram2d(data["et"][mask & signal], data["score"][mask & signal], bins=[xedges, yedges])
    frac = np.divide(sig, total, out=np.full_like(total, np.nan, dtype="float64"), where=total >= 8)
    return total, frac


def add_density_contours(ax, total, xedges, yedges) -> None:
    positive = total[total > 0]
    if positive.size < 8:
        return
    levels = np.unique(np.nanpercentile(positive, [62, 82, 94]))
    xx = 0.5 * (xedges[:-1] + xedges[1:])
    yy = 0.5 * (yedges[:-1] + yedges[1:])
    ax.contour(xx, yy, total.T, levels=levels, colors="black", linewidths=0.8, alpha=0.55)


def smooth(hist: np.ndarray, passes: int = 2) -> np.ndarray:
    kernel = np.array([1.0, 2.0, 1.0], dtype="float64")
    kernel /= kernel.sum()
    out = hist.astype("float64", copy=True)
    for _ in range(passes):
        out = np.apply_along_axis(lambda v: np.convolve(v, kernel, mode="same"), 0, out)
        out = np.apply_along_axis(lambda v: np.convolve(v, kernel, mode="same"), 1, out)
    return out


def add_class_contours(ax, data: dict[str, np.ndarray], mask: np.ndarray, xedges, yedges) -> None:
    signal = data["is_signal"].astype(bool)
    xx = 0.5 * (xedges[:-1] + xedges[1:])
    yy = 0.5 * (yedges[:-1] + yedges[1:])
    for class_mask, color in [(signal, RED), (~signal, BLUE)]:
        hist, _, _ = np.histogram2d(data["et"][mask & class_mask], data["score"][mask & class_mask], bins=[xedges, yedges])
        hist = smooth(hist, passes=2)
        if hist.sum() <= 0:
            continue
        hist = hist / hist.sum()
        positive = hist[hist > 0]
        if positive.size < 8:
            continue
        levels = np.unique(np.nanpercentile(positive, [72, 88, 96]))
        ax.contour(xx, yy, hist.T, levels=levels, colors=color, linewidths=2.0, alpha=0.95)


def panel_stats(data: dict[str, np.ndarray], mask: np.ndarray) -> dict[str, float | int]:
    signal = data["is_signal"].astype(bool)
    n = int(mask.sum())
    if n == 0:
        return {"entries": 0, "truth_fraction": float("nan"), "median_score": float("nan")}
    return {
        "entries": n,
        "truth_fraction": float(signal[mask].mean()),
        "median_score": float(np.nanmedian(data["score"][mask])),
    }


def draw_clean_busy_grid(data: dict[str, np.ndarray]) -> tuple[Path, dict]:
    thresholds = iso_thresholds_r03()
    xedges = np.linspace(PT_RANGE[0], PT_RANGE[1], 76)
    yedges = np.linspace(0.0, 1.0, 72)

    fig, axes = plt.subplots(2, 2, figsize=(W / DPI, H / DPI), dpi=DPI, sharex=True, sharey=True)
    fig.subplots_adjust(left=0.080, right=0.875, top=0.835, bottom=0.125, wspace=0.075, hspace=0.135)
    fig.text(0.050, 0.955, r"Split the cone first: the BDT score map becomes readable", fontsize=25.2, fontweight="bold", ha="left", va="top")
    fig.text(
        0.050,
        0.905,
        r"Rows are centrality; columns are default isolation state. Every panel uses cluster $E_T$ vs BDT score; green curves are the local WP80 threshold.",
        fontsize=13.8,
        color=MUTED,
        ha="left",
        va="top",
    )

    stats: dict[str, dict] = {}
    mesh = None
    for row, (clo, chi, cent_label) in enumerate(CENT_ROWS):
        iso_cut = thresholds[(clo, chi)]
        cent_mask = (data["cent"] >= clo) & (data["cent"] < chi)
        states = [
            ("clean cone", data["eiso"] < iso_cut, rf"passes isolation: $E_T^{{iso}}<{iso_cut:.2f}$ GeV"),
            ("busy cone", data["eiso"] >= iso_cut, rf"fails isolation: $E_T^{{iso}}\geq{iso_cut:.2f}$ GeV"),
        ]
        curve = np.asarray(wp80_curve(data, (clo, chi)), dtype=float)
        for col, (state_label, state_mask, state_subtitle) in enumerate(states):
            ax = axes[row, col]
            mask = cent_mask & state_mask
            total, frac = histograms(data, mask, xedges, yedges)
            mesh = ax.pcolormesh(xedges, yedges, frac.T, cmap="RdYlBu_r", norm=Normalize(0, 1), shading="auto")
            add_density_contours(ax, total, xedges, yedges)
            if curve.size:
                ax.plot(curve[:, 0], curve[:, 1], color=GREEN, lw=3.2)
            ax.set_xlim(*PT_RANGE)
            ax.set_ylim(0.0, 1.0)
            ax.grid(True, color=GRID, lw=0.7)
            ax.tick_params(labelsize=11.4)
            if row == 0:
                ax.set_title(state_label, fontsize=18.0, fontweight="bold", pad=18)
                ax.text(
                    0.5,
                    1.018,
                    state_subtitle,
                    transform=ax.transAxes,
                    ha="center",
                    va="bottom",
                    fontsize=10.9,
                    color=ORANGE,
                    fontweight="bold",
                )
            if col == 0:
                ax.set_ylabel(f"{cent_label}\nBDT score", fontsize=14.2, fontweight="bold")
            if row == 1:
                ax.set_xlabel(r"cluster $E_T$ [GeV]", fontsize=14.0)
            ax.text(
                0.030,
                0.935,
                rf"$f_\gamma={panel_stats(data, mask)['truth_fraction']:.2f}$",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=12.4,
                color=INK,
                bbox={"facecolor": "white", "edgecolor": "#d7dee8", "boxstyle": "round,pad=0.25", "alpha": 0.90},
            )
            stats[f"{cent_label}_{state_label.replace(' ', '_')}"] = {
                **panel_stats(data, mask),
                "isolation_threshold_gev": iso_cut,
            }

    cax = fig.add_axes([0.902, 0.190, 0.018, 0.585])
    fig.colorbar(mesh, cax=cax, label="truth-photon fraction")
    fig.text(0.770, 0.112, "black contours = candidate density", fontsize=11.8, color=MUTED, ha="right", va="center")
    fig.text(0.882, 0.112, "green curve = BDT WP80", fontsize=11.8, color=GREEN, ha="right", va="center", fontweight="bold")

    out = OUT_DIR / "01_clean_busy_cone_et_score_truth_fraction_grid.png"
    fig.savefig(out)
    plt.close(fig)
    return out, stats


def draw_minimal_binary_grid(data: dict[str, np.ndarray]) -> tuple[Path, dict]:
    """A simpler occupancy-only version: where do passing candidates live?"""
    thresholds = iso_thresholds_r03()
    xedges = np.linspace(PT_RANGE[0], PT_RANGE[1], 52)
    yedges = np.linspace(0.0, 1.0, 52)
    fig, axes = plt.subplots(2, 2, figsize=(W / DPI, H / DPI), dpi=DPI, sharex=True, sharey=True)
    fig.subplots_adjust(left=0.080, right=0.875, top=0.835, bottom=0.125, wspace=0.075, hspace=0.135)
    fig.text(0.050, 0.955, r"Clean cone versus busy cone: where the BDT keeps photons", fontsize=25.2, fontweight="bold", ha="left", va="top")
    fig.text(0.050, 0.905, r"Same axes in every panel. Color is candidate density; green curve is the local WP80 score cut.", fontsize=13.8, color=MUTED, ha="left", va="top")
    vmax = 1.0
    hists = []
    masks = []
    labels = []
    for clo, chi, cent_label in CENT_ROWS:
        iso_cut = thresholds[(clo, chi)]
        cent_mask = (data["cent"] >= clo) & (data["cent"] < chi)
        for state_label, state_mask in [("clean cone", data["eiso"] < iso_cut), ("busy cone", data["eiso"] >= iso_cut)]:
            mask = cent_mask & state_mask
            hist, _, _ = np.histogram2d(data["et"][mask], data["score"][mask], bins=[xedges, yedges])
            hists.append(hist)
            masks.append((mask, clo, chi))
            labels.append((cent_label, state_label, iso_cut))
            vmax = max(vmax, float(np.nanmax(hist)))
    mesh = None
    for ax, hist, (mask, clo, chi), (cent_label, state_label, iso_cut) in zip(axes.ravel(), hists, masks, labels):
        mesh = ax.pcolormesh(xedges, yedges, hist.T, cmap="magma", shading="auto", vmin=0, vmax=vmax)
        curve = np.asarray(wp80_curve(data, (clo, chi)), dtype=float)
        if curve.size:
            ax.plot(curve[:, 0], curve[:, 1], color=GREEN, lw=3.2)
        ax.set_title(f"{cent_label} | {state_label}", fontsize=16.0, fontweight="bold", pad=8)
        ax.set_xlim(*PT_RANGE)
        ax.set_ylim(0, 1)
        ax.grid(True, color=GRID, lw=0.7)
    axes[0, 0].set_ylabel("BDT score", fontsize=14.0)
    axes[1, 0].set_ylabel("BDT score", fontsize=14.0)
    axes[1, 0].set_xlabel(r"cluster $E_T$ [GeV]", fontsize=14.0)
    axes[1, 1].set_xlabel(r"cluster $E_T$ [GeV]", fontsize=14.0)
    cax = fig.add_axes([0.902, 0.190, 0.018, 0.585])
    fig.colorbar(mesh, cax=cax, label="candidates / bin")
    out = OUT_DIR / "02_clean_busy_cone_et_score_density_grid.png"
    fig.savefig(out)
    plt.close(fig)
    return out, {"vmax": vmax}


def draw_class_contour_grid(data: dict[str, np.ndarray]) -> tuple[Path, dict]:
    thresholds = iso_thresholds_r03()
    xedges = np.linspace(PT_RANGE[0], PT_RANGE[1], 74)
    yedges = np.linspace(0.0, 1.0, 70)
    fig, axes = plt.subplots(2, 2, figsize=(W / DPI, H / DPI), dpi=DPI, sharex=True, sharey=True)
    fig.subplots_adjust(left=0.082, right=0.955, top=0.770, bottom=0.125, wspace=0.075, hspace=0.125)
    fig.text(0.050, 0.955, r"Clean cone vs busy cone: the BDT separation becomes visible", fontsize=25.2, fontweight="bold", ha="left", va="top")
    fig.text(
        0.050,
        0.905,
        r"Clean cone = below the isolation cut; busy cone = above it. Red/blue shapes show where each class lives; green is the $E_T$-dependent BDT cut.",
        fontsize=13.6,
        color=MUTED,
        ha="left",
        va="top",
    )
    handles = [
        Line2D([0], [0], color=RED, lw=2.5, label="truth photons"),
        Line2D([0], [0], color=BLUE, lw=2.5, label="inclusive jets"),
        Line2D([0], [0], color=GREEN, lw=3.0, label="BDT WP80 curve"),
        Patch(facecolor="#999999", alpha=0.28, edgecolor="none", label="candidate density"),
    ]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.540, 0.865), ncol=4, frameon=False, fontsize=12.7)
    fig.text(0.300, 0.800, "clean cone", fontsize=18.5, fontweight="bold", ha="center", va="center")
    fig.text(0.725, 0.800, "busy cone", fontsize=18.5, fontweight="bold", ha="center", va="center")
    stats: dict[str, dict] = {}
    vmax = 1.0
    panel_payload = []
    for clo, chi, cent_label in CENT_ROWS:
        iso_cut = thresholds[(clo, chi)]
        cent_mask = (data["cent"] >= clo) & (data["cent"] < chi)
        for state_label, state_mask in [("clean", data["eiso"] < iso_cut), ("busy", data["eiso"] >= iso_cut)]:
            mask = cent_mask & state_mask
            total, _, _ = np.histogram2d(data["et"][mask], data["score"][mask], bins=[xedges, yedges])
            vmax = max(vmax, float(np.nanpercentile(total[total > 0], 99)) if np.any(total > 0) else 1.0)
            panel_payload.append((mask, total, clo, chi, cent_label, state_label, iso_cut))
    mesh = None
    for ax, (mask, total, clo, chi, cent_label, state_label, iso_cut) in zip(axes.ravel(), panel_payload):
        mesh = ax.pcolormesh(xedges, yedges, total.T, cmap="Greys", shading="auto", vmin=0, vmax=vmax, alpha=0.48)
        add_class_contours(ax, data, mask, xedges, yedges)
        curve = np.asarray(wp80_curve(data, (clo, chi)), dtype=float)
        if curve.size:
            ax.plot(curve[:, 0], curve[:, 1], color=GREEN, lw=3.0, zorder=5)
        ax.set_xlim(*PT_RANGE)
        ax.set_ylim(0.0, 1.0)
        ax.grid(True, color=GRID, lw=0.65)
        ax.tick_params(labelsize=11.5)
        if state_label == "clean":
            ax.set_ylabel(f"{cent_label}\nBDT score", fontsize=14.0, fontweight="bold")
        stats[f"{cent_label}_{state_label}"] = {
            **panel_stats(data, mask),
            "isolation_threshold_gev": iso_cut,
        }
    fig.text(0.515, 0.055, r"cluster $E_T$ [GeV]", fontsize=15.5, ha="center", va="center")
    out = OUT_DIR / "03_clean_busy_cone_class_contours_grid.png"
    fig.savefig(out)
    plt.close(fig)
    return out, stats


def make_contact(paths: list[Path]) -> Path:
    thumbs = []
    for path in paths:
        img = Image.open(path).convert("RGB")
        img.thumbnail((980, 550), Image.Resampling.LANCZOS)
        thumbs.append((path, img.copy()))
    sheet = Image.new("RGB", (1020 * len(thumbs), 620), "white")
    draw = ImageDraw.Draw(sheet)
    try:
        font = ImageFont.truetype("Arial.ttf", 28)
    except OSError:
        font = ImageFont.load_default()
    for idx, (path, img) in enumerate(thumbs):
        x = idx * 1020 + 20
        sheet.paste(img, (x, 52))
        draw.text((x, 14), path.name, fill=(20, 24, 34), font=font)
    out = OUT_DIR / "00_contact_sheet_clean_busy_grid.png"
    sheet.save(out)
    return out


def main() -> None:
    setup()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    data = load_arrays()
    p1, s1 = draw_clean_busy_grid(data)
    p2, s2 = draw_minimal_binary_grid(data)
    p3, s3 = draw_class_contour_grid(data)
    contact = make_contact([p1, p2, p3])
    manifest = {
        "schema": "THE11_CLEAN_BUSY_CONE_GRID_V1",
        "outputs": {
            "contact_sheet": str(contact),
            "truth_fraction_grid": str(p1),
            "density_grid": str(p2),
            "class_contour_grid": str(p3),
        },
        "sources": {
            "row_level_cache_manifest": str(MANIFEST),
            "isolation_wp_csv": str(ISO_FINE_CSV),
        },
        "columns": {"score": SCORE, "isolation": EISO},
        "notes": [
            "Rows are 0-10% and 10-20% centrality to keep the ET-dependent view inside the important 0-20% bin.",
            "Columns split default reco_eiso by the local R=0.3 WP90 isolation threshold for that centrality row.",
            "Green WP80 curves are computed from the same local row-level score cache by centrality row and ET bin.",
            "This is a local default-isolation visualization, not a row-level R=0.4 density plot.",
            "Google Slides was not mutated.",
        ],
        "stats": {
            "truth_fraction_grid": s1,
            "density_grid": s2,
            "class_contour_grid": s3,
        },
    }
    manifest_path = OUT_DIR / "clean_busy_grid_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
