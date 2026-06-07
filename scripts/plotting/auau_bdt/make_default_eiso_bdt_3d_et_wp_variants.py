#!/usr/bin/env python3
"""Make 3D/ET-aware BDT-isolation variants for THE-11 iteration."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LogNorm, Normalize  # noqa: E402
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401, E402
from PIL import Image, ImageDraw, ImageFont  # noqa: E402


REPO = Path(__file__).resolve().parents[3]
REPORT = REPO / "dataOutput/auauTightBDTValidation/model_validation_condor_20260511_194832"
MANIFEST = REPORT / "score_caches.local.list"
BDT_WP_CSV = REPO / (
    "dataOutput/auauTightBDTValidation/THE8_branchA_ladder_scorecache_fullstat_20260527/"
    "the41_centdep_wp_slides_20260604/the41_centdep_bdt_wp_summary.csv"
)
ISO_FINE_CSV = REPO / (
    "dataOutput/auauMLDiagnosticRuns/global_etcent_inclusive3_sixpack_20260516_135439/"
    "slideReady/isolation_wp_checks/"
    "auau_sliding_isolation_wp_r03_r04_eff70_80_90_fine7_no_numbers_coefficients.csv"
)
OUT_DIR = REPORT / "slideReady/default_eiso_vs_bdt_score_20260604/plot_variants_3d_et_wp_20260604"

SCORE = "score_ptFine_cent7"
EISO = "reco_eiso"
PT_RANGE = (15.0, 35.0)
CENT_020 = (0.0, 20.0)
CENT_020_BINS = [(0.0, 10.0, "0-10%"), (10.0, 20.0, "10-20%")]
ET_BINS = [
    (15.0, 17.0, "15-17"),
    (17.0, 19.0, "17-19"),
    (19.0, 21.0, "19-21"),
    (21.0, 23.0, "21-23"),
    (23.0, 25.0, "23-25"),
    (25.0, 27.0, "25-27"),
    (27.0, 30.0, "27-30"),
    (30.0, 35.0, "30-35"),
]

W, H, DPI = 2560, 1440, 200
INK = "#121826"
MUTED = "#5c6678"
GRID = "#dfe5ee"
RED = "#c94747"
BLUE = "#2869a6"
GREEN = "#188a55"
ORANGE = "#d97706"
PURPLE = "#7057a8"


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


def iso_thresholds(cone: str = "R0.3") -> dict[tuple[float, float], float]:
    out: dict[tuple[float, float], float] = {}
    with ISO_FINE_CSV.open() as handle:
        for row in csv.DictReader(handle):
            if row["cone"] == cone and abs(float(row["efficiency"]) - 0.90) < 1e-9:
                out[(float(row["cent_min"]), float(row["cent_max"]))] = float(row["threshold_gev"])
    return out


def configured_bdt_wp80() -> dict[tuple[float, float], list[tuple[float, float]]]:
    out: dict[tuple[float, float], list[tuple[float, float]]] = {}
    with BDT_WP_CSV.open() as handle:
        for row in csv.DictReader(handle):
            if row["row_type"] != "et_bin_point" or row["wp_label"] != "WP80":
                continue
            key = (float(row["centrality_min"]), float(row["centrality_max"]))
            et_center = 0.5 * (float(row["et_min"]) + float(row["et_max"]))
            out.setdefault(key, []).append((et_center, float(row["threshold"])))
    for values in out.values():
        values.sort()
    return out


def local_wp80_curves(data: dict[str, np.ndarray]) -> dict[tuple[float, float], list[tuple[float, float]]]:
    signal = data["is_signal"].astype(bool)
    out: dict[tuple[float, float], list[tuple[float, float]]] = {}
    for clo, chi, _ in CENT_020_BINS:
        values = []
        for elo, ehi, _ in ET_BINS:
            mask = signal & (data["cent"] >= clo) & (data["cent"] < chi) & (data["et"] >= elo) & (data["et"] < ehi)
            if mask.sum() < 50:
                continue
            values.append((0.5 * (elo + ehi), float(np.nanpercentile(data["score"][mask], 20.0))))
        out[(clo, chi)] = values
    return out


def sample(mask: np.ndarray, max_points: int, seed: int) -> np.ndarray:
    idx = np.flatnonzero(mask)
    if len(idx) <= max_points:
        return idx
    rng = np.random.default_rng(seed)
    return np.sort(rng.choice(idx, size=max_points, replace=False))


def draw_header(fig, title: str, subtitle: str) -> None:
    fig.text(0.045, 0.950, title, fontsize=25.5, fontweight="bold", ha="left", va="top")
    fig.text(0.045, 0.902, subtitle, fontsize=13.4, color=MUTED, ha="left", va="top")


def save_3d_data_view(data: dict[str, np.ndarray], local_wp: dict, iso_r03: dict) -> tuple[Path, dict]:
    mask020 = (data["cent"] >= CENT_020[0]) & (data["cent"] < CENT_020[1])
    signal = data["is_signal"].astype(bool)
    sig_idx = sample(mask020 & signal, 22_000, 20260604)
    bkg_idx = sample(mask020 & ~signal, 10_000, 20260605)
    xlim = tuple(np.nanpercentile(data["eiso"][mask020], [0.5, 99.2]))
    xlim = (float(xlim[0] - 0.8), float(xlim[1] + 0.8))
    ylim = PT_RANGE
    zlim = (0.0, 1.0)

    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)
    draw_header(
        fig,
        r"3D view: the BDT working point is an $E_T$-dependent surface",
        r"0-20% centrality, local row-level default $reco\_eiso$ cache. Green sheets = same-cache WP80 curves; orange walls = isolation WP90.",
    )
    ax = fig.add_axes([0.035, 0.080, 0.770, 0.770], projection="3d")
    ax.scatter(data["eiso"][bkg_idx], data["et"][bkg_idx], data["score"][bkg_idx], s=3.5, c=BLUE, alpha=0.16, depthshade=False, label="inclusive jets")
    ax.scatter(data["eiso"][sig_idx], data["et"][sig_idx], data["score"][sig_idx], s=3.5, c=RED, alpha=0.12, depthshade=False, label="truth photons")

    xsheet = np.linspace(xlim[0], xlim[1], 2)
    for clo, chi, label in CENT_020_BINS:
        curve = np.asarray(local_wp[(clo, chi)], dtype=float)
        et = curve[:, 0]
        thr = curve[:, 1]
        xx, yy = np.meshgrid(xsheet, et)
        zz = np.tile(thr[:, None], (1, 2))
        ax.plot_surface(xx, yy, zz, color=GREEN, alpha=0.16, linewidth=0, shade=False)
        ax.plot(np.full_like(et, xlim[1] - 0.8), et, thr, color=GREEN, lw=3.0, label=f"BDT WP80 {label}")

        iso = iso_r03[(clo, chi)]
        yy2, zz2 = np.meshgrid(np.linspace(ylim[0], ylim[1], 2), np.linspace(zlim[0], zlim[1], 2))
        xx2 = np.full_like(yy2, iso)
        ax.plot_surface(xx2, yy2, zz2, color=ORANGE, alpha=0.12, linewidth=0, shade=False)
        ax.plot([iso, iso], [ylim[0], ylim[1]], [0.02, 0.02], color=ORANGE, lw=3.0, label=f"iso WP90 {label}")

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_zlim(*zlim)
    ax.set_xlabel(r"default reco $E_T^{iso}$ [GeV]", labelpad=12, fontsize=12.2)
    ax.set_ylabel(r"cluster $E_T$ [GeV]", labelpad=12, fontsize=12.2)
    ax.set_zlabel("BDT score", labelpad=10, fontsize=12.2)
    ax.view_init(elev=22, azim=-56)
    ax.grid(True)
    ax.legend(loc="upper left", bbox_to_anchor=(0.01, 0.95), frameon=True, fontsize=10.6)
    fig.text(
        0.785,
        0.720,
        "Geometry check:\n"
        "green = BDT WP80\n"
        r"surface in $E_T$"
        "\n\n"
        "orange = isolation\n"
        "WP90 boundary\n\n"
        "red/blue points show\n"
        "where candidates sit.",
        fontsize=14.2,
        ha="left",
        va="top",
        color=INK,
        linespacing=1.28,
    )
    out = OUT_DIR / "01_3d_eiso_et_score_020_with_wp_surfaces.png"
    fig.savefig(out)
    plt.close(fig)
    return out, {"sampled_signal": int(len(sig_idx)), "sampled_background": int(len(bkg_idx))}


def save_et_score_map(data: dict[str, np.ndarray], local_wp: dict) -> tuple[Path, dict]:
    mask020 = (data["cent"] >= CENT_020[0]) & (data["cent"] < CENT_020[1])
    signal = data["is_signal"].astype(bool)
    xedges = np.linspace(PT_RANGE[0], PT_RANGE[1], 68)
    yedges = np.linspace(0.0, 1.0, 70)
    counts, _, _ = np.histogram2d(data["et"][mask020], data["score"][mask020], bins=[xedges, yedges])
    sig_counts, _, _ = np.histogram2d(data["et"][mask020 & signal], data["score"][mask020 & signal], bins=[xedges, yedges])
    prob = np.divide(sig_counts, counts, out=np.full_like(counts, np.nan), where=counts >= 8)

    fig, ax = plt.subplots(figsize=(W / DPI, H / DPI), dpi=DPI)
    fig.subplots_adjust(left=0.090, right=0.875, top=0.830, bottom=0.140)
    draw_header(
        fig,
        r"Corrected view: BDT WP80 moves with $E_T$",
        r"0-20% centrality. Color = truth-photon fraction; black contours = total density; green curves are WP80 in the two 10% centrality slices.",
    )
    mesh = ax.pcolormesh(xedges, yedges, prob.T, cmap="RdYlBu_r", norm=Normalize(0, 1), shading="auto")
    positive = counts[counts > 0]
    if positive.size:
        levels = np.unique(np.nanpercentile(positive, [55, 75, 90, 97]))
        xx = 0.5 * (xedges[:-1] + xedges[1:])
        yy = 0.5 * (yedges[:-1] + yedges[1:])
        ax.contour(xx, yy, counts.T, levels=levels, colors="black", linewidths=0.85, alpha=0.62)

    styles = [("-", "0-10% WP80"), ("--", "10-20% WP80")]
    for (clo, chi, _), (ls, label) in zip(CENT_020_BINS, styles):
        curve = np.asarray(local_wp[(clo, chi)], dtype=float)
        ax.plot(curve[:, 0], curve[:, 1], color=GREEN, lw=3.2, ls=ls, label=label)

    ax.text(
        0.03,
        0.94,
        "This is the missing piece:\n"
        r"the cut is a curve in $E_T$,"
        "\nnot one horizontal score.",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=13.3,
        color=INK,
        linespacing=1.25,
        bbox={"facecolor": "white", "edgecolor": "#d7dee8", "boxstyle": "round,pad=0.35", "alpha": 0.94},
    )
    ax.set_xlim(*PT_RANGE)
    ax.set_ylim(0.0, 1.0)
    ax.set_xlabel(r"cluster $E_T$ [GeV]", fontsize=15.8)
    ax.set_ylabel("BDT score", fontsize=15.8)
    ax.tick_params(labelsize=12.6)
    ax.grid(True, color=GRID, lw=0.7)
    ax.legend(loc="lower right", frameon=True, facecolor="white", edgecolor="#d7dee8", fontsize=12.0)
    cax = fig.add_axes([0.900, 0.205, 0.018, 0.560])
    fig.colorbar(mesh, cax=cax, label="truth-photon fraction")
    out = OUT_DIR / "02_et_score_truth_fraction_020_with_wp_curves.png"
    fig.savefig(out)
    plt.close(fig)
    return out, {"entries_020": int(mask020.sum())}


def save_configured_surface(config_wp: dict) -> tuple[Path, dict]:
    points = []
    for (clo, chi), values in config_wp.items():
        c = 0.5 * (clo + chi)
        for et, thr in values:
            points.append((c, et, thr))
    arr = np.asarray(points, dtype=float)
    cvals = np.unique(arr[:, 0])
    evals = np.unique(arr[:, 1])
    z = np.full((len(evals), len(cvals)), np.nan)
    for c, et, thr in arr:
        z[np.where(evals == et)[0][0], np.where(cvals == c)[0][0]] = thr
    cc, ee = np.meshgrid(cvals, evals)

    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)
    draw_header(
        fig,
        r"Configured BDT WP80 is a surface in centrality and $E_T$",
        "This is the working-point object itself, before adding isolation or candidate density.",
    )
    ax = fig.add_axes([0.080, 0.100, 0.760, 0.740], projection="3d")
    surf = ax.plot_surface(cc, ee, z, cmap="viridis", edgecolor="white", linewidth=0.65, alpha=0.92)
    ax.scatter(arr[:, 0], arr[:, 1], arr[:, 2], c=arr[:, 2], cmap="viridis", s=45, edgecolor=INK, linewidth=0.35)
    ax.set_xlabel("centrality [%]", labelpad=12, fontsize=12.5)
    ax.set_ylabel(r"cluster $E_T$ [GeV]", labelpad=12, fontsize=12.5)
    ax.set_zlabel("BDT WP80 threshold", labelpad=10, fontsize=12.5)
    ax.view_init(elev=26, azim=-132)
    cax = fig.add_axes([0.870, 0.255, 0.020, 0.430])
    fig.colorbar(surf, cax=cax, label="BDT threshold")
    out = OUT_DIR / "03_configured_bdt_wp80_cent_et_surface.png"
    fig.savefig(out)
    plt.close(fig)
    return out, {"points": int(len(points))}


def make_contact(paths: list[Path]) -> Path:
    thumbs = []
    for path in paths:
        img = Image.open(path).convert("RGB")
        img.thumbnail((780, 440), Image.Resampling.LANCZOS)
        thumbs.append((path, img.copy()))
    sheet = Image.new("RGB", (820 * len(thumbs), 500), "white")
    draw = ImageDraw.Draw(sheet)
    try:
        font = ImageFont.truetype("Arial.ttf", 24)
    except OSError:
        font = ImageFont.load_default()
    for idx, (path, img) in enumerate(thumbs):
        x = idx * 820 + 20
        sheet.paste(img, (x, 42))
        draw.text((x, 12), path.name, fill=(20, 24, 34), font=font)
    out = OUT_DIR / "00_contact_sheet_3d_et_wp_variants.png"
    sheet.save(out)
    return out


def main() -> None:
    setup()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    data = load_arrays()
    data020 = {key: val[(data["cent"] >= CENT_020[0]) & (data["cent"] < CENT_020[1])] for key, val in data.items()}
    local_wp = local_wp80_curves(data020)
    # local_wp is computed after slicing to 0-20, so centrality tests still work on 0-20 data values.
    iso_r03 = iso_thresholds("R0.3")
    config_wp = configured_bdt_wp80()
    p1, s1 = save_3d_data_view(data, local_wp, iso_r03)
    p2, s2 = save_et_score_map(data, local_wp)
    p3, s3 = save_configured_surface(config_wp)
    contact = make_contact([p1, p2, p3])
    manifest = {
        "schema": "THE11_3D_ET_AWARE_BDT_ISO_VARIANTS_V1",
        "outputs": {
            "contact_sheet": str(contact),
            "three_d_data_view": str(p1),
            "et_score_wp_curves": str(p2),
            "configured_bdt_wp80_surface": str(p3),
        },
        "sources": {
            "row_level_cache_manifest": str(MANIFEST),
            "configured_bdt_wp_csv": str(BDT_WP_CSV),
            "isolation_wp_csv": str(ISO_FINE_CSV),
        },
        "columns": {"score": SCORE, "isolation": EISO},
        "notes": [
            "The data-cloud plots use local row-level default reco_eiso and score_ptFine_cent7 caches.",
            "The WP80 curves on the data-cloud plots are computed from the same local row-level cache by ET and centrality, so the overlay is internally consistent.",
            "The configured surface plot uses the THE-41 WP80 summary CSV and is shown separately to avoid mixing score definitions.",
            "Google Slides was not mutated.",
        ],
        "stats": {
            "three_d_data_view": s1,
            "et_score_wp_curves": s2,
            "configured_bdt_wp80_surface": s3,
        },
    }
    manifest_path = OUT_DIR / "three_d_et_wp_variants_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
