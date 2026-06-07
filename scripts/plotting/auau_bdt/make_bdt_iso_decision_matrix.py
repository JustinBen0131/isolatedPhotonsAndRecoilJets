#!/usr/bin/env python3
"""Make cut-aware decision matrices for BDT/isolation understanding."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
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
OUT_DIR = REPORT / "slideReady/default_eiso_vs_bdt_score_20260604/plot_variants_decision_matrix_20260604"

SCORE = "score_ptFine_cent7"
EISO = "reco_eiso"
CENT_BINS = [
    (0.0, 10.0, "0-10%"),
    (10.0, 20.0, "10-20%"),
    (20.0, 30.0, "20-30%"),
    (30.0, 40.0, "30-40%"),
    (40.0, 50.0, "40-50%"),
    (50.0, 60.0, "50-60%"),
    (60.0, 80.0, "60-80%"),
]
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
MUTED = "#59616f"
GRID = "#e5eaf1"
GREEN = "#168a55"
ORANGE = "#d97706"
PURPLE = "#7057a8"
GRAY = "#6b7280"
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


def cmap_to(color: str) -> LinearSegmentedColormap:
    return LinearSegmentedColormap.from_list("to_color", ["#ffffff", color], N=256)


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
        & (arrays["et"] >= ET_BINS[0][0])
        & (arrays["et"] < ET_BINS[-1][1])
        & (arrays["cent"] >= CENT_BINS[0][0])
        & (arrays["cent"] < CENT_BINS[-1][1])
    )
    return {key: val[selected] for key, val in arrays.items()}


def isolation_thresholds(cone: str) -> dict[tuple[float, float], float]:
    out: dict[tuple[float, float], float] = {}
    with ISO_FINE_CSV.open() as handle:
        for row in csv.DictReader(handle):
            if row["cone"] == cone and abs(float(row["efficiency"]) - 0.90) < 1e-9:
                out[(float(row["cent_min"]), float(row["cent_max"]))] = float(row["threshold_gev"])
    return out


def configured_bdt_wp80() -> np.ndarray:
    out = np.full((len(CENT_BINS), len(ET_BINS)), np.nan)
    with BDT_WP_CSV.open() as handle:
        for row in csv.DictReader(handle):
            if row["row_type"] != "et_bin_point" or row["wp_label"] != "WP80":
                continue
            cent = (float(row["centrality_min"]), float(row["centrality_max"]))
            et = (float(row["et_min"]), float(row["et_max"]))
            ci = next((i for i, (lo, hi, _) in enumerate(CENT_BINS) if lo == cent[0] and hi == cent[1]), None)
            ei = next((i for i, (lo, hi, _) in enumerate(ET_BINS) if lo == et[0] and hi == et[1]), None)
            if ci is not None and ei is not None:
                out[ci, ei] = float(row["threshold"])
    return out


def local_bdt_wp80(data: dict[str, np.ndarray]) -> np.ndarray:
    signal = data["is_signal"].astype(bool)
    out = np.full((len(CENT_BINS), len(ET_BINS)), np.nan)
    for ci, (clo, chi, _) in enumerate(CENT_BINS):
        for ei, (elo, ehi, _) in enumerate(ET_BINS):
            mask = signal & (data["cent"] >= clo) & (data["cent"] < chi) & (data["et"] >= elo) & (data["et"] < ehi)
            if mask.sum() >= 40:
                out[ci, ei] = float(np.nanpercentile(data["score"][mask], 20.0))
    return out


def compute_decision_matrices(data: dict[str, np.ndarray]) -> dict:
    iso = isolation_thresholds("R0.3")
    bdt = configured_bdt_wp80()
    labels = [
        ("accepted", "accepted: clean cone + BDT pass", GREEN),
        ("clean_bdt_fail", "clean cone, BDT fail", ORANGE),
        ("busy_bdt_pass", "busy cone, BDT pass", PURPLE),
        ("rejected", "rejected: busy cone + BDT fail", GRAY),
    ]
    frac = {key: np.full((len(CENT_BINS), len(ET_BINS)), np.nan) for key, _, _ in labels}
    purity = {key: np.full((len(CENT_BINS), len(ET_BINS)), np.nan) for key, _, _ in labels}
    entries = np.zeros((len(CENT_BINS), len(ET_BINS)), dtype=int)
    stats: dict[str, dict] = {}
    signal = data["is_signal"].astype(bool)
    for ci, (clo, chi, clabel) in enumerate(CENT_BINS):
        iso_cut = iso[(clo, chi)]
        for ei, (elo, ehi, elabel) in enumerate(ET_BINS):
            base = (data["cent"] >= clo) & (data["cent"] < chi) & (data["et"] >= elo) & (data["et"] < ehi)
            n = int(base.sum())
            entries[ci, ei] = n
            if n == 0 or not np.isfinite(bdt[ci, ei]):
                continue
            clean = data["eiso"] < iso_cut
            bdt_pass = data["score"] >= bdt[ci, ei]
            masks = {
                "accepted": base & clean & bdt_pass,
                "clean_bdt_fail": base & clean & ~bdt_pass,
                "busy_bdt_pass": base & ~clean & bdt_pass,
                "rejected": base & ~clean & ~bdt_pass,
            }
            cell_key = f"{clabel}_{elabel}"
            stats[cell_key] = {
                "entries": n,
                "iso_r03_wp90_gev": iso_cut,
                "bdt_wp80_configured": float(bdt[ci, ei]),
            }
            for key, mask in masks.items():
                count = int(mask.sum())
                frac[key][ci, ei] = count / n
                purity[key][ci, ei] = float(signal[mask].mean()) if count else float("nan")
                stats[cell_key][f"{key}_fraction"] = float(frac[key][ci, ei])
                stats[cell_key][f"{key}_truth_fraction"] = float(purity[key][ci, ei]) if count else None
    return {"frac": frac, "purity": purity, "entries": entries, "bdt": bdt, "labels": labels, "stats": stats}


def annotate_cells(ax, values: np.ndarray, fontsize: float = 9.2) -> None:
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            val = values[i, j]
            if not np.isfinite(val):
                continue
            color = "white" if val >= 0.55 else INK
            ax.text(j, i, f"{100 * val:.0f}%", ha="center", va="center", fontsize=fontsize, color=color, fontweight="bold")


def style_matrix_axis(ax, show_y: bool, title: str) -> None:
    ax.set_title(title, fontsize=15.4, fontweight="bold", pad=8)
    ax.set_xticks(np.arange(len(ET_BINS)))
    ax.set_xticklabels([label for _, _, label in ET_BINS], fontsize=9.8)
    ax.set_yticks(np.arange(len(CENT_BINS)))
    ax.set_yticklabels([label for _, _, label in CENT_BINS] if show_y else [], fontsize=10.6)
    ax.set_xticks(np.arange(-0.5, len(ET_BINS), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(CENT_BINS), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.25)
    ax.tick_params(which="minor", bottom=False, left=False)
    for spine in ax.spines.values():
        spine.set_color("#d5dce7")


def save_decision_fraction_map(data: dict[str, np.ndarray]) -> tuple[Path, dict]:
    matrices = compute_decision_matrices(data)
    fig, axes = plt.subplots(2, 2, figsize=(W / DPI, H / DPI), dpi=DPI)
    fig.subplots_adjust(left=0.070, right=0.972, top=0.790, bottom=0.155, wspace=0.085, hspace=0.305)
    fig.text(0.045, 0.955, r"Where candidates go after the two photon-ID cuts", fontsize=24.8, fontweight="bold", ha="left", va="top")
    fig.text(
        0.045,
        0.902,
        r"Rows are centrality; columns are cluster $E_T$. Each number is the fraction of candidates in that bin assigned to the cut outcome.",
        fontsize=13.6,
        color=MUTED,
        ha="left",
        va="top",
    )
    for ax, (key, title, color), show_y in zip(axes.ravel(), matrices["labels"], [True, False, True, False]):
        ax.imshow(matrices["frac"][key], vmin=0, vmax=1, cmap=cmap_to(color), aspect="auto")
        style_matrix_axis(ax, show_y, title)
        annotate_cells(ax, matrices["frac"][key])
    axes[0, 0].set_ylabel("centrality", fontsize=12.4)
    axes[1, 0].set_ylabel("centrality", fontsize=12.4)
    fig.text(0.520, 0.085, r"cluster $E_T$ bin [GeV]", fontsize=13.0, ha="center", va="center")
    fig.text(
        0.070,
        0.050,
        r"Cut definitions in this row cache: clean cone = reco $E_T^{iso}<I_{90}^{R=0.3}(c)$; BDT pass = configured THE-41 WP80 threshold in the same centrality and $E_T$ bin.",
        fontsize=11.0,
        color=MUTED,
    )
    out = OUT_DIR / "01_decision_outcome_fraction_matrix.png"
    fig.savefig(out)
    plt.close(fig)
    return out, matrices


def save_purity_map(matrices: dict) -> Path:
    fig, axes = plt.subplots(2, 2, figsize=(W / DPI, H / DPI), dpi=DPI)
    fig.subplots_adjust(left=0.070, right=0.972, top=0.790, bottom=0.135, wspace=0.085, hspace=0.305)
    fig.text(0.045, 0.955, "Same bins, now colored by truth-photon fraction", fontsize=24.8, fontweight="bold", ha="left", va="top")
    fig.text(
        0.045,
        0.902,
        "This separates occupancy from class content: the accepted region is truth-rich, while busy-cone outcomes expose the inclusive-tail leakage.",
        fontsize=13.6,
        color=MUTED,
        ha="left",
        va="top",
    )
    for ax, (key, title, _), show_y in zip(axes.ravel(), matrices["labels"], [True, False, True, False]):
        ax.imshow(matrices["purity"][key], vmin=0, vmax=1, cmap="RdYlBu_r", aspect="auto")
        style_matrix_axis(ax, show_y, title)
        annotate_cells(ax, matrices["purity"][key])
    axes[0, 0].set_ylabel("centrality", fontsize=12.4)
    axes[1, 0].set_ylabel("centrality", fontsize=12.4)
    fig.text(0.520, 0.075, r"cluster $E_T$ bin [GeV]", fontsize=13.0, ha="center", va="center")
    out = OUT_DIR / "02_decision_outcome_truth_fraction_matrix.png"
    fig.savefig(out)
    plt.close(fig)
    return out


def save_threshold_reference() -> tuple[Path, dict]:
    bdt = configured_bdt_wp80()
    iso_r04 = isolation_thresholds("R0.4")
    iso_r03 = isolation_thresholds("R0.3")
    iso04 = np.tile(np.asarray([iso_r04[(lo, hi)] for lo, hi, _ in CENT_BINS])[:, None], (1, len(ET_BINS)))
    iso03 = np.tile(np.asarray([iso_r03[(lo, hi)] for lo, hi, _ in CENT_BINS])[:, None], (1, len(ET_BINS)))

    fig, axes = plt.subplots(1, 3, figsize=(W / DPI, H / DPI), dpi=DPI)
    fig.subplots_adjust(left=0.055, right=0.965, top=0.800, bottom=0.130, wspace=0.115)
    fig.text(0.045, 0.955, "Thresholds behind the decision map", fontsize=24.8, fontweight="bold", ha="left", va="top")
    fig.text(
        0.045,
        0.902,
        r"These are the cut values, organized on the same centrality x $E_T$ grid. BDT WP80 varies with both axes; isolation WP90 varies with centrality.",
        fontsize=13.6,
        color=MUTED,
        ha="left",
        va="top",
    )
    panels = [
        (bdt, "configured BDT WP80", "Greens", lambda x: f"{x:.3f}"),
        (iso04, r"R=0.4 isolation WP90 [GeV]", "Oranges", lambda x: f"{x:.1f}"),
        (iso03, r"default-cache R=0.3 WP90 [GeV]", "Blues", lambda x: f"{x:.1f}"),
    ]
    for idx, (ax, (data, title, cmap, fmt)) in enumerate(zip(axes, panels)):
        im = ax.imshow(data, aspect="auto", cmap=cmap)
        style_matrix_axis(ax, idx == 0, title)
        norm = im.norm
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                val = data[i, j]
                if np.isfinite(val):
                    text_color = "white" if norm(val) > 0.68 else INK
                    ax.text(j, i, fmt(val), ha="center", va="center", fontsize=8.4, color=text_color, fontweight="bold")
        ax.set_xlabel(r"cluster $E_T$ bin [GeV]", fontsize=11.3)
    axes[0].set_ylabel("centrality", fontsize=12.0)
    out = OUT_DIR / "03_threshold_reference_cent_et_matrix.png"
    fig.savefig(out)
    plt.close(fig)
    return out, {
        "configured_bdt_wp80_min": float(np.nanmin(bdt)),
        "configured_bdt_wp80_max": float(np.nanmax(bdt)),
        "r04_iso_wp90_min": float(np.nanmin(iso04)),
        "r04_iso_wp90_max": float(np.nanmax(iso04)),
        "r03_iso_wp90_min": float(np.nanmin(iso03)),
        "r03_iso_wp90_max": float(np.nanmax(iso03)),
    }


def make_contact(paths: list[Path]) -> Path:
    thumbs = []
    for path in paths:
        img = Image.open(path).convert("RGB")
        img.thumbnail((760, 430), Image.Resampling.LANCZOS)
        thumbs.append((path, img.copy()))
    sheet = Image.new("RGB", (800 * len(thumbs), 500), "white")
    draw = ImageDraw.Draw(sheet)
    try:
        font = ImageFont.truetype("Arial.ttf", 23)
    except OSError:
        font = ImageFont.load_default()
    for idx, (path, img) in enumerate(thumbs):
        x = idx * 800 + 18
        sheet.paste(img, (x, 48))
        draw.text((x, 12), path.name, fill=(20, 24, 34), font=font)
    out = OUT_DIR / "00_contact_sheet_decision_matrix.png"
    sheet.save(out)
    return out


def main() -> None:
    setup()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    data = load_arrays()
    p1, matrices = save_decision_fraction_map(data)
    p2 = save_purity_map(matrices)
    p3, threshold_stats = save_threshold_reference()
    contact = make_contact([p1, p2, p3])
    manifest = {
        "schema": "THE11_BDT_ISO_DECISION_MATRIX_V1",
        "outputs": {
            "contact_sheet": str(contact),
            "decision_fraction_matrix": str(p1),
            "truth_fraction_matrix": str(p2),
            "threshold_reference_matrix": str(p3),
        },
        "sources": {
            "row_level_cache_manifest": str(MANIFEST),
            "configured_bdt_wp80_csv": str(BDT_WP_CSV),
            "isolation_wp_csv": str(ISO_FINE_CSV),
        },
        "columns": {"score": SCORE, "isolation": EISO},
        "notes": [
            "The decision matrices use row-level default reco_eiso, exact R=0.3 WP90 isolation thresholds, and configured THE-41 BDT WP80 thresholds on the same centrality x ET grid.",
            "The threshold reference matrix separately shows configured BDT WP80 plus R=0.4 and R=0.3 isolation WP90 thresholds.",
            "The local row-level cache does not contain reco_eiso_r40, so R=0.4 candidate outcome fractions cannot be computed from local rows in this run.",
            "Google Slides was not mutated.",
        ],
        "stats": {
            "decision_cells": matrices["stats"],
            "threshold_reference": threshold_stats,
        },
    }
    manifest_path = OUT_DIR / "decision_matrix_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
