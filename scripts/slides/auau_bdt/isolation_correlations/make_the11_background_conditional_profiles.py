#!/usr/bin/env python3
"""Build THE-11 background-only BDT/isolation conditional profile diagnostics."""

from __future__ import annotations

import argparse
import csv
import json
import math
import textwrap
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np


REPO = Path(__file__).resolve().parents[4]
DEFAULT_OUT_DIR = (
    REPO
    / "dataOutput/auauTightBDTValidation/the11_weighted_bdt_iso_fullstat_basev3e_20260606/"
    "backgroundConditionalProfiles"
)
DEFAULT_REMOTE_REPORT = Path(
    "/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_eiso_cone_raw_20260518_2220/"
    "reports/model_validation_condor_the11_weighted_bdt_iso_fullstat_basev3e_20260606_1645"
)

SCORE_COLUMN = "score_centInput_pt1535"
ISO_COLUMN_R30 = "reco_eiso_r30"
ISO_COLUMN_R40 = "reco_eiso_r40"
CENT_BINS = [
    (0.0, 20.0, "0-20%", "cent_0_20"),
    (20.0, 50.0, "20-50%", "cent_20_50"),
    (50.0, 80.0, "50-80%", "cent_50_80"),
]
ET_BINS = [
    (15.0, 35.0, "15-35 GeV", "et_15_35"),
    (15.0, 20.0, "15-20 GeV", "et_15_20"),
    (20.0, 25.0, "20-25 GeV", "et_20_25"),
    (25.0, 30.0, "25-30 GeV", "et_25_30"),
    (30.0, 35.0, "30-35 GeV", "et_30_35"),
]
BDT_BANDS = [
    (0.00, 0.20, "0.0-0.2", "#5d6673"),
    (0.20, 0.40, "0.2-0.4", "#2f6f9f"),
    (0.40, 0.60, "0.4-0.6", "#16827a"),
    (0.60, 0.80, "0.6-0.8", "#4f8f2f"),
    (0.80, 0.90, "0.8-0.9", "#b18b00"),
    (0.90, 0.95, "0.9-0.95", "#d06a00"),
    (0.95, 1.01, "0.95-1.0", "#b12632"),
]
EISO_BANDS = [
    (-np.inf, 0.0, r"$E_T^{iso}<0$", "#1764a8"),
    (0.0, 3.0, r"$0<E_T^{iso}<3$", "#16827a"),
    (3.0, 6.0, r"$3<E_T^{iso}<6$", "#8a8d1e"),
    (6.0, 10.0, r"$6<E_T^{iso}<10$", "#c16b00"),
    (10.0, np.inf, r"$E_T^{iso}>10$", "#a32035"),
]
EISO_HIST_BINS = np.linspace(-10.0, 18.0, 57)
SCORE_HIST_BINS = np.linspace(0.0, 1.0, 51)
FRACTION_EISO_BINS = np.asarray([-10.0, -5.0, 0.0, 3.0, 6.0, 10.0, 14.0, 18.0])
HIGH_BDT_THRESHOLDS = [0.8, 0.9]
INTEGRATED_ET_KEY = "et_15_35"

W, H, DPI = 2560, 1440, 200
WHITE = "#ffffff"
INK = "#151515"
MUTED = "#596271"
GRID = "#dfe5ec"


def finite_float(value: float) -> float | None:
    if value is None:
        return None
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def weighted_quantile(values: np.ndarray, weights: np.ndarray, quantile: float) -> float:
    mask = np.isfinite(values) & np.isfinite(weights) & (weights > 0.0)
    if not np.any(mask):
        return float("nan")
    v = values[mask].astype("float64", copy=False)
    w = weights[mask].astype("float64", copy=False)
    order = np.argsort(v)
    v = v[order]
    w = w[order]
    cdf = np.cumsum(w)
    target = quantile * cdf[-1]
    return float(np.interp(target, cdf, v))


def weighted_pearson(x: np.ndarray, y: np.ndarray, weights: np.ndarray) -> float:
    mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(weights) & (weights > 0.0)
    if np.count_nonzero(mask) < 3:
        return float("nan")
    x = x[mask].astype("float64", copy=False)
    y = y[mask].astype("float64", copy=False)
    w = weights[mask].astype("float64", copy=False)
    wsum = np.sum(w)
    mx = np.sum(w * x) / wsum
    my = np.sum(w * y) / wsum
    vx = np.sum(w * (x - mx) ** 2) / wsum
    vy = np.sum(w * (y - my) ** 2) / wsum
    if vx <= 0.0 or vy <= 0.0:
        return float("nan")
    cov = np.sum(w * (x - mx) * (y - my)) / wsum
    return float(cov / math.sqrt(vx * vy))


def rank_ordinal(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype="float64")
    ranks[order] = np.arange(len(values), dtype="float64")
    return ranks


def spearman_rho(x: np.ndarray, y: np.ndarray, max_points: int = 250_000) -> float:
    mask = np.isfinite(x) & np.isfinite(y)
    if np.count_nonzero(mask) < 3:
        return float("nan")
    x = x[mask]
    y = y[mask]
    if len(x) > max_points:
        idx = np.linspace(0, len(x) - 1, max_points).astype("int64")
        x = x[idx]
        y = y[idx]
    return float(np.corrcoef(rank_ordinal(x), rank_ordinal(y))[0, 1])


def normalized_density(values: np.ndarray, weights: np.ndarray, bins: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    counts, _ = np.histogram(values, bins=bins, weights=weights)
    widths = np.diff(bins)
    area = float(np.sum(counts))
    if area <= 0.0:
        return np.zeros_like(counts, dtype="float64"), counts.astype("float64")
    return counts.astype("float64") / (area * widths), counts.astype("float64")


def hist_probability(values: np.ndarray, weights: np.ndarray, bins: np.ndarray) -> np.ndarray:
    counts, _ = np.histogram(values, bins=bins, weights=weights)
    total = float(np.sum(counts))
    if total <= 0.0:
        return np.zeros_like(counts, dtype="float64")
    return counts.astype("float64") / total


def js_distance(p: np.ndarray, q: np.ndarray) -> float:
    p = np.asarray(p, dtype="float64")
    q = np.asarray(q, dtype="float64")
    if np.sum(p) <= 0.0 or np.sum(q) <= 0.0:
        return float("nan")
    p = p / np.sum(p)
    q = q / np.sum(q)
    m = 0.5 * (p + q)

    def kl(a: np.ndarray, b: np.ndarray) -> float:
        mask = a > 0.0
        return float(np.sum(a[mask] * np.log(a[mask] / b[mask])))

    return float(math.sqrt(max(0.0, 0.5 * kl(p, m) + 0.5 * kl(q, m))))


def read_cache_paths(report_dir: Path, cache_manifest: Path | None = None) -> list[Path]:
    manifest = cache_manifest or (report_dir / "score_caches.list")
    if manifest.exists():
        paths = [Path(line.strip()) for line in manifest.read_text().splitlines() if line.strip()]
    else:
        paths = sorted((report_dir / "score_caches").glob("score_cache_*.npz"))
    if not paths:
        raise FileNotFoundError(f"No score caches found from {manifest}")
    return paths


def empty_cells() -> dict[tuple[str, str], dict[str, list[np.ndarray]]]:
    cells: dict[tuple[str, str], dict[str, list[np.ndarray]]] = {}
    for _, _, _, et_key in ET_BINS:
        for _, _, _, cent_key in CENT_BINS:
            cells[(et_key, cent_key)] = {
                "bkg_eiso": [],
                "bkg_score": [],
                "bkg_weight": [],
                "sig_score": [],
                "sig_weight": [],
            }
    return cells


def summarize_cell(eiso: np.ndarray, score: np.ndarray, weight: np.ndarray, sig_score: np.ndarray, sig_weight: np.ndarray) -> dict[str, Any]:
    weighted_n = float(np.sum(weight))
    wp80 = weighted_quantile(sig_score, sig_weight, 0.20)
    thresholds = HIGH_BDT_THRESHOLDS + [wp80]
    threshold_labels = ["BDT > 0.8", "BDT > 0.9", "BDT > WP80"]

    bdt_profiles: dict[str, Any] = {}
    bdt_probabilities: dict[str, np.ndarray] = {}
    for lo, hi, label, color in BDT_BANDS:
        mask = (score >= lo) & (score < hi)
        density, counts = normalized_density(eiso[mask], weight[mask], EISO_HIST_BINS)
        prob = hist_probability(eiso[mask], weight[mask], EISO_HIST_BINS)
        bdt_probabilities[label] = prob
        bdt_profiles[label] = {
            "range": [lo, min(hi, 1.0)],
            "color": color,
            "raw_n": int(np.count_nonzero(mask)),
            "weighted_n": float(np.sum(weight[mask])),
            "density": density.tolist(),
            "counts": counts.tolist(),
            "median_eiso": finite_float(weighted_quantile(eiso[mask], weight[mask], 0.50)),
            "q25_eiso": finite_float(weighted_quantile(eiso[mask], weight[mask], 0.25)),
            "q75_eiso": finite_float(weighted_quantile(eiso[mask], weight[mask], 0.75)),
        }

    eiso_profiles: dict[str, Any] = {}
    eiso_probabilities: dict[str, np.ndarray] = {}
    for lo, hi, label, color in EISO_BANDS:
        mask = (eiso >= lo) & (eiso < hi)
        density, counts = normalized_density(score[mask], weight[mask], SCORE_HIST_BINS)
        prob = hist_probability(score[mask], weight[mask], SCORE_HIST_BINS)
        eiso_probabilities[label] = prob
        eiso_profiles[label] = {
            "range": [finite_float(lo), finite_float(hi)],
            "color": color,
            "raw_n": int(np.count_nonzero(mask)),
            "weighted_n": float(np.sum(weight[mask])),
            "density": density.tolist(),
            "counts": counts.tolist(),
            "median_score": finite_float(weighted_quantile(score[mask], weight[mask], 0.50)),
            "q25_score": finite_float(weighted_quantile(score[mask], weight[mask], 0.25)),
            "q75_score": finite_float(weighted_quantile(score[mask], weight[mask], 0.75)),
        }

    fraction_bins: list[dict[str, Any]] = []
    centers = 0.5 * (FRACTION_EISO_BINS[:-1] + FRACTION_EISO_BINS[1:])
    slopes: dict[str, float | None] = {}
    for bin_idx, (lo, hi) in enumerate(zip(FRACTION_EISO_BINS[:-1], FRACTION_EISO_BINS[1:])):
        band = (eiso >= lo) & (eiso < hi)
        denom_w = float(np.sum(weight[band]))
        denom_w2 = float(np.sum(weight[band] ** 2))
        row: dict[str, Any] = {
            "eiso_lo": float(lo),
            "eiso_hi": float(hi),
            "center": float(centers[bin_idx]),
            "raw_n": int(np.count_nonzero(band)),
            "weighted_n": denom_w,
            "thresholds": {},
        }
        for threshold, threshold_label in zip(thresholds, threshold_labels):
            numerator = band & (score >= threshold)
            num_w = float(np.sum(weight[numerator]))
            frac = num_w / denom_w if denom_w > 0.0 else float("nan")
            n_eff = (denom_w * denom_w / denom_w2) if denom_w2 > 0.0 else 0.0
            err = math.sqrt(max(0.0, frac * (1.0 - frac)) / n_eff) if n_eff > 0.0 and math.isfinite(frac) else float("nan")
            row["thresholds"][threshold_label] = {
                "threshold": finite_float(threshold),
                "fraction": finite_float(frac),
                "error": finite_float(err),
                "numerator_weight": num_w,
            }
        fraction_bins.append(row)

    for threshold_label in threshold_labels:
        x_vals = []
        y_vals = []
        w_vals = []
        for row in fraction_bins:
            frac = row["thresholds"][threshold_label]["fraction"]
            if frac is not None and row["weighted_n"] > 0.0:
                x_vals.append(row["center"])
                y_vals.append(frac)
                w_vals.append(row["weighted_n"])
        if len(x_vals) >= 3:
            try:
                slopes[threshold_label] = float(np.polyfit(np.asarray(x_vals), np.asarray(y_vals), 1, w=np.sqrt(w_vals))[0])
            except np.linalg.LinAlgError:
                slopes[threshold_label] = None
        else:
            slopes[threshold_label] = None

    high_reference = "0.95-1.0"
    for candidate in ("0.95-1.0", "0.9-0.95", "0.8-0.9", "0.6-0.8"):
        if np.sum(bdt_probabilities[candidate]) > 0.0:
            high_reference = candidate
            break
    eiso_js_low_high = js_distance(bdt_probabilities["0.2-0.4"], bdt_probabilities["0.95-1.0"])
    eiso_js_low_high_populated = js_distance(bdt_probabilities["0.2-0.4"], bdt_probabilities[high_reference])
    bdt_js_clean_busy = js_distance(eiso_probabilities[r"$E_T^{iso}<0$"], eiso_probabilities[r"$E_T^{iso}>10$"])

    return {
        "raw_n": int(len(eiso)),
        "weighted_n": weighted_n,
        "wp80_signal_threshold_same_cache": finite_float(wp80),
        "median_eiso": finite_float(weighted_quantile(eiso, weight, 0.50)),
        "q25_eiso": finite_float(weighted_quantile(eiso, weight, 0.25)),
        "q75_eiso": finite_float(weighted_quantile(eiso, weight, 0.75)),
        "median_score": finite_float(weighted_quantile(score, weight, 0.50)),
        "spearman_score_eiso": finite_float(spearman_rho(score, eiso)),
        "weighted_pearson_score_eiso": finite_float(weighted_pearson(score, eiso, weight)),
        "eiso_js_0p2_0p4_vs_0p95_1p0": finite_float(eiso_js_low_high),
        "eiso_js_0p2_0p4_vs_highest_populated_bdt": finite_float(eiso_js_low_high_populated),
        "eiso_js_high_reference_bdt_band": high_reference,
        "bdt_js_eiso_lt0_vs_gt10": finite_float(bdt_js_clean_busy),
        "fraction_slope_per_gev": slopes,
        "bdt_band_profiles": bdt_profiles,
        "eiso_band_profiles": eiso_profiles,
        "high_bdt_fraction_vs_eiso": fraction_bins,
    }


def reduce_score_caches(report_dir: Path, cache_manifest: Path | None, iso_column: str) -> dict[str, Any]:
    cache_paths = read_cache_paths(report_dir, cache_manifest)
    cells = empty_cells()
    required = ["is_signal", "cluster_Et", "centrality", iso_column, SCORE_COLUMN, "event_weight"]
    cache_rows: list[dict[str, Any]] = []

    for cache_path in cache_paths:
        with np.load(cache_path, allow_pickle=True) as data:
            missing = [name for name in required if name not in data.files]
            if missing:
                raise KeyError(f"{cache_path} missing required columns: {missing}")
            is_signal = data["is_signal"].astype("int8", copy=False)
            et = data["cluster_Et"].astype("float32", copy=False)
            cent = data["centrality"].astype("float32", copy=False)
            eiso = data[iso_column].astype("float32", copy=False)
            score = data[SCORE_COLUMN].astype("float32", copy=False)
            weight = data["event_weight"].astype("float64", copy=False)

            finite = (
                np.isfinite(et)
                & np.isfinite(cent)
                & np.isfinite(eiso)
                & np.isfinite(score)
                & np.isfinite(weight)
                & (weight > 0.0)
                & (et >= 15.0)
                & (et < 35.0)
                & (cent >= 0.0)
                & (cent < 80.0)
            )
            cache_rows.append(
                {
                    "cache": str(cache_path),
                    "rows_total": int(len(et)),
                    "rows_selected_15_35_cent0_80": int(np.count_nonzero(finite)),
                    "background_selected": int(np.count_nonzero(finite & (is_signal == 0))),
                    "truth_selected": int(np.count_nonzero(finite & (is_signal == 1))),
                }
            )

            for et_lo, et_hi, _, et_key in ET_BINS:
                et_mask = finite & (et >= et_lo) & (et < et_hi)
                for cent_lo, cent_hi, _, cent_key in CENT_BINS:
                    mask = et_mask & (cent >= cent_lo) & (cent < cent_hi)
                    if not np.any(mask):
                        continue
                    bkg = mask & (is_signal == 0)
                    sig = mask & (is_signal == 1)
                    cell = cells[(et_key, cent_key)]
                    if np.any(bkg):
                        cell["bkg_eiso"].append(eiso[bkg])
                        cell["bkg_score"].append(score[bkg])
                        cell["bkg_weight"].append(weight[bkg])
                    if np.any(sig):
                        cell["sig_score"].append(score[sig])
                        cell["sig_weight"].append(weight[sig])

    panels: dict[str, Any] = {}
    for _, _, et_label, et_key in ET_BINS:
        panels[et_key] = {"label": et_label, "centrality": {}}
        for _, _, cent_label, cent_key in CENT_BINS:
            cell = cells[(et_key, cent_key)]
            eiso = np.concatenate(cell["bkg_eiso"]) if cell["bkg_eiso"] else np.asarray([], dtype="float32")
            score = np.concatenate(cell["bkg_score"]) if cell["bkg_score"] else np.asarray([], dtype="float32")
            weight = np.concatenate(cell["bkg_weight"]) if cell["bkg_weight"] else np.asarray([], dtype="float64")
            sig_score = np.concatenate(cell["sig_score"]) if cell["sig_score"] else np.asarray([], dtype="float32")
            sig_weight = np.concatenate(cell["sig_weight"]) if cell["sig_weight"] else np.asarray([], dtype="float64")
            panels[et_key]["centrality"][cent_key] = {
                "label": cent_label,
                **summarize_cell(eiso, score, weight, sig_score, sig_weight),
            }

    return {
        "metadata": {
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "report_dir": str(report_dir),
            "cache_manifest": str(cache_manifest or report_dir / "score_caches.list"),
            "cache_count": len(cache_paths),
            "score_column": SCORE_COLUMN,
            "iso_column": iso_column,
            "selection": "inclusive background only for profiles; 15 <= cluster_Et < 35 GeV; 0 <= centrality < 80; event_weight > 0",
            "normalization": "Each 1D curve is normalized to unit area inside its own centrality and ET slice.",
            "wp80_definition": "BDT > WP80 uses the same-cache weighted 20th percentile of truth-photon scores in that ET and centrality cell.",
        },
        "centrality_bins": [{"lo": lo, "hi": hi, "label": label, "key": key} for lo, hi, label, key in CENT_BINS],
        "et_bins": [{"lo": lo, "hi": hi, "label": label, "key": key} for lo, hi, label, key in ET_BINS],
        "bdt_bands": [{"lo": lo, "hi": min(hi, 1.0), "label": label, "color": color} for lo, hi, label, color in BDT_BANDS],
        "eiso_bands": [{"lo": finite_float(lo), "hi": finite_float(hi), "label": label, "color": color} for lo, hi, label, color in EISO_BANDS],
        "eiso_hist_bins": EISO_HIST_BINS.tolist(),
        "score_hist_bins": SCORE_HIST_BINS.tolist(),
        "fraction_eiso_bins": FRACTION_EISO_BINS.tolist(),
        "cache_rows": cache_rows,
        "panels": panels,
    }


def load_payload(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def setup_style() -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager

    available = {f.name for f in font_manager.fontManager.ttflist}
    for family in ("Times New Roman", "Times", "DejaVu Serif"):
        if family in available:
            plt.rcParams["font.family"] = family
            break
    plt.rcParams.update(
        {
            "figure.facecolor": WHITE,
            "savefig.facecolor": WHITE,
            "axes.facecolor": WHITE,
            "axes.edgecolor": INK,
            "axes.linewidth": 1.0,
            "axes.labelcolor": INK,
            "xtick.color": INK,
            "ytick.color": INK,
            "text.color": INK,
            "mathtext.default": "regular",
        }
    )


def safe_positive(values: np.ndarray, floor: float = 1.0e-5) -> np.ndarray:
    out = np.asarray(values, dtype="float64")
    out = np.where(out > 0.0, out, np.nan)
    return np.where(out < floor, np.nan, out)


def cone_radius_text(payload: dict[str, Any]) -> str:
    iso_column = payload.get("metadata", {}).get("iso_column", ISO_COLUMN_R30)
    return "0.4" if iso_column == ISO_COLUMN_R40 else "0.3"


def add_common_panel_text(ax: Any, cell: dict[str, Any]) -> None:
    wp80 = cell.get("wp80_signal_threshold_same_cache")
    rho = cell.get("spearman_score_eiso")
    text = f"N={cell['raw_n']:,}\nWP80={wp80:.3f}\nSpearman={rho:.2f}" if wp80 is not None and rho is not None else f"N={cell['raw_n']:,}"
    ax.text(
        0.035,
        0.965,
        text,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=12,
        color=INK,
        bbox=dict(boxstyle="round,pad=0.25", facecolor="#ffffff", edgecolor="#cfd7e3", linewidth=0.8, alpha=0.92),
    )


def plot_eiso_given_bdt(payload: dict[str, Any], out_path: Path, et_key: str = INTEGRATED_ET_KEY, *, grid: bool = False) -> None:
    import matplotlib.pyplot as plt

    setup_style()
    cone = cone_radius_text(payload)
    eiso_bins = np.asarray(payload["eiso_hist_bins"], dtype="float64")
    x = 0.5 * (eiso_bins[:-1] + eiso_bins[1:])
    if grid:
        et_keys = [row["key"] for row in payload["et_bins"] if row["key"] != INTEGRATED_ET_KEY]
        fig, axes = plt.subplots(len(et_keys), len(CENT_BINS), figsize=(16, 12), dpi=200, sharex=True, sharey=True)
        title = r"Inclusive background: isolation shape in BDT-score bands, split by $E_T^{cluster}$"
    else:
        et_keys = [et_key]
        fig, axes = plt.subplots(1, len(CENT_BINS), figsize=(W / DPI, H / DPI), dpi=DPI, sharex=True, sharey=True)
        axes = np.asarray([axes])
        title = r"Inclusive background: where high-BDT candidates live in isolation"

    for row_idx, current_et_key in enumerate(et_keys):
        et_label = payload["panels"][current_et_key]["label"]
        for col_idx, (_, _, cent_label, cent_key) in enumerate(CENT_BINS):
            ax = axes[row_idx, col_idx]
            cell = payload["panels"][current_et_key]["centrality"][cent_key]
            for band in payload["bdt_bands"]:
                prof = cell["bdt_band_profiles"][band["label"]]
                y = safe_positive(np.asarray(prof["density"], dtype="float64"))
                total_band_n = sum(
                    payload["panels"][key]["centrality"][cent["key"]]["bdt_band_profiles"][band["label"]]["raw_n"]
                    for key in et_keys
                    for cent in payload["centrality_bins"]
                )
                legend_text = f"{band['label']} (empty)" if total_band_n == 0 else band["label"]
                label = legend_text if row_idx == 0 and col_idx == 0 else None
                ax.step(x, y, where="mid", color=band["color"], linewidth=1.8, label=label)
            ax.set_yscale("log")
            ax.set_xlim(eiso_bins[0], eiso_bins[-1])
            ax.set_ylim(7.0e-4, 0.55)
            ax.grid(True, color=GRID, linewidth=0.7, alpha=0.85)
            ax.axvline(0.0, color="#8a8f99", linewidth=1.0, linestyle="--", alpha=0.75)
            if row_idx == 0:
                ax.set_title(cent_label, fontsize=18 if not grid else 13, fontweight="bold", pad=8)
            if col_idx == 0:
                ax.set_ylabel((et_label + "\n") + r"unit-area density", fontsize=14 if grid else 17)
            if row_idx == len(et_keys) - 1:
                ax.set_xlabel(rf"reco $E_T^{{iso}}$, $\Delta R<{cone}$ [GeV]", fontsize=14 if grid else 17)
            ax.tick_params(labelsize=11 if grid else 14)
            add_common_panel_text(ax, cell)

    fig.suptitle(title, fontsize=22 if not grid else 19, fontweight="bold", y=0.985 if not grid else 0.972)
    subtitle = (
        r"Unit-area $p(E_T^{iso}\mid$ BDT band) for inclusive background in each panel."
        if not grid
        else r"Control view: each curve is normalized within its own centrality and $E_T^{cluster}$ slice."
    )
    fig.text(0.5, 0.925 if not grid else 0.935, subtitle, ha="center", va="center", fontsize=14 if not grid else 12, color=MUTED)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=7, frameon=False, fontsize=13 if not grid else 10, bbox_to_anchor=(0.5, 0.018))
    fig.tight_layout(rect=[0.045, 0.075, 0.99, 0.885 if not grid else 0.905])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=DPI)
    plt.close(fig)


def plot_bdt_given_eiso(payload: dict[str, Any], out_path: Path, et_key: str = INTEGRATED_ET_KEY, *, grid: bool = False) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import FancyBboxPatch

    setup_style()
    score_bins = np.asarray(payload["score_hist_bins"], dtype="float64")
    x = 0.5 * (score_bins[:-1] + score_bins[1:])
    if not grid:
        cone = cone_radius_text(payload)
        fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)
        fig.patch.set_facecolor(WHITE)
        title = "Background BDT shape depends on isolation sideband"
        fig.text(0.055, 0.962, title, ha="left", va="top", fontsize=30, fontweight="bold", color=INK)

        for y_rule, lw in ((0.885, 1.25), (0.746, 1.25)):
            fig.add_artist(
                Line2D(
                    [0.135, 0.865],
                    [y_rule, y_rule],
                    transform=fig.transFigure,
                    color="#cbd3dc",
                    linewidth=lw,
                    solid_capstyle="butt",
                    zorder=0.18,
                )
            )

        fig.text(
            0.5,
            0.860,
            r"Raw reconstructed isolation regions overlaid, $\Delta R<"
            + cone
            + r"$",
            ha="center",
            va="center",
            fontsize=17.0,
            color=INK,
        )

        def display_eiso_label(band_label: str) -> str:
            label_map = {
                r"$E_T^{iso}<0$": r"$E_T^{iso}<0$ GeV",
                r"$0<E_T^{iso}<3$": r"$0<E_T^{iso}<3$ GeV",
                r"$3<E_T^{iso}<6$": r"$3<E_T^{iso}<6$ GeV",
                r"$6<E_T^{iso}<10$": r"$6<E_T^{iso}<10$ GeV",
                r"$E_T^{iso}>10$": r"$E_T^{iso}>10$ GeV",
            }
            return label_map.get(band_label, band_label)

        axes = []
        for idx, (_, _, cent_label, cent_key) in enumerate(CENT_BINS):
            ax = fig.add_axes([0.070 + idx * 0.305, 0.333, 0.260, 0.342])
            axes.append(ax)
            cell = payload["panels"][et_key]["centrality"][cent_key]
            for band in payload["eiso_bands"]:
                prof = cell["eiso_band_profiles"][band["label"]]
                y = safe_positive(np.asarray(prof["density"], dtype="float64"))
                ax.step(x, y, where="mid", color=band["color"], linewidth=2.15)
            ax.set_yscale("log")
            ax.set_xlim(0.0, 1.0)
            ax.set_ylim(1.0e-2, 35.0)
            ax.grid(True, color=GRID, linewidth=0.8, alpha=0.9)
            ax.set_title(cent_label + " centrality", fontsize=17, fontweight="bold", pad=8)
            ax.tick_params(labelsize=12.5)
            if idx == 0:
                ax.set_ylabel("unit-area density", fontsize=15, labelpad=8)
            else:
                ax.tick_params(labelleft=False)
            stat_text = (
                f"N={cell['raw_n']:,}\n"
                + f"median BDT={cell['median_score']:.3f}\n"
                + f"Spearman={cell['spearman_score_eiso']:.2f}"
            )
            ax.text(
                0.035,
                0.055,
                stat_text,
                transform=ax.transAxes,
                ha="left",
                va="bottom",
                fontsize=12.2,
                color=INK,
                bbox=dict(
                    boxstyle="round,pad=0.25",
                    facecolor="#ffffff",
                    edgecolor="#cfd7e3",
                    linewidth=0.8,
                    alpha=0.94,
                ),
            )

        legend_items = [(display_eiso_label(band["label"]), band["color"]) for band in payload["eiso_bands"]]
        legend_positions = [(0.232, 0.818), (0.500, 0.818), (0.768, 0.818), (0.378, 0.779), (0.625, 0.779)]
        for (legend_label, legend_color), (legend_x, legend_y) in zip(legend_items, legend_positions):
            fig.add_artist(
                Line2D(
                    [legend_x - 0.095, legend_x - 0.040],
                    [legend_y, legend_y],
                    transform=fig.transFigure,
                    color=legend_color,
                    linewidth=2.8,
                    solid_capstyle="butt",
                    zorder=3.0,
                )
            )
            fig.text(
                legend_x - 0.025,
                legend_y,
                legend_label,
                ha="left",
                va="center",
                fontsize=16.0,
                color=INK,
            )

        fig.text(0.5, 0.275, "BDT score", ha="center", va="center", fontsize=16, color=INK)

        left_card_x, right_card_x = 0.052, 0.515
        card_y, card_w, card_h = 0.055, 0.420, 0.176
        left_text_x, right_text_x = 0.072, 0.535
        heading_y, header_rule_y = 0.211, 0.178
        left_card = FancyBboxPatch(
            (left_card_x, card_y),
            card_w,
            card_h,
            boxstyle="round,pad=0.010,rounding_size=0.008",
            transform=fig.transFigure,
            facecolor="#fffaf1",
            edgecolor="#bf7a11",
            linewidth=1.05,
            zorder=0.1,
        )
        right_card = FancyBboxPatch(
            (right_card_x, card_y),
            card_w,
            card_h,
            boxstyle="round,pad=0.010,rounding_size=0.008",
            transform=fig.transFigure,
            facecolor="#f8fafc",
            edgecolor="#66707e",
            linewidth=1.05,
            zorder=0.1,
        )
        fig.add_artist(left_card)
        fig.add_artist(right_card)
        fig.add_artist(
            Line2D(
                [left_card_x + 0.018, left_card_x + card_w - 0.018],
                [header_rule_y, header_rule_y],
                transform=fig.transFigure,
                color="#e1c491",
                linewidth=0.95,
                zorder=0.25,
            )
        )
        fig.add_artist(
            Line2D(
                [right_card_x + 0.018, right_card_x + card_w - 0.018],
                [header_rule_y, header_rule_y],
                transform=fig.transFigure,
                color="#d2d8e0",
                linewidth=0.95,
                zorder=0.25,
            )
        )

        def label_body(x: float, y: float, label: str, body: str, offset: float, size: float) -> None:
            fig.text(x, y, label, ha="left", va="top", fontsize=size, color=INK, fontweight="bold")
            fig.text(x + offset, y, body, ha="left", va="top", fontsize=size, color=INK)

        fig.text(left_text_x, heading_y, "What is plotted", ha="left", va="top", fontsize=16.7, fontweight="bold", color="#a96d12")
        label_body(left_text_x, 0.169, "Sample:", "inclusive background, cluster ET = 15-35 GeV.", 0.063, 13.2)
        label_body(left_text_x, 0.137, "Curves:", r"unit-area $p(\mathrm{BDT})$ after slicing by isolation.", 0.062, 13.2)
        label_body(left_text_x, 0.101, "Reading:", "shapes hint at subtle BDT-isolation correlation.", 0.057, 12.55)
        fig.text(right_text_x, heading_y, "What it means", ha="left", va="top", fontsize=16.7, fontweight="bold", color="#5e6876")
        label_body(right_text_x, 0.169, "Observed:", "positive-isolation slices lose the high-BDT shoulder.", 0.078, 13.0)
        label_body(right_text_x, 0.137, "Physics:", "busy cone activity makes background less photon-like.", 0.066, 13.0)
        label_body(right_text_x, 0.101, "ABCD concern:", "isolation and BDT score should be independent.", 0.103, 11.65)
        fig.text(
            right_text_x + 0.014,
            0.072,
            "• But...we see a change in shape with varied iso slices.",
            ha="left",
            va="top",
            fontsize=11.65,
            color=INK,
        )

        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=DPI)
        plt.close(fig)
        return

    if grid:
        et_keys = [row["key"] for row in payload["et_bins"] if row["key"] != INTEGRATED_ET_KEY]
        fig, axes = plt.subplots(len(et_keys), len(CENT_BINS), figsize=(16, 12), dpi=200, sharex=True, sharey=True)
        title = r"Inclusive background: BDT-score shape in isolation regions, split by $E_T^{cluster}$"
    else:
        et_keys = [et_key]
        fig, axes = plt.subplots(1, len(CENT_BINS), figsize=(W / DPI, H / DPI), dpi=DPI, sharex=True, sharey=True)
        axes = np.asarray([axes])
        title = r"Inclusive background: BDT score changes with isolation region"

    for row_idx, current_et_key in enumerate(et_keys):
        et_label = payload["panels"][current_et_key]["label"]
        for col_idx, (_, _, cent_label, cent_key) in enumerate(CENT_BINS):
            ax = axes[row_idx, col_idx]
            cell = payload["panels"][current_et_key]["centrality"][cent_key]
            for band in payload["eiso_bands"]:
                prof = cell["eiso_band_profiles"][band["label"]]
                y = safe_positive(np.asarray(prof["density"], dtype="float64"))
                label = band["label"] if row_idx == 0 and col_idx == 0 else None
                ax.step(x, y, where="mid", color=band["color"], linewidth=2.0, label=label)
            ax.set_yscale("log")
            ax.set_xlim(0.0, 1.0)
            ax.set_ylim(1.0e-2, 35.0)
            ax.grid(True, color=GRID, linewidth=0.7, alpha=0.85)
            if row_idx == 0:
                ax.set_title(cent_label, fontsize=18 if not grid else 13, fontweight="bold", pad=8)
            if col_idx == 0:
                ax.set_ylabel((et_label + "\n") + "unit-area density", fontsize=14 if grid else 17)
            if row_idx == len(et_keys) - 1:
                ax.set_xlabel("BDT score", fontsize=14 if grid else 17)
            ax.tick_params(labelsize=11 if grid else 14)
            add_common_panel_text(ax, cell)

    fig.suptitle(title, fontsize=22 if not grid else 19, fontweight="bold", y=0.985 if not grid else 0.972)
    subtitle = (
        r"Unit-area $p(\mathrm{BDT}\mid E_T^{iso}$ region); curve motion tests ABCD independence."
        if not grid
        else r"Control view: curve motion with isolation is checked separately inside each $E_T^{cluster}$ slice."
    )
    fig.text(0.5, 0.925 if not grid else 0.935, subtitle, ha="center", va="center", fontsize=14 if not grid else 12, color=MUTED)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=5, frameon=False, fontsize=13 if not grid else 10, bbox_to_anchor=(0.5, 0.018))
    fig.tight_layout(rect=[0.045, 0.075, 0.99, 0.885 if not grid else 0.905])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=DPI)
    plt.close(fig)


def plot_class_split_bdt_given_eiso(npz_path: Path, summary_path: Path, out_path: Path) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import FancyBboxPatch

    setup_style()

    with np.load(npz_path, allow_pickle=True) as data:
        hist = data["hist"].astype("float64", copy=False)
        eiso_bins = data["x_bins"].astype("float64", copy=False)
        score_bins = data["y_bins"].astype("float64", copy=False)
        class_labels = [str(v) for v in data["class_labels"]]
        centrality_labels = [str(v) for v in data["centrality_labels"]]
    summary = json.loads(summary_path.read_text())

    score_centers = 0.5 * (score_bins[:-1] + score_bins[1:])
    score_widths = np.diff(score_bins)

    def band_counts(class_idx: int, cent_idx: int, lo: float, hi: float) -> np.ndarray:
        centers = 0.5 * (eiso_bins[:-1] + eiso_bins[1:])
        mask = (centers >= lo) & (centers < hi)
        return np.sum(hist[class_idx, cent_idx, mask, :], axis=0)

    def density_from_counts(counts: np.ndarray) -> np.ndarray:
        total = float(np.sum(counts))
        if total <= 0.0:
            return np.full_like(counts, np.nan, dtype="float64")
        density = counts / (total * score_widths)
        return np.where(density > 0.0, density, np.nan)

    def display_eiso_label(band_label: str) -> str:
        label_map = {
            r"$E_T^{iso}<0$": r"$E_T^{iso}<0$ GeV",
            r"$0<E_T^{iso}<3$": r"$0<E_T^{iso}<3$ GeV",
            r"$3<E_T^{iso}<6$": r"$3<E_T^{iso}<6$ GeV",
            r"$6<E_T^{iso}<10$": r"$6<E_T^{iso}<10$ GeV",
            r"$E_T^{iso}>10$": r"$E_T^{iso}>10$ GeV",
        }
        return label_map.get(band_label, band_label)

    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)
    fig.patch.set_facecolor(WHITE)
    fig.text(
        0.055,
        0.962,
        "Isolation slicing exposes class-dependent BDT behavior",
        ha="left",
        va="top",
        fontsize=29,
        fontweight="bold",
        color=INK,
    )
    fig.text(
        0.055,
        0.912,
        r"Truth photons and inclusive background are shown in separate rows; curves are unit-area $p(\mathrm{BDT})$ in raw isolation regions.",
        ha="left",
        va="top",
        fontsize=14.2,
        color=MUTED,
    )
    for y_rule in (0.858, 0.790):
        fig.add_artist(
            Line2D(
                [0.180, 0.885],
                [y_rule, y_rule],
                transform=fig.transFigure,
                color="#cbd3dc",
                linewidth=1.1,
                solid_capstyle="butt",
            )
        )
    legend_items = [(display_eiso_label(label), color) for _, _, label, color in EISO_BANDS]
    legend_positions = [(0.255, 0.835), (0.505, 0.835), (0.755, 0.835), (0.385, 0.803), (0.620, 0.803)]
    for (legend_label, legend_color), (legend_x, legend_y) in zip(legend_items, legend_positions):
        fig.add_artist(
            Line2D(
                [legend_x - 0.080, legend_x - 0.032],
                [legend_y, legend_y],
                transform=fig.transFigure,
                color=legend_color,
                linewidth=2.6,
                solid_capstyle="butt",
            )
        )
        fig.text(legend_x - 0.020, legend_y, legend_label, ha="left", va="center", fontsize=13.4, color=INK)

    class_specs = [
        ("truth photons", "truth photons", "#b12632", 0.535),
        ("inclusive background", "inclusive background", "#1764a8", 0.245),
    ]
    axes = []
    for row_idx, (class_key, row_label, row_color, y0) in enumerate(class_specs):
        class_idx = class_labels.index(class_key)
        fig.text(
            0.038,
            y0 + 0.125,
            row_label,
            ha="center",
            va="center",
            rotation=90,
            fontsize=13.8,
            fontweight="bold",
            color=row_color,
        )
        for cent_idx, cent_label in enumerate(centrality_labels):
            ax = fig.add_axes([0.085 + cent_idx * 0.300, y0, 0.250, 0.215])
            axes.append(ax)
            for lo, hi, band_label, color in EISO_BANDS:
                counts = band_counts(class_idx, cent_idx, lo, hi)
                ax.step(score_centers, density_from_counts(counts), where="mid", color=color, linewidth=1.85)
            ax.set_yscale("log")
            ax.set_xlim(0.0, 1.0)
            ax.set_ylim(8.0e-3, 45.0)
            ax.grid(True, color=GRID, linewidth=0.65, alpha=0.85)
            if row_idx == 0:
                ax.set_title(cent_label + " centrality", fontsize=15.0, fontweight="bold", pad=6)
                ax.tick_params(labelbottom=False)
            else:
                ax.set_xlabel("BDT score", fontsize=13.2, labelpad=4)
            if cent_idx == 0:
                ax.set_ylabel("")
            else:
                ax.tick_params(labelleft=False)
            ax.tick_params(labelsize=10.6)
            stats = summary["summaries"][class_key][cent_label]
            stat_text = (
                f"N={int(stats['n']):,}\n"
                + f"median BDT={float(stats['median_score']):.3f}\n"
                + f"Pearson={float(stats['weighted_pearson_score_eiso']):.2f}"
            )
            ax.text(
                0.035,
                0.060,
                stat_text,
                transform=ax.transAxes,
                ha="left",
                va="bottom",
                fontsize=8.8,
                color=INK,
                bbox=dict(
                    boxstyle="round,pad=0.22",
                    facecolor="#ffffff",
                    edgecolor="#d4dce6",
                    linewidth=0.7,
                    alpha=0.94,
                ),
            )

    left_card = FancyBboxPatch(
        (0.075, 0.045),
        0.405,
        0.082,
        boxstyle="round,pad=0.010,rounding_size=0.008",
        transform=fig.transFigure,
        facecolor="#fffaf1",
        edgecolor="#bf7a11",
        linewidth=1.0,
        zorder=0.1,
    )
    right_card = FancyBboxPatch(
        (0.520, 0.045),
        0.405,
        0.082,
        boxstyle="round,pad=0.010,rounding_size=0.008",
        transform=fig.transFigure,
        facecolor="#f8fafc",
        edgecolor="#66707e",
        linewidth=1.0,
        zorder=0.1,
    )
    fig.add_artist(left_card)
    fig.add_artist(right_card)
    fig.text(0.092, 0.112, "Signal row:", ha="left", va="top", fontsize=12.6, fontweight="bold", color="#a96d12")
    fig.text(0.175, 0.112, "truth photons remain high-score across isolation slices.", ha="left", va="top", fontsize=12.6, color=INK)
    fig.text(0.092, 0.077, "Background row:", ha="left", va="top", fontsize=12.6, fontweight="bold", color="#a96d12")
    fig.text(0.205, 0.077, "positive isolation pulls the BDT shape downward.", ha="left", va="top", fontsize=12.6, color=INK)
    fig.text(0.537, 0.112, "ABCD concern:", ha="left", va="top", fontsize=12.6, fontweight="bold", color="#5e6876")
    fig.text(0.645, 0.112, "background isolation and BDT score should be independent.", ha="left", va="top", fontsize=12.6, color=INK)
    fig.text(
        0.537,
        0.077,
        "The visible background-row motion is the closure test target.",
        ha="left",
        va="top",
        fontsize=12.6,
        color=INK,
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=DPI)
    plt.close(fig)


def plot_high_bdt_fraction(payload: dict[str, Any], out_path: Path, et_key: str = INTEGRATED_ET_KEY, *, grid: bool = False) -> None:
    import matplotlib.pyplot as plt

    setup_style()
    cone = cone_radius_text(payload)
    colors = {"BDT > 0.8": "#1764a8", "BDT > 0.9": "#16827a", "BDT > WP80": "#a32035"}
    if grid:
        et_keys = [row["key"] for row in payload["et_bins"] if row["key"] != INTEGRATED_ET_KEY]
        fig, axes = plt.subplots(len(et_keys), len(CENT_BINS), figsize=(16, 12), dpi=200, sharex=True, sharey=True)
        title = r"Inclusive background high-BDT fraction versus isolation, split by $E_T^{cluster}$"
    else:
        et_keys = [et_key]
        fig, axes = plt.subplots(1, len(CENT_BINS), figsize=(W / DPI, H / DPI), dpi=DPI, sharex=True, sharey=True)
        axes = np.asarray([axes])
        title = r"Inclusive background high-BDT fraction versus isolation"

    for row_idx, current_et_key in enumerate(et_keys):
        et_label = payload["panels"][current_et_key]["label"]
        for col_idx, (_, _, cent_label, cent_key) in enumerate(CENT_BINS):
            ax = axes[row_idx, col_idx]
            cell = payload["panels"][current_et_key]["centrality"][cent_key]
            rows = cell["high_bdt_fraction_vs_eiso"]
            centers = np.asarray([row["center"] for row in rows], dtype="float64")
            xerr = np.asarray([(row["eiso_hi"] - row["eiso_lo"]) / 2.0 for row in rows], dtype="float64")
            for label, color in colors.items():
                y = np.asarray([np.nan if row["thresholds"][label]["fraction"] is None else row["thresholds"][label]["fraction"] for row in rows])
                yerr = np.asarray([np.nan if row["thresholds"][label]["error"] is None else row["thresholds"][label]["error"] for row in rows])
                plot_label = label if row_idx == 0 and col_idx == 0 else None
                ax.errorbar(centers, y, yerr=yerr, xerr=xerr, color=color, marker="o", linewidth=2.0, markersize=4.8, capsize=2.5, label=plot_label)
            ax.set_yscale("log")
            ax.set_ylim(5.0e-5, 1.05)
            ax.set_xlim(FRACTION_EISO_BINS[0], FRACTION_EISO_BINS[-1])
            ax.grid(True, color=GRID, linewidth=0.7, alpha=0.85)
            ax.axvline(0.0, color="#8a8f99", linewidth=1.0, linestyle="--", alpha=0.75)
            if row_idx == 0:
                ax.set_title(cent_label, fontsize=18 if not grid else 13, fontweight="bold", pad=8)
            if col_idx == 0:
                ax.set_ylabel((et_label + "\n") + "background fraction", fontsize=14 if grid else 17)
            if row_idx == len(et_keys) - 1:
                ax.set_xlabel(rf"reco $E_T^{{iso}}$ bin, $\Delta R<{cone}$ [GeV]", fontsize=14 if grid else 17)
            ax.tick_params(labelsize=11 if grid else 14)
            add_common_panel_text(ax, cell)

    fig.suptitle(title, fontsize=22 if not grid else 19, fontweight="bold", y=0.985 if not grid else 0.972)
    subtitle = (
        "Flat curves support independence; falling curves mean high-BDT background is sculpted toward lower isolation."
        if not grid
        else "Control view: the same fraction test repeated in cluster-energy slices."
    )
    fig.text(0.5, 0.925 if not grid else 0.935, subtitle, ha="center", va="center", fontsize=14 if not grid else 12, color=MUTED)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False, fontsize=13 if not grid else 10, bbox_to_anchor=(0.5, 0.018))
    fig.tight_layout(rect=[0.045, 0.075, 0.99, 0.885 if not grid else 0.905])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=DPI)
    plt.close(fig)


def write_metrics(payload: dict[str, Any], csv_path: Path) -> None:
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "et_bin",
        "centrality",
        "raw_n",
        "weighted_n",
        "wp80_signal_threshold_same_cache",
        "median_eiso",
        "median_score",
        "spearman_score_eiso",
        "weighted_pearson_score_eiso",
        "eiso_js_0p2_0p4_vs_0p95_1p0",
        "eiso_js_0p2_0p4_vs_highest_populated_bdt",
        "eiso_js_high_reference_bdt_band",
        "bdt_js_eiso_lt0_vs_gt10",
        "slope_bdt_gt_0p8_per_gev",
        "slope_bdt_gt_0p9_per_gev",
        "slope_bdt_gt_wp80_per_gev",
    ]
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for et in payload["et_bins"]:
            for cent in payload["centrality_bins"]:
                cell = payload["panels"][et["key"]]["centrality"][cent["key"]]
                slopes = cell["fraction_slope_per_gev"]
                writer.writerow(
                    {
                        "et_bin": et["label"],
                        "centrality": cent["label"],
                        "raw_n": cell["raw_n"],
                        "weighted_n": cell["weighted_n"],
                        "wp80_signal_threshold_same_cache": cell["wp80_signal_threshold_same_cache"],
                        "median_eiso": cell["median_eiso"],
                        "median_score": cell["median_score"],
                        "spearman_score_eiso": cell["spearman_score_eiso"],
                        "weighted_pearson_score_eiso": cell["weighted_pearson_score_eiso"],
                        "eiso_js_0p2_0p4_vs_0p95_1p0": cell["eiso_js_0p2_0p4_vs_0p95_1p0"],
                        "eiso_js_0p2_0p4_vs_highest_populated_bdt": cell[
                            "eiso_js_0p2_0p4_vs_highest_populated_bdt"
                        ],
                        "eiso_js_high_reference_bdt_band": cell["eiso_js_high_reference_bdt_band"],
                        "bdt_js_eiso_lt0_vs_gt10": cell["bdt_js_eiso_lt0_vs_gt10"],
                        "slope_bdt_gt_0p8_per_gev": slopes["BDT > 0.8"],
                        "slope_bdt_gt_0p9_per_gev": slopes["BDT > 0.9"],
                        "slope_bdt_gt_wp80_per_gev": slopes["BDT > WP80"],
                    }
                )


def write_interpretation(payload: dict[str, Any], note_path: Path, outputs: dict[str, str]) -> None:
    def fmt(value: Any, digits: int = 3) -> str:
        if value is None:
            return "n/a"
        try:
            value = float(value)
        except (TypeError, ValueError):
            return "n/a"
        if not math.isfinite(value):
            return "n/a"
        return f"{value:.{digits}f}"

    def clean_values(values: list[Any]) -> list[float]:
        cleaned = []
        for value in values:
            try:
                value = float(value)
            except (TypeError, ValueError):
                continue
            if math.isfinite(value):
                cleaned.append(value)
        return cleaned

    def mean_or_nan(values: list[Any]) -> float:
        cleaned = clean_values(values)
        return float(np.mean(cleaned)) if cleaned else float("nan")

    integrated = payload["panels"][INTEGRATED_ET_KEY]["centrality"]
    lines = [
        "# THE-11 background-only conditional BDT/isolation profiles",
        "",
        "## What was made",
        "",
        "- Sample: inclusive/background candidates only from the existing source-preserving full-stat Au+Au validation score caches.",
        "- Selection: `15 <= cluster_Et < 35 GeV`, centrality `0-20%`, `20-50%`, `50-80%`, and "
        + f"`{payload['metadata']['iso_column']}`.",
        "- Weights: current `event_weight` from the score caches.",
        "- WP80 threshold: same-cache weighted 20th percentile of truth-photon BDT scores in each centrality and ET cell.",
        "- Normalization: every 1D curve is unit-area normalized inside its own panel.",
        "",
        "## Integrated 15-35 GeV read",
        "",
    ]
    for _, _, cent_label, cent_key in CENT_BINS:
        cell = integrated[cent_key]
        slopes = cell["fraction_slope_per_gev"]
        lines.append(
            "- "
            + f"{cent_label}: N={cell['raw_n']:,}, median Eiso={fmt(cell['median_eiso'], 2)} GeV, "
            + f"median BDT={fmt(cell['median_score'])}, Spearman={fmt(cell['spearman_score_eiso'])}, "
            + "JS[p(Eiso|0.2-0.4), "
            + f"p(Eiso|{cell['eiso_js_high_reference_bdt_band']})]="
            + f"{fmt(cell['eiso_js_0p2_0p4_vs_highest_populated_bdt'])}, "
            + f"slope f(BDT>WP80)={fmt(slopes['BDT > WP80'], 4)} per GeV."
        )
    rho_values = [integrated[cent_key]["spearman_score_eiso"] for _, _, _, cent_key in CENT_BINS]
    js_values = [integrated[cent_key]["eiso_js_0p2_0p4_vs_highest_populated_bdt"] for _, _, _, cent_key in CENT_BINS]
    slope_values = [integrated[cent_key]["fraction_slope_per_gev"]["BDT > WP80"] for _, _, _, cent_key in CENT_BINS]
    mean_rho = mean_or_nan(rho_values)
    mean_js = mean_or_nan(js_values)
    mean_slope = mean_or_nan(slope_values)
    if mean_rho < -0.15 or mean_js > 0.15 or mean_slope < -0.002:
        conclusion = (
            "The background BDT score and isolation are not independent in this sample: "
            "high-BDT background is preferentially concentrated at lower isolation, "
            "while positive-isolation background shifts toward lower BDT score. "
            "That is exactly the ABCD-relevant closure issue: the method can still be usable, "
            "but it needs explicit closure or transfer-factor validation."
        )
    else:
        conclusion = (
            "The integrated profiles show only a weak conditional motion between BDT score and isolation. "
            "That would be closer to the ABCD independence target, but the ET-sliced controls still need to be checked."
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            textwrap.fill(conclusion, width=100),
            "",
            "## Output files",
            "",
        ]
    )
    for label, path in outputs.items():
        lines.append(f"- {label}: `{path}`")
    note_path.write_text("\n".join(lines) + "\n")


def plot_all(payload_path: Path, out_dir: Path) -> dict[str, str]:
    payload = load_payload(payload_path)
    out_dir.mkdir(parents=True, exist_ok=True)
    class_split_npz = out_dir.parent / "slideReady" / "weighted_bdt_iso_fullstat_summary_r30.npz"
    class_split_json = out_dir.parent / "slideReady" / "weighted_bdt_iso_fullstat_summary_r30.json"
    outputs = {
        "p(Eiso | BDT band), integrated ET": str(out_dir / "the11_background_eiso_given_bdt_integrated_cent.png"),
        "p(BDT | Eiso band), integrated ET": str(out_dir / "the11_background_bdt_given_eiso_integrated_cent.png"),
        "p(BDT | Eiso band), signal/background rows": str(out_dir / "the11_signal_background_bdt_given_eiso_integrated_cent.png"),
        "high-BDT fraction vs Eiso, integrated ET": str(out_dir / "the11_background_highbdt_fraction_vs_eiso_integrated_cent.png"),
        "p(Eiso | BDT band), ET control grid": str(out_dir / "the11_background_eiso_given_bdt_et_control_grid.png"),
        "p(BDT | Eiso band), ET control grid": str(out_dir / "the11_background_bdt_given_eiso_et_control_grid.png"),
        "high-BDT fraction vs Eiso, ET control grid": str(out_dir / "the11_background_highbdt_fraction_vs_eiso_et_control_grid.png"),
        "metrics CSV": str(out_dir / "the11_background_conditional_profiles_metrics.csv"),
        "interpretation note": str(out_dir / "the11_background_conditional_profiles_interpretation.md"),
    }
    plot_eiso_given_bdt(payload, Path(outputs["p(Eiso | BDT band), integrated ET"]))
    plot_bdt_given_eiso(payload, Path(outputs["p(BDT | Eiso band), integrated ET"]))
    if class_split_npz.exists() and class_split_json.exists():
        plot_class_split_bdt_given_eiso(
            class_split_npz,
            class_split_json,
            Path(outputs["p(BDT | Eiso band), signal/background rows"]),
        )
    plot_high_bdt_fraction(payload, Path(outputs["high-BDT fraction vs Eiso, integrated ET"]))
    plot_eiso_given_bdt(payload, Path(outputs["p(Eiso | BDT band), ET control grid"]), grid=True)
    plot_bdt_given_eiso(payload, Path(outputs["p(BDT | Eiso band), ET control grid"]), grid=True)
    plot_high_bdt_fraction(payload, Path(outputs["high-BDT fraction vs Eiso, ET control grid"]), grid=True)
    write_metrics(payload, Path(outputs["metrics CSV"]))
    write_interpretation(payload, Path(outputs["interpretation note"]), outputs)
    return outputs


def write_manifest(payload_path: Path, out_dir: Path, outputs: dict[str, str]) -> Path:
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": str(Path(__file__).resolve()),
        "payload": str(payload_path),
        "outputs": outputs,
        "google_slides_mutation": False,
        "qa_expectation": "Visual inspect primary PNGs before calling slide-ready.",
    }
    manifest_path = out_dir / "the11_background_conditional_profiles_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest_path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=["reduce", "plot"], default="plot")
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REMOTE_REPORT)
    parser.add_argument("--cache-manifest", type=Path, default=None)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--payload", type=Path, default=None)
    parser.add_argument("--iso-column", choices=[ISO_COLUMN_R30, ISO_COLUMN_R40], default=ISO_COLUMN_R30)
    args = parser.parse_args()

    if args.mode == "reduce":
        payload = reduce_score_caches(args.report_dir, args.cache_manifest, args.iso_column)
        print("BEGIN_THE11_BACKGROUND_PROFILE_JSON")
        print(json.dumps(payload, sort_keys=True))
        print("END_THE11_BACKGROUND_PROFILE_JSON")
        return 0

    payload_path = args.payload or (args.out_dir / "the11_background_conditional_profiles_summary_r30.json")
    outputs = plot_all(payload_path, args.out_dir)
    manifest_path = write_manifest(payload_path, args.out_dir, outputs)
    print(json.dumps({"manifest": str(manifest_path), "outputs": outputs}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
