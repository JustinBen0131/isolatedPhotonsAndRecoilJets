#!/usr/bin/env python3
"""Overlay PPG12 Fig. 7 trigger diagnostics from SDCC and current pp output."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import beta
import uproot


CAMPAIGN = "the76_ppg12_parity_full_20260701_003024"
DEFAULT_CURRENT_ROOT = Path(
    "dataOutput/ppg12Parity"
    f"/{CAMPAIGN}/final_roots/pp_data_hierarchical_v1/"
    "RecoilJets_pp_ALL_PERIOD_COMBINED.root"
)
DEFAULT_SDCC_CSV = Path(
    "dataOutput/ppg12PhotonYield/ppg12_photon_yield_v1_data_20260620/"
    "shower_shape_reference_validation/fig7_trigger/"
    "ppg12_fig7_sdcc_trigger_arrays.csv"
)
DEFAULT_SDCC_BIT30_PLATEAU_CSV = Path(
    "dataOutput/ppg12PhotonYield/ppg12_photon_yield_v1_data_20260620/"
    "shower_shape_reference_validation/fig7_trigger/"
    "ppg12_fig7_sdcc_bit30_plateau_trigger_root_asymm_from_numden.csv"
)
DEFAULT_OUTDIR = Path(
    "dataOutput/ppg12Parity"
    f"/{CAMPAIGN}/data_trigger_fig7"
)

CUR_DIR = "PPG12_trigger_fig7"
SPECTRA = [
    ("scaled10", "scaled[10]", "h_ppg12_fig7_cluster_et_scaled10", "den_scaled10", "black", "o"),
    ("live29", "scaled[10] & live[29]", "h_ppg12_fig7_cluster_et_scaled10_live29", "num_live29", "#1f77b4", "s"),
    ("live30", "scaled[10] & live[30]", "h_ppg12_fig7_cluster_et_scaled10_live30", "num_live30", "#d62728", "o"),
    ("live31", "scaled[10] & live[31]", "h_ppg12_fig7_cluster_et_scaled10_live31", "num_live31", "#2ca02c", "^"),
]


def read_sdcc_rows(path: Path) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    with path.open() as handle:
        for row in csv.DictReader(handle):
            rows.append({key: float(value) for key, value in row.items()})
    if not rows:
        raise RuntimeError(f"No rows read from {path}")
    return rows


def hist_values(root_path: Path, hist_name: str) -> tuple[np.ndarray, np.ndarray]:
    with uproot.open(root_path) as handle:
        hist = handle[f"{CUR_DIR}/{hist_name}"]
        values, edges = hist.to_numpy(flow=False)
    return values.astype(float), edges.astype(float)


def hist_values_first(
    root_path: Path,
    candidates: list[str],
) -> tuple[str, np.ndarray, np.ndarray]:
    last_error: Exception | None = None
    for hist_name in candidates:
        try:
            values, edges = hist_values(root_path, hist_name)
            return hist_name, values, edges
        except uproot.exceptions.KeyInFileError as exc:
            last_error = exc
    names = ", ".join(f"{CUR_DIR}/{name}" for name in candidates)
    raise RuntimeError(f"None of the candidate histograms exists: {names}") from last_error


def values_on_reference_edges(
    values: np.ndarray,
    edges: np.ndarray,
    reference_rows: list[dict[str, float]],
    *,
    label: str,
) -> np.ndarray:
    selected: list[float] = []
    for row in reference_rows:
        lo = row["xlow"]
        hi = row["xhigh"]
        matches = np.where(np.isclose(edges[:-1], lo) & np.isclose(edges[1:], hi))[0]
        if len(matches) != 1:
            raise RuntimeError(
                f"{label}: expected exactly one current bin matching SDCC edge "
                f"[{lo}, {hi}], found {len(matches)}"
            )
        selected.append(float(values[matches[0]]))
    return np.asarray(selected, dtype=float)


def poisson_err(y: np.ndarray) -> np.ndarray:
    return np.sqrt(np.clip(y, 0.0, None))


def ratio_err(a: np.ndarray, ea: np.ndarray, b: np.ndarray, eb: np.ndarray) -> np.ndarray:
    ratio = np.divide(a, b, out=np.full_like(a, np.nan, dtype=float), where=b > 0)
    rel2 = np.zeros_like(a, dtype=float)
    rel2 += np.divide(ea, a, out=np.zeros_like(a, dtype=float), where=a > 0) ** 2
    rel2 += np.divide(eb, b, out=np.zeros_like(b, dtype=float), where=b > 0) ** 2
    return ratio * np.sqrt(rel2)


def ratio_err_asymm(
    a: np.ndarray,
    ea_low: np.ndarray,
    ea_high: np.ndarray,
    b: np.ndarray,
    eb_low: np.ndarray,
    eb_high: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    ratio = np.divide(a, b, out=np.full_like(a, np.nan, dtype=float), where=b > 0)
    rel_low2 = np.divide(ea_low, a, out=np.zeros_like(a, dtype=float), where=a > 0) ** 2
    rel_low2 += np.divide(eb_high, b, out=np.zeros_like(b, dtype=float), where=b > 0) ** 2
    rel_high2 = np.divide(ea_high, a, out=np.zeros_like(a, dtype=float), where=a > 0) ** 2
    rel_high2 += np.divide(eb_low, b, out=np.zeros_like(b, dtype=float), where=b > 0) ** 2
    return ratio * np.sqrt(rel_low2), ratio * np.sqrt(rel_high2)


def binomial_eff(pass_counts: np.ndarray, total_counts: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    eff = np.divide(
        pass_counts,
        total_counts,
        out=np.full_like(pass_counts, np.nan, dtype=float),
        where=total_counts > 0,
    )
    err = np.sqrt(np.divide(eff * (1.0 - eff), total_counts, out=np.zeros_like(eff), where=total_counts > 0))
    return eff, err


def binomial_eff_cp(
    pass_counts: np.ndarray,
    total_counts: np.ndarray,
    confidence: float = 0.682689492137,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Match ROOT TEfficiency::kFCP 1-sigma Clopper-Pearson intervals."""
    eff = np.divide(
        pass_counts,
        total_counts,
        out=np.full_like(pass_counts, np.nan, dtype=float),
        where=total_counts > 0,
    )
    alpha = 1.0 - confidence
    low = np.full_like(eff, np.nan, dtype=float)
    high = np.full_like(eff, np.nan, dtype=float)
    for i, (passed, total) in enumerate(zip(pass_counts, total_counts)):
        k = int(round(float(passed)))
        n = int(round(float(total)))
        if n <= 0:
            continue
        low[i] = 0.0 if k == 0 else beta.ppf(alpha / 2.0, k, n - k + 1)
        high[i] = 1.0 if k == n else beta.ppf(1.0 - alpha / 2.0, k + 1, n - k)
    err_low = np.maximum(0.0, eff - low)
    err_high = np.maximum(0.0, high - eff)
    err_sym = 0.5 * (err_low + err_high)
    return eff, err_low, err_high, err_sym


def finish_axes(ax_top, ax_bot, xlabel: str, ratio_ylabel: str, ymin: float | None = None, ymax: float | None = None) -> None:
    for ax in (ax_top, ax_bot):
        ax.tick_params(which="both", direction="in", top=True, right=True, length=6)
        ax.tick_params(which="minor", length=3)
        ax.minorticks_on()
    ax_bot.axhline(1.0, color="0.35", lw=1.1, ls=(0, (4, 4)))
    ax_bot.set_xlabel(xlabel, fontsize=18)
    ax_bot.set_ylabel(ratio_ylabel, fontsize=15)
    if ymin is not None and ymax is not None:
        ax_bot.set_ylim(ymin, ymax)
    ax_top.set_xlim(5.0, 15.0)
    ax_bot.set_xlim(5.0, 15.0)


def draw_header(ax, text: str, *, x: float = 0.07, y: float = 0.92) -> None:
    box = {"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 1.5}
    ax.text(x, y, r"$\it{\bf{sPHENIX}}$ Internal", transform=ax.transAxes, fontsize=15.5, bbox=box)
    ax.text(x, y - 0.075, r"$p{+}p\ \sqrt{s}=200\ \mathrm{GeV},\ 64.37\ \mathrm{pb}^{-1}$",
            transform=ax.transAxes, fontsize=12.5, bbox=box)
    ax.text(x, y - 0.145, r"$|\eta^\gamma|<0.7$", transform=ax.transAxes, fontsize=12.5, bbox=box)
    ax.text(x, y - 0.215, text, transform=ax.transAxes, fontsize=11.5, bbox=box)


def plot_shape_overlay(root_path: Path, sdcc_rows: list[dict[str, float]], outdir: Path) -> dict:
    x = np.array([row["xcenter"] for row in sdcc_rows], dtype=float)
    width = np.array([row["width"] for row in sdcc_rows], dtype=float)
    fig, (ax, rax) = plt.subplots(
        2,
        1,
        figsize=(7.2, 7.8),
        dpi=180,
        sharex=True,
        gridspec_kw={"height_ratios": [3.15, 1.0], "hspace": 0.035},
    )

    csv_rows: list[dict[str, float | str]] = []
    summary: dict[str, dict[str, float]] = {}
    for key, label, hist_name, sdcc_col, color, marker in SPECTRA:
        cur_values, cur_edges = hist_values(root_path, hist_name)
        if not np.allclose(cur_edges[: len(x) + 1], [sdcc_rows[0]["xlow"], *[row["xhigh"] for row in sdcc_rows]]):
            raise RuntimeError(f"Binning mismatch for {hist_name}")
        current = cur_values[: len(x)]
        sdcc = np.array([row[sdcc_col] for row in sdcc_rows], dtype=float)
        current_density = current / width
        sdcc_density = sdcc / width
        current_norm = current_density / np.sum(current_density * width)
        sdcc_norm = sdcc_density / np.sum(sdcc_density * width)
        current_err = poisson_err(current) / width / np.sum(current_density * width)
        sdcc_err = poisson_err(sdcc) / width / np.sum(sdcc_density * width)
        ratio = np.divide(current_norm, sdcc_norm, out=np.full_like(current_norm, np.nan), where=sdcc_norm > 0)
        ratio_unc = ratio_err(current_norm, current_err, sdcc_norm, sdcc_err)

        ax.errorbar(x, sdcc_norm, yerr=sdcc_err, xerr=width / 2.0, fmt=marker, ms=4.6,
                    color=color, mfc="white", mec=color, lw=1.0, capsize=0, alpha=0.95)
        ax.errorbar(x, current_norm, yerr=current_err, xerr=width / 2.0, fmt=marker, ms=4.1,
                    color=color, mfc=color, mec=color, lw=1.0, capsize=0, alpha=0.9)
        rax.errorbar(x, ratio, yerr=ratio_unc, xerr=width / 2.0, fmt=marker, ms=4.2,
                     color=color, mfc=color, mec=color, lw=1.0, capsize=0)

        finite = np.isfinite(ratio)
        summary[key] = {
            "current_counts": float(np.sum(current)),
            "sdcc_counts": float(np.sum(sdcc)),
            "max_abs_ratio_deviation": float(np.nanmax(np.abs(ratio[finite] - 1.0))) if np.any(finite) else math.nan,
        }
        for i in range(len(x)):
            csv_rows.append({
                "curve": key,
                "xlow": sdcc_rows[i]["xlow"],
                "xhigh": sdcc_rows[i]["xhigh"],
                "xcenter": x[i],
                "sdcc_norm_density": sdcc_norm[i],
                "current_norm_density": current_norm[i],
                "current_over_sdcc": ratio[i],
                "current_over_sdcc_err": ratio_unc[i],
                "sdcc_counts": sdcc[i],
                "current_counts": current[i],
            })

    ax.set_yscale("log")
    ax.set_ylabel("shape-normalized clusters / GeV", fontsize=16)
    ax.set_ylim(2e-4, 3.5)
    draw_header(ax, "PPG12 Fig. 7 trigger spectra", x=0.07, y=0.31)
    source_handles = [
        plt.Line2D([0], [0], marker="o", color="black", mfc="white", mec="black", lw=0, label="PPG12 SDCC"),
        plt.Line2D([0], [0], marker="o", color="black", mfc="black", mec="black", lw=0, label="Current pp output"),
    ]
    curve_handles = [
        plt.Line2D([0], [0], marker=marker, color=color, mfc=color, mec=color, lw=0, label=label)
        for _, label, _, _, color, marker in SPECTRA
    ]
    leg1 = ax.legend(handles=source_handles, loc="upper left", bbox_to_anchor=(1.02, 1.0), fontsize=11.5,
                     frameon=False, title="source", title_fontsize=12, borderaxespad=0.0)
    ax.add_artist(leg1)
    ax.legend(handles=curve_handles, loc="upper left", bbox_to_anchor=(1.02, 0.74), fontsize=10.5,
              frameon=False, title="trigger sample", title_fontsize=11.5, borderaxespad=0.0)
    finish_axes(ax, rax, r"Leading $E_T^{cluster}$ [GeV]", "Current / SDCC", 0.0, 5.6)

    out_png = outdir / "fig7_trigger_spectra_sdcc_vs_current_shape_ratio.png"
    out_csv = outdir / "fig7_trigger_spectra_sdcc_vs_current_shape_ratio.csv"
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)
    with out_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(csv_rows[0].keys()))
        writer.writeheader()
        writer.writerows(csv_rows)
    return {"png": str(out_png), "csv": str(out_csv), "summary": summary}


def plot_efficiency_overlay(root_path: Path, sdcc_rows: list[dict[str, float]], outdir: Path) -> dict:
    x = np.array([row["xcenter"] for row in sdcc_rows], dtype=float)
    width = np.array([row["width"] for row in sdcc_rows], dtype=float)
    sdcc_den = np.array([row["den_scaled10"] for row in sdcc_rows], dtype=float)
    current_den, _ = hist_values(root_path, "h_ppg12_fig7_cluster_et_scaled10")
    current_den = current_den[: len(x)]

    fig, (ax, rax) = plt.subplots(
        2,
        1,
        figsize=(7.2, 7.2),
        dpi=180,
        sharex=True,
        gridspec_kw={"height_ratios": [2.75, 1.0], "hspace": 0.04},
    )
    csv_rows: list[dict[str, float | str]] = []
    summary: dict[str, dict[str, float]] = {}
    for key, label, hist_name, sdcc_col, color, marker in SPECTRA[1:]:
        current_pass, _ = hist_values(root_path, hist_name)
        current_pass = current_pass[: len(x)]
        sdcc_pass = np.array([row[sdcc_col] for row in sdcc_rows], dtype=float)
        cur_eff, cur_err_low, cur_err_high, cur_err = binomial_eff_cp(current_pass, current_den)
        sdcc_eff, sdcc_err_low, sdcc_err_high, sdcc_err = binomial_eff_cp(sdcc_pass, sdcc_den)
        ratio = np.divide(cur_eff, sdcc_eff, out=np.full_like(cur_eff, np.nan), where=sdcc_eff > 0)
        ratio_err_low, ratio_err_high = ratio_err_asymm(
            cur_eff,
            cur_err_low,
            cur_err_high,
            sdcc_eff,
            sdcc_err_low,
            sdcc_err_high,
        )

        ax.errorbar(x, sdcc_eff, yerr=[sdcc_err_low, sdcc_err_high], xerr=width / 2.0, fmt=marker, ms=4.8,
                    color=color, mfc="white", mec=color, lw=1.0, capsize=0)
        ax.errorbar(x, cur_eff, yerr=[cur_err_low, cur_err_high], xerr=width / 2.0, fmt=marker, ms=4.2,
                    color=color, mfc=color, mec=color, lw=1.0, capsize=0, alpha=0.9)
        rax.errorbar(x, ratio, yerr=[ratio_err_low, ratio_err_high], xerr=width / 2.0, fmt=marker, ms=4.2,
                     color=color, mfc=color, mec=color, lw=1.0, capsize=0)
        finite = np.isfinite(ratio)
        summary[key] = {
            "max_abs_ratio_deviation": float(np.nanmax(np.abs(ratio[finite] - 1.0))) if np.any(finite) else math.nan,
            "median_ratio": float(np.nanmedian(ratio[finite])) if np.any(finite) else math.nan,
        }
        for i in range(len(x)):
            csv_rows.append({
                "curve": key,
                "xlow": sdcc_rows[i]["xlow"],
                "xhigh": sdcc_rows[i]["xhigh"],
                "xcenter": x[i],
                "sdcc_eff": sdcc_eff[i],
                "sdcc_eff_err_low": sdcc_err_low[i],
                "sdcc_eff_err_high": sdcc_err_high[i],
                "sdcc_eff_err_sym": sdcc_err[i],
                "current_eff": cur_eff[i],
                "current_eff_err_low": cur_err_low[i],
                "current_eff_err_high": cur_err_high[i],
                "current_eff_err_sym": cur_err[i],
                "current_over_sdcc": ratio[i],
                "current_over_sdcc_err_low": ratio_err_low[i],
                "current_over_sdcc_err_high": ratio_err_high[i],
                "sdcc_den": sdcc_den[i],
                "current_den": current_den[i],
                "statistical_model": "TEfficiency::kFCP Clopper-Pearson 68.268949%; ratio errors propagated asymmetrically",
            })

    ax.set_ylabel(r"Trigger efficiency vs scaled[10]", fontsize=16)
    ax.set_ylim(0.55, 1.08)
    draw_header(ax, "live-bit efficiency, same Fig. 7 bins", x=0.07, y=0.38)
    source_handles = [
        plt.Line2D([0], [0], marker="o", color="black", mfc="white", mec="black", lw=0, label="PPG12 SDCC"),
        plt.Line2D([0], [0], marker="o", color="black", mfc="black", mec="black", lw=0, label="Current pp output"),
    ]
    curve_handles = [
        plt.Line2D([0], [0], marker=marker, color=color, mfc=color, mec=color, lw=0, label=label.replace("scaled[10] & ", ""))
        for _, label, _, _, color, marker in SPECTRA[1:]
    ]
    leg1 = ax.legend(handles=curve_handles, loc="upper left", bbox_to_anchor=(1.02, 0.76), fontsize=12,
                     frameon=False, title="numerator", title_fontsize=12, borderaxespad=0.0)
    ax.add_artist(leg1)
    ax.legend(handles=source_handles, loc="upper left", bbox_to_anchor=(1.02, 0.43), fontsize=12,
              frameon=False, title="source", title_fontsize=12, borderaxespad=0.0)
    finish_axes(ax, rax, r"Leading $E_T^{cluster}$ [GeV]", "Current / SDCC", 0.82, 1.12)
    out_png = outdir / "fig7_trigger_efficiency_sdcc_vs_current_ratio.png"
    out_csv = outdir / "fig7_trigger_efficiency_sdcc_vs_current_ratio.csv"
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)
    with out_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(csv_rows[0].keys()))
        writer.writeheader()
        writer.writerows(csv_rows)
    return {"png": str(out_png), "csv": str(out_csv), "summary": summary}


def plot_bit30_plateau_overlay(root_path: Path, sdcc_rows: list[dict[str, float]], outdir: Path) -> dict:
    x = np.array([row["xcenter"] for row in sdcc_rows], dtype=float)
    sdcc_eff = np.array([row["eff"] for row in sdcc_rows], dtype=float)
    sdcc_cp_err_low = np.array(
        [row["eylow"] if "eylow" in row else row["eff_error"] for row in sdcc_rows],
        dtype=float,
    )
    sdcc_cp_err_high = np.array(
        [row["eyhigh"] if "eyhigh" in row else row["eff_error"] for row in sdcc_rows],
        dtype=float,
    )
    sdcc_cp_err_sym = 0.5 * (sdcc_cp_err_low + sdcc_cp_err_high)
    sdcc_th1_err = np.array(
        [row["th1_error"] if "th1_error" in row else 0.5 * (lo + hi)
         for row, lo, hi in zip(sdcc_rows, sdcc_cp_err_low, sdcc_cp_err_high)],
        dtype=float,
    )

    total_name, current_total_all, current_edges = hist_values_first(
        root_path,
        [
            "h_ppg12_fig7_eff_total_live10_bit30",
            "h_ppg12_fig7_eff_total_bit30",
        ],
    )
    pass_name, current_pass_all, pass_edges = hist_values_first(
        root_path,
        [
            "h_ppg12_fig7_eff_pass_live10_bit30",
            "h_ppg12_fig7_eff_pass_bit30",
        ],
    )
    if ("live10" in total_name) != ("live10" in pass_name):
        raise RuntimeError(
            "Mixed live[10] and scaled[10] bit30 pass/total histogram families: "
            f"{total_name}, {pass_name}"
        )
    current_contract = (
        "live[10] denominator (PPG12 Fig.7 turn-on contract)"
        if "live10" in total_name
        else "scaled[10] denominator (legacy diagnostic fallback; not exact Fig.7 turn-on)"
    )
    current_contract_short = (
        "live[10] denominator\n(PPG12 Fig.7 turn-on contract)"
        if "live10" in total_name
        else "scaled[10] denominator\n(legacy fallback; not exact Fig.7)"
    )
    if not np.allclose(current_edges, pass_edges):
        raise RuntimeError("Current bit30 pass/total histogram edges do not match")

    current_total = values_on_reference_edges(
        current_total_all, current_edges, sdcc_rows, label=total_name
    )
    current_pass = values_on_reference_edges(
        current_pass_all, pass_edges, sdcc_rows, label=pass_name
    )
    current_eff, current_err_low, current_err_high, current_err = binomial_eff_cp(current_pass, current_total)
    ratio = np.divide(current_eff, sdcc_eff, out=np.full_like(current_eff, np.nan), where=sdcc_eff > 0)
    ratio_err_low, ratio_err_high = ratio_err_asymm(
        current_eff,
        current_err_low,
        current_err_high,
        sdcc_eff,
        sdcc_cp_err_low,
        sdcc_cp_err_high,
    )

    finite = np.isfinite(ratio)
    max_idx = int(np.nanargmax(np.abs(ratio[finite] - 1.0))) if np.any(finite) else -1
    finite_indices = np.where(finite)[0]
    max_global_idx = int(finite_indices[max_idx]) if max_idx >= 0 else -1
    max_dev = float(abs(ratio[max_global_idx] - 1.0)) if max_global_idx >= 0 else math.nan

    fig, (ax, rax) = plt.subplots(
        2,
        1,
        figsize=(6.2, 6.7),
        dpi=180,
        sharex=True,
        gridspec_kw={"height_ratios": [2.85, 1.0], "hspace": 0.04},
    )

    red = "#f05a3a"
    current_color = "black"
    ax.errorbar(
        x,
        sdcc_eff,
        yerr=[sdcc_cp_err_low, sdcc_cp_err_high],
        fmt="D",
        ms=3.4,
        color=red,
        mfc=red,
        mec=red,
        lw=1.0,
        capsize=2.0,
        capthick=1.0,
        label="PPG12 SDCC",
    )
    ax.errorbar(
        x,
        current_eff,
        yerr=[current_err_low, current_err_high],
        fmt="o",
        ms=4.2,
        color=current_color,
        mfc=current_color,
        mec=current_color,
        lw=1.0,
        capsize=2.0,
        capthick=1.0,
        alpha=0.92,
        label="Current pp output",
    )
    fit_x = np.linspace(7.0, 30.0, 300)
    fit_y = 0.996634 * np.exp(-np.exp(-(fit_x - 4.10489) / 0.666567))
    ax.plot(fit_x, fit_y, color=red, lw=1.4, alpha=0.65, label="PPG12 Gumbel fit")
    rax.errorbar(
        x,
        ratio,
        yerr=[ratio_err_low, ratio_err_high],
        fmt="o",
        ms=4.0,
        color=current_color,
        mfc=current_color,
        mec=current_color,
        lw=1.0,
        capsize=2.0,
        capthick=1.0,
    )

    for panel in (ax, rax):
        panel.tick_params(which="both", direction="in", top=True, right=True, length=6)
        panel.tick_params(which="minor", length=3)
        panel.minorticks_on()

    ax.set_xlim(7.0, 30.0)
    ax.set_ylim(0.985, 1.005)
    rax.set_ylim(0.84, 1.08)
    rax.set_xlim(7.0, 30.0)
    rax.axhline(1.0, color="0.35", lw=1.1, ls=(0, (4, 4)))
    ax.set_ylabel("Photon 4 GeV / MBD N&S efficiency", fontsize=14)
    rax.set_ylabel("Current / SDCC", fontsize=13)
    rax.set_xlabel(r"Leading $E_T^{cluster}$ [GeV]", fontsize=16)

    box = {"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 1.0}
    ax.text(0.06, 0.255, r"$\it{\bf{sPHENIX}}$ Internal", transform=ax.transAxes, fontsize=12.0, bbox=box)
    ax.text(0.06, 0.185, r"$p{+}p\ \sqrt{s}=200\ \mathrm{GeV},\ 64.37\ \mathrm{pb}^{-1}$",
            transform=ax.transAxes, fontsize=9.6, bbox=box)
    ax.text(0.06, 0.125, r"$|\eta^\gamma|<0.7$", transform=ax.transAxes, fontsize=9.6, bbox=box)
    ax.text(0.06, 0.065, r"Photon 4 GeV trigger bit 30", transform=ax.transAxes, fontsize=9.0, bbox=box)
    ax.legend(
        loc="upper left",
        bbox_to_anchor=(0.08, 0.98),
        ncol=1,
        frameon=True,
        fontsize=10.4,
        title="Photon 4 GeV bit30 efficiency",
        title_fontsize=10.4,
        facecolor="white",
        edgecolor="white",
        framealpha=0.78,
        borderaxespad=0.4,
        handlelength=1.8,
        labelspacing=0.36,
    )

    csv_rows: list[dict[str, float]] = []
    for i, row in enumerate(sdcc_rows):
        csv_rows.append({
            "xlow": row["xlow"],
            "xhigh": row["xhigh"],
            "xcenter": x[i],
            "width": row["width"],
            "sdcc_eff": sdcc_eff[i],
            "sdcc_eff_err_low": sdcc_cp_err_low[i],
            "sdcc_eff_err_high": sdcc_cp_err_high[i],
            "sdcc_eff_err_sym": sdcc_cp_err_sym[i],
            "sdcc_cp_err_low": sdcc_cp_err_low[i],
            "sdcc_cp_err_high": sdcc_cp_err_high[i],
            "sdcc_cp_err_sym": sdcc_cp_err_sym[i],
            "sdcc_th1_error": sdcc_th1_err[i],
            "current_pass": current_pass[i],
            "current_total": current_total[i],
            "current_eff": current_eff[i],
            "current_eff_err_low": current_err_low[i],
            "current_eff_err_high": current_err_high[i],
            "current_eff_err_sym": current_err[i],
            "current_cp_err_low": current_err_low[i],
            "current_cp_err_high": current_err_high[i],
            "current_cp_err_sym": current_err[i],
            "current_over_sdcc": ratio[i],
            "current_over_sdcc_err_low": ratio_err_low[i],
            "current_over_sdcc_err_high": ratio_err_high[i],
            "current_total_hist": total_name,
            "current_pass_hist": pass_name,
            "current_contract": current_contract,
            "statistical_model": "Displayed bars match the PPG12 IAN Fig.7 plateau visual convention: TEfficiency::kFCP / Clopper-Pearson asymmetric intervals; TH1/binomial error is retained separately.",
        })

    out_png = outdir / "fig7_bit30_plateau_sdcc_vs_current_ratio.png"
    out_csv = outdir / "fig7_bit30_plateau_sdcc_vs_current_ratio.csv"
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)
    with out_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(csv_rows[0].keys()))
        writer.writeheader()
        writer.writerows(csv_rows)
    return {
        "png": str(out_png),
        "csv": str(out_csv),
        "summary": {
            "sdcc_rows": len(sdcc_rows),
            "current_total": float(current_total.sum()),
            "current_pass": float(current_pass.sum()),
            "max_abs_ratio_deviation": max_dev,
            "max_deviation_xcenter": float(x[max_global_idx]) if max_global_idx >= 0 else math.nan,
            "current_total_hist": total_name,
            "current_pass_hist": pass_name,
            "current_contract": current_contract,
            "display_statistical_model": "PPG12 IAN Fig.7 plateau visual convention: TEfficiency::kFCP Clopper-Pearson asymmetric intervals",
            "alternate_statistical_model": "TH1/binomial symmetric error is retained in CSV as sdcc_th1_error but is not the displayed IAN-style convention",
            "display_axis_range": "Top pad y range 0.985-1.005 to match the PPG12 IAN plateau zoom.",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--current-root", type=Path, default=DEFAULT_CURRENT_ROOT)
    parser.add_argument("--sdcc-csv", type=Path, default=DEFAULT_SDCC_CSV)
    parser.add_argument("--sdcc-bit30-plateau-csv", type=Path, default=DEFAULT_SDCC_BIT30_PLATEAU_CSV)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    sdcc_rows = read_sdcc_rows(args.sdcc_csv)
    sdcc_bit30_plateau_rows = read_sdcc_rows(args.sdcc_bit30_plateau_csv)
    shape = plot_shape_overlay(args.current_root, sdcc_rows, args.outdir)
    eff = plot_efficiency_overlay(args.current_root, sdcc_rows, args.outdir)
    bit30_plateau = plot_bit30_plateau_overlay(args.current_root, sdcc_bit30_plateau_rows, args.outdir)
    manifest = {
        "current_root": str(args.current_root.resolve()),
        "sdcc_csv": str(args.sdcc_csv.resolve()),
        "sdcc_bit30_plateau_csv": str(args.sdcc_bit30_plateau_csv.resolve()),
        "output_dir": str(args.outdir.resolve()),
        "comparison": "PPG12 Fig.7 SDCC source arrays vs current pp RecoilJets output",
        "normalization": {
            "spectra": "Each curve shape-normalized to unit area over shared SDCC bins before ratio.",
            "efficiency": "Per-bin live bit / scaled[10]; no exposure normalization.",
            "bit30_plateau": "Per-bin Photon 4 GeV bit30 trigger efficiency. Uses live[10] denominator when the current ROOT contains the exact PPG12 turn-on objects; otherwise falls back to the older scaled[10] diagnostic and labels that exception.",
        },
        "statistical_errors": {
            "ppg12_source": "ppg12codeGit/plotting/plot_trigger_bit30.C sets TEfficiency::kFCP for bit29/30/31. Drawn efficiency points use Clopper-Pearson 68.268949% asymmetric intervals.",
            "ppg12_fits": "PPG12 converts TEfficiency low/up errors to a symmetric TH1 fit error as 0.5*(err_low+err_high), with 1/sqrt(total) only as a zero-error fallback.",
            "fig7_multibit_overlay": "The multibit Fig.7 efficiency overlay uses Clopper-Pearson intervals for the TEfficiency-style live[29/30/31] cross-check and propagates low/up intervals asymmetrically in the ratio.",
            "bit30_plateau_display": "The RHS bit30 plateau overlay follows the PPG12 IAN Fig.7 visual: asymmetric Clopper-Pearson intervals from the extracted eylow/eyhigh columns. The TH1/binomial column is retained in the CSV only as an alternate audit value.",
        },
        "artifacts": {
            "shape_overlay": shape,
            "efficiency_overlay": eff,
            "bit30_plateau_overlay": bit30_plateau,
        },
    }
    manifest_path = args.outdir / "fig7_trigger_sdcc_vs_current_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
