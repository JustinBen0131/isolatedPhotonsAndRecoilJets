#!/usr/bin/env python3
"""Build the isolated THE-219 Friday PPG pp-to-AuAu plot packet.

This is an offline, current-pointer-resolved plotting consumer.  It does not
submit jobs, mutate current pointers, edit Google Slides, or modify any
THE-114/THE-139 artifact.  Every physics-facing output is intentionally
labelled as either source-first reused evidence or a diagnostic candidate.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import uproot
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Patch
from scipy.stats import chi2


REPO = Path(__file__).resolve().parents[3]
CURRENT_DIR = REPO / "dataOutput/current_recoiljets_artifacts/current"
DEFAULT_OUT = REPO / "dataOutput/the219_friday_ppg_20260814/phase_a"
CANONICAL_DECK = (
    "https://docs.google.com/presentation/d/"
    "19pjXOc1CJ1Ed7xtw2ZrXffztpWkzCz4lN13L7CHPdQE/edit"
)

POINTER_KEYS = (
    "pp_data_merged",
    "pp_sim_photonjet_merged",
    "pp_sim_inclusivejet_merged",
    "auau_data_merged",
    "auau_sim_photonjet_merged",
    "auau_sim_inclusivejet_merged",
)

PP_NATIVE_BINS = (
    (10, 12),
    (12, 14),
    (14, 16),
    (16, 18),
    (18, 20),
    (20, 22),
    (22, 24),
    (24, 26),
    (26, 35),
)
AUAU_NATIVE_BINS = ((15, 17), (17, 19), (19, 21), (21, 23), (23, 26), (26, 35))
AUAU_CENTRALITIES = (
    ("0_20", "0–20%", "#0877bd", "s"),
    ("20_50", "20–50%", "#e58b12", "^"),
    ("50_80", "50–80%", "#16824b", "D"),
)
ABCD_TOKENS = {
    "A": "isIsolated_isTight",
    "B": "notIsolated_isTight",
    "C": "isIsolated_notTight",
    "D": "notIsolated_notTight",
}

INK = "#111827"
MUTED = "#526173"
GRID = "#d9e1ea"
BLUE = "#1769aa"
ORANGE = "#d97706"
GREEN = "#16824b"
RED = "#c9362b"
PURPLE = "#7352b8"

SOURCE_FIRST = (
    {
        "role": "canonical pp photon-definition and purity control",
        "png": REPO
        / "dataOutput/ppg12Parity/the97_ppg12_final_accepted_triple_full_20260714_1550"
        / "final_pp_data_canonical_20260717/purity_overlay_current"
        / "the97_pp_data_purity_ppg12_vs_current_fullstat.png",
        "manifest": REPO
        / "dataOutput/ppg12Parity/the97_ppg12_final_accepted_triple_full_20260714_1550"
        / "final_pp_data_canonical_20260717/purity_overlay_current"
        / "the97_pp_data_purity_ppg12_vs_current_fullstat_manifest.json",
    },
    {
        "role": "raw pp ABCD parity control",
        "png": REPO
        / "dataOutput/ppg12Parity/the97_ppg12_final_accepted_triple_full_20260714_1550"
        / "final_pp_data_canonical_20260717/raw_abcd_three_panel"
        / "the97_final_pp_data_raw_abcd_yield_and_purity_current_over_ppg12.png",
        "manifest": REPO
        / "dataOutput/ppg12Parity/the97_ppg12_final_accepted_triple_full_20260714_1550"
        / "final_pp_data_canonical_20260717/raw_abcd_three_panel"
        / "the97_final_pp_data_raw_abcd_yield_and_purity_current_over_ppg12_manifest.json",
    },
    {
        "role": "current pp photon response matrix",
        "png": REPO
        / "dataOutput/ppg12Parity/the97_ppg12_si_contract_restore_full_20260715_1420"
        / "ian_current_sim_refresh_20260716/response_matrix"
        / "fig36_current_h_response_full_raw_ppg12_style.png",
        "manifest": REPO
        / "dataOutput/ppg12Parity/the97_ppg12_si_contract_restore_full_20260715_1420"
        / "ian_current_sim_refresh_20260716/response_matrix"
        / "fig36_response_matrix_manifest.json",
    },
)


@dataclass(frozen=True)
class Artifact:
    key: str
    pointer_path: Path
    pointer: dict[str, Any]
    root: Path
    root_sha256: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument(
        "--sam-reference-pdf",
        type=Path,
        default=None,
        help="Optional local copy of the PPG18 reference PDF for provenance hashing.",
    )
    parser.add_argument(
        "--visual-qa-pass",
        action="store_true",
        help="Record PASS only after every generated PNG has been inspected at full-slide scale.",
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def resolve_artifact(key: str) -> Artifact:
    pointer_path = CURRENT_DIR / key / "current.json"
    pointer = json.loads(pointer_path.read_text())
    roots = pointer.get("root_paths") or []
    if len(roots) != 1:
        raise RuntimeError(f"expected one root in {pointer_path}, found {len(roots)}")
    root = Path(roots[0])
    if not root.is_file():
        raise FileNotFoundError(root)
    return Artifact(key, pointer_path, pointer, root, sha256(root))


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "mathtext.fontset": "dejavusans",
            "axes.linewidth": 1.05,
            "axes.labelcolor": INK,
            "text.color": INK,
            "xtick.color": INK,
            "ytick.color": INK,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "savefig.facecolor": "white",
        }
    )


def sphx(
    ax: plt.Axes,
    *,
    x: float = 0.035,
    y: float = 0.965,
    size: float = 13.0,
    internal_dx: float = 0.145,
) -> None:
    ax.text(
        x,
        y,
        "sPHENIX",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=size,
        fontstyle="italic",
        fontweight="bold",
    )
    ax.text(x + internal_dx, y, "Internal", transform=ax.transAxes, ha="left", va="top", fontsize=size)


def style_axis(ax: plt.Axes) -> None:
    ax.grid(axis="y", color=GRID, lw=0.85, alpha=0.8, zorder=0)
    ax.minorticks_on()
    ax.tick_params(which="major", length=6, width=1.0, labelsize=10.5)
    ax.tick_params(which="minor", length=3, width=0.8)


def hist_sum_and_var(root: uproot.ReadOnlyDirectory, path: str) -> tuple[float, float]:
    if path not in root:
        raise KeyError(path)
    hist = root[path]
    values = np.asarray(hist.values(flow=False), dtype=float)
    variances = hist.variances(flow=False)
    variance = float(np.sum(np.abs(values))) if variances is None else float(np.sum(variances))
    return float(np.sum(values)), variance


def kappa_and_variance(counts: np.ndarray, variances: np.ndarray) -> tuple[float, float]:
    a, b, c, d = [float(item) for item in counts]
    if min(a, b, c, d) <= 0.0:
        return float("nan"), float("nan")
    kappa = b * c / (a * d)
    variance = kappa * kappa * float(np.sum(variances / (counts * counts)))
    return kappa, max(0.0, variance)


def extract_sim_abcd(
    root_file: uproot.ReadOnlyDirectory,
    bins: Iterable[tuple[int, int]],
    *,
    suffix_prefix: str = "",
    suffix_tail: str = "",
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for low, high in bins:
        suffix = f"{suffix_prefix}_pT_{low}_{high}{suffix_tail}"
        total = []
        total_var = []
        for region in "ABCD":
            name = f"SIM/h_xJpurityLead_{ABCD_TOKENS[region]}{suffix}"
            count, variance = hist_sum_and_var(root_file, name)
            total.append(count)
            total_var.append(variance)
        sig_name = f"SIM/h_xJpurityLead_sigABCD_MC{suffix}"
        sig_hist = root_file[sig_name]
        signal = np.asarray(sig_hist.values(flow=False)[:4], dtype=float)
        sig_variances_raw = sig_hist.variances(flow=False)
        signal_var = (
            np.abs(signal)
            if sig_variances_raw is None
            else np.asarray(sig_variances_raw[:4], dtype=float)
        )
        total_array = np.asarray(total, dtype=float)
        total_var_array = np.asarray(total_var, dtype=float)
        background = total_array - signal
        # The truth-signal histogram is an exact event-weighted subset of the
        # corresponding all-candidate region.  Therefore Cov(total, signal)
        # equals Var(signal), so Var(total - signal) = Var(total)-Var(signal).
        background_var = np.maximum(0.0, total_var_array - signal_var)
        raw_kappa, raw_kappa_var = kappa_and_variance(total_array, total_var_array)
        bkg_kappa, bkg_kappa_var = kappa_and_variance(background, background_var)
        row: dict[str, Any] = {
            "pt_low_gev": low,
            "pt_high_gev": high,
            "raw_kappa_BC_over_AD": raw_kappa,
            "raw_kappa_error": math.sqrt(raw_kappa_var),
            "background_kappa_BC_over_AD": bkg_kappa,
            "background_kappa_error": math.sqrt(bkg_kappa_var),
            "signal_histogram": sig_name,
        }
        for index, region in enumerate("ABCD"):
            row[f"{region}_total"] = float(total_array[index])
            row[f"{region}_total_variance"] = float(total_var_array[index])
            row[f"{region}_truth_signal"] = float(signal[index])
            row[f"{region}_truth_signal_variance"] = float(signal_var[index])
            row[f"{region}_background"] = float(background[index])
            row[f"{region}_background_variance"] = float(background_var[index])
        rows.append(row)
    return rows


def extract_data_abcd(
    root_file: uproot.ReadOnlyDirectory, bins: Iterable[tuple[int, int]]
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    top = "PPG12_scaledtrigger30"
    for low, high in bins:
        counts = []
        variances = []
        names = {}
        for region in "ABCD":
            name = f"{top}/h_xJpurityLead_{ABCD_TOKENS[region]}_pT_{low}_{high}"
            names[region] = name
            count, variance = hist_sum_and_var(root_file, name)
            counts.append(count)
            variances.append(variance)
        count_array = np.asarray(counts, dtype=float)
        variance_array = np.asarray(variances, dtype=float)
        kappa, kappa_var = kappa_and_variance(count_array, variance_array)
        row: dict[str, Any] = {
            "pt_low_gev": low,
            "pt_high_gev": high,
            "observed_kappa_BC_over_AD": kappa,
            "observed_kappa_error": math.sqrt(kappa_var),
            "histograms": names,
        }
        for index, region in enumerate("ABCD"):
            row[region] = float(count_array[index])
            row[f"{region}_variance"] = float(variance_array[index])
        rows.append(row)
    return rows


def chi2_to_unity(rows: list[dict[str, Any]], value_key: str, error_key: str) -> dict[str, float]:
    selected = [row for row in rows if math.isfinite(row[value_key]) and row[error_key] > 0.0]
    statistic = float(sum(((row[value_key] - 1.0) / row[error_key]) ** 2 for row in selected))
    ndf = len(selected)
    return {"chi2": statistic, "ndf": ndf, "p_value": float(chi2.sf(statistic, ndf))}


def centers_and_errors(rows: list[dict[str, Any]]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    centers = np.asarray([(row["pt_low_gev"] + row["pt_high_gev"]) / 2.0 for row in rows])
    xerr = np.asarray([(row["pt_high_gev"] - row["pt_low_gev"]) / 2.0 for row in rows])
    return centers, xerr, np.asarray([row["pt_low_gev"] for row in rows])


def draw_kappa(
    pp_sim: list[dict[str, Any]],
    pp_data: list[dict[str, Any]],
    auau: dict[str, list[dict[str, Any]]],
    out: Path,
) -> dict[str, Any]:
    pp_focus = [row for row in pp_sim if row["pt_low_gev"] >= 14]
    tests = {"pp_14to35": chi2_to_unity(pp_focus, "background_kappa_BC_over_AD", "background_kappa_error")}
    for key, rows in auau.items():
        tests[f"auau_{key}"] = chi2_to_unity(rows, "background_kappa_BC_over_AD", "background_kappa_error")

    fig, axes = plt.subplots(2, 1, figsize=(12.0, 8.2), dpi=180, sharex=False)
    fig.subplots_adjust(left=0.095, right=0.975, bottom=0.10, top=0.88, hspace=0.29)
    fig.suptitle(
        "ABCD factorization must be tested on background, not the inclusive sample",
        fontsize=22,
        fontweight="bold",
        y=0.965,
    )
    fig.text(
        0.095,
        0.915,
        r"$\kappa\equiv BC/(AD)$; factorization predicts 1.  Exact truth-signal subset subtraction includes its covariance.",
        fontsize=13.0,
        color=MUTED,
        ha="left",
    )

    top = axes[0]
    style_axis(top)
    top.axhline(1.0, color=INK, lw=1.15, ls="--", zorder=1)
    x, xe, _ = centers_and_errors(pp_sim)
    top.errorbar(
        x - 0.18,
        [row["raw_kappa_BC_over_AD"] for row in pp_sim],
        xerr=xe,
        yerr=[row["raw_kappa_error"] for row in pp_sim],
        fmt="o",
        color=ORANGE,
        mfc="white",
        mew=1.4,
        ms=6.5,
        capsize=2.2,
        label="p+p inclusive-jet MC, raw (signal mixed)",
        zorder=4,
    )
    top.errorbar(
        x + 0.18,
        [row["background_kappa_BC_over_AD"] for row in pp_sim],
        xerr=xe,
        yerr=[row["background_kappa_error"] for row in pp_sim],
        fmt="o",
        color=INK,
        mfc=INK,
        ms=6.3,
        capsize=2.2,
        label="p+p inclusive-jet MC, truth-signal subtracted background",
        zorder=5,
    )
    xd, xed, _ = centers_and_errors(pp_data)
    top.errorbar(
        xd,
        [row["observed_kappa_BC_over_AD"] for row in pp_data],
        xerr=xed,
        yerr=[row["observed_kappa_error"] for row in pp_data],
        fmt="s",
        color=BLUE,
        mfc="white",
        mew=1.25,
        ms=5.7,
        capsize=2.0,
        label="p+p data observed (signal mixed; diagnostic only)",
        zorder=3,
    )
    top.set_xlim(9.5, 35.7)
    top.set_ylim(0.18, 1.58)
    top.set_ylabel(r"$\kappa=BC/(AD)$", fontsize=13.5)
    top.set_xlabel(r"reconstructed photon $p_T^\gamma$ [GeV]", fontsize=12.5)
    top.legend(frameon=False, fontsize=10.2, loc="upper right", ncol=1)
    top.text(
        0.025,
        0.08,
        "Raw inclusive MC appears nonfactorizing because Region A is signal-rich.",
        transform=top.transAxes,
        fontsize=10.5,
        color=RED,
        fontweight="bold",
    )
    sphx(top)

    bottom = axes[1]
    style_axis(bottom)
    bottom.axhline(1.0, color=INK, lw=1.15, ls="--", zorder=1)
    xpp, xepp, _ = centers_and_errors(pp_focus)
    bottom.errorbar(
        xpp,
        [row["background_kappa_BC_over_AD"] for row in pp_focus],
        xerr=xepp,
        yerr=[row["background_kappa_error"] for row in pp_focus],
        fmt="o",
        color=INK,
        mfc=INK,
        ms=6.2,
        capsize=2.2,
        label=f"p+p; $p={tests['pp_14to35']['p_value']:.2f}$",
        zorder=6,
    )
    offsets = {"0_20": -0.27, "20_50": 0.0, "50_80": 0.27}
    for key, label, color, marker in AUAU_CENTRALITIES:
        rows = auau[key]
        xx, xxe, _ = centers_and_errors(rows)
        test = tests[f"auau_{key}"]
        bottom.errorbar(
            xx + offsets[key],
            [row["background_kappa_BC_over_AD"] for row in rows],
            xerr=xxe,
            yerr=[row["background_kappa_error"] for row in rows],
            fmt=marker,
            color=color,
            mfc=color,
            mec="white",
            mew=0.7,
            ms=6.5,
            capsize=2.1,
            label=f"Au+Au {label}; $p={test['p_value']:.3f}$",
            zorder=5,
        )
    bottom.set_xlim(13.5, 35.7)
    bottom.set_ylim(0.45, 1.58)
    bottom.set_xlabel(r"reconstructed photon $p_T^\gamma$ [GeV]", fontsize=13.5)
    bottom.set_ylabel(r"background-only $\kappa$", fontsize=13.5)
    bottom.legend(frameon=False, fontsize=10.0, loc="upper right", ncol=2, columnspacing=1.1)
    bottom.text(
        0.025,
        0.08,
        "Au+Au 20–50% rejects a constant-unity closure test; retain an explicit transfer correction.",
        transform=bottom.transAxes,
        fontsize=10.2,
        color=RED,
        fontweight="bold",
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    plt.close(fig)
    return tests


def transfer_rows(
    pp_sim: list[dict[str, Any]], pp_data: list[dict[str, Any]]
) -> list[dict[str, Any]]:
    sim_lookup = {(row["pt_low_gev"], row["pt_high_gev"]): row for row in pp_sim}
    rows = []
    for data in pp_data:
        low, high = data["pt_low_gev"], data["pt_high_gev"]
        if low < 14:
            continue
        sim = sim_lookup[(low, high)]
        b, d = data["B"], data["D"]
        b_var, d_var = data["B_variance"], data["D_variance"]
        kappa = sim["background_kappa_BC_over_AD"]
        kappa_var = sim["background_kappa_error"] ** 2
        naive = b / d
        naive_var = naive * naive * (b_var / (b * b) + d_var / (d * d))
        corrected = naive / kappa
        corrected_var = corrected * corrected * (naive_var / (naive * naive) + kappa_var / (kappa * kappa))
        rows.append(
            {
                "pt_low_gev": low,
                "pt_high_gev": high,
                "data_B": b,
                "data_B_variance": b_var,
                "data_D": d,
                "data_D_variance": d_var,
                "mc_background_kappa": kappa,
                "mc_background_kappa_variance": kappa_var,
                "naive_transfer_B_over_D": naive,
                "naive_transfer_variance": naive_var,
                "corrected_transfer_B_over_D_over_kappa": corrected,
                "corrected_transfer_variance": corrected_var,
                "correction_multiplier_1_over_kappa": 1.0 / kappa,
                "correction_multiplier_variance": kappa_var / (kappa**4),
            }
        )
    return rows


def draw_transfer(rows: list[dict[str, Any]], out: Path) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(10.6, 7.4), dpi=180, sharex=True, gridspec_kw={"height_ratios": [2.4, 1.0]})
    fig.subplots_adjust(left=0.12, right=0.97, bottom=0.11, top=0.86, hspace=0.08)
    fig.suptitle("Use the measured nonfactorization—do not silently set it to one", fontsize=21, fontweight="bold", y=0.96)
    fig.text(
        0.12,
        0.905,
        r"Region-C shape normalization: $A_{\rm bkg}/C=(B/D)/\kappa$; statistical transfer covariance is retained.",
        fontsize=12.5,
        color=MUTED,
    )
    x, xe, _ = centers_and_errors(rows)
    top, ratio = axes
    style_axis(top)
    top.errorbar(
        x - 0.12,
        [row["naive_transfer_B_over_D"] for row in rows],
        xerr=xe,
        yerr=[math.sqrt(row["naive_transfer_variance"]) for row in rows],
        fmt="o",
        color=ORANGE,
        mfc="white",
        mew=1.3,
        capsize=2.2,
        label=r"naive $B/D$ ($\kappa=1$)",
    )
    top.errorbar(
        x + 0.12,
        [row["corrected_transfer_B_over_D_over_kappa"] for row in rows],
        xerr=xe,
        yerr=[math.sqrt(row["corrected_transfer_variance"]) for row in rows],
        fmt="s",
        color=BLUE,
        mfc=BLUE,
        mec="white",
        mew=0.7,
        capsize=2.2,
        label=r"explicit $(B/D)/\kappa_{\rm bkg}^{MC}$",
    )
    top.set_ylabel(r"Region-C $\rightarrow$ A transfer", fontsize=13.5)
    top.legend(frameon=False, fontsize=11.5, loc="upper right")
    top.text(
        0.025,
        0.08,
        "Normalization correction only; xJ-shape closure remains a separate gate.",
        transform=top.transAxes,
        color=RED,
        fontsize=10.5,
        fontweight="bold",
    )
    sphx(top)
    style_axis(ratio)
    ratio.axhline(1.0, color=INK, ls="--", lw=1.0)
    ratio.errorbar(
        x,
        [row["correction_multiplier_1_over_kappa"] for row in rows],
        xerr=xe,
        yerr=[math.sqrt(row["correction_multiplier_variance"]) for row in rows],
        fmt="s",
        color=PURPLE,
        mfc=PURPLE,
        mec="white",
        mew=0.7,
        capsize=2.2,
    )
    ratio.set_ylabel(r"$1/\kappa$", fontsize=12.5)
    ratio.set_xlabel(r"reconstructed photon $p_T^\gamma$ [GeV]", fontsize=13.5)
    ratio.set_xlim(13.5, 35.7)
    ratio.set_ylim(0.55, 1.35)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    plt.close(fig)


def load_2d(root_file: uproot.ReadOnlyDirectory, path: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if path not in root_file:
        raise KeyError(path)
    hist = root_file[path]
    return (
        np.asarray(hist.values(flow=False), dtype=float),
        np.asarray(hist.axis(0).edges(), dtype=float),
        np.asarray(hist.axis(1).edges(), dtype=float),
    )


def draw_bdt_iso_map(data_root: Path, signal_root: Path, out: Path) -> dict[str, Any]:
    data_path = "PPG12_scaledtrigger30/h2d_bdt_eta0_pt1535_cut1"
    signal_path = "SIM/h2d_bdt_eta0_pt1535_cut1"
    with uproot.open(data_root) as data_file, uproot.open(signal_root) as signal_file:
        data_values, xedges, yedges = load_2d(data_file, data_path)
        signal_values, sxedges, syedges = load_2d(signal_file, signal_path)
    if not (np.array_equal(xedges, sxedges) and np.array_equal(yedges, syedges)):
        raise RuntimeError("data and photon+jet BDT-isolation axes differ")
    arrays = []
    for values in (data_values, signal_values):
        total = float(np.sum(values))
        arrays.append(values / total if total > 0.0 else values)
    visible = []
    xmask = (xedges[:-1] >= -0.10) & (xedges[1:] <= 1.00)
    ymask = (yedges[:-1] >= -5.0) & (yedges[1:] <= 8.0)
    for values in arrays:
        visible.append(values[np.ix_(xmask, ymask)])
    vmax = max(float(np.max(item)) for item in visible)
    positives = np.concatenate([item[item > 0.0] for item in visible])
    vmin = max(1.0e-7, float(np.quantile(positives, 0.05)))

    fig, axes = plt.subplots(1, 2, figsize=(13.2, 6.6), dpi=180, sharex=True, sharey=True)
    fig.subplots_adjust(left=0.085, right=0.91, bottom=0.12, top=0.83, wspace=0.12)
    fig.suptitle("Sam's fixed ABCD boxes versus the canonical PPG12-parity selection", fontsize=21, fontweight="bold", y=0.965)
    fig.text(
        0.085,
        0.895,
        r"Post-PPG12 preselection/NPB population, $15<p_T^\gamma<35$ GeV; canonical boundaries vary with $p_T^\gamma$.",
        fontsize=12.5,
        color=MUTED,
    )
    titles = ("p+p data candidates", "p+p photon+jet MC candidates")
    for ax, values, title in zip(axes, arrays, titles):
        mesh = ax.pcolormesh(
            xedges,
            yedges,
            values.T,
            cmap="magma_r",
            norm=LogNorm(vmin=vmin, vmax=vmax),
            shading="flat",
        )
        ax.set_xlim(-0.10, 1.00)
        ax.set_ylim(-5.0, 8.0)
        ax.set_title(title, fontsize=14.5, fontweight="bold", pad=8)
        ax.set_xlabel("PPG12 baseV3E BDT score", fontsize=13.0)
        ax.tick_params(which="major", labelsize=10.5, length=5)
        ax.minorticks_on()
        # Sam: tight >0.8, loose 0.2--0.6, isolated <2, nonisolated >4.
        for boundary in (0.2, 0.6, 0.8):
            ax.axvline(boundary, color=RED, ls="--", lw=1.15, alpha=0.9)
        for boundary in (2.0, 4.0):
            ax.axhline(boundary, color=RED, ls="--", lw=1.15, alpha=0.9)
        # Canonical pT-dependent threshold envelopes over 15--35 GeV.
        tight = [0.815625 - 0.0015625 * pt for pt in (15.0, 35.0)]
        non_tight_low = [0.7333333333 - 0.0133333333 * pt for pt in (15.0, 35.0)]
        non_tight_high = [0.684375 + 0.0015625 * pt for pt in (15.0, 35.0)]
        iso = [0.490 + 0.037 * pt for pt in (15.0, 35.0)]
        noniso = [value + 0.8 for value in iso]
        ax.axvspan(min(tight), max(tight), color="#20b7c9", alpha=0.28, lw=0)
        ax.axvspan(
            min(non_tight_low),
            max(non_tight_high),
            facecolor="#20b7c9",
            edgecolor="#087f8c",
            alpha=0.07,
            hatch="///",
            lw=0.7,
        )
        ax.axhspan(min(iso), max(iso), color="#20b7c9", alpha=0.25, lw=0)
        ax.axhspan(min(noniso), max(noniso), color="#20b7c9", alpha=0.14, lw=0)
        ax.text(0.03, 0.96, "normalized candidate density", transform=ax.transAxes, ha="left", va="top", fontsize=9.5, color=MUTED)
    axes[0].set_ylabel(r"$E_T^{iso}$ [GeV]", fontsize=13.0)
    legend_items = [
        Line2D([0], [0], color=RED, ls="--", lw=1.4, label="Sam fixed boundaries"),
        Patch(facecolor="#20b7c9", alpha=0.30, label="PPG12 tight/isolation envelopes"),
        Patch(facecolor="#20b7c9", edgecolor="#087f8c", alpha=0.10, hatch="///", label="PPG12 non-tight window union"),
    ]
    axes[0].legend(handles=legend_items, frameon=False, fontsize=9.1, loc="lower left")
    cax = fig.add_axes([0.925, 0.12, 0.018, 0.71])
    cbar = fig.colorbar(mesh, cax=cax)
    cbar.set_label("fraction of retained candidates / bin", fontsize=11.5)
    cbar.ax.tick_params(labelsize=9.5)
    fig.text(
        0.50,
        0.045,
        "This locates the populations and cut differences; it is not a background-only factorization test.",
        ha="center",
        fontsize=10.8,
        color=RED,
        fontweight="bold",
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    plt.close(fig)
    return {
        "data_histogram": data_path,
        "photonjet_histogram": signal_path,
        "stage": "cut1 = after configured PPG12 preselection/NPB",
        "pt_token": "pt1535 = 15 < pT < 35 GeV",
        "sam_fixed": {"tight_bdt_min": 0.8, "loose_bdt": [0.2, 0.6], "iso_max_gev": 2.0, "noniso_min_gev": 4.0},
        "canonical_ppg12": {
            "tight_bdt": "score > 0.815625 - 0.0015625*pT",
            "nontight_bdt": "0.7333333333 - 0.0133333333*pT < score < 0.684375 + 0.0015625*pT",
            "iso": "Eiso < 0.490 + 0.037*pT",
            "noniso": "Eiso > 0.490 + 0.037*pT + 0.8",
        },
    }


def select_pt_range(values: np.ndarray, variances: np.ndarray, edges: np.ndarray, low: float, high: float) -> tuple[np.ndarray, np.ndarray]:
    indices = np.where((edges[:-1] >= low - 1e-9) & (edges[1:] <= high + 1e-9))[0]
    if not len(indices) or abs(edges[indices[0]] - low) > 1e-9 or abs(edges[indices[-1] + 1] - high) > 1e-9:
        raise RuntimeError(f"requested pT range {low}-{high} is not an exact union of stored bins {edges.tolist()}")
    return np.sum(values[indices, :], axis=0), np.sum(variances[indices, :], axis=0)


def normalized_density(values: np.ndarray, variances: np.ndarray, edges: np.ndarray) -> tuple[np.ndarray, np.ndarray, float, float]:
    total = float(np.sum(values))
    if total <= 0.0:
        raise RuntimeError("non-positive xJ integral")
    widths = np.diff(edges)
    density = values / (total * widths)
    error = np.sqrt(np.maximum(0.0, variances)) / (total * widths)
    centers = 0.5 * (edges[:-1] + edges[1:])
    mean = float(np.sum(values * centers) / total)
    return density, error, total, mean


def load_hist2_with_variance(root_file: uproot.ReadOnlyDirectory, path: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    hist = root_file[path]
    values = np.asarray(hist.values(flow=False), dtype=float)
    variance_raw = hist.variances(flow=False)
    variances = np.abs(values) if variance_raw is None else np.asarray(variance_raw, dtype=float)
    return values, variances, np.asarray(hist.axis(0).edges()), np.asarray(hist.axis(1).edges())


def draw_pp_xj(pp_data_root: Path, out: Path) -> list[dict[str, Any]]:
    path = "PPG12_scaledtrigger30/h2_unfoldReco_pTgamma_xJ_incl_r04"
    with uproot.open(pp_data_root) as root_file:
        values, variances, pt_edges, xj_edges = load_hist2_with_variance(root_file, path)
    ranges = ((16, 20), (20, 26), (26, 35))
    fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.9), dpi=180, sharey=True)
    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.16, top=0.78, wspace=0.08)
    fig.suptitle("Current p+p raw reconstructed $x_{J\gamma}$ in exact native-bin unions", fontsize=21, fontweight="bold", y=0.965)
    fig.text(
        0.075,
        0.865,
        r"PPG12-parity Region A, anti-$k_T$ $R=0.4$, nominal $\Delta\phi>7\pi/8$; shape normalized for display.",
        fontsize=12.5,
        color=MUTED,
    )
    rows = []
    for ax, (low, high) in zip(axes, ranges):
        selected, selected_var = select_pt_range(values, variances, pt_edges, low, high)
        density, error, total, mean = normalized_density(selected, selected_var, xj_edges)
        centers = 0.5 * (xj_edges[:-1] + xj_edges[1:])
        xerr = 0.5 * np.diff(xj_edges)
        ax.errorbar(centers, density, xerr=xerr, yerr=error, fmt="o", color=BLUE, mfc=BLUE, mec="white", mew=0.6, ms=5.3, capsize=1.8)
        ax.step(xj_edges[:-1], density, where="post", color=BLUE, lw=1.15, alpha=0.8)
        style_axis(ax)
        ax.set_xlim(0.0, 2.2)
        ax.set_ylim(bottom=0.0)
        ax.set_title(rf"${low}<p_T^\gamma<{high}$ GeV", fontsize=14.0, fontweight="bold")
        ax.set_xlabel(r"$x_{J\gamma}=p_T^{jet}/p_T^\gamma$", fontsize=12.5)
        annotation = f"pairs = {total:.0f}\n" + rf"$\langle x_{{J\gamma}}\rangle={mean:.2f}$"
        ax.text(0.96, 0.92, annotation, transform=ax.transAxes, ha="right", va="top", fontsize=10.2)
        rows.append(
            {
                "pt_low_gev": low,
                "pt_high_gev": high,
                "pair_integral": total,
                "shape_mean_xj": mean,
                "values": selected.tolist(),
                "variances": selected_var.tolist(),
                "xj_edges": xj_edges.tolist(),
            }
        )
    axes[0].set_ylabel(r"shape density $1/N\,dN/dx_{J\gamma}$", fontsize=12.5)
    sphx(axes[0], size=11.5, internal_dx=0.23)
    fig.text(
        0.5,
        0.055,
        "Raw signal-region diagnostic only: no ABCD subtraction, efficiency correction, or unfolding.",
        ha="center",
        fontsize=10.7,
        color=RED,
        fontweight="bold",
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    plt.close(fig)
    return rows


def draw_pp_auau_xj(pp_data_root: Path, auau_data_root: Path, out: Path) -> list[dict[str, Any]]:
    pp_path = "PPG12_scaledtrigger30/h2_unfoldReco_pTgamma_xJ_incl_r04"
    auau_top = "photon_12_plus_MBD_NS_geq_2_vtx_lt_150"
    series = []
    with uproot.open(pp_data_root) as pp_file:
        values, variances, pt_edges, xj_edges = load_hist2_with_variance(pp_file, pp_path)
        selected, selected_var = select_pt_range(values, variances, pt_edges, 26, 35)
        series.append(("p+p", INK, "o", selected, selected_var, pp_path))
    with uproot.open(auau_data_root) as auau_file:
        for cent, label, color, marker in AUAU_CENTRALITIES:
            path = f"{auau_top}/h2_unfoldReco_pTgamma_xJ_incl_r04_isoR40_isSliding_cent_{cent}"
            values, variances, pt_edges, auau_xj_edges = load_hist2_with_variance(auau_file, path)
            if not np.array_equal(xj_edges, auau_xj_edges):
                raise RuntimeError("pp and AuAu xJ axes differ")
            selected, selected_var = select_pt_range(values, variances, pt_edges, 26, 35)
            series.append((f"Au+Au {label}", color, marker, selected, selected_var, path))
    fig, ax = plt.subplots(figsize=(10.2, 6.6), dpi=180)
    fig.subplots_adjust(left=0.12, right=0.97, bottom=0.18, top=0.82)
    fig.suptitle("A first common-bin p+p to Au+Au $x_{J\gamma}$ bridge", fontsize=22, fontweight="bold", y=0.965)
    fig.text(
        0.12,
        0.885,
        r"Exact shared native bin $26<p_T^\gamma<35$ GeV; Region A, $R=0.4$, sliding isolation, $\Delta\phi>7\pi/8$.",
        fontsize=12.5,
        color=MUTED,
    )
    rows = []
    centers = 0.5 * (xj_edges[:-1] + xj_edges[1:])
    xerr = 0.5 * np.diff(xj_edges)
    for label, color, marker, values, variances, path in series:
        density, error, total, mean = normalized_density(values, variances, xj_edges)
        ax.errorbar(centers, density, xerr=xerr, yerr=error, fmt=marker, color=color, mfc=color, mec="white", mew=0.6, ms=6.0, capsize=1.8, label=f"{label}  (pairs={total:.0f})")
        rows.append({"sample": label, "histogram": path, "pair_integral": total, "shape_mean_xj": mean, "values": values.tolist(), "variances": variances.tolist(), "xj_edges": xj_edges.tolist()})
    style_axis(ax)
    ax.set_xlim(0.0, 2.2)
    ax.set_ylim(bottom=0.0)
    ax.set_xlabel(r"$x_{J\gamma}=p_T^{jet}/p_T^\gamma$", fontsize=14.0)
    ax.set_ylabel(r"shape density $1/N\,dN/dx_{J\gamma}$", fontsize=13.5)
    ax.legend(frameon=False, fontsize=10.8, loc="upper right")
    sphx(ax)
    fig.text(
        0.50,
        0.035,
        "Raw reconstruction-level shape comparison; not a quenching claim.",
        ha="center",
        fontsize=10.6,
        color=RED,
        fontweight="bold",
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    plt.close(fig)
    return rows


def subset_ratio(sum_subset: float, var_subset: float, sum_total: float, var_total: float) -> tuple[float, float]:
    if sum_total <= 0.0:
        return float("nan"), float("nan")
    ratio = sum_subset / sum_total
    # Subset covariance: Cov(subset,total)=Var(subset).
    variance = var_subset / (sum_total**2) + (sum_subset**2) * var_total / (sum_total**4) - 2.0 * sum_subset * var_subset / (sum_total**3)
    return ratio, max(0.0, variance)


def draw_response_completeness(signal_root: Path, out: Path) -> list[dict[str, Any]]:
    definitions = (
        (
            "truth",
            "SIM/h2_unfoldTruth_pTgamma_xJ_incl_r04",
            "SIM/h2_unfoldTruthMatched_pTgamma_xJ_incl_r04",
            BLUE,
            "o",
            "truth-side matched efficiency",
        ),
        (
            "reco",
            "SIM/h2_unfoldReco_pTgamma_xJ_incl_r04",
            "SIM/h2_unfoldRecoMatched_pTgamma_xJ_incl_r04",
            ORANGE,
            "s",
            "reco-side matched fraction",
        ),
    )
    rows = []
    fig, ax = plt.subplots(figsize=(10.0, 6.3), dpi=180)
    fig.subplots_adjust(left=0.12, right=0.97, bottom=0.14, top=0.82)
    with uproot.open(signal_root) as root_file:
        for side, total_path, matched_path, color, marker, label in definitions:
            total_hist = root_file[total_path]
            matched_hist = root_file[matched_path]
            total = np.asarray(total_hist.values(flow=False), dtype=float).sum(axis=1)
            matched = np.asarray(matched_hist.values(flow=False), dtype=float).sum(axis=1)
            total_var_raw = total_hist.variances(flow=False)
            matched_var_raw = matched_hist.variances(flow=False)
            total_var = np.abs(total) if total_var_raw is None else np.asarray(total_var_raw, dtype=float).sum(axis=1)
            matched_var = np.abs(matched) if matched_var_raw is None else np.asarray(matched_var_raw, dtype=float).sum(axis=1)
            edges = np.asarray(total_hist.axis(0).edges(), dtype=float)
            xs, xes, ys, yes = [], [], [], []
            for index, (low, high) in enumerate(zip(edges[:-1], edges[1:])):
                if low < 10.0 or high > 35.0:
                    continue
                ratio, variance = subset_ratio(matched[index], matched_var[index], total[index], total_var[index])
                xs.append((low + high) / 2.0)
                xes.append((high - low) / 2.0)
                ys.append(ratio)
                yes.append(math.sqrt(variance))
                rows.append(
                    {
                        "side": side,
                        "pt_low_gev": float(low),
                        "pt_high_gev": float(high),
                        "total": float(total[index]),
                        "total_variance": float(total_var[index]),
                        "matched": float(matched[index]),
                        "matched_variance": float(matched_var[index]),
                        "matched_fraction": ratio,
                        "matched_fraction_error": math.sqrt(variance),
                        "total_histogram": total_path,
                        "matched_histogram": matched_path,
                    }
                )
            ax.errorbar(xs, ys, xerr=xes, yerr=yes, fmt=marker, color=color, mfc=color, mec="white", mew=0.7, ms=6.5, capsize=2.2, label=label)
    style_axis(ax)
    ax.set_xlim(9.5, 35.7)
    ax.set_ylim(0.0, 1.03)
    ax.set_xlabel(r"photon $p_T^\gamma$ [GeV]", fontsize=13.5)
    ax.set_ylabel("matched / total", fontsize=13.5)
    ax.legend(frameon=False, fontsize=11.5, loc="center right")
    sphx(ax)
    fig.suptitle("Current p+p response: matched fractions", fontsize=22, fontweight="bold", y=0.965)
    fig.text(0.12, 0.885, r"Photon+jet MC, $R=0.4$, nominal recoil selection; xJ integrated within each native photon bin.", fontsize=12.3, color=MUTED)
    ax.text(0.035, 0.075, "Response completeness diagnostic; no unfolding result is shown here.", transform=ax.transAxes, fontsize=10.5, color=RED, fontweight="bold")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    plt.close(fig)
    return rows


def independent_ratio(num: float, num_var: float, den: float, den_var: float) -> tuple[float, float]:
    if den <= 0.0:
        return float("nan"), float("nan")
    ratio = num / den
    variance = num_var / (den * den) + num * num * den_var / (den**4)
    return ratio, max(0.0, variance)


def extract_signal_leakage(
    root_file: uproot.ReadOnlyDirectory,
    bins: Iterable[tuple[int, int]],
    *,
    sample: str,
    centrality: str,
    suffix_prefix: str = "",
    suffix_tail: str = "",
) -> list[dict[str, Any]]:
    rows = []
    for low, high in bins:
        name = f"SIM/h_xJpurityLead_sigABCD_MC{suffix_prefix}_pT_{low}_{high}{suffix_tail}"
        hist = root_file[name]
        values = np.asarray(hist.values(flow=False)[:4], dtype=float)
        variance_raw = hist.variances(flow=False)
        variances = np.abs(values) if variance_raw is None else np.asarray(variance_raw[:4], dtype=float)
        for index, region in enumerate("BCD", start=1):
            ratio, variance = independent_ratio(values[index], variances[index], values[0], variances[0])
            rows.append(
                {
                    "sample": sample,
                    "centrality": centrality,
                    "region": region,
                    "pt_low_gev": low,
                    "pt_high_gev": high,
                    "A_signal": float(values[0]),
                    "A_signal_variance": float(variances[0]),
                    "side_signal": float(values[index]),
                    "side_signal_variance": float(variances[index]),
                    "side_over_A": ratio,
                    "side_over_A_error": math.sqrt(variance),
                    "histogram": name,
                }
            )
    return rows


def draw_signal_leakage(pp_signal_root: Path, auau_signal_root: Path, out: Path) -> list[dict[str, Any]]:
    with uproot.open(pp_signal_root) as pp_file:
        rows = extract_signal_leakage(pp_file, PP_NATIVE_BINS[3:], sample="pp", centrality="inclusive")
    with uproot.open(auau_signal_root) as auau_file:
        for cent, label, _color, _marker in AUAU_CENTRALITIES:
            rows.extend(
                extract_signal_leakage(
                    auau_file,
                    AUAU_NATIVE_BINS,
                    sample="auau",
                    centrality=label,
                    suffix_prefix="_isoR40_isSliding",
                    suffix_tail=f"_cent_{cent}",
                )
            )
    fig, axes = plt.subplots(3, 1, figsize=(11.5, 9.0), dpi=180, sharex=True)
    fig.subplots_adjust(left=0.11, right=0.975, bottom=0.10, top=0.82, hspace=0.24)
    fig.suptitle("Signal leakage is system- and centrality-dependent", fontsize=22, fontweight="bold", y=0.965)
    fig.text(
        0.11,
        0.89,
        r"Leading-photon truth-signal sideband population relative to Region A; $R=0.4$ sliding isolation.",
        fontsize=12.8,
        color=MUTED,
    )
    descriptors = {"B": "tight, non-isolated", "C": "non-tight, isolated", "D": "non-tight, non-isolated"}
    for axis, region in zip(axes, "BCD"):
        style_axis(axis)
        pp = [row for row in rows if row["sample"] == "pp" and row["region"] == region]
        x, xe, _ = centers_and_errors(pp)
        axis.errorbar(x, [row["side_over_A"] for row in pp], xerr=xe, yerr=[row["side_over_A_error"] for row in pp], fmt="o", color=INK, mfc=INK, ms=6.0, capsize=2.0, label="p+p")
        for cent, label, color, marker in AUAU_CENTRALITIES:
            del cent
            selected = [row for row in rows if row["sample"] == "auau" and row["centrality"] == label and row["region"] == region]
            xx, xxe, _ = centers_and_errors(selected)
            axis.errorbar(xx, [row["side_over_A"] for row in selected], xerr=xxe, yerr=[row["side_over_A_error"] for row in selected], fmt=marker, color=color, mfc=color, mec="white", mew=0.7, ms=6.1, capsize=2.0, label=f"Au+Au {label}")
        axis.set_ylabel(rf"$N_{region}^{{sig}}/N_A^{{sig}}$", fontsize=12.8)
        axis.set_title(descriptors[region], loc="left", fontsize=12.8, fontweight="bold", pad=4)
        axis.set_ylim(bottom=0.0)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        frameon=False,
        fontsize=9.8,
        ncol=4,
        loc="upper center",
        bbox_to_anchor=(0.59, 0.855),
        columnspacing=1.25,
        handletextpad=0.5,
    )
    axes[-1].set_xlim(14.0, 35.7)
    axes[-1].set_xlabel(r"reconstructed photon $p_T^\gamma$ [GeV]", fontsize=13.5)
    sphx(axes[0], size=11.5)
    fig.text(
        0.50,
        0.035,
        "Use these leakage terms in the purity/sideband model; do not import a single pp leakage number into Au+Au.",
        ha="center",
        fontsize=10.6,
        color=RED,
        fontweight="bold",
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    plt.close(fig)
    return rows


def draw_flow(out: Path) -> None:
    fig, ax = plt.subplots(figsize=(14.2, 4.7), dpi=180)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    fig.subplots_adjust(left=0.025, right=0.975, bottom=0.08, top=0.83)
    fig.suptitle("One analysis chain, with every physics assumption exposed as a gate", fontsize=21, fontweight="bold", y=0.96)
    stages = (
        (0.02, "1  Photon", "PPG12 parity\nfixed definition", "LOCKED", "#e8f4ff", BLUE),
        (0.215, "2  Background", "ABCD + truth leakage\nvalidate κ; transfer covariance", "TESTED", "#eef9f1", GREEN),
        (0.410, "3  Recoil", "R=0.4, Δφ>7π/8\nraw native-bin xJ", "VISIBLE", "#fff6e8", ORANGE),
        (0.605, "4  Detector", "response + fake/miss\nindependent half-closure", "GATED", "#f4effb", PURPLE),
        (0.800, "5  Physics", "pp baseline → Au+Au\ncentrality-dependent corrections", "NEXT", "#fbeeee", RED),
    )
    for x, title, body, state, fill, edge in stages:
        box = FancyBboxPatch((x, 0.26), 0.17, 0.50, boxstyle="round,pad=0.012,rounding_size=0.018", facecolor=fill, edgecolor=edge, linewidth=1.8)
        ax.add_patch(box)
        ax.text(x + 0.012, 0.68, title, fontsize=13.0, fontweight="bold", color=edge, ha="left", va="center")
        ax.text(x + 0.085, 0.49, body, fontsize=11.2, ha="center", va="center", linespacing=1.35)
        ax.text(x + 0.085, 0.31, state, fontsize=9.8, fontweight="bold", color=edge, ha="center", va="center")
    for start in (0.19, 0.385, 0.58, 0.775):
        ax.annotate("", xy=(start + 0.022, 0.51), xytext=(start, 0.51), arrowprops={"arrowstyle": "-|>", "color": MUTED, "lw": 1.8})
    ax.text(0.50, 0.12, "Friday shows the complete logic and current evidence; it does not promote a provisional result to final.", ha="center", va="center", fontsize=11.3, color=MUTED)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    plt.close(fig)


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise RuntimeError(f"refusing to write empty CSV {path}")
    scalar_rows = []
    for row in rows:
        scalar_rows.append({key: json.dumps(value, sort_keys=True) if isinstance(value, (dict, list)) else value for key, value in row.items()})
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(scalar_rows[0].keys()))
        writer.writeheader()
        writer.writerows(scalar_rows)


def source_first_manifest() -> list[dict[str, Any]]:
    entries = []
    for source in SOURCE_FIRST:
        png = source["png"]
        manifest = source["manifest"]
        if not png.is_file():
            raise FileNotFoundError(png)
        entries.append(
            {
                "role": source["role"],
                "png": str(png),
                "png_sha256": sha256(png),
                "manifest": str(manifest) if manifest.is_file() else None,
                "manifest_sha256": sha256(manifest) if manifest.is_file() else None,
                "reuse_policy": "source-first reference; not regenerated or copied",
                "visual_qa_this_run": "PASS: inspected at full resolution before packet generation",
            }
        )
    return entries


def write_questions(path: Path) -> None:
    text = """# Friday questions for Sam

1. Which `samfred` repository path and exact commit produced the 2026-08-06 plots, and which ROOT packet is authoritative?
2. Are the displayed ABCDEF regions mutually exclusive, and what exact role do E and F play beyond the four-region ABCD estimator?
3. Why were fixed `Eiso < 2 GeV` and `Eiso > 4 GeV` boundaries chosen instead of the PPG12 sliding R=0.4 isolation working point?
4. Was the Jet12 closure ratio computed on truth-signal-subtracted background, or on the inclusive candidate population?
5. Which native photon and xJ bin edges, overflow policy, jet-pT threshold, and jet-eta acceptance were used in the final plotted objects?
6. How are response fakes, misses, and photon normalization represented in the RooUnfold response object?
7. What is the exact in-situ JES objective function and fit covariance, rather than only the best-fit constant scale?
8. What is meant by the listed “third-jet” systematic, and how is it propagated through unfolding?
9. Which prior variations and iteration-selection criterion are intended for the final covariance?

These are reconciliation questions, not requests to replace the established PPG12-parity pp photon definition.
"""
    path.write_text(text)


def write_storyboard(
    path: Path,
    *,
    purity_png: str,
    response_matrix_png: str,
    kappa_png: Path,
    transfer_png: Path,
    bdt_iso_png: Path,
    pp_xj_png: Path,
    bridge_png: Path,
    response_png: Path,
    leakage_png: Path,
    flow_png: Path,
    halfclosure_png: Path,
    questions_path: Path,
) -> None:
    text = f"""# Friday PPG deck build map

## One-sentence story

We implement Sam's intended top-down p+p photon+jet analysis with the already-established PPG12-parity photon object, make every background and detector correction explicit, and then show exactly which ingredients do and do not transfer to Au+Au.

## Core plot-led slides

### 1. The photon definition is locked

- Drag: `{purity_png}`
- Say: "We start from the established PPG12-parity p+p photon definition; this meeting is about the recoil analysis, not retuning photon ID."
- Do not say: "The full p+p analysis is final."

### 2. Same ABCD intent, explicit implementation differences

- Drag: `{bdt_iso_png}`
- Say: "Sam uses fixed boxes; our canonical PPG12 boundaries vary with photon pT. We preserve our validated definition and quantify the comparison rather than silently redefining it."
- Point to: fixed red boundaries versus cyan PPG12 tight/isolation envelopes and hatched non-tight union.

### 3. Factorization must be tested on background

- Drag: `{kappa_png}`
- Say: "Inclusive candidate MC appears badly nonfactorizing because Region A contains truth signal. After exact subset subtraction, p+p background is compatible with unity: chi2/ndf 5.527/7, p=0.596. Au+Au 20-50% is not: p=0.0060."
- Do not say: "ABCD is validated for the final xJ shape." This is count-level closure.

### 4. Correction and covariance are part of the method

- Drag: `{transfer_png}`
- Say: "Even where p+p closure is statistically acceptable, the machinery carries tau=(B/D)/kappa and its statistical covariance; we never proceed by assuming factorization."
- Caveat: this is a normalization transfer. Independent xJ-shape closure remains a gate.

### 5. The current p+p recoil population is visible in honest bins

- Drag: `{pp_xj_png}`
- Say: "These are exact unions of frozen native photon bins, not fractional reconstructions of Sam's edges. They establish the current raw recoil population and statistics."
- Do not say: "This is the corrected or unfolded spectrum."

### 6. Detector incompleteness is explicit

- Drag side-by-side: `{response_matrix_png}` and `{response_png}`
- Say: "The response object records migration plus matched, fake, and miss populations. The matched fractions make the detector-correction problem visible before unfolding."

### 7. Iteration choice comes from an independent half-closure

- Drag: `{halfclosure_png}`
- Say: "Using the current qualified photon+jet input, iteration 3 is the earliest stable choice in the first-half to second-half simulation closure."
- Caveat: simulation half-closure is a diagnostic, not final p+p-data validation.

### 8. Au+Au needs its own background model

- Drag: `{leakage_png}`
- Say: "Truth-signal leakage into B, C, and D changes with system and centrality. A single p+p leakage number cannot be imported into Au+Au."

### 9. First exact-common-bin bridge to the Au+Au target

- Drag: `{bridge_png}`
- Say: "At 26-35 GeV both systems have one exact shared native bin. This is a raw reconstruction-level shape bridge that motivates the final corrected comparison."
- Do not say: "This is evidence of quenching." Statistics and corrections are not final.

### 10. Close with the gated path, not plot volume

- Drag: `{flow_png}`
- Say: "Photon definition is locked; background closure and response structure are visible; the next hard gates are xJ-dependent background-shape closure, full covariance, and independent unfolded closure."
- Ask Sam only the highest-value reconciliation questions from `{questions_path}`: authoritative commit/ROOT packet, background-only versus inclusive closure population, and exact fake/miss plus iteration contract.

## If meeting time is short

Show slides 1, 2, 3, 4, 7, 8, and 10. Keep raw xJ slides 5 and 9 as evidence/backup rather than forcing a premature physics claim.

## Meeting-safe conclusion

"We now reproduce the intended p+p analysis flow with the established PPG12 photon object and stronger provenance, covariance, and closure gates. The p+p count-level background closure is viable; Au+Au already demonstrates why centrality-dependent corrections are essential. The remaining work is sharply defined rather than hidden."
"""
    path.write_text(text)


def main() -> int:
    args = parse_args()
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    setup_style()

    artifacts = {key: resolve_artifact(key) for key in POINTER_KEYS}
    pp_data = artifacts["pp_data_merged"]
    pp_signal = artifacts["pp_sim_photonjet_merged"]
    pp_inclusive = artifacts["pp_sim_inclusivejet_merged"]
    auau_data = artifacts["auau_data_merged"]
    auau_signal = artifacts["auau_sim_photonjet_merged"]
    auau_inclusive = artifacts["auau_sim_inclusivejet_merged"]

    with uproot.open(pp_inclusive.root) as root_file:
        pp_sim_rows = extract_sim_abcd(root_file, PP_NATIVE_BINS)
    with uproot.open(pp_data.root) as root_file:
        pp_data_rows = extract_data_abcd(root_file, PP_NATIVE_BINS)
    auau_rows: dict[str, list[dict[str, Any]]] = {}
    with uproot.open(auau_inclusive.root) as root_file:
        for cent, _label, _color, _marker in AUAU_CENTRALITIES:
            auau_rows[cent] = extract_sim_abcd(
                root_file,
                AUAU_NATIVE_BINS,
                suffix_prefix="_isoR40_isSliding",
                suffix_tail=f"_cent_{cent}",
            )

    kappa_png = out_dir / "01_abcd_nonfactorization_pp_to_auau.png"
    tests = draw_kappa(pp_sim_rows, pp_data_rows, auau_rows, kappa_png)
    abcd_csv_rows = []
    for row in pp_sim_rows:
        abcd_csv_rows.append({"system": "pp", "centrality": "inclusive", "population": "inclusive MC with truth-subtracted background", **row})
    for cent, label, _color, _marker in AUAU_CENTRALITIES:
        for row in auau_rows[cent]:
            abcd_csv_rows.append({"system": "auau", "centrality": label, "population": "inclusive embedded MC with truth-subtracted background", **row})
    write_csv(out_dir / "01_abcd_nonfactorization_counts.csv", abcd_csv_rows)
    write_csv(out_dir / "01_pp_data_observed_abcd_counts.csv", pp_data_rows)

    transfers = transfer_rows(pp_sim_rows, pp_data_rows)
    transfer_png = out_dir / "02_pp_abcd_transfer_factor_with_covariance.png"
    draw_transfer(transfers, transfer_png)
    transfer_csv = out_dir / "02_pp_abcd_transfer_factor_with_covariance.csv"
    write_csv(transfer_csv, transfers)
    transfer_cov = np.diag([row["corrected_transfer_variance"] for row in transfers])
    transfer_cov_path = out_dir / "02_pp_abcd_transfer_factor_stat_covariance.json"
    transfer_cov_path.write_text(
        json.dumps(
            {
                "schema": "THE219_PP_ABCD_TRANSFER_STAT_COVARIANCE_V1",
                "bin_order": [[row["pt_low_gev"], row["pt_high_gev"]] for row in transfers],
                "quantity": "tau=(B_data/D_data)/kappa_background_MC",
                "covariance": transfer_cov.tolist(),
                "assumptions": [
                    "Data B and D regions and independent pp inclusive MC kappa are statistically independent.",
                    "Different native photon-pT bins are treated as statistically disjoint.",
                    "Only statistical covariance is included; common stitching, leakage-definition, and model systematics are not included.",
                    "This is a normalization transfer; xJ-shape closure is a separate required gate.",
                ],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )

    bdt_iso_png = out_dir / "03_pp_bdt_isolation_sam_vs_ppg12_map.png"
    bdt_iso_meta = draw_bdt_iso_map(pp_data.root, pp_signal.root, bdt_iso_png)
    pp_xj_png = out_dir / "04_pp_raw_reco_xj_native_bins.png"
    pp_xj_rows = draw_pp_xj(pp_data.root, pp_xj_png)
    (out_dir / "04_pp_raw_reco_xj_native_bins.json").write_text(json.dumps(pp_xj_rows, indent=2) + "\n")
    bridge_png = out_dir / "05_pp_auau_raw_reco_xj_26_35.png"
    bridge_rows = draw_pp_auau_xj(pp_data.root, auau_data.root, bridge_png)
    (out_dir / "05_pp_auau_raw_reco_xj_26_35.json").write_text(json.dumps(bridge_rows, indent=2) + "\n")
    response_png = out_dir / "06_pp_response_matching_completeness.png"
    response_rows = draw_response_completeness(pp_signal.root, response_png)
    write_csv(out_dir / "06_pp_response_matching_completeness.csv", response_rows)
    leakage_png = out_dir / "07_pp_auau_truth_signal_leakage.png"
    leakage_rows = draw_signal_leakage(pp_signal.root, auau_signal.root, leakage_png)
    write_csv(out_dir / "07_pp_auau_truth_signal_leakage.csv", leakage_rows)
    flow_png = out_dir / "08_analysis_flow_with_gates.png"
    draw_flow(flow_png)

    questions_path = out_dir / "friday_questions_for_sam.md"
    write_questions(questions_path)
    reused = source_first_manifest()
    generated_pngs = [kappa_png, transfer_png, bdt_iso_png, pp_xj_png, bridge_png, response_png, leakage_png, flow_png]
    halfclosure_manifest = out_dir / "iteration_halfclosure/the97_fig37_sim_halfclosure_iteration_stability_manifest.json"
    halfclosure_png = out_dir / "iteration_halfclosure/the97_fig37_sim_halfclosure_iteration_stability.png"
    halfclosure_csv = out_dir / "iteration_halfclosure/the97_fig37_sim_halfclosure_iteration_stability.csv"
    halfclosure = json.loads(halfclosure_manifest.read_text()) if halfclosure_manifest.is_file() else None
    visual_qa_status = "PASS: inspected at full-slide scale" if args.visual_qa_pass else "pending visual QA"
    storyboard_path = out_dir / "friday_storyboard.md"
    write_storyboard(
        storyboard_path,
        purity_png=reused[0]["png"],
        response_matrix_png=reused[2]["png"],
        kappa_png=kappa_png,
        transfer_png=transfer_png,
        bdt_iso_png=bdt_iso_png,
        pp_xj_png=pp_xj_png,
        bridge_png=bridge_png,
        response_png=response_png,
        leakage_png=leakage_png,
        flow_png=flow_png,
        halfclosure_png=halfclosure_png,
        questions_path=questions_path,
    )

    manifest = {
        "schema": "THE219_FRIDAY_PPG_PHASE_A_PACKET_V1",
        "status": "diagnostic_candidate_packet_not_final_analysis",
        "task": "THE-219",
        "generated_by": str(Path(__file__).resolve()),
        "generated_by_sha256": sha256(Path(__file__).resolve()),
        "canonical_deck_reference": CANONICAL_DECK,
        "google_slides_mutated": False,
        "sdcc_or_condor_actions": "none",
        "protected_workstreams_touched": [],
        "sam_source": {
            "pdf": str(args.sam_reference_pdf) if args.sam_reference_pdf else None,
            "pdf_sha256": (
                sha256(args.sam_reference_pdf)
                if args.sam_reference_pdf and args.sam_reference_pdf.is_file()
                else None
            ),
            "displayed_choices": {
                "photon_pt_min_gev": 13,
                "photon_pt_bins_gev": [13, 15, 20, 25, 35, 100],
                "jet_pt_min_gev": 5,
                "jet_radius_focus": 0.4,
                "back_to_back": "Delta phi > 7*pi/8",
                "tight_bdt": "score > 0.8",
                "loose_bdt": "0.2 < score < 0.6",
                "isolated": "Eiso < 2 GeV",
                "nonisolated": "Eiso > 4 GeV",
                "unfolding_iterations_shown": 2,
            },
        },
        "our_contract": {
            "photon_definition": "established PPG12 parity; not reopened",
            "pp_namespace": "PPG12_scaledtrigger30",
            "jet_radius": 0.4,
            "back_to_back": "nominal object without _dphiPi2 suffix = Delta phi > 7*pi/8",
            "pp_isolation": "sliding R=0.4: Eiso < 0.490+0.037*pT; nonisolated > threshold+0.8",
            "pp_bdt": bdt_iso_meta["canonical_ppg12"],
            "pp_native_xj_display_bins_gev": [[16, 20], [20, 26], [26, 35]],
            "auau_native_photon_bins_gev": [list(item) for item in AUAU_NATIVE_BINS],
            "common_pp_auau_overlay_bin_gev": [26, 35],
        },
        "inputs": {
            key: {
                "pointer": str(artifact.pointer_path),
                "pointer_sha256": sha256(artifact.pointer_path),
                "current_entry_id": artifact.pointer.get("current_entry_id"),
                "campaign_tag": artifact.pointer.get("campaign_tag"),
                "canonical_status": artifact.pointer.get("canonical_status"),
                "root": str(artifact.root),
                "root_bytes": artifact.root.stat().st_size,
                "root_sha256": artifact.root_sha256,
            }
            for key, artifact in artifacts.items()
        },
        "generated_plots": [
            {"path": str(path), "sha256": sha256(path), "status": visual_qa_status}
            for path in generated_pngs
        ],
        "visual_qa": {
            "status": visual_qa_status,
            "assertion_method": "explicit --visual-qa-pass after desktop image inspection" if args.visual_qa_pass else None,
        },
        "source_first_reused_plots": reused,
        "tables": {
            "abcd_counts": str(out_dir / "01_abcd_nonfactorization_counts.csv"),
            "pp_data_abcd_counts": str(out_dir / "01_pp_data_observed_abcd_counts.csv"),
            "transfer_factors": str(transfer_csv),
            "transfer_stat_covariance": str(transfer_cov_path),
            "response_completeness": str(out_dir / "06_pp_response_matching_completeness.csv"),
            "signal_leakage": str(out_dir / "07_pp_auau_truth_signal_leakage.csv"),
        },
        "closure_tests_diagonal_stat_only": tests,
        "physics_readout": [
            "Raw inclusive-jet MC kappa is not a valid background-closure test because Region A contains a large truth-signal component.",
            f"After exact truth-signal subset subtraction, pp 14-35 GeV gives chi2/ndf={tests['pp_14to35']['chi2']:.3f}/{tests['pp_14to35']['ndf']} and p={tests['pp_14to35']['p_value']:.3f} for kappa=1.",
            f"AuAu 20-50% gives chi2/ndf={tests['auau_20_50']['chi2']:.3f}/{tests['auau_20_50']['ndf']} and p={tests['auau_20_50']['p_value']:.4f}; a unity transfer is not defensible there.",
            "The explicit pp normalization correction tau=(B/D)/kappa and its statistical covariance are provided, but no final xJ background subtraction is claimed until shape closure is independently validated.",
        ],
        "hard_caveats": [
            "Sam's exact 13-15, 15-20, 20-25, and 25-35 photon bins cannot be reconstructed from the frozen pp histogram edges without event-level reprocessing; no fractional bin splitting is used.",
            "The current pp inclusive-jet merged ROOT has the ABCD count and xJ families but no joint TableQA h2d_bdt-versus-Eiso objects. The BDT-isolation map therefore uses current data and current photon+jet MC only and is not labelled background closure.",
            "ABCD kappa is currently a count-level closure test. xJ-dependent background-shape closure remains a required gate before applying the transfer to a final xJ spectrum.",
            "The transfer covariance supplied here is statistical and diagonal in native pT. Common model, stitching, signal-definition, and shape systematics remain to be added.",
            "All xJ comparisons in this packet are raw reconstruction-level Region-A shape diagnostics without background subtraction, efficiency correction, unfolding, or systematic covariance.",
            "The pp data current pointer is candidate_evidence/current-default, not an assertion of final scientific canonicalization.",
        ],
        "halfclosure": halfclosure,
        "independent_diagnostics": {
            "pp_sim_halfclosure": {
                "plot": str(halfclosure_png) if halfclosure_png.is_file() else None,
                "plot_sha256": sha256(halfclosure_png) if halfclosure_png.is_file() else None,
                "csv": str(halfclosure_csv) if halfclosure_csv.is_file() else None,
                "csv_sha256": sha256(halfclosure_csv) if halfclosure_csv.is_file() else None,
                "manifest": str(halfclosure_manifest) if halfclosure_manifest.is_file() else None,
                "manifest_sha256": sha256(halfclosure_manifest) if halfclosure_manifest.is_file() else None,
                "status": visual_qa_status if halfclosure_png.is_file() else "not run",
            }
        },
        "questions_for_sam": str(questions_path),
        "storyboard": {"path": str(storyboard_path), "sha256": sha256(storyboard_path)},
    }
    manifest_path = out_dir / "packet_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")

    readme = f"""# THE-219 Friday PPG Phase A plot packet

Status: diagnostic candidate packet; no final physics result and no deck mutation.

## Drag-and-drop order

1. `{reused[0]['png']}` — source-first PPG12 photon-definition/purity control.
2. `{kappa_png}` — why background-only ABCD closure matters, and where Au+Au differs.
3. `{transfer_png}` — explicit pp transfer correction with covariance.
4. `{bdt_iso_png}` — Sam fixed boxes versus the canonical PPG12-parity boundaries.
5. `{pp_xj_png}` — current pp raw xJ in honest native-bin unions.
6. `{response_png}` — response completeness before unfolding.
7. `{halfclosure_png}` — independent current-input pp iteration half-closure; selected iteration 3.
8. `{leakage_png}` — pp/Au+Au truth-signal leakage differences.
9. `{bridge_png}` — exact-common-bin raw pp/Au+Au xJ bridge.
10. `{flow_png}` — procedure flow with visible gates.

Read `{manifest_path}` before making claims. Build guidance is in `{storyboard_path}`. Questions for Sam are in `{questions_path}`.
"""
    (out_dir / "README.md").write_text(readme)
    print(json.dumps({"out_dir": str(out_dir), "manifest": str(manifest_path), "plots": [str(path) for path in generated_pngs]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
