#!/usr/bin/env python3
"""Build pp-versus-AuAu Region-A yield and photon-purity slide candidates.

The pp points use the current full-stat THE97 PPG12-parity data and
photon+jet artifacts.  The AuAu points use the completed THE100 unrestricted
BDT-complement data and photon embedding, always from the stored sliding-R=0.4
ABCD views.  Both systems use the PPG12 toy estimator for raw and
signal-leakage-corrected photon purity.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import ROOT
import uproot
from matplotlib.lines import Line2D


THIS_FILE = Path(__file__).resolve()
REPO = next((parent for parent in THIS_FILE.parents if (parent / "AGENTS.md").exists()), THIS_FILE.parents[3])
PP_HELPERS = REPO / "scripts/plotting/pp_currentian"
if str(PP_HELPERS) not in sys.path:
    sys.path.insert(0, str(PP_HELPERS))

from make_ppg12_fig3_purity_sim_sdcc_vs_current_overlay import ppg12_toy_estimate  # noqa: E402


ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

OUTDIR = REPO / "dataOutput/slides/the45_jstg_20260720/the100_pp_regionA_purity_comparison"

THE100_DATA = REPO / (
    "InputFiles/the100_auau_dualview_20260714/complement/data/"
    "RecoilJets_auau_ALL_preselectionNewPPG12_tightAuAuCentInputBase3x3BDT_"
    "nonTightAuAuBDTComplement_baseVariant.root"
)
THE100_SIGNAL = REPO / (
    "InputFiles/the100_auau_dualview_20260714/complement/signal/"
    "RecoilJets_embeddedPhoton12plus20_MERGED.root"
)

PP_BINS = ((16, 18), (18, 20), (20, 22), (22, 24), (24, 26), (26, 35))
AUAU_BINS = ((15, 17), (17, 19), (19, 21), (21, 23), (23, 26), (26, 35))
# The paired yield/purity slides use the two bookend centralities; the raw
# purity overlay also includes the completed middle-centrality view.
CENTRALITIES = (("0_20", "0--20%"), ("50_80", "50--80%"))
OVERLAY_CENTRALITIES = (("0_20", "0--20%"), ("20_50", "20--50%"), ("50_80", "50--80%"))

PP_DATA_TOP = "PPG12_scaledtrigger30"
PP_SIM_TOP = "SIM"
AUAU_DATA_TOP = "photon_12_plus_MBD_NS_geq_2_vtx_lt_150"
AUAU_SIM_TOP = "SIM"

INK = "#142235"
MUTED = "#52657A"
GRID = "#D7E0EA"
RAW = "#111827"
CORRECTED = "#2563EB"
PP_RED = "#D94B43"
AUAU_CENTRALITY_COLORS = {
    "0_20": "#C53B32",
    "20_50": "#2563EB",
    "50_80": "#15803D",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def resolve_current(sample: str) -> tuple[Path, Path, dict[str, object]]:
    pointer = REPO / "dataOutput/current_recoiljets_artifacts/current" / sample / "current.json"
    payload = json.loads(pointer.read_text())
    roots = payload.get("root_paths")
    if not isinstance(roots, list) or len(roots) != 1:
        raise RuntimeError(f"{pointer} must resolve exactly one ROOT")
    root = Path(str(roots[0]))
    if not root.is_file():
        raise FileNotFoundError(root)
    return root, pointer, payload


def hist_stats(handle: uproot.ReadOnlyDirectory, key: str) -> tuple[float, float]:
    hist = handle[key]
    values = np.asarray(hist.values(flow=False), dtype=float)
    variances = hist.variances(flow=False)
    if variances is None:
        variances = np.clip(values, 0.0, None)
    return float(np.sum(values)), float(np.sum(np.asarray(variances, dtype=float)))


def effective_count(value: float, variance: float, label: str) -> float:
    if value <= 0.0 or variance <= 0.0:
        raise RuntimeError(f"nonpositive effective count for {label}: value={value}, variance={variance}")
    return value * value / variance


def leakage_fraction(num: float, num_var: float, den: float, den_var: float, label: str) -> tuple[float, float]:
    if den <= 0.0:
        raise RuntimeError(f"nonpositive signal-region denominator for {label}")
    value = num / den
    variance = num_var / (den * den) + (num * num * den_var) / (den**4)
    return value, math.sqrt(max(variance, 0.0))


def region_key(system: str, top: str, region: str, lo: int, hi: int, cent: str | None) -> str:
    if system == "pp":
        return f"{top}/h_Eiso_ABCD_{region}_pT_{lo}_{hi}"
    assert cent is not None
    return f"{top}/h_Eiso_ABCD_{region}_isoR40_isSliding_pT_{lo}_{hi}_cent_{cent}"


def build_system_points(
    *,
    system: str,
    data_root: Path,
    signal_root: Path,
    bins: tuple[tuple[int, int], ...],
    data_top: str,
    sim_top: str,
    cent: str | None,
    rng: ROOT.TRandom3,
) -> list[dict[str, float | str | bool]]:
    points: list[dict[str, float | str | bool]] = []
    with uproot.open(data_root) as data, uproot.open(signal_root) as signal:
        for index, (lo, hi) in enumerate(bins):
            data_stats: dict[str, tuple[float, float]] = {}
            signal_stats: dict[str, tuple[float, float]] = {}
            data_keys: dict[str, str] = {}
            signal_keys: dict[str, str] = {}
            for region in "ABCD":
                data_key = region_key(system, data_top, region, lo, hi, cent)
                signal_key = region_key(system, sim_top, region, lo, hi, cent)
                data_keys[region] = data_key
                signal_keys[region] = signal_key
                data_stats[region] = hist_stats(data, data_key)
                signal_stats[region] = hist_stats(signal, signal_key)

            values = tuple(data_stats[region][0] for region in "ABCD")
            counts = tuple(
                effective_count(data_stats[region][0], data_stats[region][1], f"{system} {cent} {lo}-{hi} {region}")
                for region in "ABCD"
            )
            signal_a, signal_a_var = signal_stats["A"]
            leakage_with_error = tuple(
                leakage_fraction(
                    signal_stats[region][0],
                    signal_stats[region][1],
                    signal_a,
                    signal_a_var,
                    f"{system} {cent} {lo}-{hi} {region}/A",
                )
                for region in "BCD"
            )
            leakage = tuple(item[0] for item in leakage_with_error)
            leakage_errors = tuple(item[1] for item in leakage_with_error)
            raw, raw_error, corrected, corrected_error, diagnostics = ppg12_toy_estimate(
                rng,
                values,
                counts,
                leakage,
                leakage_errors,
                f"the45_{system}_{cent or 'all'}_{lo}_{hi}_{index}",
            )
            width = float(hi - lo)
            a_value, a_variance = data_stats["A"]
            points.append(
                {
                    "system": system,
                    "centrality": cent or "all",
                    "pt_lo": float(lo),
                    "pt_hi": float(hi),
                    "pt_center": 0.5 * (lo + hi),
                    "pt_half_width": 0.5 * width,
                    "A": a_value,
                    "A_error": math.sqrt(max(a_variance, 0.0)),
                    "A_per_GeV": a_value / width,
                    "A_per_GeV_error": math.sqrt(max(a_variance, 0.0)) / width,
                    "B": values[1],
                    "C": values[2],
                    "D": values[3],
                    "raw_purity": raw,
                    "raw_purity_error": raw_error,
                    "corrected_purity": corrected,
                    "corrected_purity_error": corrected_error,
                    "fB": leakage[0],
                    "fC": leakage[1],
                    "fD": leakage[2],
                    "fB_error": leakage_errors[0],
                    "fC_error": leakage_errors[1],
                    "fD_error": leakage_errors[2],
                    "data_A_key": data_keys["A"],
                    "data_B_key": data_keys["B"],
                    "data_C_key": data_keys["C"],
                    "data_D_key": data_keys["D"],
                    "signal_A_key": signal_keys["A"],
                    "signal_B_key": signal_keys["B"],
                    "signal_C_key": signal_keys["C"],
                    "signal_D_key": signal_keys["D"],
                    **diagnostics,
                }
            )
    return points


def arrays(points: list[dict[str, float | str | bool]], key: str) -> np.ndarray:
    return np.asarray([float(point[key]) for point in points], dtype=float)


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
            "axes.linewidth": 1.20,
            "axes.edgecolor": INK,
            "axes.labelcolor": INK,
            "xtick.color": INK,
            "ytick.color": INK,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )


def errorbar(ax: plt.Axes, points: list[dict[str, float | str | bool]], y: str, yerr: str, **kwargs: object) -> None:
    ax.errorbar(
        arrays(points, "pt_center"),
        arrays(points, y),
        xerr=arrays(points, "pt_half_width"),
        yerr=arrays(points, yerr),
        linestyle="none",
        capsize=3.0,
        elinewidth=1.25,
        **kwargs,
    )


def draw_slide(
    cent_key: str,
    cent_label: str,
    pp_points: list[dict[str, float | str | bool]],
    auau_points: list[dict[str, float | str | bool]],
) -> Path:
    setup_style()
    output = OUTDIR / f"the100_pp_auau_{cent_key}_regionA_yield_and_purity_slide.png"
    fig = plt.figure(figsize=(16, 9), dpi=160, facecolor="white")
    gs = fig.add_gridspec(
        2,
        1,
        left=0.090,
        right=0.972,
        # Lowered from 0.805 now that the caption moved up under the title, and
        # extended to the bottom margin the caption used to occupy.
        top=0.780,
        bottom=0.088,
        hspace=0.28,
        height_ratios=(0.95, 1.05),
    )
    ax_yield = fig.add_subplot(gs[0])
    ax_purity = fig.add_subplot(gs[1], sharex=ax_yield)

    fig.text(
        0.052,
        0.945,
        f"{cent_label.replace('--', '–')} Au+Au and p+p: Region-A yield and photon purity",
        ha="left",
        va="top",
        fontsize=34.5,
        fontweight="bold",
        color=INK,
    )
    # Shared JSTG subtitle-arrowhead convention: DejaVu Sans "▶" in #2468A8
    # with a 0.021-figure-width hanging indent.
    fig.text(0.054, 0.862, "▶", ha="left", va="center", fontsize=15.0, color="#2468A8", fontfamily="DejaVu Sans")
    # Deck convention: bold "Label:" then regular body.  The body is placed one
    # space after the label's measured right edge so the two read as one line.
    subtitle_label = fig.text(
        0.075, 0.862, "Top panel:", fontsize=20.0, color=INK, ha="left", va="center", fontweight="bold"
    )
    fig.canvas.draw()
    label_right = subtitle_label.get_window_extent(renderer=fig.canvas.get_renderer()).x1
    body_x = fig.transFigure.inverted().transform((label_right, 0))[0] + 0.006
    fig.text(
        body_x,
        0.862,
        r"counts divided by $E_T$ bin width only; compare the purity curves, not the yield offset.",
        fontsize=20.0,
        color=INK,
        ha="left",
        va="center",
    )

    errorbar(
        ax_yield,
        auau_points,
        "A_per_GeV",
        "A_per_GeV_error",
        marker="o",
        markersize=8.2,
        markerfacecolor=RAW,
        markeredgecolor="white",
        markeredgewidth=0.8,
        color=RAW,
        ecolor="#46505F",
        label=f"Au+Au {cent_label.replace('--', '–')}",
    )
    errorbar(
        ax_yield,
        pp_points,
        "A_per_GeV",
        "A_per_GeV_error",
        marker="s",
        markersize=8.6,
        markerfacecolor="white",
        markeredgecolor=PP_RED,
        markeredgewidth=2.0,
        color=PP_RED,
        ecolor="#EE817B",
        label="p+p PPG12 selection",
    )
    ax_yield.set_yscale("log")
    positive = np.concatenate((arrays(auau_points, "A_per_GeV"), arrays(pp_points, "A_per_GeV")))
    ax_yield.set_ylim(max(float(np.min(positive[positive > 0])) * 0.45, 1.0e-2), float(np.max(positive)) * 2.4)
    ax_yield.set_ylabel(r"Region-A candidates / GeV", fontsize=17.0)
    ax_yield.set_title("Raw tight-and-isolated Region-A candidates", fontsize=19.0, fontweight="bold", pad=9, color=INK)
    ax_yield.grid(True, which="major", color=GRID, linewidth=0.9, alpha=0.90)
    ax_yield.grid(True, which="minor", axis="y", color="#E8EEF5", linewidth=0.55, alpha=0.65)
    ax_yield.legend(loc="upper right", frameon=False, fontsize=20.0, handlelength=1.5, handletextpad=0.6, labelspacing=0.55)
    ax_yield.text(0.018, 0.625, r"$\it{\bf{sPHENIX}}$ Internal", transform=ax_yield.transAxes, ha="left", va="top", fontsize=14.0)
    ax_yield.text(0.018, 0.515, r"$\sqrt{s_{NN}}=200$ GeV", transform=ax_yield.transAxes, ha="left", va="top", fontsize=12.7, color=MUTED)
    plt.setp(ax_yield.get_xticklabels(), visible=False)

    errorbar(
        ax_purity,
        auau_points,
        "raw_purity",
        "raw_purity_error",
        marker="o",
        markersize=7.6,
        markerfacecolor=RAW,
        markeredgecolor="white",
        markeredgewidth=0.7,
        color=RAW,
        ecolor="#46505F",
        label="Au+Au raw",
    )
    errorbar(
        ax_purity,
        auau_points,
        "corrected_purity",
        "corrected_purity_error",
        marker="o",
        markersize=8.2,
        markerfacecolor=CORRECTED,
        markeredgecolor="white",
        markeredgewidth=0.7,
        color=CORRECTED,
        ecolor="#5D8EF0",
        label="Au+Au leakage-corrected",
    )
    errorbar(
        ax_purity,
        pp_points,
        "raw_purity",
        "raw_purity_error",
        marker="s",
        markersize=7.6,
        markerfacecolor="white",
        markeredgecolor=RAW,
        markeredgewidth=1.7,
        color=RAW,
        ecolor="#677080",
        label="p+p raw",
    )
    errorbar(
        ax_purity,
        pp_points,
        "corrected_purity",
        "corrected_purity_error",
        marker="s",
        markersize=8.2,
        markerfacecolor="white",
        markeredgecolor=CORRECTED,
        markeredgewidth=1.9,
        color=CORRECTED,
        ecolor="#7DA4F3",
        label="p+p leakage-corrected",
    )
    all_low: list[float] = []
    all_high: list[float] = []
    for point_set in (pp_points, auau_points):
        for value_key, error_key in (("raw_purity", "raw_purity_error"), ("corrected_purity", "corrected_purity_error")):
            values = arrays(point_set, value_key)
            errors = arrays(point_set, error_key)
            all_low.extend((values - errors).tolist())
            all_high.extend((values + errors).tolist())
    purity_low = min(-0.04, min(all_low) - 0.04)
    # Keep clear headroom above the purity points so the two-column legend
    # remains visually separate from the plotted measurements.
    purity_high = max(1.10, max(all_high) + 0.05)
    # Headroom to 1.60.  The highest drawn purity point is 0.9402 including its
    # error bar (pp corrected, 24-26 GeV); at 1.42 the enlarged two-row legend
    # came within a few pixels of it, so the axis is opened further.
    ax_purity.set_ylim(max(-0.35, purity_low), 1.60)
    ax_purity.set_xlim(14.2, 35.8)
    ax_purity.set_ylabel("Photon purity", fontsize=17.0)
    ax_purity.set_xlabel(r"reconstructed photon $E_T$ [GeV]", fontsize=17.0, labelpad=8)
    ax_purity.set_title("Raw and signal-leakage-corrected ABCD purity", fontsize=19.0, fontweight="bold", pad=9, color=INK)
    ax_purity.grid(True, which="major", color=GRID, linewidth=0.9, alpha=0.90)
    ax_purity.axhline(1.0, color="#A8B3C1", linestyle=(0, (4, 4)), linewidth=1.0, zorder=0)
    ax_purity.legend(loc="upper left", bbox_to_anchor=(0.015, 0.985), ncol=2, frameon=False, fontsize=20.0, handlelength=1.5, columnspacing=2.4, handletextpad=0.6, labelspacing=0.55)
    for ax in (ax_yield, ax_purity):
        ax.tick_params(which="both", labelsize=13.0, length=6)
        ax.tick_params(which="minor", length=3)
        ax.minorticks_on()

    fig.savefig(output, dpi=160, facecolor="white")
    plt.close(fig)
    return output


def draw_raw_purity_overlay(
    pp_points: list[dict[str, float | str | bool]],
    auau_points_by_centrality: dict[str, list[dict[str, float | str | bool]]],
) -> Path:
    """Draw one direct raw-ABCD comparison across the completed systems."""
    setup_style()
    output = OUTDIR / "the100_pp_auau_raw_purity_all_centralities_slide.png"
    fig = plt.figure(figsize=(16, 9), dpi=160, facecolor="white")
    ax = fig.add_axes((0.095, 0.205, 0.875, 0.605))

    fig.text(
        0.052,
        0.945,
        "Raw ABCD photon purity across p+p and Au+Au centrality",
        ha="left",
        va="top",
        fontsize=30.5,
        fontweight="bold",
        color=INK,
    )
    fig.text(
        0.052,
        0.858,
        "Completed production baseline: the same tight, sliding-isolation Region-A/B/C/D definition is evaluated in each system.",
        ha="left",
        va="top",
        fontsize=16.0,
        color=MUTED,
    )

    all_low: list[float] = []
    all_high: list[float] = []
    for cent_key, cent_label in OVERLAY_CENTRALITIES:
        points = auau_points_by_centrality[cent_key]
        color = AUAU_CENTRALITY_COLORS[cent_key]
        errorbar(
            ax,
            points,
            "raw_purity",
            "raw_purity_error",
            marker="o",
            markersize=8.2,
            markerfacecolor=color,
            markeredgecolor="white",
            markeredgewidth=0.8,
            color=color,
            ecolor=color,
            label=f"Au+Au {cent_label.replace('--', '–')}",
        )
        values = arrays(points, "raw_purity")
        errors = arrays(points, "raw_purity_error")
        all_low.extend((values - errors).tolist())
        all_high.extend((values + errors).tolist())

    errorbar(
        ax,
        pp_points,
        "raw_purity",
        "raw_purity_error",
        marker="s",
        markersize=8.4,
        markerfacecolor="white",
        markeredgecolor=RAW,
        markeredgewidth=2.0,
        color=RAW,
        ecolor="#677080",
        label="p+p PPG12 selection",
    )
    pp_values = arrays(pp_points, "raw_purity")
    pp_errors = arrays(pp_points, "raw_purity_error")
    all_low.extend((pp_values - pp_errors).tolist())
    all_high.extend((pp_values + pp_errors).tolist())

    purity_low = min(-0.04, min(all_low) - 0.04)
    purity_high = max(1.10, max(all_high) + 0.05)
    ax.set_ylim(max(-0.15, purity_low), min(1.45, purity_high))
    ax.set_xlim(14.2, 35.8)
    ax.set_xlabel(r"reconstructed photon $E_T$ [GeV]", fontsize=20.0, labelpad=10)
    ax.set_ylabel("Raw ABCD photon purity", fontsize=20.0)
    ax.grid(True, which="major", color=GRID, linewidth=0.9, alpha=0.90)
    ax.axhline(1.0, color="#A8B3C1", linestyle=(0, (4, 4)), linewidth=1.0, zorder=0)
    ax.legend(
        loc="upper left",
        bbox_to_anchor=(0.012, 0.985),
        ncol=2,
        frameon=False,
        fontsize=15.5,
        handlelength=1.25,
        columnspacing=1.9,
    )
    ax.text(0.985, 0.070, r"$\it{\bf{sPHENIX}}$ Internal", transform=ax.transAxes, ha="right", va="bottom", fontsize=16.0)
    ax.text(0.985, 0.017, r"$\sqrt{s_{NN}}=200$ GeV", transform=ax.transAxes, ha="right", va="bottom", fontsize=14.2, color=MUTED)
    ax.tick_params(which="both", labelsize=15.0, length=7)
    ax.tick_params(which="minor", length=3.5)
    ax.minorticks_on()

    fig.text(0.052, 0.095, r"$\blacktriangleright$", fontsize=17.0, color=INK, ha="left", va="center")
    fig.text(
        0.078,
        0.095,
        "Filled circles are Au+Au centrality intervals; open squares are the full-stat pp PPG12-parity reference."
        "  Points retain their native photon-energy bins.",
        fontsize=16.0,
        color=INK,
        ha="left",
        va="center",
    )
    fig.savefig(output, dpi=160, facecolor="white")
    plt.close(fig)
    return output


def write_csv(path: Path, *point_sets: list[dict[str, float | str | bool]]) -> None:
    rows = [point for point_set in point_sets for point in point_set]
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    pp_data, pp_data_pointer, pp_data_payload = resolve_current("pp_data_merged")
    pp_signal, pp_signal_pointer, pp_signal_payload = resolve_current("pp_sim_photonjet_merged")
    rng = ROOT.TRandom3(42)
    pp_points = build_system_points(
        system="pp",
        data_root=pp_data,
        signal_root=pp_signal,
        bins=PP_BINS,
        data_top=PP_DATA_TOP,
        sim_top=PP_SIM_TOP,
        cent=None,
        rng=rng,
    )

    outputs: dict[str, dict[str, str]] = {}
    all_auau_points: dict[str, list[dict[str, float | str | bool]]] = {}
    for cent_key, cent_label in CENTRALITIES:
        auau_points = build_system_points(
            system="auau",
            data_root=THE100_DATA,
            signal_root=THE100_SIGNAL,
            bins=AUAU_BINS,
            data_top=AUAU_DATA_TOP,
            sim_top=AUAU_SIM_TOP,
            cent=cent_key,
            rng=rng,
        )
        all_auau_points[cent_key] = auau_points
        csv_path = OUTDIR / f"the100_pp_auau_{cent_key}_regionA_yield_and_purity_points.csv"
        write_csv(csv_path, pp_points, auau_points)
        png = draw_slide(cent_key, cent_label, pp_points, auau_points)
        outputs[cent_key] = {
            "png": str(png),
            "png_sha256": sha256(png),
            "points_csv": str(csv_path),
            "points_csv_sha256": sha256(csv_path),
        }

    for cent_key, _ in OVERLAY_CENTRALITIES:
        if cent_key in all_auau_points:
            continue
        all_auau_points[cent_key] = build_system_points(
            system="auau",
            data_root=THE100_DATA,
            signal_root=THE100_SIGNAL,
            bins=AUAU_BINS,
            data_top=AUAU_DATA_TOP,
            sim_top=AUAU_SIM_TOP,
            cent=cent_key,
            rng=rng,
        )

    raw_overlay_csv = OUTDIR / "the100_pp_auau_raw_purity_all_centralities_points.csv"
    write_csv(raw_overlay_csv, pp_points, *(all_auau_points[key] for key, _ in OVERLAY_CENTRALITIES))
    raw_overlay_png = draw_raw_purity_overlay(pp_points, all_auau_points)
    outputs["raw_purity_all_centralities"] = {
        "png": str(raw_overlay_png),
        "png_sha256": sha256(raw_overlay_png),
        "points_csv": str(raw_overlay_csv),
        "points_csv_sha256": sha256(raw_overlay_csv),
    }

    manifest = {
        "schema": "THE45_THE100_PP_AUAU_REGIONA_PURITY_SLIDES_V1",
        "reference_slide": {
            "deck_id": "18qRtYLb3UHs_ClmgS4EGv2z9VODVHHhZs0YgeSFgjHQ",
            "slide_object_id": "g3efaeb2c01d_0_30",
            "matched_contract": "stacked Region-A yield and raw/leakage-corrected ABCD purity comparison",
        },
        "pp": {
            "data_pointer": str(pp_data_pointer),
            "data_pointer_payload": pp_data_payload,
            "data_root": str(pp_data),
            "data_sha256": sha256(pp_data),
            "signal_pointer": str(pp_signal_pointer),
            "signal_pointer_payload": pp_signal_payload,
            "signal_root": str(pp_signal),
            "signal_sha256": sha256(pp_signal),
            "trigger_namespace": PP_DATA_TOP,
            "pt_bins": list(PP_BINS),
        },
        "auau": {
            "campaign": "the100_auau_dualview_20260714",
            "status": "completed historical production baseline",
            "view": "unrestricted BDT complement",
            "data_root": str(THE100_DATA),
            "data_sha256": sha256(THE100_DATA),
            "signal_root": str(THE100_SIGNAL),
            "signal_sha256": sha256(THE100_SIGNAL),
            "data_namespace": AUAU_DATA_TOP,
            "isolation_view": "isoR40_isSliding",
            "pt_bins": list(AUAU_BINS),
            "centralities": [item[0] for item in OVERLAY_CENTRALITIES],
        },
        "estimator": {
            "yield": "Region-A integral divided by native ET-bin width",
            "raw_purity": "PPG12 toy estimator applied with zero leakage coefficients",
            "corrected_purity": "PPG12 quadratic physical-root toy estimator using photon-embedding B/A, C/A, and D/A leakage fractions",
            "rng": "ROOT.TRandom3(42), one sequential stream",
            "toys_per_point": 20000,
        },
        "outputs": outputs,
        "slides_mutated": False,
        "scientific_boundary": "THE100 is the completed pre-THE111 AuAu production baseline. The slides are valid historical-production comparisons and do not promote THE111 or replace the pending corrected triplet.",
    }
    manifest_path = OUTDIR / "the100_pp_regionA_purity_comparison_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")

    speaker = OUTDIR / "the100_pp_regionA_purity_comparison_speaker_script.md"
    speaker.write_text(
        "# Region-A yield and photon-purity comparison\n\n"
        "The upper panel compares the raw tight-and-isolated Region-A candidate yield in Au+Au with the current full-stat pp PPG12-parity reference. The points are divided by their native photon-energy bin widths, but they are not luminosity or exposure normalized, so the upper panel is a population and shape diagnostic rather than an absolute system-to-system yield comparison.\n\n"
        "The lower panel is the direct physics comparison. Black points are the raw ABCD purity and blue points include the prompt-photon leakage correction. Filled circles are Au+Au and open squares are pp. The two slides use the same pp reference; only the Au+Au centrality changes from 0-20 percent to 50-80 percent.\n\n"
        "The Au+Au input is the completed THE100 sliding-isolation production baseline. It predates the pending corrected THE111-scored production and is presented as the current completed-production comparison, not as a promotion decision.\n"
        "\nThe centrality-overlay slide isolates the raw ABCD estimator: it overlays the three completed Au+Au centrality intervals with the full-stat pp reference. The pp and Au+Au points retain their native photon-energy bins, so this is a common-estimator comparison rather than a bin-by-bin ratio.\n"
    )

    layout = {
        "canvas": [0, 0, 2560, 1440],
        "title_axis_x": 133,
        "title_axis_tolerance_px": 16,
        "minimum_audience_font_px": 26,
        "minimum_title_font_px": 60,
        "nodes": [
            {"kind": "text", "name": "slide title", "role": "title", "bbox": [133, 44, 2460, 137], "font_px": 68, "title_axis_align": "left"},
            {"kind": "plot", "name": "Region-A yield panel", "role": "plot", "bbox": [133, 260, 2488, 740], "title_axis_align": "left"},
            {"kind": "plot", "name": "photon purity panel", "role": "plot", "bbox": [133, 790, 2488, 1260], "title_axis_align": "left"},
            {"kind": "text", "name": "interpretation note", "role": "caption", "bbox": [133, 1300, 2460, 1380], "font_px": 34, "title_axis_align": "left"},
        ],
    }
    (OUTDIR / "the100_pp_regionA_purity_comparison_layout_nodes.json").write_text(json.dumps(layout, indent=2) + "\n")

    raw_overlay_layout = {
        "canvas": [0, 0, 2560, 1440],
        "title_axis_x": 133,
        "title_axis_tolerance_px": 16,
        "minimum_audience_font_px": 26,
        "minimum_title_font_px": 60,
        "nodes": [
            {"kind": "text", "name": "slide title", "role": "title", "bbox": [133, 44, 2460, 137], "font_px": 68, "title_axis_align": "left"},
            {"kind": "text", "name": "subtitle", "role": "caption", "bbox": [133, 170, 2450, 220], "font_px": 36, "title_axis_align": "left"},
            # This bounding box includes the y-axis title/ticks as well as the
            # data canvas, so the complete plot group aligns to the title axis.
            {"kind": "plot", "name": "raw-purity overlay", "role": "plot", "bbox": [133, 295, 2480, 1165], "title_axis_align": "left"},
            {"kind": "text", "name": "interpretation note", "role": "caption", "bbox": [133, 1275, 2460, 1360], "font_px": 36, "title_axis_align": "left"},
        ],
    }
    (OUTDIR / "the100_pp_auau_raw_purity_all_centralities_layout_nodes.json").write_text(
        json.dumps(raw_overlay_layout, indent=2) + "\n"
    )

    print(json.dumps({"outputs": outputs, "manifest": str(manifest_path), "speaker": str(speaker)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
