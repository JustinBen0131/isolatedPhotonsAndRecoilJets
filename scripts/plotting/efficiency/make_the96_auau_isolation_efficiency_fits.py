#!/usr/bin/env python3
"""Derive THE-96 Au+Au isolation-efficiency cutoffs and centrality fits.

The input is a merged embedded-photon RecoilJets ROOT file. For each photon-pT
and 5% centrality cell, the script reads the truth-matched reconstructed
isolation spectrum and finds the Eiso cut retaining 70%, 80%, and 90% of the
signal. It then fits a constant across pT in each centrality bin and a line
across the resulting centrality-dependent constants.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import uproot


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_INPUT = (
    REPO_ROOT
    / "InputFiles/the96_auau_iso_baseline5pct_20260709_1335"
    / "expanded_5to40_baseline_simembedded"
    / "RecoilJets_embeddedPhoton12plus20_MERGED.root"
)
DEFAULT_OUTDIR = (
    REPO_ROOT
    / "dataOutput/auau/isolation_baseline5pct_the96_20260709"
    / "r03_baseline_nominal_15to35"
)

CENTRALITY_BINS = [(lo, lo + 5) for lo in range(0, 80, 5)]
EXPANDED_PT_BINS = [
    (5, 8),
    (8, 10),
    (10, 12),
    (12, 14),
    (14, 15),
    (15, 17),
    (17, 19),
    (19, 21),
    (21, 23),
    (23, 26),
    (26, 30),
    (30, 35),
    (35, 40),
]
NOMINAL_PT_GROUPS = [
    ((15, 17), [(15, 17)]),
    ((17, 19), [(17, 19)]),
    ((19, 21), [(19, 21)]),
    ((21, 23), [(21, 23)]),
    ((23, 26), [(23, 26)]),
    ((26, 35), [(26, 30), (30, 35)]),
]
EFFICIENCIES = (0.70, 0.80, 0.90)
COLORS = {0.70: "#0072B2", 0.80: "#009E73", 0.90: "#D55E00"}
MARKERS = {0.70: "o", 0.80: "s", 0.90: "^"}


@dataclass(frozen=True)
class HistogramData:
    values: np.ndarray
    variances: np.ndarray
    edges: np.ndarray
    underflow: float
    overflow: float
    underflow_variance: float
    overflow_variance: float
    sources: tuple[str, ...]


def histogram_name(cone: str, pt_bin: tuple[int, int], cent_bin: tuple[int, int]) -> str:
    return (
        f"h_EisoReco_truthSigMatched_iso{cone}_"
        f"pT_{pt_bin[0]}_{pt_bin[1]}_cent_{cent_bin[0]}_{cent_bin[1]}"
    )


def read_histogram(directory, name: str) -> HistogramData:
    if name not in directory:
        raise KeyError(f"Missing ROOT histogram: SIM/{name}")
    hist = directory[name]
    values = np.asarray(hist.values(flow=False), dtype=float)
    variances = hist.variances(flow=False)
    if variances is None:
        variances = np.maximum(values, 0.0)
    variances = np.asarray(variances, dtype=float)
    edges = np.asarray(hist.axis().edges(), dtype=float)
    flow = np.asarray(hist.values(flow=True), dtype=float)
    flow_variances = hist.variances(flow=True)
    if flow_variances is None:
        flow_variances = np.maximum(flow, 0.0)
    flow_variances = np.asarray(flow_variances, dtype=float)
    return HistogramData(
        values=values,
        variances=variances,
        edges=edges,
        underflow=float(flow[0]),
        overflow=float(flow[-1]),
        underflow_variance=float(max(flow_variances[0], 0.0)),
        overflow_variance=float(max(flow_variances[-1], 0.0)),
        sources=(f"SIM/{name}",),
    )


def add_histograms(items: Iterable[HistogramData]) -> HistogramData:
    data = list(items)
    if not data:
        raise ValueError("Cannot combine an empty histogram list")
    edges = data[0].edges
    for item in data[1:]:
        if not np.array_equal(item.edges, edges):
            raise ValueError("Histogram bin edges differ within a requested pT group")
    return HistogramData(
        values=np.sum([item.values for item in data], axis=0),
        variances=np.sum([item.variances for item in data], axis=0),
        edges=edges.copy(),
        underflow=float(sum(item.underflow for item in data)),
        overflow=float(sum(item.overflow for item in data)),
        underflow_variance=float(sum(item.underflow_variance for item in data)),
        overflow_variance=float(sum(item.overflow_variance for item in data)),
        sources=tuple(source for item in data for source in item.sources),
    )


def quantile_from_histogram(data: HistogramData, probability: float) -> dict[str, float]:
    values = np.asarray(data.values, dtype=float)
    negative_tolerance = 1e-10 * max(1.0, float(np.max(np.abs(values))))
    if np.any(values < -negative_tolerance):
        raise ValueError("Isolation histogram has negative weighted bins; CDF quantile is undefined")
    values = np.maximum(values, 0.0)
    in_range_total = float(np.sum(values))
    underflow = max(float(data.underflow), 0.0)
    overflow = max(float(data.overflow), 0.0)
    total = in_range_total + underflow + overflow
    if not math.isfinite(total) or total <= 0.0:
        raise ValueError("Isolation histogram has no positive in-range weight")

    cumulative = np.cumsum(values)
    absolute_target = probability * total
    if absolute_target <= underflow:
        raise ValueError(
            f"Requested quantile {probability:.3f} lies in histogram underflow below {data.edges[0]:g} GeV"
        )
    if absolute_target > underflow + in_range_total:
        raise ValueError(
            f"Requested quantile {probability:.3f} lies in histogram overflow above {data.edges[-1]:g} GeV"
        )
    target = absolute_target - underflow
    index = int(np.searchsorted(cumulative, target, side="left"))
    index = min(max(index, 0), len(values) - 1)
    previous = float(cumulative[index - 1]) if index > 0 else 0.0
    bin_weight = float(values[index])
    fraction = 0.5 if bin_weight <= 0.0 else (target - previous) / bin_weight
    fraction = min(max(float(fraction), 0.0), 1.0)
    width = float(data.edges[index + 1] - data.edges[index])
    threshold = float(data.edges[index] + fraction * width)

    sumw2 = float(
        np.sum(np.maximum(data.variances, 0.0))
        + data.underflow_variance
        + data.overflow_variance
    )
    effective_entries = total * total / sumw2 if sumw2 > 0.0 else total
    density = bin_weight / (total * width) if width > 0.0 else 0.0
    if effective_entries > 0.0 and density > 0.0:
        probability_error = math.sqrt(probability * (1.0 - probability) / effective_entries)
        statistical_error = probability_error / density
    else:
        statistical_error = width
    threshold_error = max(0.5 * width, statistical_error)
    flow_fraction = (underflow + overflow) / max(total, 1e-12)
    return {
        "threshold_gev": threshold,
        "threshold_error_gev": float(threshold_error),
        "sumw": total,
        "sumw2": sumw2,
        "effective_entries": float(effective_entries),
        "underflow": data.underflow,
        "overflow": data.overflow,
        "flow_fraction": float(flow_fraction),
    }


def constant_fit(y: np.ndarray, yerr: np.ndarray) -> dict[str, float]:
    valid = np.isfinite(y) & np.isfinite(yerr) & (yerr > 0.0)
    if int(np.count_nonzero(valid)) < 2:
        raise ValueError("Constant fit needs at least two valid points")
    y = y[valid]
    yerr = yerr[valid]
    weights = 1.0 / np.square(yerr)
    value = float(np.sum(weights * y) / np.sum(weights))
    error = float(math.sqrt(1.0 / np.sum(weights)))
    chi2 = float(np.sum(np.square((y - value) / yerr)))
    return {
        "value_gev": value,
        "error_gev": error,
        "chi2": chi2,
        "ndf": int(y.size - 1),
        "n_points": int(y.size),
    }


def linear_fit(x: np.ndarray, y: np.ndarray, yerr: np.ndarray) -> dict:
    valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(yerr) & (yerr > 0.0)
    if int(np.count_nonzero(valid)) < 3:
        raise ValueError("Linear fit needs at least three valid points")
    x = x[valid]
    y = y[valid]
    yerr = yerr[valid]
    design = np.column_stack([np.ones_like(x), x])
    weights = 1.0 / np.square(yerr)
    normal = design.T @ (weights[:, None] * design)
    covariance = np.linalg.inv(normal)
    parameters = covariance @ (design.T @ (weights * y))
    prediction = design @ parameters
    chi2 = float(np.sum(np.square((y - prediction) / yerr)))
    return {
        "intercept_gev": float(parameters[0]),
        "slope_gev_per_percent": float(parameters[1]),
        "intercept_error_gev": float(math.sqrt(covariance[0, 0])),
        "slope_error_gev_per_percent": float(math.sqrt(covariance[1, 1])),
        "covariance": covariance.tolist(),
        "chi2": chi2,
        "ndf": int(y.size - 2),
        "n_points": int(y.size),
    }


def derive_points(root_path: Path, cone: str, pt_scheme: str) -> tuple[list[dict], list[dict], list[dict]]:
    if pt_scheme == "nominal":
        pt_groups = NOMINAL_PT_GROUPS
    else:
        pt_groups = [(pt_bin, [pt_bin]) for pt_bin in EXPANDED_PT_BINS]

    points: list[dict] = []
    with uproot.open(root_path) as root_file:
        if "SIM" not in root_file:
            raise KeyError(f"{root_path} does not contain a SIM directory")
        directory = root_file["SIM"]
        for cent_bin in CENTRALITY_BINS:
            for target_pt, source_pt_bins in pt_groups:
                direct_name = histogram_name(cone, target_pt, cent_bin)
                if target_pt == (26, 35) and direct_name in directory:
                    pieces = [read_histogram(directory, direct_name)]
                else:
                    pieces = [
                        read_histogram(directory, histogram_name(cone, source_pt, cent_bin))
                        for source_pt in source_pt_bins
                    ]
                combined = add_histograms(pieces)
                for efficiency in EFFICIENCIES:
                    stats = quantile_from_histogram(combined, efficiency)
                    points.append(
                        {
                            "efficiency": efficiency,
                            "cent_min": cent_bin[0],
                            "cent_max": cent_bin[1],
                            "cent_center": 0.5 * (cent_bin[0] + cent_bin[1]),
                            "pt_min": target_pt[0],
                            "pt_max": target_pt[1],
                            "pt_center": 0.5 * (target_pt[0] + target_pt[1]),
                            "source_histograms": " | ".join(combined.sources),
                            **stats,
                        }
                    )

    flat_fits: list[dict] = []
    for cent_bin in CENTRALITY_BINS:
        for efficiency in EFFICIENCIES:
            selected = [
                row
                for row in points
                if row["cent_min"] == cent_bin[0]
                and row["cent_max"] == cent_bin[1]
                and row["efficiency"] == efficiency
            ]
            fit = constant_fit(
                np.asarray([row["threshold_gev"] for row in selected], dtype=float),
                np.asarray([row["threshold_error_gev"] for row in selected], dtype=float),
            )
            flat_fits.append(
                {
                    "efficiency": efficiency,
                    "cent_min": cent_bin[0],
                    "cent_max": cent_bin[1],
                    "cent_center": 0.5 * (cent_bin[0] + cent_bin[1]),
                    **fit,
                }
            )

    centrality_fits: list[dict] = []
    for efficiency in EFFICIENCIES:
        selected = [row for row in flat_fits if row["efficiency"] == efficiency]
        fit = linear_fit(
            np.asarray([row["cent_center"] for row in selected], dtype=float),
            np.asarray([row["value_gev"] for row in selected], dtype=float),
            np.asarray([row["error_gev"] for row in selected], dtype=float),
        )
        centrality_fits.append({"efficiency": efficiency, **fit})
    return points, flat_fits, centrality_fits


def apply_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.linewidth": 1.1,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "legend.frameon": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def plot_centrality_grid(
    points: list[dict],
    flat_fits: list[dict],
    output: Path,
    cone_label: str,
    algorithm_label: str,
) -> None:
    apply_style()
    fig, axes = plt.subplots(4, 4, figsize=(17.0, 12.4), sharex=True, sharey=True)
    axes_flat = axes.ravel()
    thresholds = np.asarray([row["threshold_gev"] for row in points], dtype=float)
    errors = np.asarray([row["threshold_error_gev"] for row in points], dtype=float)
    y_min = float(np.nanmin(thresholds - errors))
    y_max = float(np.nanmax(thresholds + errors))
    padding = max(0.45, 0.08 * (y_max - y_min))
    x_line = np.asarray([15.0, 35.0])

    for axis, cent_bin in zip(axes_flat, CENTRALITY_BINS):
        for efficiency in EFFICIENCIES:
            selected = [
                row
                for row in points
                if row["cent_min"] == cent_bin[0]
                and row["cent_max"] == cent_bin[1]
                and row["efficiency"] == efficiency
            ]
            selected.sort(key=lambda row: row["pt_center"])
            fit = next(
                row
                for row in flat_fits
                if row["cent_min"] == cent_bin[0]
                and row["cent_max"] == cent_bin[1]
                and row["efficiency"] == efficiency
            )
            axis.errorbar(
                [row["pt_center"] for row in selected],
                [row["threshold_gev"] for row in selected],
                yerr=[row["threshold_error_gev"] for row in selected],
                fmt=MARKERS[efficiency],
                color=COLORS[efficiency],
                markersize=4.4,
                markeredgewidth=0.7,
                capsize=1.8,
                linewidth=1.0,
                label=f"{int(100 * efficiency)}%",
                zorder=3,
            )
            axis.plot(
                x_line,
                np.full_like(x_line, fit["value_gev"]),
                color=COLORS[efficiency],
                linewidth=1.35,
                linestyle="--",
                zorder=2,
            )
        axis.set_title(f"{cent_bin[0]}-{cent_bin[1]}%", fontsize=11.5, fontweight="bold")
        axis.set_xlim(14.5, 35.5)
        axis.set_ylim(y_min - padding, y_max + padding)
        axis.grid(True, color="#d9d9d9", linewidth=0.65, alpha=0.7)
        axis.minorticks_on()
        axis.tick_params(labelsize=9)

    for axis in axes[-1, :]:
        axis.set_xlabel(r"$p_T^\gamma$ [GeV]", fontsize=11)
    for axis in axes[:, 0]:
        axis.set_ylabel(r"$E_T^{\mathrm{iso,cut}}$ [GeV]", fontsize=11)

    handles, labels = axes_flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, fontsize=11, bbox_to_anchor=(0.5, 0.946))
    fig.suptitle(
        "Au+Au embedded-photon isolation efficiency thresholds",
        fontsize=17,
        fontweight="bold",
        y=0.992,
    )
    fig.text(
        0.5,
        0.962,
        f"{cone_label} {algorithm_label} | nominal 15-35 GeV | constant fits within each 5% centrality bin",
        ha="center",
        va="top",
        fontsize=11.5,
    )
    fig.text(0.052, 0.958, r"$\it{\bf{sPHENIX}}$ Internal", fontsize=12.2, ha="left", va="top")
    fig.tight_layout(rect=[0.035, 0.035, 0.995, 0.925], h_pad=1.15, w_pad=0.75)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=190)
    plt.close(fig)


def plot_centrality_summary(
    flat_fits: list[dict],
    centrality_fits: list[dict],
    output: Path,
    cone_label: str,
    algorithm_label: str,
) -> None:
    apply_style()
    fig, axis = plt.subplots(figsize=(11.2, 7.6))
    x_fit = np.linspace(0.0, 80.0, 240)
    equation_lines: list[str] = []
    all_y: list[float] = []
    all_err: list[float] = []

    for efficiency in EFFICIENCIES:
        selected = [row for row in flat_fits if row["efficiency"] == efficiency]
        selected.sort(key=lambda row: row["cent_center"])
        fit = next(row for row in centrality_fits if row["efficiency"] == efficiency)
        x = np.asarray([row["cent_center"] for row in selected], dtype=float)
        y = np.asarray([row["value_gev"] for row in selected], dtype=float)
        yerr = np.asarray([row["error_gev"] for row in selected], dtype=float)
        all_y.extend(y.tolist())
        all_err.extend(yerr.tolist())
        axis.errorbar(
            x,
            y,
            yerr=yerr,
            fmt=MARKERS[efficiency],
            color=COLORS[efficiency],
            markersize=6.5,
            capsize=2.4,
            linewidth=1.1,
            label=f"{int(100 * efficiency)}% efficiency",
            zorder=3,
        )
        prediction = fit["intercept_gev"] + fit["slope_gev_per_percent"] * x_fit
        covariance = np.asarray(fit["covariance"], dtype=float)
        design = np.column_stack([np.ones_like(x_fit), x_fit])
        prediction_error = np.sqrt(np.einsum("ij,jk,ik->i", design, covariance, design))
        axis.plot(x_fit, prediction, color=COLORS[efficiency], linewidth=2.0, zorder=2)
        axis.fill_between(
            x_fit,
            prediction - prediction_error,
            prediction + prediction_error,
            color=COLORS[efficiency],
            alpha=0.12,
            linewidth=0,
            zorder=1,
        )
        equation_lines.append(
            f"{int(100 * efficiency)}%: $E_T^{{iso,cut}}$ = {fit['intercept_gev']:.3f} "
            f"{fit['slope_gev_per_percent']:+.4f}$C$  "
            f"($\\chi^2$/ndf={fit['chi2']:.1f}/{fit['ndf']})"
        )

    y = np.asarray(all_y, dtype=float)
    yerr = np.asarray(all_err, dtype=float)
    padding = max(0.45, 0.10 * float(np.max(y + yerr) - np.min(y - yerr)))
    axis.set_xlim(0.0, 80.0)
    axis.set_ylim(float(np.min(y - yerr) - padding), float(np.max(y + yerr) + padding))
    axis.set_xlabel("Centrality [%]", fontsize=13)
    axis.set_ylabel(r"Flat-fit $E_T^{\mathrm{iso,cut}}$ [GeV]", fontsize=13)
    axis.grid(True, color="#d7d7d7", linewidth=0.75, alpha=0.75)
    axis.minorticks_on()
    axis.tick_params(labelsize=11)
    axis.legend(loc="upper right", fontsize=11)
    axis.text(
        0.035,
        0.055,
        "\n".join(equation_lines),
        transform=axis.transAxes,
        fontsize=10.5,
        va="bottom",
        ha="left",
        bbox={"facecolor": "white", "edgecolor": "#bdbdbd", "alpha": 0.92, "boxstyle": "square,pad=0.45"},
    )
    fig.suptitle(
        "Centrality-dependent Au+Au isolation working points",
        fontsize=17,
        fontweight="bold",
        y=0.985,
    )
    axis.set_title(
        f"{cone_label} {algorithm_label} | nominal 15-35 GeV flat fits",
        fontsize=12,
        pad=12,
    )
    fig.text(0.085, 0.952, r"$\it{\bf{sPHENIX}}$ Internal", fontsize=12.5, ha="left", va="top")
    fig.tight_layout(rect=[0.035, 0.035, 0.995, 0.935])
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=190)
    plt.close(fig)


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"No rows available for {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--cone", choices=("R30", "R40"), default="R30")
    parser.add_argument("--algorithm", choices=("baseline", "topocluster"), default="baseline")
    parser.add_argument("--pt-scheme", choices=("nominal", "expanded"), default="nominal")
    args = parser.parse_args()

    input_path = args.input.resolve()
    output_dir = args.outdir.resolve()
    if not input_path.is_file():
        raise SystemExit(f"Input ROOT does not exist: {input_path}")

    points, flat_fits, centrality_fits = derive_points(input_path, args.cone, args.pt_scheme)
    max_flow_fraction = max(row["flow_fraction"] for row in points)
    cone_label = "R=0.3" if args.cone == "R30" else "R=0.4"
    algorithm_label = "tower-cone baseline (no topocluster)" if args.algorithm == "baseline" else "topocluster isolation"
    stem = f"the96_auau_iso_{args.cone.lower()}_{args.algorithm}_{args.pt_scheme}"
    grid_png = output_dir / f"{stem}_5pct_flat_fit_grid.png"
    summary_png = output_dir / f"{stem}_centrality_linear_fits.png"
    points_csv = output_dir / f"{stem}_cutoff_points.csv"
    flat_csv = output_dir / f"{stem}_flat_fits.csv"
    linear_csv = output_dir / f"{stem}_linear_fits.csv"
    manifest_json = output_dir / f"{stem}_manifest.json"

    plot_centrality_grid(points, flat_fits, grid_png, cone_label, algorithm_label)
    plot_centrality_summary(flat_fits, centrality_fits, summary_png, cone_label, algorithm_label)
    write_csv(points_csv, points)
    write_csv(flat_csv, flat_fits)
    write_csv(linear_csv, centrality_fits)
    manifest = {
        "campaign": "THE-96",
        "campaign_tag": "the96_auau_iso_baseline5pct_20260709_1335",
        "input_root": str(input_path),
        "root_directory": "SIM",
        "histogram_family": f"h_EisoReco_truthSigMatched_iso{args.cone}_pT_*_cent_*",
        "algorithm": args.algorithm,
        "cone": cone_label,
        "pt_scheme": args.pt_scheme,
        "centrality_bins": [list(item) for item in CENTRALITY_BINS],
        "efficiencies": list(EFFICIENCIES),
        "cutoff_point_count": len(points),
        "flat_fit_count": len(flat_fits),
        "linear_fit_count": len(centrality_fits),
        "max_flow_fraction": max_flow_fraction,
        "outputs": {
            "grid_png": str(grid_png),
            "summary_png": str(summary_png),
            "cutoff_points_csv": str(points_csv),
            "flat_fits_csv": str(flat_csv),
            "linear_fits_csv": str(linear_csv),
        },
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_json.write_text(json.dumps(manifest, indent=2) + "\n")

    print(f"input={input_path}")
    print(f"cutoff_points={len(points)} flat_fits={len(flat_fits)} linear_fits={len(centrality_fits)}")
    print(f"max_flow_fraction={max_flow_fraction:.6g}")
    print(f"grid_png={grid_png}")
    print(f"summary_png={summary_png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
