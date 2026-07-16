#!/usr/bin/env python3
"""Plot and audit embedded AuAu MBD/sEPD low-calorimeter diagnostics."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm, TwoSlopeNorm
import numpy as np
import uproot


HIST_PREFIX = "SIM/"
CENTRALITY_GROUPS = ((0.0, 20.0), (20.0, 50.0), (50.0, 80.0))
LOW_BAND_GROUPS = ((0.0, 20.0), (20.0, 40.0), (40.0, 55.0))
CLASS_INDEX = {"all": 0, "low": 1, "main": 2}
MB_DECISIONS = ("mbPass", "mbFail", "mbMissing")

# ROOT kBird-like palette used for audience-facing calorimeter-density maps.
KBIRD_LIKE = LinearSegmentedColormap.from_list(
    "root_kbird_like",
    [
        (0.00, "#00004f"),
        (0.18, "#0033a0"),
        (0.36, "#1178bd"),
        (0.55, "#2fb7b5"),
        (0.74, "#8bd3a8"),
        (0.90, "#e7ef8a"),
        (1.00, "#ffffcc"),
    ],
)


@dataclass(frozen=True)
class Hist:
    values: np.ndarray
    variances: np.ndarray
    edges: tuple[np.ndarray, ...]


@dataclass(frozen=True)
class Sample:
    key: str
    label: str
    path: Path
    root: uproot.ReadOnlyFile


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--signal-root", required=True, type=Path)
    parser.add_argument("--inclusive-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--dpi", type=int, default=180)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_hist(sample: Sample, name: str, *, required: bool = True) -> Hist | None:
    key = HIST_PREFIX + name
    if key not in sample.root:
        if required:
            raise KeyError(f"{sample.path}: missing {key}")
        return None
    obj = sample.root[key]
    values = np.asarray(obj.values(flow=False), dtype=float)
    raw_variances = obj.variances(flow=False)
    variances = np.zeros_like(values) if raw_variances is None else np.asarray(raw_variances, dtype=float)
    edges = tuple(np.asarray(axis.edges(flow=False), dtype=float) for axis in obj.axes)
    return Hist(values, variances, edges)


def zero_like(hist: Hist) -> Hist:
    return Hist(np.zeros_like(hist.values), np.zeros_like(hist.variances), hist.edges)


def read_minbias_splits(sample: Sample, base_name: str) -> dict[str, Hist]:
    reference = read_hist(sample, base_name)
    splits: dict[str, Hist] = {}
    for suffix in MB_DECISIONS:
        hist = read_hist(sample, f"{base_name}_{suffix}", required=suffix != "mbMissing")
        splits[suffix] = hist if hist is not None else zero_like(reference)
    return splits


def centers(edges: np.ndarray) -> np.ndarray:
    return 0.5 * (edges[:-1] + edges[1:])


def select_bins(edges: np.ndarray, low: float, high: float) -> np.ndarray:
    c = centers(edges)
    return np.flatnonzero((c >= low) & (c < high))


def positive_limits(values: np.ndarray) -> tuple[float, float] | None:
    positive = values[np.isfinite(values) & (values > 0.0)]
    if positive.size == 0:
        return None
    return float(np.min(positive)), float(np.max(positive))


def short_label(sample: Sample) -> str:
    return sample.label.removeprefix("Embedded ")


def draw_density(ax: plt.Axes, hist2: np.ndarray, xedges: np.ndarray, yedges: np.ndarray) -> None:
    limits = positive_limits(hist2)
    if limits is None:
        ax.text(0.5, 0.5, "No entries", ha="center", va="center", transform=ax.transAxes)
        return
    lo, hi = limits
    norm = LogNorm(vmin=max(lo, hi * 1.0e-5), vmax=hi)
    ax.pcolormesh(xedges, yedges, hist2.T, shading="auto", cmap="viridis", norm=norm)


def shared_log_norm(densities: Iterable[np.ndarray]) -> LogNorm:
    positive = np.concatenate(
        [values[np.isfinite(values) & (values > 0.0)] for values in densities]
    )
    if positive.size == 0:
        return LogNorm(vmin=1.0, vmax=10.0)
    high = float(np.max(positive))
    return LogNorm(vmin=max(float(np.min(positive)), high * 1.0e-6), vmax=high)


def style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 11,
            "axes.linewidth": 1.0,
            "axes.grid": False,
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


def save(fig: plt.Figure, path: Path, dpi: int) -> None:
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def low_band_threshold(centrality: np.ndarray) -> np.ndarray:
    return (
        -8.484848484848325e-05 * centrality * centrality
        - 0.006181818181818246 * centrality
        + 2.5830757575757572
    )


def plot_low_band(samples: Iterable[Sample], output: Path, dpi: int) -> None:
    samples = tuple(samples)
    fig, axes = plt.subplots(1, len(samples), figsize=(12.6, 5.1), sharex=True, sharey=True)
    for ax, sample in zip(np.atleast_1d(axes), samples):
        hist = read_hist(sample, "h3_pmtDiag_logTotalCaloVsMbdCharge")
        density = np.sum(hist.values, axis=0)
        draw_density(ax, density.T, hist.edges[2], hist.edges[1])
        c = np.linspace(0.0, 55.0, 400)
        ax.plot(c, low_band_threshold(c), color="#d62728", linewidth=2.0, label="diagnostic boundary")
        ax.set_title(sample.label, fontweight="bold")
        ax.set_xlabel("Centrality [%]")
        ax.legend(frameon=False, loc="lower left")
    axes[0].set_ylabel(r"$\log_{10}(\max(0,E_{calo})+1)$")
    fig.suptitle("Signed good-tower calorimeter sum versus centrality", fontsize=17, fontweight="bold")
    fig.text(0.5, 0.92, "Au+Au embedded simulation; boundary is a diagnostic tag, not an event cut", ha="center")
    fig.subplots_adjust(top=0.82, wspace=0.08)
    save(fig, output / "low_calo_band_signal_vs_inclusive.png", dpi)


def total_calo_vs_centrality(sample: Sample, suffix: str | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    name = "h3_pmtDiag_totalCaloEnergyVsMbdCharge"
    if suffix:
        name = f"{name}_{suffix}"
    hist = read_hist(sample, name)
    density_energy_cent = np.sum(hist.values, axis=0)
    return density_energy_cent.T, hist.edges[2], hist.edges[1]


def plot_zero_based_total_calo(samples: Iterable[Sample], output: Path, dpi: int) -> None:
    samples = tuple(samples)
    payload = [(sample, *total_calo_vs_centrality(sample)) for sample in samples]
    norm = shared_log_norm(row[1] for row in payload)
    fig, axes = plt.subplots(1, len(samples), figsize=(13.4, 5.6), sharex=True, sharey=True)
    image = None
    for ax, (sample, density, cent_edges, energy_edges) in zip(np.atleast_1d(axes), payload):
        image = ax.pcolormesh(cent_edges, energy_edges, density.T, shading="auto", cmap="viridis", norm=norm)
        ax.set_title(sample.label, fontsize=13, fontweight="bold")
        ax.set_xlabel("Centrality percentile [%]")
        ax.set_xlim(0.0, 80.0)
        ax.set_ylim(0.0, min(3000.0, float(energy_edges[-1])))
    axes[0].set_ylabel(r"$E_{calo}^{total}$ [GeV]")
    if image is not None:
        fig.colorbar(image, ax=np.atleast_1d(axes).tolist(), pad=0.02, label="Weighted event count / bin")
    fig.suptitle("Centrality versus total calorimeter energy", fontsize=18, fontweight="bold")
    fig.text(0.5, 0.92, "Au+Au embedded simulation; signed finite good-tower sum; vertical range begins at zero", ha="center")
    fig.subplots_adjust(top=0.83, wspace=0.08, right=0.88)
    save(fig, output / "total_calo_vs_centrality_zero_based_all.png", dpi)

    comparison_rows = (
        (None, "No MinimumBiasClassifier requirement"),
        ("mbPass", "MinimumBiasClassifier pass"),
    )
    comparison_payload = [
        total_calo_vs_centrality(sample, suffix)
        for suffix, _ in comparison_rows
        for sample in samples
    ]
    comparison_norm = shared_log_norm(density for density, _, _ in comparison_payload)
    fig, axes = plt.subplots(2, len(samples), figsize=(13.0, 10.0), sharex=True, sharey=True)
    image = None
    for row, (suffix, row_label) in enumerate(comparison_rows):
        for column, sample in enumerate(samples):
            density, cent_edges, energy_edges = total_calo_vs_centrality(sample, suffix)
            ax = axes[row, column]
            image = ax.pcolormesh(
                cent_edges,
                energy_edges,
                density.T,
                shading="auto",
                cmap=KBIRD_LIKE,
                norm=comparison_norm,
            )
            if row == 0:
                ax.set_title(short_label(sample), fontsize=14, fontweight="bold")
            ax.set_xlim(0.0, 80.0)
            ax.set_ylim(0.0, min(3000.0, float(energy_edges[-1])))
            if row == len(comparison_rows) - 1:
                ax.set_xlabel("Centrality percentile [%]")
            if column == 0:
                ax.set_ylabel(r"$E_{calo}^{total}$ [GeV]")
                ax.text(
                    -0.31,
                    0.5,
                    row_label,
                    transform=ax.transAxes,
                    rotation=90,
                    ha="center",
                    va="center",
                    fontsize=12,
                    fontweight="bold",
                )
    if image is not None:
        fig.colorbar(image, ax=axes.ravel().tolist(), pad=0.02, label="Weighted event count / bin")
    fig.suptitle("Centrality versus total calorimeter energy", fontsize=19, fontweight="bold", y=0.975)
    fig.text(0.16, 0.895, r"$\it{\bf{sPHENIX}}$ Internal", ha="left", fontsize=14)
    fig.text(
        0.88,
        0.895,
        "Au+Au embedded simulation; signed finite good-tower sum; common logarithmic color scale",
        ha="right",
        fontsize=10,
    )
    fig.subplots_adjust(left=0.16, top=0.84, bottom=0.10, hspace=0.10, wspace=0.08, right=0.80)
    save(fig, output / "total_calo_vs_centrality_minbias_2x2_kbird.png", dpi)

    decisions = ((None, "All evaluated"), ("mbPass", "MinimumBiasClassifier pass"), ("mbFail", "MinimumBiasClassifier fail"))
    split_payload = []
    for sample in samples:
        for suffix, label in decisions:
            split_payload.append((sample, label, *total_calo_vs_centrality(sample, suffix)))
    split_norm = shared_log_norm(row[2] for row in split_payload)
    fig, axes = plt.subplots(len(samples), len(decisions), figsize=(16.0, 9.2), sharex=True, sharey=True)
    image = None
    for ax, (sample, label, density, cent_edges, energy_edges) in zip(axes.flat, split_payload):
        image = ax.pcolormesh(cent_edges, energy_edges, density.T, shading="auto", cmap="viridis", norm=split_norm)
        ax.set_title(f"{short_label(sample)}\n{label}", fontsize=12, fontweight="bold")
        ax.set_xlim(0.0, 80.0)
        ax.set_ylim(0.0, min(3000.0, float(energy_edges[-1])))
        ax.set_xlabel("Centrality percentile [%]")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$E_{calo}^{total}$ [GeV]")
    if image is not None:
        fig.colorbar(image, ax=axes.ravel().tolist(), pad=0.015, label="Weighted event count / bin")
    fig.suptitle("Total calorimeter energy by minimum-bias classification", fontsize=18, fontweight="bold")
    fig.text(0.5, 0.94, "Au+Au embedded simulation; common axes and color scale; PMT diagnostics filled before classifier rejection", ha="center")
    fig.subplots_adjust(top=0.87, hspace=0.25, wspace=0.08, right=0.90)
    save(fig, output / "total_calo_vs_centrality_zero_based_minbias_split.png", dpi)


def plot_mbd_calo_correlations(samples: Iterable[Sample], output: Path, dpi: int) -> None:
    for sample in samples:
        hist = read_hist(sample, "h3_pmtDiag_logTotalCaloVsMbdCharge")
        fig, axes = plt.subplots(1, 3, figsize=(15.2, 4.8), sharex=True, sharey=True)
        for ax, (low, high) in zip(axes, CENTRALITY_GROUPS):
            idx = select_bins(hist.edges[2], low, high)
            density = np.sum(hist.values[:, :, idx], axis=2)
            draw_density(ax, density, hist.edges[0], hist.edges[1])
            ax.set_title(f"{low:.0f}-{high:.0f}%")
            ax.set_xlabel(r"$\Sigma Q_{MBD}$")
        axes[0].set_ylabel(r"$\log_{10}(\max(0,E_{calo})+1)$")
        fig.suptitle(f"MBD charge and calorimeter-energy correlation: {sample.label}", fontsize=16, fontweight="bold")
        fig.text(0.5, 0.91, "Au+Au embedded simulation; signed finite good-tower sums", ha="center")
        fig.subplots_adjust(top=0.80, wspace=0.08)
        save(fig, output / f"mbd_vs_total_calo_by_centrality_{sample.key}.png", dpi)


def plot_minbias_classifier_comparison(samples: Iterable[Sample], output: Path, dpi: int) -> None:
    colors = {"signal": "#1f77b4", "inclusive": "#d62728"}
    fig, ax = plt.subplots(figsize=(10.2, 6.2))
    for sample in samples:
        decisions = read_minbias_splits(sample, "h3_pmtDiag_mbdNFiredTotalByClass")
        passed = decisions["mbPass"]
        failed = decisions["mbFail"]
        missing = decisions["mbMissing"]
        pass_y = np.sum(passed.values[:, :, CLASS_INDEX["all"]], axis=0)
        fail_y = np.sum(failed.values[:, :, CLASS_INDEX["all"]], axis=0)
        missing_y = np.sum(missing.values[:, :, CLASS_INDEX["all"]], axis=0)
        total = pass_y + fail_y + missing_y
        fraction = np.divide(pass_y, total, out=np.full_like(pass_y, np.nan), where=total > 0.0)
        x = centers(passed.edges[1])
        ax.plot(x, fraction, marker="o", linewidth=1.8, color=colors[sample.key], label=sample.label)
    ax.set_xlabel("Centrality [%]")
    ax.set_ylabel("MinimumBiasClassifier pass fraction")
    ax.set_ylim(0.0, 1.05)
    ax.set_xlim(0.0, 80.0)
    ax.legend(frameon=False)
    fig.suptitle("Embedded minimum-bias classifier acceptance", fontsize=17, fontweight="bold")
    fig.text(0.5, 0.91, "Pass fraction before photon or jet selection; |z| < 30 cm", ha="center")
    fig.subplots_adjust(top=0.82)
    save(fig, output / "minimum_bias_classifier_pass_fraction_by_centrality.png", dpi)

    for sample in samples:
        fig, axes = plt.subplots(2, 3, figsize=(15.2, 8.9), sharex=True, sharey=True)
        for row, (suffix, label) in enumerate((("mbPass", "classifier pass"), ("mbFail", "classifier fail"))):
            hist = read_hist(sample, f"h3_pmtDiag_logTotalCaloVsMbdCharge_{suffix}")
            for col, cent_range in enumerate(CENTRALITY_GROUPS):
                idx = select_bins(hist.edges[2], *cent_range)
                density = np.sum(hist.values[:, :, idx], axis=2)
                draw_density(axes[row, col], density, hist.edges[0], hist.edges[1])
                axes[row, col].set_title(
                    f"{label}, {cent_range[0]:.0f}-{cent_range[1]:.0f}%", fontsize=11
                )
                axes[row, col].set_xlabel(r"$\Sigma Q_{MBD}$")
                if col == 0:
                    axes[row, col].set_ylabel(r"$\log_{10}(\max(0,E_{calo})+1)$")
        fig.suptitle(f"Minimum-bias classifier comparison: {sample.label}", fontsize=17, fontweight="bold")
        fig.text(0.5, 0.94, "PMT diagnostics are filled before classifier-fail events are rejected", ha="center")
        fig.subplots_adjust(top=0.87, hspace=0.25, wspace=0.10)
        save(fig, output / f"minimum_bias_classifier_pass_fail_calo_{sample.key}.png", dpi)


def plot_old_style_correlations(samples: Iterable[Sample], output: Path, dpi: int) -> None:
    definitions = (
        ("emcal", "h3_pmtDiag_emcalEnergyVsMbdCharge", "EMCal signed-good tower sum [GeV]"),
        ("ihcal", "h3_pmtDiag_ihcalEnergyVsMbdCharge", "IHCal signed-good tower sum [GeV]"),
        ("ohcal", "h3_pmtDiag_ohcalEnergyVsMbdCharge", "OHCal signed-good tower sum [GeV]"),
        ("total", "h3_pmtDiag_totalCaloEnergyVsMbdCharge", "Total signed-good tower sum [GeV]"),
    )
    for short, name, ylabel in definitions:
        fig, axes = plt.subplots(2, 3, figsize=(15.2, 9.2), sharex=True, sharey=True)
        for row, sample in enumerate(samples):
            hist = read_hist(sample, name)
            for col, (low, high) in enumerate(CENTRALITY_GROUPS):
                ax = axes[row, col]
                idx = select_bins(hist.edges[2], low, high)
                density = np.sum(hist.values[:, :, idx], axis=2)
                draw_density(ax, density, hist.edges[0], hist.edges[1])
                ax.set_title(f"{short_label(sample)}, {low:.0f}-{high:.0f}%", fontsize=11)
                ax.set_xlabel(r"$\Sigma Q_{MBD}$")
                if col == 0:
                    ax.set_ylabel(ylabel)
        fig.suptitle(f"MBD charge versus {short.upper()} energy", fontsize=17, fontweight="bold")
        fig.text(0.5, 0.94, "Au+Au embedded simulation; historical correlation remade with current signed good-tower sums", ha="center")
        fig.subplots_adjust(top=0.88, hspace=0.24, wspace=0.10)
        save(fig, output / f"old_style_mbd_vs_{short}_by_centrality.png", dpi)


def class_distribution(hist: Hist, cent_range: tuple[float, float], cls: str) -> tuple[np.ndarray, np.ndarray]:
    idx = select_bins(hist.edges[1], *cent_range)
    values = np.sum(hist.values[:, idx, CLASS_INDEX[cls]], axis=1)
    variances = np.sum(hist.variances[:, idx, CLASS_INDEX[cls]], axis=1)
    return values, variances


def plot_class_distributions(samples: Iterable[Sample], output: Path, dpi: int) -> None:
    definitions = (
        ("mbd_fired_pmt", "h3_pmtDiag_mbdNFiredTotalByClass", "Fired MBD PMTs"),
        ("mbd_total_charge", "h3_pmtDiag_mbdChargeTotalByClass", r"$\Sigma Q_{MBD}$"),
        ("mbd_charge_asymmetry", "h3_pmtDiag_mbdChargeAsymmetryByClass", r"$(Q_N-Q_S)/(Q_N+Q_S)$"),
    )
    for short, name, xlabel in definitions:
        fig, axes = plt.subplots(2, 3, figsize=(15.2, 8.9), sharex=True, sharey=False)
        for row, sample in enumerate(samples):
            hist = read_hist(sample, name)
            x = centers(hist.edges[0])
            width = np.diff(hist.edges[0])
            for col, cent_range in enumerate(LOW_BAND_GROUPS):
                ax = axes[row, col]
                for cls, color in (("low", "#d62728"), ("main", "#1f77b4")):
                    values, variances = class_distribution(hist, cent_range, cls)
                    integral = float(np.sum(values))
                    if integral <= 0.0:
                        continue
                    density = values / integral / width
                    error = np.sqrt(np.maximum(variances, 0.0)) / integral / width
                    ax.step(x, density, where="mid", color=color, linewidth=1.8, label=cls)
                    ax.fill_between(x, np.maximum(0.0, density - error), density + error,
                                    step="mid", color=color, alpha=0.18)
                ax.set_title(f"{short_label(sample)}, {cent_range[0]:.0f}-{cent_range[1]:.0f}%", fontsize=11)
                ax.set_xlabel(xlabel)
                if col == 0:
                    ax.set_ylabel("Weighted density")
                ax.legend(frameon=False)
        fig.suptitle(f"Low-band and main-band {xlabel} distributions", fontsize=17, fontweight="bold")
        fig.text(0.5, 0.94, "Au+Au embedded simulation; event classes are defined only for centrality below 55%", ha="center")
        fig.subplots_adjust(top=0.87, hspace=0.28, wspace=0.18)
        save(fig, output / f"{short}_low_vs_main.png", dpi)


def geometry(sample: Sample) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    geometry_keys = (
        HIST_PREFIX + "h_pmtDiag_mbdGeometryEntries",
        HIST_PREFIX + "h_pmtDiag_mbdGeometryXSum",
        HIST_PREFIX + "h_pmtDiag_mbdGeometryYSum",
    )
    if not all(key in sample.root for key in geometry_keys):
        return np.arange(128, dtype=float), np.full(128, np.nan), np.zeros(128, dtype=bool)
    n = read_hist(sample, "h_pmtDiag_mbdGeometryEntries").values
    xsum = read_hist(sample, "h_pmtDiag_mbdGeometryXSum").values
    ysum = read_hist(sample, "h_pmtDiag_mbdGeometryYSum").values
    good = n > 0.0
    if np.count_nonzero(good) < 100:
        return np.arange(128, dtype=float), np.full(128, np.nan), np.zeros(128, dtype=bool)
    x = np.full_like(n, np.nan)
    y = np.full_like(n, np.nan)
    x[good] = xsum[good] / n[good]
    y[good] = ysum[good] / n[good]
    return x, y, good


def channel_metric(sample: Sample, cent_range: tuple[float, float], metric: str) -> np.ndarray:
    occ = read_hist(sample, "h3_pmtDiag_mbdPmtOccupancyByChannel")
    charge = read_hist(sample, "h3_pmtDiag_mbdPmtChargeByChannel")
    fired = read_hist(sample, "h3_pmtDiag_mbdNFiredTotalByClass")
    idx = select_bins(occ.edges[1], *cent_range)
    results = []
    for cls in ("low", "main"):
        class_idx = CLASS_INDEX[cls]
        hits = np.sum(occ.values[:, idx, class_idx], axis=1)
        events = float(np.sum(fired.values[:, idx, class_idx]))
        if metric == "occupancy":
            results.append(hits / events if events > 0.0 else np.full(128, np.nan))
        else:
            qsum = np.sum(charge.values[:, idx, class_idx], axis=1)
            results.append(np.divide(qsum, hits, out=np.full(128, np.nan), where=hits > 0.0))
    return np.divide(results[0], results[1], out=np.full(128, np.nan), where=np.asarray(results[1]) > 0.0)


def plot_channel_maps(samples: Iterable[Sample], output: Path, dpi: int) -> None:
    for metric, title in (("occupancy", "Fired-PMT occupancy"), ("mean_charge", "Mean positive PMT charge")):
        fig, axes = plt.subplots(2, 3, figsize=(14.2, 8.8))
        all_ratios = []
        payload = []
        geometry_available = True
        for sample in samples:
            x, y, good = geometry(sample)
            geometry_available &= bool(np.count_nonzero(good) >= 100)
            for cent_range in LOW_BAND_GROUPS:
                ratio = channel_metric(sample, cent_range, metric)
                finite = np.isfinite(ratio) & (ratio > 0.0)
                if np.count_nonzero(good) >= 100:
                    finite &= good
                all_ratios.extend(np.log2(ratio[finite]).tolist())
                payload.append((sample, cent_range, x, y, finite, ratio))
        vmax = max(0.25, float(np.nanpercentile(np.abs(all_ratios), 95))) if all_ratios else 1.0
        norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
        image = None
        for ax, (sample, cent_range, x, y, finite, ratio) in zip(axes.flat, payload):
            if geometry_available:
                image = ax.scatter(x[finite], y[finite], c=np.log2(ratio[finite]), cmap="coolwarm",
                                   norm=norm, marker="h", s=250, linewidths=0.5, edgecolors="black")
                ax.set_aspect("equal")
                ax.set_xlabel("MBD x [cm]")
                ax.set_ylabel("MBD y [cm]")
            else:
                channel = np.arange(ratio.size)
                values = np.log2(ratio)
                ax.plot(channel[finite], values[finite], marker="o", markersize=2.8,
                        linewidth=0.8, color="#1f77b4")
                ax.axhline(0.0, color="0.35", linestyle="--", linewidth=1.0)
                ax.axvline(63.5, color="0.55", linestyle=":", linewidth=1.0)
                ax.set_xlim(-0.5, 127.5)
                ax.set_ylim(-vmax, vmax)
                ax.set_xlabel("MBD PMT channel (arm boundary at 64)")
                ax.set_ylabel(r"$\log_2$(low/main)")
            ax.set_title(f"{short_label(sample)}, {cent_range[0]:.0f}-{cent_range[1]:.0f}%", fontsize=11)
        if geometry_available and image is not None:
            colorbar_axis = fig.add_axes([0.895, 0.16, 0.018, 0.67])
            fig.colorbar(image, cax=colorbar_axis, label=r"$\log_2$(low/main)")
        kind = "channel map" if geometry_available else "channel profile"
        fig.suptitle(f"MBD {kind}: {title} ratio", fontsize=17, fontweight="bold")
        subtitle = "Au+Au embedded simulation; low-band divided by main-band response"
        if not geometry_available:
            subtitle += "; weighted output does not retain PMT x-y metadata"
        fig.text(0.5, 0.94, subtitle, ha="center")
        fig.subplots_adjust(top=0.87, hspace=0.30, wspace=0.24,
                            right=0.86 if geometry_available else 0.96)
        suffix = "maps" if geometry_available else "channels"
        save(fig, output / f"mbd_pmt_{metric}_low_over_main_{suffix}.png", dpi)


def plot_arm_correlations(samples: Iterable[Sample], output: Path, dpi: int) -> None:
    for short, name, xlabel, ylabel in (
        ("fired", "h3_pmtDiag_mbdNFiredSouthVsNorth", r"$N_{fired}^{South}$", r"$N_{fired}^{North}$"),
        ("charge", "h3_pmtDiag_mbdChargeSouthVsNorth", r"$Q_{South}$", r"$Q_{North}$"),
    ):
        fig, axes = plt.subplots(2, 3, figsize=(14.5, 8.9), sharex=True, sharey=True)
        for row, sample in enumerate(samples):
            hist = read_hist(sample, name)
            for col, cent_range in enumerate(CENTRALITY_GROUPS):
                idx = select_bins(hist.edges[2], *cent_range)
                density = np.sum(hist.values[:, :, idx], axis=2)
                draw_density(axes[row, col], density, hist.edges[0], hist.edges[1])
                axes[row, col].set_title(f"{short_label(sample)}, {cent_range[0]:.0f}-{cent_range[1]:.0f}%", fontsize=11)
                axes[row, col].set_xlabel(xlabel)
                if col == 0:
                    axes[row, col].set_ylabel(ylabel)
        fig.suptitle(f"MBD South-North {short} correlation", fontsize=17, fontweight="bold")
        fig.text(0.5, 0.94, "Au+Au embedded simulation", ha="center")
        fig.subplots_adjust(top=0.87, hspace=0.25, wspace=0.10)
        save(fig, output / f"mbd_south_north_{short}_by_centrality.png", dpi)


def plot_sepd(samples: Iterable[Sample], output: Path, dpi: int) -> list[str]:
    available = [s for s in samples if HIST_PREFIX + "h3_pmtDiag_mbdChargeVsSepdCharge" in s.root]
    outputs: list[str] = []
    if len(available) != len(tuple(samples)):
        return outputs
    fig, axes = plt.subplots(2, 3, figsize=(14.8, 8.9), sharex=True, sharey=True)
    for row, sample in enumerate(available):
        hist = read_hist(sample, "h3_pmtDiag_mbdChargeVsSepdCharge")
        for col, cent_range in enumerate(CENTRALITY_GROUPS):
            idx = select_bins(hist.edges[2], *cent_range)
            density = np.sum(hist.values[:, :, idx], axis=2)
            draw_density(axes[row, col], density, hist.edges[0], hist.edges[1])
            axes[row, col].set_title(f"{short_label(sample)}, {cent_range[0]:.0f}-{cent_range[1]:.0f}%", fontsize=11)
            axes[row, col].set_xlabel(r"$\Sigma Q_{MBD}$")
            if col == 0:
                axes[row, col].set_ylabel(r"$\Sigma E_{sEPD}$")
    fig.suptitle("MBD and sEPD correlation", fontsize=17, fontweight="bold")
    fig.text(0.5, 0.94, "Au+Au embedded simulation; calibrated positive detector response", ha="center")
    fig.subplots_adjust(top=0.87, hspace=0.25, wspace=0.10)
    name = "mbd_vs_sepd_by_centrality.png"
    save(fig, output / name, dpi)
    outputs.append(name)
    return outputs


def class_weight(hist: Hist, cent_index: int, class_name: str) -> float:
    return float(np.sum(hist.values[:, cent_index, CLASS_INDEX[class_name]]))


def class_mean(hist: Hist, cent_index: int, class_name: str) -> float:
    values = hist.values[:, cent_index, CLASS_INDEX[class_name]]
    denominator = float(np.sum(values))
    if denominator <= 0.0:
        return math.nan
    return float(np.sum(values * centers(hist.edges[0])) / denominator)


def build_cause_rows(samples: Iterable[Sample]) -> list[dict]:
    rows: list[dict] = []
    for sample in samples:
        fired = read_hist(sample, "h3_pmtDiag_mbdNFiredTotalByClass")
        charge = read_hist(sample, "h3_pmtDiag_mbdChargeTotalByClass")
        by_decision = read_minbias_splits(sample, "h3_pmtDiag_mbdNFiredTotalByClass")
        for index, (low, high) in enumerate(zip(fired.edges[1][:-1], fired.edges[1][1:])):
            if low >= 55.0:
                continue
            all_weight = class_weight(fired, index, "all")
            row = {
                "sample": sample.key,
                "cent_low": float(low),
                "cent_high": float(high),
                "all_weight": all_weight,
                "low_fraction_all": class_weight(fired, index, "low") / all_weight if all_weight > 0.0 else math.nan,
                "mean_fired_pmt_low": class_mean(fired, index, "low"),
                "mean_fired_pmt_main": class_mean(fired, index, "main"),
                "mean_mbd_charge_low": class_mean(charge, index, "low"),
                "mean_mbd_charge_main": class_mean(charge, index, "main"),
            }
            for suffix in MB_DECISIONS:
                decision_hist = by_decision[suffix]
                decision_all = class_weight(decision_hist, index, "all")
                decision_low = class_weight(decision_hist, index, "low")
                row[f"{suffix}_all_weight"] = decision_all
                row[f"low_fraction_{suffix}"] = decision_low / decision_all if decision_all > 0.0 else math.nan
            for class_name in ("low", "main"):
                weights = {
                    suffix: class_weight(hist, index, class_name)
                    for suffix, hist in by_decision.items()
                }
                denominator = sum(weights.values())
                row[f"classifier_pass_fraction_{class_name}"] = (
                    weights["mbPass"] / denominator if denominator > 0.0 else math.nan
                )
            rows.append(row)
    return rows


def finite_ratio(numerator: float, denominator: float) -> float:
    if not math.isfinite(numerator) or not math.isfinite(denominator) or denominator == 0.0:
        return math.nan
    return numerator / denominator


def write_cause_summary(samples: Iterable[Sample], output: Path) -> tuple[list[dict], dict]:
    samples = tuple(samples)
    rows = build_cause_rows(samples)
    with (output / "low_calo_cause_by_sample_centrality.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    summary: dict = {
        "schema": "AUAU_EMBEDDED_LOW_CALO_CAUSE_SUMMARY_V1",
        "centrality_scope_percent": [0, 55],
        "samples": {},
        "interpretation_contract": {
            "classifier_test": "compare low-band fractions and pass fractions for low/main classes",
            "global_mbd_test": "compare mean fired PMTs and total MBD charge for low/main classes",
            "localized_pmt_test": "inspect low/main occupancy and mean-charge channel maps",
            "causality_limit": "diagnostic association only; these histograms do not prove detector-level causation",
        },
    }
    for sample in samples:
        fired = read_hist(sample, "h3_pmtDiag_mbdNFiredTotalByClass")
        charge = read_hist(sample, "h3_pmtDiag_mbdChargeTotalByClass")
        cent_idx = select_bins(fired.edges[1], 0.0, 55.0)
        decision_hists = read_minbias_splits(sample, "h3_pmtDiag_mbdNFiredTotalByClass")

        def summed_weight(hist: Hist, class_name: str) -> float:
            return float(np.sum(hist.values[:, cent_idx, CLASS_INDEX[class_name]]))

        def summed_mean(hist: Hist, class_name: str) -> float:
            values = np.sum(hist.values[:, cent_idx, CLASS_INDEX[class_name]], axis=1)
            denominator = float(np.sum(values))
            return float(np.sum(values * centers(hist.edges[0])) / denominator) if denominator > 0.0 else math.nan

        all_events = summed_weight(fired, "all")
        low_events = summed_weight(fired, "low")
        decision_totals = {suffix: summed_weight(hist, "all") for suffix, hist in decision_hists.items()}
        decision_lows = {suffix: summed_weight(hist, "low") for suffix, hist in decision_hists.items()}
        classifier_by_class = {}
        for class_name in ("low", "main"):
            weights = {suffix: summed_weight(hist, class_name) for suffix, hist in decision_hists.items()}
            denominator = sum(weights.values())
            classifier_by_class[class_name] = {
                "pass_fraction": weights["mbPass"] / denominator if denominator > 0.0 else math.nan,
                "weights": weights,
            }
        channel_summaries = {}
        for cent_range in LOW_BAND_GROUPS:
            key = f"{cent_range[0]:.0f}-{cent_range[1]:.0f}"
            channel_summaries[key] = {}
            for metric in ("occupancy", "mean_charge"):
                ratio = channel_metric(sample, cent_range, metric)
                finite = ratio[np.isfinite(ratio) & (ratio > 0.0)]
                channel_summaries[key][metric] = {
                    "channels": int(finite.size),
                    "median_low_over_main": float(np.median(finite)) if finite.size else math.nan,
                    "p16_low_over_main": float(np.percentile(finite, 16)) if finite.size else math.nan,
                    "p84_low_over_main": float(np.percentile(finite, 84)) if finite.size else math.nan,
                    "fraction_outside_0p8_to_1p2": float(np.mean((finite < 0.8) | (finite > 1.2))) if finite.size else math.nan,
                }
        summary["samples"][sample.key] = {
            "all_weight": all_events,
            "low_weight": low_events,
            "low_fraction_all": low_events / all_events if all_events > 0.0 else math.nan,
            "low_fraction_classifier_pass": finite_ratio(decision_lows["mbPass"], decision_totals["mbPass"]),
            "low_fraction_classifier_fail": finite_ratio(decision_lows["mbFail"], decision_totals["mbFail"]),
            "classifier_by_class": classifier_by_class,
            "mean_fired_pmt_low": summed_mean(fired, "low"),
            "mean_fired_pmt_main": summed_mean(fired, "main"),
            "mean_mbd_charge_low": summed_mean(charge, "low"),
            "mean_mbd_charge_main": summed_mean(charge, "main"),
            "channel_summaries": channel_summaries,
        }
    with (output / "low_calo_cause_summary.json").open("w") as stream:
        json.dump(summary, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return rows, summary


def plot_cause_summary(rows: list[dict], samples: Iterable[Sample], output: Path, dpi: int) -> None:
    samples = tuple(samples)
    colors = {"all": "#111111", "mbPass": "#1f77b4", "mbFail": "#d62728", "low": "#d62728", "main": "#1f77b4"}
    fig, axes = plt.subplots(len(samples), 3, figsize=(15.8, 8.7), sharex=True)
    for row_index, sample in enumerate(samples):
        sample_rows = [row for row in rows if row["sample"] == sample.key]
        x = np.asarray([0.5 * (row["cent_low"] + row["cent_high"]) for row in sample_rows])
        ax = axes[row_index, 0]
        for key, label in (("all", "all evaluated"), ("mbPass", "classifier pass"), ("mbFail", "classifier fail")):
            y_key = "low_fraction_all" if key == "all" else f"low_fraction_{key}"
            ax.plot(x, [row[y_key] for row in sample_rows], marker="o", linewidth=1.7, color=colors[key], label=label)
        ax.set_ylabel("Low-calo-band fraction")
        ax.set_ylim(0.0, 1.05)
        ax.legend(frameon=False, fontsize=9)

        ax = axes[row_index, 1]
        for class_name, label in (("low", "low-calo band"), ("main", "main band")):
            ax.plot(x, [row[f"classifier_pass_fraction_{class_name}"] for row in sample_rows], marker="o",
                    linewidth=1.7, color=colors[class_name], label=label)
        ax.set_ylabel("Classifier pass fraction")
        ax.set_ylim(0.0, 1.05)
        ax.legend(frameon=False, fontsize=9)

        ax = axes[row_index, 2]
        fired_ratio = [finite_ratio(row["mean_fired_pmt_low"], row["mean_fired_pmt_main"]) for row in sample_rows]
        charge_ratio = [finite_ratio(row["mean_mbd_charge_low"], row["mean_mbd_charge_main"]) for row in sample_rows]
        ax.plot(x, fired_ratio, marker="o", linewidth=1.7, color="#2ca02c", label="mean fired PMTs")
        ax.plot(x, charge_ratio, marker="s", linewidth=1.7, color="#9467bd", label="mean MBD charge")
        ax.axhline(1.0, color="0.35", linestyle="--", linewidth=1.0)
        ax.set_ylabel("Low band / main band")
        ax.set_ylim(bottom=0.0)
        ax.legend(frameon=False, fontsize=9)

        for column, title in enumerate(("Band population", "Minimum-bias acceptance", "MBD activity")):
            axes[row_index, column].set_title(f"{short_label(sample)}: {title}", fontsize=11, fontweight="bold")
            axes[row_index, column].set_xlabel("Centrality percentile [%]")
            axes[row_index, column].set_xlim(0.0, 55.0)
    fig.suptitle("Low-calorimeter-band diagnostic summary", fontsize=18, fontweight="bold")
    fig.text(0.5, 0.94, "Au+Au embedded simulation; low-band boundary is diagnostic only; classifier and MBD quantities are evaluated before rejection", ha="center")
    fig.subplots_adjust(top=0.87, hspace=0.30, wspace=0.25)
    save(fig, output / "low_calo_band_cause_summary.png", dpi)


def audit_and_write(samples: Iterable[Sample], output: Path) -> dict:
    manifest: dict = {
        "schema": "AUAU_EMBEDDED_PMT_LOW_CALO_DIAGNOSTICS_V2",
        "cuts": {
            "vertex_abs_cm": 30,
            "centrality_percent": [0, 80],
            "low_band_classification_max_percent": 55,
            "low_band_formula": "log10(max(0,Etotal)+1) < -8.484848484848325e-05*c^2 - 0.006181818181818246*c + 2.5830757575757572",
            "mbd_fired_definition": "finite calibrated q > 0",
            "calorimeter_sum": "finite signed TowerInfo energy; get_isGood required",
            "classification_role": "PMT diagnostics retain pass/fail events before the gate; downstream physics output requires MinimumBiasInfo::isAuAuMinimumBias() in the gated campaign",
            "embedded_minimum_bias_classifier_required": True,
        },
        "weighting": "producer event vertex/centrality weight once, then canonical sample stitch weight once",
        "sources": {
            "signal": ["run28_embeddedPhoton12", "run28_embeddedPhoton20"],
            "inclusive": ["run28_embeddedJet12", "run28_embeddedJet20", "run28_embeddedJet30", "run28_embeddedJet40"],
        },
        "samples": {},
        "qa": {},
    }
    fraction_rows = []
    channel_rows = []
    minbias_rows = []
    for sample in samples:
        audit = read_hist(sample, "h_pmtDiag_audit")
        fired = read_hist(sample, "h3_pmtDiag_mbdNFiredTotalByClass")
        occ = read_hist(sample, "h3_pmtDiag_mbdPmtOccupancyByChannel")
        charge = read_hist(sample, "h3_pmtDiag_mbdPmtChargeByChannel")
        x, y, good_geometry = geometry(sample)
        sample_info = {
            "path": str(sample.path.resolve()),
            "sha256": sha256(sample.path),
            "bytes": sample.path.stat().st_size,
            "audit_bin_values": audit.values.tolist(),
            "sepd_objects_present": HIST_PREFIX + "h3_pmtDiag_mbdChargeVsSepdCharge" in sample.root,
        }
        manifest["samples"][sample.key] = sample_info
        cent_edges = fired.edges[1]
        closure_ok = True
        bounded_ok = True
        for i, (low, high) in enumerate(zip(cent_edges[:-1], cent_edges[1:])):
            if low >= 55.0:
                continue
            all_events = float(np.sum(fired.values[:, i, CLASS_INDEX["all"]]))
            low_events = float(np.sum(fired.values[:, i, CLASS_INDEX["low"]]))
            main_events = float(np.sum(fired.values[:, i, CLASS_INDEX["main"]]))
            fraction = low_events / all_events if all_events > 0.0 else math.nan
            closure = all_events - low_events - main_events
            closure_ok &= abs(closure) <= max(1.0e-8, 1.0e-8 * abs(all_events))
            bounded_ok &= (not math.isfinite(fraction)) or (0.0 <= fraction <= 1.0)
            fraction_rows.append(
                {
                    "sample": sample.key,
                    "cent_low": low,
                    "cent_high": high,
                    "all_weight": all_events,
                    "low_weight": low_events,
                    "main_weight": main_events,
                    "low_fraction": fraction,
                    "closure_residual": closure,
                }
            )
        status_closure = True
        status_closure_max_abs = 0.0
        for name in (
            "h3_pmtDiag_mbdPmtOccupancyByChannel",
            "h3_pmtDiag_mbdPmtChargeByChannel",
            "h3_pmtDiag_mbdNFiredTotalByClass",
            "h3_pmtDiag_mbdChargeTotalByClass",
            "h3_pmtDiag_mbdChargeAsymmetryByClass",
            "h3_pmtDiag_logTotalCaloVsMbdCharge",
            "h3_pmtDiag_totalCaloEnergyVsMbdCharge",
        ):
            inclusive_hist = read_hist(sample, name)
            split_hists = read_minbias_splits(sample, name)
            split_sum = sum((split_hists[suffix].values for suffix in MB_DECISIONS),
                            np.zeros_like(inclusive_hist.values))
            residual = inclusive_hist.values - split_sum
            status_closure_max_abs = max(status_closure_max_abs, float(np.max(np.abs(residual))))
            status_closure &= bool(np.allclose(inclusive_hist.values, split_sum, rtol=2.0e-6, atol=1.0e-5))

        mb_status = read_minbias_splits(sample, "h3_pmtDiag_mbdNFiredTotalByClass")
        for i, (low, high) in enumerate(zip(cent_edges[:-1], cent_edges[1:])):
            status_weights = {
                suffix: float(np.sum(hist.values[:, i, CLASS_INDEX["all"]]))
                for suffix, hist in mb_status.items()
            }
            denominator = sum(status_weights.values())
            minbias_rows.append({
                "sample": sample.key,
                "cent_low": low,
                "cent_high": high,
                "pass_weight": status_weights["mbPass"],
                "fail_weight": status_weights["mbFail"],
                "missing_weight": status_weights["mbMissing"],
                "pass_fraction": status_weights["mbPass"] / denominator if denominator > 0.0 else math.nan,
            })

        manifest["qa"][sample.key] = {
            "class_closure": bool(closure_ok),
            "fractions_bounded": bool(bounded_ok),
            "geometry_available": bool(np.count_nonzero(good_geometry) >= 100),
            "geometry_channels": int(np.count_nonzero(good_geometry)),
            "minimum_bias_status_closure": bool(status_closure),
            "minimum_bias_status_closure_max_abs": status_closure_max_abs,
        }
        for cent_range in LOW_BAND_GROUPS:
            idx = select_bins(occ.edges[1], *cent_range)
            for channel in range(occ.values.shape[0]):
                row = {"sample": sample.key, "cent_low": cent_range[0], "cent_high": cent_range[1], "channel": channel,
                       "x_cm": x[channel], "y_cm": y[channel]}
                for cls in ("low", "main"):
                    class_idx = CLASS_INDEX[cls]
                    events = float(np.sum(fired.values[:, idx, class_idx]))
                    hits = float(np.sum(occ.values[channel, idx, class_idx]))
                    qsum = float(np.sum(charge.values[channel, idx, class_idx]))
                    row[f"{cls}_events"] = events
                    row[f"{cls}_occupancy"] = hits / events if events > 0.0 else math.nan
                    row[f"{cls}_mean_charge"] = qsum / hits if hits > 0.0 else math.nan
                channel_rows.append(row)
    if not all(v["class_closure"] and v["fractions_bounded"]
               and v["minimum_bias_status_closure"]
               for v in manifest["qa"].values()):
        raise RuntimeError(f"numerical QA failed: {manifest['qa']}")
    with (output / "low_band_fraction_by_sample_centrality.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fraction_rows[0].keys())
        writer.writeheader()
        writer.writerows(fraction_rows)
    with (output / "mbd_pmt_channel_metrics.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=channel_rows[0].keys())
        writer.writeheader()
        writer.writerows(channel_rows)
    with (output / "minimum_bias_classifier_by_sample_centrality.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=minbias_rows[0].keys())
        writer.writeheader()
        writer.writerows(minbias_rows)
    with (output / "manifest.json").open("w") as stream:
        json.dump(manifest, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return manifest


def main() -> int:
    args = parse_args()
    style()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for path in (args.signal_root, args.inclusive_root):
        if not path.is_file() or path.stat().st_size < 1024:
            raise FileNotFoundError(f"missing or tiny ROOT input: {path}")
    samples = (
        Sample("signal", "Embedded photon+jet 12+20", args.signal_root, uproot.open(args.signal_root)),
        Sample("inclusive", "Embedded inclusive+jet 12+20+30+40", args.inclusive_root, uproot.open(args.inclusive_root)),
    )
    audit_and_write(samples, args.output_dir)
    cause_rows, _ = write_cause_summary(samples, args.output_dir)
    plot_low_band(samples, args.output_dir, args.dpi)
    plot_zero_based_total_calo(samples, args.output_dir, args.dpi)
    plot_cause_summary(cause_rows, samples, args.output_dir, args.dpi)
    plot_mbd_calo_correlations(samples, args.output_dir, args.dpi)
    plot_minbias_classifier_comparison(samples, args.output_dir, args.dpi)
    plot_old_style_correlations(samples, args.output_dir, args.dpi)
    plot_class_distributions(samples, args.output_dir, args.dpi)
    plot_channel_maps(samples, args.output_dir, args.dpi)
    plot_arm_correlations(samples, args.output_dir, args.dpi)
    plot_sepd(samples, args.output_dir, args.dpi)
    pngs = sorted(path.name for path in args.output_dir.glob("*.png"))
    print(json.dumps({"output_dir": str(args.output_dir.resolve()), "pngs": pngs, "count": len(pngs)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
