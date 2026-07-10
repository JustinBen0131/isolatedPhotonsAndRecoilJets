#!/usr/bin/env python3
"""Plot AuAu embedded-photon truth-isolation diagnostics from RecoilJets output."""

from __future__ import annotations

import argparse
import csv
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import uproot


CENTRALITIES = (
    ("cent_0_20", "0-20%"),
    ("cent_20_50", "20-50%"),
    ("cent_50_80", "50-80%"),
)
PT_INTERVALS = (
    ("pT_10_15", "[10, 15) GeV", "o", "#d94f9f"),
    ("pT_15_20", "[15, 20) GeV", "s", "#2f8f2f"),
    ("pT_25_30", "[25, 30) GeV", "^", "#2878d0"),
)
SPECTRUM_STYLES = {
    "total": ("Total", "o", "#d6279f"),
    "direct": ("Direct", "s", "#477a1f"),
    "fragmentation": ("Fragmentation", "^", "#2474e5"),
}
CUTOFFS = np.arange(1.0, 21.0, 1.0)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path, help="Merged photon12+20 ROOT file")
    parser.add_argument("--outdir", required=True, type=Path, help="Output directory")
    return parser.parse_args()


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.size": 11,
            "axes.labelsize": 12,
            "axes.titlesize": 13,
            "legend.fontsize": 9.5,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "axes.linewidth": 1.1,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def sph_label(ax: plt.Axes, y: float = 0.97) -> None:
    ax.text(
        0.035,
        y,
        r"$\bf{\it{sPHENIX}}$ Internal",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=11,
    )


def object_lookup(root_file: uproot.ReadOnlyDirectory) -> dict[str, str]:
    lookup: dict[str, str] = {}
    for key in root_file.keys(recursive=True, cycle=False):
        base = key.rsplit("/", 1)[-1]
        lookup.setdefault(base, key)
    return lookup


def get_hist(
    root_file: uproot.ReadOnlyDirectory,
    lookup: dict[str, str],
    name: str,
) -> Any:
    if name not in lookup:
        raise KeyError(f"missing ROOT object: {name}")
    return root_file[lookup[name]]


def values_and_variances(hist: Any, flow: bool = False) -> tuple[np.ndarray, np.ndarray]:
    values = np.asarray(hist.values(flow=flow), dtype=float)
    variances_raw = hist.variances(flow=flow)
    if variances_raw is None:
        raise ValueError(f"histogram {hist.name} has no Sumw2/variance information")
    variances = np.asarray(variances_raw, dtype=float)
    if values.shape != variances.shape:
        raise ValueError(f"value/variance shape mismatch for {hist.name}")
    if not np.all(np.isfinite(values)) or not np.all(np.isfinite(variances)):
        raise ValueError(f"non-finite content in {hist.name}")
    if np.any(variances < 0.0):
        raise ValueError(f"negative variance in {hist.name}")
    return values, variances


def subset_fraction(
    selected: np.ndarray,
    selected_var: np.ndarray,
    rejected: np.ndarray,
    rejected_var: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    denominator = selected + rejected
    fraction = np.full_like(denominator, np.nan, dtype=float)
    error = np.full_like(denominator, np.nan, dtype=float)
    valid = denominator > 0.0
    fraction[valid] = selected[valid] / denominator[valid]
    variance = np.zeros_like(denominator, dtype=float)
    variance[valid] = (
        rejected[valid] ** 2 * selected_var[valid]
        + selected[valid] ** 2 * rejected_var[valid]
    ) / denominator[valid] ** 4
    error[valid] = np.sqrt(np.maximum(variance[valid], 0.0))
    return fraction, error


def cumulative_fraction(hist: Any) -> tuple[np.ndarray, np.ndarray]:
    values, variances = values_and_variances(hist, flow=True)
    if values.size < 3:
        raise ValueError(f"unexpected flow-bin layout for {hist.name}")
    edges = np.asarray(hist.axis().edges(flow=False), dtype=float)
    normal_values = values[1:-1]
    normal_variances = variances[1:-1]
    upper_edges = edges[1:]
    denominator = float(np.sum(values))
    denominator_var = float(np.sum(variances))
    if denominator <= 0.0:
        raise ValueError(f"empty denominator for {hist.name}")

    fractions: list[float] = []
    errors: list[float] = []
    for cutoff in CUTOFFS:
        passed_mask = upper_edges <= cutoff + 1.0e-9
        passed = float(values[0] + np.sum(normal_values[passed_mask]))
        passed_var = float(variances[0] + np.sum(normal_variances[passed_mask]))
        failed = denominator - passed
        failed_var = max(denominator_var - passed_var, 0.0)
        fraction, error = subset_fraction(
            np.asarray([passed]),
            np.asarray([passed_var]),
            np.asarray([failed]),
            np.asarray([failed_var]),
        )
        fractions.append(float(fraction[0]))
        errors.append(float(error[0]))

    fraction_array = np.asarray(fractions)
    error_array = np.asarray(errors)
    if not np.all(np.isfinite(fraction_array)) or not np.all(np.isfinite(error_array)):
        raise ValueError(f"non-finite cumulative fraction for {hist.name}")
    if np.any(fraction_array < -1.0e-10) or np.any(fraction_array > 1.0 + 1.0e-10):
        raise ValueError(f"fraction outside [0,1] for {hist.name}")
    if np.any(np.diff(fraction_array) < -1.0e-10):
        raise ValueError(f"non-monotonic cumulative fraction for {hist.name}")
    return fraction_array, error_array


def fraction_data(
    root_file: uproot.ReadOnlyDirectory,
    lookup: dict[str, str],
    photon_class: str,
) -> tuple[dict[str, dict[str, tuple[np.ndarray, np.ndarray]]], list[str]]:
    data: dict[str, dict[str, tuple[np.ndarray, np.ndarray]]] = {}
    names: list[str] = []
    for cent_tag, _ in CENTRALITIES:
        data[cent_tag] = {}
        for pt_tag, _, _, _ in PT_INTERVALS:
            name = f"h_auauTruthIso_{photon_class}_{pt_tag}_{cent_tag}"
            names.append(name)
            data[cent_tag][pt_tag] = cumulative_fraction(get_hist(root_file, lookup, name))
    return data, names


def spectrum_data(
    root_file: uproot.ReadOnlyDirectory,
    lookup: dict[str, str],
) -> tuple[dict[str, dict[str, tuple[np.ndarray, np.ndarray]]], list[str]]:
    data: dict[str, dict[str, tuple[np.ndarray, np.ndarray]]] = {}
    names: list[str] = []
    for cent_tag, _ in CENTRALITIES:
        data[cent_tag] = {}
        for photon_class in SPECTRUM_STYLES:
            name = f"h_auauTruthPt_{photon_class}_iso4_{cent_tag}"
            names.append(name)
            hist = get_hist(root_file, lookup, name)
            values, variances = values_and_variances(hist)
            edges = np.asarray(hist.axis().edges(flow=False), dtype=float)
            if values.size != 25 or not np.allclose(edges, np.arange(10.0, 36.0, 1.0)):
                raise ValueError(f"wrong pp pT binning for {name}")
            data[cent_tag][photon_class] = values, variances

        total = data[cent_tag]["total"][0]
        direct = data[cent_tag]["direct"][0]
        fragmentation = data[cent_tag]["fragmentation"][0]
        if not np.allclose(total, direct + fragmentation, rtol=1.0e-9, atol=1.0e-9):
            difference = float(np.max(np.abs(total - direct - fragmentation)))
            raise ValueError(f"total != direct + fragmentation for {cent_tag}; max diff={difference}")
    return data, names


def format_fraction_axis(
    ax: plt.Axes,
    title: str,
    photon_class: str,
    ylabel: bool = True,
) -> None:
    ax.set_title(title, pad=9)
    ax.set_xlim(0.5, 20.5)
    ax.set_xticks(np.arange(2, 21, 2))
    ax.set_ylim(0.90 if photon_class == "direct" else 0.0, 1.015)
    ax.set_xlabel(r"Truth $E_{T}^{iso}$ cutoff [GeV]")
    if ylabel:
        ax.set_ylabel("Fraction passing cutoff")
    ax.grid(axis="y", color="#d9d9d9", linewidth=0.7, alpha=0.75)


def draw_fraction_panel(
    ax: plt.Axes,
    series: dict[str, tuple[np.ndarray, np.ndarray]],
    photon_class: str,
    cent_label: str,
    ylabel: bool = True,
) -> None:
    for pt_tag, pt_label, marker, color in PT_INTERVALS:
        fractions, errors = series[pt_tag]
        label = pt_label
        if pt_tag == "pT_10_15":
            label += " (coverage starts at 12 GeV)"
        ax.errorbar(
            CUTOFFS,
            fractions,
            yerr=errors,
            linestyle="none",
            marker=marker,
            markersize=4.5,
            capsize=1.8,
            color=color,
            label=label,
        )
    format_fraction_axis(
        ax,
        f"AuAu embedded photon simulation, {cent_label}\n{photon_class.capitalize()} photons",
        photon_class,
        ylabel=ylabel,
    )
    sph_label(ax, y=0.24 if photon_class == "direct" else 0.97)
    ax.text(
        0.035,
        0.14 if photon_class == "direct" else 0.865,
        (
            r"$|\eta^\gamma|<0.7$, $|z_{vtx}|<30$ cm" + "\n" + r"$R=0.3$"
            if photon_class == "direct"
            else r"$|\eta^\gamma|<0.7$, $|z_{vtx}|<30$ cm, $R=0.3$"
        ),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9.5,
    )
    ax.legend(loc="lower right", frameon=False)


def draw_spectrum_panel(
    top: plt.Axes,
    ratio: plt.Axes,
    series: dict[str, tuple[np.ndarray, np.ndarray]],
    cent_label: str,
    ylabel: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    centers = np.arange(10.5, 35.0, 1.0)
    for photon_class, (label, marker, color) in SPECTRUM_STYLES.items():
        values, variances = series[photon_class]
        top.errorbar(
            centers,
            values,
            yerr=np.sqrt(variances),
            linestyle="none",
            marker=marker,
            markersize=4.2,
            capsize=1.5,
            color=color,
            label=label,
        )

    direct, direct_var = series["direct"]
    fragmentation, fragmentation_var = series["fragmentation"]
    fraction, error = subset_fraction(direct, direct_var, fragmentation, fragmentation_var)
    if np.any((fraction[np.isfinite(fraction)] < 0.0) | (fraction[np.isfinite(fraction)] > 1.0)):
        raise ValueError(f"direct/total outside [0,1] for {cent_label}")

    top.axvspan(10.0, 12.0, color="#f2c94c", alpha=0.22, linewidth=0.0)
    ratio.axvspan(10.0, 12.0, color="#f2c94c", alpha=0.22, linewidth=0.0)
    top.set_yscale("log")
    positive = np.concatenate([series[key][0][series[key][0] > 0.0] for key in SPECTRUM_STYLES])
    if positive.size:
        top.set_ylim(max(float(np.min(positive)) * 0.45, 1.0e-8), float(np.max(positive)) * 2.5)
    top.set_xlim(10.0, 35.0)
    top.set_title(f"AuAu embedded photon simulation, {cent_label}", pad=9)
    if ylabel:
        top.set_ylabel("Weighted photons / GeV")
        ratio.set_ylabel("Direct / total")
    top.legend(loc="upper right", frameon=False)
    top.grid(axis="y", which="both", color="#e0e0e0", linewidth=0.6, alpha=0.65)
    top.tick_params(labelbottom=False)
    sph_label(top)
    top.text(
        0.34,
        0.96,
        r"$E_{T}^{iso}<4$ GeV, $R=0.3$" + "\n" + r"$|\eta^\gamma|<0.7$, $|z_{vtx}|<30$ cm",
        transform=top.transAxes,
        ha="left",
        va="top",
        fontsize=9.2,
    )
    top.text(
        0.08,
        0.09,
        "10-12 GeV: source coverage limited",
        transform=top.transAxes,
        fontsize=8.8,
        color="#6c5600",
    )

    ratio.errorbar(
        centers,
        fraction,
        yerr=error,
        linestyle="none",
        marker="s",
        markersize=3.8,
        capsize=1.5,
        color=SPECTRUM_STYLES["direct"][2],
    )
    ratio.set_xlim(10.0, 35.0)
    ratio.set_ylim(0.0, 1.02)
    ratio.set_xlabel(r"Truth photon $p_T$ [GeV]")
    ratio.grid(axis="y", color="#d9d9d9", linewidth=0.7, alpha=0.75)
    return fraction, error


def write_fraction_plots(
    outdir: Path,
    photon_class: str,
    data: dict[str, dict[str, tuple[np.ndarray, np.ndarray]]],
) -> list[Path]:
    outputs: list[Path] = []
    for cent_tag, cent_label in CENTRALITIES:
        fig, ax = plt.subplots(figsize=(7.2, 5.6), constrained_layout=True)
        draw_fraction_panel(ax, data[cent_tag], photon_class, cent_label)
        path = outdir / f"auau_embedded_truthiso_{photon_class}_{cent_tag}.png"
        fig.savefig(path, dpi=180)
        plt.close(fig)
        outputs.append(path)

    fig, axes = plt.subplots(1, 3, figsize=(17.0, 5.25), constrained_layout=True)
    for index, ((cent_tag, cent_label), ax) in enumerate(zip(CENTRALITIES, axes)):
        draw_fraction_panel(
            ax,
            data[cent_tag],
            photon_class,
            cent_label,
            ylabel=index == 0,
        )
    summary = outdir / f"auau_embedded_truthiso_{photon_class}_centrality_summary.png"
    fig.savefig(summary, dpi=180)
    plt.close(fig)
    outputs.append(summary)
    return outputs


def write_spectrum_plots(
    outdir: Path,
    data: dict[str, dict[str, tuple[np.ndarray, np.ndarray]]],
) -> tuple[list[Path], dict[str, tuple[np.ndarray, np.ndarray]]]:
    outputs: list[Path] = []
    ratios: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for cent_tag, cent_label in CENTRALITIES:
        fig = plt.figure(figsize=(7.2, 6.8), constrained_layout=True)
        grid = fig.add_gridspec(2, 1, height_ratios=(3.1, 1.0), hspace=0.03)
        top = fig.add_subplot(grid[0])
        ratio = fig.add_subplot(grid[1], sharex=top)
        ratios[cent_tag] = draw_spectrum_panel(top, ratio, data[cent_tag], cent_label)
        path = outdir / f"auau_embedded_truthpt_spectrum_{cent_tag}.png"
        fig.savefig(path, dpi=180)
        plt.close(fig)
        outputs.append(path)

    fig = plt.figure(figsize=(17.0, 6.5), constrained_layout=True)
    grid = fig.add_gridspec(2, 3, height_ratios=(3.1, 1.0), hspace=0.03, wspace=0.08)
    for index, (cent_tag, cent_label) in enumerate(CENTRALITIES):
        top = fig.add_subplot(grid[0, index])
        ratio = fig.add_subplot(grid[1, index], sharex=top)
        draw_spectrum_panel(top, ratio, data[cent_tag], cent_label, ylabel=index == 0)
    summary = outdir / "auau_embedded_truthpt_spectrum_centrality_summary.png"
    fig.savefig(summary, dpi=180)
    plt.close(fig)
    outputs.append(summary)
    return outputs, ratios


def csv_rows(
    fraction_sets: dict[str, dict[str, dict[str, tuple[np.ndarray, np.ndarray]]]],
    spectra: dict[str, dict[str, tuple[np.ndarray, np.ndarray]]],
    ratios: dict[str, tuple[np.ndarray, np.ndarray]],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for photon_class, class_data in fraction_sets.items():
        for cent_tag, cent_label in CENTRALITIES:
            for pt_tag, pt_label, _, _ in PT_INTERVALS:
                values, errors = class_data[cent_tag][pt_tag]
                for cutoff, value, error in zip(CUTOFFS, values, errors):
                    rows.append(
                        {
                            "figure_family": "truth_iso_fraction",
                            "centrality": cent_label,
                            "photon_class": photon_class,
                            "series": pt_label,
                            "x": cutoff,
                            "x_unit": "GeV",
                            "value": value,
                            "stat_error": error,
                            "coverage_limited": pt_tag == "pT_10_15",
                        }
                    )

    centers = np.arange(10.5, 35.0, 1.0)
    for cent_tag, cent_label in CENTRALITIES:
        for photon_class in SPECTRUM_STYLES:
            values, variances = spectra[cent_tag][photon_class]
            for center, value, variance in zip(centers, values, variances):
                rows.append(
                    {
                        "figure_family": "truth_pt_spectrum",
                        "centrality": cent_label,
                        "photon_class": photon_class,
                        "series": photon_class,
                        "x": center,
                        "x_unit": "GeV",
                        "value": value,
                        "stat_error": math.sqrt(max(float(variance), 0.0)),
                        "coverage_limited": center < 12.0,
                    }
                )
        values, errors = ratios[cent_tag]
        for center, value, error in zip(centers, values, errors):
            rows.append(
                {
                    "figure_family": "truth_pt_direct_fraction",
                    "centrality": cent_label,
                    "photon_class": "direct_over_total",
                    "series": "direct/total",
                    "x": center,
                    "x_unit": "GeV",
                    "value": value,
                    "stat_error": error,
                    "coverage_limited": center < 12.0,
                }
            )
    return rows


def main() -> int:
    args = parse_args()
    root_path = args.root.expanduser().resolve()
    outdir = args.outdir.expanduser().resolve()
    if not root_path.is_file():
        raise FileNotFoundError(root_path)
    outdir.mkdir(parents=True, exist_ok=True)
    setup_style()

    with uproot.open(root_path) as root_file:
        lookup = object_lookup(root_file)
        direct, direct_names = fraction_data(root_file, lookup, "direct")
        fragmentation, fragmentation_names = fraction_data(root_file, lookup, "fragmentation")
        spectra, spectrum_names = spectrum_data(root_file, lookup)
        audit_name = "h_auauTruthIsoDiag_audit"
        audit = get_hist(root_file, lookup, audit_name)
        audit_values, _ = values_and_variances(audit)
        if audit_values.size != 16:
            raise ValueError("truth-isolation audit histogram must have 16 bins")
        if not np.isclose(audit_values[0], np.sum(audit_values[13:16]), rtol=1.0e-9, atol=1.0e-9):
            raise ValueError("accepted-event audit count does not equal centrality-bin event sum")
        if audit_values[2] <= 0.0 or audit_values[3] <= 0.0:
            raise ValueError("direct or fragmentation audit population is zero")

    pngs: list[Path] = []
    pngs.extend(write_fraction_plots(outdir, "direct", direct))
    pngs.extend(write_fraction_plots(outdir, "fragmentation", fragmentation))
    spectrum_pngs, ratios = write_spectrum_plots(outdir, spectra)
    pngs.extend(spectrum_pngs)

    rows = csv_rows(
        {"direct": direct, "fragmentation": fragmentation},
        spectra,
        ratios,
    )
    csv_path = outdir / "auau_embedded_truthiso_diagnostics_values.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    object_names = direct_names + fragmentation_names + spectrum_names + [audit_name]
    manifest = {
        "schema_version": 1,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "input_root": str(root_path),
        "selection": {
            "sample_scope": ["run28_embeddedPhoton12", "run28_embeddedPhoton20"],
            "centrality_percent": [[0, 20], [20, 50], [50, 80]],
            "truth_photon_eta_abs_max": 0.7,
            "vertex_z_abs_max_cm": 30.0,
            "truth_isolation_cone_r": 0.3,
            "truth_photon_removal_dr": 0.001,
            "truth_pt_intervals_gev": [[10, 15], [15, 20], [25, 30]],
            "truth_pt_spectrum_edges_gev": list(range(10, 36)),
            "truth_pt_spectrum_isolation_max_gev": 4.0,
            "fraction_cutoffs_gev": CUTOFFS.tolist(),
        },
        "weighting": {
            "histogram_contract": "RecoilJets RJMCWeighting weighted TH1F with Sumw2",
            "merge_contract": "photon12 and photon20 source weights applied by established embedded weighting path",
        },
        "coverage_limitation": {
            "range_gev": [10, 12],
            "reason": "Available focused signal production begins with run28_embeddedPhoton12.",
            "affected_fraction_series": "[10,15) GeV",
        },
        "root_objects": object_names,
        "validation": {
            "all_objects_present": True,
            "fractions_finite_bounded_monotonic": True,
            "truth_pt_total_equals_direct_plus_fragmentation": True,
            "accepted_events_equal_centrality_event_sum": True,
            "direct_and_fragmentation_nonzero": True,
        },
        "products": {
            "png": [str(path) for path in pngs],
            "csv": str(csv_path),
        },
        "comparison_policy": "AuAu embedded-only; no PPG12 input, overlay, normalization, or comparison ratio.",
    }
    manifest_path = outdir / "auau_embedded_truthiso_diagnostics_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    if len(pngs) != 12:
        raise RuntimeError(f"expected 12 PNG products, created {len(pngs)}")
    print(f"wrote {len(pngs)} PNG files")
    print(f"wrote {csv_path}")
    print(f"wrote {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
