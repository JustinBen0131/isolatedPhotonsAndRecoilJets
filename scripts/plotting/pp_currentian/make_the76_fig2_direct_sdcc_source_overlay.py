#!/usr/bin/env python3
"""Replot PPG12 Fig. 2 direct-photon truth-isolation points from SDCC ROOT.

The SDCC CSV input should be extracted from:
  /sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_noiso.root
  h_direct_pT_truth_isoET_0

using the same integration convention as ppg12codeGit/efficiencytool/FindTruthETCut.C.
This helper draws both an IAN-style PPG12 source reproduction and an
SDCC-vs-current RecoilJets overlay/ratio.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import uproot


REPO = Path(__file__).resolve().parents[3]
CURRENT_ARTIFACT_POINTER = (
    REPO
    / "dataOutput/current_recoiljets_artifacts/current/pp_sim_photonjet_merged/current.json"
)
DEFAULT_OUTDIR = (
    REPO
    / "dataOutput/ppg12Parity/the76_ppg12_fig24_photonjet_fix_20260702_014217"
    / "fig2_truth_iso_fraction_sdcc_direct"
)
DEFAULT_SDCC_CSV = DEFAULT_OUTDIR / "ppg12_fig2_direct_points_from_sdcc_root.csv"


@dataclass(frozen=True)
class Series:
    key: str
    label: str
    hist_label: str
    color: str
    marker: str


SERIES = [
    Series("pt10_15", r"Direct $\gamma$ $p_T$ [10, 15]", "10_15", "#d948a2", "o"),
    Series("pt15_20", r"Direct $\gamma$ $p_T$ [15, 20]", "15_20", "#2ca02c", "s"),
    Series("pt25_30", r"Direct $\gamma$ $p_T$ [25, 30]", "25_30", "#1f77ff", "^"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sdcc-csv", type=Path, default=DEFAULT_SDCC_CSV)
    parser.add_argument("--current-root", type=Path)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    return parser.parse_args()


def resolve_current_root() -> Path:
    payload = json.loads(CURRENT_ARTIFACT_POINTER.read_text())
    roots = payload.get("root_paths") or []
    if len(roots) != 1:
        raise RuntimeError(f"{CURRENT_ARTIFACT_POINTER} must contain exactly one root_paths entry")
    return Path(roots[0])


def load_sdcc_points(path: Path) -> dict[str, list[dict[str, float | str]]]:
    rows: dict[str, list[dict[str, float | str]]] = {series.key: [] for series in SERIES}
    with path.open() as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            key = row["series"]
            if key not in rows:
                continue
            rows[key].append(
                {
                    "source_root": row["source_root"],
                    "hist_name": row["hist_name"],
                    "pt_low": float(row["pt_low"]),
                    "pt_high": float(row["pt_high"]),
                    "point_index": int(row["point_index"]),
                    "draw_x_gev": float(row["draw_x_gev"]),
                    "cutoff_upper_gev": float(row["cutoff_upper_gev"]),
                    "fraction": float(row["fraction"]),
                }
            )
    for key, values in rows.items():
        if not values:
            raise RuntimeError(f"Missing SDCC Fig.2 direct series in {path}: {key}")
        values.sort(key=lambda item: int(item["point_index"]))
    return rows


def root_like_find_bin(edges: np.ndarray, value: float) -> int:
    nbins = len(edges) - 1
    if value < edges[0]:
        return 0
    if value >= edges[-1]:
        return nbins + 1
    return int(np.searchsorted(edges, value, side="right"))


def cumulative_fraction_like_ppg12(hist: uproot.behaviors.TH1.Histogram, cutoff: float) -> float:
    values, edges = hist.to_numpy(flow=False)
    values = np.asarray(values, dtype=float)
    total = float(np.sum(values))
    if total <= 0.0:
        return float("nan")
    ybinlow = root_like_find_bin(edges, -1.0)
    ybinhigh = root_like_find_bin(edges, cutoff)
    ybinlow = max(1, min(len(values), ybinlow))
    ybinhigh = max(1, min(len(values), ybinhigh))
    if ybinhigh < ybinlow:
        return float("nan")
    numer = float(np.sum(values[ybinlow - 1 : ybinhigh]))
    return numer / total


def load_current_points(
    root_path: Path, sdcc_points: dict[str, list[dict[str, float | str]]]
) -> dict[str, list[dict[str, float | str]]]:
    current: dict[str, list[dict[str, float | str]]] = {}
    with uproot.open(root_path) as root_file:
        for series in SERIES:
            hist_name = (
                f"SIM/h_ppg12_fig2_truthIso_direct_pT_{series.hist_label}_eta07_vz30_r03"
            )
            if hist_name not in root_file:
                raise KeyError(f"Missing current-analysis Fig.2 direct histogram: {hist_name}")
            hist = root_file[hist_name]
            current[series.key] = [
                {
                    "draw_x_gev": float(row["draw_x_gev"]),
                    "cutoff_upper_gev": float(row["cutoff_upper_gev"]),
                    "fraction": cumulative_fraction_like_ppg12(hist, float(row["cutoff_upper_gev"])),
                    "hist_name": hist_name,
                }
                for row in sdcc_points[series.key]
            ]
    return current


def configure_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
            "mathtext.fontset": "dejavusans",
            "axes.linewidth": 1.0,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )


def style_axis_like_ian(ax: plt.Axes) -> None:
    ax.set_xlim(0.0, 20.0)
    ax.set_ylim(0.90, 1.05)
    ax.set_xlabel(r"$E_T^{iso}$ Cutoff [GeV]", fontsize=14)
    ax.set_ylabel("Fraction of Events", fontsize=14)
    ax.set_xticks(np.arange(0, 21, 2))
    ax.set_xticks(np.arange(0, 20.5, 1), minor=True)
    ax.set_yticks(np.arange(0.90, 1.051, 0.02))
    ax.minorticks_on()
    ax.tick_params(which="major", length=6, width=1.0, labelsize=11)
    ax.tick_params(which="minor", length=3, width=0.8)


def draw_source_reproduction(path: Path, sdcc_points: dict[str, list[dict[str, float | str]]]) -> None:
    configure_style()
    fig, ax = plt.subplots(figsize=(5.65, 4.75), dpi=180)
    fig.subplots_adjust(left=0.16, right=0.985, bottom=0.15, top=0.92)
    style_axis_like_ian(ax)
    for series in SERIES:
        rows = sdcc_points[series.key]
        x = [float(row["draw_x_gev"]) for row in rows]
        y = [float(row["fraction"]) for row in rows]
        ax.scatter(
            x,
            y,
            marker=series.marker,
            s=13,
            color=series.color,
            linewidths=0.6,
            label=series.label,
            zorder=3,
        )
    ax.text(0.02, 1.045, r"$\bf{\it{sPHENIX}}$ Internal", transform=ax.transAxes, fontsize=13.5)
    ax.text(0.02, 1.000, "Photon+Jet Samples", transform=ax.transAxes, fontsize=12.0)
    ax.text(0.10, 0.925, r"Pythia, $\sqrt{s}$=200 GeV", transform=ax.transAxes, fontsize=11.0)
    ax.text(0.10, 0.875, r"vtx $|z|<30$ cm, $|\eta^\gamma|<0.7$", transform=ax.transAxes, fontsize=11.0)
    ax.text(0.50, 0.20, "PPG12 SDCC ROOT", transform=ax.transAxes, fontsize=9.5, ha="center")
    ax.legend(loc="upper right", frameon=False, fontsize=8.7, handletextpad=0.3, borderpad=0.2)
    fig.savefig(path)
    plt.close(fig)


def draw_overlay_ratio(
    path: Path,
    sdcc_points: dict[str, list[dict[str, float | str]]],
    current_points: dict[str, list[dict[str, float | str]]],
) -> None:
    configure_style()
    fig = plt.figure(figsize=(6.15, 6.35), dpi=180)
    gs = fig.add_gridspec(2, 1, height_ratios=[3.2, 1.0], hspace=0.04)
    ax = fig.add_subplot(gs[0])
    rax = fig.add_subplot(gs[1], sharex=ax)
    fig.subplots_adjust(left=0.155, right=0.985, bottom=0.12, top=0.93)
    style_axis_like_ian(ax)
    ax.tick_params(labelbottom=False)
    for series in SERIES:
        srows = sdcc_points[series.key]
        crows = current_points[series.key]
        x = np.asarray([float(row["draw_x_gev"]) for row in srows], dtype=float)
        y_sdcc = np.asarray([float(row["fraction"]) for row in srows], dtype=float)
        y_current = np.asarray([float(row["fraction"]) for row in crows], dtype=float)
        ratio = np.divide(y_sdcc, y_current, out=np.full_like(y_sdcc, np.nan), where=y_current > 0.0)
        ax.scatter(
            x,
            y_sdcc,
            marker=series.marker,
            s=22,
            facecolors="none",
            edgecolors=series.color,
            linewidths=1.2,
            zorder=4,
        )
        ax.scatter(
            x,
            y_current,
            marker=series.marker,
            s=14,
            color=series.color,
            linewidths=0.6,
            zorder=3,
        )
        rax.scatter(x, ratio, marker=series.marker, s=15, color=series.color, linewidths=0.6)
    ax.set_xlabel("")
    ax.text(0.02, 0.985, r"$\bf{\it{sPHENIX}}$ Internal", transform=ax.transAxes, fontsize=14.0, va="top")
    ax.text(0.02, 0.935, "Photon+Jet Samples", transform=ax.transAxes, fontsize=12.5, va="top")
    ax.text(0.10, 0.855, r"Pythia, $\sqrt{s}$=200 GeV", transform=ax.transAxes, fontsize=11.5)
    ax.text(0.10, 0.805, r"vtx $|z|<30$ cm, $|\eta^\gamma|<0.7$", transform=ax.transAxes, fontsize=11.5)
    sample_handles = [
        plt.Line2D(
            [0],
            [0],
            marker=series.marker,
            color=series.color,
            markerfacecolor=series.color,
            markersize=4.8,
            linestyle="None",
            label=series.label,
        )
        for series in SERIES
    ]
    source_handles = [
        plt.Line2D([0], [0], marker="o", color="black", markerfacecolor="none", markersize=5.5, linestyle="None", label="PPG12 SDCC"),
        plt.Line2D([0], [0], marker="o", color="black", markerfacecolor="black", markersize=4.5, linestyle="None", label="Current output"),
    ]
    leg1 = ax.legend(handles=sample_handles, loc="upper right", bbox_to_anchor=(0.99, 0.995), frameon=False, fontsize=10.8, handletextpad=0.35, borderpad=0.2)
    ax.add_artist(leg1)
    ax.legend(handles=source_handles, loc="lower right", bbox_to_anchor=(0.99, 0.02), frameon=False, fontsize=10.8, handletextpad=0.35, borderpad=0.2)

    rax.axhline(1.0, color="0.45", linestyle=(0, (4, 4)), linewidth=1.0)
    rax.set_xlim(0.0, 20.0)
    rax.set_ylim(0.988, 1.012)
    rax.set_xlabel(r"$E_T^{iso}$ Cutoff [GeV]", fontsize=14)
    rax.set_ylabel("SDCC / current", fontsize=11.5)
    rax.set_xticks(np.arange(0, 21, 2))
    rax.set_xticks(np.arange(0, 20.5, 1), minor=True)
    rax.minorticks_on()
    rax.tick_params(which="major", length=6, width=1.0, labelsize=10)
    rax.tick_params(which="minor", length=3, width=0.8)
    fig.savefig(path)
    plt.close(fig)


def write_points_csv(
    path: Path,
    sdcc_points: dict[str, list[dict[str, float | str]]],
    current_points: dict[str, list[dict[str, float | str]]],
) -> dict[str, dict[str, float]]:
    ratio_ranges: dict[str, dict[str, float]] = {}
    with path.open("w", newline="") as handle:
        fields = [
            "series",
            "label",
            "point_index",
            "draw_x_gev",
            "cutoff_upper_gev",
            "sdcc_fraction",
            "current_fraction",
            "sdcc_over_current",
            "sdcc_source_root",
            "sdcc_hist_name",
            "current_hist_name",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for series in SERIES:
            ratios: list[float] = []
            for srow, crow in zip(sdcc_points[series.key], current_points[series.key]):
                sdcc = float(srow["fraction"])
                current = float(crow["fraction"])
                ratio = sdcc / current if current > 0.0 else float("nan")
                ratios.append(ratio)
                writer.writerow(
                    {
                        "series": series.key,
                        "label": series.label,
                        "point_index": int(srow["point_index"]),
                        "draw_x_gev": float(srow["draw_x_gev"]),
                        "cutoff_upper_gev": float(srow["cutoff_upper_gev"]),
                        "sdcc_fraction": sdcc,
                        "current_fraction": current,
                        "sdcc_over_current": ratio,
                        "sdcc_source_root": srow["source_root"],
                        "sdcc_hist_name": srow["hist_name"],
                        "current_hist_name": crow["hist_name"],
                    }
                )
            finite = [value for value in ratios if math.isfinite(value)]
            ratio_ranges[series.key] = {
                "min": min(finite),
                "max": max(finite),
                "max_abs_delta_from_unity": max(abs(value - 1.0) for value in finite),
            }
    return ratio_ranges


def main() -> None:
    args = parse_args()
    current_root = args.current_root or resolve_current_root()
    args.outdir.mkdir(parents=True, exist_ok=True)
    sdcc_points = load_sdcc_points(args.sdcc_csv)
    current_points = load_current_points(current_root, sdcc_points)

    source_png = args.outdir / "ppg12_fig2_direct_sdcc_root_reproduction_ian_style.png"
    overlay_png = args.outdir / "ppg12_fig2_direct_sdcc_vs_current_overlay_ratio.png"
    points_csv = args.outdir / "ppg12_fig2_direct_sdcc_vs_current_points.csv"
    manifest_path = args.outdir / "ppg12_fig2_direct_sdcc_vs_current_manifest.json"

    draw_source_reproduction(source_png, sdcc_points)
    draw_overlay_ratio(overlay_png, sdcc_points, current_points)
    ratio_ranges = write_points_csv(points_csv, sdcc_points, current_points)

    first = next(iter(sdcc_points.values()))[0]
    manifest = {
        "artifact": "THE-76 PPG12 Fig.2 direct-photon truth-isolation SDCC source overlay",
        "ppg12_macro": "ppg12codeGit/efficiencytool/FindTruthETCut.C",
        "ppg12_remote_source_root": first["source_root"],
        "ppg12_histogram": first["hist_name"],
        "current_root": str(current_root),
        "current_histogram_pattern": "SIM/h_ppg12_fig2_truthIso_direct_pT_{10_15,15_20,25_30}_eta07_vz30_r03",
        "sdcc_points_csv": str(args.sdcc_csv),
        "points_csv": str(points_csv),
        "source_reproduction_png": str(source_png),
        "overlay_ratio_png": str(overlay_png),
        "fraction_definition": "Direct PPG12 Fig.2 convention from FindTruthETCut.C: pT windows [10,15], [15,20], [25,30], isolow=-1, cutoff bin upper edges 1..20 GeV; normal TH2 bins only.",
        "ratio_definition": "PPG12 SDCC ROOT fraction / current RecoilJets fraction",
        "ratio_ranges": ratio_ranges,
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
