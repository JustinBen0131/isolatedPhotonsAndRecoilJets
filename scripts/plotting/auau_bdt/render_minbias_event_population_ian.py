#!/usr/bin/env python3
"""Render plot-only MinimumBiasClassifier event-population diagnostics.

The script projects the registered three-dimensional diagnostic histograms
over MBD charge, combines the photon+jet and inclusive-jet source groups, and
renders the unrestricted and minimum-bias-classifier-pass populations with an
identical logarithmic color scale. It does not rerun reconstruction or alter
the stored histogram contents.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.ticker import LogFormatterMathtext, LogLocator, MultipleLocator
import numpy as np
import uproot


SCRIPT_REPO = Path(__file__).resolve().parents[3]
REPO = (
    SCRIPT_REPO
    if (SCRIPT_REPO / "InputFiles/t95pmtmb_v3_20260714").exists()
    else Path("/Users/patsfan753/Desktop/ThesisAnalysis")
)
DEFAULT_INCLUSIVE = (
    REPO / "InputFiles/t95pmtmb_v3_20260714/axis_safe/inclusive_axis_safe.root"
)
DEFAULT_SIGNAL = (
    REPO / "InputFiles/t95pmtmb_v3_20260714/axis_safe/signal_axis_safe.root"
)
DEFAULT_OUTPUT_DIR = (
    REPO
    / "dataOutput/auauTightBDTValidation/THE95_EmbeddedPmtLowCaloCause_20260714"
    / "ian_standalone"
)
HIST_ALL = "SIM/h3_pmtDiag_totalCaloEnergyVsMbdCharge"
HIST_PASS = "SIM/h3_pmtDiag_totalCaloEnergyVsMbdCharge_mbPass"

KBIRD = LinearSegmentedColormap.from_list(
    "sphenix_kbird",
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
KBIRD.set_bad("white")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inclusive-root", type=Path, default=DEFAULT_INCLUSIVE)
    parser.add_argument("--signal-root", type=Path, default=DEFAULT_SIGNAL)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def project(root_path: Path, key: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with uproot.open(root_path) as root:
        if key not in root:
            raise KeyError(f"missing histogram {key} in {root_path}")
        hist = root[key]
        values = np.asarray(hist.values(flow=False), dtype=float)
        if values.ndim != 3:
            raise ValueError(f"{key} is not three-dimensional: {values.shape}")
        # Producer axes: MBD charge, total calorimeter energy, centrality.
        density = np.sum(values, axis=0).T
        energy_edges = np.asarray(hist.axes[1].edges(flow=False), dtype=float)
        centrality_edges = np.asarray(hist.axes[2].edges(flow=False), dtype=float)
    return density, centrality_edges, energy_edges


def read_combined(
    inclusive_root: Path, signal_root: Path, key: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    inclusive, centrality_edges, energy_edges = project(inclusive_root, key)
    signal, signal_centrality_edges, signal_energy_edges = project(signal_root, key)
    if not np.array_equal(centrality_edges, signal_centrality_edges):
        raise ValueError(f"centrality axes differ for {key}")
    if not np.array_equal(energy_edges, signal_energy_edges):
        raise ValueError(f"energy axes differ for {key}")
    return inclusive + signal, centrality_edges, energy_edges


def configure_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "mathtext.fontset": "dejavusans",
            "axes.linewidth": 1.0,
            "axes.labelsize": 10.5,
            "axes.titlesize": 10.5,
            "xtick.labelsize": 9.0,
            "ytick.labelsize": 9.0,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "xtick.major.size": 5.0,
            "ytick.major.size": 5.0,
            "xtick.minor.size": 2.8,
            "ytick.minor.size": 2.8,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )


def render_panel(
    density: np.ndarray,
    centrality_edges: np.ndarray,
    energy_edges: np.ndarray,
    norm: LogNorm,
    title: str,
    output: Path,
) -> None:
    masked = np.ma.masked_less_equal(density, 0.0)
    fig, ax = plt.subplots(figsize=(3.70, 3.42), dpi=300)
    fig.subplots_adjust(left=0.175, right=0.825, bottom=0.165, top=0.865)

    mesh = ax.pcolormesh(
        centrality_edges,
        energy_edges,
        masked.T,
        cmap=KBIRD,
        norm=norm,
        shading="flat",
        rasterized=True,
    )
    ax.set_xlim(0.0, 80.0)
    ax.set_ylim(0.0, 3000.0)
    ax.set_xlabel("Centrality percentile [%]")
    ax.set_ylabel(r"$E_{\mathrm{calo}}^{\mathrm{total}}$ [GeV]")
    ax.set_title(title, loc="left", fontweight="bold", pad=7.0)
    ax.xaxis.set_major_locator(MultipleLocator(20.0))
    ax.xaxis.set_minor_locator(MultipleLocator(5.0))
    ax.yaxis.set_major_locator(MultipleLocator(500.0))
    ax.yaxis.set_minor_locator(MultipleLocator(100.0))

    ax.text(
        0.035,
        0.965,
        (
            r"$\it{\bf{sPHENIX}}$ Internal"
            "\n"
            "Au+Au embedded simulation"
            "\n"
            r"$\sqrt{s_{NN}}=200$ GeV"
        ),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=7.5,
        linespacing=1.25,
        color="black",
        bbox={
            "boxstyle": "square,pad=0.28",
            "facecolor": "white",
            "edgecolor": "none",
            "alpha": 0.86,
        },
    )

    cax = fig.add_axes([0.845, 0.165, 0.036, 0.700])
    colorbar = fig.colorbar(mesh, cax=cax)
    colorbar.set_label("Weighted events / bin", fontsize=8.3, labelpad=5.0)
    colorbar.ax.tick_params(labelsize=7.6, direction="in", length=3.8)
    colorbar.locator = LogLocator(base=10.0, numticks=8)
    colorbar.formatter = LogFormatterMathtext(base=10.0)
    colorbar.update_ticks()

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300)
    plt.close(fig)


def main() -> int:
    args = parse_args()
    inclusive_root = args.inclusive_root.resolve()
    signal_root = args.signal_root.resolve()
    output_dir = args.output_dir.resolve()
    for path in (inclusive_root, signal_root):
        if not path.exists():
            raise FileNotFoundError(path)

    before, centrality_edges, energy_edges = read_combined(
        inclusive_root, signal_root, HIST_ALL
    )
    after, pass_centrality_edges, pass_energy_edges = read_combined(
        inclusive_root, signal_root, HIST_PASS
    )
    if not np.array_equal(centrality_edges, pass_centrality_edges):
        raise ValueError("before/pass centrality axes differ")
    if not np.array_equal(energy_edges, pass_energy_edges):
        raise ValueError("before/pass energy axes differ")

    positive = np.concatenate(
        [array[np.isfinite(array) & (array > 0.0)] for array in (before, after)]
    )
    if positive.size == 0:
        raise ValueError("both projected populations are empty")
    high = float(np.max(positive))
    low = max(float(np.min(positive)), high * 1.0e-6)
    norm = LogNorm(vmin=low, vmax=high)

    configure_style()
    before_output = output_dir / "minbias_event_population_before.png"
    pass_output = output_dir / "minbias_event_population_pass.png"
    render_panel(
        before,
        centrality_edges,
        energy_edges,
        norm,
        "(a) All evaluated embedded events",
        before_output,
    )
    render_panel(
        after,
        centrality_edges,
        energy_edges,
        norm,
        "(b) Minimum-bias-classifier pass",
        pass_output,
    )

    summary = {
        "schema": "ian_minbias_event_population_render_v1",
        "input_roots": {
            "inclusive": {"path": str(inclusive_root), "sha256": sha256(inclusive_root)},
            "signal": {"path": str(signal_root), "sha256": sha256(signal_root)},
        },
        "histograms": {"before": HIST_ALL, "pass": HIST_PASS},
        "projection": "sum MBD-charge axis; display centrality versus total calorimeter energy",
        "combined_weighted_counts": {
            "before": float(np.sum(before)),
            "pass": float(np.sum(after)),
        },
        "shared_log_norm": {"vmin": low, "vmax": high},
        "display_range": {
            "centrality_percent": [0.0, 80.0],
            "energy_gev": [0.0, 3000.0],
        },
        "outputs": {
            before_output.name: sha256(before_output),
            pass_output.name: sha256(pass_output),
        },
    }
    (output_dir / "render_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(before_output)
    print(pass_output)
    print(output_dir / "render_summary.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
