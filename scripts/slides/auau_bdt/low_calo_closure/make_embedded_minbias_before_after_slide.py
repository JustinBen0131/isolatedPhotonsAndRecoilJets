#!/usr/bin/env python3
"""Build the audience-facing embedded MinBiasClassifier before/after slide."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm
import numpy as np
import uproot


REPO = Path(__file__).resolve().parents[4]
DEFAULT_ROOT = REPO / "InputFiles/t95pmtmb_v3_20260714/axis_safe/inclusive_axis_safe.root"
DEFAULT_SIGNAL_ROOT = REPO / "InputFiles/t95pmtmb_v3_20260714/axis_safe/signal_axis_safe.root"
DEFAULT_OUTPUT = (
    REPO
    / "dataOutput/auauTightBDTValidation/THE95_EmbeddedPmtLowCaloCause_20260714/slides"
    / "embedded_minbias_before_after_calo_centrality.png"
)

WIDTH = 2560
HEIGHT = 1440
DPI = 160
HIST_ALL = "SIM/h3_pmtDiag_totalCaloEnergyVsMbdCharge"
HIST_PASS = "SIM/h3_pmtDiag_totalCaloEnergyVsMbdCharge_mbPass"

INK = "#111827"
MUTED = "#475569"
ACCENT = "#b42318"
# Shared JSTG subtitle-arrowhead convention (THE-111 slides 11/12): literal
# DejaVu Sans "▶" in #2468A8 with a 0.021-figure-width hanging indent.
BLUE_BULLET = "#2468A8"
BULLET_INDENT = 0.021
KBIRD = LinearSegmentedColormap.from_list(
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--signal-root", type=Path, default=DEFAULT_SIGNAL_ROOT)
    parser.add_argument("--sample-mode", choices=("inclusive", "combined"), default="inclusive")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_density(root: uproot.ReadOnlyFile, key: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if key not in root:
        raise KeyError(f"missing required histogram: {key}")
    hist = root[key]
    values = np.asarray(hist.values(flow=False), dtype=float)
    if values.ndim != 3:
        raise ValueError(f"{key} must be three-dimensional, found shape {values.shape}")
    # Producer axes are MBD charge, total calorimeter energy, centrality.
    density = np.sum(values, axis=0).T
    energy_edges = np.asarray(hist.axes[1].edges(flow=False), dtype=float)
    centrality_edges = np.asarray(hist.axes[2].edges(flow=False), dtype=float)
    return density, centrality_edges, energy_edges


def shared_norm(*densities: np.ndarray) -> LogNorm:
    positive = np.concatenate(
        [density[np.isfinite(density) & (density > 0.0)] for density in densities]
    )
    if positive.size == 0:
        raise ValueError("both classifier views are empty")
    high = float(np.max(positive))
    low = max(float(np.min(positive)), high * 1.0e-6)
    return LogNorm(vmin=low, vmax=high)


def configure_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "axes.linewidth": 1.6,
            "axes.grid": False,
            "xtick.major.width": 1.5,
            "ytick.major.width": 1.5,
            "xtick.major.size": 8,
            "ytick.major.size": 8,
            "savefig.facecolor": "white",
        }
    )


def add_internal_label(ax: plt.Axes, *, right_aligned: bool = False) -> None:
    color = INK if right_aligned else "white"
    effect = [] if right_aligned else [path_effects.withStroke(linewidth=1.8, foreground="black")]
    x = 0.965 if right_aligned else 0.035
    ha = "right" if right_aligned else "left"
    ax.text(
        x,
        0.955,
        r"$\it{\bf{sPHENIX}}$ Internal",
        transform=ax.transAxes,
        color=color,
        fontsize=15.5,
        va="top",
        ha=ha,
        path_effects=effect,
    )
    ax.text(
        x,
        0.905,
        r"embedded PYTHIA $\sqrt{s_{NN}}=200$ GeV",
        transform=ax.transAxes,
        color=color,
        fontsize=13.5,
        va="top",
        ha=ha,
        path_effects=effect,
    )


def draw_panel(
    ax: plt.Axes,
    density: np.ndarray,
    centrality_edges: np.ndarray,
    energy_edges: np.ndarray,
    norm: LogNorm,
    title: str,
    right_aligned_internal_label: bool = False,
    internal_label: bool = True,
) -> matplotlib.collections.QuadMesh:
    mesh = ax.pcolormesh(
        centrality_edges,
        energy_edges,
        density.T,
        shading="auto",
        cmap=KBIRD,
        norm=norm,
        rasterized=True,
    )
    ax.set_xlim(0.0, 80.0)
    ax.set_ylim(0.0, min(3000.0, float(energy_edges[-1])))
    ax.set_title(title, fontsize=25, fontweight="bold", pad=16, color=INK)
    ax.set_xlabel("Centrality percentile [%]", fontsize=22, labelpad=10)
    ax.tick_params(axis="both", which="major", labelsize=18)
    if internal_label:
        add_internal_label(ax, right_aligned=right_aligned_internal_label)
    return mesh


def write_layout(path: Path, sample_line: str) -> None:
    payload = {
        "schema": "slide_layout_nodes_v1",
        "title_axis_x": 110,
        "minimum_audience_font_px": 43,
        "minimum_title_font_px": 74,
        "minimum_plot_annotation_font_px": 32,
        "nodes": [
            {
                "name": "title",
                "kind": "text",
                "role": "title",
                "text": "Minimum-bias classification removes the low-energy band",
                "font_px": 82,
                "bbox": [110, 54, 2440, 144],
                "title_anchor": True,
            },
            {
                "name": "sample_line",
                "kind": "text",
                "role": "body",
                "text": sample_line,
                "font_px": 45,
                "bbox": [110, 162, 2440, 218],
                "title_axis_align": "left",
            },
            {
                "name": "selection_readout",
                "kind": "text",
                "role": "body",
                "text": "Simulation requirement: keep events when MinimumBiasInfo::isAuAuMinimumBias() returns true.",
                "font_px": 45,
                "bbox": [110, 225, 2440, 286],
                "title_axis_align": "left",
            },
            {
                "name": "left_panel",
                "kind": "plot",
                "role": "plot",
                "bbox": [112, 333, 1235, 1365],
                "title_axis_align": "left",
            },
            {
                "name": "right_panel",
                "kind": "plot",
                "role": "plot",
                "bbox": [1310, 333, 2360, 1365],
            },
        ],
    }
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def write_script(path: Path, sample_description: str) -> None:
    path.write_text(
        "\n".join(
            [
                "# Embedded minimum-bias classifier comparison",
                "",
                f"This slide uses {sample_description}.",
                "The left and right panels come from the same evaluated event population and use identical axes and a shared logarithmic color scale.",
                "On the left, no embedded minimum-bias requirement is applied, and the anomalous near-zero calorimeter-energy band extends across centrality.",
                "For the canonical simulation selection on the right, we require `MinimumBiasInfo::isAuAuMinimumBias()` before downstream analysis and training-tree filling.",
                "That requirement removes the anomalous low-energy population while preserving the main centrality-dependent calorimeter-energy band.",
                "This identifies the feature as a non-minimum-bias embedded underlying-event population; it does not by itself establish the upstream detector or production cause.",
                "",
            ]
        ),
        encoding="utf-8",
    )


def main() -> int:
    args = parse_args()
    input_root = args.input_root.resolve()
    signal_root = args.signal_root.resolve()
    output = args.output.resolve()
    if not input_root.exists():
        raise FileNotFoundError(input_root)
    if args.sample_mode == "combined" and not signal_root.exists():
        raise FileNotFoundError(signal_root)

    with uproot.open(input_root) as root:
        before, centrality_edges, energy_edges = read_density(root, HIST_ALL)
        after, pass_centrality_edges, pass_energy_edges = read_density(root, HIST_PASS)
    if not np.array_equal(centrality_edges, pass_centrality_edges):
        raise ValueError("centrality axes differ between classifier views")
    if not np.array_equal(energy_edges, pass_energy_edges):
        raise ValueError("energy axes differ between classifier views")

    if args.sample_mode == "combined":
        with uproot.open(signal_root) as root:
            signal_before, signal_centrality_edges, signal_energy_edges = read_density(root, HIST_ALL)
            signal_after, signal_pass_centrality_edges, signal_pass_energy_edges = read_density(root, HIST_PASS)
        for label, lhs, rhs in (
            ("all-view centrality", centrality_edges, signal_centrality_edges),
            ("all-view energy", energy_edges, signal_energy_edges),
            ("pass-view centrality", centrality_edges, signal_pass_centrality_edges),
            ("pass-view energy", energy_edges, signal_pass_energy_edges),
        ):
            if not np.array_equal(lhs, rhs):
                raise ValueError(f"signal and inclusive {label} axes differ")
        before = before + signal_before
        after = after + signal_after
        sample_line = "Au+Au embedded training samples: photon+jet 12+20 + inclusive+jet 12+20+30+40"
        sample_description = "the canonically weighted six-sample embedded training population: Photon12+20 and Jet12+20+30+40"
        sample_names = [
            "run28_embeddedPhoton12",
            "run28_embeddedPhoton20",
            "run28_embeddedJet12",
            "run28_embeddedJet20",
            "run28_embeddedJet30",
            "run28_embeddedJet40",
        ]
    else:
        sample_line = "Au+Au embedded training sample: inclusive+jet 12+20+30+40"
        sample_description = "the canonical four-sample embedded inclusive+jet background"
        sample_names = [
            "run28_embeddedJet12",
            "run28_embeddedJet20",
            "run28_embeddedJet30",
            "run28_embeddedJet40",
        ]

    configure_style()
    norm = shared_norm(before, after)
    fig = plt.figure(figsize=(16, 9), dpi=DPI, facecolor="white")

    fig.text(
        0.043,
        0.955,
        "BDT training in Au+Au: pre-training cut in simulation",
        ha="left",
        va="top",
        fontsize=38,
        fontweight="bold",
        color=INK,
    )
    # Blue arrowhead marks the sample line; the requirement line below sits
    # flush with the title and the arrow, with no marker of its own.
    fig.text(
        0.043,
        0.875,
        "▶",
        ha="left",
        va="center",
        fontsize=18,
        color=BLUE_BULLET,
        fontfamily="DejaVu Sans",
    )
    fig.text(
        0.043 + BULLET_INDENT,
        0.875,
        sample_line,
        ha="left",
        va="center",
        fontsize=20.5,
        color=MUTED,
    )
    # Sample line sits at 0.875; drop the requirement line a little further
    # below it for breathing room without moving the line above.
    requirement_y = 0.826
    requirement_label = fig.text(
        0.043,
        requirement_y,
        "Simulation requirement:",
        ha="left",
        va="center",
        fontsize=20.5,
        fontweight="bold",
        color=INK,
    )
    # Place the sentence one space after the bold label's measured right edge
    # so the two read as a single line at any font size.
    fig.canvas.draw()
    label_right = requirement_label.get_window_extent(renderer=fig.canvas.get_renderer()).x1
    space = fig.transFigure.inverted().transform((label_right, 0))[0] + 0.007
    fig.text(
        space,
        requirement_y,
        "Keep events when MinimumBiasInfo::isAuAuMinimumBias() returns true.",
        ha="left",
        va="center",
        fontsize=19.5,
        color=INK,
    )

    # Left edge pulled in from 0.072 so the rotated y-axis title clears the
    # canvas edge; width trimmed to keep the right edge (and the gap to the
    # second panel) where it was.
    left = fig.add_axes([0.088, 0.100, 0.377, 0.650])
    right = fig.add_axes([0.510, 0.100, 0.375, 0.650])
    cax = fig.add_axes([0.902, 0.100, 0.014, 0.650])

    mesh = draw_panel(
        left,
        before,
        centrality_edges,
        energy_edges,
        norm,
        "Before: no classifier requirement",
        # One sPHENIX Internal / beam-energy stamp per slide is enough, and it
        # reads cleanly on the right panel's light corner.  The left panel's
        # white-on-dark stroked variant also rendered a stroke artifact across
        # the bold-italic "sPHENIX" mathtext.
        internal_label=False,
    )
    draw_panel(
        right,
        after,
        centrality_edges,
        energy_edges,
        norm,
        "After: classifier pass",
        right_aligned_internal_label=True,
    )
    left.set_ylabel(r"$E_{calo}^{total}$ [GeV]", fontsize=22, labelpad=2)
    right.tick_params(labelleft=False)
    colorbar = fig.colorbar(mesh, cax=cax)
    colorbar.set_label("Weighted event count / bin", fontsize=17, labelpad=8)
    colorbar.ax.tick_params(labelsize=16)

    left.annotate(
        "low-energy band",
        xy=(28.0, 70.0),
        xytext=(16.0, 440.0),
        fontsize=17,
        fontweight="bold",
        color="#7f1d1d",
        arrowprops={"arrowstyle": "-|>", "color": "#b42318", "lw": 2.5},
        path_effects=[path_effects.withStroke(linewidth=3.5, foreground="white")],
    )

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=DPI, facecolor="white")
    plt.close(fig)

    layout = output.with_suffix(".layout_nodes.json")
    manifest = output.with_suffix(".manifest.json")
    speaker = output.with_name(output.stem + "_speaker_script.md")
    audit = output.with_suffix(".audit.json")
    write_layout(layout, sample_line)
    write_script(speaker, sample_description)
    manifest.write_text(
        json.dumps(
            {
                "schema": "AUAU_EMBEDDED_MINBIAS_BEFORE_AFTER_SLIDE_V1",
                "output_png": str(output),
                "sample_mode": args.sample_mode,
                "input_roots": {
                    "inclusive": {"path": str(input_root), "sha256": sha256(input_root)},
                    "signal": (
                        {"path": str(signal_root), "sha256": sha256(signal_root)}
                        if args.sample_mode == "combined"
                        else None
                    ),
                },
                "samples": sample_names,
                "before_histogram": HIST_ALL,
                "after_histogram": HIST_PASS,
                "classifier_call": "MinimumBiasInfo::isAuAuMinimumBias()",
                "calorimeter_sum": "finite signed TowerInfo energy with get_isGood required",
                "centrality_percent": [0, 80],
                "energy_axis_gev": [0, min(3000.0, float(energy_edges[-1]))],
                "color_scale": "shared logarithmic ROOT-kBird-like scale",
                "weighting": "producer event vertex/centrality weight once, then canonical sample stitch weight once",
                "layout_nodes": str(layout),
                "speaker_script": str(speaker),
                "audit": str(audit),
                "png_dimensions": [WIDTH, HEIGHT],
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    audit_cmd = [
        sys.executable,
        str(REPO / "scripts/slides/common/post_render_slide_audit.py"),
        "--png",
        str(output),
        "--layout-nodes",
        str(layout),
        "--output",
        str(audit),
        "--require-pass",
    ]
    subprocess.run(audit_cmd, cwd=REPO, check=True)
    print(output)
    print(speaker)
    print(manifest)
    print(audit)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
