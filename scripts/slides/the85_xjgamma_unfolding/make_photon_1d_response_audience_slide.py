#!/usr/bin/env python3
"""Build the THE-85 photon pT response slide from the matrix manifest."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.patches import FancyBboxPatch


REPO = Path(__file__).resolve().parents[3]
DEFAULT_MANIFEST = (
    REPO
    / "dataOutput/the85_auau_xjgamma_unfolding_push/slides/"
    "slide02_photon_1d_response_audience_v9_no_uf_manifest.json"
)
DEFAULT_OUTPUT = (
    REPO
    / "dataOutput/the85_auau_xjgamma_unfolding_push/slides/"
    "slide02_photon_1d_response_audience_v29_generated_binlabels.png"
)


TEXT = "#111827"
TEAL = "#007A78"
CARD_EDGE = "#BCD5EB"
CARD_FILL = "#F4F8FC"


def kbird_colormap() -> LinearSegmentedColormap:
    """Approximate ROOT kBird with enough contrast for log response matrices."""

    colors = [
        (0.150, 0.130, 0.500),
        (0.050, 0.330, 0.800),
        (0.000, 0.620, 0.790),
        (0.160, 0.710, 0.610),
        (0.560, 0.730, 0.430),
        (0.930, 0.740, 0.250),
        (1.000, 0.910, 0.000),
    ]
    cmap = LinearSegmentedColormap.from_list("kBirdApprox", colors, N=256)
    cmap.set_bad("white")
    return cmap


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def clean_label(label: str) -> str:
    if label == "15-16*":
        return "14-16"
    if label == "OF":
        return "OF:\n35-40"
    return label


def short_count(value: float) -> str:
    return f"{value / 1_000_000:.1f}M"


def add_card(fig: plt.Figure, x: float, y: float, w: float, h: float, text: str) -> None:
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.010,rounding_size=0.006",
        transform=fig.transFigure,
        facecolor=CARD_FILL,
        edgecolor=CARD_EDGE,
        linewidth=1.4,
        zorder=2,
    )
    fig.patches.append(patch)
    fig.text(
        x + w / 2,
        y + h / 2,
        text,
        ha="center",
        va="center",
        fontsize=19,
        fontweight="bold",
        fontfamily="serif",
        color=TEXT,
        zorder=3,
    )


def main() -> int:
    args = parse_args()
    manifest_path = args.manifest.resolve()
    output = args.output.resolve()
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

    panels = manifest["inputs"]
    cmap = kbird_colormap()
    norm = LogNorm(vmin=1, vmax=1e7)

    fig = plt.figure(figsize=(16, 9), dpi=160, facecolor="white")

    fig.text(
        0.061,
        0.930,
        r"Photon $p_T^{\gamma}$ response for 15-35 GeV $x_{J\gamma}$ normalization",
        ha="left",
        va="center",
        fontsize=32,
        fontfamily="serif",
        fontweight="bold",
        color=TEXT,
    )
    fig.text(
        0.061,
        0.869,
        r"Leading-photon $N_{\gamma}(p_T^{\gamma})$ response before normalizing the unfolded $x_{J\gamma}$ spectrum.",
        ha="left",
        va="center",
        fontsize=21,
        fontfamily="serif",
        color=TEXT,
    )

    bullet_y = [0.812, 0.744, 0.676, 0.608]
    bullet_text = [
        r"Maps reconstructed leading-photon $p_T^{\gamma}$ to matched truth $p_T^{\gamma}$ for the photon normalization.",
        r"Panels are restricted to the 15-35 GeV analysis range used in the unfolded $x_{J\gamma}$ comparison.",
        r"This response normalizes the $x_{J\gamma}$ yield per photon; it is separate from inclusive-photon cross-section QA.",
        r"Current first bin is 14-16 GeV; future response passes should start at 15 GeV.",
    ]
    for y, text in zip(bullet_y, bullet_text):
        fig.text(0.060, y, "\u25b6", ha="left", va="center", fontsize=25, fontfamily="DejaVu Sans", color=TEAL)
        fig.text(0.092, y, text, ha="left", va="center", fontsize=19, fontfamily="serif", color=TEXT)

    lefts = [0.090, 0.392, 0.694]
    width = 0.215
    bottom = 0.104
    height = 0.300
    card_y = 0.475
    card_h = 0.050
    card_w = 0.200
    image = None

    for idx, (left, panel) in enumerate(zip(lefts, panels)):
        add_card(fig, left + (width - card_w) / 2, card_y, card_w, card_h, f"Counts: {short_count(panel['physics_sum_cells_15to35'])}")
        ax = fig.add_axes([left, bottom, width, height])
        matrix = np.array(panel["matrix_truth_x_reco_y_displayed"], dtype=float)
        masked = np.ma.masked_less_equal(matrix, 0.0)
        image = ax.imshow(masked, origin="lower", aspect="auto", cmap=cmap, norm=norm, interpolation="nearest")

        x_labels = [clean_label(x) for x in panel["y_category_labels"]]
        y_labels = [clean_label(y) for y in panel["x_category_labels"]]
        ax.set_xticks(np.arange(len(x_labels)))
        ax.set_xticklabels(x_labels, rotation=45, ha="right", rotation_mode="anchor", fontsize=13, fontfamily="serif")
        ax.set_yticks(np.arange(len(y_labels)))
        ax.set_yticklabels(y_labels, fontsize=12.2, fontfamily="serif", linespacing=0.88)
        ax.tick_params(axis="both", direction="in", top=False, right=False, width=1.1, length=5, pad=2)
        ax.set_xlabel(r"reco $p_T^{\gamma}$ [GeV]", fontsize=12.2, fontfamily="serif", labelpad=0)
        if idx == 0:
            ax.set_ylabel(r"truth $p_T^{\gamma}$ [GeV]", fontsize=18, fontfamily="serif", labelpad=10)
        else:
            ax.set_ylabel("")
        ax.set_title(panel["label"], fontsize=24, fontfamily="serif", fontweight="bold", color=TEXT, pad=6)
        for spine in ax.spines.values():
            spine.set_linewidth(1.2)

    cax = fig.add_axes([0.928, bottom, 0.014, height])
    cb = fig.colorbar(image, cax=cax)
    cb.set_label("matched photons", fontsize=15, fontfamily="serif", rotation=90, labelpad=9)
    cb.ax.tick_params(labelsize=12.5, width=1.0, length=4)
    for label in cb.ax.get_yticklabels():
        label.set_fontfamily("serif")

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=160, bbox_inches=None)
    plt.close(fig)

    out_manifest = output.with_name(output.stem + "_manifest.json")
    out_manifest.write_text(
        json.dumps(
            {
                "source_manifest": str(manifest_path),
                "output_png": str(output),
                "render_contract": "Generated from matrix arrays; no raster label patching.",
                "label_changes": {
                    "low_edge_bin": "14-16",
                    "truth_overflow_label": "OF: 35-40",
                    "overflow_support_bin_gev": "35-40",
                },
                "boundary_note": (
                    "Stored response bin is 14-16 GeV while the nominal BDT/training analysis "
                    "window is 15-35 GeV; future response production should book a 15 GeV lower edge."
                ),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    print(output)
    print(out_manifest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
