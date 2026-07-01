#!/usr/bin/env python3
"""Retouch text on the THE-85 photon-response slide.

The v12 slide image is the source of truth for the plot panels. This helper
leaves the count cards, axes, matrices, and colorbar as the original raster and
redraws the title/subtitle/bullets with clean notation and spacing.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches


REPO = Path(__file__).resolve().parents[3]
DEFAULT_SRC = (
    REPO
    / "dataOutput/the85_auau_xjgamma_unfolding_push/slides/"
    "slide02_photon_1d_response_audience_v12_root_kbird.png"
)
DEFAULT_OUT = (
    REPO
    / "dataOutput/the85_auau_xjgamma_unfolding_push/slides/"
    "slide02_photon_1d_response_audience_v20_binlabels.png"
)
DEFAULT_MANIFEST = DEFAULT_OUT.with_name(DEFAULT_OUT.stem + "_manifest.json")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SRC)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    source = args.source.resolve()
    output = args.output.resolve()
    manifest = args.manifest.resolve()

    if not source.exists():
        raise FileNotFoundError(source)

    img = plt.imread(source)
    height, width = img.shape[:2]
    if (width, height) != (2560, 1440):
        raise RuntimeError(f"Expected 2560x1440 source, got {width}x{height}")

    # Exact-size canvas so the source raster is not rescaled.
    dpi = 160
    fig = plt.figure(figsize=(width / dpi, height / dpi), dpi=dpi)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(img, extent=[0, 1, 0, 1], interpolation="none", aspect="auto")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    # Clear the title/subtitle band and the bullet band. Plot panels/count cards
    # start below this, so the original matrices are not rescaled or redrawn.
    title_clear_x0 = 0.045
    title_clear_y0 = 0.846
    title_clear_w = 0.920
    title_clear_h = 0.112
    ax.add_patch(
        patches.Rectangle(
            (title_clear_x0, title_clear_y0),
            title_clear_w,
            title_clear_h,
            transform=ax.transAxes,
            facecolor="white",
            edgecolor="none",
            zorder=2,
        )
    )

    clear_x0 = 0.045
    clear_y0 = 0.6285
    clear_w = 0.920
    clear_h = 0.2292
    ax.add_patch(
        patches.Rectangle(
            (clear_x0, clear_y0),
            clear_w,
            clear_h,
            transform=ax.transAxes,
            facecolor="white",
            edgecolor="none",
            zorder=2,
        )
    )

    bullet_color = "#007A78"
    text_color = "#111827"
    fig.text(
        0.061,
        0.930,
        r"Photon $p_T^{\gamma}$ response for 15-35 GeV $x_{J\gamma}$ normalization",
        fontsize=32,
        fontfamily="serif",
        fontweight="bold",
        color=text_color,
        ha="left",
        va="center",
        zorder=3,
    )
    fig.text(
        0.061,
        0.869,
        r"Leading-photon $N_{\gamma}(p_T^{\gamma})$ response before normalizing the unfolded $x_{J\gamma}$ spectrum.",
        fontsize=21,
        fontfamily="serif",
        color=text_color,
        ha="left",
        va="center",
        zorder=3,
    )

    bullet_x = 0.060
    text_x = 0.092
    # Even centers across the whitespace between subtitle and plot columns.
    bullet_y = [0.808, 0.736, 0.664]
    bullet_text = [
        r"Maps reconstructed leading-photon $p_T^{\gamma}$ to matched truth $p_T^{\gamma}$ for the photon normalization.",
        r"Panels are restricted to the 15-35 GeV analysis range used in the unfolded $x_{J\gamma}$ comparison.",
        r"This response normalizes the $x_{J\gamma}$ yield per photon; it is separate from inclusive-photon cross-section QA.",
    ]

    for y, text in zip(bullet_y, bullet_text):
        fig.text(
            bullet_x,
            y,
            "\u25b6",
            fontsize=25,
            fontfamily="DejaVu Sans",
            color=bullet_color,
            ha="left",
            va="center",
            zorder=3,
        )
        fig.text(
            text_x,
            y,
            text,
            fontsize=19,
            fontfamily="serif",
            color=text_color,
            ha="left",
            va="center",
            zorder=3,
        )

    # The source raster has the plot panels in the right place. Keep them fixed
    # and redraw only the x-axis titles inside the existing lower safe margin.
    for x in (0.205, 0.505, 0.805):
        fig.text(
            x,
            0.012,
            r"reco $p_T^{\gamma}$ [GeV]",
            ha="center",
            va="bottom",
            fontsize=13,
            fontfamily="serif",
            color=text_color,
            zorder=4,
        )

    def add_px_rect(x0: float, y0: float, w_px: float, h_px: float) -> None:
        ax.add_patch(
            patches.Rectangle(
                (x0 / width, 1.0 - (y0 + h_px) / height),
                w_px / width,
                h_px / height,
                transform=ax.transAxes,
                facecolor="white",
                edgecolor="none",
                zorder=5,
            )
        )

    def fig_text_px(x_px: float, y_px: float, text: str, **kwargs) -> None:
        fig.text(
            x_px / width,
            1.0 - y_px / height,
            text,
            **kwargs,
        )

    # Retouch bin labels without touching the plot panels. The plotted response
    # keeps the 35-40 GeV truth bin as an overflow-support bin.
    # Exact y-axis "OF" text locations in the source raster.
    y_tick_anchors = [220, 970, 1686]
    for x_anchor in y_tick_anchors:
        add_px_rect(x_anchor - 62, 778, 62, 58)
        fig_text_px(
            x_anchor,
            809,
            "OF: 35-40",
            ha="right",
            va="center",
            fontsize=10.5,
            fontfamily="serif",
            color=text_color,
            zorder=6,
        )

    # Remove only the star glyphs from the 15-16* bin labels. This preserves
    # the original ROOT tick-label placement and avoids plot-cell masking.
    star_masks = [
        # y-axis low-edge labels
        (205, 1192, 22, 32),
        (955, 1192, 22, 32),
        (1669, 1192, 22, 32),
        # rotated x-axis low-edge labels
        (247, 1245, 24, 28),
        (997, 1245, 24, 28),
        (1711, 1245, 24, 28),
    ]
    for x0, y0, w_px, h_px in star_masks:
        add_px_rect(x0, y0, w_px, h_px)

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=dpi)
    plt.close(fig)

    manifest.parent.mkdir(parents=True, exist_ok=True)
    manifest.write_text(
        json.dumps(
            {
                "source": str(source),
                "output": str(output),
                "operation": "retouched title, subtitle, bullet band, and x-axis titles; plot panels and count cards preserved from v12 raster",
                "redrawn_xlabel": "reco p_T^gamma [GeV]",
                "bin_label_retouch": {
                    "removed_low_edge_star": "15-16* -> 15-16",
                    "truth_overflow_label": "OF: 35-40",
                    "overflow_support_bin_gev": "35-40",
                },
                "title_clear_band_axes": {
                    "x0": title_clear_x0,
                    "y0": title_clear_y0,
                    "width": title_clear_w,
                    "height": title_clear_h,
                },
                "clear_band_axes": {
                    "x0": clear_x0,
                    "y0": clear_y0,
                    "width": clear_w,
                    "height": clear_h,
                },
                "bullet_y_axes": bullet_y,
                "unchanged_reference": "slide02_photon_1d_response_audience_v12_root_kbird.png",
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    print(output)
    print(manifest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
