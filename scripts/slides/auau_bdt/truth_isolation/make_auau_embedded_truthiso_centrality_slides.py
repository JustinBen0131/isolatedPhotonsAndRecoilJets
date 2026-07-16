#!/usr/bin/env python3
"""Build three audience-facing AuAu truth-isolation slide PNGs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from PIL import Image, ImageOps


ROOT = Path(__file__).resolve().parents[4]
DEFAULT_INPUT = ROOT / "dataOutput/auau/embedded_truthiso_the88a_20260709"
DEFAULT_OUTPUT = DEFAULT_INPUT / "slides"

WIDTH = 2560
HEIGHT = 1440
DPI = 200
INK = "#111827"
MUTED = "#475569"


SLIDES = (
    {
        "stem": "slide_01_direct_photon_truth_isolation_by_centrality",
        "title": "Direct-photon truth-isolation efficiency by centrality",
        "subtitle": (
            "Au+Au embedded photon+jet 12 + 20 simulation\n"
            "Cumulative fraction passing each truth-level $R=0.3$ isolation cutoff."
        ),
        "note": (
            "Panels compare 0-20%, 20-50%, and 50-80%; curves are split by truth-photon $p_T$."
        ),
        "plot": "auau_embedded_truthiso_direct_centrality_summary.png",
    },
    {
        "stem": "slide_02_fragmentation_photon_truth_isolation_by_centrality",
        "title": "Fragmentation-photon truth-isolation by centrality",
        "subtitle": (
            "Au+Au embedded photon+jet 12 + 20 simulation\n"
            "Cumulative fraction passing each truth-level $R=0.3$ isolation cutoff."
        ),
        "note": (
            "Panels compare 0-20%, 20-50%, and 50-80%; curves are split by truth-photon $p_T$."
        ),
        "plot": "auau_embedded_truthiso_fragmentation_centrality_summary.png",
    },
    {
        "stem": "slide_03_truth_photon_pt_after_isolation_by_centrality",
        "title": "Truth-photon $p_T$ spectra by centrality",
        "subtitle": (
            "Au+Au embedded photon+jet 12 + 20 simulation\n"
            "15-35 GeV after truth-level $E_T^{iso}<4$ GeV in $R=0.3$."
        ),
        "note": (
            "Upper panels show weighted total, direct, and fragmentation spectra; lower panels show direct / total."
        ),
        "plot": "auau_embedded_truthpt_spectrum_centrality_summary.png",
    },
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def contained_axes(fig: plt.Figure, image_path: Path) -> tuple[plt.Axes, list[float]]:
    image = Image.open(image_path)
    image_aspect = image.width / image.height
    box_left, box_bottom, box_width, box_height = 0.035, 0.145, 0.93, 0.665
    box_aspect = (box_width * WIDTH) / (box_height * HEIGHT)

    if image_aspect >= box_aspect:
        width = box_width
        height = width * WIDTH / image_aspect / HEIGHT
        left = box_left
        bottom = box_bottom + 0.5 * (box_height - height)
    else:
        height = box_height
        width = height * HEIGHT * image_aspect / WIDTH
        left = box_left + 0.5 * (box_width - width)
        bottom = box_bottom

    plot_top_px = (1.0 - bottom - height) * HEIGHT
    minimum_plot_top_px = 318.0
    if plot_top_px < minimum_plot_top_px:
        bottom -= (minimum_plot_top_px - plot_top_px) / HEIGHT

    axes = fig.add_axes([left, bottom, width, height])
    axes.imshow(mpimg.imread(image_path))
    axes.axis("off")
    bbox = [left * WIDTH, (1.0 - bottom - height) * HEIGHT,
            (left + width) * WIDTH, (1.0 - bottom) * HEIGHT]
    return axes, bbox


def layout_payload(
    spec: dict[str, str],
    plot_bbox: list[float],
    subtitle_bbox: list[float],
) -> dict:
    return {
        "canvas": [0, 0, WIDTH, HEIGHT],
        "title_axis_x": 128,
        "title_axis_tolerance_px": 12,
        "minimum_audience_font_px": 48,
        "minimum_title_font_px": 96,
        "nodes": [
            {
                "kind": "text",
                "name": "slide title",
                "role": "title",
                "text": spec["title"],
                "bbox": [128, 58, 2440, 150],
                "font_px": 100,
                "title_anchor": True,
            },
            {
                "kind": "text",
                "name": "sample and calculation description",
                "role": "audience",
                "text": spec["subtitle"],
                "bbox": subtitle_bbox,
                "font_px": 64,
                "title_axis_align": "left",
            },
            {
                "kind": "image",
                "name": "centrality comparison plot",
                "role": "evidence",
                "bbox": plot_bbox,
            },
            {
                "kind": "text",
                "name": "plot reading note",
                "role": "audience",
                "text": spec["note"],
                "bbox": [132, 1250, 2428, 1380],
                "font_px": 48,
                "title_axis_align": "left",
            },
        ],
    }


def render_slide(spec: dict[str, str], input_dir: Path, output_dir: Path) -> dict:
    plot_path = input_dir / spec["plot"]
    if not plot_path.is_file():
        raise FileNotFoundError(plot_path)

    fig = plt.figure(figsize=(WIDTH / DPI, HEIGHT / DPI), dpi=DPI)
    fig.patch.set_facecolor("white")
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
        }
    )

    fig.text(
        0.05, 0.958, spec["title"], ha="left", va="top",
        fontsize=36, fontweight="bold", color=INK,
    )
    _, plot_bbox = contained_axes(fig, plot_path)
    subtitle_height = 112.0
    subtitle_top = 0.5 * (150.0 + plot_bbox[1] - subtitle_height)
    subtitle_bbox = [132.0, subtitle_top, 2428.0, subtitle_top + subtitle_height]
    fig.text(
        0.052, 1.0 - subtitle_top / HEIGHT, spec["subtitle"],
        ha="left", va="top", fontsize=23, color=INK, linespacing=1.08,
    )
    fig.text(
        0.052, 0.050, spec["note"], ha="left", va="bottom",
        fontsize=17, color=MUTED, linespacing=1.12,
    )

    png_path = output_dir / f"{spec['stem']}.png"
    layout_path = output_dir / f"{spec['stem']}.layout.json"
    fig.savefig(png_path, dpi=DPI, facecolor="white")
    plt.close(fig)
    layout_path.write_text(
        json.dumps(layout_payload(spec, plot_bbox, subtitle_bbox), indent=2) + "\n",
        encoding="utf-8",
    )
    return {
        "png": str(png_path),
        "layout": str(layout_path),
        "source_plot": str(plot_path),
        "title": spec["title"],
    }


def write_montage(png_paths: list[Path], output_path: Path) -> None:
    thumb_width = 960
    thumb_height = 540
    montage = Image.new("RGB", (thumb_width, thumb_height * len(png_paths)), "white")
    for index, path in enumerate(png_paths):
        image = Image.open(path).convert("RGB")
        thumb = ImageOps.contain(image, (thumb_width, thumb_height))
        montage.paste(thumb, ((thumb_width - thumb.width) // 2, index * thumb_height))
    montage.save(output_path)


def main() -> int:
    args = parse_args()
    input_dir = args.input_dir.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    rendered = [render_slide(spec, input_dir, output_dir) for spec in SLIDES]
    montage_path = output_dir / "auau_truthiso_centrality_slides_montage.png"
    write_montage([Path(item["png"]) for item in rendered], montage_path)

    notes_path = output_dir / "speaker_notes.txt"
    notes_path.write_text(
        "Slide 1\n"
        "Direct photons pass the truth-isolation requirement with very high efficiency. "
        "The cutoff dependence and the three centrality selections are nearly identical.\n\n"
        "Slide 2\n"
        "Fragmentation photons have a broader truth-isolation distribution. The passing "
        "fraction rises with the cutoff and remains consistent across centrality.\n\n"
        "Slide 3\n"
        "After the 4 GeV truth-isolation requirement, direct photons dominate the weighted "
        "spectrum. The lower panels show the direct share of the selected spectrum.\n",
        encoding="utf-8",
    )

    manifest = {
        "schema_version": 1,
        "audience_facing": True,
        "visible_task_ids": False,
        "centrality_percent": [[0, 20], [20, 50], [50, 80]],
        "slides": rendered,
        "montage": str(montage_path),
        "speaker_notes": str(notes_path),
        "source_values": str(input_dir / "auau_embedded_truthiso_diagnostics_values.csv"),
        "source_manifest": str(input_dir / "auau_embedded_truthiso_diagnostics_manifest.json"),
    }
    manifest_path = output_dir / "auau_truthiso_centrality_slides_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
