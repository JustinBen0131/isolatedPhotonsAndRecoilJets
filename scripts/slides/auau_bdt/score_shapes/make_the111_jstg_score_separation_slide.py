#!/usr/bin/env python3
"""Render the validated combined-correction Au+Au BDT holdout as a JSTG slide.

The scientific panels are cropped from the SHA-pinned THE-111 validation
artifact.  This renderer changes only the slide composition; it does not
recompute scores or metrics.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from PIL import Image


THIS_FILE = Path(__file__).resolve()
REPO = (
    Path(os.environ["THESIS_ANALYSIS_REPO"]).resolve()
    if os.environ.get("THESIS_ANALYSIS_REPO")
    else next((p for p in THIS_FILE.parents if (p / "AGENTS.md").exists() or (p / ".git").exists()))
)
SCRIPTS_DIR = REPO / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.append(str(SCRIPTS_DIR))

from slides.common.slide_defaults import SLIDE_DPI, slide_figsize  # noqa: E402


DEFAULT_SOURCE_DIR = (
    REPO
    / "dataOutput/auauTightBDTValidation"
    / "the111_combined_corrected_shower_ppg12_labels_20260719_1618/paired"
)
DEFAULT_SOURCE_PNG = DEFAULT_SOURCE_DIR / "the111_combined_corrected_shower_ppg12_labels_score_separation_2x3.png"
DEFAULT_VALIDATION = DEFAULT_SOURCE_DIR / "the111_combined_corrected_shower_ppg12_labels_validation.json"
DEFAULT_OUTDIR = REPO / "dataOutput/slides/the45_jstg_20260720/slide21_bdt_score_regeneration/the111_holdout"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-png", type=Path, default=DEFAULT_SOURCE_PNG)
    parser.add_argument("--validation-json", type=Path, default=DEFAULT_VALIDATION)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    return parser.parse_args()


def panel_crops(image: Image.Image) -> list[Image.Image]:
    """Return the three validated PPG12-source-role centrality panels."""
    width, height = image.size
    # Fractions are tied to the stable 2x3 validation canvas and preserve each
    # panel title, metric box, axes, and WP80 line without raster resampling.
    boxes = [
        (0.020, 0.565, 0.350, 0.955),
        (0.365, 0.565, 0.680, 0.955),
        (0.695, 0.565, 1.000, 0.955),
    ]
    return [
        image.crop((int(x0 * width), int(y0 * height), int(x1 * width), int(y1 * height)))
        for x0, y0, x1, y1 in boxes
    ]


def main() -> int:
    args = parse_args()
    for path in (args.source_png, args.validation_json):
        if not path.exists():
            raise FileNotFoundError(path)

    payload = json.loads(args.validation_json.read_text())
    if payload.get("status") != "VALIDATED_SIMULATION_ONLY_NO_PROMOTION":
        raise RuntimeError(f"Unexpected validation status: {payload.get('status')}")
    candidate = payload["native"]["ppg12_variant"]
    bins = candidate["centrality_bins"]
    if len(bins) != 3:
        raise RuntimeError(f"Expected three centrality bins, found {len(bins)}")

    source = Image.open(args.source_png).convert("RGB")
    crops = panel_crops(source)

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "figure.facecolor": "white",
        }
    )
    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI, facecolor="white")

    fig.text(
        0.052,
        0.966,
        "Combined corrected Au+Au photon-ID BDT",
        ha="left",
        va="top",
        fontsize=24.5,
        fontweight="bold",
        color="#172033",
    )
    fig.text(
        0.052,
        0.905,
        r"Held-out score separation; $15<E_T^{\gamma}<35$ GeV, $|\eta|<0.7$; corrected shower inputs and PPG12 source-role labels",
        ha="left",
        va="top",
        fontsize=13.7,
        color="#536171",
    )

    banner = FancyBboxPatch(
        (0.052, 0.800),
        0.896,
        0.046,
        boxstyle="round,pad=0.006,rounding_size=0.008",
        transform=fig.transFigure,
        facecolor="#EAF4EE",
        edgecolor="#8BBE9A",
        linewidth=0.8,
    )
    fig.add_artist(banner)
    fig.text(
        0.070,
        0.822,
        "Higher scores are signal-like; dashed lines mark the independently derived 80% signal-efficiency working point.",
        ha="left",
        va="center",
        fontsize=12.1,
        color="#205B36",
    )

    positions = [(0.035, 0.235, 0.315, 0.525), (0.345, 0.235, 0.315, 0.525), (0.655, 0.235, 0.315, 0.525)]
    for crop, position in zip(crops, positions, strict=True):
        ax = fig.add_axes(position)
        ax.imshow(crop)
        ax.axis("off")

    callout = FancyBboxPatch(
        (0.052, 0.055),
        0.896,
        0.115,
        boxstyle="round,pad=0.012,rounding_size=0.012",
        transform=fig.transFigure,
        facecolor="#F5F8FC",
        edgecolor="#A9B8CA",
        linewidth=0.9,
    )
    fig.add_artist(callout)
    auc = [row["weighted_auc"] for row in bins]
    bkg = [row["wp80"]["background_acceptance"] for row in bins]
    fig.text(
        0.072,
        0.130,
        "Validation result",
        ha="left",
        va="center",
        fontsize=14.3,
        fontweight="bold",
        color="#173B63",
    )
    fig.text(
        0.205,
        0.130,
        "Weighted AUC: " + "  /  ".join(f"{value:.3f}" for value in auc),
        ha="left",
        va="center",
        fontsize=12.6,
        color="#172033",
    )
    fig.text(
        0.515,
        0.130,
        "WP80 background acceptance: " + "  /  ".join(f"{100.0 * value:.1f}%" for value in bkg),
        ha="left",
        va="center",
        fontsize=12.6,
        color="#172033",
    )
    fig.text(
        0.072,
        0.085,
        "Simulation holdout validation only. The matched data, photon-embedding, and inclusive-embedding application is the next production gate.",
        ha="left",
        va="center",
        fontsize=11.0,
        color="#6B4A1D",
    )

    args.outdir.mkdir(parents=True, exist_ok=True)
    stem = "the111_combined_corrected_auau_bdt_score_separation_jstg"
    png = args.outdir / f"{stem}.png"
    manifest = args.outdir / f"{stem}_manifest.json"
    speaker = args.outdir / f"{stem}_speaker_script.md"
    fig.savefig(png, dpi=SLIDE_DPI)
    plt.close(fig)

    manifest_payload = {
        "schema": "THE111_JSTG_SCORE_SEPARATION_SLIDE_V1",
        "status": "JSTG_PNG_CANDIDATE_SIMULATION_ONLY",
        "source_png": str(args.source_png),
        "source_png_sha256": sha256(args.source_png),
        "validation_json": str(args.validation_json),
        "validation_json_sha256": sha256(args.validation_json),
        "candidate_model_tmva_sha256": payload["inputs"]["ppg12_variant"]["sha256"]["tmva"],
        "validation_status": payload["status"],
        "centrality_metrics": bins,
        "rendering": {
            "scientific_content": "lossless crops of the validated PPG12-source-role row",
            "recomputed_scores": False,
            "canvas_pixels": [2560, 1440],
        },
        "output_png": str(png),
        "output_png_sha256": sha256(png),
        "speaker_script": str(speaker),
    }
    manifest.write_text(json.dumps(manifest_payload, indent=2, sort_keys=True) + "\n")
    speaker.write_text(
        "\n".join(
            [
                "# Speaker script",
                "",
                "This is the combined corrected Au+Au photon-identification BDT candidate evaluated on its held-out simulation rows.",
                "The shower inputs use the corrected common TowerInfo construction, and the labels use the PPG12 source-role contract: prompt candidates from photon sources are signal and non-prompt candidates from jet sources are background.",
                "Green is signal and orange is background. The dashed line in each centrality panel is the threshold retaining 80 percent of the weighted signal.",
                f"The weighted AUC values are {auc[0]:.3f}, {auc[1]:.3f}, and {auc[2]:.3f} from central to peripheral events.",
                f"At WP80 the weighted background acceptance falls from {100*bkg[0]:.1f} percent in 0 to 20 percent centrality to {100*bkg[2]:.1f} percent in 50 to 80 percent centrality.",
                "This is simulation validation, not data closure. The matched data and embedded production pass is the next gate.",
                "",
            ]
        )
    )
    print(png)
    print(manifest)
    print(speaker)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
