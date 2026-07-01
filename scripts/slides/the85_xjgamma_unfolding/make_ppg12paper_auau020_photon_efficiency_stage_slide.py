#!/usr/bin/env python3
"""Compose the THE-85 photon-efficiency slide with the PPG12 paper panel."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image, ImageChops


REPO = Path(__file__).resolve().parents[3]
DEFAULT_PPG12_SCREENSHOT = Path(
    "/var/folders/l3/f02nw86n5cn0tpf_zstf0ypr0000gn/T/TemporaryItems/"
    "NSIRD_screencaptureui_iYMBwd/Screenshot 2026-06-29 at 11.02.25\u202fAM.png"
)
DEFAULT_AUAU_PNG = (
    REPO
    / "dataOutput/the85_auau_xjgamma_unfolding_push/efficiency_stage/"
    / "auau020_photon_efficiency_stage_overlay.png"
)
DEFAULT_OUT = (
    REPO
    / "dataOutput/the85_auau_xjgamma_unfolding_push/efficiency_stage/"
    / "pp_vs_auau020_photon_efficiency_stage_1x2_slide.png"
)
DEFAULT_MANIFEST = DEFAULT_OUT.with_name("pp_vs_auau020_photon_efficiency_stage_1x2_manifest.json")

TITLE_COLOR = "#111827"
BODY_COLOR = "#334155"
SLIDE_FONT = "Times New Roman"


def crop_white_border(path: Path, pad: int = 10) -> Image.Image:
    image = Image.open(path).convert("RGBA")
    background = Image.new("RGBA", image.size, (255, 255, 255, 255))
    diff = ImageChops.difference(image, background).convert("L")
    bbox = diff.point(lambda px: 255 if px > 8 else 0).getbbox()
    if not bbox:
        return image
    left = max(0, bbox[0] - pad)
    upper = max(0, bbox[1] - pad)
    right = min(image.width, bbox[2] + pad)
    lower = min(image.height, bbox[3] + pad)
    return image.crop((left, upper, right, lower))


def paste_fit_bottom(canvas: Image.Image, image: Image.Image, box: tuple[int, int, int, int]) -> dict:
    x0, y0, x1, y1 = box
    max_w = x1 - x0
    max_h = y1 - y0
    scale = min(max_w / image.width, max_h / image.height)
    new_size = (int(round(image.width * scale)), int(round(image.height * scale)))
    resized = image.resize(new_size, Image.Resampling.LANCZOS)
    paste_x = x0 + (max_w - new_size[0]) // 2
    paste_y = y1 - new_size[1]
    canvas.alpha_composite(resized, (paste_x, paste_y))
    return {
        "box": list(box),
        "source_size": [image.width, image.height],
        "placed_size": list(new_size),
        "placed_xy": [paste_x, paste_y],
    }


def render_header() -> Image.Image:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": [SLIDE_FONT, "Times", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
        }
    )
    fig = plt.figure(figsize=(16, 9), dpi=160)
    fig.patch.set_facecolor("none")
    fig.text(
        0.055,
        0.955,
        r"Photon efficiency stages, $15<E_T^\gamma<35$ GeV",
        ha="left",
        va="top",
        fontsize=37,
        fontweight="bold",
        fontfamily=SLIDE_FONT,
        color=TITLE_COLOR,
    )
    fig.text(
        0.055,
        0.845,
        "Left is the PPG12 paper reference panel; right is the current Au+Au 0-20% signal-MC stage output.",
        ha="left",
        va="center",
        fontsize=20,
        fontfamily=SLIDE_FONT,
        color=BODY_COLOR,
    )
    fig.text(
        0.255,
        0.780,
        r"$p{+}p$ PPG12 paper reference",
        ha="center",
        va="center",
        fontsize=25,
        fontweight="bold",
        fontfamily=SLIDE_FONT,
        color=TITLE_COLOR,
    )
    fig.text(
        0.745,
        0.780,
        r"Au+Au embedded $\gamma$ MC, 0-20%",
        ha="center",
        va="center",
        fontsize=25,
        fontweight="bold",
        fontfamily=SLIDE_FONT,
        color=TITLE_COLOR,
    )
    fig.canvas.draw()
    image = Image.fromarray(np.asarray(fig.canvas.buffer_rgba()))
    plt.close(fig)
    return image


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ppg12-screenshot", type=Path, default=DEFAULT_PPG12_SCREENSHOT)
    parser.add_argument("--auau-png", type=Path, default=DEFAULT_AUAU_PNG)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    args = parser.parse_args()

    if not args.ppg12_screenshot.exists():
        raise FileNotFoundError(args.ppg12_screenshot)
    if not args.auau_png.exists():
        raise FileNotFoundError(args.auau_png)

    canvas = Image.new("RGBA", (2560, 1440), (255, 255, 255, 255))
    canvas.alpha_composite(render_header(), (0, 0))

    pp_img = crop_white_border(args.ppg12_screenshot, pad=8)
    auau_img = crop_white_border(args.auau_png, pad=16)

    # Keep both plot images bottom-aligned while leaving clear air below the
    # panel labels. The PPG12 screenshot has a near-square crop, so a matching
    # lower box avoids it climbing into the label region.
    placement = {
        "ppg12_paper": paste_fit_bottom(canvas, pp_img, (82, 345, 1258, 1360)),
        "auau_0_20": paste_fit_bottom(canvas, auau_img, (1316, 255, 2468, 1360)),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    canvas.convert("RGB").save(args.out, quality=95)
    manifest = {
        "slide_png": str(args.out),
        "layout": "1x2; PPG12 paper screenshot on left; AuAu current output on right; plots bottom-aligned",
        "ppg12_left": {
            "source": "user-provided screenshot crop from PPG12 paper",
            "source_png": str(args.ppg12_screenshot),
            "role": "reference panel only; not our pp output",
        },
        "auau_right": {
            "source_png": str(args.auau_png),
            "source_manifest": str(args.auau_png.with_name("auau020_photon_efficiency_stage_overlay_manifest.json")),
        },
        "placement": placement,
        "next_step_after_approval": "Run the pp RecoilJets campaign with same-definition truth-binned efficiency-stage objects, then replace the PPG12 reference screenshot with our pp output.",
    }
    args.manifest.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(args.out)
    print(args.manifest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
