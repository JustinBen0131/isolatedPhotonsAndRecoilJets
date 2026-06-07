#!/usr/bin/env python3
"""Build data-first THE-32 low-calo views without the threshold overlay."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFilter, ImageFont


W = 2560
H = 1440
OUTDIR = Path("dataOutput/auauTightBDTValidation/THE32_lowCaloDiagnosticClosure_20260603")
SOURCE_PNG = OUTDIR / "the32_low_calo_pathology_slide_v1.png"
SOURCE_JSON = OUTDIR / "the32_low_calo_pathology_slide_v1.json"

# Inner data rectangle from the validated main panel. This intentionally excludes
# the original threshold legend, axis labels, colorbar, and side panels.
DATA_BOX = (150, 315, 1655, 980)

INK = (24, 34, 48)
MUTED = (71, 84, 103)
GRID = (218, 224, 232)
BLUE = (33, 112, 181)
RED = (184, 35, 24)
YELLOW = (255, 247, 214)
YELLOW_BORDER = (227, 179, 65)
BORDER = (198, 208, 222)
GRAY = (247, 249, 252)


def font(size: int, *, bold: bool = False, italic: bool = False) -> ImageFont.FreeTypeFont:
    if bold and italic:
        names = ["Times New Roman Bold Italic.ttf", "Arial Bold Italic.ttf"]
    elif bold:
        names = ["Times New Roman Bold.ttf", "Arial Bold.ttf"]
    elif italic:
        names = ["Times New Roman Italic.ttf", "Arial Italic.ttf"]
    else:
        names = ["Times New Roman.ttf", "Arial.ttf"]
    for name in names:
        path = Path("/System/Library/Fonts/Supplemental") / name
        if path.exists():
            return ImageFont.truetype(str(path), size=size)
    return ImageFont.load_default(size=size)


F_TITLE = font(58, bold=True)
F_SUB = font(34)
F_AXIS = font(31)
F_TICK = font(26)
F_LABEL = font(32, bold=True)
F_BODY = font(31)
F_BODY_BOLD = font(34, bold=True)
F_SMALL = font(25)


def draw_wrapped(draw: ImageDraw.ImageDraw, xy, text: str, *, max_width: int, text_font, fill=INK, spacing=8) -> int:
    words = text.split()
    lines: list[str] = []
    current = ""
    for word in words:
        trial = word if not current else f"{current} {word}"
        if draw.textbbox((0, 0), trial, font=text_font)[2] <= max_width:
            current = trial
        else:
            if current:
                lines.append(current)
            current = word
    if current:
        lines.append(current)
    x, y = xy
    for line in lines:
        draw.text((x, y), line, font=text_font, fill=fill)
        y += text_font.size + spacing
    return y


def box(draw: ImageDraw.ImageDraw, xy, *, fill, outline=BORDER, width=3, radius=22) -> None:
    draw.rounded_rectangle(xy, radius=radius, fill=fill, outline=outline, width=width)


def source_layers() -> tuple[Image.Image, Image.Image, dict]:
    source = Image.open(SOURCE_PNG).convert("RGB")
    crop = source.crop(DATA_BOX)
    arr = np.asarray(crop).astype("int16")
    r, g, b = arr[:, :, 0], arr[:, :, 1], arr[:, :, 2]

    # Keep only the plotted data colors, not the black threshold/axes/legend.
    red_mask = (r > 120) & (r > g + 30) & (r > b + 30)
    blue_mask = (b > 100) & (b > r + 12) & (g > r - 8)
    blue_mask &= ~red_mask

    # Remove non-data furniture that lived inside the original axes: the
    # sPHENIX label, the original legend box, and antialiased border remnants.
    h, w = blue_mask.shape
    furniture = np.zeros((h, w), dtype=bool)
    furniture[:115, :760] = True
    furniture[520:650, :900] = True
    furniture[:25, :] = True
    furniture[-35:, :] = True
    furniture[:, :18] = True
    furniture[:, -18:] = True
    red_mask &= ~furniture
    blue_mask &= ~furniture

    blue_rgba = np.zeros((*arr.shape[:2], 4), dtype="uint8")
    blue_rgba[blue_mask, :3] = arr[blue_mask].clip(0, 255).astype("uint8")
    blue_rgba[blue_mask, 3] = 205

    red_rgba = np.zeros((*arr.shape[:2], 4), dtype="uint8")
    red_rgba[red_mask, :3] = RED
    red_rgba[red_mask, 3] = 210

    red_layer = Image.fromarray(red_rgba).filter(ImageFilter.MaxFilter(3))
    blue_layer = Image.fromarray(blue_rgba)

    pathology = json.loads(SOURCE_JSON.read_text())
    rows = pathology["event_counts"]
    total = sum(int(row["event_total"]) for row in rows)
    rejected = sum(int(row["event_rejected"]) for row in rows)
    return blue_layer, red_layer, {"event_total": total, "event_rejected": rejected}


def title(draw: ImageDraw.ImageDraw, main: str, sub: str) -> None:
    draw.text((100, 65), main, font=F_TITLE, fill=INK)
    draw_wrapped(draw, (104, 150), sub, max_width=2200, text_font=F_SUB, fill=MUTED)


def draw_axes(
    canvas: Image.Image,
    draw: ImageDraw.ImageDraw,
    *,
    plot_box: tuple[int, int, int, int],
    x_min: float,
    x_max: float,
    y_min: float = 2.65,
    y_max: float = 3.50,
) -> None:
    x0, y0, x1, y1 = plot_box
    draw.rectangle((x0, y0, x1, y1), fill="white", outline=INK, width=3)
    for xtick in np.arange(x_min, x_max + 1e-6, 10):
        px = int(x0 + (xtick - x_min) / (x_max - x_min) * (x1 - x0))
        draw.line((px, y0, px, y1), fill=GRID, width=2)
        draw.text((px - 16, y1 + 18), f"{int(xtick)}", font=F_TICK, fill=INK)
    for ytick in [2.8, 3.0, 3.2, 3.4]:
        py = int(y1 - (ytick - y_min) / (y_max - y_min) * (y1 - y0))
        draw.line((x0, py, x1, py), fill=GRID, width=2)
        draw.text((x0 - 72, py - 16), f"{ytick:.1f}", font=F_TICK, fill=INK)
    draw.text((x0 + (x1 - x0) // 2 - 155, y1 + 64), "Centrality percentile", font=F_AXIS, fill=INK)
    ylabel = "log10(CEMC + IHCal + OHCal + 1)"
    label = Image.new("RGBA", (700, 60), (255, 255, 255, 0))
    ld = ImageDraw.Draw(label)
    ld.text((0, 0), ylabel, font=F_AXIS, fill=INK)
    rotated = label.rotate(90, expand=True)
    canvas.paste(rotated, (x0 - 135, y0 + (y1 - y0) // 2 - rotated.height // 2), rotated)


def compose_data_layer(
    canvas: Image.Image,
    blue_layer: Image.Image,
    red_layer: Image.Image,
    *,
    plot_box: tuple[int, int, int, int],
    x_min: float,
    x_max: float,
    red_boost: bool = True,
) -> None:
    left, top, right, bottom = DATA_BOX
    src_w = right - left
    x0 = int((x_min / 80.0) * src_w)
    x1 = int((x_max / 80.0) * src_w)
    blue = blue_layer.crop((x0, 0, x1, blue_layer.height))
    red = red_layer.crop((x0, 0, x1, red_layer.height))
    w = plot_box[2] - plot_box[0]
    h = plot_box[3] - plot_box[1]
    blue = blue.resize((w, h), Image.Resampling.LANCZOS)
    red = red.resize((w, h), Image.Resampling.LANCZOS)
    if red_boost:
        red = red.filter(ImageFilter.MaxFilter(3))
    canvas.alpha_composite(blue, (plot_box[0], plot_box[1]))
    canvas.alpha_composite(red, (plot_box[0], plot_box[1]))


def add_legend(draw: ImageDraw.ImageDraw, x: int, y: int) -> None:
    box(draw, (x, y, x + 690, y + 112), fill="white", outline=BORDER, width=2, radius=14)
    draw.rectangle((x + 28, y + 28, x + 70, y + 62), fill=BLUE)
    draw.text((x + 88, y + 20), "normal event band", font=F_SMALL, fill=INK)
    draw.ellipse((x + 375, y + 26, x + 415, y + 66), fill=RED)
    draw.text((x + 435, y + 20), "events removed", font=F_SMALL, fill=INK)


def variant_full_data(blue_layer: Image.Image, red_layer: Image.Image, counts: dict) -> Path:
    canvas = Image.new("RGBA", (W, H), "white")
    draw = ImageDraw.Draw(canvas)
    title(
        draw,
        "Data view first: the low-calo events separate below the normal band",
        "No threshold curve is drawn here. The red points are the events removed by the cut; the blue density is the rest of the sample.",
    )
    plot_box = (290, 250, 2320, 1035)
    box(draw, (plot_box[0] - 35, plot_box[1] - 35, plot_box[2] + 35, plot_box[3] + 105), fill=GRAY, outline=BORDER)
    draw_axes(canvas, draw, plot_box=plot_box, x_min=0.0, x_max=80.0)
    compose_data_layer(canvas, blue_layer, red_layer, plot_box=plot_box, x_min=0.0, x_max=80.0)
    add_legend(draw, 1520, 305)
    box(draw, (180, 1190, 2380, 1328), fill=YELLOW, outline=YELLOW_BORDER)
    draw.text((230, 1228), "What this shows:", font=F_BODY_BOLD, fill=INK)
    draw_wrapped(
        draw,
        (570, 1232),
        f"The red events form a separate low-calorimeter-energy tail underneath the main blue population. That tail is {counts['event_rejected']:,}/{counts['event_total']:,} events in the diagnostic sample.",
        max_width=1780,
        text_font=F_BODY,
    )
    out = OUTDIR / "the32_low_calo_data_view_no_threshold_v1.png"
    canvas.convert("RGB").save(out)
    return out


def variant_central_zoom(blue_layer: Image.Image, red_layer: Image.Image, counts: dict) -> Path:
    canvas = Image.new("RGBA", (W, H), "white")
    draw = ImageDraw.Draw(canvas)
    title(
        draw,
        "Zoom on central events: the removed population is visibly below the band",
        "This is the clearest raw-data view of the separation: blue is the normal event band, red is the low-calo tail.",
    )
    plot_box = (310, 255, 2310, 1040)
    box(draw, (plot_box[0] - 35, plot_box[1] - 35, plot_box[2] + 35, plot_box[3] + 105), fill=GRAY, outline=BORDER)
    draw_axes(canvas, draw, plot_box=plot_box, x_min=0.0, x_max=35.0)
    compose_data_layer(canvas, blue_layer, red_layer, plot_box=plot_box, x_min=0.0, x_max=35.0, red_boost=True)
    add_legend(draw, 1450, 315)
    box(draw, (190, 1190, 2370, 1328), fill=YELLOW, outline=YELLOW_BORDER)
    draw.text((240, 1228), "Read it this way:", font=F_BODY_BOLD, fill=INK)
    draw_wrapped(
        draw,
        (595, 1232),
        "For the same centrality range, the red points sit at much lower total calo energy than the blue band. Those are not a new physics class; they are event-quality outliers.",
        max_width=1740,
        text_font=F_BODY,
    )
    out = OUTDIR / "the32_low_calo_data_view_central_zoom_v1.png"
    canvas.convert("RGB").save(out)
    return out


def variant_removed_only(blue_layer: Image.Image, red_layer: Image.Image, counts: dict) -> Path:
    canvas = Image.new("RGBA", (W, H), "white")
    draw = ImageDraw.Draw(canvas)
    title(
        draw,
        "Removed events only: the cut targets the low-calo tail",
        "This view suppresses the blue density so the removed population is easy to inspect.",
    )
    plot_box = (300, 260, 2315, 1040)
    box(draw, (plot_box[0] - 35, plot_box[1] - 35, plot_box[2] + 35, plot_box[3] + 105), fill=GRAY, outline=BORDER)
    draw_axes(canvas, draw, plot_box=plot_box, x_min=0.0, x_max=80.0)
    # Put a faint blue context layer underneath, then emphasize red.
    faint_blue = blue_layer.copy()
    alpha = faint_blue.getchannel("A")
    alpha = alpha.point(lambda v: int(v * 0.22))
    faint_blue.putalpha(alpha)
    compose_data_layer(canvas, faint_blue, red_layer, plot_box=plot_box, x_min=0.0, x_max=80.0, red_boost=True)
    add_legend(draw, 1510, 315)
    box(draw, (185, 1190, 2375, 1328), fill=YELLOW, outline=YELLOW_BORDER)
    draw.text((235, 1228), "Use this to see the cut:", font=F_BODY_BOLD, fill=INK)
    draw_wrapped(
        draw,
        (690, 1232),
        "The red points cluster at the bottom of the total-calo axis across centralities. This is the population removed before training.",
        max_width=1650,
        text_font=F_BODY,
    )
    out = OUTDIR / "the32_low_calo_data_view_removed_only_v1.png"
    canvas.convert("RGB").save(out)
    return out


def main() -> int:
    blue_layer, red_layer, counts = source_layers()
    outputs = [
        variant_full_data(blue_layer, red_layer, counts),
        variant_central_zoom(blue_layer, red_layer, counts),
        variant_removed_only(blue_layer, red_layer, counts),
    ]
    note = OUTDIR / "the32_low_calo_data_views_v1.md"
    note.write_text(
        "\n".join(
            [
                "# THE-32 low-calo data-first views",
                "",
                "These candidates use the validated original diagnostic plot as the data source, but remove the threshold line, legend, and side panels.",
                "The point is to first show the visual separation in the data: red events sit below the blue normal event band in total calorimeter energy at fixed centrality.",
                "",
            ]
        )
    )
    manifest = OUTDIR / "the32_low_calo_data_views_v1.json"
    manifest.write_text(
        json.dumps(
            {
                "schema": "THE32_LOW_CALO_DATA_VIEWS_V1",
                "source_png": str(SOURCE_PNG),
                "source_json": str(SOURCE_JSON),
                "source_data_box_px": list(DATA_BOX),
                "method": "Extract real blue/red plotted data pixels from the validated pathology plot, omit threshold curve and side panels, redraw clean axes.",
                "event_total": counts["event_total"],
                "event_rejected": counts["event_rejected"],
                "event_rejected_fraction": counts["event_rejected"] / counts["event_total"],
                "outputs": [str(path) for path in outputs],
                "speaker_note": str(note),
                "google_slides_mutated": False,
            },
            indent=2,
            sort_keys=True,
        )
    )
    for path in outputs:
        print(path)
    print(note)
    print(manifest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
