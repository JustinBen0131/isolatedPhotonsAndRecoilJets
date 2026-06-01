#!/usr/bin/env python3
"""Generate standalone photon-detector motif PNGs for the HP2026 title slide.

The assets are graphic-only transparent PNGs intended to sit below the
``on behalf of the sPHENIX Collaboration`` line on the title slide.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from PIL import Image, ImageDraw, ImageFilter, ImageFont


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_WORKSPACE = ROOT / "outputs/manual-20260601-hp2026-title/presentations/hp2026-title-slide"
W, H = 1320, 320
S = 3

NAVY = (19, 41, 75)
INK = (18, 22, 28)
MUTED = (76, 84, 98)
SOFT = (229, 234, 241)
BLUE = (30, 143, 214)
TEAL = (16, 92, 120)
GOLD = (246, 185, 35)
ORANGE = (255, 95, 5)
MAROON = (179, 14, 57)
WHITE = (255, 255, 255)

FONT_DIR = Path("/System/Library/Fonts/Supplemental")
TIMES_BOLD = FONT_DIR / "Times New Roman Bold.ttf"


def srgb(color: tuple[int, int, int], alpha: int = 255) -> tuple[int, int, int, int]:
    return (*color, alpha)


def p(x: float, y: float) -> tuple[int, int]:
    return int(round(x * S)), int(round(y * S))


def b(x0: float, y0: float, x1: float, y1: float) -> tuple[int, int, int, int]:
    return (*p(x0, y0), *p(x1, y1))


def scaled_points(points: list[tuple[float, float]]) -> list[tuple[int, int]]:
    return [p(x, y) for x, y in points]


def canvas() -> Image.Image:
    return Image.new("RGBA", (W * S, H * S), (255, 255, 255, 0))


def downsample(img: Image.Image) -> Image.Image:
    return img.resize((W, H), Image.Resampling.LANCZOS)


def cubic(
    p0: tuple[float, float],
    p1: tuple[float, float],
    p2: tuple[float, float],
    p3: tuple[float, float],
    n: int = 96,
) -> list[tuple[float, float]]:
    pts = []
    for i in range(n):
        t = i / (n - 1)
        u = 1 - t
        x = u**3 * p0[0] + 3 * u**2 * t * p1[0] + 3 * u * t**2 * p2[0] + t**3 * p3[0]
        y = u**3 * p0[1] + 3 * u**2 * t * p1[1] + 3 * u * t**2 * p2[1] + t**3 * p3[1]
        pts.append((x, y))
    return pts


def sine_path(
    start: tuple[float, float],
    end: tuple[float, float],
    amplitude: float,
    cycles: float,
    n: int = 150,
) -> list[tuple[float, float]]:
    dx = end[0] - start[0]
    dy = end[1] - start[1]
    length = max(1.0, math.hypot(dx, dy))
    nx, ny = -dy / length, dx / length
    pts = []
    for i in range(n):
        t = i / (n - 1)
        wave = amplitude * math.sin(2 * math.pi * cycles * t)
        pts.append((start[0] + dx * t + nx * wave, start[1] + dy * t + ny * wave))
    return pts


def glow_line(
    img: Image.Image,
    points: list[tuple[float, float]],
    color: tuple[int, int, int],
    core_width: int = 5,
    glow_width: int = 28,
    alpha: int = 230,
) -> None:
    layer = Image.new("RGBA", img.size, (255, 255, 255, 0))
    draw = ImageDraw.Draw(layer)
    for width, a in ((glow_width, 24), (int(glow_width * 0.58), 44), (int(glow_width * 0.30), 72)):
        draw.line(scaled_points(points), fill=srgb(color, a), width=max(1, width * S), joint="curve")
    draw.line(scaled_points(points), fill=srgb(color, alpha), width=core_width * S, joint="curve")
    draw.line(scaled_points(points), fill=srgb(WHITE, 150), width=max(1, int(core_width * 0.45 * S)), joint="curve")
    img.alpha_composite(layer.filter(ImageFilter.GaussianBlur(0.25 * S)))


def dash_line(
    draw: ImageDraw.ImageDraw,
    start: tuple[float, float],
    end: tuple[float, float],
    fill: tuple[int, int, int, int],
    width: int = 2,
    dash: float = 15,
    gap: float = 9,
) -> None:
    dx, dy = end[0] - start[0], end[1] - start[1]
    length = math.hypot(dx, dy)
    if length <= 0:
        return
    ux, uy = dx / length, dy / length
    dist = 0.0
    while dist < length:
        a = dist
        z = min(length, dist + dash)
        draw.line(
            [p(start[0] + ux * a, start[1] + uy * a), p(start[0] + ux * z, start[1] + uy * z)],
            fill=fill,
            width=width * S,
        )
        dist += dash + gap


def draw_collision(draw: ImageDraw.ImageDraw, cx: float, cy: float, radius: float = 9) -> None:
    for angle, color in ((0, GOLD), (70, BLUE), (140, ORANGE), (215, TEAL), (292, MAROON)):
        rad = math.radians(angle)
        x0 = cx + math.cos(rad) * radius * 0.35
        y0 = cy + math.sin(rad) * radius * 0.35
        x1 = cx + math.cos(rad) * radius * 2.15
        y1 = cy + math.sin(rad) * radius * 2.15
        draw.line([p(x0, y0), p(x1, y1)], fill=srgb(color, 210), width=max(1, int(2.3 * S)))
    draw.ellipse(b(cx - radius, cy - radius, cx + radius, cy + radius), fill=srgb(WHITE, 245), outline=srgb(NAVY, 170), width=2 * S)
    draw.ellipse(b(cx - radius * 0.45, cy - radius * 0.45, cx + radius * 0.45, cy + radius * 0.45), fill=srgb(GOLD, 235))


def detector_layers(draw: ImageDraw.ImageDraw, cx: float, cy: float, radii: list[float], squash: float = 0.70) -> None:
    colors = [srgb(BLUE, 120), srgb(TEAL, 105), srgb(NAVY, 90), srgb(MAROON, 72)]
    for i, r in enumerate(radii):
        box = b(cx - r, cy - r * squash, cx + r, cy + r * squash)
        draw.arc(box, 203, 518, fill=colors[i % len(colors)], width=max(1, int((3.2 - i * 0.15) * S)))
        for a in (226, 275, 324, 374, 430, 486):
            rad = math.radians(a)
            x = cx + math.cos(rad) * r
            y = cy + math.sin(rad) * r * squash
            draw.ellipse(b(x - 2.4, y - 2.4, x + 2.4, y + 2.4), fill=srgb(colors[i % len(colors)][:3], 115))


def emcal_grid(
    draw: ImageDraw.ImageDraw,
    x: float,
    y: float,
    cols: int,
    rows: int,
    tw: float,
    th: float,
    active: tuple[int, int],
    angle_shift: float = 0,
) -> None:
    for row in range(rows):
        for col in range(cols):
            xx = x + col * (tw + 5) + row * angle_shift
            yy = y + row * (th + 5)
            dcol = abs(col - active[0]) + abs(row - active[1])
            if dcol == 0:
                fill = srgb(GOLD, 236)
                outline = srgb(ORANGE, 210)
            elif dcol == 1:
                fill = srgb((255, 217, 108), 150)
                outline = srgb(GOLD, 110)
            elif (row + col) % 3 == 0:
                fill = srgb(BLUE, 52)
                outline = srgb(BLUE, 92)
            else:
                fill = srgb(SOFT, 96)
                outline = srgb(NAVY, 42)
            draw.rounded_rectangle(b(xx, yy, xx + tw, yy + th), radius=4 * S, fill=fill, outline=outline, width=max(1, S))


def faint_grid(draw: ImageDraw.ImageDraw) -> None:
    for x in range(40, W, 80):
        draw.line([p(x, 26), p(x, H - 28)], fill=srgb(NAVY, 13), width=S)
    for y in range(44, H, 64):
        draw.line([p(28, y), p(W - 28, y)], fill=srgb(NAVY, 10), width=S)


def draw_error_point(draw: ImageDraw.ImageDraw, x: float, y: float, ey: float, color: tuple[int, int, int]) -> None:
    draw.line([p(x, y - ey), p(x, y + ey)], fill=srgb(color, 150), width=2 * S)
    draw.line([p(x - 8, y - ey), p(x + 8, y - ey)], fill=srgb(color, 120), width=2 * S)
    draw.line([p(x - 8, y + ey), p(x + 8, y + ey)], fill=srgb(color, 120), width=2 * S)
    draw.ellipse(b(x - 6, y - 6, x + 6, y + 6), fill=srgb(color, 220), outline=srgb(WHITE, 200), width=2 * S)


def variant_a() -> Image.Image:
    img = canvas()
    draw = ImageDraw.Draw(img)
    detector_layers(draw, 205, 160, [50, 82, 118, 154])
    draw_collision(draw, 205, 160, 8)
    path = cubic((212, 160), (420, 114), (755, 108), (1017, 130), 132)
    glow_line(img, path, GOLD, core_width=5, glow_width=30)
    for i, t in enumerate((0.28, 0.43, 0.58, 0.73)):
        x, y = path[int(t * (len(path) - 1))]
        draw.ellipse(b(x - 5, y - 5, x + 5, y + 5), fill=srgb(GOLD if i % 2 else ORANGE, 190))
    emcal_grid(draw, 1010, 58, 5, 5, 45, 36, active=(1, 2), angle_shift=3.2)
    draw.arc(b(947, 62, 1184, 254), 290, 70, fill=srgb(GOLD, 130), width=3 * S)
    draw.arc(b(961, 78, 1170, 238), 290, 70, fill=srgb(BLUE, 58), width=2 * S)
    return downsample(img)


def variant_b() -> Image.Image:
    img = canvas()
    layer = Image.new("RGBA", img.size, (255, 255, 255, 0))
    draw = ImageDraw.Draw(layer)
    detector_layers(draw, 190, 170, [42, 72, 106, 138], squash=0.82)
    draw_collision(draw, 190, 170, 8)
    for pts, color in [
        (cubic((190, 170), (310, 72), (445, 64), (600, 82), 72), TEAL),
        (cubic((190, 170), (320, 240), (455, 262), (612, 232), 72), MAROON),
        (cubic((190, 170), (330, 154), (420, 195), (525, 196), 60), BLUE),
    ]:
        draw.line(scaled_points(pts), fill=srgb(color, 82), width=3 * S, joint="curve")
    cone_top, cone_bottom = (1016, 94), (1016, 246)
    draw.polygon([p(202, 170), p(*cone_top), p(*cone_bottom)], fill=srgb(GOLD, 15))
    dash_line(draw, (202, 170), cone_top, srgb(GOLD, 120), width=2)
    dash_line(draw, (202, 170), cone_bottom, srgb(GOLD, 120), width=2)
    img.alpha_composite(layer)
    path = sine_path((206, 170), (1007, 170), 5.5, 5.5)
    glow_line(img, path, GOLD, core_width=5, glow_width=25)
    draw = ImageDraw.Draw(img)
    emcal_grid(draw, 1003, 70, 5, 5, 45, 35, active=(0, 2), angle_shift=1.4)
    draw.ellipse(b(951, 103, 1120, 238), outline=srgb(GOLD, 108), width=3 * S)
    draw.ellipse(b(970, 119, 1102, 222), outline=srgb(BLUE, 45), width=2 * S)
    return downsample(img)


def variant_c() -> Image.Image:
    img = canvas()
    draw = ImageDraw.Draw(img)
    faint_grid(draw)
    cx, cy = 610, 164
    for r, color, width in [
        (92, BLUE, 4),
        (142, TEAL, 4),
        (194, NAVY, 3),
        (246, MAROON, 3),
    ]:
        draw.arc(b(cx - r * 1.65, cy - r * 0.74, cx + r * 1.65, cy + r * 0.74), 196, 344, fill=srgb(color, 106), width=width * S)
        draw.arc(b(cx - r * 1.65, cy - r * 0.74, cx + r * 1.65, cy + r * 0.74), 16, 164, fill=srgb(color, 46), width=max(1, (width - 1) * S))
    for a in range(210, 334, 24):
        rad = math.radians(a)
        x0 = cx + math.cos(rad) * 90 * 1.65
        y0 = cy + math.sin(rad) * 90 * 0.74
        x1 = cx + math.cos(rad) * 250 * 1.65
        y1 = cy + math.sin(rad) * 250 * 0.74
        draw.line([p(x0, y0), p(x1, y1)], fill=srgb(NAVY, 25), width=S)
    draw_collision(draw, cx, cy, 7)
    path = sine_path((118, 222), (1132, 102), 11, 7.0)
    glow_line(img, path, GOLD, core_width=4, glow_width=28)
    for t, color in ((0.45, GOLD), (0.55, ORANGE), (0.67, GOLD), (0.79, BLUE)):
        x, y = path[int(t * (len(path) - 1))]
        draw.ellipse(b(x - 7, y - 7, x + 7, y + 7), fill=srgb(color, 170), outline=srgb(WHITE, 160), width=2 * S)
    emcal_grid(draw, 1088, 64, 4, 5, 42, 33, active=(1, 1), angle_shift=2.4)
    return downsample(img)


def variant_d() -> Image.Image:
    img = canvas()
    draw = ImageDraw.Draw(img)
    detector_layers(draw, 174, 165, [44, 78, 114, 150], squash=0.72)
    draw_collision(draw, 174, 165, 8)
    path = cubic((184, 165), (360, 114), (530, 126), (674, 152), 96)
    glow_line(img, path, GOLD, core_width=5, glow_width=27)
    emcal_grid(draw, 628, 92, 4, 4, 39, 34, active=(1, 1), angle_shift=1.5)
    draw.arc(b(585, 78, 830, 240), 306, 58, fill=srgb(GOLD, 108), width=3 * S)
    draw.line([p(842, 162), p(1016, 162)], fill=srgb(NAVY, 56), width=3 * S)
    for x in (885, 928, 972):
        draw.ellipse(b(x - 4, 158, x + 4, 166), fill=srgb(NAVY, 76))
    plot_x0, plot_y0 = 1038, 218
    draw.line([p(plot_x0, plot_y0), p(1248, plot_y0)], fill=srgb(MUTED, 95), width=2 * S)
    draw.line([p(plot_x0, plot_y0), p(plot_x0, 82)], fill=srgb(MUTED, 95), width=2 * S)
    pts = [(1078, 190, 16), (1122, 160, 13), (1170, 127, 11), (1220, 101, 10)]
    draw.line([p(x, y) for x, y, _ in pts], fill=srgb(MAROON, 150), width=3 * S)
    for x, y, ey in pts:
        draw_error_point(draw, x, y, ey, MAROON)
    return downsample(img)


def variant_e() -> Image.Image:
    img = canvas()
    draw = ImageDraw.Draw(img)
    for x in (140, 390, 640, 890, 1140):
        draw.line([p(x, 68), p(x, 252)], fill=srgb(NAVY, 16), width=S)
    detector_layers(draw, 245, 160, [48, 84, 122, 160], squash=0.76)
    draw_collision(draw, 245, 160, 8)
    path = cubic((255, 160), (432, 160), (705, 60), (1078, 160), 140)
    glow_line(img, path, GOLD, core_width=5, glow_width=30)
    for t in (0.36, 0.50, 0.64):
        x, y = path[int(t * (len(path) - 1))]
        draw.arc(b(x - 28, y - 28, x + 28, y + 28), 220, 500, fill=srgb(BLUE, 72), width=2 * S)
    emcal_grid(draw, 1068, 89, 4, 4, 43, 35, active=(0, 2), angle_shift=2.0)
    draw.ellipse(b(1028, 84, 1204, 236), outline=srgb(GOLD, 88), width=3 * S)
    return downsample(img)


VARIANTS = {
    "A_layered_detector_ray": ("Layered detector ray", variant_a),
    "B_isolation_cone_cluster": ("Isolation cone cluster", variant_b),
    "C_blueprint_wave": ("Blueprint photon wave", variant_c),
    "D_detector_to_measurement": ("Detector to measurement", variant_d),
    "E_minimal_prompt_arc": ("Minimal prompt arc", variant_e),
}


def write_contact_sheet(files: list[tuple[str, str, Path]], out_path: Path) -> None:
    pad = 46
    thumb_w, thumb_h = 660, 160
    sheet = Image.new("RGB", (2 * thumb_w + 3 * pad, 3 * thumb_h + 4 * pad), (251, 252, 254))
    draw = ImageDraw.Draw(sheet)
    try:
        label_font = ImageFont.truetype(str(TIMES_BOLD), 22)
    except OSError:
        label_font = ImageFont.load_default()
    for i, (key, label, path) in enumerate(files):
        row, col = divmod(i, 2)
        x = pad + col * (thumb_w + pad)
        y = pad + row * (thumb_h + pad)
        draw.rounded_rectangle((x - 14, y - 14, x + thumb_w + 14, y + thumb_h + 34), radius=16, fill=(255, 255, 255), outline=(226, 231, 237), width=2)
        asset = Image.open(path).convert("RGBA").resize((thumb_w, thumb_h), Image.Resampling.LANCZOS)
        bg = Image.new("RGBA", asset.size, (255, 255, 255, 255))
        bg.alpha_composite(asset)
        sheet.paste(bg.convert("RGB"), (x, y))
        draw.text((x, y + thumb_h + 8), f"{key[0]}: {label}", font=label_font, fill=NAVY)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    sheet.save(out_path, quality=95)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workspace", type=Path, default=DEFAULT_WORKSPACE)
    parser.add_argument("--only", choices=["all", *VARIANTS.keys()], default="all")
    args = parser.parse_args()

    out_dir = args.workspace / "output/photon_detector_motif_variants"
    out_dir.mkdir(parents=True, exist_ok=True)

    rendered: list[tuple[str, str, Path]] = []
    for key, (label, render) in VARIANTS.items():
        if args.only not in ("all", key):
            continue
        path = out_dir / f"hp2026_photon_detector_motif_{key}.png"
        render().save(path)
        rendered.append((key, label, path))
        print(path)

    if args.only == "all":
        contact = out_dir / "hp2026_photon_detector_motif_contact_sheet.png"
        write_contact_sheet(rendered, contact)
        print(contact)

    manifest = {
        "title": "HP2026 title-slide photon detector motif variants",
        "size_px": [W, H],
        "background": "transparent",
        "intended_placement": "below the Justin/on-behalf block on the C3 HP2026 title slide",
        "outputs": [str(path) for _, _, path in rendered],
        "contact_sheet": str(out_dir / "hp2026_photon_detector_motif_contact_sheet.png") if args.only == "all" else None,
        "notes": [
            "Graphic-only transparent PNGs; no text, logos, slide number, or provenance footer in the individual motif assets.",
            "Designed to be scaled down and placed under the collaboration line without covering title text or the Vanderbilt image.",
            "No Google Slides mutation.",
        ],
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
