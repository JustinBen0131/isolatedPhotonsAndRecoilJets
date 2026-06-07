#!/usr/bin/env python3
"""Build a collaboration-facing stitching closure slide for embedded samples."""

from __future__ import annotations

import json
import math
from io import BytesIO
from pathlib import Path
from textwrap import wrap

import numpy as np
from PIL import Image, ImageDraw, ImageFont


REPO = Path(__file__).resolve().parents[4]
OUT = REPO / "dataOutput/slides/wp_gammajets_6_1_26/stitching_closure_20260604"

PHOTON_PLOT = (
    REPO
    / "dataOutput/stitchDiagnostics/focus21_photon_variants_20260521"
    / "photon12plus20_bounded21_blair_takeaway.png"
)
PHOTON_CSV = (
    REPO
    / "dataOutput/stitchDiagnostics/focus21_photon_variants_20260521"
    / "photon_bounded21_extended_bins.csv"
)
PHOTON_SUMMARY = (
    REPO
    / "dataOutput/stitchDiagnostics/focus21_photon_variants_20260521"
    / "photon12plus20_bounded21_slide4_replacement_summary.json"
)
INCLUSIVE_PLOT_FULLSLIDE = (
    REPO
    / "dataOutput/stitchDiagnostics/jet40_slide6_spectrum_20260526"
    / "embeddedInclusiveJet12plus20plus30plus40_slide6_replacement.png"
)
INCLUSIVE_SUMMARY = (
    REPO
    / "dataOutput/stitchDiagnostics/jet40_slide6_spectrum_20260526"
    / "embeddedInclusiveJet12plus20plus30plus40_slide6_replacement_summary.json"
)
INCLUSIVE_CSV = (
    REPO
    / "dataOutput/stitchDiagnostics/jet40_slide6_spectrum_20260526"
    / "inclusive4_stitch_bins.csv"
)

W, H = 2560, 1440
FONT_DIR = Path("/System/Library/Fonts/Supplemental")
TIMES = FONT_DIR / "Times New Roman.ttf"
TIMES_BOLD = FONT_DIR / "Times New Roman Bold.ttf"
TIMES_ITALIC = FONT_DIR / "Times New Roman Italic.ttf"
TIMES_BOLD_ITALIC = FONT_DIR / "Times New Roman Bold Italic.ttf"

INK = (18, 24, 38)
MUTED = (78, 86, 101)
LINE = (208, 216, 228)
PANEL = (248, 250, 252)
BLUE = (33, 76, 195)
RED = (210, 35, 42)
ORANGE = (238, 114, 18)
MAGENTA = (204, 43, 168)
GREEN = (35, 150, 55)
SOFT_BLUE = (239, 245, 255)
SOFT_GREEN = (238, 251, 244)
SOFT_AMBER = (255, 249, 235)
SIGNAL_ACCENT = (26, 111, 116)
SIGNAL_SOFT = (235, 248, 248)
SIGNAL_LINE = (174, 214, 216)
BACKGROUND_ACCENT = (92, 78, 152)
BACKGROUND_SOFT = (244, 241, 250)
BACKGROUND_LINE = (204, 198, 230)
RESULT_ACCENT = (42, 91, 117)
RESULT_SOFT = (241, 246, 249)


def font(size: int, *, bold: bool = False, italic: bool = False) -> ImageFont.FreeTypeFont:
    path = TIMES
    if bold and italic:
        path = TIMES_BOLD_ITALIC
    elif bold:
        path = TIMES_BOLD
    elif italic:
        path = TIMES_ITALIC
    try:
        return ImageFont.truetype(str(path), size)
    except OSError:
        return ImageFont.load_default()


F = {
    "title": font(60, bold=True),
    "subtitle": font(30),
    "panel_title": font(39, bold=True),
    "body": font(30),
    "body_bold": font(30, bold=True),
    "band_body": font(30),
    "band_bold": font(33, bold=True),
    "small": font(27),
    "small_bold": font(27, bold=True),
    "tiny": font(24),
    "tiny_bold": font(24, bold=True),
    "card_title": font(44, bold=True),
    "card_body": font(37),
    "card_body_bold": font(37, bold=True),
    "row_name": font(31, bold=True),
    "row": font(28),
    "row_bold": font(28, bold=True),
    "table_header": font(26, bold=True),
    "table_body": font(29),
    "table_body_small": font(27),
    "sub": font(18),
    "sub_bold": font(18, bold=True),
    "row_sub_bold": font(21, bold=True),
}

def load_json(path: Path) -> dict:
    return json.loads(path.read_text())


def read_csv_points(path: Path, variant: str) -> list[dict[str, float]]:
    points: list[dict[str, float]] = []
    header_seen = False
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or line.startswith("["):
            continue
        if line.startswith("variant,"):
            header_seen = True
            continue
        if not header_seen:
            continue
        tok = line.split(",")
        if len(tok) != 7 or tok[0] != variant:
            continue
        lo = float(tok[2])
        hi = float(tok[3])
        points.append(
            {
                "x": float(tok[4]),
                "ex": 0.5 * (hi - lo),
                "y": float(tok[5]),
                "ey": float(tok[6]),
            }
        )
    return points


def rebin_points(points: list[dict[str, float]], x_min: float, x_max: float, width: float) -> list[dict[str, float]]:
    out: list[dict[str, float]] = []
    lo = x_min
    while lo < x_max - 1.0e-9:
        hi = lo + width
        y = 0.0
        err2 = 0.0
        for p in points:
            plo = p["x"] - p["ex"]
            phi = p["x"] + p["ex"]
            if plo >= lo - 1.0e-9 and phi <= hi + 1.0e-9:
                y += p["y"]
                err2 += p["ey"] * p["ey"]
        if y > 0.0:
            out.append({"x": 0.5 * (lo + hi), "ex": 0.5 * width, "y": y, "ey": math.sqrt(err2)})
        lo = hi
    return out


def fit_modified_power_law(points: list[dict[str, float]]):
    xs = np.array([p["x"] for p in points if p["y"] > 0.0], dtype=float)
    ys = np.array([p["y"] for p in points if p["y"] > 0.0], dtype=float)
    logx = np.log(xs)
    l_inv = np.log(1.0 / xs)
    design = np.column_stack([np.ones_like(xs), l_inv, logx * l_inv, xs * l_inv])
    coeff, *_ = np.linalg.lstsq(design, np.log(ys), rcond=None)

    def eval_y(x: float) -> float:
        lx = math.log(x)
        li = math.log(1.0 / x)
        logy = coeff[0] + coeff[1] * li + coeff[2] * lx * li + coeff[3] * x * li
        return math.exp(logy) if math.isfinite(logy) else 0.0

    return eval_y


def render_stitch_plot(
    *,
    title: str,
    subtitle: str,
    points: list[dict[str, float]],
    pieces: list[dict],
    x_min: float,
    x_max: float,
    y_min: float,
    y_max: float,
    ratio_ylim: tuple[float, float],
    x_label: str,
    y_label: str,
    ratio_label: str,
    boundaries: list[float],
    jump_text: str,
    accent: tuple[int, int, int],
    sample_label: str,
) -> Image.Image:
    import matplotlib.pyplot as plt
    from matplotlib.ticker import LogLocator, NullFormatter

    plt.rcParams.update(
        {
            "font.family": "Times New Roman",
            "axes.linewidth": 1.2,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )
    fit = fit_modified_power_law(points)
    ratio_values: list[float] = []
    for piece in pieces:
        vals = [p for p in points if piece["lo"] <= p["x"] < piece["hi"] and p["y"] > 0.0]
        if not vals:
            continue
        for p in vals:
            ref = fit(float(p["x"]))
            if ref > 0.0 and math.isfinite(ref):
                ratio_values.append(p["y"] / ref)
    ratio_rms = 0.0
    if ratio_values:
        ratios = np.array(ratio_values, dtype=float)
        ratio_rms = float(np.sqrt(np.mean(np.square(ratios - 1.0))))

    fig = plt.figure(figsize=(9.8, 6.0), dpi=100, facecolor="white")
    gs = fig.add_gridspec(
        2,
        1,
        height_ratios=[3.6, 1.15],
        left=0.115,
        right=0.965,
        bottom=0.125,
        top=0.875,
        hspace=0.05,
    )
    ax = fig.add_subplot(gs[0])
    rax = fig.add_subplot(gs[1], sharex=ax)

    label_box = {
        "boxstyle": "round,pad=0.28",
        "facecolor": "white",
        "edgecolor": "0.82",
        "alpha": 0.92,
        "linewidth": 0.8,
    }
    ax.text(
        0.965,
        0.955,
        rf"$\it{{\bf{{sPHENIX}}}}$ Internal" + f"\n{sample_label}",
        transform=ax.transAxes,
        fontsize=14.0,
        ha="right",
        va="top",
        linespacing=1.05,
        color=np.array(INK) / 255.0,
        bbox=label_box,
        zorder=6,
    )
    ax.text(
        0.965,
        0.735,
        title,
        transform=ax.transAxes,
        fontsize=15.0,
        ha="right",
        va="top",
        fontweight="bold",
        color=np.array(INK) / 255.0,
        bbox=label_box,
        zorder=6,
    )

    xs_fit = np.linspace(x_min, x_max, 500)
    ys_fit = np.array([fit(float(x)) for x in xs_fit])
    ax.plot(xs_fit, ys_fit, color="black", linestyle="--", linewidth=1.8, label="Modified power-law fit", zorder=1)
    if ratio_rms > 0.0:
        band_lo = max(ratio_ylim[0], 1.0 - ratio_rms)
        band_hi = min(ratio_ylim[1], 1.0 + ratio_rms)
        rax.axhspan(
            band_lo,
            band_hi,
            color=np.array(accent) / 255.0,
            alpha=0.12,
            linewidth=0,
            zorder=0,
        )

    for piece in pieces:
        vals = [p for p in points if piece["lo"] <= p["x"] < piece["hi"] and p["y"] > 0.0]
        if not vals:
            continue
        marker = {"circle": "o", "square": "s", "triangle": "^", "down": "v"}[piece["marker"]]
        color = np.array(piece["color"]) / 255.0
        x = np.array([p["x"] for p in vals])
        y = np.array([p["y"] for p in vals])
        yerr = np.array([p["ey"] for p in vals])
        ax.errorbar(
            x,
            y,
            yerr=yerr,
            fmt=marker,
            markersize=3.7,
            linestyle="none",
            elinewidth=0.9,
            capsize=1.6,
            color=color,
            label=piece["label"],
            zorder=3,
        )
        ref = np.array([fit(float(xx)) for xx in x])
        ratio = y / ref
        ratio_err = yerr / ref
        rax.errorbar(
            x,
            ratio,
            yerr=ratio_err,
            fmt=marker,
            markersize=3.7,
            linestyle="none",
            elinewidth=0.9,
            capsize=1.6,
            color=color,
            zorder=3,
        )

    ax.set_yscale("log")
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    rax.set_ylim(*ratio_ylim)
    rax.axhline(1.0, color="0.55", linestyle=":", linewidth=1.0, zorder=1)
    ax.set_ylabel(y_label, fontsize=15.0)
    rax.set_ylabel(ratio_label, fontsize=13.5)
    rax.set_xlabel(x_label, fontsize=15.0, loc="right")
    ax.tick_params(labelbottom=False, labelsize=13, length=6)
    rax.tick_params(labelsize=13, length=6)
    ax.yaxis.set_major_locator(LogLocator(base=10))
    ax.yaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10) * 0.1))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.grid(axis="y", which="major", linestyle=":", linewidth=0.7, color="0.86")
    rax.grid(axis="y", which="major", linestyle=":", linewidth=0.7, color="0.86")
    ax.text(
        0.035,
        0.105,
        subtitle,
        transform=ax.transAxes,
        fontsize=18.0,
        color=np.array(INK) / 255.0,
        va="bottom",
        linespacing=1.18,
        bbox={
            "boxstyle": "round,pad=0.30",
            "facecolor": "white",
            "edgecolor": "0.82",
            "alpha": 0.88,
            "linewidth": 0.8,
        },
    )
    rax.text(
        0.72,
        0.83,
        f"RMS band = {100.0 * ratio_rms:.1f}%",
        transform=rax.transAxes,
        fontsize=13.0,
        color=np.array(INK) / 255.0,
        fontweight="normal",
    )

    handles, labels = ax.get_legend_handles_labels()
    if handles:
        fit_handle = handles[0]
        fit_label = labels[0]
        data_handles = handles[1:]
        data_labels = labels[1:]
        ax.legend(
            data_handles + [fit_handle],
            data_labels + [fit_label],
            loc="upper right",
            bbox_to_anchor=(0.985, 0.93),
            fontsize=15.0,
            frameon=False,
            handlelength=1.9,
            labelspacing=0.32,
            borderaxespad=0.25,
            markerscale=1.1,
        )

    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=100, facecolor="white")
    plt.close(fig)
    buf.seek(0)
    return Image.open(buf).convert("RGB")


def build_standardized_plots(photon: dict, inclusive: dict) -> tuple[Image.Image, Image.Image, Path, Path]:
    photon_points_raw = read_csv_points(PHOTON_CSV, "bounded_12_21_ge21")
    photon_points = [p for p in photon_points_raw if 12.0 <= p["x"] <= 40.0 and p["y"] > 0.0]
    inclusive_points = [p for p in read_csv_points(INCLUSIVE_CSV, "inclusive4") if 12.0 <= p["x"] <= 50.0 and p["y"] > 0.0]

    photon_plot = render_stitch_plot(
        title="PhotonJet12+20 bounded 21 GeV ownership",
        subtitle="PhotonJet12: 12-21 GeV\nPhotonJet20: >=21 GeV",
        points=photon_points,
        pieces=[
            {"label": "PhotonJet12 stitched", "lo": 12.0, "hi": 21.0, "color": BLUE, "marker": "circle"},
            {"label": "PhotonJet20 stitched", "lo": 21.0, "hi": 40.1, "color": RED, "marker": "square"},
        ],
        x_min=12.0,
        x_max=40.0,
        y_min=2.0e3,
        y_max=1.6e8,
        ratio_ylim=(0.90, 1.10),
        x_label=r"Truth photon filter $p_T$ [GeV]",
        y_label="Weighted stitched entries / bin",
        ratio_label="stitched / fit",
        boundaries=[21.0],
        jump_text=f"smooth handoff: jump = {photon['jump_21_over_19_21']:.3f}",
        accent=SIGNAL_ACCENT,
        sample_label="Embedded Photon+Jet 12+20",
    )

    inclusive_plot = render_stitch_plot(
        title="Inclusive Jet12+20+30+40 ownership",
        subtitle="Jet12: 12-21 GeV     Jet20: 21-31 GeV\nJet30: 31-41 GeV     Jet40: >=41 GeV",
        points=inclusive_points,
        pieces=[
            {"label": "Jet12 stitched", "lo": 12.0, "hi": 21.0, "color": BLUE, "marker": "circle"},
            {"label": "Jet20 stitched", "lo": 21.0, "hi": 31.0, "color": ORANGE, "marker": "square"},
            {"label": "Jet30 stitched", "lo": 31.0, "hi": 41.0, "color": MAGENTA, "marker": "triangle"},
            {"label": "Jet40 stitched", "lo": 41.0, "hi": 50.1, "color": GREEN, "marker": "down"},
        ],
        x_min=12.0,
        x_max=50.0,
        y_min=1.0,
        y_max=1.2e6,
        ratio_ylim=(0.88, 1.10),
        x_label=r"Truth jet filter $p_T$ [GeV]",
        y_label="Weighted stitched entries / bin",
        ratio_label="stitched / fit",
        boundaries=[21.0, 31.0, 41.0],
        jump_text=(
            f"boundary jumps: {inclusive['jump_21_over_19_21']:.3f}, "
            f"{inclusive['jump_31_over_29_31']:.3f}, {inclusive['jump_41_over_39_41']:.3f}"
        ),
        accent=BACKGROUND_ACCENT,
        sample_label="Embedded Inclusive Jet 12+20+30+40",
    )

    photon_png = OUT / "standardized_photonjet12plus20_stitch_plot.png"
    inclusive_png = OUT / "standardized_inclusivejet12plus20plus30plus40_stitch_plot.png"
    photon_plot.save(photon_png)
    inclusive_plot.save(inclusive_png)
    return photon_plot, inclusive_plot, photon_png, inclusive_png


def rounded(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], fill, outline=LINE, width=2, radius=16) -> None:
    draw.rounded_rectangle(box, radius=radius, fill=fill, outline=outline, width=width)


def text(draw: ImageDraw.ImageDraw, xy: tuple[int, int], msg: str, fnt, fill=INK, anchor=None) -> None:
    draw.text(xy, msg, font=fnt, fill=fill, anchor=anchor)


def text_advance(draw: ImageDraw.ImageDraw, x: int, y: int, msg: str, fnt, fill=INK) -> int:
    draw.text((x, y), msg, font=fnt, fill=fill)
    bbox = draw.textbbox((x, y), msg, font=fnt)
    return bbox[2]


def rich_line(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    parts: list[tuple[str, ImageFont.FreeTypeFont, tuple[int, int, int]]],
) -> int:
    cursor = x
    max_h = 0
    for msg, fnt, fill in parts:
        bbox = draw.textbbox((cursor, y), msg, font=fnt)
        draw.text((cursor, y), msg, font=fnt, fill=fill)
        cursor += bbox[2] - bbox[0]
        max_h = max(max_h, bbox[3] - bbox[1])
    return y + max_h


def sigma_sub(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    sub: str,
    fnt: ImageFont.FreeTypeFont,
    sub_fnt: ImageFont.FreeTypeFont,
    fill=INK,
) -> int:
    draw.text((x, y), "σ", font=fnt, fill=fill)
    sigma_box = draw.textbbox((x, y), "σ", font=fnt)
    sub_x = sigma_box[2] + 1
    sub_y = y + int(fnt.size * 0.48)
    draw.text((sub_x, sub_y), sub, font=sub_fnt, fill=fill)
    sub_box = draw.textbbox((sub_x, sub_y), sub, font=sub_fnt)
    return sub_box[2]


def sigma_equals(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    sub: str,
    value: str,
    fnt: ImageFont.FreeTypeFont,
    sub_fnt: ImageFont.FreeTypeFont,
    fill=INK,
) -> int:
    cursor = sigma_sub(draw, x, y, sub, fnt, sub_fnt, fill)
    draw.text((cursor + 5, y), f" = {value}", font=fnt, fill=fill)
    bbox = draw.textbbox((cursor + 5, y), f" = {value}", font=fnt)
    return bbox[2]


def wrapped(draw: ImageDraw.ImageDraw, xy: tuple[int, int], msg: str, fnt, fill=INK, width_chars=68, spacing=7) -> int:
    x, y = xy
    lines: list[str] = []
    for para in msg.split("\n"):
        lines.extend(wrap(para, width=width_chars) or [""])
    for line in lines:
        draw.text((x, y), line, font=fnt, fill=fill)
        y += fnt.size + spacing
    return y


def wrapped_pixels(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    msg: str,
    fnt,
    *,
    fill=INK,
    max_width: int,
    spacing: int = 5,
) -> int:
    x, y = xy
    for para in msg.split("\n"):
        words = para.split()
        lines: list[str] = []
        line = ""
        for word in words:
            candidate = f"{line} {word}".strip()
            bbox = draw.textbbox((0, 0), candidate, font=fnt)
            if line and bbox[2] - bbox[0] > max_width:
                lines.append(line)
                line = word
            else:
                line = candidate
        lines.append(line)
        for line in lines:
            draw.text((x, y), line, font=fnt, fill=fill)
            y += fnt.size + spacing
    return y


def paste_contained(
    canvas: Image.Image,
    img: Image.Image,
    box: tuple[int, int, int, int],
    *,
    fill=(255, 255, 255),
) -> tuple[int, int, int, int]:
    x0, y0, x1, y1 = box
    bw, bh = x1 - x0, y1 - y0
    scale = min(bw / img.width, bh / img.height)
    nw, nh = int(img.width * scale), int(img.height * scale)
    resized = img.resize((nw, nh), Image.Resampling.LANCZOS)
    bg = Image.new("RGB", (bw, bh), fill)
    ox, oy = (bw - nw) // 2, (bh - nh) // 2
    bg.paste(resized, (ox, oy))
    canvas.paste(bg, (x0, y0))
    return (x0 + ox, y0 + oy, x0 + ox + nw, y0 + oy + nh)


def add_plot_panel(
    canvas: Image.Image,
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    plot: Image.Image,
    *,
    accent: tuple[int, int, int],
    outline: tuple[int, int, int],
    crop: tuple[int, int, int, int] | None = None,
) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, (255, 255, 255), outline=outline, width=3, radius=12)
    plot_box = (x0 + 14, y0 + 14, x1 - 14, y1 - 14)
    plot_img = plot.crop(crop) if crop else plot
    paste_contained(canvas, plot_img.convert("RGB"), plot_box)


def fmt_weight(value: float, decimals: int = 2) -> str:
    if value >= 100:
        return f"{value:.0f}x"
    if value >= 10:
        return f"{value:.1f}x"
    return f"{value:.{decimals}f}x"


def draw_weight_table(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    title_msg: str,
    rows: list[dict],
    accent,
    soft,
    anchor_label: str,
) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, soft, outline=accent, width=3, radius=12)
    draw.rectangle((x0 + 22, y0 + 22, x0 + 36, y0 + 72), fill=accent)
    text(draw, (x0 + 58, y0 + 15), title_msg, F["panel_title"], INK)
    header_y = y0 + 93
    cols = [x0 + 34, x0 + 300, x0 + 540, x0 + 805]
    headers = ["sample", "owned window", "σ_eff [pb]", f"rel. w / {anchor_label}"]
    for col, header in zip(cols, headers):
        text(draw, (col, header_y), header, F["table_header"], MUTED)
    draw.line((x0 + 30, header_y + 38, x1 - 30, header_y + 38), fill=accent, width=2)
    row_start = header_y + 55
    row_font = F["table_body"] if len(rows) <= 3 else F["table_body_small"]
    row_step = min(56, max(42, int((y1 - row_start - 24) / max(1, len(rows)))))
    for i, row in enumerate(rows):
        y = row_start + i * row_step
        if i % 2 == 0:
            draw.rounded_rectangle((x0 + 24, y - 7, x1 - 24, y + row_step - 9), radius=7, fill=(255, 255, 255))
        for col, key in zip(cols, ["sample", "window", "sigma", "rel"]):
            fill = row.get("color", INK) if key == "sample" else INK
            text(draw, (col, y), row[key], row_font, fill)


def build() -> tuple[Path, Path, Path]:
    OUT.mkdir(parents=True, exist_ok=True)
    photon = load_json(PHOTON_SUMMARY)
    inclusive = load_json(INCLUSIVE_SUMMARY)

    canvas = Image.new("RGB", (W, H), "white")
    draw = ImageDraw.Draw(canvas)

    title = "embedded stitch-weight closure audit"
    text(draw, (72, 46), title, F["title"], INK)
    subtitle = "Same closure logic as pp: identify the generator slice, count the effective cross section, then stitch by the owned truth-pT window."
    text(draw, (76, 126), subtitle, F["subtitle"], MUTED)

    photon_img, inclusive_img, photon_plot_png, inclusive_plot_png = build_standardized_plots(photon, inclusive)

    add_plot_panel(
        canvas,
        draw,
        (70, 174, 1248, 890),
        photon_img,
        accent=SIGNAL_ACCENT,
        outline=SIGNAL_LINE,
    )
    add_plot_panel(
        canvas,
        draw,
        (1322, 174, 2490, 890),
        inclusive_img,
        accent=BACKGROUND_ACCENT,
        outline=BACKGROUND_LINE,
    )

    draw_weight_table(
        draw,
        (70, 930, 1248, 1268),
        "PhotonJet 12+20 embedded weights",
        [
            {
                "sample": "PhotonJet12",
                "window": "12-21",
                "sigma": f"{photon['photonjet12_sigma_eff_pb']:.1f}",
                "rel": fmt_weight(photon["photonjet12_merge_scale"]),
                "color": BLUE,
            },
            {
                "sample": "PhotonJet20",
                "window": ">=21",
                "sigma": f"{photon['photonjet20_sigma_eff_pb']:.2f}",
                "rel": "1x",
                "color": RED,
            },
        ],
        SIGNAL_ACCENT,
        SIGNAL_SOFT,
        "PhotonJet20",
    )
    draw_weight_table(
        draw,
        (1322, 930, 2490, 1268),
        "Inclusive Jet 12+20+30+40 embedded weights",
        [
            {
                "sample": "Jet12",
                "window": "12-21",
                "sigma": f"{inclusive['jet12_sigma_eff_pb']:.3g}",
                "rel": fmt_weight(inclusive["jet12_relative_weight_to_jet40"], 0),
                "color": BLUE,
            },
            {
                "sample": "Jet20",
                "window": "21-31",
                "sigma": f"{inclusive['jet20_sigma_eff_pb']:.4g}",
                "rel": fmt_weight(inclusive["jet20_relative_weight_to_jet40"], 1),
                "color": ORANGE,
            },
            {
                "sample": "Jet30",
                "window": "31-41",
                "sigma": f"{inclusive['jet30_sigma_eff_pb']:.1f}",
                "rel": fmt_weight(inclusive["jet30_relative_weight_to_jet40"], 2),
                "color": MAGENTA,
            },
            {
                "sample": "Jet40",
                "window": ">=41",
                "sigma": f"{inclusive['jet40_sigma_eff_pb']:.1f}",
                "rel": "1x",
                "color": GREEN,
            },
        ],
        BACKGROUND_ACCENT,
        BACKGROUND_SOFT,
        "Jet40",
    )

    note_box = (70, 1287, 2490, 1407)
    rounded(draw, note_box, SOFT_BLUE, outline=(147, 197, 253), width=2, radius=10)
    left_x, right_x = 96, 1540
    text(draw, (left_x, 1301), "Closure result", F["band_bold"], SIGNAL_ACCENT)
    result = (
        f"Photon handoff = {photon['jump_21_over_19_21']:.3f}; inclusive handoffs = "
        f"{inclusive['jump_21_over_19_21']:.3f}, {inclusive['jump_31_over_29_31']:.3f}, "
        f"{inclusive['jump_41_over_39_41']:.3f}. Counted σ_eff weights in owned filter windows."
    )
    wrapped_pixels(draw, (left_x, 1343), result, F["band_body"], fill=INK, max_width=1350, spacing=3)
    text(draw, (right_x, 1301), "Interpretation", F["band_bold"], BACKGROUND_ACCENT)
    wrapped_pixels(
        draw,
        (right_x, 1343),
        "Signal and inclusive-background embedded samples close with the same ownership logic used for pp.",
        F["band_body"],
        fill=INK,
        max_width=850,
        spacing=3,
    )

    png = OUT / "embedded_signal_inclusive_stitching_closure_slide.png"
    script = OUT / "embedded_signal_inclusive_stitching_closure_slide_script.md"
    manifest = OUT / "embedded_signal_inclusive_stitching_closure_slide_manifest.json"
    canvas.save(png)

    script.write_text(
        """# WP GammaJets Stitching Closure Slide Script - Embedded Signal And Background

This slide is meant to close the loop on the embedded stitching inputs. The point is not that we have many separate samples; the point is that each sample has a clear ownership window, and the merged spectrum is continuous once those windows and weights are applied.

On the left is the embedded photon-plus-jet signal side. PhotonJet12 owns truth-filter photons from 12 to 21 GeV, and PhotonJet20 owns the region above 21 GeV. The effective cross sections, sigma_eff, come from counting how many generator events pass the photon production filter in those windows. With those weights applied, the handoff is smooth at about the one percent level.

On the right is the embedded inclusive-jet background side. The same idea is extended to four samples: Jet12, Jet20, Jet30, and Jet40, with ownership windows at 12-21, 21-31, 31-41, and above 41 GeV. Jet40 is the relative-weight anchor. The boundary checks at 21, 31, and 41 GeV are all near unity, so the combined spectrum is continuous rather than double counted or gapped.

The practical takeaway is that the signal and inclusive-background embedded inputs are now stitched with explicit event ownership and counted effective cross sections. That gives us a controlled foundation for the downstream photon-ID and working-point comparisons.

The dashed curve is the modified power-law fit used only as a smooth closure reference: A*(1/x)^(b + c ln(x) + d x). The binned stitched points are not connected to each other.
""",
        encoding="utf-8",
    )

    manifest.write_text(
        json.dumps(
            {
                "png": str(png),
                "script": str(script),
                "google_slides_mutation": False,
                "source_plots": {
                    "standardized_photon": str(photon_plot_png),
                    "standardized_inclusive": str(inclusive_plot_png),
                    "legacy_photon_root_macro_png": str(PHOTON_PLOT),
                    "legacy_inclusive_root_macro_png": str(INCLUSIVE_PLOT_FULLSLIDE),
                },
                "source_csvs": {
                    "photon": str(PHOTON_CSV),
                    "inclusive": str(INCLUSIVE_CSV),
                },
                "source_summaries": {
                    "photon": str(PHOTON_SUMMARY),
                    "inclusive": str(INCLUSIVE_SUMMARY),
                },
                "photon": {
                    "windows": {
                        "PhotonJet12": photon["photonjet12_gate"],
                        "PhotonJet20": photon["photonjet20_gate"],
                    },
                    "sigma_eff_pb": {
                        "PhotonJet12": photon["photonjet12_sigma_eff_pb"],
                        "PhotonJet20": photon["photonjet20_sigma_eff_pb"],
                    },
                    "relative_weights": {
                        "PhotonJet12": photon["photonjet12_merge_scale"],
                        "PhotonJet20": 1.0,
                    },
                    "boundary_jump": photon["jump_21_over_19_21"],
                },
                "inclusive": {
                    "windows": {
                        "Jet12": inclusive["jet12_gate"],
                        "Jet20": inclusive["jet20_gate"],
                        "Jet30": inclusive["jet30_gate"],
                        "Jet40": inclusive["jet40_gate"],
                    },
                    "sigma_eff_pb": {
                        "Jet12": inclusive["jet12_sigma_eff_pb"],
                        "Jet20": inclusive["jet20_sigma_eff_pb"],
                        "Jet30": inclusive["jet30_sigma_eff_pb"],
                        "Jet40": inclusive["jet40_sigma_eff_pb"],
                    },
                    "relative_weights_to_jet40": {
                        "Jet12": inclusive["jet12_relative_weight_to_jet40"],
                        "Jet20": inclusive["jet20_relative_weight_to_jet40"],
                        "Jet30": inclusive["jet30_relative_weight_to_jet40"],
                        "Jet40": inclusive["jet40_relative_weight_to_jet40"],
                    },
                    "boundary_jumps": {
                        "21": inclusive["jump_21_over_19_21"],
                        "31": inclusive["jump_31_over_29_31"],
                        "41": inclusive["jump_41_over_39_41"],
                    },
                },
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    return png, script, manifest


if __name__ == "__main__":
    paths = build()
    for path in paths:
        print(path)
