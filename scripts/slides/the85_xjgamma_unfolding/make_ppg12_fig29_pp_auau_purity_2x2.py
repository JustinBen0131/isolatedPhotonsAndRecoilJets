#!/usr/bin/env python3
"""Build a PPG12 Fig. 29 style pp-vs-AuAu purity slide candidate."""

from __future__ import annotations

import csv
import json
import math
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image, ImageDraw, ImageFont
from scipy.optimize import curve_fit
from scipy.special import erf


REPO = Path(__file__).resolve().parents[3]
OUT_DIR = REPO / "dataOutput/the85_auau_xjgamma_unfolding_push/slides/fig29_pp_auau_purity_2x2"
OUT_DIR.mkdir(parents=True, exist_ok=True)

IAN_PDF = REPO / "usefulDocs/ppg12/current/PPG12_CURRENT_IAN_2026-05-21_v4.pdf"
PDFTOPPM = Path("/Users/patsfan753/.cache/codex-runtimes/codex-primary-runtime/dependencies/bin/pdftoppm")
PAGE_RENDER = OUT_DIR / "ppg12_current_ian_page39-039.png"
PP_RAW_CROP = OUT_DIR / "ppg12_fig29_pp_raw_corrected_panel.png"
PP_FIT_CROP = OUT_DIR / "ppg12_fig29_pp_fit_panel.png"

POINTS_CSV = REPO / (
    "dataOutput/the85_auau_xjgamma_unfolding_push/slides/"
    "slide02_purity_leakage_corrected_available_bins_1x3_v3_points.csv"
)

OUT_PNG = OUT_DIR / "ppg12_fig29_style_pp_vs_auau020_purity_2x2.png"
OUT_MANIFEST = OUT_DIR / "ppg12_fig29_style_pp_vs_auau020_purity_2x2_manifest.json"
OUT_SCRIPT = OUT_DIR / "ppg12_fig29_style_pp_vs_auau020_purity_2x2_speaker_script.md"


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    candidates = [
        "/System/Library/Fonts/Supplemental/Times New Roman Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Times New Roman.ttf",
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/Library/Fonts/Arial Bold.ttf" if bold else "/Library/Fonts/Arial.ttf",
    ]
    for path in candidates:
        p = Path(path)
        if p.exists():
            return ImageFont.truetype(str(p), size=size)
    return ImageFont.load_default()


def ensure_pp_crops() -> None:
    if not PAGE_RENDER.exists() or PAGE_RENDER.stat().st_size < 100_000:
        prefix = OUT_DIR / "ppg12_current_ian_page39"
        subprocess.run(
            [
                str(PDFTOPPM),
                "-png",
                "-f",
                "39",
                "-l",
                "39",
                "-r",
                "220",
                str(IAN_PDF),
                str(prefix),
            ],
            check=True,
        )

    page = Image.open(PAGE_RENDER).convert("RGB")
    # Crops are from current IAN page 39 at 220 dpi, Figure 29 only.
    page.crop((285, 245, 920, 895)).save(PP_RAW_CROP)
    page.crop((920, 245, 1610, 895)).save(PP_FIT_CROP)


def trim_white_margins(img: Image.Image, tolerance: int = 248, pad: int = 8) -> Image.Image:
    rgb = img.convert("RGB")
    arr = np.asarray(rgb)
    mask = np.any(arr < tolerance, axis=2)
    if not mask.any():
        return rgb
    ys, xs = np.where(mask)
    left = max(int(xs.min()) - pad, 0)
    top = max(int(ys.min()) - pad, 0)
    right = min(int(xs.max()) + pad + 1, rgb.width)
    bottom = min(int(ys.max()) + pad + 1, rgb.height)
    return rgb.crop((left, top, right, bottom))


def read_auau020_points() -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    with POINTS_CSV.open() as handle:
        for row in csv.DictReader(handle):
            if row["label"] != "Au+Au 0-20%":
                continue
            parsed: dict[str, float] = {}
            for key, value in row.items():
                try:
                    parsed[key] = float(value)
                except ValueError:
                    parsed[key] = math.nan
            rows.append(parsed)
    if not rows:
        raise RuntimeError(f"no Au+Au 0-20% rows in {POINTS_CSV}")
    return rows


def pade11(x: np.ndarray, a: float, b: float, c: float) -> np.ndarray:
    return (a + b * x) / (1.0 + c * x)


def erf_alt(x: np.ndarray, a: float, b: float, c: float, d: float) -> np.ndarray:
    return a + b * 0.5 * (1.0 + erf((x - c) / np.maximum(d, 1.0e-6)))


def fit_curves(rows: list[dict[str, float]]) -> dict[str, object]:
    x = np.array([r["pt_mid"] for r in rows], dtype=float)
    y = np.array([r["corrected_purity"] for r in rows], dtype=float)
    ey = np.array([r["corrected_purity_err"] for r in rows], dtype=float)
    mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(ey) & (ey > 0.0) & (ey < 1.0)
    x = x[mask]
    y = y[mask]
    ey = ey[mask]

    meta: dict[str, object] = {"fit_points": [[float(a), float(b), float(c)] for a, b, c in zip(x, y, ey)]}
    xx = np.linspace(10.0, 35.0, 240)
    meta["x_grid"] = xx.tolist()

    try:
        popt, _ = curve_fit(
            pade11,
            x,
            y,
            sigma=ey,
            absolute_sigma=True,
            p0=[0.35, 0.03, 0.02],
            bounds=([-2.0, -1.0, -0.09], [2.0, 1.0, 0.2]),
            maxfev=200000,
        )
        resid = (y - pade11(x, *popt)) / ey
        meta["pade11"] = {
            "ok": True,
            "params": [float(v) for v in popt],
            "chi2": float(np.sum(resid * resid)),
            "ndf": int(len(x) - len(popt)),
            "y_grid": pade11(xx, *popt).clip(0.0, 1.2).tolist(),
        }
    except Exception as exc:  # pragma: no cover - diagnostic output
        meta["pade11"] = {"ok": False, "error": str(exc), "y_grid": [math.nan for _ in xx]}

    try:
        popt, _ = curve_fit(
            erf_alt,
            x,
            y,
            sigma=ey,
            absolute_sigma=True,
            p0=[0.35, 0.4, 18.0, 8.0],
            bounds=([0.0, 0.0, 5.0, 0.5], [1.2, 1.2, 40.0, 40.0]),
            maxfev=200000,
        )
        resid = (y - erf_alt(x, *popt)) / ey
        meta["erf_alt"] = {
            "ok": True,
            "params": [float(v) for v in popt],
            "chi2": float(np.sum(resid * resid)),
            "ndf": int(len(x) - len(popt)),
            "y_grid": erf_alt(xx, *popt).clip(0.0, 1.2).tolist(),
        }
    except Exception as exc:  # pragma: no cover - diagnostic output
        meta["erf_alt"] = {"ok": False, "error": str(exc), "y_grid": [math.nan for _ in xx]}

    return meta


def apply_common_style(ax) -> None:
    ax.set_xlim(10.0, 35.0)
    ax.set_ylim(0.0, 1.2)
    ax.minorticks_on()
    ax.tick_params(which="both", direction="in", top=True, right=True, labelsize=18, length=7)
    ax.tick_params(which="minor", length=4)
    ax.set_xlabel(r"$E_T^{\gamma}$ [GeV]", fontsize=22, loc="right")
    ax.set_ylabel("Purity", fontsize=23)
    for spine in ax.spines.values():
        spine.set_linewidth(1.2)


def save_auau_raw_panel(rows: list[dict[str, float]]) -> Path:
    out = OUT_DIR / "auau020_raw_vs_leakage_corrected_panel.png"
    x = np.array([r["pt_mid"] for r in rows], dtype=float)
    ex = np.array([r["pt_width"] / 2.0 for r in rows], dtype=float)
    raw = np.array([r["raw_purity"] for r in rows], dtype=float)
    raw_e = np.array([r["raw_purity_err"] for r in rows], dtype=float)
    corr = np.array([r["corrected_purity"] for r in rows], dtype=float)
    corr_e = np.array([r["corrected_purity_err"] for r in rows], dtype=float)

    plt.rcParams.update({"font.family": "DejaVu Sans", "mathtext.fontset": "dejavuserif"})
    fig, ax = plt.subplots(figsize=(6.4, 5.9), dpi=210)
    fig.patch.set_facecolor("white")
    apply_common_style(ax)
    ax.errorbar(
        x,
        corr,
        xerr=ex,
        yerr=corr_e,
        fmt="o",
        ms=7.5,
        mfc="white",
        mec="#144dff",
        ecolor="#144dff",
        elinewidth=1.15,
        capsize=2.5,
        label="w/ sig. leak. corr.",
        zorder=3,
    )
    ax.errorbar(
        x,
        raw,
        xerr=ex,
        yerr=raw_e,
        fmt="o",
        ms=7.0,
        mfc="black",
        mec="black",
        ecolor="black",
        elinewidth=1.15,
        capsize=2.5,
        label="w/o sig. leak. corr.",
        zorder=2,
    )
    ax.text(0.07, 0.95, r"$\bf{\it{sPHENIX}}$ Internal", transform=ax.transAxes, ha="left", va="top", fontsize=19)
    ax.text(0.07, 0.875, r"Au+Au $\sqrt{s_{NN}}=200$ GeV", transform=ax.transAxes, ha="left", va="top", fontsize=18)
    ax.text(0.07, 0.805, r"0--20%, default BDT", transform=ax.transAxes, ha="left", va="top", fontsize=18)
    ax.text(0.73, 0.95, r"$|\eta^\gamma|<0.7$", transform=ax.transAxes, ha="left", va="top", fontsize=19)
    ax.legend(loc="lower left", bbox_to_anchor=(0.05, 0.08), frameon=False, fontsize=17, handlelength=1.35)
    fig.subplots_adjust(left=0.16, right=0.985, top=0.985, bottom=0.16)
    fig.savefig(out, dpi=210)
    plt.close(fig)
    return out


def save_auau_fit_panel(rows: list[dict[str, float]], fit: dict[str, object]) -> Path:
    out = OUT_DIR / "auau020_leakage_corrected_fit_panel.png"
    x = np.array([r["pt_mid"] for r in rows], dtype=float)
    ex = np.array([r["pt_width"] / 2.0 for r in rows], dtype=float)
    corr = np.array([r["corrected_purity"] for r in rows], dtype=float)
    corr_e = np.array([r["corrected_purity_err"] for r in rows], dtype=float)
    xx = np.array(fit["x_grid"], dtype=float)

    fig, ax = plt.subplots(figsize=(6.4, 5.9), dpi=210)
    fig.patch.set_facecolor("white")
    apply_common_style(ax)
    ax.errorbar(
        x,
        corr,
        xerr=ex,
        yerr=corr_e,
        fmt="o",
        ms=6.6,
        mfc="black",
        mec="black",
        ecolor="black",
        elinewidth=1.15,
        capsize=2.5,
        label="Data purity (leakage-corr.)",
        zorder=3,
    )
    pade = fit.get("pade11", {})
    if isinstance(pade, dict) and pade.get("ok"):
        yy = np.array(pade["y_grid"], dtype=float)
        chi = float(pade["chi2"])
        ndf = int(pade["ndf"])
        ax.plot(xx, yy, color="#d62728", lw=2.0, label=rf"Padé[1/1] fit, $\chi^2$={chi:.1f}/{ndf}")
    erf_fit = fit.get("erf_alt", {})
    if isinstance(erf_fit, dict) and erf_fit.get("ok"):
        yy = np.array(erf_fit["y_grid"], dtype=float)
        chi = float(erf_fit["chi2"])
        ndf = int(erf_fit["ndf"])
        ax.plot(xx, yy, color="#0066cc", lw=2.0, ls=":", label=rf"Erf alt., $\chi^2$={chi:.1f}/{ndf}")
    ax.text(0.07, 0.95, r"$\bf{\it{sPHENIX}}$ Internal", transform=ax.transAxes, ha="left", va="top", fontsize=18)
    ax.text(0.07, 0.875, r"Au+Au 0--20%, $\sqrt{s_{NN}}=200$ GeV", transform=ax.transAxes, ha="left", va="top", fontsize=17)
    ax.text(0.07, 0.805, r"$|\eta^\gamma|<0.7$", transform=ax.transAxes, ha="left", va="top", fontsize=17)
    ax.legend(loc="lower right", bbox_to_anchor=(0.99, 0.055), frameon=False, fontsize=15.6, handlelength=2.0)
    fig.subplots_adjust(left=0.16, right=0.985, top=0.985, bottom=0.16)
    fig.savefig(out, dpi=210)
    plt.close(fig)
    return out


def fit_image(img: Image.Image, box: tuple[int, int], fill=(255, 255, 255)) -> Image.Image:
    w, h = box
    canvas = Image.new("RGBA", (w, h), fill + (255,))
    src = trim_white_margins(img).convert("RGBA")
    scale = min(w / src.width, h / src.height)
    scaled = src.resize((max(1, int(src.width * scale)), max(1, int(src.height * scale))), Image.Resampling.LANCZOS)
    canvas.alpha_composite(scaled, ((w - scaled.width) // 2, (h - scaled.height) // 2))
    return canvas


def compose_slide(auau_raw: Path, auau_fit: Path, fit_meta: dict[str, object], rows: list[dict[str, float]]) -> None:
    slide = Image.new("RGB", (2560, 1440), "white")
    draw = ImageDraw.Draw(slide)

    title_font = font(64, bold=True)
    subtitle_font = font(35)
    header_font = font(39, bold=True)
    note_font = font(25)

    title = "ABCD purity: pp reference and Au+Au 0-20%"
    draw.text((80, 38), title, font=title_font, fill=(10, 19, 35))
    subtitle = "Top row: raw sideband estimate and signal-leakage correction.  Bottom row: fits to leakage-corrected purity."
    draw.text((82, 112), subtitle, font=subtitle_font, fill=(38, 47, 65))

    col_w, row_h = 950, 555
    x_pp, x_auau = 205, 1400
    y_top, y_bot = 235, 820

    headers = [
        (x_pp, "p+p reference (PPG12 IAN Fig. 29)"),
        (x_auau, "Au+Au 0-20% (current default BDT output)"),
    ]
    for x, text in headers:
        bbox = draw.textbbox((0, 0), text, font=header_font)
        draw.text((x + (col_w - (bbox[2] - bbox[0])) // 2, 174), text, font=header_font, fill=(12, 25, 45))

    panels = [
        (PP_RAW_CROP, x_pp, y_top),
        (auau_raw, x_auau, y_top),
        (PP_FIT_CROP, x_pp, y_bot),
        (auau_fit, x_auau, y_bot),
    ]
    for path, x, y in panels:
        panel = fit_image(Image.open(path).convert("RGB"), (col_w, row_h))
        slide.paste(panel.convert("RGB"), (x, y), panel)

    note = "Au+Au bins are current-output bins overlapping 15 < E_T < 35 GeV; the fit row is a guide to current statistics."
    bbox = draw.textbbox((0, 0), note, font=note_font)
    draw.text(((2560 - (bbox[2] - bbox[0])) // 2, 1372), note, font=note_font, fill=(74, 85, 101))

    slide.save(OUT_PNG)

    manifest = {
        "schema": "THE85_FIG29_STYLE_PP_AUAU020_PURITY_2X2_V1",
        "output_png": str(OUT_PNG),
        "speaker_script": str(OUT_SCRIPT),
        "pp_reference": {
            "source_pdf": str(IAN_PDF),
            "page": 39,
            "figure": 29,
            "raw_corrected_crop": str(PP_RAW_CROP),
            "fit_crop": str(PP_FIT_CROP),
            "note": "Literal crops from the current PPG12 IAN figure 29.",
        },
        "auau_current": {
            "points_csv": str(POINTS_CSV),
            "selected_label": "Au+Au 0-20%",
            "raw_definition": "P_raw=max(0,A-B*C/D)/A from DATA ABCD counts.",
            "corrected_definition": "P_corr=S_A/A where S_A solves S_A=A-(B-fB*S_A)*(C-fC*S_A)/(D-fD*S_A).",
            "fit_note": "Fits are to the shown leakage-corrected AuAu 0-20% points; this is a current-statistics guide, not a final purity parameterization.",
            "points": rows,
            "fit_meta": fit_meta,
        },
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")

    OUT_SCRIPT.write_text(
        "# Speaker script: pp and Au+Au ABCD purity\n\n"
        "This slide puts the PPG12 Figure 29 pp purity reference next to the current Au+Au 0-20% result in the same visual language. "
        "The top row compares the raw ABCD sideband estimate to the signal-leakage-corrected purity. "
        "The bottom row shows the leakage-corrected points with fit forms. "
        "For Au+Au the fit should be read as a guide to the current 0-20% statistics, especially above about 24 GeV where the errors are large.\n",
        encoding="utf-8",
    )


def main() -> None:
    ensure_pp_crops()
    rows = read_auau020_points()
    fit_meta = fit_curves(rows)
    auau_raw = save_auau_raw_panel(rows)
    auau_fit = save_auau_fit_panel(rows, fit_meta)
    compose_slide(auau_raw, auau_fit, fit_meta, rows)
    print(OUT_PNG)
    print(OUT_MANIFEST)
    print(OUT_SCRIPT)


if __name__ == "__main__":
    main()
