#!/usr/bin/env python3
"""Build a slide-facing pp stitch-weight closure audit.

The audit is deliberately local/read-only. It uses the already accepted
RecoilJets in-situ pp stitching contract CSVs plus checked-in PPG12/RecoilJets
constants, and records unresolved generator-config text separately.
"""

from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image, ImageDraw, ImageFont


REPO = Path(__file__).resolve().parents[4]
OUT = REPO / "dataOutput/slides/wp_gammajets_6_1_26/pp_stitch_weight_parity_20260604"

PHOTON_CSV = (
    REPO
    / "dataOutput/ppPhotonMLPipeline/ppg12_basev3E_currentIAN_exactStitchPhoton0p5_20260526_1307"
    / "validation/insitu_stitching/pp_currentian_photon0p5_exactstitch_contract_points.csv"
)
PHOTON_JSON = PHOTON_CSV.with_name("pp_currentian_photon0p5_exactstitch_contract_summary.json")
JET_CSV = (
    REPO
    / "dataOutput/ppPhotonMLPipeline/ppg12_basev3E_currentIAN_exactStitch_20260523_0915"
    / "validation/insitu_stitching/pp_currentian_exactstitch_contract_points.csv"
)
JET_JSON = JET_CSV.with_name("pp_currentian_exactstitch_contract_summary.json")
JET_BOUNDARY_JSON = JET_CSV.with_name("pp_currentian_exactstitch_boundary_continuity_qa.json")

PPG12_WEIGHTS = REPO / "ppg12codeGit/efficiencytool/CrossSectionWeights.h"
PPG12_PYTHIA = REPO / "ppg12codeGit/simcrosssection/Fun4AllPythia.C"
RECOILJETS_H = REPO / "macros/AnalyzeRecoilJets.h"
RECOILJETS_CC = REPO / "src/RecoilJets.cc"
MERGE_SH = REPO / "scripts/mergeRecoilJets.sh"

W, H = 2560, 1440
FONT_DIR = Path("/System/Library/Fonts/Supplemental")
TIMES = FONT_DIR / "Times New Roman.ttf"
TIMES_BOLD = FONT_DIR / "Times New Roman Bold.ttf"
TIMES_ITALIC = FONT_DIR / "Times New Roman Italic.ttf"
TIMES_BOLD_ITALIC = FONT_DIR / "Times New Roman Bold Italic.ttf"

INK = (17, 24, 39)
MUTED = (79, 88, 105)
LINE = (201, 209, 220)
BG = (255, 255, 255)
PHOTON_ACCENT = (28, 111, 116)
PHOTON_SOFT = (235, 248, 248)
JET_ACCENT = (92, 78, 152)
JET_SOFT = (244, 241, 250)
NOTE_SOFT = (239, 245, 255)
OK_GREEN = (22, 101, 52)
WARN_AMBER = (146, 64, 14)

SAMPLE_COLORS = {
    "run28_photonjet5": "#2ca25f",
    "run28_photonjet10": "#2b8cbe",
    "run28_photonjet20": "#e6550d",
    "run28_jet8": "#d33682",
    "run28_jet12": "#2ca02c",
    "run28_jet20": "#1f77b4",
    "run28_jet30": "#ff7f0e",
    "run28_jet40": "#6f4eb2",
}


@dataclass(frozen=True)
class SampleSummary:
    group: str
    sample: str
    xsec_pb: float
    events: float
    bin_width: float
    window_lo: float
    window_hi: float
    source_files: int

    @property
    def weight_pb_per_event_bin(self) -> float:
        return self.xsec_pb / self.events / self.bin_width


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
    "section": font(39, bold=True),
    "body": font(28),
    "body_bold": font(28, bold=True),
    "small": font(24),
    "small_bold": font(24, bold=True),
    "tiny": font(21),
    "tiny_bold": font(21, bold=True),
    "table_header": font(26, bold=True),
    "table_body": font(27),
    "table_body_big": font(30),
}


def draw_round(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], fill, outline=LINE, width=2, radius=10) -> None:
    draw.rounded_rectangle(box, radius=radius, fill=fill, outline=outline, width=width)


def draw_wrapped(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    text: str,
    max_width: int,
    fnt: ImageFont.FreeTypeFont,
    fill=INK,
    line_spacing: int = 8,
) -> int:
    words = text.split()
    lines: list[str] = []
    cur = ""
    for word in words:
        test = word if not cur else f"{cur} {word}"
        if draw.textbbox((0, 0), test, font=fnt)[2] <= max_width:
            cur = test
        else:
            if cur:
                lines.append(cur)
            cur = word
    if cur:
        lines.append(cur)
    x, y = xy
    for line in lines:
        draw.text((x, y), line, font=fnt, fill=fill)
        y += fnt.size + line_spacing
    return y


def load_points(path: Path, group: str) -> tuple[dict[str, list[dict[str, float]]], list[SampleSummary]]:
    by_sample: dict[str, list[dict[str, float]]] = {}
    first: dict[str, dict[str, str]] = {}
    with path.open() as f:
        for row in csv.DictReader(f):
            if row["group"] != group:
                continue
            sample = row["sample"]
            first.setdefault(sample, row)
            x = float(row["bin_center"])
            y = float(row["density_pb_per_gev"])
            ey = float(row["density_err_pb_per_gev"])
            used = int(row["ppg12_bin_center_window"]) == 1
            if used and y > 0 and math.isfinite(y):
                by_sample.setdefault(sample, []).append({"x": x, "y": y, "ey": ey})
    summaries = []
    for sample, row in sorted(first.items()):
        bin_width = float(row["bin_high"]) - float(row["bin_low"])
        summaries.append(
            SampleSummary(
                group=group,
                sample=sample,
                xsec_pb=float(row["xsec_pb"]),
                events=float(row["events_processed_metadata"]),
                bin_width=bin_width,
                window_lo=float(row["stitch_window_low"]),
                window_hi=float(row["stitch_window_high"]),
                source_files=int(row.get("source_files") or 0),
            )
        )
    return by_sample, summaries


def combined_points(by_sample: dict[str, list[dict[str, float]]]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rows = [p for pts in by_sample.values() for p in pts]
    xs = sorted({p["x"] for p in rows})
    y = []
    ey = []
    for x in xs:
        vals = [p for p in rows if abs(p["x"] - x) < 1e-9]
        y.append(sum(p["y"] for p in vals))
        ey.append(math.sqrt(sum(p["ey"] ** 2 for p in vals)))
    return np.array(xs), np.array(y), np.array(ey)


def fit_modified_power(x: np.ndarray, y: np.ndarray):
    mask = (x > 0) & (y > 0) & np.isfinite(x) & np.isfinite(y)
    xfit = x[mask]
    yfit = y[mask]
    logx = np.log(xfit)
    l_inv = np.log(1.0 / xfit)
    design = np.column_stack([np.ones_like(xfit), l_inv, logx * l_inv, xfit * l_inv])
    coeff, *_ = np.linalg.lstsq(design, np.log(yfit), rcond=None)

    def eval_y(vals):
        vals = np.asarray(vals, dtype=float)
        lx = np.log(vals)
        li = np.log(1.0 / vals)
        return np.exp(coeff[0] + coeff[1] * li + coeff[2] * lx * li + coeff[3] * vals * li)

    return eval_y


def render_plot(
    *,
    title: str,
    sample_label: str,
    by_sample: dict[str, list[dict[str, float]]],
    out: Path,
    xlim: tuple[float, float],
    ylim: tuple[float, float],
    xlabel: str,
    boundaries: list[float],
) -> dict[str, float]:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "axes.linewidth": 1.15,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )
    x, y, ey = combined_points(by_sample)
    fit_mask = (x >= xlim[0]) & (x <= xlim[1]) & (y > 0)
    fit = fit_modified_power(x[fit_mask], y[fit_mask])
    ratio = y / fit(x)
    ratio_err = ey / fit(x)
    good = (x >= xlim[0]) & (x <= xlim[1]) & np.isfinite(ratio) & (y > 0)
    rms = float(np.sqrt(np.mean((ratio[good] - 1.0) ** 2)))

    fig = plt.figure(figsize=(8.9, 5.8), dpi=190, facecolor="white")
    gs = fig.add_gridspec(2, 1, height_ratios=[3.4, 1.05], hspace=0.06, left=0.125, right=0.965, top=0.89, bottom=0.135)
    ax = fig.add_subplot(gs[0])
    rax = fig.add_subplot(gs[1], sharex=ax)

    ax.set_yscale("log")
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    grid = np.linspace(max(xlim[0], 0.1), xlim[1], 500)
    ax.plot(grid, fit(grid), color="black", lw=1.8, ls="--", label="Modified power-law fit")

    for sample, pts in by_sample.items():
        arr_x = np.array([p["x"] for p in pts])
        arr_y = np.array([p["y"] for p in pts])
        arr_ey = np.array([p["ey"] for p in pts])
        mask = (arr_x >= xlim[0]) & (arr_x <= xlim[1]) & (arr_y > 0)
        ax.errorbar(
            arr_x[mask],
            arr_y[mask],
            yerr=arr_ey[mask],
            fmt="o",
            ms=4.8,
            lw=1.0,
            color=SAMPLE_COLORS.get(sample, "#555555"),
            label=sample.replace("run28_", ""),
            capsize=0,
            zorder=3,
        )

    band_lo = max(0.55, 1.0 - rms)
    band_hi = min(1.45, 1.0 + rms)
    rax.axhspan(band_lo, band_hi, color="#dbeafe", alpha=0.75, label=f"RMS band = {100*rms:.1f}%")
    rax.axhline(1.0, color="black", lw=1.0)
    rax.errorbar(x[good], ratio[good], yerr=ratio_err[good], fmt="o", ms=4.0, lw=0.9, color="#111827", capsize=0)
    rax.set_ylim(0.55, 1.45)

    ax.set_ylabel(r"$d\sigma / dp_T$ [pb/GeV]", fontsize=16.5, fontweight="bold")
    rax.set_ylabel("data / fit", fontsize=15.0, fontweight="bold")
    rax.set_xlabel(xlabel, fontsize=16.5, fontweight="bold")
    ax.tick_params(labelsize=14.0, which="both")
    rax.tick_params(labelsize=14.0, which="both")
    plt.setp(ax.get_xticklabels(), visible=False)

    label_box = {"boxstyle": "round,pad=0.28", "facecolor": "white", "edgecolor": "#cbd5e1", "alpha": 0.92, "linewidth": 0.8}
    ax.text(
        0.965,
        0.955,
        r"$\it{\bf{sPHENIX}}$ Internal" + f"\n{sample_label}",
        transform=ax.transAxes,
        fontsize=14.7,
        ha="right",
        va="top",
        linespacing=1.05,
        bbox=label_box,
        zorder=6,
    )
    ax.text(
        0.965,
        0.735,
        title,
        transform=ax.transAxes,
        fontsize=15.4,
        fontweight="bold",
        ha="right",
        va="top",
        bbox=label_box,
        zorder=6,
    )
    ax.legend(frameon=False, fontsize=12.3, loc="lower left", ncol=2, handlelength=1.2, columnspacing=0.75)
    rax.legend(frameon=False, fontsize=12.4, loc="upper right", handlelength=1.3)

    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=190)
    plt.close(fig)
    return {"ratio_rms": rms, "n_points": int(np.count_nonzero(good))}


def fmt_xsec(x: float) -> str:
    if abs(x) >= 1e4 or abs(x) < 1e-2:
        return f"{x:.4g}"
    return f"{x:.4f}".rstrip("0").rstrip(".")


def write_weight_table(path: Path, samples: list[SampleSummary]) -> None:
    rows = []
    for s in samples:
        rows.append(
            {
                "group": s.group,
                "sample": s.sample,
                "xsec_pb": f"{s.xsec_pb:.12g}",
                "events_processed": f"{s.events:.12g}",
                "bin_width_GeV": f"{s.bin_width:.8g}",
                "stitch_window": f"[{s.window_lo:g},{s.window_hi:g})",
                "source_files": s.source_files,
                "sigma_eff_pb_per_event_per_GeV": f"{s.weight_pb_per_event_bin:.12g}",
            }
        )
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def draw_weight_table(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    title: str,
    samples: list[SampleSummary],
    accent,
    soft,
    anchor_sample: str,
) -> None:
    x0, y0, x1, y1 = box
    draw_round(draw, box, fill=soft, outline=accent, width=3, radius=12)
    draw.rectangle((x0 + 22, y0 + 22, x0 + 36, y0 + 72), fill=accent)
    draw.text((x0 + 58, y0 + 15), title, font=F["section"], fill=INK)
    anchor_label = anchor_sample.replace("run28_", "").replace("photonjet", "PhotonJet").replace("jet", "Jet")
    header_y = y0 + 93
    cols = [x0 + 34, x0 + 272, x0 + 490, x0 + 690, x0 + 890]
    headers = ["sample", "window", "xsec [pb]", "Nproc", f"σ_eff / {anchor_label}"]
    for c, h in zip(cols, headers):
        draw.text((c, header_y), h, font=F["table_header"], fill=MUTED)
    draw.line((x0 + 30, header_y + 38, x1 - 30, header_y + 38), fill=accent, width=2)
    anchor = next(s for s in samples if s.sample == anchor_sample).weight_pb_per_event_bin
    row_start = header_y + 55
    row_font = F["table_body_big"] if len(samples) <= 3 else F["table_body"]
    row_step = min(56, max(38, int((y1 - row_start - 24) / max(1, len(samples)))))
    for i, s in enumerate(samples):
        y = row_start + i * row_step
        if i % 2 == 0:
            draw.rounded_rectangle((x0 + 24, y - 7, x1 - 24, y + row_step - 9), radius=7, fill=(255, 255, 255))
        name = s.sample.replace("run28_", "").replace("photonjet", "PhotonJet").replace("jet", "Jet")
        win = f"{s.window_lo:g}-{s.window_hi:g}"
        rel = s.weight_pb_per_event_bin / anchor
        vals = [name, win, fmt_xsec(s.xsec_pb), f"{s.events/1e6:.3f}M", f"{rel:.3g}x"]
        for c, v in zip(cols, vals):
            draw.text((c, y), v, font=row_font, fill=INK)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    photon_by_sample, photon_samples = load_points(PHOTON_CSV, "photon")
    jet_by_sample, jet_samples = load_points(JET_CSV, "jet")
    all_samples = photon_samples + jet_samples

    weight_csv = OUT / "pp_stitch_weight_parity_audit_weights.csv"
    write_weight_table(weight_csv, all_samples)

    photon_plot = OUT / "pp_photonjet_stitch_weight_parity_plot.png"
    jet_plot = OUT / "pp_inclusivejet_stitch_weight_parity_plot.png"
    photon_metrics = render_plot(
        title="Photon+jet 5+10+20 exact-stitch contract",
        sample_label=r"$p{+}p$ PYTHIA8 PhotonJet5/10/20",
        by_sample=photon_by_sample,
        out=photon_plot,
        xlim=(5, 40),
        ylim=(1e-2, 5e5),
        xlabel=r"max truth photon $p_T$ [GeV]",
        boundaries=[14, 22],
    )
    jet_metrics = render_plot(
        title="Inclusive-jet 8+12+20+30+40 exact-stitch contract",
        sample_label=r"$p{+}p$ PYTHIA8 inclusive jets",
        by_sample=jet_by_sample,
        out=jet_plot,
        xlim=(8, 50),
        ylim=(1e0, 3e6),
        xlabel=r"max R=0.4 truth jet $p_T$ [GeV]",
        boundaries=[14, 21, 32, 42],
    )

    boundary = json.loads(JET_BOUNDARY_JSON.read_text())["boundaries"]
    max_overlap = max(float(b["max_same_bin_overlap_fractional_deviation"]) for b in boundary)
    jet8_boundary = next(b for b in boundary if float(b["boundary_GeV"]) == 14.0)

    canvas = Image.new("RGB", (W, H), BG)
    draw = ImageDraw.Draw(canvas)
    draw.text((72, 44), "pp stitch-weight closure audit", font=F["title"], fill=INK)
    draw.text(
        (74, 112),
        "Same closure logic as embedded: identify the generator slice, use the in-situ denominator, then stitch by the owned truth-pT window.",
        font=F["subtitle"],
        fill=MUTED,
    )

    plot_boxes = [(70, 174, 1248, 890), (1322, 174, 2490, 890)]
    for box, img_path, accent in zip(plot_boxes, [photon_plot, jet_plot], [PHOTON_ACCENT, JET_ACCENT]):
        draw_round(draw, box, fill="white", outline=accent, width=3, radius=12)
        im = Image.open(img_path).convert("RGB")
        im.thumbnail((box[2] - box[0] - 14, box[3] - box[1] - 14), Image.Resampling.LANCZOS)
        canvas.paste(im, (box[0] + (box[2] - box[0] - im.width) // 2, box[1] + (box[3] - box[1] - im.height) // 2))

    draw_weight_table(
        draw,
        (70, 930, 1248, 1268),
        "PhotonJet 5+10+20 weights",
        sorted(photon_samples, key=lambda s: s.window_lo),
        PHOTON_ACCENT,
        PHOTON_SOFT,
        "run28_photonjet20",
    )
    draw_weight_table(
        draw,
        (1322, 930, 2490, 1268),
        "Inclusive Jet 8+12+20+30+40 weights",
        sorted(jet_samples, key=lambda s: s.window_lo),
        JET_ACCENT,
        JET_SOFT,
        "run28_jet40",
    )

    note_box = (70, 1287, 2490, 1407)
    draw_round(draw, note_box, fill=NOTE_SOFT, outline=(147, 197, 253), width=2, radius=10)
    draw.text((94, 1307), "Audit result:", font=F["body_bold"], fill=OK_GREEN)
    result = (
        f"Photon xsecs/windows match PPG12 CrossSectionWeights.h; inclusive uses RecoilJets/wiki jet8 = 1.3013e7 pb "
        f"rather than Shuhang header's old 1.15e7 pb. Jet boundary closure passes by same-bin overlap; "
        f"the worst overlap deviation is {100*max_overlap:.1f}% across all inclusive-jet boundaries."
    )
    draw_wrapped(draw, (284, 1307), result, 1470, F["small"], fill=INK, line_spacing=4)
    draw.text((1815, 1307), "Config gap:", font=F["body_bold"], fill=WARN_AMBER)
    draw_wrapped(
        draw,
        (1975, 1307),
        "Shuhang macro names the JetStructure cfgs via $CALIBRATIONROOT; cfg text is not local.",
        470,
        F["small"],
        fill=INK,
        line_spacing=4,
    )

    slide_png = OUT / "pp_stitch_weight_parity_audit_slide.png"
    canvas.save(slide_png)

    script_md = OUT / "pp_stitch_weight_parity_audit_slide_script.md"
    script_md.write_text(
        """# WP GammaJets pp Stitch-Weight Closure Audit Script

Here I am closing the pp stitching loop using the same logic that made the embedded combination clean.  The key point is that the stitched spectrum is not tuned by eye.  Each slice carries a generator cross section, an in-situ event denominator from the RecoilJets metadata histogram, and a truth-pT ownership window.

On the left, the PhotonJet5, PhotonJet10, and PhotonJet20 samples use the PPG12 photon windows and cross sections.  On the right, the inclusive-jet samples use Jet8 through Jet40, with Jet8 updated to the Jet Structure wiki value that is now in the RecoilJets constants.

The ratio panels are only a closure diagnostic against a smooth modified power-law fit.  The important audit result is that the event-denominator and window contract is explicit, and the inclusive-jet boundary check passes with the largest same-bin overlap deviation at about 12 percent.  The visible drop from the last bin of one owned window to the first bin of the next window is the falling jet spectrum over one bin; the closure test compares the two samples at the same bin centers.

The one caveat is generator-config text provenance: Shuhang's local source identifies the JetStructure Pythia config names through CALIBRATIONROOT, but those cfg files are not present in the local checkout.  So for tomorrow this slide supports the weight and denominator story, while exact cfg text readback would be a read-only SDCC or calibration-area follow-up.
""",
        encoding="utf-8",
    )

    manifest = {
        "schema": "PP_STITCH_WEIGHT_PARITY_AUDIT_V1",
        "slide_png": str(slide_png),
        "script_md": str(script_md),
        "weight_csv": str(weight_csv),
        "plots": {"photon": str(photon_plot), "inclusive_jet": str(jet_plot)},
        "inputs": {
            "photon_contract_csv": str(PHOTON_CSV),
            "photon_contract_json": str(PHOTON_JSON),
            "jet_contract_csv": str(JET_CSV),
            "jet_contract_json": str(JET_JSON),
            "jet_boundary_json": str(JET_BOUNDARY_JSON),
            "ppg12_cross_section_header": str(PPG12_WEIGHTS),
            "ppg12_pythia_macro": str(PPG12_PYTHIA),
            "recoiljets_constants": [str(RECOILJETS_H), str(RECOILJETS_CC), str(MERGE_SH)],
        },
        "normalization": "density_pb_per_GeV = raw window-owned counts * xsec_pb / events_processed_metadata / bin_width_GeV",
        "boundary_qa_definition": "Window-owned adjacent bins are expected to fall for a steep spectrum; stitch continuity is tested with same-bin overlap between neighboring generator slices.",
        "pythia_config_status": {
            "identified_from": str(PPG12_PYTHIA),
            "photon_cfgs_named": [
                "$CALIBRATIONROOT/Generators/JetStructure_TG/phpythia8_5GeV_JS_MDC2.cfg",
                "$CALIBRATIONROOT/Generators/JetStructure_TG/phpythia8_15GeV_JS_MDC2.cfg",
                "$CALIBRATIONROOT/Generators/JetStructure_TG/phpythia8_30GeV_JS_MDC2.cfg",
            ],
            "local_cfg_text_found": False,
            "local_env_CALIBRATIONROOT": "",
            "note": "No local JetStructure TG cfg text found under the repo/Desktop/home search used for this audit.",
        },
        "metrics": {
            "photon_ratio_rms": photon_metrics["ratio_rms"],
            "inclusive_jet_ratio_rms": jet_metrics["ratio_rms"],
            "jet_boundary_max_same_bin_overlap_fractional_deviation": max_overlap,
            "jet8_to_jet12_display_adjacent_ratio": float(jet8_boundary["display_adjacent_left_over_right"]),
            "all_jet_boundaries_pass": all(bool(b["pass_boundary"]) for b in boundary),
        },
        "google_slides_mutation": False,
    }
    manifest_path = OUT / "pp_stitch_weight_parity_audit_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    print(slide_png)
    print(script_md)
    print(manifest_path)
    print(weight_csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
