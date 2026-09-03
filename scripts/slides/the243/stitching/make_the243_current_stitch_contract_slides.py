#!/usr/bin/env python3
"""Render THE-243 current-production stitch contract slide candidates.

These slides audit the accepted schema-10 sample identity and exact bound
normalization source.  Au+Au embedded inclusive jets consume the canonical
``sigma_eff/Npass`` source-stitch artifact; pp and photon-signal rows retain
their accepted catalog contracts.  They do not infer a fine-binned spectral
closure from compact response products.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
import sys

import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator
from PIL import Image, ImageDraw, ImageFont


REPO = Path(__file__).resolve().parents[4]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.data_prep.recoiljets.auau_embedded_inclusive_schema10_weighting import (  # noqa: E402
    load_source_stitch_artifact,
)

W, H = 2560, 1440
FONT_DIR = Path("/System/Library/Fonts/Supplemental")
TIMES = FONT_DIR / "Times New Roman.ttf"
TIMES_BOLD = FONT_DIR / "Times New Roman Bold.ttf"
TIMES_ITALIC = FONT_DIR / "Times New Roman Italic.ttf"

INK = (17, 24, 39)
MUTED = (78, 88, 104)
LINE = (201, 210, 222)
TEAL = (24, 112, 116)
TEAL_SOFT = (235, 248, 248)
VIOLET = (91, 77, 151)
VIOLET_SOFT = (244, 241, 250)
BLUE_SOFT = (239, 246, 255)
GREEN = (24, 112, 62)
AMBER = (153, 78, 0)

COLORS = {
    "pp_photon5": "#269f67",
    "pp_photon10": "#2488bd",
    "pp_photon20": "#e56518",
    "pp_jet8": "#d33682",
    "pp_jet12": "#2ca02c",
    "pp_jet20": "#1f77b4",
    "pp_jet30": "#ff7f0e",
    "pp_jet40": "#6f4eb2",
    "auau_photon12": "#2369bd",
    "auau_photon20": "#d83b36",
    "auau_jet12": "#2369bd",
    "auau_jet20": "#ee7b1a",
    "auau_jet30": "#c53fa9",
    "auau_jet40": "#29954f",
}


@dataclass(frozen=True)
class Slice:
    sample_id: str
    lo: float
    hi: float | None
    xsec_pb: float
    normalization_denominator_events: int
    weight: float
    campaign: str


WINDOWS = {
    "pp": {
        "photon": [("pp_photon5", 5, 14), ("pp_photon10", 14, 22), ("pp_photon20", 22, 40)],
        "jet": [
            ("pp_jet8", 9, 14),
            ("pp_jet12", 14, 21),
            ("pp_jet20", 21, 32),
            ("pp_jet30", 32, 42),
            ("pp_jet40", 42, 50),
        ],
    },
    "auau": {
        "photon": [("auau_photon12", 12, 21), ("auau_photon20", 21, 40)],
        "jet": [
            ("auau_jet12", 12, 21),
            ("auau_jet20", 21, 31),
            ("auau_jet30", 31, 41),
            ("auau_jet40", 41, 50),
        ],
    },
}


def font(size: int, *, bold: bool = False, italic: bool = False) -> ImageFont.FreeTypeFont:
    path = TIMES_BOLD if bold else TIMES_ITALIC if italic else TIMES
    try:
        return ImageFont.truetype(str(path), size)
    except OSError:
        return ImageFont.load_default()


F = {
    "title": font(60, bold=True),
    "subtitle": font(30),
    "panel": font(37, bold=True),
    "body": font(28),
    "body_bold": font(28, bold=True),
    "small": font(24),
    "small_bold": font(24, bold=True),
    "tiny": font(21),
    "tiny_bold": font(21, bold=True),
    "foot": font(25),
    "foot_bold": font(25, bold=True),
}


def rounded(draw: ImageDraw.ImageDraw, box, fill, outline=LINE, width=2, radius=12):
    draw.rounded_rectangle(box, radius=radius, fill=fill, outline=outline, width=width)


def load_slices(
    system: str,
    group: str,
    *,
    catalog_path: Path,
    source_stitch_assembly: Path,
    source_stitch_receipt: Path,
) -> list[Slice]:
    payload = json.loads(catalog_path.read_text())
    samples = payload["samples"]
    source_stitch = (
        load_source_stitch_artifact(source_stitch_assembly, source_stitch_receipt)
        if system == "auau" and group == "jet"
        else None
    )
    out: list[Slice] = []
    for sample_id, lo, hi in WINDOWS[system][group]:
        row = samples[sample_id]
        if source_stitch is not None:
            canonical = source_stitch.sample(sample_id)
            row_xsec = canonical.ownership_effective_cross_section_pb
            row_events = canonical.normalization_denominator_events
            row_weight = canonical.stitching_weight_pb_per_owned_event
            row_lo = canonical.ownership_low_gev
            row_hi = canonical.ownership_high_gev
        else:
            row_xsec = float(row["cross_section_pb"])
            row_events = int(row["generated_events"])
            row_weight = float(row["cross_section_weight_pb_per_event"])
            row_lo = float(lo)
            row_hi = float(hi)
        out.append(
            Slice(
                sample_id=sample_id,
                lo=row_lo,
                hi=row_hi,
                xsec_pb=row_xsec,
                normalization_denominator_events=row_events,
                weight=row_weight,
                campaign=str(row["production_campaign_tag"]),
            )
        )
    return out


def render_contract_plot(slices: list[Slice], *, title: str, out: Path) -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "axes.linewidth": 1.2,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )
    fig, ax = plt.subplots(figsize=(9.4, 4.65), dpi=180, facecolor="white")
    fig.subplots_adjust(left=0.13, right=0.97, bottom=0.19, top=0.82)
    ref = slices[-1].weight
    for row in slices:
        display_hi = 50.0 if row.hi is None else row.hi
        y = row.weight / ref
        c = COLORS[row.sample_id]
        ax.plot([row.lo, display_hi], [y, y], color=c, lw=9, solid_capstyle="butt")
        ax.scatter([0.5 * (row.lo + display_hi)], [y], s=70, color=c, edgecolor="white", lw=1.0, zorder=3)
        ax.text(
            0.5 * (row.lo + display_hi),
            y * 1.22,
            row.sample_id.replace("pp_", "").replace("auau_", ""),
            color=c,
            ha="center",
            va="bottom",
            fontsize=13.5,
            fontweight="bold",
        )
    for row in slices[:-1]:
        assert row.hi is not None
        ax.axvline(row.hi, color="#aeb8c6", lw=1.0, ls=(0, (3, 3)), zorder=0)
    ax.set_yscale("log")
    ax.yaxis.set_major_locator(LogLocator(base=10))
    ax.set_xlim(min(r.lo for r in slices) - 1, max(50.0 if r.hi is None else r.hi for r in slices) + 1)
    ymin = min(r.weight / ref for r in slices) / 2.8
    ymax = max(r.weight / ref for r in slices) * 4.5
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel("owned generator truth-$p_T$ window [GeV]", fontsize=16)
    ax.set_ylabel("per-event weight / highest-slice weight", fontsize=16)
    ax.grid(axis="y", which="major", color="#d9e0ea", lw=0.8, alpha=0.8)
    ax.set_title(title, loc="left", fontsize=19, fontweight="bold", pad=14)
    ax.text(
        0.985,
        1.075,
        r"$\it{\bf{sPHENIX}}$ Internal",
        transform=ax.transAxes,
        ha="right",
        va="center",
        fontsize=15,
    )
    ax.text(
        0.985,
        0.97,
        "exact current catalog weights\n(no spectral-shape inference)",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=12.5,
        color="#4f5869",
        bbox={"boxstyle": "round,pad=0.3", "fc": "white", "ec": "#ccd4df", "alpha": 0.94},
    )
    fig.savefig(out, dpi=180, facecolor="white")
    plt.close(fig)


def paste_fit(canvas: Image.Image, panel: Image.Image, box: tuple[int, int, int, int]) -> None:
    x0, y0, x1, y1 = box
    panel = panel.copy()
    panel.thumbnail((x1 - x0, y1 - y0), Image.Resampling.LANCZOS)
    x = x0 + (x1 - x0 - panel.width) // 2
    y = y0 + (y1 - y0 - panel.height) // 2
    canvas.paste(panel, (x, y))


def fmt_xsec(value: float) -> str:
    if value >= 1.0e5 or value < 0.01:
        return f"{value:.3e}"
    if value >= 100:
        return f"{value:,.1f}"
    return f"{value:.4g}"


def draw_table(draw: ImageDraw.ImageDraw, box, slices: list[Slice], *, accent, soft) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, soft, accent, 2)
    draw.rectangle((x0 + 20, y0 + 18, x0 + 34, y0 + 70), fill=accent)
    label = " + ".join(r.sample_id.split("_")[-1].replace("photon", "Photon").replace("jet", "Jet") for r in slices)
    draw.text((x0 + 52, y0 + 18), f"{label} weights", font=F["panel"], fill=INK)
    uses_npass = any(row.sample_id.startswith("auau_jet") for row in slices)
    headers = ["sample", "owned window", "xsec [pb]", "Npass" if uses_npass else "Ngen", "weight [pb/event]"]
    widths = [0.18, 0.19, 0.19, 0.17, 0.27]
    inner_x0, inner_x1 = x0 + 26, x1 - 26
    xs = [inner_x0]
    for frac in widths[:-1]:
        xs.append(xs[-1] + int((inner_x1 - inner_x0) * frac))
    hy = y0 + 94
    for x, h in zip(xs, headers):
        draw.text((x, hy), h, font=F["tiny_bold"], fill=MUTED)
    draw.line((inner_x0, hy + 33, inner_x1, hy + 33), fill=accent, width=2)
    row_y = hy + 44
    row_h = 48
    for idx, row in enumerate(slices):
        if idx % 2 == 0:
            draw.rounded_rectangle((inner_x0 - 8, row_y - 5, inner_x1 + 4, row_y + 37), 7, fill=(255, 255, 255))
        vals = [
            row.sample_id.replace("pp_", "").replace("auau_", ""),
            f">={row.lo:g}" if row.hi is None else f"{row.lo:g}-{row.hi:g}",
            fmt_xsec(row.xsec_pb),
            f"{row.normalization_denominator_events/1e6:.3f}M",
            f"{row.weight:.4g}",
        ]
        for x, value in zip(xs, vals):
            draw.text((x, row_y), value, font=F["small"], fill=INK)
        row_y += row_h


def render_slide(
    system: str,
    *,
    catalog_path: Path,
    output_dir: Path,
    source_stitch_assembly: Path,
    source_stitch_receipt: Path,
) -> Path:
    photon = load_slices(
        system,
        "photon",
        catalog_path=catalog_path,
        source_stitch_assembly=source_stitch_assembly,
        source_stitch_receipt=source_stitch_receipt,
    )
    jets = load_slices(
        system,
        "jet",
        catalog_path=catalog_path,
        source_stitch_assembly=source_stitch_assembly,
        source_stitch_receipt=source_stitch_receipt,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    tag = "pp" if system == "pp" else "auau"
    pplot = output_dir / f"current_{tag}_photon_stitch_contract_plot.png"
    jplot = output_dir / f"current_{tag}_jet_stitch_contract_plot.png"
    render_contract_plot(photon, title="Photon-triggered generator slices", out=pplot)
    render_contract_plot(jets, title="Inclusive-jet generator slices", out=jplot)

    canvas = Image.new("RGB", (W, H), "white")
    draw = ImageDraw.Draw(canvas)
    title = "pp stitch-weight contract — current THE-121/122 production" if system == "pp" else "embedded stitch-weight contract — current THE-121/122 production"
    subtitle = (
        "Exact sample normalization and ownership used by the schema-10 response inputs feeding THE-243."
    )
    draw.text((74, 42), title, font=F["title"], fill=INK)
    draw.text((76, 112), subtitle, font=F["subtitle"], fill=MUTED)

    left = (72, 174, 1263, 858)
    right = (1325, 174, 2488, 858)
    rounded(draw, left, (255, 255, 255), TEAL, 3)
    rounded(draw, right, (255, 255, 255), VIOLET, 3)
    paste_fit(canvas, Image.open(pplot), (92, 196, 1243, 840))
    paste_fit(canvas, Image.open(jplot), (1345, 196, 2468, 840))

    draw_table(draw, (72, 900, 1263, 1282), photon, accent=TEAL, soft=TEAL_SOFT)
    draw_table(draw, (1325, 900, 2488, 1282), jets, accent=VIOLET, soft=VIOLET_SOFT)

    rounded(draw, (72, 1310, 2488, 1406), BLUE_SOFT, (155, 194, 246), 2)
    draw.text((98, 1328), "Certified now:", font=F["foot_bold"], fill=GREEN)
    draw.text(
        (268, 1328),
        "14 isolated sample contracts · 20,006 exact schema-10 inputs · cross sections, denominators, weights and source provenance bound.",
        font=F["foot"],
        fill=INK,
    )
    draw.text((98, 1362), "Open display item:", font=F["foot_bold"], fill=AMBER)
    draw.text(
        (298, 1362),
        "fine-binned max-truth-pT spectral closure was not retained in the compact response; this slide does not infer it from coarse bins.",
        font=F["foot"],
        fill=INK,
    )

    out = output_dir / f"S05_{tag}_current_the121_the122_stitch_contract_candidate.png"
    canvas.save(out)
    manifest = {
        "schema": "THE243CurrentStitchContractSlideCandidateV1",
        "system": system,
        "status": "CURRENT_PRODUCTION_CONTRACT__NOT_FINE_BINNED_SPECTRAL_CLOSURE",
        "catalog_path": str(catalog_path),
        "catalog_sample_count": 14,
        "catalog_input_root_count": 20006,
        "auau_embedded_inclusive_source_stitch": (
            {
                "assembly_path": str(source_stitch_assembly),
                "assembly_sha256": load_source_stitch_artifact(
                    source_stitch_assembly, source_stitch_receipt
                ).assembly_sha256,
                "receipt_path": str(source_stitch_receipt),
            }
            if system == "auau"
            else None
        ),
        "photon_samples": [r.__dict__ for r in photon],
        "jet_samples": [r.__dict__ for r in jets],
        "output_png": str(out),
    }
    (output_dir / f"S05_{tag}_current_the121_the122_stitch_contract_candidate_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--source-stitch-assembly", type=Path, required=True)
    parser.add_argument("--source-stitch-receipt", type=Path, required=True)
    args = parser.parse_args()
    for required in (
        args.catalog,
        args.source_stitch_assembly,
        args.source_stitch_receipt,
    ):
        if not required.is_file():
            raise FileNotFoundError(required)
    for system in ("pp", "auau"):
        print(
            render_slide(
                system,
                catalog_path=args.catalog,
                output_dir=args.output_dir,
                source_stitch_assembly=args.source_stitch_assembly,
                source_stitch_receipt=args.source_stitch_receipt,
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
