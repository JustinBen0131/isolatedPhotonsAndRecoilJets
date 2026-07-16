#!/usr/bin/env python3
"""Build the THE-88 nominal AuAu signal-leakage centrality slide."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from matplotlib.patches import Patch


THIS_FILE = Path(__file__).resolve()
REPO = next((p for p in THIS_FILE.parents if (p / "AGENTS.md").exists()), THIS_FILE.parents[4])
SCRIPTS = REPO / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.append(str(SCRIPTS))

from slides.common.slide_defaults import (  # noqa: E402
    SLIDE_DPI,
    SLIDE_HEIGHT_PX,
    SLIDE_WIDTH_PX,
    slide_figsize,
)


AUAU_ROOT = (
    REPO
    / "InputFiles/the88_bounded_nontight_20260708/simembedded"
    / "preselectionNewPPG12_tightAuAuCentInputBase3x3BDT_nonTightAuAuBDTSideband_baseVariant"
    / "photonJet12and20merged_SIM/RecoilJets_embeddedPhoton12plus20_MERGED.root"
)
OUT_DIR = (
    REPO
    / "dataOutput/auau/the88_auau_data_centrality_fixed_20260713"
    / "signal_leakage_auau_centrality"
)
STEM = "the88_auau_signal_leakage_bcd_r40_cent_sliding_by_centrality_slide"
OUT_PNG = OUT_DIR / f"{STEM}.png"
OUT_CSV = OUT_DIR / f"{STEM}.csv"
OUT_JSON = OUT_DIR / f"{STEM}.manifest.json"
OUT_LAYOUT = OUT_DIR / f"{STEM}.layout_nodes.json"
OUT_AUDIT = OUT_DIR / f"{STEM}.audit.json"
OUT_SCRIPT = OUT_DIR / f"{STEM}_speaker_script.md"

TOPDIR = "SIM"
AUAU_ISO_TAG = "isoR40_isSliding"

AUAU_BINS = [(15, 17), (17, 19), (19, 21), (21, 23), (23, 26), (26, 35)]

SAMPLES = [
    ("auau_0_20", "Au+Au 0-20%", "#0284c7", "0_20"),
    ("auau_20_50", "Au+Au 20-50%", "#f59e0b", "20_50"),
    ("auau_50_80", "Au+Au 50-80%", "#15803d", "50_80"),
]

REGIONS = [
    ("B", 2, r"$f_B=N_B^{sig}/N_A^{sig}$", "Region B: tight, non-isolated"),
    ("C", 3, r"$f_C=N_C^{sig}/N_A^{sig}$", "Region C: non-tight, isolated"),
    ("D", 4, r"$f_D=N_D^{sig}/N_A^{sig}$", "Region D: non-tight, non-isolated"),
]

TIMES_DIR = Path("/System/Library/Fonts/Supplemental")
TIMES_FONTS = [
    TIMES_DIR / "Times New Roman.ttf",
    TIMES_DIR / "Times New Roman Bold.ttf",
    TIMES_DIR / "Times New Roman Italic.ttf",
    TIMES_DIR / "Times New Roman Bold Italic.ttf",
]

INK = "#111827"
MUTED = "#475569"
GRID = "#dbe4ef"
PANEL_FILL = "#fbfdff"
EDGE = "#c8d3e2"


@dataclass(frozen=True)
class Count:
    value: float
    error: float


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def configure_fonts() -> None:
    for font in TIMES_FONTS:
        if font.exists():
            font_manager.fontManager.addfont(str(font))
    plt.rcParams.update(
        {
            "font.family": "Times New Roman",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "custom",
            "mathtext.rm": "Times New Roman",
            "mathtext.it": "Times New Roman:italic",
            "mathtext.bf": "Times New Roman:bold",
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "axes.unicode_minus": False,
        }
    )


def open_root(path: Path):
    import ROOT

    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.AddDirectory(False)
    root_file = ROOT.TFile.Open(str(path), "READ")
    if not root_file or root_file.IsZombie() or root_file.TestBit(ROOT.TFile.kRecovered):
        raise OSError(f"Could not open clean ROOT file: {path}")
    return root_file


def get_count(hist, bin_index: int) -> Count:
    value = float(hist.GetBinContent(bin_index))
    error = float(hist.GetBinError(bin_index))
    if error <= 0.0 and value > 0.0:
        error = math.sqrt(value)
    return Count(value=value, error=error)


def ratio(num: Count, den: Count) -> tuple[float, float]:
    if den.value <= 0.0:
        raise ValueError("Region-A truth-signal denominator is not positive")
    value = num.value / den.value
    if num.value <= 0.0:
        return value, num.error / den.value
    rel2 = (num.error / num.value) ** 2 + (den.error / den.value) ** 2
    return value, abs(value) * math.sqrt(max(0.0, rel2))


def hist_path(cent_key: str, lo: int, hi: int) -> str:
    return f"{TOPDIR}/h_sigABCD_MC_{AUAU_ISO_TAG}_pT_{lo}_{hi}_cent_{cent_key}"


def collect_rows() -> list[dict]:
    auau_file = open_root(AUAU_ROOT)
    rows: list[dict] = []
    try:
        for pt_index, (lo, hi) in enumerate(AUAU_BINS):
            for sample_key, sample_label, _color, cent_key in SAMPLES:
                path = hist_path(cent_key, lo, hi)
                hist = auau_file.Get(path)
                if not hist:
                    raise KeyError(f"Missing histogram: {path}")
                a_count = get_count(hist, 1)
                for region, bin_index, ratio_label, description in REGIONS:
                    side_count = get_count(hist, bin_index)
                    value, error = ratio(side_count, a_count)
                    rows.append(
                        {
                            "pt_index": pt_index,
                            "sample_key": sample_key,
                            "sample_label": sample_label,
                            "centrality": cent_key,
                            "pt_lo_gev": lo,
                            "pt_hi_gev": hi,
                            "region": region,
                            "region_description": description,
                            "ratio_label": ratio_label,
                            "n_a_signal": a_count.value,
                            "n_side_signal": side_count.value,
                            "leakage_fraction": value,
                            "leakage_fraction_error": error,
                            "histogram": path,
                        }
                    )
    finally:
        auau_file.Close()
    return rows


def matrix_for(rows: list[dict], region: str) -> tuple[np.ndarray, np.ndarray]:
    values = np.full((len(SAMPLES), len(AUAU_BINS)), np.nan, dtype=float)
    errors = np.full_like(values, np.nan)
    lookup = {(row["sample_key"], row["pt_index"], row["region"]): row for row in rows}
    for sample_index, (sample_key, _label, _color, _cent) in enumerate(SAMPLES):
        for pt_index in range(len(AUAU_BINS)):
            row = lookup[(sample_key, pt_index, region)]
            values[sample_index, pt_index] = row["leakage_fraction"]
            errors[sample_index, pt_index] = row["leakage_fraction_error"]
    return values, errors


def style_axis(ax: plt.Axes, title: str, ylabel: str, show_x: bool) -> None:
    ax.set_facecolor(PANEL_FILL)
    ax.set_title(title, loc="left", fontsize=19.5, fontweight="bold", pad=8, color=INK)
    ax.set_ylabel(ylabel, fontsize=16.3, color=INK, labelpad=10)
    ax.grid(axis="y", color=GRID, linewidth=1.15, zorder=1)
    ax.tick_params(axis="y", labelsize=15.0, colors=INK, length=5, width=1.0)
    ax.tick_params(axis="x", labelsize=14.5, colors=INK, length=0, pad=7, labelbottom=show_x)
    for spine in ax.spines.values():
        spine.set_edgecolor(EDGE)
        spine.set_linewidth(1.1)


def draw_bars(ax: plt.Axes, rows: list[dict], region: str, title: str, ylabel: str, show_x: bool) -> None:
    values, errors = matrix_for(rows, region)
    style_axis(ax, title, ylabel, show_x)
    x = np.arange(len(AUAU_BINS), dtype=float)
    width = 0.22
    offsets = (np.arange(len(SAMPLES)) - 1.0) * width

    for pt_index in range(len(AUAU_BINS)):
        if pt_index % 2 == 0:
            ax.axvspan(pt_index - 0.5, pt_index + 0.5, color="#eef4fb", alpha=0.68, zorder=0)

    for sample_index, (_key, label, color, _cent) in enumerate(SAMPLES):
        ax.bar(
            x + offsets[sample_index],
            values[sample_index],
            width=width * 0.92,
            color=color,
            edgecolor="white",
            linewidth=0.8,
            label=label,
            zorder=3,
        )
        ax.errorbar(
            x + offsets[sample_index],
            values[sample_index],
            yerr=errors[sample_index],
            fmt="none",
            ecolor="#0f172a",
            elinewidth=0.8,
            capsize=2.2,
            zorder=4,
        )

    top = float(np.nanmax(values + errors)) * 1.25
    floors = {"B": 0.15, "C": 0.35, "D": 0.05}
    ax.set_ylim(0.0, max(top, floors[region]))
    ax.set_xlim(-0.55, len(AUAU_BINS) - 0.45)
    ax.set_xticks(x)
    labels = [f"{lo}-{hi}" for lo, hi in AUAU_BINS]
    ax.set_xticklabels(labels, fontweight="bold")
    if show_x:
        ax.set_xlabel(
            r"cluster $E_T$ interval [GeV]",
            fontsize=17.5,
            fontweight="bold",
            labelpad=9,
        )


def write_csv(rows: list[dict]) -> None:
    with OUT_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_layout() -> None:
    payload = {
        "schema": "slide_layout_nodes_v1",
        "slide_size": [SLIDE_WIDTH_PX, SLIDE_HEIGHT_PX],
        "title_axis_x": 126,
        "minimum_title_font_px": 72,
        "minimum_audience_font_px": 31,
        "nodes": [
            {
                "kind": "text",
                "name": "title",
                "role": "title",
                "title_anchor": True,
                "text": "Truth-signal sideband leakage by Au+Au centrality",
                "bbox": [126, 62, 2380, 145],
                "font_px": 100,
            },
            {
                "kind": "text",
                "name": "sample definition",
                "role": "audience",
                "title_axis_align": "left",
                "text": "Embedded photon+jet 12 + 20; historical nominal tight and bounded non-tight BDT definitions",
                "bbox": [126, 180, 2400, 224],
                "font_px": 34,
            },
            {
                "kind": "text",
                "name": "isolation definition",
                "role": "audience",
                "title_axis_align": "left",
                "text": "Nominal reco R=0.4, Eiso < 7.57 - 0.0658*c GeV, c = event centrality percentile",
                "bbox": [126, 236, 2400, 278],
                "font_px": 31,
            },
            {
                "kind": "panel",
                "name": "region B panel",
                "bbox": [256, 360, 2420, 605],
            },
            {
                "kind": "panel",
                "name": "region C panel",
                "bbox": [256, 675, 2420, 920],
            },
            {
                "kind": "panel",
                "name": "region D panel",
                "bbox": [256, 990, 2420, 1235],
            },
        ],
    }
    OUT_LAYOUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def write_speaker_script(rows: list[dict]) -> None:
    region_ranges: dict[str, dict[str, tuple[float, float]]] = {}
    for region, _bin_index, _ratio, _description in REGIONS:
        region_ranges[region] = {}
        for sample_key, sample_label, _color, _cent in SAMPLES:
            values = [
                row["leakage_fraction"]
                for row in rows
                if row["region"] == region and row["sample_key"] == sample_key
            ]
            region_ranges[region][sample_label] = (min(values), max(values))

    all_b = [value for value_range in region_ranges["B"].values() for value in value_range]
    all_c = [value for value_range in region_ranges["C"].values() for value in value_range]
    all_d = [value for value_range in region_ranges["D"].values() for value in value_range]

    text = f"""# Speaker script

This slide compares truth-signal contamination in the three ABCD sideband regions with the signal in Region A. The bars show N B, N C, or N D divided by N A, so lower values mean less prompt-signal leakage into that sideband.

The input is the merged Au+Au embedded photon+jet 12 and 20 GeV production. The three colors show 0-20, 20-50, and 50-80 percent centrality. All three use the nominal R equals 0.4 reconstructed isolation requirement, with E T iso less than 7.57 minus 0.0658 times the event centrality percentile.

The x axis now uses only the native Au+Au cluster transverse-energy intervals, from 15 to 35 GeV. There is no pp comparison and no mixed binning on this slide.

Across the displayed intervals, Region B leakage spans {min(all_b):.3f} to {max(all_b):.3f}, Region C spans {min(all_c):.3f} to {max(all_c):.3f}, and Region D spans {min(all_d):.3f} to {max(all_d):.3f}. Region C is the largest prompt-leakage term, while Region D remains the smallest. The error bars are statistical uncertainties from the stored weighted histograms, and each panel uses its own labeled y-axis scale.

The photon identification uses the historical pre-MinBiasClassifier AuAuCentInputBase3x3BDT tight working point and the bounded AuAuBDTSideband non-tight definition. That training provenance is the only exception relative to the future retrained embedded baseline.
"""
    OUT_SCRIPT.write_text(text, encoding="utf-8")


def write_manifest(rows: list[dict]) -> None:
    manifest = {
        "schema_version": 1,
        "artifact_type": "audience_facing_full_slide_png",
        "output_png": str(OUT_PNG),
        "csv": str(OUT_CSV),
        "layout_nodes": str(OUT_LAYOUT),
        "audit": str(OUT_AUDIT),
        "speaker_script": str(OUT_SCRIPT),
        "sources": {
            "auau": {
                "path": str(AUAU_ROOT),
                "sha256": sha256(AUAU_ROOT),
                "samples": ["run28_embeddedPhoton12", "run28_embeddedPhoton20"],
                "centrality_bins_percent": [[0, 20], [20, 50], [50, 80]],
                "isolation": "nominal reco R=0.4 with Eiso < 7.57 - 0.0658*centrality_percentile GeV",
                "histogram_family": f"SIM/h_sigABCD_MC_{AUAU_ISO_TAG}_pT_<lo>_<hi>_cent_<cent>",
                "embedded_minbias_classifier_required": False,
            },
        },
        "leakage_definition": "f_X = N_sig(X) / N_sig(A), X in {B,C,D}",
        "abcd_regions": {
            "A": "tight and isolated",
            "B": "tight and non-isolated",
            "C": "non-tight and isolated",
            "D": "non-tight and non-isolated",
        },
        "native_auau_bins_gev": AUAU_BINS,
        "binning_note": "Only native AuAu intervals are shown; no pp comparison, fractional rebinning, or interpolation.",
        "photon_id_contract": {
            "tight": "AuAuCentInputBase3x3BDT centlinear threshold 0.5387310379 + 0.0011102647*centrality",
            "nontight": "AuAuBDTSideband bounded relative-to-tight offsets [-0.20,-0.03]",
        },
        "status_note": "Nominal AuAu reconstructed isolation and historical pre-MinBiasClassifier BDT selection.",
        "row_count": len(rows),
        "png_dimensions": [SLIDE_WIDTH_PX, SLIDE_HEIGHT_PX],
        "google_slides_mutated": False,
    }
    OUT_JSON.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def draw_slide(rows: list[dict]) -> None:
    configure_fonts()
    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI, facecolor="white")

    fig.text(
        0.049,
        0.954,
        "Truth-signal sideband leakage by Au+Au centrality",
        ha="left",
        va="top",
        fontsize=34.5,
        fontweight="bold",
        color=INK,
    )
    fig.text(
        0.049,
        0.875,
        r"Embedded photon+jet 12 + 20; historical AuAuCentInputBase3x3BDT tight + bounded AuAuBDTSideband non-tight",
        ha="left",
        va="top",
        fontsize=17.5,
        color=MUTED,
    )
    fig.text(
        0.049,
        0.837,
        r"Nominal reco isolation: $R=0.4$, $E_T^{iso}<7.57-0.0658c$ GeV, $c=$ event centrality percentile",
        ha="left",
        va="top",
        fontsize=15.7,
        color=MUTED,
    )

    positions = [
        [0.100, 0.580, 0.845, 0.170],
        [0.100, 0.360, 0.845, 0.170],
        [0.100, 0.140, 0.845, 0.170],
    ]
    axes: list[plt.Axes] = []
    for index, ((region, _bin, ylabel, description), position) in enumerate(zip(REGIONS, positions)):
        ax = fig.add_axes(position)
        draw_bars(ax, rows, region, description, ylabel, show_x=(index == 2))
        axes.append(ax)

    axes[0].text(
        0.012,
        0.955,
        r"$\bf{\it{sPHENIX}}$ Internal",
        transform=axes[0].transAxes,
        ha="left",
        va="top",
        fontsize=17.5,
        color=INK,
    )
    handles = [Patch(facecolor=color, edgecolor="none", label=label) for _key, label, color, _cent in SAMPLES]
    fig.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.70, 0.805),
        frameon=False,
        ncol=3,
        fontsize=14.5,
        columnspacing=1.25,
        handlelength=1.40,
        handletextpad=0.45,
    )

    fig.savefig(OUT_PNG, dpi=SLIDE_DPI, facecolor="white")
    plt.close(fig)


def run_audit() -> None:
    command = [
        sys.executable,
        str(REPO / "scripts/slides/common/post_render_slide_audit.py"),
        "--png",
        str(OUT_PNG),
        "--layout-nodes",
        str(OUT_LAYOUT),
        "--output",
        str(OUT_AUDIT),
        "--require-pass",
    ]
    subprocess.run(command, cwd=REPO, check=True)


def main() -> int:
    if not AUAU_ROOT.is_file():
        raise FileNotFoundError(AUAU_ROOT)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = collect_rows()
    if len(rows) != len(AUAU_BINS) * len(SAMPLES) * len(REGIONS):
        raise RuntimeError(f"Unexpected output row count: {len(rows)}")
    if not all(math.isfinite(row["leakage_fraction"]) and row["leakage_fraction"] >= 0.0 for row in rows):
        raise RuntimeError("Non-finite or negative leakage fraction found")
    write_csv(rows)
    draw_slide(rows)
    write_layout()
    write_speaker_script(rows)
    write_manifest(rows)
    run_audit()
    for path in (OUT_PNG, OUT_CSV, OUT_JSON, OUT_LAYOUT, OUT_AUDIT, OUT_SCRIPT):
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
